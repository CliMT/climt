from sympl import Stepper, get_constant, initialize_numpy_arrays_with_properties
import numpy as np
# from scipy.interpolate import CubicSpline
from scipy import sparse
from scipy.sparse.linalg import spsolve


class IceSheet(Stepper):
    """
    1-d snow-ice energy balance model.
    """

    input_properties = {
        'downwelling_longwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'downwelling_shortwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'upwelling_longwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'upwelling_shortwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'surface_upward_latent_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_upward_sensible_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'land_ice_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'sea_ice_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'surface_snow_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'area_type': {
            'dims': ['*'],
            'units': 'dimensionless',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'snow_and_ice_temperature': {
            'dims': ['ice_interface_levels', '*'],
            'units': 'degK',
        },
        'sea_surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'soil_surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'height_on_ice_interface_levels': {
            'dims': ['ice_interface_levels', '*'],
            'units': 'm',
        },
    }

    output_properties = {
        'land_ice_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'sea_ice_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'surface_snow_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'snow_and_ice_temperature': {
            'dims': ['ice_interface_levels', '*'],
            'units': 'degK',
        },
        'height_on_ice_interface_levels': {
            'dims': ['ice_interface_levels', '*'],
            'units': 'm',
        },
        'sea_surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
    }

    diagnostic_properties = {
        'upward_heat_flux_at_ground_level_in_soil': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'heat_flux_into_sea_water_due_to_sea_ice': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_downward_heat_flux_in_sea_ice': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_albedo_for_direct_shortwave': {
            'dims': ['*'],
            'units': 'dimensionless'
        },
        'surface_albedo_for_diffuse_shortwave': {
            'dims': ['*'],
            'units': 'dimensionless'
        },
    }

    def __init__(self, maximum_snow_ice_height=10, **kwargs):
        """
        Args:
            maximum_snow_ice_height (float):
                The maximum combined height of snow and ice handled by the model in :math:`m`.
            levels (int):
                The number of levels on which temperature must be output.
        """
        self._max_height = maximum_snow_ice_height
        self._update_constants()
        self._epsilon = 1e-5
        super(IceSheet, self).__init__(**kwargs)

    def _update_constants(self):
        self._Kice = get_constant('thermal_conductivity_of_solid_phase_as_ice', 'W/m/degK')
        self._Ksnow = get_constant('thermal_conductivity_of_solid_phase_as_snow', 'W/m/degK')
        self._rho_ice = get_constant('density_of_solid_phase_as_ice', 'kg/m^3')
        self._C_ice = get_constant('heat_capacity_of_solid_phase_as_ice', 'J/kg/degK')
        self._rho_snow = get_constant('density_of_solid_phase_as_snow', 'kg/m^3')
        self._C_snow = get_constant('heat_capacity_of_solid_phase_as_snow', 'J/kg/degK')
        self._Lf = get_constant('latent_heat_of_fusion', 'J/kg')
        self._melting_temperature = get_constant('freezing_temperature_of_liquid_phase', 'degK')

    def array_call(self, raw_state, timestep):
        self._update_constants()

        num_cols = raw_state['area_type'].shape[0]

        net_heat_flux = (
            raw_state['downwelling_shortwave_flux_in_air'][:, 0] +
            raw_state['downwelling_longwave_flux_in_air'][:, 0] -
            raw_state['upwelling_shortwave_flux_in_air'][:, 0] -
            raw_state['upwelling_longwave_flux_in_air'][:, 0] -
            raw_state['surface_upward_sensible_heat_flux'] -
            raw_state['surface_upward_latent_heat_flux']
        )

        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, raw_state, self.input_properties
        )

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )

        # Copy input values
        outputs['surface_temperature'][:] = raw_state['surface_temperature']
        outputs['sea_surface_temperature'][:] = raw_state['sea_surface_temperature']
        outputs['land_ice_thickness'][:] = raw_state['land_ice_thickness']
        outputs['sea_ice_thickness'][:] = raw_state['sea_ice_thickness']
        outputs['surface_snow_thickness'][:] = raw_state['surface_snow_thickness']
        outputs['snow_and_ice_temperature'][:] = raw_state['snow_and_ice_temperature']

        for col in range(num_cols):
            area_type = raw_state['area_type'][col].astype(str)
            total_height = 0.
            surface_temperature = raw_state['snow_and_ice_temperature'][:, col][-1]
            soil_surface_temperature = None

            if area_type == 'land_ice':
                total_height = raw_state['land_ice_thickness'][col] \
                    + raw_state['surface_snow_thickness'][col]
                soil_surface_temperature = raw_state['soil_surface_temperature'][col]
            elif area_type == 'sea_ice':
                if raw_state['sea_ice_thickness'][col] == 0:
                    # No sea ice, so skip calculation
                    continue
                total_height = raw_state['sea_ice_thickness'][col] \
                    + raw_state['surface_snow_thickness'][col]
            elif area_type == 'land':
                total_height = raw_state['surface_snow_thickness'][col]
                soil_surface_temperature = raw_state['soil_surface_temperature'][col]
            if total_height > self._max_height:
                raise ValueError("Total height exceeds maximum value of {} m.".format(self._max_height))

            if total_height < self._epsilon:  # Some epsilon
                continue

            snow_height_fraction = raw_state['surface_snow_thickness'][col] / total_height

            temp_profile = raw_state['snow_and_ice_temperature'][:, col]
            num_layers = temp_profile.shape[0]
            dz = float(total_height / num_layers)

            snow_level = int((1 - snow_height_fraction)*num_layers) - 1
            levels = np.arange(num_layers - 1)

            # Create vertically varying profiles
            rho_snow_ice = self._rho_ice*np.ones(num_layers - 1)
            rho_snow_ice[levels > snow_level] = self._rho_snow

            heat_capacity_snow_ice = self._C_ice*np.ones(num_layers - 1)
            heat_capacity_snow_ice[levels > snow_level] = self._C_snow

            kappa_snow_ice = self._Kice*np.ones(num_layers - 1)
            kappa_snow_ice[levels > snow_level] = self._Ksnow

            check_melting = True
            if surface_temperature < self._melting_temperature - self._epsilon:
                check_melting = False

            new_temperature = self.calculate_new_ice_temperature(
                rho_snow_ice, heat_capacity_snow_ice,
                kappa_snow_ice, temp_profile,
                timestep.total_seconds(), dz,
                num_layers,
                surface_temperature,
                net_heat_flux[col],
                soil_surface_temperature)

            # Cool down from melting temperature if
            # heat flux through ice is greater than
            # net forcing heat flux at surface
            heat_flux_through_ice = ((new_temperature[-1] - new_temperature[-2]) *
                                     (kappa_snow_ice[-1] + kappa_snow_ice[-2])*0.5/dz)

            if temp_profile[-1] > self._melting_temperature - self._epsilon:
                if heat_flux_through_ice > net_heat_flux[col]:
                    surface_temperature = temp_profile[-1] - 10*self._epsilon
                    temp_profile[-1] = temp_profile[-1] - 10*self._epsilon

                    new_temperature = self.calculate_new_ice_temperature(
                        rho_snow_ice, heat_capacity_snow_ice,
                        kappa_snow_ice, temp_profile,
                        timestep.total_seconds(), dz,
                        num_layers,
                        surface_temperature,
                        net_heat_flux[col],
                        soil_surface_temperature)

                    check_melting = False

            # Energy balance for lower surface of snow/ice
            if area_type == 'sea_ice':
                # TODO Add ocean heat flux parameterization
                # At sea surface
                heat_flux_to_sea_water = (new_temperature[1] - new_temperature[0])*(
                    kappa_snow_ice[0] + kappa_snow_ice[1])*0.5/dz
                heat_flux_to_sea_water = round(heat_flux_to_sea_water, 6)

                # If heat_flux_to_sea_water is positive, flux of heat into water
                # an impossible situation which means ice is above freezing point.
                assert heat_flux_to_sea_water <= 0

                height_of_growing_ice = -(heat_flux_to_sea_water*timestep.total_seconds() /
                                          (rho_snow_ice[0]*self._Lf))

                outputs['sea_ice_thickness'][col] += height_of_growing_ice
                diagnostics['heat_flux_into_sea_water_due_to_sea_ice'][col]\
                    = heat_flux_to_sea_water

            elif area_type in ['land_ice', 'land']:
                # At land surface
                heat_flux_to_land = (new_temperature[0] - new_temperature[1]) * kappa_snow_ice[0] / dz

                diagnostics['upward_heat_flux_at_ground_level_in_soil'][col] \
                    = heat_flux_to_land

                height_of_growing_ice = 0

            else:
                continue

            # Energy balance at atmosphere surface
            heat_flux_through_ice = ((new_temperature[-1] - new_temperature[-2]) *
                                     (kappa_snow_ice[-1] + kappa_snow_ice[-2])*0.5/dz)
            diagnostics['surface_downward_heat_flux_in_sea_ice'][col] = heat_flux_through_ice

            height_of_melting_ice = 0
            # Surface is melting
            if check_melting:
                energy_to_melt_ice = (net_heat_flux[col] - heat_flux_through_ice)*timestep.total_seconds()
                energy_to_melt_ice = round(energy_to_melt_ice, 6)
                if energy_to_melt_ice < 0:
                    print(net_heat_flux[col], heat_flux_through_ice)
                    assert False

                height_of_melting_ice = (energy_to_melt_ice / (rho_snow_ice[-1]*self._Lf))

                if height_of_melting_ice > raw_state['surface_snow_thickness'][col]:

                    outputs['sea_ice_thickness'][col] -= (
                        height_of_melting_ice - raw_state['surface_snow_thickness'][col])
                    outputs['surface_snow_thickness'][col] = 0

                else:
                    outputs['surface_snow_thickness'][col] -= height_of_melting_ice

            total_height += (height_of_growing_ice + height_of_melting_ice)

            outputs['snow_and_ice_temperature'][:, col] = new_temperature

            outputs['surface_temperature'][col] = new_temperature[-1]
            outputs['height_on_ice_interface_levels'][:, col] = np.linspace(
                0,
                total_height,
                outputs['height_on_ice_interface_levels'].shape[0],
                endpoint=True
            )

            if outputs['surface_snow_thickness'][col] > 0:
                diagnostics['surface_albedo_for_direct_shortwave'][col] = 0.8
                diagnostics['surface_albedo_for_diffuse_shortwave'][col] = 0.8
            elif area_type == 'sea_ice' and outputs['sea_ice_thickness'][col] > 0:
                diagnostics['surface_albedo_for_direct_shortwave'][col] = 0.5
                diagnostics['surface_albedo_for_diffuse_shortwave'][col] = 0.5
            elif area_type == 'sea_ice' and outputs['sea_ice_thickness'][col] > 0:
                diagnostics['surface_albedo_for_direct_shortwave'][col] = 0.5
                diagnostics['surface_albedo_for_diffuse_shortwave'][col] = 0.5

            if height_of_melting_ice > 0:
                diagnostics['surface_albedo_for_direct_shortwave'][col] = 0.2
                diagnostics['surface_albedo_for_diffuse_shortwave'][col] = 0.2

        return diagnostics, outputs

    def calculate_new_ice_temperature(self, rho, specific_heat, kappa,
                                      temp_profile, dt, dz,
                                      num_layers, surf_temperature, net_flux,
                                      soil_temperature=None):

        r = np.zeros(num_layers)
        a_sub = np.zeros(num_layers)
        a_sup = np.zeros(num_layers)

        K_interface = 0.5*(kappa[:-1] + kappa[1:])
        K_mid = kappa

        heat_capacity = rho * specific_heat
        heat_capacity_int = 0.5*(heat_capacity[:-1] + heat_capacity[1:])

        mu_inv_int = dt / (heat_capacity_int * 2 * dz * dz)

        r[1:-1] = K_interface*mu_inv_int

        dp = (1 + 2*r)
        dm = (1 - 2*r)

        a_sub[:-2] = -mu_inv_int*K_mid[:-1]
        a_sup[2:] = -mu_inv_int*K_mid[1:]

        mat_lhs = sparse.spdiags([a_sub, dp, a_sup], [-1, 0, 1], num_layers, num_layers, format='csc')

        mat_rhs = sparse.spdiags([-a_sub, dm, -a_sup], [-1, 0, 1], num_layers, num_layers, format='csc')

        rhs = mat_rhs * temp_profile

        # Set flux condition if temperature is below melting point,
        # and dirichlet condition above melting point
        if surf_temperature < self._melting_temperature - self._epsilon:
            mat_lhs[-1, -1] = -1
            mat_lhs[-1, -2] = 1
            rhs[-1] = -net_flux*dz/K_mid[-1]
        else:
            mat_lhs[-1, -1] = 1
            mat_lhs[-1, -2] = 0
            rhs[-1] = self._melting_temperature

        mat_lhs[0, 0] = 1
        mat_lhs[0, 1] = 0
        if soil_temperature is None:
            rhs[0] = self._melting_temperature
        else:
            rhs[0] = soil_temperature

        return spsolve(mat_lhs, rhs)
