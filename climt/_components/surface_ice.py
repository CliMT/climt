from climt import ClimtImplicit
from sympl import replace_none_with_default
import numpy as np
from scipy.interpolate import CubicSpline
from scipy import sparse
from scipy.sparse.linalg import spsolve


class IceSheet(ClimtImplicit):
    """
    1-d snow-ice energy balance model.
    """

    _climt_inputs = {
        'downwelling_longwave_flux_in_air': 'W m^-2',
        'downwelling_shortwave_flux_in_air': 'W m^-2',
        'upwelling_longwave_flux_in_air': 'W m^-2',
        'upwelling_shortwave_flux_in_air': 'W m^-2',
        'surface_upward_latent_heat_flux': 'W m^-2',
        'surface_upward_sensible_heat_flux': 'W m^-2',
        'land_ice_thickness': 'm',
        'sea_ice_thickness': 'm',
        'surface_snow_thickness': 'm',
        'area_type': 'dimensionless',
        'surface_temperature': 'degK',
        'snow_ice_temperature_poly': 'dimensionless',
        'sea_surface_temperature': 'degK',
        'soil_surface_temperature': 'degK'
    }

    _climt_outputs = {
        'land_ice_thickness': 'm',
        'sea_ice_thickness': 'm',
        'surface_snow_thickness': 'm',
        'surface_temperature': 'degK',
        'snow_ice_temperature': 'degK',
        'sea_surface_temperature': 'degK',
    }

    _climt_diagnostics = {
        'upward_heat_flux_at_ground_level_in_soil': 'W m^-2',
        'heat_flux_into_sea_water_due_to_sea_ice': 'W m^-2',
        'snow_ice_temperature_poly': 'dimensionless',
    }

    quantity_descriptions = {
        'snow_ice_temperature': {
            'dims': ['x', 'y', 'ice_vertical_levels'],
            'units': 'degK',
            'default_value': 260.,
        },
    }

    extra_dimensions = {'ice_vertical_levels': np.arange(10)}

    def __init__(self,
                 vertical_resolution=0.1,
                 maximum_snow_ice_height=10,
                 number_vertical_levels=30,
                 thermal_conductivity_of_ice=2.22,
                 thermal_conductivity_of_snow=0.2,
                 density_of_ice=916.7,
                 density_of_snow=100.,
                 specific_heat_of_ice=2108.,
                 specific_heat_of_snow=2108.,
                 latent_heat_of_fusion=333550.,
                 melting_point_of_ice=273.):
        """

        Args:
            vertical_resolution(float):
                The vertical resolution of the model in :math:`m`.

            maximum_snow_ice_height(float):
                The maximum combined height of snow and ice handled by the model in :math:`m`.

            number_vertical_levels(float):
                The number of levels on which temperature must be output.

            thermal_conductivity_of_ice (float, optional):
                The thermal conductivity of ice in :math:`W m^{-1} K^{-1}`.

            thermal_conductivity_of_snow (float, optional):
                The thermal conductivity of ice in :math:`W m^{-1} K^{-1}`.

            density_of_ice (float, optional):
                The density of ice in :math:`kg m^{-3}`.

            density_of_snow (float, optional):
                The density of snow in :math:`kg m^{-3}`.

            specific_heat_of_ice (float, optional):
                The specific heat of ice in :math:`J kg^{-1}`.

            specific_heat_of_snow (float, optional):
                The specific heat of snow in :math:`J kg^{-1}`.

            latent_heat_of_fusion (float, optional):
                The latent heat of fusion of ice in :math:`J kg^{-1}`

            melting_point_of_ice (float, optional):
                The melting point of ice in :math:`degK`.
        """

        self._dz = vertical_resolution

        self._max_height = maximum_snow_ice_height

        self._output_levels = number_vertical_levels
        self.extra_dimensions['ice_vertical_levels'] = np.arange(self._output_levels)

        self.quantity_descriptions['default_value'] = CubicSpline(
            np.linspace(0, maximum_snow_ice_height, 50), 271.*np.ones(50))

        self._Kice = replace_none_with_default(
            'thermal_conductivity_of_ice', thermal_conductivity_of_ice)

        self._Ksnow = replace_none_with_default(
            'thermal_conductivity_of_snow', thermal_conductivity_of_snow)

        self._rho_ice = replace_none_with_default(
            'density_of_ice', density_of_ice)

        self._C_ice = replace_none_with_default(
            'specific_heat_of_ice', specific_heat_of_ice)

        self._rho_snow = replace_none_with_default(
            'density_of_snow', density_of_snow)

        self._C_snow = replace_none_with_default(
            'specific_heat_of_snow', specific_heat_of_snow)

        self._Lf = replace_none_with_default(
            'latent_heat_of_fusion_of_ice', latent_heat_of_fusion)

        self._temp_melt = replace_none_with_default(
            'melting_point_of_ice', melting_point_of_ice)

    def __call__(self, state, timestep):
        """
        Calculate new ice sheet height.

        Args:
            state (dict):
                The state dictionary.

            time_step (timedelta):
                The model timestep.

        Returns:
            new_state(dict), diagnostics(dict):
                * The new state calculated
                * Any diagnostics

        """
        # TODO use land and sea surface temperature for surface temperature calculations
        input_arrays = self.get_numpy_arrays_from_state('_climt_inputs', state)

        num_lons, num_lats = input_arrays['area_type'].shape

        net_heat_flux = (input_arrays['downwelling_shortwave_flux_in_air'][:, :, 0] +
                         input_arrays['downwelling_longwave_flux_in_air'][:, :, 0] -
                         input_arrays['upwelling_shortwave_flux_in_air'][:, :, 0] -
                         input_arrays['upwelling_longwave_flux_in_air'][:, :, 0] -
                         input_arrays['surface_upward_sensible_heat_flux'] -
                         input_arrays['surface_upward_latent_heat_flux'])

        output_dict = self.create_state_dict_for('_climt_outputs', state)
        output_arrays = self.get_numpy_arrays_from_state('_climt_outputs', output_dict)

        diagnostic_dict = self.create_state_dict_for('_climt_diagnostics', state)
        diagnostic_arrays = self.get_numpy_arrays_from_state('_climt_diagnostics', diagnostic_dict)

        # Copy input values
        output_arrays['surface_temperature'][:] = input_arrays['surface_temperature']
        output_arrays['sea_surface_temperature'][:] = input_arrays['sea_surface_temperature']
        output_arrays['land_ice_thickness'][:] = input_arrays['land_ice_thickness']
        output_arrays['sea_ice_thickness'][:] = input_arrays['sea_ice_thickness']
        output_arrays['surface_snow_thickness'][:] = input_arrays['surface_snow_thickness']

        for lon in range(num_lons):
            for lat in range(num_lats):

                    area_type = input_arrays['area_type'][lon, lat].astype(str)
                    total_height = 0.
                    surface_temperature = input_arrays['surface_temperature'][lon, lat]
                    soil_surface_temperature = None

                    if area_type == 'land_ice':
                        total_height = input_arrays['land_ice_thickness'][lon, lat] \
                            + input_arrays['surface_snow_thickness'][lon, lat]
                        soil_surface_temperature = input_arrays['soil_surface_temperature'][lon, lat]
                    elif area_type == 'sea_ice':
                        if input_arrays['sea_ice_thickness'][lon, lat] == 0:
                            # No sea ice, so skip calculation
                            continue
                        total_height = input_arrays['sea_ice_thickness'][lon, lat] \
                            + input_arrays['surface_snow_thickness'][lon, lat]
                    elif area_type == 'land':
                        total_height = input_arrays['surface_snow_thickness'][lon, lat]
                        soil_surface_temperature = input_arrays['soil_surface_temperature'][lon, lat]
                    if total_height > self._max_height:
                        raise ValueError("Total height exceeds maximum value of {} m.".format(self._max_height))

                    if total_height < 1e-8:  # Some epsilon
                        continue

                    # print('total_height', total_height)
                    snow_height_fraction = input_arrays['surface_snow_thickness'][lon, lat] / total_height

                    num_layers = int(total_height / self._dz)

                    # Create temperature profile from cubic spline
                    vertical_coordinates = np.linspace(0, total_height, num_layers)
                    temp_profile = input_arrays['snow_ice_temperature_poly'][lon, lat](vertical_coordinates)

                    snow_level = int((1 - snow_height_fraction)*num_layers)
                    levels = np.arange(num_layers)
                    # print('snow_level', snow_level)
                    # print('num_layers', num_layers)

                    # Create vertically varying profiles
                    rho_snow_ice = self._rho_ice*np.ones(num_layers)
                    rho_snow_ice[levels > snow_level] = self._rho_snow

                    heat_capacity_snow_ice = self._C_ice*np.ones(num_layers)
                    heat_capacity_snow_ice[levels > snow_level] = self._C_snow

                    kappa_snow_ice = self._Kice*np.ones(num_layers)
                    kappa_snow_ice[levels > snow_level] = self._Ksnow

                    check_melting = True
                    if surface_temperature < self._temp_melt:
                        check_melting = False

                    # Calculate new temp_profile based using implicit method
                    new_temp = self.calculate_new_ice_temperature(
                        rho_snow_ice, heat_capacity_snow_ice,
                        kappa_snow_ice, temp_profile,
                        timestep.total_seconds(), num_layers,
                        surface_temperature,
                        net_heat_flux[lon, lat],
                        soil_surface_temperature)

                    # print('net_heat_flux', net_heat_flux[lon, lat])
                    # print('SH', input_arrays['surface_upward_sensible_heat_flux'][lon, lat])
                    # print('LH', input_arrays['surface_upward_latent_heat_flux'][lon, lat])

                    # Energy balance for lower surface of snow/ice
                    if area_type == 'sea_ice':
                        # TODO Add ocean heat flux parameterization
                        # At sea surface
                        heat_flux_to_sea_water = (new_temp[1] - new_temp[0])*kappa_snow_ice[0]/self._dz
                        # print('heat_flux_to_sea_water', heat_flux_to_sea_water)

                        # If heat_flux_to_sea_water is positive, flux of heat into water
                        # an impossible situation which means ice is above freezing point.
                        assert heat_flux_to_sea_water <= 0

                        height_of_growing_ice = -(heat_flux_to_sea_water*timestep.total_seconds() /
                                                  (rho_snow_ice[0]*self._Lf))

                        output_arrays['sea_ice_thickness'][lon, lat] += height_of_growing_ice
                        diagnostic_arrays['heat_flux_into_sea_water_due_to_sea_ice'][lon, lat]\
                            = heat_flux_to_sea_water

                        # print('height_of_growing_ice', height_of_growing_ice)

                    elif area_type in ['land_ice', 'land']:
                        # At land surface
                        heat_flux_to_land = (new_temp[0] - new_temp[1]) * kappa_snow_ice[0] / self._dz
                        # print('heat_flux_to_land', heat_flux_to_land)

                        diagnostic_arrays['upward_heat_flux_at_ground_level_in_soil'][lon, lat] \
                            = heat_flux_to_land

                        height_of_growing_ice = 0

                    else:
                        continue

                    # Energy balance at atmosphere surface
                    heat_flux_to_atmosphere = ((new_temp[-1] - new_temp[-2]) *
                                               (kappa_snow_ice[-1] + kappa_snow_ice[-2])*0.5/self._dz)
                    # print('heat_flux_to_atmosphere', heat_flux_to_atmosphere)

                    height_of_melting_ice = 0
                    # Surface is melting
                    if check_melting:
                        energy_to_melt_ice = (net_heat_flux[lon, lat] + heat_flux_to_atmosphere)

                        height_of_melting_ice = (energy_to_melt_ice*timestep.total_seconds() /
                                                 (rho_snow_ice[-1]*self._Lf))
                        # print('height_of_melting_ice', height_of_melting_ice)

                        if height_of_melting_ice > input_arrays['surface_snow_thickness'][lon, lat]:

                            output_arrays['sea_ice_thickness'][lon, lat] -= (
                                height_of_melting_ice - input_arrays['surface_snow_thickness'][lon, lat])
                            output_arrays['surface_snow_thickness'][lon, lat] = 0

                        else:
                            output_arrays['surface_snow_thickness'][lon, lat] -= height_of_melting_ice

                    total_height += (height_of_growing_ice + height_of_melting_ice)

                    diagnostic_arrays['snow_ice_temperature_poly'][lon, lat] = CubicSpline(
                        np.linspace(0, total_height, num_layers), new_temp)

                    output_arrays['snow_ice_temperature'][lon, lat] = \
                        diagnostic_arrays['snow_ice_temperature_poly'][lon, lat](
                            np.linspace(0, total_height, self._output_levels))

                    output_arrays['surface_temperature'][lon, lat] = new_temp[-1]

        # print()
        # print()
        return output_dict, diagnostic_dict

    def calculate_new_ice_temperature(self, rho, specific_heat, kappa,
                                      temp_profile, dt,
                                      num_layers, surf_temp, net_flux,
                                      soil_temperature=None):

        r = np.zeros(num_layers)
        a_sub = np.zeros(num_layers)
        a_sup = np.zeros(num_layers)

        K_bar = 0.25*(kappa[2:] + kappa[:-2]) + 0.5 * kappa[1:-1]
        K_mid = 0.5*(kappa[1:] + kappa[:-1])

        mu_inv = dt / (rho * specific_heat * 2 * self._dz * self._dz)

        r[1:-1] = K_bar*mu_inv[1:-1]

        dp = (1 + 2*r)
        dm = (1 - 2*r)

        a_sub[:-1] = -mu_inv[1:]*K_mid
        a_sup[1:] = -mu_inv[:-1]*K_mid

        mat_lhs = sparse.spdiags([a_sub, dp, a_sup], [-1, 0, 1], num_layers, num_layers, format='csc')

        mat_rhs = sparse.spdiags([-a_sub, dm, -a_sup], [-1, 0, 1], num_layers, num_layers, format='csc')

        rhs = mat_rhs * temp_profile

        # Set flux condition if temperature is below melting point,
        # and dirichlet condition above melting point
        if surf_temp < self._temp_melt:
            mat_lhs[-1, -1] = -1
            mat_lhs[-1, -2] = 1
            rhs[-1] = -net_flux*self._dz/K_mid[-1]
        else:
            mat_lhs[-1, -1] = 1
            mat_lhs[-1, -2] = 0
            rhs[-1] = self._temp_melt

        mat_lhs[0, 0] = 1
        mat_lhs[0, 1] = 0
        if soil_temperature is None:
            rhs[0] = self._temp_melt
        else:
            rhs[0] = soil_temperature

        return spsolve(mat_lhs, rhs)
