from sympl import DataArray
from sympl import replace_none_with_default
from sympl import get_numpy_array
from sympl import combine_dimensions
from climt import ClimtImplicit, ClimtImplicitPrognostic
import numpy as np
import copy
try:
    from . import _simple_physics as phys
except ImportError:
    print('Import failed. Simple Physics will not be available!')


class SimplePhysics(ClimtImplicit):
    """
    Interface to the simple physics package.

    Reed and Jablonowski 2012:
    title = {Idealized tropical cyclone simulations of intermediate complexity: a test case for {AGCMs}}
    journal = {Journal of Advances in Modeling Earth Systems}

    """

    _climt_inputs = {
        'air_temperature': 'degK',
        'air_pressure': 'mbar',
        'air_pressure_on_interface_levels': 'mbar',
        'surface_air_pressure': 'mbar',
        'surface_temperature': 'degK',
        'specific_humidity': 'g/g',
        'northward_wind': 'm s^-1',
        'eastward_wind': 'm s^-1',
        'latitude': 'degrees_north',
        'surface_specific_humidity': 'g/g',
    }

    _climt_diagnostics = {
        'stratiform_precipitation_rate': 'm s^-1',
        'surface_upward_latent_heat_flux': 'W m^-2',
        'surface_upward_sensible_heat_flux': 'W m^-2',
    }

    _climt_outputs = {
        'air_temperature': 'degK',
        # 'air_pressure': 'mbar',
        # 'air_pressure_on_interface_levels': 'mbar',
        # 'surface_air_pressure': 'mbar',
        # 'surface_temperature': 'degK',
        'specific_humidity': 'g/g',
        'northward_wind': 'm s^-1',
        'eastward_wind': 'm s^-1',
    }

    def __init__(
            self,
            simulate_cyclone=False,
            large_scale_condensation=True,
            boundary_layer=True,
            surface_fluxes=True,
            use_external_surface_temperature=True,
            use_external_surface_specific_humidity=False,
            top_of_boundary_layer=85000.0,
            boundary_layer_influence_height=20000.0,
            drag_coefficient_heat_fluxes=0.0011,
            base_momentum_drag_coefficient=0.0007,
            wind_dependent_momentum_drag_coefficient=0.000065,
            maximum_momentum_drag_coefficient=0.002,
            gravitational_acceleration=None,
            gas_constant_dry_air=None,
            specific_heat_dry_air=None,
            latent_heat_condensation=None,
            gas_constant_condensible=None,
            planetary_radius=None,
            planetary_rotation_rate=None,
            density_condensed_phase=1000.):
        """

        Args:

            simulate_cyclone (bool):
                Option indicating whether the package must
                simulate a tropical cyclone. This was the original test case this
                physics package was used for.
                Default value is False.

            large_scale_condensation (bool):
                Option indicating whether the package
                must add moisture and heating tendencies due to large scale condensation.
                Default value is True.

            boundary_layer (bool):
                Option indicating whether the package must simulate
                a simple boundary layer. **It is recommended that this option remain True
                unless another boundary layer component is being used**.
                Default value is True.

            surface_fluxes (bool):
                Option indicating whether the package must calculate
                surface fluxes. **It is recommended that this option remain True unless the
                fluxes are being calculated by another component**.
                Default value is True.

            use_external_surface_temperature (bool):
                Option indicating whether the package
                must use surface temperature available in the model state.
                If False, an internally generated surface temperature is used.
                Default value is True.

            top_of_boundary_layer (float):
                The nominal top of the boundary layer in :math:`Pa`.

            boundary_layer_influence_height (float):
                The decay of the influence of the boundary layer above
                :code:`top_of_boundary_layer` in :math:`Pa`. The influence
                reduces to :math:`1/e` times the boundary layer value at
                a pressure given by :code:`top_of_boundary_layer+boundary_layer_influence_height`.

            drag_coefficient_heat_fluxes (float):
                The wind speed independent drag coefficient for latent and sensible
                heat fluxes.

            base_momentum_drag_coefficient (float):
                The minimum drag coefficient for winds.

            wind_dependent_momentum_drag_coefficient (float):
                The part of the momentum drag coefficient that depends on the surface wind
                speed. The total drag coefficient is given by
                :code:`base_momentum_drag_coefficient + wind_dependent_momentum_drag_coefficient*u_base`,
                where :code:`u_base` is the surface wind speed.

            maximum_momentum_drag_coefficient (float):
                This drag coefficient is used for surface wind speeds exceeding :math:`20 m/s`.

            gravitational_acceleration (float):
                value of acceleration due to gravity in :math:`m s^{-1}`.
                Default value from climt.default_constants is used if None.

            gas_constant_dry_air (float):
                The gas constant of dry air in :math:`J K^{-1} kg^{-1}`.
                Default value from climt.default_constants is used if None.

            specific_heat_dry_air (float):
                The specific heat of dry air in :math:`J K^{-1} kg^{-1}`.
                Default value from climt.default_constants is used if None.

            latent_heat_condensation (float):
                The latent heat of condensation of the condensible in :math:`J kg^{-1}`.
                Default value from climt.default_constants (for water vapor) is used if None.

            gas_constant_condensible (float):
                The gas constant of the condensible substance in :math:`J {K^-1} kg^{-1}`.
                Default value from climt.default_constants (for water vapor) is used if None.

            planetary_radius (float):
                The radius of the planet in :math:`m`.
                Default value from climt.default_constants (for Earth) is used if None.

            planetary_rotation_rate (float):
                The rotation rate of the planet in :math:`s^{-1}`.
                Default value from climt.default_constants (for Earth) is used if None.

        """

        self._cyclone = simulate_cyclone
        self._lsc = large_scale_condensation
        self._pbl = boundary_layer
        self._surface_flux = surface_fluxes
        self._use_ext_ts = use_external_surface_temperature
        self._use_ext_qsurf = use_external_surface_specific_humidity
        self._return_tend = False

        phys.init_simple_physics(self._cyclone, self._lsc, self._pbl,
                                 self._surface_flux, self._use_ext_ts,
                                 self._use_ext_qsurf)

        self._g = replace_none_with_default(
            'gravitational_acceleration', gravitational_acceleration)

        self._Cpd = replace_none_with_default(
            'heat_capacity_of_dry_air_at_constant_pressure',
            specific_heat_dry_air)

        self._Rair = replace_none_with_default(
            'gas_constant_of_dry_air',
            gas_constant_dry_air)

        self._Rcond = replace_none_with_default(
            'gas_constant_of_water_vapor',
            gas_constant_condensible)

        self._radius = replace_none_with_default(
            'planetary_radius',
            planetary_radius)

        self._Omega = replace_none_with_default(
            'planetary_rotation_rate',
            planetary_rotation_rate)

        self._Lv = replace_none_with_default(
            'latent_heat_of_vaporization_of_water',
            latent_heat_condensation)

        self._Ct = drag_coefficient_heat_fluxes

        self._pbl_top = top_of_boundary_layer

        self._delta_pbl = boundary_layer_influence_height

        self._Cd0 = base_momentum_drag_coefficient

        self._Cd1 = wind_dependent_momentum_drag_coefficient

        self._Cm = maximum_momentum_drag_coefficient

        self._rho_condensible = replace_none_with_default(
            'density_of_liquid_water',
            density_condensed_phase)

        phys.set_physical_constants(self._g, self._Cpd, self._Rair, self._Lv,
                                    self._Rcond, self._radius, self._Omega,
                                    self._rho_condensible, self._pbl_top,
                                    self._delta_pbl, self._Ct, self._Cd0,
                                    self._Cd1, self._Cm)

    def __call__(self, state, timestep):
        '''
        Calculate surface and boundary layer tendencies.

        Args:
            state (dict):
                The model state dictionary

            timestep (timedelta):
                The model timestep

        Returns:
            state (dict), diagnostics(dict) :

            * The updated model state.
            * diagnostics for Simple Physics
        '''

        U = get_numpy_array(state['eastward_wind'].to_units('m/s'), ['x', 'y', 'z'])
        V = get_numpy_array(state['northward_wind'].to_units('m/s'), ['x', 'y', 'z'])
        P = get_numpy_array(state['air_pressure'].to_units('Pa'), ['x', 'y', 'z'])
        Pint = get_numpy_array(state['air_pressure_on_interface_levels'].to_units('Pa'), ['x', 'y', 'z'])
        T = get_numpy_array(state['air_temperature'].to_units('degK'), ['x', 'y', 'z'])
        q = get_numpy_array(state['specific_humidity'].to_units('g/g'), ['x', 'y', 'z'])
        q_surface = get_numpy_array(state['surface_specific_humidity'].to_units('g/g'), ['x', 'y'])

        Ts = get_numpy_array(state['surface_temperature'].to_units('degK'), ['x', 'y'])
        Ps = get_numpy_array(state['surface_air_pressure'].to_units('Pa'), ['x', 'y'])
        lats = get_numpy_array(state['latitude'].to_units('degrees_north'), ['x', 'y'])
        lats = np.asfortranarray(lats, dtype=np.double)

        if lats.shape[0] == 1:  # 1-d array only
            num_longitudes = Ts.shape[0]
            lat_list = [lats[0, :] for i in range(num_longitudes)]
            lats = np.asfortranarray(np.stack(lat_list, axis=0))

        dims_mid = combine_dimensions([
            state['surface_temperature'],
            state['air_temperature'],
            state['eastward_wind'],
            state['northward_wind']],
            ['x', 'y', 'z'])

        (t_out, u_out, v_out,
         q_out, precip_out,
         sensible_heat_flux, latent_heat_flux) = phys.get_new_state(
             U, V, T, P,
             Pint, q, Ps, Ts, q_surface,
             lats, timestep.total_seconds())

        if self._return_tend:
            t_out = (t_out - T)/timestep.total_seconds()
            u_out = (u_out - U)/timestep.total_seconds()
            v_out = (v_out - V)/timestep.total_seconds()
            q_out = (q_out - q)/timestep.total_seconds()
            return (t_out, u_out, v_out, q_out, precip_out,
                    sensible_heat_flux, latent_heat_flux)

        new_state = {
            'eastward_wind': DataArray(
                u_out, dims=dims_mid, attrs=state['eastward_wind'].attrs
            ),
            'northward_wind': DataArray(
                v_out, dims=dims_mid, attrs=state['northward_wind'].attrs
            ),
            'air_temperature': DataArray(
                t_out, dims=dims_mid, attrs=state['air_temperature'].attrs
            ),
            'specific_humidity': DataArray(
                q_out, dims=dims_mid, attrs=state['specific_humidity'].attrs
            ),
        }

        diagnostics = self.create_state_dict_for('_climt_diagnostics', state)

        diagnostics['stratiform_precipitation_rate'].values[:] = precip_out
        diagnostics['surface_upward_sensible_heat_flux'].values[:] = sensible_heat_flux
        diagnostics['surface_upward_latent_heat_flux'].values[:] = latent_heat_flux

        return diagnostics, new_state

    def get_prognostic_version(self, model_time_step):
        """
        Get a component which returns tendencies instead of new model state.

        This is mainly to allow SimplePhysics to be used with the spectral
        dynamical core.

        Args:
            model_time_step (timedelta):
                The model time step used.

        Returns:
            SimplePhysicsWithTendencies(ClimtImplicitPrognostic):
                Version of SimplePhysics which returns tendencies.

        """

        new_self = copy.copy(self)
        new_self._return_tend = True
        return SimplePhysicsWithTendencies(new_self, model_time_step)


class SimplePhysicsWithTendencies(ClimtImplicitPrognostic):
    """
    Version of simple physics that returns tendencies
    """

    def __init__(self, simple_physics, timestep):

        self._implicit_simple_physics = simple_physics
        self.current_time_step = timestep

        self._climt_inputs = simple_physics._climt_inputs
        self._climt_diagnostics = simple_physics._climt_diagnostics
        self._climt_tendencies = simple_physics._climt_outputs
        self._climt_tendencies['air_temperature'] = 'degK/s'
        self._climt_tendencies['eastward_wind'] = 'm s^-2'
        self._climt_tendencies['northward_wind'] = 'm s^-2'
        self._climt_tendencies['specific_humidity'] = 'g g^-1 s^-1'

    def __call__(self, state):
        (t_out, u_out, v_out, q_out, precip_out,
         sensible_heat_flux, latent_heat_flux) =\
            self._implicit_simple_physics(state,
                                          self.current_time_step)

        tendencies = self.create_state_dict_for('_climt_tendencies', state)
        tendencies['air_temperature'].values = t_out
        tendencies['eastward_wind'].values = u_out
        tendencies['northward_wind'].values = v_out
        tendencies['specific_humidity'].values = q_out

        diagnostics = self.create_state_dict_for('_climt_diagnostics', state)

        diagnostics['stratiform_precipitation_rate'].values = precip_out
        diagnostics['surface_upward_sensible_heat_flux'].values = sensible_heat_flux
        diagnostics['surface_upward_latent_heat_flux'].values = latent_heat_flux

        return tendencies, diagnostics
