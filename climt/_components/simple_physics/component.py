from sympl import Implicit
from sympl import DataArray
from sympl import replace_none_with_default
from sympl import get_numpy_array
from sympl import combine_dimensions
from  . import _simple_physics as phys


class SimplePhysics(Implicit):
    """
    Interface to the simple physics package

    Reed and Jablonowski 2012:
    title = {Idealized tropical cyclone simulations of intermediate complexity: a test case for {AGCMs}}
    journal = {Journal of Advances in Modeling Earth Systems}

    """

    inputs = (
        'air_temperature', 'air_pressure', 'air_pressure_on_interface_levels',
        'surface_pressure', 'surface_temperature', 'specific_humidity',
        'northward_wind', 'eastward_wind', 'latitude'
    )

    diagnostic_outputs = ()

    tendency_outputs = (
        'air_temperature', 'air_pressure', 'air_pressure_on_interface_levels',
        'surface_pressure', 'surface_temperature', 'specific_humidity',
        'northward_wind', 'eastward_wind', 'latitude'
    )

    def __init__(
            self,
            simulate_cyclone=False,
            large_scale_condensation=True,
            boundary_layer=True,
            surface_fluxes=True,
            use_external_surface_temperature=True,
            gravitational_acceleration=None,
            gas_constant_dry_air=None,
            specific_heat_dry_air=None,
            latent_heat_condensation=None,
            gas_constant_condensible=None,
            planetary_radius=None,
            planetary_rotation_rate=None):
        """

        Args:
            simulate_cyclone (bool): Option indicating whether the package must
                simulate a tropical cyclone. This was the original test case this
                physics package was used for.
                Default value is False.

            large_scale_condensation (bool): Option indicating whether the package
                must add moisture and heating tendencies due to large scale condensation.
                Default value is True.

            boundary_layer (bool): Option indicating whether the package must simulate
                a simple boundary layer. **It is recommended that this option remain True
                unless another boundary layer component is being used**.
                Default value is True.

            surface_fluxes (bool): Option indicating whether the package must calculate
                surface fluxes. **It is recommended that this option remain True unless the
                fluxes are being calculated by another component**.
                Default value is True.

            use_external_surface_temperature (bool): Option indicating whether the package
                must use surface temperature available in the model state.
                If False, an internally generated surface temperature is used.
                Default value is True.

            gravitational_acceleration (float): value of acceleration due to gravity.
                Default value from climt.default_constants is used if None.

            gas_constant_dry_air (float): The gas constant of dry air.
                Default value from climt.default_constants is used if None.

            specific_heat_dry_air (float): The specific heat of dry air.
                Default value from climt.default_constants is used if None.

            latent_heat_condensation (float): The latent heat of condensation of the condensible.
                Default value from climt.default_constants (for water vapor) is used if None.

            gas_constant_condensible (float): The gas constant of the condensible substance.
                Default value from climt.default_constants (for water vapor) is used if None.

            planetary_radius (float): The radius of the planet.
                Default value from climt.default_constants (for Earth) is used if None.

            planetary_rotation_rate (float): The rotation rate of the planet.
                Default value from climt.default_constants (for Earth) is used if None.


        """

        self._cyclone = simulate_cyclone
        self._lsc = large_scale_condensation
        self._pbl = boundary_layer
        self._surface_flux = surface_fluxes
        self._use_ext_ts = use_external_surface_temperature

        phys.init_simple_physics(self._cyclone, self._lsc, self._pbl,
                                 self._surface_flux, self._use_ext_ts)

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

        phys.set_physical_constants(self._g, self._Cpd, self._Rair, self._Lv,
                                    self._Rcond, self._radius, self._Omega)

    def __call__(self, state, timestep):
        '''
        Returns the tendencies for a simplified moist physics simulation.

        Args:
            state (dict) : The model state dictionary

            timestep (float) : The model timestep

        Returns:
            state (dict) : The updated model state.

            diagnostics (dict) : empty for this component (No diagnostics)
        '''

        U = get_numpy_array(state['eastward_wind'].to_units('m/s'), ['x', 'y', 'z'])
        V = get_numpy_array(state['northward_wind'].to_units('m/s'), ['x', 'y', 'z'])
        P = get_numpy_array(state['air_pressure'].to_units('Pa'), ['x', 'y', 'z'])
        Pint = get_numpy_array(state['air_pressure_on_interface_levels'].to_units('Pa'), ['x', 'y', 'z'])
        T = get_numpy_array(state['air_temperature'].to_units('degK'), ['x', 'y', 'z'])
        q = get_numpy_array(state['specific_humidity'].to_units('kg/kg'), ['x', 'y', 'z'])

        Ts = get_numpy_array(state['surface_temperature'].to_units('degK'), ['x', 'y'])
        Ps = get_numpy_array(state['surface_pressure'].to_units('Pa'), ['x', 'y'])
        lats = state['latitude'].values

        dims_mid = combine_dimensions([
            state['surface_temperature'],
            state['air_temperature'],
            state['eastward_wind'],
            state['northward_wind']],
            ['x', 'y', 'z'])

        dims_surf = combine_dimensions([
            state['surface_pressure'],
            state['surface_temperature']],
            ['x', 'y'])

        (t_out, u_out, v_out,
         q_out, precip_out) = phys.get_new_state(
             U, V, T, P,
             Pint, q, Ps, Ts,
             lats, timestep.seconds)

        new_state = {
            'eastward_wind': DataArray(
                u_out, dims=dims_mid, attrs=state['eastward_wind'].attrs
            ).squeeze(),
            'northward_wind': DataArray(
                v_out, dims=dims_mid, attrs=state['northward_wind'].attrs
            ).squeeze(),
            'air_temperature': DataArray(
                t_out, dims=dims_mid, attrs=state['air_temperature'].attrs
            ).squeeze(),
            'specific_humidity': DataArray(
                q_out, dims=dims_mid, attrs=state['specific_humidity'].attrs
            ).squeeze(),
            'precipitation': DataArray(
                precip_out, dims=dims_surf, attrs={'units': 'm s^-1'}
            ).squeeze(),
        }

        state.update(new_state)

        return {}, state
