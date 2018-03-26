from sympl import DataArray
from sympl import get_numpy_array
from sympl import combine_dimensions
from ..._core import ClimtImplicit, get_constant
import numpy as np
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
            maximum_momentum_drag_coefficient=0.002):
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

        """

        self._cyclone = simulate_cyclone
        self._lsc = large_scale_condensation
        self._pbl = boundary_layer
        self._surface_flux = surface_fluxes
        self._use_ext_ts = use_external_surface_temperature
        self._use_ext_qsurf = use_external_surface_specific_humidity

        phys.init_simple_physics(self._cyclone, self._lsc, self._pbl,
                                 self._surface_flux, self._use_ext_ts,
                                 self._use_ext_qsurf)

        self._g = get_constant('gravitational_acceleration', 'm/s^2')

        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')

        self._Rair = get_constant('gas_constant_of_dry_air', 'J/kg/degK')

        self._Rcond = get_constant('gas_constant_of_vapor_phase', 'J/kg/degK')

        self._radius = get_constant('planetary_radius', 'm')

        self._Omega = get_constant('planetary_rotation_rate', 's^-1')

        self._Lv = get_constant('latent_heat_of_condensation', 'J/kg')

        self._Ct = drag_coefficient_heat_fluxes

        self._pbl_top = top_of_boundary_layer

        self._delta_pbl = boundary_layer_influence_height

        self._Cd0 = base_momentum_drag_coefficient

        self._Cd1 = wind_dependent_momentum_drag_coefficient

        self._Cm = maximum_momentum_drag_coefficient

        self._rho_condensible = get_constant('density_of_liquid_water', 'kg/m^3')

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

        if 'latitude' not in state.keys():
            raise IndexError(
                'Simple Physics: State must contain a quantity called latitude')
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

        latent_heat_flux[latent_heat_flux < 0] = 0

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
