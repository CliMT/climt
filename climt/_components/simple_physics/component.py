from ..._core import ensure_contiguous_state
from sympl import Stepper, get_constant
import logging
try:
    from . import _simple_physics as phys
except ImportError as error:
    logging.warning(
        'Import failed. Simple Physics is likely not compiled and will not be'
        'available.'
    )
    print(error)


class SimplePhysics(Stepper):
    """
    Interface to the simple physics package.

    Reed and Jablonowski 2012:
    title = {Idealized tropical cyclone simulations of intermediate complexity: a test case for {AGCMs}}
    journal = {Journal of Advances in Modeling Earth Systems}

    """

    input_properties = {
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK',
        },
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa',
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'Pa',
        },
        'surface_air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'kg/kg',
        },
        'northward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'eastward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'surface_specific_humidity': {
            'dims': ['*'],
            'units': 'kg/kg',
        },
        'latitude': {
            'dims': ['*'],
            'units': 'degrees_north',
        }
    }

    diagnostic_properties = {
        'stratiform_precipitation_rate': {
            'dims': ['*'],
            'units': 'm s^-1',
        },
        'surface_upward_latent_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_upward_sensible_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
    }

    output_properties = {
        'air_temperature': {'units': 'degK'},
        'specific_humidity': {'units': 'kg/kg'},
        'northward_wind': {'units': 'm s^-1'},
        'eastward_wind': {'units': 'm s^-1'},
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
            **kwargs):
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

        self._Ct = drag_coefficient_heat_fluxes
        self._pbl_top = top_of_boundary_layer
        self._delta_pbl = boundary_layer_influence_height
        self._Cd0 = base_momentum_drag_coefficient
        self._Cd1 = wind_dependent_momentum_drag_coefficient
        self._Cm = maximum_momentum_drag_coefficient
        self._set_fortran_constants()
        super(SimplePhysics, self).__init__(**kwargs)

    def _set_fortran_constants(self):
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
        self._Rair = get_constant('gas_constant_of_dry_air', 'J/kg/degK')
        self._Rcond = get_constant('gas_constant_of_vapor_phase', 'J/kg/degK')
        self._radius = get_constant('planetary_radius', 'm')
        self._Omega = get_constant('planetary_rotation_rate', 's^-1')
        self._Lv = get_constant('latent_heat_of_condensation', 'J/kg')
        self._rho_condensible = get_constant('density_of_liquid_water', 'kg/m^3')
        phys.set_physical_constants(self._g, self._Cpd, self._Rair, self._Lv,
                                    self._Rcond, self._radius, self._Omega,
                                    self._rho_condensible, self._pbl_top,
                                    self._delta_pbl, self._Ct, self._Cd0,
                                    self._Cd1, self._Cm)

    @ensure_contiguous_state
    def array_call(self, state, timestep):
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
        self._set_fortran_constants()
        (t_out, u_out, v_out, q_out, precip_out,
         sensible_heat_flux, latent_heat_flux) = phys.get_new_state(
            state['eastward_wind'],
            state['northward_wind'],
            state['air_temperature'],
            state['air_pressure'],
            state['air_pressure_on_interface_levels'],
            state['specific_humidity'],
            state['surface_air_pressure'],
            state['surface_temperature'],
            state['surface_specific_humidity'],
            state['latitude'],
            timestep.total_seconds()
        )
        latent_heat_flux[latent_heat_flux < 0] = 0
        new_state = {
            'eastward_wind': u_out,
            'northward_wind': v_out,
            'air_temperature': t_out,
            'specific_humidity': q_out,
        }
        diagnostics = {
            'stratiform_precipitation_rate': precip_out,
            'surface_upward_sensible_heat_flux': sensible_heat_flux,
            'surface_upward_latent_heat_flux': latent_heat_flux,
        }
        return diagnostics, new_state
