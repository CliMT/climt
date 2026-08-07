from typing import NamedTuple

import numpy as np
from sympl import (
    Stepper,
    get_constant,
    initialize_numpy_arrays_with_properties,
)

from ..._core.backend import jit_compile, prange
from ..._core.tridiagonal import solve_tridiagonal


class BoundaryLayerParams(NamedTuple):
    Rd: float
    Cp: float
    g: float
    k: float
    z0: float
    fb: float
    P0: float
    Ric: float
    Lv: float


@jit_compile
def _richardson_diffusivity(Ri_a, u_fric, C, z, params):
    """Surface-layer diffusion coefficient K_b (thesis Eqn 2.8).

    Uses the *surface-layer* Richardson number ``Ri_a`` in the multiplier
    (NOT the local Ri of canonical Frierson 2006 Eqn 2.6), so K_b is
    continuous at ``Ri_a == 0`` and decays 1 -> 0 as ``Ri_a -> Ric``.
    """
    base = params.k * u_fric * np.sqrt(C) * z
    if Ri_a <= 0:
        return base
    return base / (
        1.0
        + Ri_a / params.Ric * np.log(z / params.z0) / (1.0 - Ri_a / params.Ric)
    )


@jit_compile
def _diffuse_profile(
    profile, p, p_int, rho, diff, dt, g, surface_exchange, surface_source
):
    """Implicit vertical diffusion of ``profile`` with a surface boundary term.

    Assembles the tridiagonal system and solves it with the shared Thomas
    solver. ``rho``/``diff`` have length ``num_layers - 1``.

    The lowest layer carries the surface boundary condition. With
    ``dp0 = p_int[0] - p_int[1]`` the layer-0 mass:

    * ``surface_exchange`` is the implicit bulk-exchange coefficient
      ``beta = g * rho_s * C * |v| * dt / dp0``, added to ``diag[0]``. It makes
      an unknown-dependent flux ``rho_s C |v| (X_s - X_0)`` backward-Euler.
    * ``surface_source`` is the explicit right-hand-side addition: ``beta * X_s``
      for that bulk flux, or ``g * dt * F / dp0`` for a prescribed flux ``F``
      (positive upward into the atmosphere).

    Both zero reproduces the no-flux solve.
    """
    num_layers = profile.shape[0]
    diag_m = np.zeros(num_layers)
    diag_p = np.zeros(num_layers)
    diag = np.zeros(num_layers)
    for i in range(num_layers):
        if i != 0:
            diag_m[i] = (
                g * g * rho[i - 1] * rho[i - 1] * diff[i - 1] * dt
                / (p[i - 1] - p[i])
                / (p_int[i] - p_int[i + 1])
            )
        if i != num_layers - 1:
            diag_p[i] = (
                g * g * rho[i] * rho[i] * diff[i] * dt
                / (p[i] - p[i + 1])
                / (p_int[i] - p_int[i + 1])
            )
        diag[i] = 1.0 + diag_m[i] + diag_p[i]
    diag[0] = diag[0] + surface_exchange
    rhs = profile.copy()
    rhs[0] = rhs[0] + surface_source
    return solve_tridiagonal(-diag_m[1:], diag, -diag_p[:-1], rhs)


@jit_compile(parallel=True)
def _boundary_layer_kernel(
    air_temperature, surface_temperature, air_pressure, air_pressure_int,
    surface_pressure, specific_humidity, surface_humidity, northward_wind,
    eastward_wind, sensible_heat_flux, latent_heat_flux,
    new_air_temperature, new_specific_humidity, new_northward_wind,
    new_eastward_wind, north_wind_stress, east_wind_stress, boundary_height,
    applied_sensible_flux, applied_latent_flux,
    params, num_cols, timestep, flux_mode,
):
    Rd = params.Rd
    Cp = params.Cp
    g = params.g
    k = params.k
    z0 = params.z0
    fb = params.fb
    P0 = params.P0
    Ri_c = params.Ric
    Lv = params.Lv

    for col in prange(num_cols):
        col_T = air_temperature[:, col]
        col_Ts = surface_temperature[col]
        col_p = air_pressure[:, col]
        col_p_int = air_pressure_int[:, col]
        col_ps = surface_pressure[col]
        col_q = specific_humidity[:, col]
        col_qs = surface_humidity[col]
        col_v = northward_wind[:, col]
        col_u = eastward_wind[:, col]

        col_v_int = 0.5 * (col_v[1:] + col_v[:-1])
        col_u_int = 0.5 * (col_u[1:] + col_u[:-1])
        col_T_int = 0.5 * (col_T[1:] + col_T[:-1])
        col_q_int = 0.5 * (col_q[1:] + col_q[:-1])
        col_rho = col_p_int[1:-1] / (
            Rd * (1.0 + 0.608 * col_q_int) * col_T_int
        )

        n = col_T_int.shape[0]

        pot_virt_temp = (
            col_T_int
            * np.power(P0 / col_p_int[1:-1], Rd / Cp)
            * (1.0 + 0.608 * col_q_int)
        )
        pot_virt_temp_surf = (
            col_Ts * np.power(P0 / col_ps, Rd / Cp) * (1.0 + 0.608 * col_qs)
        )

        z = np.zeros(n)
        z[0] = (
            Rd * (1.0 + 0.608 * col_q[0]) * col_T[0] / g
        ) * np.log(col_ps / col_p_int[1:-1][0])
        for i in range(1, n):
            z[i] = z[i - 1] + (
                Rd * (1.0 + 0.608 * col_q[i]) * col_T[i] / g
            ) * np.log(col_p_int[1:-1][i - 1] / col_p_int[1:-1][i])

        wind_int = np.sqrt(col_v_int * col_v_int + col_u_int * col_u_int)
        for i in range(wind_int.shape[0]):
            if wind_int[i] < 1.0:
                wind_int[i] = 1.0

        Ri_a = (
            g * z[0] * (pot_virt_temp[0] - pot_virt_temp_surf)
            / (pot_virt_temp_surf * wind_int[0] * wind_int[0])
        )
        if Ri_a < 0:
            C = k * k * np.power(np.log(z[0] / z0), -2)
        elif Ri_a < Ri_c:
            C = k * k * np.power(np.log(z[0] / z0), -2) * np.power(
                1.0 - Ri_a / Ri_c, 2
            )
        else:
            C = 0.0

        count = 0
        Rich = np.zeros(n)
        for i in range(n):
            Rich[i] = (
                g * z[i] * (pot_virt_temp[i] - pot_virt_temp[0])
                / (pot_virt_temp[0] * wind_int[i] * wind_int[i])
            )
            if Rich[i] > Ri_c:
                count = i + 1
                break
        h = z[count - 1]
        boundary_height[col] = h

        u_fric = wind_int[0]

        dp0 = col_p_int[0] - col_p_int[1]
        bulk_conductance = col_rho[0] * C * wind_int[0]
        beta = g * bulk_conductance * timestep / dp0

        if flux_mode == 1:
            scalar_exchange = beta
            source_T = beta * col_Ts
            source_q = beta * col_qs
        elif flux_mode == 2:
            scalar_exchange = 0.0
            source_T = g * timestep * sensible_heat_flux[col] / (Cp * dp0)
            source_q = g * timestep * latent_heat_flux[col] / (Lv * dp0)
        else:
            scalar_exchange = 0.0
            source_T = 0.0
            source_q = 0.0

        if flux_mode == 0:
            wind_exchange = 0.0
        else:
            wind_exchange = beta

        diff = np.zeros(n)
        for i in range(count):
            if z[i] < fb * h:
                diff[i] = _richardson_diffusivity(Ri_a, u_fric, C, z[i], params)
            else:
                diff[i] = (
                    _richardson_diffusivity(Ri_a, u_fric, C, fb * h, params)
                    * z[i] / (h * fb)
                    * np.power(1.0 - (z[i] - fb * h) / ((1.0 - fb) * h), 2)
                )

        new_air_temperature[:, col] = _diffuse_profile(
            col_T, col_p, col_p_int, col_rho, diff, timestep, g,
            scalar_exchange, source_T,
        )
        new_specific_humidity[:, col] = _diffuse_profile(
            col_q, col_p, col_p_int, col_rho, diff, timestep, g,
            scalar_exchange, source_q,
        )
        new_northward_wind[:, col] = _diffuse_profile(
            col_v, col_p, col_p_int, col_rho, diff, timestep, g,
            wind_exchange, 0.0,
        )
        new_eastward_wind[:, col] = _diffuse_profile(
            col_u, col_p, col_p_int, col_rho, diff, timestep, g,
            wind_exchange, 0.0,
        )

        applied_sensible_flux[col] = Cp * bulk_conductance * (
            col_Ts - new_air_temperature[0, col]
        )
        applied_latent_flux[col] = Lv * bulk_conductance * (
            col_qs - new_specific_humidity[0, col]
        )
        # Post-solve layer-0 winds, so the stress equals the momentum the
        # column actually lost. The surface value for momentum is no-slip
        # (X_s == 0), hence no (0 - u) difference here -- just the drag.
        north_wind_stress[col] = (
            bulk_conductance * new_northward_wind[0, col]
        )
        east_wind_stress[col] = bulk_conductance * new_eastward_wind[0, col]


_FLUX_MODES = {None: 0, 'bulk': 1, 'external': 2}


class SimpleBoundaryLayer(Stepper):
    """A simple boundary-layer scheme that diffuses heat, moisture and
    momentum upward from the lowest model level.

    This is the surface-flux / boundary-layer formulation of
    Frierson, Held & Zurita-Gotor (2006), with the surface-layer diffusion
    coefficient modified (thesis Eqn 2.8) to use the surface-layer Richardson
    number in its multiplier, making it continuous at ``Ri_a == 0``.

    The component assumes a surface-flux component has already applied the
    surface fluxes at the lowest model level; it then diffuses the resulting
    profiles using diffusion coefficients from a simplified Monin-Obukhov
    theory with a K-profile capped by a critical Richardson number.
    """

    input_properties = {
        'air_temperature': {'dims': ['mid_levels', '*'], 'units': 'degK'},
        'specific_humidity': {'dims': ['mid_levels', '*'], 'units': 'kg/kg'},
        'air_pressure': {'dims': ['mid_levels', '*'], 'units': 'Pa'},
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'], 'units': 'Pa',
        },
        'northward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
        'eastward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
        'surface_air_pressure': {'dims': ['*'], 'units': 'Pa'},
        'surface_temperature': {'dims': ['*'], 'units': 'degK'},
        'surface_specific_humidity': {'dims': ['*'], 'units': 'kg/kg'},
    }

    output_properties = {
        'air_temperature': {'dims': ['mid_levels', '*'], 'units': 'degK'},
        'specific_humidity': {'dims': ['mid_levels', '*'], 'units': 'kg/kg'},
        'northward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
        'eastward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
    }

    diagnostic_properties = {
        'northward_wind_stress': {'dims': ['*'], 'units': 'Pa'},
        'eastward_wind_stress': {'dims': ['*'], 'units': 'Pa'},
        'boundary_layer_height': {'dims': ['*'], 'units': 'm'},
    }

    def __init__(self, surface_fluxes='bulk', von_karman_constant=0.4,
                 roughness_length=0.0000321, specific_fraction=0.1,
                 reference_pressure=100000, critical_richardson_number=1,
                 **kwargs):
        """
        Args:
            surface_fluxes: how surface fluxes enter the lowest model level.

                * ``'bulk'`` (default): the component computes Frierson (2006)
                  bulk fluxes of heat, moisture and momentum from its own
                  exchange coefficient and applies them implicitly. It then
                  reports the applied heat and moisture fluxes as diagnostics.
                * ``'external'``: ``surface_upward_sensible_heat_flux`` and
                  ``surface_upward_latent_heat_flux`` become required inputs and
                  are applied as prescribed fluxes. Momentum still uses the
                  internal bulk drag.
                * ``None``: no surface exchange at all. The diffusion has
                  no-flux boundaries and conserves every column integral. A
                  separate surface-flux component must supply the fluxes.

                Pairing ``'bulk'`` with a component that already applies
                surface fluxes -- ``SimplePhysics(surface_fluxes=True)``, which
                is its default -- applies them twice.
            von_karman_constant: von Karman constant k.
            roughness_length: surface roughness length z0 (m).
            specific_fraction: surface-layer fraction fb of the boundary
                layer depth.
            reference_pressure: reference pressure P0 (Pa) for potential
                temperature.
            critical_richardson_number: critical Richardson number Ric that
                caps the diffusion and sets the boundary-layer top.
        """
        if surface_fluxes not in _FLUX_MODES:
            raise ValueError(
                "surface_fluxes must be 'bulk', 'external' or None, got %r"
                % (surface_fluxes,)
            )
        self._surface_fluxes = surface_fluxes
        self._flux_mode = _FLUX_MODES[surface_fluxes]
        self._k = von_karman_constant
        self._z0 = roughness_length
        self._fb = specific_fraction
        self._P0 = reference_pressure
        self._Ric = critical_richardson_number

        if surface_fluxes == 'bulk':
            # Instance-level override: only bulk mode reports the fluxes it
            # applied. In 'external' mode these names are inputs, so they must
            # not also be diagnostics.
            self.diagnostic_properties = dict(self.diagnostic_properties)
            self.diagnostic_properties['surface_upward_sensible_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }
            self.diagnostic_properties['surface_upward_latent_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }
        elif surface_fluxes == 'external':
            # Instance-level override: only external mode consumes fluxes
            # computed by a separate surface component. The fluxes were
            # computed with that component's own drag coefficient (e.g.
            # BucketHydrology's bulk_coefficient=0.0011), which differs from
            # this component's Monin-Obukhov C. That inconsistency is inherent
            # to the modular split and is not reconciled here.
            self.input_properties = dict(self.input_properties)
            self.input_properties['surface_upward_sensible_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }
            self.input_properties['surface_upward_latent_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }

        self._update_constants()
        super(SimpleBoundaryLayer, self).__init__(**kwargs)

    def _update_constants(self):
        self._Rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        self._Cp = get_constant(
            'heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1'
        )
        self._g = get_constant('gravitational_acceleration', 'm s^-2')
        self._Lv = get_constant('latent_heat_of_condensation', 'J kg^-1')

    def array_call(self, state, timestep):
        """Diffuse temperature, humidity and wind profiles for each column."""
        num_cols = state['air_temperature'].shape[1]

        new_state = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties
        )
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        params = BoundaryLayerParams(
            Rd=self._Rd, Cp=self._Cp, g=self._g, k=self._k, z0=self._z0,
            fb=self._fb, P0=self._P0, Ric=self._Ric, Lv=self._Lv,
        )

        # numba needs real arrays for every branch, so unused slots are zeros.
        zeros = np.zeros(num_cols)
        if self._flux_mode == 2:
            sensible_flux = state['surface_upward_sensible_heat_flux']
            latent_flux = state['surface_upward_latent_heat_flux']
        else:
            sensible_flux = zeros
            latent_flux = zeros

        if self._flux_mode == 1:
            applied_sensible = diagnostics['surface_upward_sensible_heat_flux']
            applied_latent = diagnostics['surface_upward_latent_heat_flux']
        else:
            applied_sensible = np.zeros(num_cols)
            applied_latent = np.zeros(num_cols)

        _boundary_layer_kernel(
            state['air_temperature'],
            state['surface_temperature'],
            state['air_pressure'],
            state['air_pressure_on_interface_levels'],
            state['surface_air_pressure'],
            state['specific_humidity'],
            state['surface_specific_humidity'],
            state['northward_wind'],
            state['eastward_wind'],
            sensible_flux,
            latent_flux,
            new_state['air_temperature'],
            new_state['specific_humidity'],
            new_state['northward_wind'],
            new_state['eastward_wind'],
            diagnostics['northward_wind_stress'],
            diagnostics['eastward_wind_stress'],
            diagnostics['boundary_layer_height'],
            applied_sensible,
            applied_latent,
            params,
            num_cols,
            timestep.total_seconds(),
            self._flux_mode,
        )

        return diagnostics, new_state
