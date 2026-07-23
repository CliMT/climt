"""Subsurface heat + moisture + ice transport (BEST conduction section).

Implicit tridiagonal for the diffusive (density) terms; explicit Euler for
the freeze/melt source Gamma. Neumann bottom BC. Node 0 = bottom, node n-1
= surface. The linear solve uses the shared jitable Thomas solver
(:func:`climt._core.tridiagonal.solve_tridiagonal`) rather than
``scipy.sparse``/``spsolve``, so the whole numeric step is numba-jitable.
SecondBEST still does not depend on the ocean/ice ``snow_ice_column`` core.
"""
import numpy as np
from sympl import get_constant

from ...._core.backend import jit_compile
from ...._core.tridiagonal import solve_tridiagonal


class SubsurfaceTransport:
    def __call__(self, profiles, surface_flux_bc, timestep, dz):
        raise NotImplementedError


@jit_compile
def _best_subsurface_step(T, X_w, X_i, surface_flux_bc, dt, dz,
                          kappa, cv, Tf, Lf, rho_w):
    """Jitable implicit-diffusion + explicit freeze/melt subsurface step.

    ``sparse.diags`` used the standard banded layout (sub-diagonal =
    ``off[1:]``, super-diagonal = ``off[:-1]``), which for the constant
    ``off = -rr`` is just ``-rr`` on both off-diagonals; the two Neumann
    boundary rows only alter the first/last main-diagonal entries.
    """
    n = T.shape[0]
    rr = kappa * dt / (cv * dz * dz)

    lower = np.full(n - 1, -rr)
    upper = np.full(n - 1, -rr)
    main = np.full(n, 1.0 + 2.0 * rr)
    # Neumann bottom (node 0) and surface (node n-1)
    main[0] = 1.0 + rr
    main[-1] = 1.0 + rr

    rhs = T.copy()
    rhs[-1] = T[-1] + surface_flux_bc * dt / (cv * dz)
    T_diff = solve_tridiagonal(lower, main, upper, rhs)

    # --- freeze/melt source Gamma (explicit) ---
    Gamma = (cv / Lf) * (Tf - T_diff) / dt   # kg/m^3/s of ice
    # limit by available water (freezing) or ice (melting)
    max_freeze = rho_w * X_w / dt
    max_melt = rho_w * X_i / dt
    Gamma = np.minimum(np.maximum(Gamma, -max_melt), max_freeze)

    X_i_new = X_i + Gamma * dt / rho_w
    X_w_new = X_w - Gamma * dt / rho_w
    # latent heat release/absorption adjusts temperature back
    T_new = T_diff + Lf * Gamma * dt / cv
    if surface_flux_bc <= 0:
        T_new = np.minimum(T_new, Tf)
    return T_new, np.maximum(X_w_new, 0.0), np.maximum(X_i_new, 0.0)


class BestSubsurfaceTransport(SubsurfaceTransport):
    def __init__(self, thermal_conductivity=2.0, volumetric_heat_capacity=2.0e6):
        self._kappa = thermal_conductivity
        self._cv = volumetric_heat_capacity

    def __call__(self, profiles, surface_flux_bc, timestep, dz):
        T = np.array(profiles["T"], dtype=float)
        X_w = np.array(profiles["X_w"], dtype=float)
        X_i = np.array(profiles["X_i"], dtype=float)
        Tf = get_constant("freezing_temperature_of_liquid_phase", "degK")
        Lf = get_constant("latent_heat_of_fusion", "J/kg")
        rho_w = get_constant("density_of_liquid_water", "kg/m^3")

        T_new, X_w_new, X_i_new = _best_subsurface_step(
            T, X_w, X_i, float(surface_flux_bc), float(timestep), float(dz),
            self._kappa, self._cv, Tf, Lf, rho_w,
        )
        return {"T": T_new, "X_w": X_w_new, "X_i": X_i_new}
