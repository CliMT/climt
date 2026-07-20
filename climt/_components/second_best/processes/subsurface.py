"""Subsurface heat + moisture + ice transport (BEST conduction section).

Implicit tridiagonal for the diffusive (density) terms; explicit Euler for
the freeze/melt source Gamma. Neumann bottom BC. Node 0 = bottom, node n-1
= surface. This solver is self-contained (its own scipy sparse tridiagonal)
so SecondBEST does not depend on the ocean/ice plan's snow_ice_column core.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from sympl import get_constant


class SubsurfaceTransport:
    def __call__(self, profiles, surface_flux_bc, timestep, dz):
        raise NotImplementedError


class BestSubsurfaceTransport(SubsurfaceTransport):
    def __init__(self, thermal_conductivity=2.0, volumetric_heat_capacity=2.0e6):
        self._kappa = thermal_conductivity
        self._cv = volumetric_heat_capacity

    def __call__(self, profiles, surface_flux_bc, timestep, dz):
        T = np.array(profiles["T"], dtype=float)
        X_w = np.array(profiles["X_w"], dtype=float)
        X_i = np.array(profiles["X_i"], dtype=float)
        n = T.shape[0]
        dt = timestep
        Tf = get_constant("freezing_temperature_of_liquid_phase", "degK")
        Lf = get_constant("latent_heat_of_fusion", "J/kg")
        rho_w = get_constant("density_of_liquid_water", "kg/m^3")

        # --- implicit heat diffusion ---
        rr = self._kappa * dt / (self._cv * dz * dz)
        main = np.full(n, 1 + 2 * rr)
        off = np.full(n, -rr)
        # Build in LIL format so boundary-row overwrites below are efficient;
        # convert to CSC only once assembly is complete, right before spsolve.
        A = sparse.diags([off[1:], main, off[:-1]], [-1, 0, 1], format="lil")
        # Neumann bottom (node 0): dT/dz = 0
        A[0, 0] = 1 + rr; A[0, 1] = -rr
        rhs = T.copy()
        # surface flux Neumann (node n-1): kappa dT/dz = -flux
        A[-1, -1] = 1 + rr; A[-1, -2] = -rr
        rhs[-1] = T[-1] + surface_flux_bc * dt / (self._cv * dz)
        A = A.tocsc()
        T_diff = spsolve(A, rhs)

        # --- freeze/melt source Gamma (explicit) ---
        # potential ice production rate: drives T toward Tf, converting phase
        Gamma = (self._cv / Lf) * (Tf - T_diff) / dt   # kg/m^3/s of ice
        # limit by available water (freezing) or ice (melting)
        max_freeze = rho_w * X_w / dt
        max_melt = rho_w * X_i / dt
        Gamma = np.clip(Gamma, -max_melt, max_freeze)

        X_i_new = X_i + Gamma * dt / rho_w
        X_w_new = X_w - Gamma * dt / rho_w
        # latent heat release/absorption adjusts temperature back
        T_new = T_diff + Lf * Gamma * dt / self._cv
        T_new = np.minimum(T_new, Tf) if surface_flux_bc <= 0 else T_new
        return {"T": T_new, "X_w": np.clip(X_w_new, 0.0, None),
                "X_i": np.clip(X_i_new, 0.0, None)}
