"""Surface-layer drag process (BEST Section 6)."""
import numpy as np
from sympl import get_constant


class SurfaceLayer:
    def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
        raise NotImplementedError


class BestSurfaceLayer(SurfaceLayer):
    def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
        kappa = get_constant("von_karman_constant", "dimensionless")
        g = get_constant("gravitational_acceleration", "m/s^2")
        U = np.maximum(wind_speed, 1e-3)
        C_DN = (kappa / (np.log(z_mid) - np.log(z0))) ** 2
        zeta = np.exp(-kappa / np.sqrt(C_DN))
        # Ri < 0 unstable (T_surf > T_air)
        Ri = -(g * z_mid / (T_surf * U * U)) * (T_surf - T_air)
        eps = 1.0 if area_type in ("sea", "sea_ice") else 0.01
        if np.ndim(Ri) == 0:
            Ri = float(Ri)
        unstable = Ri < 0.0
        # Momentum
        C_Dm = np.where(
            unstable,
            C_DN * (1 - 8 * Ri / (1 + 56.768 * C_DN * np.sqrt(np.abs(Ri) / zeta))),
            C_DN * ((1 - 4 * eps * Ri) ** 2) / (1 + 8 * (1 - eps) * Ri),
        )
        # Heat / scalars
        C_Dh = np.where(
            unstable,
            C_DN * (1 - 12 * Ri / (1 + 41.801 * C_DN * np.sqrt(np.abs(Ri) / zeta))),
            C_DN * ((1 - 4 * eps * Ri) / (1 + (6 - 4 * eps) * Ri)) ** 2,
        )
        return {"C_Dm": float(C_Dm), "C_Dh": float(C_Dh),
                "C_DN": float(C_DN), "Ri": float(Ri)}
