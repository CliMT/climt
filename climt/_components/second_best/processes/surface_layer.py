"""Surface-layer drag process (BEST Section 6)."""
import numpy as np
from sympl import get_constant


class SurfaceLayer:
    def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
        raise NotImplementedError

    def interpolate_to_height(self, drag, z0, z_mid, z_target,
                              surface_value, level_value, kind):
        """Diagnose a value at ``z_target`` between the surface and the lowest
        model level, consistent with the stability-dependent drag in ``drag``.

        The Monin-Obukhov stability functions are recovered from the bulk
        transfer coefficients this process returned (``C_Dm``/``C_Dh`` relative
        to the neutral ``C_DN``) and assumed to vary linearly with height across
        the surface layer, so the diagnosed profile is consistent with the
        surface fluxes rather than a neutral log-law. Reduces exactly to the
        neutral log-law when ``C_Dm == C_Dh == C_DN``.

        ``kind='wind'``: ``level_value`` is the lowest-level wind *speed*; the
        surface value is 0. Returns the wind speed at ``z_target``.
        ``kind='scalar'``: interpolates between ``surface_value`` (at the
        roughness height) and ``level_value`` (at ``z_mid``).
        """
        kappa = get_constant("von_karman_constant", "dimensionless")
        ln_mid = np.log(z_mid / z0)
        ln_tgt = np.log(z_target / z0)
        C_Dm = drag["C_Dm"]
        C_Dh = drag["C_Dh"]
        frac = z_target / z_mid
        if kind == "wind":
            # psi_m recovered from momentum drag; u* profile.
            psi_m = ln_mid - kappa / np.sqrt(C_Dm)
            denom = ln_mid - psi_m                      # == kappa/sqrt(C_Dm) > 0
            num = max(ln_tgt - psi_m * frac, 0.0)
            return level_value * num / denom
        # psi_h recovered from heat drag; scalar profile.
        psi_h = ln_mid - kappa * np.sqrt(C_Dm) / C_Dh
        denom = ln_mid - psi_h                          # == kappa*sqrt(C_Dm)/C_Dh > 0
        weight = (ln_tgt - psi_h * frac) / denom
        weight = min(max(weight, 0.0), 1.0)
        return surface_value + (level_value - surface_value) * weight


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
