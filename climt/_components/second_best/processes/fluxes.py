"""Surface flux process (BEST Section 8): bulk SHF/LHF with beta wetness."""
import numpy as np
from sympl import get_constant


class SurfaceFluxes:
    def __call__(self, drag, atmos, soil, soil_props, timestep):
        raise NotImplementedError


class BestSurfaceFluxes(SurfaceFluxes):
    def __call__(self, drag, atmos, soil, soil_props, timestep):
        Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK")
        Lv = get_constant("latent_heat_of_vaporization", "J/kg")
        Lf = get_constant("latent_heat_of_fusion", "J/kg")
        Li = Lv + Lf  # sublimation
        rho = atmos["air_density"]
        U = atmos["wind_speed"]
        C_Dh = drag["C_Dh"]
        dT = soil["surface_temperature"] - atmos["air_temperature"]
        shf = rho * Cpd * U * C_Dh * dT

        # wetness factor beta_u (soil): frozen sublimation term + exfiltration ratio
        W_Lu = soil["W_Lu"]
        W_Fu = soil["W_Fu"]
        c_u = C_Dh * U                         # soil conductance
        dq = soil["saturation_specific_humidity"] - atmos["air_specific_humidity"]
        E_usP = rho * c_u * dq                 # potential evaporation rate
        # exfiltration limit E_usmax (Eq. in Section 8) — clamp Theta to [0.01, 1]
        B = 4.0                                # placeholder overwritten below if given
        beta_u = np.clip(W_Lu, 0.0, 1.0)       # simplified when soil_props lacks B/K_H0
        if "B" in soil_props:
            B = soil_props["B"]
            K_H0 = soil_props["K_H0"]
            Theta = np.clip((W_Lu - 0.01) / max(1.0 - W_Fu, 1e-6), 1e-3, 1.0)
            rho_w = get_constant("density_of_liquid_water", "kg/m^3")
            Xv = soil_props["porosity"]
            psi0 = soil_props["psi_0"]
            K_HD = (-4 * K_H0 * B * psi0 * rho_w * Xv * (1 - W_Fu)) / (np.pi * timestep)
            E_usmax = K_HD * Theta ** (0.5 * B + 2) - K_H0 * Theta ** (2 * B + 3)
            frozen_term = (W_Fu * Lv / Li) if Li > 0 else 0.0
            ratio = np.clip(E_usmax / E_usP, 0.0, 1.0) if abs(E_usP) > 1e-12 else 0.0
            beta_u = np.clip(frozen_term + ratio, 0.0, 1.0)

        evaporation = beta_u * E_usP / rho     # m/s equivalent (per unit density)
        lhf = Lv * rho * evaporation
        u, v = atmos["u"], atmos["v"]
        momentum_flux = -rho * drag["C_Dm"] * U * np.array([u, v])
        return {"sensible_heat_flux": float(shf), "latent_heat_flux": float(lhf),
                "momentum_flux": momentum_flux, "evaporation": float(evaporation)}
