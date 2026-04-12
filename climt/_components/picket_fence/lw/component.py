import numpy as np
from sympl import TendencyComponent, get_constant

from ..common import compute_column_amount, compute_heating_rate, njit, prange
from ..optics.parmentier import (
    compute_thermal_opacities,
    load_parmentier_coefficients,
    lookup_ratio_coefficients,
)
from .kernels import lw_transport


class PicketFenceLongwave(TendencyComponent):
    def __init__(
        self,
        optics="parmentier",
        table=None,
        coefficients="solar_composition",
        rosseland_mean_fit="freedman2014",
        **kwargs,
    ):
        self._optics_mode = optics
        self._num_bands = 2 if optics == "parmentier" else None

        if optics == "parmentier":
            self._coefficients = load_parmentier_coefficients(coefficients)
            self._num_bands = 2
        elif optics == "correlated_k":
            raise NotImplementedError("Correlated-k mode not yet implemented")
        else:
            raise ValueError(f"Unknown optics mode: {optics}")

        super(PicketFenceLongwave, self).__init__(**kwargs)

    @property
    def input_properties(self):
        props = {
            "air_temperature": {
                "dims": ["mid_levels", "*"],
                "units": "degK",
                "alias": "T",
            },
            "air_pressure": {
                "dims": ["mid_levels", "*"],
                "units": "Pa",
                "alias": "p",
            },
            "air_pressure_on_interface_levels": {
                "dims": ["interface_levels", "*"],
                "units": "Pa",
                "alias": "p_int",
            },
            "surface_temperature": {
                "dims": ["*"],
                "units": "degK",
                "alias": "T_surf",
            },
        }
        if self._optics_mode == "parmentier":
            props["irradiation_temperature"] = {
                "dims": ["*"],
                "units": "degK",
                "alias": "T_irr",
            }
            props["internal_temperature"] = {
                "dims": ["*"],
                "units": "degK",
                "alias": "T_int",
            }
        return props

    @property
    def tendency_properties(self):
        return {
            "air_temperature": {"units": "degK s^-1"},
        }

    @property
    def diagnostic_properties(self):
        return {
            "upwelling_longwave_flux_in_air": {
                "dims": ["interface_levels", "*"],
                "units": "W m^-2",
            },
            "downwelling_longwave_flux_in_air": {
                "dims": ["interface_levels", "*"],
                "units": "W m^-2",
            },
            "upwelling_longwave_flux_in_air_per_band": {
                "dims": ["interface_levels", "*", "num_longwave_bands"],
                "units": "W m^-2",
            },
            "downwelling_longwave_flux_in_air_per_band": {
                "dims": ["interface_levels", "*", "num_longwave_bands"],
                "units": "W m^-2",
            },
            "longwave_heating_rate": {
                "dims": ["mid_levels", "*"],
                "units": "degK day^-1",
            },
        }

    @property
    def num_longwave_bands(self):
        return self._num_bands

    def array_call(self, state):
        T = np.asarray(getattr(state["T"], "data", state["T"]))
        p = np.asarray(getattr(state["p"], "data", state["p"]))
        p_int = np.asarray(getattr(state["p_int"], "data", state["p_int"]))
        T_surf = np.asarray(getattr(state["T_surf"], "data", state["T_surf"]))

        orig_shape_T = T.shape
        orig_shape_pint = p_int.shape
        nlev = T.shape[0]

        T_flat = T.reshape(nlev, -1)
        p_flat = p.reshape(nlev, -1)
        p_int_flat = p_int.reshape(nlev + 1, -1)
        T_surf_flat = T_surf.reshape(-1)
        ncol = T_flat.shape[1]

        sigma = get_constant("stefan_boltzmann_constant", "W/m^2/K^4")
        g = get_constant("gravitational_acceleration", "m/s^2")
        cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")

        if self._optics_mode == "parmentier":
            T_irr = np.asarray(getattr(state["T_irr"], "data", state["T_irr"]))
            T_int = np.asarray(getattr(state["T_int"], "data", state["T_int"]))
            T_irr_flat = T_irr.reshape(-1)
            T_int_flat = T_int.reshape(-1)
            tau, planck_src, surf_src = self._parmentier_optics(
                T_flat,
                p_flat,
                p_int_flat,
                T_surf_flat,
                T_irr_flat,
                T_int_flat,
                sigma,
                g,
            )
        else:
            raise NotImplementedError

        nband = tau.shape[0]
        ngpt = tau.shape[1]
        emissivity = np.ones((nband, ncol))
        weights = np.ones((nband, ngpt))

        up_band, down_band, up_broad, down_broad = lw_transport(
            T_flat,
            T_surf_flat,
            tau,
            planck_src,
            surf_src,
            emissivity,
            weights,
            sigma,
        )

        net_flux = up_broad - down_broad
        heating_rate = compute_heating_rate(net_flux, p_int_flat, g, cpd)

        # Reshape outputs
        tendency = heating_rate.reshape(orig_shape_T)
        up_broad_out = up_broad.reshape(orig_shape_pint)
        down_broad_out = down_broad.reshape(orig_shape_pint)
        heating_rate_kday = heating_rate.reshape(orig_shape_T) * 86400.0

        # Per-band: (nband, nlev+1, ncol) -> (nlev+1, ncol, nband) -> reshape
        up_band_out = np.moveaxis(up_band, 0, -1)  # (nlev+1, ncol, nband)
        down_band_out = np.moveaxis(down_band, 0, -1)
        target_band_shape = orig_shape_pint + (nband,)
        up_band_out = up_band_out.reshape(target_band_shape)
        down_band_out = down_band_out.reshape(target_band_shape)

        return (
            {"T": tendency},
            {
                "upwelling_longwave_flux_in_air": up_broad_out,
                "downwelling_longwave_flux_in_air": down_broad_out,
                "upwelling_longwave_flux_in_air_per_band": up_band_out,
                "downwelling_longwave_flux_in_air_per_band": down_band_out,
                "longwave_heating_rate": heating_rate_kday,
            },
        )

    def _parmentier_optics(self, T, p, p_int, T_surf, T_irr, T_int, sigma, g):
        """Compute optical depths and sources for Parmentier mode."""
        nlev, ncol = T.shape
        nband = 2
        ngpt = 1

        tau = np.zeros((nband, ngpt, nlev, ncol))
        planck_src = np.zeros((nband, ngpt, nlev, ncol))
        surf_src = np.zeros((nband, ngpt, ncol))

        for i in range(ncol):
            # Compute T_eff per column from irradiation and internal temperatures
            # Following Lee et al. (2021) Eq. 20:
            #   T_eff^4 = T_int^4 + (1 - A_B) * mu_* * T_irr^4
            # For now, assume A_B=0 (Bond albedo) and mu_*=1/4 (global average).
            # These can be refined once the SW solver provides a proper Bond albedo.
            A_B = 0.0
            mu_star = 0.25
            T_eff = (T_int[i] ** 4 + (1.0 - A_B) * mu_star * T_irr[i] ** 4) ** 0.25
            # Floor at 100 K to avoid degenerate coefficient lookups
            T_eff = max(T_eff, 100.0)

            gv1, gv2, gv3, beta, gamma_P, R = lookup_ratio_coefficients(
                self._coefficients, T_eff
            )

            for k in range(nlev):
                # Simplified Rosseland mean opacity (placeholder)
                # A proper Freedman 2014 fit would use T[k,i] and p[k,i]
                kappa_R = 1e-4  # m^2/kg — Earth-like placeholder
                kappa_1, kappa_2 = compute_thermal_opacities(kappa_R, gamma_P, beta, R)

                # Layer mass: dp / g
                dp = abs(p_int[k + 1, i] - p_int[k, i])
                mass = dp / g

                tau[0, 0, k, i] = kappa_1 * mass
                tau[1, 0, k, i] = kappa_2 * mass

                planck_val = sigma * T[k, i] ** 4
                planck_src[0, 0, k, i] = beta * planck_val
                planck_src[1, 0, k, i] = (1.0 - beta) * planck_val

            surf_planck = sigma * T_surf[i] ** 4
            surf_src[0, 0, i] = beta * surf_planck
            surf_src[1, 0, i] = (1.0 - beta) * surf_planck

        return tau, planck_src, surf_src
