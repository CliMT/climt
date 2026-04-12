# climt/_components/picket_fence/sw/component.py
import numpy as np
from sympl import TendencyComponent, get_constant

from ..common import compute_column_amount, compute_heating_rate, njit, prange
from ..optics.parmentier import load_parmentier_coefficients, lookup_ratio_coefficients
from .kernels import sw_two_stream


class PicketFenceShortwave(TendencyComponent):
    def __init__(
        self,
        optics="parmentier",
        table=None,
        coefficients="solar_composition",
        stellar_spectrum="sun",
        rosseland_mean_fit="freedman2014",
        **kwargs,
    ):
        self._optics_mode = optics

        if optics == "parmentier":
            self._coefficients = load_parmentier_coefficients(coefficients)
            self._num_bands = 3
            # Solar flux per band is computed dynamically from T_irr in array_call.
            # For Earth (non-irradiated) fallback, we store a default.
            self._default_solar_flux_per_band = np.array([1361.0 / 3.0] * 3)
        elif optics == "correlated_k":
            raise NotImplementedError("Correlated-k SW mode not yet implemented")
        else:
            raise ValueError(f"Unknown optics mode: {optics}")

        super(PicketFenceShortwave, self).__init__(**kwargs)

    @property
    def input_properties(self):
        props = {
            "air_temperature": {
                "dims": ["mid_levels", "*"],
                "units": "degK",
                "alias": "T",
            },
            "air_pressure": {"dims": ["mid_levels", "*"], "units": "Pa", "alias": "p"},
            "air_pressure_on_interface_levels": {
                "dims": ["interface_levels", "*"],
                "units": "Pa",
                "alias": "p_int",
            },
            "surface_temperature": {"dims": ["*"], "units": "degK", "alias": "T_surf"},
            "zenith_angle": {"dims": ["*"], "units": "radians", "alias": "zenith"},
            "surface_albedo_for_direct_shortwave": {
                "dims": ["*"],
                "units": "dimensionless",
                "alias": "albedo",
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
        return {"air_temperature": {"units": "degK s^-1"}}

    @property
    def diagnostic_properties(self):
        return {
            "upwelling_shortwave_flux_in_air": {
                "dims": ["interface_levels", "*"],
                "units": "W m^-2",
            },
            "downwelling_shortwave_flux_in_air": {
                "dims": ["interface_levels", "*"],
                "units": "W m^-2",
            },
            "upwelling_shortwave_flux_in_air_per_band": {
                "dims": ["interface_levels", "*", "num_shortwave_bands"],
                "units": "W m^-2",
            },
            "downwelling_shortwave_flux_in_air_per_band": {
                "dims": ["interface_levels", "*", "num_shortwave_bands"],
                "units": "W m^-2",
            },
            "shortwave_heating_rate": {
                "dims": ["mid_levels", "*"],
                "units": "degK day^-1",
            },
        }

    @property
    def num_shortwave_bands(self):
        return self._num_bands

    def array_call(self, state):
        T = np.asarray(getattr(state["T"], "data", state["T"]))
        p = np.asarray(getattr(state["p"], "data", state["p"]))
        p_int = np.asarray(getattr(state["p_int"], "data", state["p_int"]))
        zenith = np.asarray(getattr(state["zenith"], "data", state["zenith"]))
        albedo = np.asarray(getattr(state["albedo"], "data", state["albedo"]))

        orig_shape_T = T.shape
        orig_shape_pint = p_int.shape
        nlev = T.shape[0]

        T_flat = T.reshape(nlev, -1)
        p_flat = p.reshape(nlev, -1)
        p_int_flat = p_int.reshape(nlev + 1, -1)
        zenith_flat = zenith.reshape(-1)
        albedo_flat = albedo.reshape(-1)
        ncol = T_flat.shape[1]

        sigma = get_constant("stefan_boltzmann_constant", "W/m^2/K^4")
        g_const = get_constant("gravitational_acceleration", "m/s^2")
        cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")

        nband = self._num_bands
        ngpt = 1  # Parmentier mode: 1 g-point per band

        if self._optics_mode == "parmentier":
            T_irr = np.asarray(getattr(state["T_irr"], "data", state["T_irr"]))
            T_int = np.asarray(getattr(state["T_int"], "data", state["T_int"]))
            T_irr_flat = T_irr.reshape(-1)
            T_int_flat = T_int.reshape(-1)
            tau, ssa, asym = self._parmentier_sw_optics(
                T_flat, p_flat, p_int_flat, T_irr_flat, T_int_flat, g_const
            )
            # Compute TOA stellar flux from T_irr (Parmentier et al. 2015, Eq. 2):
            #   F_0 = sigma * T_irr^4, split equally across 3 visible bands.
            # For Earth-like (T_irr=0), fall back to stored default (1361 W/m^2).
            T_irr_max = T_irr_flat.max()
            if T_irr_max > 0:
                F0 = sigma * T_irr_max**4
                solar_flux_per_band = np.array([F0 / 3.0] * 3)
            else:
                solar_flux_per_band = self._default_solar_flux_per_band
        else:
            raise NotImplementedError

        solar_flux = solar_flux_per_band.reshape(nband, 1) * np.ones((nband, ngpt))
        weights = np.ones((nband, ngpt))

        up_band, down_band, up_broad, down_broad = sw_two_stream(
            tau, ssa, asym, zenith_flat, albedo_flat, solar_flux, weights
        )

        net_flux = up_broad - down_broad
        heating_rate = compute_heating_rate(net_flux, p_int_flat, g_const, cpd)

        tendency = heating_rate.reshape(orig_shape_T)
        up_broad_out = up_broad.reshape(orig_shape_pint)
        down_broad_out = down_broad.reshape(orig_shape_pint)
        heating_kday = heating_rate.reshape(orig_shape_T) * 86400.0

        up_band_out = np.moveaxis(up_band, 0, -1).reshape(orig_shape_pint + (nband,))
        down_band_out = np.moveaxis(down_band, 0, -1).reshape(
            orig_shape_pint + (nband,)
        )

        return (
            {"T": tendency},
            {
                "upwelling_shortwave_flux_in_air": up_broad_out,
                "downwelling_shortwave_flux_in_air": down_broad_out,
                "upwelling_shortwave_flux_in_air_per_band": up_band_out,
                "downwelling_shortwave_flux_in_air_per_band": down_band_out,
                "shortwave_heating_rate": heating_kday,
            },
        )

    def _parmentier_sw_optics(self, T, p, p_int, T_irr, T_int, g):
        """Compute SW optical depths for Parmentier mode (3 visible bands)."""
        nlev, ncol = T.shape
        nband = 3
        ngpt = 1

        tau = np.zeros((nband, ngpt, nlev, ncol))
        ssa = np.zeros((nband, ngpt, nlev, ncol))  # pure absorption for v1
        asym = np.zeros((nband, ngpt, nlev, ncol))

        for i in range(ncol):
            # Compute T_eff per column (same formula as LW, Lee et al. 2021 Eq. 20)
            A_B = 0.0
            mu_star = 0.25
            T_eff = (T_int[i] ** 4 + (1.0 - A_B) * mu_star * T_irr[i] ** 4) ** 0.25
            T_eff = max(T_eff, 100.0)
            gv1, gv2, gv3, beta, gamma_P, R = lookup_ratio_coefficients(
                self._coefficients, T_eff
            )
            gamma_vs = [gv1, gv2, gv3]

            for k in range(nlev):
                kappa_R = 1e-4  # placeholder Rosseland mean
                dp = abs(p_int[k + 1, i] - p_int[k, i])
                mass = dp / g

                for b in range(3):
                    kappa_v = gamma_vs[b] * kappa_R
                    tau[b, 0, k, i] = kappa_v * mass

        return tau, ssa, asym
