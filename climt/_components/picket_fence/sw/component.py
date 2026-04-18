# climt/_components/picket_fence/sw/component.py
import numpy as np
from sympl import TendencyComponent, get_constant

from ..common import (
    MOLAR_MASS,
    MOLAR_MASS_DRY_AIR,
    compute_column_amount,
    compute_heating_rate,
    njit,
    prange,
)
from ..optics.parmentier import (
    compute_rosseland_mean_opacity,
    load_freedman2014_coefficients,
    load_parmentier_coefficients,
    lookup_ratio_coefficients,
)
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
            self._freedman_coeffs = load_freedman2014_coefficients()
            self._num_bands = 3
            # Solar flux per band is computed dynamically from T_irr in array_call.
            # For Earth (non-irradiated) fallback, we store a default.
            # Compute per-band solar flux from stellar spectrum.
            # Parmentier mode: 3 equal-width visible bands spanning the full spectrum.
            from ..optics.stellar import (
                integrate_spectrum_over_bands,
                load_stellar_spectrum,
            )

            try:
                spec = load_stellar_spectrum(stellar_spectrum)
                wn = spec["wavenumber"]
                wn_lo, wn_hi = wn.min(), wn.max()
                band_width = (wn_hi - wn_lo) / 3.0
                band_limits = np.array(
                    [
                        [wn_lo, wn_lo + band_width],
                        [wn_lo + band_width, wn_lo + 2 * band_width],
                        [wn_lo + 2 * band_width, wn_hi],
                    ]
                )
                self._solar_flux_per_band = integrate_spectrum_over_bands(
                    spec, band_limits
                )
            except (FileNotFoundError, KeyError):
                self._solar_flux_per_band = np.array([1361.0 / 3.0] * 3)
        elif optics == "correlated_k":
            from ..optics.correlated_k import load_k_table

            self._table = load_k_table(table)
            self._num_bands = self._table["k_coefficients"].shape[1]
            self._num_gpts = self._table["k_coefficients"].shape[2]
            self._gas_names = list(self._table["gas_names"])
            self._solar_source = self._table["solar_source_per_gpoint"]
            self._rayleigh = self._table.get("rayleigh_coefficient", None)
        else:
            raise ValueError(f"Unknown optics mode: {optics}")

        from climt._core.initialization import set_num_shortwave_bands

        set_num_shortwave_bands(self._num_bands)
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
            "flux_adjustment_for_earth_sun_distance": {
                "dims": ["*"],
                "units": "dimensionless",
                "alias": "earth_sun_factor",
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
        elif self._optics_mode == "correlated_k":
            _GAS_CF_NAME = {
                "h2o": "specific_humidity",
                "co2": "mole_fraction_of_carbon_dioxide_in_air",
            }
            _GAS_UNITS = {"h2o": "kg/kg"}
            for gas in self._gas_names:
                cf_name = _GAS_CF_NAME.get(gas, f"mole_fraction_of_{gas}_in_air")
                units = _GAS_UNITS.get(gas, "mole/mole")
                props[cf_name] = {
                    "dims": ["mid_levels", "*"],
                    "units": units,
                    "alias": gas,
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
            "shortwave_optical_depth_per_band": {
                "dims": ["mid_levels", "*", "num_shortwave_bands"],
                "units": "dimensionless",
            },
            "shortwave_heating_rate_per_band": {
                "dims": ["mid_levels", "*", "num_shortwave_bands"],
                "units": "degK day^-1",
            },
        }

    @property
    def num_shortwave_bands(self):
        return self._num_bands

    def array_call(self, state):
        T = state["T"]
        p = state["p"]
        p_int = state["p_int"]

        orig_shape_T = T.shape
        orig_shape_pint = p_int.shape
        nlev = T.shape[0]

        T_flat = T.reshape(nlev, -1)
        p_flat = p.reshape(nlev, -1)
        p_int_flat = p_int.reshape(nlev + 1, -1)
        zenith_flat = state["zenith"].reshape(-1)
        albedo_flat = state["albedo"].reshape(-1)
        ncol = T_flat.shape[1]

        g_const = get_constant("gravitational_acceleration", "m/s^2")
        cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")

        if self._optics_mode == "parmentier":
            sigma = get_constant("stefan_boltzmann_constant", "W/m^2/K^4")
            nband = self._num_bands
            ngpt = 1

            T_irr_flat = state["T_irr"].reshape(-1)
            T_int_flat = state["T_int"].reshape(-1)
            tau, ssa, asym = self._parmentier_sw_optics(
                T_flat, p_flat, p_int_flat, T_irr_flat, T_int_flat, g_const
            )
            T_irr_max = T_irr_flat.max()
            if T_irr_max > 0:
                F0 = sigma * T_irr_max**4
                solar_flux_per_band = np.array([F0 / 3.0] * 3)
            else:
                solar_flux_per_band = self._solar_flux_per_band

            earth_sun_factor = float(state["earth_sun_factor"].reshape(-1)[0])
            solar_flux = (
                solar_flux_per_band.reshape(nband, 1)
                * np.ones((nband, ngpt))
                * earth_sun_factor
            )
            weights = np.ones((nband, ngpt))

        elif self._optics_mode == "correlated_k":
            ngas = len(self._gas_names)
            gas_amounts = np.zeros((ngas, nlev, ncol))
            for ig, gas in enumerate(self._gas_names):
                q_gas_flat = state[gas].reshape(nlev, -1)
                if gas != "h2o":
                    M_gas = MOLAR_MASS.get(gas, MOLAR_MASS_DRY_AIR)
                    q_gas_flat = q_gas_flat * (M_gas / MOLAR_MASS_DRY_AIR)
                gas_amounts[ig, :, :] = compute_column_amount(
                    q_gas_flat, p_int_flat, g_const
                )

            from ..optics.correlated_k import compute_ck_optical_depth

            result = compute_ck_optical_depth(self._table, T_flat, p_flat, gas_amounts)
            if isinstance(result, tuple):
                tau_abs, weights = result
            else:
                tau_abs = result
                weights = self._table["gpoint_weights"]

            nband = tau_abs.shape[0]
            ngpt = tau_abs.shape[1]

            ssa = np.zeros((nband, ngpt, nlev, ncol))
            asym = np.zeros((nband, ngpt, nlev, ncol))
            tau = tau_abs.copy()

            if self._rayleigh is not None:
                for b in range(nband):
                    for k in range(nlev):
                        dp = abs(p_int_flat[k + 1, :] - p_int_flat[k, :])
                        tau_ray = self._rayleigh[b] * dp / g_const
                        for gp in range(ngpt):
                            tau_total = tau_abs[b, gp, k, :] + tau_ray
                            ssa[b, gp, k, :] = np.where(
                                tau_total > 0,
                                tau_ray / tau_total,
                                0.0,
                            )
                            tau[b, gp, k, :] = tau_total

            earth_sun_factor = float(state["earth_sun_factor"].reshape(-1)[0])
            if ngpt == self._solar_source.shape[1]:
                solar_flux = self._solar_source * earth_sun_factor
            else:
                # TODO: ESFT expands g-points as the outer product across gases.
                # The solar source should be replicated per original g-point index,
                # not cycled with %. Replace with a proper mapping once ESFT solar
                # distribution is specified.
                ngpt_orig = self._solar_source.shape[1]
                solar_flux = np.zeros((nband, ngpt))
                for b in range(nband):
                    for idx in range(ngpt):
                        g_idx = idx % ngpt_orig
                        solar_flux[b, idx] = (
                            self._solar_source[b, g_idx] * earth_sun_factor
                        )
        else:
            raise NotImplementedError

        up_band, down_band, up_broad, down_broad = sw_two_stream(
            tau, ssa, asym, zenith_flat, albedo_flat, solar_flux, weights
        )

        net_flux = up_broad - down_broad
        heating_rate = compute_heating_rate(net_flux, p_int_flat, g_const, cpd)

        # Per-band optical depth: sum over g-points (weighted)
        tau_band = np.zeros((nband, nlev, ncol))
        for b in range(nband):
            for g_pt in range(ngpt):
                tau_band[b] += weights[b, g_pt] * tau[b, g_pt]
        tau_band_out = np.moveaxis(tau_band, 0, -1).reshape(orig_shape_T + (nband,))

        # Per-band heating rate from per-band net flux divergence
        hr_band = np.zeros((nband, nlev, ncol))
        for b in range(nband):
            net_band = up_band[b] - down_band[b]
            hr_band[b] = compute_heating_rate(net_band, p_int_flat, g_const, cpd) * 86400.0
        hr_band_out = np.moveaxis(hr_band, 0, -1).reshape(orig_shape_T + (nband,))

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
                "shortwave_optical_depth_per_band": tau_band_out,
                "shortwave_heating_rate_per_band": hr_band_out,
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
                kappa_R = compute_rosseland_mean_opacity(
                    T[k, i], p[k, i], self._freedman_coeffs
                )
                dp = abs(p_int[k + 1, i] - p_int[k, i])
                mass = dp / g

                for b in range(3):
                    kappa_v = gamma_vs[b] * kappa_R
                    tau[b, 0, k, i] = kappa_v * mass

        return tau, ssa, asym
