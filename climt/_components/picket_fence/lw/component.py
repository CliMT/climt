import numpy as np
from sympl import TendencyComponent, get_constant

from ..common import (
    MOLAR_MASS,
    MOLAR_MASS_DRY_AIR,
    compute_column_amount,
    compute_heating_rate,
)
from ..optics.parmentier import (
    compute_rosseland_mean_opacity,
    compute_thermal_opacities,
    load_freedman2014_coefficients,
    load_parmentier_coefficients,
    lookup_ratio_coefficients,
)
from .kernels import lw_transport, planck_sources_kernel


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
        self._has_co2_axis = False
        self._num_bands = 2 if optics == "parmentier" else None

        if optics == "parmentier":
            self._coefficients = load_parmentier_coefficients(coefficients)
            self._freedman_coeffs = load_freedman2014_coefficients()
            self._num_bands = 2
        elif optics == "correlated_k":
            from ..optics.correlated_k import load_k_table

            self._table = load_k_table(table)
            self._num_bands = self._table["k_coefficients"].shape[1]
            self._num_gpts = self._table["k_coefficients"].shape[2]
            self._gas_names = list(self._table["gas_names"])
            _has_h2o_axis = "h2o_vmr_grid" in self._table
            self._has_co2_axis = "co2_vmr_grid" in self._table
            self._fully_premixed = (self._gas_names == ["effective"]) and not _has_h2o_axis
            self._premixed_bg = (
                (self._gas_names == ["effective"] and _has_h2o_axis)
                or str(self._table.get("background_is_premixed", np.array("")))
                .lower() == "true"
            )
        else:
            raise ValueError(f"Unknown optics mode: {optics}")

        self._diagnostics_level = kwargs.pop("diagnostics_level", 0)
        from climt._core.initialization import set_num_longwave_bands
        set_num_longwave_bands(self._num_bands)
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
            "surface_longwave_emissivity": {
                "dims": ["num_longwave_bands", "*"],
                "units": "dimensionless",
                "alias": "emissivity",
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
            if self._premixed_bg:
                # Background is pre-mixed; H2O VMR is the only live axis.
                props["specific_humidity"] = {
                    "dims": ["mid_levels", "*"],
                    "units": "kg/kg",
                    "alias": "h2o",
                }
                if self._has_co2_axis:
                    props["mole_fraction_of_carbon_dioxide_in_air"] = {
                        "dims": ["mid_levels", "*"],
                        "units": "mole/mole",
                        "alias": "co2",
                    }
            elif not self._fully_premixed:
                # H2O enters as specific humidity (kg/kg); all other gases are
                # mole fractions (mol/mol) per CF conventions.
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
        # Cloud optical depth per band (spectrally-grey: same value across all
        # g-points within a band). Uses num_longwave_bands set at scheme init.
        props["longwave_optical_thickness_due_to_cloud"] = {
            "dims": ["mid_levels", "*", "num_longwave_bands"],
            "units": "dimensionless",
            "alias": "tau_cloud_lw",
        }
        return props

    @property
    def tendency_properties(self):
        return {
            "air_temperature": {"units": "degK s^-1"},
        }

    @property
    def diagnostic_properties(self):
        props = {
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
            "air_temperature_tendency_from_longwave": {
                "dims": ["mid_levels", "*"],
                "units": "degK day^-1",
            },
            "longwave_optical_depth_per_band": {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "dimensionless",
            },
            "longwave_transmittance_per_band": {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "dimensionless",
            },
            "air_temperature_tendency_from_longwave_per_band": {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "degK day^-1",
            },
        }
        if self._diagnostics_level >= 1:
            props["lw_layer_transmittance"] = {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "dimensionless",
            }
            props["lw_up_per_gpoint"] = {
                "dims": ["interface_levels", "*", "num_longwave_bands"],
                "units": "W m^-2",
            }
            props["lw_down_per_gpoint"] = {
                "dims": ["interface_levels", "*", "num_longwave_bands"],
                "units": "W m^-2",
            }
        return props

    @property
    def num_longwave_bands(self):
        return self._num_bands

    def array_call(self, state):
        T = state["T"]
        p = state["p"]
        p_int = state["p_int"]
        T_surf = state["T_surf"]

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
            T_irr_flat = state["T_irr"].reshape(-1)
            T_int_flat = state["T_int"].reshape(-1)
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
            weights = np.ones((tau.shape[0], tau.shape[1]))
        elif self._optics_mode == "correlated_k":
            ngas = len(self._gas_names)
            gas_amounts = np.zeros((ngas, nlev, ncol))
            h2o_vmr = None
            co2_vmr = None
            if self._fully_premixed:
                # k is m²/kg-of-air (all gases baked in); gas_amount = air column.
                one = np.ones((nlev, ncol))
                gas_amounts[0, :, :] = compute_column_amount(one, p_int_flat, g)
            elif self._premixed_bg:
                # Table k is m^2/kg-of-AIR with H2O as the runtime X-axis.
                # gas_amount is therefore column mass of air, not of H2O.
                # H2O enters only via specific_humidity → mole-fraction lookup
                # for the trilinear (T, logP, logX) k interpolation.
                q_h2o = state["h2o"].reshape(nlev, -1)  # kg/kg (specific humidity)
                one = np.ones_like(q_h2o)
                gas_amounts[0, :, :] = compute_column_amount(one, p_int_flat, g)
                M_H2O = MOLAR_MASS["h2o"] / MOLAR_MASS_DRY_AIR  # ratio (dimensionless)
                # x_H2O = (q/M_H2O) / (q/M_H2O + (1-q)/M_dry); express as
                # x = q / (q + (1-q) * M_H2O/M_dry) in mass-ratio form.
                h2o_vmr = q_h2o / np.maximum(
                    q_h2o + (1.0 - q_h2o) * M_H2O, 1e-30
                )
                if self._has_co2_axis:
                    # mole_fraction_of_carbon_dioxide_in_air IS already a VMR.
                    co2_vmr = state["co2"].reshape(nlev, -1)
            else:
                for ig, gas in enumerate(self._gas_names):
                    q_gas_flat = state[gas].reshape(nlev, -1)
                    if gas != "h2o":
                        M_gas = MOLAR_MASS.get(gas, MOLAR_MASS_DRY_AIR)
                        q_gas_flat = q_gas_flat * (M_gas / MOLAR_MASS_DRY_AIR)
                    gas_amounts[ig, :, :] = compute_column_amount(
                        q_gas_flat, p_int_flat, g
                    )

            tau, planck_src, surf_src, weights = self._correlated_k_optics(
                T_flat, p_flat, gas_amounts, T_surf_flat, sigma,
                h2o_vmr=h2o_vmr, co2_vmr=co2_vmr,
            )
        else:
            raise NotImplementedError

        nband = tau.shape[0]
        ngpt = tau.shape[1]

        emissivity = state["emissivity"].reshape(nband, ncol)

        # Add cloud optical depth per band (shape: nlev, *horiz, nband).
        # Applied uniformly across g-points within each band.
        tau_cloud_flat = state["tau_cloud_lw"].reshape(nlev, ncol, nband)  # (nlev, ncol, nband)
        # Rearrange to (nband, nlev, ncol) then broadcast over ngpt
        tau = tau + tau_cloud_flat.transpose(2, 0, 1)[:, np.newaxis, :, :]

        lw_result = lw_transport(
            T_flat,
            T_surf_flat,
            tau,
            planck_src,
            surf_src,
            emissivity,
            weights,
            sigma,
            diagnostics_level=self._diagnostics_level,
        )
        if self._diagnostics_level > 0:
            up_band, down_band, up_broad, down_broad, kernel_diag = lw_result
        else:
            up_band, down_band, up_broad, down_broad = lw_result

        net_flux = up_broad - down_broad
        heating_rate = compute_heating_rate(net_flux, p_int_flat, g, cpd)

        # Per-band optical depth: sum over g-points (weighted)
        # tau shape: (nband, ngpt, nlev, ncol)
        D = 1.66  # diffusivity factor
        tau_band = np.zeros((nband, nlev, ncol))
        for b in range(nband):
            for g_pt in range(ngpt):
                tau_band[b] += weights[b, g_pt] * tau[b, g_pt]
        tau_band_out = np.moveaxis(tau_band, 0, -1).reshape(orig_shape_T + (nband,))

        trans_band = np.exp(-D * tau_band)
        trans_band_out = np.moveaxis(trans_band, 0, -1).reshape(orig_shape_T + (nband,))

        # Per-band heating rate from per-band net flux divergence
        hr_band = np.zeros((nband, nlev, ncol))
        for b in range(nband):
            net_band = up_band[b] - down_band[b]  # (nlev+1, ncol)
            hr_band[b] = compute_heating_rate(net_band, p_int_flat, g, cpd) * 86400.0
        hr_band_out = np.moveaxis(hr_band, 0, -1).reshape(orig_shape_T + (nband,))

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

        diagnostics = {
            "upwelling_longwave_flux_in_air": up_broad_out,
            "downwelling_longwave_flux_in_air": down_broad_out,
            "upwelling_longwave_flux_in_air_per_band": up_band_out,
            "downwelling_longwave_flux_in_air_per_band": down_band_out,
            "air_temperature_tendency_from_longwave": heating_rate_kday,
            "longwave_optical_depth_per_band": tau_band_out,
            "longwave_transmittance_per_band": trans_band_out,
            "air_temperature_tendency_from_longwave_per_band": hr_band_out,
        }
        if self._diagnostics_level > 0:
            w_sum = weights.sum(axis=1)  # (nband,)

            def _avg_layer(arr):
                avg = np.einsum("bgnc,bg->bnc", arr, weights) / w_sum[:, None, None]
                return np.moveaxis(avg, 0, -1).reshape(orig_shape_T + (nband,))

            def _avg_iface(arr):
                avg = np.einsum("bgnc,bg->bnc", arr, weights) / w_sum[:, None, None]
                return np.moveaxis(avg, 0, -1).reshape(orig_shape_pint + (nband,))

            diagnostics["lw_layer_transmittance"] = _avg_layer(kernel_diag["transmittance"])
            diagnostics["lw_up_per_gpoint"] = _avg_iface(kernel_diag["up_per_gpoint"])
            diagnostics["lw_down_per_gpoint"] = _avg_iface(kernel_diag["down_per_gpoint"])

        return ({"T": tendency}, diagnostics)

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
                kappa_R = compute_rosseland_mean_opacity(
                    T[k, i], p[k, i], self._freedman_coeffs
                )
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

    def _correlated_k_optics(self, T, p, gas_amounts, T_surf, sigma,
                             h2o_vmr=None, co2_vmr=None):
        from ..optics.correlated_k import compute_ck_optical_depth

        nlev, ncol = T.shape

        result = compute_ck_optical_depth(
            self._table, T, p, gas_amounts, h2o_vmr=h2o_vmr, co2_vmr=co2_vmr
        )
        if isinstance(result, tuple):
            # ESFT mode: returns (tau, combined_weights)
            tau_raw, weights = result
        else:
            # Additive mode: returns tau only
            tau_raw = result
            weights = self._table["gpoint_weights"]

        nband = tau_raw.shape[0]
        ngpt = tau_raw.shape[1]

        planck_frac = self._table["planck_fraction"]  # (nband_orig, ngpt_orig, nT)
        T_grid = self._table["temperature_grid"]
        nband_orig = planck_frac.shape[0]
        ngpt_orig = planck_frac.shape[1]

        overlap = str(self._table.get("overlap_method", np.array("additive")))
        ngas = self._table["k_coefficients"].shape[0]

        planck_src = np.zeros((nband, ngpt, nlev, ncol))
        surf_src = np.zeros((nband, ngpt, ncol))

        is_esft = (overlap == "esft" and ngas > 1)
        planck_sources_kernel(
            np.ascontiguousarray(planck_frac),
            np.ascontiguousarray(T_grid, dtype=np.float64),
            np.ascontiguousarray(T, dtype=np.float64),
            np.ascontiguousarray(T_surf, dtype=np.float64),
            float(sigma),
            nband, ngpt, is_esft, ngpt_orig, nband_orig,
            planck_src, surf_src,
        )

        return tau_raw, planck_src, surf_src, weights
