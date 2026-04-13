import numpy as np
from sympl import TendencyComponent, get_constant

from ..common import compute_column_amount, compute_heating_rate, njit, prange
from ..optics.parmentier import (
    compute_rosseland_mean_opacity,
    compute_thermal_opacities,
    load_freedman2014_coefficients,
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
            self._freedman_coeffs = load_freedman2014_coefficients()
            self._num_bands = 2
        elif optics == "correlated_k":
            from ..optics.correlated_k import load_k_table

            self._table = load_k_table(table)
            self._num_bands = self._table["k_coefficients"].shape[1]
            self._num_gpts = self._table["k_coefficients"].shape[2]
            self._gas_names = list(self._table["gas_names"])
        else:
            raise ValueError(f"Unknown optics mode: {optics}")

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
            # H2O enters as specific humidity (kg/kg); all other gases are mole
            # fractions (mol/mol) per CF conventions.
            _GAS_CF_NAME = {
                "h2o": "specific_humidity",
                "co2": "mole_fraction_of_carbon_dioxide_in_air",
            }
            _GAS_UNITS = {"h2o": "kg/kg"}  # default for others: mole/mole
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
            "alias": "tau_cloud",
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
            # Molar masses (g/mol) for unit conversion from mole fraction to kg/kg.
            # H2O is already kg/kg (specific humidity); others are mole fractions.
            _MOLAR_MASS = {
                "h2o": 18.015,
                "co2": 44.010,
                "o3": 47.998,
                "ch4": 16.043,
                "n2o": 44.013,
                "o2": 31.998,
            }
            _M_AIR = 28.970  # g/mol dry air

            ngas = len(self._gas_names)
            gas_amounts = np.zeros((ngas, nlev, ncol))
            for ig, gas in enumerate(self._gas_names):
                q_gas_flat = state[gas].reshape(nlev, -1)
                # Convert mole fraction → mass mixing ratio (kg/kg) for non-H2O gases
                if gas != "h2o":
                    M_gas = _MOLAR_MASS.get(gas, _M_AIR)  # unknown gas: assume ~air
                    q_gas_flat = q_gas_flat * (M_gas / _M_AIR)
                gas_amounts[ig, :, :] = compute_column_amount(q_gas_flat, p_int_flat, g)

            tau, planck_src, surf_src, weights = self._correlated_k_optics(
                T_flat, p_flat, gas_amounts, T_surf_flat, sigma
            )
        else:
            raise NotImplementedError

        nband = tau.shape[0]
        ngpt = tau.shape[1]

        emissivity = state["emissivity"].reshape(nband, ncol)

        # Add cloud optical depth per band (shape: nlev, *horiz, nband).
        # Applied uniformly across g-points within each band.
        tau_cloud_flat = state["tau_cloud"].reshape(nlev, ncol, nband)  # (nlev, ncol, nband)
        # Rearrange to (nband, nlev, ncol) then broadcast over ngpt
        tau = tau + tau_cloud_flat.transpose(2, 0, 1)[:, np.newaxis, :, :]

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

    def _correlated_k_optics(self, T, p, gas_amounts, T_surf, sigma):
        from ..optics.correlated_k import compute_ck_optical_depth

        nlev, ncol = T.shape

        result = compute_ck_optical_depth(self._table, T, p, gas_amounts)
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
        nT = len(T_grid)
        nband_orig = planck_frac.shape[0]
        ngpt_orig = planck_frac.shape[1]

        overlap = str(self._table.get("overlap_method", np.array("additive")))
        ngas = self._table["k_coefficients"].shape[0]

        planck_src = np.zeros((nband, ngpt, nlev, ncol))
        surf_src = np.zeros((nband, ngpt, ncol))

        for icol in range(ncol):
            # Surface source
            T_s = T_surf[icol]
            iTs = np.searchsorted(T_grid, T_s) - 1
            iTs = max(0, min(iTs, nT - 2))
            fTs = (T_s - T_grid[iTs]) / (T_grid[iTs + 1] - T_grid[iTs])
            fTs = max(0.0, min(1.0, fTs))

            surf_planck = sigma * T_s**4
            for ib in range(nband):
                for igp in range(ngpt):
                    if overlap == "esft" and ngas > 1:
                        g_idx_orig = igp % ngpt_orig
                    else:
                        g_idx_orig = igp
                    ib_orig = min(ib, nband_orig - 1)
                    frac_s = (
                        planck_frac[ib_orig, g_idx_orig, iTs] * (1 - fTs)
                        + planck_frac[ib_orig, g_idx_orig, iTs + 1] * fTs
                    )
                    surf_src[ib, igp, icol] = frac_s * surf_planck

            # Layer source
            for k in range(nlev):
                T_l = T[k, icol]
                iTl = np.searchsorted(T_grid, T_l) - 1
                iTl = max(0, min(iTl, nT - 2))
                fTl = (T_l - T_grid[iTl]) / (T_grid[iTl + 1] - T_grid[iTl])
                fTl = max(0.0, min(1.0, fTl))

                layer_planck = sigma * T_l**4
                for ib in range(nband):
                    for igp in range(ngpt):
                        if overlap == "esft" and ngas > 1:
                            g_idx_orig = igp % ngpt_orig
                        else:
                            g_idx_orig = igp
                        ib_orig = min(ib, nband_orig - 1)
                        frac_l = (
                            planck_frac[ib_orig, g_idx_orig, iTl] * (1 - fTl)
                            + planck_frac[ib_orig, g_idx_orig, iTl + 1] * fTl
                        )
                        planck_src[ib, igp, k, icol] = frac_l * layer_planck

        return tau_raw, planck_src, surf_src, weights
