import numpy as np
from sympl import ImplicitTendencyComponent, get_constant

from ..._core import bolton_q_sat


class EmanuelConvectionPython(ImplicitTendencyComponent):
    """
    Pure Python implementation of the Emanuel Convection scheme (version 4.3c).
    """

    input_properties = {
        "air_temperature": {"dims": ["mid_levels", "*"], "units": "degK"},
        "specific_humidity": {"dims": ["mid_levels", "*"], "units": "g/g"},
        "air_pressure": {"dims": ["mid_levels", "*"], "units": "mbar"},
        "air_pressure_on_interface_levels": {
            "dims": ["interface_levels", "*"],
            "units": "mbar",
        },
        "eastward_wind": {"dims": ["mid_levels", "*"], "units": "m s^-1"},
        "northward_wind": {"dims": ["mid_levels", "*"], "units": "m s^-1"},
    }

    tendency_properties = {
        "air_temperature": {"units": "degK s^-1"},
        "specific_humidity": {"units": "g g^-1 s^-1"},
        "eastward_wind": {"units": "m s^-2"},
        "northward_wind": {"units": "m s^-2"},
    }

    diagnostic_properties = {
        "convective_precipitation_rate": {"dims": ["*"], "units": "mm day^-1"},
    }

    def __init__(
        self,
        elcrit=0.0011,
        tlcrit=-55.0,
        entp=1.5,
        sigd=0.05,
        sigs=0.12,
        omtrain=50.0,
        omtsnow=5.5,
        coeffr=1.0,
        coeffs=0.8,
        cu=0.7,
        beta=10.0,
        dtmax=0.9,
        alpha=0.1,
        damp=0.1,
        ipbl=0,
        minorig=1,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.params = {
            "elcrit": elcrit,
            "tlcrit": tlcrit,
            "entp": entp,
            "sigd": sigd,
            "sigs": sigs,
            "omtrain": omtrain,
            "omtsnow": omtsnow,
            "coeffr": coeffr,
            "coeffs": coeffs,
            "cu": cu,
            "beta": beta,
            "dtmax": dtmax,
            "alpha": alpha,
            "damp": damp,
            "ipbl": ipbl,
            "minorig": minorig,
        }
        self._update_constants()
        self.cbmf = None  # Persistent cloud base mass flux

    def _update_constants(self):
        self.CPD = get_constant(
            "heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"
        )
        self.CPV = get_constant("heat_capacity_of_vapor_phase", "J/kg/degK")
        self.CL = 2500.0  # Specific heat of liquid water (Note: Artificially small value from Emanuel Scheme)
        self.RV = get_constant("gas_constant_of_vapor_phase", "J/kg/degK")
        self.RD = get_constant("gas_constant_of_dry_air", "J/kg/degK")
        self.LV0 = get_constant("latent_heat_of_condensation", "J/kg")
        self.G = get_constant("gravitational_acceleration", "m/s^2")
        self.ROWL = 1000.0  # Density of liquid water
        self.DELT0 = 300.0  # Reference time scale for relaxation

    def array_call(self, state, timestep):
        """
        Step the Emanuel Convection scheme.
        """
        self._update_constants()
        nlev, ncol = state["air_temperature"].shape

        # Initialize persistent state
        if self.cbmf is None or self.cbmf.shape[0] != ncol:
            self.cbmf = np.zeros(ncol)

        # 1. Thermodynamics
        p = state["air_pressure"]
        ph = state["air_pressure_on_interface_levels"]
        t = state["air_temperature"]
        q = state["specific_humidity"]

        # Calculate Saturation Specific Humidity using Bolton's formula (as in Fortran)
        qs = bolton_q_sat(t, p * 100.0, self.RD, self.RV)

        eps = self.RD / self.RV
        epsi = 1.0 / eps

        # Virtual Temperature
        tv = t * (1.0 + q * epsi - q)

        # Geopotential (hydrostatic)
        gz = np.zeros((nlev + 1, ncol))
        for i in range(1, nlev + 1):
            tvx = tv[i - 1]
            tvy = tv[i - 2] if i > 1 else tv[i - 1]
            p_prev = p[i - 2] if i > 1 else ph[0]
            gz[i] = (
                gz[i - 1]
                + 0.5 * self.RD * (tvx + tvy) * (p_prev - p[i - 1]) / ph[i - 1]
            )

        # Moist Static Energy (h_m) and Enthalpy (h)
        lv = self.LV0 - (self.CL - self.CPV) * (t - 273.15)
        cpn = self.CPD * (1.0 - q) + self.CPV * q
        h = t * cpn + gz[:nlev]
        hm = (self.CPD * (1.0 - q) + self.CL * q) * (t - t[0]) + lv * q + gz[:nlev]

        # 2. Find Parcel Origin Level (NK)
        nk = self._find_parcel_origin(hm, self.params["minorig"], nlev, ncol)

        # 3. Calculate LCL
        plcl = self._calculate_lcl(p, t, q, qs, nk, ncol)

        # 4. Find Cloud Base (ICB)
        icb = self._find_icb(p, plcl, nk, nlev, ncol)

        # 5. Parcel Lifting (TLIFT)
        tp, qp, clw, tvp = self._tlift(state, nk, icb, tv, gz)

        # Relaxation rate and damping (normalized by timestep)
        damps = self.params["damp"] * timestep.total_seconds() / self.params["beta"]
        alpha_eff = 0.1 * self.params["alpha"]

        # Buoyancy at LCL
        dtma = np.zeros(ncol)
        for c in range(ncol):
            idx = icb[c]
            if idx < nlev:
                dtma[c] = tvp[idx, c] - tv[idx, c] + self.params["dtmax"]
            else:
                dtma[c] = 0.0

        # Update persistent CBMF
        self.cbmf = (1.0 - damps) * self.cbmf + alpha_eff * dtma
        self.cbmf = np.maximum(self.cbmf, 0.0)

        # --- Simplified Mixing and Tendency (Preserved for now) ---
        m_rates = np.zeros_like(t)
        for i in range(1, nlev):
            b_curr = tvp[i] - tv[i]
            b_prev = tvp[i - 1] - tv[i - 1]
            db = np.abs(b_curr - b_prev)
            dp = np.abs(p[i] - p[i - 1])
            m_rates[i] = db + self.params["entp"] * 0.02 * dp

        dbosum = np.sum(m_rates, axis=0)
        dbosum = np.maximum(dbosum, 1e-20)

        m_flux = np.zeros_like(t)
        for i in range(nlev):
            m_flux[i] = self.cbmf * m_rates[i] / dbosum

        ep = np.zeros_like(t)
        elcrit = self.params["elcrit"]
        tlcrit = self.params["tlcrit"]
        epmax = 0.999

        for i in range(nlev):
            tca = tp[i] - 273.15
            elacrit = np.where(tca >= 0.0, elcrit, elcrit * (1.0 - tca / tlcrit))
            elacrit = np.maximum(elacrit, 0.0)
            ep[i] = epmax * (1.0 - elacrit / np.maximum(clw[i], 1e-8))
            ep[i] = np.clip(ep[i], 0.0, epmax)

        mp = np.zeros_like(t)
        evap = np.zeros_like(t)
        water = np.zeros_like(t)

        sigd = self.params["sigd"]
        sigs = self.params["sigs"]
        ginv = 1.0 / self.G

        inb_idx = nlev - 1

        for i in range(inb_idx, -1, -1):
            wd_train = self.G * ep[i] * m_flux[i] * clw[i]
            is_rain = t[i] > 273.15
            coeff = np.where(is_rain, self.params["coeffr"], self.params["coeffs"])
            wt_vel = np.where(is_rain, self.params["omtrain"], self.params["omtsnow"])

            # Use calculated qs
            qs_local = qs[i]
            qsm = 0.5 * (q[i] + qp[i])
            afac = (
                coeff
                * p[i]
                * np.maximum(qs_local - qsm, 0.0)
                / (1.0e4 + 2.0e3 * p[i] * qs_local)
            )

            b6 = 100.0 * (ph[i] - ph[i + 1]) * sigs * afac / wt_vel
            water_prev = water[i + 1] if i < inb_idx else 0.0
            c6 = (water_prev * wt_vel + wd_train / sigd) / wt_vel
            revap = 0.5 * (-b6 + np.sqrt(b6 * b6 + 4.0 * c6))

            evap[i] = sigs * afac * revap
            water[i] = revap * revap

            if i > 0:
                dhdp = np.maximum((h[i] - h[i - 1]) / (p[i - 1] - p[i]), 10.0)
                mp[i] = 100.0 * ginv * lv[i] * sigd * evap[i] / dhdp

        precip = wt_vel * sigd * water[0] * 3600.0 * 24000.0 / (self.ROWL * self.G)

        ft = np.zeros_like(t)
        fq = np.zeros_like(q)
        fu = np.zeros_like(state["eastward_wind"])
        fv = np.zeros_like(state["northward_wind"])

        for i in range(nlev):
            dp = ph[i] - ph[i + 1]
            dpinv = 0.01 / dp

            ft[i] = self.G * dpinv * m_flux[i] * (tp[i] - t[i])
            fq[i] = self.G * dpinv * m_flux[i] * (qp[i] - q[i])

            ft[i] += self.G * dpinv * mp[i] * (t[i] - t[i])
            fq[i] += self.G * dpinv * mp[i] * (qp[i] - q[i])

        for c in range(ncol):
            column_enthalpy_tendency = np.sum(
                (self.CPD * ft[:, c] + lv[:, c] * fq[:, c]) * (ph[:-1, c] - ph[1:, c])
            )
            correction = column_enthalpy_tendency / np.sum(ph[:-1, c] - ph[1:, c])
            ft[:, c] -= correction / self.CPD

        return {
            "air_temperature": ft,
            "specific_humidity": fq,
            "eastward_wind": fu,
            "northward_wind": fv,
        }, {
            "convective_precipitation_rate": precip,
        }

        return {
            "air_temperature": ft,
            "specific_humidity": fq,
            "eastward_wind": fu,
            "northward_wind": fv,
        }, {
            "convective_precipitation_rate": precip,
        }

    def _tlift(self, state, nk, icb, tv, gz):
        """
        Calculates lifted parcel properties using an iterative Newton-Raphson solver.
        """
        nlev, ncol = state["air_temperature"].shape
        p = state["air_pressure"]
        t = state["air_temperature"]
        q = state["specific_humidity"]
        qs = state.get("saturation_specific_humidity", q)  # Fallback to environment q

        tp = np.zeros_like(t)
        qp = np.zeros_like(q)
        clw = np.zeros_like(q)
        tvp = np.zeros_like(t)

        eps = self.RD / self.RV
        epsi = 1.0 / eps
        cpvmcl = self.CL - self.CPV

        # Vectorized scalar properties at parcel origin (NK)
        rows = np.arange(ncol)
        t_nk = t[nk, rows]
        q_nk = q[nk, rows]
        gz_nk = gz[nk, rows]

        # Static Energy at source (Eq. from Fortran)
        ah0 = (
            (self.CPD * (1.0 - q_nk) + self.CL * q_nk) * t_nk
            + q_nk * (self.LV0 - cpvmcl * (t_nk - 273.15))
            + gz_nk
        )

        cpp = self.CPD * (1.0 - q_nk) + q_nk * self.CPV
        cpinv = 1.0 / cpp

        # Loop over all levels
        for i in range(nlev):
            # Mask for levels below ICB (Dry Adiabatic)
            # Since i is fixed, we check which columns have icb > i

            # Dry Adiabatic (i < icb)
            is_dry = i < icb

            if np.any(is_dry):
                # TPK(I) = T(NK) - (GZ(I)-GZ(NK))*CPINV
                tp_dry = t_nk - (gz[i] - gz_nk) * cpinv
                # TVP(I) = TPK(I)*(1.+Q(NK)*EPSI)
                tvp_dry = tp_dry * (1.0 + q_nk * epsi)

                # Apply to dry columns
                tp[i] = np.where(is_dry, tp_dry, tp[i])
                qp[i] = np.where(is_dry, q_nk, qp[i])  # Conserved q
                clw[i] = np.where(is_dry, 0.0, clw[i])
                tvp[i] = np.where(is_dry, tvp_dry, tvp[i])

            # Moist Adiabatic (i >= icb)
            is_moist = ~is_dry
            if np.any(is_moist):
                # Trial values
                tg = t[i]
                qg = q[i]  # Initial guess is environment

                # Iterate 2 times
                for _ in range(2):
                    alv = self.LV0 - cpvmcl * (tg - 273.15)
                    s = 1.0 / (self.CPD + alv * alv * qg / (self.RV * tg * tg))
                    ahg = (
                        self.CPD * tg
                        + (self.CL - self.CPD) * q_nk * tg
                        + alv * qg
                        + gz[i]
                    )

                    tg = tg + s * (ah0 - ahg)
                    tg = np.maximum(tg, 35.0)

                    tc = tg - 273.15
                    denom = 243.5 + tc
                    # Saturation vapor pressure
                    es = np.where(
                        tc >= 0.0,
                        6.112 * np.exp(17.67 * tc / denom),
                        np.exp(23.33086 - 6111.72784 / tg + 0.15215 * np.log(tg)),
                    )

                    qg = eps * es / (p[i] - es * (1.0 - eps))

                # Final values for moist
                tp_moist = (
                    ah0 - (self.CL - self.CPD) * q_nk * tg - gz[i] - alv * qg
                ) / self.CPD
                qp_moist = qg
                clw_moist = np.maximum(0.0, q_nk - qg)
                rg = qg / (1.0 - q_nk)
                tvp_moist = tp_moist * (1.0 + rg * epsi)

                tp[i] = np.where(is_moist, tp_moist, tp[i])
                qp[i] = np.where(is_moist, qp_moist, qp[i])
                clw[i] = np.where(is_moist, clw_moist, clw[i])
                tvp[i] = np.where(is_moist, tvp_moist, tvp[i])

        return tp, qp, clw, tvp

    def _find_parcel_origin(self, hm, minorig, nlev, ncol):
        """Find the level of maximum moist static energy (NK) below the level of minimum MSE."""
        nk = np.zeros(ncol, dtype=int)

        # 1. Find level of minimum MSE (IHMIN)
        # Search from minorig upwards to find local minimum
        ahmin = np.full(ncol, 1.0e12)
        ihmin = np.full(ncol, nlev - 1, dtype=int)

        for i in range(minorig, nlev):
            is_new_min = (hm[i] < ahmin) & (hm[i] < hm[i - 1])
            ahmin = np.where(is_new_min, hm[i], ahmin)
            ihmin = np.where(is_new_min, i, ihmin)

        # 2. Find Max MSE below IHMIN
        ahmax = np.zeros(ncol)
        nk[:] = minorig

        for i in range(minorig, nlev):
            is_below_min = i <= ihmin
            is_new_max = (hm[i] > ahmax) & is_below_min

            ahmax = np.where(is_new_max, hm[i], ahmax)
            nk = np.where(is_new_max, i, nk)

        return nk

    def _calculate_lcl(self, p, t, q, qs, nk, ncol):
        """Calculate LCL using Bolton's approximation."""
        rows = np.arange(ncol)
        t_nk = t[nk, rows]
        q_nk = q[nk, rows]
        qs_nk = qs[nk, rows]
        p_nk = p[nk, rows]

        rh = q_nk / qs_nk
        chi = t_nk / (1669.0 - 122.0 * rh - t_nk)
        plcl = p_nk * (rh**chi)

        return plcl

    def _find_icb(self, p, plcl, nk, nlev, ncol):
        """Find the first level above LCL (ICB)."""
        icb = np.full(ncol, nlev - 1, dtype=int)

        for i in range(nlev):
            is_above_origin = i > nk
            is_above_lcl = p[i] < plcl

            # Update to minimum index satisfying condition
            update_mask = is_above_origin & is_above_lcl
            icb = np.where(update_mask, np.minimum(icb, i), icb)

        return icb
