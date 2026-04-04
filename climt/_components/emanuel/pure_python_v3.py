# -*- coding: utf-8 -*-
from typing import NamedTuple

import numpy as np
from sympl import (
    ImplicitTendencyComponent,
    # get_constant,
    # initialize_numpy_arrays_with_properties,
)

from ..._core import ensure_contiguous_state
from ..._core.backend import jit_compile, prange

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(x, **kwargs):
        return x

try:
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp

    HAS_JAX = True
except ImportError:
    HAS_JAX = False


class EmanuelParams(NamedTuple):
    IPBL: int
    MINORIG: int
    ELCRIT: float
    TLCRIT: float
    ENTP: float
    SIGD: float
    SIGS: float
    OMTRAIN: float
    OMTSNOW: float
    COEFFR: float
    COEFFS: float
    CU: float
    BETA: float
    DTMAX: float
    ALPHA: float
    DAMP: float
    CPD: float
    CPV: float
    CL: float
    RV: float
    RD: float
    LV0: float
    G: float
    ROWL: float
    DELT0: float


class EmanuelConvectionPythonV3(ImplicitTendencyComponent):
    input_properties = {
        "air_temperature": {"dims": ["*", "mid_levels"], "units": "degK"},
        "specific_humidity": {"dims": ["*", "mid_levels"], "units": "kg/kg"},
        "eastward_wind": {"dims": ["*", "mid_levels"], "units": "m s^-1"},
        "northward_wind": {"dims": ["*", "mid_levels"], "units": "m s^-1"},
        "air_pressure": {"dims": ["*", "mid_levels"], "units": "mbar"},
        "air_pressure_on_interface_levels": {
            "dims": ["*", "interface_levels"],
            "units": "mbar",
        },
        "cloud_base_mass_flux": {"dims": ["*"], "units": "kg m^-2 s^-1"},
    }

    diagnostic_properties = {
        "convective_state": {
            "dims": ["*"],
            "units": "dimensionless",
            "dtype": np.int32,
        },
        "convective_precipitation_rate": {"dims": ["*"], "units": "mm day^-1"},
        "convective_downdraft_velocity_scale": {"dims": ["*"], "units": "m s^-1"},
        "convective_downdraft_temperature_scale": {"dims": ["*"], "units": "degK"},
        "convective_downdraft_specific_humidity_scale": {
            "dims": ["*"],
            "units": "kg/kg",
        },
        "cloud_base_mass_flux": {"dims": ["*"], "units": "kg m^-2 s^-1"},
        "atmosphere_convective_available_potential_energy": {
            "dims": ["*"],
            "units": "J kg^-1",
        },
        "air_temperature_tendency_from_convection": {
            "dims": ["*", "mid_levels"],
            "units": "degK day^-1",
        },
    }

    tendency_properties = {
        "air_temperature": {"units": "degK s^-1"},
        "specific_humidity": {"units": "kg/kg s^-1"},
        "eastward_wind": {"units": "m s^-2"},
        "northward_wind": {"units": "m s^-2"},
    }

    def __init__(self, **kwargs):
        self.IPBL = 0
        self.MINORIG = 1
        self.ELCRIT = 0.0011
        self.TLCRIT = -55.0
        self.ENTP = 1.5
        self.SIGD = 0.05
        self.SIGS = 0.12
        self.OMTRAIN = 50.0
        self.OMTSNOW = 5.5
        self.COEFFR = 1.0
        self.COEFFS = 0.8
        self.CU = 0.7
        self.BETA = 10.0
        self.DTMAX = 0.9
        self.ALPHA = 0.1
        self.DAMP = 0.1
        self.CPD = 1005.7
        self.CPV = 1870.0
        self.CL = 2500.0
        self.RV = 461.5
        self.RD = 287.04
        self.LV0 = 2.501e6
        self.G = 9.8
        self.ROWL = 1000.0
        self.DELT0 = 300.0
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
        self._params = EmanuelParams(
            IPBL=int(self.IPBL),
            MINORIG=int(self.MINORIG),
            ELCRIT=float(self.ELCRIT),
            TLCRIT=float(self.TLCRIT),
            ENTP=float(self.ENTP),
            SIGD=float(self.SIGD),
            SIGS=float(self.SIGS),
            OMTRAIN=float(self.OMTRAIN),
            OMTSNOW=float(self.OMTSNOW),
            COEFFR=float(self.COEFFR),
            COEFFS=float(self.COEFFS),
            CU=float(self.CU),
            BETA=float(self.BETA),
            DTMAX=float(self.DTMAX),
            ALPHA=float(self.ALPHA),
            DAMP=float(self.DAMP),
            CPD=float(self.CPD),
            CPV=float(self.CPV),
            CL=float(self.CL),
            RV=float(self.RV),
            RD=float(self.RD),
            LV0=float(self.LV0),
            G=float(self.G),
            ROWL=float(self.ROWL),
            DELT0=float(self.DELT0),
        )
        super(EmanuelConvectionPythonV3, self).__init__(**kwargs)

    @ensure_contiguous_state
    def array_call(self, state, timestep):
        t = state["air_temperature"]
        q = state["specific_humidity"]
        u = state["eastward_wind"]
        v = state["northward_wind"]
        p = state["air_pressure"]
        ph = state["air_pressure_on_interface_levels"]
        transposed = False
        if t.shape[0] != ph.shape[0] - 1:
            t = t.T
            q = q.T
            u = u.T
            v = v.T
            p = p.T
            ph = ph.T
            transposed = True
        nlev, ncol = t.shape
        delt = timestep.total_seconds()
        _use_jax = HAS_JAX and isinstance(t, jax.Array)
        if _use_jax:
            qs = _bolton_q_sat_jax(t, p * 100, self.RD, self.RV)
            cbmf_val = state.get("cloud_base_mass_flux", None)
            cbmf = jnp.zeros(ncol) if cbmf_val is None else jnp.array(cbmf_val)
            results = _jax_vectorized_convect(
                t, q, qs, u, v, p, ph,
                nlev, nlev - 3, 0, delt, cbmf, self._params,
            )
        else:
            from climt._core import bolton_q_sat
            qs = bolton_q_sat(t, p * 100, self.RD, self.RV)
            cbmf = state.get("cloud_base_mass_flux", np.zeros(ncol)).copy()
            ntra = 0
            tra = np.zeros((nlev, 1))
            tra_vector = np.broadcast_to(tra[:, :, np.newaxis], (nlev, 1, ncol))
            results = _numpy_vectorized_convect(
                t, q, qs, u, v, p, ph,
                nlev, nlev - 3, ntra, delt, cbmf, tra_vector, self._params,
            )
        ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results
        if transposed:
            ft = ft.T
            fq = fq.T
            fu = fu.T
            fv = fv.T
        tendencies = {
            "air_temperature": ft,
            "specific_humidity": fq,
            "eastward_wind": fu,
            "northward_wind": fv,
        }
        diagnostics = {
            "convective_state": iflag,
            "convective_precipitation_rate": precip,
            "convective_downdraft_velocity_scale": wd,
            "convective_downdraft_temperature_scale": tprime,
            "convective_downdraft_specific_humidity_scale": qprime,
            "cloud_base_mass_flux": cbmf_new,
            "atmosphere_convective_available_potential_energy": outcape,
        }
        diagnostics["air_temperature_tendency_from_convection"] = ft * 86400.0
        return tendencies, diagnostics


@njit
def _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params):
    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS
    AH0 = (
        (params.CPD * (1.0 - Q[NK]) + params.CL * Q[NK]) * T[NK]
        + Q[NK] * (params.LV0 - CPVMCL * (T[NK] - 273.15))
        + GZ[NK]
    )
    CPP = params.CPD * (1.0 - Q[NK]) + Q[NK] * params.CPV
    CPINV = 1.0 / CPP
    if KK == 1:
        for i in range(ICB):
            CLW[i] = 0.0
        for i in range(NK, ICB):
            TPK[i] = T[NK] - (GZ[i] - GZ[NK]) * CPINV
            TVP[i] = TPK[i] * (1.0 + Q[NK] * EPSI)
    NST = NL if KK == 2 else ICB
    NSB = ICB + 1 if KK == 2 else ICB
    for i in range(NSB, NST + 1):
        TG = T[i]
        QG = QS[i]
        ALV = params.LV0 - CPVMCL * (T[i] - 273.15)
        for j in range(2):
            S = 1.0 / (params.CPD + ALV * ALV * QG / (params.RV * T[i] * T[i]))
            AHG = (
                params.CPD * TG
                + (params.CL - params.CPD) * Q[NK] * T[i]
                + ALV * QG
                + GZ[i]
            )
            TG = max(TG + S * (AH0 - AHG), 35.0)
            TC = TG - 273.15
            DENOM = 243.5 + TC
            ES = (
                6.112 * np.exp(17.67 * TC / DENOM)
                if TC >= 0.0
                else np.exp(23.33086 - 6111.72784 / TG + 0.15215 * np.log(TG))
            )
            QG = EPS * ES / (P[i] - ES * (1.0 - EPS))
        TPK[i] = (
            AH0 - (params.CL - params.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG
        ) / params.CPD
        CLW[i] = max(0.0, Q[NK] - QG)
        TVP[i] = TPK[i] * (1.0 + (QG / (1.0 - Q[NK])) * EPSI)
    return TVP, TPK, CLW


@njit
def _convect_functional_np(
    T_in,
    Q_in,
    QS_in,
    U_in,
    V_in,
    P_in,
    PH_in,
    ND,
    NL,
    NTRA,
    DELT,
    CBMF_in,
    TRA_in,
    params,
):
    T = T_in.copy()
    Q = Q_in.copy()
    QS = QS_in.copy()
    U = U_in.copy()
    V = V_in.copy()
    P = P_in.copy()
    PH = PH_in.copy()
    TRA = TRA_in.copy()
    CBMF = CBMF_in
    FT = np.zeros(ND)
    FQ = np.zeros(ND)
    FU = np.zeros(ND)
    FV = np.zeros(ND)
    FTRA = np.zeros((ND, max(1, NTRA)))
    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS
    GINV = 1.0 / params.G
    DELTI = 1.0 / DELT
    PRECIP = 0.0
    WD = 0.0
    TPRIME = 0.0
    QPRIME = 0.0
    IFLAG = 0
    OUTCAPE = 0.0
    TH = np.zeros(NL + 1)
    for i in range(NL + 1):
        RDCP = (params.RD * (1.0 - Q[i]) + Q[i] * params.RV) / (
            params.CPD * (1.0 - Q[i]) + Q[i] * params.CPV
        )
        TH[i] = T[i] * (1000.0 / P[i]) ** RDCP
    GZ = np.zeros(ND + 1)
    CPN = np.zeros(ND + 1)
    H = np.zeros(ND + 1)
    LV = np.zeros(ND + 1)
    HM = np.zeros(ND + 1)
    TV = np.zeros(ND + 1)
    GZ[0] = 0.0
    CPN[0] = params.CPD * (1.0 - Q[0]) + Q[0] * params.CPV
    H[0] = T[0] * CPN[0]
    LV[0] = params.LV0 - CPVMCL * (T[0] - 273.15)
    HM[0] = LV[0] * Q[0]
    TV[0] = T[0] * (1.0 + Q[0] * EPSI - Q[0])
    AHMIN = 1.0e12
    IHMIN = NL
    for i in range(1, NL + 1):
        TVX = T[i] * (1.0 + Q[i] * EPSI - Q[i])
        TVY = T[i - 1] * (1.0 + Q[i - 1] * EPSI - Q[i - 1])
        GZ[i] = GZ[i - 1] + 0.5 * params.RD * (TVX + TVY) * (P[i - 1] - P[i]) / PH[i]
        CPN[i] = params.CPD * (1.0 - Q[i]) + params.CPV * Q[i]
        H[i] = T[i] * CPN[i] + GZ[i]
        LV[i] = params.LV0 - CPVMCL * (T[i] - 273.15)
        HM[i] = (
            (params.CPD * (1.0 - Q[i]) + params.CL * Q[i]) * (T[i] - T[0])
            + LV[i] * Q[i]
            + GZ[i]
        )
        TV[i] = T[i] * (1.0 + Q[i] * EPSI - Q[i])
        if (i + 1) >= params.MINORIG and HM[i] < AHMIN and HM[i] < HM[i - 1]:
            AHMIN = HM[i]
            IHMIN = i
    IHMIN = min(IHMIN, NL - 1)
    AHMAX = 0.0
    NK = 0
    for i in range(params.MINORIG - 1, IHMIN + 1):
        if HM[i] > AHMAX:
            NK = i
            AHMAX = HM[i]
    if T[NK] < 250.0 or Q[NK] <= 0.0 or IHMIN == (NL - 1):
        CBMF = 0.0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
    RH = Q[NK] / QS[NK]
    CHI = T[NK] / (1669.0 - 122.0 * RH - T[NK])
    PLCL = P[NK] * (RH**CHI)
    if PLCL < 200.0 or PLCL >= 2000.0:
        IFLAG = 2
        CBMF = 0.0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
    ICB = NL - 1
    for i in range(NK + 1, NL + 1):
        if P[i] < PLCL:
            ICB = min(ICB, i)
    if ICB >= (NL - 1):
        IFLAG = 3
        CBMF = 0.0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
    TVP = np.zeros(ND)
    TP = np.zeros(ND)
    CLW = np.zeros(ND)
    TVP, TP, CLW = _tlift_functional_np(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params
    )
    for i in range(NK, ICB + 1):
        TVP[i] -= TP[i] * Q[NK]
    if CBMF == 0.0 and TVP[ICB] <= (TV[ICB] - params.DTMAX):
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
    if IFLAG != 4:
        IFLAG = 1
    TVP, TP, CLW = _tlift_functional_np(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params
    )
    EP = np.zeros(ND)
    SIGP = np.zeros(ND)
    for i in range(NK + 1):
        EP[i] = 0.0
        SIGP[i] = params.SIGS
    for i in range(NK + 1, NL + 1):
        TCA = TP[i] - 273.15
        ELACRIT = (
            params.ELCRIT if TCA >= 0.0 else params.ELCRIT * (1.0 - TCA / params.TLCRIT)
        )
        ELACRIT = max(ELACRIT, 0.0)
        EPMAX = 0.999
        EP[i] = EPMAX * (1.0 - ELACRIT / max(CLW[i], 1.0e-8))
        EP[i] = max(min(EP[i], EPMAX), 0.0)
        SIGP[i] = params.SIGS
    for i in range(ICB + 1, NL + 1):
        TVP[i] -= TP[i] * Q[NK]
    HP = H.copy()
    NENT = np.zeros(ND + 1, dtype=np.int32)
    WATER = np.zeros(ND + 1)
    EVAP = np.zeros(ND + 1)
    WT = np.full(ND + 1, params.OMTSNOW)
    MP = np.zeros(ND + 1)
    M = np.zeros(ND + 1)
    LVCP = np.zeros(ND + 1)
    for i in range(NL + 1):
        LVCP[i] = LV[i] / CPN[i]
    QENT = np.zeros((ND + 1, ND + 1))
    ELIJ = np.zeros((ND + 1, ND + 1))
    MENT = np.zeros((ND + 1, ND + 1))
    SIJ = np.zeros((ND + 1, ND + 1))
    UENT = np.zeros((ND + 1, ND + 1))
    VENT = np.zeros((ND + 1, ND + 1))
    TRAENT = np.zeros((ND + 1, ND + 1, max(1, NTRA)))
    for i in range(NL + 1):
        for j in range(NL + 1):
            QENT[i, j] = Q[j]
            UENT[i, j] = U[j]
            VENT[i, j] = V[j]
            for k in range(NTRA):
                TRAENT[i, j, k] = TRA[j, k]
    QP = np.zeros(ND + 1)
    UP = np.zeros(ND + 1)
    VP = np.zeros(ND + 1)
    TRAP = np.zeros((ND + 1, max(1, NTRA)))
    QP[0] = Q[0]
    UP[0] = U[0]
    VP[0] = V[0]
    for i in range(NTRA):
        TRAP[0, i] = TRA[0, i]
    for i in range(1, NL + 1):
        QP[i] = Q[i - 1]
        UP[i] = U[i - 1]
        VP[i] = V[i - 1]
        for j in range(NTRA):
            TRAP[i, j] = TRA[i - 1, j]
    CAPE = 0.0
    CAPEM = 0.0
    INB = ICB + 1
    INB1 = INB
    BYP = 0.0
    for i in range(ICB + 1, NL):
        BY = (TVP[i] - TV[i]) * (PH[i] - PH[i + 1]) / P[i]
        CAPE += BY
        if BY >= 0.0:
            INB1 = i + 1
        if CAPE > 0.0:
            INB = i + 1
            BYP = (TVP[i + 1] - TV[i + 1]) * (PH[i + 1] - PH[i + 2]) / P[i + 1]
            CAPEM = CAPE
    INB = max(INB, INB1)
    CAPE = CAPEM + BYP
    DEFRAC = max(CAPEM - CAPE, 0.001)
    FRAC = min(max(-CAPE / DEFRAC, 0.0), 1.0)
    OUTCAPE = CAPE
    for i in range(ICB, INB + 1):
        HP[i] = H[NK] + (LV[i] + (params.CPD - params.CPV) * T[i]) * EP[i] * CLW[i]
    TVPPLCL = TVP[ICB - 1] - params.RD * TVP[ICB - 1] * (P[ICB - 1] - PLCL) / (
        CPN[ICB - 1] * P[ICB - 1]
    )
    TVAPLCL = TV[ICB] + (TVP[ICB] - TVP[ICB + 1]) * (PLCL - P[ICB]) / (
        P[ICB] - P[ICB + 1]
    )
    DTPBL = 0.0
    for i in range(NK, ICB):
        DTPBL += (TVP[i] - TV[i]) * (PH[i] - PH[i + 1])
    DTPBL /= PH[NK] - PH[ICB]
    DTMA = TVPPLCL - TVAPLCL + params.DTMAX + DTPBL
    CBMFOLD = CBMF
    DAMPS = params.DAMP * DELT / params.DELT0
    CBMF = (1.0 - DAMPS) * CBMF + 0.1 * params.ALPHA * DTMA
    CBMF = max(CBMF, 0.0)
    if CBMF == 0.0 and CBMFOLD == 0.0:
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
    M[ICB] = 0.0
    DBOSUM = 0.0
    for i in range(ICB + 1, INB + 1):
        k = min(i, INB1)
        DBO = abs(TV[k] - TVP[k]) + params.ENTP * 0.02 * (PH[k] - PH[k + 1])
        DBOSUM += DBO
        M[i] = CBMF * DBO
    if DBOSUM > 0:
        for i in range(ICB + 1, INB + 1):
            M[i] /= DBOSUM
    for i in range(ICB + 1, INB + 1):
        QTI = Q[NK] - EP[i] * CLW[i]
        for j in range(ICB, INB + 1):
            BF2 = 1.0 + LV[j] * LV[j] * QS[j] / (params.RV * T[j] * T[j] * params.CPD)
            ANUM = H[j] - HP[i] + (params.CPV - params.CPD) * T[j] * (QTI - Q[j])
            DENOM = H[i] - HP[i] + (params.CPD - params.CPV) * (Q[i] - QTI) * T[j]
            DEI = DENOM
            if abs(DEI) < 0.01:
                DEI = 0.01
            SIJ[i, j] = ANUM / DEI
            SIJ[i, i] = 1.0
            ALTEM = (SIJ[i, j] * Q[i] + (1.0 - SIJ[i, j]) * QTI - QS[j]) / BF2
            CWAT = CLW[j] * (1.0 - EP[j])
            if (SIJ[i, j] < 0.0 or SIJ[i, j] > 1.0 or ALTEM > CWAT) and j > i:
                ANUM -= LV[j] * (QTI - QS[j] - CWAT * BF2)
                DENOM += LV[j] * (Q[i] - QTI)
                if abs(DENOM) < 0.01:
                    DENOM = 0.01
                SIJ[i, j] = ANUM / DENOM
                ALTEM = (
                    SIJ[i, j] * Q[i]
                    + (1.0 - SIJ[i, j]) * QTI
                    - QS[j]
                    - (BF2 - 1.0) * CWAT
                )
            if 0.0 < SIJ[i, j] < 0.9:
                QENT[i, j] = SIJ[i, j] * Q[i] + (1.0 - SIJ[i, j]) * QTI
                UENT[i, j] = SIJ[i, j] * U[i] + (1.0 - SIJ[i, j]) * U[NK]
                VENT[i, j] = SIJ[i, j] * V[i] + (1.0 - SIJ[i, j]) * V[NK]
                for k in range(NTRA):
                    TRAENT[i, j, k] = (
                        SIJ[i, j] * TRA[i, k] + (1.0 - SIJ[i, j]) * TRA[NK, k]
                    )
                ELIJ[i, j] = max(0.0, ALTEM)
                MENT[i, j] = M[i] / (1.0 - SIJ[i, j])
                NENT[i] += 1
            SIJ[i, j] = min(max(SIJ[i, j], 0.0), 1.0)
        if NENT[i] == 0:
            MENT[i, i] = M[i]
            QENT[i, i] = Q[NK] - EP[i] * CLW[i]
            UENT[i, i] = U[NK]
            VENT[i, i] = V[NK]
            for k in range(NTRA):
                TRAENT[i, i, k] = TRA[NK, k]
            ELIJ[i, i] = CLW[i]
            SIJ[i, i] = 1.0
    SIJ[INB, INB] = 1.0
    for i in range(ICB + 1, INB + 1):
        if NENT[i] != 0:
            QP1 = Q[NK] - EP[i] * CLW[i]
            ANUM = H[i] - HP[i] - LV[i] * (QP1 - QS[i])
            DENOM = H[i] - HP[i] + LV[i] * (Q[i] - QP1)
            if abs(DENOM) < 0.01:
                DENOM = 0.01
            SCRIT = ANUM / DENOM
            ALT = QP1 - QS[i] + SCRIT * (Q[i] - QP1)
            if ALT < 0.0:
                SCRIT = 1.0
            SCRIT = max(SCRIT, 0.0)
            ASIJ = 0.0
            SMIN = 1.0
            for j in range(ICB, INB + 1):
                if 0.0 < SIJ[i, j] < 0.9:
                    if j > i:
                        SMID = min(SIJ[i, j], SCRIT)
                        SJMAX = SMID
                        SJMIN = SMID
                        if SMID < SMIN and SIJ[i, j + 1] < SMID:
                            SMIN = SMID
                            SJMAX = min(min(SIJ[i, j + 1], SIJ[i, j]), SCRIT)
                            SJMIN = min(max(SIJ[i, j - 1], SIJ[i, j]), SCRIT)
                    else:
                        SJMAX = max(SIJ[i, j + 1], SCRIT)
                        SMID = max(SIJ[i, j], SCRIT)
                        SJMIN = max(SIJ[i, j - 1] if j > 0 else 0.0, SCRIT)
                    ASIJ += (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (
                        PH[j] - PH[j + 1]
                    )
                    MENT[i, j] *= (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (
                        PH[j] - PH[j + 1]
                    )
            ASIJ = max(1.0e-21, ASIJ)
            for j in range(ICB, INB + 1):
                MENT[i, j] /= ASIJ
            s = 0.0
            for j in range(ICB, INB + 1):
                s += MENT[i, j]
            if s < 1.0e-18:
                NENT[i] = 0
                MENT[i, i] = M[i]
                QENT[i, i] = Q[NK] - EP[i] * CLW[i]
                UENT[i, i] = U[NK]
                VENT[i, i] = V[NK]
                for k in range(NTRA):
                    TRAENT[i, i, k] = TRA[NK, k]
                ELIJ[i, i] = CLW[i]
                SIJ[i, i] = 1.0
    if EP[INB] >= 0.0001:
        JTT = 0
        for i in range(INB, -1, -1):
            WDTRAIN = params.G * EP[i] * M[i] * CLW[i]
            if i > 0:
                for j in range(i):
                    WDTRAIN += (
                        params.G
                        * max(0.0, ELIJ[j, i] - (1.0 - EP[i]) * CLW[i])
                        * MENT[j, i]
                    )
            COEFF = params.COEFFR if T[i] > 273.0 else params.COEFFS
            WT[i] = params.OMTRAIN if T[i] > 273.0 else params.OMTSNOW
            AFAC = max(
                COEFF
                * PH[i]
                * (QS[i] - 0.5 * (Q[i] + QP[i + 1]))
                / (1.0e4 + 2.0e3 * PH[i] * QS[i]),
                0.0,
            )
            SIGT = min(max(SIGP[i], 0.0), 1.0)
            B6 = 100.0 * (PH[i] - PH[i + 1]) * SIGT * AFAC / WT[i]
            C6 = (WATER[i + 1] * WT[i + 1] + WDTRAIN / params.SIGD) / WT[i]
            REVAP = 0.5 * (-B6 + np.sqrt(B6 * B6 + 4.0 * C6))
            EVAP[i] = SIGT * AFAC * REVAP
            WATER[i] = REVAP * REVAP
            if i > 0:
                DHDP = max((H[i] - H[i - 1]) / (P[i - 1] - P[i]), 10.0)
                MP[i] = 100.0 * GINV * LV[i] * params.SIGD * EVAP[i] / DHDP
                FAC = 20.0 / (PH[i - 1] - PH[i])
                MP[i] = (FAC * MP[i + 1] + MP[i]) / (1.0 + FAC)
                if P[i] > (0.949 * P[0]):
                    JTT = max(JTT, i)
                    MP[i] = MP[JTT] * (P[0] - P[i]) / (P[0] - P[JTT])
            if i != INB:
                QSTM = QS[0] if i == 0 else QS[i - 1]
                if MP[i] > MP[i + 1]:
                    RAT = MP[i + 1] / MP[i]
                    QP[i] = (
                        QP[i + 1] * RAT
                        + Q[i] * (1.0 - RAT)
                        + 100.0
                        * GINV
                        * params.SIGD
                        * (PH[i] - PH[i + 1])
                        * (EVAP[i] / MP[i])
                    )
                    UP[i] = UP[i + 1] * RAT + U[i] * (1.0 - RAT)
                    VP[i] = VP[i + 1] * RAT + V[i] * (1.0 - RAT)
                    for k in range(NTRA):
                        TRAP[i, k] = TRAP[i + 1, k] * RAT + TRAP[i, k] * (1.0 - RAT)
                elif MP[i + 1] > 0.0:
                    QP[i] = (
                        GZ[i + 1]
                        - GZ[i]
                        + QP[i + 1] * (LV[i + 1] + T[i + 1] * (params.CL - params.CPD))
                        + params.CPD * (T[i + 1] - T[i])
                    ) / (LV[i] + T[i] * (params.CL - params.CPD))
                    UP[i] = UP[i + 1]
                    VP[i] = VP[i + 1]
                    for k in range(NTRA):
                        TRAP[i, k] = TRAP[i + 1, k]
                QP[i] = max(min(QP[i], QSTM), 0.0)
        PRECIP += (
            WT[0] * params.SIGD * WATER[0] * 3600.0 * 24000.0 / (params.ROWL * params.G)
        )
    WD = params.BETA * abs(MP[ICB]) * 0.01 * params.RD * T[ICB] / (params.SIGD * P[ICB])
    QPRIME = 0.5 * (QP[0] - Q[0])
    TPRIME = params.LV0 * QPRIME / params.CPD
    DPINV = 0.01 / (PH[0] - PH[1])
    AM = 0.0
    if NK == 0:
        for ii in range(1, INB + 1):
            AM += M[ii]
    if (2.0 * params.G * DPINV * AM) >= DELTI:
        IFLAG = 4
    FT[0] += (
        params.G * DPINV * AM * (T[1] - T[0] + (GZ[1] - GZ[0]) / CPN[0])
        - LVCP[0] * params.SIGD * EVAP[0]
        + params.SIGD
        * WT[1]
        * (params.CL - params.CPD)
        * WATER[1]
        * (T[1] - T[0])
        * DPINV
        / CPN[0]
    )
    FQ[0] += (
        params.G * MP[1] * (QP[1] - Q[0]) * DPINV
        + params.SIGD * EVAP[0]
        + params.G * AM * (Q[1] - Q[0]) * DPINV
    )
    FU[0] += params.G * DPINV * (MP[1] * (UP[1] - U[0]) + AM * (U[1] - U[0]))
    FV[0] += params.G * DPINV * (MP[1] * (VP[1] - V[0]) + AM * (V[1] - V[0]))
    for j in range(NTRA):
        FTRA[0, j] += (
            params.G
            * DPINV
            * (MP[1] * (TRAP[1, j] - TRA[0, j]) + AM * (TRA[1, j] - TRA[0, j]))
        )
    for j in range(1, INB + 1):
        FQ[0] += params.G * DPINV * MENT[j, 0] * (QENT[j, 0] - Q[0])
        FU[0] += params.G * DPINV * MENT[j, 0] * (UENT[j, 0] - U[0])
        FV[0] += params.G * DPINV * MENT[j, 0] * (VENT[j, 0] - V[0])
        for k in range(NTRA):
            FTRA[0, k] += params.G * DPINV * MENT[j, 0] * (TRAENT[j, 0, k] - TRA[0, k])
    for i in range(1, INB + 1):
        DPINV = 0.01 / (PH[i] - PH[i + 1])
        CPINV = 1.0 / CPN[i]
        AMP1 = 0.0
        if i >= NK:
            for ii in range(i + 1, INB + 2):
                AMP1 += M[ii]
        for k in range(i + 1):
            for jj in range(i + 1, INB + 2):
                AMP1 += MENT[k, jj]
        if (2.0 * params.G * DPINV * AMP1) >= DELTI:
            IFLAG = 4
        AD = 0.0
        for k in range(i):
            for jj in range(i, INB + 1):
                AD += MENT[jj, k]
        FT[i] += (
            params.G
            * DPINV
            * (
                AMP1 * (T[i + 1] - T[i] + (GZ[i + 1] - GZ[i]) * CPINV)
                - AD * (T[i] - T[i - 1] + (GZ[i] - GZ[i - 1]) * CPINV)
            )
            - params.SIGD * LVCP[i] * EVAP[i]
        )
        FT[i] += (
            params.G
            * DPINV
            * MENT[i, i]
            * (HP[i] - H[i] + T[i] * (params.CPV - params.CPD) * (Q[i] - QENT[i, i]))
            * CPINV
        )
        FT[i] += (
            params.SIGD
            * WT[i + 1]
            * (params.CL - params.CPD)
            * WATER[i + 1]
            * (T[i + 1] - T[i])
            * DPINV
            * CPINV
        )
        FQ[i] += params.G * DPINV * (AMP1 * (Q[i + 1] - Q[i]) - AD * (Q[i] - Q[i - 1]))
        FU[i] += params.G * DPINV * (AMP1 * (U[i + 1] - U[i]) - AD * (U[i] - U[i - 1]))
        FV[i] += params.G * DPINV * (AMP1 * (V[i + 1] - V[i]) - AD * (V[i] - V[i - 1]))
        for k in range(NTRA):
            FTRA[i, k] += (
                params.G
                * DPINV
                * (
                    AMP1 * (TRA[i + 1, k] - TRA[i, k])
                    - AD * (TRA[i, k] - TRA[i - 1, k])
                )
            )
        for k in range(i):
            AWAT = max(ELIJ[k, i] - (1.0 - EP[i]) * CLW[i], 0.0)
            FQ[i] += params.G * DPINV * MENT[k, i] * (QENT[k, i] - AWAT - Q[i])
            FU[i] += params.G * DPINV * MENT[k, i] * (UENT[k, i] - U[i])
            FV[i] += params.G * DPINV * MENT[k, i] * (VENT[k, i] - V[i])
            for j in range(NTRA):
                FTRA[i, j] += (
                    params.G * DPINV * MENT[k, i] * (TRAENT[k, i, j] - TRA[i, j])
                )
        for k in range(i, INB + 1):
            FQ[i] += params.G * DPINV * MENT[k, i] * (QENT[k, i] - Q[i])
            FU[i] += params.G * DPINV * MENT[k, i] * (UENT[k, i] - U[i])
            FV[i] += params.G * DPINV * MENT[k, i] * (VENT[k, i] - V[i])
            for j in range(NTRA):
                FTRA[i, j] += (
                    params.G * DPINV * MENT[k, i] * (TRAENT[k, i, j] - TRA[i, j])
                )
        FQ[i] += (
            params.SIGD * EVAP[i]
            + params.G
            * (MP[i + 1] * (QP[i + 1] - Q[i]) - MP[i] * (QP[i] - Q[i - 1]))
            * DPINV
        )
        FU[i] += (
            params.G
            * (MP[i + 1] * (UP[i + 1] - U[i]) - MP[i] * (UP[i] - U[i - 1]))
            * DPINV
        )
        FV[i] += (
            params.G
            * (MP[i + 1] * (VP[i + 1] - V[i]) - MP[i] * (VP[i] - V[i - 1]))
            * DPINV
        )
        for j in range(NTRA):
            FTRA[i, j] += (
                params.G
                * DPINV
                * (
                    MP[i + 1] * (TRAP[i + 1, j] - TRA[i, j])
                    - MP[i] * (TRAP[i, j] - TRA[i - 1, j])
                )
            )
    FQOLD = FQ[INB]
    FQ[INB] *= 1.0 - FRAC
    FQ[INB - 1] += (
        FRAC
        * FQOLD
        * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
        * LV[INB]
        / LV[INB - 1]
    )
    FTOLD = FT[INB]
    FT[INB] *= 1.0 - FRAC
    FT[INB - 1] += (
        FRAC
        * FTOLD
        * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
        * CPN[INB]
        / CPN[INB - 1]
    )
    FUOLD = FU[INB]
    FU[INB] *= 1.0 - FRAC
    FU[INB - 1] += FRAC * FUOLD * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
    FVOLD = FV[INB]
    FV[INB] *= 1.0 - FRAC
    FV[INB - 1] += FRAC * FVOLD * ((PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB]))
    for k in range(NTRA):
        FTRAOLD = FTRA[INB, k]
        FTRA[INB, k] *= 1.0 - FRAC
        FTRA[INB - 1, k] += (
            FRAC * FTRAOLD * (PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB])
        )
    ENTS = 0.0
    UAV = 0.0
    VAV = 0.0
    DENOM2 = PH[0] - PH[INB + 1]
    for i in range(INB + 1):
        dp = PH[i] - PH[i + 1]
        ENTS += (CPN[i] * FT[i] + LV[i] * FQ[i]) * dp
        UAV += FU[i] * dp
        VAV += FV[i] * dp
    ENTS /= DENOM2
    UAV /= DENOM2
    VAV /= DENOM2
    for i in range(INB + 1):
        FT[i] -= ENTS / CPN[i]
        FU[i] = (1.0 - params.CU) * (FU[i] - UAV)
        FV[i] = (1.0 - params.CU) * (FV[i] - VAV)
    for k in range(NTRA):
        TRAAV = 0.0
        for i in range(INB + 1):
            TRAAV += FTRA[i, k] * (PH[i] - PH[i + 1])
        TRAAV /= DENOM2
        for i in range(INB + 1):
            FTRA[i, k] -= TRAAV
    return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG


@jit_compile(backend=np, parallel=True)
def _numpy_vectorized_convect(
    t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, tra, params
):
    nlev, ncol = t.shape
    ft = np.zeros(t.shape)
    fq = np.zeros(q.shape)
    fu = np.zeros(u.shape)
    fv = np.zeros(v.shape)
    precip = np.zeros(ncol)
    wd = np.zeros(ncol)
    tprime = np.zeros(ncol)
    qprime = np.zeros(ncol)
    cbmf_new = np.zeros(ncol)
    outcape = np.zeros(ncol)
    iflag = np.zeros(ncol, dtype=np.int32)
    for i in prange(ncol):
        res = _convect_functional_np(
            t[:, i],
            q[:, i],
            qs[:, i],
            u[:, i],
            v[:, i],
            p[:, i],
            ph[:, i],
            ND,
            NL,
            NTRA,
            DELT,
            cbmf[i],
            tra[:, :, i],
            params,
        )
        (
            ft[:, i],
            fq[:, i],
            fu[:, i],
            fv[:, i],
            precip[i],
            wd[i],
            tprime[i],
            qprime[i],
            cbmf_new[i],
            outcape[i],
            iflag[i],
        ) = res
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag


# ---------------------------------------------------------------------------
# JAX path (Phase 2: Python loops, not JIT-able)
# ---------------------------------------------------------------------------


def _bolton_q_sat_jax(T, p, Rd, Rh2O):
    es = 611.2 * jnp.exp(17.67 * (T - 273.15) / (T - 29.65))
    epsilon = Rd / Rh2O
    return epsilon * es / (p - (1 - epsilon) * es)


def _tlift_functional_jax(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params):
    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS
    qNK = float(Q[NK])
    tNK = float(T[NK])
    gzNK = float(GZ[NK])
    AH0 = (
        (params.CPD * (1.0 - qNK) + params.CL * qNK) * tNK
        + qNK * (params.LV0 - CPVMCL * (tNK - 273.15))
        + gzNK
    )
    CPP = params.CPD * (1.0 - qNK) + qNK * params.CPV
    CPINV = 1.0 / CPP
    if KK == 1:
        for i in range(ICB):
            CLW = CLW.at[i].set(0.0)
        for i in range(NK, ICB):
            tpk_i = tNK - (float(GZ[i]) - gzNK) * CPINV
            TPK = TPK.at[i].set(tpk_i)
            TVP = TVP.at[i].set(tpk_i * (1.0 + qNK * EPSI))
    NST = NL if KK == 2 else ICB
    NSB = ICB + 1 if KK == 2 else ICB
    for i in range(NSB, NST + 1):
        TG = float(T[i])
        QG = float(QS[i])
        ALV = params.LV0 - CPVMCL * (float(T[i]) - 273.15)
        ti = float(T[i])
        gzi = float(GZ[i])
        pi = float(P[i])
        for j in range(2):
            S = 1.0 / (params.CPD + ALV * ALV * QG / (params.RV * ti * ti))
            AHG = (
                params.CPD * TG
                + (params.CL - params.CPD) * qNK * ti
                + ALV * QG
                + gzi
            )
            TG = max(TG + S * (AH0 - AHG), 35.0)
            TC = TG - 273.15
            DENOM = 243.5 + TC
            if TC >= 0.0:
                ES = 6.112 * float(jnp.exp(17.67 * TC / DENOM))
            else:
                ES = float(jnp.exp(23.33086 - 6111.72784 / TG + 0.15215 * jnp.log(TG)))
            QG = EPS * ES / (pi - ES * (1.0 - EPS))
        tpk_i = (AH0 - (params.CL - params.CPD) * qNK * ti - gzi - ALV * QG) / params.CPD
        TPK = TPK.at[i].set(tpk_i)
        CLW = CLW.at[i].set(max(0.0, qNK - QG))
        TVP = TVP.at[i].set(tpk_i * (1.0 + (QG / (1.0 - qNK)) * EPSI))
    return TVP, TPK, CLW


def _convect_functional_jax(
    T_in, Q_in, QS_in, U_in, V_in, P_in, PH_in,
    ND, NL, NTRA, DELT, CBMF_in, params,
):
    T = jnp.array(T_in)
    Q = jnp.array(Q_in)
    QS = jnp.array(QS_in)
    U = jnp.array(U_in)
    V = jnp.array(V_in)
    P = jnp.array(P_in)
    PH = jnp.array(PH_in)
    CBMF = float(CBMF_in)

    FT = jnp.zeros(ND)
    FQ = jnp.zeros(ND)
    FU = jnp.zeros(ND)
    FV = jnp.zeros(ND)
    FTRA = jnp.zeros((ND, 1))

    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS
    GINV = 1.0 / params.G
    DELTI = 1.0 / DELT
    PRECIP = 0.0
    WD = 0.0
    TPRIME = 0.0
    QPRIME = 0.0
    IFLAG = 0
    OUTCAPE = 0.0

    # --- Block 2: Thermodynamic profiles ---
    TH = jnp.zeros(NL + 1)
    for i in range(NL + 1):
        qi = float(Q[i])
        RDCP = (params.RD * (1.0 - qi) + qi * params.RV) / (
            params.CPD * (1.0 - qi) + qi * params.CPV
        )
        TH = TH.at[i].set(float(T[i]) * (1000.0 / float(P[i])) ** RDCP)

    GZ = jnp.zeros(ND + 1)
    CPN = jnp.zeros(ND + 1)
    H = jnp.zeros(ND + 1)
    LV = jnp.zeros(ND + 1)
    HM = jnp.zeros(ND + 1)
    TV = jnp.zeros(ND + 1)

    q0 = float(Q[0])
    t0 = float(T[0])
    cpn0 = params.CPD * (1.0 - q0) + q0 * params.CPV
    CPN = CPN.at[0].set(cpn0)
    H = H.at[0].set(t0 * cpn0)
    lv0 = params.LV0 - CPVMCL * (t0 - 273.15)
    LV = LV.at[0].set(lv0)
    HM = HM.at[0].set(lv0 * q0)
    TV = TV.at[0].set(t0 * (1.0 + q0 * EPSI - q0))

    AHMIN = 1.0e12
    IHMIN = NL
    for i in range(1, NL + 1):
        qi = float(Q[i])
        ti = float(T[i])
        qim1 = float(Q[i - 1])
        tim1 = float(T[i - 1])
        TVX = ti * (1.0 + qi * EPSI - qi)
        TVY = tim1 * (1.0 + qim1 * EPSI - qim1)
        gz_i = float(GZ[i - 1]) + 0.5 * params.RD * (TVX + TVY) * (float(P[i - 1]) - float(P[i])) / float(PH[i])
        GZ = GZ.at[i].set(gz_i)
        cpn_i = params.CPD * (1.0 - qi) + params.CPV * qi
        CPN = CPN.at[i].set(cpn_i)
        H = H.at[i].set(ti * cpn_i + gz_i)
        lv_i = params.LV0 - CPVMCL * (ti - 273.15)
        LV = LV.at[i].set(lv_i)
        hm_i = (
            (params.CPD * (1.0 - qi) + params.CL * qi) * (ti - t0)
            + lv_i * qi
            + gz_i
        )
        HM = HM.at[i].set(hm_i)
        TV = TV.at[i].set(TVX)
        if (i + 1) >= params.MINORIG and hm_i < AHMIN and hm_i < float(HM[i - 1]):
            AHMIN = hm_i
            IHMIN = i

    IHMIN = min(IHMIN, NL - 1)

    # --- Block 3: Launch level NK ---
    AHMAX = 0.0
    NK = 0
    for i in range(params.MINORIG - 1, IHMIN + 1):
        hm_i = float(HM[i])
        if hm_i > AHMAX:
            NK = i
            AHMAX = hm_i

    # --- Block 4: Early exits ---
    if float(T[NK]) < 250.0 or float(Q[NK]) <= 0.0 or IHMIN == (NL - 1):
        CBMF = 0.0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    RH = float(Q[NK]) / float(QS[NK])
    CHI = float(T[NK]) / (1669.0 - 122.0 * RH - float(T[NK]))
    PLCL = float(P[NK]) * (RH ** CHI)

    if PLCL < 200.0 or PLCL >= 2000.0:
        IFLAG = 2
        CBMF = 0.0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    ICB = NL - 1
    for i in range(NK + 1, NL + 1):
        if float(P[i]) < PLCL:
            ICB = min(ICB, i)

    if ICB >= (NL - 1):
        IFLAG = 3
        CBMF = 0.0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    # --- Block 5: TLIFT ---
    TVP = jnp.zeros(ND)
    TP = jnp.zeros(ND)
    CLW = jnp.zeros(ND)
    TVP, TP, CLW = _tlift_functional_jax(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params
    )
    for i in range(NK, ICB + 1):
        TVP = TVP.at[i].set(float(TVP[i]) - float(TP[i]) * float(Q[NK]))

    if CBMF == 0.0 and float(TVP[ICB]) <= (float(TV[ICB]) - params.DTMAX):
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    if IFLAG != 4:
        IFLAG = 1

    TVP, TP, CLW = _tlift_functional_jax(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params
    )

    # --- Block 6: EP/SIGP ---
    EP = jnp.zeros(ND)
    SIGP = jnp.zeros(ND)
    for i in range(NK + 1):
        EP = EP.at[i].set(0.0)
        SIGP = SIGP.at[i].set(params.SIGS)

    for i in range(NK + 1, NL + 1):
        TCA = float(TP[i]) - 273.15
        ELACRIT = params.ELCRIT if TCA >= 0.0 else params.ELCRIT * (1.0 - TCA / params.TLCRIT)
        ELACRIT = max(ELACRIT, 0.0)
        EPMAX = 0.999
        ep_i = EPMAX * (1.0 - ELACRIT / max(float(CLW[i]), 1.0e-8))
        ep_i = max(min(ep_i, EPMAX), 0.0)
        EP = EP.at[i].set(ep_i)
        SIGP = SIGP.at[i].set(params.SIGS)

    for i in range(ICB + 1, NL + 1):
        TVP = TVP.at[i].set(float(TVP[i]) - float(TP[i]) * float(Q[NK]))

    # --- Block 7: Setup matrices ---
    HP = jnp.array(H)
    NENT = jnp.zeros(ND + 1, dtype=jnp.int32)
    WATER = jnp.zeros(ND + 1)
    EVAP = jnp.zeros(ND + 1)
    WT = jnp.full(ND + 1, params.OMTSNOW)
    MP = jnp.zeros(ND + 1)
    M = jnp.zeros(ND + 1)
    LVCP = jnp.zeros(ND + 1)
    for i in range(NL + 1):
        LVCP = LVCP.at[i].set(float(LV[i]) / float(CPN[i]))

    QENT = jnp.zeros((ND + 1, ND + 1))
    ELIJ = jnp.zeros((ND + 1, ND + 1))
    MENT = jnp.zeros((ND + 1, ND + 1))
    SIJ = jnp.zeros((ND + 1, ND + 1))
    UENT = jnp.zeros((ND + 1, ND + 1))
    VENT = jnp.zeros((ND + 1, ND + 1))
    for i in range(NL + 1):
        for j in range(NL + 1):
            QENT = QENT.at[i, j].set(float(Q[j]))
            UENT = UENT.at[i, j].set(float(U[j]))
            VENT = VENT.at[i, j].set(float(V[j]))

    QP = jnp.zeros(ND + 1)
    UP = jnp.zeros(ND + 1)
    VP = jnp.zeros(ND + 1)
    QP = QP.at[0].set(float(Q[0]))
    UP = UP.at[0].set(float(U[0]))
    VP = VP.at[0].set(float(V[0]))
    for i in range(1, NL + 1):
        QP = QP.at[i].set(float(Q[i - 1]))
        UP = UP.at[i].set(float(U[i - 1]))
        VP = VP.at[i].set(float(V[i - 1]))

    # --- Block 8: CAPE loop → INB ---
    CAPE = 0.0
    CAPEM = 0.0
    INB = ICB + 1
    INB1 = INB
    BYP = 0.0
    for i in range(ICB + 1, NL):
        BY = (float(TVP[i]) - float(TV[i])) * (float(PH[i]) - float(PH[i + 1])) / float(P[i])
        CAPE += BY
        if BY >= 0.0:
            INB1 = i + 1
        if CAPE > 0.0:
            INB = i + 1
            BYP = (float(TVP[i + 1]) - float(TV[i + 1])) * (float(PH[i + 1]) - float(PH[i + 2])) / float(P[i + 1])
            CAPEM = CAPE
    INB = max(INB, INB1)
    CAPE = CAPEM + BYP
    DEFRAC = max(CAPEM - CAPE, 0.001)
    FRAC = min(max(-CAPE / DEFRAC, 0.0), 1.0)
    OUTCAPE = CAPE

    # --- Block 9: HP update, CBMF ---
    for i in range(ICB, INB + 1):
        HP = HP.at[i].set(
            float(H[NK]) + (float(LV[i]) + (params.CPD - params.CPV) * float(T[i])) * float(EP[i]) * float(CLW[i])
        )

    TVPPLCL = float(TVP[ICB - 1]) - params.RD * float(TVP[ICB - 1]) * (float(P[ICB - 1]) - PLCL) / (
        float(CPN[ICB - 1]) * float(P[ICB - 1])
    )
    TVAPLCL = float(TV[ICB]) + (float(TVP[ICB]) - float(TVP[ICB + 1])) * (PLCL - float(P[ICB])) / (
        float(P[ICB]) - float(P[ICB + 1])
    )
    DTPBL = 0.0
    for i in range(NK, ICB):
        DTPBL += (float(TVP[i]) - float(TV[i])) * (float(PH[i]) - float(PH[i + 1]))
    DTPBL /= float(PH[NK]) - float(PH[ICB])
    DTMA = TVPPLCL - TVAPLCL + params.DTMAX + DTPBL
    CBMFOLD = CBMF
    DAMPS = params.DAMP * DELT / params.DELT0
    CBMF = (1.0 - DAMPS) * CBMF + 0.1 * params.ALPHA * DTMA
    CBMF = max(CBMF, 0.0)

    if CBMF == 0.0 and CBMFOLD == 0.0:
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    # --- Block 10: Updraft mass flux M ---
    M = M.at[ICB].set(0.0)
    DBOSUM = 0.0
    for i in range(ICB + 1, INB + 1):
        k = min(i, INB1)
        DBO = abs(float(TV[k]) - float(TVP[k])) + params.ENTP * 0.02 * (float(PH[k]) - float(PH[k + 1]))
        DBOSUM += DBO
        M = M.at[i].set(CBMF * DBO)
    if DBOSUM > 0:
        for i in range(ICB + 1, INB + 1):
            M = M.at[i].set(float(M[i]) / DBOSUM)

    # --- Block 11: Entrainment/detrainment (SIJ/MENT) ---
    for i in range(ICB + 1, INB + 1):
        QTI = float(Q[NK]) - float(EP[i]) * float(CLW[i])
        for j in range(ICB, INB + 1):
            BF2 = 1.0 + float(LV[j]) * float(LV[j]) * float(QS[j]) / (params.RV * float(T[j]) * float(T[j]) * params.CPD)
            ANUM = float(H[j]) - float(HP[i]) + (params.CPV - params.CPD) * float(T[j]) * (QTI - float(Q[j]))
            DENOM = float(H[i]) - float(HP[i]) + (params.CPD - params.CPV) * (float(Q[i]) - QTI) * float(T[j])
            DEI = DENOM
            if abs(DEI) < 0.01:
                DEI = 0.01
            SIJ = SIJ.at[i, j].set(ANUM / DEI)
            SIJ = SIJ.at[i, i].set(1.0)
            sij_ij = float(SIJ[i, j])  # re-read: 1.0 when j==i
            ALTEM = (sij_ij * float(Q[i]) + (1.0 - sij_ij) * QTI - float(QS[j])) / BF2
            CWAT = float(CLW[j]) * (1.0 - float(EP[j]))
            if (sij_ij < 0.0 or sij_ij > 1.0 or ALTEM > CWAT) and j > i:
                ANUM -= float(LV[j]) * (QTI - float(QS[j]) - CWAT * BF2)
                DENOM += float(LV[j]) * (float(Q[i]) - QTI)
                if abs(DENOM) < 0.01:
                    DENOM = 0.01
                SIJ = SIJ.at[i, j].set(ANUM / DENOM)
                sij_ij = float(SIJ[i, j])
                ALTEM = (
                    sij_ij * float(Q[i])
                    + (1.0 - sij_ij) * QTI
                    - float(QS[j])
                    - (BF2 - 1.0) * CWAT
                )
            if 0.0 < sij_ij < 0.9:
                QENT = QENT.at[i, j].set(sij_ij * float(Q[i]) + (1.0 - sij_ij) * QTI)
                UENT = UENT.at[i, j].set(sij_ij * float(U[i]) + (1.0 - sij_ij) * float(U[NK]))
                VENT = VENT.at[i, j].set(sij_ij * float(V[i]) + (1.0 - sij_ij) * float(V[NK]))
                ELIJ = ELIJ.at[i, j].set(max(0.0, ALTEM))
                MENT = MENT.at[i, j].set(float(M[i]) / (1.0 - sij_ij))
                NENT = NENT.at[i].add(1)
            SIJ = SIJ.at[i, j].set(min(max(sij_ij, 0.0), 1.0))

        if int(NENT[i]) == 0:
            MENT = MENT.at[i, i].set(float(M[i]))
            QENT = QENT.at[i, i].set(float(Q[NK]) - float(EP[i]) * float(CLW[i]))
            UENT = UENT.at[i, i].set(float(U[NK]))
            VENT = VENT.at[i, i].set(float(V[NK]))
            ELIJ = ELIJ.at[i, i].set(float(CLW[i]))
            SIJ = SIJ.at[i, i].set(1.0)

    SIJ = SIJ.at[INB, INB].set(1.0)

    for i in range(ICB + 1, INB + 1):
        if int(NENT[i]) != 0:
            QP1 = float(Q[NK]) - float(EP[i]) * float(CLW[i])
            ANUM = float(H[i]) - float(HP[i]) - float(LV[i]) * (QP1 - float(QS[i]))
            DENOM = float(H[i]) - float(HP[i]) + float(LV[i]) * (float(Q[i]) - QP1)
            if abs(DENOM) < 0.01:
                DENOM = 0.01
            SCRIT = ANUM / DENOM
            ALT = QP1 - float(QS[i]) + SCRIT * (float(Q[i]) - QP1)
            if ALT < 0.0:
                SCRIT = 1.0
            SCRIT = max(SCRIT, 0.0)
            ASIJ = 0.0
            SMIN = 1.0
            for j in range(ICB, INB + 1):
                sij_ij = float(SIJ[i, j])
                if 0.0 < sij_ij < 0.9:
                    if j > i:
                        SMID = min(sij_ij, SCRIT)
                        SJMAX = SMID
                        SJMIN = SMID
                        if SMID < SMIN and float(SIJ[i, j + 1]) < SMID:
                            SMIN = SMID
                            SJMAX = min(min(float(SIJ[i, j + 1]), sij_ij), SCRIT)
                            SJMIN = min(max(float(SIJ[i, j - 1]), sij_ij), SCRIT)
                    else:
                        SJMAX = max(float(SIJ[i, j + 1]), SCRIT)
                        SMID = max(sij_ij, SCRIT)
                        SJMIN = max(float(SIJ[i, j - 1]) if j > 0 else 0.0, SCRIT)
                    ASIJ += (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (
                        float(PH[j]) - float(PH[j + 1])
                    )
                    MENT = MENT.at[i, j].set(
                        float(MENT[i, j]) * (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (
                            float(PH[j]) - float(PH[j + 1])
                        )
                    )
            ASIJ = max(1.0e-21, ASIJ)
            for j in range(ICB, INB + 1):
                MENT = MENT.at[i, j].set(float(MENT[i, j]) / ASIJ)
            s = 0.0
            for j in range(ICB, INB + 1):
                s += float(MENT[i, j])
            if s < 1.0e-18:
                NENT = NENT.at[i].set(0)
                MENT = MENT.at[i, i].set(float(M[i]))
                QENT = QENT.at[i, i].set(float(Q[NK]) - float(EP[i]) * float(CLW[i]))
                UENT = UENT.at[i, i].set(float(U[NK]))
                VENT = VENT.at[i, i].set(float(V[NK]))
                ELIJ = ELIJ.at[i, i].set(float(CLW[i]))
                SIJ = SIJ.at[i, i].set(1.0)

    # --- Block 12: Downdraft ---
    if float(EP[INB]) >= 0.0001:
        JTT = 0
        for i in range(INB, -1, -1):
            WDTRAIN = params.G * float(EP[i]) * float(M[i]) * float(CLW[i])
            if i > 0:
                for j in range(i):
                    WDTRAIN += (
                        params.G
                        * max(0.0, float(ELIJ[j, i]) - (1.0 - float(EP[i])) * float(CLW[i]))
                        * float(MENT[j, i])
                    )
            COEFF = params.COEFFR if float(T[i]) > 273.0 else params.COEFFS
            WT = WT.at[i].set(params.OMTRAIN if float(T[i]) > 273.0 else params.OMTSNOW)
            AFAC = max(
                COEFF
                * float(PH[i])
                * (float(QS[i]) - 0.5 * (float(Q[i]) + float(QP[i + 1])))
                / (1.0e4 + 2.0e3 * float(PH[i]) * float(QS[i])),
                0.0,
            )
            SIGT = min(max(float(SIGP[i]), 0.0), 1.0)
            B6 = 100.0 * (float(PH[i]) - float(PH[i + 1])) * SIGT * AFAC / float(WT[i])
            C6 = (float(WATER[i + 1]) * float(WT[i + 1]) + WDTRAIN / params.SIGD) / float(WT[i])
            REVAP = 0.5 * (-B6 + float(jnp.sqrt(B6 * B6 + 4.0 * C6)))
            EVAP = EVAP.at[i].set(SIGT * AFAC * REVAP)
            WATER = WATER.at[i].set(REVAP * REVAP)

            if i > 0:
                DHDP = max((float(H[i]) - float(H[i - 1])) / (float(P[i - 1]) - float(P[i])), 10.0)
                mp_i = 100.0 * GINV * float(LV[i]) * params.SIGD * float(EVAP[i]) / DHDP
                FAC = 20.0 / (float(PH[i - 1]) - float(PH[i]))
                mp_i = (FAC * float(MP[i + 1]) + mp_i) / (1.0 + FAC)
                MP = MP.at[i].set(mp_i)
                if float(P[i]) > (0.949 * float(P[0])):
                    JTT = max(JTT, i)
                    mp_i = float(MP[JTT]) * (float(P[0]) - float(P[i])) / (float(P[0]) - float(P[JTT]))
                    MP = MP.at[i].set(mp_i)

            if i != INB:
                QSTM = float(QS[0]) if i == 0 else float(QS[i - 1])
                if float(MP[i]) > float(MP[i + 1]):
                    RAT = float(MP[i + 1]) / float(MP[i])
                    QP = QP.at[i].set(
                        float(QP[i + 1]) * RAT
                        + float(Q[i]) * (1.0 - RAT)
                        + 100.0 * GINV * params.SIGD
                        * (float(PH[i]) - float(PH[i + 1]))
                        * (float(EVAP[i]) / float(MP[i]))
                    )
                    UP = UP.at[i].set(float(UP[i + 1]) * RAT + float(U[i]) * (1.0 - RAT))
                    VP = VP.at[i].set(float(VP[i + 1]) * RAT + float(V[i]) * (1.0 - RAT))
                elif float(MP[i + 1]) > 0.0:
                    QP = QP.at[i].set(
                        (float(GZ[i + 1]) - float(GZ[i])
                         + float(QP[i + 1]) * (float(LV[i + 1]) + float(T[i + 1]) * (params.CL - params.CPD))
                         + params.CPD * (float(T[i + 1]) - float(T[i]))
                        ) / (float(LV[i]) + float(T[i]) * (params.CL - params.CPD))
                    )
                    UP = UP.at[i].set(float(UP[i + 1]))
                    VP = VP.at[i].set(float(VP[i + 1]))
                QP = QP.at[i].set(max(min(float(QP[i]), QSTM), 0.0))

        PRECIP += (
            float(WT[0]) * params.SIGD * float(WATER[0]) * 3600.0 * 24000.0 / (params.ROWL * params.G)
        )

    # --- Blocks 13–15: Tendency accumulation ---
    WD = params.BETA * abs(float(MP[ICB])) * 0.01 * params.RD * float(T[ICB]) / (params.SIGD * float(P[ICB]))
    QPRIME = 0.5 * (float(QP[0]) - float(Q[0]))
    TPRIME = params.LV0 * QPRIME / params.CPD

    DPINV = 0.01 / (float(PH[0]) - float(PH[1]))
    AM = 0.0
    if NK == 0:
        for ii in range(1, INB + 1):
            AM += float(M[ii])

    if (2.0 * params.G * DPINV * AM) >= DELTI:
        IFLAG = 4

    FT = FT.at[0].add(
        params.G * DPINV * AM * (float(T[1]) - float(T[0]) + (float(GZ[1]) - float(GZ[0])) / float(CPN[0]))
        - float(LVCP[0]) * params.SIGD * float(EVAP[0])
        + params.SIGD * float(WT[1]) * (params.CL - params.CPD) * float(WATER[1])
        * (float(T[1]) - float(T[0])) * DPINV / float(CPN[0])
    )
    FQ = FQ.at[0].add(
        params.G * float(MP[1]) * (float(QP[1]) - float(Q[0])) * DPINV
        + params.SIGD * float(EVAP[0])
        + params.G * AM * (float(Q[1]) - float(Q[0])) * DPINV
    )
    FU = FU.at[0].add(params.G * DPINV * (float(MP[1]) * (float(UP[1]) - float(U[0])) + AM * (float(U[1]) - float(U[0]))))
    FV = FV.at[0].add(params.G * DPINV * (float(MP[1]) * (float(VP[1]) - float(V[0])) + AM * (float(V[1]) - float(V[0]))))

    for j in range(1, INB + 1):
        FQ = FQ.at[0].add(params.G * DPINV * float(MENT[j, 0]) * (float(QENT[j, 0]) - float(Q[0])))
        FU = FU.at[0].add(params.G * DPINV * float(MENT[j, 0]) * (float(UENT[j, 0]) - float(U[0])))
        FV = FV.at[0].add(params.G * DPINV * float(MENT[j, 0]) * (float(VENT[j, 0]) - float(V[0])))

    for i in range(1, INB + 1):
        DPINV = 0.01 / (float(PH[i]) - float(PH[i + 1]))
        CPINV = 1.0 / float(CPN[i])
        AMP1 = 0.0
        if i >= NK:
            for ii in range(i + 1, INB + 2):
                AMP1 += float(M[ii])
        for k in range(i + 1):
            for jj in range(i + 1, INB + 2):
                AMP1 += float(MENT[k, jj])
        if (2.0 * params.G * DPINV * AMP1) >= DELTI:
            IFLAG = 4
        AD = 0.0
        for k in range(i):
            for jj in range(i, INB + 1):
                AD += float(MENT[jj, k])

        FT = FT.at[i].add(
            params.G * DPINV * (
                AMP1 * (float(T[i + 1]) - float(T[i]) + (float(GZ[i + 1]) - float(GZ[i])) * CPINV)
                - AD * (float(T[i]) - float(T[i - 1]) + (float(GZ[i]) - float(GZ[i - 1])) * CPINV)
            )
            - params.SIGD * float(LVCP[i]) * float(EVAP[i])
        )
        FT = FT.at[i].add(
            params.G * DPINV * float(MENT[i, i])
            * (float(HP[i]) - float(H[i]) + float(T[i]) * (params.CPV - params.CPD) * (float(Q[i]) - float(QENT[i, i])))
            * CPINV
        )
        FT = FT.at[i].add(
            params.SIGD * float(WT[i + 1]) * (params.CL - params.CPD)
            * float(WATER[i + 1]) * (float(T[i + 1]) - float(T[i]))
            * DPINV * CPINV
        )
        FQ = FQ.at[i].add(params.G * DPINV * (AMP1 * (float(Q[i + 1]) - float(Q[i])) - AD * (float(Q[i]) - float(Q[i - 1]))))
        FU = FU.at[i].add(params.G * DPINV * (AMP1 * (float(U[i + 1]) - float(U[i])) - AD * (float(U[i]) - float(U[i - 1]))))
        FV = FV.at[i].add(params.G * DPINV * (AMP1 * (float(V[i + 1]) - float(V[i])) - AD * (float(V[i]) - float(V[i - 1]))))

        for k in range(i):
            AWAT = max(float(ELIJ[k, i]) - (1.0 - float(EP[i])) * float(CLW[i]), 0.0)
            FQ = FQ.at[i].add(params.G * DPINV * float(MENT[k, i]) * (float(QENT[k, i]) - AWAT - float(Q[i])))
            FU = FU.at[i].add(params.G * DPINV * float(MENT[k, i]) * (float(UENT[k, i]) - float(U[i])))
            FV = FV.at[i].add(params.G * DPINV * float(MENT[k, i]) * (float(VENT[k, i]) - float(V[i])))

        for k in range(i, INB + 1):
            FQ = FQ.at[i].add(params.G * DPINV * float(MENT[k, i]) * (float(QENT[k, i]) - float(Q[i])))
            FU = FU.at[i].add(params.G * DPINV * float(MENT[k, i]) * (float(UENT[k, i]) - float(U[i])))
            FV = FV.at[i].add(params.G * DPINV * float(MENT[k, i]) * (float(VENT[k, i]) - float(V[i])))

        FQ = FQ.at[i].add(
            params.SIGD * float(EVAP[i])
            + params.G * (float(MP[i + 1]) * (float(QP[i + 1]) - float(Q[i])) - float(MP[i]) * (float(QP[i]) - float(Q[i - 1])))
            * DPINV
        )
        FU = FU.at[i].add(
            params.G * (float(MP[i + 1]) * (float(UP[i + 1]) - float(U[i])) - float(MP[i]) * (float(UP[i]) - float(U[i - 1])))
            * DPINV
        )
        FV = FV.at[i].add(
            params.G * (float(MP[i + 1]) * (float(VP[i + 1]) - float(V[i])) - float(MP[i]) * (float(VP[i]) - float(V[i - 1])))
            * DPINV
        )

    # --- Block 14: FRAC smoothing at INB ---
    FQOLD = float(FQ[INB])
    FQ = FQ.at[INB].set(FQOLD * (1.0 - FRAC))
    FQ = FQ.at[INB - 1].add(
        FRAC * FQOLD * ((float(PH[INB]) - float(PH[INB + 1])) / (float(PH[INB - 1]) - float(PH[INB])))
        * float(LV[INB]) / float(LV[INB - 1])
    )
    FTOLD = float(FT[INB])
    FT = FT.at[INB].set(FTOLD * (1.0 - FRAC))
    FT = FT.at[INB - 1].add(
        FRAC * FTOLD * ((float(PH[INB]) - float(PH[INB + 1])) / (float(PH[INB - 1]) - float(PH[INB])))
        * float(CPN[INB]) / float(CPN[INB - 1])
    )
    FUOLD = float(FU[INB])
    FU = FU.at[INB].set(FUOLD * (1.0 - FRAC))
    FU = FU.at[INB - 1].add(
        FRAC * FUOLD * ((float(PH[INB]) - float(PH[INB + 1])) / (float(PH[INB - 1]) - float(PH[INB])))
    )
    FVOLD = float(FV[INB])
    FV = FV.at[INB].set(FVOLD * (1.0 - FRAC))
    FV = FV.at[INB - 1].add(
        FRAC * FVOLD * ((float(PH[INB]) - float(PH[INB + 1])) / (float(PH[INB - 1]) - float(PH[INB])))
    )

    # --- Block 15: Conservation correction ---
    ENTS = 0.0
    UAV = 0.0
    VAV = 0.0
    DENOM2 = float(PH[0]) - float(PH[INB + 1])
    for i in range(INB + 1):
        dp = float(PH[i]) - float(PH[i + 1])
        ENTS += (float(CPN[i]) * float(FT[i]) + float(LV[i]) * float(FQ[i])) * dp
        UAV += float(FU[i]) * dp
        VAV += float(FV[i]) * dp
    ENTS /= DENOM2
    UAV /= DENOM2
    VAV /= DENOM2
    for i in range(INB + 1):
        FT = FT.at[i].set(float(FT[i]) - ENTS / float(CPN[i]))
        FU = FU.at[i].set((1.0 - params.CU) * (float(FU[i]) - UAV))
        FV = FV.at[i].set((1.0 - params.CU) * (float(FV[i]) - VAV))

    return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG


def _jax_vectorized_convect(
    t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, params,
):
    nlev, ncol = t.shape
    ft_list, fq_list, fu_list, fv_list = [], [], [], []
    precip_list, wd_list, tprime_list, qprime_list = [], [], [], []
    cbmf_list, outcape_list, iflag_list = [], [], []
    for i in range(ncol):
        res = _convect_functional_jax(
            t[:, i], q[:, i], qs[:, i], u[:, i], v[:, i],
            p[:, i], ph[:, i],
            ND, NL, NTRA, DELT, cbmf[i], params,
        )
        ft_list.append(res[0])
        fq_list.append(res[1])
        fu_list.append(res[2])
        fv_list.append(res[3])
        precip_list.append(res[4])
        wd_list.append(res[5])
        tprime_list.append(res[6])
        qprime_list.append(res[7])
        cbmf_list.append(res[8])
        outcape_list.append(res[9])
        iflag_list.append(res[10])
    ft = jnp.stack(ft_list, axis=1)
    fq = jnp.stack(fq_list, axis=1)
    fu = jnp.stack(fu_list, axis=1)
    fv = jnp.stack(fv_list, axis=1)
    precip = jnp.array(precip_list)
    wd = jnp.array(wd_list)
    tprime = jnp.array(tprime_list)
    qprime = jnp.array(qprime_list)
    cbmf_new = jnp.array(cbmf_list)
    outcape = jnp.array(outcape_list)
    iflag = jnp.array(iflag_list, dtype=jnp.int32)
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
