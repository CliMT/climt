# -*- coding: utf-8 -*-
import functools
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
    # KK is always a Python literal (1 or 2) — static control flow on KK is fine.
    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS
    qNK = Q[NK]
    tNK = T[NK]
    gzNK = GZ[NK]
    AH0 = (
        (params.CPD * (1.0 - qNK) + params.CL * qNK) * tNK
        + qNK * (params.LV0 - CPVMCL * (tNK - 273.15))
        + gzNK
    )
    CPP = params.CPD * (1.0 - qNK) + qNK * params.CPV
    CPINV = 1.0 / CPP
    if KK == 1:
        # ICB and NK are traced → use full range with masks
        for i in range(NL + 1):
            CLW = CLW.at[i].set(jnp.where(i < ICB, 0.0, CLW[i]))
        for i in range(NL + 1):
            tpk_i = tNK - (GZ[i] - gzNK) * CPINV
            in_range = (i >= NK) & (i < ICB)
            TPK = TPK.at[i].set(jnp.where(in_range, tpk_i, TPK[i]))
            TVP = TVP.at[i].set(jnp.where(in_range, tpk_i * (1.0 + qNK * EPSI), TVP[i]))
    # Outer saturation-adjustment loop.
    # KK==1: original range(ICB, ICB+1) → mask i == ICB
    # KK==2: original range(ICB+1, NL+1) → mask i > ICB
    # ICB is traced, so use full static range with jnp.where guards.
    for i in range(NL + 1):
        if KK == 1:
            loop_mask = (i == ICB)
        else:
            loop_mask = (i > ICB)
        TG = T[i]
        QG = QS[i]
        ALV = params.LV0 - CPVMCL * (T[i] - 273.15)
        ti = T[i]
        gzi = GZ[i]
        pi = P[i]
        for j in range(2):
            S = 1.0 / (params.CPD + ALV * ALV * QG / (params.RV * ti * ti))
            AHG = (
                params.CPD * TG
                + (params.CL - params.CPD) * qNK * ti
                + ALV * QG
                + gzi
            )
            TG = jnp.maximum(TG + S * (AH0 - AHG), 35.0)
            TC = TG - 273.15
            DENOM = 243.5 + TC
            ES = jnp.where(
                TC >= 0.0,
                6.112 * jnp.exp(17.67 * TC / DENOM),
                jnp.exp(23.33086 - 6111.72784 / TG + 0.15215 * jnp.log(TG)),
            )
            QG = EPS * ES / (pi - ES * (1.0 - EPS))
        tpk_i = (AH0 - (params.CL - params.CPD) * qNK * ti - gzi - ALV * QG) / params.CPD
        TPK = TPK.at[i].set(jnp.where(loop_mask, tpk_i, TPK[i]))
        CLW = CLW.at[i].set(jnp.where(loop_mask, jnp.maximum(0.0, qNK - QG), CLW[i]))
        TVP = TVP.at[i].set(jnp.where(loop_mask, tpk_i * (1.0 + (QG / (1.0 - qNK)) * EPSI), TVP[i]))
    return TVP, TPK, CLW


@functools.partial(jax.jit, static_argnums=(7, 8, 9))
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
    CBMF = CBMF_in

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
        qi = Q[i]
        RDCP = (params.RD * (1.0 - qi) + qi * params.RV) / (
            params.CPD * (1.0 - qi) + qi * params.CPV
        )
        TH = TH.at[i].set(T[i] * (1000.0 / P[i]) ** RDCP)

    GZ = jnp.zeros(ND + 1)
    CPN = jnp.zeros(ND + 1)
    H = jnp.zeros(ND + 1)
    LV = jnp.zeros(ND + 1)
    HM = jnp.zeros(ND + 1)
    TV = jnp.zeros(ND + 1)

    q0 = Q[0]
    t0 = T[0]
    cpn0 = params.CPD * (1.0 - q0) + q0 * params.CPV
    CPN = CPN.at[0].set(cpn0)
    H = H.at[0].set(t0 * cpn0)
    lv0 = params.LV0 - CPVMCL * (t0 - 273.15)
    LV = LV.at[0].set(lv0)
    HM = HM.at[0].set(lv0 * q0)
    TV = TV.at[0].set(t0 * (1.0 + q0 * EPSI - q0))

    for i in range(1, NL + 1):
        qi = Q[i]
        ti = T[i]
        qim1 = Q[i - 1]
        tim1 = T[i - 1]
        TVX = ti * (1.0 + qi * EPSI - qi)
        TVY = tim1 * (1.0 + qim1 * EPSI - qim1)
        gz_i = GZ[i - 1] + 0.5 * params.RD * (TVX + TVY) * (P[i - 1] - P[i]) / PH[i]
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

    # Vectorized IHMIN: find index where HM is local minimum and i+1 >= MINORIG
    indices = jnp.arange(ND + 1)
    hm_prev = jnp.roll(HM, 1).at[0].set(jnp.inf)  # HM[-1] sentinel = inf so HM[0] < HM[-1] never true at i=0
    ihmin_valid = (indices >= 1) & ((indices + 1) >= params.MINORIG) & (indices <= NL) & (HM < hm_prev)
    ihmin_vals = jnp.where(ihmin_valid, HM, jnp.inf)
    IHMIN = jnp.argmin(ihmin_vals)
    # If no valid entry found (all inf), default to NL
    IHMIN = jnp.where(jnp.all(~ihmin_valid), NL, IHMIN)
    IHMIN = jnp.minimum(IHMIN, NL - 1)

    # --- Block 3: Launch level NK ---
    # Vectorized NK: argmax of HM in [MINORIG-1, IHMIN]
    nk_valid = (indices >= (params.MINORIG - 1)) & (indices <= IHMIN)
    nk_vals = jnp.where(nk_valid, HM, -jnp.inf)
    NK = jnp.argmax(nk_vals)

    # --- Block 4: Early exits (accumulated as `active` flag for JIT) ---
    active = jnp.ones((), dtype=jnp.bool_)

    # Exit condition 1: T[NK] < 250 or Q[NK] <= 0 or IHMIN == NL-1
    exit1 = (T[NK] < 250.0) | (Q[NK] <= 0.0) | (IHMIN == (NL - 1))
    active = active & ~exit1

    RH = Q[NK] / jnp.maximum(QS[NK], 1.0e-30)
    CHI = T[NK] / jnp.maximum(1669.0 - 122.0 * RH - T[NK], 1.0e-30)
    PLCL = P[NK] * (RH ** CHI)

    # Exit condition 2: PLCL out of range
    exit2 = (PLCL < 200.0) | (PLCL >= 2000.0)
    IFLAG = jnp.where(exit2 & active, 2, IFLAG)
    active = active & ~exit2

    # Vectorized ICB: first index > NK where P < PLCL, default NL-1
    # P has length ND; use indices[:ND] for comparison
    p_indices = jnp.arange(ND)
    icb_candidates = (p_indices > NK) & (p_indices <= NL) & (P < PLCL)
    ICB = jnp.where(jnp.any(icb_candidates), jnp.argmax(icb_candidates), NL - 1)

    # Exit condition 3: ICB >= NL-1
    exit3 = ICB >= (NL - 1)
    IFLAG = jnp.where(exit3 & active, 3, IFLAG)
    active = active & ~exit3

    # --- Block 5: TLIFT ---
    TVP = jnp.zeros(ND)
    TP = jnp.zeros(ND)
    CLW = jnp.zeros(ND)
    TVP, TP, CLW = _tlift_functional_jax(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params
    )
    # Variable-bound loop [NK, ICB+1) → full range with mask
    for i in range(NL + 1):
        in_range = (i >= NK) & (i <= ICB)
        TVP = TVP.at[i].set(jnp.where(in_range, TVP[i] - TP[i] * Q[NK], TVP[i]))

    # Early return: CBMF == 0 and TVP[ICB] too cold → fold into active flag
    exit4 = (CBMF == 0.0) & (TVP[ICB] <= (TV[ICB] - params.DTMAX))
    active = active & ~exit4

    IFLAG = jnp.where(IFLAG != 4, 1, IFLAG)

    TVP, TP, CLW = _tlift_functional_jax(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params
    )

    # --- Block 6: EP/SIGP ---
    # Full-range computation with jnp.where masks instead of traced-bound loops
    EP = jnp.zeros(ND)
    SIGP = jnp.full(ND, params.SIGS)

    for i in range(NL + 1):
        in_nk_range = (i > NK)
        TCA = TP[i] - 273.15
        ELACRIT = jnp.where(TCA >= 0.0, params.ELCRIT,
                            params.ELCRIT * (1.0 - TCA / params.TLCRIT))
        ELACRIT = jnp.maximum(ELACRIT, 0.0)
        EPMAX = 0.999
        ep_i = EPMAX * (1.0 - ELACRIT / jnp.maximum(CLW[i], 1.0e-8))
        ep_i = jnp.maximum(jnp.minimum(ep_i, EPMAX), 0.0)
        EP = EP.at[i].set(jnp.where(in_nk_range, ep_i, 0.0))

    # TVP adjustment for [ICB+1, NL+1) → full range with mask
    for i in range(NL + 1):
        in_range = (i > ICB)
        TVP = TVP.at[i].set(jnp.where(in_range, TVP[i] - TP[i] * Q[NK], TVP[i]))

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
        LVCP = LVCP.at[i].set(LV[i] / CPN[i])

    ELIJ = jnp.zeros((ND + 1, ND + 1))
    MENT = jnp.zeros((ND + 1, ND + 1))
    SIJ = jnp.zeros((ND + 1, ND + 1))
    # Broadcast Q, U, V into (ND+1, ND+1) matrices: QENT[i,j] = Q[j] for all i
    Q_pad_init = jnp.concatenate([Q, jnp.zeros(ND + 1 - ND)])  # pad to ND+1
    U_pad_init = jnp.concatenate([U, jnp.zeros(ND + 1 - ND)])
    V_pad_init = jnp.concatenate([V, jnp.zeros(ND + 1 - ND)])
    QENT = jnp.broadcast_to(Q_pad_init[None, :], (ND + 1, ND + 1)).copy()
    UENT = jnp.broadcast_to(U_pad_init[None, :], (ND + 1, ND + 1)).copy()
    VENT = jnp.broadcast_to(V_pad_init[None, :], (ND + 1, ND + 1)).copy()

    # QP[0] = Q[0], QP[i] = Q[i-1] for i>=1 to NL; zero beyond NL
    QP = jnp.concatenate([Q[:1], Q[:NL], jnp.zeros(ND - NL)])
    UP = jnp.concatenate([U[:1], U[:NL], jnp.zeros(ND - NL)])
    VP = jnp.concatenate([V[:1], V[:NL], jnp.zeros(ND - NL)])

    # --- Block 8: CAPE loop → INB (JIT-compatible vectorized) ---
    # Compute buoyancy BY for all levels
    BY_all = (TVP[:ND] - TV[:ND]) * (PH[:ND] - PH[1:ND + 1]) / P
    cape_lvl = jnp.arange(ND)
    cape_mask = (cape_lvl > ICB) & (cape_lvl < NL)
    BY_masked = jnp.where(cape_mask, BY_all, 0.0)
    CAPE_cum = jnp.cumsum(BY_masked)

    # INB1: (last level in range where BY >= 0) + 1, default ICB+1
    inb1_cand = cape_mask & (BY_all >= 0.0)
    INB1 = jnp.where(jnp.any(inb1_cand),
                      jnp.max(jnp.where(inb1_cand, cape_lvl, -1)) + 1,
                      ICB + 1)

    # INB from CAPE: (last level where cumulative CAPE > 0) + 1, default ICB+1
    inb_cand = cape_mask & (CAPE_cum > 0.0)
    last_pos_i = jnp.max(jnp.where(inb_cand, cape_lvl, -1))
    has_pos_cape = jnp.any(inb_cand)
    INB_cape = jnp.where(has_pos_cape, last_pos_i + 1, ICB + 1)

    # CAPEM = cumulative CAPE at last positive level
    CAPEM = jnp.where(has_pos_cape, CAPE_cum[last_pos_i], 0.0)
    # BYP = buoyancy one level beyond last positive CAPE level
    byp_idx = jnp.minimum(last_pos_i + 1, ND - 1)
    BYP = jnp.where(has_pos_cape, BY_all[byp_idx], 0.0)

    INB = jnp.maximum(INB_cape, INB1)
    CAPE = CAPEM + BYP
    DEFRAC = jnp.maximum(CAPEM - CAPE, 0.001)
    FRAC = jnp.minimum(jnp.maximum(-CAPE / DEFRAC, 0.0), 1.0)
    OUTCAPE = CAPE

    # --- Block 9: HP update, CBMF (JIT-compatible) ---
    # HP update: full range with mask for [ICB, INB+1)
    hp_lvl = jnp.arange(ND + 1)
    hp_mask = (hp_lvl >= ICB) & (hp_lvl <= INB)
    # Pad T, EP, CLW from ND to ND+1 for uniform indexing
    _T_hp = jnp.concatenate([T, jnp.zeros(1)])
    _EP_hp = jnp.concatenate([EP, jnp.zeros(1)])
    _CLW_hp = jnp.concatenate([CLW, jnp.zeros(1)])
    HP_new = H[NK] + (LV + (params.CPD - params.CPV) * _T_hp) * _EP_hp * _CLW_hp
    HP = jnp.where(hp_mask, HP_new, HP)

    TVPPLCL = TVP[ICB - 1] - params.RD * TVP[ICB - 1] * (P[ICB - 1] - PLCL) / (
        CPN[ICB - 1] * P[ICB - 1]
    )
    TVAPLCL = TV[ICB] + (TVP[ICB] - TVP[ICB + 1]) * (PLCL - P[ICB]) / (
        P[ICB] - P[ICB + 1]
    )
    # DTPBL: sum over [NK, ICB) with mask
    dtpbl_mask = (hp_lvl[:ND] >= NK) & (hp_lvl[:ND] < ICB)
    DTPBL = jnp.sum(jnp.where(dtpbl_mask, (TVP[:ND] - TV[:ND]) * (PH[:ND] - PH[1:ND + 1]), 0.0))
    DTPBL = DTPBL / jnp.maximum(PH[NK] - PH[ICB], 1e-30)
    DTMA = TVPPLCL - TVAPLCL + params.DTMAX + DTPBL
    CBMFOLD = CBMF
    DAMPS = params.DAMP * DELT / params.DELT0
    CBMF = (1.0 - DAMPS) * CBMF + 0.1 * params.ALPHA * DTMA
    CBMF = jnp.maximum(CBMF, 0.0)

    # Early return → active flag
    exit5 = (CBMF == 0.0) & (CBMFOLD == 0.0)
    active = active & ~exit5

    # --- Block 10: Updraft mass flux M (differentiable) ---
    M = M.at[ICB].set(0.0)
    lvl = jnp.arange(ND + 1)
    k_idx = jnp.minimum(lvl, INB1)
    DBO_all = jnp.abs(TV[k_idx] - TVP[k_idx]) + params.ENTP * 0.02 * (PH[k_idx] - PH[k_idx + 1])
    mass_mask = (lvl >= ICB + 1) & (lvl <= INB)
    DBO_masked = jnp.where(mass_mask, DBO_all, 0.0)
    DBOSUM = jnp.sum(DBO_masked)
    M = jnp.where(mass_mask, CBMF * DBO_masked / jnp.maximum(DBOSUM, 1e-30), 0.0)

    # --- Block 11: Entrainment/detrainment (SIJ/MENT) (differentiable) ---
    # Pad ND-sized arrays to ND+1 for 2D broadcasting
    _z1 = jnp.zeros(1)
    T_pad = jnp.concatenate([T, T[-1:]])  # safe non-zero sentinel for gradient
    Q_pad = jnp.concatenate([Q, _z1])
    QS_pad = jnp.concatenate([QS, _z1])
    U_pad = jnp.concatenate([U, _z1])
    V_pad = jnp.concatenate([V, _z1])
    EP_pad = jnp.concatenate([EP, _z1])
    CLW_pad = jnp.concatenate([CLW, _z1])

    # 2D index arrays: ii = updraft level, jj = environment level
    ii = jnp.arange(ND + 1)[:, None]  # (ND+1, 1)
    jj = jnp.arange(ND + 1)[None, :]  # (1, ND+1)

    # Active masks
    i_active = (ii >= ICB + 1) & (ii <= INB)
    j_active = (jj >= ICB) & (jj <= INB)
    ij_active = i_active & j_active

    # QTI[i] = Q[NK] - EP[i] * CLW[i]
    QTI_vec = Q_pad[NK] - EP_pad * CLW_pad  # (ND+1,)

    # BF2[j]
    BF2_vec = 1.0 + LV * LV * QS_pad / (params.RV * T_pad * T_pad * params.CPD)  # (ND+1,)

    # First-pass SIJ computation for all (i,j) pairs
    # ANUM[i,j] = H[j] - HP[i] + (CPV - CPD) * T[j] * (QTI[i] - Q[j])
    ANUM_2d = H[None, :] - HP[:, None] + (params.CPV - params.CPD) * T_pad[None, :] * (QTI_vec[:, None] - Q_pad[None, :])
    # DENOM[i,j] = H[i] - HP[i] + (CPD - CPV) * (Q[i] - QTI[i]) * T[j]
    DENOM_2d = H[:, None] - HP[:, None] + (params.CPD - params.CPV) * (Q_pad[:, None] - QTI_vec[:, None]) * T_pad[None, :]
    DEI_2d = jnp.where(jnp.abs(DENOM_2d) < 0.01, 0.01, DENOM_2d)
    SIJ_first = ANUM_2d / DEI_2d

    # For diagonal (j==i), SIJ = 1.0 — but we compute off-diagonal first, then set diagonal
    diag_mask = (ii == jj)
    SIJ_first = jnp.where(diag_mask, 1.0, SIJ_first)

    # ALTEM and CWAT for first pass
    ALTEM_first = (SIJ_first * Q_pad[:, None] + (1.0 - SIJ_first) * QTI_vec[:, None] - QS_pad[None, :]) / BF2_vec[None, :]
    CWAT_2d = CLW_pad[None, :] * (1.0 - EP_pad[None, :])

    # Recalculation condition: (SIJ < 0 or SIJ > 1 or ALTEM > CWAT) and j > i
    recalc_cond = ((SIJ_first < 0.0) | (SIJ_first > 1.0) | (ALTEM_first > CWAT_2d)) & (jj > ii)

    # Second-pass ANUM/DENOM for recalculation
    ANUM_recalc = ANUM_2d - LV[None, :] * (QTI_vec[:, None] - QS_pad[None, :] - CWAT_2d * BF2_vec[None, :])
    DENOM_recalc = DENOM_2d + LV[None, :] * (Q_pad[:, None] - QTI_vec[:, None])
    DENOM_recalc = jnp.where(jnp.abs(DENOM_recalc) < 0.01, 0.01, DENOM_recalc)
    SIJ_second = ANUM_recalc / DENOM_recalc

    # Pick recalculated or first-pass SIJ
    SIJ_raw = jnp.where(recalc_cond & ~diag_mask, SIJ_second, SIJ_first)

    # Recompute ALTEM with final SIJ
    ALTEM_final_norec = (SIJ_raw * Q_pad[:, None] + (1.0 - SIJ_raw) * QTI_vec[:, None] - QS_pad[None, :]) / BF2_vec[None, :]
    ALTEM_final_rec = (SIJ_raw * Q_pad[:, None] + (1.0 - SIJ_raw) * QTI_vec[:, None] - QS_pad[None, :] - (BF2_vec[None, :] - 1.0) * CWAT_2d)
    ALTEM_final = jnp.where(recalc_cond & ~diag_mask, ALTEM_final_rec, ALTEM_final_norec)

    # Entrainment condition: 0 < SIJ < 0.9 (and not diagonal)
    ent_cond = (SIJ_raw > 0.0) & (SIJ_raw < 0.9) & ~diag_mask & ij_active

    # Compute entrained quantities
    QENT_ij = SIJ_raw * Q_pad[:, None] + (1.0 - SIJ_raw) * QTI_vec[:, None]
    UENT_ij = SIJ_raw * U_pad[:, None] + (1.0 - SIJ_raw) * U_pad[NK]
    VENT_ij = SIJ_raw * V_pad[:, None] + (1.0 - SIJ_raw) * V_pad[NK]
    ELIJ_ij = jnp.maximum(0.0, ALTEM_final)
    MENT_ij = M[:, None] / jnp.maximum(1.0 - SIJ_raw, 1e-30)

    # Apply entrainment condition
    QENT = jnp.where(ent_cond, QENT_ij, QENT)
    UENT = jnp.where(ent_cond, UENT_ij, UENT)
    VENT = jnp.where(ent_cond, VENT_ij, VENT)
    ELIJ = jnp.where(ent_cond, ELIJ_ij, ELIJ)
    MENT = jnp.where(ent_cond, MENT_ij, MENT)

    # Count entraining levels per updraft level i
    NENT_count = jnp.sum(ent_cond.astype(jnp.int32), axis=1)  # (ND+1,)

    # Clamp SIJ to [0, 1] for non-diagonal active entries
    SIJ_clamped = jnp.clip(SIJ_raw, 0.0, 1.0)
    SIJ = jnp.where(ij_active & ~diag_mask, SIJ_clamped, SIJ)
    # Set diagonal to 1.0 for active i levels
    SIJ = jnp.where(diag_mask & i_active, 1.0, SIJ)

    # Fallback for levels with no entraining levels (NENT==0) — vectorized
    no_ent_mask = (NENT_count == 0) & (lvl >= ICB + 1) & (lvl <= INB)  # (ND+1,)
    diag_no_ent = diag_mask & no_ent_mask[:, None]  # (ND+1, ND+1) — diagonal of no-ent levels
    MENT = jnp.where(diag_no_ent, M[:, None], MENT)
    QENT = jnp.where(diag_no_ent, (Q_pad[NK] - EP_pad * CLW_pad)[:, None], QENT)
    UENT = jnp.where(diag_no_ent, U_pad[NK], UENT)
    VENT = jnp.where(diag_no_ent, V_pad[NK], VENT)
    ELIJ = jnp.where(diag_no_ent, CLW_pad[:, None], ELIJ)
    SIJ = jnp.where(diag_no_ent, 1.0, SIJ)

    SIJ = SIJ.at[INB, INB].set(1.0)

    # --- Block 11 second half: SCRIT normalization (JIT-compatible) ---
    # Static-range loops with masks (ICB, INB are traced)
    for i in range(NL + 1):
        i_active = (i >= ICB + 1) & (i <= INB) & (NENT_count[i] > 0)

        QP1 = Q_pad[NK] - EP_pad[i] * CLW_pad[i]
        ANUM_s = H[i] - HP[i] - LV[i] * (QP1 - QS_pad[i])
        DENOM_s = H[i] - HP[i] + LV[i] * (Q_pad[i] - QP1)
        DENOM_s = jnp.where(jnp.abs(DENOM_s) < 0.01, 0.01, DENOM_s)
        SCRIT = ANUM_s / DENOM_s
        ALT = QP1 - QS_pad[i] + SCRIT * (Q_pad[i] - QP1)
        SCRIT = jnp.where(ALT < 0.0, 1.0, SCRIT)
        SCRIT = jnp.maximum(SCRIT, 0.0)

        ASIJ_acc = 0.0
        SMIN = 1.0
        weights = jnp.zeros(ND + 1)
        for j in range(NL + 1):
            j_valid = (j >= ICB) & (j <= INB)
            sij_ij = SIJ[i, j]
            is_ent_j = (sij_ij > 0.0) & (sij_ij < 0.9) & j_valid
            if j > i:
                SMID = jnp.minimum(sij_ij, SCRIT)
                SJMAX = SMID
                SJMIN = SMID
                sij_next = SIJ[i, min(j + 1, ND)]
                sij_prev = SIJ[i, max(j - 1, 0)]
                update_minmax = (SMID < SMIN) & (sij_next < SMID)
                SJMAX = jnp.where(update_minmax, jnp.minimum(jnp.minimum(sij_next, sij_ij), SCRIT), SJMAX)
                SJMIN = jnp.where(update_minmax, jnp.minimum(jnp.maximum(sij_prev, sij_ij), SCRIT), SJMIN)
                SMIN = jnp.where(update_minmax & is_ent_j, SMID, SMIN)
            else:
                sij_next = SIJ[i, min(j + 1, ND)]
                sij_prev = SIJ[i, max(j - 1, 0)] if j > 0 else 0.0
                SJMAX = jnp.maximum(sij_next, SCRIT)
                SMID = jnp.maximum(sij_ij, SCRIT)
                SJMIN = jnp.maximum(sij_prev, SCRIT)

            w = (jnp.abs(SJMAX - SMID) + jnp.abs(SJMIN - SMID)) * (PH[j] - PH[j + 1])
            w = jnp.where(is_ent_j, w, 0.0)
            weights = weights.at[j].set(w)
            ASIJ_acc = ASIJ_acc + w

        ASIJ_acc = jnp.maximum(ASIJ_acc, 1.0e-21)
        MENT_row = MENT[i, :] * weights / ASIJ_acc
        ment_sum = jnp.sum(MENT_row)

        # If sum of MENT is tiny, fall back to diagonal
        tiny = (ment_sum < 1.0e-18)
        MENT_fallback = jnp.zeros(ND + 1).at[i].set(M[i])
        QENT_fallback = QENT[i, :].at[i].set(Q_pad[NK] - EP_pad[i] * CLW_pad[i])
        UENT_fallback = UENT[i, :].at[i].set(U_pad[NK])
        VENT_fallback = VENT[i, :].at[i].set(V_pad[NK])
        ELIJ_fallback = ELIJ[i, :].at[i].set(CLW_pad[i])
        SIJ_fallback = SIJ[i, :].at[i].set(1.0)

        new_ment_row = jnp.where(tiny, MENT_fallback, MENT_row)
        new_qent_row = jnp.where(tiny, QENT_fallback, QENT[i, :])
        new_uent_row = jnp.where(tiny, UENT_fallback, UENT[i, :])
        new_vent_row = jnp.where(tiny, VENT_fallback, VENT[i, :])
        new_elij_row = jnp.where(tiny, ELIJ_fallback, ELIJ[i, :])
        new_sij_row = jnp.where(tiny, SIJ_fallback, SIJ[i, :])

        # Only update if i_active (has entraining levels and in valid range)
        MENT = MENT.at[i, :].set(jnp.where(i_active, new_ment_row, MENT[i, :]))
        QENT = QENT.at[i, :].set(jnp.where(i_active, new_qent_row, QENT[i, :]))
        UENT = UENT.at[i, :].set(jnp.where(i_active, new_uent_row, UENT[i, :]))
        VENT = VENT.at[i, :].set(jnp.where(i_active, new_vent_row, VENT[i, :]))
        ELIJ = ELIJ.at[i, :].set(jnp.where(i_active, new_elij_row, ELIJ[i, :]))
        SIJ = SIJ.at[i, :].set(jnp.where(i_active, new_sij_row, SIJ[i, :]))

    # --- Block 12: Downdraft (differentiable) ---
    # Pre-compute WDTRAIN for all levels as a vectorized operation.
    # WDTRAIN[i] = G*EP[i]*M[i]*CLW[i] + sum_{j<i} G*max(0, ELIJ[j,i]-(1-EP[i])*CLW[i])*MENT[j,i]
    WDTRAIN_base = params.G * EP_pad * M * CLW_pad  # (ND+1,)
    elij_term = ELIJ - (1.0 - EP_pad[None, :]) * CLW_pad[None, :]  # (ND+1, ND+1)
    wdtrain_contrib = params.G * jnp.maximum(0.0, elij_term) * MENT  # (ND+1, ND+1)
    j_lt_i = jnp.arange(ND + 1)[:, None] < jnp.arange(ND + 1)[None, :]
    WDTRAIN_sum = jnp.sum(jnp.where(j_lt_i, wdtrain_contrib, 0.0), axis=0)  # (ND+1,)
    lvl = jnp.arange(ND + 1)
    WDTRAIN_all = WDTRAIN_base + jnp.where(lvl > 0, WDTRAIN_sum, 0.0)

    # Pre-compute per-level quantities for the downdraft
    SIGP_pad = jnp.concatenate([SIGP, _z1])
    WT_dd = jnp.where(T_pad > 273.0, params.OMTRAIN, params.OMTSNOW)  # (ND+1,)
    COEFF_dd = jnp.where(T_pad > 273.0, params.COEFFR, params.COEFFS)  # (ND+1,)
    SIGT_dd = jnp.clip(SIGP_pad, 0.0, 1.0)  # (ND+1,)

    def downdraft_step(carry, i_rev):
        """Scan from INB down to 0 (fixed-length for JIT; invalid steps are no-ops)."""
        (WATER_above, WT_above, MP_above, QP_above, UP_above, VP_above,
         JTT, MP_JTT, P_JTT) = carry
        i_raw = INB - i_rev
        valid = (i_raw >= 0) & (i_raw <= INB)
        i = jnp.maximum(i_raw, 0)  # safe index for array access

        # AFAC depends on QP_above (= QP[i+1] from previous step)
        SIGT_i = SIGT_dd[i]
        COEFF_i = COEFF_dd[i]
        WT_i = WT_dd[i]
        AFAC = jnp.maximum(
            COEFF_i * PH[i] * (QS_pad[i] - 0.5 * (Q_pad[i] + QP_above))
            / (1.0e4 + 2.0e3 * PH[i] * QS_pad[i]),
            0.0,
        )
        B6 = 100.0 * (PH[i] - PH[i + 1]) * SIGT_i * AFAC / WT_i
        C6 = (WATER_above * WT_above + WDTRAIN_all[i] / params.SIGD) / WT_i
        REVAP = 0.5 * (-B6 + jnp.sqrt(jnp.maximum(B6 * B6 + 4.0 * C6, 1e-30)))
        EVAP_i = SIGT_i * AFAC * REVAP
        WATER_i = REVAP * REVAP

        # MP[i]: smoothed downdraft mass flux (only when i > 0)
        im1 = jnp.maximum(i - 1, 0)
        DHDP = jnp.maximum((H[i] - H[im1]) / jnp.maximum(P[im1] - P[i], 1e-30), 10.0)
        MP_raw = 100.0 * GINV * LV[i] * params.SIGD * EVAP_i / DHDP
        FAC = 20.0 / jnp.maximum(PH[im1] - PH[i], 1e-30)
        MP_smooth = (FAC * MP_above + MP_raw) / (1.0 + FAC)

        # Boundary-layer correction: JTT tracks first BL level
        in_bl = (i > 0) & (P[i] > 0.949 * P[0])
        first_bl = in_bl & (JTT == 0)
        JTT_new = jnp.where(first_bl, jnp.int32(i), JTT)
        MP_JTT_new = jnp.where(first_bl, MP_smooth, MP_JTT)
        P_JTT_new = jnp.where(first_bl, P[i], P_JTT)
        MP_bl = MP_JTT_new * (P[0] - P[i]) / jnp.maximum(P[0] - P_JTT_new, 1e-30)
        MP_i = jnp.where(i > 0, jnp.where(in_bl, MP_bl, MP_smooth), 0.0)

        # QP/UP/VP updates (when i != INB)
        not_inb = (i != INB)
        mp_increasing = MP_i > MP_above
        safe_MP_i = jnp.where(MP_i > 0, MP_i, 1.0)  # gradient-safe denom
        RAT = jnp.where(mp_increasing & (MP_i > 0), MP_above / safe_MP_i, 1.0)

        ip1 = jnp.minimum(i + 1, ND - 1)
        QP_inc = (
            QP_above * RAT + Q_pad[i] * (1.0 - RAT)
            + 100.0 * GINV * params.SIGD
            * (PH[i] - PH[i + 1]) * EVAP_i / safe_MP_i
        )
        QP_dec = (
            (GZ[i + 1] - GZ[i]
             + QP_above * (LV[i + 1] + T_pad[ip1] * (params.CL - params.CPD))
             + params.CPD * (T_pad[ip1] - T_pad[i]))
            / jnp.maximum(LV[i] + T_pad[i] * (params.CL - params.CPD), 1e-30)
        )
        QP_i = jnp.where(
            not_inb & mp_increasing, QP_inc,
            jnp.where(not_inb & (MP_above > 0.0), QP_dec, QP[i]),
        )
        QSTM_i = jnp.where(i == 0, QS_pad[0], QS_pad[im1])
        QP_i = jnp.where(not_inb, jnp.clip(QP_i, 0.0, QSTM_i), QP[i])

        UP_i = jnp.where(
            not_inb & mp_increasing,
            UP_above * RAT + U_pad[i] * (1.0 - RAT),
            jnp.where(not_inb & (MP_above > 0.0), UP_above, UP[i])
        )
        VP_i = jnp.where(
            not_inb & mp_increasing,
            VP_above * RAT + V_pad[i] * (1.0 - RAT),
            jnp.where(not_inb & (MP_above > 0.0), VP_above, VP[i])
        )

        # Gate carry and outputs by validity (no-op for out-of-range steps)
        new_carry = (
            jnp.where(valid, WATER_i, WATER_above),
            jnp.where(valid, WT_i, WT_above),
            jnp.where(valid, MP_i, MP_above),
            jnp.where(valid, QP_i, QP_above),
            jnp.where(valid, UP_i, UP_above),
            jnp.where(valid, VP_i, VP_above),
            jnp.where(valid, JTT_new, JTT),
            jnp.where(valid, MP_JTT_new, MP_JTT),
            jnp.where(valid, P_JTT_new, P_JTT),
        )
        outputs = (
            jnp.where(valid, WATER_i, 0.0),
            jnp.where(valid, WT_i, 0.0),
            jnp.where(valid, MP_i, 0.0),
            jnp.where(valid, EVAP_i, 0.0),
            jnp.where(valid, QP_i, 0.0),
            jnp.where(valid, UP_i, 0.0),
            jnp.where(valid, VP_i, 0.0),
        )
        return new_carry, outputs

    ep_active = EP_pad[INB] >= 0.0001
    # Always run the scan (unconditional for JIT); gate results by ep_active
    inb_p1 = jnp.minimum(INB + 1, ND)
    init_carry = (
        jnp.array(0.0),               # WATER[INB+1] = 0
        jnp.array(params.OMTSNOW),    # WT[INB+1] initial
        jnp.array(0.0),               # MP[INB+1] = 0
        QP[inb_p1],                   # QP[INB+1]
        UP[inb_p1],                   # UP[INB+1]
        VP[inb_p1],                   # VP[INB+1]
        jnp.array(0, dtype=jnp.int32),  # JTT
        jnp.array(0.0),               # MP_JTT
        P[0],                          # P_JTT (safe default)
    )
    _, scan_out = jax.lax.scan(downdraft_step, init_carry, jnp.arange(NL + 1))
    WATER_sc, WT_sc, MP_sc, EVAP_sc, QP_sc, UP_sc, VP_sc = scan_out

    # Scatter scan results into full arrays (static loop with mask)
    for k in range(NL + 1):
        level_raw = INB - k
        level = jnp.maximum(level_raw, 0)  # safe index
        k_valid = (level_raw >= 0) & ep_active
        WATER = WATER.at[level].set(jnp.where(k_valid, WATER_sc[k], WATER[level]))
        WT = WT.at[level].set(jnp.where(k_valid, WT_sc[k], WT[level]))
        MP = MP.at[level].set(jnp.where(k_valid, MP_sc[k], MP[level]))
        EVAP = EVAP.at[level].set(jnp.where(k_valid, EVAP_sc[k], EVAP[level]))
        QP = QP.at[level].set(jnp.where(k_valid, QP_sc[k], QP[level]))
        UP = UP.at[level].set(jnp.where(k_valid, UP_sc[k], UP[level]))
        VP = VP.at[level].set(jnp.where(k_valid, VP_sc[k], VP[level]))

    PRECIP = jnp.where(ep_active,
        PRECIP + WT_sc[INB] * params.SIGD * WATER_sc[INB]
        * 3600.0 * 24000.0 / (params.ROWL * params.G),
        PRECIP)

    # --- Blocks 13–15: Tendency accumulation (differentiable) ---
    WD = params.BETA * jnp.abs(MP[ICB]) * 0.01 * params.RD * T_pad[ICB] / (params.SIGD * P[ICB])
    QPRIME = 0.5 * (QP[0] - Q_pad[0])
    TPRIME = params.LV0 * QPRIME / params.CPD

    # --- Block 13: FT, FQ, FU, FV computation ---
    # Precompute DPINV for all ND levels: 0.01 / (PH[i] - PH[i+1])
    DPINV_all = 0.01 / (PH[:ND] - PH[1:ND + 1])  # (ND,)
    CPINV_all = 1.0 / jnp.maximum(CPN[:ND], 1e-30)  # (ND,) safe for levels > NL

    # --- Surface level (i=0) ---
    DPINV_0 = DPINV_all[0]
    # AM = sum of M[1:INB+1] when NK==0
    lvl_m = jnp.arange(ND + 1)
    AM = jnp.where(NK == 0, jnp.sum(jnp.where((lvl_m >= 1) & (lvl_m <= INB), M, 0.0)), 0.0)

    IFLAG = jnp.where((2.0 * params.G * DPINV_0 * AM) >= DELTI, 4, IFLAG)

    FT = FT.at[0].add(
        params.G * DPINV_0 * AM * (T_pad[1] - T_pad[0] + (GZ[1] - GZ[0]) / CPN[0])
        - LVCP[0] * params.SIGD * EVAP[0]
        + params.SIGD * WT[1] * (params.CL - params.CPD) * WATER[1]
        * (T_pad[1] - T_pad[0]) * DPINV_0 / CPN[0]
    )
    FQ = FQ.at[0].add(
        params.G * MP[1] * (QP[1] - Q_pad[0]) * DPINV_0
        + params.SIGD * EVAP[0]
        + params.G * AM * (Q_pad[1] - Q_pad[0]) * DPINV_0
    )
    FU = FU.at[0].add(params.G * DPINV_0 * (MP[1] * (UP[1] - U_pad[0]) + AM * (U_pad[1] - U_pad[0])))
    FV = FV.at[0].add(params.G * DPINV_0 * (MP[1] * (VP[1] - V_pad[0]) + AM * (V_pad[1] - V_pad[0])))

    # MENT column sums for surface level: sum over j=1..INB of MENT[j,0]*(XENT[j,0]-X[0])
    j_mask_surf = (lvl_m >= 1) & (lvl_m <= INB)  # (ND+1,)
    FQ = FQ.at[0].add(params.G * DPINV_0 * jnp.sum(
        jnp.where(j_mask_surf, MENT[:, 0] * (QENT[:, 0] - Q_pad[0]), 0.0)))
    FU = FU.at[0].add(params.G * DPINV_0 * jnp.sum(
        jnp.where(j_mask_surf, MENT[:, 0] * (UENT[:, 0] - U_pad[0]), 0.0)))
    FV = FV.at[0].add(params.G * DPINV_0 * jnp.sum(
        jnp.where(j_mask_surf, MENT[:, 0] * (VENT[:, 0] - V_pad[0]), 0.0)))

    # --- Interior levels (i=1..INB) ---
    # Precompute AMP1[i] and AD[i] for all levels i=1..INB as arrays.
    #
    # AMP1[i] = sum of M[ii] for ii in [i+1, INB+1] (when i >= NK)
    #         + sum of MENT[k, jj] for k in [0..i], jj in [i+1, INB+1]
    #
    # AD[i] = sum of MENT[jj, k] for k in [0..i-1], jj in [i, INB]

    for i in range(1, ND):
        i_active = (i <= INB)  # mask: only compute when i <= INB
        DPINV_i = DPINV_all[i]
        CPINV_i = CPINV_all[i]

        # AMP1: mass flux part
        AMP1_mass = jnp.where(i >= NK,
            jnp.sum(jnp.where((lvl_m >= i + 1) & (lvl_m <= INB + 1), M, 0.0)),
            0.0)
        # AMP1: entrainment part — sum MENT[k, jj] for k<=i, jj>i, jj<=INB+1
        k_le_i = (lvl_m[:, None] <= i)           # (ND+1, 1) broadcast
        jj_gt_i = (lvl_m[None, :] >= i + 1) & (lvl_m[None, :] <= INB + 1)  # (1, ND+1)
        AMP1_ent = jnp.sum(jnp.where(k_le_i & jj_gt_i, MENT, 0.0))
        AMP1 = AMP1_mass + AMP1_ent

        IFLAG = jnp.where(i_active & ((2.0 * params.G * DPINV_i * AMP1) >= DELTI), 4, IFLAG)

        # AD: sum MENT[jj, k] for k < i, jj in [i, INB]
        k_lt_i = (lvl_m[None, :] < i)            # columns k < i
        jj_ge_i = (lvl_m[:, None] >= i) & (lvl_m[:, None] <= INB)  # rows jj >= i, jj <= INB
        AD = jnp.sum(jnp.where(jj_ge_i & k_lt_i, MENT, 0.0))

        # FT: advection + evaporation + self-entrainment + water loading
        ft_adv = (
            params.G * DPINV_i * (
                AMP1 * (T_pad[i + 1] - T_pad[i] + (GZ[i + 1] - GZ[i]) * CPINV_i)
                - AD * (T_pad[i] - T_pad[i - 1] + (GZ[i] - GZ[i - 1]) * CPINV_i)
            )
            - params.SIGD * LVCP[i] * EVAP[i]
        )
        ft_self = (
            params.G * DPINV_i * MENT[i, i]
            * (HP[i] - H[i] + T_pad[i] * (params.CPV - params.CPD) * (Q_pad[i] - QENT[i, i]))
            * CPINV_i
        )
        ft_water = (
            params.SIGD * WT[i + 1] * (params.CL - params.CPD)
            * WATER[i + 1] * (T_pad[i + 1] - T_pad[i])
            * DPINV_i * CPINV_i
        )
        FT = FT.at[i].add(jnp.where(i_active, ft_adv + ft_self + ft_water, 0.0))

        # FQ, FU, FV: advection terms
        fq_adv = params.G * DPINV_i * (AMP1 * (Q_pad[i + 1] - Q_pad[i]) - AD * (Q_pad[i] - Q_pad[i - 1]))
        fu_adv = params.G * DPINV_i * (AMP1 * (U_pad[i + 1] - U_pad[i]) - AD * (U_pad[i] - U_pad[i - 1]))
        fv_adv = params.G * DPINV_i * (AMP1 * (V_pad[i + 1] - V_pad[i]) - AD * (V_pad[i] - V_pad[i - 1]))

        # MENT contributions for k < i (with AWAT water term in FQ)
        k_below = (lvl_m < i)  # (ND+1,)
        AWAT_vec = jnp.maximum(ELIJ[:, i] - (1.0 - EP_pad[i]) * CLW_pad[i], 0.0)  # (ND+1,)
        fq_ment_below = params.G * DPINV_i * jnp.sum(
            jnp.where(k_below, MENT[:, i] * (QENT[:, i] - AWAT_vec - Q_pad[i]), 0.0))
        fu_ment_below = params.G * DPINV_i * jnp.sum(
            jnp.where(k_below, MENT[:, i] * (UENT[:, i] - U_pad[i]), 0.0))
        fv_ment_below = params.G * DPINV_i * jnp.sum(
            jnp.where(k_below, MENT[:, i] * (VENT[:, i] - V_pad[i]), 0.0))

        # MENT contributions for k >= i (no AWAT)
        k_above = (lvl_m >= i) & (lvl_m <= INB)  # (ND+1,)
        fq_ment_above = params.G * DPINV_i * jnp.sum(
            jnp.where(k_above, MENT[:, i] * (QENT[:, i] - Q_pad[i]), 0.0))
        fu_ment_above = params.G * DPINV_i * jnp.sum(
            jnp.where(k_above, MENT[:, i] * (UENT[:, i] - U_pad[i]), 0.0))
        fv_ment_above = params.G * DPINV_i * jnp.sum(
            jnp.where(k_above, MENT[:, i] * (VENT[:, i] - V_pad[i]), 0.0))

        # Downdraft mass flux contributions
        fq_down = (
            params.SIGD * EVAP[i]
            + params.G * (MP[i + 1] * (QP[i + 1] - Q_pad[i]) - MP[i] * (QP[i] - Q_pad[i - 1]))
            * DPINV_i
        )
        fu_down = (
            params.G * (MP[i + 1] * (UP[i + 1] - U_pad[i]) - MP[i] * (UP[i] - U_pad[i - 1]))
            * DPINV_i
        )
        fv_down = (
            params.G * (MP[i + 1] * (VP[i + 1] - V_pad[i]) - MP[i] * (VP[i] - V_pad[i - 1]))
            * DPINV_i
        )

        FQ = FQ.at[i].add(jnp.where(i_active, fq_adv + fq_ment_below + fq_ment_above + fq_down, 0.0))
        FU = FU.at[i].add(jnp.where(i_active, fu_adv + fu_ment_below + fu_ment_above + fu_down, 0.0))
        FV = FV.at[i].add(jnp.where(i_active, fv_adv + fv_ment_below + fv_ment_above + fv_down, 0.0))

    # --- Block 14: FRAC smoothing at INB (differentiable) ---
    ph_ratio = (PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB])
    FQOLD = FQ[INB]
    FQ = FQ.at[INB].set(FQOLD * (1.0 - FRAC))
    FQ = FQ.at[INB - 1].add(FRAC * FQOLD * ph_ratio * LV[INB] / LV[INB - 1])
    FTOLD = FT[INB]
    FT = FT.at[INB].set(FTOLD * (1.0 - FRAC))
    FT = FT.at[INB - 1].add(FRAC * FTOLD * ph_ratio * CPN[INB] / CPN[INB - 1])
    FUOLD = FU[INB]
    FU = FU.at[INB].set(FUOLD * (1.0 - FRAC))
    FU = FU.at[INB - 1].add(FRAC * FUOLD * ph_ratio)
    FVOLD = FV[INB]
    FV = FV.at[INB].set(FVOLD * (1.0 - FRAC))
    FV = FV.at[INB - 1].add(FRAC * FVOLD * ph_ratio)

    # --- Block 15: Conservation correction (differentiable) ---
    active_levels = jnp.arange(ND) <= INB  # (ND,) mask for levels 0..INB
    dp = PH[:ND] - PH[1:ND + 1]  # (ND,)
    DENOM2 = PH[0] - PH[INB + 1]
    ENTS = jnp.sum(jnp.where(active_levels, (CPN[:ND] * FT + LV[:ND] * FQ) * dp, 0.0)) / DENOM2
    UAV = jnp.sum(jnp.where(active_levels, FU * dp, 0.0)) / DENOM2
    VAV = jnp.sum(jnp.where(active_levels, FV * dp, 0.0)) / DENOM2
    FT = jnp.where(active_levels, FT - ENTS * CPINV_all, FT)
    FU = jnp.where(active_levels, (1.0 - params.CU) * (FU - UAV), FU)
    FV = jnp.where(active_levels, (1.0 - params.CU) * (FV - VAV), FV)

    # Gate all outputs by `active` flag — zero out if any early exit was triggered
    zero = 0.0
    FT = jnp.where(active, FT, zero)
    FQ = jnp.where(active, FQ, zero)
    FU = jnp.where(active, FU, zero)
    FV = jnp.where(active, FV, zero)
    PRECIP = jnp.where(active, PRECIP, zero)
    WD = jnp.where(active, WD, zero)
    TPRIME = jnp.where(active, TPRIME, zero)
    QPRIME = jnp.where(active, QPRIME, zero)
    CBMF = jnp.where(active, CBMF, zero)
    OUTCAPE = jnp.where(active, OUTCAPE, zero)

    return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG


def _jax_vectorized_convect(
    t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, params,
):
    def column_call(T, Q, QS, U, V, P, PH, CBMF):
        return _convect_functional_jax(T, Q, QS, U, V, P, PH, ND, NL, NTRA, DELT, CBMF, params)

    # vmap over columns: input arrays are (nlev, ncol), vmap over axis 1 → per-column calls
    # cbmf is (ncol,), vmap over axis 0
    results = jax.vmap(column_call, in_axes=(1, 1, 1, 1, 1, 1, 1, 0))(
        t, q, qs, u, v, p, ph, cbmf)

    # vmap returns (ncol, nlev) for 1D outputs → transpose to (nlev, ncol)
    ft, fq, fu, fv = results[0].T, results[1].T, results[2].T, results[3].T
    precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results[4:]
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
