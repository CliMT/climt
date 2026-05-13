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
from ..._core.condensibles import CondensibleParams, get_condensible_params, _lcl_pressure, _sat_vap_pressure

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(x, **kwargs):
        return x


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
    RD: float
    G: float
    DELT0: float


class EmanuelConvectionPython(ImplicitTendencyComponent):
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
        self.RD = 287.04
        self.G = 9.8
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
            RD=float(self.RD),
            G=float(self.G),
            DELT0=float(self.DELT0),
        )
        self._cond = get_condensible_params()
        super(EmanuelConvectionPython, self).__init__(**kwargs)

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
        from climt._core import bolton_q_sat

        qs = bolton_q_sat(t, p * 100, self.RD, self.RV)
        cbmf = state.get("cloud_base_mass_flux", np.zeros(ncol)).copy()
        ntra = 0
        tra = np.zeros((nlev, 1))
        delt = timestep.total_seconds()
        tra_vector = np.broadcast_to(tra[:, :, np.newaxis], (nlev, 1, ncol))
        results = _numpy_vectorized_convect(
            t,
            q,
            qs,
            u,
            v,
            p,
            ph,
            nlev,
            nlev - 3,
            ntra,
            delt,
            cbmf,
            tra_vector,
            self._params,
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
def _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params, cond):
    CPVMCL = cond.CL - cond.CPV
    EPS = params.RD / cond.RV
    EPSI = 1.0 / EPS
    AH0 = (
        (params.CPD * (1.0 - Q[NK]) + cond.CL * Q[NK]) * T[NK]
        + Q[NK] * (cond.LV0 - CPVMCL * (T[NK] - cond.T_freeze))
        + GZ[NK]
    )
    CPP = params.CPD * (1.0 - Q[NK]) + Q[NK] * cond.CPV
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
        ALV = cond.LV0 - CPVMCL * (T[i] - cond.T_freeze)
        for j in range(2):
            S = 1.0 / (params.CPD + ALV * ALV * QG / (cond.RV * T[i] * T[i]))
            AHG = (
                params.CPD * TG
                + (cond.CL - params.CPD) * Q[NK] * T[i]
                + ALV * QG
                + GZ[i]
            )
            TG = max(TG + S * (AH0 - AHG), 35.0)
            ES = _sat_vap_pressure(TG, cond.species_id)
            QG = EPS * ES / (P[i] - ES * (1.0 - EPS))
        TPK[i] = (
            AH0 - (cond.CL - params.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG
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
    CPVMCL = cond.CL - cond.CPV
    EPS = params.RD / cond.RV
    EPSI = 1.0 / EPS
    GINV = 1.0 / params.G
    DELTI = 1.0 / DELT
    PRECIP = 0.0
    WD = 0.0
    TPRIME = 0.0
    QPRIME = 0.0
    IFLAG = 0
    OUTCAPE = 0.0
    # TH (potential temperature) is computed but never used in this function.
    # It is a carry-over from the original Fortran code. Commented out rather
    # than deleted to preserve the original structure. The 1000 hPa reference
    # pressure is Earth-specific; generalisation deferred to a future pass.
    # TH = np.zeros(NL + 1)
    # for i in range(NL + 1):
    #     RDCP = (params.RD * (1.0 - Q[i]) + Q[i] * params.RV) / (
    #         params.CPD * (1.0 - Q[i]) + Q[i] * params.CPV
    #     )
    #     TH[i] = T[i] * (1000.0 / P[i]) ** RDCP
    GZ = np.zeros(ND + 1)
    CPN = np.zeros(ND + 1)
    H = np.zeros(ND + 1)
    LV = np.zeros(ND + 1)
    HM = np.zeros(ND + 1)
    TV = np.zeros(ND + 1)
    GZ[0] = 0.0
    CPN[0] = params.CPD * (1.0 - Q[0]) + Q[0] * cond.CPV
    H[0] = T[0] * CPN[0]
    LV[0] = cond.LV0 - CPVMCL * (T[0] - cond.T_freeze)
    HM[0] = LV[0] * Q[0]
    TV[0] = T[0] * (1.0 + Q[0] * EPSI - Q[0])
    AHMIN = 1.0e12
    IHMIN = NL
    for i in range(1, NL + 1):
        TVX = T[i] * (1.0 + Q[i] * EPSI - Q[i])
        TVY = T[i - 1] * (1.0 + Q[i - 1] * EPSI - Q[i - 1])
        GZ[i] = GZ[i - 1] + 0.5 * params.RD * (TVX + TVY) * (P[i - 1] - P[i]) / PH[i]
        CPN[i] = params.CPD * (1.0 - Q[i]) + cond.CPV * Q[i]
        H[i] = T[i] * CPN[i] + GZ[i]
        LV[i] = cond.LV0 - CPVMCL * (T[i] - cond.T_freeze)
        HM[i] = (
            (params.CPD * (1.0 - Q[i]) + cond.CL * Q[i]) * (T[i] - T[0])
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
    PLCL = _lcl_pressure(P[NK], RH, T[NK], cond)
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
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params, cond
    )
    for i in range(NK, ICB + 1):
        TVP[i] -= TP[i] * Q[NK]
    if CBMF == 0.0 and TVP[ICB] <= (TV[ICB] - params.DTMAX):
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
    if IFLAG != 4:
        IFLAG = 1
    TVP, TP, CLW = _tlift_functional_np(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params, cond
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
        HP[i] = H[NK] + (LV[i] + (params.CPD - cond.CPV) * T[i]) * EP[i] * CLW[i]
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
            BF2 = 1.0 + LV[j] * LV[j] * QS[j] / (cond.RV * T[j] * T[j] * params.CPD)
            ANUM = H[j] - HP[i] + (cond.CPV - params.CPD) * T[j] * (QTI - Q[j])
            DENOM = H[i] - HP[i] + (params.CPD - cond.CPV) * (Q[i] - QTI) * T[j]
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
            COEFF = params.COEFFR if T[i] > cond.T_freeze else params.COEFFS
            WT[i] = params.OMTRAIN if T[i] > cond.T_freeze else params.OMTSNOW
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
                        + QP[i + 1] * (LV[i + 1] + T[i + 1] * (cond.CL - params.CPD))
                        + params.CPD * (T[i + 1] - T[i])
                    ) / (LV[i] + T[i] * (cond.CL - params.CPD))
                    UP[i] = UP[i + 1]
                    VP[i] = VP[i + 1]
                    for k in range(NTRA):
                        TRAP[i, k] = TRAP[i + 1, k]
                QP[i] = max(min(QP[i], QSTM), 0.0)
        PRECIP += (
            WT[0] * params.SIGD * WATER[0] * 3600.0 * 24000.0 / (cond.ROWL * params.G)
        )
    WD = params.BETA * abs(MP[ICB]) * 0.01 * params.RD * T[ICB] / (params.SIGD * P[ICB])
    QPRIME = 0.5 * (QP[0] - Q[0])
    TPRIME = cond.LV0 * QPRIME / params.CPD
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
            * (HP[i] - H[i] + T[i] * (cond.CPV - params.CPD) * (Q[i] - QENT[i, i]))
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


EmanuelConvectionPythonV3 = EmanuelConvectionPython
