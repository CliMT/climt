# -*- coding: utf-8 -*-
import numpy as np
from typing import NamedTuple
from ..._core.backend import jit_compile, prange


def set_item(arr, idx, val):
    arr[idx] = val
    return arr

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x: x

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

class EmanuelConvectionPythonV2(object):
    def __init__(self, **kwargs):
        # Default parameters from Fortran code
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

        # Update with kwargs
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
        
        # Create a NamedTuple for Numba compatibility
        self._params = EmanuelParams(
            IPBL=int(self.IPBL), MINORIG=int(self.MINORIG), ELCRIT=float(self.ELCRIT),
            TLCRIT=float(self.TLCRIT), ENTP=float(self.ENTP), SIGD=float(self.SIGD),
            SIGS=float(self.SIGS), OMTRAIN=float(self.OMTRAIN), OMTSNOW=float(self.OMTSNOW),
            COEFFR=float(self.COEFFR), COEFFS=float(self.COEFFS), CU=float(self.CU),
            BETA=float(self.BETA), DTMAX=float(self.DTMAX), ALPHA=float(self.ALPHA),
            DAMP=float(self.DAMP), CPD=float(self.CPD), CPV=float(self.CPV),
            CL=float(self.CL), RV=float(self.RV), RD=float(self.RD),
            LV0=float(self.LV0), G=float(self.G), ROWL=float(self.ROWL),
            DELT0=float(self.DELT0)
        )

    def array_call(self, state, timestep):
        t = state['air_temperature']
        q = state['specific_humidity']
        u = state['eastward_wind']
        v = state['northward_wind']
        p = state['air_pressure']
        ph = state['air_pressure_on_interface_levels']

        nlev, ncol = t.shape

        from climt._core import bolton_q_sat
        qs = bolton_q_sat(t, p * 100, self.RD, self.RV)

        cbmf = state.get('cloud_base_mass_flux', np.zeros(ncol)).copy()
        ntra = 0; tra = np.zeros((nlev, 1)); delt = timestep.total_seconds()

        # Call vectorized function
        tra_vector = np.broadcast_to(tra[:, :, np.newaxis], (nlev, 1, ncol))

        results = _numpy_vectorized_convect(
            t, q, qs, u, v, p, ph, nlev, nlev-3, ntra, delt, cbmf, tra_vector, self._params,
            _convect_functional_numba, _tlift_functional_numba
        )
        
        ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results

        tendencies = {'air_temperature': ft, 'specific_humidity': fq, 'eastward_wind': fu, 'northward_wind': fv}
        diagnostics = {
            'convective_state': iflag, 'convective_precipitation_rate': precip,
            'convective_downdraft_velocity_scale': wd, 'convective_downdraft_temperature_scale': tprime,
            'convective_downdraft_specific_humidity_scale': qprime, 'cloud_base_mass_flux': cbmf_new,
            'atmosphere_convective_available_potential_energy': outcape
        }
        return tendencies, diagnostics

def _tlift_functional(xp, P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params):
    # ND, NL, NK, KK, ICB are 0-based indices
    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS

    AH0 = (params.CPD * (1. - Q[NK]) + params.CL * Q[NK]) * T[NK] + \
          Q[NK] * (params.LV0 - CPVMCL * (T[NK] - 273.15)) + GZ[NK]
    CPP = params.CPD * (1. - Q[NK]) + Q[NK] * params.CPV
    CPINV = 1.0 / CPP

    if KK == 1:
        for i in range(ICB):
            CLW = set_item(CLW, i, 0.0)
        for i in range(NK, ICB):
            TPK = set_item(TPK, i, T[NK] - (GZ[i] - GZ[NK]) * CPINV)
            TVP = set_item(TVP, i, TPK[i] * (1. + Q[NK] * EPSI))

    NST = ICB
    NSB = ICB
    if KK == 2:
        NST = NL
        NSB = ICB + 1

    for i in range(NSB, NST + 1):
        TG = T[i]
        QG = QS[i]
        ALV = params.LV0 - CPVMCL * (T[i] - 273.15)
        for j in range(2):
            S = params.CPD + ALV * ALV * QG / (params.RV * T[i] * T[i])
            S = 1.0 / S
            AHG = params.CPD * TG + (params.CL - params.CPD) * Q[NK] * T[i] + ALV * QG + GZ[i]
            TG = TG + S * (AH0 - AHG)
            TG = xp.maximum(TG, 35.0)
            TC = TG - 273.15
            DENOM = 243.5 + TC
            if TC >= 0.0:
                ES = 6.112 * xp.exp(17.67 * TC / DENOM)
            else:
                ES = xp.exp(23.33086 - 6111.72784 / TG + 0.15215 * xp.log(TG))
            QG = EPS * ES / (P[i] - ES * (1. - EPS))

        TPK = set_item(TPK, i, (AH0 - (params.CL - params.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG) / params.CPD)
        CLW = set_item(CLW, i, Q[NK] - QG)
        CLW = set_item(CLW, i, xp.maximum(0.0, CLW[i]))
        RG = QG / (1. - Q[NK])
        TVP = set_item(TVP, i, TPK[i] * (1. + RG * EPSI))

    return TVP, TPK, CLW

def _convect_functional(xp, tlift_func, T_in, Q_in, QS_in, U_in, V_in, P_in, PH_in, ND, NL, NTRA, DELT, CBMF_in, TRA_in, params):
    T = T_in.copy()
    Q = Q_in.copy()
    QS = QS_in.copy()
    U = U_in.copy()
    V = V_in.copy()
    P = P_in.copy()
    PH = PH_in.copy()
    TRA = TRA_in.copy()
    CBMF = CBMF_in

    FT = xp.zeros(ND)
    FQ = xp.zeros(ND)
    FU = xp.zeros(ND)
    FV = xp.zeros(ND)
    FTRA = xp.zeros((ND, xp.maximum(1, NTRA)))

    CPVMCL = params.CL - params.CPV
    EPS = params.RD / params.RV
    EPSI = 1.0 / EPS
    GINV = 1.0 / params.G
    DELTI = 1.0 / DELT

    TH = xp.zeros(NL + 1)
    for i in range(NL + 1):
        RDCP = (params.RD * (1. - Q[i]) + Q[i] * params.RV) / (params.CPD * (1. - Q[i]) + Q[i] * params.CPV)
        TH = set_item(TH, i, T[i] * (1000.0 / P[i])**RDCP)

    PRECIP = 0.0
    WD = 0.0
    TPRIME = 0.0
    QPRIME = 0.0
    IFLAG = 0
    OUTCAPE = 0.0

    GZ = xp.zeros(ND + 1)
    CPN = xp.zeros(ND + 1)
    H = xp.zeros(ND + 1)
    LV = xp.zeros(ND + 1)
    HM = xp.zeros(ND + 1)
    TV = xp.zeros(ND + 1)

    GZ = set_item(GZ, 0, 0.0)
    CPN = set_item(CPN, 0, params.CPD * (1. - Q[0]) + Q[0] * params.CPV)
    H = set_item(H, 0, T[0] * CPN[0])
    LV = set_item(LV, 0, params.LV0 - CPVMCL * (T[0] - 273.15))
    HM = set_item(HM, 0, LV[0] * Q[0])
    TV = set_item(TV, 0, T[0] * (1. + Q[0] * EPSI - Q[0]))

    AHMIN = 1.0E12
    IHMIN = NL

    for i in range(1, NL + 1):
        TVX = T[i] * (1. + Q[i] * EPSI - Q[i])
        TVY = T[i-1] * (1. + Q[i-1] * EPSI - Q[i-1])
        GZ = set_item(GZ, i, GZ[i-1] + 0.5 * params.RD * (TVX + TVY) * (P[i-1] - P[i]) / PH[i])
        CPN = set_item(CPN, i, params.CPD * (1. - Q[i]) + params.CPV * Q[i])
        H = set_item(H, i, T[i] * CPN[i] + GZ[i])
        LV = set_item(LV, i, params.LV0 - CPVMCL * (T[i] - 273.15))
        HM = set_item(HM, i, (params.CPD * (1. - Q[i]) + params.CL * Q[i]) * (T[i] - T[0]) + LV[i] * Q[i] + GZ[i])
        TV = set_item(TV, i, T[i] * (1. + Q[i] * EPSI - Q[i]))
        if (i+1) >= params.MINORIG and HM[i] < AHMIN and HM[i] < HM[i-1]:
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
        IFLAG = 0
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

    TVP = xp.zeros(ND)
    TP = xp.zeros(ND)
    CLW = xp.zeros(ND)
    TVP, TP, CLW = tlift_func(xp, P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params)
    for i in range(NK, ICB + 1):
        TVP = set_item(TVP, i, TVP[i] - TP[i] * Q[NK])

    if CBMF == 0.0 and TVP[ICB] <= (TV[ICB] - params.DTMAX):
        IFLAG = 0
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG
            
    if IFLAG != 4:
        IFLAG = 1

    TVP, TP, CLW = tlift_func(xp, P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params)

    EP = xp.zeros(ND)
    SIGP = xp.zeros(ND)
    for i in range(NK + 1):
        EP = set_item(EP, i, 0.0)
        SIGP = set_item(SIGP, i, params.SIGS)
    for i in range(NK + 1, NL + 1):
        TCA = TP[i] - 273.15
        ELACRIT = params.ELCRIT if TCA >= 0.0 else params.ELCRIT * (1.0 - TCA / params.TLCRIT)
        ELACRIT = xp.maximum(ELACRIT, 0.0)
        EPMAX = 0.999
        EP = set_item(EP, i, EPMAX * (1.0 - ELACRIT / xp.maximum(CLW[i], 1.0E-8)))
        EP = set_item(EP, i, xp.maximum(xp.minimum(EP[i], EPMAX), 0.0))
        SIGP = set_item(SIGP, i, params.SIGS)

    for i in range(ICB + 1, NL + 1):
        TVP = set_item(TVP, i, TVP[i] - TP[i] * Q[NK])

    HP = H.copy()
    NENT = xp.zeros(ND + 1, dtype=xp.int32)
    WATER = xp.zeros(ND + 1)
    EVAP = xp.zeros(ND + 1)
    WT = xp.full(ND + 1, params.OMTSNOW)
    MP = xp.zeros(ND + 1)
    M = xp.zeros(ND + 1)
    LVCP = xp.zeros(ND + 1)
    LVCP = set_item(LVCP, slice(0, NL+1), LV[:NL+1] / CPN[:NL+1])

    QENT = xp.zeros((ND + 1, ND + 1))
    ELIJ = xp.zeros((ND + 1, ND + 1))
    MENT = xp.zeros((ND + 1, ND + 1))
    SIJ = xp.zeros((ND + 1, ND + 1))
    UENT = xp.zeros((ND + 1, ND + 1))
    VENT = xp.zeros((ND + 1, ND + 1))
    TRAENT = xp.zeros((ND + 1, ND + 1, xp.maximum(1, NTRA)))

    for i in range(NL + 1):
        for j in range(NL + 1):
            QENT = set_item(QENT, (i, j), Q[j])
            UENT = set_item(UENT, (i, j), U[j])
            VENT = set_item(VENT, (i, j), V[j])
            for k in range(NTRA):
                TRAENT = set_item(TRAENT, (i, j, k), TRA[j, k])

    QP = xp.zeros(ND + 1)
    UP = xp.zeros(ND + 1)
    VP = xp.zeros(ND + 1)
    TRAP = xp.zeros((ND + 1, xp.maximum(1, NTRA)))

    QP = set_item(QP, 0, Q[0]); UP = set_item(UP, 0, U[0]); VP = set_item(VP, 0, V[0])
    for i in range(NTRA): TRAP = set_item(TRAP, (0, i), TRA[0, i])
    for i in range(1, NL + 1):
        QP = set_item(QP, i, Q[i-1]); UP = set_item(UP, i, U[i-1]); VP = set_item(VP, i, V[i-1])
        for j in range(NTRA): TRAP = set_item(TRAP, (i, j), TRA[i-1, j])

    CAPE = 0.0; CAPEM = 0.0; INB = ICB + 1; INB1 = INB; BYP = 0.0
    for i in range(ICB + 1, NL):
        BY = (TVP[i] - TV[i]) * (PH[i] - PH[i+1]) / P[i]
        CAPE += BY
        if BY >= 0.0: INB1 = i + 1
        if CAPE > 0.0:
            INB = i + 1
            BYP = (TVP[i+1] - TV[i+1]) * (PH[i+1] - PH[i+2]) / P[i+1]
            CAPEM = CAPE
    INB = max(INB, INB1)
    CAPE = CAPEM + BYP
    DEFRAC = xp.maximum(CAPEM - CAPE, 0.001)
    FRAC = xp.minimum(xp.maximum(-CAPE / DEFRAC, 0.0), 1.0)
    OUTCAPE = CAPE

    for i in range(ICB, INB + 1):
        HP = set_item(HP, i, H[NK] + (LV[i] + (params.CPD - params.CPV) * T[i]) * EP[i] * CLW[i])

    TVPPLCL = TVP[ICB-1] - params.RD * TVP[ICB-1] * (P[ICB-1] - PLCL) / (CPN[ICB-1] * P[ICB-1])
    TVAPLCL = TV[ICB] + (TVP[ICB] - TVP[ICB+1]) * (PLCL - P[ICB]) / (P[ICB] - P[ICB+1])
    DTPBL = 0.0
    for i in range(NK, ICB):
        DTPBL += (TVP[i] - TV[i]) * (PH[i] - PH[i+1])
    DTPBL /= (PH[NK] - PH[ICB])
    DTMA = TVPPLCL - TVAPLCL + params.DTMAX + DTPBL

    CBMFOLD = CBMF
    DAMPS = params.DAMP * DELT / params.DELT0
    CBMF = (1. - DAMPS) * CBMF + 0.1 * params.ALPHA * DTMA
    CBMF = xp.maximum(CBMF, 0.0)

    if CBMF == 0.0 and CBMFOLD == 0.0:
        return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

    M = set_item(M, ICB, 0.0); DBOSUM = 0.0
    for i in range(ICB + 1, INB + 1):
        k = min(i, INB1)
        DBO = xp.abs(TV[k] - TVP[k]) + params.ENTP * 0.02 * (PH[k] - PH[k+1])
        DBOSUM += DBO
        M = set_item(M, i, CBMF * DBO)
    if DBOSUM > 0:
        for i in range(ICB + 1, INB + 1): M = set_item(M, i, M[i] / DBOSUM)

    for i in range(ICB + 1, INB + 1):
        QTI = Q[NK] - EP[i] * CLW[i]
        for j in range(ICB, INB + 1):
            BF2 = 1. + LV[j] * LV[j] * QS[j] / (params.RV * T[j] * T[j] * params.CPD)
            ANUM = H[j] - HP[i] + (params.CPV - params.CPD) * T[j] * (QTI - Q[j])
            DENOM = H[i] - HP[i] + (params.CPD - params.CPV) * (Q[i] - QTI) * T[j]
            DEI = DENOM
            if xp.abs(DEI) < 0.01: DEI = 0.01
            SIJ = set_item(SIJ, (i, j), ANUM / DEI)
            SIJ = set_item(SIJ, (i, i), 1.0)
            ALTEM = (SIJ[i, j] * Q[i] + (1. - SIJ[i, j]) * QTI - QS[j]) / BF2
            CWAT = CLW[j] * (1. - EP[j])
            if (SIJ[i, j] < 0.0 or SIJ[i, j] > 1.0 or ALTEM > CWAT) and j > i:
                ANUM -= LV[j] * (QTI - QS[j] - CWAT * BF2)
                DENOM += LV[j] * (Q[i] - QTI)
                if xp.abs(DENOM) < 0.01: DENOM = 0.01
                SIJ = set_item(SIJ, (i, j), ANUM / DENOM)
                ALTEM = SIJ[i, j] * Q[i] + (1. - SIJ[i, j]) * QTI - QS[j] - (BF2 - 1.) * CWAT

            if 0.0 < SIJ[i, j] < 0.9:
                QENT = set_item(QENT, (i, j), SIJ[i, j] * Q[i] + (1. - SIJ[i, j]) * QTI)
                UENT = set_item(UENT, (i, j), SIJ[i, j] * U[i] + (1. - SIJ[i, j]) * U[NK])
                VENT = set_item(VENT, (i, j), SIJ[i, j] * V[i] + (1. - SIJ[i, j]) * V[NK])
                for k in range(NTRA): TRAENT = set_item(TRAENT, (i, j, k), SIJ[i, j] * TRA[i, k] + (1. - SIJ[i, j]) * TRA[NK, k])
                ELIJ = set_item(ELIJ, (i, j), xp.maximum(0.0, ALTEM))
                MENT = set_item(MENT, (i, j), M[i] / (1. - SIJ[i, j]))
                NENT = set_item(NENT, i, NENT[i] + 1)
            SIJ = set_item(SIJ, (i, j), xp.minimum(xp.maximum(SIJ[i, j], 0.0), 1.0))
        if NENT[i] == 0:
            MENT = set_item(MENT, (i, i), M[i]); QENT = set_item(QENT, (i, i), Q[NK] - EP[i] * CLW[i]); UENT = set_item(UENT, (i, i), U[NK]); VENT = set_item(VENT, (i, i), V[NK])
            for k in range(NTRA): TRAENT = set_item(TRAENT, (i, i, k), TRA[NK, k])
            ELIJ = set_item(ELIJ, (i, i), CLW[i]); SIJ = set_item(SIJ, (i, i), 1.0)
    SIJ = set_item(SIJ, (INB, INB), 1.0)

    for i in range(ICB + 1, INB + 1):
        if NENT[i] != 0:
            QP1 = Q[NK] - EP[i] * CLW[i]
            ANUM = H[i] - HP[i] - LV[i] * (QP1 - QS[i])
            DENOM = H[i] - HP[i] + LV[i] * (Q[i] - QP1)
            if xp.abs(DENOM) < 0.01: DENOM = 0.01
            SCRIT = ANUM / DENOM
            ALT = QP1 - QS[i] + SCRIT * (Q[i] - QP1)
            if ALT < 0.0: SCRIT = 1.0
            SCRIT = xp.maximum(SCRIT, 0.0)
            ASIJ = 0.0; SMIN = 1.0
            for j in range(ICB, INB + 1):
                if 0.0 < SIJ[i, j] < 0.9:
                    if j > i:
                        SMID = xp.minimum(SIJ[i, j], SCRIT)
                        SJMAX = SMID; SJMIN = SMID
                        if SMID < SMIN and SIJ[i, j+1] < SMID:
                            SMIN = SMID; SJMAX = xp.minimum(xp.minimum(SIJ[i, j+1], SIJ[i, j]), SCRIT)
                            SJMIN = xp.minimum(xp.maximum(SIJ[i, j-1], SIJ[i, j]), SCRIT)
                    else:
                        SJMAX = xp.maximum(SIJ[i, j+1], SCRIT); SMID = xp.maximum(SIJ[i, j], SCRIT)
                        SJMIN = xp.maximum(SIJ[i, j-1] if j > 0 else 0.0, SCRIT)
                    ASIJ += (xp.abs(SJMAX - SMID) + xp.abs(SJMIN - SMID)) * (PH[j] - PH[j+1])
                    MENT = set_item(MENT, (i, j), MENT[i, j] * (xp.abs(SJMAX - SMID) + xp.abs(SJMIN - SMID)) * (PH[j] - PH[j+1]))
            ASIJ = xp.maximum(1.0E-21, ASIJ)
            for j in range(ICB, INB + 1): MENT = set_item(MENT, (i, j), MENT[i, j] / ASIJ)
            if xp.sum(MENT[i, ICB:INB+1]) < 1.0E-18:
                NENT = set_item(NENT, i, 0); MENT = set_item(MENT, (i, i), M[i]); QENT = set_item(QENT, (i, i), Q[NK] - EP[i] * CLW[i]); UENT = set_item(UENT, (i, i), U[NK]); VENT = set_item(VENT, (i, i), V[NK])
                for k in range(NTRA): TRAENT = set_item(TRAENT, (i, i, k), TRA[NK, k])
                ELIJ = set_item(ELIJ, (i, i), CLW[i]); SIJ = set_item(SIJ, (i, i), 1.0)

    if EP[INB] >= 0.0001:
        JTT = 0
        for i in range(INB, -1, -1):
            WDTRAIN = params.G * EP[i] * M[i] * CLW[i]
            if i > 0:
                for j in range(i):
                    WDTRAIN += params.G * xp.maximum(0.0, ELIJ[j, i] - (1. - EP[i]) * CLW[i]) * MENT[j, i]
            COEFF = params.COEFFR if T[i] > 273.0 else params.COEFFS
            WT = set_item(WT, i, params.OMTRAIN if T[i] > 273.0 else params.OMTSNOW)
            AFAC = xp.maximum(COEFF * PH[i] * (QS[i] - 0.5 * (Q[i] + QP[i+1])) / (1.0E4 + 2.0E3 * PH[i] * QS[i]), 0.0)
            SIGT = xp.minimum(xp.maximum(SIGP[i], 0.0), 1.0)
            B6 = 100. * (PH[i] - PH[i+1]) * SIGT * AFAC / WT[i]
            C6 = (WATER[i+1] * WT[i+1] + WDTRAIN / params.SIGD) / WT[i]
            REVAP = 0.5 * (-B6 + xp.sqrt(B6 * B6 + 4. * C6))
            EVAP = set_item(EVAP, i, SIGT * AFAC * REVAP)
            WATER = set_item(WATER, i, REVAP * REVAP)
            if i > 0:
                DHDP = xp.maximum((H[i] - H[i-1]) / (P[i-1] - P[i]), 10.0)
                MP = set_item(MP, i, 100. * GINV * LV[i] * params.SIGD * EVAP[i] / DHDP)
                FAC = 20.0 / (PH[i-1] - PH[i])
                MP = set_item(MP, i, (FAC * MP[i+1] + MP[i]) / (1. + FAC))
                if P[i] > (0.949 * P[0]):
                    JTT = max(JTT, i)
                    MP = set_item(MP, i, MP[JTT] * (P[0] - P[i]) / (P[0] - P[JTT]))
            if i != INB:
                QSTM = QS[0] if i == 0 else QS[i-1]
                if MP[i] > MP[i+1]:
                    RAT = MP[i+1] / MP[i]
                    QP = set_item(QP, i, QP[i+1] * RAT + Q[i] * (1.0 - RAT) + 100. * GINV * params.SIGD * (PH[i] - PH[i+1]) * (EVAP[i] / MP[i]))
                    UP = set_item(UP, i, UP[i+1] * RAT + U[i] * (1. - RAT)); VP = set_item(VP, i, VP[i+1] * RAT + V[i] * (1. - RAT))
                    for k in range(NTRA): TRAP = set_item(TRAP, (i, k), TRAP[i+1, k] * RAT + TRAP[i, k] * (1. - RAT))
                elif MP[i+1] > 0.0:
                    QP = set_item(QP, i, (GZ[i+1] - GZ[i] + QP[i+1] * (LV[i+1] + T[i+1] * (params.CL - params.CPD)) + params.CPD * (T[i+1] - T[i])) / (LV[i] + T[i] * (params.CL - params.CPD)))
                    UP = set_item(UP, i, UP[i+1]); VP = set_item(VP, i, VP[i+1])
                    for k in range(NTRA): TRAP = set_item(TRAP, (i, k), TRAP[i+1, k])
                QP = set_item(QP, i, xp.maximum(xp.minimum(QP[i], QSTM), 0.0))
        PRECIP += WT[0] * params.SIGD * WATER[0] * 3600. * 24000. / (params.ROWL * params.G)

    WD = params.BETA * xp.abs(MP[ICB]) * 0.01 * params.RD * T[ICB] / (params.SIGD * P[ICB])
    QPRIME = 0.5 * (QP[0] - Q[0]); TPRIME = params.LV0 * QPRIME / params.CPD
    DPINV = 0.01 / (PH[0] - PH[1]); AM = 0.0
    if NK == 0: AM = xp.sum(M[1:INB+1])
    if (2. * params.G * DPINV * AM) >= DELTI: IFLAG = 4
    FT = set_item(FT, 0, FT[0] + params.G * DPINV * AM * (T[1] - T[0] + (GZ[1] - GZ[0]) / CPN[0]) - LVCP[0] * params.SIGD * EVAP[0] + params.SIGD * WT[1] * (params.CL - params.CPD) * WATER[1] * (T[1] - T[0]) * DPINV / CPN[0])
    FQ = set_item(FQ, 0, FQ[0] + params.G * MP[1] * (QP[1] - Q[0]) * DPINV + params.SIGD * EVAP[0] + params.G * AM * (Q[1] - Q[0]) * DPINV)
    FU = set_item(FU, 0, FU[0] + params.G * DPINV * (MP[1] * (UP[1] - U[0]) + AM * (U[1] - U[0])))
    FV = set_item(FV, 0, FV[0] + params.G * DPINV * (MP[1] * (VP[1] - V[0]) + AM * (V[1] - V[0])))
    for j in range(NTRA): FTRA = set_item(FTRA, (0, j), FTRA[0, j] + params.G * DPINV * (MP[1] * (TRAP[1, j] - TRA[0, j]) + AM * (TRA[1, j] - TRA[0, j])))
    for j in range(1, INB + 1):
        FQ = set_item(FQ, 0, FQ[0] + params.G * DPINV * MENT[j, 0] * (QENT[j, 0] - Q[0]))
        FU = set_item(FU, 0, FU[0] + params.G * DPINV * MENT[j, 0] * (UENT[j, 0] - U[0]))
        FV = set_item(FV, 0, FV[0] + params.G * DPINV * MENT[j, 0] * (VENT[j, 0] - V[0]))
        for k in range(NTRA): FTRA = set_item(FTRA, (0, k), FTRA[0, k] + params.G * DPINV * MENT[j, 0] * (TRAENT[j, 0, k] - TRA[0, k]))

    for i in range(1, INB + 1):
        DPINV = 0.01 / (PH[i] - PH[i+1]); CPINV = 1.0 / CPN[i]
        AMP1 = xp.sum(M[i+1:INB+2]) if i >= NK else 0.0
        for k in range(i + 1): AMP1 += xp.sum(MENT[k, i+1:INB+2])
        if (2. * params.G * DPINV * AMP1) >= DELTI: IFLAG = 4
        AD = 0.0
        for k in range(i): AD += xp.sum(MENT[i:INB+1, k])
        FT = set_item(FT, i, FT[i] + params.G * DPINV * (AMP1 * (T[i+1] - T[i] + (GZ[i+1] - GZ[i]) * CPINV) - AD * (T[i] - T[i-1] + (GZ[i] - GZ[i-1]) * CPINV)) - params.SIGD * LVCP[i] * EVAP[i])
        FT = set_item(FT, i, FT[i] + params.G * DPINV * MENT[i, i] * (HP[i] - H[i] + T[i] * (params.CPV - params.CPD) * (Q[i] - QENT[i, i])) * CPINV)
        FT = set_item(FT, i, FT[i] + params.SIGD * WT[i+1] * (params.CL - params.CPD) * WATER[i+1] * (T[i+1] - T[i]) * DPINV * CPINV)
        FQ = set_item(FQ, i, FQ[i] + params.G * DPINV * (AMP1 * (Q[i+1] - Q[i]) - AD * (Q[i] - Q[i-1])))
        FU = set_item(FU, i, FU[i] + params.G * DPINV * (AMP1 * (U[i+1] - U[i]) - AD * (U[i] - U[i-1])))
        FV = set_item(FV, i, FV[i] + params.G * DPINV * (AMP1 * (V[i+1] - V[i]) - AD * (V[i] - V[i-1])))
        for k in range(NTRA): FTRA = set_item(FTRA, (i, k), FTRA[i, k] + params.G * DPINV * (AMP1 * (TRA[i+1, k] - TRA[i, k]) - AD * (TRA[i, k] - TRA[i-1, k])))
        for k in range(i):
            AWAT = xp.maximum(ELIJ[k, i] - (1. - EP[i]) * CLW[i], 0.0)
            FQ = set_item(FQ, i, FQ[i] + params.G * DPINV * MENT[k, i] * (QENT[k, i] - AWAT - Q[i]))
            FU = set_item(FU, i, FU[i] + params.G * DPINV * MENT[k, i] * (UENT[k, i] - U[i])); FV = set_item(FV, i, FV[i] + params.G * DPINV * MENT[k, i] * (VENT[k, i] - V[i]))
            for j in range(NTRA): FTRA = set_item(FTRA, (i, j), FTRA[i, j] + params.G * DPINV * MENT[k, i] * (TRAENT[k, i, j] - TRA[i, j]))
        for k in range(i, INB + 1):
            FQ = set_item(FQ, i, FQ[i] + params.G * DPINV * MENT[k, i] * (QENT[k, i] - Q[i]))
            FU = set_item(FU, i, FU[i] + params.G * DPINV * MENT[k, i] * (UENT[k, i] - U[i])); FV = set_item(FV, i, FV[i] + params.G * DPINV * MENT[k, i] * (VENT[k, i] - V[i]))
            for j in range(NTRA): FTRA = set_item(FTRA, (i, j), FTRA[i, j] + params.G * DPINV * MENT[k, i] * (TRAENT[k, i, j] - TRA[i, j]))
        FQ = set_item(FQ, i, FQ[i] + params.SIGD * EVAP[i] + params.G * (MP[i+1] * (QP[i+1] - Q[i]) - MP[i] * (QP[i] - Q[i-1])) * DPINV)
        FU = set_item(FU, i, FU[i] + params.G * (MP[i+1] * (UP[i+1] - U[i]) - MP[i] * (UP[i] - U[i-1])) * DPINV)
        FV = set_item(FV, i, FV[i] + params.G * (MP[i+1] * (VP[i+1] - V[i]) - MP[i] * (VP[i] - V[i-1])) * DPINV)
        for j in range(NTRA): FTRA = set_item(FTRA, (i, j), FTRA[i, j] + params.G * DPINV * (MP[i+1] * (TRAP[i+1, j] - TRA[i, j]) - MP[i] * (TRAP[i, j] - TRA[i-1, j])))

    FQOLD = FQ[INB]; FQ = set_item(FQ, INB, FQ[INB] * (1. - FRAC))
    FQ = set_item(FQ, INB-1, FQ[INB-1] + FRAC * FQOLD * ((PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])) * LV[INB] / LV[INB-1])
    FTOLD = FT[INB]; FT = set_item(FT, INB, FT[INB] * (1. - FRAC))
    FT = set_item(FT, INB-1, FT[INB-1] + FRAC * FTOLD * ((PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])) * CPN[INB] / CPN[INB-1])
    FUOLD = FU[INB]; FU = set_item(FU, INB, FU[INB] * (1. - FRAC))
    FU = set_item(FU, INB-1, FU[INB-1] + FRAC * FUOLD * ((PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])))
    FVOLD = FV[INB]; FV = set_item(FV, INB, FV[INB] * (1. - FRAC))
    FV = set_item(FV, INB-1, FV[INB-1] + FRAC * FVOLD * ((PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])))
    for k in range(NTRA):
        FTRAOLD = FTRA[INB, k]; FTRA = set_item(FTRA, (INB, k), FTRA[INB, k] * (1. - FRAC))
        FTRA = set_item(FTRA, (INB-1, k), FTRA[INB-1, k] + FRAC * FTRAOLD * (PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB]))

    ENTS = xp.sum((CPN[:INB+1] * FT[:INB+1] + LV[:INB+1] * FQ[:INB+1]) * (PH[:INB+1] - PH[1:INB+2])) / (PH[0] - PH[INB+1])
    UAV = xp.sum(FU[:INB+1] * (PH[:INB+1] - PH[1:INB+2])) / (PH[0] - PH[INB+1])
    VAV = xp.sum(FV[:INB+1] * (PH[:INB+1] - PH[1:INB+2])) / (PH[0] - PH[INB+1])
    for i in range(INB + 1):
        FT = set_item(FT, i, FT[i] - ENTS / CPN[i])
        FU = set_item(FU, i, (1. - params.CU) * (FU[i] - UAV))
        FV = set_item(FV, i, (1. - params.CU) * (FV[i] - VAV))
    for k in range(NTRA):
        TRAAV = xp.sum(FTRA[:INB+1, k] * (PH[:INB+1] - PH[1:INB+2])) / (PH[0] - PH[INB+1])
        FTRA = set_item(FTRA, (slice(0, INB+1), k), FTRA[:INB+1, k] - TRAAV)

    return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

if HAS_NUMBA:
    _tlift_functional_numba = njit(_tlift_functional)
    _convect_functional_numba = njit(_convect_functional)
else:
    _tlift_functional_numba = _tlift_functional
    _convect_functional_numba = _convect_functional

@jit_compile(backend=np)
def _numpy_vectorized_convect(t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, tra, params, convect_func, tlift_func):
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
        # Explicitly call numba version
        res = convect_func(
            np, tlift_func, t[:, i], q[:, i], qs[:, i], u[:, i], v[:, i],
            p[:, i], ph[:, i], ND, NL, NTRA, DELT, cbmf[i], tra[:, :, i], params
        )
        (ft[:, i], fq[:, i], fu[:, i], fv[:, i], precip[i], 
         wd[i], tprime[i], qprime[i], cbmf_new[i], outcape[i], iflag[i]) = res
            
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
