# -*- coding: utf-8 -*-
import numpy as np
from typing import NamedTuple
from ..._core.backend import get_array_namespace, set_item, jit_compile, vectorize_component, prange

try:
    from numba import njit, extending
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x: x

class EmanuelParams(NamedTuple):
    IPBL: int; MINORIG: int; ELCRIT: float; TLCRIT: float; ENTP: float; SIGD: float; SIGS: float; OMTRAIN: float; OMTSNOW: float; COEFFR: float; COEFFS: float; CU: float; BETA: float; DTMAX: float; ALPHA: float; DAMP: float; CPD: float; CPV: float; CL: float; RV: float; RD: float; LV0: float; G: float; ROWL: float; DELT0: float

class EmanuelConvectionPythonV3(object):
    def __init__(self, **kwargs):
        self.IPBL = 0; self.MINORIG = 1; self.ELCRIT = 0.0011; self.TLCRIT = -55.0; self.ENTP = 1.5; self.SIGD = 0.05; self.SIGS = 0.12; self.OMTRAIN = 50.0; self.OMTSNOW = 5.5; self.COEFFR = 1.0; self.COEFFS = 0.8; self.CU = 0.7; self.BETA = 10.0; self.DTMAX = 0.9; self.ALPHA = 0.1; self.DAMP = 0.1; self.CPD = 1005.7; self.CPV = 1870.0; self.CL = 2500.0; self.RV = 461.5; self.RD = 287.04; self.LV0 = 2.501e6; self.G = 9.8; self.ROWL = 1000.0; self.DELT0 = 300.0
        for key, value in kwargs.items():
            if hasattr(self, key): setattr(self, key, value)
        self._params = EmanuelParams(IPBL=int(self.IPBL), MINORIG=int(self.MINORIG), ELCRIT=float(self.ELCRIT), TLCRIT=float(self.TLCRIT), ENTP=float(self.ENTP), SIGD=float(self.SIGD), SIGS=float(self.SIGS), OMTRAIN=float(self.OMTRAIN), OMTSNOW=float(self.OMTSNOW), COEFFR=float(self.COEFFR), COEFFS=float(self.COEFFS), CU=float(self.CU), BETA=float(self.BETA), DTMAX=float(self.DTMAX), ALPHA=float(self.ALPHA), DAMP=float(self.DAMP), CPD=float(self.CPD), CPV=float(self.CPV), CL=float(self.CL), RV=float(self.RV), RD=float(self.RD), LV0=float(self.LV0), G=float(self.G), ROWL=float(self.ROWL), DELT0=float(self.DELT0))

    def array_call(self, state, timestep):
        t = state['air_temperature']; q = state['specific_humidity']; u = state['eastward_wind']; v = state['northward_wind']; p = state['air_pressure']; ph = state['air_pressure_on_interface_levels']
        nlev, ncol = t.shape; xp = get_array_namespace(t); from climt._core import bolton_q_sat
        qs = bolton_q_sat(t, p * 100, self.RD, self.RV); cbmf = state.get('cloud_base_mass_flux', xp.zeros(ncol)).copy()
        ntra = 0; tra = xp.zeros((nlev, 1)); delt = timestep.total_seconds(); tra_vector = xp.broadcast_to(tra[:, :, xp.newaxis], (nlev, 1, ncol))
        if xp is np:
            results = _numpy_vectorized_convect(t, q, qs, u, v, p, ph, nlev, nlev-3, ntra, delt, cbmf, tra_vector, self._params)
        else:
            import jax
            @jax.jit
            def vectorized_jax_call(t, q, qs, u, v, p, ph, cbmf, tra_vector):
                def jax_column_wrapper(T, Q, QS, U, V, P, PH, CBMF, TRA):
                    return _convect_functional_jax(T, Q, QS, U, V, P, PH, nlev, nlev-3, ntra, delt, CBMF, TRA, self._params)
                return jax.vmap(jax_column_wrapper, in_axes=(1, 1, 1, 1, 1, 1, 1, 0, 2))(t, q, qs, u, v, p, ph, cbmf, tra_vector)
            results = vectorized_jax_call(t, q, qs, u, v, p, ph, cbmf, tra_vector)
        ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results
        tendencies = {'air_temperature': ft, 'specific_humidity': fq, 'eastward_wind': fu, 'northward_wind': fv}
        diagnostics = {'convective_state': iflag, 'convective_precipitation_rate': precip, 'convective_downdraft_velocity_scale': wd, 'convective_downdraft_temperature_scale': tprime, 'convective_downdraft_specific_humidity_scale': qprime, 'cloud_base_mass_flux': cbmf_new, 'atmosphere_convective_available_potential_energy': outcape}
        return tendencies, diagnostics

@njit
def _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params):
    CPVMCL = params.CL - params.CPV; EPS = params.RD / params.RV; EPSI = 1.0 / EPS
    AH0 = (params.CPD * (1. - Q[NK]) + params.CL * Q[NK]) * T[NK] + Q[NK] * (params.LV0 - CPVMCL * (T[NK] - 273.15)) + GZ[NK]
    CPP = params.CPD * (1. - Q[NK]) + Q[NK] * params.CPV; CPINV = 1.0 / CPP
    if KK == 1:
        for i in range(ICB): CLW[i] = 0.0
        for i in range(NK, ICB):
            TPK[i] = T[NK] - (GZ[i] - GZ[NK]) * CPINV
            TVP[i] = TPK[i] * (1. + Q[NK] * EPSI)
    NST = NL if KK == 2 else ICB; NSB = ICB + 1 if KK == 2 else ICB
    for i in range(NSB, NST + 1):
        TG = T[i]; QG = QS[i]; ALV = params.LV0 - CPVMCL * (T[i] - 273.15)
        for j in range(2):
            S = 1.0 / (params.CPD + ALV * ALV * QG / (params.RV * T[i] * T[i]))
            AHG = params.CPD * TG + (params.CL - params.CPD) * Q[NK] * T[i] + ALV * QG + GZ[i]
            TG = max(TG + S * (AH0 - AHG), 35.0); TC = TG - 273.15; DENOM = 243.5 + TC
            ES = 6.112 * np.exp(17.67 * TC / DENOM) if TC >= 0.0 else np.exp(23.33086 - 6111.72784 / TG + 0.15215 * np.log(TG))
            QG = EPS * ES / (P[i] - ES * (1. - EPS))
        TPK[i] = (AH0 - (params.CL - params.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG) / params.CPD
        CLW[i] = max(0.0, Q[NK] - QG); TVP[i] = TPK[i] * (1. + (QG / (1. - Q[NK])) * EPSI)
    return TVP, TPK, CLW

def _tlift_functional_jax(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params):
    import jax.numpy as jnp
    CPVMCL = params.CL - params.CPV; EPS = params.RD / params.RV; EPSI = 1.0 / EPS
    AH0 = (params.CPD * (1. - Q[NK]) + params.CL * Q[NK]) * T[NK] + Q[NK] * (params.LV0 - CPVMCL * (T[NK] - 273.15)) + GZ[NK]
    CPP = params.CPD * (1. - Q[NK]) + Q[NK] * params.CPV; CPINV = 1.0 / CPP
    lvl = jnp.arange(ND); mask_kk1 = (KK == 1)
    CLW = jnp.where((lvl < ICB) & mask_kk1, 0.0, CLW)
    TPK_new = T[NK] - (GZ[:ND] - GZ[NK]) * CPINV; TVP_new = TPK_new * (1. + Q[NK] * EPSI)
    mask_in_range = (lvl >= NK) & (lvl < ICB) & mask_kk1
    TPK = jnp.where(mask_in_range, TPK_new, TPK); TVP = jnp.where(mask_in_range, TVP_new, TVP)
    NST = jnp.where(KK == 2, NL, ICB); NSB = jnp.where(KK == 2, ICB + 1, ICB)
    ALV_all = params.LV0 - CPVMCL * (T - 273.15); TG_it = T; QG_it = QS
    for _ in range(2):
        S = 1.0 / (params.CPD + ALV_all * ALV_all * QG_it / (params.RV * T * T))
        AHG = params.CPD * TG_it + (params.CL - params.CPD) * Q[NK] * T + ALV_all * QG_it + GZ[:ND]
        TG_it = jnp.maximum(TG_it + S * (AH0 - AHG), 35.0); TC = TG_it - 273.15; DENOM = 243.5 + TC
        ES = jnp.where(TC >= 0.0, 6.112 * jnp.exp(17.67 * TC / DENOM), jnp.exp(23.33086 - 6111.72784 / TG_it + 0.15215 * jnp.log(TG_it)))
        QG_it = EPS * ES / (P - ES * (1. - EPS))
    tpk_vec = (AH0 - (params.CL - params.CPD) * Q[NK] * T - GZ[:ND] - ALV_all * QG_it) / params.CPD
    clw_vec = jnp.maximum(0.0, Q[NK] - QG_it); tvp_vec = tpk_vec * (1. + (QG_it / (1. - Q[NK])) * EPSI)
    is_active = (lvl >= NSB) & (lvl <= NST)
    TPK = jnp.where(is_active, tpk_vec, TPK); CLW = jnp.where(is_active, clw_vec, CLW); TVP = jnp.where(is_active, tvp_vec, TVP)
    return TVP, TPK, CLW

@njit
def _convect_functional_np(T_in, Q_in, QS_in, U_in, V_in, P_in, PH_in, ND, NL, NTRA, DELT, CBMF_in, TRA_in, params):
    T = T_in.copy(); Q = Q_in.copy(); QS = QS_in.copy(); U = U_in.copy(); V = V_in.copy(); P = P_in.copy(); PH = PH_in.copy(); TRA = TRA_in.copy(); CBMF = CBMF_in
    FT = np.zeros(ND); FQ = np.zeros(ND); FU = np.zeros(ND); FV = np.zeros(ND); FTRA = np.zeros((ND, max(1, NTRA)))
    CPVMCL = params.CL - params.CPV; EPS = params.RD / params.RV; EPSI = 1.0 / EPS; GINV = 1.0 / params.G; DELTI = 1.0 / DELT
    TH = np.zeros(NL + 1)
    for i in range(NL + 1): RDCP = (params.RD * (1. - Q[i]) + Q[i] * params.RV) / (params.CPD * (1. - Q[i]) + Q[i] * params.CPV); TH[i] = T[i] * (1000.0 / P[i])**RDCP
    GZ = np.zeros(ND + 1); CPN = np.zeros(ND + 1); H = np.zeros(ND + 1); LV = np.zeros(ND + 1); HM = np.zeros(ND + 1); TV = np.zeros(ND + 1)
    GZ[0] = 0.0; CPN[0] = params.CPD * (1. - Q[0]) + Q[0] * params.CPV; H[0] = T[0] * CPN[0]; LV[0] = params.LV0 - CPVMCL * (T[0] - 273.15); HM[0] = LV[0] * Q[0]; TV[0] = T[0] * (1. + Q[0] * EPSI - Q[0])
    AHMIN = 1.0E12; IHMIN = 0
    for i in range(1, NL + 1):
        TVX = T[i] * (1. + Q[i] * EPSI - Q[i]); TVY = T[i-1] * (1. + Q[i-1] * EPSI - Q[i-1]); GZ[i] = GZ[i-1] + 0.5 * params.RD * (TVX + TVY) * (P[i-1] - P[i]) / PH[i]
        CPN[i] = params.CPD * (1. - Q[i]) + params.CPV * Q[i]; H[i] = T[i] * CPN[i] + GZ[i]; LV[i] = params.LV0 - CPVMCL * (T[i] - 273.15)
        HM[i] = (params.CPD * (1. - Q[i]) + params.CL * Q[i]) * (T[i] - T[0]) + LV[i] * Q[i] + GZ[i]; TV[i] = T[i] * (1. + Q[i] * EPSI - Q[i])
        if (i+1) >= params.MINORIG and HM[i] < AHMIN and HM[i] < HM[i-1]: AHMIN = HM[i]; IHMIN = i
    IHMIN = min(IHMIN, NL - 1); AHMAX = -1.0E12; NK = 0
    for i in range(params.MINORIG - 1, IHMIN + 1):
        if HM[i] > AHMAX: NK = i; AHMAX = HM[i]
    if T[NK] < 250.0 or Q[NK] <= 0.0 or IHMIN == (NL - 1): return FT, FQ, FU, FV, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0
    RH = Q[NK] / QS[NK]; CHI = T[NK] / (1669.0 - 122.0 * RH - T[NK]); PLCL = P[NK] * (RH**CHI)
    if PLCL < 200.0 or PLCL >= 2000.0: return FT, FQ, FU, FV, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2
    ICB = NL - 1
    for i in range(NK + 1, NL + 1):
        if P[i] < PLCL: ICB = i; break
    if ICB >= (NL - 1): return FT, FQ, FU, FV, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3
    TVP = np.zeros(ND); TP = np.zeros(ND); CLW = np.zeros(ND); TVP, TP, CLW = _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params)
    for i in range(NK, ICB + 1): TVP[i] -= TP[i] * Q[NK]
    if CBMF == 0.0 and TVP[ICB] <= (TV[ICB] - params.DTMAX): return FT, FQ, FU, FV, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0
    IFLAG = 1; TVP, TP, CLW = _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params)
    EP = np.zeros(ND); SIGP = np.zeros(ND)
    for i in range(NK + 1, NL + 1):
        TCA = TP[i] - 273.15; ELACRIT = params.ELCRIT if TCA >= 0.0 else params.ELCRIT * (1.0 - TCA / params.TLCRIT); EP[i] = 0.999 * (1.0 - max(ELACRIT, 0.0) / max(CLW[i], 1.0E-8)); EP[i] = max(min(EP[i], 0.999), 0.0); SIGP[i] = params.SIGS
    for i in range(ICB + 1, NL + 1): TVP[i] -= TP[i] * Q[NK]
    HP = H.copy(); LVCP = LV / CPN; QP = np.zeros(ND + 1); UP = np.zeros(ND + 1); VP = np.zeros(ND + 1); QP[0] = Q[0]; UP[0] = U[0]; VP[0] = V[0]
    for i in range(1, NL + 1): QP[i] = Q[i-1]; UP[i] = U[i-1]; VP[i] = V[i-1]
    CAPE = 0.0; CAPEM = 0.0; INB = ICB + 1; INB1 = INB; BYP = 0.0
    for i in range(ICB + 1, NL):
        BY = (TVP[i] - TV[i]) * (PH[i] - PH[i+1]) / P[i]; CAPE += BY
        if BY >= 0.0: INB1 = i + 1
        if CAPE > 0.0: INB = i + 1; BYP = (TVP[i+1] - TV[i+1]) * (PH[i+1] - PH[i+2]) / P[i+1]; CAPEM = CAPE
    INB = max(INB, INB1); DEFRAC = max(CAPEM - CAPE, 0.001); FRAC = min(max(-CAPE / DEFRAC, 0.0), 1.0); OUTCAPE = CAPEM + BYP
    for i in range(ICB, INB + 1): HP[i] = H[NK] + (LV[i] + (params.CPD - params.CPV) * T[i]) * EP[i] * CLW[i]
    TVPPLCL = TVP[ICB-1] - params.RD * TVP[ICB-1] * (P[ICB-1] - PLCL) / (CPN[ICB-1] * P[ICB-1]); TVAPLCL = TV[ICB] + (TVP[ICB] - TVP[ICB+1]) * (PLCL - P[ICB]) / (P[ICB] - P[ICB+1])
    DTPBL = 0.0
    for i in range(NK, ICB): DTPBL += (TVP[i] - TV[i]) * (PH[i] - PH[i+1])
    DTPBL /= max(PH[NK] - PH[ICB], 1.0); DTMA = TVPPLCL - TVAPLCL + params.DTMAX + DTPBL
    CBMF = max((1. - params.DAMP * DELT / params.DELT0) * CBMF + 0.1 * params.ALPHA * DTMA, 0.0)
    if CBMF == 0.0: return FT, FQ, FU, FV, 0.0, 0.0, 0.0, 0.0, 0.0, OUTCAPE, IFLAG
    # Mixing matrices
    QENT = np.zeros((ND + 1, ND + 1)); ELIJ = np.zeros((ND + 1, ND + 1)); MENT = np.zeros((ND + 1, ND + 1)); SIJ = np.zeros((ND + 1, ND + 1)); UENT = np.zeros((ND + 1, ND + 1)); VENT = np.zeros((ND + 1, ND + 1)); NENT = np.zeros(ND + 1, dtype=np.int32); M = np.zeros(ND + 1)
    DBOSUM = 0.0
    for i in range(ICB + 1, INB + 1):
        k = min(i, INB1); DBO = abs(TV[k] - TVP[k]) + params.ENTP * 0.02 * (PH[k] - PH[k+1]); DBOSUM += DBO; M[i] = CBMF * DBO
    if DBOSUM > 0:
        for i in range(ICB + 1, INB + 1): M[i] /= DBOSUM
    for i in range(ICB + 1, INB + 1):
        QTI = Q[NK] - EP[i] * CLW[i]
        for j in range(ICB, INB + 1):
            BF2 = 1. + LV[j] * LV[j] * QS[j] / (params.RV * T[j] * T[j] * params.CPD); ANUM = H[j] - HP[i] + (params.CPV - params.CPD) * T[j] * (QTI - Q[j]); DENOM = H[i] - HP[i] + (params.CPD - params.CPV) * (Q[i] - QTI) * T[j]; DEI = max(abs(DENOM), 0.01); SIJ[i, j] = ANUM / DEI; SIJ[i, i] = 1.0; ALTEM = (SIJ[i, j] * Q[i] + (1. - SIJ[i, j]) * QTI - QS[j]) / BF2; CWAT = CLW[j] * (1. - EP[j])
            if (SIJ[i, j] < 0.0 or SIJ[i, j] > 1.0 or ALTEM > CWAT) and j > i:
                ANUM -= LV[j] * (QTI - QS[j] - CWAT * BF2); DENOM += LV[j] * (Q[i] - QTI); DEI = max(abs(DENOM), 0.01); SIJ[i, j] = ANUM / DEI; ALTEM = SIJ[i, j] * Q[i] + (1. - SIJ[i, j]) * QTI - QS[j] - (BF2 - 1.) * CWAT
            if 0.0 < SIJ[i, j] < 0.9:
                QENT[i, j] = SIJ[i, j] * Q[i] + (1. - SIJ[i, j]) * QTI; UENT[i, j] = SIJ[i, j] * U[i] + (1. - SIJ[i, j]) * U[NK]; VENT[i, j] = SIJ[i, j] * V[i] + (1. - SIJ[i, j]) * V[NK]; ELIJ[i, j] = max(0.0, ALTEM); MENT[i, j] = M[i] / (1. - SIJ[i, j]); NENT[i] += 1
            SIJ[i, j] = min(max(SIJ[i, j], 0.0), 1.0)
        if NENT[i] == 0: MENT[i, i] = M[i]; QENT[i, i] = Q[NK] - EP[i] * CLW[i]; UENT[i, i] = U[NK]; VENT[i, i] = V[NK]; ELIJ[i, i] = CLW[i]; SIJ[i, i] = 1.0
    for i in range(ICB + 1, INB + 1):
        if NENT[i] != 0:
            QP1 = Q[NK] - EP[i] * CLW[i]; ANUM = H[i] - HP[i] - LV[i] * (QP1 - QS[i]); DENOM = H[i] - HP[i] + LV[i] * (Q[i] - QP1); DEI = max(abs(DENOM), 0.01); SCRIT = ANUM / DEI; ALT = QP1 - QS[i] + SCRIT * (Q[i] - QP1); SCRIT = max(SCRIT, 0.0); ASIJ = 0.0
            for j in range(ICB, INB + 1):
                if 0.0 < SIJ[i, j] < 0.9:
                    if j > i: SMID = min(SIJ[i, j], SCRIT); SJMAX = SMID; SJMIN = SMID
                    else: SJMAX = max(SIJ[i, j+1], SCRIT); SMID = max(SIJ[i, j], SCRIT); SJMIN = max(SIJ[i, j-1] if j > 0 else 0.0, SCRIT)
                    ASIJ += (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (PH[j] - PH[j+1]); MENT[i, j] *= (abs(SJMAX - SMID) + abs(SJMIN - SMID)) * (PH[j] - PH[j+1])
            ASIJ = max(ASIJ, 1.0E-21)
            for j in range(ICB, INB + 1): MENT[i, j] /= ASIJ
            if MENT[i, ICB:INB+1].sum() < 1.0E-18: NENT[i] = 0; MENT[i, i] = M[i]; QENT[i, i] = Q[NK] - EP[i] * CLW[i]; UENT[i, i] = U[NK]; VENT[i, i] = V[NK]; ELIJ[i, i] = CLW[i]; SIJ[i, i] = 1.0
    WATER = np.zeros(ND + 1); EVAP = np.zeros(ND + 1); WT = np.full(ND + 1, params.OMTSNOW); MP = np.zeros(ND + 1); PRECIP = 0.0
    if EP[INB] >= 0.0001:
        for i in range(INB, -1, -1):
            WDTRAIN = params.G * EP[i] * M[i] * CLW[i]
            if i > 0:
                for j in range(i): WDTRAIN += params.G * max(0.0, ELIJ[j, i] - (1. - EP[i]) * CLW[i]) * MENT[j, i]
            AFAC = max(params.COEFFR * PH[i] * (QS[i] - 0.5 * (Q[i] + (QP[i+1] if i < ND else Q[i]))) / (1.0E4 + 2.0E3 * PH[i] * QS[i]), 0.0) if i < ND else 0.0
            REVAP = 0.5 * (-(100.*(PH[i]-PH[i+1])*params.SIGS*AFAC/WT[i]) + np.sqrt((100.*(PH[i]-PH[i+1])*params.SIGS*AFAC/WT[i])**2 + 4.*((WATER[i+1] if i < ND else 0.0)*WT[i+1] + WDTRAIN / params.SIGD)/WT[i])) if i < ND else 0.0
            EVAP[i] = params.SIGS * AFAC * REVAP; WATER[i] = REVAP * REVAP
            if i > 0:
                DHDP = max((H[i] - H[i-1]) / (P[i-1] - P[i]), 10.0); MP[i] = 100. * GINV * LV[i] * params.SIGD * EVAP[i] / DHDP
                FAC = 20.0 / (PH[i-1] - PH[i]); MP[i] = ((FAC * MP[i+1] if i < ND else 0.0) + MP[i]) / (1. + FAC)
            if i != INB:
                QSTM = QS[max(0, i-1)]; RAT = MP[i+1] / MP[i] if (MP[i] > MP[i+1] and MP[i] > 0) else 1.0
                if MP[i] > MP[i+1]: QP[i] = QP[i+1] * RAT + Q[i] * (1.0 - RAT) + 100. * GINV * params.SIGD * (PH[i] - PH[i+1]) * (EVAP[i] / MP[i])
                elif MP[i+1] > 0.0: QP[i] = ((GZ[i+1] if i < ND else GZ[i]) - GZ[i] + (QP[i+1] if i < ND else QP[i]) * ((LV[i+1] if i < ND else LV[i]) + (T[i+1] if i < ND else T[i]) * (params.CL - params.CPD)) + params.CPD * ((T[i+1] if i < ND else T[i]) - T[i])) / (LV[i] + T[i] * (params.CL - params.CPD))
                QP[i] = max(min(QP[i], QSTM), 0.0)
        PRECIP = WT[0] * params.SIGD * WATER[0] * 3600. * 24000. / (params.ROWL * params.G)
    WD = params.BETA * abs(MP[ICB]) * 0.01 * params.RD * T[ICB] / (params.SIGD * P[ICB]); QPRIME = 0.5 * (QP[0] - Q[0]); TPRIME = params.LV0 * QPRIME / params.CPD; DPINV = 0.01 / (PH[0] - PH[1]); AM = M[1:INB+1].sum() if NK == 0 else 0.0
    FT[0] = params.G * DPINV * AM * (T[1] - T[0] + (GZ[1] - GZ[0]) / CPN[0]) - LVCP[0] * params.SIGD * EVAP[0] + params.SIGD * WT[1] * (params.CL - params.CPD) * WATER[1] * (T[1] - T[0]) * DPINV / CPN[0]
    FQ[0] = params.G * MP[1] * (QP[1] - Q[0]) * DPINV + params.SIGD * EVAP[0] + params.G * AM * (Q[1] - Q[0]) * DPINV
    for i in range(1, INB + 1):
        DPINV = 0.01 / (PH[i] - PH[i+1]); CPINV = 1.0 / CPN[i]; AMP1 = M[i+1:INB+2].sum() if i >= NK else 0.0
        for k in range(i + 1): AMP1 += MENT[k, i+1:INB+2].sum()
        AD = MENT[i:INB+1, :i].sum(); FT[i] = params.G * DPINV * (AMP1 * (T[i+1] - T[i] + (GZ[i+1] - GZ[i]) * CPINV) - AD * (T[i] - T[i-1] + (GZ[i] - GZ[i-1]) * CPINV)) - params.SIGD * LVCP[i] * EVAP[i]; FT[i] += params.G * DPINV * MENT[i, i] * (HP[i] - H[i] + T[i] * (params.CPV - params.CPD) * (Q[i] - QENT[i, i])) * CPINV; FT[i] += params.SIGD * WT[i+1] * (params.CL - params.CPD) * WATER[i+1] * (T[i+1] - T[i]) * DPINV * CPINV
        FQ[i] = params.G * DPINV * (AMP1 * (Q[i+1] - Q[i]) - AD * (Q[i] - Q[i-1]))
        for k in range(i): FQ[i] += params.G * DPINV * MENT[k, i] * (QENT[k, i] - max(0.0, ELIJ[k, i] - (1. - EP[i]) * CLW[i]) - Q[i])
        for k in range(i, INB + 1): FQ[i] += params.G * DPINV * MENT[k, i] * (QENT[k, i] - Q[i])
        FQ[i] += params.SIGD * EVAP[i] + params.G * (MP[i+1] * (QP[i+1] - Q[i]) - MP[i] * (QP[i] - Q[i-1])) * DPINV
    FQOLD = FQ[INB]; FQ[INB] *= (1. - FRAC); FQ[INB-1] += FRAC * FQOLD * ((PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])) * LV[INB] / LV[INB-1]; FTOLD = FT[INB]; FT[INB] *= (1. - FRAC); FT[INB-1] += FRAC * FTOLD * ((PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])) * CPN[INB] / CPN[INB-1]
    ENTS = (CPN[:INB+1] * FT[:INB+1] + LV[:INB+1] * FQ[:INB+1]).dot(PH[:INB+1] - PH[1:INB+2]) / (PH[0] - PH[INB+1])
    for i in range(INB + 1): FT[i] -= ENTS / CPN[i]
    return FT, FQ, np.zeros(ND), np.zeros(ND), PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG

def _convect_functional_jax(T_in, Q_in, QS_in, U_in, V_in, P_in, PH_in, ND, NL, NTRA, DELT, CBMF_in, TRA_in, params):
    import jax; import jax.numpy as jnp
    T = T_in; Q = Q_in; QS = QS_in; U = U_in; V = V_in; P = P_in; PH = PH_in; TRA = TRA_in; CBMF = CBMF_in
    FT = jnp.zeros(ND); FQ = jnp.zeros(ND); FTRA = jnp.zeros((ND, max(1, NTRA))); CPVMCL = params.CL - params.CPV; EPS = params.RD / params.RV; EPSI = 1.0 / EPS; GINV = 1.0 / params.G; DELTI = 1.0 / DELT; RDCP = (params.RD * (1. - Q) + Q * params.RV) / (params.CPD * (1. - Q) + Q * params.CPV); TH = T * (1000.0 / P)**RDCP; GZ_acc = jnp.zeros(ND + 1); TV_all = T * (1. + Q * EPSI - Q)
    def gz_step(carry, i): gz_curr = carry + 0.5 * params.RD * (TV_all[i] + TV_all[i-1]) * (P[i-1] - P[i]) / PH[i]; return gz_curr, gz_curr
    _, gz_scan = jax.lax.scan(gz_step, 0.0, jnp.arange(1, NL + 1)); GZ = jnp.concatenate([jnp.zeros(1), gz_scan, jnp.zeros(ND-NL)]); CPN = params.CPD * (1. - Q) + params.CPV * Q; H = T * CPN + GZ[:ND]; LV = params.LV0 - CPVMCL * (T - 273.15); HM = (params.CPD * (1. - Q) + params.CL * Q) * (T - T[0]) + LV * Q + GZ[:ND]; TV = T * (1. + Q * EPSI - Q); AHMIN = 1.0E12; IHMIN = jnp.array(0, dtype=jnp.int32); lvl = jnp.arange(ND); lvl_h = jnp.arange(ND + 1)
    for i in range(1, NL + 1): cond = ((i+1) >= params.MINORIG) & (HM[i] < AHMIN) & (HM[i] < HM[i-1]); AHMIN = jnp.where(cond, HM[i], AHMIN); IHMIN = jnp.where(cond, i, IHMIN)
    nk_mask = (lvl >= params.MINORIG - 1) & (lvl <= jnp.minimum(IHMIN, NL - 1)); NK = jnp.argmax(jnp.where(nk_mask, HM, -1.0E15)); active = (T[NK] >= 250.0) & (Q[NK] > 0.0) & (IHMIN != (NL - 1)); RH = Q[NK] / QS[NK]; CHI = T[NK] / (1669.0 - 122.0 * RH - T[NK]); PLCL = P[NK] * (RH**CHI); active = active & (PLCL >= 200.0) & (PLCL < 2000.0); icb_mask = (lvl >= NK + 1) & (lvl <= NL) & (P < PLCL); ICB = jnp.where(jnp.any(icb_mask), jnp.argmax(icb_mask), NL - 1); active = active & (ICB < (NL - 1)); TVP = jnp.zeros(ND); TP = jnp.zeros(ND); CLW = jnp.zeros(ND); TVP, TP, CLW = _tlift_functional_jax(P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params)
    TVP = jnp.where((lvl >= NK) & (lvl <= ICB), TVP - TP * Q[NK], TVP); tvp_icb = jnp.sum(jnp.where(lvl == ICB, TVP, 0.0)); tv_icb = jnp.sum(jnp.where(lvl == ICB, TV, 0.0)); active = active & ~((CBMF == 0.0) & (tvp_icb <= (tv_icb - params.DTMAX))); IFLAG = jnp.where(active, jnp.array(1, dtype=jnp.int32), jnp.array(0, dtype=jnp.int32)); TVP, TP, CLW = _tlift_functional_jax(P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params)
    TCA = TP - 273.15; ELACRIT = jnp.where(TCA >= 0.0, params.ELCRIT, params.ELCRIT * (1.0 - TCA / params.TLCRIT)); ep_val = 0.999 * (1.0 - jnp.maximum(ELACRIT, 0.0) / jnp.maximum(CLW, 1.0E-8)); EP = jnp.where((lvl > NK) & (lvl <= NL), jnp.maximum(jnp.minimum(ep_val, 0.999), 0.0), 0.0); TVP = jnp.where((lvl > ICB) & (lvl <= NL), TVP - TP * Q[NK], TVP); HP = H.copy(); QP = jnp.zeros(ND + 1); QP = QP.at[0].set(Q[0])
    for i in range(1, NL + 1): QP = QP.at[i].set(Q[i-1])
    def cape_step(carry, i): cape_acc, inb_acc, inb1_acc, capem_acc, byp_acc = carry; is_range = (i >= ICB + 1) & (i < NL); BY = (TVP[i] - TV[i]) * (PH[i] - PH[i+1]) / P[i]; cape_next = cape_acc + jnp.where(is_range, BY, 0.0); inb1_next = jnp.where(is_range & (BY >= 0.0), i + 1, inb1_acc); cond_inb = is_range & (cape_next > 0.0); inb_next = jnp.where(cond_inb, i + 1, inb_acc); capem_next = jnp.where(cond_inb, cape_next, capem_acc); byp_next = jnp.where(cond_inb, (jnp.where(lvl == i+1, TVP, 0.0).sum() - jnp.where(lvl == i+1, TV, 0.0).sum()) * (PH[i+1] - PH[i+2]) / P[i+1], byp_acc); return (cape_next, inb_next, inb1_next, capem_next, byp_next), None
    (CAPE_fin, INB, INB1, CAPEM, BYP), _ = jax.lax.scan(cape_step, (0.0, ICB + 1, ICB + 1, 0.0, 0.0), jnp.arange(ND)); INB = jnp.maximum(INB, INB1); DEFRAC = jnp.maximum(CAPEM - (CAPEM + BYP), 0.001); FRAC = jnp.minimum(jnp.maximum(-(CAPEM + BYP) / DEFRAC, 0.0), 1.0); OUTCAPE = CAPEM + BYP; hp_val = H[NK] + (LV + (params.CPD - params.CPV) * T) * EP * CLW; HP = HP.at[:ND].set(jnp.where((lvl >= ICB) & (lvl <= INB), hp_val, HP[:ND])); TVPPLCL = jnp.sum(jnp.where(lvl == ICB-1, TVP - params.RD * TVP * (P - PLCL) / (CPN * P), 0.0)); TVAPLCL = tv_icb + jnp.sum(jnp.where(lvl == ICB, (TVP - jnp.roll(TVP, -1)) * (PLCL - P) / (P - jnp.roll(P, -1)), 0.0)); DTPBL = jnp.sum(jnp.where((lvl >= NK) & (lvl < ICB), (TVP - TV) * (PH[:ND] - PH[1:]), 0.0)); pbl_depth = jnp.maximum(jnp.sum(jnp.where(lvl_h == NK, PH, 0.0)) - jnp.sum(jnp.where(lvl_h == ICB, PH, 0.0)), 1.0); DTMA = TVPPLCL - TVAPLCL + params.DTMAX + DTPBL / pbl_depth; CBMF = jnp.maximum((1. - params.DAMP * DELT / params.DELT0) * CBMF + 0.1 * params.ALPHA * DTMA, 0.0); mask_final = active & (CBMF > 0.0)
    return (jnp.where(mask_final, FT, 0.0), jnp.where(mask_final, FQ, 0.0), jnp.zeros(ND), jnp.zeros(ND), jnp.array(0.0), jnp.array(0.0), jnp.array(0.0), jnp.array(0.0), jnp.where(active, CBMF, 0.0), jnp.where(mask_final, OUTCAPE, 0.0), IFLAG)

@jit_compile(backend=np, parallel=True)
def _numpy_vectorized_convect(t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, tra, params):
    nlev, ncol = t.shape; ft = np.zeros(t.shape); fq = np.zeros(q.shape); fu = np.zeros(u.shape); fv = np.zeros(v.shape); precip = np.zeros(ncol); wd = np.zeros(ncol); tprime = np.zeros(ncol); qprime = np.zeros(ncol); cbmf_new = np.zeros(ncol); outcape = np.zeros(ncol); iflag = np.zeros(ncol, dtype=np.int32)
    for i in prange(ncol):
        res = _convect_functional_np(t[:, i], q[:, i], qs[:, i], u[:, i], v[:, i], p[:, i], ph[:, i], ND, NL, NTRA, DELT, cbmf[i], tra[:, :, i], params)
        (ft[:, i], fq[:, i], fu[:, i], fv[:, i], precip[i], wd[i], tprime[i], qprime[i], cbmf_new[i], outcape[i], iflag[i]) = res
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
