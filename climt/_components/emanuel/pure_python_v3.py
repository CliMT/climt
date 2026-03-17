# -*- coding: utf-8 -*-
import numpy as np
from typing import NamedTuple
from ..._core.backend import jit_compile, prange
from sympl import (
    ImplicitTendencyComponent,
    get_constant,
    initialize_numpy_arrays_with_properties,
)
from ..._core import ensure_contiguous_state

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x: x

class EmanuelParams(NamedTuple):
    IPBL: int; MINORIG: int; ELCRIT: float; TLCRIT: float; ENTP: float; SIGD: float; SIGS: float; OMTRAIN: float; OMTSNOW: float; COEFFR: float; COEFFS: float; CU: float; BETA: float; DTMAX: float; ALPHA: float; DAMP: float; CPD: float; CPV: float; CL: float; RV: float; RD: float; LV0: float; G: float; ROWL: float; DELT0: float

class EmanuelConvectionPythonV3(ImplicitTendencyComponent):
    input_properties = {
        "air_temperature": {"dims": ["*", "mid_levels"], "units": "degK"},
        "specific_humidity": {"dims": ["*", "mid_levels"], "units": "kg/kg"},
        "eastward_wind": {"dims": ["*", "mid_levels"], "units": "m s^-1"},
        "northward_wind": {"dims": ["*", "mid_levels"], "units": "m s^-1"},
        "air_pressure": {"dims": ["*", "mid_levels"], "units": "mbar"},
        "air_pressure_on_interface_levels": {"dims": ["*", "interface_levels"], "units": "mbar"},
        "cloud_base_mass_flux": {"dims": ["*"], "units": "kg m^-2 s^-1"},
    }

    diagnostic_properties = {
        "convective_state": {"dims": ["*"], "units": "dimensionless", "dtype": np.int32},
        "convective_precipitation_rate": {"dims": ["*"], "units": "mm day^-1"},
        "convective_downdraft_velocity_scale": {"dims": ["*"], "units": "m s^-1"},
        "convective_downdraft_temperature_scale": {"dims": ["*"], "units": "degK"},
        "convective_downdraft_specific_humidity_scale": {"dims": ["*"], "units": "kg/kg"},
        "cloud_base_mass_flux": {"dims": ["*"], "units": "kg m^-2 s^-1"},
        "atmosphere_convective_available_potential_energy": {"dims": ["*"], "units": "J kg^-1"},
        "air_temperature_tendency_from_convection": {"dims": ["*", "mid_levels"], "units": "degK day^-1"},
    }

    tendency_properties = {
        "air_temperature": {"units": "degK s^-1"},
        "specific_humidity": {"units": "kg/kg s^-1"},
        "eastward_wind": {"units": "m s^-2"},
        "northward_wind": {"units": "m s^-2"},
    }

    def __init__(self, **kwargs):
        self.IPBL = 0; self.MINORIG = 1; self.ELCRIT = 0.0011; self.TLCRIT = -55.0; self.ENTP = 1.5; self.SIGD = 0.05; self.SIGS = 0.12; self.OMTRAIN = 50.0; self.OMTSNOW = 5.5; self.COEFFR = 1.0; self.COEFFS = 0.8; self.CU = 0.7; self.BETA = 10.0; self.DTMAX = 0.9; self.ALPHA = 0.1; self.DAMP = 0.1; self.CPD = 1005.7; self.CPV = 1870.0; self.CL = 2500.0; self.RV = 461.5; self.RD = 287.04; self.LV0 = 2.501e6; self.G = 9.8; self.ROWL = 1000.0; self.DELT0 = 300.0
        for key, value in kwargs.items():
            if hasattr(self, key): setattr(self, key, value)
        self._params = EmanuelParams(IPBL=int(self.IPBL), MINORIG=int(self.MINORIG), ELCRIT=float(self.ELCRIT), TLCRIT=float(self.TLCRIT), ENTP=float(self.ENTP), SIGD=float(self.SIGD), SIGS=float(self.SIGS), OMTRAIN=float(self.OMTRAIN), OMTSNOW=float(self.OMTSNOW), COEFFR=float(self.COEFFR), COEFFS=float(self.COEFFS), CU=float(self.CU), BETA=float(self.BETA), DTMAX=float(self.DTMAX), ALPHA=float(self.ALPHA), DAMP=float(self.DAMP), CPD=float(self.CPD), CPV=float(self.CPV), CL=float(self.CL), RV=float(self.RV), RD=float(self.RD), LV0=float(self.LV0), G=float(self.G), ROWL=float(self.ROWL), DELT0=float(self.DELT0))
        super(EmanuelConvectionPythonV3, self).__init__(**kwargs)

    @ensure_contiguous_state
    def array_call(self, state, timestep):
        t = state['air_temperature']; q = state['specific_humidity']; u = state['eastward_wind']; v = state['northward_wind']; p = state['air_pressure']; ph = state['air_pressure_on_interface_levels']
        # Handle sympl (ncol, nlev) -> Emanuel (nlev, ncol)
        # Note: if array_call is called directly with (nlev, ncol) this might fail.
        # But sympl calls always provide (ncol, nlev).
        transposed = False
        if t.shape[0] != ph.shape[0]-1: # Heuristic to detect if transpose is needed
             t = t.T; q = q.T; u = u.T; v = v.T; p = p.T; ph = ph.T
             transposed = True
        
        nlev, ncol = t.shape; from climt._core import bolton_q_sat
        qs = bolton_q_sat(t, p * 100, self.RD, self.RV); cbmf = state.get('cloud_base_mass_flux', np.zeros(ncol)).copy()
        ntra = 0; tra = np.zeros((nlev, 1)); delt = timestep.total_seconds(); tra_vector = np.broadcast_to(tra[:, :, np.newaxis], (nlev, 1, ncol))
        results = _numpy_vectorized_convect(t, q, qs, u, v, p, ph, nlev, nlev-3, ntra, delt, cbmf, tra_vector, self._params)
        
        ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results
        # Transpose back if we transposed in
        if transposed:
            ft = ft.T; fq = fq.T; fu = fu.T; fv = fv.T
        tendencies = {'air_temperature': ft, 'specific_humidity': fq, 'eastward_wind': fu, 'northward_wind': fv}
        diagnostics = {'convective_state': iflag, 'convective_precipitation_rate': precip, 'convective_downdraft_velocity_scale': wd, 'convective_downdraft_temperature_scale': tprime, 'convective_downdraft_specific_humidity_scale': qprime, 'cloud_base_mass_flux': cbmf_new, 'atmosphere_convective_available_potential_energy': outcape}
        diagnostics['air_temperature_tendency_from_convection'] = ft * 86400.0
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
    return FT, FQ, np.zeros(ND), np.zeros(ND), 0.0, 0.0, 0.0, 0.0, CBMF, OUTCAPE, IFLAG


@jit_compile(backend=np, parallel=True)
def _numpy_vectorized_convect(t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, tra, params):
    nlev, ncol = t.shape; ft = np.zeros(t.shape); fq = np.zeros(q.shape); fu = np.zeros(u.shape); fv = np.zeros(v.shape); precip = np.zeros(ncol); wd = np.zeros(ncol); tprime = np.zeros(ncol); qprime = np.zeros(ncol); cbmf_new = np.zeros(ncol); outcape = np.zeros(ncol); iflag = np.zeros(ncol, dtype=np.int32)
    for i in prange(ncol):
        res = _convect_functional_np(t[:, i], q[:, i], qs[:, i], u[:, i], v[:, i], p[:, i], ph[:, i], ND, NL, NTRA, DELT, cbmf[i], tra[:, :, i], params)
        (ft[:, i], fq[:, i], fu[:, i], fv[:, i], precip[i], wd[i], tprime[i], qprime[i], cbmf_new[i], outcape[i], iflag[i]) = res
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
