# -*- coding: utf-8 -*-
"""Differentiable JAX kernel for the Emanuel (1991) convection scheme.

Imported lazily from ``EmanuelConvectionPython.array_call`` when JAX arrays
are supplied so ``import climt`` never imports ``jax``.
"""
import functools

import jax
import jax.numpy as jnp


def _sat_vap_pressure_jax(T, species_id):
    """JAX mirror of condensibles._sat_vap_pressure (hPa). species_id is static."""
    if species_id == 0:  # H2O
        TC = T - 273.15
        # Safe operands in both branches so reverse-mode AD never evaluates
        # 1/(243.5 + TC) <= 0 or log(T <= 0) on the unselected side of jnp.where.
        TC_warm = jnp.where(TC >= 0.0, TC, 0.0)
        T_cold = jnp.where(TC < 0.0, jnp.maximum(T, 1.0), 273.15)
        return jnp.where(
            TC >= 0.0,
            6.112 * jnp.exp(17.67 * TC_warm / (243.5 + TC_warm)),
            jnp.exp(23.33086 - 6111.72784 / T_cold + 0.15215 * jnp.log(T_cold)),
        )
    elif species_id == 1:  # CH4
        T_safe = jnp.maximum(T, 1.0)
        return 1013.25 * jnp.exp((5.1e5 / 518.3) * (1.0 / 111.65 - 1.0 / T_safe))
    else:  # CO2
        T_safe = jnp.maximum(T, 1.0)
        return 1013.25 * jnp.exp((5.71e5 / 188.9) * (1.0 / 194.7 - 1.0 / T_safe))


def compute_qs_jax(T, P, cond, RD):
    """JAX mirror of condensibles.compute_qs (P in hPa). Vectorized over the array T."""
    EPS = RD / cond.RV
    es = _sat_vap_pressure_jax(T, cond.species_id)
    return EPS * es / (P - (1.0 - EPS) * es)


def _tlift_functional_jax(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params, cond):
    # KK is always a Python literal (1 or 2) — static control flow on KK is fine.
    CPVMCL = cond.CL - cond.CPV
    EPS = params.RD / cond.RV
    EPSI = 1.0 / EPS
    qNK = Q[NK]
    tNK = T[NK]
    gzNK = GZ[NK]
    AH0 = (
        (params.CPD * (1.0 - qNK) + cond.CL * qNK) * tNK
        + qNK * (cond.LV0 - CPVMCL * (tNK - cond.T_freeze))
        + gzNK
    )
    CPP = params.CPD * (1.0 - qNK) + qNK * cond.CPV
    CPINV = 1.0 / CPP

    lvl = jnp.arange(NL + 1)
    T_nl = T[:NL + 1]
    QS_nl = QS[:NL + 1]
    P_nl = P[:NL + 1]
    GZ_nl = GZ[:NL + 1]
    TVP_nl = TVP[:NL + 1]
    TPK_nl = TPK[:NL + 1]
    CLW_nl = CLW[:NL + 1]

    if KK == 1:
        CLW_nl = jnp.where(lvl < ICB, 0.0, CLW_nl)
        tpk_dry = tNK - (GZ_nl - gzNK) * CPINV
        in_dry = (lvl >= NK) & (lvl < ICB)
        TPK_nl = jnp.where(in_dry, tpk_dry, TPK_nl)
        TVP_nl = jnp.where(in_dry, tpk_dry * (1.0 + qNK * EPSI), TVP_nl)

    # Outer saturation-adjustment over levels 0..NL (2 fixed Newton iterations unrolled).
    loop_mask = (lvl == ICB) if KK == 1 else (lvl > ICB)
    TG = T_nl
    QG = QS_nl
    ALV = cond.LV0 - CPVMCL * (T_nl - cond.T_freeze)
    for _ in range(2):
        S = 1.0 / (params.CPD + ALV * ALV * QG / (cond.RV * T_nl * T_nl))
        AHG = (
            params.CPD * TG
            + (cond.CL - params.CPD) * qNK * T_nl
            + ALV * QG
            + GZ_nl
        )
        TG = jnp.maximum(TG + S * (AH0 - AHG), 35.0)
        ES = _sat_vap_pressure_jax(TG, cond.species_id)
        QG = EPS * ES / jnp.maximum(P_nl - ES * (1.0 - EPS), 1.0e-30)
    tpk_sat = (AH0 - (cond.CL - params.CPD) * qNK * T_nl - GZ_nl - ALV * QG) / params.CPD
    TPK_nl = jnp.where(loop_mask, tpk_sat, TPK_nl)
    CLW_nl = jnp.where(loop_mask, jnp.maximum(0.0, qNK - QG), CLW_nl)
    TVP_nl = jnp.where(
        loop_mask,
        tpk_sat * (1.0 + (QG / jnp.maximum(1.0 - qNK, 1.0e-30)) * EPSI),
        TVP_nl,
    )

    TVP = jnp.concatenate([TVP_nl, TVP[NL + 1:]])
    TPK = jnp.concatenate([TPK_nl, TPK[NL + 1:]])
    CLW = jnp.concatenate([CLW_nl, CLW[NL + 1:]])
    return TVP, TPK, CLW


@functools.partial(jax.jit, static_argnums=(7, 8, 9, 13))
def _convect_functional_jax(
    T_in, Q_in, QS_in, U_in, V_in, P_in, PH_in,
    ND, NL, NTRA, DELT, CBMF_in, params, cond,
):
    T = jnp.array(T_in)
    Q = jnp.array(Q_in)
    QS = jnp.array(QS_in)
    U = jnp.array(U_in)
    V = jnp.array(V_in)
    P = jnp.array(P_in)
    PH = jnp.array(PH_in)
    CBMF = CBMF_in

    CPVMCL = cond.CL - cond.CPV
    EPS = params.RD / cond.RV
    EPSI = 1.0 / EPS
    GINV = 1.0 / params.G
    DELTI = 1.0 / DELT
    PRECIP = 0.0
    WD = 0.0
    TPRIME = 0.0
    QPRIME = 0.0
    IFLAG = jnp.int32(0)
    OUTCAPE = 0.0

    # --- Block 2: Thermodynamic profiles ---
    T_nl = T[:NL + 1]
    Q_nl = Q[:NL + 1]
    P_nl = P[:NL + 1]

    TV_nl = T_nl * (1.0 + Q_nl * EPSI - Q_nl)
    dgz = 0.5 * params.RD * (TV_nl[1:] + TV_nl[:-1]) * (P_nl[:-1] - P_nl[1:]) / PH[1:NL + 1]
    GZ_nl = jnp.concatenate([jnp.zeros(1), jnp.cumsum(dgz)])
    CPN_nl = params.CPD * (1.0 - Q_nl) + cond.CPV * Q_nl
    H_nl = T_nl * CPN_nl + GZ_nl
    LV_nl = cond.LV0 - CPVMCL * (T_nl - cond.T_freeze)
    HM_nl = (
        (params.CPD * (1.0 - Q_nl) + cond.CL * Q_nl) * (T_nl - T_nl[0])
        + LV_nl * Q_nl
        + GZ_nl
    )
    pad_nl_to_nd1 = jnp.zeros(ND - NL)
    GZ = jnp.concatenate([GZ_nl, pad_nl_to_nd1])
    CPN = jnp.concatenate([CPN_nl, pad_nl_to_nd1])
    H = jnp.concatenate([H_nl, pad_nl_to_nd1])
    LV = jnp.concatenate([LV_nl, pad_nl_to_nd1])
    HM = jnp.concatenate([HM_nl, pad_nl_to_nd1])
    TV = jnp.concatenate([TV_nl, pad_nl_to_nd1])

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
    RH_safe = jnp.maximum(RH, 1.0e-30)
    if cond.species_id == 0:
        CHI = T[NK] / jnp.maximum(1669.0 - 122.0 * RH - T[NK], 1.0e-30)
        PLCL = jnp.where(Q[NK] > 0.0, P[NK] * (RH_safe ** CHI), 0.0)
    else:
        chi = cond.RV * T[NK] / cond.LV0
        PLCL = jnp.where(Q[NK] > 0.0, P[NK] * (RH_safe ** chi), 0.0)

    # Exit condition 2: PLCL out of range
    exit2 = (PLCL < 200.0) | (PLCL >= 2000.0)
    IFLAG = jnp.where(exit2 & active, jnp.int32(2), IFLAG)
    active = active & ~exit2

    # Vectorized ICB: first index > NK where P < PLCL, default NL-1
    # P has length ND; use indices[:ND] for comparison
    p_indices = jnp.arange(ND)
    icb_candidates = (p_indices > NK) & (p_indices <= NL) & (P < PLCL)
    ICB = jnp.where(jnp.any(icb_candidates), jnp.argmax(icb_candidates), NL - 1)

    # Exit condition 3: ICB >= NL-1
    exit3 = ICB >= (NL - 1)
    IFLAG = jnp.where(exit3 & active, jnp.int32(3), IFLAG)
    active = active & ~exit3

    # --- Block 5: TLIFT ---
    TVP = jnp.zeros(ND)
    TP = jnp.zeros(ND)
    CLW = jnp.zeros(ND)
    TVP, TP, CLW = _tlift_functional_jax(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params, cond
    )
    lvl_nd = jnp.arange(ND)
    in_nk_icb = (lvl_nd <= NL) & (lvl_nd >= NK) & (lvl_nd <= ICB)
    TVP = jnp.where(in_nk_icb, TVP - TP * Q[NK], TVP)

    # Early return: CBMF == 0 and TVP[ICB] too cold → fold into active flag
    exit4 = (CBMF == 0.0) & (TVP[ICB] <= (TV[ICB] - params.DTMAX))
    active = active & ~exit4

    IFLAG = jnp.where(active & (IFLAG != 4), jnp.int32(1), IFLAG)

    TVP, TP, CLW = _tlift_functional_jax(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 2, TVP, TP, CLW, params, cond
    )

    # --- Block 6: EP/SIGP ---
    SIGP = jnp.full(ND, params.SIGS)
    TCA = TP[:NL + 1] - cond.T_freeze
    ELACRIT = jnp.where(
        TCA >= 0.0,
        params.ELCRIT,
        params.ELCRIT * (1.0 - TCA / params.TLCRIT),
    )
    ELACRIT = jnp.maximum(ELACRIT, 0.0)
    EPMAX = 0.999
    ep_nl = EPMAX * (1.0 - ELACRIT / jnp.maximum(CLW[:NL + 1], 1.0e-8))
    ep_nl = jnp.maximum(jnp.minimum(ep_nl, EPMAX), 0.0)
    EP = jnp.concatenate([
        jnp.where(jnp.arange(NL + 1) > NK, ep_nl, 0.0),
        jnp.zeros(ND - NL - 1),
    ])

    in_above_icb = (lvl_nd <= NL) & (lvl_nd > ICB)
    TVP = jnp.where(in_above_icb, TVP - TP * Q[NK], TVP)

    # --- Block 7: Setup matrices ---
    HP = jnp.array(H)
    WATER = jnp.zeros(ND + 1)
    EVAP = jnp.zeros(ND + 1)
    WT = jnp.full(ND + 1, params.OMTSNOW)
    MP = jnp.zeros(ND + 1)
    M = jnp.zeros(ND + 1)
    LVCP = jnp.concatenate([LV[:NL + 1] / CPN[:NL + 1], jnp.zeros(ND - NL)])

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
    INB1 = jnp.where(
        jnp.any(inb1_cand),
        jnp.max(jnp.where(inb1_cand, cape_lvl, -1)) + 1,
        ICB + 1,
    )

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
    OUTCAPE = jnp.where(active, CAPE, 0.0)

    # --- Block 9: HP update, CBMF (JIT-compatible) ---
    # HP update: full range with mask for [ICB, INB+1)
    hp_lvl = jnp.arange(ND + 1)
    hp_mask = (hp_lvl >= ICB) & (hp_lvl <= INB)
    # Pad T, EP, CLW from ND to ND+1 for uniform indexing
    _T_hp = jnp.concatenate([T, jnp.zeros(1)])
    _EP_hp = jnp.concatenate([EP, jnp.zeros(1)])
    _CLW_hp = jnp.concatenate([CLW, jnp.zeros(1)])
    HP_new = H[NK] + (LV + (params.CPD - cond.CPV) * _T_hp) * _EP_hp * _CLW_hp
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
    BF2_vec = 1.0 + LV * LV * QS_pad / (cond.RV * T_pad * T_pad * params.CPD)  # (ND+1,)

    # First-pass SIJ computation for all (i,j) pairs
    # ANUM[i,j] = H[j] - HP[i] + (CPV - CPD) * T[j] * (QTI[i] - Q[j])
    ANUM_2d = H[None, :] - HP[:, None] + (cond.CPV - params.CPD) * T_pad[None, :] * (QTI_vec[:, None] - Q_pad[None, :])
    # DENOM[i,j] = H[i] - HP[i] + (CPD - CPV) * (Q[i] - QTI[i]) * T[j]
    DENOM_2d = H[:, None] - HP[:, None] + (params.CPD - cond.CPV) * (Q_pad[:, None] - QTI_vec[:, None]) * T_pad[None, :]
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
    i_sc = jnp.arange(NL + 1)
    j_sc = jnp.arange(NL + 1)
    i_active_scrit = (i_sc >= ICB + 1) & (i_sc <= INB) & (NENT_count[:NL + 1] > 0)

    QP1 = Q_pad[NK] - EP_pad[:NL + 1] * CLW_pad[:NL + 1]
    ANUM_s = H[:NL + 1] - HP[:NL + 1] - LV[:NL + 1] * (QP1 - QS_pad[:NL + 1])
    DENOM_s = H[:NL + 1] - HP[:NL + 1] + LV[:NL + 1] * (Q_pad[:NL + 1] - QP1)
    DENOM_s = jnp.where(jnp.abs(DENOM_s) < 0.01, 0.01, DENOM_s)
    SCRIT = ANUM_s / DENOM_s
    ALT = QP1 - QS_pad[:NL + 1] + SCRIT * (Q_pad[:NL + 1] - QP1)
    SCRIT = jnp.maximum(jnp.where(ALT < 0.0, 1.0, SCRIT), 0.0)

    sij_sub = SIJ[:NL + 1, :NL + 1]
    sij_next = SIJ[:NL + 1, 1:NL + 2]
    sij_prev_gt = jnp.concatenate([SIJ[:NL + 1, :1], SIJ[:NL + 1, :NL]], axis=1)
    sij_prev_le = jnp.concatenate([jnp.zeros((NL + 1, 1)), SIJ[:NL + 1, :NL]], axis=1)

    j_valid = (j_sc[None, :] >= ICB) & (j_sc[None, :] <= INB)
    is_ent = (sij_sub > 0.0) & (sij_sub < 0.9) & j_valid
    gt_mask = j_sc[None, :] > i_sc[:, None]
    scrit_col = SCRIT[:, None]

    smid_gt = jnp.minimum(sij_sub, scrit_col)
    smin_cand = jnp.where(gt_mask & (sij_next < smid_gt) & is_ent, smid_gt, 1.0)
    smin_incl = jax.lax.associative_scan(jnp.minimum, smin_cand, axis=1)
    smin_before = jnp.concatenate([jnp.ones((NL + 1, 1)), smin_incl[:, :-1]], axis=1)
    update_minmax = (smid_gt < smin_before) & (sij_next < smid_gt)
    sjmax_gt = jnp.where(update_minmax, jnp.minimum(jnp.minimum(sij_next, sij_sub), scrit_col), smid_gt)
    sjmin_gt = jnp.where(update_minmax, jnp.minimum(jnp.maximum(sij_prev_gt, sij_sub), scrit_col), smid_gt)

    smid_le = jnp.maximum(sij_sub, scrit_col)
    sjmax_le = jnp.maximum(sij_next, scrit_col)
    sjmin_le = jnp.maximum(sij_prev_le, scrit_col)

    SMID = jnp.where(gt_mask, smid_gt, smid_le)
    SJMAX = jnp.where(gt_mask, sjmax_gt, sjmax_le)
    SJMIN = jnp.where(gt_mask, sjmin_gt, sjmin_le)

    dp_sc = PH[:NL + 1] - PH[1:NL + 2]
    w_sub = jnp.where(
        is_ent,
        (jnp.abs(SJMAX - SMID) + jnp.abs(SJMIN - SMID)) * dp_sc[None, :],
        0.0,
    )
    weights = jnp.concatenate([w_sub, jnp.zeros((NL + 1, ND - NL))], axis=1)
    ASIJ_acc = jnp.maximum(jnp.sum(weights, axis=1), 1.0e-21)
    MENT_rows = MENT[:NL + 1, :] * weights / ASIJ_acc[:, None]
    ment_sum = jnp.sum(MENT_rows, axis=1)

    tiny = (ment_sum < 1.0e-18)[:, None]
    diag_nl = i_sc[:, None] == jnp.arange(ND + 1)[None, :]
    MENT_fb = jnp.where(diag_nl, M[:NL + 1, None], 0.0)
    QENT_fb = jnp.where(diag_nl, QP1[:, None], QENT[:NL + 1, :])
    UENT_fb = jnp.where(diag_nl, U_pad[NK], UENT[:NL + 1, :])
    VENT_fb = jnp.where(diag_nl, V_pad[NK], VENT[:NL + 1, :])
    ELIJ_fb = jnp.where(diag_nl, CLW_pad[:NL + 1, None], ELIJ[:NL + 1, :])
    SIJ_fb = jnp.where(diag_nl, 1.0, SIJ[:NL + 1, :])

    act_row = i_active_scrit[:, None]
    MENT = jnp.concatenate(
        [jnp.where(act_row, jnp.where(tiny, MENT_fb, MENT_rows), MENT[:NL + 1, :]), MENT[NL + 1:, :]],
        axis=0,
    )
    QENT = jnp.concatenate(
        [jnp.where(act_row & tiny, QENT_fb, QENT[:NL + 1, :]), QENT[NL + 1:, :]],
        axis=0,
    )
    UENT = jnp.concatenate(
        [jnp.where(act_row & tiny, UENT_fb, UENT[:NL + 1, :]), UENT[NL + 1:, :]],
        axis=0,
    )
    VENT = jnp.concatenate(
        [jnp.where(act_row & tiny, VENT_fb, VENT[:NL + 1, :]), VENT[NL + 1:, :]],
        axis=0,
    )
    ELIJ = jnp.concatenate(
        [jnp.where(act_row & tiny, ELIJ_fb, ELIJ[:NL + 1, :]), ELIJ[NL + 1:, :]],
        axis=0,
    )
    SIJ = jnp.concatenate(
        [jnp.where(act_row & tiny, SIJ_fb, SIJ[:NL + 1, :]), SIJ[NL + 1:, :]],
        axis=0,
    )

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
    WT_dd = jnp.where(T_pad > cond.T_freeze, params.OMTRAIN, params.OMTSNOW)  # (ND+1,)
    COEFF_dd = jnp.where(T_pad > cond.T_freeze, params.COEFFR, params.COEFFS)  # (ND+1,)
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
             + QP_above * (LV[i + 1] + T_pad[ip1] * (cond.CL - params.CPD))
             + params.CPD * (T_pad[ip1] - T_pad[i]))
            / jnp.maximum(LV[i] + T_pad[i] * (cond.CL - params.CPD), 1e-30)
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

    # Gather downdraft scan outputs directly by level (scan step k = INB - level)
    k_for_lvl = INB - lvl
    k_safe = jnp.clip(k_for_lvl, 0, NL)
    scatter_valid = ep_active & (k_for_lvl >= 0) & (k_for_lvl <= NL)
    WATER = jnp.where(scatter_valid, WATER_sc[k_safe], WATER)
    WT = jnp.where(scatter_valid, WT_sc[k_safe], WT)
    MP = jnp.where(scatter_valid, MP_sc[k_safe], MP)
    EVAP = jnp.where(scatter_valid, EVAP_sc[k_safe], EVAP)
    QP = jnp.where(scatter_valid, QP_sc[k_safe], QP)
    UP = jnp.where(scatter_valid, UP_sc[k_safe], UP)
    VP = jnp.where(scatter_valid, VP_sc[k_safe], VP)

    inb_sc_idx = jnp.clip(INB, 0, NL)
    PRECIP = jnp.where(
        ep_active & (INB <= NL),
        PRECIP + WT_sc[inb_sc_idx] * params.SIGD * WATER_sc[inb_sc_idx]
        * 3600.0 * 24000.0 / (cond.ROWL * params.G),
        PRECIP,
    )

    # --- Blocks 13–15: Tendency accumulation (differentiable) ---
    WD = params.BETA * jnp.abs(MP[ICB]) * 0.01 * params.RD * T_pad[ICB] / (params.SIGD * P[ICB])
    QPRIME = 0.5 * (QP[0] - Q_pad[0])
    TPRIME = cond.LV0 * QPRIME / params.CPD

    # --- Block 13: FT, FQ, FU, FV computation ---
    DPINV_all = 0.01 / (PH[:ND] - PH[1:ND + 1])
    CPINV_all = 1.0 / jnp.maximum(CPN[:ND], 1e-30)

    # --- Surface level (i=0) ---
    DPINV_0 = DPINV_all[0]
    lvl_m = jnp.arange(ND + 1)
    AM = jnp.where(NK == 0, jnp.sum(jnp.where((lvl_m >= 1) & (lvl_m <= INB), M, 0.0)), 0.0)

    IFLAG = jnp.where(active & ((2.0 * params.G * DPINV_0 * AM) >= DELTI), jnp.int32(4), IFLAG)

    ft_0 = (
        params.G * DPINV_0 * AM * (T_pad[1] - T_pad[0] + (GZ[1] - GZ[0]) / CPN[0])
        - LVCP[0] * params.SIGD * EVAP[0]
        + params.SIGD * WT[1] * (cond.CL - params.CPD) * WATER[1]
        * (T_pad[1] - T_pad[0]) * DPINV_0 / CPN[0]
    )
    j_mask_surf = (lvl_m >= 1) & (lvl_m <= INB)
    fq_0 = (
        params.G * MP[1] * (QP[1] - Q_pad[0]) * DPINV_0
        + params.SIGD * EVAP[0]
        + params.G * AM * (Q_pad[1] - Q_pad[0]) * DPINV_0
        + params.G * DPINV_0 * jnp.sum(
            jnp.where(j_mask_surf, MENT[:, 0] * (QENT[:, 0] - Q_pad[0]), 0.0)
        )
    )
    fu_0 = (
        params.G * DPINV_0 * (MP[1] * (UP[1] - U_pad[0]) + AM * (U_pad[1] - U_pad[0]))
        + params.G * DPINV_0 * jnp.sum(
            jnp.where(j_mask_surf, MENT[:, 0] * (UENT[:, 0] - U_pad[0]), 0.0)
        )
    )
    fv_0 = (
        params.G * DPINV_0 * (MP[1] * (VP[1] - V_pad[0]) + AM * (V_pad[1] - V_pad[0]))
        + params.G * DPINV_0 * jnp.sum(
            jnp.where(j_mask_surf, MENT[:, 0] * (VENT[:, 0] - V_pad[0]), 0.0)
        )
    )

    # --- Interior levels (i=1..ND-1) ---
    i_vec = jnp.arange(1, ND)
    i_act = i_vec <= INB
    DPINV_int = DPINV_all[1:ND]
    CPINV_int = CPINV_all[1:ND]

    m_above_mask = (lvl_m[None, :] >= i_vec[:, None] + 1) & (lvl_m[None, :] <= INB + 1)
    AMP1_mass = jnp.where(
        i_vec >= NK,
        jnp.sum(jnp.where(m_above_mask, M[None, :], 0.0), axis=1),
        0.0,
    )
    k_le_i = lvl_m[None, :, None] <= i_vec[:, None, None]
    jj_gt_i = (lvl_m[None, None, :] >= i_vec[:, None, None] + 1) & (lvl_m[None, None, :] <= INB + 1)
    AMP1_ent = jnp.sum(jnp.where(k_le_i & jj_gt_i, MENT[None, :, :], 0.0), axis=(1, 2))
    AMP1 = AMP1_mass + AMP1_ent

    IFLAG = jnp.where(active & jnp.any(i_act & ((2.0 * params.G * DPINV_int * AMP1) >= DELTI)), jnp.int32(4), IFLAG)

    jj_ge_i = (lvl_m[None, :, None] >= i_vec[:, None, None]) & (lvl_m[None, :, None] <= INB)
    k_lt_i = lvl_m[None, None, :] < i_vec[:, None, None]
    AD = jnp.sum(jnp.where(jj_ge_i & k_lt_i, MENT[None, :, :], 0.0), axis=(1, 2))

    T_i = T_pad[1:ND]
    T_ip1 = T_pad[2:ND + 1]
    T_im1 = T_pad[:ND - 1]
    GZ_i = GZ[1:ND]
    GZ_ip1 = GZ[2:ND + 1]
    GZ_im1 = GZ[:ND - 1]
    Q_i = Q_pad[1:ND]
    Q_ip1 = Q_pad[2:ND + 1]
    Q_im1 = Q_pad[:ND - 1]
    U_i = U_pad[1:ND]
    U_ip1 = U_pad[2:ND + 1]
    U_im1 = U_pad[:ND - 1]
    V_i = V_pad[1:ND]
    V_ip1 = V_pad[2:ND + 1]
    V_im1 = V_pad[:ND - 1]

    MENT_diag = jnp.diag(MENT)[1:ND]
    QENT_diag = jnp.diag(QENT)[1:ND]

    ft_adv = (
        params.G * DPINV_int * (
            AMP1 * (T_ip1 - T_i + (GZ_ip1 - GZ_i) * CPINV_int)
            - AD * (T_i - T_im1 + (GZ_i - GZ_im1) * CPINV_int)
        )
        - params.SIGD * LVCP[1:ND] * EVAP[1:ND]
    )
    ft_self = (
        params.G * DPINV_int * MENT_diag
        * (HP[1:ND] - H[1:ND] + T_i * (cond.CPV - params.CPD) * (Q_i - QENT_diag))
        * CPINV_int
    )
    ft_water = (
        params.SIGD * WT[2:ND + 1] * (cond.CL - params.CPD)
        * WATER[2:ND + 1] * (T_ip1 - T_i)
        * DPINV_int * CPINV_int
    )
    ft_int = jnp.where(i_act, ft_adv + ft_self + ft_water, 0.0)

    fq_adv = params.G * DPINV_int * (AMP1 * (Q_ip1 - Q_i) - AD * (Q_i - Q_im1))
    fu_adv = params.G * DPINV_int * (AMP1 * (U_ip1 - U_i) - AD * (U_i - U_im1))
    fv_adv = params.G * DPINV_int * (AMP1 * (V_ip1 - V_i) - AD * (V_i - V_im1))

    k_below = lvl_m[:, None] < i_vec[None, :]
    k_above = (lvl_m[:, None] >= i_vec[None, :]) & (lvl_m[:, None] <= INB)
    MENT_cols = MENT[:, 1:ND]
    QENT_cols = QENT[:, 1:ND]
    UENT_cols = UENT[:, 1:ND]
    VENT_cols = VENT[:, 1:ND]
    AWAT = jnp.maximum(ELIJ[:, 1:ND] - (1.0 - EP_pad[None, 1:ND]) * CLW_pad[None, 1:ND], 0.0)

    fq_ment_below = params.G * DPINV_int * jnp.sum(
        jnp.where(k_below, MENT_cols * (QENT_cols - AWAT - Q_i[None, :]), 0.0), axis=0
    )
    fu_ment_below = params.G * DPINV_int * jnp.sum(
        jnp.where(k_below, MENT_cols * (UENT_cols - U_i[None, :]), 0.0), axis=0
    )
    fv_ment_below = params.G * DPINV_int * jnp.sum(
        jnp.where(k_below, MENT_cols * (VENT_cols - V_i[None, :]), 0.0), axis=0
    )

    fq_ment_above = params.G * DPINV_int * jnp.sum(
        jnp.where(k_above, MENT_cols * (QENT_cols - Q_i[None, :]), 0.0), axis=0
    )
    fu_ment_above = params.G * DPINV_int * jnp.sum(
        jnp.where(k_above, MENT_cols * (UENT_cols - U_i[None, :]), 0.0), axis=0
    )
    fv_ment_above = params.G * DPINV_int * jnp.sum(
        jnp.where(k_above, MENT_cols * (VENT_cols - V_i[None, :]), 0.0), axis=0
    )

    MP_i = MP[1:ND]
    MP_ip1 = MP[2:ND + 1]
    fq_down = (
        params.SIGD * EVAP[1:ND]
        + params.G * (MP_ip1 * (QP[2:ND + 1] - Q_i) - MP_i * (QP[1:ND] - Q_im1)) * DPINV_int
    )
    fu_down = (
        params.G * (MP_ip1 * (UP[2:ND + 1] - U_i) - MP_i * (UP[1:ND] - U_im1)) * DPINV_int
    )
    fv_down = (
        params.G * (MP_ip1 * (VP[2:ND + 1] - V_i) - MP_i * (VP[1:ND] - V_im1)) * DPINV_int
    )

    fq_int = jnp.where(i_act, fq_adv + fq_ment_below + fq_ment_above + fq_down, 0.0)
    fu_int = jnp.where(i_act, fu_adv + fu_ment_below + fu_ment_above + fu_down, 0.0)
    fv_int = jnp.where(i_act, fv_adv + fv_ment_below + fv_ment_above + fv_down, 0.0)

    FT = jnp.concatenate([jnp.atleast_1d(ft_0), ft_int])
    FQ = jnp.concatenate([jnp.atleast_1d(fq_0), fq_int])
    FU = jnp.concatenate([jnp.atleast_1d(fu_0), fu_int])
    FV = jnp.concatenate([jnp.atleast_1d(fv_0), fv_int])

    # --- Block 14: FRAC smoothing at INB (differentiable) ---
    ph_ratio = (PH[INB] - PH[INB + 1]) / (PH[INB - 1] - PH[INB])
    FQOLD = FQ[INB]
    FQ = FQ.at[INB].set(FQOLD * (1.0 - FRAC))
    FQ = FQ.at[INB - 1].add(FRAC * FQOLD * ph_ratio * LV[INB] / jnp.maximum(LV[INB - 1], 1.0e-30))
    FTOLD = FT[INB]
    FT = FT.at[INB].set(FTOLD * (1.0 - FRAC))
    FT = FT.at[INB - 1].add(FRAC * FTOLD * ph_ratio * CPN[INB] / jnp.maximum(CPN[INB - 1], 1.0e-30))
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

    return FT, FQ, FU, FV, PRECIP, WD, TPRIME, QPRIME, CBMF, OUTCAPE, IFLAG


def emanuel_convect_jax(
    t, q, qs, u, v, p, ph, nlev, nl, ntra, delt, cbmf, params, cond,
):
    def column_call(T, Q, QS, U, V, P, PH, CBMF):
        return _convect_functional_jax(T, Q, QS, U, V, P, PH, nlev, nl, ntra, delt, CBMF, params, cond)

    # vmap over columns: input arrays are (nlev, ncol), vmap over axis 1 → per-column calls
    # cbmf is (ncol,), vmap over axis 0
    results = jax.vmap(column_call, in_axes=(1, 1, 1, 1, 1, 1, 1, 0))(
        t, q, qs, u, v, p, ph, cbmf
    )

    # vmap returns (ncol, nlev) for 1D outputs → transpose to (nlev, ncol)
    ft, fq, fu, fv = results[0].T, results[1].T, results[2].T, results[3].T
    precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results[4:]
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
