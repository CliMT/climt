"""Experiment #11 — PF table band-mean k vs LBL at WARM nodes (230/290 K).

The RCE test (exp #10) left a tropospheric warm bias: PF is +5 to +16 K warmer
than RRTMG from the surface to ~50 hPa (peak ~+16 K at 600-700 hPa), shared by
all PF tables. exp #4 attributed this to a warm-T CO2 over-opacity (PF/LBL ~1.8x
at 288 K), but exps #7/#8 only validated k-fidelity at the COLD nodes (110/170
K). This script closes that gap: compare the PF table band-mean k DIRECTLY at the
warm T-nodes (230, 290 K) at tropospheric p-nodes (1e4, 1e5 Pa) against LBL.

LBL band-means come from debug_data/chaverot_vs_lbl_lbl.npz (exp #9 envelope;
grid-converged at p=1e4 and 1e5 Pa per that entry's caveat). No interpolation:
we read PF k at exact (T, p) nodes and g-integrate with the table's own weights.

Run: python scripts/experiments/pf_warm_node_fidelity.py   (climt env is fine)
"""
from __future__ import annotations

import os
import numpy as np
import netCDF4 as nc

_HERE = os.path.dirname(__file__)
_CK = os.path.join(_HERE, "..", "..", "climt", "_data", "picket_fence", "correlated_k")
_LBL = os.path.join(_HERE, "..", "..", "debug_data", "chaverot_vs_lbl_lbl.npz")

TABLES = {
    "Chaverot 6-band (co2refined_gl)": "earth_low_res_lw_co2refined_gl",
    "linepyline 6-band (fix)":         "earth_low_res_lw_co2refined_linepyline",
}
# (band edge tuple) -> LBL npz key
BAND_KEY = {
    (500.0, 630.0): "band_500_630",
    (630.0, 700.0): "band_630_700",
    (700.0, 800.0): "band_700_800",
}
T_TARGETS = [230.0, 290.0]
P_TARGETS = [1.0e4, 1.0e5]   # Pa — converged LBL nodes (strat control + troposphere)


def band_mean_pf(table_path, T_tgt, p_tgt_pa):
    """PF table band-mean k (m^2/kg) at the exact (T, p) nodes, dry, per band."""
    with nc.Dataset(table_path) as ds:
        kc = ds["k_coefficients"][...]        # (1, nband, ngpt, nT, nP, nX)
        gw = ds["gpoint_weights"][...]         # (nband, ngpt)
        be = ds["band_wavenumber_limits"][...]
        Tg = np.asarray(ds["temperature_grid"][...])
        lpg = np.asarray(ds["pressure_grid_log"][...])
        # pressure grid stored as natural log of Pa (ln) -> convert to Pa.
        pg = np.exp(lpg)
    iT = int(np.argmin(np.abs(Tg - T_tgt)))
    iP = int(np.argmin(np.abs(pg - p_tgt_pa)))
    iX = 0  # dry (smallest H2O vmr node); LBL is CO2-only
    out = {}
    for b in range(be.shape[0]):
        edge = (float(be[b, 0]), float(be[b, 1]))
        if edge not in BAND_KEY:
            continue
        kg = kc[0, b, :, iT, iP, iX]
        w = gw[b]
        out[edge] = float(np.sum(kg * w) / np.sum(w))
    return out, Tg[iT], pg[iP]


def main():
    L = np.load(_LBL)
    Tn, pn = L["T_nodes"], L["p_nodes"]

    for name, base in TABLES.items():
        path = os.path.join(_CK, base + ".nc")
        print(f"\n=== {name} ===")
        print(f"{'band':14s} {'T':>5s} {'p_hPa':>7s} {'PF k':>11s} {'LBL k':>11s} {'PF/LBL':>8s}")
        for T_tgt in T_TARGETS:
            for p_tgt in P_TARGETS:
                pf, Tnode, pnode = band_mean_pf(path, T_tgt, p_tgt)
                iTl = int(np.argmin(np.abs(Tn - T_tgt)))
                iPl = int(np.argmin(np.abs(pn - p_tgt)))
                for edge, key in BAND_KEY.items():
                    if edge not in pf:
                        continue
                    lbl = float(L[key][iTl, iPl])
                    r = pf[edge] / lbl if lbl > 0 else np.inf
                    tag = "core" if edge == (630.0, 700.0) else "wing"
                    print(f"{edge[0]:.0f}-{edge[1]:.0f} ({tag:4s}) {Tnode:5.0f} "
                          f"{pnode/100:7.1f} {pf[edge]:11.3e} {lbl:11.3e} {r:8.2f}")


if __name__ == "__main__":
    main()
