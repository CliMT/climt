"""Experiment T2 — isolate PF table node-fidelity from interpolation error.

Experiment #4 found the shipped PF wing band (500-630 cm^-1) is 4.15x more
opaque than LBL at the strat hotspot (p=97 hPa, T=147 K). Experiment T1
showed linear-in-k T-interpolation over the coarse 110<->170 K Chaverot
bracket explains only ~1.5x of that. This script isolates the remaining
factor by comparing the PF table k DIRECTLY at its T-nodes (110, 170 K) at
the p=1e4 Pa node (no p-interp, no T-interp) against LBL band-mean k at the
same (T, p).

  ratio_node = PF_node(T) / LBL(T)        -> pure table/source fidelity
  ratio_147  = PF_linterp(147) / LBL(147) -> node fidelity + interp curvature

If ratio_node ~ 1 at both nodes, the whole bias is interpolation (loader
fix). If ratio_node >> 1, the bias is baked into the table (bin_down/GL
rebin or the Chaverot source CDF) and needs a regen / better source.

LBL band-mean k comes from debug_data/lbl_kT_curvature.npz (written by
lbl_kT_curvature.py). Run this from the `radiation` env (netCDF4 + numpy).
"""
from __future__ import annotations

import os
import sys

import netCDF4 as nc
import numpy as np

_CK_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "climt", "_data",
    "picket_fence", "correlated_k")
# Optional CLI arg: table basename (no .nc) or full path. Defaults to the
# shipped Chaverot co2refined table.
_arg = sys.argv[1] if len(sys.argv) > 1 else "earth_low_res_lw_co2refined_gl"
PF_TABLE = _arg if _arg.endswith(".nc") else os.path.join(_CK_DIR, _arg + ".nc")
LBL_NPZ = os.path.join(
    os.path.dirname(__file__), "..", "..", "debug_data", "lbl_kT_curvature.npz",
)

P_NODE = 1.0e4          # Pa — exact Chaverot p-node (~100 hPa)
T_LO, T_HI, T_TGT = 110.0, 170.0, 147.0
X_H2O_TARGET = 1e-6     # dry — LBL was CO2-only

# PF band edges -> LBL npz key
BAND_TO_KEY = {
    (500.0, 630.0): "k_wing_500_630",
    (630.0, 700.0): "k_core_630_700",
    (700.0, 800.0): "k_wing_700_800",
    (1800.0, 3250.0): "k_v3_1800_3250",
}


def main():
    with nc.Dataset(PF_TABLE) as ds:
        kc = ds["k_coefficients"][...]       # (1, nband, ngpt, nT, nP, nX)
        gw = ds["gpoint_weights"][...]        # (nband, ngpt)
        be = ds["band_wavenumber_limits"][...]
        Tg = ds["temperature_grid"][...]
        lpg = ds["pressure_grid_log"][...]
        xg = ds["h2o_vmr_grid"][...]
    # g-integrate to band-mean k(b, T, P, X)
    k_band = np.einsum("bgtpx,bg->btpx", kc[0], gw)

    iT_lo = int(np.argmin(np.abs(Tg - T_LO)))
    iT_hi = int(np.argmin(np.abs(Tg - T_HI)))
    iP = int(np.argmin(np.abs(np.exp(lpg) - P_NODE)))
    iX = int(np.argmin(np.abs(xg - X_H2O_TARGET)))
    assert abs(Tg[iT_lo] - T_LO) < 1 and abs(Tg[iT_hi] - T_HI) < 1
    assert abs(np.exp(lpg)[iP] - P_NODE) < 1.0

    w = (T_TGT - T_LO) / (T_HI - T_LO)

    lbl = np.load(LBL_NPZ)
    T_lbl = lbl["T_grid"]
    jl = {110.0: int(np.argmin(np.abs(T_lbl - 110.0))),
          147.0: int(np.argmin(np.abs(T_lbl - 147.0))),
          170.0: int(np.argmin(np.abs(T_lbl - 170.0)))}

    print(f"PF table: {os.path.basename(PF_TABLE)}")
    print(f"p node = {np.exp(lpg)[iP]:.0f} Pa, T nodes = {Tg[iT_lo]:.0f}/"
          f"{Tg[iT_hi]:.0f} K, X = {xg[iX]:.1e}\n")
    print(f"{'band':<14} {'PF@110':>9} {'LBL@110':>9} {'node110':>8}  "
          f"{'PF@170':>9} {'LBL@170':>9} {'node170':>8}  "
          f"{'PFinterp147':>11} {'LBL@147':>9} {'r147':>6}")

    for b in range(be.shape[0]):
        edges = (float(be[b, 0]), float(be[b, 1]))
        key = BAND_TO_KEY.get(edges)
        if key is None:
            continue
        pf110 = float(k_band[b, iT_lo, iP, iX])
        pf170 = float(k_band[b, iT_hi, iP, iX])
        pf147_lin = (1 - w) * pf110 + w * pf170
        l110 = float(lbl[key][jl[110.0]])
        l147 = float(lbl[key][jl[147.0]])
        l170 = float(lbl[key][jl[170.0]])
        r110 = pf110 / l110 if l110 > 0 else np.nan
        r170 = pf170 / l170 if l170 > 0 else np.nan
        r147 = pf147_lin / l147 if l147 > 0 else np.nan
        name = f"{int(edges[0])}-{int(edges[1])}"
        print(f"{name:<14} {pf110:>9.2e} {l110:>9.2e} {r110:>8.2f}  "
              f"{pf170:>9.2e} {l170:>9.2e} {r170:>8.2f}  "
              f"{pf147_lin:>11.2e} {l147:>9.2e} {r147:>6.2f}")

    print("\nInterpretation:")
    print("  node110/node170 ~ 1  -> table faithful at nodes; bias is interp")
    print("  node110/node170 >> 1 -> bias baked into table (rebin/source)")


if __name__ == "__main__":
    main()
