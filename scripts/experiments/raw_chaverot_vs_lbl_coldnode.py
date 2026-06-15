"""Experiment #6 — raw Chaverot R=500 source vs LBL at the cold strat node.

T2 showed the shipped PF table is 14.5x too opaque in the 500-630 cm^-1
wing at the 110 K node (p=1e4 Pa), and that this is a TABLE NODE error, not
interpolation. This script asks the final question: is the error in the
*raw Chaverot source* CDF, or introduced by our bin_down + GL 20->8 rebin?

It reads the native R=500 Chaverot Ktable5d via exo_k (same loader the
generator uses), computes the uniform-nu band-mean k over the diagnostic
CO2 bands at the exact (T, p) nodes, converts m^2/molecule -> m^2/kg with
the same factor the generator applies (N_A / M_air), and compares to LBL
band-mean k from lbl_kT_curvature.npz.

  native/LBL ~ 1  -> source faithful; bias is our rebin pipeline
  native/LBL >> 1 -> bias is baked into the Chaverot source CDF

Run from the `radiation` env (needs exo_k + numpy).
"""
from __future__ import annotations

import os

import numpy as np
import exo_k

SRC = ("/Users/joymonteiro/github/chaverot/data/correlated-k_tables/"
       "Earth_376ppmCO2_R500_10-30k.ktable.SI.h5")
LBL_NPZ = os.path.join(
    os.path.dirname(__file__), "..", "..", "debug_data", "lbl_kT_curvature.npz")

N_A = 6.02214076e23
M_AIR = 0.028964            # kg/mol — generator default --mixture-molar-mass
CONV = N_A / M_AIR          # m^2/molecule -> m^2/kg (generate_picket_fence_tables.py:359)

P_NODE = 1.0e4
X_DRY = 1e-6
BANDS = {
    "500-630": (500.0, 630.0, "k_wing_500_630"),
    "630-700": (630.0, 700.0, "k_core_630_700"),
    "700-800": (700.0, 800.0, "k_wing_700_800"),
}


def main():
    kt = exo_k.Ktable5d(filename=SRC)
    wns = np.asarray(kt.wns)
    bw = np.diff(np.asarray(kt.wnedges))
    w = np.asarray(kt.weights)
    p = np.asarray(kt.pgrid)
    x = np.asarray(kt.xgrid)
    ip = int(np.argmin(np.abs(p - P_NODE)))
    ix = int(np.argmin(np.abs(x - X_DRY)))

    lbl = np.load(LBL_NPZ)
    Tl = lbl["T_grid"]

    def lblval(key, T):
        return float(lbl[key][int(np.argmin(np.abs(Tl - T)))])

    def native_mean(iT, lo, hi):
        sel = (wns >= lo) & (wns < hi)
        binmean = (kt.kdata[ip, iT, ix, :, :][sel] * w[None, :]).sum(axis=1)
        return np.average(binmean, weights=bw[sel]) * CONV

    # Chaverot T nodes: index 1 = 110 K, index 2 = 170 K.
    print(f"kdata_unit={kt.kdata_unit}  CONV={CONV:.3e} m2/kg per m2/molec")
    print(f"p node={p[ip]:.0f} Pa, x={x[ix]:.1e} (dry)\n")
    print(f"{'band':<9}{'T':>5}{'nativeR500':>13}{'LBL':>12}{'native/LBL':>12}")
    for T, iT in [(110.0, 1), (170.0, 2)]:
        for nm, (lo, hi, key) in BANDS.items():
            nv = native_mean(iT, lo, hi)
            l = lblval(key, T)
            r = nv / l if l > 0 else np.nan
            print(f"{nm:<9}{T:>5.0f}{nv:>13.3e}{l:>12.3e}{r:>12.2f}")
        print()
    print("native/LBL >> 1 in 500-630 at 110 K -> Chaverot SOURCE CDF defect")


if __name__ == "__main__":
    main()
