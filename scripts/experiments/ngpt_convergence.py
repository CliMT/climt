"""Experiment #8 — does the GL band-mean converge to LBL as ngpt grows?

Exp #7 left the corrected linepyline table reading band-mean k at 0.4-0.9x LBL
(via the shipped ngpt=8 GL quadrature), even though the per-nu kappa is exactly
1.00x LBL. The claim was: that deficit is a Gauss-Legendre quadrature artifact
(8 nodes under-sample the skewed k-CDF's high-k tail), NOT a kappa error, and it
must vanish as ngpt increases. This script tests that claim directly.

It samples CO2 kappa ONCE at the diagnostic node (T, p=1e4 Pa, dry 376 ppm) with
the validated recipe (binning=False, dnu=0.1, pseudovoigt), then runs the REAL
pipeline function kappa_to_k_coeffs at ngpt = 2,4,8,16,32,64,128 and computes the
g-integrated band-mean sum(k_i * w_i) for each. Convergence target = the plain
arithmetic mean of the sampled kappa over the band (which == the LBL band-mean,
since LBL band-mean is the same uniform-nu mean and kappa is 1.00x LBL).

Run from the linepyline env (needs linepyline + the pure-numpy cork_table_builder):
    /Users/joymonteiro/miniconda3/envs/linepyline/bin/python \
        scripts/experiments/ngpt_convergence.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import xarray as xr

import linepyline as lpl

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ))
from cork_table_builder.k_distribution import kappa_to_k_coeffs

CO2_VMR = 376e-6
P_FIXED = 1.0e4
P_PROF = np.array([1.0e3, P_FIXED, 1.0e5])
J = 1
NU_MIN, NU_MAX, DNU = 10.0, 3250.0, 0.1
NGPTS = [2, 4, 8, 16, 32, 64, 128]
BANDS = {"wing_500_630": (500.0, 630.0),
         "core_630_700": (630.0, 700.0),
         "wing_700_800": (700.0, 800.0),
         "v3_1800_3250": (1800.0, 3250.0)}

_MM = {"CO2": 44.01, "air": 28.97}
MEAN_MW = CO2_VMR * _MM["CO2"] + (1.0 - CO2_VMR) * _MM["air"]
MASS_FRAC = (CO2_VMR * _MM["CO2"]) / MEAN_MW

LBL_NPZ = os.path.join(os.path.dirname(__file__), "..", "..",
                       "debug_data", "lbl_kT_curvature.npz")
LBL_KEY = {"wing_500_630": "k_wing_500_630", "core_630_700": "k_core_630_700",
           "wing_700_800": "k_wing_700_800", "v3_1800_3250": "k_v3_1800_3250"}


def main():
    lbl = np.load(LBL_NPZ)
    Tl = lbl["T_grid"]

    def lblval(key, T):
        return float(lbl[key][int(np.argmin(np.abs(Tl - T)))])

    rtm = lpl.rtm(background_gas="air", use_numba=True)
    band_edges = np.array([10.0, 500.0, 630.0, 700.0, 800.0, 1800.0, 3250.0])

    for T in (110.0, 170.0):
        T_arr = np.full(len(P_PROF), T)
        kda = rtm.get_kappa_hitran(
            "CO2", NU_MIN, NU_MAX, DNU, P_PROF, T_arr, p_self=P_PROF * CO2_VMR,
            line_shape="pseudovoigt", binning=False, remove_plinth=False,
        )
        nu = np.asarray(kda["nu"].values)
        kappa_row = np.maximum(kda.transpose("p", "nu").values[J], 0.0) * MASS_FRAC

        # band index lookup in band_edges for the 4 diagnostic bands
        print(f"=== T = {T:.0f} K, p = {P_FIXED:.0f} Pa "
              f"(kappa binning=False dnu=0.1, * CO2 mass frac) ===")
        hdr = f"{'band':<14}{'LBL':>11}{'arith_mean':>12}" + \
              "".join(f"{('g='+str(n)):>9}" for n in NGPTS)
        print(hdr)
        for name, (lo, hi) in BANDS.items():
            sel = (nu >= lo) & (nu < hi)
            arith = float(np.mean(kappa_row[sel]))
            lblv = lblval(LBL_KEY[name], T)
            # run the real pipeline at each ngpt; kappa shaped (1, nNu)
            row = []
            for ng in NGPTS:
                kc, w = kappa_to_k_coeffs(kappa_row[None, :], nu, band_edges, ng)
                # find which band index matches (lo, hi)
                ib = int(np.argmin(np.abs(band_edges[:-1] - lo)))
                bandmean = float((kc[0, ib] * w[ib]).sum())
                row.append(bandmean)
            print(f"{name:<14}{lblv:>11.3e}{arith:>12.3e}" +
                  "".join(f"{v:>9.3e}" for v in row))
        # also show as ratio to LBL for readability
        print(f"{'  (ratio/LBL)':<14}{'1.000':>11}", end="")
        print(f"{'':>12}", end="")
        for name in BANDS:  # noqa: just spacing; ratios printed next block
            break
        print()
        for name, (lo, hi) in BANDS.items():
            sel = (nu >= lo) & (nu < hi)
            lblv = lblval(LBL_KEY[name], T)
            ratios = []
            for ng in NGPTS:
                kc, w = kappa_to_k_coeffs(kappa_row[None, :], nu, band_edges, ng)
                ib = int(np.argmin(np.abs(band_edges[:-1] - lo)))
                bandmean = float((kc[0, ib] * w[ib]).sum())
                ratios.append(bandmean / lblv if lblv > 0 else np.nan)
            print(f"{name+' /LBL':<14}{'':>11}{'':>12}" +
                  "".join(f"{r:>9.2f}" for r in ratios))
        print()

    print("Expect each band's g-integrated band-mean to climb monotonically "
          "toward arith_mean (== LBL) as ngpt grows -> confirms the 0.4-0.9x "
          "deficit at ngpt=8 is pure GL quadrature truncation, not a kappa error.")


if __name__ == "__main__":
    main()
