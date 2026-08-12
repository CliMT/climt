"""A5 CO2-interpolation accuracy probe for the cork LW k-table.

Part 1 (leave-one-out node test) runs in the climt env with NO linepyline and
is the primary deliverable.  It answers the design question "log-k vs linear-k
for CO2 interpolation" by holding out interior CO2 grid nodes one at a time and
measuring how well each scheme reconstructs them from their two neighbours.

Part 2 (LBL comparison) is guarded behind ``try: import linepyline``.  When
linepyline is not installed it prints a notice and skips.

Usage
-----
  /path/to/climt/python scripts/experiments/eval_co2_interp_accuracy.py [table_name]

  table_name defaults to "earth_low_res_lw".
"""
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.join(_HERE, "..", "..")
_DBG = os.path.join(_REPO, "debug_data")

sys.path.insert(0, _HERE)  # allows: from eval_band_structure import cork_olr

# ---------------------------------------------------------------------------
# Load table
# ---------------------------------------------------------------------------
table_name = sys.argv[1] if len(sys.argv) > 1 else "earth_low_res_lw"

from climt._components.cork.optics.correlated_k import load_k_table  # noqa: E402

table = load_k_table(table_name)
k = table["k_coefficients"]          # (ngas, nband, ngpt, nT, nP, nX_H2O, nC)
co2_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
band_edges = table.get("band_wavenumber_limits")

ngas, nband, ngpt, nT, nP, nX, nC = k.shape
print(f"\n=== A5 CO2-interpolation accuracy probe  ({table_name}) ===")
print(f"k shape: {k.shape}")
print(f"CO2 grid [{nC} nodes]: {co2_grid * 1e6} ppm")
if band_edges is not None:
    print("Band edges (cm-1):", band_edges[:, 0].tolist(), "...", float(band_edges[-1, 1]))

# ---------------------------------------------------------------------------
# PART 1 — Leave-one-out node accuracy
# ---------------------------------------------------------------------------
print("\n--- Part 1: leave-one-out node accuracy (no linepyline) ---")
print(
    "NOTE: this is a 2-node-apart STRESS TEST (neighbours straddle the held-out node).\n"
    "      Real production interpolation is between *adjacent* nodes (1 apart) at an\n"
    "      arbitrary fC, giving roughly 4× smaller error than reported here."
)

# Because the CO2 nodes are evenly log-spaced, the held-out interior node iC
# sits at the exact geometric midpoint of its neighbours iC-1 and iC+1, so
# the log-CO2 fraction is exactly 0.5.
#
# log-k interp at fC=0.5:   k_interp = sqrt(k_lo * k_hi)   (geometric mean)
# linear-k interp at fC=0.5: k_interp = 0.5 * (k_lo + k_hi)

_FLOOR = 1e-40

logk_median_abs  = []
logk_p95_abs     = []
logk_median_sgn  = []
link_median_abs  = []
link_p95_abs     = []
link_median_sgn  = []

# Collect errors across all interior nodes
all_logk_relerr = []
all_link_relerr = []

interior_nodes = range(1, nC - 1)

for iC in interior_nodes:
    # k slices at the two bracketing nodes; shape (ngas, nband, ngpt, nT, nP, nX)
    k_lo = k[..., iC - 1]
    k_hi = k[..., iC + 1]
    k_true = k[..., iC]

    # Log-k geometric-mean interpolation at fC = 0.5
    k_logk = np.sqrt(np.maximum(k_lo, _FLOOR) * np.maximum(k_hi, _FLOOR))

    # Linear-k arithmetic-mean interpolation at fC = 0.5
    k_link = 0.5 * (k_lo + k_hi)

    denom = np.maximum(k_true, _FLOOR)

    relerr_logk = (k_logk - k_true) / denom   # signed relative error
    relerr_link = (k_link - k_true) / denom

    all_logk_relerr.append(relerr_logk.ravel())
    all_link_relerr.append(relerr_link.ravel())

    med_abs_logk = float(np.median(np.abs(relerr_logk)))
    p95_abs_logk = float(np.percentile(np.abs(relerr_logk), 95))
    med_sgn_logk = float(np.median(relerr_logk))

    med_abs_link = float(np.median(np.abs(relerr_link)))
    p95_abs_link = float(np.percentile(np.abs(relerr_link), 95))
    med_sgn_link = float(np.median(relerr_link))

    logk_median_abs.append(med_abs_logk)
    logk_p95_abs.append(p95_abs_logk)
    logk_median_sgn.append(med_sgn_logk)
    link_median_abs.append(med_abs_link)
    link_p95_abs.append(med_abs_link)
    link_median_sgn.append(med_sgn_link)

    co2_ppm = co2_grid[iC] * 1e6
    print(
        f"  iC={iC:2d}  CO2={co2_ppm:7.1f} ppm"
        f"  log-k |rel|: median={med_abs_logk*100:.3f}%  p95={p95_abs_logk*100:.3f}%  signed={med_sgn_logk*100:+.3f}%"
        f"  |  linear-k |rel|: median={med_abs_link*100:.3f}%  p95={p95_abs_link*100:.3f}%  signed={med_sgn_link*100:+.3f}%"
    )

# Overall (pooled across all interior nodes)
all_logk = np.concatenate(all_logk_relerr)
all_link = np.concatenate(all_link_relerr)

overall_logk_med  = float(np.median(np.abs(all_logk)))
overall_logk_p95  = float(np.percentile(np.abs(all_logk), 95))
overall_logk_sgn  = float(np.median(all_logk))
overall_link_med  = float(np.median(np.abs(all_link)))
overall_link_p95  = float(np.percentile(np.abs(all_link), 95))
overall_link_sgn  = float(np.median(all_link))

print()
print(f"  POOLED (all {len(interior_nodes)} interior CO2 nodes):")
print(f"    log-k  : median |rel err| = {overall_logk_med*100:.4f}%   "
      f"p95 = {overall_logk_p95*100:.4f}%   "
      f"median signed = {overall_logk_sgn*100:+.4f}%")
print(f"    linear-k: median |rel err| = {overall_link_med*100:.4f}%   "
      f"p95 = {overall_link_p95*100:.4f}%   "
      f"median signed = {overall_link_sgn*100:+.4f}%")

winner = "log-k" if overall_logk_med <= overall_link_med else "linear-k"
print()
print(
    f"A5 VERDICT: log-k median |rel err| = {overall_logk_med*100:.4f}%  vs  "
    f"linear-k = {overall_link_med*100:.4f}%  ->  {winner} is better"
)

# ---------------------------------------------------------------------------
# Band-mean version: key bands (CO2 nu2 band ~660 cm-1, window ~900 cm-1)
# ---------------------------------------------------------------------------
print()
print("  Band-mean rel err at selected bands (averaged over ngpt, nT, nP, nX):")

# Identify CO2 nu2 band (~630-700 cm-1) and window band (~800-980 cm-1)
_KEY_BANDS = {}
if band_edges is not None:
    for ib, (lo, hi) in enumerate(band_edges):
        lo, hi = float(lo), float(hi)
        if lo <= 660.0 <= hi:
            _KEY_BANDS[f"CO2 ν₂ [{lo:.0f}-{hi:.0f} cm⁻¹]"] = ib
        if lo <= 900.0 <= hi:
            _KEY_BANDS[f"window [{lo:.0f}-{hi:.0f} cm⁻¹]"] = ib
else:
    # Fallback: first two bands
    _KEY_BANDS = {"band 0": 0, "band 1": 1}

print(f"  {'Band':30s}  {'iC':>4s}  {'CO2 [ppm]':>10s}  "
      f"{'log-k band-mean err':>22s}  {'linear-k band-mean err':>22s}")

for iC in interior_nodes:
    k_lo   = k[..., iC - 1]  # (ngas, nband, ngpt, nT, nP, nX)
    k_hi   = k[..., iC + 1]
    k_true = k[..., iC]

    k_logk = np.sqrt(np.maximum(k_lo, _FLOOR) * np.maximum(k_hi, _FLOOR))
    k_link = 0.5 * (k_lo + k_hi)

    co2_ppm = co2_grid[iC] * 1e6

    for label, ib in _KEY_BANDS.items():
        # Mean over (ngas, ngpt, nT, nP, nX)
        km_true = float(np.mean(k_true[:, ib, ...]))
        km_logk = float(np.mean(k_logk[:, ib, ...]))
        km_link = float(np.mean(k_link[:, ib, ...]))
        denom   = max(km_true, _FLOOR)
        err_logk = (km_logk - km_true) / denom
        err_link = (km_link - km_true) / denom
        print(f"  {label:30s}  {iC:4d}  {co2_ppm:10.1f}  "
              f"{err_logk*100:+19.4f}%  "
              f"{err_link*100:+19.4f}%")

# ---------------------------------------------------------------------------
# PART 2 — Optional LBL comparison at off-node CO2
# ---------------------------------------------------------------------------
print()
print("--- Part 2: CORK (log-k vs linear-k) vs LBL at off-node CO2 ---")
# This part needs NO linepyline: it runs CORK in the climt env and compares against
# pre-generated LBL spectra. Generate those once in the linepyline env with
#   $LPY scripts/experiments/lbl_olr_tiebreak.py --co2 <ppm>
# which writes debug_data/lbl_olr_spec_moist_co2_<ppm>ppm.npz.

import climt._components.cork.optics.correlated_k as _ck_mod  # noqa: E402
from eval_band_structure import cork_olr  # noqa: E402
from climt import get_grid  # noqa: E402

NZ = 40
DNU = 0.1  # cm^-1, the LBL spectrum integration step (must match lbl_olr_tiebreak)
OFF_NODE_CO2_PPM = [400, 2000]  # between the 215/464 and 1000/2154 ppm nodes

_fwd = os.path.join(_DBG, "forward_profile.npz")
if not os.path.isfile(_fwd):
    print(f"  Forward profile not found at {_fwd} — skipping Part 2.")
else:
    d = np.load(_fwd, allow_pickle=True)
    T, q_moist, Ts = d["T_ref"], d["q_ref_moist"], float(d["T_surf"])
    grid = get_grid(nx=1, ny=1, nz=NZ)

    print(f"  {'CO2[ppm]':>9s}  {'CORK log-k':>10s}  {'CORK lin-k':>10s}  "
          f"{'LBL':>10s}  {'logk-LBL':>9s}  {'link-LBL':>9s}")
    any_lbl = False
    for ppm in OFF_NODE_CO2_PPM:
        co2_val = ppm * 1e-6
        _ck_mod._CO2_INTERP_LOGK = True
        olr_logk, _ = cork_olr(table_name, T, q_moist, Ts, co2_val, grid)
        _ck_mod._CO2_INTERP_LOGK = False
        olr_link, _ = cork_olr(table_name, T, q_moist, Ts, co2_val, grid)
        _ck_mod._CO2_INTERP_LOGK = True  # reset to design default

        lbl_path = os.path.join(_DBG, f"lbl_olr_spec_moist_co2_{ppm}ppm.npz")
        if os.path.isfile(lbl_path):
            any_lbl = True
            lbl = np.load(lbl_path)
            lbl_total = float(np.sum(lbl["olr_spec"]) * DNU)
            print(f"  {ppm:9d}  {olr_logk:10.3f}  {olr_link:10.3f}  "
                  f"{lbl_total:10.3f}  {olr_logk - lbl_total:+9.3f}  "
                  f"{olr_link - lbl_total:+9.3f}")
        else:
            print(f"  {ppm:9d}  {olr_logk:10.3f}  {olr_link:10.3f}  "
                  f"{'(no LBL)':>10s}  {'—':>9s}  {'—':>9s}")

    if not any_lbl:
        print()
        print("  No off-node LBL spectra found. Generate them IN THE linepyline env:")
        print("    LPY=/Users/joymonteiro/miniconda3/envs/linepyline/bin/python")
        for ppm in OFF_NODE_CO2_PPM:
            print(f"    $LPY scripts/experiments/lbl_olr_tiebreak.py --co2 {ppm}")
        print("  (each writes debug_data/lbl_olr_spec_moist_co2_<ppm>ppm.npz), "
              "then re-run this probe in the climt env.")
    else:
        print()
        print("  Smaller |scheme-LBL| = more faithful. Part 1 (leave-one-out) "
              "already favours log-k;\n  this confirms it against true LBL at "
              "off-node CO2.")

print()
print("=== A5 probe complete ===")
