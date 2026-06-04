"""Experiment #23 (climt env) — fast band-structure evaluator vs LBL.

Inner-loop tool for designing a PF LW band structure that captures single-column
OLR. Computes PF per-band TOA OLR on the FIXED exp #16 profile (from
forward_profile.npz) and compares to the line-by-line truth by re-binning the
saved LBL spectrum (lbl_olr_spec_{moist,dry}.npz, computed once at D=1.66) onto
the candidate table's own band edges. No linepyline re-run per iteration.

GOAL METRIC: |LBL - PF| total OLR, moist & dry. Target < 5 W/m^2 (= RRTMG's own
agreement with LBL: 4.1 moist / 7.3 dry). Per-band LBL-PF localizes residual
over-trapping.

Usage: /Users/joymonteiro/miniconda3/envs/climt/bin/python \
         scripts/experiments/eval_band_structure.py <table_name>
"""
import os
import sys
import numpy as np
import xarray as xr

from climt import get_default_state, get_grid
from climt._components.picket_fence import PicketFenceLongwave

NZ = 40
DNU = 0.1
_HERE = os.path.dirname(__file__)
_DBG = os.path.join(_HERE, "..", "..", "debug_data")
CKDIR = os.path.join("climt", "_data", "picket_fence", "correlated_k")

# LBL totals (exp #21, fixed profile, D=1.66). RRTMG totals for reference.
LBL_TOTAL = {"moist": 246.76, "dry": 347.89}
RRTMG_TOTAL = {"moist": 242.62, "dry": 340.62}


def pf_olr(table, T, q, Ts, co2, grid):
    pf = PicketFenceLongwave(optics="correlated_k", table=table)
    state = get_default_state([pf], grid_state=grid)
    state["air_temperature"].values[:, 0, 0] = T
    state["specific_humidity"].values[:, 0, 0] = q
    state["surface_temperature"].values[:] = Ts
    if "surface_longwave_emissivity" in state:
        state["surface_longwave_emissivity"].values[:] = 1.0
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2
    if "mole_fraction_of_ozone_in_air" in state:
        state["mole_fraction_of_ozone_in_air"].values[:] = 0.0
    _t, diag = pf(state)
    total = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    per_band = diag["upwelling_longwave_flux_in_air_per_band"].values[-1, 0, 0, :]
    return total, np.asarray(per_band)


def main():
    table = sys.argv[1]
    lbl_suffix = ""
    co2_override = None
    for a in sys.argv[2:]:
        if a.startswith("co2="):
            co2_override = float(a.split("=", 1)[1])
        elif not a.startswith("co2="):
            lbl_suffix = a
    d = np.load(os.path.join(_DBG, "forward_profile.npz"), allow_pickle=True)
    p_mid, T, q_moist, Ts, co2 = (d["p_mid_Pa"], d["T_ref"],
                                  d["q_ref_moist"], float(d["T_surf"]),
                                  float(d["CO2"]))
    if co2_override is not None:
        co2 = co2_override
    grid = get_grid(nx=1, ny=1, nz=NZ)

    with xr.open_dataset(os.path.join(CKDIR, table + ".nc")) as ds:
        edges = ds["band_wavenumber_limits"].values  # (nband, 2)

    print(f"=== {table}  ({len(edges)} bands) ===")
    print(f"{'case':6s} {'RRTMG':>8s} {'PF':>8s} {'LBL':>8s} "
          f"{'LBL-PF':>8s} {'LBL-RRTMG':>10s}")
    for tag, q in (("moist", q_moist), ("dry", np.zeros_like(q_moist))):
        total_pf, per_band_pf = pf_olr(table, T, q, Ts, co2, grid)
        lbl = np.load(os.path.join(_DBG, f"lbl_olr_spec_{tag}{lbl_suffix}.npz"))
        nu, olr_spec = lbl["nu"], lbl["olr_spec"]
        lbl_total = float(np.sum(olr_spec) * DNU)
        rr = RRTMG_TOTAL[tag] if lbl_suffix == "" else float("nan")
        print(f"{tag:6s} {rr:8.2f} {total_pf:8.2f} "
              f"{lbl_total:8.2f} {lbl_total-total_pf:+8.2f} "
              f"{lbl_total-rr:+10.2f}")
        # per-band LBL-PF
        print(f"   {'band [cm^-1]':16s} {'PF':>8s} {'LBL':>8s} {'LBL-PF':>8s}")
        for (lo, hi), pf_b in zip(edges, per_band_pf):
            sel = (nu >= lo) & (nu < hi)
            lbl_b = float(np.sum(olr_spec[sel]) * DNU)
            flag = "  <<" if abs(lbl_b - pf_b) > 4.0 else ""
            print(f"   {f'{lo:.0f}-{hi:.0f}':16s} {pf_b:8.2f} {lbl_b:8.2f} "
                  f"{lbl_b - pf_b:+8.2f}{flag}")
        print()


if __name__ == "__main__":
    main()
