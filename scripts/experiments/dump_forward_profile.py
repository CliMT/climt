"""Experiment #21 step 1 (climt env) — dump the exp #16 fixed profile + PF/RRTMG OLR.

Reproduces the identical (T, p, q) profile used by lw_forward_pf_vs_rrtmg.py and
records RRTMG and PF (winsplit table) OLR for both moist and dry cases, then saves
everything to debug_data/forward_profile.npz so the linepyline-env tiebreak script
(lbl_olr_tiebreak.py) can run the LBL on the SAME profile.

Run: /Users/joymonteiro/miniconda3/envs/climt/bin/python scripts/experiments/dump_forward_profile.py
"""
import os
import re
import subprocess
import sys
import numpy as np

from climt import RRTMGLongwave, get_default_state, get_grid

CO2 = 376e-6
NZ = 40
T_SURF = 288.0
PF_TABLE = "earth_low_res_lw_co2refined_linepyline_winsplit"
_HERE = os.path.dirname(__file__)
_FWD = os.path.join(_HERE, "lw_forward_pf_vs_rrtmg.py")


def reference_profile(p_hpa):
    p_s = p_hpa.max()
    T = np.maximum(200.0, 288.0 * (p_hpa / p_s) ** 0.18)
    q = 0.015 * (p_hpa / p_s) ** 3.0
    return T, q


def forward_olr(dry):
    """Run lw_forward_pf_vs_rrtmg.py in its own process (one RRTMG per process —
    a second RRTMG instance in the same process returns garbage) and parse the
    'OLR (TOA up)  RRTMG  PF  PF-RRTMG' line."""
    cmd = [sys.executable, _FWD, "--table", PF_TABLE]
    if dry:
        cmd.append("--dry")
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        if line.strip().startswith("OLR (TOA up)"):
            nums = re.findall(r"-?\d+\.\d+", line)
            return float(nums[0]), float(nums[1])  # RRTMG, PF
    raise RuntimeError("could not parse OLR from forward script output")


def main():
    grid = get_grid(nx=1, ny=1, nz=NZ)
    rr = RRTMGLongwave()
    rr_state = get_default_state([rr], grid_state=grid)
    p_hpa = rr_state["air_pressure"].values[:, 0, 0] / 100.0
    p_int_hpa = rr_state["air_pressure_on_interface_levels"].values[:, 0, 0] / 100.0
    T_ref, q_ref_moist = reference_profile(p_hpa)

    results = {}
    for tag, dry in (("moist", False), ("dry", True)):
        olr_rr, olr_pf = forward_olr(dry)
        results[f"olr_rrtmg_{tag}"] = olr_rr
        results[f"olr_pf_{tag}"] = olr_pf
        print(f"[{tag}] RRTMG OLR = {olr_rr:.2f}   PF OLR = {olr_pf:.2f}   "
              f"PF-RRTMG = {olr_pf - olr_rr:+.2f}")

    out = os.path.abspath(os.path.join(_HERE, "..", "..",
                                       "debug_data", "forward_profile.npz"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(
        out,
        p_mid_Pa=p_hpa * 100.0,
        p_int_Pa=p_int_hpa * 100.0,
        T_ref=T_ref,
        q_ref_moist=q_ref_moist,
        T_surf=T_SURF,
        CO2=CO2,
        pf_table=PF_TABLE,
        **results,
    )
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
