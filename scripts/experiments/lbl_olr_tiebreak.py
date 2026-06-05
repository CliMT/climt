"""Experiment #21 step 2 (linepyline env) — LBL OLR tiebreak: PF vs RRTMG vs LBL.

The PF-vs-RRTMG OLR gap (exp #16) could mean PF over-traps OR RRTMG under-traps.
linepyline (line-by-line, the same reference that built the PF k-table) settles it:
run the LBL on the IDENTICAL fixed profile (from dump_forward_profile.py) and see
whether the true OLR sits with PF or with RRTMG.

To make this apples-to-apples we set the linepyline two-stream diffusivity D=1.66
(PF's value; rte-rrtmgp uses 1.64) — so the ONLY difference between LBL and PF is
full spectral resolution vs PF's 7-band/8-g-point correlated-k. The OLR is also
binned onto PF's band edges to localize any per-band over-trapping.

Run: /Users/joymonteiro/miniconda3/envs/linepyline/bin/python scripts/experiments/lbl_olr_tiebreak.py
"""
import argparse
import os
import numpy as np
import xarray as xr
import linepyline as lpl

# PF winsplit band edges (cm^-1); LBL OLR is binned onto these.
BAND_EDGES = [10.0, 500.0, 630.0, 700.0, 800.0, 1250.0, 1800.0, 3250.0]
NU_MIN, NU_MAX, DNU = 10.0, 3250.0, 0.1
D_PF = 1.66  # match PF's diffusivity so only spectral resolution differs
INCLUDE_CONTINUUM = True  # overridden by --no-continuum
SUFFIX = ""               # output filename suffix


def lbl_olr(rtm, p_pa, T, Ts, q, co2):
    """Return (total_olr, per_band_olr[list], nu, olr_spec) over [NU_MIN, NU_MAX]."""
    p_da = xr.DataArray(p_pa, dims=("p",), coords={"p": p_pa})
    T_da = xr.DataArray(T, dims=("p",), coords={"p": p_pa})
    q_da = xr.DataArray(q, dims=("p",), coords={"p": p_pa})
    ds = rtm.radiative_transfer(
        NU_MIN, NU_MAX, DNU, p_da, float(p_pa[-1]), T_da, float(Ts),
        q=q_da, absorbers={"CO2": co2},
        D=D_PF, line_shape="pseudovoigt", binning=False,
        include_mtckd_continuum=INCLUDE_CONTINUUM,
    )
    nu = ds["nu"].values
    olr_spec = np.asarray(ds["olr"].values)  # W/m^2/cm^-1 at TOA
    total = float(np.sum(olr_spec) * DNU)
    bands = []
    for lo, hi in zip(BAND_EDGES[:-1], BAND_EDGES[1:]):
        sel = (nu >= lo) & (nu < hi)
        bands.append(float(np.sum(olr_spec[sel]) * DNU))
    return total, bands, nu, olr_spec


def main():
    global INCLUDE_CONTINUUM, SUFFIX
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-continuum", action="store_true", default=False,
                    help="disable H2O MT_CKD continuum in the LBL (diagnostic)")
    ap.add_argument("--co2", type=float, default=None,
                    help="CO2 in ppm; overrides the profile CO2 (376) and tags "
                         "the output file lbl_olr_spec_<tag>_co2_<ppm>ppm.npz "
                         "(for the A5 off-node CO2-interpolation check)")
    a = ap.parse_args()
    if a.no_continuum:
        INCLUDE_CONTINUUM = False
        SUFFIX = "_nocont"

    npz = os.path.join(os.path.dirname(__file__), "..", "..",
                       "debug_data", "forward_profile.npz")
    d = np.load(npz, allow_pickle=True)

    p_mid = d["p_mid_Pa"]
    T = d["T_ref"]
    q_moist = d["q_ref_moist"]
    Ts = float(d["T_surf"])
    co2 = float(d["CO2"])
    if a.co2 is not None:
        co2 = a.co2 * 1e-6
        SUFFIX += f"_co2_{int(round(a.co2))}ppm"

    # linepyline wants p ascending (TOA first, surface last).
    order = np.argsort(p_mid)
    p_s = p_mid[order]
    T_s = T[order]
    q_ms = q_moist[order]

    rtm = lpl.rtm(background_gas="air", use_numba=True)

    print(f"LBL OLR tiebreak — fixed profile, D={D_PF} (PF-matched), "
          f"CO2={int(co2*1e6)} ppm, O3=0, Ts={Ts} K\n")

    for tag, q in (("MOIST", q_ms), ("DRY", np.zeros_like(q_ms))):
        # Reference PF/RRTMG OLR are baked at the profile CO2 (376 ppm); they are
        # meaningless under --co2 override, so fall back to NaN if absent/stale.
        _k_pf, _k_rr = f"olr_pf_{tag.lower()}", f"olr_rrtmg_{tag.lower()}"
        olr_pf = float(d[_k_pf]) if (a.co2 is None and _k_pf in d.files) else float("nan")
        olr_rr = float(d[_k_rr]) if (a.co2 is None and _k_rr in d.files) else float("nan")
        total, bands, nu, olr_spec = lbl_olr(rtm, p_s, T_s, Ts, q, co2)

        print(f"===== {tag} =====")
        print(f"  {'TOTAL OLR (10-3250 cm^-1)':32s}  RRTMG={olr_rr:7.2f}   "
              f"PF={olr_pf:7.2f}   LBL={total:7.2f}")
        print(f"  {'LBL-PF':32s}  {total-olr_pf:+7.2f}     "
              f"LBL-RRTMG {total-olr_rr:+7.2f}")
        nearer = "PF" if abs(total - olr_pf) < abs(total - olr_rr) else "RRTMG"
        print(f"  -> LBL is closer to {nearer}\n")
        print(f"  {'band [cm^-1]':18s} {'PF OLR':>9s} {'LBL OLR':>9s} {'LBL-PF':>9s}")
        for (lo, hi), bolr in zip(zip(BAND_EDGES[:-1], BAND_EDGES[1:]), bands):
            print(f"  {f'{lo:.0f}-{hi:.0f}':18s} {'':>9s} {bolr:9.2f}")
        print()

        # save spectra for optional plotting
        out = os.path.join(os.path.dirname(__file__), "..", "..",
                           "debug_data", f"lbl_olr_spec_{tag.lower()}{SUFFIX}.npz")
        np.savez(out, nu=nu, olr_spec=olr_spec, total=total,
                 bands=np.array(bands), olr_pf=olr_pf, olr_rrtmg=olr_rr)


if __name__ == "__main__":
    main()
