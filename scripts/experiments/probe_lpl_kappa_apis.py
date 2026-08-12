"""Pin the ~4.85x discrepancy between the two linepyline kappa paths.

The new linepyline-built CORK table reads ~5x too opaque in EVERY band vs LBL
(even the core, which should be ~1.0). Two code paths produce kappa:

  * LBL reference (lbl_kT_curvature.py): rtm.radiative_transfer(...)["kappa"]
    with line_shape="pseudovoigt"  -> core@110 ~ 0.110 m^2/kg (matches Chaverot)
  * table-builder (kappa_sampling.sample_kappa_grid): rtm.get_kappa_hitran(...)
    * mass_fraction, default line_shape="lorentz"  -> core@110 ~ 0.533 (4.85x)

This script computes band-mean kappa at ONE node (T=110 & 170 K, p=1e4 Pa,
dry CO2 376 ppm) several ways so we can attribute the discrepancy to:
  (a) line_shape (lorentz vs pseudovoigt),
  (b) the mass_fraction multiplier (units of get_kappa_hitran), and
  (c) anything else radiative_transfer does that get_kappa_hitran does not.

Run from the linepyline env:
    /Users/joymonteiro/miniconda3/envs/linepyline/bin/python \
        scripts/experiments/probe_lpl_kappa_apis.py
"""
from __future__ import annotations

import numpy as np
import xarray as xr

import linepyline as lpl

CO2_VMR = 376e-6
P_FIXED = 1.0e4
P_PROF = np.array([1.0e3, P_FIXED, 1.0e5])
J = 1                      # middle row == P_FIXED
NU_MIN, NU_MAX, DNU = 10.0, 3250.0, 0.1
BANDS = {"wing_500_630": (500.0, 630.0),
         "core_630_700": (630.0, 700.0),
         "wing_700_800": (700.0, 800.0),
         "v3_1800_3250": (1800.0, 3250.0)}

# mass fraction of CO2 in dry air, exactly as kappa_sampling._mass_fraction
_MM = {"CO2": 44.01, "air": 28.97}
MEAN_MW = CO2_VMR * _MM["CO2"] + (1.0 - CO2_VMR) * _MM["air"]
MASS_FRAC = (CO2_VMR * _MM["CO2"]) / MEAN_MW


def bandmean(nu, row, lo, hi):
    sel = (nu >= lo) & (nu < hi)
    return float(np.mean(row[sel]))


def main():
    print(f"CO2 mass fraction in dry air = {MASS_FRAC:.4e}  "
          f"(mean MW = {MEAN_MW:.3f} g/mol)\n")
    rtm = lpl.rtm(background_gas="air", use_numba=True)

    for T in (110.0, 170.0):
        p_da = xr.DataArray(P_PROF, dims=("p",), coords={"p": P_PROF})
        T_da = xr.DataArray(np.full_like(P_PROF, T), dims=("p",),
                            coords={"p": P_PROF})

        # --- (A) radiative_transfer kappa (LBL reference) -----------------
        ds = rtm.radiative_transfer(
            NU_MIN, NU_MAX, DNU, p_da, float(P_PROF[-1]), T_da, float(T),
            absorbers={"CO2": CO2_VMR}, line_shape="pseudovoigt",
        )
        k_rt = ds["kappa"].values
        if ds["kappa"].dims[0] != "p":
            k_rt = k_rt.T
        good = np.isfinite(k_rt).all(axis=1)
        k_rt_row = k_rt[good][J]
        nu_rt = ds["nu"].values

        # --- (B,C,D) get_kappa_hitran: binning True/False, dnu 0.1/0.5 -----
        T_arr = np.full(len(P_PROF), T)
        p_self = P_PROF * CO2_VMR
        variants = {
            "bin_dnu0.1": dict(binning=True, dnu=0.1),
            "nobin_dnu0.1": dict(binning=False, dnu=0.1),
            "bin_dnu0.5": dict(binning=True, dnu=0.5),
            "nobin_dnu0.5": dict(binning=False, dnu=0.5),
        }
        rows = {}
        for tag, kw in variants.items():
            kda = rtm.get_kappa_hitran(
                "CO2", NU_MIN, NU_MAX, kw["dnu"], P_PROF, T_arr, p_self=p_self,
                line_shape="pseudovoigt", binning=kw["binning"],
                remove_plinth=False,
            )
            arr = kda.transpose("p", "nu").values
            rows[tag] = (np.asarray(kda["nu"].values), np.maximum(arr[J], 0.0))

        print(f"=== T = {T:.0f} K, p = {P_FIXED:.0f} Pa "
              f"(RT=radiative_transfer pvoigt, *MF applied to hitran) ===")
        print(f"{'band':<14} {'RT':>11}", end="")
        for tag in variants:
            print(f" {tag+'*MF':>16}", end="")
        print(f" {'ratios vs RT':>28}")
        for name, (lo, hi) in BANDS.items():
            rt = bandmean(nu_rt, k_rt_row, lo, hi)
            vals = {tag: bandmean(rows[tag][0], rows[tag][1], lo, hi) * MASS_FRAC
                    for tag in variants}
            print(f"{name:<14} {rt:>11.3e}", end="")
            for tag in variants:
                print(f" {vals[tag]:>16.3e}", end="")
            ratios = "  ".join(
                f"{tag.split('_')[0]}{tag.split('dnu')[1]}={vals[tag]/rt:.2f}"
                if rt > 0 else "nan" for tag in variants)
            print(f"   {ratios}")
        print()

    print("Reading: if hitran_lor*mf / RT ~ 1, the builder is consistent and "
          "the table opacity bug is elsewhere.\n"
          "If RT ~ hitran*mf only WITHOUT the mass_frac multiplier, then "
          "get_kappa_hitran already returns m^2/kg-mixture and the builder's "
          "extra *mass_frac is the (T-independent) bug.")


if __name__ == "__main__":
    main()
