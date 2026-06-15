"""Experiment T1 — LBL k(T) curvature check at a fixed strat pressure.

Runs linepyline LBL at p≈100 hPa (1e4 Pa) across a fine cold-T grid and
computes uniform-nu band-mean kappa(T) over the diagnostic CO2 bands. Then,
treating the Chaverot source T-grid (which brackets the strat with the
110 K and 170 K nodes), compares two interpolation estimates of k(147 K)
against the LBL-true k(147 K):

  * linear-in-k  : what the climt loader currently does
  * log-in-k     : a candidate cheap fix (loader-only, no table regen)

If linear-in-k overestimates k(147) substantially but log-k is close to
truth, the cold-strat bias is (at least partly) an interpolation artifact
of the coarse Chaverot T-grid, fixable in the loader.

Usage (linepyline env):
    /Users/joymonteiro/miniconda3/envs/linepyline/bin/python \
        scripts/experiments/lbl_kT_curvature.py
"""
from __future__ import annotations

import os

import numpy as np
import xarray as xr

import linepyline as lpl

CO2_VMR = 376e-6
P_FIXED = 1.0e4  # Pa  (~100 hPa, the strat-bias hotspot layer)

# Fine temperature grid: brackets the Chaverot 110/170 K nodes and samples
# the strat target (147 K) plus colder layers where the bias explodes.
T_GRID = np.array([94.0, 110.0, 123.0, 147.0, 170.0, 184.0, 230.0])

# Chaverot source nodes that bracket the strat target T=147 K.
T_LO, T_HI, T_TARGET = 110.0, 170.0, 147.0

# Diagnostic CO2 bands (cm^-1): the two 15-um wings, the core, and v3.
BANDS = {
    "wing_500_630": (500.0, 630.0),
    "core_630_700": (630.0, 700.0),
    "wing_700_800": (700.0, 800.0),
    "v3_1800_3250": (1800.0, 3250.0),
}

NU_MIN = 10.0
NU_MAX = 3250.0
DNU = 0.1


def band_mean_kappa(nu, kappa_row, lo, hi):
    """Uniform-nu mean of kappa over [lo, hi)."""
    sel = (nu >= lo) & (nu < hi)
    return float(np.mean(kappa_row[sel]))


def main():
    rtm = lpl.rtm(background_gas="air", use_numba=True)

    # linepyline requires a multi-level profile (len(p) > 1). Use a 3-level
    # isothermal column bracketing the target pressure and extract the
    # p == P_FIXED row so pressure broadening is correct at that level.
    p_prof = np.array([1.0e3, P_FIXED, 1.0e5])
    j_fixed = 1

    nu_ref = None
    # kband[band][i] = band-mean kappa at T_GRID[i]
    kband = {name: np.zeros(len(T_GRID)) for name in BANDS}

    for i, T in enumerate(T_GRID):
        p_da = xr.DataArray(p_prof, dims=("p",), coords={"p": p_prof})
        T_da = xr.DataArray(np.full_like(p_prof, T), dims=("p",),
                            coords={"p": p_prof})
        ds = rtm.radiative_transfer(
            NU_MIN, NU_MAX, DNU, p_da, float(p_prof[-1]), T_da, float(T),
            absorbers={"CO2": CO2_VMR}, line_shape="pseudovoigt",
        )
        kappa = ds["kappa"].values
        if ds["kappa"].dims[0] != "p":
            kappa = kappa.T
        good = np.isfinite(kappa).all(axis=1)
        kappa_row = kappa[good][j_fixed]
        nu = ds["nu"].values
        if nu_ref is None:
            nu_ref = nu
        for name, (lo, hi) in BANDS.items():
            kband[name][i] = band_mean_kappa(nu, kappa_row, lo, hi)

    iT = {float(t): j for j, t in enumerate(T_GRID)}
    j_lo, j_hi, j_tgt = iT[T_LO], iT[T_HI], iT[T_TARGET]
    w = (T_TARGET - T_LO) / (T_HI - T_LO)  # linear weight toward T_HI

    print(f"p = {P_FIXED:.0f} Pa (~{P_FIXED/100:.0f} hPa), CO2 = {CO2_VMR*1e6:.0f} ppm")
    print(f"interp target T={T_TARGET} K between nodes {T_LO}/{T_HI} K (w={w:.3f})\n")
    print(f"{'band':<14} {'k_true':>10} {'lin_interp':>11} {'log_interp':>11} "
          f"{'lin/true':>9} {'log/true':>9}")
    rows = []
    for name in BANDS:
        k = kband[name]
        k_true = k[j_tgt]
        k_lo, k_hi = k[j_lo], k[j_hi]
        lin = (1 - w) * k_lo + w * k_hi
        # log-in-k (guard against zeros)
        if k_lo > 0 and k_hi > 0:
            log = np.exp((1 - w) * np.log(k_lo) + w * np.log(k_hi))
        else:
            log = lin
        r_lin = lin / k_true if k_true > 0 else np.nan
        r_log = log / k_true if k_true > 0 else np.nan
        print(f"{name:<14} {k_true:>10.3e} {lin:>11.3e} {log:>11.3e} "
              f"{r_lin:>9.3f} {r_log:>9.3f}")
        rows.append((name, k_true, lin, log, r_lin, r_log))

    print("\nfull k(T) table (band-mean kappa, m^2/kg or model units):")
    header = "T[K]    " + "  ".join(f"{n:>12}" for n in BANDS)
    print(header)
    for i, T in enumerate(T_GRID):
        vals = "  ".join(f"{kband[n][i]:>12.4e}" for n in BANDS)
        print(f"{T:6.1f}  {vals}")

    # Save raw arrays for downstream plotting / re-analysis.
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "debug_data",
    )
    os.makedirs(out_dir, exist_ok=True)
    np.savez(
        os.path.join(out_dir, "lbl_kT_curvature.npz"),
        T_grid=T_GRID,
        p_fixed=P_FIXED,
        **{f"k_{n}": kband[n] for n in BANDS},
    )
    print(f"\nwrote {os.path.join(out_dir, 'lbl_kT_curvature.npz')}")


if __name__ == "__main__":
    main()
