"""Experiment #13 — PF table H2O-VMR-axis fidelity vs LBL (lines + MT_CKD).

All prior k-fidelity work (#1-#11) validated the CO2 bands at the DRY node
(X_H2O=1e-6) only. This validates the water-vapour axis the moist RCE depends
on: it compares the PF table band-mean k against the linepyline line-by-line
band-mean kappa across the full H2O VMR axis, at tropospheric (T, p) nodes.

LBL reference = `pf_table_builder.kappa_sampling.sample_kappa_grid` with the
validated recipe (binning=False, dnu=0.1, pseudovoigt) — the exact kappa the
earth linepyline table was built from (exp #7). We sample it TWICE, with the
MT_CKD continuum ON and OFF, so we can:
  * confirm the table tracks the continuum-ON reference (continuum is in the
    table, as the code path in generate_pf_tables_linepyline implies), and
  * quantify how much the continuum adds per band/X (the window 800-1800 is the
    region it matters for moist columns).

This is a SOURCE-faithful node-fidelity check (kappa recipe was validated vs
linepyline.radiative_transfer for CO2 in exp #7). A fully independent
radiative_transfer cross-check of the H2O+continuum kappa remains a follow-up.

Run (linepyline env, has linepyline + MT_CKD + netCDF4):
  /Users/joymonteiro/miniconda3/envs/linepyline/bin/python \
      scripts/experiments/pf_h2o_axis_fidelity.py [--dnu 0.1] [--T 230 290] [--p 1e4 1e5]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import netCDF4 as nc

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.join(_HERE, "..", "..")
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_table_builder.kappa_sampling import sample_kappa_grid  # noqa: E402

_CK = os.path.join(_REPO, "climt", "_data", "picket_fence", "correlated_k")

# earth scenario composition (matches generate_pf_tables_linepyline.py SCENARIOS["earth"])
BACKGROUND = "air"
ABSORBERS = {"CO2": 376e-6}
H2O_VMR_GRID = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0])
NU_LO, NU_HI = 10.0, 3250.0

TABLES = {
    "linepyline 6band": "earth_low_res_lw_co2refined_linepyline",
    "shipped 4band":    "earth_low_res_lw_8gpt",
}


def pf_band_means(table_path, T_tgt, p_tgt_pa, X_targets):
    """PF table band-mean k (m^2/kg) at exact (T, p, X) nodes -> dict[band_edge] = array over X."""
    with nc.Dataset(table_path) as ds:
        kc = ds["k_coefficients"][...]          # (1, nband, ngpt, nT, nP, nX)
        gw = ds["gpoint_weights"][...]           # (nband, ngpt)
        be = ds["band_wavenumber_limits"][...]
        Tg = np.asarray(ds["temperature_grid"][...])
        pg = np.exp(np.asarray(ds["pressure_grid_log"][...]))   # ln(Pa) -> Pa
        xg = np.asarray(ds["h2o_vmr_grid"][...])
    iT = int(np.argmin(np.abs(Tg - T_tgt)))
    iP = int(np.argmin(np.abs(pg - p_tgt_pa)))
    edges = [(float(be[b, 0]), float(be[b, 1])) for b in range(be.shape[0])]
    out = {}
    for b, edge in enumerate(edges):
        vals = []
        for X in X_targets:
            iX = int(np.argmin(np.abs(xg - X)))
            vals.append(float(np.sum(kc[0, b, :, iT, iP, iX] * gw[b]) / np.sum(gw[b])))
        out[edge] = np.array(vals)
    return out, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dnu", type=float, default=0.1)
    ap.add_argument("--T", type=float, nargs="+", default=[230.0, 290.0])
    ap.add_argument("--p", type=float, nargs="+", default=[1e4, 1e5])
    args = ap.parse_args()

    T_grid = np.array(args.T)
    p_grid = np.array(args.p)
    nu_grid = np.arange(NU_LO, NU_HI + args.dnu / 2, args.dnu)
    print(f"Sampling LBL kappa: nT={len(T_grid)} nP={len(p_grid)} "
          f"nX={len(H2O_VMR_GRID)} nNu={len(nu_grid)} (dnu={args.dnu})  [continuum ON & OFF]")

    common = dict(background_gas=BACKGROUND, absorbers=ABSORBERS,
                  h2o_vmr_grid=H2O_VMR_GRID, T_grid=T_grid, p_grid=p_grid,
                  nu_grid=nu_grid, line_shape="pseudovoigt", binning=False)
    kappa_on = sample_kappa_grid(**common, include_mtckd_continuum=True)
    kappa_off = sample_kappa_grid(**common, include_mtckd_continuum=False)
    # shapes: (nT, nP, nX, nNu)

    def band_mean_lbl(kappa, edge):
        m = (nu_grid >= edge[0]) & (nu_grid < edge[1])
        return kappa[:, :, :, m].mean(axis=3)   # (nT, nP, nX)

    for tname, base in TABLES.items():
        path = os.path.join(_CK, base + ".nc")
        print(f"\n{'='*78}\n{tname}\n{'='*78}")
        for iT, Tv in enumerate(T_grid):
            for iP, pv in enumerate(p_grid):
                pfm, edges = pf_band_means(path, Tv, pv, H2O_VMR_GRID)
                print(f"\n-- T={Tv:.0f} K, p={pv/100:.0f} hPa --")
                print(f"{'band':12s} {'X_H2O':>8s} {'PF k':>10s} {'LBL_on':>10s} "
                      f"{'PF/LBL':>7s} {'cont x':>7s}")
                for edge in edges:
                    lon = band_mean_lbl(kappa_on, edge)[iT, iP]    # (nX,)
                    loff = band_mean_lbl(kappa_off, edge)[iT, iP]
                    pfk = pfm[edge]
                    for ix, X in enumerate(H2O_VMR_GRID):
                        r = pfk[ix] / lon[ix] if lon[ix] > 0 else np.inf
                        contx = lon[ix] / loff[ix] if loff[ix] > 0 else np.inf
                        print(f"{edge[0]:.0f}-{edge[1]:.0f}".ljust(12) +
                              f" {X:8.0e} {pfk[ix]:10.3e} {lon[ix]:10.3e} "
                              f"{r:7.2f} {contx:7.2f}")


if __name__ == "__main__":
    main()
