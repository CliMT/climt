"""Plot correlated-k k(g) curves from a cork LW table — teaching figure.

Demonstrates how a correlated-k distribution looks, and how it FAILS when a band
lumps two spectrally distinct populations (the transparent 8-12 um window + the
opaque H2O nu2 band in the 800-1800 cm^-1 band of the Earth LW table). See
docs/notes/corr-k-explainer-seed.md and experiment #16 in
docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md.

k(g) = the band's absorption coefficients sorted from weakest to strongest, with
g in [0,1] the cumulative fraction of the band. Panel A: the 800-1800 band dry
vs moist (the moist curve develops a low-g transparent "shelf" + a high-g opaque
"cliff" = two populations in one band). Panel B: that lumped band vs spectrally
coherent bands (single smooth curves) at the same humidity.

Run: /Users/joymonteiro/miniconda3/envs/climt/bin/python scripts/experiments/plot_kg_curves.py
"""
import argparse
import os
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TABLE = "climt/_data/cork/correlated_k/earth_low_res_lw_co2refined_linepyline.nc"
T_TARGET, P_TARGET = 290.0, 1e5      # K, Pa
WINDOW_BAND = 4                       # index of 800-1800 cm^-1 in the 6-band table


def gcoord(w):
    """Centred cumulative-g coordinate for a set of g-point weights."""
    return np.cumsum(w) - 0.5 * w


def main():
    ap = argparse.ArgumentParser(description="Plot correlated-k k(g) curves")
    ap.add_argument("--out", default=os.path.join("debug_data", "kg_curves_window_band.png"),
                    help="output PNG path (default: debug_data/kg_curves_window_band.png)")
    args = ap.parse_args()

    with nc.Dataset(TABLE) as ds:
        k = ds["k_coefficients"][...]                       # (1,nband,ngpt,nT,nP,nX)
        gw = ds["gpoint_weights"][...]                       # (nband,ngpt)
        be = ds["band_wavenumber_limits"][...]
        Tg = np.asarray(ds["temperature_grid"][...])
        pg = np.exp(np.asarray(ds["pressure_grid_log"][...]))  # ln(Pa) -> Pa
        xg = np.asarray(ds["h2o_vmr_grid"][...])
    iT = int(np.argmin(np.abs(Tg - T_TARGET)))
    iP = int(np.argmin(np.abs(pg - P_TARGET)))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A: window band, dry -> moist
    g = gcoord(gw[WINDOW_BAND])
    shown = [(0, "X=1e-6 (dry)"), (3, "X=1e-3"), (5, "X=1e-1 (moist)")]
    cmap = plt.get_cmap("viridis", len(shown))
    print(f"Band {WINDOW_BAND} = {be[WINDOW_BAND,0]:.0f}-{be[WINDOW_BAND,1]:.0f} cm^-1, "
          f"T={Tg[iT]:.0f} K, p={pg[iP]/100:.0f} hPa")
    print(f"{'g':>8s} " + " ".join(f"{lab:>14s}" for _, lab in shown))
    for i in range(len(g)):
        vals = " ".join(f"{k[0,WINDOW_BAND,i,iT,iP,ix]:14.2e}" for ix, _ in shown)
        print(f"{g[i]:8.3f} {vals}")
    for j, (ix, lab) in enumerate(shown):
        axA.plot(g, k[0, WINDOW_BAND, :, iT, iP, ix], "-o", color=cmap(j), label=lab)
    axA.set_yscale("log")
    axA.set_xlabel("g  (cumulative fraction of band)")
    axA.set_ylabel(r"k  [m$^2$/kg]")
    axA.set_title(f"Band {WINDOW_BAND}: 800–1800 cm$^{{-1}}$ (window + H$_2$O $\\nu_2$, LUMPED)")
    axA.grid(alpha=0.3, which="both")
    axA.legend()

    # Panel B: lumped band vs spectrally coherent bands, at moist X
    for Bc, lab, c in [(4, "800–1800 (lumped)", "#d62728"),
                       (2, "630–700 (CO2 core, coherent)", "#1f77b4"),
                       (0, "10–500 (H2O rot)", "#2ca02c")]:
        axB.plot(gcoord(gw[Bc]), k[0, Bc, :, iT, iP, 5], "-o", color=c, label=lab)
    axB.set_yscale("log")
    axB.set_xlabel("g")
    axB.set_ylabel(r"k  [m$^2$/kg]")
    axB.set_title("k(g) at X=1e-1: lumped (two-piece) vs coherent (smooth)")
    axB.grid(alpha=0.3, which="both")
    axB.legend()

    fig.suptitle("Correlated-k k(g) curves — linepyline Earth LW table (T=290 K, p=1000 hPa)")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = args.out
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
