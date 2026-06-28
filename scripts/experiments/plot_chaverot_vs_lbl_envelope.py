"""Plot native Chaverot R=500 vs LBL across the (T, p) envelope.

Reads debug_data/chaverot_vs_lbl_native.npz and chaverot_vs_lbl_lbl.npz (written
by the two extraction scripts) and produces two figures in debug_data/:

  chaverot_vs_lbl_spectra.png  -- kappa(nu) overlay, 480-820 cm^-1, panel grid
                                  rows=T, cols=p. Native R=500 (step) vs LBL.
                                  The 500-630 far-wing band is shaded.
  chaverot_vs_lbl_ratio.png    -- native/LBL band-mean ratio heatmaps for the
                                  three diagnostic 15-um bands across (T, p).

Run from the radiation env (has matplotlib + numpy):
    /Users/joymonteiro/miniconda3/envs/radiation/bin/python \
        scripts/experiments/plot_chaverot_vs_lbl_envelope.py
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

DBG = os.path.join(os.path.dirname(__file__), "..", "..", "debug_data")
NAT = np.load(os.path.join(DBG, "chaverot_vs_lbl_native.npz"))
LBL = np.load(os.path.join(DBG, "chaverot_vs_lbl_lbl.npz"))

# Only p=1e4 and p=1e5 Pa are grid-converged in the LBL reference: at lower p
# the CO2 lines are Doppler-narrow (~1e-3 cm^-1) and dnu=0.1 under-samples them,
# so those band-means are unreliable on OUR side (not source defects). Restrict
# to the converged, radiatively-relevant strat (1e4) + troposphere (1e5) columns.
# (verified via dnu 0.1->0.005 convergence: p=1e4 stable to <2%; p<=1e3 not.)
P_KEEP = [1.0e4, 1.0e5]
T = NAT["T_nodes"]
P_ALL = list(NAT["p_nodes"])
PCOL = [P_ALL.index(p) for p in P_KEEP]
P = np.array(P_KEEP)
BANDS = ["500_630", "630_700", "700_800"]
BAND_LABEL = {"500_630": "500-630", "630_700": "630-700", "700_800": "700-800"}


def plot_spectra():
    nT, nP = len(T), len(P)
    fig, axes = plt.subplots(nT, nP, figsize=(4.0 * nP, 3.0 * nT),
                             sharex=True, squeeze=False)
    wn_n = NAT["wns_spec"]
    wn_l = LBL["wns_spec"]
    for a in range(nT):
        for b in range(nP):
            ax = axes[a][b]
            bb = PCOL[b]
            ax.axvspan(500, 630, color="0.88", zorder=0)
            ax.step(wn_n, NAT["spec"][a, bb], where="mid", lw=1.1,
                    color="C3", label="Chaverot R=500")
            ax.plot(wn_l, LBL["spec"][a, bb], lw=0.7, color="C0",
                    label="LBL (linepyline)")
            ax.set_yscale("log")
            ax.set_xlim(480, 820)
            ax.set_ylim(1e-9, 1e2)
            ax.set_title(f"T={T[a]:.0f} K, p={P[b]:.0f} Pa", fontsize=9)
            if a == nT - 1:
                ax.set_xlabel("wavenumber (cm$^{-1}$)")
            if b == 0:
                ax.set_ylabel("$\\kappa$ (m$^2$/kg)")
            if a == 0 and b == 0:
                ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("CO$_2$ mass absorption: Chaverot R=500 source vs line-by-line "
                 "(dry, 376 ppm)\n500-630 cm$^{-1}$ far wing shaded",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(DBG, "chaverot_vs_lbl_spectra.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


def plot_ratio():
    fig, axes = plt.subplots(1, len(BANDS), figsize=(5.2 * len(BANDS), 4.2),
                             squeeze=False)
    for j, band in enumerate(BANDS):
        ax = axes[0][j]
        ratio = (NAT[f"band_{band}"] / LBL[f"band_{band}"])[:, PCOL]
        im = ax.imshow(ratio, origin="upper", aspect="auto",
                       norm=LogNorm(vmin=0.1, vmax=60), cmap="RdBu_r")
        ax.set_xticks(range(len(P)))
        ax.set_xticklabels([f"{p:.0f}" for p in P], rotation=45, ha="right")
        ax.set_yticks(range(len(T)))
        ax.set_yticklabels([f"{t:.0f}" for t in T])
        ax.set_xlabel("p (Pa)")
        if j == 0:
            ax.set_ylabel("T (K)")
        ax.set_title(f"{BAND_LABEL[band]} cm$^{{-1}}$\nnative / LBL")
        for a in range(len(T)):
            for b in range(len(P)):
                v = ratio[a, b]
                ax.text(b, a, f"{v:.1f}", ha="center", va="center",
                        fontsize=8,
                        color="white" if (v > 8 or v < 0.3) else "black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label="native / LBL")
    fig.suptitle("Chaverot R=500 / line-by-line band-mean $\\kappa$ ratio "
                 "(>1 = source too opaque)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(DBG, "chaverot_vs_lbl_ratio.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    plot_spectra()
    plot_ratio()
