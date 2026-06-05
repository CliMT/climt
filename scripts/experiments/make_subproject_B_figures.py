# scripts/experiments/make_subproject_B_figures.py
"""Single entry point for the picket-fence Experiments-post hero figures.

    --figure kg          window-band k(g) shelf+cliff (dry vs moist)
    --figure lbl         PF vs RRTMG vs line-by-line OLR overlay
    --figure rce         before/after moist RCE temperature + tropopause markers
    --figure throughput  PF vs RRTMG throughput bar (us/col)

Reads existing committed data (debug_data/*.npz) and the .npz produced by the
RCE/bench --save hooks; writes one PNG to --out. Wired into the experiment's
sources.yml so `make experiments` regenerates stale figures only.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEBUG = os.path.join(REPO, "debug_data")


def fig_kg(out):
    """Delegate to plot_kg_curves.py --out (it owns the k(g) construction)."""
    subprocess.run(
        [sys.executable, os.path.join(HERE, "plot_kg_curves.py"), "--out", out],
        cwd=REPO, check=True)


def fig_lbl(out):
    """PF/RRTMG scalar OLR vs line-by-line spectrum, moist + dry at 400 ppm."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, kind in zip(axes, ("moist", "dry")):
        d = np.load(os.path.join(DEBUG, f"lbl_olr_spec_{kind}_co2_400ppm.npz"))
        ax.plot(d["nu"], d["olr_spec"], lw=0.4, color="#444",
                label="line-by-line")
        ax.set_title(f"{kind}: LBL={float(d['total']):.1f}  "
                     f"PF={float(d['olr_pf']):.1f}  "
                     f"RRTMG={float(d['olr_rrtmg']):.1f} W/m$^2$")
        ax.set_xlabel("wavenumber (cm$^{-1}$)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("spectral OLR (W/m$^2$/cm$^{-1}$)")
    fig.suptitle("Line-by-line truth vs PF and RRTMG (CO$_2$=400 ppm)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=140)
    plt.close(fig)


def _load_rce(npz_path):
    d = np.load(npz_path)
    cols = sorted({k.split("__")[0] for k in d.files})
    return {c: {k.split("__")[1]: d[k] for k in d.files if k.startswith(c + "__")}
            for c in cols}


def fig_rce(out, before_npz, after_npz):
    """Before/after moist RCE temperature profiles with tropopause markers."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    from tropopause import find_tropopause
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    for ax, npz, title in (
        (axes[0], before_npz, "before: 4-band default"),
        (axes[1], after_npz, "after: 14-band CO$_2$ table"),
    ):
        cols = _load_rce(npz)
        for label, s in cols.items():
            ax.plot(s["T"], s["p_hpa"], "-o", ms=3, label=label)
            tp = find_tropopause(s["T"], s["p_hpa"])
            ax.plot(np.interp(tp["p_curvature"], s["p_hpa"][::-1], s["T"][::-1]),
                    tp["p_curvature"], "*", ms=16, mec="k")
        ax.set_yscale("log"); ax.grid(alpha=0.3)
        ax.set_title(title); ax.set_xlabel("Temperature (K)")
        ax.legend(fontsize=8)
    # Invert ONCE: the axes share y (sharey=True), so a single invert flips both.
    # Inverting per-axis in the loop would double-invert back to upright.
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Pressure (hPa)")
    fig.suptitle("Moist RCE: PF vs RRTMG (★ = θ-curvature tropopause)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=140)
    plt.close(fig)


def fig_throughput(out, npz):
    d = np.load(npz)
    fig, ax = plt.subplots(figsize=(5, 5))
    vals = [float(d["rrtmg_us_per_col"]), float(d["pf_us_per_col"])]
    ax.bar(["RRTMG-LW", "PF-LW\n(14b×8g, njit)"], vals,
           color=["#888", "#1f77b4"])
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.1f}", ha="center", va="bottom")
    ax.set_ylabel("µs / column"); ax.set_title("Longwave throughput (lower is better)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", required=True,
                    choices=["kg", "lbl", "rce", "throughput"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--before-npz")
    ap.add_argument("--after-npz")
    ap.add_argument("--npz")
    a = ap.parse_args()
    if a.figure == "kg":
        fig_kg(a.out)
    elif a.figure == "lbl":
        fig_lbl(a.out)
    elif a.figure == "rce":
        fig_rce(a.out, a.before_npz, a.after_npz)
    elif a.figure == "throughput":
        fig_throughput(a.out, a.npz)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
