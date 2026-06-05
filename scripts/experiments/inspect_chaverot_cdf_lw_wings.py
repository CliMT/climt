"""Diagnostic: plot Chaverot R=500 LW source CDF k(g) for wing vs core CO2 bands.

Tests the hypothesis from the 2026-05-18 picket-fence experiment log:
the cold-stratosphere bias is upstream of the picket-fence integrator —
the Chaverot R=500 source k-distribution itself lacks moderate-strength
CO2-wing line opacity at low pressure that warms RRTMG's stratosphere.

If true, then at low p (~10-100 hPa) and fixed T (~250 K) the wing-band
k(g) curves (500-630 and 700-820 cm^-1) should collapse to values orders
of magnitude below the 630-700 core, with little or no moderate-k
plateau in g-space — i.e. the wings carry only the very weakest line
tails, not a meaningful continuum of moderate absorbers.

Run from the `radiation` env (the one with exo_k installed) — NOT climt:
    conda activate radiation
    python examples/inspect_chaverot_cdf_lw_wings.py
"""
import os

import matplotlib.pyplot as plt
import numpy as np

import exo_k

CHAVEROT_DIR = "/Users/joymonteiro/github/chaverot/data/correlated-k_tables"
SOURCE = "Earth_376ppmCO2_R500_10-30k.ktable.SI.h5"
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "debug_data")

# Three target bands: two wings flanking the CO2 core
BAND_EDGES = np.array([500.0, 630.0, 700.0, 820.0])
BAND_LABELS = ["500-630 (wing)", "630-700 (CO2 core)", "700-820 (wing)"]
BAND_COLORS = ["C0", "C3", "C2"]

# Fixed atmospheric state for the slice
T_TARGET = 250.0  # K
# Representative low-stratosphere H2O VMR (dry). Chaverot grid usually includes
# something close to 1e-6. We snap to the nearest grid value.
X_H2O_TARGET = 1e-6
# Pressures to overlay (Pa). Spans low-strat through upper-trop.
P_TARGETS_PA = [1e2, 1e3, 1e4, 1e5, 1e6]  # 1, 10, 100, 1000, 10000 hPa


def _nearest_idx(grid, value):
    return int(np.argmin(np.abs(np.asarray(grid) - value)))


def main():
    src = os.path.join(CHAVEROT_DIR, SOURCE)
    print(f"Loading {src}")
    kt = exo_k.Ktable5d(filename=src)
    print(f"  shape (P, T, X, nu, g) = "
          f"({kt.Np}, {kt.Nt}, {kt.Nx}, {kt.Nw}, {kt.Ng})")
    print(f"  P (Pa) range: {kt.pgrid.min():.3g} .. {kt.pgrid.max():.3g}")
    print(f"  T (K)  range: {kt.tgrid.min():.3g} .. {kt.tgrid.max():.3g}")
    print(f"  X_H2O grid (sample): {np.asarray(kt.xgrid)[:5]} ... "
          f"{np.asarray(kt.xgrid)[-3:]}")
    print(f"  Native nu_edges: {kt.wnedges.min():.1f} .. "
          f"{kt.wnedges.max():.1f} cm^-1 across {kt.Nw} native bands")

    # Re-bin to our 3 diagnostic bands. This averages over the native R=500
    # nodes inside each band but keeps the full native g-grid.
    kt_binned = kt.copy()
    kt_binned.bin_down(BAND_EDGES)
    print(f"  after bin_down -> Nw={kt_binned.Nw}, Ng={kt_binned.Ng}")

    # Axis indices for the slice
    it = _nearest_idx(kt_binned.tgrid, T_TARGET)
    ix = _nearest_idx(kt_binned.xgrid, X_H2O_TARGET)
    T_actual = float(np.asarray(kt_binned.tgrid)[it])
    X_actual = float(np.asarray(kt_binned.xgrid)[ix])
    print(f"  slicing at T={T_actual:.2f} K (idx {it}), "
          f"X_H2O={X_actual:.3g} (idx {ix})")

    g = np.asarray(kt_binned.ggrid)  # (Ng,)
    k = kt_binned.kdata  # (Np, Nt, Nx, Nw, Ng), units m^2/molecule
    p_grid_pa = np.asarray(kt_binned.pgrid)

    # One panel per pressure, all three bands overlaid.
    npan = len(P_TARGETS_PA)
    fig, axes = plt.subplots(1, npan, figsize=(4.2 * npan, 4.0), sharey=True)
    if npan == 1:
        axes = [axes]

    for ax, p_target in zip(axes, P_TARGETS_PA):
        ip = _nearest_idx(p_grid_pa, p_target)
        p_actual_pa = float(p_grid_pa[ip])
        for b, (label, color) in enumerate(zip(BAND_LABELS, BAND_COLORS)):
            kline = k[ip, it, ix, b, :]
            kline = np.maximum(kline, 1e-60)  # safe log
            ax.semilogy(g, kline, color=color, lw=1.6, label=label)
        ax.set_title(f"p = {p_actual_pa/100:.1f} hPa")
        ax.set_xlabel("g (cumulative)")
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("k  [m$^2$/molecule]")
    axes[-1].legend(loc="lower right", fontsize=8)
    fig.suptitle(
        f"Chaverot R=500 LW source k(g): wing vs core\n"
        f"T = {T_actual:.0f} K, X_H2O = {X_actual:.1e}, "
        f"source = {SOURCE}",
        fontsize=10,
    )
    fig.tight_layout()
    out_png = os.path.join(OUT_DIR, "chaverot_cdf_lw_wings_kg.png")
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(out_png, dpi=140)
    print(f"  wrote {out_png}")

    # Quantitative summary: band-mean k (linear average over g, weighted by
    # uniform g-weights since native g-grid is already on the unit interval)
    # and ratio wing/core at each diagnostic pressure.
    print()
    print("Band-mean k(g) [m^2/molecule] at T = "
          f"{T_actual:.0f} K, X_H2O = {X_actual:.1e}:")
    print(f"{'p [hPa]':>10s}  "
          f"{'wing 500-630':>14s}  {'core 630-700':>14s}  "
          f"{'wing 700-820':>14s}  {'wing/core (left)':>18s}  "
          f"{'wing/core (right)':>18s}")
    gw = np.asarray(kt_binned.weights)
    for p_target in P_TARGETS_PA:
        ip = _nearest_idx(p_grid_pa, p_target)
        p_actual_pa = float(p_grid_pa[ip])
        means = []
        for b in range(3):
            kline = k[ip, it, ix, b, :]
            means.append(float(np.sum(kline * gw)))
        wL, core, wR = means
        ratL = wL / core if core > 0 else np.nan
        ratR = wR / core if core > 0 else np.nan
        print(f"{p_actual_pa/100:10.2f}  "
              f"{wL:14.3e}  {core:14.3e}  {wR:14.3e}  "
              f"{ratL:18.3e}  {ratR:18.3e}")


if __name__ == "__main__":
    main()
