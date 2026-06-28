"""Compare Chaverot R=500 vs linepyline LBL opacities in the CO2 15 um band.

Tests whether the wing/core opacity contrast observed in
`inspect_chaverot_cdf_lw_wings.py` (~0.5%) reflects the true CO2 spectrum or
an artifact of the R=500 source CDF.

For each diagnostic band [500-630], [630-700], [700-820] cm^-1, computes:
  - Chaverot full-resolution band-mean k:
      * walks native R=500 sub-bands inside the window
      * per sub-band: integrates k(g) over g with native quadrature weights
      * bandwidth-weighted average across sub-bands in the window
      * linearly interpolates in T between the two grid nodes bracketing T_target
      * units converted m^2/molecule -> m^2/kg via N_A / M_air
  - LBL band-mean kappa: uniform-dnu mean of kappa(nu) over the same window.

Both are reported as m^2/kg-air at matching (p, T). Wing/core ratios within
each dataset and the LBL/Chaverot ratio per band quantify whether the
Chaverot CDF is faithful to LBL in the wings.

Run from the `radiation` env (has exo_k + h5py).
"""
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np

import exo_k

CHAVEROT = ("/Users/joymonteiro/github/chaverot/data/correlated-k_tables/"
            "Earth_376ppmCO2_R500_10-30k.ktable.SI.h5")
LBL_NC = os.path.join(os.path.dirname(__file__), "lbl_k_iso.nc")

T_TARGET = 250.0           # LBL was generated at T_ISO = 250 K
X_H2O_TARGET = 1e-6        # dry slice (LBL has CO2 only)
P_TARGETS_PA = [1e3, 1e4, 1e5]   # 10, 100, 1000 hPa
BANDS_CM = [(500.0, 630.0), (630.0, 700.0), (700.0, 820.0)]
BAND_LABELS = ["wing 500-630", "core 630-700", "wing 700-820"]
BAND_COLORS = ["C0", "C3", "C2"]

N_AVOGADRO = 6.02214076e23
M_AIR = 0.028964  # kg/mol — same as generate_cork_tables.py default
CHAVEROT_M2_PER_MOLEC_TO_M2_PER_KG = N_AVOGADRO / M_AIR


def _bracket_indices(grid, value):
    """Return (i_lo, i_hi, frac) so value = grid[lo] + frac * (grid[hi] - grid[lo])."""
    g = np.asarray(grid)
    if value <= g[0]:
        return 0, 0, 0.0
    if value >= g[-1]:
        return len(g) - 1, len(g) - 1, 0.0
    hi = int(np.searchsorted(g, value))
    lo = hi - 1
    frac = (value - g[lo]) / (g[hi] - g[lo])
    return lo, hi, float(frac)


def chaverot_band_mean_k_m2_per_kg(kt, ip, T_target, ix, wn_lo, wn_hi):
    """Bandwidth-weighted, g-integrated mean k inside [wn_lo, wn_hi], at
    grid pressure index ip, X index ix, and arbitrary T (linear interp).
    Returns m^2/kg-air.
    """
    wnedges = np.asarray(kt.wnedges)
    wns = np.asarray(kt.wns)
    weights = np.asarray(kt.weights)
    kdata = kt.kdata  # (Np, Nt, Nx, Nw, Ng)
    tgrid = np.asarray(kt.tgrid)

    # native sub-bands whose centers fall in window
    centers = 0.5 * (wnedges[:-1] + wnedges[1:])
    band_widths = np.diff(wnedges)
    mask = (centers >= wn_lo) & (centers < wn_hi)
    if not mask.any():
        return np.nan

    # T interpolation
    it_lo, it_hi, t_frac = _bracket_indices(tgrid, T_target)

    # k(p, T_target, X, sub_band) = g-integral over g
    k_sub_lo = np.sum(kdata[ip, it_lo, ix, mask, :] * weights, axis=-1)
    k_sub_hi = np.sum(kdata[ip, it_hi, ix, mask, :] * weights, axis=-1)
    k_sub = k_sub_lo + t_frac * (k_sub_hi - k_sub_lo)  # m^2/molec, per sub-band

    # bandwidth-weighted mean across sub-bands
    w = band_widths[mask]
    k_band_m2_molec = float(np.sum(k_sub * w) / np.sum(w))
    return k_band_m2_molec * CHAVEROT_M2_PER_MOLEC_TO_M2_PER_KG


def chaverot_native_k_spectrum_m2_per_kg(kt, ip, T_target, ix, wn_lo, wn_hi):
    """k vs band-center for every native sub-band overlapping [wn_lo, wn_hi]."""
    wnedges = np.asarray(kt.wnedges)
    weights = np.asarray(kt.weights)
    kdata = kt.kdata
    tgrid = np.asarray(kt.tgrid)
    centers = 0.5 * (wnedges[:-1] + wnedges[1:])
    mask = (centers >= wn_lo) & (centers < wn_hi)
    it_lo, it_hi, t_frac = _bracket_indices(tgrid, T_target)
    k_lo = np.sum(kdata[ip, it_lo, ix, mask, :] * weights, axis=-1)
    k_hi = np.sum(kdata[ip, it_hi, ix, mask, :] * weights, axis=-1)
    k = (k_lo + t_frac * (k_hi - k_lo)) * CHAVEROT_M2_PER_MOLEC_TO_M2_PER_KG
    return centers[mask], k


def lbl_band_mean_kappa(nu, kappa_at_p, wn_lo, wn_hi):
    """Uniform-dnu mean of kappa over [wn_lo, wn_hi]. nu in cm^-1, kappa in m^2/kg."""
    sel = (nu >= wn_lo) & (nu < wn_hi)
    if not sel.any():
        return np.nan
    return float(np.mean(kappa_at_p[sel]))


def _nearest_idx(grid, value):
    return int(np.argmin(np.abs(np.asarray(grid) - value)))


def main():
    print(f"Chaverot:  {CHAVEROT}")
    print(f"LBL:       {LBL_NC}")
    kt = exo_k.Ktable5d(filename=CHAVEROT)

    with h5py.File(LBL_NC, "r") as f:
        nu_lbl = np.array(f["nu"])            # cm^-1
        p_lbl = np.array(f["p"])              # Pa
        T_lbl = np.array(f["T"])              # K, per level
        kappa_lbl = np.array(f["kappa"])      # (Np_lbl, Nnu)
        Ts_K = float(f.attrs["Ts_K"][0])
        CO2_vmr = float(f.attrs["CO2_vmr"][0])
        line_shape = bytes(f.attrs["line_shape"]).decode()
        dnu = float(f.attrs["dnu_cm-1"][0])

    print(f"  LBL: T_iso={Ts_K:.0f} K, CO2={CO2_vmr*1e6:.0f} ppm, "
          f"shape={line_shape}, dnu={dnu} cm^-1")
    print(f"  Chaverot: T_grid {np.asarray(kt.tgrid)} K")
    print(f"  Chaverot: p_grid {np.asarray(kt.pgrid)} Pa")
    print(f"  Chaverot: X_H2O_grid {np.asarray(kt.xgrid)}")

    ix = _nearest_idx(kt.xgrid, X_H2O_TARGET)
    X_actual = float(np.asarray(kt.xgrid)[ix])
    print(f"  Chaverot X_H2O slice: {X_actual:.2e}  (target {X_H2O_TARGET:.0e})")
    print()

    # ── Tables ───────────────────────────────────────────────────────────────
    print(f"Band-mean k [m^2/kg-air] at T = {T_TARGET:.0f} K, "
          f"Chaverot interpolated in T:")
    print(f"{'p [hPa]':>9s}  {'band':<16s}  "
          f"{'Chaverot':>11s}  {'LBL':>11s}  {'LBL/Chav':>9s}  "
          f"{'Chav w/c':>9s}  {'LBL w/c':>9s}")
    for p_target in P_TARGETS_PA:
        ip = _nearest_idx(kt.pgrid, p_target)
        p_actual = float(np.asarray(kt.pgrid)[ip])
        ip_lbl = _nearest_idx(p_lbl, p_target)
        p_actual_lbl = float(p_lbl[ip_lbl])

        chav = []
        lbl = []
        for (wn_lo, wn_hi) in BANDS_CM:
            chav.append(chaverot_band_mean_k_m2_per_kg(
                kt, ip, T_TARGET, ix, wn_lo, wn_hi))
            lbl.append(lbl_band_mean_kappa(
                nu_lbl, kappa_lbl[ip_lbl, :], wn_lo, wn_hi))
        chav = np.array(chav)
        lbl = np.array(lbl)
        chav_wc_L = chav[0] / chav[1] if chav[1] > 0 else np.nan
        chav_wc_R = chav[2] / chav[1] if chav[1] > 0 else np.nan
        lbl_wc_L = lbl[0] / lbl[1] if lbl[1] > 0 else np.nan
        lbl_wc_R = lbl[2] / lbl[1] if lbl[1] > 0 else np.nan

        for i, label in enumerate(BAND_LABELS):
            wc_chav = (chav_wc_L if i == 0 else
                       chav_wc_R if i == 2 else 1.0)
            wc_lbl = (lbl_wc_L if i == 0 else
                      lbl_wc_R if i == 2 else 1.0)
            ratio = lbl[i] / chav[i] if chav[i] > 0 else np.nan
            print(f"{p_actual/100:>9.2f}  {label:<16s}  "
                  f"{chav[i]:>11.3e}  {lbl[i]:>11.3e}  "
                  f"{ratio:>9.3f}  "
                  f"{wc_chav:>9.3e}  {wc_lbl:>9.3e}")
        if abs(p_actual - p_actual_lbl) / p_actual > 0.1:
            print(f"          (NOTE: Chaverot p={p_actual/100:.2f} hPa snapped "
                  f"to LBL p={p_actual_lbl/100:.2f} hPa)")
        print()

    # ── Plot k(nu) overlay at one diagnostic pressure ─────────────────────
    p_plot = 1e4  # 100 hPa
    ip = _nearest_idx(kt.pgrid, p_plot)
    ip_lbl = _nearest_idx(p_lbl, p_plot)
    p_actual = float(np.asarray(kt.pgrid)[ip]) / 100.0
    p_actual_lbl = float(p_lbl[ip_lbl]) / 100.0

    fig, ax = plt.subplots(figsize=(11, 5))
    nu_lo, nu_hi = 480.0, 840.0
    sel_lbl = (nu_lbl >= nu_lo) & (nu_lbl <= nu_hi)
    # bin LBL down to 1 cm^-1 for readability
    bin_edges = np.arange(nu_lo, nu_hi + 1.0, 1.0)
    bin_mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    k_lbl_binned = np.zeros_like(bin_mids)
    for i in range(len(bin_mids)):
        s = (nu_lbl >= bin_edges[i]) & (nu_lbl < bin_edges[i+1])
        k_lbl_binned[i] = np.mean(kappa_lbl[ip_lbl, s]) if s.any() else np.nan
    ax.semilogy(bin_mids, k_lbl_binned, color="black", lw=0.8,
                label=f"LBL (linepyline), 1 cm$^{{-1}}$ mean, T=250 K", zorder=2)

    for (wn_lo, wn_hi), label, color in zip(BANDS_CM, BAND_LABELS, BAND_COLORS):
        nu_c, k_c = chaverot_native_k_spectrum_m2_per_kg(
            kt, ip, T_TARGET, ix, wn_lo, wn_hi)
        ax.semilogy(nu_c, k_c, "o-", color=color, ms=2.5, lw=0.8,
                    label=f"Chaverot R=500, {label}", zorder=3)
        ax.axvspan(wn_lo, wn_hi, color=color, alpha=0.06)

    ax.set_xlabel(r"wavenumber [cm$^{-1}$]")
    ax.set_ylabel(r"k  [m$^2$/kg-air]")
    ax.set_title(
        f"CO$_2$ 15 µm band: Chaverot R=500 (g-mean per native band) vs LBL\n"
        f"p ~ {p_actual:.1f} hPa, T = {T_TARGET:.0f} K, "
        f"CO$_2$ = 376 ppm, dry"
    )
    ax.set_xlim(nu_lo, nu_hi)
    ax.legend(loc="lower center", ncol=2, fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    out_dir = os.path.join(os.path.dirname(__file__), "..", "..", "debug_data")
    os.makedirs(out_dir, exist_ok=True)
    out_png = os.path.join(out_dir, "chaverot_vs_lbl_wings_kspec.png")
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
