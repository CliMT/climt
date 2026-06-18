"""Compare CORK k(p, T) to LBL kappa(p, T) across a real T(p) profile.

For each CORK k-table of interest, walks the LBL profile (60 levels of dry
adiabat from `lbl_k_adiabat.nc`) and compares:
  - CORK: g-integrated band-mean k at each layer's (p, T), bilinearly
    interpolated from the table's (log p, T) grid at X_H2O = 1e-6.
  - LBL: uniform-dnu mean of kappa(nu) inside the CORK band's wavenumber
    window, at the same layer.

Both in m^2/kg-air. Plots per-band k(p) overlay + ratio panel, and prints
worst-case ratio per band.

Tests whether the spot-check agreement from experiment #2 (T=250 K,
p=100 hPa, three bands) holds across all CORK bands and all stratospheric
+ tropospheric layers. If yes, the source k-table is exonerated as a
contributor to the cold-strat bias and the next thing to instrument is
the CORK integrator. If no, we've localized which bands/layers diverge.

Run from `radiation` env (needs netCDF4 + h5py + matplotlib).
"""
import os

import h5py
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

LBL_NC = os.path.join(os.path.dirname(__file__), "lbl_k_adiabat.nc")
CORK_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                      "climt", "_data", "cork", "correlated_k")
CORK_TABLES = [
    "earth_low_res_lw_co2refined_gl",
    "earth_low_res_lw_10band_gl",
    "earth_low_res_lw_8gpt_twostretch",  # the "4-band GL shipped" baseline
]
X_H2O_TARGET = 1e-6  # dry — LBL was CO2-only


def load_cork_band_mean_k(table_path):
    """Return (band_edges, T_grid, log_p_grid, x_grid, k_band(b, T, P, X)) in m^2/kg."""
    with nc.Dataset(table_path) as ds:
        kc = ds["k_coefficients"][...]      # (1, nband, ngpt, nT, nP, nX)
        gw = ds["gpoint_weights"][...]       # (nband, ngpt)
        be = ds["band_wavenumber_limits"][...]   # (nband, 2)
        T_grid = ds["temperature_grid"][...]
        log_p_grid = ds["pressure_grid_log"][...]
        x_grid = ds["h2o_vmr_grid"][...]
    # g-integrate: k_band(b, T, P, X) = sum_g k(0, b, g, T, P, X) * gw(b, g)
    k_band = np.einsum("bgtpx,bg->btpx", kc[0], gw)
    return be, T_grid, log_p_grid, x_grid, k_band


def cork_k_at_layer(band_edges, T_grid, log_p_grid, x_grid, k_band,
                  p_pa, T_K, x_h2o):
    """Bilinear interp in (log p, T) at fixed X_H2O index for all bands.
    Returns k(nband,) in m^2/kg-air.
    """
    ix = int(np.argmin(np.abs(x_grid - x_h2o)))
    lp = np.log(p_pa)
    T_clip = np.clip(T_K, T_grid.min(), T_grid.max())
    lp_clip = np.clip(lp, log_p_grid.min(), log_p_grid.max())

    # T bracket
    it = int(np.searchsorted(T_grid, T_clip)) - 1
    it = max(0, min(it, len(T_grid) - 2))
    t_frac = (T_clip - T_grid[it]) / (T_grid[it+1] - T_grid[it])

    # log p bracket — log_p_grid may be ascending or descending; ensure ascending
    asc = log_p_grid[0] < log_p_grid[-1]
    lpg = log_p_grid if asc else log_p_grid[::-1]
    ip = int(np.searchsorted(lpg, lp_clip)) - 1
    ip = max(0, min(ip, len(lpg) - 2))
    p_frac = (lp_clip - lpg[ip]) / (lpg[ip+1] - lpg[ip])
    if not asc:
        ip = len(log_p_grid) - 2 - ip
        p_frac = 1.0 - p_frac

    # k_band shape: (nband, nT, nP, nX) — bilinear in (T, P)
    k00 = k_band[:, it,   ip,   ix]
    k01 = k_band[:, it,   ip+1, ix]
    k10 = k_band[:, it+1, ip,   ix]
    k11 = k_band[:, it+1, ip+1, ix]
    k = ((1 - t_frac) * (1 - p_frac) * k00
         + (1 - t_frac) * p_frac     * k01
         + t_frac       * (1 - p_frac) * k10
         + t_frac       * p_frac     * k11)
    return k  # (nband,)


def lbl_band_mean(nu, kappa_p, wn_lo, wn_hi):
    sel = (nu >= wn_lo) & (nu < wn_hi)
    if not sel.any():
        return np.nan
    return float(np.mean(kappa_p[sel]))


def main():
    print(f"LBL: {LBL_NC}")
    with h5py.File(LBL_NC, "r") as f:
        nu_lbl = np.array(f["nu"])
        p_lbl = np.array(f["p"])
        T_lbl = np.array(f["T"])
        kappa_lbl = np.array(f["kappa"])
        Ts = float(f.attrs["Ts_K"][0])
        co2 = float(f.attrs["CO2_vmr"][0])
    print(f"  T_surf={Ts:.0f} K, CO2={co2*1e6:.0f} ppm, "
          f"60 levels, {len(nu_lbl)} nu nodes")

    for table in CORK_TABLES:
        path = os.path.join(CORK_DIR, table + ".nc")
        if not os.path.exists(path):
            print(f"missing {path}, skipping")
            continue
        print()
        print(f"=== {table} ===")
        be, Tg, lpg, xg, k_band = load_cork_band_mean_k(path)
        nband = be.shape[0]
        print(f"  bands ({nband}):", be.tolist())
        print(f"  T_grid [{Tg.min():.0f}, {Tg.max():.0f}] ({len(Tg)})")
        print(f"  p_grid [{np.exp(lpg).min():.3g}, {np.exp(lpg).max():.3g}] Pa "
              f"({len(lpg)})")

        # Per-layer k for CORK and LBL
        k_cork = np.zeros((len(p_lbl), nband))
        k_lbl = np.zeros_like(k_cork)
        for k, (p, T) in enumerate(zip(p_lbl, T_lbl)):
            k_cork[k, :] = cork_k_at_layer(be, Tg, lpg, xg, k_band,
                                       p, T, X_H2O_TARGET)
            for b in range(nband):
                k_lbl[k, b] = lbl_band_mean(nu_lbl, kappa_lbl[k, :],
                                            be[b, 0], be[b, 1])
        ratio = np.where(k_lbl > 0, k_cork / k_lbl, np.nan)

        # CO2-dominated bands only — LBL is CO2-only, so other bands are
        # apples-to-oranges (CORK includes H2O baseline, LBL doesn't).
        CO2_BANDS = {(500, 630), (500, 800), (500, 820),
                     (630, 700), (700, 800), (700, 820),
                     (1800, 2380), (1800, 3250), (2380, 3250)}
        co2_mask = np.array(
            [(int(be[b, 0]), int(be[b, 1])) in CO2_BANDS for b in range(nband)]
        )

        # Restrict comparison to layers where T is in CORK table's T_grid
        # (no extrapolation below T_grid.min() = 50 K, and stay >= 100 K so we
        # avoid weird extrapolation effects near T_FLOOR=80 K in the LBL run).
        T_min_valid = 100.0
        valid = T_lbl >= T_min_valid
        print(f"  Restricting to layers with T >= {T_min_valid} K: "
              f"{valid.sum()}/{len(T_lbl)} layers, "
              f"p range {p_lbl[valid].min()/100:.2f} - {p_lbl[valid].max()/100:.2f} hPa")

        # Per-band ratio summary in valid T range, CO2-dominated bands only
        print(f"  {'band':<22s}  {'ratio (valid T): min / median / max':<40s}")
        for b in range(nband):
            if not co2_mask[b]:
                continue
            r = ratio[valid, b]
            r_finite = r[np.isfinite(r) & (r > 0)]
            if r_finite.size == 0:
                continue
            print(f"  {int(be[b,0]):4d}-{int(be[b,1]):4d} cm^-1 (CO2)  "
                  f"{r_finite.min():6.3f} / {np.median(r_finite):6.3f} / "
                  f"{r_finite.max():6.3f}")

        # Spot-check at meaningful column layers
        spot_p_hpa = [1000, 500, 200, 100, 50, 20]
        print()
        print(f"  Per-layer CORK vs LBL band-mean k (CO2-dominated bands), m^2/kg-air:")
        print(f"  {'p [hPa]':>9s}  {'T [K]':>6s}  "
              + "  ".join(f"{int(be[b,0])}-{int(be[b,1])}" for b in range(nband)
                          if co2_mask[b]))
        for target_p in spot_p_hpa:
            k_idx = int(np.argmin(np.abs(p_lbl/100 - target_p)))
            row_cork = "  ".join(f"{k_cork[k_idx, b]:.2e}" for b in range(nband)
                                if co2_mask[b])
            row_lbl = "  ".join(f"{k_lbl[k_idx, b]:.2e}" for b in range(nband)
                                 if co2_mask[b])
            row_rat = "  ".join(f"{ratio[k_idx, b]:.2f}" for b in range(nband)
                                 if co2_mask[b])
            label_w = max(len(f"{int(be[b,0])}-{int(be[b,1])}") for b in range(nband)
                          if co2_mask[b])
            print(f"  {p_lbl[k_idx]/100:9.2f}  {T_lbl[k_idx]:6.1f}  CORK:  {row_cork}")
            print(f"  {'':>9s}  {'':>6s}  LBL: {row_lbl}")
            print(f"  {'':>9s}  {'':>6s}  R:   {row_rat}")

        # Plot
        ncol = 3
        nrow = (nband + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow),
                                 sharey=True, squeeze=False)
        p_hpa = p_lbl / 100.0
        for b in range(nband):
            ax = axes[b // ncol, b % ncol]
            ax.semilogx(k_cork[:, b], p_hpa, "C0-", label="CORK")
            ax.semilogx(k_lbl[:, b], p_hpa, "k--", label="LBL", lw=1)
            ax.set_yscale("log")
            ax.invert_yaxis()
            ax.set_xlabel("k [m²/kg]")
            ax.set_title(f"band {b}: {int(be[b,0])}-{int(be[b,1])} cm⁻¹")
            ax.grid(True, alpha=0.3, which="both")
            ax.legend(fontsize=8)
            # ratio inset
            ax2 = ax.twiny()
            ax2.plot(ratio[:, b], p_hpa, "C3-", alpha=0.5)
            ax2.axvline(1.0, color="C3", ls=":", lw=0.8, alpha=0.5)
            ax2.set_xscale("log")
            ax2.set_xlabel("CORK/LBL ratio", color="C3", fontsize=8)
            ax2.tick_params(axis='x', colors="C3", labelsize=7)
            ax2.set_xlim(0.05, 20)
        for col in range(ncol):
            axes[0, col].set_ylim(p_hpa.max(), p_hpa.min())
        axes[0, 0].set_ylabel("p [hPa]")
        fig.suptitle(
            f"CORK k vs LBL kappa on dry adiabat (T_sfc=288 K, CO2=376 ppm)\n"
            f"table: {table}",
            fontsize=10,
        )
        fig.tight_layout()
        out_dir = os.path.join(os.path.dirname(__file__), "..", "..", "debug_data")
        os.makedirs(out_dir, exist_ok=True)
        out = os.path.join(out_dir, f"cork_vs_lbl_k_profile__{table}.png")
        fig.savefig(out, dpi=130)
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
