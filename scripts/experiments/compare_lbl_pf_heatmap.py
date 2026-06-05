"""Heatmaps of PF/LBL band-mean k ratio over (band, pressure) for multiple
profiles and PF tables.

Reads:
  - examples/lbl_k_iso.nc     (from lbl_band_mean_k.py --profile iso)
  - examples/lbl_k_adiabat.nc (from lbl_band_mean_k.py --profile adiabat)
  - all candidate PF tables under climt/_data/picket_fence/correlated_k/

Output:
  - examples/lbl_pf_heatmap.png

Run in the climt env:
    conda run -n climt python examples/compare_lbl_pf_heatmap.py
"""
import glob
import os

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr
from scipy.io import netcdf_file

H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
K_BOLT = 1.380649e-23
TABLE_DIR = "climt/_data/picket_fence/correlated_k"
LBL_FILES = {
    "iso (T=250 K)":     "examples/lbl_k_iso.nc",
    "dry adiabat (Ts=288 K)": "examples/lbl_k_adiabat.nc",
}


def planck_wavenumber(nu_cm, T):
    nu_m = nu_cm * 100.0
    c1 = 2.0 * H_PLANCK * C_LIGHT**2
    c2 = H_PLANCK * C_LIGHT / K_BOLT
    return c1 * nu_m**3 / (np.exp(c2 * nu_m / T) - 1.0)


def lbl_band_mean(nu, k_at_p, T, lo, hi):
    m = (nu >= lo) & (nu < hi)
    if not np.any(m):
        return np.nan
    B = planck_wavenumber(nu[m], T)
    return np.trapezoid(B * k_at_p[m], nu[m]) / np.trapezoid(B, nu[m])


def load_pf(nc_path):
    nc = netcdf_file(nc_path, "r", mmap=False)
    out = {
        "band_lims": nc.variables["band_wavenumber_limits"][:].copy(),
        "w":         nc.variables["gpoint_weights"][:].copy(),
        "k":         nc.variables["k_coefficients"][:].copy(),  # (gas, b, g, T, P, X)
        "T_grid":    nc.variables["temperature_grid"][:].copy(),
        "p_log":     nc.variables["pressure_grid_log"][:].copy(),
        "x_grid":    nc.variables["h2o_vmr_grid"][:].copy(),
    }
    nc.close()
    return out


def bracket(grid, v):
    i = np.searchsorted(grid, v) - 1
    i = max(0, min(i, len(grid) - 2))
    f = (v - grid[i]) / (grid[i + 1] - grid[i])
    f = max(0.0, min(1.0, f))
    return i, f


def pf_k_at(table, T, p_pa, h2o):
    iT, fT = bracket(table["T_grid"], T)
    iP, fP = bracket(table["p_log"], np.log(p_pa))
    iX, fX = bracket(table["x_grid"], h2o)
    sub = table["k"][0, :, :, iT:iT+2, iP:iP+2, iX:iX+2]
    k_interp = np.einsum(
        "bgtpx,t,p,x->bg", sub,
        np.array([1 - fT, fT]), np.array([1 - fP, fP]), np.array([1 - fX, fX]),
    )
    return (k_interp * table["w"]).sum(axis=1)  # (band,)


def build_ratio(table, nu, k_lbl_pT, p_pa, T_p, h2o):
    """Return ratio (n_band, n_p)."""
    band_lims = table["band_lims"]
    nb = band_lims.shape[0]
    np_ = len(p_pa)
    out = np.full((nb, np_), np.nan)
    for j, (p, Tj) in enumerate(zip(p_pa, T_p)):
        k_pf = pf_k_at(table, Tj, p, h2o)
        for i, (lo, hi) in enumerate(band_lims):
            k_l = lbl_band_mean(nu, k_lbl_pT[j], Tj, lo, hi)
            if k_l > 0:
                out[i, j] = k_pf[i] / k_l
    return out


def main():
    tables = []
    for pat in ("earth_low_res_lw_co2refined*.nc",
                "earth_low_res_lw_10band*.nc",
                "earth_low_res_lw_*tightcore*.nc"):
        tables += sorted(glob.glob(os.path.join(TABLE_DIR, pat)))
    tables = list(dict.fromkeys(tables))
    if not tables:
        raise SystemExit("no PF tables matched")

    profiles = {}
    for label, path in LBL_FILES.items():
        if not os.path.exists(path):
            raise SystemExit(f"missing {path} — run lbl_band_mean_k.py first")
        ds = xr.open_dataset(path)
        k = ds["kappa"].values
        if ds["kappa"].dims[0] != "p":
            k = k.T
        p_pa = ds["p"].values
        T_p = ds["T"].values
        print(f"LBL {label}: {len(p_pa)} levels, "
              f"p=[{p_pa.min():.4g}, {p_pa.max():.4g}] Pa, "
              f"T=[{T_p.min():.1f}, {T_p.max():.1f}] K, "
              f"nu=[{ds['nu'].values.min():.3g}, {ds['nu'].values.max():.3g}] cm⁻¹, "
              f"kappa=[{k.min():.3g}, {k.max():.3g}]")
        profiles[label] = {
            "nu": ds["nu"].values,
            "p":  p_pa,
            "T":  T_p,
            "k":  k,
        }

    nrows = len(profiles)
    ncols = len(tables)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(3.6 * ncols, 4.2 * nrows),
        squeeze=False, sharey=True,
    )
    norm = mcolors.LogNorm(vmin=0.05, vmax=20.0)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("lightgray")

    for r, (plabel, pdata) in enumerate(profiles.items()):
        for c, tpath in enumerate(tables):
            ax = axes[r, c]
            table = load_pf(tpath)
            ratio = build_ratio(
                table, pdata["nu"], pdata["k"], pdata["p"], pdata["T"], h2o=1e-6,
            )
            nb = ratio.shape[0]
            band_centers = 0.5 * (table["band_lims"][:, 0] + table["band_lims"][:, 1])
            finite = np.isfinite(ratio)
            print(f"  {os.path.basename(tpath)}: "
                  f"ratio shape={ratio.shape}, "
                  f"finite={finite.sum()}/{ratio.size}, "
                  f"range=[{np.nanmin(ratio):.3g}, {np.nanmax(ratio):.3g}], "
                  f"median={np.nanmedian(ratio):.3g}")
            # Build edges: x = band index half-step, y = log-p midpoints
            x_edges = np.arange(nb + 1) - 0.5
            p_hpa = pdata["p"] / 100.0
            log_p = np.log(p_hpa)
            y_edges = np.empty(len(p_hpa) + 1)
            y_edges[1:-1] = np.exp(0.5 * (log_p[:-1] + log_p[1:]))
            y_edges[0] = np.exp(log_p[0] - 0.5 * (log_p[1] - log_p[0]))
            y_edges[-1] = np.exp(log_p[-1] + 0.5 * (log_p[-1] - log_p[-2]))
            masked = np.ma.masked_invalid(ratio.T)
            ax.set_facecolor("lightgray")
            im = ax.pcolormesh(
                x_edges, y_edges, masked,
                cmap=cmap, norm=norm, shading="flat",
            )
            ax.set_yscale("log")
            ax.set_ylim(y_edges[-1], y_edges[0])  # inverted: surface at bottom
            ax.set_xticks(np.arange(nb))
            ax.set_xticklabels([f"{c:.0f}" for c in band_centers],
                               rotation=45, fontsize=7)
            if r == 0:
                ax.set_title(os.path.basename(tpath).replace(
                    "earth_low_res_lw_", "").replace(".nc", ""), fontsize=8)
            if c == 0:
                ax.set_ylabel(f"{plabel}\np (hPa)")
            ax.set_xlabel("band center (cm⁻¹)")

    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.012, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("PF / LBL  (1 = perfect, <1 PF under, >1 PF over)")
    fig.suptitle("Planck-weighted band-mean k: PF vs LBL  (376 ppm CO2, dry)",
                 fontsize=12, y=0.995)
    out = "examples/lbl_pf_heatmap.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
