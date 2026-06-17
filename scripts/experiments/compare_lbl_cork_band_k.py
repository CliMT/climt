"""Compare Planck-weighted band-mean k from LBL vs CORK correlated-k tables.

Reads:
  - examples/lbl_k_iso250.nc  (from lbl_band_mean_k.py)
  - all earth_low_res_lw_co2refined*.nc tables shipped under
    climt/_data/cork/correlated_k/

For each table:
  k_CORK(band)  = Σ_gpt w_gpt × k_coefficients(band, gpt; T=250 K, p=P_REF, q=0)
  k_LBL(band) = ∫_band B(ν, 250 K) k_LBL(ν) dν  /  ∫_band B(ν, 250 K) dν
              (k_LBL(ν) at the same p as CORK, interpolated onto the LBL p-grid)

Prints a table per .nc with: band edges, k_CORK, k_LBL, ratio.

Run in the climt env:
    conda run -n climt python examples/compare_lbl_cork_band_k.py
"""
import glob
import os

import numpy as np
import xarray as xr
from scipy.io import netcdf_file

T_REF = 250.0
P_REF_HPA = 10.0           # reference pressure for headline comparison
H2O_VMR_REF = 0.0           # dry
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
K_BOLT = 1.380649e-23
TABLE_DIR = "climt/_data/cork/correlated_k"
LBL_FILE = "examples/lbl_k_iso250.nc"


def planck_wavenumber(nu_cm, T):
    """Planck radiance B(ν, T) per cm⁻¹ (arbitrary units OK — we normalize)."""
    nu_m = nu_cm * 100.0  # cm⁻¹ → m⁻¹
    c1 = 2.0 * H_PLANCK * C_LIGHT**2
    c2 = H_PLANCK * C_LIGHT / K_BOLT
    return c1 * nu_m**3 / (np.exp(c2 * nu_m / T) - 1.0)


def lbl_band_mean(nu, k_at_p, T, lo, hi):
    """Planck-weighted band mean of k_LBL over [lo, hi] cm⁻¹."""
    m = (nu >= lo) & (nu < hi)
    if not np.any(m):
        return np.nan
    B = planck_wavenumber(nu[m], T)
    return np.trapz(B * k_at_p[m], nu[m]) / np.trapz(B, nu[m])


def cork_band_mean(nc_path, T_ref, p_ref_pa, h2o_ref):
    """Return (band_edges (nband,2), k_cork (nband,), p_grid_pa, T_grid, h2o_grid)."""
    nc = netcdf_file(nc_path, "r", mmap=False)
    band_lims = nc.variables["band_wavenumber_limits"][:].copy()
    w = nc.variables["gpoint_weights"][:].copy()        # (band, gpt)
    k = nc.variables["k_coefficients"][:].copy()        # (gas, band, gpt, T, P, X)
    T_grid = nc.variables["temperature_grid"][:].copy()
    p_log = nc.variables["pressure_grid_log"][:].copy()
    x_grid = nc.variables["h2o_vmr_grid"][:].copy()
    nc.close()

    # Interpolate (T, log p, X) trilinearly. h2o_ref is 0 → first node if grid starts at 0.
    log_p = np.log(p_ref_pa)
    # Find brackets
    def bracket(grid, v):
        i = np.searchsorted(grid, v) - 1
        i = max(0, min(i, len(grid) - 2))
        f = (v - grid[i]) / (grid[i + 1] - grid[i])
        f = max(0.0, min(1.0, f))
        return i, f

    iT, fT = bracket(T_grid, T_ref)
    iP, fP = bracket(p_log, log_p)
    iX, fX = bracket(x_grid, h2o_ref)

    # k[gas=0, band, gpt, iT:iT+2, iP:iP+2, iX:iX+2]
    sub = k[0, :, :, iT:iT+2, iP:iP+2, iX:iX+2]
    # weights along T, P, X
    w_T = np.array([1 - fT, fT])
    w_P = np.array([1 - fP, fP])
    w_X = np.array([1 - fX, fX])
    k_interp = np.einsum("bgtpx,t,p,x->bg", sub, w_T, w_P, w_X)
    k_cork = (k_interp * w).sum(axis=1)  # (band,)
    return band_lims, k_cork


def main():
    lbl = xr.open_dataset(LBL_FILE)
    k_lbl = lbl["kappa"].values  # (p, nu) or (nu, p)
    # Normalize axes: want (p, nu)
    dims = lbl["kappa"].dims
    if dims[0] != "p":
        k_lbl = k_lbl.T
    nu = lbl["nu"].values
    p_lbl = lbl["p"].values  # Pa

    # Interpolate k_LBL to P_REF in log-p
    log_p_lbl = np.log(p_lbl)
    log_p_ref = np.log(P_REF_HPA * 100.0)
    j = np.searchsorted(log_p_lbl, log_p_ref) - 1
    j = max(0, min(j, len(log_p_lbl) - 2))
    fp = (log_p_ref - log_p_lbl[j]) / (log_p_lbl[j+1] - log_p_lbl[j])
    k_at_p = (1 - fp) * k_lbl[j] + fp * k_lbl[j + 1]

    tables = sorted(glob.glob(os.path.join(TABLE_DIR, "earth_low_res_lw_co2refined*.nc")))
    tables += sorted(glob.glob(os.path.join(TABLE_DIR, "earth_low_res_lw_10band*.nc")))
    tables += sorted(glob.glob(os.path.join(TABLE_DIR, "earth_low_res_lw_*tightcore*.nc")))
    tables = list(dict.fromkeys(tables))
    if not tables:
        raise SystemExit(f"no co2refined tables found in {TABLE_DIR}")

    print(f"\nLBL vs CORK band-mean k  @ T={T_REF} K, p={P_REF_HPA} hPa, dry, 376 ppm CO2")
    print("=" * 78)
    for tpath in tables:
        name = os.path.basename(tpath)
        band_lims, k_cork = cork_band_mean(tpath, T_REF, P_REF_HPA * 100.0, H2O_VMR_REF)
        print(f"\n{name}")
        print(f"  {'band (cm⁻¹)':>18s}  {'k_CORK (m²/kg)':>14s}  {'k_LBL (m²/kg)':>14s}  {'CORK/LBL':>8s}")
        for (lo, hi), kp in zip(band_lims, k_cork):
            kl = lbl_band_mean(nu, k_at_p, T_REF, lo, hi)
            ratio = kp / kl if kl > 0 else np.nan
            print(f"  {lo:7.1f} – {hi:7.1f}  {kp:14.4e}  {kl:14.4e}  {ratio:8.3f}")


if __name__ == "__main__":
    main()
