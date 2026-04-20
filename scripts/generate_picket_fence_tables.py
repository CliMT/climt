"""Convert a Chaverot (Exo_k-format) correlated-k table into climt's picket-fence
netCDF format and re-bin to picket-fence resolution.

Usage:
    python scripts/generate_picket_fence_tables.py \\
        --input /path/to/chaverot_table.h5 \\
        --output climt/_data/picket_fence/correlated_k/earth_low_res_lw.nc \\
        --kind lw \\
        --bands 10,500,1250,2500,3250 \\
        --ngpt 2 \\
        --vmr 0.01

Notes:
    * Requires `exo_k` (GPL-3), not a climt runtime dep. Run from a dev env.
    * Chaverot tables (Zenodo 16795590, CC-BY 4.0) are Ktable5d objects: the
      x-axis is H2O VMR. Use --vmr to select the target H2O mixing ratio
      (log-linear interpolation).  Typical values:
        Earth present-day  ~0.01  (1 % mean tropospheric H2O)
        Mars               ~1e-5  (very dry)
        Venus              ~3e-5  (20-30 ppm H2O)
    * `bands` argument is N+1 band edges in wavenumber (cm^-1). For LW on Earth,
      something like 10, 500, 1250, 2500 cm^-1 gives a window/non-window split.
"""
import argparse
import os

import numpy as np


def _collapse_vmr(ktable5d, vmr):
    """Interpolate a Ktable5d in log-VMR space to produce a (Np,Nt,Nw,Ng) array."""
    xg = ktable5d.xgrid
    log_xg = np.log(xg)
    log_x = np.log(vmr)
    # Clamp to grid edges
    log_x = np.clip(log_x, log_xg[0], log_xg[-1])
    ix = np.searchsorted(log_xg, log_x) - 1
    ix = int(np.clip(ix, 0, len(xg) - 2))
    f = (log_x - log_xg[ix]) / (log_xg[ix + 1] - log_xg[ix])
    # kdata shape: (Np, Nt, Nx, Nw, Ng)
    k0 = ktable5d.kdata[:, :, ix, :, :]
    k1 = ktable5d.kdata[:, :, ix + 1, :, :]
    return np.exp(np.log(np.maximum(k0, 1e-300)) * (1 - f)
                  + np.log(np.maximum(k1, 1e-300)) * f)


def _rebin_to_bands(ktable, band_edges_wn, ngpt, vmr=None):
    """Re-bin an Exo_k Ktable/Ktable5d to new bands with `ngpt` g-points per band.

    Uses Exo_k's built-in wavenumber-bin averaging and then re-quadratures the
    within-band k-distribution to `ngpt` g-points using Gauss-Legendre weights.
    """
    k_binned = ktable.copy()
    k_binned.bin_down(np.asarray(band_edges_wn))

    nT = k_binned.Nt
    nP = k_binned.Np
    nband = k_binned.Nw
    original_g = k_binned.ggrid  # cumulative g-points, shape (Ng,)

    is_5d = hasattr(k_binned, 'Nx')
    if is_5d:
        if vmr is None:
            raise ValueError("--vmr is required for Ktable5d (Chaverot) tables")
        kdata_4d = _collapse_vmr(k_binned, vmr)  # (Np, Nt, Nw, Ng)
    else:
        kdata_4d = k_binned.kdata  # (Np, Nt, Nw, Ng)

    xi, wi = np.polynomial.legendre.leggauss(ngpt)
    xi = 0.5 * (xi + 1.0)
    wi = 0.5 * wi

    k_out = np.zeros((nband, ngpt, nT, nP))
    for b in range(nband):
        for t in range(nT):
            for p in range(nP):
                kline = kdata_4d[p, t, b, :]
                k_out[b, :, t, p] = np.interp(xi, original_g, kline)

    return {
        "k": k_out,
        "g_weights": np.broadcast_to(wi, (nband, ngpt)).copy(),
        "T_grid": np.asarray(k_binned.tgrid),
        "log_p_grid": np.log(np.asarray(k_binned.pgrid)),
        "wn_edges": np.asarray(band_edges_wn),
    }


def _write_climt_netcdf(out_path, rebinned, kind, gas_names=("effective",),
                        rayleigh=None, solar_per_gpt=None):
    from scipy.io import netcdf_file

    nband, ngpt, nT, nP = rebinned["k"].shape
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("gas", len(gas_names))
        nc.createDimension("bounds", 2)

        k = nc.createVariable("k_coefficients", "f4",
                              ("gas", "band", "gpoint", "temperature", "pressure"))
        k[:] = rebinned["k"][np.newaxis]

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = rebinned["g_weights"]

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = rebinned["T_grid"]

        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = rebinned["log_p_grid"]

        edges = rebinned["wn_edges"]
        limits = np.column_stack([edges[:-1], edges[1:]])
        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = limits

        if kind == "lw":
            pf = nc.createVariable("planck_fraction", "f4",
                                   ("band", "gpoint", "temperature"))
            pf[:] = 1.0 / ngpt
        else:
            sf = nc.createVariable("solar_source_per_gpoint", "f4",
                                   ("band", "gpoint"))
            sf[:] = (solar_per_gpt if solar_per_gpt is not None
                     else (1361.0 / nband / ngpt))
            r = nc.createVariable("rayleigh_coefficient", "f4", ("band",))
            r[:] = rayleigh if rayleigh is not None else 0.0

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = "additive"
        nc.resolution = "low"
        nc.source = "chaverot_zenodo_16795590_rebinned"


def main():
    import exo_k

    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Chaverot Exo_k table (.h5 or .nc)")
    ap.add_argument("--output", required=True, help="climt .nc path")
    ap.add_argument("--kind", choices=("lw", "sw"), required=True)
    ap.add_argument("--bands", required=True,
                    help="comma-separated band edges in wavenumber (cm^-1)")
    ap.add_argument("--ngpt", type=int, default=2)
    ap.add_argument("--vmr", type=float, default=None,
                    help="H2O VMR to select from a Ktable5d (required for Chaverot tables)")
    args = ap.parse_args()

    edges = [float(s) for s in args.bands.split(",")]
    try:
        ktable = exo_k.Ktable(filename=args.input)
    except ValueError:
        ktable = exo_k.Ktable5d(filename=args.input)

    rebinned = _rebin_to_bands(ktable, edges, args.ngpt, vmr=args.vmr)
    _write_climt_netcdf(args.output, rebinned, args.kind)
    print(f"wrote {args.output}  shape={rebinned['k'].shape}")


if __name__ == "__main__":
    main()
