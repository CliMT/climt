"""Convert a Chaverot (Exo_k-format) correlated-k table into climt's picket-fence
netCDF format and re-bin to picket-fence resolution.

Usage:
    python scripts/generate_picket_fence_tables.py \\
        --input /path/to/chaverot_table.h5 \\
        --output climt/_data/picket_fence/correlated_k/earth_low_res_lw.nc \\
        --kind lw \\
        --bands 10,500,1250,2500,3250  \\
        --ngpt 2

Notes:
    * Requires `exo_k` (GPL-3), not a climt runtime dep. Run from a dev env.
    * Chaverot tables (Zenodo 16795590, CC-BY 4.0) are pre-mixed. Each output
      table therefore corresponds to one atmospheric scenario (e.g. present-day
      Earth) rather than a per-gas table.
    * `bands` argument is N+1 band edges in wavenumber (cm^-1). For LW on Earth,
      something like 10, 500, 1250, 2500 cm^-1 gives a window/non-window split.
"""
import argparse
import os

import numpy as np


def _rebin_to_bands(ktable, band_edges_wn, ngpt):
    """Re-bin an Exo_k Ktable to new bands with `ngpt` g-points per band.

    Uses Exo_k's built-in wavenumber-bin averaging and then re-quadratures the
    within-band k-distribution to `ngpt` g-points using Gauss-Legendre weights.
    """
    k_new = ktable.copy()
    k_new.bin_down(np.asarray(band_edges_wn))
    nT = k_new.Nt
    nP = k_new.Np
    nband = k_new.Nw
    k_out = np.zeros((nband, ngpt, nT, nP))
    xi, wi = np.polynomial.legendre.leggauss(ngpt)
    xi = 0.5 * (xi + 1.0)
    wi = 0.5 * wi
    for b in range(nband):
        original_g = k_new.ggrid
        for t in range(nT):
            for p in range(nP):
                kline = k_new.kdata[p, t, b, :]
                k_out[b, :, t, p] = np.interp(xi, original_g, kline)
    return {
        "k": k_out,
        "g_weights": np.broadcast_to(wi, (nband, ngpt)).copy(),
        "T_grid": np.asarray(k_new.tgrid),
        "log_p_grid": np.log(np.asarray(k_new.pgrid)),
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
    args = ap.parse_args()

    edges = [float(s) for s in args.bands.split(",")]
    ktable = exo_k.Ktable(filename=args.input)
    rebinned = _rebin_to_bands(ktable, edges, args.ngpt)
    _write_climt_netcdf(args.output, rebinned, args.kind)
    print(f"wrote {args.output}  shape={rebinned['k'].shape}")


if __name__ == "__main__":
    main()
