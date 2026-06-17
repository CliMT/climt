"""Generate a single-band, single-g-point, constant-opacity cork
correlated-k table for grey-limit unit testing.

This script has no external dependencies beyond numpy and scipy and can be
run on any machine without the Chaverot source data.

    python scripts/generate_single_band_unit_table.py \\
        --output climt/_data/cork/correlated_k/single_band_unit_lw.nc \\
        --k 1e-2
"""
import argparse
import os

import numpy as np


def write_single_band_unit_lw(out_path, k_value=1e-2,
                              band_wn=(10.0, 3250.0)):
    """Write a single-band, single-g-point LW k-table with constant opacity.

    Args:
        out_path: output netCDF path
        k_value: constant mass absorption coefficient (m^2/kg)
        band_wn: (low, high) wavenumber edges in cm^-1
    """
    from scipy.io import netcdf_file

    T_grid = np.linspace(100.0, 400.0, 8).astype("f4")
    log_p_grid = np.log(np.logspace(-1, 6, 8)).astype("f4")
    nband, ngpt, nT, nP = 1, 1, T_grid.size, log_p_grid.size

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("gas", 1)
        nc.createDimension("bounds", 2)

        k = nc.createVariable("k_coefficients", "f4",
                              ("gas", "band", "gpoint", "temperature", "pressure"))
        k[:] = np.full((1, nband, ngpt, nT, nP), k_value, dtype="f4")

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = np.ones((nband, ngpt), dtype="f4")

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = T_grid
        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = log_p_grid

        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = np.asarray([band_wn], dtype="f4")

        cork = nc.createVariable("planck_fraction", "f4",
                               ("band", "gpoint", "temperature"))
        cork[:] = np.ones((nband, ngpt, nT), dtype="f4")

        nc.gas_names = "co2"
        nc.overlap_method = "additive"
        nc.resolution = "unit"
        nc.source = "synthetic_unit_fixture"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--k", type=float, default=1e-1,
                    help="Constant mass absorption coefficient, m^2/kg")
    args = ap.parse_args()
    write_single_band_unit_lw(args.output, k_value=args.k)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
