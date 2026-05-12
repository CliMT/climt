"""Write picket-fence correlated-k netCDF tables in climt's schema."""
from __future__ import annotations

import os

import numpy as np
from scipy.io import netcdf_file


def _ensure_parent(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)


def write_lw_table(
    out_path,
    *,
    k_coefficients,
    gpoint_weights,
    T_grid,
    log_p_grid,
    band_edges,
    planck_fraction,
    h2o_vmr_grid=None,
    gas_names=("effective",),
    overlap_method="additive",
    resolution="low",
    source="linepyline",
):
    """Write LW k-table.

    k_coefficients shape:
        (ngas, nband, ngpt, nT, nP) or (ngas, nband, ngpt, nT, nP, nX) if H2O axis.
    """
    _ensure_parent(out_path)
    has_x = (h2o_vmr_grid is not None)
    if has_x:
        ngas, nband, ngpt, nT, nP, nX = k_coefficients.shape
    else:
        ngas, nband, ngpt, nT, nP = k_coefficients.shape
    edges = np.asarray(band_edges)
    limits = np.column_stack([edges[:-1], edges[1:]])

    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("gas", ngas)
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("bounds", 2)
        if has_x:
            nc.createDimension("h2o_vmr", nX)
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure", "h2o_vmr"))
        else:
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure"))
        v[:] = k_coefficients.astype("f4")

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = gpoint_weights.astype("f4")

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = np.asarray(T_grid, dtype="f4")

        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = np.asarray(log_p_grid, dtype="f4")

        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = limits.astype("f4")

        pf = nc.createVariable("planck_fraction", "f4",
                               ("band", "gpoint", "temperature"))
        pf[:] = planck_fraction.astype("f4")

        if has_x:
            xv = nc.createVariable("h2o_vmr_grid", "f4", ("h2o_vmr",))
            xv[:] = np.asarray(h2o_vmr_grid, dtype="f4")

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = overlap_method
        nc.resolution = resolution
        nc.source = source


def write_sw_table(
    out_path, *,
    k_coefficients, gpoint_weights, T_grid, log_p_grid, band_edges,
    solar_source_per_gpoint, rayleigh_coefficient,
    h2o_vmr_grid=None,
    gas_names=("effective",),
    overlap_method="additive", resolution="low", source="linepyline",
):
    _ensure_parent(out_path)
    has_x = (h2o_vmr_grid is not None)
    if has_x:
        ngas, nband, ngpt, nT, nP, nX = k_coefficients.shape
    else:
        ngas, nband, ngpt, nT, nP = k_coefficients.shape
    edges = np.asarray(band_edges)
    limits = np.column_stack([edges[:-1], edges[1:]])

    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("gas", ngas)
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("bounds", 2)
        if has_x:
            nc.createDimension("h2o_vmr", nX)
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure", "h2o_vmr"))
        else:
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure"))
        v[:] = k_coefficients.astype("f4")

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = gpoint_weights.astype("f4")

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = np.asarray(T_grid, dtype="f4")

        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = np.asarray(log_p_grid, dtype="f4")

        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = limits.astype("f4")

        ss = nc.createVariable("solar_source_per_gpoint", "f4", ("band", "gpoint"))
        ss[:] = solar_source_per_gpoint.astype("f4")

        rc = nc.createVariable("rayleigh_coefficient", "f4", ("band",))
        rc[:] = np.asarray(rayleigh_coefficient, dtype="f4")

        if has_x:
            xv = nc.createVariable("h2o_vmr_grid", "f4", ("h2o_vmr",))
            xv[:] = np.asarray(h2o_vmr_grid, dtype="f4")

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = overlap_method
        nc.resolution = resolution
        nc.source = source
