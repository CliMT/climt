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
    co2_vmr_grid=None,
    gas_names=("effective",),
    overlap_method="additive",
    resolution="low",
    source="linepyline",
    continuum_kappa=None,
):
    """Write LW k-table.

    k_coefficients shape:
        (ngas, nband, ngpt, nT, nP)                  no live-X axis
        (ngas, nband, ngpt, nT, nP, nX_H2O)          H2O axis
        (ngas, nband, ngpt, nT, nP, nX_H2O, nX_CO2)  H2O + CO2 axes

    co2_vmr_grid (optional): (nX_CO2,) log-spaced CO2 volume mixing ratios.
        Requires h2o_vmr_grid. Adds a trailing CO2 dimension to k_coefficients.

    continuum_kappa (optional): band-grey H2O continuum mass-absorption
        (m^2/kg-air), shape (nband, nT, nP, nX_H2O). H2O-only and CO2-independent
        (no CO2 axis) — decouples its quadratic-in-H2O scaling from the line CDF
        (cf. RRTMG's separate tauself/taufor).
    """
    _ensure_parent(out_path)
    has_x = (h2o_vmr_grid is not None)
    has_co2 = (co2_vmr_grid is not None)
    if has_co2 and not has_x:
        raise ValueError("co2_vmr_grid requires an h2o_vmr axis")
    if has_co2:
        ngas, nband, ngpt, nT, nP, nX, nC = k_coefficients.shape
    elif has_x:
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
        if has_co2:
            nc.createDimension("co2_vmr", nC)
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure", "h2o_vmr", "co2_vmr"))
        elif has_x:
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

        if has_co2:
            cvg = nc.createVariable("co2_vmr_grid", "f4", ("co2_vmr",))
            cvg[:] = np.asarray(co2_vmr_grid, dtype="f4")

        if continuum_kappa is not None:
            if not has_x:
                raise ValueError("continuum_kappa requires an h2o_vmr axis")
            ck = np.asarray(continuum_kappa)
            if ck.shape != (nband, nT, nP, nX):
                raise ValueError(
                    f"continuum_kappa shape {ck.shape} != "
                    f"(nband,nT,nP,nX)={(nband, nT, nP, nX)}")
            cv = nc.createVariable(
                "continuum_kappa", "f4",
                ("band", "temperature", "pressure", "h2o_vmr"))
            cv[:] = ck.astype("f4")

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
