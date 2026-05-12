"""Sample linepyline line-by-line mass absorption coefficients on a (T, p[, X_H2O]) grid.

This is the slow step in table generation. Each (T, p[, X]) point requires one
call to linepyline.rtm.get_kappa_hitran (and MT_CKD continuum for H2O) per
absorber. With `binning=True` and ~0.5 cm^-1 wavenumber resolution, a 10x15
(T, p) grid over 10-30000 cm^-1 takes ~5-15 min per absorber on a modern laptop.
"""
from __future__ import annotations

import numpy as np
import xarray as xr


def sample_kappa_grid(
    *,
    background_gas,
    absorbers: dict,
    T_grid: np.ndarray,
    p_grid: np.ndarray,
    nu_grid: np.ndarray,
    h2o_vmr_grid: np.ndarray | None = None,
    line_shape: str = "lorentz",
    binning: bool = True,
    include_mtckd_continuum: bool = True,
    surface_gravity: float = 9.81,
):
    """Compute mass absorption coefficient kappa(T, p[, X_H2O], nu).

    Args:
        background_gas: "air", "N2", or None (passed to linepyline.rtm).
        absorbers: dict mapping HITRAN species (H2O, CO2, O3, CH4, NH3) to
            volume mixing ratio. Values are scalars — homogeneous columns only.
            For Hab1-style runs with variable H2O, pass H2O=0.0 here and supply
            h2o_vmr_grid; H2O kappa is computed separately and added with the
            grid VMR.
        T_grid: (nT,) temperatures in K (ascending).
        p_grid: (nP,) pressures in Pa (ascending; log-spaced recommended).
        nu_grid: (nNu,) wavenumber grid in cm^-1 (linepyline expects uniform).
        h2o_vmr_grid: optional (nX,) H2O VMRs. If given, output has an extra
            axis and H2O kappa is recomputed per VMR.
        line_shape, binning, include_mtckd_continuum: passed to linepyline.
        surface_gravity: passed to linepyline.rtm (does not affect kappa, only
            optical-depth integration, but linepyline requires it at __init__).

    Returns:
        kappa array, shape (nT, nP, nNu) or (nT, nP, nX, nNu) if h2o_vmr_grid
        is given. Units: m^2/kg of *total moist atmosphere* (consistent with
        how climt's correlated_k.compute_ck_optical_depth multiplies it by
        column air mass).
    """
    import linepyline

    rtm = linepyline.rtm(background_gas=background_gas,
                         surface_gravity=surface_gravity)

    nu_min, nu_max = float(nu_grid[0]), float(nu_grid[-1])
    dnu = float(np.diff(nu_grid).mean())
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)

    # Pre-compute kappa(T, p, nu) for each non-H2O absorber.
    # linepyline.get_kappa_hitran accepts 1D arrays of (T, p) and returns
    # kappa(p, nu) for a single T. We call it once per T.
    non_h2o = {s: v for s, v in absorbers.items() if s != "H2O"}
    base_kappa = np.zeros((nT, nP, nNu))
    for iT, Tval in enumerate(T_grid):
        T_arr = np.full(nP, Tval)
        for species, vmr in non_h2o.items():
            p_self = p_grid * vmr  # partial pressure (Pa)
            k_da = rtm.get_kappa_hitran(
                species, nu_min, nu_max, dnu,
                p_grid, T_arr, p_self=p_self,
                line_shape=line_shape, binning=binning,
                remove_plinth=False,
            )
            # rtm output is xarray (p, nu); interp wavenumber onto target grid.
            # Clamp to 0: linepyline can produce small negatives at band edges
            # due to Lorentz-wing extrapolation beyond the requested nu range.
            k_arr = np.maximum(k_da.transpose("p", "nu").interp(nu=nu_grid).values, 0.0)
            # Weight by mass mixing ratio of this species relative to total moist gas
            # For non-trace species we use molecular masses; linepyline's internal
            # mass weighting is done in get_optical_depth, so here we multiply by
            # the dimensionless mass fraction.
            mass_frac = _mass_fraction(species, vmr, background_gas, non_h2o)
            base_kappa[iT] += k_arr * mass_frac

    h2o_vmr_used = absorbers.get("H2O", 0.0)
    if h2o_vmr_grid is None:
        if h2o_vmr_used > 0.0:
            base_kappa += _h2o_kappa(
                rtm, nu_grid, p_grid, T_grid, h2o_vmr_used,
                line_shape, binning, include_mtckd_continuum,
                background_gas, non_h2o,
            )
        return base_kappa

    nX = len(h2o_vmr_grid)
    out = np.zeros((nT, nP, nX, nNu))
    for iX, X in enumerate(h2o_vmr_grid):
        h2o_contrib = _h2o_kappa(
            rtm, nu_grid, p_grid, T_grid, float(X),
            line_shape, binning, include_mtckd_continuum,
            background_gas, non_h2o,
        )
        out[:, :, iX, :] = base_kappa + h2o_contrib
    return out


def _h2o_kappa(rtm, nu_grid, p_grid, T_grid, X_h2o,
               line_shape, binning, include_mtckd_continuum,
               background_gas, other_absorbers):
    """H2O kappa (lines + optional MT_CKD continuum), shape (nT, nP, nNu)."""
    nu_min, nu_max = float(nu_grid[0]), float(nu_grid[-1])
    dnu = float(np.diff(nu_grid).mean())
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)
    out = np.zeros((nT, nP, nNu))
    for iT, Tval in enumerate(T_grid):
        T_arr = np.full(nP, Tval)
        p_self = p_grid * X_h2o
        k_lines = rtm.get_kappa_hitran(
            "H2O", nu_min, nu_max, dnu,
            p_grid, T_arr, p_self=p_self,
            line_shape=line_shape, binning=binning,
            remove_plinth=include_mtckd_continuum,
        )
        k = np.maximum(k_lines.transpose("p", "nu").interp(nu=nu_grid).values, 0.0)
        if include_mtckd_continuum:
            k_cont = rtm.get_kappa_mtckd(nu_min, nu_max, dnu,
                                         p_grid, T_arr, p_self)
            k = k + np.maximum(k_cont.transpose("p", "nu").interp(nu=nu_grid).values, 0.0)
        mass_frac = _mass_fraction(
            "H2O", X_h2o, background_gas,
            {**other_absorbers, "H2O": X_h2o},
        )
        out[iT] = k * mass_frac
    return out


# HITRAN molar masses (g/mol). Background gas masses follow phys.gases entries.
_MOLAR_MASS = {"H2O": 18.015, "CO2": 44.01, "O3": 47.998, "CH4": 16.04,
               "NH3": 17.031, "air": 28.97, "N2": 28.014}


def _mass_fraction(species, vmr, background_gas, all_absorbers):
    """Mass fraction of `species` in a moist gas mixture."""
    mean_mw = 0.0
    f_tot = 0.0
    for s, v in all_absorbers.items():
        mean_mw += v * _MOLAR_MASS[s]
        f_tot += v
    if background_gas is not None:
        mean_mw += (1.0 - f_tot) * _MOLAR_MASS[background_gas]
    return (vmr * _MOLAR_MASS[species]) / mean_mw
