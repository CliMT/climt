# climt/_components/picket_fence/optics/correlated_k.py
import os

import importlib_resources
import numpy as np

from ..common import njit, prange


def load_k_table(name_or_path):
    """Load a correlated-k table.

    Args:
        name_or_path: table name (e.g., "test_2band_lw") or path to .npz file

    Returns:
        dict-like npz object
    """
    if os.path.isfile(name_or_path):
        return np.load(name_or_path, allow_pickle=True)

    data_path = importlib_resources.files(
        "climt._data.picket_fence.correlated_k"
    ).joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(data_path) as f:
        return np.load(f, allow_pickle=True)


def interpolate_k(table, T, p):
    """Bilinear interpolation of k-coefficients in (log p, T) space.

    Args:
        table: loaded k-table
        T: (ncol,) temperature values, K
        p: (ncol,) pressure values, Pa

    Returns:
        k_interp: (ngas, nband, ngpt, ncol)
    """
    k = table["k_coefficients"]  # (ngas, nband, ngpt, nT, nP)
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]
    ngas, nband, ngpt, nT, nP = k.shape
    ncol = len(T)

    result = np.zeros((ngas, nband, ngpt, ncol))
    log_p = np.log(np.maximum(p, 1.0))

    for col in range(ncol):
        # Find T indices
        iT = np.searchsorted(T_grid, T[col]) - 1
        iT = max(0, min(iT, nT - 2))
        fT = (T[col] - T_grid[iT]) / (T_grid[iT + 1] - T_grid[iT])
        fT = max(0.0, min(1.0, fT))

        # Find log(p) indices
        iP = np.searchsorted(p_grid_log, log_p[col]) - 1
        iP = max(0, min(iP, nP - 2))
        fP = (log_p[col] - p_grid_log[iP]) / (p_grid_log[iP + 1] - p_grid_log[iP])
        fP = max(0.0, min(1.0, fP))

        # Bilinear interpolation
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    v00 = k[ig, ib, igp, iT, iP]
                    v10 = k[ig, ib, igp, iT + 1, iP]
                    v01 = k[ig, ib, igp, iT, iP + 1]
                    v11 = k[ig, ib, igp, iT + 1, iP + 1]
                    result[ig, ib, igp, col] = (
                        v00 * (1 - fT) * (1 - fP)
                        + v10 * fT * (1 - fP)
                        + v01 * (1 - fT) * fP
                        + v11 * fT * fP
                    )

    return result


def compute_esft_weights(gpoint_weights, ngas):
    """Compute ESFT combined g-point weights for multiple gases.

    For N gases each with G g-points, produces G^N combined weights per band,
    where the combined weight is the product of individual gas weights.

    Args:
        gpoint_weights: (nband, ngpt) per-gas g-point weights (same for all gases)
        ngas: number of gases

    Returns:
        combined_weights: (nband, ngpt^ngas) combined weights
    """
    nband, ngpt = gpoint_weights.shape
    ngpt_combined = ngpt ** ngas

    combined = np.zeros((nband, ngpt_combined))

    for b in range(nband):
        w = gpoint_weights[b]
        for idx in range(ngpt_combined):
            weight = 1.0
            remainder = idx
            for gas in range(ngas):
                g_idx = remainder % ngpt
                remainder //= ngpt
                weight *= w[g_idx]
            combined[b, idx] = weight

    return combined


def compute_ck_optical_depth(table, T, p, gas_amounts):
    """Compute optical depths from correlated-k table.

    Supports both additive and ESFT overlap methods.

    Args:
        table: loaded k-table
        T: (nlev, ncol) temperature, K
        p: (nlev, ncol) pressure, Pa
        gas_amounts: (ngas, nlev, ncol) column amount per gas, kg/m^2

    Returns:
        For additive overlap:
            tau: (nband, ngpt, nlev, ncol) optical depth per layer
        For ESFT overlap:
            tuple of (tau, weights) where
            tau: (nband, ngpt^ngas, nlev, ncol) optical depth per combined g-point
            weights: (nband, ngpt^ngas) combined g-point weights
    """
    overlap = str(table.get("overlap_method", np.array("additive")))

    if overlap == "esft":
        return _compute_ck_optical_depth_esft(table, T, p, gas_amounts)
    else:
        return _compute_ck_optical_depth_additive(table, T, p, gas_amounts)


def _compute_ck_optical_depth_additive(table, T, p, gas_amounts):
    """Additive overlap: sum optical depths from all gases on the same g-point grid."""
    k_data = table["k_coefficients"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    tau = np.zeros((nband, ngpt, nlev, ncol))

    for k_lev in range(nlev):
        k_interp = interpolate_k(table, T[k_lev, :], p[k_lev, :])
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    for icol in range(ncol):
                        tau[ib, igp, k_lev, icol] += (
                            k_interp[ig, ib, igp, icol] * gas_amounts[ig, k_lev, icol]
                        )

    return tau


def _compute_ck_optical_depth_esft(table, T, p, gas_amounts):
    """ESFT overlap: outer product of g-points across gases.

    For combined g-point index idx, decode per-gas g-point indices by
    treating idx as a mixed-radix number (each digit base ngpt).
    """
    k_data = table["k_coefficients"]  # (ngas, nband, ngpt, nT, nP)
    gpoint_weights = table["gpoint_weights"]  # (nband, ngpt)
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    ngpt_combined = ngpt ** ngas
    tau = np.zeros((nband, ngpt_combined, nlev, ncol))

    # Precompute combined weights
    combined_weights = compute_esft_weights(gpoint_weights, ngas)

    for k_lev in range(nlev):
        k_interp = interpolate_k(table, T[k_lev, :], p[k_lev, :])
        # k_interp: (ngas, nband, ngpt, ncol)

        for ib in range(nband):
            for idx in range(ngpt_combined):
                # Decode combined index to per-gas g-point indices
                remainder = idx
                for ig in range(ngas):
                    g_idx = remainder % ngpt
                    remainder //= ngpt
                    for icol in range(ncol):
                        tau[ib, idx, k_lev, icol] += (
                            k_interp[ig, ib, g_idx, icol] * gas_amounts[ig, k_lev, icol]
                        )

    return tau, combined_weights
