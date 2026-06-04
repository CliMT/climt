# climt/_components/picket_fence/optics/correlated_k.py
import os

import importlib_resources
import numpy as np

from ..common import njit, prange


_NETCDF_VARS = (
    "k_coefficients",
    "gpoint_weights",
    "temperature_grid",
    "pressure_grid_log",
    "h2o_vmr_grid",
    "co2_vmr_grid",
    "band_wavenumber_limits",
    "planck_fraction",
    "solar_source_per_gpoint",
    "rayleigh_coefficient",
    "continuum_kappa",
)


def _load_netcdf_table(path):
    from scipy.io import netcdf_file

    with netcdf_file(path, "r", mmap=False) as nc:
        out = {}
        for name in _NETCDF_VARS:
            if name in nc.variables:
                arr = np.asarray(nc.variables[name][:]).copy()
                # scipy.io.netcdf writes big-endian by default; numba kernels
                # reject non-native dtypes, so normalize to native byte order.
                if arr.dtype.byteorder not in ("=", "|"):
                    arr = arr.astype(arr.dtype.newbyteorder("="))
                out[name] = arr

        if "gas_names" in nc.variables:
            raw = nc.variables["gas_names"][:]
            out["gas_names"] = np.asarray(
                [_decode(x) for x in np.atleast_1d(raw)]
            )
        else:
            gn = getattr(nc, "gas_names", None)
            if gn is not None:
                out["gas_names"] = np.asarray(
                    [s.strip() for s in _decode(gn).split(",") if s.strip()]
                )

        for attr in ("overlap_method", "resolution", "background_is_premixed"):
            val = getattr(nc, attr, None)
            if val is not None:
                out[attr] = np.asarray(_decode(val))

        return out


def _decode(x):
    if isinstance(x, bytes):
        return x.decode("utf-8")
    if isinstance(x, np.ndarray) and x.dtype.kind == "S":
        return x.tobytes().decode("utf-8").rstrip("\x00")
    return str(x)


def load_k_table(name_or_path):
    """Load a correlated-k table.

    Args:
        name_or_path: table name (e.g., ``"test_2band_lw"``) or path to a
            ``.npz`` / ``.nc`` file. Shipped tables under
            ``climt/_data/picket_fence/correlated_k/`` are resolved by name,
            preferring ``.nc`` (design-spec format) then falling back to
            ``.npz`` (legacy unit fixtures).

    Returns:
        dict-like object exposing ``k_coefficients``, ``gpoint_weights``,
        ``temperature_grid``, ``pressure_grid_log`` and (LW/SW-specific)
        auxiliary arrays.
    """
    if os.path.isfile(name_or_path):
        if name_or_path.endswith(".nc"):
            return _load_netcdf_table(name_or_path)
        return np.load(name_or_path, allow_pickle=True)

    pkg = importlib_resources.files("climt._data.picket_fence.correlated_k")
    nc_path = pkg.joinpath(f"{name_or_path}.nc")
    with importlib_resources.as_file(nc_path) as f:
        if os.path.isfile(f):
            return _load_netcdf_table(str(f))

    npz_path = pkg.joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(npz_path) as f:
        return np.load(f, allow_pickle=True)


def interpolate_k(table, T, p, h2o_vmr=None):
    """Interpolate k-coefficients in (T, log p[, log X_H2O]) space.

    Args:
        table: loaded k-table. May have ``k_coefficients`` shaped
            ``(ngas, nband, ngpt, nT, nP)`` (bilinear, back-compat) or
            ``(ngas, nband, ngpt, nT, nP, nX)`` when an H2O VMR axis is
            present (trilinear).
        T: (ncol,) temperature values, K
        p: (ncol,) pressure values, Pa
        h2o_vmr: (ncol,) H2O mole fraction per column. Required if the
            table has an H2O axis; ignored otherwise.

    Returns:
        k_interp: (ngas, nband, ngpt, ncol)
    """
    k = table["k_coefficients"]
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]
    has_x_axis = k.ndim == 6

    ncol = len(T)
    log_p = np.log(np.maximum(p, 1.0))

    if has_x_axis:
        ngas, nband, ngpt, nT, nP, nX = k.shape
        x_grid = table["h2o_vmr_grid"]
        log_x_grid = np.log(np.maximum(np.asarray(x_grid, dtype=np.float64), 1e-30))
        if h2o_vmr is None:
            raise ValueError(
                "k-table has an h2o_vmr_grid axis but h2o_vmr was not provided"
            )
        # Clamp to grid edges and take log for log-linear interpolation in X.
        x_clamped = np.clip(h2o_vmr, float(x_grid[0]), float(x_grid[-1]))
        log_x = np.log(np.maximum(x_clamped, 1e-30))
    else:
        ngas, nband, ngpt, nT, nP = k.shape
        nX = None

    result = np.zeros((ngas, nband, ngpt, ncol))

    for col in range(ncol):
        iT = np.searchsorted(T_grid, T[col]) - 1
        iT = max(0, min(iT, nT - 2))
        fT = (T[col] - T_grid[iT]) / (T_grid[iT + 1] - T_grid[iT])
        fT = max(0.0, min(1.0, fT))

        iP = np.searchsorted(p_grid_log, log_p[col]) - 1
        iP = max(0, min(iP, nP - 2))
        fP = (log_p[col] - p_grid_log[iP]) / (p_grid_log[iP + 1] - p_grid_log[iP])
        fP = max(0.0, min(1.0, fP))

        if has_x_axis:
            iX = np.searchsorted(log_x_grid, log_x[col]) - 1
            iX = max(0, min(iX, nX - 2))
            fX = (log_x[col] - log_x_grid[iX]) / (
                log_x_grid[iX + 1] - log_x_grid[iX]
            )
            fX = max(0.0, min(1.0, fX))

            for ig in range(ngas):
                for ib in range(nband):
                    for igp in range(ngpt):
                        # Trilinear: (T, P, X). Two bilinear interpolations
                        # on adjacent X slices, then linear in X.
                        base = k[ig, ib, igp]  # (nT, nP, nX)
                        v000 = base[iT, iP, iX]
                        v100 = base[iT + 1, iP, iX]
                        v010 = base[iT, iP + 1, iX]
                        v110 = base[iT + 1, iP + 1, iX]
                        v001 = base[iT, iP, iX + 1]
                        v101 = base[iT + 1, iP, iX + 1]
                        v011 = base[iT, iP + 1, iX + 1]
                        v111 = base[iT + 1, iP + 1, iX + 1]
                        x0 = (v000 * (1 - fT) * (1 - fP)
                              + v100 * fT * (1 - fP)
                              + v010 * (1 - fT) * fP
                              + v110 * fT * fP)
                        x1 = (v001 * (1 - fT) * (1 - fP)
                              + v101 * fT * (1 - fP)
                              + v011 * (1 - fT) * fP
                              + v111 * fT * fP)
                        result[ig, ib, igp, col] = x0 * (1 - fX) + x1 * fX
        else:
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


def compute_ck_optical_depth(table, T, p, gas_amounts, h2o_vmr=None):
    """Compute optical depths from correlated-k table.

    Supports both additive and ESFT overlap methods. When the table has an
    ``h2o_vmr_grid`` axis, ``h2o_vmr`` (nlev, ncol) must be provided so
    interpolate_k can do trilinear (T, log p, log X_H2O) lookup.

    Args:
        table: loaded k-table
        T: (nlev, ncol) temperature, K
        p: (nlev, ncol) pressure, Pa
        gas_amounts: (ngas, nlev, ncol) column amount per gas, kg/m^2.
            For premixed-background tables this should be the column mass
            of air (the "gas" is bulk air with H2O as the live axis).
        h2o_vmr: (nlev, ncol) H2O mole fraction per level/column.

    Returns:
        For additive overlap:
            tau: (nband, ngpt, nlev, ncol) optical depth per layer
        For ESFT overlap:
            tuple of (tau, weights)
    """
    overlap = str(table.get("overlap_method", np.array("additive")))

    if overlap == "esft":
        return _compute_ck_optical_depth_esft(table, T, p, gas_amounts, h2o_vmr)
    else:
        return _compute_ck_optical_depth_additive(table, T, p, gas_amounts, h2o_vmr)


def interpolate_continuum(table, T, p, h2o_vmr):
    """Interpolate the band-grey H2O continuum mass-absorption (m^2/kg-air).

    Returns (nband, ncol) or None if the table has no decoupled continuum.

    Interpolation is **log-linear in the continuum value** over (T, log p,
    log X_H2O). This is essential: the H2O self-continuum scales ~X^2, and
    linear-in-value interpolation across the widely log-spaced X nodes badly
    OVER-estimates a convex power law between nodes (e.g. ~6x at X~0.024 between
    the 1e-2 and 1e-1 nodes). log-space interpolation is exact for any power law
    in X — this is the whole point of decoupling the continuum from the
    line k-distribution (cf. RRTMG's analytic selffac ~ H2O^2 scaling).
    """
    if "continuum_kappa" not in table:
        return None
    cont = np.asarray(table["continuum_kappa"])
    if cont.ndim != 4:  # needs the (nband, nT, nP, nX) H2O-axis form
        return None
    if h2o_vmr is None:
        return None

    nband, nT, nP, nX = cont.shape
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]
    x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
    log_x_grid = np.log(np.maximum(x_grid, 1e-30))

    _FLOOR = 1e-40
    log_cont = np.log(np.maximum(cont, _FLOOR))  # (nband, nT, nP, nX)

    ncol = len(T)
    log_p = np.log(np.maximum(p, 1.0))
    x_clamped = np.clip(h2o_vmr, float(x_grid[0]), float(x_grid[-1]))
    log_x = np.log(np.maximum(x_clamped, 1e-30))

    out = np.zeros((nband, ncol))
    for col in range(ncol):
        iT = max(0, min(np.searchsorted(T_grid, T[col]) - 1, nT - 2))
        fT = (T[col] - T_grid[iT]) / (T_grid[iT + 1] - T_grid[iT])
        fT = min(1.0, max(0.0, fT))

        iP = max(0, min(np.searchsorted(p_grid_log, log_p[col]) - 1, nP - 2))
        fP = (log_p[col] - p_grid_log[iP]) / (p_grid_log[iP + 1] - p_grid_log[iP])
        fP = min(1.0, max(0.0, fP))

        iX = max(0, min(np.searchsorted(log_x_grid, log_x[col]) - 1, nX - 2))
        fX = (log_x[col] - log_x_grid[iX]) / (log_x_grid[iX + 1] - log_x_grid[iX])
        fX = min(1.0, max(0.0, fX))

        for ib in range(nband):
            b = log_cont[ib]
            x0 = (b[iT, iP, iX] * (1 - fT) * (1 - fP)
                  + b[iT + 1, iP, iX] * fT * (1 - fP)
                  + b[iT, iP + 1, iX] * (1 - fT) * fP
                  + b[iT + 1, iP + 1, iX] * fT * fP)
            x1 = (b[iT, iP, iX + 1] * (1 - fT) * (1 - fP)
                  + b[iT + 1, iP, iX + 1] * fT * (1 - fP)
                  + b[iT, iP + 1, iX + 1] * (1 - fT) * fP
                  + b[iT + 1, iP + 1, iX + 1] * fT * fP)
            out[ib, col] = np.exp(x0 * (1 - fX) + x1 * fX)
    return out


def _compute_ck_optical_depth_additive(table, T, p, gas_amounts, h2o_vmr=None):
    """Additive overlap: sum optical depths from all gases on the same g-point grid."""
    k_data = table["k_coefficients"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    tau = np.zeros((nband, ngpt, nlev, ncol))
    has_continuum = (
        "continuum_kappa" in table
        and np.asarray(table["continuum_kappa"]).ndim == 4
        and h2o_vmr is not None
    )

    for k_lev in range(nlev):
        x_lev = h2o_vmr[k_lev, :] if h2o_vmr is not None else None
        k_interp = interpolate_k(table, T[k_lev, :], p[k_lev, :], h2o_vmr=x_lev)
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    for icol in range(ncol):
                        tau[ib, igp, k_lev, icol] += (
                            k_interp[ig, ib, igp, icol] * gas_amounts[ig, k_lev, icol]
                        )
        if has_continuum:
            # Band-grey continuum: same air-column mass scaling as the premixed
            # line k (gas 0), added equally to every g-point in the band.
            cont = interpolate_continuum(table, T[k_lev, :], p[k_lev, :], x_lev)
            if cont is not None:
                for ib in range(nband):
                    for igp in range(ngpt):
                        for icol in range(ncol):
                            tau[ib, igp, k_lev, icol] += (
                                cont[ib, icol] * gas_amounts[0, k_lev, icol]
                            )

    return tau


def _compute_ck_optical_depth_esft(table, T, p, gas_amounts, h2o_vmr=None):
    """ESFT overlap: outer product of g-points across gases."""
    k_data = table["k_coefficients"]
    gpoint_weights = table["gpoint_weights"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    ngpt_combined = ngpt ** ngas
    tau = np.zeros((nband, ngpt_combined, nlev, ncol))

    combined_weights = compute_esft_weights(gpoint_weights, ngas)

    for k_lev in range(nlev):
        x_lev = h2o_vmr[k_lev, :] if h2o_vmr is not None else None
        k_interp = interpolate_k(table, T[k_lev, :], p[k_lev, :], h2o_vmr=x_lev)

        for ib in range(nband):
            for idx in range(ngpt_combined):
                remainder = idx
                for ig in range(ngas):
                    g_idx = remainder % ngpt
                    remainder //= ngpt
                    for icol in range(ncol):
                        tau[ib, idx, k_lev, icol] += (
                            k_interp[ig, ib, g_idx, icol] * gas_amounts[ig, k_lev, icol]
                        )

    return tau, combined_weights
