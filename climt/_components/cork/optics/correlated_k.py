# climt/_components/cork/optics/correlated_k.py
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

# CO2 runtime axis interpolation: True = geometric (log-k vs log-X_CO2),
# the design default (CO2 band-mean k is convex/saturating in amount, so
# linear-in-value over log-spaced nodes over-estimates it). Flip to False to
# evaluate linear-in-k during the A5 CO2-accuracy check.
_CO2_INTERP_LOGK = True


@njit
def _ck_bracket(grid, v):
    n = grid.shape[0]
    i = np.searchsorted(grid, v) - 1
    if i < 0:
        i = 0
    elif i > n - 2:
        i = n - 2
    f = (v - grid[i]) / (grid[i + 1] - grid[i])
    if f < 0.0:
        f = 0.0
    elif f > 1.0:
        f = 1.0
    return i, f


@njit
def _ck_txx7(k, ig, ib, igp, iT, fT, iP, fP, iX, fX, iC):
    v000 = k[ig, ib, igp, iT,     iP,     iX,     iC]
    v100 = k[ig, ib, igp, iT + 1, iP,     iX,     iC]
    v010 = k[ig, ib, igp, iT,     iP + 1, iX,     iC]
    v110 = k[ig, ib, igp, iT + 1, iP + 1, iX,     iC]
    v001 = k[ig, ib, igp, iT,     iP,     iX + 1, iC]
    v101 = k[ig, ib, igp, iT + 1, iP,     iX + 1, iC]
    v011 = k[ig, ib, igp, iT,     iP + 1, iX + 1, iC]
    v111 = k[ig, ib, igp, iT + 1, iP + 1, iX + 1, iC]
    x0 = (v000 * (1.0 - fT) * (1.0 - fP) + v100 * fT * (1.0 - fP)
          + v010 * (1.0 - fT) * fP + v110 * fT * fP)
    x1 = (v001 * (1.0 - fT) * (1.0 - fP) + v101 * fT * (1.0 - fP)
          + v011 * (1.0 - fT) * fP + v111 * fT * fP)
    return x0 * (1.0 - fX) + x1 * fX


@njit
def _ck_txx_cont(log_cont, ib, iT, fT, iP, fP, iX, fX):
    b000 = log_cont[ib, iT,     iP,     iX]
    b100 = log_cont[ib, iT + 1, iP,     iX]
    b010 = log_cont[ib, iT,     iP + 1, iX]
    b110 = log_cont[ib, iT + 1, iP + 1, iX]
    b001 = log_cont[ib, iT,     iP,     iX + 1]
    b101 = log_cont[ib, iT + 1, iP,     iX + 1]
    b011 = log_cont[ib, iT,     iP + 1, iX + 1]
    b111 = log_cont[ib, iT + 1, iP + 1, iX + 1]
    x0 = (b000 * (1.0 - fT) * (1.0 - fP) + b100 * fT * (1.0 - fP)
          + b010 * (1.0 - fT) * fP + b110 * fT * fP)
    x1 = (b001 * (1.0 - fT) * (1.0 - fP) + b101 * fT * (1.0 - fP)
          + b011 * (1.0 - fT) * fP + b111 * fT * fP)
    return x0 * (1.0 - fX) + x1 * fX


@njit(parallel=True)
def _ck_tau_additive_co2_kernel(
    k, T_grid, p_grid_log, log_x_grid, log_c_grid,
    Tarr, log_p, log_x, log_c, gas_amounts,
    has_cont, log_cont, co2_logk, tau,
):
    ngas, nband, ngpt = k.shape[0], k.shape[1], k.shape[2]
    nlev = Tarr.shape[0]
    ncol = Tarr.shape[1]
    FLOOR = 1e-40
    for i in prange(ncol):
        for kk in range(nlev):
            iT, fT = _ck_bracket(T_grid, Tarr[kk, i])
            iP, fP = _ck_bracket(p_grid_log, log_p[kk, i])
            iX, fX = _ck_bracket(log_x_grid, log_x[kk, i])
            iC, fC = _ck_bracket(log_c_grid, log_c[kk, i])
            for ib in range(nband):
                if has_cont:
                    lc = _ck_txx_cont(log_cont, ib, iT, fT, iP, fP, iX, fX)
                    cont_val = np.exp(lc)
                else:
                    cont_val = 0.0
                for igp in range(ngpt):
                    acc = 0.0
                    for ig in range(ngas):
                        c0 = _ck_txx7(k, ig, ib, igp, iT, fT, iP, fP, iX, fX, iC)
                        c1 = _ck_txx7(k, ig, ib, igp, iT, fT, iP, fP, iX, fX, iC + 1)
                        if co2_logk:
                            l0 = np.log(c0 if c0 > FLOOR else FLOOR)
                            l1 = np.log(c1 if c1 > FLOOR else FLOOR)
                            kv = np.exp(l0 * (1.0 - fC) + l1 * fC)
                        else:
                            kv = c0 * (1.0 - fC) + c1 * fC
                        acc += kv * gas_amounts[ig, kk, i]
                    if has_cont:
                        acc += cont_val * gas_amounts[0, kk, i]
                    tau[ib, igp, kk, i] = acc


def _load_netcdf_table(path):
    try:
        from scipy.io import netcdf_file
    except ImportError as exc:
        raise ImportError(
            "Reading .nc k-tables requires scipy. Install scipy, or use a .npz "
            "table (the shipped demo tables are available as .npz)."
        ) from exc

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
            ``climt/_data/cork/correlated_k/`` are resolved by name,
            preferring ``.npz`` (numpy-native, scipy-free) then falling back
            to ``.nc`` (design-spec format, requires scipy).

    Returns:
        dict-like object exposing ``k_coefficients``, ``gpoint_weights``,
        ``temperature_grid``, ``pressure_grid_log`` and (LW/SW-specific)
        auxiliary arrays.
    """
    if os.path.isfile(name_or_path):
        if name_or_path.endswith(".nc"):
            return _load_netcdf_table(name_or_path)
        return np.load(name_or_path, allow_pickle=True)

    pkg = importlib_resources.files("climt._data.cork.correlated_k")
    npz_path = pkg.joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(npz_path) as f:
        if os.path.isfile(f):
            return np.load(f, allow_pickle=True)

    nc_path = pkg.joinpath(f"{name_or_path}.nc")
    with importlib_resources.as_file(nc_path) as f:
        if os.path.isfile(f):
            return _load_netcdf_table(str(f))

    raise FileNotFoundError(f"No k-table named {name_or_path!r} (.npz or .nc)")


def interpolate_k(table, T, p, h2o_vmr=None, co2_vmr=None):
    """Interpolate k-coefficients in (T, log p[, log X_H2O[, log X_CO2]]) space.

    Args:
        table: loaded k-table. ``k_coefficients`` may be shaped
            ``(ngas, nband, ngpt, nT, nP)`` (bilinear),
            ``(..., nP, nX)`` (trilinear, H2O axis), or
            ``(..., nP, nX, nC)`` (quadrilinear, H2O + CO2 axes).
        T: (ncol,) temperature values, K
        p: (ncol,) pressure values, Pa
        h2o_vmr: (ncol,) H2O mole fraction per column. Required if the table
            has an H2O axis.
        co2_vmr: (ncol,) CO2 mole fraction per column. Required if the table
            has a CO2 axis.

    Returns:
        k_interp: (ngas, nband, ngpt, ncol)
    """
    k = table["k_coefficients"]
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]
    has_co2_axis = k.ndim == 7
    has_x_axis = k.ndim >= 6

    ncol = len(T)
    log_p = np.log(np.maximum(p, 1.0))

    nX = nC = None
    if has_co2_axis:
        ngas, nband, ngpt, nT, nP, nX, nC = k.shape
    elif has_x_axis:
        ngas, nband, ngpt, nT, nP, nX = k.shape
    else:
        ngas, nband, ngpt, nT, nP = k.shape

    if has_x_axis:
        x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
        log_x_grid = np.log(np.maximum(x_grid, 1e-30))
        if h2o_vmr is None:
            raise ValueError(
                "k-table has an h2o_vmr_grid axis but h2o_vmr was not provided"
            )
        x_clamped = np.clip(h2o_vmr, float(x_grid[0]), float(x_grid[-1]))
        log_x = np.log(np.maximum(x_clamped, 1e-30))

    if has_co2_axis:
        c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
        log_c_grid = np.log(np.maximum(c_grid, 1e-30))
        if co2_vmr is None:
            raise ValueError(
                "k-table has a co2_vmr_grid axis but co2_vmr was not provided"
            )
        c_clamped = np.clip(co2_vmr, float(c_grid[0]), float(c_grid[-1]))
        log_c = np.log(np.maximum(c_clamped, 1e-30))

    result = np.zeros((ngas, nband, ngpt, ncol))

    for col in range(ncol):
        iT = max(0, min(np.searchsorted(T_grid, T[col]) - 1, nT - 2))
        fT = (T[col] - T_grid[iT]) / (T_grid[iT + 1] - T_grid[iT])
        fT = max(0.0, min(1.0, fT))

        iP = max(0, min(np.searchsorted(p_grid_log, log_p[col]) - 1, nP - 2))
        fP = (log_p[col] - p_grid_log[iP]) / (p_grid_log[iP + 1] - p_grid_log[iP])
        fP = max(0.0, min(1.0, fP))

        if has_x_axis:
            iX = max(0, min(np.searchsorted(log_x_grid, log_x[col]) - 1, nX - 2))
            fX = (log_x[col] - log_x_grid[iX]) / (
                log_x_grid[iX + 1] - log_x_grid[iX])
            fX = max(0.0, min(1.0, fX))

        if has_co2_axis:
            iC = max(0, min(np.searchsorted(log_c_grid, log_c[col]) - 1, nC - 2))
            fC = (log_c[col] - log_c_grid[iC]) / (
                log_c_grid[iC + 1] - log_c_grid[iC])
            fC = max(0.0, min(1.0, fC))

        def _txx(base):
            # Trilinear (T, P, X_H2O), linear-in-k. base shape (nT, nP, nX).
            v000 = base[iT, iP, iX]
            v100 = base[iT + 1, iP, iX]
            v010 = base[iT, iP + 1, iX]
            v110 = base[iT + 1, iP + 1, iX]
            v001 = base[iT, iP, iX + 1]
            v101 = base[iT + 1, iP, iX + 1]
            v011 = base[iT, iP + 1, iX + 1]
            v111 = base[iT + 1, iP + 1, iX + 1]
            x0 = (v000 * (1 - fT) * (1 - fP) + v100 * fT * (1 - fP)
                  + v010 * (1 - fT) * fP + v110 * fT * fP)
            x1 = (v001 * (1 - fT) * (1 - fP) + v101 * fT * (1 - fP)
                  + v011 * (1 - fT) * fP + v111 * fT * fP)
            return x0 * (1 - fX) + x1 * fX

        def _tp(base):
            # Bilinear (T, P), linear-in-k. base shape (nT, nP).
            return (base[iT, iP] * (1 - fT) * (1 - fP)
                    + base[iT + 1, iP] * fT * (1 - fP)
                    + base[iT, iP + 1] * (1 - fT) * fP
                    + base[iT + 1, iP + 1] * fT * fP)

        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    if has_co2_axis:
                        c0 = _txx(k[ig, ib, igp, :, :, :, iC])
                        c1 = _txx(k[ig, ib, igp, :, :, :, iC + 1])
                        if _CO2_INTERP_LOGK:
                            lc0 = np.log(max(c0, 1e-40))
                            lc1 = np.log(max(c1, 1e-40))
                            result[ig, ib, igp, col] = np.exp(
                                lc0 * (1 - fC) + lc1 * fC)
                        else:
                            result[ig, ib, igp, col] = c0 * (1 - fC) + c1 * fC
                    elif has_x_axis:
                        result[ig, ib, igp, col] = _txx(k[ig, ib, igp])
                    else:
                        result[ig, ib, igp, col] = _tp(k[ig, ib, igp])

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


def compute_ck_optical_depth(table, T, p, gas_amounts, h2o_vmr=None,
                             co2_vmr=None):
    """Compute optical depths from correlated-k table.

    Supports both additive and ESFT overlap methods. When the table has an
    ``h2o_vmr_grid`` axis, ``h2o_vmr`` (nlev, ncol) must be provided so
    interpolate_k can do trilinear (T, log p, log X_H2O) lookup. When the
    table has a ``co2_vmr_grid`` axis, ``co2_vmr`` (nlev, ncol) must be
    provided for quadrilinear lookup.

    Args:
        table: loaded k-table
        T: (nlev, ncol) temperature, K
        p: (nlev, ncol) pressure, Pa
        gas_amounts: (ngas, nlev, ncol) column amount per gas, kg/m^2.
            For premixed-background tables this should be the column mass
            of air (the "gas" is bulk air with H2O as the live axis).
        h2o_vmr: (nlev, ncol) H2O mole fraction per level/column.
        co2_vmr: (nlev, ncol) CO2 mole fraction per level/column, required
            when the table has a co2_vmr_grid axis.

    Returns:
        For additive overlap:
            tau: (nband, ngpt, nlev, ncol) optical depth per layer
        For ESFT overlap:
            tuple of (tau, weights)
    """
    overlap = str(table.get("overlap_method", np.array("additive")))

    if overlap == "esft":
        return _compute_ck_optical_depth_esft(
            table, T, p, gas_amounts, h2o_vmr, co2_vmr)
    else:
        return _compute_ck_optical_depth_additive(
            table, T, p, gas_amounts, h2o_vmr, co2_vmr)


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


def _compute_ck_optical_depth_additive(table, T, p, gas_amounts, h2o_vmr=None,
                                       co2_vmr=None):
    """Additive overlap: sum optical depths from all gases on the same g-point grid."""
    k_data = table["k_coefficients"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    if k_data.ndim == 7:
        # Fast compiled path for the CO2-axis (premixed-bg) table — the default.
        if h2o_vmr is None:
            raise ValueError(
                "k-table has an h2o_vmr_grid axis but h2o_vmr was not provided")
        if co2_vmr is None:
            raise ValueError(
                "k-table has a co2_vmr_grid axis but co2_vmr was not provided")
        return _additive_co2_fast(table, T, p, gas_amounts, h2o_vmr, co2_vmr)

    tau = np.zeros((nband, ngpt, nlev, ncol))
    has_continuum = (
        "continuum_kappa" in table
        and np.asarray(table["continuum_kappa"]).ndim == 4
        and h2o_vmr is not None
    )

    for k_lev in range(nlev):
        x_lev = h2o_vmr[k_lev, :] if h2o_vmr is not None else None
        c_lev = co2_vmr[k_lev, :] if co2_vmr is not None else None
        k_interp = interpolate_k(
            table, T[k_lev, :], p[k_lev, :], h2o_vmr=x_lev, co2_vmr=c_lev)
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


def _additive_co2_fast(table, T, p, gas_amounts, h2o_vmr, co2_vmr):
    """Host wrapper: prep arrays and call the njit CO2-axis tau kernel."""
    k = np.ascontiguousarray(table["k_coefficients"])
    T_grid = np.asarray(table["temperature_grid"], dtype=np.float64)
    p_grid_log = np.asarray(table["pressure_grid_log"], dtype=np.float64)
    x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
    c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
    log_x_grid = np.log(np.maximum(x_grid, 1e-30))
    log_c_grid = np.log(np.maximum(c_grid, 1e-30))

    nband, ngpt = k.shape[1], k.shape[2]
    nlev, ncol = T.shape

    Tarr = np.ascontiguousarray(T, dtype=np.float64)
    log_p = np.log(np.maximum(p, 1.0))
    x_clamped = np.clip(h2o_vmr, float(x_grid[0]), float(x_grid[-1]))
    log_x = np.log(np.maximum(x_clamped, 1e-30))
    c_clamped = np.clip(co2_vmr, float(c_grid[0]), float(c_grid[-1]))
    log_c = np.log(np.maximum(c_clamped, 1e-30))
    gas64 = np.ascontiguousarray(gas_amounts, dtype=np.float64)

    has_cont = ("continuum_kappa" in table
                and np.asarray(table["continuum_kappa"]).ndim == 4)
    if has_cont:
        cont = np.asarray(table["continuum_kappa"], dtype=np.float64)
        log_cont = np.log(np.maximum(cont, 1e-40))
    else:
        log_cont = np.zeros((nband, 1, 1, 1))

    tau = np.zeros((nband, ngpt, nlev, ncol))
    _ck_tau_additive_co2_kernel(
        k, T_grid, p_grid_log, log_x_grid, log_c_grid,
        Tarr, log_p, log_x, log_c, gas64,
        has_cont, log_cont, bool(_CO2_INTERP_LOGK), tau,
    )
    return tau


def _compute_ck_optical_depth_esft(table, T, p, gas_amounts, h2o_vmr=None,
                                   co2_vmr=None):
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
        c_lev = co2_vmr[k_lev, :] if co2_vmr is not None else None
        k_interp = interpolate_k(
            table, T[k_lev, :], p[k_lev, :], h2o_vmr=x_lev, co2_vmr=c_lev)

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
