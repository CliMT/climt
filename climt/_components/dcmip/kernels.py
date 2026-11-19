"""numba port of the DCMIP initial conditions for tests 4 and 5.

Line-by-line ports of ``test4_baroclinic_wave``
(``fortran/_lib/dcmip/dcmip_initial_conditions_test_4_v3.f90``) and
``test5_tropical_cyclone`` (``..._test_5_v1.f90``), evaluated with pressure
as the vertical coordinate (``zcoords = 0``), one point at a time over levels
and columns as the Cython wrapper ``_dcmip.pyx`` did. All Fortran literals are
double precision (``d0``), so no single-precision constants are needed.
"""

import math

import numpy as np

from ..._core.backend import HAS_NUMBA

if HAS_NUMBA:
    from numba import njit

    def _jit(func):
        return njit(cache=True)(func)
else:
    def _jit(func):
        return func


# ---------------------------------------------------------------------------
# Test 4: baroclinic wave (Jablonowski and Williamson 2006)
# ---------------------------------------------------------------------------
_PI = 4.0 * math.atan(1.0)
_A4 = 6371220.0
_RD4 = 287.0
_G4 = 9.80616
_OMEGA4 = 7.29212e-5
_P0_4 = 100000.0
_DEG2RAD = _PI / 180.0
_ETA_TROPO = 0.2
_U0 = 35.0
_T0_4 = 288.0
_ETA0 = 0.252
_RADIUS = 10.0
_PERT_AMP = 1.0
_PERT_LON = 20.0
_PERT_LAT = 40.0
_ETA_SFC = 1.0
_DELTA_T = 480000.0
_GAMMA4 = 0.005
_A_OMEGA = _A4 * _OMEGA4
_EXPONENT4 = _RD4 * _GAMMA4 / _G4
_Q0_4 = 0.021
_LAT_HW = 2.0 * _PI / 9.0
_P_HW = 34000.0


@_jit
def _t_mean(eta):
    if eta >= _ETA_TROPO:
        return _T0_4 * eta ** _EXPONENT4
    return _T0_4 * eta ** _EXPONENT4 + _DELTA_T * (_ETA_TROPO - eta) ** 5


@_jit
def _t_deviation(lat, eta):
    factor = eta * _PI * _U0 / _RD4
    phi_vertical = (eta - _ETA0) * 0.5 * _PI
    return factor * 1.5 * math.sin(phi_vertical) * (
        math.cos(phi_vertical)) ** 0.5 * (
        (-2.0 * math.sin(lat) ** 6 * (math.cos(lat) ** 2 + 1.0 / 3.0)
         + 10.0 / 63.0) * _U0 * math.cos(phi_vertical) ** 1.5
        + (8.0 / 5.0 * math.cos(lat) ** 3 * (math.sin(lat) ** 2 + 2.0 / 3.0)
           - _PI / 4.0) * _A_OMEGA * 0.5)


@_jit
def _surface_geopotential(lat):
    cos_tmp = _U0 * math.cos((_ETA_SFC - _ETA0) * _PI * 0.5) ** 1.5
    return ((-2.0 * math.sin(lat) ** 6 * (math.cos(lat) ** 2 + 1.0 / 3.0)
             + 10.0 / 63.0) * cos_tmp
            + (8.0 / 5.0 * math.cos(lat) ** 3
               * (math.sin(lat) ** 2 + 2.0 / 3.0)
               - _PI / 4.0) * _A_OMEGA) * cos_tmp


@_jit
def _u_wind(lon, lat, eta, perturb):
    phi_vertical = (eta - _ETA0) * 0.5 * _PI
    u = (math.cos(phi_vertical) ** 1.5 * 4.0 * _U0 * math.sin(lat) ** 2
         * math.cos(lat) ** 2)
    if perturb:
        perturb_lon = _PERT_LON * _DEG2RAD
        perturb_lat = _PERT_LAT * _DEG2RAD
        sin_tmp = math.sin(perturb_lat) * math.sin(lat)
        cos_tmp = math.cos(perturb_lat) * math.cos(lat)
        r = math.acos(sin_tmp + cos_tmp * math.cos(lon - perturb_lon))
        u_perturb = _PERT_AMP * math.exp(-(r * _RADIUS) ** 2)
        u = u_perturb + u
    return u


@_jit
def _baroclinic_wave_point(lon, lat, p, perturb, moist):
    ps = _P0_4
    eta = p / ps
    u = _u_wind(lon, lat, eta, perturb)
    v = 0.0
    t = _t_mean(eta) + _t_deviation(lat, eta)
    phis = _surface_geopotential(lat)
    if moist:
        q = (_Q0_4 * math.exp(-(lat / _LAT_HW) ** 4)
             * math.exp(-((eta - 1.0) * _P0_4 / _P_HW) ** 2))
        t = t / (1.0 + 0.608 * q)
    else:
        q = 0.0
    return u, v, t, q, ps, phis


@_jit
def _baroclinic_wave_columns(pressure, longitude, latitude, perturb, moist):
    nlev, ncol = pressure.shape
    u = np.empty((nlev, ncol))
    v = np.empty((nlev, ncol))
    t = np.empty((nlev, ncol))
    q = np.empty((nlev, ncol))
    ps = np.empty(ncol)
    phis = np.empty(ncol)
    for i in range(ncol):
        for k in range(nlev):
            (u[k, i], v[k, i], t[k, i], q[k, i], ps[i],
             phis[i]) = _baroclinic_wave_point(
                longitude[i], latitude[i], pressure[k, i], perturb, moist)
    return u, v, t, q, ps, phis


def baroclinic_wave_ics(pressure, longitude, latitude, perturb, moist):
    """pressure (nlev, ncol) Pa; longitude, latitude (ncol,) radians;
    perturb, moist bool.

    Returns (u, v, t, q, surface_pressure, surface_geopotential):
    (nlev, ncol) x4, (ncol,) x2."""
    return _baroclinic_wave_columns(
        np.ascontiguousarray(pressure, dtype=np.float64),
        np.ascontiguousarray(longitude, dtype=np.float64),
        np.ascontiguousarray(latitude, dtype=np.float64),
        bool(perturb), bool(moist))


# ---------------------------------------------------------------------------
# Test 5: tropical cyclone (Reed and Jablonowski 2011)
# ---------------------------------------------------------------------------
_A5 = 6371220.0
_RD5 = 287.0
_G5 = 9.80616
_OMEGA5 = 7.292115e-5
_CONVERT = 180.0 / _PI
_RP = 282000.0
_DP = 1115.0
_ZP = 7000.0
_Q0_5 = 0.021
_GAMMA5 = 0.007
_TS0 = 302.15
_P00 = 101500.0
_CEN_LAT = 10.0
_CEN_LON = 180.0
_ZQ1 = 3000.0
_ZQ2 = 8000.0
_EXPPR = 1.5
_EXPPZ = 2.0
_ZTROP = 15000.0
_QTROP = 1.0e-11
_RFPI = 1000000.0
_CONST_TV = 0.608
_DELTAZ = 2.0e-13
_EPSILON = 1.0e-25
_EXPONENT5 = _RD5 * _GAMMA5 / _G5
_T0_5 = _TS0 * (1.0 + _CONST_TV * _Q0_5)
_TTROP = _T0_5 - _GAMMA5 * _ZTROP


@_jit
def _fpi_f(phere, gr, zhere):
    return phere - (_P00 - _DP * math.exp(-(gr / _RP) ** _EXPPR)
                    * math.exp(-(zhere / _ZP) ** _EXPPZ)) * (
        (_T0_5 - _GAMMA5 * zhere) / _T0_5) ** (_G5 / (_RD5 * _GAMMA5))


@_jit
def _fpi_dfdz(gr, zhere):
    return (-_EXPPZ * zhere * _DP * math.exp(-(gr / _RP) ** _EXPPR)
            * math.exp(-(zhere / _ZP) ** _EXPPZ) / (_ZP * _ZP)
            * ((_T0_5 - _GAMMA5 * zhere) / _T0_5) ** (_G5 / (_RD5 * _GAMMA5))
            + _G5 / (_RD5 * _T0_5)
            * (_P00 - _DP * math.exp(-(gr / _RP) ** _EXPPR)
               * math.exp(-(zhere / _ZP) ** _EXPPZ))
            * ((_T0_5 - _GAMMA5 * zhere) / _T0_5)
            ** (_G5 / (_RD5 * _GAMMA5) - 1.0))


@_jit
def _tropical_cyclone_point(lon, lat, p):
    f = 2.0 * _OMEGA5 * math.sin(_CEN_LAT / _CONVERT)
    gr = _A5 * math.acos(
        math.sin(_CEN_LAT / _CONVERT) * math.sin(lat)
        + (math.cos(_CEN_LAT / _CONVERT) * math.cos(lat)
           * math.cos(lon - _CEN_LON / _CONVERT)))
    ps = _P00 - _DP * math.exp(-(gr / _RP) ** _EXPPR)

    height = (_T0_5 / _GAMMA5) * (1.0 - (p / ps) ** _EXPONENT5)
    if gr < _RFPI:
        # Fixed-point iteration of the Fortran's goto loop; like the
        # Fortran, it gives up (without the message) after 20 iterations.
        zhere = height
        n = 1
        while True:
            n = n + 1
            zn = zhere - _fpi_f(p, gr, zhere) / _fpi_dfdz(gr, zhere)
            if n > 20:
                break
            elif abs(zn - zhere) / abs(zn) > _DELTAZ:
                zhere = zn
            else:
                break
        height = zn

    d1 = (math.sin(_CEN_LAT / _CONVERT) * math.cos(lat)
          - math.cos(_CEN_LAT / _CONVERT) * math.sin(lat)
          * math.cos(lon - _CEN_LON / _CONVERT))
    d2 = math.cos(_CEN_LAT / _CONVERT) * math.sin(lon - _CEN_LON / _CONVERT)
    d = max(_EPSILON, math.sqrt(d1 ** 2.0 + d2 ** 2.0))
    ufac = d1 / d
    vfac = d2 / d
    if height > _ZTROP:
        u = 0.0
        v = 0.0
    else:
        speed = (-f * gr / 2.0 + math.sqrt(
            (f * gr / 2.0) ** 2.0
            - _EXPPR * (gr / _RP) ** _EXPPR * _RD5 * (_T0_5 - _GAMMA5 * height)
            / (_EXPPZ * height * _RD5 * (_T0_5 - _GAMMA5 * height)
               / (_G5 * _ZP ** _EXPPZ)
               + (1.0 - _P00 / _DP * math.exp((gr / _RP) ** _EXPPR)
                  * math.exp((height / _ZP) ** _EXPPZ)))))
        v = vfac * speed
        u = ufac * speed

    if height > _ZTROP:
        q = _QTROP
    else:
        q = _Q0_5 * math.exp(-height / _ZQ1) * math.exp(
            -(height / _ZQ2) ** _EXPPZ)

    if height > _ZTROP:
        t = _TTROP
    else:
        t = ((_T0_5 - _GAMMA5 * height) / (1.0 + _CONST_TV * q)
             / (1.0 + _EXPPZ * _RD5 * (_T0_5 - _GAMMA5 * height) * height
                / (_G5 * _ZP ** _EXPPZ
                   * (1.0 - _P00 / _DP * math.exp((gr / _RP) ** _EXPPR)
                      * math.exp((height / _ZP) ** _EXPPZ)))))
    phis = 0.0
    return u, v, t, q, ps, phis


@_jit
def _tropical_cyclone_columns(pressure, longitude, latitude):
    nlev, ncol = pressure.shape
    u = np.empty((nlev, ncol))
    v = np.empty((nlev, ncol))
    t = np.empty((nlev, ncol))
    q = np.empty((nlev, ncol))
    ps = np.empty(ncol)
    phis = np.empty(ncol)
    for i in range(ncol):
        for k in range(nlev):
            (u[k, i], v[k, i], t[k, i], q[k, i], ps[i],
             phis[i]) = _tropical_cyclone_point(
                longitude[i], latitude[i], pressure[k, i])
    return u, v, t, q, ps, phis


def tropical_cyclone_ics(pressure, longitude, latitude):
    """Same shapes and return order as :func:`baroclinic_wave_ics`; the
    tropical cyclone test is always moist."""
    return _tropical_cyclone_columns(
        np.ascontiguousarray(pressure, dtype=np.float64),
        np.ascontiguousarray(longitude, dtype=np.float64),
        np.ascontiguousarray(latitude, dtype=np.float64))
