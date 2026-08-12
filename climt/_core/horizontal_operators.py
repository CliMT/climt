"""Spherical horizontal differential operators on a lat/lon grid (degrees).

Fields and coordinates are 2-D arrays shaped (nlat, nlon). Returns zeros
when a dimension is too small to difference (single column / degenerate).
"""
import numpy as np
from sympl import get_constant

_RADIUS = None


def _earth_radius():
    global _RADIUS
    if _RADIUS is None:
        try:
            _RADIUS = get_constant("planetary_radius", "m")
        except Exception:
            _RADIUS = 6.371e6
    return _RADIUS


def _grads(field, lat, lon):
    a = _earth_radius()
    latr = np.deg2rad(lat)
    lonr = np.deg2rad(lon)
    if field.shape[0] < 3 or field.shape[1] < 3:
        z = np.zeros_like(field, dtype=float)
        return z, z
    dfdlat = np.gradient(field, axis=0) / np.gradient(latr, axis=0)
    dfdlon = np.gradient(field, axis=1) / np.gradient(lonr, axis=1)
    df_dy = dfdlat / a
    df_dx = dfdlon / (a * np.cos(latr))
    return df_dx, df_dy


def divergence(fx, fy, lat, lon):
    dfx_dx, _ = _grads(fx, lat, lon)
    _, dfy_dy = _grads(fy, lat, lon)
    return dfx_dx + dfy_dy


def curl_z(fx, fy, lat, lon):
    dfy_dx, _ = _grads(fy, lat, lon)
    _, dfx_dy = _grads(fx, lat, lon)
    return dfy_dx - dfx_dy
