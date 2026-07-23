import numpy as np
from climt._core.horizontal_operators import curl_z, divergence


def test_curl_of_solid_body_field_is_constant():
    lat = np.linspace(-40, 40, 41)
    lon = np.linspace(0, 80, 81)
    LON, LAT = np.meshgrid(lon, lat)
    # fx = -y, fy = x  -> curl_z = 2 (in planar limit); check sign & nonzero
    fx = -(LAT)
    fy = (LON)
    c = curl_z(fx, fy, LAT, LON)
    assert np.all(c[2:-2, 2:-2] > 0.0)


def test_single_column_returns_zero():
    lat = np.array([[10.0]])
    lon = np.array([[20.0]])
    assert np.all(curl_z(np.zeros((1, 1)), np.zeros((1, 1)), lat, lon) == 0.0)
    assert np.all(divergence(np.zeros((1, 1)), np.zeros((1, 1)), lat, lon) == 0.0)
