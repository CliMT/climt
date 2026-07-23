# tests/test_topography.py
"""Smoke tests for the bundled Earth orography + land-ice-thickness file.

Provides ``surface_geopotential`` (read by the GFS dynamical core) and
``land_ice_thickness`` (read by SeaIce/LandIce), both on the 2-degree grid of
the bundled land mask, derived from ETOPO 2022.
"""
import os

import numpy as np
import xarray as xr

DATA = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "topography",
    "earth_topography_2deg.nc")
MASK = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "land_mask",
    "earth_landmask_2deg.nc")
G = 9.80665


def test_bundled_topography_layout():
    # scipy engine -> must be netCDF-3 classic, like the other boundary files
    ds = xr.open_dataset(DATA, engine="scipy")
    assert set(ds.data_vars) == {"surface_geopotential", "land_ice_thickness"}
    for v in ds.data_vars:
        assert ds[v].dims == ("lat", "lon")
        assert ds[v].shape == (90, 180)
        assert np.isfinite(ds[v].values).all()
    assert ds["surface_geopotential"].attrs.get("units") == "m^2 s^-2"
    assert ds["land_ice_thickness"].attrs.get("units") == "m"


def test_topography_shares_the_mask_grid():
    ds = xr.open_dataset(DATA, engine="scipy")
    mask = xr.open_dataset(MASK, engine="scipy")
    np.testing.assert_array_equal(ds["lat"].values, mask["lat"].values)
    np.testing.assert_array_equal(ds["lon"].values, mask["lon"].values)


def test_orography_and_ice_are_physical():
    ds = xr.open_dataset(DATA, engine="scipy")
    mask = xr.open_dataset(MASK, engine="scipy")
    is_sea = mask["area_type_code"].values == 0
    orog = ds["surface_geopotential"].values / G
    ice = ds["land_ice_thickness"].values

    # non-negative; sea cells carry no orography or ice
    assert orog.min() >= 0.0 and ice.min() >= 0.0
    assert np.allclose(orog[is_sea], 0.0)
    assert np.allclose(ice[is_sea], 0.0)

    # gross-geography sanity: highlands and the two ice sheets
    def at(field, lat, lon):
        return float(ds[field].sel(lat=lat, lon=lon, method="nearest"))

    assert at("surface_geopotential", 28, 85) / G > 3000.0   # Tibet/Himalaya
    assert at("land_ice_thickness", 72, 320) > 1500.0        # Greenland
    assert at("land_ice_thickness", -80, 90) > 2000.0        # East Antarctica
    assert at("land_ice_thickness", 0, 200) == 0.0           # mid-Pacific: none
    # ice caps out near the real ~4 km Antarctic maximum, not runaway
    assert ice.max() < 5000.0
