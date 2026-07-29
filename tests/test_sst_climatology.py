# tests/test_sst_climatology.py
"""Smoke tests for the bundled 12-month climatological SST used by DataOcean."""
import datetime
import os

import numpy as np
import xarray as xr

from climt import DataOcean, LandMask, get_default_state, get_grid

DATA = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "data_ocean",
    "earth_sst_climatology_1deg.nc",
)


def test_bundled_sst_layout():
    # Must open with the scipy engine (DataOcean uses it -> netCDF-3 classic).
    ds = xr.open_dataset(DATA, engine="scipy")
    assert "tos" in ds.data_vars
    tos = ds["tos"]
    assert tos.dims == ("time", "lat", "lon")
    assert tos.shape == (12, 180, 360)
    assert tos.attrs.get("units") == "degC"
    # coords well-formed: ascending lat, 0-360 lon
    assert float(ds["lat"].min()) < -80.0 and float(ds["lat"].max()) > 80.0
    assert float(ds["lon"].min()) >= 0.0 and float(ds["lon"].max()) < 360.0
    # physical Celsius range, no NaNs in the bundled field
    assert np.isfinite(tos.values).all()
    assert -5.0 <= float(tos.min()) and float(tos.max()) <= 45.0


def test_data_ocean_consumes_bundled_sst():
    ocean = DataOcean(DATA)  # default sst_variable="tos", Celsius auto-detected
    mask = LandMask()
    state = get_default_state(
        [mask, ocean], grid_state=get_grid(nx=64, ny=32, nz=10))
    state.update(mask(state))  # set area_type before prescribing SST
    state["time"] = datetime.datetime(2000, 7, 15)
    diag = ocean(state)
    sst = np.asarray(diag["sea_surface_temperature"])
    # full sea coverage, physical Kelvin range (freezing seawater .. warm pool)
    assert np.isfinite(sst).all()
    assert 270.0 < float(sst.min()) and float(sst.max()) < 310.0


def test_climatology_loops_across_years():
    # DataOcean interpolates by day-of-year, so the field is year-agnostic and
    # loops indefinitely; seasons differ but the same month in different years
    # is identical.
    ocean = DataOcean(DATA)
    mask = LandMask()
    state = get_default_state(
        [mask, ocean], grid_state=get_grid(nx=64, ny=32, nz=10))
    state.update(mask(state))

    state["time"] = datetime.datetime(1850, 1, 15)
    jan_1850 = np.asarray(ocean(state)["sea_surface_temperature"])
    state["time"] = datetime.datetime(2400, 1, 15)
    jan_2400 = np.asarray(ocean(state)["sea_surface_temperature"])
    state["time"] = datetime.datetime(2000, 7, 15)
    jul = np.asarray(ocean(state)["sea_surface_temperature"])

    np.testing.assert_allclose(jan_1850, jan_2400)     # loops across years
    assert not np.allclose(jan_1850, jul)              # seasonal cycle present
