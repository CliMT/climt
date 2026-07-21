# tests/test_land_mask.py
import os
import numpy as np
import xarray as xr

DATA = os.path.join(
    os.path.dirname(__file__), "..", "climt", "_data", "land_mask",
    "earth_landmask_2deg.nc",
)

def test_bundled_mask_gross_geography():
    ds = xr.open_dataset(DATA, engine="scipy")
    code = ds["area_type_code"]
    # mid-Pacific is sea (0)
    assert int(code.sel(lat=0.0, lon=200.0, method="nearest")) == 0
    # central Asia is land (1)
    assert int(code.sel(lat=45.0, lon=90.0, method="nearest")) == 1
    # central Antarctica is land_ice (2)
    assert int(code.sel(lat=-85.0, lon=100.0, method="nearest")) == 2
    # only valid codes present
    assert set(np.unique(code.values)).issubset({0, 1, 2})
