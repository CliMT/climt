# tests/test_land_mask.py
import os
import numpy as np
import pytest
import xarray as xr
import climt
from climt import get_grid, get_default_state

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


def _synthetic_mask(tmp_path):
    lat = np.arange(-89.0, 90.0, 2.0)
    lon = np.arange(1.0, 360.0, 2.0)
    code = np.zeros((lat.size, lon.size), dtype=np.int8)  # all sea
    code[lat > 0.0, :] = 1                                # northern hemi land
    da = xr.DataArray(code, coords={"lat": lat, "lon": lon},
                      dims=("lat", "lon"), name="area_type_code")
    da.attrs["flag_meanings"] = "sea land land_ice"
    path = str(tmp_path / "synth_mask.nc")
    da.to_dataset().to_netcdf(path, engine="scipy")
    return path


def test_land_mask_sets_area_type_from_geography(tmp_path):
    from climt import LandMask
    mask = LandMask(mask_dataset=_synthetic_mask(tmp_path))
    state = get_default_state([mask], grid_state=get_grid(nx=32, ny=16, nz=10))
    diag = mask(state)
    at = diag["area_type"].values.astype(str)
    lat = state["latitude"].values
    # northern cells -> land, southern -> sea
    assert np.all(at[lat > 0.0] == "land")
    assert np.all(at[lat < 0.0] == "sea")
    # only valid categories
    assert set(np.unique(at)).issubset({"sea", "land", "land_ice"})


@pytest.mark.parametrize("nx,ny", [(8, 4), (64, 32), (144, 90)])
def test_land_mask_arbitrary_resolution_all_assigned(nx, ny):
    from climt import LandMask
    mask = LandMask()  # bundled default
    state = get_default_state([mask], grid_state=get_grid(nx=nx, ny=ny, nz=10))
    at = mask(state)["area_type"].values.astype(str)
    assert at.size == state["latitude"].values.size
    assert np.all(np.isin(at, ["sea", "land", "land_ice"]))  # every cell valid
    assert not np.any(at == "")                              # none left blank


G = 9.80665


def test_land_mask_loads_topography_by_default():
    from climt import LandMask
    mask = LandMask()  # load_topography=True by default
    assert set(mask.diagnostic_properties) == {
        "area_type", "surface_geopotential", "land_ice_thickness"}
    state = get_default_state(
        [mask], grid_state=get_grid(nx=128, ny=62, nz=10))
    diag = mask(state)
    area = diag["area_type"].values.astype(str)
    orog = diag["surface_geopotential"].values / G
    ice = diag["land_ice_thickness"].values

    # non-negative everywhere, finite, and flat/ice-free over ocean
    assert np.isfinite(orog).all() and np.isfinite(ice).all()
    assert orog.min() >= 0.0 and ice.min() >= 0.0
    assert np.allclose(orog[area == "sea"], 0.0)
    assert np.allclose(ice[area == "sea"], 0.0)
    # gross geography: real highlands and the two ice sheets show up
    assert orog.max() > 3000.0          # Himalaya/Antarctic plateau
    assert 2000.0 < ice.max() < 5000.0  # Antarctica ~4 km, not runaway
    assert np.any(ice[area == "land_ice"] > 1000.0)


def test_land_mask_topography_can_be_disabled():
    from climt import LandMask
    mask = LandMask(load_topography=False)
    assert set(mask.diagnostic_properties) == {"area_type"}
    state = get_default_state(
        [mask], grid_state=get_grid(nx=32, ny=16, nz=10))
    diag = mask(state)
    assert "surface_geopotential" not in diag
    assert "land_ice_thickness" not in diag
