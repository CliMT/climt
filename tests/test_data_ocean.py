import numpy as np
import pytest
import xarray as xr
from datetime import datetime
from climt._components.data_ocean._sst_interpolation import (
    mid_month_values, interp_time,
)
from climt import get_default_state, get_grid


def test_mid_month_preserves_monthly_means():
    rng = np.random.default_rng(1)
    means = 285.0 + rng.normal(0, 3, size=12)
    mm = mid_month_values(means)
    # Independent geometric reconstruction: month m spans from the midpoint of
    # (mm[m-1],mm[m]) to the midpoint of (mm[m],mm[m+1]); average the two-segment
    # piecewise-linear curve (boundary -> peak mm[m] -> boundary) over the month.
    recon = np.empty(12)
    for m in range(12):
        a = 0.5 * (mm[(m - 1) % 12] + mm[m])   # start-of-month boundary value
        b = 0.5 * (mm[m] + mm[(m + 1) % 12])   # end-of-month boundary value
        recon[m] = 0.5 * (0.5 * (a + mm[m]) + 0.5 * (mm[m] + b))
    assert np.allclose(recon, means, atol=1e-6)


def test_mid_month_preserves_means_numerically():
    rng = np.random.default_rng(2)
    means = 285.0 + rng.normal(0, 3, size=12)
    mm = mid_month_values(means)

    # Independent numerical-integration check: sample the piecewise-linear
    # interpolant through mid-month values (mid-months at fractions
    # (m+0.5)/12 of the year) on a fine uniform grid, and numerically
    # average it over each equal month-interval [m/12, (m+1)/12].
    n = 12
    mid_fracs = (np.arange(n) + 0.5) / n
    # Extend cyclically so interpolation covers the full [0, 1] range.
    ext_fracs = np.concatenate(
        [[mid_fracs[-1] - 1.0], mid_fracs, [mid_fracs[0] + 1.0]]
    )
    ext_vals = np.concatenate([[mm[-1]], mm, [mm[0]]])

    n_samples = 20000
    x = np.linspace(0.0, 1.0, n_samples, endpoint=False) + 0.5 / n_samples
    interp_vals = np.interp(x, ext_fracs, ext_vals)

    recon = np.empty(n)
    for m in range(n):
        lo, hi = m / n, (m + 1) / n
        mask = (x >= lo) & (x < hi)
        recon[m] = interp_vals[mask].mean()

    assert np.allclose(recon, means, atol=1e-2)


def test_interp_time_within_bounds():
    means = np.linspace(280, 300, 12)
    mm = mid_month_values(means)
    for month in range(1, 13):
        val = interp_time(mm, datetime(2000, month, 15, 12))
        assert mm.min() - 1e-6 <= val <= mm.max() + 1e-6


def _make_sst_file(tmp_path):
    lat = np.arange(-88.0, 90.0, 4.0)
    lon = np.arange(2.0, 360.0, 4.0)
    time = np.arange(12)  # month index
    # SST = 300 - (lat/3)^2 scaled, constant in time for a simple check
    base = 300.0 - (lat[:, None] / 3.0) ** 2 * 0.01 + 0 * lon[None, :]
    data = np.repeat(base[None, :, :], 12, axis=0)
    ds = xr.Dataset(
        {"tos": (("time", "lat", "lon"), data)},
        coords={"time": time, "lat": lat, "lon": lon},
    )
    ds["tos"].attrs["units"] = "K"
    path = str(tmp_path / "sst.nc")
    ds.to_netcdf(path, engine="scipy")
    return path


def test_data_ocean_prescribes_sea_only(tmp_path):
    from climt import DataOcean
    ocean = DataOcean(_make_sst_file(tmp_path), sst_variable="tos")
    state = get_default_state([ocean], grid_state=get_grid(nx=32, ny=16, nz=10))
    lat = state["latitude"].values
    state["area_type"].values[:] = "sea"
    state["area_type"].values[lat > 45.0] = "land"
    land_before = state["surface_temperature"].values[lat > 45.0].copy()
    state["time"] = datetime(2000, 7, 15, 12)
    diag = ocean(state)
    sst = diag["sea_surface_temperature"].values
    # every sea cell finite
    assert np.all(np.isfinite(sst[lat <= 45.0]))
    # land cells passed through unchanged in surface_temperature
    assert np.allclose(diag["surface_temperature"].values[lat > 45.0], land_before)


@pytest.mark.parametrize("nx,ny", [(8, 4), (64, 40)])
def test_data_ocean_arbitrary_resolution(tmp_path, nx, ny):
    from climt import DataOcean
    ocean = DataOcean(_make_sst_file(tmp_path), sst_variable="tos")
    state = get_default_state([ocean], grid_state=get_grid(nx=nx, ny=ny, nz=10))
    state["area_type"].values[:] = "sea"
    state["time"] = datetime(2000, 1, 15, 12)
    sst = ocean(state)["sea_surface_temperature"].values
    assert np.all(np.isfinite(sst))
    assert sst.min() > 200.0 and sst.max() < 350.0  # bounded/physical


def test_data_ocean_full_sea_coverage_offset_mask(tmp_path):
    from climt import DataOcean
    # source with a land hole; model sea cells over the hole must still be filled
    lat = np.arange(-88.0, 90.0, 4.0)
    lon = np.arange(2.0, 360.0, 4.0)
    data = np.full((12, lat.size, lon.size), 290.0)
    data[:, lat.size // 2, lon.size // 2] = np.nan  # masked source point
    ds = xr.Dataset({"tos": (("time", "lat", "lon"), data)},
                    coords={"time": np.arange(12), "lat": lat, "lon": lon})
    ds["tos"].attrs["units"] = "K"
    path = str(tmp_path / "holey.nc")
    ds.to_netcdf(path, engine="scipy")
    ocean = DataOcean(path, sst_variable="tos")
    state = get_default_state([ocean], grid_state=get_grid(nx=90, ny=45, nz=10))
    state["area_type"].values[:] = "sea"
    state["time"] = datetime(2000, 6, 15, 12)
    sst = ocean(state)["sea_surface_temperature"].values
    assert np.all(np.isfinite(sst))   # nearest-valid fill covered the hole
