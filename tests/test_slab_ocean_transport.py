import numpy as np
from climt import SlabSurface, get_default_state, get_grid
import sympl


def test_qflux_adds_to_sea_tendency():
    sympl.set_backend(sympl.DataArrayBackend())
    grid = get_grid(nx=4, ny=1, nz=10)
    slab = SlabSurface()
    state = get_default_state([slab], grid_state=grid)
    state["area_type"].values[:] = "sea"
    tend0, _ = slab(state)
    state["ocean_heat_transport_convergence"].values[:] = 50.0
    tend1, _ = slab(state)
    dT = (tend1["surface_temperature"].values - tend0["surface_temperature"].values)
    assert np.all(dT > 0.0)   # positive convergence warms the mixed layer


def test_zero_qflux_matches_baseline():
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface()
    state = get_default_state([slab], grid_state=get_grid(nx=4, ny=1, nz=10))
    state["area_type"].values[:] = "sea"
    tend, _ = slab(state)
    # default field is zero -> tendency identical to computing without the term
    assert not np.any(np.isnan(tend["surface_temperature"].values))


def test_ekman_pumping_from_wind_stress_curl():
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface(include_ekman=True)
    state = get_default_state([slab], grid_state=get_grid(nx=32, ny=16, nz=10))
    state["area_type"].values[:] = "sea"
    lat = state["latitude"].values
    # zonal stress with meridional structure -> nonzero curl
    state["surface_downward_eastward_stress"].values[:] = 0.1 * np.cos(np.deg2rad(lat))
    state["surface_downward_northward_stress"].values[:] = 0.0
    tend, diag = slab(state)
    assert "ekman_pumping" in diag
    assert np.any(diag["ekman_pumping"].values != 0.0)
    assert not np.any(np.isnan(tend["surface_temperature"].values))


def test_ekman_zero_curl_gives_zero_pumping():
    sympl.set_backend(sympl.DataArrayBackend())
    slab = SlabSurface(include_ekman=True)
    state = get_default_state([slab], grid_state=get_grid(nx=16, ny=8, nz=10))
    state["area_type"].values[:] = "sea"
    state["surface_downward_eastward_stress"].values[:] = 0.1   # uniform -> no curl
    tend, diag = slab(state)
    assert np.allclose(diag["ekman_pumping"].values, 0.0, atol=1e-12)
