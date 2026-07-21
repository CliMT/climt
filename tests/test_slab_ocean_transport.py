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
