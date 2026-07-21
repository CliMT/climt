from datetime import timedelta

import numpy as np

from climt import SeaIce, get_default_state, get_grid


def _sea_ice_state():
    ice = SeaIce()
    state = get_default_state([ice], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["area_type"].values[:] = "sea_ice"
    state["sea_ice_thickness"].values[:] = 1.0
    state["surface_snow_thickness"].values[:] = 0.0
    state["snow_and_ice_temperature"].values[:] = 260.0
    return ice, state


def test_sea_ice_grows_when_cold():
    ice, state = _sea_ice_state()
    # strong upward (cooling) surface flux, no ocean heat flux
    state["downwelling_longwave_flux_in_air"].values[:] = 100.0
    state["upwelling_longwave_flux_in_air"].values[:] = 320.0
    diag, new = ice(state, timedelta(seconds=3600))
    assert new["sea_ice_thickness"].values[0] >= state["sea_ice_thickness"].values[0]
    assert not np.any(np.isnan(new["snow_and_ice_temperature"].values))


def test_sea_ice_thickness_never_negative_under_strong_melt():
    ice, state = _sea_ice_state()
    state["sea_ice_thickness"].values[:] = 0.05
    state["snow_and_ice_temperature"].values[:] = 273.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 800.0
    state["downwelling_longwave_flux_in_air"].values[:] = 400.0
    diag, new = ice(state, timedelta(seconds=6 * 3600))
    assert new["sea_ice_thickness"].values[0] >= 0.0  # defect 6 fix


def test_sea_ice_basal_ocean_heat_flux_direction():
    # Regression coverage for the basal Flux boundary condition (Task 3's
    # bottom-Flux code path was previously untested). Use near-neutral top
    # forcing so the basal ocean heat flux dominates the thickness response.
    #
    # Sign convention (verified numerically against solve_column directly,
    # see task-4-report.md): "heat_flux_into_sea_water_due_to_sea_ice" is
    # positive when heat leaves the ice and enters the ocean below (which
    # this simple Stefan-style bookkeeping treats as basal melt, shrinking
    # sea_ice_thickness), and negative when heat is drawn out of the ocean
    # into the ice base (basal freezing, growing sea_ice_thickness). A
    # positive stored value should therefore thin the ice relative to a
    # zero-flux baseline, and a negative stored value should thicken it.
    ice_base, state_base = _sea_ice_state()
    state_base["downwelling_longwave_flux_in_air"].values[:] = 200.0
    state_base["upwelling_longwave_flux_in_air"].values[:] = 200.0
    state_base["heat_flux_into_sea_water_due_to_sea_ice"].values[:] = 0.0
    _, new_base = ice_base(state_base, timedelta(seconds=3600))

    ice_melt, state_melt = _sea_ice_state()
    state_melt["downwelling_longwave_flux_in_air"].values[:] = 200.0
    state_melt["upwelling_longwave_flux_in_air"].values[:] = 200.0
    state_melt["heat_flux_into_sea_water_due_to_sea_ice"].values[:] = 500.0
    _, new_melt = ice_melt(state_melt, timedelta(seconds=3600))

    ice_grow, state_grow = _sea_ice_state()
    state_grow["downwelling_longwave_flux_in_air"].values[:] = 200.0
    state_grow["upwelling_longwave_flux_in_air"].values[:] = 200.0
    state_grow["heat_flux_into_sea_water_due_to_sea_ice"].values[:] = -500.0
    _, new_grow = ice_grow(state_grow, timedelta(seconds=3600))

    assert new_melt["sea_ice_thickness"].values[0] < new_base["sea_ice_thickness"].values[0]
    assert new_grow["sea_ice_thickness"].values[0] > new_base["sea_ice_thickness"].values[0]
    assert not np.any(np.isnan(new_melt["snow_and_ice_temperature"].values))
    assert not np.any(np.isnan(new_grow["snow_and_ice_temperature"].values))
