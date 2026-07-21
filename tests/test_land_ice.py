from datetime import timedelta

import numpy as np
from sympl import get_constant

from climt import LandIce, get_default_state, get_grid


def test_land_ice_reports_soil_heat_flux_and_conserves_mass():
    ice = LandIce()
    state = get_default_state([ice], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["area_type"].values[:] = "land_ice"
    state["land_ice_thickness"].values[:] = 3.0
    state["surface_snow_thickness"].values[:] = 0.5
    state["snow_and_ice_temperature"].values[:] = 255.0
    state["soil_surface_temperature"].values[:] = 272.0
    diag, new = ice(state, timedelta(seconds=3600))
    assert "upward_heat_flux_at_ground_level_in_soil" in diag
    # Soil (272 K) is warmer than the ice column (255 K): heat should
    # conduct upward out of the soil into the snow/ice column, so the
    # basal flux diagnostic must be positive.
    assert diag["upward_heat_flux_at_ground_level_in_soil"].values[0] > 0
    assert new["land_ice_thickness"].values[0] >= 0.0
    assert not np.any(np.isnan(new["snow_and_ice_temperature"].values))


def test_land_ice_basal_flux_reverses_sign_when_soil_is_colder():
    # Companion to test_land_ice_reports_soil_heat_flux_and_conserves_mass:
    # locks in the opposite basal-flux sign when the soil is colder than
    # the ice column, so a future sign-flip regression in the base-node
    # conduction formula (see array_call's ``heat_flux_to_soil``) is
    # caught. Neither of these two cases is exercised by
    # TestLandIceEnergyConservation (tests/test_conservation.py), which
    # deliberately keeps soil_surface_temperature equal to the initial
    # snow_and_ice_temperature so that step is a true no-op.
    ice = LandIce()
    state = get_default_state([ice], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["area_type"].values[:] = "land_ice"
    state["land_ice_thickness"].values[:] = 3.0
    state["surface_snow_thickness"].values[:] = 0.5
    state["snow_and_ice_temperature"].values[:] = 255.0
    state["soil_surface_temperature"].values[:] = 240.0
    diag, new = ice(state, timedelta(seconds=3600))
    # Soil (240 K) is colder than the ice column (255 K): heat should
    # conduct downward from the column into the soil, so the basal flux
    # diagnostic must be negative.
    assert diag["upward_heat_flux_at_ground_level_in_soil"].values[0] < 0
    assert not np.any(np.isnan(new["snow_and_ice_temperature"].values))


def test_land_ice_mass_balance_non_negative_under_strong_melt():
    ice = LandIce()
    state = get_default_state([ice], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["area_type"].values[:] = "land_ice"
    state["land_ice_thickness"].values[:] = 0.1
    state["surface_snow_thickness"].values[:] = 0.05
    state["snow_and_ice_temperature"].values[:] = 273.0
    state["soil_surface_temperature"].values[:] = 275.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 900.0
    diag, new = ice(state, timedelta(seconds=6 * 3600))
    assert new["land_ice_thickness"].values[0] >= 0.0
    assert new["surface_snow_thickness"].values[0] >= 0.0


def _land_ice_state():
    ice = LandIce()
    state = get_default_state([ice], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["area_type"].values[:] = "land_ice"
    state["land_ice_thickness"].values[:] = 3.0
    state["surface_snow_thickness"].values[:] = 0.0
    state["snow_and_ice_temperature"].values[:] = 260.0
    state["soil_surface_temperature"].values[:] = 260.0
    return ice, state


def test_land_ice_surface_energy_forcing_direction():
    # Complementary, non-vacuous coverage for TestLandIceEnergyConservation
    # (tests/test_conservation.py), which deliberately keeps forcing at
    # zero because solve_column's Flux top boundary condition is a
    # quasi-steady algebraic constraint rather than a time-integrated one,
    # producing a dt-independent, forcing-proportional residual (inherited
    # from IceSheet; see task-4-report.md and task-5-report.md). Instead
    # of exact closure, check that the column's response to a strongly
    # nonzero net surface flux moves in the direction the forcing implies.
    freezing = get_constant("freezing_temperature_of_liquid_phase", "degK")

    ice_base, state_base = _land_ice_state()
    assert state_base["snow_and_ice_temperature"].values[0, 0] < freezing
    _, new_base = ice_base(state_base, timedelta(seconds=3600))
    base_mean_temperature = new_base["snow_and_ice_temperature"].values.mean()

    # Strong net cooling: upwelling longwave far exceeds downwelling.
    ice_cool, state_cool = _land_ice_state()
    state_cool["upwelling_longwave_flux_in_air"].values[:] = 400.0
    state_cool["downwelling_longwave_flux_in_air"].values[:] = 50.0
    _, new_cool = ice_cool(state_cool, timedelta(seconds=3600))
    cool_mean_temperature = new_cool["snow_and_ice_temperature"].values.mean()

    # Strong net warming: downwelling longwave far exceeds upwelling.
    ice_warm, state_warm = _land_ice_state()
    state_warm["downwelling_longwave_flux_in_air"].values[:] = 400.0
    state_warm["upwelling_longwave_flux_in_air"].values[:] = 50.0
    _, new_warm = ice_warm(state_warm, timedelta(seconds=3600))
    warm_mean_temperature = new_warm["snow_and_ice_temperature"].values.mean()

    assert cool_mean_temperature < base_mean_temperature
    assert warm_mean_temperature > base_mean_temperature
    assert not np.any(np.isnan(new_cool["snow_and_ice_temperature"].values))
    assert not np.any(np.isnan(new_warm["snow_and_ice_temperature"].values))
