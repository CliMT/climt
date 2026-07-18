import numpy as np
from datetime import timedelta
from climt import IceSheet, get_default_state


def _run_once():
    ice = IceSheet()
    state = get_default_state([ice])
    # ensure a non-trivial ice column so calculate_new_ice_temperature runs
    state["area_type"].values[:] = "sea_ice"
    state["sea_ice_thickness"].values[:] = 2.0
    state["snow_and_ice_temperature"].values[:] = 260.0
    diag, new = ice(state, timedelta(hours=1))
    return new["snow_and_ice_temperature"].values.copy()


def test_icesheet_matches_golden():
    golden = np.load("tests/_golden/icesheet_snow_ice_temperature.npy")
    np.testing.assert_allclose(_run_once(), golden, rtol=1e-9, atol=1e-9)
