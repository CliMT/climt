import warnings
from datetime import timedelta

import numpy as np

from climt import IceSheet, get_default_state, get_grid


def test_icesheet_warns_and_dispatches():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ice = IceSheet()
        assert any(issubclass(x.category, DeprecationWarning) for x in w)
    state = get_default_state([ice], grid_state=get_grid(nx=2, ny=1, nz=10))
    # get_default_state returns arrays shaped (lat, lon) = (1, 2) here, not
    # a flat length-2 array, so columns are indexed as [0, 0] / [0, 1].
    state["area_type"].values[0, 0] = "sea_ice"
    state["area_type"].values[0, 1] = "land_ice"
    state["sea_ice_thickness"].values[0, 0] = 1.0
    state["land_ice_thickness"].values[0, 1] = 2.0
    diag, new = ice(state, timedelta(seconds=600))
    assert not np.any(np.isnan(new["snow_and_ice_temperature"].values))
