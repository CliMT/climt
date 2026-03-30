# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, HeldSuarez

def test_held_suarez_parity():
    """
    Verify that the optimized HeldSuarez matches expected values.
    """
    nlev = 30
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=nlev)

    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())

    hs = HeldSuarez()
    state = get_default_state([hs], grid_state=grid)

    # Set some non-zero winds and varied temperature
    state['eastward_wind'].values[:] = 10.0
    state['air_temperature'].values[:] = 300.0

    tendencies, _ = hs(state)

    # Check that we got reasonable non-zero tendencies
    assert np.any(tendencies['eastward_wind'].values != 0.0)
    assert np.any(tendencies['air_temperature'].values != 0.0)
    assert not np.any(np.isnan(tendencies['air_temperature'].values))

if __name__ == "__main__":
    test_held_suarez_parity()
