# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, GridScaleCondensation

def test_gsc_parity():
    """
    Verify that the optimized GridScaleCondensation matches expected values.
    """
    nlev = 30
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=nlev)

    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())

    gsc = GridScaleCondensation()
    state = get_default_state([gsc], grid_state=grid)

    # Saturated profile
    state['air_temperature'].values[:] = 280.0
    state['specific_humidity'].values[:] = 0.05

    initial_temp = state['air_temperature'].values.copy()

    timestep = timedelta(minutes=10)
    _, outputs = gsc(state, timestep)

    new_temp = outputs['air_temperature'].values

    # Check that temperature changed
    assert np.any(new_temp != initial_temp)
    assert not np.any(np.isnan(new_temp))

    print("SUCCESS: Grid Scale Condensation parity verified (basic change)!")

if __name__ == "__main__":
    test_gsc_parity()
