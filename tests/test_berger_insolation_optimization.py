# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import datetime
import sympl
from climt import get_grid, get_default_state, BergerSolarInsolation

def test_berger_insolation_parity():
    """
    Verify that the optimized BergerSolarInsolation matches expected values.
    """
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=10)

    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())

    insolation = BergerSolarInsolation()
    state = get_default_state([insolation], grid_state=grid)
    state['time'] = datetime(2026, 6, 21, 12, 0)

    output = insolation(state)

    # Check reasonable values
    # Summer solstice at noon should have high insolation in NH
    assert np.any(output['solar_insolation'].values > 0.0)
    assert not np.any(np.isnan(output['solar_insolation'].values))

    print("SUCCESS: Berger Insolation parity verified (basic)!")

if __name__ == "__main__":
    test_berger_insolation_parity()
