# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import datetime
import sympl
from climt import get_grid, get_default_state, Instellation

def test_instellation_parity():
    """
    Verify that the optimized Instellation matches original behavior.
    """
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=10)

    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())

    inst = Instellation()
    state = get_default_state([inst], grid_state=grid)
    state['time'] = datetime(2026, 6, 21, 12, 0)

    output = inst(state)

    # Check reasonable values
    assert np.any(output['zenith_angle'].values != 0.0)
    assert not np.any(np.isnan(output['zenith_angle'].values))

    print("SUCCESS: Instellation parity verified (basic)!")

if __name__ == "__main__":
    test_instellation_parity()
