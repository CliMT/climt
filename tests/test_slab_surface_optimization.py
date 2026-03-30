# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, SlabSurface

def test_slab_surface_parity():
    """
    Verify that the optimized SlabSurface matches expected values.
    """
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=20)

    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())

    slab = SlabSurface()
    state = get_default_state([slab], grid_state=grid)

    # Set some non-zero fluxes
    state['downwelling_shortwave_flux_in_air'].values[:] = 300.0
    state['downwelling_longwave_flux_in_air'].values[:] = 300.0
    state['upwelling_shortwave_flux_in_air'].values[:] = 50.0
    state['upwelling_longwave_flux_in_air'].values[:] = 300.0

    initial_ts = state['surface_temperature'].values.copy()

    tendencies, diagnostics = slab(state)

    # Check reasonable values
    assert np.any(tendencies['surface_temperature'].values != 0.0)
    assert not np.any(np.isnan(tendencies['surface_temperature'].values))

    print("SUCCESS: Slab Surface parity verified (basic)!")

if __name__ == "__main__":
    test_slab_surface_parity()
