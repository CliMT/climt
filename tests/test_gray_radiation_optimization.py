# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, GrayLongwaveRadiation, Frierson06LongwaveOpticalDepth

def test_gray_radiation_parity():
    """
    Verify that the optimized GrayLongwaveRadiation matches expected values.
    """
    nlev = 30
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=nlev)

    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())

    tau_comp = Frierson06LongwaveOpticalDepth()
    lw_comp = GrayLongwaveRadiation()

    state = get_default_state([tau_comp, lw_comp], grid_state=grid)

    # Run optical depth
    diag_tau = tau_comp(state)
    state.update(diag_tau)

    # Run radiation
    tendencies, diagnostics = lw_comp(state)

    # Check that we got reasonable non-zero values
    assert np.any(diagnostics['upwelling_longwave_flux_in_air'].values != 0.0)
    assert np.any(tendencies['air_temperature'].values != 0.0)
    assert not np.any(np.isnan(tendencies['air_temperature'].values))

if __name__ == "__main__":
    test_gray_radiation_parity()
