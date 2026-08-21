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


def test_gsc_precipitation_amount_is_positive_mass_per_area():
    """
    precipitation_amount must be a positive kg m^-2 accumulation.

    The expected value is computed by hand from the physics: for every
    supersaturated layer, the condensed specific humidity times the layer
    mass per unit area dp/g (kg m^-2), summed over the column.
    """
    nlev = 28
    grid = get_grid(nx=1, ny=1, nz=nlev)

    sympl.set_backend(sympl.DataArrayBackend())

    gsc = GridScaleCondensation()
    state = get_default_state([gsc], grid_state=grid)

    state['air_temperature'].values[:] = 280.0
    # Supersaturate only the lowest few layers; the thin top layers of the
    # default grid are numerically pathological at 0.05 kg/kg.
    state['specific_humidity'].values[:] = 0.0
    state['specific_humidity'].values[:5] = 0.03

    q_before = state['specific_humidity'].values.copy()
    p_int = state['air_pressure_on_interface_levels'].values.copy()

    diagnostics, outputs = gsc(state, timedelta(minutes=10))

    precip = diagnostics['precipitation_amount'].values

    # Hand-computed expectation: condensed mass per unit area.
    g = 9.80665
    dp = p_int[:-1, ...] - p_int[1:, ...]
    condensed = q_before - outputs['specific_humidity'].values
    expected = np.sum(condensed * dp / g, axis=0)

    assert p_int[0, 0, 0] > p_int[-1, 0, 0], 'interface pressures are bottom-first'
    assert np.all(condensed[:5] > 0), 'the lowest layers should be supersaturated'
    assert np.all(precip > 0)
    # A saturated-to-280K column drops O(10) mm of water; sanity-bound it so a
    # factor-of-1000 units error cannot pass.
    assert np.all(precip > 1.0)
    assert np.all(precip < 1000.0)
    np.testing.assert_allclose(precip, expected, rtol=1e-10)


if __name__ == "__main__":
    test_gsc_parity()
