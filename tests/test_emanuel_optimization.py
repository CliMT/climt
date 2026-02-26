# -*- coding: utf-8 -*-
import numpy as np
import pytest
from datetime import timedelta
from climt import get_grid, get_default_state, EmanuelConvection
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v2 import EmanuelConvectionPythonV2
from climt._core.util import numpy_version_of

def create_test_state(nlev, ncol, moisture_type='moist'):
    """Create a variety of atmospheric profiles for testing."""
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    
    # Base temperature profile (isothermal)
    state['air_temperature'].values[:] = 300.0
    
    q = np.zeros_like(state['specific_humidity'].values)
    if moisture_type == 'moist':
        # High surface moisture to force convection
        q[0:5, :, :] = 0.020
    elif moisture_type == 'dry':
        # Very dry profile
        q[:] = 1e-6
    elif moisture_type == 'unstable':
        # Lapse rate that encourages instability
        for i in range(nlev):
            state['air_temperature'].values[i, :, :] = 310.0 - (i * 2.0)
        q[0:10, :, :] = 0.015
        
    state['specific_humidity'].values[:] = q
    return state

@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("moisture_type", ['moist', 'dry', 'unstable'])
def test_emanuel_v2_parity(ncol, moisture_type):
    """
    Formally test parity between Original Python and V2 Abstraction.
    """
    nlev = 30
    state = create_test_state(nlev, ncol, moisture_type)
    timestep = timedelta(minutes=15)
    
    raw_state = numpy_version_of(state)
    
    # Prepare state for array_call (nlev, ncol)
    def prepare_input(arr):
        # arr is (nlev, ny, nx) from get_default_state + get_grid
        # We want (nlev, ncol) where ncol = ny*nx
        shape = arr.shape
        return arr.reshape(shape[0], -1)

    python_state = {
        'air_temperature': prepare_input(raw_state['air_temperature']),
        'specific_humidity': prepare_input(raw_state['specific_humidity']),
        'air_pressure': prepare_input(raw_state['air_pressure'] / 100.0),
        'air_pressure_on_interface_levels': prepare_input(raw_state['air_pressure_on_interface_levels'] / 100.0),
        'eastward_wind': prepare_input(raw_state['eastward_wind']),
        'northward_wind': prepare_input(raw_state['northward_wind']),
    }

    conv_orig = EmanuelConvectionPython()
    conv_v2 = EmanuelConvectionPythonV2()

    # Get results
    tend_orig, diag_orig = conv_orig.array_call(python_state, timestep)
    tend_v2, diag_v2 = conv_v2.array_call(python_state, timestep)

    # Check Tendencies
    for key in tend_orig:
        np.testing.assert_allclose(
            tend_orig[key], 
            tend_v2[key], 
            atol=1e-12, 
            rtol=1e-12,
            err_msg=f"Tendency parity failed for {key} (ncol={ncol}, type={moisture_type})"
        )

    # Check Diagnostics
    for key in diag_orig:
        # iflag (convective_state) might be slightly different in type but should match values
        np.testing.assert_allclose(
            diag_orig[key], 
            diag_v2[key], 
            atol=1e-12, 
            rtol=1e-12,
            err_msg=f"Diagnostic parity failed for {key} (ncol={ncol}, type={moisture_type})"
        )
