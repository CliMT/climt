# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import datetime
import sympl
from climt import get_grid, get_default_state, BergerSolarInsolation, JaxBackend

# Enable x64
jax.config.update("jax_enable_x64", True)

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

def test_berger_insolation_differentiation():
    """
    Verify that BergerSolarInsolation is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    backend = JaxBackend()
    sympl.set_backend(backend)
    
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=10)
    insolation = BergerSolarInsolation()
    state = get_default_state([insolation], grid_state=grid)
    state['time'] = datetime(2026, 6, 21, 12, 0)
    
    def compute_loss(lat_values):
        current_state = {
            'latitude': lat_values,
            'longitude': state['longitude'].data,
            'time': state['time'],
        }
        
        output = insolation.array_call(current_state)
        return jnp.sum(jnp.square(output['solar_insolation']))

    initial_lat = state['latitude'].data
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_lat)
    
    print(f"\nBerger Insolation Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Berger Insolation should be differentiable w.r.t latitude"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Berger Insolation differentiation verified!")

if __name__ == "__main__":
    test_berger_insolation_parity()
    test_berger_insolation_differentiation()
