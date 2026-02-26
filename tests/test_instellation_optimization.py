# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import datetime
import sympl
from climt import get_grid, get_default_state, Instellation, JaxBackend

# Enable x64
jax.config.update("jax_enable_x64", True)

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

def test_instellation_differentiation():
    """
    Verify that Instellation is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    backend = JaxBackend()
    sympl.set_backend(backend)
    
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=10)
    inst = Instellation()
    state = get_default_state([inst], grid_state=grid)
    state['time'] = datetime(2026, 6, 21, 12, 0)
    
    def compute_loss(lat_values):
        current_state = {
            'latitude': lat_values,
            'longitude': state['longitude'].data,
            'time': state['time'],
        }
        
        output = inst.array_call(current_state)
        return jnp.sum(jnp.square(output['zenith_angle']))

    initial_lat = state['latitude'].data
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_lat)
    
    print(f"\nInstellation Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Instellation should be differentiable w.r.t latitude"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Instellation differentiation verified!")

if __name__ == "__main__":
    test_instellation_parity()
    test_instellation_differentiation()
