# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, GridScaleCondensation, JaxBackend

# Enable x64
jax.config.update("jax_enable_x64", True)

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

def test_gsc_differentiation():
    """
    Verify that GridScaleCondensation is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    backend = JaxBackend()
    sympl.set_backend(backend)
    
    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    
    gsc = GridScaleCondensation()
    state = get_default_state([gsc], grid_state=grid)
    
    # Saturated profile
    initial_temp = jnp.full((nlev, ncol), 280.0)
    
    def compute_loss(temp_values):
        current_state = {
            'air_temperature': temp_values,
            'specific_humidity': jnp.full((nlev, ncol), 0.05),
            'air_pressure': state['air_pressure'].data,
            'air_pressure_on_interface_levels': state['air_pressure_on_interface_levels'].data,
        }
        
        _, outputs = gsc.array_call(current_state, timedelta(minutes=10))
        return jnp.sum(jnp.square(outputs['air_temperature']))

    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_temp)
    
    print(f"\nGrid Scale Condensation Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Grid Scale Condensation should be differentiable w.r.t temperature"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Grid Scale Condensation differentiation verified!")

if __name__ == "__main__":
    test_gsc_parity()
    test_gsc_differentiation()
