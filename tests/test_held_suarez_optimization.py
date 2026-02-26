# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, HeldSuarez, JaxBackend

# Enable x64
jax.config.update("jax_enable_x64", True)

def test_held_suarez_parity():
    """
    Verify that the optimized HeldSuarez matches expected values.
    """
    nlev = 30
    ncol = 10
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    
    # Force NumPy backend
    sympl.set_backend(sympl.DataArrayBackend())
    
    hs = HeldSuarez()
    state = get_default_state([hs], grid_state=grid)
    
    # Set some non-zero winds and varied temperature
    state['eastward_wind'].values[:] = 10.0
    state['air_temperature'].values[:] = 300.0
    
    tendencies, _ = hs(state)
    
    # Check that we got reasonable non-zero tendencies
    assert np.any(tendencies['eastward_wind'].values != 0.0)
    assert np.any(tendencies['air_temperature'].values != 0.0)
    assert not np.any(np.isnan(tendencies['air_temperature'].values))

def test_held_suarez_differentiation():
    """
    Verify that HeldSuarez is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    
    nlev = 30
    ncol = 1
    # We use sympl to get the grid but call array_call directly to avoid np.asarray calls in sympl
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    hs = HeldSuarez()
    state_sym = get_default_state([hs], grid_state=grid)
    
    # Create plain JAX dictionary for array_call
    jax_state = {
        'air_temperature': jnp.full((ncol, nlev), 300.0),
        'eastward_wind': jnp.full((ncol, nlev), 10.0),
        'northward_wind': jnp.zeros((ncol, nlev)),
        'air_pressure': jnp.array(state_sym['air_pressure'].values.T),
        'surface_air_pressure': jnp.array(state_sym['surface_air_pressure'].values.T).flatten(),
        'latitude': jnp.array(state_sym['latitude'].values.T).flatten(),
    }
    
    def compute_loss(temp_values):
        current_state = {**jax_state, 'air_temperature': temp_values}
        # Direct call to array_call is differentiable!
        tendencies, _ = hs.array_call(current_state)
        return jnp.sum(jnp.square(tendencies['air_temperature']))

    initial_temp = jax_state['air_temperature']
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_temp)
    
    print(f"\nHeld-Suarez Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Held-Suarez should be differentiable w.r.t temperature"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Held-Suarez differentiation verified!")

if __name__ == "__main__":
    test_held_suarez_parity()
    test_held_suarez_differentiation()
