# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, EmanuelConvectionPythonV3, JaxBackend

# 1. Enable x64 for scientific precision
jax.config.update("jax_enable_x64", True)

def test_emanuel_jax_differentiation_direct():
    """
    Test that we can compute gradients through the array_call interface
    of EmanuelConvectionPythonV3. 
    """
    # Force CPU
    jax.config.update('jax_platform_name', 'cpu')
    
    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    conv = EmanuelConvectionPythonV3()
    timestep = timedelta(minutes=15)
    state_sympl = get_default_state([conv], grid_state=grid)
    
    # We need a profile that has significant SENSITIVITY to temperature
    # A tiny change in T should change the buoyancy and thus the tendencies.
    
    # Create an unstable profile
    temp = jnp.linspace(310, 250, nlev).reshape(nlev, ncol)
    q = jnp.zeros((nlev, ncol)).at[0:10].set(0.015)
    
    jax_state = {
        'air_temperature': temp,
        'specific_humidity': q,
        'air_pressure': jnp.array(state_sympl['air_pressure'].values / 100.0).reshape(nlev, ncol),
        'air_pressure_on_interface_levels': jnp.array(state_sympl['air_pressure_on_interface_levels'].values / 100.0).reshape(nlev+1, ncol),
        'eastward_wind': jnp.zeros((nlev, ncol)),
        'northward_wind': jnp.zeros((nlev, ncol)),
    }

    # 3. Define a differentiable function through array_call
    def compute_loss(temp_values):
        current_state = {**jax_state, 'air_temperature': temp_values}
        tendencies, diagnostics = conv.array_call(current_state, timestep)
        # Sum of squares of tendencies to ensure we pick up ANY change
        return jnp.sum(jnp.square(tendencies['air_temperature']))

    # 4. Compute Gradient
    initial_temp = jax_state['air_temperature']
    grad_func = jax.grad(compute_loss)
    
    print("\nComputing gradient through array_call interface...")
    gradient = grad_func(initial_temp)
    
    print(f"Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    if jnp.any(gradient != 0.0):
        idx = jnp.argmax(jnp.abs(gradient))
        print(f"Max gradient found at index {idx}: {gradient.flatten()[idx]}")
    
    # Verification
    assert jnp.any(gradient != 0.0), "Gradient is zero. Even with unstable profile, check if logic is truly differentiable."
    assert not jnp.any(jnp.isnan(gradient)), "Gradient contains NaNs"
    print("SUCCESS: Gradients computed successfully through array_call!")

if __name__ == "__main__":
    test_emanuel_jax_differentiation_direct()
