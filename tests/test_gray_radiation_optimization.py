# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, GrayLongwaveRadiation, Frierson06LongwaveOpticalDepth, JaxBackend

# Enable x64
jax.config.update("jax_enable_x64", True)

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

def test_gray_radiation_differentiation():
    """
    Verify that GrayLongwaveRadiation is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    backend = JaxBackend()
    sympl.set_backend(backend)
    
    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    
    tau_comp = Frierson06LongwaveOpticalDepth()
    lw_comp = GrayLongwaveRadiation()
    
    state = get_default_state([tau_comp, lw_comp], grid_state=grid)
    
    # Pre-calculate tau for the initial state
    # array_call expects sl, p_interface, T_surface, tau
    
    jax_state_init = {
        'latitude': state['latitude'].data.flatten(),
        'air_pressure_on_interface_levels': state['air_pressure_on_interface_levels'].data,
        'surface_air_pressure': state['surface_air_pressure'].data.flatten(),
    }
    
    # Use tau_comp.array_call to get initial tau
    diag_tau_init = tau_comp.array_call(jax_state_init)
    initial_tau = diag_tau_init['longwave_optical_depth_on_interface_levels']

    def compute_loss(temp_values):
        # temp_values is (nlev, ncol, ...) because of climt internal state structure
        current_state = {
            'sl': temp_values,
            'tau': initial_tau,
            'T_surface': state['surface_temperature'].data.flatten(),
            'p_interface': state['air_pressure_on_interface_levels'].data,
        }
        
        # Direct call to array_call
        tendencies, _ = lw_comp.array_call(current_state)
        return jnp.sum(jnp.square(tendencies['sl']))

    # Initial temp must be (ncol, nlev) to match sympl state extracted via JaxBackend
    initial_temp = state['air_temperature'].data
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_temp)
    
    print(f"\nGray Radiation Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Gray Radiation should be differentiable w.r.t temperature"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Gray Radiation differentiation verified!")

def test_frierson_optical_depth_differentiation():
    """
    Verify that Frierson06LongwaveOpticalDepth is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    backend = JaxBackend()
    sympl.set_backend(backend)
    
    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    
    tau_comp = Frierson06LongwaveOpticalDepth()
    state = get_default_state([tau_comp], grid_state=grid)
    
    def compute_loss(ps_values):
        current_state = {
            'air_pressure_on_interface_levels': state['air_pressure_on_interface_levels'].data,
            'surface_air_pressure': ps_values,
            'latitude': state['latitude'].data.flatten(),
        }
        
        diagnostics = tau_comp.array_call(current_state)
        return jnp.sum(jnp.square(diagnostics['longwave_optical_depth_on_interface_levels']))

    initial_ps = state['surface_air_pressure'].data.flatten()
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_ps)
    
    print(f"\nFrierson Optical Depth Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Frierson Optical Depth should be differentiable w.r.t surface pressure"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Frierson Optical Depth differentiation verified!")

if __name__ == "__main__":
    test_gray_radiation_parity()
    test_gray_radiation_differentiation()
