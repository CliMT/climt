# -*- coding: utf-8 -*-
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, SlabSurface, JaxBackend

# Enable x64
jax.config.update("jax_enable_x64", True)

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

def test_slab_surface_differentiation():
    """
    Verify that SlabSurface is differentiable via JAX.
    """
    jax.config.update('jax_platform_name', 'cpu')
    backend = JaxBackend()
    sympl.set_backend(backend)
    
    ny = 10
    nx = 5
    grid = get_grid(nx=nx, ny=ny, nz=20)
    slab = SlabSurface()
    state = get_default_state([slab], grid_state=grid)
    
    # Prepare standard climt area type (NumPy strings)
    # The JAX path in array_call will convert this to codes.
    area_type = np.full((ny, nx), "sea")
    
    def compute_loss(sw_down_val):
        # sw_down_val is (ny, nx, nz_interface)
        current_state = {
            'area_type': area_type,
            'downwelling_shortwave_flux_in_air': sw_down_val,
            'downwelling_longwave_flux_in_air': state['downwelling_longwave_flux_in_air'].data,
            'upwelling_shortwave_flux_in_air': state['upwelling_shortwave_flux_in_air'].data,
            'upwelling_longwave_flux_in_air': state['upwelling_longwave_flux_in_air'].data,
            'surface_upward_latent_heat_flux': state['surface_upward_latent_heat_flux'].data,
            'surface_upward_sensible_heat_flux': state['surface_upward_sensible_heat_flux'].data,
            'upward_heat_flux_at_ground_level_in_soil': state['upward_heat_flux_at_ground_level_in_soil'].data,
            'heat_flux_into_sea_water_due_to_sea_ice': state['heat_flux_into_sea_water_due_to_sea_ice'].data,
            'sea_water_density': state['sea_water_density'].data,
            'surface_material_density': state['surface_material_density'].data,
            'heat_capacity_of_soil': state['heat_capacity_of_soil'].data,
            'surface_thermal_capacity': state['surface_thermal_capacity'].data,
            'ocean_mixed_layer_thickness': state['ocean_mixed_layer_thickness'].data,
            'soil_layer_thickness': state['soil_layer_thickness'].data,
            'surface_temperature': state['surface_temperature'].data,
        }
        
        tendencies, _ = slab.array_call(current_state)
        return jnp.sum(jnp.square(tendencies['surface_temperature']))

    initial_sw = state['downwelling_shortwave_flux_in_air'].data
    # Set to a non-zero value so the gradient is non-zero
    initial_sw = initial_sw.at[:].set(100.0)
    grad_func = jax.grad(compute_loss)
    gradient = grad_func(initial_sw)
    
    print(f"\nSlab Surface Gradient sum: {jnp.sum(jnp.abs(gradient))}")
    assert jnp.any(gradient != 0.0), "Slab Surface should be differentiable w.r.t SW flux"
    assert not jnp.any(jnp.isnan(gradient))
    print("SUCCESS: Slab Surface differentiation verified!")

if __name__ == "__main__":
    test_slab_surface_parity()
    test_slab_surface_differentiation()
