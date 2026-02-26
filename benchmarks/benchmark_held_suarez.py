# -*- coding: utf-8 -*-
import time
import numpy as np
from datetime import timedelta
import sympl
from climt import get_grid, get_default_state, HeldSuarez, JaxBackend, UnytBackend

try:
    import jax
    import jax.numpy as jnp
    DEVICE = jax.devices()[0]
    PLATFORM = DEVICE.platform.upper()
    if PLATFORM == "METAL":
        jax.config.update("jax_enable_x64", False)
    else:
        jax.config.update("jax_enable_x64", True)
    HAS_JAX = True
except ImportError:
    HAS_JAX = False
    PLATFORM = "NONE"

def benchmark_held_suarez(ncol=1000, n_iter=100):
    nlev = 30
    print(f"\nBenchmarking Held-Suarez with {ncol} columns, {nlev} levels...")
    print(f"JAX Platform: {PLATFORM}")

    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    hs = HeldSuarez()

    # 1. Numba Performance (using UnytBackend)
    unyt_backend = UnytBackend()
    sympl.set_backend(unyt_backend)
    # Re-call get_grid to ensure it uses the new backend if it creates quantities
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    state_unyt = get_default_state([hs], grid_state=grid)
    
    # Warmup Numba
    hs(state_unyt)
    
    start = time.time()
    for _ in range(n_iter):
        hs(state_unyt)
    numba_time = (time.time() - start) / float(n_iter)
    print(f"V2 Numba JIT (Unyt CPU):    {numba_time:.6f} s per call")

    # 2. JAX JIT Performance
    if HAS_JAX:
        jax_backend = JaxBackend()
        sympl.set_backend(jax_backend)
        grid = get_grid(nx=ncol, ny=1, nz=nlev)
        jax_state = get_default_state([hs], grid_state=grid)
        
        # Warmup JAX
        hs(jax_state)
        
        start = time.time()
        for _ in range(n_iter):
            res_tend, _ = hs(jax_state)
            # block_until_ready for accurate timing on GPU
            res_tend['air_temperature'].values.block_until_ready()
        jax_time = (time.time() - start) / float(n_iter)
        print(f"V2 JAX JIT ({PLATFORM}):       {jax_time:.6f} s per call")
    else:
        jax_time = None

    print("-" * 40)
    if jax_time:
        print(f"JAX ({PLATFORM}) is {numba_time/jax_time:.1f}x faster than Numba (Unyt CPU)")

if __name__ == "__main__":
    benchmark_held_suarez(ncol=1000, n_iter=100)
    benchmark_held_suarez(ncol=10000, n_iter=50)
    benchmark_held_suarez(ncol=100000, n_iter=10)
