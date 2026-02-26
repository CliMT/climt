# -*- coding: utf-8 -*-
import time
import numpy as np
from datetime import timedelta
from climt import get_grid, get_default_state, EmanuelConvection
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v2 import EmanuelConvectionPythonV2
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPythonV3
from climt._core.util import numpy_version_of

try:
    import jax
    import jax.numpy as jnp
    HAS_JAX = True
    DEVICE = jax.devices()[0]
    PLATFORM = DEVICE.platform.upper()
    # Disable x64 for METAL as it is not supported
    if PLATFORM == "METAL":
        jax.config.update("jax_enable_x64", False)
        print("Disabling JAX x64 for METAL (GPU) backend...")
    else:
        jax.config.update("jax_enable_x64", True)
except ImportError:
    HAS_JAX = False
    PLATFORM = "NONE"

def benchmark_emanuel(ncol=10, n_iter=10):
    nlev = 30
    print(f"\nBenchmarking Emanuel Convection with {ncol} columns, {nlev} levels...")
    print(f"JAX Platform: {PLATFORM}")

    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state['air_temperature'].values[:] = 300.0
    state['specific_humidity'].values[0:5, :, :] = 0.020
    timestep = timedelta(minutes=15)
    
    conv_fortran = EmanuelConvection()
    conv_python_orig = EmanuelConvectionPython()
    conv_python_v3 = EmanuelConvectionPythonV3()

    # 1. Fortran
    start = time.time()
    for _ in range(n_iter): conv_fortran(state, timestep)
    fortran_time = (time.time() - start) / float(n_iter)
    print(f"Fortran Version:         {fortran_time:.4f} s per call")

    # 2. Original Python
    raw_state = numpy_version_of(state)
    def prepare_input(arr): return arr.reshape(arr.shape[0], -1)
    python_state = {
        'air_temperature': prepare_input(raw_state['air_temperature']),
        'specific_humidity': prepare_input(raw_state['specific_humidity']),
        'air_pressure': prepare_input(raw_state['air_pressure'] / 100.0),
        'air_pressure_on_interface_levels': prepare_input(raw_state['air_pressure_on_interface_levels'] / 100.0),
        'eastward_wind': prepare_input(raw_state['eastward_wind']),
        'northward_wind': prepare_input(raw_state['northward_wind']),
    }
    start = time.time()
    for _ in range(n_iter): conv_python_orig.array_call(python_state, timestep)
    orig_time = (time.time() - start) / float(n_iter)
    print(f"Original Python:         {orig_time:.4f} s per call")

    # 3. Numba
    print(f"Warming up V3 Numba...")
    conv_python_v3.array_call(python_state, timestep)
    start = time.time()
    for _ in range(n_iter): conv_python_v3.array_call(python_state, timestep)
    v3_numba_time = (time.time() - start) / float(n_iter)
    print(f"V3 Python (Numba JIT):   {v3_numba_time:.4f} s per call")

    # 4. JAX
    if HAS_JAX:
        print(f"Warming up V3 JAX ({PLATFORM})...")
        # Ensure float32 for Metal
        dtype = jnp.float32 if PLATFORM == "METAL" else jnp.float64
        jax_state = {k: jnp.array(v, dtype=dtype) for k, v in python_state.items()}
        conv_python_v3.array_call(jax_state, timestep)
        
        start = time.time()
        for _ in range(n_iter):
            res_tend, _ = conv_python_v3.array_call(jax_state, timestep)
            res_tend['air_temperature'].block_until_ready()
        v3_jax_time = (time.time() - start) / float(n_iter)
        print(f"V3 Python (JAX JIT):     {v3_jax_time:.4f} s per call")
    else:
        v3_jax_time = None

    print("-" * 40)
    print(f"V3 (Numba) is {orig_time/v3_numba_time:.2f}x faster than Original Python")
    if v3_jax_time:
        print(f"V3 (JAX) is {orig_time/v3_jax_time:.2f}x faster than Original Python")

if __name__ == "__main__":
    benchmark_emanuel(ncol=10, n_iter=10)
    benchmark_emanuel(ncol=1000, n_iter=5)
