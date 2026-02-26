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
    jax.config.update('jax_platform_name', 'cpu')
    jax.config.update("jax_enable_x64", True)
    HAS_JAX = True
except ImportError:
    HAS_JAX = False

def benchmark_emanuel():
    nlev = 30
    ncol = 10 # Reduced for faster debugging
    print(f"Benchmarking Emanuel Convection with {ncol} columns, {nlev} levels...")

    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state['air_temperature'].values[:] = 300.0
    state['specific_humidity'].values[0:5, :, :] = 0.020
    timestep = timedelta(minutes=15)
    
    conv_fortran = EmanuelConvection()
    conv_python_orig = EmanuelConvectionPython()
    conv_python_v3 = EmanuelConvectionPythonV3()

    start = time.time()
    for _ in range(10): conv_fortran(state, timestep)
    fortran_time = (time.time() - start) / 10.0
    print(f"Fortran Version:         {fortran_time:.4f} s per call")

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
    for _ in range(10): conv_python_orig.array_call(python_state, timestep)
    orig_time = (time.time() - start) / 10.0
    print(f"Original Python:         {orig_time:.4f} s per call")

    print("Warming up V3 Python (Numba JIT)...")
    conv_python_v3.array_call(python_state, timestep)
    start = time.time(); n_iter = 10
    for _ in range(n_iter): conv_python_v3.array_call(python_state, timestep)
    v3_numba_time = (time.time() - start) / float(n_iter)
    print(f"V3 Python (Numba JIT):   {v3_numba_time:.4f} s per call")

    if HAS_JAX:
        print("Warming up V3 Python (JAX JIT - Compilation might take a minute)...")
        jax_state = {k: jnp.array(v) for k, v in python_state.items()}
        s = time.time()
        conv_python_v3.array_call(jax_state, timestep)
        print(f"JAX Warmup/Compilation took {time.time() - s:.2f} s")
        
        start = time.time()
        jax_iter = 5
        for _ in range(jax_iter):
            res_tend, _ = conv_python_v3.array_call(jax_state, timestep)
            res_tend['air_temperature'].block_until_ready()
        v3_jax_time = (time.time() - start) / float(jax_iter)
        print(f"V3 Python (JAX JIT):     {v3_jax_time:.4f} s per call")
    else:
        v3_jax_time = None

    print("-" * 40)
    print(f"Fortran is {orig_time/fortran_time:.1f}x faster than Original Python")
    print(f"V3 (Numba) is {orig_time/v3_numba_time:.2f}x faster than Original Python")
    if v3_jax_time:
        print(f"V3 (JAX) is {orig_time/v3_jax_time:.2f}x faster than Original Python")

if __name__ == "__main__":
    benchmark_emanuel()
