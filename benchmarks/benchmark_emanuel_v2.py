# -*- coding: utf-8 -*-
import time
import numpy as np
from datetime import timedelta
from climt import get_grid, get_default_state, EmanuelConvection
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v2 import EmanuelConvectionPythonV2
from climt._core.util import numpy_version_of

def benchmark_emanuel():
    nlev = 30
    # Use a larger number of columns to see vectorization/parallelization effects
    ncol = 100
    print(f"Benchmarking Emanuel Convection with {ncol} columns, {nlev} levels...")

    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state['air_temperature'].values[:] = 300.0
    state['specific_humidity'].values[0:5, :, :] = 0.020
    
    timestep = timedelta(minutes=15)
    
    # Components
    conv_fortran = EmanuelConvection()
    conv_python_orig = EmanuelConvectionPython()
    conv_python_v2 = EmanuelConvectionPythonV2()

    # 1. Fortran Benchmark
    start = time.time()
    for _ in range(10):
        conv_fortran(state, timestep)
    fortran_time = (time.time() - start) / 10.0
    print(f"Fortran Version:         {fortran_time:.4f} s per call")

    # Prepare raw state for Python versions
    raw_state = numpy_version_of(state)
    def prepare_input(arr):
        shape = arr.shape
        return arr.reshape(shape[0], -1)

    python_state = {
        'air_temperature': prepare_input(raw_state['air_temperature']),
        'specific_humidity': prepare_input(raw_state['specific_humidity']),
        'air_pressure': prepare_input(raw_state['air_pressure'] / 100.0),
        'air_pressure_on_interface_levels': prepare_input(raw_state['air_pressure_on_interface_levels'] / 100.0),
        'eastward_wind': prepare_input(raw_state['eastward_wind']),
        'northward_wind': prepare_input(raw_state['northward_wind']),
    }

    # 2. Original Python Benchmark
    start = time.time()
    for _ in range(10):
        conv_python_orig.array_call(python_state, timestep)
    orig_time = (time.time() - start) / 10.0
    print(f"Original Python:         {orig_time:.4f} s per call")

    # 3. V2 Python Benchmark (JIT)
    print("Warming up V2 Python (JIT)...")
    conv_python_v2.array_call(python_state, timestep)
    
    start = time.time()
    n_iter = 50
    for _ in range(n_iter):
        conv_python_v2.array_call(python_state, timestep)
    v2_time = (time.time() - start) / float(n_iter)
    print(f"V2 Python (JIT):         {v2_time:.4f} s per call")

    # Calculate speedups
    print("-" * 40)
    print(f"Fortran is {orig_time/fortran_time:.1f}x faster than Original Python")
    print(f"V2 is {orig_time/v2_time:.2f}x faster than Original Python")

if __name__ == "__main__":
    benchmark_emanuel()
