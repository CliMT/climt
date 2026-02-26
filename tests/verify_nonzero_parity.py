# -*- coding: utf-8 -*-
import numpy as np
from datetime import timedelta
from climt import get_grid, get_default_state, EmanuelConvection
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v2 import EmanuelConvectionPythonV2
from climt._core.util import numpy_version_of

def create_test_state(nlev, ncol, moisture_type='moist'):
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state['air_temperature'].values[:] = 300.0
    q = np.zeros_like(state['specific_humidity'].values)
    if moisture_type == 'moist':
        q[0:5, :, :] = 0.020
    elif moisture_type == 'dry':
        q[:] = 1e-6
    elif moisture_type == 'unstable':
        for i in range(nlev):
            state['air_temperature'].values[i, :, :] = 310.0 - (i * 2.0)
        q[0:10, :, :] = 0.015
    state['specific_humidity'].values[:] = q
    return state

def verify_profile(moisture_type):
    print(f"\n{'='*30} PROFILE: {moisture_type.upper()} {'='*30}")
    nlev = 30
    ncol = 1
    state = create_test_state(nlev, ncol, moisture_type)
    timestep = timedelta(minutes=15)
    raw_state = numpy_version_of(state)
    
    python_state = {
        'air_temperature': raw_state['air_temperature'].reshape(nlev, ncol),
        'specific_humidity': raw_state['specific_humidity'].reshape(nlev, ncol),
        'air_pressure': (raw_state['air_pressure'] / 100.0).reshape(nlev, ncol),
        'air_pressure_on_interface_levels': (raw_state['air_pressure_on_interface_levels'] / 100.0).reshape(nlev + 1, ncol),
        'eastward_wind': raw_state['eastward_wind'].reshape(nlev, ncol),
        'northward_wind': raw_state['northward_wind'].reshape(nlev, ncol),
    }

    conv_fortran = EmanuelConvection()
    conv_orig = EmanuelConvectionPython()
    conv_v2 = EmanuelConvectionPythonV2()

    # 1. Fortran results
    tend_fort, _ = conv_fortran(state, timestep)
    ft_fort = tend_fort['air_temperature'].values.flatten() * 86400
    fq_fort = tend_fort['specific_humidity'].values.flatten() * 86400 * 1000

    # 2. Python results
    tend_orig, _ = conv_orig.array_call(python_state, timestep)
    tend_v2, _ = conv_v2.array_call(python_state, timestep)

    ft_orig = tend_orig['air_temperature'].flatten() * 86400
    fq_orig = tend_orig['specific_humidity'].flatten() * 86400 * 1000
    ft_v2 = tend_v2['air_temperature'].flatten() * 86400
    fq_v2 = tend_v2['specific_humidity'].flatten() * 86400 * 1000
    
    p = raw_state['air_pressure'].flatten() / 100.0

    print("\n--- Temperature Tendency (K/day) ---")
    print(f"{'Lev':<4} | {'Press':<8} | {'T_Fort':<12} | {'T_Py_Orig':<12} | {'T_Py_V2':<12}")
    print("-" * 75)
    
    count_t = 0
    for i in range(nlev):
        if abs(ft_fort[i]) > 1e-8 or abs(ft_orig[i]) > 1e-8 or abs(ft_v2[i]) > 1e-8:
            print(f"{i:<4} | {p[i]:<8.2f} | {ft_fort[i]:12.8f} | {ft_orig[i]:12.8f} | {ft_v2[i]:12.8f}")
            count_t += 1
    
    if count_t == 0:
        print("No non-zero temperature tendencies.")

    print("\n--- Specific Humidity Tendency (g/kg/day) ---")
    print(f"{'Lev':<4} | {'Press':<8} | {'Q_Fort':<12} | {'Q_Py_Orig':<12} | {'Q_Py_V2':<12}")
    print("-" * 75)
    
    count_q = 0
    for i in range(nlev):
        if abs(fq_fort[i]) > 1e-8 or abs(fq_orig[i]) > 1e-8 or abs(fq_v2[i]) > 1e-8:
            print(f"{i:<4} | {p[i]:<8.2f} | {fq_fort[i]:12.8f} | {fq_orig[i]:12.8f} | {fq_v2[i]:12.8f}")
            count_q += 1
            
    if count_q == 0:
        print("No non-zero humidity tendencies.")

if __name__ == "__main__":
    for p in ['moist', 'dry', 'unstable']:
        verify_profile(p)
