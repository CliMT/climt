# -*- coding: utf-8 -*-
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

    conv_orig = EmanuelConvectionPython()
    conv_v2 = EmanuelConvectionPythonV2()
    conv_v3 = EmanuelConvectionPythonV3()

    # Python results (Numba path)
    tend_orig, diag_orig = conv_orig.array_call(python_state, timestep)
    tend_v2, diag_v2 = conv_v2.array_call(python_state, timestep)
    tend_v3, diag_v3 = conv_v3.array_call(python_state, timestep)

    # JAX results
    if HAS_JAX:
        jax_state = {k: jnp.array(v) for k, v in python_state.items()}
        print("Calling V3 JAX...")
        tend_jax, diag_jax = conv_v3.array_call(jax_state, timestep)
        ft_jax = np.array(tend_jax['air_temperature']).flatten() * 86400
        fq_jax = np.array(tend_jax['specific_humidity']).flatten() * 86400 * 1000
        cbmf_jax = float(diag_jax['cloud_base_mass_flux'][0])
    else:
        ft_jax = np.zeros(nlev)
        fq_jax = np.zeros(nlev)
        cbmf_jax = 0.0

    print(f"CBMF Comparison: Orig={diag_orig['cloud_base_mass_flux'][0]:.6f}, V2={diag_v2['cloud_base_mass_flux'][0]:.6f}, V3_Numba={diag_v3['cloud_base_mass_flux'][0]:.6f}, V3_JAX={cbmf_jax:.6f}")

if __name__ == "__main__":
    for p in ['moist', 'unstable']:
        verify_profile(p)
