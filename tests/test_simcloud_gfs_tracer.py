import sys
import numpy as np
import climt
from sympl import register_tracer, get_tracer_names
from sympl._core.tracers import reset_tracers, reset_packers
from climt import SimCloud

from gfs_dynamical_core import GFSDynamicalCore
from sympl import DataArray
from datetime import datetime, timedelta

def test_simcloud_gfs_tracer():
    print("Setting up GFS Tracer Advection and Conservation Test...")
    
    reset_tracers()
    reset_packers()
    
    register_tracer('cloud_liquid_water_mixing_ratio', 'kg/kg')
    register_tracer('cloud_ice_mixing_ratio', 'kg/kg')
    
    print(f"Registered tracers in sympl: {get_tracer_names()}")
    
    dycore = GFSDynamicalCore()
    grid = climt.get_grid(nx=32, ny=16, nz=28)
    
    state = climt.get_default_state([dycore], grid_state=grid)
    
    # wind field (eastward wind = 25 m/s * cos(lat))
    lat_deg = state['latitude'].values
    lat_rad = np.deg2rad(lat_deg)
    state['eastward_wind'][:] = 25.0 * np.cos(lat_rad)[np.newaxis, :, :]
    state['northward_wind'][:] = 0.0
    
    p = state['air_pressure'].values
    ps = state['surface_air_pressure'].values
    t_surf = 300.0
    gamma = 0.0065
    t_profile = t_surf * (p / ps[np.newaxis, :, :])**(287.0 * gamma / 9.81)
    state['air_temperature'][:] = np.maximum(200.0, t_profile)

    simc = SimCloud()
    temp = state['air_temperature'].values
    tc = temp - 273.15
    tc_clip_w = np.maximum(tc, -100.0)
    tc_clip_i = np.maximum(tc, -100.0)
    
    es_w = np.exp(34.4942 - 4924.99 / (tc_clip_w + 237.1)) / (tc_clip_w + 105.0)**1.57
    es_i = np.exp(43.4942 - 6545.8 / (tc_clip_i + 278.0)) / (tc_clip_i + 868.0)**2.0
    
    liquid_fraction = np.clip((temp - simc.params.T_min) / (simc.params.T_max - simc.params.T_min), 0.0, 1.0)
    sat_vapor_pressure = liquid_fraction * es_w + (1.0 - liquid_fraction) * es_i
    sat_specific_humidity = simc.params.epsilon * sat_vapor_pressure / (p - (1.0 - simc.params.epsilon) * sat_vapor_pressure)
    
    rh = np.where(p < 20000.0, 0.0, 0.99)  # 99% RH below 200 hPa, 0% above
    state['specific_humidity'][:] = np.maximum(rh * sat_specific_humidity, 1e-8)
    
    state['cloud_liquid_water_mixing_ratio'] = DataArray(
        1e-5 * np.cos(lat_rad)[np.newaxis, :, :] * np.ones_like(p),
        dims=state['air_pressure'].dims,
        attrs={'units': 'kg/kg'}
    )
    state['cloud_ice_mixing_ratio'] = DataArray(
        5e-6 * np.cos(lat_rad)[np.newaxis, :, :] * np.ones_like(p),
        dims=state['air_pressure'].dims,
        attrs={'units': 'kg/kg'}
    )
   
    g = dycore._g
    
    def compute_total_mass(mixing_ratio_field, state_obj):
        dp_curr = np.abs(np.diff(state_obj['air_pressure_on_interface_levels'].values, axis=0))
        lat_rad_curr = np.deg2rad(state_obj['latitude'].values)
        cos_lat_curr = np.cos(lat_rad_curr)
        w_curr = cos_lat_curr[np.newaxis, :, :] * dp_curr
        return np.sum(mixing_ratio_field * w_curr) / g

    mass_q_init = compute_total_mass(state['specific_humidity'].values, state)
    mass_ql_init = compute_total_mass(state['cloud_liquid_water_mixing_ratio'].values, state)
    mass_qi_init = compute_total_mass(state['cloud_ice_mixing_ratio'].values, state)
    
    print("\n--- Initial State ---")
    print(f"Total specific humidity mass: {mass_q_init:.6f} kg/m^2 (weighted)")
    print(f"Total cloud liquid mass:      {mass_ql_init:.6f} kg/m^2 (weighted)")
    print(f"Total cloud ice mass:         {mass_qi_init:.6f} kg/m^2 (weighted)")
    
    dt = timedelta(minutes=15)
    num_steps = 2
    
    print(f"\nRunning GFS dynamical core for {num_steps} steps (dt=15m)...")
    
    for step in range(num_steps):
        diags, state = dycore(state, timestep=dt)
        state.update(diags)
        
    # Final masses
    mass_q_final = compute_total_mass(state['specific_humidity'].values, state)
    mass_ql_final = compute_total_mass(state['cloud_liquid_water_mixing_ratio'].values, state)
    mass_qi_final = compute_total_mass(state['cloud_ice_mixing_ratio'].values, state)
    
    print("\n--- Final State ---")
    print(f"Total specific humidity mass: {mass_q_final:.6f} kg/m^2 (weighted)")
    print(f"Total cloud liquid mass:      {mass_ql_final:.6f} kg/m^2 (weighted)")
    print(f"Total cloud ice mass:         {mass_qi_final:.6f} kg/m^2 (weighted)")
    
    # Check conservation (mass should be conserved within a few percent)
    q_err = abs(mass_q_final - mass_q_init) / mass_q_init * 100 if mass_q_init > 0.0 else 0.0
    ql_err = abs(mass_ql_final - mass_ql_init) / mass_ql_init * 100 if mass_ql_init > 0.0 else 0.0
    qi_err = abs(mass_qi_final - mass_qi_init) / mass_qi_init * 100 if mass_qi_init > 0.0 else 0.0
    
    print(f"\nConservation Errors:")
    print(f"Specific humidity error: {q_err:.4f}%")
    print(f"Cloud liquid error:      {ql_err:.4f}%")
    print(f"Cloud ice error:         {qi_err:.4f}%")
    
    # Check that tracers are advected (no longer perfectly uniform)
    std_ql = np.std(state['cloud_liquid_water_mixing_ratio'].values)
    std_qi = np.std(state['cloud_ice_mixing_ratio'].values)
    print(f"Standard deviation of cloud liquid: {std_ql:.4e}")
    print(f"Standard deviation of cloud ice:    {std_qi:.4e}")
    
    # Verify bounds: no negative values
    min_q = np.min(state['specific_humidity'].values)
    min_ql = np.min(state['cloud_liquid_water_mixing_ratio'].values)
    min_qi = np.min(state['cloud_ice_mixing_ratio'].values)
    print(f"Min mixing ratios: q={min_q:.4e}, ql={min_ql:.4e}, qi={min_qi:.4e}")
    
    assert ql_err < 1.0, "Cloud liquid water mass conservation failed (> 1.0%)!"
    assert qi_err < 1.0, "Cloud ice mass conservation failed (> 1.0%)!"
    assert std_ql > 0.0, "Cloud liquid water was not advected (remained uniform)!"
    assert std_qi > 0.0, "Cloud ice was not advected (remained uniform)!"
    assert min_ql >= 0.0, "Negative cloud liquid water values encountered!"
    assert min_qi >= 0.0, "Negative cloud ice values encountered!"
    
    print("\nGFS Tracer Advection and Conservation Test PASSED!")

if __name__ == "__main__":
    try:
        test_simcloud_gfs_tracer()
    except Exception as e:
        print(f"\nGFS test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
