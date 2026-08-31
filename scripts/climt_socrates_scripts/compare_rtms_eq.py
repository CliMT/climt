import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))
spectral_base = os.path.join(repo_root, 'climt', '_lib', 'socrates', 'spectral_files_THAI')

from sympl import (
    AdamsBashforth, set_constant, NetCDFMonitor)
from climt import (
    SlabSurface, get_default_state, get_grid,
    SocratesLongwave, SocratesShortwave)

# TRAPPIST-1e Global Constants (based on Fauchez et al. 2020)
set_constant('gravitational_acceleration', 9.114, 'm/s^2')
set_constant('planetary_radius', 5797610.0, 'm')
set_constant('solar_constant', 900.0 / 4.0, 'W m^-2')
set_constant('gas_constant_of_dry_air', 296.8, 'J kg^-1 K^-1')
set_constant('heat_capacity_of_dry_air_at_constant_pressure', 1039.0, 'J/kg/K')
set_constant('zenith_angle', np.pi*0.25, 'dimensionless')

images_dir = os.path.join(script_dir, "plots")
netcdf_dir = os.path.join(images_dir, "netcdf")
os.makedirs(netcdf_dir, exist_ok=True)

# Band-resolved flux diagnostics (see socrates_runes.F90): each run_simulation()
# call points these at a case-specific tag so the SW/LW band flux time series
# for DSA Mars and TRAPPIST-1e don't get appended into the same file.
band_flux_dir = os.path.join(images_dir, "band_fluxes")
os.makedirs(band_flux_dir, exist_ok=True)
os.environ['SOCRATES_BAND_FLUX_DIR'] = band_flux_dir
# Log every 100th call to match the NetCDF snapshot cadence below, rather than
# every hourly timestep, to keep the files a manageable size.
os.environ['SOCRATES_BAND_FLUX_STRIDE'] = '100'

# Define file paths for the comparison
dsa_mars_sw_path = os.path.join(spectral_base, 'dsa_mars', 'sp_sw_42_dsa_mars_trappist1')
dsa_mars_lw_path = os.path.join(spectral_base, 'dsa_mars', 'sp_lw_17_dsa_mars')
trappist1e_sw_path = os.path.join(spectral_base, 'trappist_1e', 'sp_sw_42_trappist1e')
trappist1e_lw_path = os.path.join(spectral_base, 'trappist_1e', 'sp_lw_17_trappist1e')

# Composition cases from Fauchez et al. 2020 paper
composition_cases = {
    # 'No_Absorbers': {
    #     'title': 'No Molecular Absorbers (Pure Scattering/Grey)',
    #     'co2': 0.0,
    #     'n2': 0.0,
    #     'h2o': 0.0
    # },
    'Ben1_N2_CO2': {
        'title': 'Ben1/Hab1 (1 bar N2 + 400 ppm CO2)',
        'co2': 400e-6,
        'n2': 1.0 - 400e-6,
        'h2o': 0.0
    },
    'Ben2_CO2_only': {
        'title': 'Ben2/Hab2 (1 bar pure CO2)',
        'co2': 1.0,
        'n2': 0.0,
        'h2o': 0.0
    }
}

# Fields to store in NetCDF files
fields_to_store = [
    'air_temperature',
    'air_pressure',
    'air_pressure_on_interface_levels',
    'surface_temperature',
    'upwelling_longwave_flux_in_air',
    'downwelling_longwave_flux_in_air',
    'upwelling_shortwave_flux_in_air',
    'downwelling_shortwave_flux_in_air',
    'air_temperature_tendency_from_shortwave',
    'air_temperature_tendency_from_longwave',
    'mole_fraction_of_carbon_dioxide_in_air',
    'mole_fraction_of_nitrogen_in_air',
    'specific_humidity'
]

def run_simulation(sw_file, lw_file, sw_bands, lw_bands, case_name, file_suffix, comp_dict):
    print(f"\n==================================================")
    print(f"Running: {case_name}")
    print(f"SW file: {os.path.basename(sw_file)}")
    print(f"LW file: {os.path.basename(lw_file)}")
    print(f"==================================================")

    # Point the Fortran band-flux diagnostic at a file tag unique to this
    # case, and clear out any files left over from a previous run of this
    # script (the Fortran side appends, so stale rows would otherwise mix
    # with this run's call_index sequence).
    os.environ['SOCRATES_BAND_FLUX_TAG'] = file_suffix
    for source in ('sw', 'lw'):
        stale_path = os.path.join(band_flux_dir, f"band_fluxes_{file_suffix}_{source}.txt")
        if os.path.exists(stale_path):
            os.remove(stale_path)

    # Initialize radiation and slab surface components
    socrates_sw = SocratesShortwave(spectral_file=sw_file, n_band=sw_bands, albedo_sw=0.3)
    socrates_lw = SocratesLongwave(spectral_file=lw_file, n_band=lw_bands)
    surface = SlabSurface()
    
    time_stepper = AdamsBashforth([socrates_sw, socrates_lw, surface])
    timestep = timedelta(hours=1)
    
    # Initialize state with 60 levels
    grid = get_grid(nz=60)
    state = get_default_state([socrates_sw, socrates_lw, surface], grid_state=grid)
    
    # Set time field in state
    state['time'] = datetime(2000, 1, 1)
    
    # Anchor reference pressure to 1 bar
    p_surf = 1e5
    state['air_pressure'][:] = state['air_pressure'] * (p_surf / 101325.0)
    state['air_pressure_on_interface_levels'][:] = state['air_pressure_on_interface_levels'] * (p_surf / 101325.0)
    
    # Start from an isothermal profile of 300K
    state['air_temperature'][:] = 300.0
    state['surface_temperature'][:] = 300.0
    
    # Set composition
    state['specific_humidity'][:] = comp_dict['h2o'] if comp_dict['h2o'] > 0 else 1e-12
    state['mole_fraction_of_carbon_dioxide_in_air'][:] = comp_dict['co2']
    state['mole_fraction_of_nitrogen_in_air'][:] = comp_dict['n2']
    
    # Initialize some required surface fields
    state['ocean_mixed_layer_thickness'][:] = 1.0
    state['area_type'][:] = 'land' 
    state['surface_material_density'][:] = 2000.0 
    state['heat_capacity_of_soil'][:] = 1000.0
    state['surface_upward_latent_heat_flux'][:] = 0.0
    state['surface_upward_sensible_heat_flux'][:] = 0.0
    state['upward_heat_flux_at_ground_level_in_soil'][:] = 0.0
    state['heat_flux_into_sea_water_due_to_sea_ice'][:] = 0.0
    
    # Initialize NetCDF monitor
    nc_filename = f"{file_suffix}.nc"
    nc_filepath = os.path.join(netcdf_dir, nc_filename)
    if os.path.exists(nc_filepath):
        os.remove(nc_filepath)
        
    netcdf_monitor = NetCDFMonitor(
        nc_filepath,
        write_on_store=True,
        store_names=fields_to_store
    )
    
    n_steps = 15000
    last_100_Ts = 300.0
    
    for i in range(n_steps):
        diagnostics, state = time_stepper(state, timestep)
        state.update(diagnostics)
        state['time'] += timestep
        
        # Store state in NetCDF file every 100 steps
        if i % 100 == 0:
            netcdf_monitor.store(state)
            
            ts_val = state['surface_temperature'].values.item()
            if i > 0:
                temp_diff = np.abs(ts_val - last_100_Ts)
                # If temperature change over 100 hours is less than 5e-4 K, we have converged
                if temp_diff < 5e-4:
                    print(f"  Step {i:5d}/{n_steps:5d}: Ts = {ts_val:7.2f} K (Converged! diff={temp_diff:.2e})")
                    break
                else:
                    print(f"  Step {i:5d}/{n_steps:5d}: Ts = {ts_val:7.2f} K, diff={temp_diff:.2e}")
            else:
                print(f"  Step {i:5d}/{n_steps:5d}: Ts = {ts_val:7.2f} K")
            last_100_Ts = ts_val
            
    # Save the very final converged state to the NetCDF monitor
    netcdf_monitor.store(state)
    print(f"Saved NetCDF output to: {nc_filepath}")
            
    return state

def main():
    for comp_name, comp_dict in composition_cases.items():
        # Run DSA Mars
        state_mars = run_simulation(
            sw_file=dsa_mars_sw_path,
            lw_file=dsa_mars_lw_path,
            sw_bands=42,
            lw_bands=17,
            case_name=f"{comp_name} (DSA Mars Spectral)",
            file_suffix=f"dsa_mars_{comp_name}",
            comp_dict=comp_dict
        )
        
        # Run TRAPPIST-1e
        state_trappist = run_simulation(
            sw_file=trappist1e_sw_path,
            lw_file=trappist1e_lw_path,
            sw_bands=42,
            lw_bands=17,
            case_name=f"{comp_name} (TRAPPIST-1e Spectral)",
            file_suffix=f"trappist1e_{comp_name}",
            comp_dict=comp_dict
        )
        
        # Generate the comparison plot
        fig, (ax_temp, ax_flux) = plt.subplots(1, 2, figsize=(14, 6))
        
        # 1. Temperature Profile Plot
        p_levels_mars = state_mars['air_pressure'].values.flatten() / 100.0
        temp_mars = state_mars['air_temperature'].values.flatten()
        
        p_levels_trappist = state_trappist['air_pressure'].values.flatten() / 100.0
        temp_trappist = state_trappist['air_temperature'].values.flatten()
        
        # Plot profiles (solid blue for trappist1e, dashed magenta for dsa mars)
        ax_temp.plot(temp_trappist, p_levels_trappist, label='TRAPPIST-1e Spectral Files', color='#0066cc', linewidth=2.5)
        ax_temp.plot(temp_mars, p_levels_mars, label='DSA Mars Spectral Files', color='#cc0066', linewidth=2.5, linestyle='--')
        
        # Scatter markers for surface temperature
        ts_trappist = state_trappist['surface_temperature'].values.item()
        ts_mars = state_mars['surface_temperature'].values.item()
        ax_temp.scatter([ts_trappist], [p_levels_trappist[-1]], color='#0066cc', s=80, zorder=5)
        ax_temp.scatter([ts_mars], [p_levels_mars[-1]], color='#cc0066', s=80, zorder=5)
        
        ax_temp.set_yscale('log')
        ax_temp.set_ylim(1000.0, 0.1)
        ax_temp.axes.invert_yaxis()
        ax_temp.grid(True, which="both", ls=":")
        ax_temp.set_xlabel('Air Temperature (K)', fontsize=12)
        ax_temp.set_ylabel('Pressure (hPa)', fontsize=12)
        ax_temp.set_title(f'Temperature Profiles\nTRAPPIST-1e Ts: {ts_trappist:.1f} K | DSA Mars Ts: {ts_mars:.1f} K', fontsize=11)
        ax_temp.legend(fontsize=10, loc='best')
        
        # 2. Net Radiative Fluxes Plot
        p_interface_mars = state_mars['air_pressure_on_interface_levels'].values.flatten() / 100.0
        p_interface_trappist = state_trappist['air_pressure_on_interface_levels'].values.flatten() / 100.0
        
        # Net LW = Up - Down
        lw_net_mars = state_mars['upwelling_longwave_flux_in_air'].values.flatten() - state_mars['downwelling_longwave_flux_in_air'].values.flatten()
        lw_net_trappist = state_trappist['upwelling_longwave_flux_in_air'].values.flatten() - state_trappist['downwelling_longwave_flux_in_air'].values.flatten()
        
        # Net SW = Down - Up
        sw_net_mars = state_mars['downwelling_shortwave_flux_in_air'].values.flatten() - state_mars['upwelling_shortwave_flux_in_air'].values.flatten()
        sw_net_trappist = state_trappist['downwelling_shortwave_flux_in_air'].values.flatten() - state_trappist['upwelling_shortwave_flux_in_air'].values.flatten()
        
        # Plot Net Fluxes
        ax_flux.plot(sw_net_trappist, p_interface_trappist, label='Net SW (Down-Up) - TRAPPIST-1e', color='#1f77b4', linewidth=2)
        ax_flux.plot(sw_net_mars, p_interface_mars, label='Net SW (Down-Up) - DSA Mars', color='#aec7e8', linewidth=2, linestyle='--')
        
        ax_flux.plot(lw_net_trappist, p_interface_trappist, label='Net LW (Up-Down) - TRAPPIST-1e', color='#ff7f0e', linewidth=2)
        ax_flux.plot(lw_net_mars, p_interface_mars, label='Net LW (Up-Down) - DSA Mars', color='#ffbb78', linewidth=2, linestyle='--')
        
        ax_flux.set_yscale('log')
        ax_flux.set_ylim(1000.0, 0.1)
        ax_flux.axes.invert_yaxis()
        ax_flux.grid(True, which="both", ls=":")
        ax_flux.set_xlabel('Net Flux (W/m^2)', fontsize=12)
        ax_flux.set_ylabel('Pressure (hPa)', fontsize=12)
        ax_flux.set_title('Net Radiative Fluxes\n(SW Positive Downward, LW Positive Upward)', fontsize=11)
        ax_flux.legend(fontsize=10, loc='best')
        
        plt.suptitle(comp_dict['title'], fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        plot_path = os.path.join(images_dir, f"compare_rtms_{comp_name}.png")
        plt.savefig(plot_path, dpi=150)
        plt.close(fig)
        print(f"Saved comparison plot to {plot_path}\n")

if __name__ == '__main__':
    main()
