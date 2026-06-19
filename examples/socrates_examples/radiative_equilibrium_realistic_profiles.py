import os
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt
from sympl import AdamsBashforth, PlotFunctionMonitor
from climt import (
    SocratesShortwave, SocratesLongwave, get_default_state, get_grid
)


def plot_function(fig, state):
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(
        state['air_temperature_tendency_from_shortwave'].values.flatten(),
        state['air_pressure'].values.flatten() / 100, '-o', label='SW')
    ax.axes.invert_yaxis()
    ax.plot(
        state['air_temperature_tendency_from_longwave'].values.flatten(),
        state['air_pressure'].values.flatten() / 100, '-o', label='LW')
    ax.axes.invert_yaxis()
    ax.set_title('Heating Rates')
    ax.grid()
    ax.set_xlabel('K/day')
    ax.set_ylabel('millibar')
    ax.legend()
    ax.set_yscale('log')
    ax.set_ylim(1e3, 10.)

    ax = fig.add_subplot(1, 2, 2)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].values.flatten() / 100, '-o')
    ax.axes.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylim(1e3, 10.)
    ax.set_title('Temperature')
    ax.grid()
    ax.set_xlabel('K')
    ax.set_yticklabels([])

    plt.suptitle('Radiative Equilibrium with SOCRATES (Realistic Profiles)')


# Initialize monitor and components
monitor = PlotFunctionMonitor(plot_function)
rad_sw = SocratesShortwave(spectral_file='sp_sw_ga9')
rad_lw = SocratesLongwave(spectral_file='sp_lw_ga9')
time_stepper = AdamsBashforth([rad_sw, rad_lw])
timestep = timedelta(hours=3)

# Build default state on a 60-level grid to match realistic profiles
grid = get_grid(nx=1, ny=1, nz=60)
state = get_default_state([rad_sw, rad_lw], grid_state=grid)

# Resolve paths to the npz files in the parent directory
script_dir = os.path.dirname(os.path.abspath(__file__))
tp_path = os.path.abspath(
    os.path.join(script_dir, "..", "thermodynamic_profiles.npz")
)
mol_path = os.path.abspath(
    os.path.join(script_dir, "..", "molecule_profiles.npz")
)

print(f"Loading profiles:\n - {tp_path}\n - {mol_path}")
tp_profiles = np.load(tp_path)
mol_profiles = np.load(mol_path)

state['air_pressure'].values[:] = (
    tp_profiles['air_pressure'][:, np.newaxis, np.newaxis]
)
state['air_temperature'].values[:] = (
    tp_profiles['air_temperature'][:, np.newaxis, np.newaxis]
)
state['air_pressure_on_interface_levels'].values[:] = (
    tp_profiles['interface_pressures'][:, np.newaxis, np.newaxis]
)

state['specific_humidity'].values[:] = (
    mol_profiles['specific_humidity'][:, np.newaxis, np.newaxis] * 1e-3
)
state['mole_fraction_of_carbon_dioxide_in_air'].values[:] = (
    mol_profiles['carbon_dioxide'][:, np.newaxis, np.newaxis]
)
state['mole_fraction_of_ozone_in_air'].values[:] = (
    mol_profiles['ozone'][:, np.newaxis, np.newaxis]
)

print("Running Radiative Equilibrium (Realistic Profiles) example...")
# Run time loop
for i in range(1000):
    diagnostics, new_state = time_stepper(state, timestep)
    state.update(diagnostics)
    if i % 100 == 0:
        monitor.store(state)
        print(f"Step {i:04d} complete")
    state = new_state

print("Simulation finished. Close the plot window to exit.")
plt.show()
