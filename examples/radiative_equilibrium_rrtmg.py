from sympl import (
    AdamsBashforth, PlotFunctionMonitor)
from climt import RRTMGShortwave, get_default_state
import numpy as np
from datetime import timedelta


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['shortwave_heating_rate'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    # print(state['air_temperature'].values.flatten())
    # print(state['air_pressure'].values.flatten())
    # ax.set_yscale('log')
    ax.set_ylim(1e5, 1.)


monitor = PlotFunctionMonitor(plot_function)
radiation = RRTMGShortwave()
time_stepper = AdamsBashforth([radiation])
timestep = timedelta(hours=4)

z_information = {'label': 'vertical_level',
                 'values': np.arange(60),
                 'units': ''}

state = get_default_state([radiation], z=z_information)

tp_profiles = np.load('thermodynamic_profiles.npz')
mol_profiles = np.load('molecule_profiles.npz')

state['air_pressure'].values[0, 0, :] = tp_profiles['air_pressure']
state['air_temperature'].values[0, 0, :] = tp_profiles['air_temperature']
state['air_pressure_on_interface_levels'].values[0, 0, :] = tp_profiles['interface_pressures']

# state['specific_humidity'].values[0, 0, :] = mol_profiles['specific_humidity']
state['mole_fraction_of_carbon_dioxide_in_air'].values[0, 0, :] = mol_profiles['carbon_dioxide']
state['mole_fraction_of_ozone_in_air'].values[0, 0, :] = mol_profiles['ozone']

for i in range(8000):

    # print(i)
    diagnostics, new_state = time_stepper.__call__(state, timestep)
    state.update(diagnostics)
    if i % 200 == 0:
        monitor.store(state)
    state = new_state
