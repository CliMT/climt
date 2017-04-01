from sympl import (
    AdamsBashforth, PlotFunctionMonitor)
from climt import RRTMGLongwave, get_default_state
import numpy as np
from datetime import timedelta


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    # print(state['air_temperature'].values.flatten())
    # print(state['air_pressure'].values.flatten())
    # ax.set_yscale('log')
    ax.set_ylim(1e5, 1.)


monitor = PlotFunctionMonitor(plot_function)
radiation = RRTMGLongwave()
time_stepper = AdamsBashforth([radiation])
timestep = timedelta(hours=4)

z_information = {'label': 'vertical_level',
                 'values': np.arange(60),
                 'units': ''}

state = get_default_state([radiation], z=z_information)

# format: p_mid, t_mid, p_int, t_int
tp_profiles = np.loadtxt('tp_profiles').transpose()
# format: h2o, co2, o3
mol_profiles = np.loadtxt('mol_profiles').transpose()

state['air_pressure'].values[0, 0, :] = tp_profiles[0, :]*100
state['air_temperature'].values[0, 0, :] = tp_profiles[1, :]
state['air_pressure_on_interface_levels'].values[0, 0, 1:] = tp_profiles[2, :]*100
state['air_pressure_on_interface_levels'].values[0, 0, 0] = 101300.

state['specific_humidity'].values[0, 0, :] = mol_profiles[0, :]*(18.02/28.964)*1000
state['carbon_dioxide_mixing_ratio'].values[0, 0, :] = mol_profiles[1, :]
state['ozone_mixing_ratio'].values[0, 0, :] = mol_profiles[2, :]

for i in range(8000):

    # print(i)
    diagnostics, new_state = time_stepper.__call__(state, timestep)
    state.update(diagnostics)
    if i % 200 == 0:
        monitor.store(state)
    state = new_state
