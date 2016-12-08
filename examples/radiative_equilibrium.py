from climt import (
    DataArray, Frierson06GrayLongwaveRadiation, AdamsBashforth,
    PlotFunctionMonitor)
import numpy as np
from datetime import timedelta

pressure_axis = np.array(
    [1e5, 9e4, 8e4, 7e4, 5e4, 3e4, 1e4, 8e3, 4e3, 1e3, 7e2, 4e2, 1e2,
     7., 4., 1.])

state = {
    'air_temperature': DataArray(
        np.ones((1, 1, len(pressure_axis)))*250.,
        dims=('x', 'y', 'half_levels'),
        attrs={'units': 'degK'}),
}
constant_state = {
    'surface_temperature': DataArray(
        np.ones((1, 1))*274., dims=('x', 'y'), attrs={'units': 'degK'}),
    'surface_pressure': DataArray(
        np.ones((1, 1))*1e5, dims=('x', 'y'), attrs={'units': 'Pa'}),
    'air_pressure': DataArray(
        pressure_axis[None, None, :], dims=('x', 'y', 'half_levels'),
        attrs={'units': 'Pa'}),
    'latitude': DataArray(
        np.zeros((1,)), dims=('y',), attrs={'units': 'degrees_north'}),
}
state.update(constant_state)


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].values.flatten())
    ax.axes.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylim(1e5, 100.)

monitor = PlotFunctionMonitor(plot_function)
radiation = Frierson06GrayLongwaveRadiation()
time_stepper = AdamsBashforth([radiation])
timestep = timedelta(hours=4)
for i in range(6*7*4*8):
    print(i)
    new_state = time_stepper.step(state, timestep)
    new_state.update(constant_state)
    if i % 5 == 0:
        monitor.store(state)
    state = new_state
