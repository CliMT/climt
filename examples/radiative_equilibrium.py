from climt import DataArray, Frierson06GrayLongwaveRadiation, AdamsBashforth
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt
import time


pressure_axis = np.array(
    [1e5, 9e4, 8e4, 7e4, 5e4, 3e4, 1e4, 8e3, 4e3, 1e3, 7e2, 4e2, 1e2, 7., 4., 1.])

nz = len(pressure_axis)

state = {
    'air_temperature': DataArray(
        np.ones((1, 1, nz))*250., dims=('x', 'y', 'half_levels'),
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

plt.ion()
fig = plt.figure()
plt.plot(
    state['air_temperature'].values.flatten(),
    state['air_pressure'].values.flatten())
plt.gca().invert_yaxis()
plt.yscale('log')
plt.ylim(1e5, 100.)

radiation = Frierson06GrayLongwaveRadiation()
time_stepper = AdamsBashforth([radiation])
timestep = timedelta(hours=4)
new_state = time_stepper.step(state, timestep)
state = new_state
state.update(constant_state)
start = time.time()
for i in range(6*7*4*8):
    new_state = time_stepper.step(state, timestep)
    state = new_state
    state.update(constant_state)
    if i % 5 == 0:
        plt.clf()
        plt.plot(
            state['air_temperature'].values.flatten(),
            state['air_pressure'].values.flatten())
        plt.yscale('log')
        plt.ylim(1e5, 10.)
        plt.draw_all()
        plt.pause(1e-3)

print('Time elapsed: {} seconds'.format(time.time() - start))

