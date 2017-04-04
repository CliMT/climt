from climt import RRTMGLongwave, get_default_state
import numpy as np


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


radiation = RRTMGLongwave()

z_information = {'label': 'vertical_level',
                 'values': np.arange(30),
                 'units': ''}

y_information = {'label': 'latitude',
                 'values': np.arange(10),
                 'units': 'degrees north'}

state = get_default_state([radiation], y=y_information, z=z_information)
one_d_state = get_default_state([radiation])
# print(state['specific_humidity'])
state['specific_humidity'].values[:, :, :17] = 10.
one_d_state['specific_humidity'].values[:, :, :17] = 10.
# print(state['specific_humidity'])
tendencies, diagnostics = radiation(state)
one_d_tend, one_d_diag = radiation(one_d_state)

print(np.all(tendencies['air_temperature'][0, 5, :] == one_d_tend['air_temperature'][0, 0, :]))
