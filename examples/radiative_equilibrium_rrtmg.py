from sympl import (
    AdamsBashforth, PlotFunctionMonitor)
from climt import RRTMGShortwave, RRTMGLongwave, get_default_state, get_grid
from datetime import timedelta
import matplotlib.pyplot as plt


def plot_function(fig, state):
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(
        state['air_temperature_tendency_from_shortwave'].values.flatten(),
        state['air_pressure'].values.flatten()/100, '-o', label='SW')
    ax.axes.invert_yaxis()
    ax.plot(
        state['air_temperature_tendency_from_longwave'].values.flatten(),
        state['air_pressure'].values.flatten()/100, '-o', label='LW')
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
        state['air_pressure'].values.flatten()/100, '-o')
    ax.axes.invert_yaxis()

    ax.set_yscale('log')
    ax.set_ylim(1e3, 10.)
    ax.set_title('Temperature')
    ax.grid()
    ax.set_xlabel('K')
    ax.set_yticklabels([])

    plt.suptitle('Radiative Eq. with RRTMG')


monitor = PlotFunctionMonitor(plot_function)
rad_sw = RRTMGShortwave()
rad_lw = RRTMGLongwave()
time_stepper = AdamsBashforth([rad_sw, rad_lw])
timestep = timedelta(hours=3)

grid = get_grid(nx=1, ny=1, nz=30)
state = get_default_state([rad_sw, rad_lw], grid_state=grid)

for i in range(100000):

    diagnostics, new_state = time_stepper(state, timestep)
    state.update(diagnostics)
    if i % 2 == 0:
        monitor.store(state)
    state = new_state
