import climt
from sympl import PlotFunctionMonitor
import numpy as np
import matplotlib.pyplot as plt


def plot_function(fig, state):

    fig.set_size_inches(10, 5)

    ax = fig.add_subplot(1, 2, 1)
    state['air_temperature'].mean(dim='longitude').transpose().plot.contourf(
        ax=ax, levels=16)
    ax.set_title('Temperature')

    ax = fig.add_subplot(1, 2, 2)
    state['eastward_wind'].mean(dim='longitude').transpose().plot.contourf(
        ax=ax, levels=16)
    ax.set_title('Zonal Wind')

    plt.suptitle('Time: '+str(state['time']))


monitor = PlotFunctionMonitor(plot_function)

dycore = climt.GFSDynamicalCore(number_of_longitudes=128,
                                number_of_latitudes=64)
held_suarez = climt.HeldSuarez()

my_state = climt.get_default_state([dycore], x=dycore.grid_definition['x'],
                                   y=dycore.grid_definition['y'],
                                   mid_levels=dycore.grid_definition['mid_levels'],
                                   interface_levels=dycore.grid_definition['interface_levels'])

my_state['eastward_wind'].values[:] = np.random.randn(*my_state['eastward_wind'].shape)

dycore.prognostics = [held_suarez]

for i in range(10000):
    output, diag = dycore(my_state)
    if (my_state['time'].hour % 2 == 0 and
            my_state['time'].minute == 0):
        print('max. zonal wind: ', np.amax(my_state['eastward_wind'].values))
        monitor.store(my_state)
    my_state.update(output)
    my_state['time'] += dycore._time_step
