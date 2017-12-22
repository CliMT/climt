import climt
from sympl import PlotFunctionMonitor
from datetime import timedelta
import numpy as np


def plot_function(fig, state):

    ax = fig.add_subplot(1, 1, 1)
    state['zenith_angle'].transpose().plot.contourf(
        ax=ax, levels=16, robust=True)

    fig.suptitle('Zenith Angle at time: '+str(state['time']))
#    plt.savefig(str(state['time'])+'.png')


monitor = PlotFunctionMonitor(plot_function)

instellation = climt.Instellation()

y_grid = dict(label='latitude',
              values=np.linspace(-90, 90, 100), units='degrees_north')
x_grid = dict(label='longitude',
              values=np.linspace(0, 360, 100), units='degrees_east')

state = climt.get_default_state([instellation], x=x_grid, y=y_grid)

time_step = timedelta(hours=6)

for i in range(8000):
    diag = instellation(state)
    state.update(diag)
    monitor.store(state)
    state['time'] += time_step
