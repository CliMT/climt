import climt
from sympl import PlotFunctionMonitor
from datetime import timedelta


def plot_function(fig, state):

    ax = fig.add_subplot(1, 1, 1)
    ax.contourf(state['longitude'], state['latitude'],
                state['zenith_angle'])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    fig.suptitle('Zenith Angle at time: '+str(state['time']))


monitor = PlotFunctionMonitor(plot_function)

instellation = climt.Instellation()

state = climt.get_default_state([instellation],
                                grid_state=climt.get_grid(nx=100, ny=100,
                                                          latitude_grid='regular'))

time_step = timedelta(hours=6)

for i in range(8000):
    diag = instellation(state)
    state.update(diag)
    monitor.store(state)
    state['time'] += time_step
