import climt
from sympl import PlotFunctionMonitor, set_constant
from datetime import timedelta
import matplotlib.pyplot as plt


def plot_function(fig, state):

    ax = fig.add_subplot(1, 1, 1)
    CS = ax.contourf(state['longitude'], state['latitude'],
                     state['surface_air_pressure'].to_units('mbar'))
    plt.colorbar(CS)
    ax.set_title('Surface Pressure at: '+str(state['time']))


monitor = PlotFunctionMonitor(plot_function)

set_constant('reference_air_pressure', value=1e5, units='Pa')
dycore = climt.GFSDynamicalCore()
dcmip = climt.DcmipInitialConditions(add_perturbation=True)

grid = climt.get_grid(nx=128, ny=64, nz=20)

my_state = climt.get_default_state([dycore], grid_state=grid)

timestep = timedelta(minutes=10)

out = dcmip(my_state)

my_state.update(out)

for i in range(1000):
    diag, output = dycore(my_state, timestep)
    monitor.store(my_state)
    my_state.update(output)
    my_state['time'] += timestep
