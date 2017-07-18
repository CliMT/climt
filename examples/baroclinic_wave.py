import climt
from sympl import PlotFunctionMonitor


def plot_function(fig, state):

    ax = fig.add_subplot(1, 1, 1)
    state['surface_air_pressure'].transpose().plot.contourf(
        ax=ax, levels=16)


monitor = PlotFunctionMonitor(plot_function)

dycore = climt.GfsDynamicalCore(number_of_longitudes=198,
                                number_of_latitudes=94,
                                dry_pressure=1e5)
dcmip = climt.DcmipInitialConditions()

my_state = climt.get_default_state([dycore], x=dycore.grid_definition['x'],
                                   y=dycore.grid_definition['y'],
                                   z=dycore.grid_definition['mid_levels'])

my_state['surface_air_pressure'].values[:] = 1e5
dycore(my_state)

out = dcmip(my_state, add_perturbation=True)

my_state.update(out)

for i in range(1000):
    diag, output = dycore(my_state)
    monitor.store(output)
    my_state.update(output)
