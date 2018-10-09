from sympl import PlotFunctionMonitor
from climt import SimplePhysics, get_default_state
from datetime import timedelta


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['specific_humidity'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    ax.set_ylim(1e5, 100.)
    ax.set_xlabel('kg/kg')
    ax.set_ylabel('Pa')
    ax.set_title('Specific Humidity')


monitor = PlotFunctionMonitor(plot_function)
simple_physics = SimplePhysics()
state = get_default_state([simple_physics])
state['eastward_wind'].values[:] = 10.
timestep = timedelta(hours=1)

for i in range(600):

    diagnostics, new_state = simple_physics(state, timestep)
    print('SH Flux:', diagnostics['surface_upward_sensible_heat_flux'].values.item())
    print('LH Flux:', diagnostics['surface_upward_latent_heat_flux'].values.item())
    state.update(diagnostics)
    monitor.store(state)
    state.update(new_state)
