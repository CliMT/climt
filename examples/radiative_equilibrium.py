from sympl import (
    AdamsBashforth, PlotFunctionMonitor)
from climt import (
    Frierson06LongwaveOpticalDepth, GrayLongwaveRadiation,
    get_default_state)
from datetime import timedelta


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['air_temperature'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylim(1e5, 100.)
    ax.set_xlabel('Kelvin')
    ax.set_ylabel('Pa')
    ax.set_title('Radiative equilibrium in a Grey Gas atmosphere')
    ax.grid()


monitor = PlotFunctionMonitor(plot_function)
diagnostic = Frierson06LongwaveOpticalDepth()
radiation = GrayLongwaveRadiation()
time_stepper = AdamsBashforth([radiation])
timestep = timedelta(hours=4)

state = get_default_state([radiation, diagnostic])

for i in range(6*7*4*10):
    state.update(diagnostic(state))
    diagnostics, new_state = time_stepper(state, timestep)
    state.update(diagnostics)
    if i % 5 == 0:
        monitor.store(state)
    state.update(new_state)
