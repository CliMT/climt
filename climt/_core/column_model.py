"""A column model stepped on raw arrays: compiled components composed into one step.

``build_column_model`` compiles every component once against an example
state; ``ColumnModel.step`` then reproduces one iteration of the sympl loop
of ``examples/gmd_radiative_convective_unyt.py`` using only NumPy:

1. diagnostic components on state *n*, outputs written into the state;
2. tendency components on state *n*, tendencies of one quantity summed in
   component order;
3. the time-stepping update (Euler or Adams-Bashforth with sympl's
   coefficients and start-up at lower order) of every quantity with a
   tendency, then the tendency components' diagnostics written;
4. each stepper in order: its diagnostics, then its outputs, written;
5. prescribed quantities reset to their build-time values; time advanced.
"""
import numpy as np
import sympl

from .compile import CompileError, compile_component

__all__ = ("ColumnModel", "RunResult", "build_column_model")

_TIME_STEPPING = ("euler", "adams_bashforth")


def build_column_model(state, *, timestep, tendency_components=(), diagnostic_components=(),
                       steppers=(), time_stepping="adams_bashforth", order=3,
                       prescribed=(), scenarios=None, threads=None):
    """Compile components against ``state`` into a ``ColumnModel``.

    Parameters
    ----------
    state : dict
        Example state; its units and dimension order become the canonical form.
    timestep : datetime.timedelta or UnytTimeDelta
    tendency_components, diagnostic_components, steppers : sequence
        sympl components, applied as described in the module docstring.
    time_stepping : {"euler", "adams_bashforth"}
    order : int
        Adams-Bashforth order, 1 to 4.
    prescribed : sequence of str
        Quantities reset to their values in ``state`` after every step.
    scenarios, threads
        Not yet supported; must be None.

    Raises
    ------
    CompileError
        For any unsupported component or configuration.
    """
    if scenarios is not None:
        raise NotImplementedError("scenarios are not supported yet")
    if threads is not None:
        raise NotImplementedError("threads are not supported yet")
    if time_stepping not in _TIME_STEPPING:
        raise CompileError(
            "time_stepping {!r} is not supported; use one of {}".format(
                time_stepping, _TIME_STEPPING))
    if time_stepping == "euler":
        order = 1
    if isinstance(order, float) and order.is_integer():
        order = int(order)
    if not isinstance(order, int) or not 1 <= order <= 4:
        raise CompileError("Adams-Bashforth order must be an integer from 1 to 4, got {!r}"
                           .format(order))

    groups = (
        (diagnostic_components, ("diagnostic",), "diagnostic_components"),
        (tendency_components, ("tendency", "implicit"), "tendency_components"),
        (steppers, ("stepper",), "steppers"),
    )
    compiled = []
    for components, kinds, label in groups:
        group = []
        for component in components:
            c = compile_component(component, state)
            if c.kind not in kinds:
                raise CompileError("{}: a {} component cannot be given in {}".format(
                    type(component).__name__, c.kind, label))
            group.append(c)
        compiled.append(group)
    diagnostic_cs, tendency_cs, stepper_cs = compiled

    # A diagnostic produced by two components.
    producer = {}
    for c in diagnostic_cs + tendency_cs + stepper_cs:
        for name in c.diagnostics:
            if name in producer:
                raise CompileError(
                    "diagnostic {} is produced by both {} and {}".format(
                        name, type(producer[name].component).__name__,
                        type(c.component).__name__))
            producer[name] = c

    specs = {}
    raw = {}
    for c in diagnostic_cs + tendency_cs + stepper_cs:
        for name, array in c.extract(state).items():
            if name not in raw:
                raw[name] = array
                specs[name] = c.inputs[name]
    for c in diagnostic_cs + tendency_cs + stepper_cs:
        for name, spec in list(c.diagnostics.items()) + list(c.outputs.items()):
            if name not in raw:
                raw[name] = np.zeros(spec.shape)
                specs[name] = spec
        for name in c.tendencies:
            if name not in raw:
                raise CompileError("{}: tendency of {}, which is not in the state".format(
                    type(c.component).__name__, name))

    prescribed_values = {}
    for name in prescribed:
        if name not in raw:
            if name not in state:
                raise CompileError("prescribed quantity {} is not in the state".format(name))
            raise CompileError(
                "prescribed quantity {} is not read or written by any component".format(name))
        prescribed_values[name] = raw[name].copy()

    return ColumnModel(raw, specs, state["time"], timestep, diagnostic_cs, tendency_cs,
                       stepper_cs, order, prescribed_values)


class RunResult(object):
    """Stored steps of ``ColumnModel.run``.

    Attributes
    ----------
    times : list of datetime
    arrays : dict
        ``{name: ndarray of shape (n_stored, *canonical shape)}``.
    specs : dict
        ``{name: ArraySpec}``.
    """

    def __init__(self, times, arrays, specs):
        self.times = times
        self.arrays = arrays
        self.specs = specs

    def to_xarray(self):
        """An ``xarray.Dataset`` with dims ``("time", *spec.dims)``; units in attrs."""
        import xarray as xr
        data_vars = {
            name: (("time",) + tuple(self.specs[name].dims), array,
                   {"units": self.specs[name].units})
            for name, array in self.arrays.items()
        }
        return xr.Dataset(data_vars, coords={"time": list(self.times)})


class ColumnModel(object):
    """A model stepped on raw canonical arrays; built by ``build_column_model``.

    Attributes
    ----------
    state : dict
        Raw arrays in canonical form of every quantity any component reads or writes.
    specs : dict
        ``{name: ArraySpec}`` of ``state``.
    time : datetime.datetime
    """

    def __init__(self, state, specs, time, timestep, diagnostic_components,
                 tendency_components, steppers, order, prescribed):
        self.state = state
        self.specs = specs
        self.time = time
        self._timestep = timestep
        self._dt = timestep.total_seconds()
        self._diagnostic_components = diagnostic_components
        self._tendency_components = tendency_components
        self._steppers = steppers
        self._order = order
        self._prescribed = prescribed
        self._history = []

    def step(self):
        """Advance the model one timestep, in place."""
        state = self.state
        dt = self._dt
        time = self.time

        for c in self._diagnostic_components:
            _, diagnostics = c(state, time=time)
            state.update(diagnostics)

        tendencies = {}
        tendency_diagnostics = {}
        for c in self._tendency_components:
            if c.takes_timestep:
                tnd, diagnostics = c(state, timestep=dt, time=time)
            else:
                tnd, diagnostics = c(state, time=time)
            for name, value in tnd.items():
                if name in tendencies:
                    tendencies[name] += value
                else:
                    tendencies[name] = value.copy()
            tendency_diagnostics.update(diagnostics)

        if self._tendency_components:
            self._history.append(tendencies)
            state.update(self._update(state, self._history, dt))
            if len(self._history) == self._order:
                self._history.pop(0)
        state.update(tendency_diagnostics)

        for c in self._steppers:
            diagnostics, outputs = c(state, timestep=dt, time=time)
            state.update(diagnostics)
            state.update(outputs)

        for name, value in self._prescribed.items():
            state[name] = value.copy()
        self.time = self.time + self._timestep

    def _update(self, state, history, dt):
        # The expressions of sympl/_components/timesteppers.py, in its order.
        order = min(self._order, len(history))
        t = history
        new = {}
        for key in t[0]:
            if order == 1:
                new[key] = state[key] + t[-1][key] * dt
            elif order == 2:
                new[key] = state[key] + dt * (1.5 * t[-1][key] - 0.5 * t[-2][key])
            elif order == 3:
                new[key] = state[key] + dt * (
                    23. / 12 * t[-1][key] - 4. / 3 * t[-2][key] + 5. / 12 * t[-3][key])
            else:
                new[key] = state[key] + dt * (
                    55. / 24 * t[-1][key] - 59. / 24 * t[-2][key]
                    + 37. / 24 * t[-3][key] - 3. / 8 * t[-4][key])
        return new

    def run(self, n_steps, output_every=1, outputs=None):
        """Step ``n_steps`` times, storing every ``output_every``-th step.

        ``outputs`` names the quantities to store (default: all of ``state``).
        """
        names = list(self.state) if outputs is None else list(outputs)
        for name in names:
            if name not in self.state:
                raise KeyError("{} is not in the model state".format(name))
        times = []
        stored = {name: [] for name in names}
        for i in range(n_steps):
            self.step()
            if (i + 1) % output_every == 0:
                times.append(self.time)
                for name in names:
                    stored[name].append(np.array(self.state[name]))
        arrays = {
            name: (np.stack(values) if values
                   else np.empty((0,) + self.specs[name].shape))
            for name, values in stored.items()
        }
        return RunResult(times, arrays, {name: self.specs[name] for name in names})

    def to_state(self, scenario=None):
        """The current state as quantities of the active sympl backend."""
        if scenario is not None:
            raise NotImplementedError("scenarios are not supported yet")
        backend = sympl.get_backend()
        out = {"time": self.time}
        for name, array in self.state.items():
            spec = self.specs[name]
            out[name] = backend.create_quantity(
                np.array(array), name=name, units=spec.units, dims=spec.dims)
        return out
