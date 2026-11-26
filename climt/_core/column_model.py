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
import copy
import warnings

import numpy as np
import sympl

from .compile import CompileError, compile_component
from .parameters import with_parameters

__all__ = ("ColumnModel", "RunResult", "Scenarios", "build_column_model")

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
    scenarios : Scenarios, optional
        Many scenarios run as columns; see ``Scenarios``.
    threads
        Not yet supported; must be None.

    Raises
    ------
    CompileError
        For any unsupported component or configuration.
    """
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

    columns = None
    grouped = {}
    if scenarios is not None:
        state, columns = _scenario_state(state, scenarios.n)
        replace, grouped = _resolve_parameters(scenarios)
        diagnostic_components = [replace.get(id(c), c) for c in diagnostic_components]
        tendency_components = [replace.get(id(c), c) for c in tendency_components]
        steppers = [replace.get(id(c), c) for c in steppers]

    groups = (
        (diagnostic_components, ("diagnostic",), "diagnostic_components"),
        (tendency_components, ("tendency", "implicit"), "tendency_components"),
        (steppers, ("stepper",), "steppers"),
    )
    compiled = []
    for components, kinds, label in groups:
        group = []
        for component in components:
            if id(component) in grouped:
                c = _GroupedComponent(component, grouped[id(component)], state, columns)
            else:
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

    if scenarios is not None:
        for name, value in scenarios.inputs.items():
            if name not in raw:
                raise ValueError(
                    "scenario input {} is not read or written by any component".format(name))
            _set_columns(raw[name], specs[name].dims, columns, value, name)

    prescribed_values = {}
    for name in prescribed:
        if name not in raw:
            if name not in state:
                raise CompileError("prescribed quantity {} is not in the state".format(name))
            raise CompileError(
                "prescribed quantity {} is not read or written by any component".format(name))
        prescribed_values[name] = raw[name].copy()

    return ColumnModel(raw, specs, state["time"], timestep, diagnostic_cs, tendency_cs,
                       stepper_cs, order, prescribed_values, columns)


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
                 tendency_components, steppers, order, prescribed, columns=None):
        self.state = state
        self.specs = specs
        self.time = time
        self._timestep = timestep
        # UnytTimeDelta.total_seconds() is a unyt quantity in s: keep a float.
        self._dt = float(getattr(timestep.total_seconds(), "value", timestep.total_seconds()))
        self._diagnostic_components = diagnostic_components
        self._tendency_components = tendency_components
        self._steppers = steppers
        self._order = order
        self._prescribed = prescribed
        self._history = []
        # (horizontal dims, their lengths) when built with scenarios.
        self._columns = columns

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
        """The current state as quantities of the active sympl backend.

        With ``scenario=i``, only column ``i`` (as a one-column state).
        """
        backend = sympl.get_backend()
        out = {"time": self.time}
        for name, array in self.state.items():
            spec = self.specs[name]
            if scenario is not None:
                array = _get_column(array, spec.dims, self._columns, scenario)
            out[name] = backend.create_quantity(
                np.array(array), name=name, units=spec.units, dims=spec.dims)
        return out


class Scenarios(object):
    """Many single-column scenarios, run as the columns of one model.

    Parameters
    ----------
    n : int
        Number of scenarios.
    inputs : dict, optional
        ``{quantity name: array}`` in the example state's units, of shape
        ``(n,)`` (one value per scenario, broadcast over the quantity's other
        dims) or ``(n, *quantity shape without its horizontal dims)``.
    parameters : dict, optional
        ``{(component, name): array of shape (n,)}``; ``name`` is in
        ``type(component).per_column_parameters`` or
        ``type(component).grouped_parameters``.

    Raises
    ------
    ValueError
        For an array of the wrong length or shape; names the quantity or parameter.
    """

    def __init__(self, n, inputs=None, parameters=None):
        if not isinstance(n, (int, np.integer)) or n < 1:
            raise ValueError("the number of scenarios must be a positive integer, got {!r}"
                             .format(n))
        self.n = int(n)
        self.inputs = {}
        for name, value in (inputs or {}).items():
            value = np.asarray(value)
            if value.ndim == 0 or value.shape[0] != self.n:
                raise ValueError(
                    "scenario input {} has shape {}; expected ({}, ...)".format(
                        name, value.shape, self.n))
            self.inputs[name] = value
        self.parameters = {}
        for key, value in (parameters or {}).items():
            component, name = key
            value = np.asarray(value, dtype=np.float64)
            if value.shape != (self.n,):
                raise ValueError(
                    "scenario parameter {} of {} has shape {}; expected ({},)".format(
                        name, type(component).__name__, value.shape, self.n))
            self.parameters[key] = value


def _horizontal_dims(state):
    for name in ("latitude", "longitude"):
        if name in state:
            return tuple(sympl.get_backend().get_dims(state[name]))
    dims = set()
    for name, value in state.items():
        if name != "time":
            dims.update(sympl.get_backend().get_dims(value))
    return tuple(d for d in ("lat", "lon") if d in dims)


def _dim_lengths(state):
    backend = sympl.get_backend()
    lengths = {}
    for name, value in state.items():
        if name != "time":
            lengths.update(zip(backend.get_dims(value), backend.get_shape(value)))
    return lengths


def _map_quantities(state, dim, function):
    """A copy of ``state`` with ``function(array, axis)`` applied along ``dim``."""
    backend = sympl.get_backend()
    out = {}
    for name, value in state.items():
        dims = tuple(backend.get_dims(value)) if name != "time" else ()
        if dim not in dims:
            out[name] = value
            continue
        units = getattr(value, "attrs", {}).get("units")
        array = np.asarray(backend.get_array(
            value, name, units, list(dims), dict(zip(dims, backend.get_shape(value)))))
        out[name] = backend.create_quantity(
            function(array, dims.index(dim)), name=name, units=units, dims=dims)
    return out


def _scenario_state(state, n):
    """The state to compile on, and (horizontal dims, their lengths)."""
    horizontal = _horizontal_dims(state)
    lengths = _dim_lengths(state)
    shape = tuple(lengths[d] for d in horizontal)
    ncol = int(np.prod(shape)) if shape else 1
    if ncol == n:
        return state, (horizontal, shape)
    if ncol == 1:
        if "lon" not in horizontal:
            raise ValueError("cannot expand the state to {} scenarios: it has no lon dim"
                             .format(n))
        state = _map_quantities(state, "lon", lambda a, axis: np.repeat(a, n, axis=axis))
        shape = tuple(n if d == "lon" else 1 for d in horizontal)
        return state, (horizontal, shape)
    raise ValueError("{} scenarios, but the state has {} columns".format(n, ncol))


def _resolve_parameters(scenarios):
    """``({id: component with per-column values}, {id: {attribute: values}})``."""
    by_component = {}
    for (component, name), value in scenarios.parameters.items():
        by_component.setdefault(id(component), (component, {}))[1][name] = value
    replace = {}
    grouped = {}
    for key, (component, values) in by_component.items():
        cls = type(component)
        per_column = tuple(getattr(cls, "per_column_parameters", ()))
        group_names = dict(getattr(cls, "grouped_parameters", {}))
        column_values = {}
        for name, value in values.items():
            if name in group_names:
                grouped.setdefault(key, {})[group_names[name]] = value
            elif name in per_column:
                column_values[name] = value
            else:
                raise CompileError(
                    "{} does not declare {!r} as a per-scenario parameter; declared: {}"
                    .format(cls.__name__, name,
                            ", ".join(per_column + tuple(group_names)) or "none"))
        new = component
        if column_values:
            new = with_parameters(component, column_values, scenarios.n)
        replace[key] = new
        if key in grouped:
            grouped[id(new)] = grouped.pop(key)
    return replace, grouped


def _column_axes(dims, columns):
    horizontal, _ = columns
    return [dims.index(d) for d in horizontal if d in dims]


def _set_columns(array, dims, columns, value, name):
    """Write per-scenario ``value`` into the canonical ``array``, in place."""
    horizontal, shape = columns
    dims = tuple(dims)
    if not all(d in dims for d in horizontal):
        raise ValueError("scenario input {} has no column dims {}".format(name, horizontal))
    axes = _column_axes(dims, columns)
    view = np.moveaxis(array, axes, list(range(len(axes))))
    rest = view.shape[len(axes):]
    if value.shape == (value.shape[0],):
        value = value.reshape(shape + (1,) * len(rest))
    elif value.shape == (value.shape[0],) + rest:
        value = value.reshape(shape + rest)
    else:
        raise ValueError(
            "scenario input {} has shape {}; expected ({},) or ({}, {})".format(
                name, value.shape, value.shape[0], value.shape[0],
                ", ".join(str(r) for r in rest)))
    view[...] = value


def _get_column(array, dims, columns, scenario):
    """Column ``scenario`` of a canonical array, keeping its column dims (length 1)."""
    horizontal, shape = columns
    index = np.unravel_index(scenario, shape)
    out = array
    for d, i in zip(horizontal, index):
        if d in dims:
            axis = tuple(dims).index(d)
            out = np.take(out, [i], axis=axis)
    return np.array(out)


class _GroupedComponent(object):
    """A component called once per group of columns sharing grouped-parameter values.

    Behaves as a ``CompiledComponent`` on the full canonical state.
    """

    def __init__(self, component, values, state, columns):
        horizontal, shape = columns
        varying = [d for d, length in zip(horizontal, shape) if length > 1]
        name = type(component).__name__
        if len(varying) > 1:
            raise CompileError(
                "{}: grouped parameters need the scenarios along one horizontal dim, "
                "but they span {}".format(name, varying))
        self.dim = varying[0] if varying else horizontal[-1]
        attributes = sorted(values)
        keys = list(zip(*(values[a] for a in attributes)))
        order = []
        members = {}
        for i, key in enumerate(keys):
            if key not in members:
                members[key] = []
                order.append(key)
            members[key].append(i)
        if len(order) > 1:
            warnings.warn("{}: grouped parameters {} split the scenarios into {} groups; "
                          "the component is called once per group".format(
                              name, ", ".join(attributes), len(order)), UserWarning)
        self._groups = []
        for key in order:
            indices = np.array(members[key])
            copied = copy.copy(component)
            for attribute, value in zip(attributes, key):
                setattr(copied, attribute, float(value))
            substate = _map_quantities(
                state, self.dim, lambda a, axis: np.take(a, indices, axis=axis))
            self._groups.append((indices, compile_component(copied, substate)))
        first = copy.copy(component)
        for attribute, value in zip(attributes, order[0]):
            setattr(first, attribute, float(value))
        full = compile_component(first, state)
        self.component = component
        self.kind = full.kind
        self.takes_timestep = full.takes_timestep
        self.inputs = full.inputs
        self.tendencies = full.tendencies
        self.diagnostics = full.diagnostics
        self.outputs = full.outputs
        self.extract = full.extract
        # Specs of (first, second) of what __call__ returns.
        if self.kind == "stepper":
            self._specs = (full.diagnostics, full.outputs)
        else:
            self._specs = (full.tendencies, full.diagnostics)

    def __call__(self, raw, timestep=None, time=None):
        results = ({}, {})
        for indices, c in self._groups:
            sub = {}
            for name, spec in c.inputs.items():
                if self.dim in spec.dims:
                    sub[name] = np.take(raw[name], indices, axis=spec.dims.index(self.dim))
                else:
                    sub[name] = raw[name]
            returned = c(sub, timestep=timestep, time=time)
            for out, arrays, specs in zip(results, returned, self._specs):
                for name, array in arrays.items():
                    spec = specs[name]
                    if self.dim not in spec.dims:
                        out.setdefault(name, array)
                        continue
                    if name not in out:
                        out[name] = np.empty(spec.shape, dtype=np.asarray(array).dtype)
                    axis = spec.dims.index(self.dim)
                    np.moveaxis(out[name], axis, 0)[indices] = np.moveaxis(array, axis, 0)
        return results
