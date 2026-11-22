"""Compile a sympl component into a call on raw NumPy arrays.

``compile_component(component, state)`` resolves, once, everything sympl
resolves on every call: wildcard dimensions, unit conversions, dimension
order, broadcasting and aliases. The resulting ``CompiledComponent`` takes
and returns raw arrays in *canonical form*:

* an input keeps the units string and dimension order it has in ``state``;
* a tendency of a quantity in ``state`` has units ``"(<state units>)/s"``
  and the state's dims and shape;
* a diagnostic in ``state`` has the state's units and dims;
* a stepper output (the new value of a quantity) in ``state`` has the
  state's units and dims;
* a tendency, diagnostic or stepper output of a quantity not in ``state`` has the
  component's declared units and dims, ``"*"`` replaced by the wildcard dims
  in the order sympl matched them.

At call time only NumPy array operations run: no pint, unyt, xarray or
``DataArray`` method is touched.
"""
from datetime import timedelta
from typing import NamedTuple

import numpy as np
import sympl
from sympl._core.units import unit_registry
from sympl._core.wildcard import get_wildcard_matches_and_dim_lengths

from .unyt_backend import UnytBackend

__all__ = ("ArraySpec", "CompileError", "CompiledComponent", "compile_component")

_PROBE_DIM = "_compile_probe"

# Component base class -> kind. Checked in order.
_KINDS = (
    (sympl.Stepper, "stepper"),
    (sympl.TendencyComponent, "tendency"),
    (sympl.ImplicitTendencyComponent, "implicit"),
    (sympl.DiagnosticComponent, "diagnostic"),
)
_TAKES_TIMESTEP = ("implicit", "stepper")
_WRAPPERS = (sympl.ScalingWrapper, sympl.UpdateFrequencyWrapper, sympl.TimeDifferencingWrapper)


class CompileError(ValueError):
    """A component or state that cannot be compiled to a raw-array call."""


class ArraySpec(NamedTuple):
    units: str
    dims: tuple
    shape: tuple


def tendency_units(state_units):
    """Canonical units of the tendency of a quantity held in ``state_units``.

    Written without spaces so that both pint and climt's unyt backend parse it.
    """
    if not state_units.strip():
        return "1/s"
    return "({})/s".format(state_units)


class _InputPlan(object):
    """State array (canonical form) -> the array the component receives."""

    __slots__ = ("name", "key", "perm", "expand_shape", "scale", "offset",
                 "broadcast_shape", "flat_shape")

    def apply(self, array):
        a = array
        if self.perm is not None:
            a = a.transpose(self.perm)
        if self.expand_shape is not None:
            a = a.reshape(self.expand_shape)
        fresh = False
        if self.scale is not None:
            a = a * self.scale
            fresh = True
        if self.offset is not None:
            a = a + self.offset
            fresh = True
        if self.broadcast_shape is not None:
            a = np.broadcast_to(a, self.broadcast_shape)
            fresh = False
        # The component may write into its inputs: never hand it a view of raw.
        a = np.asarray(a, order="C") if fresh else np.array(a, order="C")
        if self.flat_shape is not None:
            a = a.reshape(self.flat_shape)
        return a


class _OutputPlan(object):
    """Array returned by the component -> canonical form.

    ``source`` indexes the pair array_call returns (for a diagnostic
    component, ``({}, diagnostics)``); ``key`` is the name in that dict.
    """

    __slots__ = ("name", "source", "key", "wildcard_axis", "wildcard_shape",
                 "perm", "expand_shape", "broadcast_shape", "scale", "offset")

    def apply(self, array):
        a = np.asarray(array)
        if self.wildcard_axis is not None:
            i = self.wildcard_axis
            a = a.reshape(a.shape[:i] + self.wildcard_shape + a.shape[i + 1:])
        if self.perm is not None:
            a = a.transpose(self.perm)
        if self.expand_shape is not None:
            a = a.reshape(self.expand_shape)
        if self.scale is not None:
            a = a * self.scale
        if self.offset is not None:
            a = a + self.offset
        if self.broadcast_shape is not None:
            a = np.array(np.broadcast_to(a, self.broadcast_shape))
        return a


class CompiledComponent(object):
    """A sympl component with units, dimensions and aliases resolved.

    Attributes
    ----------
    component
        The component it was built from.
    kind : str
        ``"tendency"``, ``"implicit"``, ``"diagnostic"`` or ``"stepper"``.
    takes_timestep : bool
        True if ``__call__`` requires ``timestep``.
    inputs, tendencies, diagnostics, outputs : dict
        ``{name: ArraySpec}`` in canonical form. ``outputs`` is empty except
        for steppers.
    """

    def __init__(self, component, kind, inputs, tendencies, diagnostics, outputs,
                 input_plans, first_plans, second_plans):
        self.component = component
        self.kind = kind
        self.takes_timestep = kind in _TAKES_TIMESTEP
        self.inputs = inputs
        self.tendencies = tendencies
        self.diagnostics = diagnostics
        self.outputs = outputs
        self._input_plans = input_plans
        self._first_plans = first_plans
        self._second_plans = second_plans
        self._array_call = component.array_call
        # (first, second) of what __call__ returns.
        if kind == "stepper":
            self._specs = (diagnostics, outputs)
            # sympl's Stepper inserts (new - old) / timestep * time_unit for
            # each output under "<name>_tendency_from_<component name>".
            self._tendency_diagnostics = (
                [(component._get_tendency_name(name), name,
                  component.time_unit_timedelta.total_seconds())
                 for name in component.output_properties]
                if component.tendencies_in_diagnostics else [])
        else:
            self._specs = (tendencies, diagnostics)
            self._tendency_diagnostics = []

    def __repr__(self):
        return "CompiledComponent({}, kind={!r})".format(
            type(self.component).__name__, self.kind)

    def extract(self, state):
        """Copies of the inputs in ``state``, in canonical form."""
        backend = sympl.get_backend()
        raw = {}
        for name, spec in self.inputs.items():
            array = np.array(
                backend.get_array(state[name], name, spec.units, list(spec.dims),
                                  dict(zip(spec.dims, spec.shape))),
                order="C")
            if array.shape != spec.shape:
                raise ValueError(
                    "{}: {} has shape {}, but the component was compiled for "
                    "shape {}".format(type(self.component).__name__, name,
                                      array.shape, spec.shape))
            raw[name] = array
        return raw

    def __call__(self, raw, timestep=None, time=None):
        """Run the component on raw canonical arrays; ``raw`` is not modified.

        Returns ``(tendencies, diagnostics)`` as dicts of raw canonical arrays;
        for a stepper ``(diagnostics, outputs)``, as sympl's Stepper does.
        ``timestep`` is in seconds and is required iff ``takes_timestep``.
        """
        if self.takes_timestep and timestep is None:
            raise TypeError(
                "{} ({} component) requires a timestep in seconds".format(
                    type(self.component).__name__, self.kind))
        state = {plan.key: plan.apply(raw[plan.name]) for plan in self._input_plans}
        if time is not None:
            state["time"] = time
        if self.kind == "diagnostic":
            returned = ({}, self._array_call(state))
        elif self.takes_timestep:
            returned = self._array_call(state, timedelta(seconds=timestep))
        else:
            returned = self._array_call(state)
        if self._tendency_diagnostics:
            diagnostics, new_state = returned
            for out_name, name, time_unit in self._tendency_diagnostics:
                diagnostics[out_name] = (
                    (new_state[name] - state[name]) / timestep * time_unit)
        first = {plan.name: plan.apply(returned[plan.source][plan.key])
                 for plan in self._first_plans}
        second = {plan.name: plan.apply(returned[plan.source][plan.key])
                  for plan in self._second_plans}
        return first, second

    def restore(self, first, second):
        """Wrap the arrays ``__call__`` returned as quantities of the active backend."""
        backend = sympl.get_backend()

        def wrap(arrays, specs):
            return {
                name: backend.create_quantity(
                    array, name=name, units=specs[name].units, dims=specs[name].dims)
                for name, array in arrays.items()
            }

        return wrap(first, self._specs[0]), wrap(second, self._specs[1])


def compile_component(component, state):
    """Resolve ``component``'s units, dimensions and aliases against ``state``.

    Parameters
    ----------
    component : sympl.TendencyComponent, sympl.ImplicitTendencyComponent, sympl.DiagnosticComponent or sympl.Stepper
        A component with ``array_call``. sympl's wrappers are refused.
    state : dict
        An example state; its units and dimension order become the canonical
        form of the compiled component's inputs.

    Returns
    -------
    CompiledComponent

    Raises
    ------
    CompileError
        If the component is not supported, an input is missing from the state,
        or a unit cannot be converted affinely.
    """
    cls = type(component).__name__
    kind = _kind_of(component, cls)
    backend = sympl.get_backend()

    input_properties = component.input_properties
    for name in input_properties:
        if name not in state:
            raise CompileError("{}: input {} is missing from the state".format(cls, name))
    try:
        wildcard_names, dim_lengths = get_wildcard_matches_and_dim_lengths(
            state, input_properties)
    except Exception as err:
        raise CompileError("{}: {}".format(cls, err)) from err

    ctx = _Context(cls, backend, state, wildcard_names, dim_lengths)

    inputs = {}
    input_plans = []
    for name, properties in input_properties.items():
        inputs[name] = ctx.state_spec(name)
        input_plans.append(ctx.input_plan(name, properties))

    tendencies = {}
    tendency_plans = []
    for name, properties in getattr(component, "tendency_properties", {}).items():
        dims = _declared_dims(cls, "tendency", name, properties, input_properties)
        key = _raw_name(name, properties, input_properties)
        spec, plan = ctx.output(name, properties["units"], dims, key, 0, tendency=True)
        tendencies[name] = spec
        tendency_plans.append(plan)

    # A stepper returns (diagnostics, new_state); the others (x, diagnostics).
    diagnostic_source = 0 if kind == "stepper" else 1
    outputs = {}
    output_plans = []
    for name, properties in getattr(component, "output_properties", {}).items():
        dims = _declared_dims(cls, "output", name, properties, input_properties)
        key = _raw_name(name, properties, input_properties)
        spec, plan = ctx.output(name, properties["units"], dims, key, 1)
        outputs[name] = spec
        output_plans.append(plan)

    diagnostics = {}
    diagnostic_plans = []
    added = set(getattr(component, "_added_diagnostic_names", ()))
    for name, properties in component.diagnostic_properties.items():
        if name in added:
            continue
        key = _raw_name(name, properties, input_properties)
        spec, plan = ctx.output(name, properties["units"], properties["dims"], key,
                                diagnostic_source)
        diagnostics[name] = spec
        diagnostic_plans.append(plan)
    if kind != "stepper" and getattr(component, "tendencies_in_diagnostics", False):
        # sympl inserts each tendency, in its tendency units, under
        # "<name>_tendency_from_<component name>".
        for name, properties in component.tendency_properties.items():
            out_name = component._get_tendency_name(name)
            out_properties = component.diagnostic_properties[out_name]
            key = _raw_name(name, properties, input_properties)
            spec, plan = ctx.output(out_name, out_properties["units"],
                                    out_properties["dims"], key, 0)
            diagnostics[out_name] = spec
            diagnostic_plans.append(plan)

    if kind == "stepper":
        first_plans, second_plans = diagnostic_plans, output_plans
    else:
        first_plans, second_plans = tendency_plans, diagnostic_plans
    return CompiledComponent(
        component, kind, inputs, tendencies, diagnostics, outputs,
        input_plans, first_plans, second_plans)


def _declared_dims(cls, what, name, properties, input_properties):
    """Declared dims of a tendency or output; sympl falls back to the input's."""
    dims = properties.get("dims")
    if dims is None:
        if name not in input_properties:
            raise CompileError(
                "{}: {} {} declares no dims and is not an input".format(cls, what, name))
        dims = input_properties[name]["dims"]
    return dims


def _kind_of(component, cls):
    if isinstance(component, _WRAPPERS):
        raise CompileError(
            "{}: sympl wrappers cannot be compiled; compile the wrapped "
            "component instead".format(cls))
    for base, kind in _KINDS:
        if isinstance(component, base):
            if not callable(getattr(component, "array_call", None)):
                raise CompileError("{}: component has no array_call".format(cls))
            if getattr(component, "uses_tracers", False):
                raise CompileError("{}: components that use tracers are not supported".format(cls))
            return kind
    raise CompileError(
        "{}: not a sympl TendencyComponent, ImplicitTendencyComponent, "
        "DiagnosticComponent or Stepper".format(cls))


def _raw_name(name, properties, input_properties):
    """The key under which array_call returns ``name`` (sympl's get_alias_or_name)."""
    if "alias" in properties:
        return properties["alias"]
    if "alias" in input_properties.get(name, {}):
        return input_properties[name]["alias"]
    return name


def _identity_to_none(perm):
    return None if perm == tuple(range(len(perm))) else perm


class _Context(object):
    """What compile_component knows about the state while building plans."""

    def __init__(self, cls, backend, state, wildcard_names, dim_lengths):
        self.cls = cls
        self.backend = backend
        self.state = state
        self.wildcard_names = wildcard_names
        self.dim_lengths = dim_lengths

    def state_spec(self, name):
        value = self.state[name]
        units = getattr(value, "attrs", {}).get("units", "")
        return ArraySpec(units, tuple(self.backend.get_dims(value)),
                         tuple(self.backend.get_shape(value)))

    def has_units(self, name):
        # The unyt backend stores non-numeric data (LandMask's area_type)
        # without units and passes it through unconverted.
        return "units" in getattr(self.state[name], "attrs", {})

    def expand_wildcard(self, name, dims):
        """``dims`` with "*" replaced by the wildcard names; and the index of "*"."""
        dims = list(dims)
        if "*" not in dims:
            return dims, None
        if self.wildcard_names is None:
            raise CompileError(
                "{}: {} has wildcard dims, but no input does".format(self.cls, name))
        i = dims.index("*")
        dims[i:i + 1] = self.wildcard_names
        return dims, i

    def affine(self, name, from_units, to_units):
        """(scale, offset) with ``to = from * scale + offset``; None for 1 and 0."""
        probe = np.array([0.0, 1.0, 2.0])
        try:
            quantity = self.backend.create_quantity(probe, name, from_units, (_PROBE_DIM,))
            f = np.asarray(
                self.backend.get_array(quantity, name, to_units, [_PROBE_DIM], {_PROBE_DIM: 3}),
                dtype=float)
        except Exception as err:
            raise CompileError(
                "{}: cannot convert {} from {!r} to {!r}: {}".format(
                    self.cls, name, from_units, to_units, err)) from err
        offset = f[0]
        scale = f[1] - f[0]
        if not abs(f[2] - (offset + 2 * scale)) <= 1e-9 * max(1.0, abs(f[2])):
            raise CompileError(
                "{}: conversion of {} from {!r} to {!r} is not affine".format(
                    self.cls, name, from_units, to_units))
        return (None if scale == 1.0 else float(scale),
                None if offset == 0.0 else float(offset))

    def has_offset(self, units):
        if isinstance(self.backend, UnytBackend):
            return self.backend._get_cached_unit(units).base_offset != 0
        # pint has no public test; to_base_units fails on sympl's degrees_north.
        return not unit_registry(units)._is_multiplicative

    def input_plan(self, name, properties):
        spec = self.state_spec(name)
        target, i = self.expand_wildcard(name, properties["dims"])
        plan = _InputPlan()
        plan.name = name
        plan.key = properties.get("alias", name)
        present = [d for d in target if d in spec.dims]
        plan.perm = _identity_to_none(tuple(spec.dims.index(d) for d in present))
        expanded = tuple(self.dim_lengths[d] if d in spec.dims else 1 for d in target)
        full = tuple(self.dim_lengths.get(d, 1) for d in target)
        plan.expand_shape = None if len(target) == len(spec.dims) else expanded
        plan.broadcast_shape = None if full == expanded else full
        if i is None:
            plan.flat_shape = None
        else:
            n = len(self.wildcard_names)
            plan.flat_shape = full[:i] + (int(np.prod(full[i:i + n])),) + full[i + n:]
        if self.has_units(name):
            plan.scale, plan.offset = self.affine(name, spec.units, properties["units"])
        else:
            plan.scale = plan.offset = None
        return plan

    def output(self, name, units, dims, key, source, tendency=False):
        """(ArraySpec, _OutputPlan) for an output declared with ``units``, ``dims``."""
        expanded, i = self.expand_wildcard(name, dims)
        plan = _OutputPlan()
        plan.name = name
        plan.source = source
        plan.key = key
        plan.wildcard_axis = i
        if i is not None:
            plan.wildcard_shape = tuple(self.dim_lengths[d] for d in self.wildcard_names)
        plan.perm = plan.expand_shape = plan.broadcast_shape = None

        if name not in self.state:
            unknown = [d for d in expanded if d not in self.dim_lengths]
            if unknown:
                raise CompileError(
                    "{}: dims {} of {} have no length in the state".format(
                        self.cls, unknown, name))
            plan.scale = plan.offset = None
            shape = tuple(self.dim_lengths[d] for d in expanded)
            return ArraySpec(units, tuple(expanded), shape), plan

        state_spec = self.state_spec(name)
        extra = [d for d in expanded if d not in state_spec.dims]
        if extra:
            raise CompileError(
                "{}: {} has dims {} in the state but is returned with dims {}; "
                "broadcast {} in the example state to dims {}".format(
                    self.cls, name, state_spec.dims, tuple(expanded), name,
                    tuple(expanded)))
        lengths = dict(zip(state_spec.dims, state_spec.shape))
        for d in expanded:
            if self.dim_lengths.get(d, lengths[d]) != lengths[d]:
                raise CompileError(
                    "{}: dim {} of {} has length {} in the state, but {} in the "
                    "inputs".format(self.cls, d, name, lengths[d], self.dim_lengths[d]))
        present = [d for d in state_spec.dims if d in expanded]
        plan.perm = _identity_to_none(tuple(expanded.index(d) for d in present))
        if len(present) != len(state_spec.dims):
            plan.expand_shape = tuple(
                lengths[d] if d in expanded else 1 for d in state_spec.dims)
            if plan.expand_shape != state_spec.shape:
                plan.broadcast_shape = state_spec.shape

        if tendency:
            try:
                offset = self.has_offset(state_spec.units)
            except Exception as err:
                raise CompileError("{}: cannot read the units {!r} of {}: {}".format(
                    self.cls, state_spec.units, name, err)) from err
            if offset:
                raise CompileError(
                    "{}: {} has a tendency but its state units {!r} have an "
                    "offset; give it in an absolute unit such as degK".format(
                        self.cls, name, state_spec.units))
            canonical = tendency_units(state_spec.units)
        elif self.has_units(name):
            canonical = state_spec.units
        else:
            plan.scale = plan.offset = None
            return state_spec, plan
        plan.scale, plan.offset = self.affine(name, units, canonical)
        return ArraySpec(canonical, state_spec.dims, state_spec.shape), plan
