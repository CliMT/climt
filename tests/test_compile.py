"""compile_component: one-time resolution of units, dimensions and aliases.

The oracle throughout is sympl's own ``component(state)`` call.
"""
import copy
from datetime import timedelta

import numpy as np
import pytest
import sympl
import unyt

import climt
from climt import _components

DT = 900.0
REQUIRED = ["GrayLongwaveRadiation", "CorkLongwaveRadiation", "HeldSuarez",
            "EmanuelConvectionPython", "Frierson06LongwaveOpticalDepth",
            "RRTMGLongwave", "RRTMGShortwave", "EmanuelConvection", "SlabSurface"]
COMPILABLE = (sympl.TendencyComponent, sympl.ImplicitTendencyComponent, sympl.DiagnosticComponent)
OTHERS = sorted(c.__name__ for c in _components.__all__
                if issubclass(c, COMPILABLE) and c.__name__ not in REQUIRED)
STEPPERS = ["SimplePhysics", "GridScaleCondensation", "SimpleBoundaryLayer", "DryConvectiveAdjustment"]
ICE_STEPPERS = ["IceSheet", "LandIce", "SeaIce"]  # an input that grows dims
OTHER_STEPPERS = sorted(c.__name__ for c in _components.__all__
                        if issubclass(c, sympl.Stepper)
                        and c.__name__ not in STEPPERS + ICE_STEPPERS)


@pytest.fixture(autouse=True)
def dataarray_backend():
    # Other test modules set UnytBackend at import time; conftest restores
    # whatever was active, so pin the default backend here.
    sympl.set_backend(sympl.DataArrayBackend())


def make(name, **kwargs):
    component = getattr(climt, name)(**kwargs)
    grid = climt.get_grid(nx=3, ny=2, nz=10)
    state = climt.get_default_state([component], grid_state=grid)
    return component, state


def sympl_call(component, state):
    state = copy.deepcopy(state)
    if isinstance(component, sympl.DiagnosticComponent):
        return {}, component(state)
    if isinstance(component, (sympl.ImplicitTendencyComponent, sympl.Stepper)):
        return component(state, timedelta(seconds=DT))
    return component(state)


def compiled_call(component, state):
    compiled = climt.compile_component(component, state)
    raw = compiled.extract(state)
    before = {name: value.copy() for name, value in raw.items()}
    kwargs = {"time": state["time"]}
    if compiled.takes_timestep:
        kwargs["timestep"] = DT
    first, second = compiled(raw, **kwargs)
    assert raw.keys() == before.keys()
    for name in raw:
        np.testing.assert_array_equal(raw[name], before[name])
    return compiled, compiled.restore(first, second)


def assert_same(expected, actual):
    assert set(expected) == set(actual)
    for name, value in actual.items():
        assert set(expected[name].dims) == set(value.dims), name
        converted = expected[name]
        if "units" in value.attrs:  # the unyt backend keeps strings without units
            converted = converted.to_units(value.attrs["units"])
        converted = converted.transpose(*value.dims)
        actual_values, expected_values = np.asarray(value.values), np.asarray(converted.values)
        if np.issubdtype(expected_values.dtype, np.number):
            np.testing.assert_allclose(actual_values, expected_values, rtol=1e-12, atol=0, err_msg=name)
        else:  # e.g. LandMask's area_type strings
            np.testing.assert_array_equal(actual_values, expected_values, err_msg=name)


def check_against_sympl(component, state):
    expected = sympl_call(component, state)
    compiled, actual = compiled_call(component, state)
    assert_same(expected[0], actual[0])
    assert_same(expected[1], actual[1])
    for name in component.input_properties:
        assert compiled.inputs[name] == climt._core.compile.ArraySpec(
            state[name].attrs.get("units", ""), tuple(state[name].dims), tuple(state[name].shape))
    for name, spec in compiled.tendencies.items():
        assert spec.dims == tuple(state[name].dims)
        assert spec.shape == tuple(state[name].shape)
    for name, spec in compiled.outputs.items():
        if name not in state:  # declared units and dims
            continue
        assert spec.units == state[name].attrs.get("units", "")
        assert spec.dims == tuple(state[name].dims)
        assert spec.shape == tuple(state[name].shape)
    return compiled


@pytest.mark.parametrize("name", REQUIRED)
def test_matches_sympl(name):
    component, state = make(name)
    check_against_sympl(component, state)


@pytest.mark.parametrize("name", REQUIRED)
def test_matches_sympl_under_the_unyt_backend(name):
    sympl.set_backend(climt.UnytBackend())
    component, state = make(name)
    check_against_sympl(component, state)


@pytest.mark.parametrize("name", OTHERS)
def test_other_components_match_sympl(name):
    try:
        component, state = make(name)
        sympl_call(component, state)
    except Exception as err:
        pytest.skip("{} cannot be built or run with defaults: {!r}".format(name, err))
    check_against_sympl(component, state)


@pytest.mark.parametrize("name", STEPPERS)
def test_stepper_matches_sympl(name):
    component, state = make(name)
    compiled = check_against_sympl(component, state)
    assert compiled.kind == "stepper"
    assert compiled.takes_timestep
    assert set(compiled.outputs) == set(component.output_properties)


@pytest.mark.parametrize("name", STEPPERS)
def test_stepper_matches_sympl_under_the_unyt_backend(name):
    sympl.set_backend(climt.UnytBackend())
    component, state = make(name)
    check_against_sympl(component, state)


@pytest.mark.parametrize("name", OTHER_STEPPERS)
def test_other_steppers_match_sympl(name):
    try:
        component, state = make(name)
        sympl_call(component, state)
    except Exception as err:
        pytest.skip("{} cannot be built or run with defaults: {!r}".format(name, err))
    check_against_sympl(component, state)


@pytest.mark.parametrize("name", ICE_STEPPERS)
def test_quantity_that_grows_dims_is_refused(name):
    component, state = make(name)
    with pytest.raises(climt.CompileError, match="height_on_ice_interface_levels") as info:
        climt.compile_component(component, state)
    assert name in str(info.value)
    assert "broadcast" in str(info.value)


@pytest.mark.parametrize("name", ICE_STEPPERS)
def test_quantity_broadcast_in_the_state_matches_sympl(name):
    component, state = make(name)
    quantity = "height_on_ice_interface_levels"
    old = state[quantity]
    values = np.broadcast_to(
        old.values[:, None, None],
        old.shape + (2, 3))  # make() uses ny=2, nx=3
    state[quantity] = sympl.DataArray(
        np.array(values), dims=("ice_interface_levels", "lat", "lon"), attrs=dict(old.attrs))
    try:
        sympl_call(component, state)
    except Exception as err:
        pytest.skip("sympl cannot run {} on the broadcast state: {!r}".format(name, err))
    check_against_sympl(component, state)


def test_stepper_units_and_dimension_order_are_resolved():
    component, state = resolved_state("SimplePhysics")
    compiled = check_against_sympl(component, state)
    assert compiled.outputs["air_temperature"].dims == ("lon", "lat", "mid_levels")


def test_stepper_tendencies_in_diagnostics_match_sympl():
    base = climt.DryConvectiveAdjustment
    isolated = type("IsolatedDCA", (base,), {
        "diagnostic_properties": copy.deepcopy(base.diagnostic_properties),
        "input_properties": copy.deepcopy(base.input_properties)})
    component = isolated(tendencies_in_diagnostics=True)
    state = climt.get_default_state([component], grid_state=climt.get_grid(nx=3, ny=2, nz=10))
    compiled = check_against_sympl(component, state)
    assert "air_temperature_tendency_from_IsolatedDCA" in compiled.diagnostics


def resolved_state(name):
    component, state = make(name)
    state["air_pressure"] = state["air_pressure"].to_units("hPa")
    state["air_pressure_on_interface_levels"] = (
        state["air_pressure_on_interface_levels"].to_units("bar"))
    state["air_temperature"] = state["air_temperature"].transpose("lon", "lat", "mid_levels")
    return component, state


@pytest.mark.parametrize("name", ["GrayLongwaveRadiation", "EmanuelConvectionPython", "CorkLongwaveRadiation"])
def test_units_and_dimension_order_are_resolved(name):
    component, state = resolved_state(name)
    compiled = check_against_sympl(component, state)
    assert compiled.inputs["air_pressure"].units == "hPa"
    assert compiled.inputs["air_pressure_on_interface_levels"].units == "bar"
    assert compiled.tendencies["air_temperature"].dims == ("lon", "lat", "mid_levels")


def test_tendencies_in_diagnostics_match_sympl():
    # sympl adds the tendency diagnostics to diagnostic_properties in place;
    # a subclass with its own dict keeps that out of GrayLongwaveRadiation.
    base = climt.GrayLongwaveRadiation
    isolated = type("IsolatedGray", (base,), {
        "diagnostic_properties": copy.deepcopy(base.diagnostic_properties)})
    component = isolated(tendencies_in_diagnostics=True)
    state = climt.get_default_state([component], grid_state=climt.get_grid(nx=3, ny=2, nz=10))
    compiled = check_against_sympl(component, state)
    assert "air_temperature_tendency_from_IsolatedGray" in compiled.diagnostics


def test_offset_unit_input_is_converted():
    component, state = make("GrayLongwaveRadiation")
    state["surface_temperature"] = state["surface_temperature"].to_units("degC")
    compiled = check_against_sympl(component, state)
    assert compiled.inputs["surface_temperature"].units == "degC"


def test_offset_unit_with_a_tendency_is_refused():
    component, state = make("GrayLongwaveRadiation")
    state["air_temperature"] = state["air_temperature"].to_units("degC")
    with pytest.raises(climt.CompileError, match="air_temperature") as info:
        climt.compile_component(component, state)
    assert "degK" in str(info.value)


def test_refusals():
    wrapped = sympl.TimeDifferencingWrapper(climt.SimplePhysics())
    with pytest.raises(climt.CompileError, match="TimeDifferencingWrapper"):
        climt.compile_component(wrapped, {})

    component, state = make("GrayLongwaveRadiation")
    del state["air_temperature"]
    with pytest.raises(climt.CompileError, match="air_temperature"):
        climt.compile_component(component, state)

    component, state = make("EmanuelConvectionPython")
    compiled = climt.compile_component(component, state)
    with pytest.raises(TypeError, match="timestep"):
        compiled(compiled.extract(state), time=state["time"])


def test_compile_error_is_a_value_error():
    assert issubclass(climt.CompileError, ValueError)


class WritesItsInput(sympl.TendencyComponent):
    input_properties = {"air_temperature": {"dims": ["*", "mid_levels"], "units": "degK", "alias": "T"}}
    tendency_properties = {"air_temperature": {"units": "degK s^-1"}}
    diagnostic_properties = {}

    def array_call(self, state):
        state["T"] += 1.0
        return {"T": state["T"] * 1e-6}, {}  # sympl keys the tendency by the input alias


def test_raw_is_not_modified_by_a_component_that_writes_its_inputs():
    component = WritesItsInput()
    state = climt.get_default_state([component], grid_state=climt.get_grid(nx=3, ny=2, nz=10))
    check_against_sympl(component, state)  # asserts raw is unchanged


# Components whose own array_call reads no sympl constant (several call
# get_constant, which converts units, inside array_call).
@pytest.mark.parametrize("name", ["EmanuelConvectionPython", "RRTMGLongwave"])
def test_call_uses_no_unit_machinery(monkeypatch, name):
    # hPa, bar and transposed inputs; RRTMGLongwave also a degC input and a
    # degK day^-1 tendency: unit conversions and transposes are all active.
    component, state = resolved_state(name)
    if "surface_temperature" in state:
        state["surface_temperature"] = state["surface_temperature"].to_units("degC")
    compiled = climt.compile_component(component, state)
    raw = compiled.extract(state)
    kwargs = {"time": state["time"], "timestep": DT}
    expected = compiled(raw, **kwargs)

    def forbidden(*args, **kwargs):
        raise AssertionError("unit machinery called during a compiled call")

    monkeypatch.setattr(sympl.DataArray, "to_units", forbidden)
    monkeypatch.setattr(unyt.unyt_array, "to", forbidden)
    tendencies, diagnostics = compiled(raw, **kwargs)
    for name in expected[0]:
        np.testing.assert_array_equal(tendencies[name], expected[0][name])
    for name in expected[1]:
        np.testing.assert_array_equal(diagnostics[name], expected[1][name])


def test_works_with_the_unyt_backend():
    sympl.set_backend(climt.UnytBackend())
    component, state = make("GrayLongwaveRadiation")
    assert isinstance(state["air_temperature"], climt.UnytStateContainer)
    expected_tendencies, _ = component(copy.deepcopy(state))
    compiled, (tendencies, _) = compiled_call(component, state)
    actual = tendencies["air_temperature"]
    assert isinstance(actual, climt.UnytStateContainer)
    assert actual.dims == state["air_temperature"].dims
    expected = expected_tendencies["air_temperature"].to_units(
        compiled.tendencies["air_temperature"].units).transpose(*actual.dims)
    np.testing.assert_allclose(actual.values, expected.values, rtol=1e-12, atol=0)
