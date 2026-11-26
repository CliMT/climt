"""build_column_model: one model step from compiled components.

The oracle is the sympl loop of examples/gmd_radiative_convective_unyt.py.
"""
import copy
from datetime import timedelta

import numpy as np
import pytest
import sympl

import climt

N_STEPS = 50
TIMESTEP = timedelta(minutes=5)


@pytest.fixture(autouse=True)
def dataarray_backend():
    sympl.set_backend(sympl.DataArrayBackend())


def set_values(state, name, value):
    state[name].values[:] = value


def gmd_column(convection="EmanuelConvection"):
    """(state, kwargs) for the GMD column of examples/gmd_radiative_convective_unyt.py."""
    convection = getattr(climt, convection)()
    tendency_components = [convection, climt.RRTMGLongwave(), climt.RRTMGShortwave(),
                           climt.SlabSurface()]
    stepper = climt.SimplePhysics()
    state = climt.get_default_state([stepper] + tendency_components)
    set_values(state, "air_temperature", 270)
    for name in ("surface_albedo_for_direct_shortwave", "surface_albedo_for_direct_near_infrared",
                 "surface_albedo_for_diffuse_shortwave"):
        set_values(state, name, 0.5)
    set_values(state, "zenith_angle", np.pi / 2.5)
    set_values(state, "surface_temperature", 300.)
    set_values(state, "ocean_mixed_layer_thickness", 5)
    set_values(state, "area_type", "sea")
    set_values(state, "eastward_wind", 3.)
    kwargs = dict(timestep=TIMESTEP, tendency_components=tendency_components,
                  steppers=[stepper], prescribed=["eastward_wind"])
    return state, kwargs


def dry_column():
    radiation = climt.GrayLongwaveRadiation()
    stepper = climt.DryConvectiveAdjustment()
    state = climt.get_default_state([radiation, stepper])
    set_values(state, "air_temperature", 300.)
    set_values(state, "surface_temperature", 300.)
    kwargs = dict(timestep=TIMESTEP, tendency_components=[radiation], steppers=[stepper],
                  time_stepping="euler")
    return state, kwargs


def sympl_run(state, kwargs, n_steps):
    """The sympl loop of the example, returning the states after each step."""
    state = copy.deepcopy(state)
    if kwargs.get("time_stepping", "adams_bashforth") == "euler":
        order = 1
    else:
        order = kwargs.get("order", 3)
    time_stepper = sympl.AdamsBashforth(*kwargs["tendency_components"], order=order)
    timestep = kwargs["timestep"]
    prescribed = {name: state[name].values.copy() for name in kwargs.get("prescribed", ())}
    states = []
    for _ in range(n_steps):
        diagnostics, state = time_stepper(state, timestep)
        state.update(diagnostics)
        for stepper in kwargs.get("steppers", ()):
            diagnostics, new_state = stepper(state, timestep)
            state.update(diagnostics)
            state.update(new_state)
        state["time"] += timestep
        for name, value in prescribed.items():
            state[name].values[:] = value
        states.append(copy.deepcopy(state))
    return states


def assert_matches(model, expected, step):
    assert model.time == expected["time"]
    for name, array in model.state.items():
        spec = model.specs[name]
        value = expected[name]
        if spec.units and "units" in value.attrs:
            value = value.to_units(spec.units)
        value = np.asarray(value.transpose(*spec.dims).values)
        if np.issubdtype(value.dtype, np.number):
            np.testing.assert_allclose(array, value, rtol=1e-12, atol=0,
                                       err_msg="{} at step {}".format(name, step))
        else:
            np.testing.assert_array_equal(array, value, err_msg=name)


def check_run(make, **overrides):
    state, kwargs = make()
    kwargs.update(overrides)
    expected = sympl_run(state, kwargs, N_STEPS)
    state, kwargs = make()  # fresh component instances
    kwargs.update(overrides)
    model = climt.build_column_model(state, **kwargs)
    for step in range(N_STEPS):
        model.step()
        assert_matches(model, expected[step], step + 1)


@pytest.mark.parametrize("convection", ["EmanuelConvection", "EmanuelConvectionPython"])
@pytest.mark.parametrize("time_stepping,order",
                         [("euler", 1), ("adams_bashforth", 2), ("adams_bashforth", 3)])
def test_steps_match_sympl(convection, time_stepping, order):
    check_run(lambda: gmd_column(convection), time_stepping=time_stepping, order=order)


def test_dry_column_matches_sympl():
    check_run(dry_column)


@pytest.mark.parametrize("bad,match", [
    ({"time_stepping": "leapfrog"}, "leapfrog"),
    ({"order": 5}, "order"),
    ({"steppers": [sympl.UpdateFrequencyWrapper(climt.SimplePhysics(), timedelta(hours=1))]},
     "UpdateFrequencyWrapper"),
])
def test_build_refusals(bad, match):
    state, kwargs = gmd_column()
    kwargs.update(bad)
    with pytest.raises(climt.CompileError, match=match):
        climt.build_column_model(state, **kwargs)


def test_duplicate_diagnostic_is_refused():
    radiation = climt.GrayLongwaveRadiation()
    other = climt.GrayLongwaveRadiation()
    state = climt.get_default_state([radiation])
    with pytest.raises(climt.CompileError,
                       match="downwelling_longwave_flux_in_air.*GrayLongwaveRadiation"):
        climt.build_column_model(state, timestep=TIMESTEP,
                                 tendency_components=[radiation, other])


def test_duplicate_between_tendency_and_diagnostic_component_is_refused():
    radiation = climt.GrayLongwaveRadiation()
    optical_depth = climt.Frierson06LongwaveOpticalDepth()
    state = climt.get_default_state([radiation, optical_depth])
    with pytest.raises(climt.CompileError, match="Frierson06LongwaveOpticalDepth"):
        climt.build_column_model(state, timestep=TIMESTEP, tendency_components=[radiation],
                                 diagnostic_components=[optical_depth,
                                                        climt.Frierson06LongwaveOpticalDepth()])


def test_run_returns_stored_steps():
    state, kwargs = dry_column()
    model = climt.build_column_model(state, **kwargs)
    result = model.run(10, output_every=5, outputs=["air_temperature"])
    shape = model.state["air_temperature"].shape
    assert len(result.times) == 2
    assert result.times[-1] == model.time
    assert result.arrays["air_temperature"].shape == (2,) + shape
    np.testing.assert_array_equal(result.arrays["air_temperature"][-1],
                                  model.state["air_temperature"])
    ds = result.to_xarray()
    assert ds["air_temperature"].dims == ("time", "mid_levels", "lat", "lon")
    assert ds["air_temperature"].attrs["units"] == "degK"


def test_to_state_round_trips():
    state, kwargs = dry_column()
    model = climt.build_column_model(state, **kwargs)
    model.step()
    restored = model.to_state()
    assert restored["time"] == model.time
    np.testing.assert_array_equal(restored["air_temperature"].values,
                                  model.state["air_temperature"])
