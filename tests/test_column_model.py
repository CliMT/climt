"""build_column_model: one model step from compiled components.

The oracle is the sympl loop of examples/gmd_radiative_convective_unyt.py.
"""
import copy
import warnings
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


def test_dry_column_matches_sympl_under_the_unyt_backend():
    sympl.set_backend(climt.UnytBackend())
    state, kwargs = dry_column()
    assert isinstance(state["air_temperature"], climt.UnytStateContainer)
    check_run(dry_column, timestep=climt.UnytTimeDelta(seconds=300))


# Scenarios: many single-column scenarios as the columns of one model.

SURFACE_TEMPERATURES = [295., 298., 301., 304.]
ALPHAS = [0.1, 0.2, 0.2, 0.3]
DRAG = [7e-4] * 4


def scenario_column():
    """The GMD column with EmanuelConvectionPython; (state, kwargs, components by role)."""
    state, kwargs = gmd_column("EmanuelConvectionPython")
    convection, _, radiation_sw, _ = kwargs["tendency_components"]
    roles = {"convection": convection, "simple_physics": kwargs["steppers"][0],
             "radiation_sw": radiation_sw}
    return state, kwargs, roles


def alone(surface_temperature=300., alpha=None, drag=None, solar=None):
    """(state, kwargs) of one scenario run alone, with scalar parameters."""
    state, kwargs, roles = scenario_column()
    set_values(state, "surface_temperature", surface_temperature)
    if alpha is not None:
        c = roles["convection"]
        c._params = c._params._replace(ALPHA=alpha)
    if drag is not None:
        c = roles["simple_physics"]
        c._params = c._params._replace(base_momentum_drag_coefficient=drag)
    if solar is not None:
        roles["radiation_sw"]._solar_const = solar
    return state, kwargs


def column(model, name, i):
    from climt._core.column_model import _get_column
    return _get_column(model.state[name], model.specs[name].dims, model._columns, i)


def assert_column_matches(model, i, expected, step):
    assert model.time == expected["time"]
    for name in model.state:
        spec = model.specs[name]
        value = expected[name]
        if spec.units and "units" in value.attrs:
            value = value.to_units(spec.units)
        value = np.asarray(value.transpose(*spec.dims).values)
        actual = column(model, name, i)
        if np.issubdtype(value.dtype, np.number):
            np.testing.assert_allclose(actual, value, rtol=1e-12, atol=0,
                                       err_msg="{} of scenario {} at step {}".format(
                                           name, i, step))
        else:
            np.testing.assert_array_equal(actual, value, err_msg=name)


def scenario_model(surface_temperatures=SURFACE_TEMPERATURES, alphas=ALPHAS, solar=None):
    state, kwargs, roles = scenario_column()
    parameters = {(roles["convection"], "ALPHA"): alphas,
                  (roles["simple_physics"], "base_momentum_drag_coefficient"): DRAG}
    if solar is not None:
        parameters = {(roles["radiation_sw"], "stellar_irradiance"): solar}
    scenarios = climt.Scenarios(
        len(surface_temperatures), inputs={"surface_temperature": surface_temperatures},
        parameters=parameters)
    return climt.build_column_model(state, scenarios=scenarios, **kwargs)


def check_scenarios(model, runs):
    expected = [sympl_run(state, kwargs, N_STEPS) for state, kwargs in runs]
    for step in range(N_STEPS):
        model.step()
        for i, states in enumerate(expected):
            assert_column_matches(model, i, states[step], step + 1)


def test_scenarios_match_sympl_runs_alone():
    model = scenario_model()
    check_scenarios(model, [alone(t, alpha=a, drag=d)
                            for t, a, d in zip(SURFACE_TEMPERATURES, ALPHAS, DRAG)])


def test_identical_scenarios_give_identical_columns():
    model = scenario_model(surface_temperatures=[295., 298., 298., 304.])
    model.run(N_STEPS)
    for name in model.state:
        np.testing.assert_array_equal(column(model, name, 1), column(model, name, 2),
                                      err_msg=name)


def test_grouped_stellar_irradiance_matches_and_warns():
    solar = [1360., 1360., 1300., 1300.]
    # Not pytest.warns: it would turn a FortranUnavailableError (pure install)
    # into "DID NOT WARN" before conftest can report it as a skip.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = scenario_model(solar=solar)
    assert any(issubclass(w.category, UserWarning) and "2 groups" in str(w.message) for w in caught)
    check_scenarios(model, [alone(t, solar=s) for t, s in zip(SURFACE_TEMPERATURES, solar)])
    # The solar constant reaches the Fortran: top-of-atmosphere insolation scales with it.
    toa = [column(model, "downwelling_shortwave_flux_in_air", i).max() for i in range(4)]
    np.testing.assert_allclose(toa[2] / toa[0], 1300. / 1360., rtol=1e-9)


def test_continuing_a_run():
    model = scenario_model()
    model.run(25)
    model.run(25)
    fresh = scenario_model()
    fresh.run(50)
    assert model.time == fresh.time
    for name in fresh.state:
        np.testing.assert_array_equal(model.state[name], fresh.state[name], err_msg=name)


def test_to_state_of_a_scenario():
    model = scenario_model()
    model.step()
    restored = model.to_state(scenario=2)
    assert restored["surface_temperature"].shape == (1, 1)
    assert restored["air_temperature"].dims == ("mid_levels", "lat", "lon")
    np.testing.assert_array_equal(restored["air_temperature"].values,
                                  column(model, "air_temperature", 2))


def test_scenario_count_must_match_the_columns():
    state, kwargs, _ = scenario_column()
    state = climt._core.column_model._scenario_state(state, 4)[0]
    with pytest.raises(ValueError, match="3 scenarios.*4 columns"):
        climt.build_column_model(state, scenarios=climt.Scenarios(3), **kwargs)


@pytest.mark.parametrize("scenarios,match", [
    (lambda c: climt.Scenarios(4, inputs={"surface_temperature": np.ones(3)}),
     "surface_temperature"),
    (lambda c: climt.Scenarios(4, inputs={"air_temperature": np.ones((4, 3))}),
     "air_temperature"),
    (lambda c: climt.Scenarios(4, parameters={(c, "ALPHA"): np.ones((4, 2))}), "ALPHA"),
    (lambda c: climt.Scenarios(4, parameters={(c, "IPBL"): np.ones(4)}), "IPBL"),
])
def test_scenario_errors(scenarios, match):
    state, kwargs, roles = scenario_column()
    with pytest.raises(ValueError, match=match):
        climt.build_column_model(state, scenarios=scenarios(roles["convection"]), **kwargs)


# Threads: numba's thread count during run.

def test_threads_default_from_columns():
    numba = pytest.importorskip("numba")
    from climt._core.parameters import PARALLEL_MIN_COLUMNS
    state, kwargs = dry_column()
    assert climt.build_column_model(state, **kwargs).threads == 1
    model = climt.build_column_model(
        state, scenarios=climt.Scenarios(PARALLEL_MIN_COLUMNS), **kwargs)
    assert model.threads == numba.config.NUMBA_NUM_THREADS


def spy_model(threads, raise_on=None):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < threads:
        pytest.skip("needs {} numba threads".format(threads))
    state, kwargs = dry_column()
    model = climt.build_column_model(state, threads=threads, **kwargs)
    seen = []
    step = model.step

    def spy():
        seen.append(numba.get_num_threads())
        if raise_on is not None and len(seen) == raise_on:
            raise RuntimeError("step {}".format(raise_on))
        step()
    model.step = spy
    return numba, model, seen


def test_threads_set_during_run_and_restored():
    numba, model, seen = spy_model(2)
    numba.set_num_threads(1)
    model.run(3)
    assert seen == [2, 2, 2]
    assert numba.get_num_threads() == 1


def test_threads_restored_after_exception():
    numba, model, seen = spy_model(2, raise_on=3)
    numba.set_num_threads(1)
    with pytest.raises(RuntimeError, match="step 3"):
        model.run(5)
    assert numba.get_num_threads() == 1


@pytest.mark.parametrize("extra", [None, 1])
def test_threads_out_of_range(extra):
    numba = pytest.importorskip("numba")
    threads = 0 if extra is None else numba.config.NUMBA_NUM_THREADS + extra
    state, kwargs = dry_column()
    with pytest.raises(ValueError, match="threads"):
        climt.build_column_model(state, threads=threads, **kwargs)


def test_with_parameters_is_exported():
    assert climt.with_parameters is climt._core.parameters.with_parameters


# Final-review fixes.

def dry_scenarios(**scenario_kwargs):
    state, kwargs = dry_column()
    return climt.build_column_model(state, scenarios=climt.Scenarios(3, **scenario_kwargs),
                                    **kwargs)


def test_scenario_input_with_units_is_converted():
    unyt = pytest.importorskip("unyt")
    celsius = dry_scenarios(
        inputs={"surface_temperature": unyt.unyt_array([20., 25., 30.], "degC")})
    kelvin = dry_scenarios(inputs={"surface_temperature": [293.15, 298.15, 303.15]})
    np.testing.assert_allclose(celsius.state["surface_temperature"],
                               kelvin.state["surface_temperature"], rtol=1e-12, atol=0)


def test_scenario_input_with_dataarray_units_is_converted():
    value = sympl.DataArray(np.array([20., 25., 30.]), dims=["s"], attrs={"units": "degC"})
    celsius = dry_scenarios(inputs={"surface_temperature": value})
    kelvin = dry_scenarios(inputs={"surface_temperature": [293.15, 298.15, 303.15]})
    np.testing.assert_allclose(celsius.state["surface_temperature"],
                               kelvin.state["surface_temperature"], rtol=1e-12, atol=0)


def test_scenario_input_with_inconvertible_units_is_refused():
    unyt = pytest.importorskip("unyt")
    with pytest.raises(ValueError, match="surface_temperature"):
        dry_scenarios(inputs={"surface_temperature": unyt.unyt_array([1., 2., 3.], "m")})


def test_scenario_parameter_with_units_is_refused():
    unyt = pytest.importorskip("unyt")
    sw = climt.RRTMGShortwave()
    with pytest.raises(ValueError, match="stellar_irradiance"):
        climt.Scenarios(3, parameters={
            (sw, "stellar_irradiance"): unyt.unyt_array([1360., 1361., 1362.], "W/m**2")})


def test_mutating_caller_array_after_build_has_no_effect():
    temperatures = np.array([290., 300., 310.])
    model = dry_scenarios(inputs={"surface_temperature": temperatures})
    before = model.state["surface_temperature"].copy()
    temperatures[:] = 0.
    model.step()
    fresh = dry_scenarios(inputs={"surface_temperature": [290., 300., 310.]})
    fresh.step()
    np.testing.assert_array_equal(model.state["air_temperature"],
                                  fresh.state["air_temperature"])
    assert before.min() == 290.


def test_scenario_arrays_are_copies():
    temperatures = np.array([290., 300., 310.])
    scenarios = climt.Scenarios(3, inputs={"surface_temperature": temperatures})
    temperatures[:] = 0.
    np.testing.assert_array_equal(scenarios.inputs["surface_temperature"], [290., 300., 310.])


def test_grouped_and_per_column_parameters_on_one_component_are_refused():
    state, kwargs, roles = scenario_column()
    sw = roles["radiation_sw"]
    cls = type(sw)
    original = cls.__dict__.get("per_column_parameters")
    cls.per_column_parameters = ("fake",)
    try:
        scenarios = climt.Scenarios(4, parameters={(sw, "stellar_irradiance"): np.ones(4),
                                                   (sw, "fake"): np.ones(4)})
        with pytest.raises(climt.CompileError, match="grouped and per-column"):
            climt.build_column_model(state, scenarios=scenarios, **kwargs)
    finally:
        if original is None:
            del cls.per_column_parameters
        else:
            cls.per_column_parameters = original


def test_slab_surface_with_ekman_refused_with_scenarios():
    state, kwargs = gmd_column()
    kwargs["tendency_components"][3] = climt.SlabSurface(include_ekman=True)
    with pytest.raises(climt.CompileError, match="SlabSurface.*Ekman"):
        climt.build_column_model(state, scenarios=climt.Scenarios(4), **kwargs)


@pytest.mark.parametrize("order", [True, False])
def test_bool_order_is_refused(order):
    state, kwargs = dry_column()
    kwargs["time_stepping"] = "adams_bashforth"
    with pytest.raises(climt.CompileError, match="order"):
        climt.build_column_model(state, order=order, **kwargs)


@pytest.mark.parametrize("args,match", [((-1,), "n_steps"), ((2, 0), "output_every")])
def test_run_argument_validation(args, match):
    state, kwargs = dry_column()
    model = climt.build_column_model(state, **kwargs)
    with pytest.raises(ValueError, match=match):
        model.run(*args)
