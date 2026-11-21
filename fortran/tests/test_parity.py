"""The cffi modules reproduce the Cython reference outputs, exactly and as fast.

Each case replays a call recorded from climt by ``generate_cython_reference.py``:
the component's set-up calls, then the call itself with the same arguments.
"""

import importlib
import importlib.util
import inspect
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from benchmark_calls import BARE_CALLS, MODULES, bare_call, time_calls_and_bare_min_of_runs  # noqa: E402
from generate_cython_reference import record_call  # noqa: E402

REQUIRED = os.environ.get("CLIMT_REQUIRE_FORTRAN") == "1"
FIXTURES = os.path.join(HERE, "fixtures")
REFERENCE = np.load(os.path.join(FIXTURES, "cython_reference.npz"))
with open(os.path.join(FIXTURES, "cython_timings.json")) as timings_file:
    CYTHON_TIMINGS = json.load(timings_file)
SPEED_TOLERANCE = 1.10


def _value(array):
    """A recorded argument: scalars (0-d arrays) as Python scalars, arrays as fresh copies."""
    return array.item() if array.ndim == 0 else array.copy()


def _arguments(prefix):
    """The recorded arguments under ``prefix``, in call order."""
    return {key[len(prefix):]: _value(REFERENCE[key]) for key in REFERENCE.files if key.startswith(prefix)}


def reference_inputs(call):
    """Arguments of the recorded call, by name, in call order."""
    return _arguments(call + "/in/")


def reference_outputs(call):
    """Array arguments of the recorded call after it, by name."""
    return _arguments(call + "/out/")


def cffi_module(call):
    """Import the module behind ``call``; skip (or fail, if required) unless its _cffi module is built."""
    name = "climt_fortran." + MODULES[call]
    if importlib.util.find_spec(name + "_cffi") is None:
        if REQUIRED:
            pytest.fail(name + "_cffi does not exist")
        pytest.skip(name + "_cffi is not built")
    return importlib.import_module(name)


def replay_set_up(call, module):
    """Make the set-up calls the component made before ``call``."""
    for function in REFERENCE[call + "/init_order"]:
        getattr(module, function)(*_arguments("%s/init/%s/" % (call, function)).values())


def run_call(call, inputs):
    """Replay the set-up calls, then call the module positionally, as climt does.

    Returns the array arguments after the call, by name.
    """
    module = cffi_module(call)
    replay_set_up(call, module)
    function = getattr(module, REFERENCE[call + "/function"].item())
    assert list(inspect.signature(function).parameters)[:len(inputs)] == list(inputs)
    function(*inputs.values())
    return {name: value for name, value in inputs.items() if isinstance(value, np.ndarray)}


# Machine-specific: the Cython ratios in cython_timings.json were recorded on
# one developer machine, so the gates mean nothing on other machines (CI).
speed_gate = pytest.mark.skipif(
    os.environ.get("CLIMT_SPEED_GATES") != "1",
    reason="speed gate against Cython timings from one machine; set CLIMT_SPEED_GATES=1 to run it")


@pytest.mark.parametrize("call", ["emanuel", "rrtmg_lw", "rrtmg_lw_mcica", "rrtmg_sw", "rrtmg_sw_mcica"])
def test_matches_cython_reference(call):
    outputs = run_call(call, reference_inputs(call))
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


def test_reference_exercises_convection():
    # Guard against a vacuous parity check: some recorded columns convect
    # (and some do not, so state left over from one column would show).
    outputs = reference_outputs("emanuel")
    assert set(outputs["conv_state"]) == {0, 1}
    assert np.any(outputs["dtemp"] != 0) and np.any(outputs["precip"] > 0)


def test_emanuel_rejects_float32():
    inputs = reference_inputs("emanuel")
    inputs["temp"] = inputs["temp"].astype(np.float32)
    with pytest.raises(TypeError, match="float64"):
        run_call("emanuel", inputs)


def test_emanuel_transposed_input_matches():
    # Inputs are read correctly, and outputs filled in place, when the
    # caller's arrays are not C-contiguous.
    inputs = {name: np.asfortranarray(value) if isinstance(value, np.ndarray) else value
              for name, value in reference_inputs("emanuel").items()}
    assert not inputs["temp"].flags.c_contiguous and not inputs["dtemp"].flags.c_contiguous
    outputs = run_call("emanuel", inputs)
    for name, expected in reference_outputs("emanuel").items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


def test_emanuel_fills_strided_outputs():
    # Outputs (and inputs) that are strided views of larger arrays.
    def strided(value):
        view = np.zeros(value.shape + (2,), value.dtype)[..., 0]
        view[...] = value
        return view

    inputs = {name: strided(value) if isinstance(value, np.ndarray) else value
              for name, value in reference_inputs("emanuel").items()}
    assert not inputs["precip"].flags.c_contiguous and not inputs["dtemp"].flags.c_contiguous
    outputs = run_call("emanuel", inputs)
    for name, expected in reference_outputs("emanuel").items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


def test_emanuel_requires_initialisation(monkeypatch):
    module = cffi_module("emanuel")
    monkeypatch.setattr(module, "_initialised", False)
    with pytest.raises(ValueError, match="not initialised"):
        module.convect(*reference_inputs("emanuel").values())


def test_emanuel_callable_from_numba():
    numba = pytest.importorskip("numba")
    from numba.core.typing import cffi_utils

    module = cffi_module("emanuel")
    compiled = importlib.import_module("climt_fortran._emanuel_convection_cffi")
    cffi_utils.register_module(compiled)
    ffi, lib = compiled.ffi, compiled.lib

    @numba.njit
    def convect(nlevs, max_conv_lev, num_tracers, dt, temp, q, qs, u, v, pmid, pint,
                conv_state, precip, vel, tscale, qscale, cbmf, cape, dtemp, dq, du, dv):
        ints = np.array([nlevs, max_conv_lev, num_tracers], dtype=np.int32)
        timestep = np.array([dt])
        tracers = np.zeros(1)
        dtracers = np.zeros(1)
        for col in range(temp.shape[0]):
            # Fortran may write to its inputs: give it copies, as the wrapper does.
            column = [temp[col].copy(), q[col].copy(), qs[col].copy(), u[col].copy(),
                      v[col].copy(), pmid[col].copy(), pint[col].copy()]
            lib.emanuel_convection(
                ffi.from_buffer(column[0]), ffi.from_buffer(column[1]), ffi.from_buffer(column[2]),
                ffi.from_buffer(column[3]), ffi.from_buffer(column[4]), ffi.from_buffer(column[5]),
                ffi.from_buffer(column[6]),
                ffi.from_buffer(ints[0:1]), ffi.from_buffer(ints[1:2]), ffi.from_buffer(ints[2:3]),
                ffi.from_buffer(timestep), ffi.from_buffer(conv_state[col:col + 1]),
                ffi.from_buffer(dtemp[col]), ffi.from_buffer(dq[col]), ffi.from_buffer(du[col]),
                ffi.from_buffer(dv[col]), ffi.from_buffer(precip[col:col + 1]),
                ffi.from_buffer(vel[col:col + 1]), ffi.from_buffer(tscale[col:col + 1]),
                ffi.from_buffer(qscale[col:col + 1]), ffi.from_buffer(cbmf[col:col + 1]),
                ffi.from_buffer(cape[col:col + 1]), ffi.from_buffer(tracers), ffi.from_buffer(dtracers))

    inputs = reference_inputs("emanuel")
    replay_set_up("emanuel", module)
    assert inputs["num_tracers"] == 0
    convect(*[value for name, value in inputs.items() if name != "ncols"])
    for name, expected in reference_outputs("emanuel").items():
        np.testing.assert_array_equal(inputs[name], expected, err_msg=name)


def _assert_not_slower_than_cython(calls):
    """The speed gate: for each call, (wrapper time / bare Fortran time), both
    timed in this run, is at most SPEED_TOLERANCE times the Cython module's
    ratio recorded in cython_timings.json. A ratio to a bare call of the same
    Fortran routines does not drift with the machine's state as absolute
    times do (see benchmark_calls)."""
    for call in calls:
        cffi_module(call)
    wrapper, bare = time_calls_and_bare_min_of_runs(runs=5, calls=calls)
    failures = []
    for call in calls:
        ratio, cython_ratio = wrapper[call] / bare[call], CYTHON_TIMINGS[call + "/ratio_to_bare"]
        if ratio > SPEED_TOLERANCE * cython_ratio:
            failures.append(
                "%s: wrapper %.2f us / bare %.2f us = %.3f > %.2f x %.3f (Cython / bare when recorded, "
                "bare then %.2f us)" % (call, wrapper[call], bare[call], ratio, SPEED_TOLERANCE, cython_ratio,
                                        CYTHON_TIMINGS[call + "/bare"]))
    assert not failures, "; ".join(failures)


@speed_gate
def test_emanuel_not_slower_than_cython():
    _assert_not_slower_than_cython(["emanuel"])


LW_CALLS = ["rrtmg_lw", "rrtmg_lw_mcica"]
# The LW arguments RRTMG writes to: fluxes, heating rates, McICA sub-columns.
LW_OUTPUTS = ("upward_longwave_flux", "downward_longwave_flux", "longwave_heating_rate",
              "up_lw_flux_clearsky", "down_lw_flux_clearsky", "lw_heating_rate_clearsky",
              "cloud_fraction_mcica", "cloud_tau_mcica", "cloud_ice_path_mcica",
              "cloud_liq_path_mcica", "cloud_ice_eff_size_mcica", "cloud_droplet_eff_radius_mcica")


def _strided(value):
    """A strided (not C-contiguous) view holding ``value``."""
    view = np.zeros(value.shape + (2,), value.dtype)[..., 0]
    view[...] = value
    return view


def test_reference_exercises_longwave():
    # Guard against a vacuous parity check: fluxes are filled, and the McICA
    # case has clouds in its sub-columns.
    for call in LW_CALLS:
        outputs = reference_outputs(call)
        assert np.all(outputs["upward_longwave_flux"] > 0), call
        assert np.any(outputs["longwave_heating_rate"] != 0), call
    assert np.any(reference_outputs("rrtmg_lw_mcica")["cloud_fraction_mcica"] > 0)


@pytest.mark.parametrize("call", LW_CALLS)
def test_rrtmg_lw_fills_non_contiguous_output(call):
    inputs = reference_inputs(call)
    for name in LW_OUTPUTS:
        if name in inputs:
            inputs[name] = _strided(inputs[name])
    assert not inputs["upward_longwave_flux"].flags.c_contiguous
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


@pytest.mark.parametrize("call", LW_CALLS)
def test_rrtmg_lw_transposed_input_matches(call):
    # Every array argument in Fortran order: inputs are read correctly, and
    # outputs filled in place.
    inputs = {name: np.asfortranarray(value) if isinstance(value, np.ndarray) else value
              for name, value in reference_inputs(call).items()}
    assert not inputs["layer_pressure"].flags.c_contiguous
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


@pytest.mark.parametrize("call", LW_CALLS)
@pytest.mark.parametrize("name", ["layer_temp", "upward_longwave_flux", "aerosol_tau"])
def test_rrtmg_lw_rejects_float32(call, name):
    inputs = reference_inputs(call)
    inputs[name] = inputs[name].astype(np.float32)
    with pytest.raises(TypeError, match="%s: expected a float64" % name):
        run_call(call, inputs)


def _array_cases(calls):
    """``(call, name)`` for every array argument of each call."""
    return [(call, name) for call in calls
            for name, value in reference_inputs(call).items() if isinstance(value, np.ndarray)]


def _trailing_axis_cases(calls):
    """``(call, name)`` for every array argument of more than one axis of each call."""
    return [(call, name) for call in calls
            for name, value in reference_inputs(call).items()
            if isinstance(value, np.ndarray) and value.ndim > 1]


@pytest.mark.parametrize("call, name", _array_cases(LW_CALLS))
def test_rrtmg_lw_rejects_wrong_shape(call, name):
    # Every array, C-contiguous and float64 but one row short: the compiled
    # kernel must refuse it before Fortran runs.
    inputs = reference_inputs(call)
    inputs[name] = inputs[name][:-1]
    assert inputs[name].flags.c_contiguous
    with pytest.raises(ValueError, match="%s: expected shape" % name):
        run_call(call, inputs)


@pytest.mark.parametrize("call, name", _trailing_axis_cases(LW_CALLS))
def test_rrtmg_lw_rejects_wrong_trailing_axis(call, name):
    # One short on the last axis (e.g. 15 of 16 bands of cloud_tau, 139 of
    # 140 g-points of cloud_tau_mcica), C-contiguous: refused before Fortran
    # writes anything.
    inputs = reference_inputs(call)
    inputs[name] = np.ascontiguousarray(inputs[name][..., :-1])
    flux = "downward_longwave_flux" if name == "upward_longwave_flux" else "upward_longwave_flux"
    inputs[flux][...] = -1
    with pytest.raises(ValueError, match="%s: expected shape" % name):
        run_call(call, inputs)
    assert np.all(inputs[flux] == -1)


@pytest.mark.parametrize("call, name, value", [
    (call, name, value) for call in LW_CALLS
    for name, value in [("num_cols", 2**32 + 1), ("num_layers", -2**31 - 1), ("rrtm_iplon", 2**31)]
    if call == "rrtmg_lw_mcica" or name != "rrtm_iplon"])
def test_rrtmg_lw_rejects_integers_outside_int32(call, name, value):
    # 1.0 converted them to int32_t, raising OverflowError; they must not wrap.
    inputs = reference_inputs(call)
    inputs[name] = value
    inputs["upward_longwave_flux"][...] = -1
    with pytest.raises(OverflowError):
        run_call(call, inputs)
    assert np.all(inputs["upward_longwave_flux"] == -1)


@pytest.mark.parametrize("call, name, value", [
    (call, name, value) for call in LW_CALLS for name in ("rrtm_iplon", "num_cols", "num_layers")
    for value in (2**31 - 1, -2**31) if call == "rrtmg_lw_mcica" or name != "rrtm_iplon"])
def test_rrtmg_lw_accepts_int32_limits(call, name, value):
    # The int32_t limits themselves pass the conversion. For the array sizes,
    # calculate_change_up_flux is set, so that no (num_cols, 1) duflx_dt
    # placeholders are allocated (2**31 - 1 rows would be 17 GB): the sizes
    # then disagree with the arrays, a ValueError naming the first array.
    # Fortran does not use rrtm_iplon.
    module = cffi_module(call)
    inputs = reference_inputs(call)
    inputs[name] = value
    if name in ("num_cols", "num_layers"):
        init = REFERENCE[call + "/init_order"][-1]
        replay_set_up(call, module)
        try:
            getattr(module, init)(**{**_arguments("%s/init/%s/" % (call, init)), "calculate_change_up_flux": 1})
            with pytest.raises(ValueError, match="layer_pressure: expected shape"):
                getattr(module, REFERENCE[call + "/function"].item())(*inputs.values())
        finally:
            replay_set_up(call, module)
    else:
        outputs = run_call(call, inputs)
        for output, expected in reference_outputs(call).items():
            np.testing.assert_array_equal(outputs[output], expected, err_msg=output)


@pytest.mark.parametrize("call", LW_CALLS)
def test_rrtmg_lw_read_only_inputs_match(call):
    # Read-only C-contiguous inputs are copied, as strided ones are.
    inputs = reference_inputs(call)
    for name, value in inputs.items():
        if isinstance(value, np.ndarray) and name not in LW_OUTPUTS:
            value.flags.writeable = False
    assert not inputs["layer_pressure"].flags.writeable and inputs["layer_pressure"].flags.c_contiguous
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


@pytest.mark.parametrize("call", ["rrtmg_lw", "rrtmg_lw_mcica", "rrtmg_sw", "rrtmg_sw_mcica"])
def test_rrtmg_masked_array_arguments_match(call):
    # numba rejects np.ma.MaskedArray; the wrapper passes plain views instead.
    inputs = {name: np.ma.asarray(value) if isinstance(value, np.ndarray) else value
              for name, value in reference_inputs(call).items()}
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(np.asarray(outputs[name]), expected, err_msg=name)


@pytest.mark.parametrize("call", LW_CALLS)
def test_rrtmg_lw_change_up_flux(call):
    # With calculate_change_up_flux set, duflx_dt and duflxc_dt
    # ((num_layers + 1, num_cols)) are required, and filled in place.
    module = cffi_module(call)
    inputs = reference_inputs(call)
    init = REFERENCE[call + "/init_order"][-1]
    arguments = _arguments("%s/init/%s/" % (call, init))
    arguments["calculate_change_up_flux"] = 1
    function = getattr(module, REFERENCE[call + "/function"].item())
    shape = inputs["upward_longwave_flux"].shape
    duflx_dt, duflxc_dt = _strided(np.zeros(shape)), np.zeros(shape)
    replay_set_up(call, module)
    try:
        getattr(module, init)(**arguments)
        with pytest.raises(TypeError, match="duflx_dt: expected a float64 array, got NoneType"):
            function(*inputs.values())
        function(*inputs.values(), duflx_dt=duflx_dt, duflxc_dt=duflxc_dt)
    finally:
        replay_set_up(call, module)
    assert np.all(duflx_dt > 0) and np.all(duflxc_dt > 0)
    expected = reference_outputs(call)["upward_longwave_flux"]
    np.testing.assert_array_equal(inputs["upward_longwave_flux"], expected)


def test_rrtmg_lw_initialise_checks_every_flag_first():
    module = cffi_module("rrtmg_lw")
    replay_set_up("rrtmg_lw", module)
    before = module._flags[module._CLOUD_OVERLAP_METHOD]
    with pytest.raises(TypeError):
        module.initialise_rrtm_radiation(1004.64, before + 1, 0, 2, 1, 1.5)
    with pytest.raises(OverflowError):
        module.initialise_rrtm_radiation_mcica(1004.64, before + 1, permute_seed=2**31)
    assert module._flags[module._CLOUD_OVERLAP_METHOD] == before


def test_rrtmg_lw_keeps_flags_rewritten_by_fortran():
    # RRTMG replaces an out-of-range cloud overlap method by 2, through the
    # pointer it is given; as in 1.0, later calls see 2.
    module = cffi_module("rrtmg_lw")
    replay_set_up("rrtmg_lw", module)
    arguments = _arguments("rrtmg_lw/init/initialise_rrtm_radiation/")
    arguments["cloud_overlap_method"] = 7
    try:
        module.initialise_rrtm_radiation(**arguments)
        module.rrtm_calculate_longwave_fluxes(*reference_inputs("rrtmg_lw").values())
        assert module._flags[module._CLOUD_OVERLAP_METHOD] == 2
    finally:
        replay_set_up("rrtmg_lw", module)


@speed_gate
def test_rrtmg_lw_not_slower_than_cython():
    _assert_not_slower_than_cython(LW_CALLS)


SW_CALLS = ["rrtmg_sw", "rrtmg_sw_mcica"]
# The SW arguments RRTMG writes to: fluxes, heating rates, McICA sub-columns.
SW_OUTPUTS = ("upward_shortwave_flux", "downward_shortwave_flux", "shortwave_heating_rate",
              "up_sw_flux_clearsky", "down_sw_flux_clearsky", "sw_heating_rate_clearsky",
              "cloud_fraction_mcica", "cloud_tau_mcica", "cloud_single_scat_alb_mcica",
              "cloud_asym_mcica", "cloud_fwd_scat_frac_mcica", "cloud_ice_path_mcica",
              "cloud_liq_path_mcica", "cloud_ice_eff_size_mcica", "cloud_droplet_eff_radius_mcica")


def _sw_init(call):
    """``(module, name, arguments)`` of the initialise call recorded for ``call``."""
    init = REFERENCE[call + "/init_order"][-1]
    return cffi_module(call), init, _arguments("%s/init/%s/" % (call, init))


def test_reference_exercises_shortwave():
    # Guard against a vacuous parity check: fluxes are filled, and the McICA
    # case has clouds in its sub-columns.
    for call in SW_CALLS:
        outputs = reference_outputs(call)
        assert np.all(outputs["downward_shortwave_flux"] > 0), call
        assert np.all(outputs["upward_shortwave_flux"] > 0), call
        assert np.any(outputs["shortwave_heating_rate"] != 0), call
    assert np.any(reference_outputs("rrtmg_sw_mcica")["cloud_fraction_mcica"] > 0)
    assert np.any(reference_outputs("rrtmg_sw_mcica")["cloud_liq_path_mcica"] > 0)


@pytest.mark.parametrize("call", SW_CALLS)
def test_rrtmg_sw_fills_non_contiguous_output(call):
    inputs = reference_inputs(call)
    for name in SW_OUTPUTS:
        if name in inputs:
            inputs[name] = _strided(inputs[name])
    assert not inputs["upward_shortwave_flux"].flags.c_contiguous
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


def _every_other(value):
    """A view of every other element of a longer 1-D array, holding ``value``."""
    view = np.zeros(2 * value.size)[::2]
    view[...] = value
    return view


@pytest.mark.parametrize("call", SW_CALLS)
def test_rrtmg_sw_transposed_input_matches(call):
    # Every array argument in Fortran order, and the 1-D ones (surface
    # temperature, albedos, zenith angle) every other element of a longer
    # array: inputs are read correctly, and outputs filled in place.
    inputs = {name: (np.asfortranarray(value) if value.ndim > 1 else _every_other(value))
              if isinstance(value, np.ndarray) else value
              for name, value in reference_inputs(call).items()}
    assert not inputs["layer_pressure"].flags.c_contiguous
    assert not inputs["sfc_direct_sw_albedo"].flags.c_contiguous
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


@pytest.mark.parametrize("call", SW_CALLS)
@pytest.mark.parametrize("one_d", ["contiguous", "strided"])
def test_rrtmg_sw_read_only_inputs_match(call, one_d):
    # Read-only C-contiguous inputs are copied, as strided ones are; also
    # with the 1-D inputs read-only views of every other element.
    inputs = reference_inputs(call)
    for name, value in inputs.items():
        if isinstance(value, np.ndarray) and name not in SW_OUTPUTS:
            if one_d == "strided" and value.ndim == 1:
                inputs[name] = value = _every_other(value)
            value.flags.writeable = False
    assert not inputs["layer_pressure"].flags.writeable and inputs["layer_pressure"].flags.c_contiguous
    assert inputs["surface_temp"].flags.c_contiguous == (one_d == "contiguous")
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


@pytest.mark.parametrize("call", SW_CALLS)
@pytest.mark.parametrize("name", ["layer_temp", "upward_shortwave_flux", "aerosol_at_55u", "cloud_tau"])
def test_rrtmg_sw_rejects_float32(call, name):
    inputs = reference_inputs(call)
    inputs[name] = inputs[name].astype(np.float32)
    with pytest.raises(TypeError, match="%s: expected a float64" % name):
        run_call(call, inputs)


@pytest.mark.parametrize("call", SW_CALLS)
def test_rrtmg_sw_rejects_read_only_output(call):
    inputs = reference_inputs(call)
    inputs["shortwave_heating_rate"].flags.writeable = False
    with pytest.raises(ValueError, match="shortwave_heating_rate: array is read-only"):
        run_call(call, inputs)


@pytest.mark.parametrize("call, name", _array_cases(SW_CALLS))
def test_rrtmg_sw_rejects_wrong_shape(call, name):
    # Every array, C-contiguous and float64 but one row short: refused
    # before Fortran runs.
    inputs = reference_inputs(call)
    inputs[name] = inputs[name][:-1]
    assert inputs[name].flags.c_contiguous
    with pytest.raises(ValueError, match="%s: expected shape" % name):
        run_call(call, inputs)


@pytest.mark.parametrize("call, name", _trailing_axis_cases(SW_CALLS))
def test_rrtmg_sw_rejects_wrong_trailing_axis(call, name):
    # One short on the last axis (e.g. 13 of 14 bands of cloud_tau, 111 of
    # 112 g-points of cloud_tau_mcica), C-contiguous: refused before Fortran
    # writes anything.
    inputs = reference_inputs(call)
    inputs[name] = np.ascontiguousarray(inputs[name][..., :-1])
    flux = "downward_shortwave_flux" if name == "upward_shortwave_flux" else "upward_shortwave_flux"
    inputs[flux][...] = -1
    with pytest.raises(ValueError, match="%s: expected shape" % name):
        run_call(call, inputs)
    assert np.all(inputs[flux] == -1)


def _sw_integer_cases(values):
    return [(call, name, value) for call in SW_CALLS
            for name in ("rrtm_iplon", "num_cols", "num_layers", "day_of_year") for value in values
            if call == "rrtmg_sw_mcica" or name != "rrtm_iplon"]


@pytest.mark.parametrize("call, name, value", _sw_integer_cases([2**31, -2**31 - 1, 2**32 + 1]))
def test_rrtmg_sw_rejects_integers_outside_int32(call, name, value):
    # 1.0 converted them to int32_t, raising OverflowError; they must not wrap.
    inputs = reference_inputs(call)
    inputs[name] = value
    inputs["upward_shortwave_flux"][...] = -1
    with pytest.raises(OverflowError):
        run_call(call, inputs)
    assert np.all(inputs["upward_shortwave_flux"] == -1)


@pytest.mark.parametrize("call, name, value", _sw_integer_cases([2**31 - 1, -2**31]))
def test_rrtmg_sw_accepts_int32_limits(call, name, value):
    # The int32_t limits themselves pass the conversion. The array sizes then
    # disagree with the arrays (a ValueError naming the first array); the
    # others reach Fortran, which does not use rrtm_iplon and takes the
    # day of year as given.
    inputs = reference_inputs(call)
    inputs[name] = value
    if name in ("num_cols", "num_layers"):
        with pytest.raises(ValueError, match="layer_pressure: expected shape"):
            run_call(call, inputs)
    else:
        outputs = run_call(call, inputs)
        assert np.all(np.isfinite(outputs["upward_shortwave_flux"]))
        assert np.all(outputs["downward_shortwave_flux"] > 0)


@pytest.mark.parametrize("call", SW_CALLS)
def test_rrtmg_sw_rejects_non_numeric_scalars(call):
    inputs = reference_inputs(call)
    inputs["solar_cycle_fraction"] = "0.5"
    with pytest.raises(TypeError):
        run_call(call, inputs)


@pytest.mark.parametrize("call", SW_CALLS)
def test_rrtmg_sw_accepts_integer_scalars(call):
    # flux_adj_sun_earth = 1 (an int) is the recorded 1.0.
    inputs = reference_inputs(call)
    assert inputs["flux_adj_sun_earth"] == 1.0
    inputs["flux_adj_sun_earth"] = 1
    outputs = run_call(call, inputs)
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)


def test_rrtmg_sw_initialise_checks_every_argument_first():
    # A bad argument raises before RRTMG or this module's state changes.
    module, init, arguments = _sw_init("rrtmg_sw_mcica")
    replay_set_up("rrtmg_sw_mcica", module)
    flags, solar = module._flags.copy(), module._solar_arrays
    bad = [("permute_seed", 2**31, OverflowError), ("cloud_overlap_method", 1.5, TypeError),
           ("solar_var_flag", -2**31 - 1, OverflowError),
           ("facular_sunspot_ampl", np.ones(2, np.float32), TypeError),
           ("facular_sunspot_ampl", np.ones(1), ValueError),
           ("solar_var_banded", np.ones(13), ValueError),
           ("solar_var_banded", np.ones((16, 1)), ValueError),
           ("solar_constant", "1367", TypeError)]
    read_only = np.ones(2)
    read_only.flags.writeable = False
    bad.append(("facular_sunspot_ampl", read_only, ValueError))
    plain = {key: value for key, value in arguments.items()
             if key not in ("permute_seed", "random_number_generator")}
    for name, value, error in bad:
        calls = [(module.initialise_rrtm_radiation_mcica, arguments)]
        if name in plain:
            calls.append((module.initialise_rrtm_radiation, plain))
        for initialise, valid in calls:
            with pytest.raises(error):
                initialise(**{**valid, name: value})
            np.testing.assert_array_equal(module._flags, flags)
            assert module._solar_arrays is solar
    module.initialise_rrtm_radiation_mcica(**{**arguments, "permute_seed": 2**31 - 1})
    assert module._flags[module._PERMUTE_SEED] == 2**31 - 1
    replay_set_up("rrtmg_sw_mcica", module)


@pytest.mark.parametrize("call", SW_CALLS)
def test_rrtmg_sw_requires_initialisation(call, monkeypatch):
    module = cffi_module(call)
    replay_set_up(call, module)
    monkeypatch.setattr(module, "_solar_arrays", None)
    function = getattr(module, REFERENCE[call + "/function"].item())
    with pytest.raises(ValueError, match="not initialised"):
        function(*reference_inputs(call).values())


@pytest.mark.parametrize("call", SW_CALLS)
@pytest.mark.parametrize("layout", ["contiguous", "strided"])
def test_rrtmg_sw_writes_sunspot_amplitudes_back(call, layout):
    # RRTMG rescales facular_sunspot_ampl (indsolvar) in place, once per
    # column, by the solar cycle fraction. 1.0 kept the caller's array and
    # passed it to every call, so the caller's array holds the rescaled
    # values and the next call starts from them.
    module, init, arguments = _sw_init(call)
    amplitudes = np.array([1.2, 0.8]) if layout == "contiguous" else _strided(np.array([1.2, 0.8]))
    banded = arguments["solar_var_banded"]
    inputs = reference_inputs(call)
    inputs["solar_cycle_fraction"] = fraction = 0.2
    function = getattr(module, REFERENCE[call + "/function"].item())
    replay_set_up(call, module)
    try:
        getattr(module, init)(**{**arguments, "facular_sunspot_ampl": amplitudes, "solar_var_flag": 1})
        weight = (fraction - 0.0229) / (0.3817 - 0.0229)
        expected = np.array([1.2, 0.8])
        for _ in range(2):
            function(*inputs.values())
            expected = 1 + weight ** inputs["num_cols"] * (expected - 1)
            np.testing.assert_allclose(amplitudes, expected, rtol=1e-12)
        np.testing.assert_array_equal(banded, 1.0)
    finally:
        replay_set_up(call, module)


@pytest.mark.parametrize("call, flag, value, rewritten", [
    ("rrtmg_sw", "cloud_overlap_method", 7, 2), ("rrtmg_sw", "aerosol_input_flag", 3, 0),
    ("rrtmg_sw_mcica", "random_number_generator", 5, 1)])
def test_rrtmg_sw_keeps_flags_rewritten_by_fortran(call, flag, value, rewritten):
    # RRTMG replaces an out-of-range overlap method by 2 and aerosol flag by
    # 0, and the McICA generator a nonzero generator choice by 1, through the
    # pointers it is given; as in 1.0, later calls see the new values.
    module, init, arguments = _sw_init(call)
    replay_set_up(call, module)
    try:
        getattr(module, init)(**{**arguments, flag: value})
        getattr(module, REFERENCE[call + "/function"].item())(*reference_inputs(call).values())
        assert module._flags[getattr(module, "_" + flag.upper())] == rewritten
    finally:
        replay_set_up(call, module)


@speed_gate
def test_rrtmg_sw_not_slower_than_cython():
    _assert_not_slower_than_cython(SW_CALLS)


@pytest.mark.parametrize("call", BARE_CALLS)
def test_bare_call_does_the_wrappers_work(call):
    # The speed gate's yardstick: the bare Fortran calls fill the same outputs,
    # bit for bit, as the wrapper given the same recorded column.
    module = cffi_module(call)
    init_calls, (function, inputs, _) = record_call(call, num_columns=1)
    run_bare, _, bare = bare_call(call)
    run_bare()
    for init_function, init_inputs, _ in init_calls:
        getattr(module, init_function)(*[_value(value) for value in init_inputs.values()])
    arguments = [_value(value) for value in inputs.values()]
    getattr(module, function)(*arguments)
    wrapper = {name: value for name, value in zip(inputs, arguments) if isinstance(value, np.ndarray)}
    assert set(wrapper) == set(bare)
    for name in wrapper:
        np.testing.assert_array_equal(bare[name], wrapper[name], err_msg=name)
    if call != "emanuel":
        flux = "upward_longwave_flux" if call.startswith("rrtmg_lw") else "upward_shortwave_flux"
        assert np.all(wrapper[flux] > 0)
