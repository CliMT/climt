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

from benchmark_calls import MODULES, time_calls_min_of_runs  # noqa: E402

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


@pytest.mark.parametrize("call", ["emanuel", "rrtmg_lw", "rrtmg_lw_mcica"])
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


def test_emanuel_not_slower_than_cython():
    import numba

    cffi_module("emanuel")
    threads = numba.get_num_threads()
    numba.set_num_threads(1)  # the gate protocol of benchmark_calls
    try:
        microseconds = time_calls_min_of_runs(runs=5, calls=["emanuel"])["emanuel"]
    finally:
        numba.set_num_threads(threads)
    limit = SPEED_TOLERANCE * CYTHON_TIMINGS["emanuel"]
    assert microseconds <= limit, "emanuel: %.1f us > %.1f us (%.2f x Cython)" % (
        microseconds, limit, SPEED_TOLERANCE)


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


def _lw_array_cases():
    return [(call, name) for call in LW_CALLS
            for name, value in reference_inputs(call).items() if isinstance(value, np.ndarray)]


@pytest.mark.parametrize("call, name", _lw_array_cases())
def test_rrtmg_lw_rejects_wrong_shape(call, name):
    # Every array, C-contiguous and float64 but one row short: the compiled
    # kernel must refuse it before Fortran runs.
    inputs = reference_inputs(call)
    inputs[name] = inputs[name][:-1]
    assert inputs[name].flags.c_contiguous
    with pytest.raises(ValueError, match="%s: expected shape" % name):
        run_call(call, inputs)


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


def test_rrtmg_lw_not_slower_than_cython():
    import numba

    cffi_module("rrtmg_lw")
    threads = numba.get_num_threads()
    numba.set_num_threads(1)  # the gate protocol of benchmark_calls
    try:
        timings = time_calls_min_of_runs(runs=5, calls=LW_CALLS)
    finally:
        numba.set_num_threads(threads)
    for call in LW_CALLS:
        limit = SPEED_TOLERANCE * CYTHON_TIMINGS[call]
        assert timings[call] <= limit, "%s: %.1f us > %.1f us (%.2f x Cython)" % (
            call, timings[call], limit, SPEED_TOLERANCE)
