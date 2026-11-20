"""Each Fortran routine receives each argument in the place the 1.0 Cython
module put it.

The recorded reference cannot see every swap: many of its inputs hold one
value (the four SW albedos are all 0.06, the aerosol arrays and ``cloud_tau``
are zero, ``cloud_props_flag = 2`` leaves the cloud optical properties
unused). So these tests run the wrappers with the Fortran routines replaced
by recorders and compare what each routine is given, argument by argument,
with the calls of the deleted ``.pyx`` modules (``_rrtmg_lw.pyx`` at
2c7e8e74, ``_rrtmg_sw.pyx`` at 28b18038), which are written out below.

The wrappers call Fortran from numba kernels, through ``ExternalFunction``
symbols. A recorder is a cffi callback with the routine's prototype,
registered under a new symbol; the test compiles a copy of the kernel whose
globals name the recorders instead of the routines (not cached, so numba's
cache never holds it), and has the wrapper use it. Integer and double
scalars are passed by pointer: they are given distinct values here, and
compared by value; arrays are compared by address.
"""

import itertools
import types

import llvmlite.binding
import numba
import numpy as np
import pytest
from numba import types as numba_types

from climt_fortran import _numba_cffi
from test_parity import REFERENCE, _arguments, cffi_module, reference_inputs, replay_set_up

_SYMBOLS = itertools.count()


def _record_fortran_calls(monkeypatch, module, kernel_name, signature, routines):
    """Make ``module.<kernel_name>`` record its Fortran calls instead of making
    them. ``routines`` maps each global of the kernel that is an
    ``ExternalFunction`` to the name of its routine in ``module.lib``.
    Returns the list the calls are appended to, as ``(routine, arguments)``,
    each argument as ``(address, value at that address)``."""
    ffi, calls, keep = module.ffi, [], []
    namespace = dict(vars(module))
    for global_name, routine in routines.items():
        def record(*pointers, routine=routine):
            calls.append((routine, [(int(ffi.cast("uintptr_t", p)), p[0]) for p in pointers]))

        callback = ffi.callback(ffi.typeof(getattr(module.lib, routine)), record)
        keep.append(callback)
        symbol = "test_argument_order_%s_%d" % (routine, next(_SYMBOLS))
        llvmlite.binding.add_symbol(symbol, int(ffi.cast("uintptr_t", callback)))
        namespace[global_name] = numba_types.ExternalFunction(symbol, getattr(module, global_name).sig)
    kernel = getattr(module, kernel_name)
    copy = types.FunctionType(kernel.__code__, namespace, kernel.__name__ + "_recording")
    compiled = numba.njit(signature)(copy)  # not cached
    monkeypatch.setitem(_numba_cffi._KERNELS, (copy, id(signature)), (signature, compiled))
    monkeypatch.setattr(module, kernel_name, copy)
    monkeypatch.setattr(module, "_recording_callbacks", keep, raising=False)
    return calls


def _assert_arguments(routine, recorded, expected, arrays):
    """``recorded`` arguments of ``routine`` against ``expected``, in order: a
    name in ``arrays`` is that array's address; a number is the value pointed to."""
    expected_pairs = [("array", arrays[e].ctypes.data) if isinstance(e, str) else ("value", e) for e in expected]
    found = [("array", address) if isinstance(e, str) else ("value", value)
             for e, (address, value) in zip(expected, recorded)]
    names = {array.ctypes.data: name for name, array in arrays.items()}
    readable = [names.get(v, v) if kind == "array" else v for kind, v in found]
    assert len(recorded) == len(expected), (routine, readable)
    wrong = [(index, e, r) for index, (e, r, f, x) in enumerate(zip(expected, readable, found, expected_pairs))
             if f != x]
    assert not wrong, "%s: (argument, expected, found): %s" % (routine, wrong)


# Longwave: integers distinct from each other and from the sizes (4 columns,
# 30 layers). Only the wrapper's flags change: Fortran never runs.
LW_FLAGS = dict(cloud_overlap_method=3, calculate_change_up_flux=1, cloud_props_flag=5, ice_props_flag=6,
                liq_droplet_flag=7)
LW_MCICA_FLAGS = dict(LW_FLAGS, permute_seed=8, random_number_generator=9)
LW_IPLON = 2
LW_STATE = ["layer_pressure", "interface_pressure", "layer_temp", "interface_temp", "surface_temp", "h2o_vmr",
            "o3_vmr", "co2_vmr", "ch4_vmr", "n2o_vmr", "o2_vmr", "cfc11_vmr", "cfc12_vmr", "cfc22_vmr",
            "ccl4_vmr", "surf_emissivity"]
LW_FLUXES = ["aerosol_tau", "upward_longwave_flux", "downward_longwave_flux", "longwave_heating_rate",
             "up_lw_flux_clearsky", "down_lw_flux_clearsky", "lw_heating_rate_clearsky", "duflx_dt", "duflxc_dt"]


def _lw_radiation_call(flags, clouds):
    """``rrtmg_lw_[no]mcica_wrapper`` as 1.0 called it."""
    return ([4, 30, flags["cloud_overlap_method"], flags["calculate_change_up_flux"]] + LW_STATE
            + [flags["cloud_props_flag"], flags["ice_props_flag"], flags["liq_droplet_flag"]] + clouds + LW_FLUXES)


LW_EXPECTED = {
    "rrtmg_lw": [("rrtmg_lw_nomcica_wrapper", _lw_radiation_call(LW_FLAGS, [
        "cloud_fraction", "cloud_tau", "cloud_ice_path", "cloud_liq_path", "cloud_ice_eff_size",
        "cloud_droplet_eff_radius"]))],
    "rrtmg_lw_mcica": [
        ("mcica_subcol_lw_wrapper", [
            LW_IPLON, 4, 30, LW_MCICA_FLAGS["cloud_overlap_method"], LW_MCICA_FLAGS["permute_seed"],
            LW_MCICA_FLAGS["random_number_generator"], "layer_pressure", "cloud_fraction", "cloud_ice_path",
            "cloud_liq_path", "cloud_ice_eff_size", "cloud_droplet_eff_radius", "cloud_tau",
            "cloud_fraction_mcica", "cloud_ice_path_mcica", "cloud_liq_path_mcica", "cloud_ice_eff_size_mcica",
            "cloud_droplet_eff_radius_mcica", "cloud_tau_mcica"]),
        ("rrtmg_lw_mcica_wrapper", _lw_radiation_call(LW_MCICA_FLAGS, [
            "cloud_fraction_mcica", "cloud_tau_mcica", "cloud_ice_path_mcica", "cloud_liq_path_mcica",
            "cloud_ice_eff_size_mcica", "cloud_droplet_eff_radius_mcica"])),
    ],
}
LW_KERNELS = {
    "rrtmg_lw": ("_nomcica_kernel", "_NOMCICA_SIGNATURE", {"_rrtmg_lw_nomcica": "rrtmg_lw_nomcica_wrapper"}),
    "rrtmg_lw_mcica": ("_mcica_kernel", "_MCICA_SIGNATURE", {"_mcica_subcol": "mcica_subcol_lw_wrapper",
                                                             "_rrtmg_lw_mcica": "rrtmg_lw_mcica_wrapper"}),
}


@pytest.mark.parametrize("call", ["rrtmg_lw", "rrtmg_lw_mcica"])
def test_rrtmg_lw_argument_order(call, monkeypatch):
    module = cffi_module(call)
    kernel, signature, routines = LW_KERNELS[call]
    init = REFERENCE[call + "/init_order"][-1]
    flags = LW_MCICA_FLAGS if call == "rrtmg_lw_mcica" else LW_FLAGS
    inputs = reference_inputs(call)
    shape = inputs["upward_longwave_flux"].shape
    inputs.update(duflx_dt=np.zeros(shape), duflxc_dt=np.zeros(shape))
    if call == "rrtmg_lw_mcica":
        inputs["rrtm_iplon"] = LW_IPLON
    replay_set_up(call, module)
    try:
        getattr(module, init)(**{**_arguments("%s/init/%s/" % (call, init)), **flags})
        calls = _record_fortran_calls(monkeypatch, module, kernel, getattr(module, signature), routines)
        getattr(module, REFERENCE[call + "/function"].item())(*inputs.values())
    finally:
        replay_set_up(call, module)
    arrays = {name: value for name, value in inputs.items() if isinstance(value, np.ndarray)}
    assert [routine for routine, _ in calls] == [routine for routine, _ in LW_EXPECTED[call]]
    for (routine, recorded), (_, expected) in zip(calls, LW_EXPECTED[call]):
        _assert_arguments(routine, recorded, expected, arrays)
