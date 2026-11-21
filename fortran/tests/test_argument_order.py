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
compared by value; the flags and the SW solar constant also by address (in
the module's ``_flags`` and ``_solar_constant``, where Fortran's rewrites
must land); arrays by address.
"""

import itertools
import types
from typing import NamedTuple

import llvmlite.binding
import numba
import numpy as np
import pytest
from numba import types as numba_types

from climt_fortran import _numba_cffi
from test_parity import REFERENCE, _arguments, cffi_module, reference_inputs, replay_set_up

_SYMBOLS = itertools.count()
# Every recording callback ever registered with llvmlite.binding.add_symbol:
# the symbols cannot be removed, so the callbacks must never be freed.
_CALLBACKS = []


class Held(NamedTuple):
    """An argument pointing to ``module.<array>[module.<index>]`` (index 0 if
    None), which holds ``value``: the module's flags and solar constant,
    passed by pointer so that Fortran's rewrites stick."""
    array: str
    index: object
    value: object


def _flag(flags, name):
    """The expected argument for flag ``name``, set to ``flags[name]``."""
    return Held("_flags", "_" + name.upper(), flags[name])


def _record_fortran_calls(monkeypatch, module, kernel_name, signature, routines):
    """Make ``module.<kernel_name>`` record its Fortran calls instead of making
    them. ``routines`` maps each global of the kernel that is an
    ``ExternalFunction`` to the name of its routine in ``module.lib``.
    Returns the list the calls are appended to, as ``(routine, arguments)``,
    each argument as ``(address, value at that address)``."""
    ffi, calls = module.ffi, []
    namespace = dict(vars(module))
    for global_name, routine in routines.items():
        def record(*pointers, routine=routine):
            calls.append((routine, [(int(ffi.cast("uintptr_t", p)), p[0]) for p in pointers]))

        callback = ffi.callback(ffi.typeof(getattr(module.lib, routine)), record)
        _CALLBACKS.append(callback)
        symbol = "test_argument_order_%s_%d" % (routine, next(_SYMBOLS))
        llvmlite.binding.add_symbol(symbol, int(ffi.cast("uintptr_t", callback)))
        namespace[global_name] = numba_types.ExternalFunction(symbol, getattr(module, global_name).sig)
    kernel = getattr(module, kernel_name)
    copy = types.FunctionType(kernel.__code__, namespace, kernel.__name__ + "_recording")
    compiled = numba.njit(signature)(copy)  # not cached
    monkeypatch.setitem(_numba_cffi._KERNELS, (copy, id(signature)), (signature, compiled))
    monkeypatch.setattr(module, kernel_name, copy)
    return calls


def _expected_pair(module, arrays, expected):
    """What an argument must be: ``("array", address)`` for a name in
    ``arrays``, ``("held", (address, value))`` for a ``Held``, else
    ``("value", expected)``, the value pointed to."""
    if isinstance(expected, str):
        return "array", arrays[expected].ctypes.data
    if isinstance(expected, Held):
        array = getattr(module, expected.array)
        index = 0 if expected.index is None else getattr(module, expected.index)
        return "held", (array.ctypes.data + array.itemsize * index, expected.value)
    return "value", expected


def _found_pair(expected, address, value):
    """The recorded argument, in the form of ``_expected_pair``."""
    if isinstance(expected, str):
        return "array", address
    if isinstance(expected, Held):
        return "held", (address, value)
    return "value", value


def _assert_arguments(module, routine, recorded, expected, arrays):
    """``recorded`` arguments of ``routine`` against ``expected``, in order: a
    name in ``arrays`` is that array's address, a ``Held`` an address in the
    module's arrays and its value, a number the value pointed to."""
    expected_pairs = [_expected_pair(module, arrays, e) for e in expected]
    found = [_found_pair(e, address, value) for e, (address, value) in zip(expected, recorded)]
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
    return ([4, 30, _flag(flags, "cloud_overlap_method"), _flag(flags, "calculate_change_up_flux")] + LW_STATE
            + [_flag(flags, "cloud_props_flag"), _flag(flags, "ice_props_flag"), _flag(flags, "liq_droplet_flag")]
            + clouds + LW_FLUXES)


LW_EXPECTED = {
    "rrtmg_lw": [("rrtmg_lw_nomcica_wrapper", _lw_radiation_call(LW_FLAGS, [
        "cloud_fraction", "cloud_tau", "cloud_ice_path", "cloud_liq_path", "cloud_ice_eff_size",
        "cloud_droplet_eff_radius"]))],
    "rrtmg_lw_mcica": [
        ("mcica_subcol_lw_wrapper", [
            LW_IPLON, 4, 30, _flag(LW_MCICA_FLAGS, "cloud_overlap_method"), _flag(LW_MCICA_FLAGS, "permute_seed"),
            _flag(LW_MCICA_FLAGS, "random_number_generator"), "layer_pressure", "cloud_fraction", "cloud_ice_path",
            "cloud_liq_path", "cloud_ice_eff_size", "cloud_droplet_eff_radius", "cloud_tau",
            "cloud_fraction_mcica", "cloud_ice_path_mcica", "cloud_liq_path_mcica", "cloud_ice_eff_size_mcica",
            "cloud_droplet_eff_radius_mcica", "cloud_tau_mcica"]),
        ("rrtmg_lw_mcica_wrapper", _lw_radiation_call(LW_MCICA_FLAGS, [
            "cloud_fraction_mcica", "cloud_tau_mcica", "cloud_ice_path_mcica", "cloud_liq_path_mcica",
            "cloud_ice_eff_size_mcica", "cloud_droplet_eff_radius_mcica"])),
    ],
}
# Shortwave: as for longwave, with distinct day of year and double scalars.
SW_FLAGS = dict(cloud_overlap_method=3, cloud_props_flag=5, ice_props_flag=6, liq_droplet_flag=7,
                aerosol_input_flag=8, solar_var_flag=9)
SW_MCICA_FLAGS = dict(SW_FLAGS, permute_seed=10, random_number_generator=11)
SW_IPLON, SW_DAY, SW_FLUX_ADJ, SW_CYCLE, SW_SOLAR_CONSTANT = 2, 17, 1.25, 0.375, 1361.5
SW_STATE = ["layer_pressure", "interface_pressure", "layer_temp", "interface_temp", "surface_temp", "h2o_vmr",
            "o3_vmr", "co2_vmr", "ch4_vmr", "n2o_vmr", "o2_vmr", "sfc_direct_sw_albedo", "sfc_diffuse_sw_albedo",
            "sfc_direct_nir_albedo", "sfc_diffuse_nir_albedo", "cos_zenith_angle"]
SW_FLUXES = ["aerosol_tau", "aerosol_single_scat_alb", "aerosol_asym", "aerosol_at_55u", "upward_shortwave_flux",
             "downward_shortwave_flux", "shortwave_heating_rate", "up_sw_flux_clearsky", "down_sw_flux_clearsky",
             "sw_heating_rate_clearsky", "solar_var_banded", "facular_sunspot_ampl", SW_CYCLE]
SW_CLOUDS = ["cloud_fraction", "cloud_tau", "cloud_single_scat_alb", "cloud_asym", "cloud_fwd_scat_frac",
             "cloud_ice_path", "cloud_liq_path", "cloud_ice_eff_size", "cloud_droplet_eff_radius"]


def _sw_radiation_call(flags, clouds):
    """``rrtmg_sw_[no]mcica_wrapper`` as 1.0 called it."""
    return ([4, 30, _flag(flags, "cloud_overlap_method"), _flag(flags, "aerosol_input_flag")] + SW_STATE
            + [SW_FLUX_ADJ, SW_DAY, Held("_solar_constant", None, SW_SOLAR_CONSTANT), _flag(flags, "solar_var_flag"),
               _flag(flags, "cloud_props_flag"), _flag(flags, "ice_props_flag"), _flag(flags, "liq_droplet_flag")]
            + clouds + SW_FLUXES)


EXPECTED = {
    **LW_EXPECTED,
    "rrtmg_sw": [("rrtmg_sw_nomcica_wrapper", _sw_radiation_call(SW_FLAGS, SW_CLOUDS))],
    "rrtmg_sw_mcica": [
        ("mcica_subcol_sw_wrapper", [
            SW_IPLON, 4, 30, _flag(SW_MCICA_FLAGS, "cloud_overlap_method"), _flag(SW_MCICA_FLAGS, "permute_seed"),
            _flag(SW_MCICA_FLAGS, "random_number_generator"), "layer_pressure", "cloud_fraction", "cloud_ice_path",
            "cloud_liq_path", "cloud_ice_eff_size", "cloud_droplet_eff_radius", "cloud_tau",
            "cloud_single_scat_alb", "cloud_asym", "cloud_fwd_scat_frac", "cloud_fraction_mcica",
            "cloud_ice_path_mcica", "cloud_liq_path_mcica", "cloud_ice_eff_size_mcica",
            "cloud_droplet_eff_radius_mcica", "cloud_tau_mcica", "cloud_single_scat_alb_mcica",
            "cloud_asym_mcica", "cloud_fwd_scat_frac_mcica"]),
        # The radiation call takes the non-McICA ice size and droplet radius.
        ("rrtmg_sw_mcica_wrapper", _sw_radiation_call(SW_MCICA_FLAGS, [
            "cloud_fraction_mcica", "cloud_tau_mcica", "cloud_single_scat_alb_mcica", "cloud_asym_mcica",
            "cloud_fwd_scat_frac_mcica", "cloud_ice_path_mcica", "cloud_liq_path_mcica", "cloud_ice_eff_size",
            "cloud_droplet_eff_radius"])),
    ],
}
# For each call: the kernel, its signature, and its ExternalFunction globals.
KERNELS = {
    "rrtmg_lw": ("_nomcica_kernel", "_NOMCICA_SIGNATURE", {"_rrtmg_lw_nomcica": "rrtmg_lw_nomcica_wrapper"}),
    "rrtmg_lw_mcica": ("_mcica_kernel", "_MCICA_SIGNATURE", {"_mcica_subcol": "mcica_subcol_lw_wrapper",
                                                             "_rrtmg_lw_mcica": "rrtmg_lw_mcica_wrapper"}),
    "rrtmg_sw": ("_nomcica_kernel", "_NOMCICA_SIGNATURE", {"_rrtmg_sw_nomcica": "rrtmg_sw_nomcica_wrapper"}),
    "rrtmg_sw_mcica": ("_mcica_kernel", "_MCICA_SIGNATURE", {"_mcica_subcol": "mcica_subcol_sw_wrapper",
                                                             "_rrtmg_sw_mcica": "rrtmg_sw_mcica_wrapper"}),
}


def _set_up(call, inputs):
    """Initialise with the distinct flags (and SW scalars) above; return the
    initialise call's arrays, by name. Changes ``inputs`` to match."""
    init = REFERENCE[call + "/init_order"][-1]
    arguments = _arguments("%s/init/%s/" % (call, init))
    mcica = call.endswith("mcica")
    if call.startswith("rrtmg_lw"):
        arguments.update(LW_MCICA_FLAGS if mcica else LW_FLAGS)
        shape = inputs["upward_longwave_flux"].shape
        inputs.update(duflx_dt=np.zeros(shape), duflxc_dt=np.zeros(shape))  # calculate_change_up_flux=1
        iplon = LW_IPLON
    else:
        arguments.update(SW_MCICA_FLAGS if mcica else SW_FLAGS, solar_constant=SW_SOLAR_CONSTANT)
        inputs.update(day_of_year=SW_DAY, flux_adj_sun_earth=SW_FLUX_ADJ, solar_cycle_fraction=SW_CYCLE)
        iplon = SW_IPLON
    if mcica:
        inputs["rrtm_iplon"] = iplon
    getattr(cffi_module(call), init)(**arguments)
    return {name: value for name, value in arguments.items() if isinstance(value, np.ndarray)}


@pytest.mark.parametrize("call", ["rrtmg_lw", "rrtmg_lw_mcica", "rrtmg_sw", "rrtmg_sw_mcica"])
def test_argument_order(call, monkeypatch):
    module = cffi_module(call)
    kernel, signature, routines = KERNELS[call]
    inputs = reference_inputs(call)
    replay_set_up(call, module)
    try:
        init_arrays = _set_up(call, inputs)
        calls = _record_fortran_calls(monkeypatch, module, kernel, getattr(module, signature), routines)
        getattr(module, REFERENCE[call + "/function"].item())(*inputs.values())
    finally:
        replay_set_up(call, module)
    arrays = {name: value for name, value in inputs.items() if isinstance(value, np.ndarray)}
    arrays.update(init_arrays)
    assert [routine for routine, _ in calls] == [routine for routine, _ in EXPECTED[call]]
    for (routine, recorded), (_, expected) in zip(calls, EXPECTED[call]):
        _assert_arguments(module, routine, recorded, expected, arrays)
