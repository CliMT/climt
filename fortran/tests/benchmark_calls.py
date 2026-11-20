"""Per-call timings of the climt-fortran calls at one column.

Each call's arguments are recorded from the climt component that makes it
(see ``generate_cython_reference.record_call``), with one column; the
component's set-up calls are replayed once, then one untimed warm-up call is
made, then only the call itself is timed.

Every array argument is restored from a pristine copy before each call,
outside the timed region, so every timed call sees the same inputs (Emanuel's
``cloud_base_mass_flux`` is read and updated in place by ``convect``). A
repeat's per-call time is the sum of its ``n`` single-call times over ``n``;
``time_calls`` returns the best repeat.

``fixtures/cython_timings.json`` holds the yardstick, from the Cython build
of the five ``.pyx`` modules at commit 73e8688c (same flags, same ``_lib``),
timed with those modules in place of the installed ones:

* ``<call>``: Cython per-call time, all five calls, the minimum over 5 runs
  of ``time_calls`` (``NUMBA_NUM_THREADS=1``, ``n=2000``, ``repeats=5``).
  Single runs varied by up to ~40% on the SW calls, hence the minimum.
* ``<call>/bare`` and ``<call>/ratio_to_bare``, for ``BARE_CALLS``: the bare
  Fortran time (``bare_call``: the same routines called through cffi with
  pointers built once) and Cython time / bare time, both minimums over the
  same 5 runs, each run timing the Cython call and then the bare one
  (``--runs 5 --save``). Recorded in a second session, in which the Cython
  times were within 1% of the ``<call>`` values.

Absolute times drift with the machine's state: an earlier baseline was 1.6
times this one a day later. So the gate is a ratio: a module passes if, for
each call, its time over the bare time, both from the same runs
(``time_calls_and_bare_min_of_runs``, 5 runs, ``NUMBA_NUM_THREADS=1``), is at
most 1.10 times ``<call>/ratio_to_bare``. The bare call does the same Fortran
work, so the ratio measures the wrapper's overhead and moves little with the
machine's speed. If the Fortran sources or compiler flags change, re-record
the file from that build: import its modules into ``sys.modules`` as
``climt_fortran._<name>``, then run ``main()`` with ``--save``.

Run: ``NUMBA_NUM_THREADS=1 python fortran/tests/benchmark_calls.py [--runs N] [--save path.json]``.
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

from generate_cython_reference import CALLS, record_call  # noqa: E402

MODULES = {
    "emanuel": "_emanuel_convection",
    "rrtmg_lw": "_rrtmg_lw",
    "rrtmg_lw_mcica": "_rrtmg_lw",
    "rrtmg_sw": "_rrtmg_sw",
    "rrtmg_sw_mcica": "_rrtmg_sw",
}
# Calls whose speed is gated, against a bare call to the same Fortran routines.
BARE_CALLS = ("emanuel", "rrtmg_lw", "rrtmg_lw_mcica")


def _best_time(call, resets, n, repeats):
    """Best per-call time of ``call()`` in microseconds: one untimed warm-up
    call, then ``repeats`` repeats of ``n`` calls, each preceded (untimed) by
    restoring every ``(array, original)`` of ``resets``."""
    call()  # warm-up, untimed
    best = float("inf")
    for _ in range(repeats):
        total = 0
        for _ in range(n):
            for array, original in resets:
                np.copyto(array, original)
            start = time.perf_counter_ns()
            call()
            total += time.perf_counter_ns() - start
        best = min(best, total)
    return best / n / 1e3


def time_calls(n=2000, repeats=5, calls=None):
    """Best per-call time, in microseconds, of each call (default ``CALLS``) at one column."""
    import importlib

    timings = {}
    for name in calls or CALLS:
        module = importlib.import_module("climt_fortran." + MODULES[name])
        init_calls, (function, inputs, _) = record_call(name, num_columns=1)
        for init_function, init_inputs, _ in init_calls:
            getattr(module, init_function)(*init_inputs.values())
        function = getattr(module, function)
        pristine = [a.item() if a.ndim == 0 else a for a in inputs.values()]
        args = [a.copy(order="K") if isinstance(a, np.ndarray) else a for a in pristine]
        resets = [(a, p) for a, p in zip(args, pristine) if isinstance(a, np.ndarray)]
        timings[name] = _best_time(lambda: function(*args), resets, n, repeats)
    return timings


def bare_call(name):
    """The Fortran routines behind call ``name`` (one of ``BARE_CALLS``), set
    up for one recorded column, as ``(call, resets, arrays)``: ``call()``
    makes the same sequence of Fortran calls the wrapper makes, through the
    cffi module, with pointers built once to contiguous copies of the recorded
    arguments; ``resets`` pairs each copy with its original, and ``arrays``
    maps each recorded array argument's name to its copy.

    The Fortran set-up is replayed through the cffi module too (a Cython
    module has its own copy of the Fortran state).
    """
    import importlib

    compiled = importlib.import_module("climt_fortran." + MODULES[name] + "_cffi")
    ffi, lib = compiled.ffi, compiled.lib
    init_calls, (_, inputs, _) = record_call(name, num_columns=1)
    init = {function: {arg: value.item() for arg, value in arguments.items()}
            for function, arguments, _ in init_calls}
    resets = []

    def array(copy, dtype=np.float64):
        """Pointer to the C-contiguous array ``copy`` (restored before each call)."""
        assert copy.dtype == dtype and copy.flags.c_contiguous
        resets.append((copy, copy.copy()))
        return ffi.from_buffer("int32_t[]" if dtype == np.int32 else "double[]", copy, True)

    def integer(value):
        return ffi.new("int32_t *", value)

    def double(value):
        return ffi.new("double *", value)

    copies = {arg: np.array(value, order="C") for arg, value in inputs.items() if value.ndim > 0}
    arrays = {arg: array(copy, copy.dtype) for arg, copy in copies.items()}
    if name == "emanuel":
        # convect: one Fortran call per column; one column here. No tracers.
        lib.init_emanuel_convection_fortran(
            integer(0), integer(init["init_emanuel_convection"]["py_least_conv_level"]),
            *[double(value) for arg, value in init["init_emanuel_convection"].items()
              if arg != "py_least_conv_level"])
        assert inputs["num_tracers"].item() == 0
        sequence = [(lib.emanuel_convection, [
            *[arrays[arg] for arg in ("temp", "q", "qs", "u", "v", "pmid", "pint")],
            integer(inputs["nlevs"].item()), integer(inputs["max_conv_lev"].item()), integer(0),
            double(inputs["dt"].item()), arrays["conv_state"],
            *[arrays[arg] for arg in ("dtemp", "dq", "du", "dv", "precip", "downdraft_vel_scale",
                                      "downdraft_temp_scale", "downdraft_q_scale",
                                      "cloud_base_mass_flux", "cape")],
            array(np.zeros(1)), array(np.zeros(1))])]
    else:
        lib.rrtmg_set_constants(*[double(value) for value in init["set_constants"].values()])
        flags = init["initialise_rrtm_radiation_mcica" if name.endswith("mcica") else "initialise_rrtm_radiation"]
        lib.rrtmg_lw_ini_wrapper(double(flags["cp_d_air"]))
        assert flags["calculate_change_up_flux"] == 0
        num_cols, num_layers = integer(inputs["num_cols"].item()), integer(inputs["num_layers"].item())
        overlap, change = integer(flags["cloud_overlap_method"]), integer(0)
        change_arrays = [array(np.zeros((inputs["num_cols"].item(), 1))) for _ in range(2)]
        common = [num_cols, num_layers, overlap, change,
                  *[arrays[arg] for arg in ("layer_pressure", "interface_pressure", "layer_temp", "interface_temp",
                                            "surface_temp", "h2o_vmr", "o3_vmr", "co2_vmr", "ch4_vmr", "n2o_vmr",
                                            "o2_vmr", "cfc11_vmr", "cfc12_vmr", "cfc22_vmr", "ccl4_vmr",
                                            "surf_emissivity")],
                  integer(flags["cloud_props_flag"]), integer(flags["ice_props_flag"]),
                  integer(flags["liq_droplet_flag"])]
        fluxes = [arrays[arg] for arg in ("aerosol_tau", "upward_longwave_flux", "downward_longwave_flux",
                                          "longwave_heating_rate", "up_lw_flux_clearsky", "down_lw_flux_clearsky",
                                          "lw_heating_rate_clearsky")] + change_arrays
        clouds = ("cloud_fraction", "cloud_tau", "cloud_ice_path", "cloud_liq_path", "cloud_ice_eff_size",
                  "cloud_droplet_eff_radius")
        if name == "rrtmg_lw":
            sequence = [(lib.rrtmg_lw_nomcica_wrapper, common + [arrays[arg] for arg in clouds] + fluxes)]
        else:
            sub_columns = [arrays[arg + "_mcica"] for arg in clouds]
            # The sub-column generator takes tau last.
            sequence = [
                (lib.mcica_subcol_lw_wrapper, [
                    integer(inputs["rrtm_iplon"].item()), num_cols, num_layers, overlap,
                    integer(flags["permute_seed"]), integer(flags["random_number_generator"]),
                    arrays["layer_pressure"],
                    *[arrays[arg] for arg in clouds[:1] + clouds[2:] + clouds[1:2]],
                    *(sub_columns[:1] + sub_columns[2:] + sub_columns[1:2])]),
                (lib.rrtmg_lw_mcica_wrapper, common + sub_columns + fluxes),
            ]

    def call():
        for routine, pointers in sequence:
            routine(*pointers)

    return call, resets, copies


def time_bare_calls(n=2000, repeats=5, calls=BARE_CALLS):
    """Best per-call time, in microseconds, of the bare Fortran calls behind
    each call, with the protocol of ``time_calls``."""
    return {name: _best_time(*bare_call(name)[:2], n, repeats) for name in calls}


def time_calls_min_of_runs(runs=5, n=2000, repeats=5, calls=None, on_run=None):
    """Per-call minimum, in microseconds, over ``runs`` separate runs of ``time_calls``.

    ``on_run(index, timings)`` is called after each run.
    """
    results = []
    for index in range(runs):
        results.append(time_calls(n, repeats, calls))
        if on_run is not None:
            on_run(index, results[-1])
    return {name: min(run[name] for run in results) for name in results[0]}


def time_calls_and_bare_min_of_runs(runs=5, n=2000, repeats=5, calls=BARE_CALLS, on_run=None):
    """``(wrapper, bare)``: per-call minima, in microseconds, over ``runs``
    runs, each timing ``calls`` through their modules and then bare. This is
    the gate statistic. ``on_run(index, timings)`` is called after each run,
    with the bare times under ``"<call>/bare"``."""
    wrapper, bare = [], []
    for index in range(runs):
        wrapper.append(time_calls(n, repeats, calls))
        bare.append(time_bare_calls(n, repeats, calls))
        if on_run is not None:
            on_run(index, {**wrapper[-1], **{name + "/bare": t for name, t in bare[-1].items()}})
    return ({name: min(run[name] for run in wrapper) for name in calls},
            {name: min(run[name] for run in bare) for name in calls})


def print_run(index, timings):
    print("run %d: %s" % (index + 1, ", ".join("%s %.1f" % item for item in timings.items())))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--save", help="write the timings and gate ratios as JSON to this path")
    parser.add_argument("-n", type=int, default=2000, help="calls per repeat")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--runs", type=int, default=5, help="separate runs; the minimum per call is kept")
    options = parser.parse_args()
    results = {}
    for index in range(options.runs):
        run = time_calls(options.n, options.repeats)
        run.update({name + "/bare": t for name, t in time_bare_calls(options.n, options.repeats).items()})
        print_run(index, run)
        for name, microseconds in run.items():
            results[name] = min(results.get(name, float("inf")), microseconds)
    for name in BARE_CALLS:
        results[name + "/ratio_to_bare"] = results[name] / results[name + "/bare"]
    for name, value in results.items():
        print("%-28s %10.3f" % (name, value))
    if options.save:
        with open(options.save, "w") as output:
            json.dump(results, output, indent=2)
            output.write("\n")


if __name__ == "__main__":
    main()
