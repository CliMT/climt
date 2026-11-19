"""Per-call timings of the climt-fortran calls at one column.

Each call's arguments are recorded from the climt component that makes it
(see ``generate_cython_reference.record_call``), with one column; the
component's set-up calls are replayed once, then only the call itself is timed.

Every array argument is restored from a pristine copy before each call,
outside the timed region, so every timed call sees the same inputs (Emanuel's
``cloud_base_mass_flux`` is read and updated in place by ``convect``). A
repeat's per-call time is the sum of its ``n`` single-call times over ``n``;
``time_calls`` returns the best repeat.

``fixtures/cython_timings.json`` is the speed gate for the cffi modules. It was
recorded as the minimum, per call, over 5 separate runs of ``time_calls``
(``NUMBA_NUM_THREADS=1``, ``n=2000``, ``repeats=5``; ``--runs 5 --save``),
because single runs varied by up to ~40% run to run on the SW calls.

Run: ``NUMBA_NUM_THREADS=1 python fortran/tests/benchmark_calls.py [--save path.json]``.
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


def time_calls(n=2000, repeats=5):
    """Best per-call time, in microseconds, of each call in ``CALLS`` at one column."""
    import importlib

    timings = {}
    for name in CALLS:
        module = importlib.import_module("climt_fortran." + MODULES[name])
        init_calls, (function, inputs, _) = record_call(name, num_columns=1)
        for init_function, init_inputs, _ in init_calls:
            getattr(module, init_function)(*init_inputs.values())
        call = getattr(module, function)
        pristine = [a.item() if a.ndim == 0 else a for a in inputs.values()]
        args = [a.copy(order="K") if isinstance(a, np.ndarray) else a for a in pristine]
        resets = [(a, p) for a, p in zip(args, pristine) if isinstance(a, np.ndarray)]
        best = float("inf")
        for _ in range(repeats):
            total = 0
            for _ in range(n):
                for array, original in resets:
                    np.copyto(array, original)
                start = time.perf_counter_ns()
                call(*args)
                total += time.perf_counter_ns() - start
            best = min(best, total)
        timings[name] = best / n / 1e3
    return timings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", help="write the timings as JSON to this path")
    parser.add_argument("-n", type=int, default=2000, help="calls per repeat")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--runs", type=int, default=1, help="separate runs of time_calls; the minimum per call is kept")
    options = parser.parse_args()
    runs = [time_calls(options.n, options.repeats) for _ in range(options.runs)]
    for index, run in enumerate(runs):
        print("run %d: %s" % (index + 1, ", ".join("%s %.1f" % item for item in run.items())))
    timings = {name: min(run[name] for run in runs) for name in CALLS}
    for name, microseconds in timings.items():
        print("%-16s %10.1f us" % (name, microseconds))
    if options.save:
        with open(options.save, "w") as output:
            json.dump(timings, output, indent=2)
            output.write("\n")


if __name__ == "__main__":
    main()
