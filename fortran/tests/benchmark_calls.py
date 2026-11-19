"""Per-call timings of the climt-fortran calls at one column.

Each call's arguments are recorded from the climt component that makes it
(see ``generate_cython_reference.record_call``), with one column; the
component's set-up calls are replayed once, then only the call itself is timed.

Run: ``NUMBA_NUM_THREADS=1 python fortran/tests/benchmark_calls.py [--save path.json]``.
"""

import argparse
import json
import os
import sys
import timeit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

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
        args = [a.item() if a.ndim == 0 else a for a in inputs.values()]
        best = min(timeit.repeat(lambda: call(*args), number=n, repeat=repeats))
        timings[name] = best / n * 1e6
    return timings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", help="write the timings as JSON to this path")
    parser.add_argument("-n", type=int, default=2000, help="calls per repeat")
    parser.add_argument("--repeats", type=int, default=5)
    options = parser.parse_args()
    timings = time_calls(options.n, options.repeats)
    for name, microseconds in timings.items():
        print("%-16s %10.1f us" % (name, microseconds))
    if options.save:
        with open(options.save, "w") as output:
            json.dump(timings, output, indent=2)
            output.write("\n")


if __name__ == "__main__":
    main()
