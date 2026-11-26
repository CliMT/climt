"""Benchmark build_column_model (approach C) against the sympl loop.

The GMD column of examples/gmd_radiative_convective_unyt.py, with
EmanuelConvection and with EmanuelConvectionPython. Times (a) the sympl loop
under UnytBackend (1 and 16 columns) and (b) ColumnModel.run (1, 16, 256 and
4096 columns), with 1 and NUMBA_NUM_THREADS threads: 20 warm-up steps, then
200 timed steps. Reports ms per step, ms per step per column and the share of
time spent outside the components' array_call.

Writes benchmarks/results/scm_build.json.

Usage: python benchmarks/scm_build.py [--steps N] [--warmup N] [--columns 1,16,...]
"""
import argparse
import json
import os
import platform
import subprocess
import time
from datetime import timedelta

import numba
import numpy as np
import sympl

import climt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "scm_build.json")
SYMPL_COLUMNS = (1, 16)


class Timer(object):
    """Wraps an array_call, accumulating the time spent in it."""

    def __init__(self, call):
        self.call = call
        self.total = 0.

    def __call__(self, *args, **kwargs):
        start = time.perf_counter()
        try:
            return self.call(*args, **kwargs)
        finally:
            self.total += time.perf_counter() - start


def set_values(quantity, value):
    if hasattr(quantity, "values"):
        quantity.values[:] = value
    else:
        quantity.data[:] = value


def make_column(convection, ncol=None):
    """(state, components, stepper) of the GMD column; ``ncol`` columns along lon."""
    conv = getattr(climt, convection)()
    tendency_components = [conv, climt.RRTMGLongwave(), climt.RRTMGShortwave(),
                           climt.SlabSurface()]
    stepper = climt.SimplePhysics()
    components = [stepper] + tendency_components
    if ncol is None:
        state = climt.get_default_state(components)
    else:
        state = climt.get_default_state(components, grid_state=climt.get_grid(nx=ncol, ny=1))
    values = {
        "air_temperature": 270., "surface_albedo_for_direct_shortwave": 0.5,
        "surface_albedo_for_direct_near_infrared": 0.5,
        "surface_albedo_for_diffuse_shortwave": 0.5, "zenith_angle": np.pi / 2.5,
        "surface_temperature": 300., "ocean_mixed_layer_thickness": 5., "area_type": "sea",
        "eastward_wind": 3.,
    }
    for name, value in values.items():
        set_values(state[name], value)
    return state, tendency_components, stepper


def bench_sympl(convection, ncol, warmup, steps):
    sympl.set_backend(climt.UnytBackend())
    try:
        state, tendency_components, stepper = make_column(convection, ncol)
        timers = []
        for c in tendency_components + [stepper]:
            c.array_call = Timer(c.array_call)
            timers.append(c.array_call)
        time_stepper = sympl.AdamsBashforth(*tendency_components, order=3)
        timestep = climt.UnytTimeDelta(timedelta(minutes=5))

        def step(state):
            diagnostics, state = time_stepper(state, timestep)
            state.update(diagnostics)
            diagnostics, new_state = stepper(state, timestep)
            state.update(diagnostics)
            state.update(new_state)
            state["time"] += timestep
            set_values(state["eastward_wind"], 3.)
            return state

        for _ in range(warmup):
            state = step(state)
        for t in timers:
            t.total = 0.
        start = time.perf_counter()
        for _ in range(steps):
            state = step(state)
        elapsed = time.perf_counter() - start
    finally:
        sympl.set_backend(sympl.DataArrayBackend())
    return elapsed, sum(t.total for t in timers)


def bench_compiled(convection, ncol, threads, warmup, steps):
    sympl.set_backend(sympl.DataArrayBackend())
    state, tendency_components, stepper = make_column(convection)
    kwargs = dict(timestep=timedelta(minutes=5), tendency_components=tendency_components,
                  steppers=[stepper], prescribed=["eastward_wind"], threads=threads)
    if ncol > 1:
        kwargs["scenarios"] = climt.Scenarios(ncol)
    model = climt.build_column_model(state, **kwargs)
    timers = []
    for c in model._tendency_components + model._steppers:
        c._array_call = Timer(c._array_call)
        timers.append(c._array_call)
    model.run(warmup, outputs=[])
    for t in timers:
        t.total = 0.
    start = time.perf_counter()
    model.run(steps, outputs=[])
    elapsed = time.perf_counter() - start
    return elapsed, sum(t.total for t in timers)


def row(convection, ncol, threads, path, elapsed, inside, steps):
    ms = 1e3 * elapsed / steps
    r = {"convection": convection, "columns": ncol, "threads": threads, "path": path,
         "ms_per_step": ms, "ms_per_step_per_column": ms / ncol,
         "outside_kernels": 1. - inside / elapsed}
    print("{convection:24s} {columns:5d} cols {threads:2d} thr {path:9s} "
          "{ms_per_step:10.3f} ms/step {ms_per_step_per_column:9.4f} ms/step/col "
          "outside {outside_kernels:.0%}".format(**r), flush=True)
    return r


def commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                                       text=True).strip()
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--warmup", type=int, default=20)
    parser.add_argument("--columns", default="1,16,256,4096")
    parser.add_argument("--out", default=OUT)
    args = parser.parse_args()
    columns = [int(c) for c in args.columns.split(",")]
    thread_counts = sorted({1, numba.config.NUMBA_NUM_THREADS})

    results = []
    for convection in ("EmanuelConvection", "EmanuelConvectionPython"):
        for ncol in columns:
            if ncol in SYMPL_COLUMNS:
                elapsed, inside = bench_sympl(convection, ncol, args.warmup, args.steps)
                results.append(row(convection, ncol, 1, "sympl", elapsed, inside, args.steps))
            for threads in thread_counts:
                elapsed, inside = bench_compiled(convection, ncol, threads, args.warmup,
                                                 args.steps)
                results.append(row(convection, ncol, threads, "compiled", elapsed, inside,
                                   args.steps))

    out = {"machine": "{} {} ({} cores)".format(platform.system(), platform.machine(),
                                                os.cpu_count()),
           "python": platform.python_version(), "numba": numba.__version__,
           "commit": commit(), "approach": "C", "steps": args.steps,
           "warmup": args.warmup, "results": results}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
