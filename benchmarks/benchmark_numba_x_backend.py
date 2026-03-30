"""
Benchmark: Numba x Backend matrix.

Measures wall-clock time for each pure-Python component across 4 configs:
  1. DataArrayBackend, no Numba
  2. DataArrayBackend, with Numba
  3. UnytBackend, no Numba
  4. UnytBackend, with Numba

Each config runs in a subprocess to ensure clean JIT state.

Usage:
    python benchmarks/benchmark_numba_x_backend.py
    python benchmarks/benchmark_numba_x_backend.py --ncol 4096 --iters 30
"""

import argparse
import json
import os
import subprocess
import sys

COMPONENTS = [
    "HeldSuarez",
    "GrayLongwaveRadiation",
    "Frierson06LongwaveOpticalDepth",
    "GridScaleCondensation",
    "DryConvectiveAdjustment",
    "BergerSolarInsolation",
    "SlabSurface",
    "Instellation",
    "BucketHydrology",
    "EmanuelConvectionPythonV3",
]

# Components that need a timestep (Steppers / ImplicitTendency)
NEEDS_TIMESTEP = {
    "GridScaleCondensation",
    "DryConvectiveAdjustment",
    "BucketHydrology",
    "EmanuelConvectionPythonV3",
}


def run_single_config(component_name, backend_name, disable_jit, ncol, nlev, iters):
    """Run a single component+backend+jit config in a subprocess, return elapsed seconds."""
    env = os.environ.copy()
    if disable_jit:
        env["NUMBA_DISABLE_JIT"] = "1"
    else:
        env.pop("NUMBA_DISABLE_JIT", None)

    script = f"""
import time, json
import numpy as np
import sympl
from datetime import timedelta, datetime
from sympl._core.backend import DataArrayBackend

import climt
from climt import (
    {component_name},
    UnytBackend,
    get_default_state,
    get_grid,
)

NCOL = {ncol}
NLEV = {nlev}
ITERS = {iters}
TIMESTEP = timedelta(minutes=10)

backend_name = "{backend_name}"
if backend_name == "unyt":
    sympl.set_backend(UnytBackend())
else:
    sympl.set_backend(DataArrayBackend())

grid = get_grid(nx=NCOL, ny=1, nz=NLEV)
comp = {component_name}()
state = get_default_state([comp], grid_state=grid)
state["time"] = datetime(2000, 6, 21, 12)

needs_ts = "{component_name}" in {repr(NEEDS_TIMESTEP)}
is_diagnostic = isinstance(comp, sympl.DiagnosticComponent)
is_stepper = isinstance(comp, sympl.Stepper)

# Warmup
try:
    if is_diagnostic:
        comp(state)
    elif is_stepper or needs_ts:
        comp(state, TIMESTEP)
    else:
        comp(state)
except Exception as e:
    print(json.dumps({{"error": str(e)}}))
    raise SystemExit(1)

# Timed iterations
t0 = time.perf_counter()
for _ in range(ITERS):
    if is_diagnostic:
        comp(state)
    elif is_stepper or needs_ts:
        comp(state, TIMESTEP)
    else:
        comp(state)
elapsed = time.perf_counter() - t0
print(json.dumps({{"elapsed": elapsed}}))
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env=env,
        timeout=300,
    )
    if result.returncode != 0:
        stderr_short = result.stderr.strip().split("\n")[-3:]
        return None, "\n".join(stderr_short)

    try:
        data = json.loads(result.stdout.strip())
        if "error" in data:
            return None, data["error"]
        return data["elapsed"], None
    except (json.JSONDecodeError, KeyError):
        return None, result.stdout.strip()[:200]


def main():
    parser = argparse.ArgumentParser(description="Numba x Backend benchmark matrix")
    parser.add_argument("--ncol", type=int, default=4096, help="Number of columns")
    parser.add_argument("--nlev", type=int, default=30, help="Number of levels")
    parser.add_argument("--iters", type=int, default=30, help="Iterations per config")
    args = parser.parse_args()

    configs = [
        ("DA+noJIT", "dataarray", True),
        ("DA+Numba", "dataarray", False),
        ("Unyt+noJIT", "unyt", True),
        ("Unyt+Numba", "unyt", False),
    ]

    header = f"{'Component':<35}"
    for label, _, _ in configs:
        header += f" | {label:>12}"
    header += f" | {'Numba spd':>10} | {'Unyt spd':>10}"
    print(f"Benchmark: {args.ncol} cols, {args.nlev} levels, {args.iters} iters")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    for comp_name in COMPONENTS:
        row = f"{comp_name:<35}"
        times = {}
        for label, backend, disable_jit in configs:
            elapsed, err = run_single_config(
                comp_name, backend, disable_jit, args.ncol, args.nlev, args.iters
            )
            if elapsed is not None:
                row += f" | {elapsed:>10.3f}s"
                times[label] = elapsed
            else:
                row += f" | {'FAIL':>11}s"
                print(f"  [{comp_name} {label}] Error: {err}", file=sys.stderr)

        # Compute speedups
        numba_speedup = ""
        if "DA+noJIT" in times and "DA+Numba" in times:
            numba_speedup = f"{times['DA+noJIT'] / times['DA+Numba']:.1f}x"

        unyt_speedup = ""
        if "DA+Numba" in times and "Unyt+Numba" in times:
            unyt_speedup = f"{times['DA+Numba'] / times['Unyt+Numba']:.1f}x"

        row += f" | {numba_speedup:>10} | {unyt_speedup:>10}"
        print(row)

    print()
    print("Numba spd = DA+noJIT / DA+Numba (kernel speedup from Numba)")
    print(
        "Unyt spd  = DA+Numba / Unyt+Numba (backend speedup from Unyt over DataArray)"
    )


if __name__ == "__main__":
    main()
