"""
Benchmark: Held-Suarez GCM with GFS dynamical core.

Runs 300 time steps of the Held-Suarez test case using UnytBackend,
with and without Numba JIT. Each config runs in a subprocess to
ensure clean JIT state.

Usage:
    python benchmarks/benchmark_held_suarez_gcm.py
    python benchmarks/benchmark_held_suarez_gcm.py --steps 100
"""
import argparse
import json
import os
import subprocess
import sys
import time


WORKER_SCRIPT = '''
import time, json, os
import numpy as np
import sympl
from datetime import timedelta

import climt
from climt import HeldSuarez, UnytBackend, get_default_state, get_grid
from gfs_dynamical_core import GFSDynamicalCore

STEPS = {steps}

sympl.set_backend(UnytBackend())

model_time_step = timedelta(seconds=600)

held_suarez = HeldSuarez()
dycore = GFSDynamicalCore([held_suarez])
grid = get_grid(nx=128, ny=62)
my_state = get_default_state([dycore], grid_state=grid)

my_state['eastward_wind'].values[:] = np.random.randn(
    *my_state['eastward_wind'].shape
)

# Warmup (1 step, includes JIT compilation)
t_warmup_start = time.perf_counter()
diag, output = dycore(my_state, model_time_step)
my_state.update(output)
my_state['time'] += model_time_step
t_warmup = time.perf_counter() - t_warmup_start

# Timed run
t0 = time.perf_counter()
for i in range(STEPS):
    diag, output = dycore(my_state, model_time_step)
    my_state.update(output)
    my_state['time'] += model_time_step
elapsed = time.perf_counter() - t0

numba_mode = "OFF" if os.environ.get("NUMBA_DISABLE_JIT") == "1" else "ON"
print(json.dumps({{
    "elapsed": elapsed,
    "warmup": t_warmup,
    "steps": STEPS,
    "per_step": elapsed / STEPS,
    "numba": numba_mode,
}}))
'''


def run_config(steps, disable_jit):
    env = os.environ.copy()
    if disable_jit:
        env["NUMBA_DISABLE_JIT"] = "1"
    else:
        env.pop("NUMBA_DISABLE_JIT", None)

    script = WORKER_SCRIPT.format(steps=steps)
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env=env,
        timeout=1800,
    )
    if result.returncode != 0:
        print(f"FAILED (exit {result.returncode}):", file=sys.stderr)
        for line in result.stderr.strip().split("\n")[-5:]:
            print(f"  {line}", file=sys.stderr)
        return None

    # GFS core may print debug lines to stdout; find the JSON line
    for line in reversed(result.stdout.strip().split("\n")):
        line = line.strip()
        if line.startswith("{"):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                continue
    print(f"No JSON in output: {result.stdout[:300]}", file=sys.stderr)
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Held-Suarez GCM benchmark (UnytBackend, Numba on/off)"
    )
    parser.add_argument("--steps", type=int, default=300, help="Timesteps to run")
    args = parser.parse_args()

    print(f"Held-Suarez GCM benchmark: {args.steps} steps, 128x62 grid, UnytBackend")
    print("=" * 70)

    for label, disable_jit in [("Numba OFF", True), ("Numba ON", False)]:
        print(f"\nRunning {label}...", flush=True)
        t_wall_start = time.perf_counter()
        data = run_config(args.steps, disable_jit)
        t_wall = time.perf_counter() - t_wall_start

        if data:
            print(f"  Warmup (1 step):  {data['warmup']:.2f}s")
            print(f"  Timed ({data['steps']} steps): {data['elapsed']:.2f}s")
            print(f"  Per step:         {data['per_step']:.4f}s")
            print(f"  Wall clock total: {t_wall:.1f}s")
        else:
            print(f"  FAILED")

    print()


if __name__ == "__main__":
    main()
