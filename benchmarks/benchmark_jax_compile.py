#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Measure JAX trace graph size, compile (first-call) time, and steady-state run time
for the Emanuel convection JAX kernel."""
import argparse
from datetime import datetime
from pathlib import Path
import platform
import statistics
import sys
import time

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

from climt._core.jax_support import count_jaxpr_equations  # noqa: E402
from tests.jax_cases import emanuel_case, emanuel_forward, emanuel_loss  # noqa: E402

RESULTS_PATH = Path(__file__).resolve().parent / "results" / "emanuel_jax_compile.txt"


def _measure(fn, arg, steady_runs=20):
    if hasattr(jax, "clear_caches"):
        jax.clear_caches()
    t0 = time.perf_counter()
    equations = count_jaxpr_equations(fn, arg)
    trace_s = time.perf_counter() - t0

    if hasattr(jax, "clear_caches"):
        jax.clear_caches()
    compiled = jax.jit(fn)
    t1 = time.perf_counter()
    jax.block_until_ready(compiled(arg))
    first_call_s = time.perf_counter() - t1

    steady_times = []
    for _ in range(steady_runs):
        ts = time.perf_counter()
        jax.block_until_ready(compiled(arg))
        steady_times.append(time.perf_counter() - ts)
    steady_median_s = statistics.median(steady_times)
    return equations, trace_s, first_call_s, steady_median_s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nlev", type=int, default=30)
    parser.add_argument("--ncol", type=int, default=1)
    parser.add_argument("--skip-grad", action="store_true")
    parser.add_argument("--label", type=str, default="")
    args = parser.parse_args()

    component, state, timestep = emanuel_case(nlev=args.nlev, ncol=args.ncol, profile="unstable", as_jax=True)
    temp = state["air_temperature"]

    fwd_fn = emanuel_forward(component, state, timestep)
    fwd_eqns, fwd_trace_s, fwd_first_s, fwd_steady_s = _measure(fwd_fn, temp)

    lines = [
        "=" * 76,
        f"label: {args.label or '(unlabeled)'}",
        f"date: {datetime.now().isoformat(timespec='seconds')}",
        f"jax_version: {jax.__version__}",
        f"machine: {platform.machine()}",
        f"nlev: {args.nlev}  ncol: {args.ncol}",
        (
            f"forward:  equations={fwd_eqns}  trace={fwd_trace_s:.3f}s  "
            f"first_call={fwd_first_s:.3f}s  steady_median={fwd_steady_s * 1e3:.3f}ms ({fwd_steady_s:.6f}s)"
        ),
    ]

    if not args.skip_grad:
        grad_fn = jax.grad(emanuel_loss(component, state, timestep))
        grad_eqns, grad_trace_s, grad_first_s, grad_steady_s = _measure(grad_fn, temp)
        lines.append(
            f"gradient: equations={grad_eqns}  trace={grad_trace_s:.3f}s  "
            f"first_call={grad_first_s:.3f}s  steady_median={grad_steady_s * 1e3:.3f}ms ({grad_steady_s:.6f}s)"
        )

    block = "\n".join(lines) + "\n"
    print(block, end="", flush=True)
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULTS_PATH.open("a", encoding="utf-8") as fh:
        fh.write(block)


if __name__ == "__main__":
    main()
