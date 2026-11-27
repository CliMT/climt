#!/usr/bin/env python
"""Measure traced graph size and compile/steady wall time for Emanuel JAX kernel."""

import argparse
import datetime
import platform
import statistics
import time
from pathlib import Path

import jax

jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

from climt._core.jax_support import count_jaxpr_equations  # noqa: E402
from tests.jax_cases import emanuel_case, emanuel_forward, emanuel_loss  # noqa: E402

RESULTS_PATH = Path(__file__).resolve().parent / "results" / "emanuel_jax_compile.txt"


def measure(function, arg, n_steady=20):
    jax.clear_caches()
    t0 = time.perf_counter()
    eqns = count_jaxpr_equations(function, arg)
    trace_s = time.perf_counter() - t0

    jax.clear_caches()
    jitted = jax.jit(function)
    t0 = time.perf_counter()
    out = jitted(arg)
    jax.block_until_ready(out)
    first_call_s = time.perf_counter() - t0

    steady_times = []
    for _ in range(n_steady):
        t0 = time.perf_counter()
        out = jitted(arg)
        jax.block_until_ready(out)
        steady_times.append(time.perf_counter() - t0)
    median_s = statistics.median(steady_times)
    return eqns, trace_s, first_call_s, median_s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nlev", type=int, default=30)
    parser.add_argument("--ncol", type=int, default=1)
    parser.add_argument("--skip-grad", action="store_true")
    parser.add_argument("--label", type=str, default="")
    args = parser.parse_args()

    component, state, timestep = emanuel_case(nlev=args.nlev, ncol=args.ncol, profile="unstable", as_jax=True)
    T0 = state["air_temperature"]

    fwd_fn = emanuel_forward(component, state, timestep)
    fwd_eqns, fwd_trace, fwd_first, fwd_med = measure(fwd_fn, T0)

    lines = [
        "=" * 72,
        f"date:    {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"label:   {args.label or '(none)'}",
        f"jax:     {jax.__version__}",
        f"machine: {platform.machine()}",
        f"shape:   nlev={args.nlev}, ncol={args.ncol}",
        (
            f"forward: eqns={fwd_eqns:6d}  trace={fwd_trace:7.3f}s  "
            f"first_call={fwd_first:7.3f}s  steady_median={fwd_med * 1e3:7.3f}ms"
        ),
    ]

    if not args.skip_grad:
        grad_fn = jax.grad(emanuel_loss(component, state, timestep))
        grd_eqns, grd_trace, grd_first, grd_med = measure(grad_fn, T0)
        lines.append(
            f"grad:    eqns={grd_eqns:6d}  trace={grd_trace:7.3f}s  "
            f"first_call={grd_first:7.3f}s  steady_median={grd_med * 1e3:7.3f}ms"
        )

    block = "\n".join(lines) + "\n"
    print(block, end="")
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with RESULTS_PATH.open("a", encoding="utf-8") as handle:
        handle.write(block)


if __name__ == "__main__":
    main()
