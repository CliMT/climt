"""Benchmark Cork-LW (consolidated njit(parallel) transport) vs RRTMG-LW.

Wall-clock per call over NCOL columns, warm JIT. Also profiles the CORK call to
split time between the optical-depth assembly (still pure-Python: interpolate_k +
additive tau loop + continuum/CO2 interpolation) and the prange transport kernel
(lw_transport), so we can see what the parallel kernel actually buys end-to-end.

Run: /Users/joymonteiro/miniconda3/envs/climt/bin/python \
        scripts/experiments/bench_cork_vs_rrtmg.py [NCOL] [NZ]
"""
import argparse
import cProfile
import io
import os
import pstats
import sys
import time

import numpy as np
import sympl

from climt import get_default_state, get_grid

# Optional NCOL/NZ are read straight from sys.argv (argparse has no positionals
# registered, so it would reject them). The isdigit guard lets flags like --save
# fall through to argparse in main(). Constraint: don't pass NCOL/NZ positionals
# in the same invocation as --save (argparse would reject the leftover tokens).
NCOL = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].lstrip("-").isdigit() else 100
NZ = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].lstrip("-").isdigit() else 30
N_WARM = 2
N_TIME = 20


def _time_component(lw, state, n):
    # Warm up (numba JIT compile on first call).
    for _ in range(N_WARM):
        lw(state)
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        lw(state)
        ts.append(time.perf_counter() - t0)
    ts = np.array(ts) * 1e3  # ms
    return ts.mean(), ts.std(), ts.min()


def _make_state(lw):
    grid = get_grid(nx=NCOL, ny=1, nz=NZ)
    state = get_default_state([lw], grid_state=grid)
    if "specific_humidity" in state:
        # Simple decreasing-with-height humidity so the optics do real work.
        q = np.linspace(8e-3, 1e-6, NZ)[:, None, None]
        state["specific_humidity"].values[:] = q
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 400e-6
    return state


def main():
    ap = argparse.ArgumentParser(description="Benchmark CORK-LW vs RRTMG-LW throughput")
    ap.add_argument("--save", default=None,
                    help="if given, write throughput scalars to this .npz "
                         "(keys: rrtmg_us_per_col, cork_us_per_col, ncol)")
    args = ap.parse_args()

    sympl.set_backend(sympl.DataArrayBackend())
    nthreads = os.environ.get("NUMBA_NUM_THREADS", "(default = all cores)")
    print(f"=== LW benchmark: NCOL={NCOL}, NZ={NZ}, "
          f"NUMBA_NUM_THREADS={nthreads} ===\n")

    # --- RRTMG-LW (Fortran) ---
    from climt import RRTMGLongwave
    rrtmg = RRTMGLongwave()
    st_r = _make_state(rrtmg)
    mean_r, std_r, min_r = _time_component(rrtmg, st_r, N_TIME)
    print(f"RRTMG-LW          : {mean_r:7.2f} ± {std_r:5.2f} ms/call  "
          f"(min {min_r:6.2f})  -> {mean_r/NCOL*1e3:6.1f} µs/col")

    # --- Cork-LW (correlated-k, 14-band, prange transport) ---
    from climt._components.cork import CorkLongwaveRadiation
    cork = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    st_p = _make_state(cork)
    mean_p, std_p, min_p = _time_component(cork, st_p, N_TIME)
    print(f"CORK-LW (14b x 8g)  : {mean_p:7.2f} ± {std_p:5.2f} ms/call  "
          f"(min {min_p:6.2f})  -> {mean_p/NCOL*1e3:6.1f} µs/col")

    print(f"\nPF / RRTMG ratio  : {mean_p/mean_r:5.2f}x")

    if args.save:
        # mean_* are ms/call; /NCOL*1e3 -> µs/col.
        np.savez(args.save,
                 rrtmg_us_per_col=np.asarray(mean_r / NCOL * 1e3),
                 cork_us_per_col=np.asarray(mean_p / NCOL * 1e3),
                 ncol=np.asarray(NCOL))
        print(f"saved throughput -> {args.save}")

    # --- Where does CORK spend its time? Profile one warm call. ---
    print("\n--- CORK call profile (cumulative top functions) ---")
    pr = cProfile.Profile()
    pr.enable()
    cork(st_p)
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats(12)
    # Trim to the informative lines.
    for line in s.getvalue().splitlines():
        if any(k in line for k in ("correlated_k", "kernels", "lw_transport",
                                   "interpolate_", "compute_ck", "array_call",
                                   "ncalls", "function calls")):
            print(line)


if __name__ == "__main__":
    main()
