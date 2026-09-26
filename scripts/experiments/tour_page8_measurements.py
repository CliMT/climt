"""Measure every number the Modelling Tour's page 8 quotes.

Page 8 (``docs/modelling-tour/08-turbulent-heat-exchange.qmd``) adds
``SimpleBoundaryLayer(surface_fluxes='bulk')`` to page 7's gray column and
holds a wind up with ``stepping.wind_relaxation``. Its column is built here
exactly as the page builds it: ``tour_gray_lw`` at diffusivity 2, nz = 28,
a 2 m slab, 240 W m^-2 absorbed at the surface, dt = 1 h.

Physics runs use numba if it is there (the step count is a physics quantity,
numba-independent); ``cost`` must be run with ``NUMBA_DISABLE_JIT=1``, the
native stand-in for the browser, or its timings are worthless.

    python scripts/experiments/tour_page8_measurements.py all
    NUMBA_DISABLE_JIT=1 python scripts/experiments/tour_page8_measurements.py cost

Why 30-day means. At equilibrium the boundary layer's depth is chosen level
by level, and every six hours or so it briefly mixes a kilometre or two deep
before collapsing back to ~550 m. The surface--air jump and the sensible heat
flux flicker with it, by about +-0.3 K and +-2 W m^-2, so a single step's
value is one sample of a sawtooth. Every number below is reported both ways:
the last step's value, and the mean over the last 30 days.
"""
import copy
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import sympl

import climt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "docs", "modelling-tour", "_tour"))
import budgets      # noqa: E402
import stepping     # noqa: E402

SOLAR = 240.0
DT = climt.UnytTimeDelta(hours=1)
COLD_STEPS = 12000          # 500 days
RESTART_STEPS = (1000, 2000, 3000)


def column(z0=1e-3, wind=5.0, boundary_layer=True, slab_depth=2.0,
           timescale_hours=24.0):
    """Page 8's column. ``wind=None``: no momentum source."""
    sympl.set_backend(climt.UnytBackend())
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table="tour_gray_lw", diffusivity_factor=2.0)
    surface = climt.SlabSurface()
    steppers = ([climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                           roughness_length=z0)]
                if boundary_layer else [])
    state = climt.get_default_state(
        [longwave, surface] + steppers,
        grid_state=climt.get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    tendencies = [longwave, surface]
    if boundary_layer and wind is not None:
        tendencies.append(stepping.wind_relaxation(
            state, wind, timescale_hours=timescale_hours))
    return tendencies, steppers, state


def lowest_wind(state):
    if "eastward_wind" not in state:       # the radiative-only column
        return 0.0
    return float(state["eastward_wind"].values[0, 0, 0])


def recorder(dt=DT):
    return stepping.Recorder(dt, jump=budgets.surface_air_jump,
                             flux=budgets.sensible_heat_flux,
                             wind=lowest_wind,
                             ts=lambda s: float(
                                 s["surface_temperature"].values.ravel()[0]))


def summarise(label, record, state, wall=None):
    line = (f"{label:34s} last: jump {record['jump'][-1]:6.2f} "
            f"SH {record['flux'][-1]:5.1f} u0 {record['wind'][-1]:5.2f} "
            f"Ts {record['ts'][-1]:7.2f} | 30 d: jump "
            f"{record.mean('jump'):6.2f} SH {record.mean('flux'):5.1f} "
            f"u0 {record.mean('wind'):5.2f} Ts {record.mean('ts'):7.2f} | "
            f"TOA {budgets.toa_imbalance(state):+.3f}")
    if wall is not None:
        line += f"  [{wall:.0f} s]"
    return line


def hold(speed, record):
    def after_step(state):
        state["eastward_wind"].values[:] = speed
        record(state)
    return after_step


# ---------------------------------------------------------- cold starts

def cold(spec):
    """One 12 000-step cold start. ``spec`` is (label, column kwargs, extra)."""
    label, kwargs, extra = spec
    tendencies, steppers, state = column(**kwargs)
    if "initial_wind" in extra:
        state["eastward_wind"].values[:] = extra["initial_wind"]
    record = recorder()
    after = hold(extra["reset"], record) if "reset" in extra else record
    start = time.time()
    try:
        return _cold_run(label, tendencies, steppers, state, record, after,
                         extra, start)
    except ValueError as error:
        return (f"{label:34s} FAILED after {len(record['days'])} steps "
                f"({len(record['days']) / 24:.1f} d): {str(error)[:60]}")


def _cold_run(label, tendencies, steppers, state, record, after, extra,
              start):
    if extra.get("euler_reset"):
        # A loop of one-step integrate(..., 1) calls, for comparison.
        # Each call builds a fresh AdamsBashforth, so this is forward Euler.
        for _ in range(COLD_STEPS):
            stepping.integrate(tendencies, steppers, state, DT, 1)
            state["eastward_wind"].values[:] = extra["euler_reset"]
            record(state)
    else:
        stepping.integrate(tendencies, steppers, state, DT, COLD_STEPS,
                           after_step=after)
    return summarise(label, record, state, time.time() - start)


COLD = [
    ("radiative only (dt 1 h)", dict(boundary_layer=False), {}),
    ("relax 2 m/s", dict(wind=2.0), {}),
    ("relax 5 m/s", dict(wind=5.0), {}),
    ("relax 10 m/s", dict(wind=10.0), {}),
    ("hard reset 5 m/s (AB3)", dict(wind=None), dict(reset=5.0)),
    ("hard reset 5 m/s (Euler loop)", dict(wind=None), dict(euler_reset=5.0)),
    ("no source, u0=0", dict(wind=None), dict(initial_wind=0.0)),
    ("no source, u0=5", dict(wind=None), dict(initial_wind=5.0)),
    ("no source, u0=10", dict(wind=None), dict(initial_wind=10.0)),
    ("z0 3.21e-5, 5 m/s", dict(z0=3.21e-5), {}),
    ("z0 1e-1, 5 m/s", dict(z0=1e-1), {}),
    ("z0 3.21e-5, calm", dict(z0=3.21e-5, wind=None), {}),
    ("z0 1e-1, calm", dict(z0=1e-1, wind=None), {}),
    ("tau 100 d, 5 m/s", dict(timescale_hours=2400.0), {}),
    # Page 8's third physics exercise. tau = dt puts the relaxation outside
    # AdamsBashforth's stability region (|dt/tau| < 6/11 on the real axis for
    # the third-order scheme), so this one is expected to blow up.
    ("tau 1 h, 5 m/s", dict(timescale_hours=1.0), {}),
    ("tau 2 h, 5 m/s", dict(timescale_hours=2.0), {}),
]


# ------------------------------------------------------ timestep choice

def timestep(hours):
    """Why page 8 steps hourly: the turbulent equilibrium depends on dt."""
    dt = climt.UnytTimeDelta(minutes=int(hours * 60))
    tendencies, steppers, state = column()
    record = recorder(dt)
    n_steps = int(round(COLD_STEPS / hours))
    stepping.integrate(tendencies, steppers, state, dt, n_steps,
                       after_step=record)
    return summarise(f"relax 5 m/s, dt {hours:g} h x {n_steps}", record,
                     state)


# ----------------------------------------------- restarts from equilibrium

def restart(_):
    """Page 8's knob cells start from the headline's equilibrium. Check that
    RESTART_STEPS from there lands where a 12 000-step cold start does."""
    tendencies, steppers, base = column()
    stepping.integrate(tendencies, steppers, base, DT, COLD_STEPS)
    lines = []
    cases = [("relax 2 m/s", dict(wind=2.0)),
             ("relax 10 m/s", dict(wind=10.0)),
             ("hard reset 5 m/s", dict(reset=5.0)),
             ("z0 3.21e-5", dict(z0=3.21e-5)),
             ("z0 1e-1", dict(z0=1e-1))]
    for label, case in cases:
        state = copy.deepcopy(base)
        longwave, surface = tendencies[:2]
        steppers_here = steppers
        if "z0" in case:
            steppers_here = [climt.SimpleBoundaryLayer(
                surface_fluxes="bulk", roughness_length=case["z0"])]
        record = recorder()
        if "reset" in case:
            here, after = [longwave, surface], hold(case["reset"], record)
        else:
            wind = case.get("wind", 5.0)
            here = [longwave, surface,
                    stepping.wind_relaxation(state, wind, initialise=False)]
            after = record
        done = 0
        for target in RESTART_STEPS:
            stepping.integrate(here, steppers_here, state, DT, target - done,
                               after_step=after)
            done = target
            lines.append(summarise(f"restart {label} +{target}", record,
                                   state))
    return "\n".join(lines)


# ------------------------------------------------------------- spin-down

def spindown(_):
    """How fast the drag removes a 10 m/s wind with no momentum source."""
    tendencies, steppers, state = column(wind=None)
    state["eastward_wind"].values[:] = 10.0
    record = recorder()
    stepping.integrate(tendencies, steppers, state, DT, 24 * 200,
                       after_step=record)
    days, wind = record["days"], record["wind"]
    out = [f"spin-down from 10 m/s, no source (lowest level):"]
    for day in (1, 2, 5, 10, 20, 30, 60, 100, 200):
        out.append(f"  day {day:3d}: {wind[np.argmin(abs(days - day))]:.4f}")
    for level in (1.0, 0.1, 0.01, 0.001):
        below = np.nonzero(np.abs(wind) < level)[0]
        out.append(f"  below {level} m/s from day "
                   f"{days[below[0]] if len(below) else float('nan'):.1f}")
    column_max = float(np.abs(state["eastward_wind"].values).max())
    out.append(f"  column max |u| at day 200: {column_max:.4f}; levels still "
               f"above 9 m/s: {int((state['eastward_wind'].values > 9).sum())}"
               f" of {state['eastward_wind'].values.shape[0]}")
    return "\n".join(out)


# ------------------------------------------------------------------ cost

def cost(_):
    """Per-step cost of each page-8 loop. Run with NUMBA_DISABLE_JIT=1."""
    out = [f"NUMBA_DISABLE_JIT={os.environ.get('NUMBA_DISABLE_JIT', '')}"]
    for label, kwargs, dt_hours in (
            ("radiative only (dt 12 h)", dict(boundary_layer=False), 12),
            ("boundary layer + relaxation", dict(), 1),
            ("boundary layer, no source", dict(wind=None), 1)):
        tendencies, steppers, state = column(**kwargs)
        timestep = climt.UnytTimeDelta(hours=dt_hours)
        record = (recorder() if kwargs.get("boundary_layer", True)
                  else None)
        stepping.integrate(tendencies, steppers, state, timestep, 5,
                           after_step=record)
        start = time.time()
        stepping.integrate(tendencies, steppers, state, timestep, 300,
                           after_step=record)
        per_step = (time.time() - start) / 300
        out.append(f"  {label:30s} {per_step * 1e3:6.2f} ms/step native "
                   f"(~{3.5 * per_step * 1e3:.0f} ms/step browser)")
    return "\n".join(out)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which == "cost":
        print(cost(None))
        return
    jobs = []
    if which in ("all", "cold"):
        jobs += [(cold, spec) for spec in COLD]
    if which == "tau":
        jobs += [(cold, spec) for spec in COLD if spec[0].startswith("tau")]
    if which in ("all", "timestep"):
        jobs += [(timestep, hours) for hours in (0.5, 12.0)]
    if which in ("all", "restart"):
        jobs.append((restart, None))
    if which in ("all", "spindown"):
        jobs.append((spindown, None))
    with Pool(int(os.environ.get("JOBS", os.cpu_count() or 1))) as pool:
        results = [pool.apply_async(func, (arg,)) for func, arg in jobs]
        for result in results:
            print(result.get(), flush=True)


if __name__ == "__main__":
    main()
