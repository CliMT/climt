"""Measure every number the Modelling Tour's page 9 quotes.

Page 9 (``docs/modelling-tour/09-moisture-and-buoyancy.qmd``) takes page 8's
boundary-layer column, switches its radiation to the 14-band
``earth_low_res_lw`` table, and gives the slab a wet surface: a
``stepping.SurfaceHumidity`` stepper holds ``surface_specific_humidity`` at a
fixed relative humidity of the surface's own temperature, so the boundary
layer has a latent heat flux to compute. ``GridScaleCondensation`` follows the
boundary layer in the stepper list. The column is built here exactly as the
page builds it: nz = 28, a 2 m slab, 240 W m^-2 absorbed at the surface, the
wind relaxed toward 5 m/s over 24 h, z0 = 1e-3.

Physics runs use numba if it is there (the step count is a physics quantity,
numba-independent); ``cost`` must be run with ``NUMBA_DISABLE_JIT=1``, the
native stand-in for the browser, or its timings are worthless.

    python scripts/experiments/tour_page9_measurements.py cold
    python scripts/experiments/tour_page9_measurements.py timestep
    python scripts/experiments/tour_page9_measurements.py short
    python scripts/experiments/tour_page9_measurements.py restart
    python scripts/experiments/tour_page9_measurements.py supersat
    NUMBA_DISABLE_JIT=1 python scripts/experiments/tour_page9_measurements.py cost

Every number is a 30-day mean unless it says otherwise. The boundary layer's
depth saw-tooths from step to step, as on page 8, and the fluxes with it.
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
import soundings    # noqa: E402
import stepping     # noqa: E402

SOLAR = 240.0
DT = climt.UnytTimeDelta(hours=1)
COLD_DAYS = 1000
CHECKPOINTS = (30, 60, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)


def column(relative_humidity=1.0, condensation=True, condense_first=False,
           slab_depth=2.0):
    """Page 9's column. Returns (tendencies, steppers, state)."""
    sympl.set_backend(climt.UnytBackend())
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                               roughness_length=1e-3)
    moist = [boundary_layer]
    if condensation:
        if condense_first:
            moist.insert(0, climt.GridScaleCondensation())
        else:
            moist.append(climt.GridScaleCondensation())
    steppers = [stepping.SurfaceHumidity(relative_humidity)] + moist
    state = climt.get_default_state(
        [longwave, surface] + steppers,
        grid_state=climt.get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    relaxation = stepping.wind_relaxation(state, 5.0, timescale_hours=24.0)
    return [longwave, surface, relaxation], steppers, state


def _scalar(state, name):
    return float(np.asarray(state[name].values, dtype=float).ravel()[0])


def _theta_v_gap(state):
    """theta_v - theta at the lowest level, K."""
    Rd = float(sympl.get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Cp = float(sympl.get_constant(
        "heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"))
    p = float(state["air_pressure"].values[0, 0, 0])
    T = float(state["air_temperature"].values[0, 0, 0])
    q = float(state["specific_humidity"].values[0, 0, 0])
    return T * (1.0e5 / p) ** (Rd / Cp) * 0.61 * q


def recorder(dt=DT):
    return stepping.Recorder(
        dt,
        ts=lambda s: _scalar(s, "surface_temperature"),
        sh=budgets.sensible_heat_flux,
        lh=lambda s: _scalar(s, "surface_upward_latent_heat_flux"),
        jump=budgets.surface_air_jump,
        blh=lambda s: _scalar(s, "boundary_layer_height"),
        q0=lambda s: float(s["specific_humidity"].values[0, 0, 0]),
        peak_rh=lambda s: float(np.max(soundings.relative_humidity(s))),
        precip=lambda s: budgets.precipitation_rate(s, dt),
        toa=budgets.toa_imbalance,
        gap=_theta_v_gap,
    )


def window(record, day, days=30.0):
    """Means of every series over the ``days`` before ``day``."""
    elapsed = record["days"]
    mask = (elapsed > day - days) & (elapsed <= day + 1e-9)
    out = {name: float(np.mean(record[name][mask]))
           for name in ("ts", "sh", "lh", "jump", "blh", "q0", "peak_rh",
                        "precip", "toa", "gap")}
    out["blh_median"] = float(np.median(record["blh"][mask]))
    out["peak_rh_max"] = float(np.max(record["peak_rh"][mask]))
    out["bowen"] = out["sh"] / out["lh"] if out["lh"] else float("nan")
    return out


def line(label, w):
    return (f"{label:26s} Ts {w['ts']:7.2f}  SH {w['sh']:5.1f}  "
            f"LH {w['lh']:5.1f}  Bowen {w['bowen']:5.3f}  "
            f"jump {w['jump']:5.2f}  BLh mean {w['blh']:5.0f} "
            f"median {w['blh_median']:5.0f} m  q0 {w['q0'] * 1e3:5.2f} g/kg  "
            f"thv-th {w['gap']:4.2f} K  peakRH {w['peak_rh'] * 100:6.1f}% "
            f"(max {w['peak_rh_max'] * 100:6.1f}%)  P {w['precip']:5.2f} "
            f"mm/d  TOA {w['toa']:+6.2f}")


# ---------------------------------------------------------- cold starts

COLD = [
    ("RH 1.0 (headline)", dict()),
    ("RH 0.4", dict(relative_humidity=0.4)),
    ("RH 0 (dry surface)", dict(relative_humidity=0.0)),
    ("RH 1.0, no condensation", dict(condensation=False)),
    ("RH 1.0, condense first", dict(condense_first=True)),
]


def cold(spec):
    label, kwargs = spec
    tendencies, steppers, state = column(**kwargs)
    record = recorder()
    start = time.time()
    n_steps = COLD_DAYS * 24
    try:
        stepping.integrate(tendencies, steppers, state, DT, n_steps,
                           after_step=record)
    except Exception as error:     # noqa: BLE001 -- report and carry on
        return [f"{label}: FAILED after {len(record['days'])} steps: "
                f"{str(error)[:300]}"]
    lines = [f"{label}  [{time.time() - start:.0f} s]"]
    for day in CHECKPOINTS:
        lines.append(line(f"  day {day:4d}", window(record, day)))
    return lines


# ------------------------------------------------------ timestep choice

def timestep(minutes):
    """The partition's dependence on dt, cold start to 1000 days."""
    dt = climt.UnytTimeDelta(minutes=minutes)
    tendencies, steppers, state = column()
    record = recorder(dt)
    n_steps = int(round(COLD_DAYS * 24 * 60 / minutes))
    stepping.integrate(tendencies, steppers, state, dt, n_steps,
                       after_step=record)
    return [line(f"dt {minutes:4d} min x {n_steps}",
                 window(record, COLD_DAYS))]


# ------------------------------------------------ what the cells print

SHORT_DAYS = 30          # the headline and knob cells' cold starts
SHORT_MEAN_DAYS = 10.0   # what their prints average over
RESTART_DAYS = 10        # the condensation cell's continuation


def _short_run(kwargs, days=SHORT_DAYS):
    tendencies, steppers, state = column(**kwargs)
    record = recorder()
    stepping.integrate(tendencies, steppers, state, DT, days * 24,
                       after_step=record)
    return tendencies, steppers, state, record


def short(label_kwargs):
    """The page's cold-start cells, exactly: SHORT_DAYS at 1 h."""
    label, kwargs = label_kwargs
    _, _, _, record = _short_run(kwargs)
    return [line(f"{label}, day {SHORT_DAYS}",
                 window(record, SHORT_DAYS, days=SHORT_MEAN_DAYS))]


SHORT = [
    ("RH 1.0 (headline cell)", dict()),
    ("RH 0.4 (knob cell)", dict(relative_humidity=0.4)),
    ("RH 1.0, condense first", dict(condense_first=True)),
    ("RH 1.0, no condensation", dict(condensation=False)),
]


def restart(_):
    """The condensation cell: from the headline cell's day-30 state, step
    RESTART_DAYS more with the sink and without it."""
    tendencies, steppers, base, _ = _short_run(dict())
    lines = []
    for label, keep_sink in (("with condensation", True),
                             ("without condensation", False)):
        state = copy.deepcopy(base)
        here = steppers if keep_sink else steppers[:2]
        record = recorder()
        stepping.integrate(tendencies, here, state, DT, RESTART_DAYS * 24,
                           after_step=record)
        w = window(record, RESTART_DAYS, days=RESTART_DAYS)
        lines.append(line(f"restart +{RESTART_DAYS} d, {label}", w))
        for day in (1, 2, 5, 10):
            mask = record["days"] <= day
            lines.append(f"    peak RH by day {day:2d}: "
                         f"{np.max(record['peak_rh'][mask]) * 100:6.1f}%  "
                         f"(level of the peak now: "
                         f"{_peak_level(state):.0f} hPa)")
    return lines


def _peak_level(state):
    relative_humidity = soundings.relative_humidity(state)
    k = int(np.argmax(relative_humidity))
    return float(state["air_pressure"].values[k, 0, 0]) / 100.0


def supersat(_):
    """No sink at all, from a cold start: how far and how fast the column
    supersaturates, and when it stops running."""
    tendencies, steppers, state = column(condensation=False)
    record = recorder()
    try:
        stepping.integrate(tendencies, steppers, state, DT, COLD_DAYS * 24,
                           after_step=record)
        ending = "ran to the end"
    except Exception as error:     # noqa: BLE001 -- report and carry on
        ending = (f"FAILED after {len(record['days'])} steps "
                  f"({len(record['days']) / 24:.1f} d): {str(error)[:300]}")
    lines = ["no condensation, cold start: " + ending]
    for day in (1, 2, 5, 10, 30, 60, 100, 200, 300, 360):
        mask = record["days"] <= day
        if not mask.any() or record["days"][-1] < day:
            break
        at = np.searchsorted(record["days"], day)
        lines.append(f"    day {day:3d}: peak RH {record['peak_rh'][at] * 100:8.1f}%"
                     f"  q0 {record['q0'][at] * 1e3:6.2f} g/kg  "
                     f"LH {record['lh'][at]:5.1f}  Ts {record['ts'][at]:6.2f}")
    return lines


# ------------------------------------------------------------------ cost

def cost(n=30, warmup=5):
    """Per-step cost of page 9's stack. Run with NUMBA_DISABLE_JIT=1."""
    tendencies, steppers, state = column()
    stepping.integrate(tendencies, steppers, state, DT, warmup)
    start = time.time()
    stepping.integrate(tendencies, steppers, state, DT, n)
    per_step = (time.time() - start) / n
    return [f"page 9 stack: {per_step * 1e3:.1f} ms/step native "
            f"(~{per_step * 3.5:.3f} s/step browser), "
            f"NUMBA_DISABLE_JIT={os.environ.get('NUMBA_DISABLE_JIT', '0')}"]


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "cold"
    if which == "cold":
        jobs, run = COLD, cold
    elif which == "timestep":
        jobs, run = [30, 60, 180, 360], timestep
    elif which == "short":
        jobs, run = SHORT, short
    elif which == "restart":
        jobs, run = [None], restart
    elif which == "supersat":
        jobs, run = [None], supersat
    elif which == "cost":
        for text in cost():
            print(text)
        return
    else:
        raise SystemExit(f"unknown measurement {which!r}")
    with Pool(min(4, len(jobs))) as pool:
        for lines in pool.imap(run, jobs):
            print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
