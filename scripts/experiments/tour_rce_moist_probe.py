"""Measure the moist RCE stack (page 12) at its real timestep, 5 min.

The main campaign (tour_rce_measurements.py) runs at NUMBA_DISABLE_JIT=1 for
browser-cost numbers and near 6 h, where the moist column is unstable. This probe
runs WITH numba at dt = 5 min -- the repo's Emanuel convention -- to get the moist
stack's steps-to-equilibrium, 2xCO2 warming, and start-independence, which the
6 h campaign could not.
"""
import json, os, sys, time
import numpy as np
import sympl, climt

sys.path.insert(0, os.path.join("scripts", "experiments"))
import tour_rce_measurements as m   # reuse column(), STACKS, steps_to_equilibrium, budgets

DT = climt.UnytTimeDelta(minutes=5)
MAX_STEPS = 400_000        # ~1400 sim-days at 5 min; a safety cap, not a target
CHECK_EVERY = 200
THRESHOLD = 0.5            # W/m^2 on |TOA imbalance|, same as the campaign


def run_to_equilibrium(state, tendencies, steppers, label):
    started = time.time()
    n_steps, imbalance, surface = m.steps_to_equilibrium(
        tendencies, steppers, state, DT,
        threshold=THRESHOLD, max_steps=MAX_STEPS, check_every=CHECK_EVERY)
    elapsed = time.time() - started
    days = None if n_steps is None else n_steps * 5.0 / 60.0 / 24.0
    print(f"  {label:12s} n_steps={n_steps} "
          f"({'--' if days is None else f'{days:.0f} d'}) "
          f"TOA={imbalance:+.3f} Tsurf={surface:.2f}K [{elapsed:.0f}s]")
    return n_steps, imbalance, surface


def main():
    sympl.set_backend(climt.UnytBackend())
    os.makedirs("debug_data", exist_ok=True)
    kwargs = m.STACKS["12-moist-rce"]
    out = {"dt_minutes": 5, "numba": "on"}

    # A. default start -> equilibrium
    print("== default start ==")
    tendencies, steppers, state = m.column(**kwargs)
    try:
        n, toa, tsurf = run_to_equilibrium(state, tendencies, steppers, "default")
    except Exception as e:
        out["default"] = {"blew_up": str(e)[:200]}
        json.dump(out, open("debug_data/tour_rce_moist_probe.json", "w"),
                  indent=2, default=float)
        print("BLEW UP on default start:", str(e)[:200]); return
    default_profile = state["air_temperature"].values[:, 0, 0].copy()
    out["default"] = dict(n_steps=n, toa=toa, surface_temperature=tsurf)

    # B. 2xCO2 from that equilibrium
    print("== 2xCO2 ==")
    before = float(state["surface_temperature"].values.ravel()[0])
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= 2.0
    n2, toa2, after = run_to_equilibrium(state, tendencies, steppers, "2xCO2")
    out["co2_doubling"] = dict(before=before, after=after, warming=after - before,
                               n_steps=n2, toa=toa2)

    # C. warm-steep start -> equilibrium, then start-independence
    print("== warm-steep start ==")
    tendencies2, steppers2, state2 = m.column(**kwargs)
    p = state2["air_pressure"].values[:, 0, 0]
    sp = float(state2["surface_air_pressure"].values.ravel()[0])
    T = 328.0 - 9.0e-3 * 7000.0 * np.log(sp / p)
    state2["air_temperature"].values[:, 0, 0] = np.maximum(T, 170.0)
    state2["surface_temperature"].values[:] = 330.0
    nw, toaw, tsurfw = run_to_equilibrium(state2, tendencies2, steppers2, "warm-steep")
    warm_profile = state2["air_temperature"].values[:, 0, 0].copy()
    max_dT = float(np.abs(default_profile - warm_profile).max())
    surf_dT = abs(tsurf - tsurfw)
    verdict = ("START-INDEPENDENT" if max_dT < 0.5
               else "PATH-DEPENDENT (or one run did not converge)")
    print(f"  max|dT|={max_dT:.3f}K surface|dT|={surf_dT:.3f}K -> {verdict}")
    out["independence"] = dict(warm_steep=dict(n_steps=nw, surface_temperature=tsurfw),
                               max_profile_difference=max_dT,
                               surface_difference=surf_dT, verdict=verdict)

    json.dump(out, open("debug_data/tour_rce_moist_probe.json", "w"),
              indent=2, default=float)
    print("\nwrote debug_data/tour_rce_moist_probe.json")


if __name__ == "__main__":
    main()
