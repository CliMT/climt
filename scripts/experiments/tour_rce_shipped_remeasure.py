"""Re-measure Task 0's sensitivity numbers against the states that shipped.

Task 0 (`tour_rce_measurements.py`, `tour_rce_moist_probe.py`) measured the
2xCO2 response and start-independence from runs stopped by a TOA-only gate.
The shipped equilibria use a stricter one -- TOA balance AND a stationary
surface temperature -- and sit elsewhere: the dry state at 266.48 K rather
than 266.60 K. (The moist state was 279.97 K rather than 286.67 K until
2026-09-27, when it was regenerated over a surface saturated at its own
temperature, with a trend gate of its own; see the generator's MOIST_GATE.)
So those numbers belong to the wrong equilibria. This script takes the configuration,
the gate and the loop from `scripts/generate_tour_equilibria.py` itself, so
there is one definition of each, and measures:

    co2-dry      page 11's 2xCO2 from rce_dry_equilibrium.npz: warming and the
                 steps it takes (page 11 runs this live, so the steps matter).
    indep-dry    a warm, steep start (+40 K, 9 K/km) under the same gate,
                 compared with rce_dry_equilibrium.npz.
    indep-moist  the same for the moist column, against rce_moist_equilibrium.npz.

    settle-moist both shipped moist states stepped on 60 000 steps: where
                 each settles, and the settled 2xCO2 warming page 12 quotes.

Page 12's 2xCO2 *state* is not made here: it ships, as
rce_moist_2xco2_equilibrium.npz, from `generate_tour_equilibria.py
--moist-2xco2`. `settle-moist` measures the warming from it.

    python scripts/experiments/tour_rce_shipped_remeasure.py co2-dry indep-dry
    python scripts/experiments/tour_rce_shipped_remeasure.py indep-moist
    python scripts/experiments/tour_rce_shipped_remeasure.py settle-moist

Each measurement writes debug_data/tour_rce_remeasure_<name>.json.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
import sympl

import climt

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
_spec = importlib.util.spec_from_file_location(
    "generate_tour_equilibria",
    os.path.join(REPO, "scripts", "generate_tour_equilibria.py"))
gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gen)


def _surface(state):
    return float(state["surface_temperature"].values.ravel()[0])


def _profile(state):
    return state["air_temperature"].values[:, 0, 0].copy()


def co2_dry():
    tendencies, steppers, state, base = gen.load_equilibrium("dry")
    before = _surface(state)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = \
        gen.CO2_DOUBLING * base["co2_ppm"] * 1e-6
    started = time.time()
    n_steps = gen.run_to_equilibrium(tendencies, steppers, state,
                                     gen.DRY_DT_HOURS, print_every=500,
                                     **gen.gate_for("dry"))
    after = _surface(state)
    result = dict(before_k=before, after_k=after, warming_k=after - before,
                  n_steps=n_steps, days=n_steps * gen.DRY_DT_HOURS / 24.0,
                  toa_imbalance=gen.budgets.toa_imbalance(state),
                  native_wall_s=time.time() - started)
    print(f"  dry 2xCO2: {before:.3f} -> {after:.3f} K "
          f"({after - before:+.3f} K) in {n_steps} steps")
    return result


def _independence(kind):
    """A warm, steep start against the shipped state, under the shipped gate.

    The same perturbed start Task 0 used: 328 K at the surface lapsing at
    9 K/km (scale height 7 km), floored at 170 K, and a 330 K surface.
    """
    moist = (kind == "moist")
    components = gen.moist_components() if moist else gen.dry_components()
    dt_hours = gen.MOIST_DT_HOURS if moist else gen.DRY_DT_HOURS
    tendencies, steppers = gen.split(components)
    state, relaxation = gen.build_state(components, moist)
    p = state["air_pressure"].values[:, 0, 0]
    surface_pressure = float(state["surface_air_pressure"].values.ravel()[0])
    T = 328.0 - 9.0e-3 * 7000.0 * np.log(surface_pressure / p)
    state["air_temperature"].values[:, 0, 0] = np.maximum(T, 170.0)
    state["surface_temperature"].values[:] = 330.0

    started = time.time()
    n_steps = gen.run_to_equilibrium(
        tendencies + [relaxation], steppers, state, dt_hours,
        print_every=500 if not moist else 10000, **gen.gate_for(kind))

    _, _, shipped, provenance = gen.load_equilibrium(kind)
    difference = np.abs(_profile(state) - _profile(shipped))
    p_hpa = shipped["air_pressure"].values[:, 0, 0] / 100.0
    worst = int(np.argmax(difference))
    result = dict(
        n_steps=n_steps, days=n_steps * dt_hours / 24.0,
        surface_k=_surface(state), shipped_surface_k=_surface(shipped),
        shipped_n_steps=provenance["n_steps"],
        surface_difference_k=abs(_surface(state) - _surface(shipped)),
        max_profile_difference_k=float(difference.max()),
        max_difference_at_hpa=float(p_hpa[worst]),
        troposphere_max_difference_k=float(difference[p_hpa > 200.0].max()),
        toa_imbalance=gen.budgets.toa_imbalance(state),
        native_wall_s=time.time() - started)
    print(f"  {kind} warm-steep: {n_steps} steps, Tsurf {_surface(state):.3f} "
          f"vs shipped {_surface(shipped):.3f} K; column max |dT| "
          f"{difference.max():.3f} K at {p_hpa[worst]:.0f} hPa; below 200 hPa "
          f"{result['troposphere_max_difference_k']:.3f} K")
    return result


def settle_moist(n_steps=60000, window_steps=30000):
    """Step both shipped moist states on, and compare where they settle.

    The moist column never reaches TOA = 0. DryConvectiveAdjustment's
    moist-cp bookkeeping adds ~1 W/m^2 the budget does not count (page 12;
    Emanuel adds ~0.03), so the column settles with a steady TOA
    imbalance (about -1.1 W/m^2 over the saturated surface; it was +0.3 over
    the old fixed-humidity one). The generator's moist gate stops on flat
    60-day trends, so the files should sit where the column settles; this
    checks that, and measures page 12's 2xCO2 warming as the difference
    between the two *settled* surface temperatures: the mean over the last
    ``window_steps`` of ``n_steps``, from each shipped file.
    """
    settled = {}
    for label, filename in (("base", "rce_moist_equilibrium.npz"),
                            ("2xco2", "rce_moist_2xco2_equilibrium.npz")):
        tendencies, steppers, state, provenance = gen.load_equilibrium(
            "moist", filename=filename)
        timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
        rows = []
        for k in range(n_steps // 50):
            gen.stepping.integrate(tendencies, steppers, state, timestep, 50)
            rows.append(((k + 1) * 50, gen.budgets.toa_imbalance(state),
                         _surface(state)))
        rows = np.array(rows)
        tail = rows[rows[:, 0] > n_steps - window_steps]
        days = tail[:, 0] * provenance["dt_hours"] / 24.0
        settled[label] = dict(
            shipped_surface_k=float(provenance.get(
                "surface_temperature_k", rows[0, 2])),
            settled_surface_k=float(tail[:, 2].mean()),
            settled_toa=float(tail[:, 1].mean()),
            settled_toa_std=float(tail[:, 1].std()),
            trend_k_per_100d=float(100.0 * np.polyfit(days, tail[:, 2], 1)[0]))
        print(f"  {label:6s} settles at {settled[label]['settled_surface_k']:.4f} K, "
              f"TOA {settled[label]['settled_toa']:+.3f} "
              f"(std {settled[label]['settled_toa_std']:.3f}), trend "
              f"{settled[label]['trend_k_per_100d']:+.4f} K/100 d", flush=True)
    warming = (settled["2xco2"]["settled_surface_k"]
               - settled["base"]["settled_surface_k"])
    print(f"  settled 2xCO2 warming: {warming:+.3f} K")
    return dict(settled, warming_k=warming, n_steps=n_steps,
                window_steps=window_steps)


MEASUREMENTS = {
    "co2-dry": co2_dry,
    "indep-dry": lambda: _independence("dry"),
    "indep-moist": lambda: _independence("moist"),
    "settle-moist": settle_moist,
}


def main():
    sympl.set_backend(climt.UnytBackend())
    out_dir = os.path.join(REPO, "debug_data")
    os.makedirs(out_dir, exist_ok=True)
    for name in sys.argv[1:] or list(MEASUREMENTS):
        print(f"== {name}", flush=True)
        result = MEASUREMENTS[name]()
        result.update(climt_version=climt.__version__)
        path = os.path.join(out_dir, f"tour_rce_remeasure_{name}.json")
        with open(path, "w") as handle:
            json.dump(result, handle, indent=2)
        print(f"  wrote {path}", flush=True)


if __name__ == "__main__":
    main()
