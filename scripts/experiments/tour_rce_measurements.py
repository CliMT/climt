"""Measure everything the Modelling Tour's RCE pages need to quote.

The spec's Task 0: no page is written until these five are known. Run under
NUMBA_DISABLE_JIT=1 for anything cost-related -- Pyodide has no numba, so a
JIT-compiled timing is off by ~20x on the 14-band radiation and worthless as a
browser estimate.

    NUMBA_DISABLE_JIT=1 conda run -n climt python \\
        scripts/experiments/tour_rce_measurements.py all

Individual measurements:

    ... tour_rce_measurements.py stability     # 1: largest stable timestep
    ... tour_rce_measurements.py convergence   # 2: steps to equilibrium
    ... tour_rce_measurements.py cost          # 3: per-step cost per page
    ... tour_rce_measurements.py independence  # 4: start-independence
    ... tour_rce_measurements.py supersat      # 5: what condensation removes

Named measurements each write their own shard,
``debug_data/tour_rce_measurements_<name>.json``, so they can be run in
parallel; ``... tour_rce_measurements.py merge`` combines the shards into
``debug_data/tour_rce_measurements.json``. ``all`` writes the combined file
directly, after every measurement, so a late crash loses nothing.
"""
import glob
import json
import os
import sys
import time

import numpy as np
import sympl

import climt

sys.path.insert(0, os.path.join("docs", "modelling-tour", "_tour"))
import budgets      # noqa: E402
import stepping     # noqa: E402

OUT_DIR = "debug_data"
SOLAR = 240.0
NZ = 28
SLAB_DEPTH = 2.0
WIND = 5.0           # m/s, the wind the relaxation holds the column at
Z0 = 1e-3            # surface roughness length, m
# Page 07's gray half must reproduce page 04's analytic profile, and that only
# works against the table page 04 calibrated: tour_gray_lw at diffusivity 2.0.
GRAY = "tour_gray_lw"
GRAY_DIFFUSIVITY = 2.0
BAND14 = "earth_low_res_lw"


def column(table, with_boundary_layer=False, with_dry_adjustment=False,
           with_condensation=False, with_emanuel=False,
           saturated_surface=None, diffusivity_factor=None,
           slab_depth=SLAB_DEPTH, nz=NZ, co2_ppm=None,
           wind=WIND, roughness_length=Z0):
    """Build one page's exact stack. Returns (tendencies, steppers, state)."""
    if diffusivity_factor is None:
        longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                               table=table)
    else:
        longwave = climt.CorkLongwaveRadiation(
            optics="correlated_k", table=table,
            diffusivity_factor=diffusivity_factor)
    surface = climt.SlabSurface()
    tendencies = [longwave, surface]
    steppers = []
    if with_boundary_layer:
        steppers.append(climt.SimpleBoundaryLayer(
            surface_fluxes="bulk", roughness_length=roughness_length))
    if with_dry_adjustment:
        steppers.append(climt.DryConvectiveAdjustment())
    if with_emanuel:
        tendencies.append(climt.EmanuelConvectionPython())
    if with_condensation:
        steppers.append(climt.GridScaleCondensation())
    if saturated_surface is None:
        saturated_surface = with_condensation or with_emanuel

    grid = climt.get_grid(nx=1, ny=1, nz=nz)
    state = climt.get_default_state(tendencies + steppers, grid_state=grid)
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    if co2_ppm is not None:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = (
            co2_ppm * 1e-6)
    if saturated_surface:
        # A saturated surface, so the latent flux is on. Page 9's knob.
        state["surface_specific_humidity"].values[:] = 0.015
    if with_boundary_layer:
        # Without this the column spins itself down to dead calm within days
        # and its surface fluxes are those of a windless planet. See
        # _tour/stepping.wind_relaxation.
        tendencies.append(stepping.wind_relaxation(state, wind))
    return tendencies, steppers, state


STACKS = {
    "07-gray":   dict(table=GRAY, diffusivity_factor=GRAY_DIFFUSIVITY),
    "07-14band": dict(table=BAND14),
    "08-bl":     dict(table=GRAY, diffusivity_factor=GRAY_DIFFUSIVITY,
                      with_boundary_layer=True),
    "09-moist-bl": dict(table=BAND14, with_boundary_layer=True,
                        with_condensation=True),
    "11-dry-rce": dict(table=BAND14, with_boundary_layer=True,
                       with_dry_adjustment=True),
    "12-moist-rce": dict(table=BAND14, with_boundary_layer=True,
                         with_dry_adjustment=True, with_condensation=True,
                         with_emanuel=True),
}


# ---------------------------------------------------------- 1: stability

def stability(hours=(0.25, 0.5, 1, 2, 3, 6, 12, 24), n_steps=200):
    """Largest timestep each stack survives n_steps without going non-finite.

    "Survives" is: no exception (Task 2's guard raises on non-finite fluxes),
    no NaN anywhere, and a surface temperature still in [150, 400] K.
    """
    results = {}
    for name, kwargs in STACKS.items():
        largest = None
        for dt_hours in hours:
            tendencies, steppers, state = column(**kwargs)
            timestep = climt.UnytTimeDelta(hours=dt_hours)
            try:
                stepping.integrate(tendencies, steppers, state, timestep,
                                   n_steps)
                T = state["air_temperature"].values
                surface = float(state["surface_temperature"].values.ravel()[0])
                ok = (np.all(np.isfinite(T)) and 150.0 < surface < 400.0)
            except Exception as error:                   # noqa: BLE001
                ok = False
                print(f"    {name} dt={dt_hours}h -> {type(error).__name__}: "
                      f"{str(error)[:70]}")
            if ok:
                largest = dt_hours
            else:
                break
        results[name] = largest
        print(f"  {name:14s} largest stable dt = {largest} h "
              f"({n_steps} steps)")
    return results


# -------------------------------------------------------- 2: convergence

def steps_to_equilibrium(tendencies, steppers, state, timestep,
                         threshold=0.5, max_steps=6000, check_every=25):
    """Steps until |TOA imbalance| stays under `threshold` W/m^2.

    Returns (n_steps, final_toa, final_surface_temperature) or
    (None, ...) if it never gets there.
    """
    for step in range(0, max_steps, check_every):
        stepping.integrate(tendencies, steppers, state, timestep, check_every)
        imbalance = budgets.toa_imbalance(state)
        if abs(imbalance) < threshold:
            return (step + check_every, imbalance,
                    float(state["surface_temperature"].values.ravel()[0]))
    return (None, budgets.toa_imbalance(state),
            float(state["surface_temperature"].values.ravel()[0]))


def convergence(dt_hours=None):
    """Steps to equilibrium from a cold start, per stack, at its stable dt.

    `dt_hours` overrides the per-stack choice; otherwise each stack uses the
    largest dt `stability()` cleared, capped at 12 h (beyond which the
    snapshots on page 7 are too coarse to show an approach).
    """
    sigma = float(sympl.get_constant("stefan_boltzmann_constant",
                                     "W/m^2/degK^4"))
    emission_temperature = (SOLAR / sigma) ** 0.25
    stable = stability()
    results = {}
    for name, kwargs in STACKS.items():
        dt = dt_hours or min(stable[name] or 1, 12)
        tendencies, steppers, state = column(**kwargs)
        started = time.time()
        n_steps, imbalance, surface = steps_to_equilibrium(
            tendencies, steppers, state, climt.UnytTimeDelta(hours=dt))
        elapsed = time.time() - started
        air = state["air_temperature"].values[:, 0, 0]
        discontinuity = surface - float(air[0])
        top = float(air[-1])
        results[name] = dict(dt_hours=dt, n_steps=n_steps,
                             toa_imbalance=imbalance,
                             surface_temperature=surface,
                             surface_air_discontinuity=discontinuity,
                             top_air_temperature=top,
                             skin_ratio=top / emission_temperature,
                             wall_seconds=elapsed)
        days = None if n_steps is None else n_steps * dt / 24.0
        print(f"  {name:14s} dt={dt:5.2f}h  {str(n_steps):>6s} steps "
              f"({'' if days is None else f'{days:.0f} d'})  "
              f"TOA {imbalance:+7.3f}  Tsurf {surface:6.2f} K  "
              f"dT(sfc-air) {discontinuity:+6.2f} K  "
              f"Ttop {top:6.2f} K "
              f"(Ttop/Te = {top / emission_temperature:.4f})  "
              f"[{elapsed:.1f} s native]")
    return results


def co2_doubling_response(name="11-dry-rce", dt_hours=1.0, threshold=0.5):
    """From equilibrium, double CO2 and re-converge: the sensitivity."""
    kwargs = STACKS[name]
    tendencies, steppers, state = column(**kwargs)
    steps_to_equilibrium(tendencies, steppers, state,
                         climt.UnytTimeDelta(hours=dt_hours),
                         threshold=threshold)
    before = float(state["surface_temperature"].values.ravel()[0])
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= 2.0
    n_steps, imbalance, after = steps_to_equilibrium(
        tendencies, steppers, state, climt.UnytTimeDelta(hours=dt_hours),
        threshold=threshold)
    print(f"  {name}: 2xCO2  {before:.2f} -> {after:.2f} K "
          f"(dT = {after - before:+.2f} K) in {n_steps} steps at "
          f"dt = {dt_hours} h")
    return dict(before=before, after=after, warming=after - before,
                n_steps=n_steps, toa_imbalance=imbalance)


# ---------------------------------------------------------------- 3: cost

def cost(n=30, warmup=5, dt_hours=1.0):
    """Per-step wall clock for each page's exact component list."""
    results = {}
    for name, kwargs in STACKS.items():
        tendencies, steppers, state = column(**kwargs)
        timestep = climt.UnytTimeDelta(hours=dt_hours)
        stepping.integrate(tendencies, steppers, state, timestep, warmup)
        started = time.time()
        stepping.integrate(tendencies, steppers, state, timestep, n)
        per_step = (time.time() - started) / n * 1000.0
        results[name] = per_step
        print(f"  {name:14s} {per_step:7.2f} ms/step native "
              f"(~{per_step * 3.5 / 1000:.3f} s/step browser)")
    return results


# ------------------------------------------------------ 4: independence

def independence(name="11-dry-rce", dt_hours=1.0, threshold=0.5,
                 tolerance=0.5):
    """Two different initial conditions must reach the same equilibrium.

    Load-bearing: shipping an equilibrium state is only honest if the
    equilibrium does not depend on how you got there, and page 12 says so in
    prose. Start A is climt's default profile; start B is 40 K warmer at the
    surface with a steeper lapse rate.
    """
    kwargs = STACKS[name]
    finals = {}
    for label in ("default", "warm-steep"):
        tendencies, steppers, state = column(**kwargs)
        # Note: the wind is NOT part of the initial condition to vary. It is
        # relaxed to a fixed target, so it is a boundary condition here, the
        # same in both runs. Varying it would test something else.
        if label == "warm-steep":
            p = state["air_pressure"].values[:, 0, 0]
            surface_pressure = float(
                state["surface_air_pressure"].values.ravel()[0])
            T = 328.0 - 9.0e-3 * 7000.0 * np.log(surface_pressure / p)
            state["air_temperature"].values[:, 0, 0] = np.maximum(T, 170.0)
            state["surface_temperature"].values[:] = 330.0
        n_steps, imbalance, surface = steps_to_equilibrium(
            tendencies, steppers, state, climt.UnytTimeDelta(hours=dt_hours),
            threshold=threshold, max_steps=12000)
        finals[label] = dict(
            n_steps=n_steps, toa_imbalance=imbalance,
            surface_temperature=surface,
            profile=state["air_temperature"].values[:, 0, 0].copy())
        print(f"  {label:11s} {str(n_steps):>6s} steps  "
              f"Tsurf {surface:6.2f} K  TOA {imbalance:+7.3f}")

    difference = np.abs(finals["default"]["profile"]
                        - finals["warm-steep"]["profile"])
    surface_difference = abs(finals["default"]["surface_temperature"]
                             - finals["warm-steep"]["surface_temperature"])
    print(f"  max |dT| through the column: {difference.max():.3f} K; "
          f"surface: {surface_difference:.3f} K")
    verdict = ("START-INDEPENDENT" if difference.max() < tolerance
               else "PATH-DEPENDENT — DO NOT SHIP A STATE")
    print(f"  -> {verdict} (tolerance {tolerance} K)")
    return dict(max_profile_difference=float(difference.max()),
                surface_difference=float(surface_difference),
                verdict=verdict,
                default=finals["default"]["surface_temperature"],
                warm_steep=finals["warm-steep"]["surface_temperature"])


# ---------------------------------------------------------- 5: supersat

def supersaturation(dt_hours=1.0, n_steps=200):
    """How much supersaturation GridScaleCondensation removes on page 9.

    Page 9's stack is the boundary layer over a saturated surface with 14-band
    radiation. Without a moisture sink the column supersaturates: the BL mixes
    moist air into colder air aloft. This runs the stack twice -- once with
    condensation, once without -- and reports the peak relative humidity each
    reaches, and the precipitation the sink produces.
    """
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Rv = float(get_constant("gas_constant_of_vapor_phase", "J/kg/degK"))
    epsilon = Rd / Rv

    def saturation_specific_humidity(T, p):
        # Bolton (1980), the same form _tour/soundings.py uses.
        e_sat = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
        return epsilon * e_sat / np.maximum(p - (1.0 - epsilon) * e_sat, 1.0)

    results = {}
    for label, with_sink in (("with condensation", True),
                             ("without condensation", False)):
        tendencies, steppers, state = column(
            BAND14, with_boundary_layer=True, saturated_surface=True,
            with_condensation=with_sink)

        peak = 0.0
        timestep = climt.UnytTimeDelta(hours=dt_hours)
        for _ in range(n_steps):
            stepping.integrate(tendencies, steppers, state, timestep, 1)
            T = state["air_temperature"].values[:, 0, 0]
            q = state["specific_humidity"].values[:, 0, 0]
            p = state["air_pressure"].values[:, 0, 0]
            relative_humidity = q / saturation_specific_humidity(T, p)
            peak = max(peak, float(np.max(relative_humidity)))
        precipitation = budgets.precipitation_rate(state, timestep)
        sensible = float(
            state["surface_upward_sensible_heat_flux"].values.ravel()[0])
        latent = float(
            state["surface_upward_latent_heat_flux"].values.ravel()[0])
        bowen = sensible / latent if latent != 0.0 else float("nan")
        results[label] = dict(peak_relative_humidity=peak,
                              precipitation_mm_day=precipitation,
                              sensible_heat_flux=sensible,
                              latent_heat_flux=latent,
                              bowen_ratio=bowen)
        print(f"  {label:22s} peak RH {peak * 100:6.1f}%   "
              f"precip {precipitation:6.3f} mm/day   "
              f"SH {sensible:7.2f}  LH {latent:7.2f}  "
              f"Bowen {bowen:6.3f}")
    return results


# -------------------------------------------------------------------- driver

MEASUREMENTS = {
    "stability": stability,
    "convergence": convergence,
    "cost": cost,
    "independence": independence,
    "independence-moist": lambda: independence(name="12-moist-rce",
                                               dt_hours=0.25),
    "supersat": supersaturation,
    "co2": co2_doubling_response,
}

STEP_ORDER = ["cost", "stability", "convergence", "co2", "independence",
              "independence-moist", "supersat"]


def _write(path, collected):
    with open(path, "w") as handle:
        json.dump(collected, handle, indent=2, default=float)


def merge():
    """Combine the per-measurement shards into one file."""
    collected = {}
    pattern = os.path.join(OUT_DIR, "tour_rce_measurements_*.json")
    shards = sorted(glob.glob(pattern))
    for shard in shards:
        stem = os.path.basename(shard)[:-len(".json")]
        name = stem[len("tour_rce_measurements_"):]
        with open(shard) as handle:
            collected[name] = json.load(handle)
        print(f"  found {shard}")
    path = os.path.join(OUT_DIR, "tour_rce_measurements.json")
    _write(path, collected)
    print(f"\nwrote {path} from {len(shards)} shard(s)")


def main():
    sympl.set_backend(climt.UnytBackend())
    os.makedirs(OUT_DIR, exist_ok=True)
    if os.environ.get("NUMBA_DISABLE_JIT") != "1":
        print("WARNING: numba JIT is ON. Cost numbers will not reflect the "
              "browser. Re-run with NUMBA_DISABLE_JIT=1.\n")

    requested = sys.argv[1:] or ["all"]
    if requested == ["merge"]:
        merge()
        return

    everything = requested == ["all"]
    names = STEP_ORDER if everything else requested
    collected = {}
    for name in names:
        print(f"\n== {name} ==")
        result = MEASUREMENTS[name]()
        collected[name] = result
        # Written after every measurement: a crash late in a two-hour campaign
        # must not lose what already ran.
        if everything:
            path = os.path.join(OUT_DIR, "tour_rce_measurements.json")
            _write(path, collected)
        else:
            path = os.path.join(OUT_DIR,
                                f"tour_rce_measurements_{name}.json")
            _write(path, result)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
