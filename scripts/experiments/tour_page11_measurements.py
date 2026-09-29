"""Measure every number the Modelling Tour's page 11 quotes.

Page 11 (``docs/modelling-tour/11-dry-rce.qmd``) loads the shipped dry
radiative-convective equilibrium, ``rce_dry_equilibrium.npz``, and perturbs
it. Its component list is ``scripts/generate_tour_equilibria.py``'s
``dry_components()``, and the wind relaxation is rebuilt on the loaded state
with ``initialise=False``, exactly as the page does it (through
``generate_tour_equilibria.load_equilibrium``).

    python scripts/experiments/tour_page11_measurements.py all
    python scripts/experiments/tour_page11_measurements.py convection radiative
    NUMBA_DISABLE_JIT=1 python scripts/experiments/tour_page11_measurements.py cost

Physics runs use numba if it is there; ``cost`` must be run with
``NUMBA_DISABLE_JIT=1``, the native stand-in for the browser.

Why 30-day means. With the boundary layer diffusing dry static energy the
shipped equilibrium is a fixed point: the surface temperature does not move in
its third decimal place from step to step, and the convecting layer runs from
the lowest level to 370 hPa on every step. Everything below is still reported
as the mean over the last 30 days (60 steps), which is what the page compares,
and which stays honest for a perturbed run that has not finished wobbling.

    convection  the shipped state stepped 120 steps: where dry adjustment acts,
                the lapse rate there, the surface temperature
    radiative   page 7's 14-band column (LW + slab, diffusivity 1.66, from
                climt's default state) stepped to TOA balance: the radiative
                equilibrium page 11 draws dashed
    co2         2x and 0.5x CO2 from the shipped state, 1000 steps each, as the
                page's knob cell runs it
    settle      the shipped state stepped 10 000 more steps: how much its top
                level still cools
    order       the craft callout: the adjustment before the boundary layer
    noadjust    code exercise 2: the shipped state with DryConvectiveAdjustment
                removed, 1000 steps
    cost        seconds per step of the page's stack, JIT off
"""
import copy
import importlib.util
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
stepping, budgets = gen.stepping, gen.budgets

DT = climt.UnytTimeDelta(hours=gen.DRY_DT_HOURS)
MONTH = 60                       # 30 days at 12 h
RADIATIVE_STEPS = 10000          # page 7's column, to TOA balance


def _constants():
    from sympl import get_constant
    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Cp = float(get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                            "J/kg/degK"))
    g = float(get_constant("gravitational_acceleration", "m/s^2"))
    return Rd, Cp, g


def theta(state):
    Rd, Cp, _ = _constants()
    T = state["air_temperature"].values[:, 0, 0]
    p = state["air_pressure"].to_units("Pa").values[:, 0, 0]
    return T * (1.0e5 / p) ** (Rd / Cp)


def lapse_rate(state):
    """-dT/dz in K/km between mid levels, hypsometric thickness."""
    Rd, _, g = _constants()
    T = state["air_temperature"].values[:, 0, 0]
    p = state["air_pressure"].to_units("Pa").values[:, 0, 0]
    dz = Rd * 0.5 * (T[:-1] + T[1:]) / g * np.log(p[:-1] / p[1:])
    return -np.diff(T) / dz * 1000.0


def neutral(state):
    """Level pairs sharing one theta: the layer dry adjustment has mixed."""
    return np.abs(np.diff(theta(state))) < 1e-6


def surface(state):
    return float(state["surface_temperature"].values.ravel()[0])


def convection(n_steps=2 * MONTH):
    tendencies, steppers, state, _ = gen.load_equilibrium("dry")
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    print("shipped state, as loaded:")
    mask = neutral(state)
    levels = np.where(mask)[0]
    print(f"  neutral between levels {levels[0]} and {levels[-1] + 1} "
          f"({p[levels[0]]:.0f} to {p[levels[-1] + 1]:.0f} hPa), lapse "
          f"{lapse_rate(state)[mask].mean():.3f} K/km; Ts {surface(state):.3f}"
          f"; lowest air {state['air_temperature'].values[0, 0, 0]:.2f} K")
    print("  lapse rate by layer (K/km):",
          np.round(lapse_rate(state), 2).tolist())
    print("  T (K):", np.round(state["air_temperature"].values[:, 0, 0],
                               2).tolist())

    adjustment = steppers[-1]
    model = sympl.AdamsBashforth(*tendencies)
    tops, bottoms, lapses, surfaces, toas = [], [], [], [], []
    for _ in range(n_steps):
        diagnostics, new_state = model(state, DT)
        state.update(new_state)
        state.update(diagnostics)
        for stepper in steppers[:-1]:
            d, s = stepper(state, DT)
            state.update(s)
            state.update(d)
        before = state["air_temperature"].values[:, 0, 0].copy()
        d, s = adjustment(state, DT)
        state.update(s)
        state.update(d)
        state["time"] += DT
        touched = np.where(np.abs(
            state["air_temperature"].values[:, 0, 0] - before) > 1e-8)[0]
        tops.append(touched[-1])
        bottoms.append(touched[0])
        lapses.append(lapse_rate(state)[neutral(state)].mean())
        surfaces.append(surface(state))
        toas.append(budgets.toa_imbalance(state))
    tops, bottoms = np.array(tops[-MONTH:]), np.array(bottoms[-MONTH:])
    print(f"last 30 days ({MONTH} steps):")
    for name, levels in (("top", tops), ("bottom", bottoms)):
        counts = {f"{p[k]:.0f} hPa": int((levels == k).sum())
                  for k in np.unique(levels)}
        print(f"  {name} of the adjusted layer: {counts}")
    print(f"  lapse rate across the neutral layer: mean "
          f"{np.mean(lapses[-MONTH:]):.3f}, range "
          f"{np.min(lapses[-MONTH:]):.3f}-{np.max(lapses[-MONTH:]):.3f} K/km")
    print(f"  surface temperature: mean {np.mean(surfaces[-MONTH:]):.3f} K, "
          f"range {np.min(surfaces[-MONTH:]):.3f}-"
          f"{np.max(surfaces[-MONTH:]):.3f}")
    print(f"  TOA imbalance: mean {np.mean(toas[-MONTH:]):+.3f} W/m^2")


def radiative_column():
    """Page 7's non-grey column, built exactly as page 7 builds it."""
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table="earth_low_res_lw",
        diffusivity_factor=1.66)
    slab = climt.SlabSurface()
    state = climt.get_default_state([longwave, slab],
                                    grid_state=climt.get_grid(nx=1, ny=1,
                                                              nz=gen.NZ))
    state["ocean_mixed_layer_thickness"].values[:] = gen.SLAB_DEPTH_M
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = gen.SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, slab], state


def radiative():
    components, state = radiative_column()
    done = 0
    for n in (900, 1400, 2000, RADIATIVE_STEPS - 4300):
        stepping.integrate(components, [], state, DT, n)
        done += n
        print(f"  page 7 column, {done:5d} steps: Ts {surface(state):.3f} K, "
              f"TOA {budgets.toa_imbalance(state):+.3f}, top "
              f"{state['air_temperature'].values[-1, 0, 0]:.2f} K", flush=True)
    T_re = state["air_temperature"].values[:, 0, 0].copy()
    print("  radiative-equilibrium T (K), bottom first:")
    print("  " + repr(np.round(T_re, 1).tolist()))

    _, _, rce, _ = gen.load_equilibrium("dry")
    T = rce["air_temperature"].values[:, 0, 0]
    p = rce["air_pressure"].values[:, 0, 0] / 100.0
    difference = T - T_re
    k = int(np.argmax(difference))
    print(f"  RCE minus RE: surface {surface(rce) - surface(state):+.2f} K; "
          f"largest {difference.max():+.2f} K at {p[k]:.0f} hPa; lowest air "
          f"{difference[0]:+.2f}; 100 hPa level ({p[22]:.0f}) "
          f"{difference[22]:+.2f}; top {difference[-1]:+.2f} K")
    print(f"  RE surface-air jump {surface(state) - T_re[0]:.2f} K; RE lapse "
          f"rate in the lowest layer {lapse_rate(state)[0]:.1f} K/km, "
          f"largest {lapse_rate(state).max():.1f} K/km")

    # Why a much warmer troposphere costs the surface so little: one
    # longwave call on hybrids of the two states, temperatures held fixed.
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table=gen.TABLE)

    def olr(column):
        return float(longwave(column)[1][
            "upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    hybrid = copy.deepcopy(rce)
    hybrid["air_temperature"].values[:, 0, 0] = T_re
    print(f"  OLR, one LW call: RCE {olr(rce):.2f}; RCE surface under RE air "
          f"{olr(hybrid):.2f} W/m^2 (air warming worth "
          f"{olr(rce) - olr(hybrid):+.2f})")
    # The split: chill the *surface* and hold the air fixed. Chilling the air
    # instead changes its temperature-dependent k-distribution and makes the
    # column more transparent (it gives 244.06, more than the whole OLR).
    sigma = 5.670374419e-8
    Ts = surface(rce)
    cold_ground = copy.deepcopy(rce)
    cold_ground["surface_temperature"].values[:] = 1.0
    air_only = olr(cold_ground)
    print(f"  air's own emission to space (surface at 1 K): {air_only:.2f}; "
          f"surface shining through {olr(rce) - air_only:.2f} W/m^2 = "
          f"{(olr(rce) - air_only) / (sigma * Ts ** 4):.0%} of sigma Ts^4")
    warm_ground = copy.deepcopy(rce)
    warm_ground["surface_temperature"].values[:] = Ts + 1.0
    direct = olr(warm_ground) - olr(rce)
    print(f"  surface +1 K, air fixed: OLR {direct:+.2f} W/m^2 = "
          f"{direct / (4 * sigma * Ts ** 3):.0%} of 4 sigma Ts^3; x "
          f"{surface(state) - Ts:.2f} K = "
          f"{(surface(state) - Ts) * direct:.2f}")
    cold_air = copy.deepcopy(rce)
    cold_air["air_temperature"].values[:] = 1.0
    print(f"  (the wrong split, air at 1 K: {olr(cold_air):.2f} W/m^2)")


def _perturbed(factor, n_steps=1000, drop_adjustment=False,
               drop_steppers=False):
    tendencies, steppers, state, provenance = gen.load_equilibrium("dry")
    if drop_steppers:
        # Without the boundary layer nothing updates the turbulent fluxes, and
        # SlabSurface would go on losing the last step's sensible heat flux.
        steppers = []
        state["surface_upward_sensible_heat_flux"].values[:] = 0.0
        state["surface_upward_latent_heat_flux"].values[:] = 0.0
    if drop_adjustment:
        steppers = [s for s in steppers
                    if not isinstance(s, climt.DryConvectiveAdjustment)]
    before = surface(state)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= factor
    record = stepping.Recorder(
        DT, surface=surface, toa=budgets.toa_imbalance,
        mixing_top=lambda s: float(
            s["boundary_layer_height"].values.ravel()[0])
        if "boundary_layer_height" in s else np.nan)
    stepping.integrate(tendencies, steppers, state, DT, n_steps,
                       after_step=record)
    for step in range(100, n_steps + 1, 100):
        print(f"    step {step:5d}  Ts {record['surface'][step - 1]:.3f} K  "
              f"TOA {record['toa'][step - 1]:+.3f}")
    return before, state, record


def co2():
    # The base, 30-day-mean: the shipped state stepped on, unperturbed.
    _, _, base_record = _perturbed(1.0, n_steps=MONTH)
    base_mean = base_record.mean("surface", days=30)
    for factor in (2.0, 0.5):
        print(f"  CO2 x{factor}:")
        before, state, record = _perturbed(factor)
        after_mean = record.mean("surface", days=30)
        print(f"  CO2 x{factor}: file {before:.3f} K -> last step "
              f"{surface(state):.3f} K ({surface(state) - before:+.3f}); "
              f"30-day means {base_mean:.3f} -> {after_mean:.3f} K "
              f"({after_mean - base_mean:+.3f} K); as the page defines it, "
              f"30-day mean minus the file, {after_mean - before:+.3f} K; "
              f"last-step TOA "
              f"{budgets.toa_imbalance(state):+.3f}, 30-day TOA "
              f"{record.mean('toa', days=30):+.3f}, surface imbalance "
              f"{budgets.surface_imbalance(state):+.3f} W/m^2", flush=True)


def noadjust():
    before, state, record = _perturbed(1.0, drop_adjustment=True)
    T = state["air_temperature"].values[:, 0, 0]
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    components, radiative_state = radiative_column()
    stepping.integrate(components, [], radiative_state, DT, RADIATIVE_STEPS)
    T_re = radiative_state["air_temperature"].values[:, 0, 0]
    print(f"  boundary layer height, last 30 days: mean "
          f"{record.mean('mixing_top'):.0f} m, max "
          f"{record['mixing_top'][-MONTH:].max():.0f} m")
    print(f"  no adjustment, 1000 steps: Ts {before:.3f} -> "
          f"{surface(state):.3f} K (30-day {record.mean('surface'):.3f}), "
          f"TOA {budgets.toa_imbalance(state):+.3f}")
    print("  lapse rate (K/km):", np.round(lapse_rate(state), 1).tolist())
    difference = T - T_re
    print(f"  minus page 7's RE: surface "
          f"{surface(state) - surface(radiative_state):+.2f} K, air "
          f"{difference.min():+.2f} to {difference.max():+.2f} K "
          f"(largest at {p[np.argmax(np.abs(difference))]:.0f} hPa)")

    # The control: take the boundary layer out as well.
    _, state, _ = _perturbed(1.0, drop_steppers=True)
    difference = state["air_temperature"].values[:, 0, 0] - T_re
    print(f"  no boundary layer either, 1000 steps: air minus page 7's RE "
          f"{difference.min():+.2f} to {difference.max():+.2f} K")


def order():
    """The craft callout: the two steppers the other way round."""
    _, _, base_record = _perturbed(1.0, n_steps=MONTH)
    tendencies, steppers, state, _ = gen.load_equilibrium("dry")
    shipped_unstable = []
    stepping.integrate(tendencies, steppers, state, DT, MONTH,
                       after_step=lambda column: shipped_unstable.append(
                           bool(np.any(np.diff(theta(column)) < -1e-6))))
    print(f"  shipped order: steps ending with theta falling somewhere: "
          f"{sum(shipped_unstable)} of {MONTH}")
    tendencies, steppers, state, _ = gen.load_equilibrium("dry")
    swapped = steppers[::-1]
    print("  steppers:", [type(s).__name__ for s in swapped])
    record = stepping.Recorder(DT, surface=surface,
                               toa=budgets.toa_imbalance)
    unstable = []

    def after(column):
        record(column)
        unstable.append(np.where(np.diff(theta(column)) < -1e-6)[0])

    stepping.integrate(tendencies, swapped, state, DT, 1000, after_step=after)
    th = theta(state)
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    print(f"  adjustment first, 1000 steps: 30-day T_s "
          f"{record.mean('surface'):.3f} K against "
          f"{base_record.mean('surface'):.3f} K as shipped; TOA "
          f"{record.mean('toa'):+.3f}")
    where = sorted({int(round(p[k])) for ks in unstable[-MONTH:] for k in ks})
    print(f"  steps ending with theta falling somewhere: "
          f"{sum(len(k) > 0 for k in unstable[-MONTH:])} of the last {MONTH}"
          f"; at {where} hPa; lowest thetas now {np.round(th[:5], 2).tolist()}")


def settle(n_steps=10000):
    """How far the shipped state still has to go: 10 000 more steps."""
    tendencies, steppers, state, _ = gen.load_equilibrium("dry")
    start = state["air_temperature"].values[:, 0, 0].copy()
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    stepping.integrate(tendencies, steppers, state, DT, n_steps)
    change = state["air_temperature"].values[:, 0, 0] - start
    print(f"  after {n_steps} more steps, dT by level (K):")
    for k in range(len(p) - 1, 17, -1):
        print(f"    {p[k]:7.1f} hPa  {change[k]:+.2f}")


def cost(n=20, warmup=3):
    tendencies, steppers, state, _ = gen.load_equilibrium("dry")
    stepping.integrate(tendencies, steppers, state, DT, warmup)
    started = time.perf_counter()
    stepping.integrate(tendencies, steppers, state, DT, n)
    per_step = (time.perf_counter() - started) / n
    jit = "off" if os.environ.get("NUMBA_DISABLE_JIT") == "1" else "ON"
    print(f"  page 11 stack, JIT {jit}: {1e3 * per_step:.1f} ms/step native, "
          f"~{3.5 * per_step:.2f} s/step browser; 1000 steps ~"
          f"{3.5 * per_step * 1000 / 60:.1f} min browser")
    # Cells 0-2: load + one adjustment call.
    started = time.perf_counter()
    tendencies, steppers, state, _ = gen.load_equilibrium("dry")
    copy.deepcopy(state)
    print(f"  load: {time.perf_counter() - started:.2f} s native")


MEASUREMENTS = dict(convection=convection, radiative=radiative, co2=co2,
                    noadjust=noadjust, order=order, settle=settle,
                    cost=cost)


def main():
    sympl.set_backend(climt.UnytBackend())
    names = sys.argv[1:] or ["all"]
    if names == ["all"]:
        names = [n for n in MEASUREMENTS if n != "cost"]
    for name in names:
        print(f"== {name}", flush=True)
        MEASUREMENTS[name]()


if __name__ == "__main__":
    main()
