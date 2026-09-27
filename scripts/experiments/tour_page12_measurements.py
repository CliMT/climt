"""Measure every number the Modelling Tour's page 12 quotes.

Page 12 (``docs/modelling-tour/12-moist-rce.qmd``) loads the shipped moist
radiative-convective equilibrium, ``rce_moist_equilibrium.npz``, and its
doubled-CO2 twin, ``rce_moist_2xco2_equilibrium.npz``. Its component list is
``scripts/generate_tour_equilibria.py``'s ``moist_components()``, and the wind
relaxation is rebuilt on the loaded state with ``initialise=False``, exactly as
the page does it (through ``generate_tour_equilibria.load_equilibrium``).

    python scripts/experiments/tour_page12_measurements.py month energy
    python scripts/experiments/tour_page12_measurements.py timestep
    NUMBA_DISABLE_JIT=1 python scripts/experiments/tour_page12_measurements.py cost

Physics runs use numba if it is there; ``cost`` must be run with
``NUMBA_DISABLE_JIT=1``, the native stand-in for the browser. Set
``NUMBA_NUM_THREADS=1`` when running several at once.

Why 30-day means. The moist column is noisy: grid-scale condensation rains in
single-step bursts (up to ~160 mm/day for one 5-minute step), and the surface
temperature wanders by ~0.05 K. One step's value is one sample of that, so
everything below is a mean, over 30 days unless it says otherwise.

    month      each shipped state stepped 30 days at 5 min: surface
               temperature, P and E, the fluxes, the lapse-rate profile
    energy     5 days: the column enthalpy (cp_dry T + Lv q) each component
               adds, W/m^2 -- where the settled TOA residual comes from
    heating    5 days: each stepper's heating by level, where the column is
               saturated, and Emanuel's CAPE and mass flux
    timestep   the shipped state stepped 60 days at 2.5, 5, 10 and 20 min
    nogsc      code exercise 2: GridScaleCondensation removed, 300 days
    nodca      DryConvectiveAdjustment removed, 300 days
    feedback   radiation-only: the 2xCO2 forcing and the OLR change from the
               temperature and from the humidity, each alone
    adiabat    moist adiabats: the lapse rate on this column's temperatures,
               and what saturated parcels average from 1000 to 500 hPa
    cost       seconds per step of the page's stack, JIT off
"""
import copy
import importlib.util
import os
import sys
import time
import warnings

import numpy as np
import sympl

import climt

REPO = os.path.join(os.path.dirname(__file__), "..", "..")
_spec = importlib.util.spec_from_file_location(
    "generate_tour_equilibria",
    os.path.join(REPO, "scripts", "generate_tour_equilibria.py"))
gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gen)
stepping, budgets, states = gen.stepping, gen.budgets, gen.states
import soundings  # noqa: E402  (on sys.path via the generator)

STEPS_PER_DAY = 288              # at 5 min
MONTH = 30 * STEPS_PER_DAY
DOUBLED = "rce_moist_2xco2_equilibrium.npz"


def _constants():
    from sympl import get_constant

    def c(name, units):
        return float(get_constant(name, units))
    return dict(
        Rd=c("gas_constant_of_dry_air", "J/kg/degK"),
        Cp=c("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"),
        g=c("gravitational_acceleration", "m/s^2"),
        Lv=c("latent_heat_of_condensation", "J/kg"),
    )


def column(state, name):
    values = np.asarray(state[name].values, dtype=float)
    return values.reshape(values.shape[0], -1)[:, 0].copy()


def scalar(state, name):
    return float(np.asarray(state[name].values).ravel()[0])


def lapse_rate(T, p):
    """-dT/dz in K/km between adjacent levels, hypsometric thickness."""
    c = _constants()
    dz = (c["Rd"] * 0.5 * (T[:-1] + T[1:]) / c["g"]) * np.log(p[:-1] / p[1:])
    return -np.diff(T) / dz * 1000.0


def lower_troposphere(p):
    """The layers below 500 hPa -- the page's lower troposphere."""
    return p[:-1] > 5.0e4


def load(filename=None):
    return gen.load_equilibrium("moist", filename=filename)


def recorder(timestep):
    return stepping.Recorder(
        timestep, profiles=("air_temperature", "specific_humidity"),
        profile_every=1,
        surface=lambda s: scalar(s, "surface_temperature"),
        precipitation=lambda s: budgets.precipitation_rate(s, timestep),
        convective=lambda s: scalar(s, "convective_precipitation_rate"),
        mass_flux=lambda s: scalar(s, "cloud_base_mass_flux"),
        evaporation=budgets.evaporation_rate,
        sensible=budgets.sensible_heat_flux,
        latent=budgets.latent_heat_flux,
        toa=budgets.toa_imbalance)


def report(label, record, p, days=30.0):
    elapsed = record["days"]
    last = elapsed > elapsed[-1] - days
    mean = {name: float(np.mean(record[name][last]))
            for name in ("surface", "precipitation", "convective", "mass_flux",
                         "evaporation", "sensible", "latent", "toa")}
    T = record["air_temperature"][last].mean(axis=0)
    q = record["specific_humidity"][last].mean(axis=0)
    lapse = lapse_rate(T, p)
    print(f"  {label}: last {days:.0f} days of {elapsed[-1]:.0f}")
    print(f"    Ts {mean['surface']:.3f} K  P {mean['precipitation']:.3f} "
          f"(Emanuel {mean['convective']:.3f}) E {mean['evaporation']:.3f} "
          f"mm/day  P-E {mean['precipitation'] - mean['evaporation']:+.3f}")
    print(f"    SH {mean['sensible']:.2f} LH {mean['latent']:.2f} W/m^2, Bowen "
          f"{mean['sensible'] / mean['latent']:.2f}; TOA {mean['toa']:+.3f}; "
          f"cloud-base mass flux {mean['mass_flux']:.2e}")
    p_mid = np.sqrt(p[:-1] * p[1:])
    band = (p_mid < 79220.0) & (p_mid > 5.0e4)   # above the adjusted layer
    print(f"    lapse below 500 hPa {lapse[lower_troposphere(p)].mean():.2f}, "
          f"792 to 500 hPa {lapse[band].mean():.2f} K/km; lowest q "
          f"{1e3 * q[0]:.2f} g/kg")
    print("    lapse by layer (K/km): "
          + " ".join(f"{x:.1f}" for x in lapse[:18]))
    return mean, T, q, lapse


def run(filename=None, minutes=5.0, days=30.0, drop=()):
    tendencies, steppers, state, _ = load(filename)
    tendencies = [c for c in tendencies if type(c).__name__ not in drop]
    steppers = [c for c in steppers if type(c).__name__ not in drop]
    timestep = climt.UnytTimeDelta(minutes=minutes)
    record = recorder(timestep)
    n_steps = int(round(days * 1440 / minutes))
    start = time.time()
    stepping.integrate(tendencies, steppers, state, timestep, n_steps,
                       after_step=record)
    print(f"  ({n_steps} steps, {time.time() - start:.0f} s)")
    return record, column(state, "air_pressure")


def month():
    """30 days of each shipped state, at the shipped 5 min."""
    base, p = run()
    base_mean, T0, q0, lapse0 = report("base", base, p)
    doubled, _ = run(DOUBLED)
    doubled_mean, T1, q1, lapse1 = report("2xCO2", doubled, p)
    print(f"  2xCO2 warming, 30-day mean minus 30-day mean: "
          f"{doubled_mean['surface'] - base_mean['surface']:+.3f} K")
    print("  level  p (hPa)  T base  dT    q base  q ratio  RH base")
    RH = q0 / soundings.saturation_specific_humidity(T0, p)
    for k in range(20):
        print(f"  {k:5d} {p[k] / 100:8.1f} {T0[k]:7.2f} {T1[k] - T0[k]:+5.2f} "
              f"{1e3 * q0[k]:6.3f} {q1[k] / q0[k] if q0[k] > 0 else np.nan:7.2f}"
              f" {100 * RH[k]:6.0f}")
    for days in (1, 2, 5):
        n = days * STEPS_PER_DAY
        blocks = len(base["precipitation"]) // n
        residual = (base["precipitation"][:blocks * n].reshape(blocks, n).mean(1)
                    - base["evaporation"][:blocks * n].reshape(blocks, n).mean(1))
        print(f"  base, {days}-day means of P - E: sd {residual.std():.3f}, "
              f"range {residual.min():+.3f} to {residual.max():+.3f} mm/day")


def _enthalpy_rate(state, tendencies):
    """Column cp_dry dT/dt + Lv dq/dt, W/m^2, of one component's tendencies."""
    c = _constants()
    mass = budgets.pressure_thickness(state) / c["g"]
    total = 0.0
    if "air_temperature" in tendencies:
        total += np.sum(mass * c["Cp"] * np.asarray(
            tendencies["air_temperature"].to_units("degK/s").values,
            dtype=float).reshape(len(mass), -1)[:, 0])
    if "specific_humidity" in tendencies:
        total += np.sum(mass * c["Lv"] * np.asarray(
            tendencies["specific_humidity"].to_units("kg/kg/s").values,
            dtype=float).reshape(len(mass), -1)[:, 0])
    return total


def energy(days=5):
    """Which component adds energy the rest of the model does not count.

    Each tendency component is evaluated on the state before the step (a
    diagnostic call: nothing is applied), and each stepper's change to
    ``budgets.column_enthalpy`` is measured as it is applied. At equilibrium
    the sources should be zero except the boundary layer's, which is the
    surface's sensible and latent heat. What is left over is what the TOA
    residual is made of.
    """
    from sympl import AdamsBashforth

    tendencies, steppers, state, _ = load()
    timestep = climt.UnytTimeDelta(minutes=5)
    seconds = 300.0
    names = ([type(c).__name__ for c in tendencies]
             + [type(c).__name__ for c in steppers])
    added = {name: 0.0 for name in names}
    extra = dict(surface_turbulent=0.0, toa=0.0, moist_cp_adjustment=0.0)
    model = AdamsBashforth(*tendencies)
    n_steps = days * STEPS_PER_DAY
    c = _constants()
    Cvap = float(sympl.get_constant("heat_capacity_of_vapor_phase", "J/kg/K"))
    for _ in range(n_steps):
        for component in tendencies:
            try:
                result = component(state)
            except TypeError:
                result = component(state, timestep)
            added[type(component).__name__] += _enthalpy_rate(state, result[0])
        diagnostics, new = model(state, timestep)
        state.update(new)
        state.update(diagnostics)
        for component in steppers:
            name = type(component).__name__
            before = budgets.column_enthalpy(state)
            T0 = column(state, "air_temperature")
            q0 = column(state, "specific_humidity")
            diagnostics, new = component(state, timestep)
            state.update(new)
            state.update(diagnostics)
            added[name] += (budgets.column_enthalpy(state) - before) / seconds
            if name == "DryConvectiveAdjustment":
                # The enthalpy the scheme itself conserves: cp of moist air.
                T1 = column(state, "air_temperature")
                q1 = column(state, "specific_humidity")
                mass = budgets.pressure_thickness(state) / c["g"]
                cp0 = c["Cp"] * (1 - q0) + Cvap * q0
                cp1 = c["Cp"] * (1 - q1) + Cvap * q1
                extra["moist_cp_adjustment"] += np.sum(
                    mass * (cp1 * T1 - cp0 * T0)) / seconds
        state["time"] += timestep
        extra["surface_turbulent"] += (budgets.sensible_heat_flux(state)
                                       + budgets.latent_heat_flux(state))
        extra["toa"] += budgets.toa_imbalance(state)
    print(f"  column enthalpy (cp_dry T + Lv q) added, mean over {days} days, "
          "W/m^2:")
    for name in names:
        print(f"    {name:28s} {added[name] / n_steps:+9.3f}")
    print(f"    {'(surface SH + LH)':28s} "
          f"{extra['surface_turbulent'] / n_steps:+9.3f}")
    print(f"    {'(adjustment, moist cp)':28s} "
          f"{extra['moist_cp_adjustment'] / n_steps:+9.4f}")
    print(f"    {'(TOA imbalance)':28s} {extra['toa'] / n_steps:+9.3f}")


def heating(days=5):
    """Where each stepper heats, how saturated each level is, and what
    Emanuel's own diagnostics say."""
    tendencies, steppers, state, _ = load()
    timestep = climt.UnytTimeDelta(minutes=5)
    names = [type(c).__name__ for c in steppers]
    heat = {name: np.zeros(28) for name in names}
    condensing = np.zeros(28)
    humidity = np.zeros(28)
    cape, flux = [], []

    def after(s):
        humidity[:] += soundings.relative_humidity(s)
        cape.append(scalar(s, "atmosphere_convective_available_potential_energy"))
        flux.append(scalar(s, "cloud_base_mass_flux"))

    class Metered:
        def __init__(self, stepper):
            self.stepper = stepper

        def __call__(self, s, dt):
            T0 = column(s, "air_temperature")
            q0 = column(s, "specific_humidity")
            diagnostics, new = self.stepper(s, dt)
            name = type(self.stepper).__name__
            if "air_temperature" in new:     # SurfaceHumidity has none
                heat[name][:] += (np.asarray(new["air_temperature"].values,
                                             dtype=float).reshape(28, -1)[:, 0]
                                  - T0)
            if name == "GridScaleCondensation":
                condensing[:] += np.asarray(
                    new["specific_humidity"].values,
                    dtype=float).reshape(28, -1)[:, 0] < q0 - 1e-12
            return diagnostics, new

    n_steps = days * STEPS_PER_DAY
    stepping.integrate(tendencies, [Metered(s) for s in steppers], state,
                       timestep, n_steps, after_step=after)
    p = column(state, "air_pressure") / 100
    print(f"  mean over {days} days; heating in K/day")
    print(f"  {'p (hPa)':>8} {'RH %':>5} " + " ".join(f"{n[:12]:>12}"
                                                     for n in names)
          + "  condensing (fraction of steps)")
    for k in range(20):
        print(f"  {p[k]:8.1f} {100 * humidity[k] / n_steps:5.0f} "
              + " ".join(f"{heat[n][k] / days:12.3f}" for n in names)
              + f"  {condensing[k] / n_steps:.3f}")
    cape, flux = np.array(cape), np.array(flux)
    print(f"  Emanuel CAPE: mean {cape.mean():.2f}, max {cape.max():.2f} J/kg; "
          f"cloud-base mass flux > 0 on {np.mean(flux > 0):.0%} of steps")


def timestep(days=60):
    """The shipped state stepped at four timesteps; the last 30 days' means,
    and the first two days' (what the page's cell 5 sees)."""
    for minutes in (2.5, 5.0, 10.0, 20.0):
        record, p = run(minutes=minutes, days=days)
        report(f"dt = {minutes:g} min", record, p)
        first = record["days"] <= 2.0 + 1e-9
        print(f"    first 2 days: Ts {record['surface'][first].mean():.3f}  "
              f"P {record['precipitation'][first].mean():.3f} (Emanuel "
              f"{record['convective'][first].mean():.3f})  E "
              f"{record['evaporation'][first].mean():.3f}  mass flux "
              f"{record['mass_flux'][first].mean():.2e}")


def _dropped(name, days=300):
    record, p = run(days=days, drop=(name,))
    for day in (2, 5, 10, 30, 100, days):
        window = (record["days"] > day - 1) & (record["days"] <= day)
        T = record["air_temperature"][window].mean(axis=0)
        lapse = lapse_rate(T, p)
        print(f"    day {day:3d}: Ts {record['surface'][window].mean():.3f}  "
              f"P {record['precipitation'][window].mean():.3f} (Emanuel "
              f"{record['convective'][window].mean():.3f})  lapse 792-589 hPa "
              + " ".join(f"{x:.1f}" for x in lapse[8:13])
              + f"  below 500 hPa {lapse[lower_troposphere(p)].mean():.2f}")
    report(f"without {name}", record, p)


def nogsc():
    _dropped("GridScaleCondensation")


def nodca():
    _dropped("DryConvectiveAdjustment")


def feedback():
    """Why the moist column warms twice as much as the dry one: the OLR change
    between the two shipped states, taken apart one field at a time with the
    others held at the base state's values -- radiation calls only."""
    _, _, base, _ = load()
    _, _, doubled, _ = load(DOUBLED)
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")

    def olr(s):
        return float(longwave(s)[1]["upwelling_longwave_flux_in_air"]
                     .values[-1, 0, 0])

    def swap(*names):
        s = copy.deepcopy(base)
        for name in names:
            s[name].values[:] = doubled[name].values
        return s

    warming = (scalar(doubled, "surface_temperature")
               - scalar(base, "surface_temperature"))
    forcing = olr(base) - olr(swap("mole_fraction_of_carbon_dioxide_in_air"))
    temperature = olr(swap("air_temperature", "surface_temperature")) - olr(base)
    vapour = olr(swap("specific_humidity")) - olr(base)
    print(f"  file-to-file warming {warming:+.3f} K")
    print(f"  forcing (CO2 doubled, nothing else)   {forcing:+.2f} W/m^2")
    print(f"  OLR change, temperatures alone        {temperature:+.2f} W/m^2 "
          f"({temperature / warming:.2f} per K)")
    print(f"  OLR change, humidity alone            {vapour:+.2f} W/m^2 "
          f"({vapour / warming:.2f} per K)")
    print(f"  sum {temperature + vapour:+.2f} against the forcing "
          f"{forcing:+.2f}")
    print(f"  warming with the humidity held fixed: "
          f"{forcing / (temperature / warming):.2f} K")
    surface = copy.deepcopy(base)
    surface["surface_temperature"].values[:] += 1.0
    print(f"  surface +1 K, air fixed: OLR {olr(surface) - olr(base):+.2f} "
          "W/m^2")


def adiabat():
    """The moist adiabat against this column."""
    _, _, state, _ = load()
    T = column(state, "air_temperature")
    p = column(state, "air_pressure")
    rate = soundings.moist_adiabatic_lapse_rate(T, p)
    print("  moist adiabatic lapse rate on the shipped column's own T:")
    for k in range(18):
        print(f"    {p[k] / 100:7.1f} hPa  T {T[k]:.2f}  {rate[k]:.2f} K/km")
    levels = np.geomspace(1.0e5, 5.0e4, 400)
    c = _constants()
    for T0 in (280.0, 283.0, 284.0, 285.0, 286.0, 290.0, 300.0):
        Tp = soundings.moist_adiabat(T0, 1.0e5, levels)
        z = np.sum(c["Rd"] * 0.5 * (Tp[1:] + Tp[:-1]) / c["g"]
                   * np.log(levels[:-1] / levels[1:]))
        print(f"  saturated parcel from {T0:.0f} K at 1000 hPa: mean lapse "
              f"1000-500 hPa {1e3 * (Tp[0] - Tp[-1]) / z:.2f} K/km; "
              f"{Tp[-1]:.1f} K at 500 hPa")


def cost(n=20, warmup=3):
    if os.environ.get("NUMBA_DISABLE_JIT") != "1":
        print("  (JIT is on: this is not the browser stand-in)")
    tendencies, steppers, state, _ = load()
    timestep = climt.UnytTimeDelta(minutes=5)
    stepping.integrate(tendencies, steppers, state, timestep, warmup)
    start = time.perf_counter()
    stepping.integrate(tendencies, steppers, state, timestep, n)
    per_step = (time.perf_counter() - start) / n
    print(f"  page 12 stack: {1e3 * per_step:.1f} ms/step native, "
          f"~{3.5 * per_step:.2f} s/step browser; 576 steps "
          f"~{576 * 3.5 * per_step / 60:.1f} min, 288 steps "
          f"~{288 * 3.5 * per_step:.0f} s")
    start = time.perf_counter()
    load()
    load(DOUBLED)
    print(f"  loading both states: {time.perf_counter() - start:.2f} s native")


def main():
    sympl.set_backend(climt.UnytBackend())
    # The ImplicitTendencyComponent warning, once per AdamsBashforth: the page
    # explains it; here it is noise.
    warnings.filterwarnings("ignore", message="Using an ImplicitTendency")
    tasks = dict(month=month, energy=energy, heating=heating,
                 timestep=timestep, nogsc=nogsc, nodca=nodca,
                 feedback=feedback, adiabat=adiabat, cost=cost)
    names = sys.argv[1:] or ["all"]
    if names == ["all"]:
        names = [n for n in tasks if n != "cost"]
    for name in names:
        print(f"== {name}", flush=True)
        tasks[name]()


if __name__ == "__main__":
    main()
