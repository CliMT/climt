"""Generate the two radiative-convective equilibrium states pages 11 and 12 ship.

Pages 11 and 12 perturb an equilibrium rather than spending thousands of
in-browser steps finding one. This produces those states.

**The defaults here are exactly what shipped.** Running this with no arguments
must reproduce the committed files to within the convergence threshold. That is
not a nicety: `scripts/generate_tour_spectrum_table.py` shipped with defaults
that differed from what it had been run with, and reconstructing the real
invocation cost an afternoon.

    conda run -n climt python scripts/generate_tour_equilibria.py
    conda run -n climt python scripts/generate_tour_equilibria.py --moist
    conda run -n climt python scripts/generate_tour_equilibria.py --out /tmp

The configuration constants below are Task 0's measured decisions (see the plan
`docs/superpowers/plans/2026-08-20-modelling-tour-rce.md`, Task 8's log): the
dry stack converges at dt = 12 h, the moist Emanuel stack at dt = 5 min, and
|TOA imbalance| < 0.5 W/m^2 is the accepted equilibrium gate for both -- the
spec author has ruled the residual upper-column spread benign, so the states
ship directly from this convergence.

Regenerate deliberately, when `tests/test_modelling_tour.py`'s residual test
says the physics moved -- never on a schedule, and never from a dependency
hash. See that test for why.
"""
import argparse
import os
import sys

import numpy as np
import sympl

import climt

sys.path.insert(0, os.path.join("docs", "modelling-tour", "_tour"))
import budgets      # noqa: E402
import states       # noqa: E402
import stepping     # noqa: E402

DATA_DIR = os.path.join("docs", "modelling-tour", "_data")

# ---- the configuration the shipped states were made with. Changing any of
# ---- these invalidates both files and every number on pages 11 and 12.
NZ = 28
SOLAR = 240.0                 # prescribed absorbed shortwave at the surface
SLAB_DEPTH_M = 2.0
CO2_PPM = 330.0
TABLE = "earth_low_res_lw"
DRY_DT_HOURS = 12.0           # Task 0: dry stack converges at 12 h
MOIST_DT_HOURS = 5.0 / 60.0   # Task 0: the moist Emanuel stack's stable dt, 5 min
CONVERGENCE_W_M2 = 0.5        # |TOA imbalance| accepted as equilibrium (Task 0)
# Convergence is TOA balance AND a surface that has stopped moving. The
# surface *flux* imbalance is not a usable gate: the one-step flux lag
# `_tour/stepping.py` documents makes it oscillate several W/m^2 step to step
# even at equilibrium, so it never sits under any small threshold. The surface
# *temperature*, by contrast, settles into a tiny limit cycle. So the second
# gate is that the mean surface temperature over one window equals the mean
# over the previous one (means cancel the oscillation whatever its period),
# which is a strict form of the drift the residual test checks.
SURFACE_STEADY_K = 0.02       # window-to-window mean surface-temperature drift
STEADY_WINDOW_STEPS = 1000    # the averaging window for that drift
MAX_STEPS = 400000            # ~1400 sim-days for the moist run; a cap, not a target
SURFACE_SPECIFIC_HUMIDITY = 0.015    # saturated surface, page 9's default
WIND_M_S = 5.0                # the wind the relaxation holds the column at
WIND_TIMESCALE_HOURS = 24.0
ROUGHNESS_LENGTH_M = 1e-3


# The wind relaxation is NOT in these lists. It is built per-state by
# `stepping.wind_relaxation`, because it has to write its target fields into
# the state it will act on -- so it cannot exist before the state does. It is
# appended to the tendency list in `build_state`, and `split` never sees it.
# The residual test rebuilds it the same way; see `components_for`.


def dry_components():
    """Page 11's physics components. The single definition; the test imports it."""
    return [
        climt.CorkLongwaveRadiation(optics="correlated_k", table=TABLE),
        climt.SlabSurface(),
        climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                  roughness_length=ROUGHNESS_LENGTH_M),
        climt.DryConvectiveAdjustment(),
    ]


def moist_components():
    """Page 12's physics components. The single definition; the test imports it."""
    return [
        climt.CorkLongwaveRadiation(optics="correlated_k", table=TABLE),
        climt.SlabSurface(),
        climt.EmanuelConvectionPython(),
        climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                  roughness_length=ROUGHNESS_LENGTH_M),
        climt.DryConvectiveAdjustment(),
        climt.GridScaleCondensation(),
    ]


def split(components):
    """(tendency components, stepper components), in the order given.

    A sympl Stepper declares ``output_properties``; a TendencyComponent and an
    ImplicitTendencyComponent do not. ``EmanuelConvectionPython`` is the latter
    and lands in the tendency list, which is where climt has always run it --
    see _tour/stepping.py and page 12 for the warning that provokes.
    """
    tendencies = [c for c in components if not hasattr(c, "output_properties")]
    steppers = [c for c in components if hasattr(c, "output_properties")]
    return tendencies, steppers


def build_state(components, moist):
    """The state, and the wind relaxation that has to be built alongside it.

    Returns ``(state, relaxation)``. The relaxation is a tendency component,
    but it cannot be constructed before the state exists, because it writes
    its equilibrium and timescale fields into that state.
    """
    state = climt.get_default_state(
        components, grid_state=climt.get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = SLAB_DEPTH_M
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_PPM * 1e-6
    if moist:
        state["surface_specific_humidity"].values[:] = \
            SURFACE_SPECIFIC_HUMIDITY
    else:
        # Page 11's column is genuinely dry: no water vapour anywhere, so its
        # 14-band radiation sees CO2 alone. The page says so, and it is why
        # page 11 equilibrates well below Earth's surface temperature.
        state["specific_humidity"].values[:] = 0.0
        state["surface_specific_humidity"].values[:] = 0.0
    relaxation = stepping.wind_relaxation(
        state, WIND_M_S, timescale_hours=WIND_TIMESCALE_HOURS)
    return state, relaxation


def run_to_equilibrium(components, moist, dt_hours, check_every=50):
    tendencies, steppers = split(components)
    state, relaxation = build_state(components, moist)
    tendencies = tendencies + [relaxation]
    timestep = climt.UnytTimeDelta(hours=dt_hours)

    window = max(1, STEADY_WINDOW_STEPS // check_every)   # checks per window
    surface_history = []
    for step in range(0, MAX_STEPS, check_every):
        stepping.integrate(tendencies, steppers, state, timestep, check_every)
        toa = budgets.toa_imbalance(state)
        surface_imbalance = budgets.surface_imbalance(state)
        surface = float(state["surface_temperature"].values.ravel()[0])
        surface_history.append(surface)

        # Trend in the surface temperature, oscillation removed: the mean over
        # the last window minus the mean over the one before it. Undefined
        # until two full windows have been seen.
        if len(surface_history) >= 2 * window:
            recent = np.mean(surface_history[-window:])
            earlier = np.mean(surface_history[-2 * window:-window])
            drift = abs(recent - earlier)
        else:
            drift = float("inf")

        print(f"  step {step + check_every:6d}  "
              f"({(step + check_every) * dt_hours / 24.0:7.1f} d)  "
              f"TOA {toa:+8.4f}  surf(flux) {surface_imbalance:+8.4f} W/m^2  "
              f"Tsurf {surface:8.4f} K  drift {drift:.4f} K")
        if abs(toa) < CONVERGENCE_W_M2 and drift < SURFACE_STEADY_K:
            return state, step + check_every
    raise SystemExit(
        f"did not converge in {MAX_STEPS} steps at dt = {dt_hours} h: need "
        f"|TOA| < {CONVERGENCE_W_M2} W/m^2 and a surface drift < "
        f"{SURFACE_STEADY_K} K; last TOA imbalance "
        f"{budgets.toa_imbalance(state):+.4f}")


def generate(kind, out_dir):
    moist = (kind == "moist")
    components = moist_components() if moist else dry_components()
    dt_hours = MOIST_DT_HOURS if moist else DRY_DT_HOURS
    print(f"{kind} equilibrium: dt = {dt_hours} h, nz = {NZ}, "
          f"slab = {SLAB_DEPTH_M} m, CO2 = {CO2_PPM} ppm")

    state, n_steps = run_to_equilibrium(components, moist, dt_hours)
    provenance = dict(
        table=TABLE, nz=NZ, dt_hours=dt_hours, n_steps=n_steps,
        slab_depth_m=SLAB_DEPTH_M, solar=SOLAR, co2_ppm=CO2_PPM,
        wind_m_s=WIND_M_S, wind_timescale_hours=WIND_TIMESCALE_HOURS,
        roughness_length_m=ROUGHNESS_LENGTH_M,
        components=[type(c).__name__ for c in components] + ["UnytRelaxation"],
        toa_imbalance=budgets.toa_imbalance(state),
        surface_imbalance=budgets.surface_imbalance(state),
        convergence_threshold_w_m2=CONVERGENCE_W_M2,
        surface_steady_k=SURFACE_STEADY_K,
        steady_window_steps=STEADY_WINDOW_STEPS,
    )
    path = os.path.join(out_dir, f"rce_{kind}_equilibrium.npz")
    states.save(path, state, provenance)
    print(states.describe({**provenance,
                           "climt_version": climt.__version__,
                           "saved_at": "(just now)"}))
    print(f"wrote {path} ({os.path.getsize(path) / 1024:.1f} kB)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry", action="store_true",
                        help="generate only the dry equilibrium")
    parser.add_argument("--moist", action="store_true",
                        help="generate only the moist equilibrium")
    parser.add_argument("--out", default=DATA_DIR,
                        help=f"output directory (default: {DATA_DIR})")
    args = parser.parse_args()

    sympl.set_backend(climt.UnytBackend())
    os.makedirs(args.out, exist_ok=True)
    wanted = [k for k, on in (("dry", args.dry), ("moist", args.moist)) if on]
    for kind in wanted or ["dry", "moist"]:
        generate(kind, args.out)


if __name__ == "__main__":
    main()
