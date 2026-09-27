"""Generate the radiative-convective equilibrium states pages 11 and 12 ship.

Pages 11 and 12 perturb an equilibrium rather than spending thousands of
in-browser steps finding one. This produces those states: the dry and the
moist equilibrium, and the moist one re-equilibrated at doubled CO2. Page 11's
2xCO2 runs live (1000 steps at 12 h, ~3 min in the browser); page 12's
cannot, because at the moist column's 5 min timestep it takes tens of
thousands of steps, so it ships too.

**The defaults here are exactly what shipped.** Running this with no arguments
must reproduce the committed files to within the convergence threshold. That is
not a nicety: `scripts/generate_tour_spectrum_table.py` shipped with defaults
that differed from what it had been run with, and reconstructing the real
invocation cost an afternoon.

    conda run -n climt python scripts/generate_tour_equilibria.py
    conda run -n climt python scripts/generate_tour_equilibria.py --moist
    conda run -n climt python scripts/generate_tour_equilibria.py --moist-2xco2
    conda run -n climt python scripts/generate_tour_equilibria.py --out /tmp

The configuration constants below are Task 0's measured decisions (see the plan
`docs/superpowers/plans/2026-08-20-modelling-tour-rce.md`, Task 8's log): the
dry stack converges at dt = 12 h and the moist Emanuel stack at dt = 5 min.
The dry gate is |TOA imbalance| < 0.5 W/m^2 with a stationary surface; the
moist column, whose surface is held saturated at its own temperature by page
9's SurfaceHumidity, has a trend gate of its own (MOIST_GATE, below, says
why). The spec author has ruled the residual upper-column spread benign, so
the states ship directly from this convergence.

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

# Absolute, so the generator imports cleanly from any working directory: the
# residual test and the experiment scripts load it with importlib, and
# stepping.SurfaceHumidity imports `soundings` from this directory at call time.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "docs", "modelling-tour", "_tour"))
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
# The moist column's gate is different, for two measured reasons (2026-09-27,
# saturated surface; the plan's Task 9 log has the numbers):
#   * it is noisy. Episodic convection swings the surface flux by +-40 W/m^2
#     and the instantaneous TOA by ~0.5 W/m^2 (1 sd), while its late drift is
#     slow, ~0.05 K per 30 days. At 5 min a 1000-step window is 3.5 days, and
#     an instantaneous |TOA| < 0.5 is passed by noise: the dry gate stopped the
#     spin-up at 194 400 steps with the 30-day mean TOA still -1.7 and falling.
#   * it never reaches TOA = 0. EmanuelConvectionPython is not fully
#     energy-conserving, so the column settles with a steady TOA imbalance --
#     about -1.1 W/m^2 over a saturated surface (the stratosphere is in
#     radiative balance and the slab has stopped moving, so this is not
#     storage). No |TOA| threshold below that is a convergence test.
# So the moist gate is on *trends* over 60-day (17 280-step) windows: the mean
# surface temperature and the mean TOA must each match the previous window's,
# to 0.01 K and 0.1 W/m^2, and the mean TOA must be inside the +-2 W/m^2 the
# scheme's residual lives in (which rules out the turning point early in the
# spin-up, where the surface peaks while TOA is still -25). Replayed on the
# measured trajectory it stops at ~296 000 steps, 0.01 K from where the column
# settles; the 30-day windows tried first stopped 0.1 K short.
DRY_GATE = dict(window_steps=STEADY_WINDOW_STEPS, surface_steady_k=SURFACE_STEADY_K,
                toa_limit_w_m2=CONVERGENCE_W_M2, toa_steady_w_m2=None,
                mean_toa=False)
MOIST_GATE = dict(window_steps=17280, surface_steady_k=0.01,
                  toa_limit_w_m2=2.0, toa_steady_w_m2=0.1, mean_toa=True)
MAX_STEPS = 400000            # ~1400 sim-days for the moist run; a cap, not a target
# The moist surface is saturated *at its own temperature*, every step, by
# page 9's stepping.SurfaceHumidity. Until 2026-09-27 this was a fixed
# surface_specific_humidity of 0.015 kg/kg written once, which at the ~280 K
# the column settled at was ~246 % of saturation -- the unphysical surface
# page 9 teaches against. See the plan's Task 9 log.
SURFACE_RELATIVE_HUMIDITY = 1.0
WIND_M_S = 5.0                # the wind the relaxation holds the column at
WIND_TIMESCALE_HOURS = 24.0
ROUGHNESS_LENGTH_M = 1e-3
CO2_DOUBLING = 2.0            # the perturbation rce_moist_2xco2_equilibrium.npz ships


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
    """Page 12's physics components. The single definition; the test imports it.

    ``SurfaceHumidity`` is the first *stepper* (``split`` keeps list order, and
    the three entries before it are tendency components), so the boundary
    layer exchanges water with the surface as the slab has just left it --
    the order page 9 wires and explains.
    """
    return [
        climt.CorkLongwaveRadiation(optics="correlated_k", table=TABLE),
        climt.SlabSurface(),
        climt.EmanuelConvectionPython(),
        stepping.SurfaceHumidity(SURFACE_RELATIVE_HUMIDITY),
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
    if not moist:
        # Page 11's column is genuinely dry: no water vapour anywhere, so its
        # 14-band radiation sees CO2 alone. The page says so, and it is why
        # page 11 equilibrates well below Earth's surface temperature.
        state["specific_humidity"].values[:] = 0.0
        state["surface_specific_humidity"].values[:] = 0.0
    relaxation = stepping.wind_relaxation(
        state, WIND_M_S, timescale_hours=WIND_TIMESCALE_HOURS)
    return state, relaxation


def run_to_equilibrium(tendencies, steppers, state, dt_hours,
                       check_every=50, print_every=50, window_steps=None,
                       surface_steady_k=None, toa_limit_w_m2=None,
                       toa_steady_w_m2=None, mean_toa=False, summary=None):
    """Step ``state`` in place until it passes the gate; return the steps taken.

    ``tendencies`` must already include the wind relaxation. The starting
    state is the caller's: a fresh ``build_state`` for the spin-ups, a loaded
    equilibrium for the 2xCO2 re-equilibration, or a deliberately different
    profile for a start-independence check.

    The gate, over consecutive windows of ``window_steps``: the mean surface
    temperature changes by less than ``surface_steady_k`` from one window to
    the next; |TOA| < ``toa_limit_w_m2``, instantaneous or (``mean_toa``) the
    window mean; and, if ``toa_steady_w_m2`` is set, the window-mean TOA also
    changes by less than that. The keyword defaults are the dry gate; pass
    ``**gate_for(kind)`` -- ``MOIST_GATE`` explains why the moist one differs.

    If ``summary`` is a dict, the final window's mean TOA and surface
    temperature are written into it as ``window_mean_toa_w_m2`` and
    ``window_mean_surface_temperature_k``.
    """
    window_steps = window_steps or STEADY_WINDOW_STEPS
    surface_steady_k = surface_steady_k or SURFACE_STEADY_K
    toa_limit_w_m2 = toa_limit_w_m2 or CONVERGENCE_W_M2
    timestep = climt.UnytTimeDelta(hours=dt_hours)

    window = max(1, window_steps // check_every)   # checks per window
    surface_history = []
    toa_history = []
    for step in range(0, MAX_STEPS, check_every):
        stepping.integrate(tendencies, steppers, state, timestep, check_every)
        toa = budgets.toa_imbalance(state)
        surface_imbalance = budgets.surface_imbalance(state)
        surface = float(state["surface_temperature"].values.ravel()[0])
        surface_history.append(surface)
        toa_history.append(toa)

        # Trends, oscillation removed: the mean over the last window minus the
        # mean over the one before it. Undefined until two full windows.
        if len(surface_history) >= 2 * window:
            drift = abs(np.mean(surface_history[-window:])
                        - np.mean(surface_history[-2 * window:-window]))
            toa_trend = abs(np.mean(toa_history[-window:])
                            - np.mean(toa_history[-2 * window:-window]))
        else:
            drift = toa_trend = float("inf")
        gated_toa = np.mean(toa_history[-window:]) if mean_toa else toa

        converged = (abs(gated_toa) < toa_limit_w_m2
                     and drift < surface_steady_k
                     and (toa_steady_w_m2 is None
                          or toa_trend < toa_steady_w_m2))
        if converged or (step + check_every) % print_every == 0:
            mean_note = (f"  TOA(mean) {gated_toa:+8.4f} (trend "
                         f"{toa_trend:.3f})" if mean_toa else "")
            print(f"  step {step + check_every:6d}  "
                  f"({(step + check_every) * dt_hours / 24.0:7.1f} d)  "
                  f"TOA {toa:+8.4f}{mean_note}  surf(flux) "
                  f"{surface_imbalance:+8.4f} "
                  f"W/m^2  Tsurf {surface:8.4f} K  drift {drift:.4f} K",
                  flush=True)
        if converged:
            if summary is not None:
                summary.update(
                    window_mean_toa_w_m2=float(np.mean(toa_history[-window:])),
                    window_mean_surface_temperature_k=float(
                        np.mean(surface_history[-window:])))
            return step + check_every
    raise SystemExit(
        f"did not converge in {MAX_STEPS} steps at dt = {dt_hours} h: need "
        f"|TOA| < {toa_limit_w_m2} W/m^2 and a surface drift < "
        f"{surface_steady_k} K; last TOA imbalance "
        f"{budgets.toa_imbalance(state):+.4f}")


def gate_for(kind):
    """The ``run_to_equilibrium`` keyword arguments for ``kind``'s gate."""
    return dict(MOIST_GATE) if kind == "moist" else dict(DRY_GATE)


def _provenance(components, dt_hours, n_steps, state, co2_ppm, gate,
                summary):
    provenance = dict(summary,
        table=TABLE, nz=NZ, dt_hours=dt_hours, n_steps=n_steps,
        slab_depth_m=SLAB_DEPTH_M, solar=SOLAR, co2_ppm=co2_ppm,
        wind_m_s=WIND_M_S, wind_timescale_hours=WIND_TIMESCALE_HOURS,
        roughness_length_m=ROUGHNESS_LENGTH_M,
        components=[type(c).__name__ for c in components] + ["UnytRelaxation"],
        toa_imbalance=budgets.toa_imbalance(state),
        surface_imbalance=budgets.surface_imbalance(state),
        convergence_threshold_w_m2=gate["toa_limit_w_m2"],
        steady_window_steps=gate["window_steps"],
        surface_steady_k=gate["surface_steady_k"],
        toa_gate=("window mean" if gate["mean_toa"] else "instantaneous"),
    )
    if gate["toa_steady_w_m2"] is not None:
        provenance["toa_steady_w_m2"] = gate["toa_steady_w_m2"]
    # A moist state records the surface relative humidity it was held at, so
    # a reader of the file can tell a SurfaceHumidity state from the old
    # fixed-surface-q ones (which lack this key).
    for component in components:
        if isinstance(component, stepping.SurfaceHumidity):
            provenance["surface_relative_humidity"] = \
                component.relative_humidity
    return provenance


def _write(path, state, provenance):
    states.save(path, state, provenance)
    print(states.describe({**provenance,
                           "climt_version": climt.__version__,
                           "saved_at": "(just now)"}))
    print(f"wrote {path} ({os.path.getsize(path) / 1024:.1f} kB)")


def generate(kind, out_dir):
    moist = (kind == "moist")
    components = moist_components() if moist else dry_components()
    dt_hours = MOIST_DT_HOURS if moist else DRY_DT_HOURS
    print(f"{kind} equilibrium: dt = {dt_hours} h, nz = {NZ}, "
          f"slab = {SLAB_DEPTH_M} m, CO2 = {CO2_PPM} ppm")

    tendencies, steppers = split(components)
    state, relaxation = build_state(components, moist)
    summary = {}
    n_steps = run_to_equilibrium(tendencies + [relaxation], steppers, state,
                                 dt_hours, print_every=1000 if moist else 50,
                                 summary=summary, **gate_for(kind))
    _write(os.path.join(out_dir, f"rce_{kind}_equilibrium.npz"), state,
           _provenance(components, dt_hours, n_steps, state, CO2_PPM,
                       gate_for(kind), summary))


def load_equilibrium(kind, data_dir=DATA_DIR, filename=None):
    """A shipped equilibrium, ready to step: ``(tendencies, steppers, state,
    provenance)``, with the wind relaxation rebuilt on the loaded state as the
    pages and the residual test do it. ``filename`` defaults to
    ``rce_<kind>_equilibrium.npz``; pass another to load, say, the 2xCO2
    state with ``kind``'s components."""
    components = moist_components() if kind == "moist" else dry_components()
    tendencies, steppers = split(components)
    filename = filename or f"rce_{kind}_equilibrium.npz"
    state, provenance = states.load(
        os.path.join(data_dir, filename), components,
        grid_state=climt.get_grid(nx=1, ny=1, nz=NZ))
    # initialise=False: keep the spun-up, drag-sheared wind profile.
    tendencies = tendencies + [stepping.wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]
    return tendencies, steppers, state, provenance


def generate_doubled(kind, out_dir):
    """Page 12's 2xCO2 experiment, run offline: double CO2 on the shipped
    equilibrium and re-converge under the same gate.

    It ships rather than running in the browser because at dt = 5 min it takes
    tens of thousands of steps, and page 12 costs ~0.19 s a step there. It
    starts from the file in ``out_dir``, so regenerate the base state first
    if that has changed.
    """
    components = moist_components() if kind == "moist" else dry_components()
    dt_hours = MOIST_DT_HOURS if kind == "moist" else DRY_DT_HOURS
    tendencies, steppers, state, base = load_equilibrium(kind, out_dir)
    before = float(state["surface_temperature"].values.ravel()[0])
    co2_ppm = CO2_DOUBLING * base["co2_ppm"]
    print(f"{kind} 2xCO2: from rce_{kind}_equilibrium.npz "
          f"(Tsurf {before:.3f} K), CO2 {base['co2_ppm']} -> {co2_ppm} ppm, "
          f"dt = {dt_hours} h")

    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2_ppm * 1e-6
    summary = {}
    n_steps = run_to_equilibrium(tendencies, steppers, state, dt_hours,
                                 print_every=STEADY_WINDOW_STEPS,
                                 summary=summary, **gate_for(kind))
    after = float(state["surface_temperature"].values.ravel()[0])
    provenance = _provenance(components, dt_hours, n_steps, state, co2_ppm,
                             gate_for(kind), summary)
    provenance.update(
        perturbed_from=f"rce_{kind}_equilibrium.npz",
        perturbed_from_saved_at=base["saved_at"],
        base_surface_temperature_k=before,
        surface_temperature_k=after,
        warming_k=after - before,
    )
    # The moist gate works on window means, and the column's surface
    # temperature wanders ~0.1 K step to step, so the warming worth quoting
    # is between the two final window means, not the two final instants.
    if ("window_mean_surface_temperature_k" in base
            and "window_mean_surface_temperature_k" in summary):
        provenance["window_mean_warming_k"] = (
            summary["window_mean_surface_temperature_k"]
            - base["window_mean_surface_temperature_k"])
    print(f"  warming {after - before:+.3f} K in {n_steps} steps "
          f"({n_steps * dt_hours / 24.0:.1f} d)")
    _write(os.path.join(out_dir, f"rce_{kind}_2xco2_equilibrium.npz"), state,
           provenance)


TARGETS = ("dry", "moist", "moist-2xco2")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry", action="store_true",
                        help="generate only the dry equilibrium")
    parser.add_argument("--moist", action="store_true",
                        help="generate only the moist equilibrium")
    parser.add_argument("--moist-2xco2", action="store_true",
                        help="generate only page 12's 2xCO2 re-equilibration, "
                             "from the moist equilibrium already in --out")
    parser.add_argument("--out", default=DATA_DIR,
                        help=f"output directory (default: {DATA_DIR})")
    args = parser.parse_args()

    sympl.set_backend(climt.UnytBackend())
    os.makedirs(args.out, exist_ok=True)
    chosen = dict(dry=args.dry, moist=args.moist,
                  **{"moist-2xco2": args.moist_2xco2})
    # In this order: the 2xCO2 run starts from the moist file just written.
    for target in [t for t in TARGETS if chosen[t]] or TARGETS:
        if target == "moist-2xco2":
            generate_doubled("moist", args.out)
        else:
            generate(target, args.out)


if __name__ == "__main__":
    main()
