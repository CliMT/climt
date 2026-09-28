# Handoff: make page 12's moist column warm enough that Emanuel fires

**Branch:** `claude/sweet-fermi-fmaktk` (1:1 with origin at the time of writing).
**Plan:** `docs/superpowers/plans/2026-08-20-modelling-tour-rce.md` (Task 15 = page 12; Task 8 log §7 and Task 9 log have the measurement history).
**Env:** `CLIMT_PURE_PYTHON=1 pip install -e . numba pytest matplotlib flake8` into a venv works (conda not needed). Always run with `NUMBA_NUM_THREADS=1`: a single column gains nothing from numba threads, and several jobs at once on 4 cores ran ~10x slower without it.

## The decision (user, 2026-09-28)

Page 12's reveal is that latent heating relaxes the lapse rate toward the moist adiabat. **It doesn't happen in the shipped column**: 9.48 K/km over 792–500 hPa against page 11's dry 9.75 K/km. Asked how to handle it, the user said: **"make the column warmer so that emanuel fires. That's the whole point."** So the job is a stack or forcing change, then regenerating the moist states and re-measuring page 12. The spec is not to be amended, and page 12 is not to ship as measured.

## Why Emanuel doesn't fire now (from `docs/modelling-tour/12-moist-rce.qmd`, "Why this column is not on its moist adiabat")

- The surface is saturated (`stepping.SurfaceHumidity(1.0)`), but relative humidity in the dry-adjusted layer is 32 % at its bottom and 83 % at its top. Surface temperature is 285.99 K.
- Emanuel's closure looks only at cloud base. Its parcel crosses the stable layer near the ground and arrives up to 0.89 K colder than the air around it, right at the scheme's `DTMAX` tolerance of 0.9 K. The closure's cloud-base buoyancy averages −0.27 K, so the mass flux stays about 1e-3 kg m⁻² s⁻¹. There is plenty of instability above that: the parcel is buoyant from 836 to 180 hPa, with about 1700 J/kg of CAPE.
- Emanuel supplies 0.03 of the 2.65 mm/day; `GridScaleCondensation` does the rest, in place.
- Without GSC, Emanuel carries all the rain and the lapse rate *does* go moist: 7.2–8.5 K/km over 792–589 hPa. But the air at 792 hPa then reaches 157 % relative humidity, so GSC stays in the stack.
- Emanuel's rain roughly halves each time the timestep doubles (the closure's push is per call). If Emanuel starts doing the raining, **re-run the timestep check** (`tour_page12_measurements.py timestep`), because that stops being harmless.

## Where I was about to start

1. **Lever:** the prescribed absorbed shortwave, `SOLAR = 240.0` in `scripts/generate_tour_equilibria.py`. It is shared by the dry and moist columns and by earlier pages, and goes into `downwelling_shortwave_flux_in_air[0]` in `build_state`. CO₂ is too weak a lever (2×CO₂ gives +2.24 K). My guess is that the column needs to be roughly 10–15 K warmer, near a tropical ~300 K, before the boundary layer is moist enough for positive cloud-base buoyancy. That is not measured.
2. **Screen before spinning up:** start from the shipped `rce_moist_equilibrium.npz` and run about 120 days (≈35 000 steps at 5 min, a few minutes each) at `SOLAR` values of 260, 280, 300 and 320. In parallel, one job per core. Over the last 30 days, record Emanuel's share of the rain (`convective_precipitation_rate` against the total), the cloud-base mass flux, the surface temperature, and the lapse rate over 792–500 hPa. `scripts/experiments/tour_page12_measurements.py` already has helpers for this: `recorder`, `report`, `run`, `closure`, `lapse_rate`, `lower_troposphere`. Reuse them.
3. **Decide the scope with the user before regenerating.** Should the higher `SOLAR` apply to the moist column only (a new `MOIST_SOLAR`), or to everything? Moist-only keeps pages 7–11 untouched. But it breaks page 12's direct surface comparison with page 11 ("surface … against page 11's"), and the water-vapour-greenhouse framing in `_data/README.md`.
4. **Then regenerate** with `python scripts/generate_tour_equilibria.py --moist` followed by `--moist-2xco2`. They use `MOIST_GATE`, a 60-day trend gate; the base took about 296 000 steps. After that, re-run every page-12 measurement (`tour_page12_measurements.py`: month, closure, heating, energy, timestep, nogsc, nodca, feedback, adiabat, cost) and `tour_rce_shipped_remeasure.py settle-moist`. Rewrite page 12's numbers and the "why not on its moist adiabat" section, update the plan's Task 15 log, Task 9 log and `_data/README.md`, and run `pytest tests/test_modelling_tour.py`. The residual tests compare against each file's recorded settled TOA residual.

## Facts to keep straight

- The moist column never reaches TOA = 0. It settles at about −1.1 W/m² over the saturated surface. The generator's comments attribute most of that to `DryConvectiveAdjustment` using a moist heat capacity (+0.97 W/m²) and a little to Emanuel (+0.03). That's why the moist gate is on trends, not on the size of the TOA imbalance.
- Current shipped numbers, which will all change: moist base 285.99 K; 2×CO₂ 288.25 K; settled warming +2.24 K. Dry: 266.48 K, 2×CO₂ +1.15 to 1.20 K (page 11 quotes +1.15), run live for 1000 steps.
- `tests/test_modelling_tour.py::test_draw_evolution_builds_the_four_panel_figure` fails in a fresh venv with matplotlib 3.11 (no "non-interactive" warning). The failure predates this work.
- The container restart on 2026-09-27 stopped five waiters ("measurement job", "timestep/restart/spindown results"). Whatever they were waiting on was never re-run. The page 12 commits (`36c8062`, `590c682`) are in the branch.

## State at context clear (2026-09-28, ~03:10 UTC)

Work had already moved past "raise SOLAR". Commit `fd5dae0` added an opt-in `SimpleBoundaryLayer(diffuse='dry_static_energy')`. It mixes Cp T + g z instead of T, which removes the ~5.5 K/km stable layer near the ground that kept Emanuel's cloud-base buoyancy negative.

**Uncommitted work in progress** (16 files): the SBL component, its tests and caches, `HISTORY.rst`, the SBL user guide, pages 07–09, `rce_dry_equilibrium.npz`, the generator, and the page 7/8/9 measurement scripts. Don't commit it blind. Check each piece against the job logs below first.

**Detached jobs that were running** (reparented to init, so they survive a context clear but not a container restart):
- `generate_tour_equilibria.py --moist`, logging to `/tmp/claude-0/-home-user-climt/7f238135-7f48-4305-85c5-07b5ccfdb1c7/scratchpad/gen_moist.log`
- `tour_page8_measurements.py all`, logging to `.../7f238135-.../scratchpad/p8.log` in the same directory
- `tour_page11_measurements.py all` (its log is somewhere in that same scratchpad)

**First thing to do:** check whether each is still alive with `pgrep -af "generate_tour|tour_page"`. If one is gone, read the end of its log. If the container restarted, re-run it with `NUMBA_NUM_THREADS=1`.

## Decision (user, 2026-09-28): SOLAR = 240 everywhere, with the DSE boundary layer

The DSE boundary layer on its own makes Emanuel fire. At SOLAR = 240, Emanuel carries all 2.84 mm/day of the rain, and the 792–500 hPa lapse rate is 8.82 K/km against the dry column's 9.76. With the same diffusion, page 9 settles at 286.2 K moist against 266.7 K dry.

A previous session had raised SOLAR to 265 on every page. **Revert that to 240 on every page.** At 265, page 9's moist column (no convection scheme) runs away: 314 K after 1000 days and still rising. Any uncommitted SOLAR = 265 edits in the working tree must go back to 240.

The detached jobs listed above had already exited when this section was written. `pgrep` found none, so read their logs and re-run at 240.
