# A Modelling Tour of the Climate System — Radiative-Convective Equilibrium Tranche

**Date:** 2026-08-20
**Status:** design, awaiting review
**Scope:** tranche 2 of the Modelling Tour. Course chapters 10–12: turbulent heat
exchange, dry convection, radiative-convective equilibrium.

## Summary

Six new pages, `07`–`12`, in `docs/modelling-tour/`, continuing the section shipped
in PR #226. Where tranche 1 made single radiation calls on profiles the reader
prescribed, this tranche **integrates in time**: it introduces `sympl` time
stepping, a slab surface, a turbulent boundary layer, dry convective adjustment and
moist convection, and arrives at radiative-convective equilibrium.

It also **absorbs and deletes** `docs/radiative-transfer/09-live-rce.qmd`, whose
gray/non-grey stepping demo is what page 07 teaches properly.

## The thesis

Tranche 1 *prescribed* Γ = 6.5 K/km and a 200 K isothermal stratosphere on every
page. Every number those six pages quote rests on a profile the reader typed in.

Tranche 2 makes the model **produce** that profile. Page 11 shows that pure
radiation plus dry convection gives ~9.8 K/km — not the 6.5 that tranche 1 assumed.
Page 12 shows that latent heating relaxes it back toward ≈6.5, which is where the
assumption came from all along. The tranche closes a loop opened on page 01.

## Decisions

| Question | Decision |
|---|---|
| How far past the notes | Chapters 10–12 only, but with an interactive surface. No clouds, no sea ice, no polar capstone. |
| Surface | `SlabSurface` at **2 m** mixed-layer depth by default, dropping to 1 m if Task 0 shows the perturbation runs are too slow. Thin on purpose: response time is linear in depth, and it is the difference between a usable page and an unusable one. |
| Shortwave | **Prescribed absorbed flux at the surface** (`SOLAR = 240` W/m², as the live-RCE page does) on every page. `CorkShortwaveRadiation` is *deliberately deferred* — see Deferred below. |
| Moisture | Dry first, then moist, twice over: dry BL → moist BL (virtual potential temperature, buoyancy) on chapter 10; dry RCE → moist RCE on chapter 12. |
| Condensation | `GridScaleCondensation` is included unconditionally from page 09 onward. |
| Moist convection | `EmanuelConvectionPython` (pure Python; the Fortran `EmanuelConvection` is absent in the browser). |
| Reaching equilibrium | **Ship precomputed equilibrium states** for pages 11 and 12; those pages load one and run short perturbation experiments. Page 07 integrates live. |
| Existing live-RCE page | **Deleted outright**, no stub. |
| Backend | `UnytBackend`, as tranche 1. |
| Page numbers | `07`–`12`, continuing the tour. Provisional until the arc is locked in implementation. |

## Audience and objectives

Same audience as tranche 1: EC2213 students, plus readers arriving cold.

**Science objectives.** By the end, a student can explain and *compute*: how a
column approaches radiative equilibrium and how long it takes; why radiative
equilibrium has a surface–air temperature discontinuity and what erases it; how
bulk turbulent fluxes are computed and what sets their partition between sensible
and latent; why moisture changes buoyancy through the virtual temperature; why an
unstable column mixes to a dry adiabat while conserving enthalpy; what lapse rate
radiative-convective equilibrium actually produces, dry and moist; and what a
climate sensitivity is when you have measured one yourself.

**Craft objectives.** A student can distinguish `TendencyComponent`,
`ImplicitTendencyComponent` and `Stepper` and knows why the distinction changes the
time loop; can build and run an `AdamsBashforth` integration with the correct
timestep type; can compose several components in a defensible order; can read an
energy budget as a convergence test; and can save, ship and reload a model state.

## Site architecture

```
docs/modelling-tour/
  07-stepping-a-column.qmd
  08-turbulent-heat-exchange.qmd
  09-moisture-and-buoyancy.qmd
  10-dry-convection.qmd
  11-dry-rce.qmd
  12-moist-rce.qmd
  _tour/
    assets.py       NEW  browser/native asset resolution, lifted out of tables.py
    stepping.py     NEW  the time loops and the snapshot figure
    states.py       NEW  save/load a sympl state as .npz, with provenance
    budgets.py      NEW  TOA and surface energy budget diagnostics
    tables.py       becomes a caller of assets.py
  _data/
    rce_dry_equilibrium.npz     NEW
    rce_moist_equilibrium.npz   NEW
scripts/
  generate_tour_equilibria.py   NEW
  serve_wheel.py                MOVED from docs/radiative-transfer/_live/
```

Pages are `format: live-html`, `engine: jupyter`, including
`docs/_includes/climt-live-boot.qmd` and declaring their `_tour` modules and any
`_data` asset under `pyodide: resources:`, exactly as pages 01–06 do.

### Deleting `radiative-transfer/09-live-rce.qmd`

Every inbound reference was checked; there are five, and none is a public link:

- `docs/_quarto.yml` — sidebar entry removed; six added under Modelling Tour. The
  comment at line 125, which explains why HTML appearance is declared top-level and
  names `09-live-rce.qmd` as one of two `live-html` families, becomes tour-only.
- `docs/modelling-tour/04-gray-equilibrium-tested.qmd:305` — the one real
  cross-link, re-pointed at page 07, where it becomes a forward reference within
  the same tour.
- `docs/modelling-tour/01-emissivity-spectrum.qmd:9` — a comment pointing at 09 for
  how to preview an unreleased wheel. That instruction moves into the header of
  `docs/_includes/climt-live-boot.qmd`, which is where anyone looking for it now is.
- `docs/superpowers/HANDOFF-in-browser-rce-demo.md` and the 2026-07-19 plan —
  historical records; left alone, as plans are.
- `radiative-transfer/index.qmd`, the site landing page and the README never linked
  it. The RT track loses a page it did not advertise.

**The include consolidation.** `docs/_includes/climt-live-setup.qmd` (227 lines) had
exactly one consumer: the live-RCE page. It is deleted, leaving `climt-live-boot.qmd` as the
site's single boot include for all twelve live pages. This undoes the deliberate
two-include duplication introduced by fix #2 of the tranche 1 whole-branch review,
which existed to stop two page families breaking each other; there is now one family.

The helpers do **not** move into the boot include, which stays minimal by design.
`integrate_with_snapshots` and `integrate_to_equilibrium` become
`_tour/stepping.py`, delivered by the `resources:` mechanism and declared only by
pages 07–12, so pages 01–06 download nothing new.
`docs/radiative-transfer/_live/rce_helpers.py` is deleted together with the
keep-two-copies-in-sync contract its docstring describes: `_tour/stepping.py` is the
single copy, imported rather than duplicated into a cell.
`tests/test_live_rce_demo.py` (217 lines) folds into `tests/test_modelling_tour.py`.

### Corrections inherited, not copied

The live-RCE page ships `DT_HOURS = 2, N_STEPS = 400` — 33 simulated days, leaving the column
**29 K above** its own analytic skin temperature, with a comment admitting it. Page
07 must not inherit those numbers. Measured on that exact configuration: `dt = 12 h`
with ~675 steps converges to **0.65 K** of the analytic skin temperature at the same
wall clock (5.1 s vs 5.3 s native). Page 07 also uses `UnytBackend` instead of
`DataArrayBackend` (3.2× per step on a gray column) and `nz = 28` instead of 18,
which roughly halves the top-of-column layer-to-layer gaps so the skin layer reads
as flat, for +16% per step.

## Page designs

Each page keeps tranche 1's anatomy: chapter anchor, visible prescribed state, the
calls, a figure, one clearly marked knob, Physics and Code exercises, and a "going
deeper" link. The craft thread stays cumulative.

**Moisture and radiation, page by page.** The arc is not monotone in moisture —
chapter 10 goes dry then moist, and chapter 12 does it again — so each page's
configuration is fixed here rather than left to the implementer:

| Page | Water vapour | LW table | Convection | Time loop |
|---|---|---|---|---|
| 07 | none (q = 0) | gray, then 14-band | none | yes, live to equilibrium |
| 08 | none (q = 0) | gray | none | yes, live |
| 09 | evolves; supplied by the surface | 14-band | none | yes, live, short |
| 10 | none (q = 0) | none — no radiation at all | dry adjustment | **no** |
| 11 | none (q = 0) | 14-band, CO₂ only | dry adjustment | yes, from shipped state |
| 12 | evolves | 14-band | Emanuel + condensation | yes, from shipped state |

Two consequences to state on the pages themselves. Page 11's column is genuinely
dry, so its 14-band radiation sees CO₂ alone and it equilibrates **well below**
Earth's surface temperature; no number on it is comparable to tranche 1's. And the
surface temperature difference between pages 11 and 12 is, to first order, the water
vapour greenhouse plus its feedback — the quantity tranche 1's page 06 drew from
prescribed profiles, now produced by a model that supplies its own water.

**Two facts verified against the code, not assumed.** `CorkLongwaveRadiation` with
`single_band_gray_lw` does not take `specific_humidity` as an input at all — its
optical depth is a function of pressure alone — so pages 07 and 08 are dry by
construction and the gray column reproduces page 04's analytic profile whatever the
reader does to humidity. That is worth saying on page 07 rather than leaving the
reader to wonder.

And climt's default `specific_humidity` is **0.0 everywhere**. The live-RCE page's
14-band column is therefore already a CO₂-only atmosphere, while its prose describes
"spectral windows letting surface radiation escape" as though water vapour were
present. Page 07 inherits the run but not the ambiguity: it states that this column
is dry, which is precisely what makes page 12's comparison worth drawing.

### 7. Stepping a column

*Chapters 9–10, and the bridge from tranche 1.*

The first page in the tour where anything moves. Builds a single column over a 1–2 m
slab with a prescribed absorbed shortwave flux at the surface, steps it with
`AdamsBashforth`, and draws the evolution: temperature, heating rate and up/down
longwave flux profiles at several times, plus OLR and surface temperature as time
series.

**The reveal.** Started from an arbitrary profile, the gray column *finds* the
analytic equilibrium that page 04 prescribed and verified. Page 04 checked that
ℋ ≈ 0 on a profile handed to it; page 07 produces that profile from nothing.
Then, changing one string — the absorption table — the 14-band column finds a
different one, and the two are compared.

- **Knob:** mixed-layer depth. 1 m and 5 m reach the *same* equilibrium at very
  different speeds — the cleanest possible demonstration that heat capacity sets
  response time and not equilibrium.
- **Craft:** `TendencyComponent` vs `Stepper`; `AdamsBashforth`; the `UnytTimeDelta`
  gotcha (a plain `datetime.timedelta` raises `UnitOperationError` under
  `UnytBackend`, because `total_seconds()` returns a bare float that does not cancel
  the `/s`); update order — the prognostic state is applied before the diagnostics,
  or the stepper carries flux diagnostics forward at their pre-step value;
  `state["time"]`.
- **Claims to test:** the converged profile matches page 04's analytic T(τ) within a
  stated tolerance; the top level sits within ~0.7 K of `Te/2^(1/4)`.

### 8. Turbulent heat exchange

*Chapter 10.*

Adds `SimpleBoundaryLayer(surface_fluxes='bulk')` to page 07's dry column. The
scheme computes Frierson (2006) bulk fluxes of heat, moisture and momentum from its
own exchange coefficient, applies them implicitly, and reports the applied fluxes as
diagnostics — so the flux the reader plots is the flux the model used.

**The reveal.** Page 04 derived radiative equilibrium's surface–air temperature
discontinuity as an analytic feature and left it there. This is the process that
erases it. The figure puts the radiative-equilibrium profile and the
radiative-plus-turbulent profile side by side near the surface, with the sensible
heat flux time series beneath.

- **Knob:** wind speed, and roughness length `z0`.
- **Craft:** running a `Stepper` inside the same loop as tendency components;
  diagnostics that are fluxes; the `'bulk'` / `'external'` / `None` modes and the
  double-counting trap that pairing `'bulk'` with a component that already applies
  surface fluxes creates.
- **Claims to test:** the surface–air discontinuity shrinks from its
  radiative-equilibrium value to a stated smaller one; with `surface_fluxes=None`
  the diffusion conserves every column integral.

### 9. Moisture changes the buoyancy

*Chapter 10, continued.*

Turns on the latent flux by giving the slab a saturated surface humidity, and
introduces the virtual potential temperature θ_v = θ(1 + 0.61q) as the quantity
buoyancy actually depends on. Radiation switches to the 14-band table, so the water
vapour the boundary layer supplies is radiatively active — the greenhouse effect
tranche 1 measured, now supplied by the surface rather than typed in.
`GridScaleCondensation` enters here as the minimum moisture sink: the boundary layer
mixes moist air into colder air aloft, and without a sink the column supersaturates.

- **Knob:** surface relative humidity, which moves the Bowen ratio.
- **Craft:** composing components that write the same variable; deriving θ_v from
  state quantities with correct units; why the order of two components that both
  touch specific humidity is a physical choice, not a formatting one.
- **Claims to test:** the latent flux exceeds the sensible flux over a saturated
  ~288 K surface; θ_v − θ ≈ 0.61 q θ to stated precision; the boundary-layer top
  moves by a stated amount relative to page 08; the magnitude of the supersaturation
  that `GridScaleCondensation` removes.

### 10. Dry convection

*Chapter 11.*

**The only page in the tranche with no time loop.** Prescribes a superadiabatic
profile, calls `DryConvectiveAdjustment` once, and shows before and after in both T
and θ. Static stability and enthalpy-conserving mixing are seen on their own, with
no convergence noise in the way — the same single-call discipline that made tranche
1 legible.

- **Knob:** the depth and strength of the initial instability.
- **Craft:** calling a `Stepper` directly — it takes a timestep and returns
  `(diagnostics, new_state)`, no tendencies anywhere; writing a conservation check
  as the test of whether you understood the scheme.
- **Claims to test:** column-integrated enthalpy conserved to machine precision
  (`tests/test_conservation.py` already asserts this, so the page's claim is guarded
  by a test that exists); θ uniform across the adjusted layers.

### 11. Dry RCE

*Chapter 12.*

Everything dry, together, to equilibrium: longwave radiation, slab surface,
boundary layer, dry convective adjustment. The page **loads a shipped equilibrium
state** rather than spinning one up, then runs perturbation experiments from it.

**The reveal.** The convecting layer sits on a ~9.8 K/km dry adiabat — not the
6.5 K/km every tranche 1 page assumed — and the tropopause **emerges** as the level
where convection stops, instead of being prescribed as a 200 K cap. Compared
against page 07's radiative equilibrium: a much cooler surface, a similar upper
atmosphere.

- **Knob:** CO₂. Doubling it and integrating to the new equilibrium is a **measured
  dry climate sensitivity** — the number tranche 1's page 05 could only compute as a
  forcing.
- **Craft:** ordering several components defensibly; the TOA and surface energy
  budgets as a convergence test rather than eyeballing a curve; loading and saving
  a model state.
- **Claims to test:** the lapse rate in the convecting layer is within a stated
  percentage of dry adiabatic; the shipped state's TOA imbalance is below threshold;
  the 2×CO₂ surface warming falls in a stated range.

### 12. Moist RCE

*Chapter 12 and beyond.*

Adds `EmanuelConvectionPython`. Loads the shipped moist equilibrium and perturbs it.

**The reveal.** Latent heating relaxes the lapse rate toward the moist adiabat — a
value near the 6.5 K/km tranche 1 assumed, so the assumption is finally *explained*
rather than replaced. Precipitation appears as a diagnostic, and the moisture budget
closes against the surface latent heat flux. The climate sensitivity differs from
page 11's dry one, and the page says why.

- **Knob:** CO₂, and surface temperature.
- **Craft:** `ImplicitTendencyComponent` — the third and last component protocol,
  and the warning it provokes (see below); precipitation and convective
  diagnostics; timestep sensitivity, demonstrated rather than asserted.
- **Claims to test:** the lapse rate lies between dry and moist adiabatic and near
  the latter through the lower troposphere; precipitation ≈ latent heat flux / L_v
  at equilibrium; the sensitivity value.

**A warning the page will print, and should explain rather than obey.**
`EmanuelConvectionPython` is an `ImplicitTendencyComponent`, so constructing an
`AdamsBashforth` around it makes sympl emit: *"Using an ImplicitTendencyComponent
in sympl TendencyStepper objects may lead to scientifically invalid results."*
**Emanuel goes inside `AdamsBashforth` anyway** — that is how climt has always run
it, and it works. The warning is generic to the component protocol, not a finding
about this scheme at these timesteps, and restructuring the loop around it would
buy nothing and cost the reader a special case to remember.

What the page owes the reader is the sentence explaining it, because the warning
*will* appear in their browser output on the first run and an unexplained warning
in a teaching page reads as a bug. One line: this is sympl saying it cannot verify
that an implicitly-formulated component is safe inside a multi-stage stepper; here
it is, and page 12's own timestep-sensitivity experiment is where you check that
claim yourself. `_tour/stepping.py` therefore keeps exactly two component
categories, tendency and stepper, as page 07 introduced them.

## Shipped equilibrium states

Two single-column states: `rce_dry_equilibrium.npz` (page 11) and
`rce_moist_equilibrium.npz` (page 12). Page 07 ships nothing — watching the approach
*is* that page, and its gray run is ~18 s in-browser.

**Where they live.** `docs/modelling-tour/_data/`, beside `earth_spectrum_lw.npz`,
using the mechanism that directory's README already documents: `project: resources:`
publishes them, `pyodide: resources:` stages them into the Pyodide filesystem, and a
`_tour` module locates them. At `nz = 28`, one column, a dozen quantities, these are
a few kB — four orders of magnitude below the table already committed, so the
"every clone carries the blob" cost accepted in tranche 1's item 19 does not
re-open.

`_tour/tables.py` already owns browser/native path resolution. That logic is lifted
into `_tour/assets.py` and `tables.py` becomes a caller — the one refactor this
tranche makes to tranche 1 code, in service of the current goal.

**Contents.** Values, dims and units per quantity, plus provenance: climt version,
table name, `nz`, `dt`, step count, slab depth, prescribed `SOLAR`, CO₂, the
component list, and the final TOA and surface imbalances. `_tour/states.py`
reconstitutes by calling `get_default_state()` for the page's components and
overwriting arrays, so a sympl change fails loudly at a named quantity instead of
silently producing a state with wrong dims.

**Staleness, without the 200-day bill.** These are **not** wired into
`build_experiments.py --check`. That dependency-hash machinery is what made a
class-attribute one-liner in `cork/lw/component.py` re-run 57 600 five-minute steps
twice, and a content hash over `cork/**/*.py` cannot distinguish a no-op from a
physics change. The guard is instead a **residual test**: load the shipped state,
step it ten steps under the page's own component list, assert the TOA imbalance
stays under threshold and the surface temperature drifts less than ~0.05 K. Seconds
to run, in CI, and it fails exactly when the physics genuinely moved — at which
point `scripts/generate_tour_equilibria.py` is re-run deliberately. Its defaults
match what shipped, which is the lesson `generate_tour_spectrum_table.py` taught by
not doing so.

**What the pages do with them.** Load, change one knob, integrate until it settles.
With a 1–2 m slab and an effective feedback near 2 W m⁻² K⁻¹ the response time is
roughly three weeks, so three e-folds is ~70 simulated days: order 1700 steps at
dt = 1 h. Those are estimates from the slab's heat capacity, not measurements; Task
0 replaces them. Each page states plainly that it is perturbing an equilibrium
computed for one specific configuration.

## Cost and stability

Per-step budget, from tranche 1's measured cost model (native with
`NUMBA_DISABLE_JIT=1`, the correct proxy for Pyodide, × ≈3.5 for the browser):

| stack | native | ≈ browser |
|---|---|---|
| page 07 gray: LW + slab | ~8 ms | ~0.03 s |
| page 07 non-grey: 14-band LW + slab | ~36 ms | ~0.13 s |
| page 11 dry RCE: + SBL + dry adjustment | ~41 ms | ~0.14 s |
| page 12 moist RCE: + Emanuel + condensation | ~44 ms | ~0.15 s |

Radiation dominates by ~30×, so every component this tranche adds is nearly free.
**Cost is set almost entirely by the timestep**, and that is the number we do not
have. Gray radiation is stable at dt = 12 h and produces NaN at 24 h. Implicit BL
diffusion and the adjustment schemes should tolerate large steps. Emanuel is the
unknown: the scheme is conventionally run at 10–20 minutes, and if that holds here,
page 12's ~70-day perturbation is ~5 000 steps ≈ 12 minutes in-browser rather than
the few minutes estimated above.

If the measurement is bad, the levers, in order: a thinner slab (response time is
linear in depth, and depth is already a knob whose lesson *is* response time); a
perturbation the column answers faster than a CO₂ doubling; then fewer levels.
**Not** a lever: degrading the physics to make a cell finish. A run that takes
minutes while the lecture continues is the accepted trade.

## Library change

One, and it belongs in this tranche because six pages invite the experiment that
triggers it. When the timestep passes stability, `cork/lw/component.py` emits
`RuntimeWarning: invalid value encountered in subtract` at
`net_band = up_band[b] - down_band[b]` and returns NaN. In the browser that appears
as a blank figure rather than an error. The component should detect the
non-finite result and raise with a message naming the timestep as the likely cause.
Small, testable, and it protects every page in the tranche.

## Task 0 — measure before designing pages

No `.qmd` is written until these are measured and recorded:

1. Largest stable timestep for `EmanuelConvectionPython` in this configuration, and
   for the full dry stack.
2. Steps to equilibrium, dry and moist, from a cold start — and from the shipped
   equilibrium after a CO₂ doubling.
3. Per-step cost of each page's exact component list, native, no numba.
4. **Start-independence:** two different initial conditions converging to the same
   equilibrium within tolerance. Load-bearing — shipping an equilibrium state is
   only honest if the equilibrium does not depend on how you got there, and page 12
   says so in prose.
5. The magnitude of the supersaturation `GridScaleCondensation` removes on page 09's
   stack, reported as a number the page quotes.

Every number the six pages quote descends from this task.

## Testing

The pattern tranche 1 established, unchanged: page computation lives in `_tour/` as
importable, natively testable Python; page cells stay thin; `{pyodide}` cells cannot
run in CI, so the helpers carry the assertions.

- `tests/test_modelling_tour.py` gains the tranche 2 helpers, and absorbs
  `tests/test_live_rce_demo.py`.
- Each page's physics claim above becomes an assertion, on the table that page
  actually uses. **This is the tranche 1 lesson that cost the most**: pages 1–3 ran
  on the 56-band table while every test ran on the 14-band one, so every number
  rewritten in that task was guarded by nothing. Tests are parametrised over the
  table the page declares.
- The residual test for both shipped equilibrium states.
- `_tour/states.py` round-trip: save a state, reload it, assert every quantity
  matches in values, dims and units.
- `_tour/assets.py` keeps `tables.py`'s existing tests — resolution from a working
  directory nowhere near the docs tree, and graceful fallback when an asset is
  absent — extended to the equilibrium states.
- Static fallback figures in `_artifacts/`, regenerated by script and committed, so
  pages degrade if Pyodide fails.

## Deferred

**`CorkShortwaveRadiation` is deliberately not used in this tranche.** It has not
been stress-tested, and this tranche is not the place to find out. Every page
prescribes an absorbed shortwave flux at the surface instead. The consequence is
that surface albedo and zenith angle are not knobs anywhere in the tour, and the
surface energy balance is one-sided in the shortwave — worth stating on the pages,
and worth a tranche of its own once the component has been exercised.

## Out of scope

- Clouds, sea ice, and the polar-column capstone.
- `SecondBEST`, `BucketHydrology`, land surface components.
- Any compiled component: `RRTMGLongwave`, `RRTMGShortwave`, the Fortran
  `EmanuelConvection`, `SimplePhysics`, `BergerSolarInsolation` are all absent in
  the browser.
- Multi-column or latitude-resolved configurations.
- Replacing `tutorial/`. Tutorial 2's subject — time integration — is absorbed here
  in context, as tranche 1 absorbed Tutorial 1's.

## Risks

| Risk | Mitigation |
|---|---|
| Emanuel's stable timestep forces a 10+ minute browser cell | Task 0 measures it first. Levers in stated order; physics is not one of them. |
| The shipped equilibrium depends on the path taken to it | Task 0 item 4 tests start-independence explicitly before either state is committed. |
| A shipped equilibrium silently goes stale after a physics change | Residual test in CI, deliberately *not* dependency hashing. |
| Deleting the live-RCE page loses a working demo | Its figure function and integration loop move to `_tour/stepping.py` intact and are covered by the tests that already cover them; page 07 renders the same figure with corrected numbers. |
| Prose numbers drift from the configuration they were measured on | Every quoted number comes from a cell the reader can run, per the tranche 1 correction; tests are parametrised over the table the page declares. |
| The moist column runs away or collapses at some knob setting | Pages state the range each knob was tested over; the NaN guard turns a blank figure into a message. |

## References

- Course notes: <https://joymonteiro.github.io/principles_planetary_climate/>,
  chapters 10–12.
- Tranche 1 spec: `docs/superpowers/specs/2026-08-12-modelling-tour-radiation-design.md`
- Tranche 1 plan, including the whole-branch review and the reviewer's Minor list:
  `docs/superpowers/plans/2026-08-12-modelling-tour-radiation.md`
- In-browser demo spec (CORS, wheel hosting):
  `docs/superpowers/specs/2026-07-19-in-browser-nongrey-rce-demo-design.md`
