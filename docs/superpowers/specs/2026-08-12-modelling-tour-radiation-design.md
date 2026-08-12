# A Modelling Tour of the Climate System — Radiation Tranche

**Date:** 2026-08-12
**Status:** design, awaiting review
**Scope:** tranche 1 of a multi-part section. Radiation only.

## Summary

Add a new top-level docs section, **A Modelling Tour of the Climate System**, that
supplies the computational half of the existing course notes at
<https://joymonteiro.github.io/principles_planetary_climate/> (EC2213, *Principles
of Planetary Climate*). Those notes derive analytically and contain no figures and
no numerical examples. The tour supplies exactly those, live in the browser.

This tranche covers course chapters 5–8 (Shell Model → Transition to a continuous
atmosphere → Basics of Radiative Transfer → Gray Gas radiative equilibrium) in six
pages. Every page prescribes a profile, makes one or a few `CorkLongwaveRadiation`
calls, and reads CORK's **per-band** diagnostics. **Nothing integrates in time.**

Two learning objectives carry equal weight: the science, and *using climt*.

## Motivation

### The gap

The course notes stop at the point where computation would begin. The Mar 9 lecture
builds the Schwarzschild equation, its integral solution, `OLR = I⁺(τ_∞)`, back
radiation `I⁻(0)`, and the weighting function `W ∝ τ* e^{−τ*}` peaking at `τ* = 1`.
It then says:

> "this is all for a particular wavelength… the OLR that we had from our shell is
> actually this integrated over all lambda and over all theta"

Students never perform that integral. CORK performs exactly that integral, and its
per-band diagnostics let students **un-do** it — decomposing the OLR they already
know into the spectral pieces that produce it.

It also closes a question left open in the Jan 29 lecture: if the shell is gray with
ε ≠ 1, the radiating level "is not defined at all, it's somewhere in between." With
per-band optical depth, the answer is that there is no single radiating height — there
is a *spectrum* of them, and it can be plotted.

### Why now

Capabilities that did not exist when the course was taught:

- `CorkLongwaveRadiation` / `CorkShortwaveRadiation` expose per-band optical depth,
  transmittance, up/down flux and heating rate
  (`climt/_components/cork/lw/component.py`, `diagnostic_properties`).
- The Earth LW table carries **runtime CO₂ and H₂O axes**, so both are live knobs.
- climt runs in the browser via Pyodide (`docs/radiative-transfer/09-live-rce.qmd`),
  so a 95-student class needs zero installs. `tutorial/tutorial.md` currently opens
  on virtualenvs and conda, which is where large classes lose people.

## Decisions

| Question | Decision |
|---|---|
| Relationship to the existing Radiative Transfer track | Parallel, cross-linked. Tour = student-facing and experiment-first; RT track stays the research-grade theory reference. |
| Scope of this spec | Radiation only, in depth. Chapters 10–12 are tranche 2. |
| Time integration | **None.** Prescribed profiles and single calls throughout. |
| Spectral resolution | Two tables: shipped 14×8 for quantitative work, new ~100-band table for spectra. |
| g-points in the new table | 4, subject to a validation gate (see Risks). |
| Where the new table lives | Static asset on the docs site, not in the wheel. |
| Opening page | Continuity with the students' own shell model, then the spectrum breaks it. |
| climt-craft teaching | A cumulative thread across all six pages, not an intro page. |
| Backend | `UnytBackend` (see Backend note). |

## Audience and objectives

Undergraduate/early-graduate students who have taken, or are following, EC2213.
Readers arriving cold from the web are a secondary audience: each page names its
chapter and links to it, so the tour is usable without the course.

**Science objectives.** By the end, a student can explain and *compute*: why a single
emissivity is a fiction; where in the spectrum the atmosphere is transparent; what
sets the radiating level and why it differs by wavelength; why the gray radiative
equilibrium profile is not an equilibrium for a real atmosphere; why CO₂ forcing is
logarithmic and comes from band wings; and why water vapour drives both the strongest
feedback and a hard limit on OLR.

**Craft objectives.** A student can build a model state, inspect a component's
input/diagnostic contracts, set and read state quantities with correct units and
dims, distinguish `mid_levels` from `interface_levels`, sweep a parameter, and
combine components.

## Site architecture

```
docs/modelling-tour/
  index.qmd                  landing: what this is, how it maps to the course
  01-emissivity-spectrum.qmd
  02-window-measured.qmd
  03-radiating-level.qmd
  04-gray-equilibrium-tested.qmd
  05-co2-knob.qmd
  06-water-vapour-limit.qmd
  _tour/
    soundings.py             prescribed-profile builders (natively testable)
    spectra.py               Planck, brightness temperature, band-plot helpers
    tables.py                fetch/locate the hi-res table (browser + native)
  _artifacts/                static fallback figures
  _data/
    earth_spectrum_lw.npz    hi-res table, served as a site asset
```

`docs/_quarto.yml` gains a navbar entry and a docked sidebar between "Radiative
Transfer" and "Experiments", plus `modelling-tour/**/*.qmd` in `project: render:`.

Pages are `format: live-html`, `engine: jupyter`, reusing
`docs/_includes/climt-live-setup.qmd` for the Pyodide boot.

## Page designs

Every page follows one anatomy:

1. **Chapter anchor** — a callout naming the course chapter and what it derived.
2. **Prescribed state** — built inline from `_tour/soundings.py`, visibly, not hidden.
3. **One or a few CORK calls.**
4. **A per-band figure.**
5. **A live knob** — one clearly marked value to change and re-run.
6. **Exercises**, split into *Physics* and *Code*.
7. **Going deeper** — link to the RT track chapter.

### 1. ε is not a number

*Chapter 5 (Shell Model), Chapter 7.*

Opens with the students' own single-ε shell model and its prediction. Then computes
per-band OLR on a realistic sounding (Tₛ = 288 K, RH = 80%, Γ = 6.5 K/km, isothermal
200 K stratosphere) and converts to brightness temperature.

The reveal: their ε is an average over a spectrum whose brightness temperature swings
from ~200 K (CO₂ core) to ~285 K (window) — verified during design; see Appendix A.
Overlaying Planck curves reproduces the classic Hanel/IRIS figure.

- **Knob:** surface temperature.
- **Craft:** components as objects; `get_grid`, `get_default_state`; a component call
  returns `(tendencies, diagnostics)`; reading `.diagnostic_properties`.

### 2. The window, measured

*Chapters 6, 7.*

Per-band optical depth and transmittance profiles. The atmospheric window stops being
an assertion and becomes a measured width. Directly answers the Jan 29 lecture's
question — "what would the window be if you had a pure nitrogen atmosphere?" — by
turning CO₂ and H₂O down to the bottom of their axes.

- **Knob:** CO₂ concentration and a humidity scale factor.
- **Craft:** the state *contract* — `.input_properties`, units, dims; setting
  `mole_fraction_of_carbon_dioxide_in_air` and `specific_humidity`.

### 3. Where photons come from

*Chapter 7.*

Builds the weighting function `W = τ* e^{−τ*}` **per band** from
`longwave_optical_depth_per_band`, cumulated from the top down, and finds each band's
`τ* = 1` level. Plots those levels on the temperature sounding, next to the brightness
temperatures from page 1 — showing that `T_b(band) ≈ T(z where τ* = 1)`.

This is the page that resolves the "radiating height is somewhere in between" problem:
it is a distribution, not a height.

- **Knob:** humidity scaling, to watch emission levels rise.
- **Craft:** `mid_levels` vs `interface_levels`; deriving quantities from diagnostics.

### 4. Gray equilibrium, tested

*Chapter 8.*

Chapter 8 derives `T(τ) = T_e[(1 + τ_∞ − τ)/2]^{1/4}`, the skin temperature
`2^{−1/4} T_e`, and the surface discontinuity, then stops. This page:

1. Prescribes that exact analytic profile.
2. Runs `CorkLongwaveRadiation` with `single_band_gray_lw` and `D = 2` and confirms
   the heating rate ℋ ≈ 0 — **the analytic solution validated by a radiation code.**
3. Hands the *same* profile to the 14-band table and shows ℋ ≠ 0 — the gray
   equilibrium is not an equilibrium for a real atmosphere, band by band.

Step 3 is the strongest single argument for non-grey radiation in the whole tour, and
it costs two calls.

- **Knob:** τ_∞.
- **Craft:** prescribing a whole profile; `TendencyComponent` vs `Stepper`; why ℋ is
  a tendency, not a state variable.

### 5. The CO₂ knob

*Chapter 4 (forcing, sensitivity, Planck response).*

Sweeps CO₂ across the table's runtime axis (10–10 000 ppm, 10 log nodes) with the
14×8 table, plotting OLR against log CO₂. Students measure the forcing themselves.
Design-time values on the standard sounding: 280 → 560 ppm gives **−3.59 W/m²**,
560 → 1120 gives −3.86 (canonical ≈ 3.7).

Then decomposes the forcing by band: the 15 µm core (630–700 cm⁻¹) barely moves
because it is saturated; the wings do the work. This is *band saturation*, computed.

- **Knob:** the sweep range, and the sounding's lapse rate.
- **Craft:** sweeping a parameter — loop, mutate state, collect diagnostics.

### 6. Water vapour, and the limit

*Chapter 4 and beyond.*

OLR as a function of Tₛ, computed twice: at **fixed relative humidity** (water vapour
responds) and at **fixed specific humidity** (it does not). The gap between the two
curves *is* the water vapour feedback, drawn rather than asserted.

Pushed to high Tₛ on saturated soundings, OLR flattens — the Simpson–Nakajima limit
and the runaway greenhouse. The table's axes support this: H₂O VMR reaches 1.0
(q ≈ 0.38 kg/kg), temperature 1000 K, pressure 100 bar. Per-band diagnostics show the
window closing, which is *why* the limit exists.

Introduces `CorkShortwaveRadiation` so the energy balance is two-sided.

- **Knob:** relative humidity.
- **Craft:** two components together; `num_longwave_bands` and the per-band axis.

## Spectral tables

Two tables with distinct jobs.

| | `earth_low_res_lw` (existing) | `earth_spectrum_lw` (new) |
|---|---|---|
| Resolution | 14 bands × 8 g-points | ~100 bands × 4 g-points |
| Ships in | the wheel | the docs site, as a static asset |
| Size | 3.0 MB | ~5.4 MB (CO₂ axis thinned 10 → 5 nodes) |
| Cost/call (browser) | ~0.12 s | ~0.42 s |
| Used by | pages 4, 5, 6 — anything quantitative or swept | pages 1, 2, 3 — single-call spectra |

**Format must be `.npz`.** The pure wheel has no scipy, and `_load_netcdf_table` in
`climt/_components/cork/optics/correlated_k.py` requires `scipy.io.netcdf_file`.

**Hosting.** The docs site is same-origin with the page, so a static asset can be
fetched straight into the Pyodide filesystem — no CORS problem (unlike GitHub release
assets, per the 2026-07-19 spec) and no wheel bloat. Only pages 1–3 download it.
`load_k_table` already accepts a filesystem path, so no library change is needed.

**Generation.** Offline, via `scripts/generate_cork_tables_linepyline.py`, which
already accepts `--bands` and `--ngpt`. Requires the linepyline conda environment.
The CO₂ axis is thinned to 5 log nodes over 10–10 000 ppm to hold the file near
5.4 MB; k varies smoothly in log CO₂, but this needs checking (see Risks).

**Band edges must be a refinement of the existing 14-band grid** — that is, all 15
existing edges (10, 250, 350, 500, 630, 700, 800, 980, 1080, 1180, 1250, 1400, 1600,
1800, 3250 cm⁻¹) appear in the new edge list, with additional edges subdividing them.
This makes validation exact: the new table's per-band OLR aggregates onto the 14-band
grid by simple summation, with no interpolation error to disentangle from real
disagreement.

## Measured cost model

All figures single-column, `nz = 18`, native with `NUMBA_DISABLE_JIT=1` — the correct
proxy for Pyodide, which has no numba. The browser factor of **≈3.5×** is calibrated
from the existing live-RCE page's measured "~2 min for 400 steps".

**All of these assume the `NpzFile` prerequisite below has landed.**

| call | native (no numba) | ≈ browser |
|---|---|---|
| CORK LW 14 × 8 | 35 ms | 0.12 s |
| CORK LW ~100 × 4 (projected) | ~120 ms | ~0.42 s |
| CORK LW 1-band gray | 7.0 ms | 0.025 s |
| CORK SW 3 × 2 | 12.7 ms | 0.045 s |
| DryConvectiveAdjustment | 2.8 ms | ~10 ms |
| SimpleBoundaryLayer | 2.4 ms | ~8 ms |
| EmanuelConvectionPython | 2.6 ms | ~9 ms |
| GridScaleCondensation | 0.7 ms | ~3 ms |

Cost is linear in `nband × ngpt` at **~0.29–0.31 ms per quadrature point** (measured
by tiling the real table to 28/56/112/224 bands). Radiation dominates everything
else by ~30×.

A 20-point sweep with the 14×8 table costs ~2.4 s in-browser, which sets the design
rule: **sweeps use 14×8; spectra use the hi-res table for a single call.**

### Browser-available components

Compiled, therefore **absent** in the browser: `RRTMGLongwave`, `RRTMGShortwave`,
`EmanuelConvection` (Fortran), `SimplePhysics`, `BergerSolarInsolation`,
`DcmipInitialConditions` (`setup.py`, `ext_modules`). Everything else is pure Python.
The tour therefore uses CORK, and tranche 2 will use `EmanuelConvectionPython` and
`SimpleBoundaryLayer` in place of their compiled equivalents.

### Backend note

Use `UnytBackend`. For the single calls in this tranche it is neutral (34.1 ms vs
35.1 ms), but it removes ~5.5 ms/step of state-container overhead when stepping, which
matters for tranche 2: gray stepping goes 7.98 → 2.52 ms/step (**3.17×**), 14-band
40.12 → 34.92 ms/step (1.15×). Adopting it now keeps both tranches consistent.

**Gotcha to document:** `sympl.AdamsBashforth` under `UnytBackend` requires
`climt.UnytTimeDelta` as the timestep. A plain `datetime.timedelta` raises
`UnitOperationError` (`degK` + `degK/s`), because `timedelta.total_seconds()` returns
a bare float and does not cancel the `/s`.

## Prerequisites

Two library changes land before the pages, both outside the tour proper.

1. **Materialise `NpzFile` in `load_k_table`.** It currently returns numpy's lazy
   `NpzFile`, which re-decompresses the requested array from its zip on *every*
   `__getitem__` — 16.5 ms for the 3 MB `k_coefficients`, several times per call.
   Returning `{k: npz[k] for k in npz.files}` gives 54.6 → 2.9 ms/call with numba
   (**19×**) and 85.1 → 35.1 ms without (2.5×). Tracked separately; this spec's cost
   model assumes it.

2. **`D` (diffusivity factor) as a constructor argument** on
   `CorkLongwaveRadiation`. Currently a module constant `DIFFUSIVITY_FACTOR = 1.66`
   in `climt/_components/cork/lw/kernels.py`. The default stays **1.66**, so existing
   behaviour and all shipped-table calibrations are unchanged; the tour passes
   `D = 2` explicitly where it reproduces the course notes' algebra. See Notation.

## Notation reconciliation

Course chapter 7 derives the two-stream equations with a diffusivity factor of **2**
(μ = 0.5, the 60° photon). CORK uses **1.66**, the Elsasser value. Left alone, the
tour would silently disagree with the notes, and page 4's ℋ ≈ 0 check would sit off
zero for reasons unrelated to the physics being taught.

Resolution: make `D` settable, use `D = 2` where the tour reproduces the notes'
algebra, then flip to 1.66 as a teaching point about hemispheric averaging.

Symbols follow the course notes throughout: `τ` measured up from the surface
(`τ = 0` ground, `τ_∞` TOA), `τ*` measured down from space, `F⁺`/`F⁻`, `F_net`,
`ℋ` for the heating rate, `B(T)`. Where CORK's variable names differ, pages state
the mapping explicitly.

## Testing

Pyodide cells cannot run in CI. This tranche follows the pattern already established
by `docs/radiative-transfer/_live/rce_helpers.py` and `tests/test_live_rce_demo.py`:

- Every page's computational core lives in `docs/modelling-tour/_tour/` as importable,
  natively testable Python. Page cells stay thin — build state, call helper, plot.
- `tests/test_modelling_tour.py` exercises each helper natively in the `climt` conda
  environment, asserting on the physics the pages claim: per-band OLR sums to the
  broadband diagnostic; brightness temperature in the window exceeds that in the CO₂
  core; the CO₂ doubling forcing falls in 3.0–4.5 W/m²; ℋ ≈ 0 for the analytic gray
  profile at `D = 2`; OLR flattens at high Tₛ.
- **Table validation gate:** `tests/test_spectrum_table.py` compares
  `earth_spectrum_lw` against `earth_low_res_lw` and against line-by-line reference
  output for broadband OLR, per-band OLR aggregated to the 14-band grid, and heating
  rate profiles. This test decides whether 4 g-points survives.
- Static fallback figures in `_artifacts/` are regenerated by a script and committed,
  so pages degrade gracefully if Pyodide fails.

## Out of scope

- **Time integration** of any kind, including the existing live-RCE page's approach.
- **Chapters 10–12** — turbulent heat exchange, dry convection, RCE. These map onto
  `SimpleBoundaryLayer`, `DryConvectiveAdjustment` and `EmanuelConvectionPython` and
  become tranche 2, where integration, sympl Steppers and `UnytTimeDelta` are taught.
- **Clouds.** `longwave_optical_thickness_due_to_cloud` is per-band and would suit
  the tour well, but it belongs with the tranche that has a surface energy balance.
- Porting the remaining `gemini_output` notebooks verbatim. NB1–5 content is absorbed
  and reorganised here; NB6–10 belong to tranche 2 or later.
- Replacing `tutorial/`. The craft thread absorbs Tutorial 1's material in context;
  Tutorial 2's subject (time integration) is explicitly tranche 2.

## Risks

| Risk | Mitigation |
|---|---|
| 4 g-points too coarse at ~100 bands, degrading fluxes or heating rates | The validation gate above runs before any page depends on the table. Fall back to 50 × 8 at the same file size if it fails. |
| Thinning the CO₂ axis 10 → 5 nodes distorts the forcing | Validate the 5-node table's CO₂ forcing against the 10-node table across 280–1120 ppm before adopting. |
| linepyline environment drift — its real API has previously differed from what plans assumed | Generation is a single offline task with a checked-in output and a recorded sha256, as in `MANIFEST.md`. It is not on the critical path: pages 1–3 fall back to the shipped 14×8 table, which Appendix A shows already reproduces the IRIS figure recognisably. The hi-res table improves those three pages; it does not enable them. |
| 5.4 MB extra download on pages 1–3 | Browser-cached after first load; pages 4–6 never fetch it; state the cost in the page callout as the existing live pages do. |
| `NpzFile` prerequisite slips | Pages still work, at 2.5× the quoted cost — 0.3 s and ~1 s per call. Unpleasant, not blocking. |

## Appendix A — design-time verification

Computed during design with the shipped 14-band table. Two separate checks:

**Consistency.** On climt's default state (`nz = 30`), per-band OLR summed to the
broadband diagnostic exactly: 445.50 W/m² both ways.

**Physics.** On the standard sounding (Tₛ = 288 K, RH = 80%, Γ = 6.5 K/km,
T_strat = 200 K, `nz = 40`):

Brightness temperatures: CO₂ core (630–700 cm⁻¹) ≈ 200 K; window (800–1180) ≈ 282–285 K
against a 288 K surface; H₂O rotational band (<500) ≈ 220–240 K; H₂O 6.3 µm ≈ 225 K.

| CO₂ | OLR (W/m²) | window 800–1180 | 15 µm core |
|---|---|---|---|
| 280 ppm | 246.71 | 91.92 | 6.61 |
| 560 ppm | 243.12 | 91.54 | 6.36 |
| 1120 ppm | 239.26 | 90.88 | 6.24 |

Table axes confirmed adequate for page 6: H₂O VMR 1e-6 … 1.0, CO₂ 10 … 10 000 ppm
(10 log nodes), T 50 … 1000 K, p 1 Pa … 100 bar.

## References

- Course notes: <https://joymonteiro.github.io/principles_planetary_climate/>
  (chapters 5–8 for this tranche; 10–12 for tranche 2).
- Lecture transcripts, Jan 29 / Feb 3 / Mar 4 / Mar 9 2026, for framing and notation.
- `docs/superpowers/specs/2026-07-19-in-browser-nongrey-rce-demo-design.md` — Pyodide
  hosting, CORS constraints, and the live-cell template this section reuses.
- `docs/superpowers/specs/2026-05-18-climt-website-design.md` — site structure and the
  notebook-placement convention.
- `climt/_data/cork/correlated_k/MANIFEST.md` — table provenance and regeneration.
- Existing course notebooks in `~/github/model_tour_climate/gemini_output` (NB1–5) and
  tutorials in `~/github/model_tour_climate/tutorial`.
