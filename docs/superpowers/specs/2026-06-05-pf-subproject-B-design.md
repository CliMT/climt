# Sub-project B design — "What it takes to build a longwave scheme that rivals RRTMG"

**Status:** design (approved via brainstorming). Next step: `superpowers:writing-plans`.
**Date:** 2026-06-05
**Seed:** `docs/superpowers/specs/2026-06-05-pf-subproject-B-seed.md`
**Parent (sub-project A):** `docs/superpowers/specs/2026-06-03-pf-co2-adjustable-hifi-earth-lw-design.md`
**Investigation log (all numbers/arc):** `docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md`

## Purpose

A student-facing pedagogical artifact telling the discovery story behind the
shipped picket-fence longwave scheme: how a warm-biased correlated-k table was
root-caused and fixed to rival RRTMG. It teaches correlated-k, its real failure
mode, and the fix — grounded in an actual debugging session — and demonstrates
"fast vs faithful." Written against the **now-shipped** capability A delivered.

## Audience

Advanced undergraduate (intro atmospheric science). Assumes basic
Beer–Lambert / blackbody emission; explains k-distributions, g-points, and the
"correlated" assumption from scratch. Gentle math, every term defined, generous
scaffolding.

## Deliverables

1. **A Sphinx page** (`docs/*.rst`, renders on ReadTheDocs) carrying the
   **complete, self-contained** discovery narrative + pre-rendered figures. A
   student reads it end-to-end without a kernel.
2. **A companion notebook** in `examples/` (runs in the `climt` conda env) that
   reproduces every figure and hard number on the page as a runnable cell — the
   "evidence you can re-run" to the page's "claims."
3. **One small code helper:** `scripts/experiments/tropopause.py` — a robust
   cross-climate tropopause locator (the only new code in B), with unit tests.
4. **One figure script:** `scripts/experiments/make_subproject_B_figures.py` —
   single entry point that regenerates all of the page's hero figures.
5. **A forward-hook** written into the Pyodide spec for the future
   JupyterLite-interactive version of the notebook (no code in B).

### Format decision (page vs notebook)

Page = the full story (standalone prose + pre-rendered figures). Notebook = the
lab bench (re-runnable evidence). Chosen because it matches the
Sphinx/ReadTheDocs home (page teaches without a kernel), keeps the static build
independent of notebook execution, and ports cleanly to a future JupyterLite
version (the prose is already standalone; the notebook just becomes the
interactive surface).

**Page → notebook contract:** every figure and every hard number on the page is
reproducible by a labeled cell in the notebook.

## Narrative spine (the page)

Discovery-driven: lead with the symptom, investigate live, hit and reject the
dead-ends, arrive at the fix. The rejected attempts are first-class teaching
beats. Sections §1–7 are the spine (the surface-temperature detective story);
§8–9 collapse into a single short "two more things this bought us" coda.

1. **The symptom.** Moist RCE with the shipped-*before* 4-band table sits
   **+11.9 K** warmer at the surface than RRTMG (history: +20.1 K), surface
   humidity runs away **3.2×** (6.49 vs 2.02 g/kg), tropopause **broken at
   583 hPa**. Frame the question: line data, or the scheme?
2. **Primer (scaffolded).** Just enough RT: band-mean transmission can't use a
   mean k (⟨exp⟩ ≠ exp⟨⟩); line-by-line is exact but costly; the k-distribution
   trick (sort k(ν) → smooth k(g)); g-points as quadrature; the "correlated"
   assumption for stacking layers. Retold from `docs/notes/corr-k-explainer-seed.md`.
3. **Isolating the culprit.** Line-by-line truth (linepyline at D=1.66) *matches*
   RRTMG → line data is fine. PF over-traps ~17 (dry) to ~36 (moist) W/m² → it's
   the **correlated-k scheme**, not the physics input. The pivot of the story.
4. **The dead-ends (first-class beats).** More g-points (8→16→32) **doesn't
   help**; finer vertical resolution **doesn't help**; wing-boost hacks rejected.
   Lesson: when adding resolution doesn't move the error, you're fighting a
   *structural* assumption, not a numerical one — a spectral-correlation floor.
5. **Root cause #1 — the continuum (dominant moist lever).** The H₂O MT_CKD
   continuum, folded into the line-sorted k-distribution and interpolated
   linearly across log-spaced H₂O nodes, over-traps ~17 W/m². **Fix: decouple
   it** — line-only k-distribution + a separate band-grey continuum term, log-X
   interpolated. The same separation RRTMG uses (`tauself`/`taufor`).
6. **Root cause #2 — band structure (secondary lever).** The k(g) figure: the
   800–1800 cm⁻¹ "window" band lumps the transparent 8–12 µm window with the
   opaque H₂O ν₂ band. Dry it's invisible (both halves transparent); moist it
   grows a low-g shelf (window) + high-g cliff (ν₂ lines), ~30× step near
   g≈0.5 — the fingerprint of two populations in one band. **Fix:** a 14-band
   partition isolating the window and refining far-IR/ν₂.
7. **The payoff + the diagnostic-trust lesson.** Before→after: moist T_sfc
   +11.9 K → **+1.9–2.1 K** (indep. verify +2.13 K), dry +0.8 K, q_sfc runaway
   gone (2.2 vs 2.0), tropopause 583 → **196 hPa co-located**. Teaching subtlety:
   fixed-profile single-column LBL−PF reads **+11 W/m² for both** the candD
   prototype and the shipped table, yet full RCE converges to ~2 K — the
   fixed-profile metric **over-states** the error a self-adjusting column
   actually incurs. A lesson about trusting the right diagnostic.
8–9. **Coda — "two more things this bought us."** (a) *Fast vs faithful:* making
   the table faithful made the Python assembly the bottleneck (467× slower than
   RRTMG); `@njit(parallel=True)` on the two hot loops → **~540× speedup,
   bit-for-bit** → faithful scheme now *faster* than RRTMG (37.9 vs 45.2 µs/col).
   Accuracy and speed weren't a trade-off — the cost was un-compiled assembly,
   not physics. (b) *The free knob:* CO₂ adjustable 10–10000 ppm via one runtime
   axis, log-k interpolated; off-node OLR-vs-LBL bias = on-node bias →
   interpolating CO₂ between nodes is **fidelity-free** (leave-one-out: log-k
   5.5% vs linear-k 30.9%).
10. **Coda / answer key.** The final 14-band recipe stated cleanly + pointers to
   the notebook.

### The answer-key recipe (for §10)

14 bands (edges cm⁻¹: `10, 250, 350, 500, 630, 700, 800, 980, 1080, 1180, 1250,
1400, 1600, 1800, 3250`), 8 g-points; line-only k-distribution + decoupled
band-grey H₂O continuum (log-X interpolated); runtime CO₂ axis (10 log nodes,
10–10000 ppm); quadrilinear interp — T (linear), log p (linear), log X_H₂O
(linear-k), log X_CO₂ (log-k). Default = `earth_low_res_lw`; "before" table
preserved as `earth_low_res_lw_4band_ngpt2_before.nc`.

## The companion notebook

`examples/` notebook (climt conda env), one section of cells per page beat:

1. **Setup.** Imports; load both tables (`earth_low_res_lw_4band_ngpt2_before.nc`,
   `earth_low_res_lw.nc`); construct `PicketFenceLongwave` + RRTMG.
2. **Reproduce the symptom.** Moist RCE before-table vs RRTMG → +11.9 K / runaway
   q / 583 hPa, tropopause located by the new θ-curvature diagnostic.
3. **k-distribution sandbox.** Build and plot a k(g) curve from a band (reader
   *constructs* the abstraction), then the dry-vs-moist window-band shelf+cliff.
4. **Isolate the culprit.** Load pre-computed linepyline LBL OLR (`.npz` in
   `debug_data/`); overlay PF vs LBL vs RRTMG.
5. **Dead-ends, runnable.** Swap g-points 8→16→32 and bump vertical resolution;
   show the forward-flux comparison barely moves — the reader feels the floor.
6. **Before→after.** Re-run moist RCE with the shipped 14-band table → +1.9 K, q
   recovered, tropopause 196 hPa. Side-by-side with cell 2.
7. **Diagnostic-trust subtlety.** Compute fixed-profile single-column LBL−PF (+11
   for both) next to converged RCE ΔT (~2 K) — the over-statement lesson as
   numbers the reader generates.
8. **Coda cells.** Throughput numbers (`bench_pf_vs_rrtmg.py`) + a CO₂-sweep cell.

The notebook is committed **pre-executed** (outputs visible on GitHub/nbviewer)
so it reads as a finished artifact without a kernel.

### JupyterLite boundary (forward-degradation)

RRTMG (Fortran) and linepyline do **not** run in Pyodide. Cells 2, 4, 6, 7
depend on them, so they are authored to read pre-computed `.npz` artifacts
(already in `debug_data/`) **as the data source**, with the live-RRTMG call shown
but guarded (`if RRTMG_AVAILABLE:`). Same notebook degrades gracefully: full live
in the climt env, pre-baked data in a future browser port. The PF scheme itself,
the k(g) sandbox (cell 3), and the CO₂ sweep (cell 8) are pure-Python and run
live in JupyterLite.

## The tropopause diagnostic (the only new code)

**What it does.** Given a column's `air_temperature` and `air_pressure`,
return the tropopause pressure, robust across dry / moist / runaway-warm climates
(the constant-θ+10K marker failed here, producing the 583 hPa artifact).

**Primary criterion — θ-curvature.** Compute θ(p); the tropopause is the kink
where θ transitions from the convective near-constant troposphere to the rapidly
increasing stratosphere — the pressure of maximum upward curvature of θ vs
log-p (peak of d²θ/d(ln p)², lowest robust peak above a floor). Directly captures
"top of convection," which the broken marker was reaching for. **This is the
headline value** for the 196/583 hPa claims.

**Cross-check — cold-point.** Also report the temperature-minimum pressure as a
secondary marker. The notebook plots both; the page teaches the caveat
(cold-point rides high/low-pressure in the moist tropics, so the two diverge
there — itself a teachable moment).

**Home.** `scripts/experiments/tropopause.py` with a pure function
`find_tropopause(temperature, pressure) -> {p_curvature, p_coldpoint}` (numpy in,
dict out, no climt-state coupling). Rationale: it's experiment/diagnostic tooling
(matches the `scripts/experiments/` convention), notebook-importable,
independently unit-testable, and pure numpy → JupyterLite-safe. **Not** a shipped
`climt` diagnostic component (that would be scope creep beyond a teaching
artifact).

**Tests (TDD)** in `tests/`, synthetic profiles with known answers:
- Clean RCE-like profile (constant-θ troposphere + warm stratosphere) →
  curvature finds the prescribed kink.
- Pathological runaway-warm profile that broke the old marker → asserts it does
  **not** return the ~583 hPa artifact and lands in the physically sensible upper
  troposphere.
- Degenerate / isothermal column → returns gracefully (documented sentinel, no
  crash).

## Figures & data assets

- **Single figure script:** `scripts/experiments/make_subproject_B_figures.py`
  regenerates **all** of the page's hero figures cleanly, with consistent
  publication styling, into `docs/radiative_transfer/figures/`. It pulls from the
  existing experiment scripts / `.npz` data internally, but there is one obvious
  place to find and re-run figure code (no scattered per-figure invocations).
- Hero figures: the k(g) shelf+cliff, the before/after RCE pair, the throughput
  bar (at minimum).
- The page embeds the **copied** figures from `docs/radiative_transfer/figures/`
  so the doc build does not depend on `debug_data/` (a scratch dir).
- The notebook reads the `.npz` artifacts already in `debug_data/` for the
  Fortran/LBL-dependent cells.

## JupyterLite forward-hook

Add a short forward-looking section to the Pyodide spec
(`docs/superpowers/specs/2026-03-29-pyodide-pure-python-support-design.md`)
recording that this notebook is a target deliverable for in-browser deployment,
and mapping the boundary: pure-Python-live cells (k(g) sandbox, CO₂ sweep, PF
scheme) vs cells needing pre-baked `.npz` (anything touching RRTMG/linepyline).
No code in B for this — just the written hand-off.

## Non-goals (scope discipline)

B is documentation + teaching + the one tropopause helper. Explicitly out of
scope (the seed warns against re-opening these):
- No band-structure re-optimization.
- No "line-side log-X" work.
- No new scheme capability.
- No shipped `climt` diagnostic component.
If the discovery story surfaces a tempting improvement, it gets **noted, not
built**.

## Testing & verification

- `tropopause.py` → TDD unit tests (clean / pathological-runaway / degenerate),
  per above.
- Notebook → executes end-to-end in the climt env
  (`jupyter nbconvert --execute`) as the acceptance check; committed
  pre-executed.
- Page → builds clean under Sphinx (no broken refs/figures).

## Source numbers (for figures/claims)

| metric | before (4-band) | after (14-band CO₂) |
|---|---|---|
| moist RCE T_sfc vs RRTMG | +11.9 K (hist. +20.1 K) | +1.9–2.1 K (indep. +2.13 K) |
| dry RCE T_sfc vs RRTMG | — | +0.8 K |
| moist q_sfc vs RRTMG (g/kg) | 6.49 vs 2.02 (3.2× runaway) | 2.2 vs 2.0 (~15%) |
| tropopause (moist) | 583 hPa (broken) | 196 hPa (co-located) |
| single-column LBL−PF | — | +11.5 moist / +11.0 dry (corr-k floor) |
| A5 log-k vs linear-k (LOO p95) | — | 5.5% vs 30.9% |
| throughput, 1000 cols | — | PF 37.9 µs/col vs RRTMG 45.2 (~540× over start) |

## Reusable artifacts (from sub-project A)

- `docs/notes/corr-k-explainer-seed.md` — worked correlated-k explainer.
- `scripts/experiments/plot_kg_curves.py` → k(g) figure.
- `scripts/experiments/bench_pf_vs_rrtmg.py` — fast-vs-faithful benchmark.
- `scripts/experiments/eval_band_structure.py` — single-column LBL−PF per band.
- `scripts/experiments/eval_co2_interp_accuracy.py` — A5 probe.
- `scripts/experiments/rce_{moist,dry}_pf_vs_rrtmg.py` — RCE (`--table`, `--co2`).
- `scripts/experiments/lbl_olr_tiebreak.py --co2 <ppm>` — LBL OLR (linepyline env).
- `debug_data/*.npz`, `debug_data/*.png` — pre-computed data and figures.
- Tables: `earth_low_res_lw_4band_ngpt2_before.nc` (before) vs
  `earth_low_res_lw.nc` (after, default).
