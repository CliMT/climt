# Seed brief — Sub-project B: "What it takes to build a longwave scheme that rivals RRTMG"

**Status:** pre-brainstorming seed (NOT a design). Start B with `superpowers:brainstorming`, then its own design → plan cycle.
**Date:** 2026-06-05
**Purpose:** Let a fresh session launch sub-project B cold without re-deriving sub-project A. This distills the narrative threads, the hard numbers, the reusable artifacts, and the open questions a B brainstorm must resolve.

## What B is (from the A spec)

A student-facing **pedagogical page + companion runnable notebook** telling the discovery story and demonstrating "fast vs faithful," written **against the now-shipped capability** A delivered. Source of record: `docs/superpowers/specs/2026-06-03-cork-co2-adjustable-hifi-earth-lw-design.md` (§"Relationship to sub-project B").

## The three narrative threads

1. **The discovery story (faithful).** The shipped 4-band Earth-LW table ran far too warm vs RRTMG. A single-column + RCE investigation root-caused the warm bias to the **correlated-k longwave scheme**, not the line data:
   - Line-by-line truth (linepyline at CORK's D=1.66) matches RRTMG; CORK over-traps by ~17 (dry) to ~36 (moist) W/m².
   - Dominant moist lever: the **H₂O MT_CKD continuum**. Folding it into the line-sorted k-distribution and interpolating its ~X² scaling linearly across log-spaced H₂O nodes over-traps ~17 W/m². **Fix: decouple it** (line-only k-distribution + a separate band-grey continuum term, **log-X interpolated**) — the same separation RRTMG uses (`tauself`/`taufor`).
   - Secondary lever: **band structure** — a 14-band partition isolating the transparent window and refining far-IR/ν₂.
   - **Neither more g-points (8→16→32) nor finer vertical resolution help** — the residual is the correlated-k spectral-correlation floor.

2. **Fast vs faithful (the performance arc).** Making the table faithful (14 bands × 8 g-pts + runtime CO₂) made the Python optical-depth/Planck assembly the bottleneck (467× slower than RRTMG). Compiling the two hot loops with `@njit(parallel=True)` (on top of the already-compiled transport kernel) gave **~540× speedup, bit-for-bit**. The faithful scheme is now **faster than RRTMG**. Lesson: accuracy and speed are not a trade-off here — the cost was un-compiled assembly loops, not the physics.

3. **The adjustable knob is free (A5).** The table is CO₂-adjustable 10–10000 ppm via a single runtime axis, interpolated geometrically (log-k vs log-X_CO₂). Validation shows **interpolating CO₂ between table nodes adds no fidelity penalty** — off-node OLR-vs-LBL bias equals the on-node bias. log-k confirmed both by a node-level leave-one-out test and against true line-by-line OLR.

## Hard numbers (for figures/claims — all in the investigation log)

| metric | old 4-band default | new 14-band CO₂ table |
|---|---|---|
| moist RCE T_sfc vs RRTMG | **+11.9 K** (history: +20.1 K for the original) | **+1.9–2.1 K** (indep. verify +2.13 K) |
| dry RCE T_sfc vs RRTMG | — | **+0.8 K** |
| moist q_sfc vs RRTMG (g/kg) | 6.49 vs 2.02 (**3.2× runaway**) | 2.2 vs 2.0 (~15%) |
| tropopause (moist) | **583 hPa (broken)** | **196 hPa (co-located)** |
| single-column LBL−CORK | — | +11.5 moist / +11.0 dry (= candD; the corr-k floor) |
| A5 log-k vs linear-k (leave-one-out p95) | — | **5.5% vs 30.9%** |
| A5 off-node vs true LBL (400/2000 ppm) | — | log-k closer; **CO₂ axis fidelity-free** |
| throughput, 1000 cols | — | **CORK 37.9 µs/col vs RRTMG 45.2** (~540× over the start) |

**Best teaching subtlety:** single-column fixed-profile LBL−CORK reads +11 W/m² for *both* the candD prototype and the shipped table, yet full RCE converges to ~2 K. The fixed-profile metric **over-states** the error a self-adjusting column actually incurs — a lesson about which diagnostic to trust.

## Table recipe (the "answer key")

14 bands (edges cm⁻¹: `10, 250, 350, 500, 630, 700, 800, 980, 1080, 1180, 1250, 1400, 1600, 1800, 3250`), 8 g-points; line-only k-distribution + **decoupled band-grey H₂O continuum (log-X interpolated)**; runtime **CO₂ axis** (10 log nodes, 10–10000 ppm); quadrilinear interp — T (linear), log p (linear), log X_H₂O (linear-k), log X_CO₂ (**log-k**). Default = `earth_low_res_lw`; "before" table preserved as `earth_low_res_lw_4band_ngpt2_before.nc`.

## Reusable artifacts (demos + figures, mostly runnable in the climt env)

- `docs/notes/corr-k-explainer-seed.md` — worked correlated-k explainer + `scripts/experiments/plot_kg_curves.py` → `debug_data/kg_curves_window_band.png` (the k(g) figure; shows why window-band lumping fails).
- `scripts/experiments/bench_cork_vs_rrtmg.py` — the fast-vs-faithful benchmark (CORK vs RRTMG, N columns).
- `scripts/experiments/eval_band_structure.py` — single-column LBL−CORK per band (`co2=` override).
- `scripts/experiments/eval_co2_interp_accuracy.py` — the A5 probe (Part 1 leave-one-out; Part 2 off-node vs LBL).
- `scripts/experiments/rce_{moist,dry}_cork_vs_rrtmg.py` — RCE (`--table`, `--co2` flags) for the headline T_sfc contrast.
- `scripts/experiments/lbl_olr_tiebreak.py --co2 <ppm>` — line-by-line OLR ground truth (linepyline env).
- "Before/after" tables: `earth_low_res_lw_4band_ngpt2_before.nc` vs `earth_low_res_lw.nc`.

## Open questions for the B brainstorm (resolve before designing)

1. **Audience & depth** — advanced undergrad / intro grad atmos-science? How much radiative-transfer background to assume?
2. **Format & home** — static page (where? `docs/`?) + companion notebook; how interactive (CO₂ slider, band-structure toggle, live benchmark)? Pyodide/JupyterLite is supported (see memory) — run-in-browser notebook is feasible.
3. **Narrative balance** — how much of the messy discovery (the rejected attempts: more g-points, finer vertical res, wing-boost hacks) vs the clean final recipe? The rejections are pedagogically valuable.
4. **The tropopause diagnostic (REAL WORK, flagged by the A spec as B's).** The current constant-θ+10 K marker **mislocates the tropopause for runaway-warm columns** (the 583 hPa artifact). The notebook needs a **robust cross-climate criterion** — WMO lapse-rate or θ-curvature. This is a genuine implementation task inside B, not just prose.
5. **Scope discipline** — B is documentation/teaching, not new scheme capability. Resist re-opening band optimization or the deferred "line-side log-X" (both A non-goals).

## Pointers

- A spec: `docs/superpowers/specs/2026-06-03-cork-co2-adjustable-hifi-earth-lw-design.md`
- A plan: `docs/superpowers/plans/2026-06-04-cork-co2-adjustable-hifi-earth-lw.md`
- Investigation log (all results, the full arc): `docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md`
- Knowledge graph: `graphify-out/` (use `graphify query "..."`).
