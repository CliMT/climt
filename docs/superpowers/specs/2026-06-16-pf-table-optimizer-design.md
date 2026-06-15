# Picket-Fence Correlated-k Table Optimizer — Design

**Date:** 2026-06-16
**Branch:** `feature/picket-fence-radiation`
**Status:** Approved design, pending implementation plan

## Problem

Building a picket-fence (PF) correlated-k LW table that matches line-by-line
(LBL) radiative transfer is currently a manual, expert-driven loop: a human
picks band edges, runs the generator, runs PicketFence on a fixed profile,
compares per-band OLR against a re-binned LBL spectrum, reasons about where the
residual lives, and iterates. The 28-experiment campaign in
`docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md` distilled a set of
hard-won heuristics for *where* to place bands and *what* actually moves the
flux error. This project encodes those heuristics into an automated optimizer
that places bands, decides how many are needed, builds the table, compares PF
against LBL, and iterates until acceptance criteria are met.

## Learnings being encoded (from the band-refinement campaign)

These are settled, do-not-re-investigate results that the optimizer must respect:

1. **Over-trapping is driven by band width / spectral inhomogeneity, not
   g-point count** (exp #22). Within a wide inhomogeneous band, a single
   g-point's k maps to different wavenumbers at different layers — the
   correlated-k correlation assumption fails. Adding g-points samples the same
   broken within-band CDF more finely and is radiatively inert. **→ ngpt is NOT
   an optimization knob.** Only narrower, spectrally-homogeneous bands reduce
   the error.
2. **The transparent window wants to be wide** (exp #18). Subdividing the
   window band makes over-trapping worse. **→ window bands are non-splittable.**
3. **Moist over-trapping is the H₂O MT_CKD continuum** (exp #24–25). Decoupling
   the continuum from the line k-distribution and interpolating it in log-X cut
   the moist OLR error roughly in half. **→ continuum decoupling is a knob,
   triggered by moist-specific error.**
4. **Splitting fixes dry, not moist** (exp #23). **→ band splitting is chosen
   from the dry per-band residual; the moist error is addressed by the
   continuum lever.**
5. **Target accuracy** is RRTMG's own agreement with LBL: |LBL − PF| total OLR
   < 5 W/m² moist & dry on a fixed column (RRTMG: 4.1 moist / 7.3 dry).

## Scope

- **v1 target:** Earth LW, built so other scenarios and SW plug in later.
- **Optimization knobs:** band edges + band count (primary); continuum
  decoupling (moist lever); T/p/X grid nodes (last-resort lever). ngpt is
  explicitly excluded.
- **Out of scope for v1:** global black-box refinement (the B approach),
  shortwave, non-Earth scenarios (the code is structured to admit them but they
  are not validated here).

## Key architectural insight: band-independent vs band-dependent work

The expensive linepyline work does **not** depend on band placement. Only two
cheap, pure-numpy steps do. This is what makes the loop fast.

| Work | Depends on | Cost | Cadence |
|---|---|---|---|
| LBL ground-truth RT (OLR + heating-rate profile) | grid/profile | high (linepyline) | once per grid, cached |
| High-res κ(ν,T,p,X) sampling | grid | high (linepyline) | once per grid, cached |
| κ → k-coefficients (re-bin to band edges) | band edges | low (numpy) | every iteration |
| PicketFence RT solve | band edges + table | low (njit) | every iteration |

**Decision (reinterpreting "full linepyline RT per iteration"):** the accept
metric is true LBL radiative transfer (flux + heating rate), NOT the re-bin
proxy used in `eval_band_structure.py`. But because LBL RT is band-independent,
it is computed **once per grid and cached**; linepyline re-runs only when the
grid knob changes. The per-band re-binned LBL spectrum is retained solely for
*residual attribution* (deciding which band to split), never as the accept
metric.

## Cross-environment orchestration

linepyline and climt live in separate conda envs. The optimizer runs as the
loop in the **climt** env (it needs PicketFence). It shells out to
`conda run -n linepyline python -m pf_optimizer.<module>` only to (re)build the
κ-cache and LBL truth, which happens once per grid — rare. This matches the
existing generate-in-linepyline / evaluate-in-climt pattern.

## Architecture

New package `scripts/pf_optimizer/`:

- **`kappa_cache.py`** (linepyline env) — sample κ(ν,T,p,X) on the scenario grid
  via `pf_table_builder.kappa_sampling`, including the decoupled-continuum
  pieces (line-only, continuum-only). Cache to disk keyed by a grid hash. This
  is the dominant linepyline cost.
- **`truth.py`** (linepyline env) — compute LBL ground-truth radiative transfer
  (spectral OLR → total OLR, and the heating-rate profile) for each validation
  column. Cache keyed by (column, dnu, line_shape, continuum settings).
- **`table_build.py`** (numpy only) — given a cached κ-cube + candidate
  `band_edges` (+ decouple flag), build k-coefficients and write a temporary
  `.nc` table. Reuses `pf_table_builder.k_distribution`,
  `planck_fraction`, `netcdf_writer`. No linepyline.
- **`evaluate.py`** (climt env) — run `PicketFenceLongwave` on the validation
  columns, return residuals vs cached truth: total OLR error (moist & dry),
  heating-rate RMSE/level, and per-band re-binned LBL−PF for attribution.
- **`diagnostics.py`** (numpy) — the correlation-breakdown split-point metric,
  per-band residual attribution, window-band detection, and the
  leave-one-node-out κ-interpolation-error metric per grid axis. Pure functions.
- **`actions.py`** (numpy) — the brain: given residuals + diagnostics, return
  the single next action and its reason. Pure functions, encode the learnings
  and guardrails.
- **`loop.py`** — orchestrator: seed → build → evaluate → decide → apply →
  repeat; stopping logic; per-iteration logging.
- **`scripts/optimize_pf_table.py`** — CLI driver (scenario, kind, seed bands,
  targets, max_bands/iters, output path).

## The search engine (hybrid greedy A + inhomogeneity-seeded C)

```
seed band structure (coarse: window + 2–3 absorbing bands)
loop:
  table   = build_table(kappa_cube, band_edges, decouple)        # numpy
  resid   = evaluate(table, validation_columns) vs cached truth  # climt
  if max(|OLR_err_moist|, |OLR_err_dry|) < OLR_TARGET
     and heating_rate_rmse < HR_TARGET:
       ACCEPT → save table + report; stop
  action  = choose_action(resid, diagnostics)                    # the brain
  if action is None:  STALL → save best-so-far + report; stop
  apply(action); record iteration
stop also on: band_count >= max_bands | iter >= max_iters
```

### Action selection (`choose_action`, priority-ordered, one per iteration)

1. **Moist-specific error → decouple continuum.** If
   `|OLR_err_moist| − |OLR_err_dry|` exceeds a threshold and continuum is not
   yet decoupled → action = enable `--decouple-continuum`. (exp #24–25.)
2. **Otherwise split the worst dry band.** Choose the band with the largest dry
   per-band `|LBL − PF|` OLR residual (splitting fixes dry — exp #23). Place the
   new edge at the **correlation-breakdown** point: the ν within the band where
   the layer-to-layer rank-correlation of sorted-κ degrades most (directly
   targets the correlated-k assumption failure from exp #22).
3. **Guardrails.** Window bands (low mean κ / near-unity transmission) are
   non-splittable (exp #18). Enforce `min_band_width` and `max_bands`. ngpt is
   never modified (exp #22).
4. **Grid refinement (last resort).** Only if splitting + continuum are
   exhausted *and* the leave-one-node-out κ-interpolation error on a T/p/X axis
   exceeds the remaining flux gap → refine that axis and trigger a linepyline
   re-cache. Rare by design.
5. **Stall.** If no action is expected to reduce the residual → return None →
   loop reports best-so-far.

### Split-point metric — correlation-breakdown

For a candidate band spanning [ν_lo, ν_hi): at each ν we have κ across the
(T, p, X) layers of a representative column. The correlated-k assumption holds
when the rank-ordering of κ across layers is consistent across ν. Compute, along
ν, the running rank-correlation of the per-layer κ vector against a band
reference ordering; the split point is the ν at which this correlation drops
most sharply (the spectral location where a single g-point can no longer
represent the band coherently). This is computed from the **cached κ-cube** —
no linepyline, no PF run needed to choose the split.

## Acceptance criteria

- `max(|OLR_err_moist|, |OLR_err_dry|) < 5 W/m²` (configurable `OLR_TARGET`).
- Heating-rate RMSE per level `< HR_TARGET` (configurable; default chosen to be
  comparable to RRTMG−LBL on the same column).
- Evaluated across the validation columns (not a single profile).

## Validation columns

Seeded from the existing fixed `forward_profile.npz` (moist + dry), plus 1–2
variants (different T_surf and/or CO₂) so band structure is not overfit to one
column. Configurable list. Each column gets its own cached LBL truth.

## Stall behavior (v1)

On stall: stop, save the best table found, and write a diagnostic report
(residual breakdown, per-band attribution, the action history and why each was
taken). The user decides the next move manually. The global black-box refiner
(approach B) is explicitly deferred — `loop.py` exposes a clean hook for it but
v1 does not implement it.

## Logging

Per the project convention, every iteration is logged (action, reason,
resulting OLR/heating-rate residuals, accept/reject) to the active picket-fence
plan's experiment log. The best table and the full iteration history are saved.

## Testing (TDD)

- **`diagnostics.py`, `actions.py`** — pure functions, unit-tested on synthetic
  κ-cubes with known inhomogeneity / known residual structure (e.g., a
  deliberately inhomogeneous band must yield a split at the planted breakdown ν;
  a moist-dominant residual must yield the continuum action).
- **`table_build.py`** — round-trips a known κ-cube to a table and checks
  k-coefficients/shape against `pf_table_builder` directly.
- **`evaluate.py`** — runs PF on a tiny column, checks residual wiring against a
  cached toy truth.
- **`loop.py`** — integration test on a tiny 2-band grid: assert the loop
  terminates, never touches ngpt, never splits the window band, and produces a
  table that loads in PicketFence.

## Open items / assumptions to confirm during implementation

- Exact `HR_TARGET` value (calibrate against RRTMG−LBL heating-rate RMSE on the
  validation columns).
- Whether the κ-cache stores the full (nT,nP,nX,nNu) cube or a compressed form;
  decide based on memory footprint of the Earth grid.
