# Design — CO₂-adjustable, RRTMG-fidelity Earth-LW picket-fence table

**Date:** 2026-06-03
**Status:** approved design (sub-project A of a two-part effort)
**Related:** investigation log `docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md` (exps #16–#28); corr-k seed `docs/notes/corr-k-explainer-seed.md`.

## Context / motivation

The shipped Earth-LW picket-fence table runs **far too warm vs RRTMG** in radiative-convective equilibrium. A long single-column + RCE investigation (see the plan above) root-caused the warm bias to the **correlated-k longwave scheme**, not the line data:

- On an identical (T,p,q) profile, line-by-line truth (linepyline at PF's diffusivity) gives OLR that matches RRTMG; PF over-traps by ~17 (dry) to ~36 (moist) W/m².
- The dominant moist lever is the **H₂O MT_CKD continuum**: folding it into the line-sorted k-distribution and interpolating its ~X² scaling linearly across log-spaced H₂O nodes over-traps ~17 W/m². Decoupling it (line-only k-distribution + a separate band-grey continuum term, **log-X interpolated**) fixes it — the same separation RRTMG uses (`tauself`/`taufor` with analytic `selffac ∝ H₂O²·p`).
- The secondary lever is **band structure**: a 14-band partition that isolates the transparent window and refines the far-IR/ν₂ regions.
- Neither more g-points (8→16→32) nor finer vertical resolution help; the residual is the correlated-k spectral-correlation floor.

Combining the two fixes (table "candD": 14 bands, decoupled+log-X continuum, ngpt=8) brought converged moist-RCE surface temperature from **+20.1 K** (current shipped 4-band) to **+2.7 K** vs RRTMG, dry to +1.0 K, with the water-vapour feedback tamed (q_sfc 1.33 vs RRTMG 0.95; old default 6.94). The fix currently lives in experiment scripts/tables; this sub-project productionizes it.

This is **sub-project A**. Sub-project **B** (a student-facing pedagogical page + companion notebook documenting "what it takes to build a longwave scheme that rivals RRTMG") follows in its own spec, written against the shipped capability A delivers.

## Goal

Ship a production Earth-LW correlated-k scheme and table that:
1. Eliminates the H₂O-continuum over-trapping (decoupled + log-X continuum).
2. Uses the refined 14-band partition.
3. Is **CO₂-adjustable at runtime over 10–10000 ppm** from a single table.
4. Becomes the default `earth_low_res_lw`, matching RRTMG RCE to a few K.

## Non-goals (future / separate)

- The matching **shortwave** table (and its CO₂ axis).
- The **line-side log-X interpolation** (lines also ~convex in X; untested, riskier — affects all tables).
- A **robust cross-climate tropopause diagnostic** (the constant-θ+10 K marker mislocates the tropopause for runaway-warm columns; WMO lapse-rate or θ-curvature criterion needed — belongs with sub-project B's notebook).
- Further **band optimization** beyond candD's 14 bands.
- Regenerating **Titan / TRAPPIST-1e** tables.

## Design

### The table

- **Band partition:** 14 bands, edges (cm⁻¹) `10, 250, 350, 500, 630, 700, 800, 980, 1080, 1180, 1250, 1400, 1600, 1800, 3250`.
- **g-points:** 8 per band (16 gave no measurable benefit).
- **H₂O continuum — decoupled:** the line k-distribution is built **without** the MT_CKD continuum; the continuum is stored separately as a band-grey field `continuum_kappa(nband, nT, nP, nX_H₂O)` and added to every g-point's optical depth at runtime, **log-X interpolated**. The continuum is H₂O-only and CO₂-independent → **no CO₂ axis on the continuum**.
- **CO₂ — runtime axis:** the line k-distribution gains a CO₂-VMR dimension:
  `k_coefficients(ngas=1, nband, ngpt, nT, nP, nX_H₂O, nX_CO₂)`.
  Each `(X_H₂O, X_CO₂)` node is a **true premixed mixture** sampled and k-distributed independently → no gas-overlap approximation. **10** log-spaced CO₂ nodes spanning 10–10000 ppm (≈ ⅓-decade spacing: 10, 21.5, 46.4, 100, 215, 464, 1000, 2154, 4642, 10000). Runtime interpolation is **quadrilinear**: T (linear), log p (linear), log X_H₂O, and **log X_CO₂**.
  - *CO₂-axis interpolation scheme:* default to **log-k vs log-X_CO₂** (geometric), since linear-in-value across log-spaced nodes over-estimates convex quantities — the exact failure that caused the continuum bug (exp #25); CO₂ band-mean k vs amount is similarly convex/saturating. The final choice (log-k vs linear-k) is **confirmed by the A5 CO₂-accuracy check**, but log-k is the starting design. This is independent of the existing H₂O-line axis (which stays linear-in-k; switching it is the deferred "line-side log-X" non-goal). The ⅓-decade node spacing also keeps interpolation error small regardless.
- **Naming / default:** regenerated as the canonical **`earth_low_res_lw`**, replacing the current 4-band/ngpt=2 file in place (the feature is new — no users depend on the old contents). The current `earth_low_res_lw_8gpt.nc` and the original `earth_low_res_lw.nc` content are kept under explicit names as "before" references for sub-project B.

### Code changes (extend the existing H₂O-axis machinery)

1. **Generation** — `scripts/generate_pf_tables_linepyline.py` + `scripts/pf_table_builder/kappa_sampling.py`:
   - Sample line kappa on the 4-D `(T, p, X_H₂O, X_CO₂)` grid (CO₂ as a swept absorber).
   - Sample the H₂O continuum once on `(T, p, X_H₂O)` (CO₂-independent), as the existing `--decouple-continuum` path already does.
   - Build the line k-distribution per `(X_H₂O, X_CO₂)` node; band-grey-average the continuum per band.
   - One-time **offline cost grows ~10×** (the CO₂ axis) — likely hours; document and run on an adequate machine.
2. **Table format** — `scripts/pf_table_builder/netcdf_writer.py`: add a `co2_vmr_grid` variable and the CO₂ dimension on `k_coefficients`. `continuum_kappa` already supported.
3. **Loader** — `climt/_components/picket_fence/optics/correlated_k.py`: read `co2_vmr_grid` and the extra k dimension.
4. **Interpolation** — `interpolate_k`: add the CO₂ axis (quadrilinear; linear-in-k vs log-X_CO₂). `interpolate_continuum` unchanged (H₂O-only).
5. **Component** — `climt/_components/picket_fence/lw/component.py`: declare `mole_fraction_of_carbon_dioxide_in_air` as a **runtime input**, convert to VMR as needed, and pass it into the optical-depth interpolation. (Today CO₂ is baked in and not read.)
6. **Backward compatibility:** tables lacking a CO₂ axis and/or `continuum_kappa` behave exactly as today (the new axes are optional and detected by presence).

### Performance pass

Consolidate the band / g-point / weight-accumulation loops — currently **pure Python** in `lw_transport`, with only the inner single-g-point sweep `@njit(prange over columns)` — into a **single `@njit(parallel=True)`** routine that loops over a flattened band (or band×column) axis and accumulates weighted fluxes inside the compiled kernel. This removes the Python-level `nband×ngpt×nlev×ncol` overhead, which matters now that band count ~3.5×'s and per-call work grows (continuum + CO₂ interpolation). Must reproduce current fluxes bit-for-bit on existing tables (regression-guarded by the test suite).

### Validation & testing (gates the default switch)

- **True baseline:** run converged RCE with the *actual* current default (`earth_low_res_lw`, ngpt=2) for the honest "before" number (the +20 K figure is for the `_8gpt` variant).
- **Single-column LBL fidelity** (`scripts/experiments/eval_band_structure.py`, extended for CO₂): confirm |LBL−PF| holds across a CO₂×H₂O grid, not just at 376 ppm.
- **CO₂-interpolation accuracy:** compare the table at off-node CO₂ (e.g. 400, 2000 ppm) against LBL to verify log-CO₂ interpolation fidelity.
- **RCE at ≥2 CO₂ values** (e.g. 280 and 1120 ppm), moist + dry, converged, vs RRTMG.
- **Test suite:** existing PF tests stay green; add unit tests for (a) loading + interpolating the CO₂ axis, (b) the decoupled continuum load/interp, (c) the consolidated njit kernel reproducing the old fluxes.

### Risks / notes

- **Generation cost** (4-D grid × LBL spectrum) is the main practical risk; mitigate by documenting the run and using an adequate machine / overnight run.
- **CO₂ log-interpolation fidelity** must be verified (band-mean CO₂ k is expected smooth in log-CO₂, but the window/wing bands should be checked).
- The continuum band-grey approximation is validated only for Earth H₂O continuum at the tested conditions; the CO₂×H₂O grid validation extends that coverage.

## Key files

- `scripts/generate_pf_tables_linepyline.py`, `scripts/pf_table_builder/kappa_sampling.py`, `scripts/pf_table_builder/netcdf_writer.py`
- `climt/_components/picket_fence/optics/correlated_k.py` (loader, `interpolate_k`, `interpolate_continuum`)
- `climt/_components/picket_fence/lw/component.py` (CO₂ input), `climt/_components/picket_fence/lw/kernels.py` (njit consolidation)
- New shipped table `climt/_data/picket_fence/correlated_k/earth_low_res_lw.nc`
- Tests under `tests/` and `tests/pf_table_builder/`
- Examples/docs referencing `earth_low_res_lw` (verify still valid; CO₂ now settable)

## Relationship to sub-project B

B (the student-facing pedagogical page + runnable notebook telling the discovery story and demonstrating fast-vs-faithful) is written **after** A lands, against the shipped scheme and final numbers. It gets its own design → plan cycle.
