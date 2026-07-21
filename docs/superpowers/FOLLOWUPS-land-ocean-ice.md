# Land-Ocean-Ice Surface Components — Issues & Follow-ups

This branch (`feature/land-ocean-ice-components`, PR #219 → `develop`) bundles
**two** feature sets implemented in successive sessions on the same branch:

1. **Land-surface** (commits `c2cce5e`..`5e7fcf4`): two-layer `BucketHydrology`
   (`num_layers=2`), the modular `SecondBEST` land model, its five pluggable
   BEST process objects, a `soil` vertical grid, and `von_karman_constant`.
   Plan: `docs/superpowers/plans/2026-07-20-land-surface-components.md`.
2. **Ocean & ice** (commits `bf1d974`..`3a7834e`): `LandMask`, the
   `_core/snow_ice_column.py` solver, `SeaIce`/`LandIce`, `IceSheet` as a
   deprecating shim, `SlabSurface` q-flux + Ekman, `_core/surface_fluxes.py`,
   and `DataOcean`.
   Plan: `docs/superpowers/plans/2026-07-20-ocean-ice-components.md`.

Every task in both halves passed an individual spec+quality review, and each
half passed a final whole-branch review. **Nothing below blocks merge** — these
are the deferred findings, deliberate deviations, and plan defects that were
found and fixed, captured explicitly so reviewers and future work have the full
picture.

Full suite on the branch HEAD: **401 passed, 2 skipped, 132 failed**. All 132
failures are pre-existing (unbuilt Cython/Fortran extensions: RRTMG, Emanuel,
Simple Physics, DCMIP) and reproduce identically on the pre-work baseline — zero
regressions from either half.

---

## A. Needs a human decision

### A1. `SecondBEST` `z_mid` reference height — RESOLVED
~~`z_mid = Rd·T_air/g` (~8–9 km scale height) produced implausibly small drag
coefficients.~~ **Fixed:** the surface-layer reference height is now the
hypsometric height of the lowest model level above the surface,
`z = (Rd·T/g)·ln(p_surface/p_lowest)` (~tens of metres), which also makes
`surface_air_pressure` (previously a declared-but-unread input) actually used.
Covered by `test_surface_layer_reference_height_uses_lowest_model_level`;
`TestSecondBEST` cached output reseeded.

### A2. Flux/Neumann boundary is a quasi-steady algebraic constraint (Important, ocean-ice)
`_core/snow_ice_column.py`'s `Flux` boundary condition is an algebraic
constraint rather than a time-integrated flux (inherited from the original
`IceSheet`). This produces a `dt`-independent, forcing-proportional
energy-conservation artifact, which makes `atol`-based forced-conservation tests
impossible for `SeaIce`/`LandIce` (they use zero-forcing no-op + directional
energy-response tests instead). The ocean-ice plan mandated "no new physics," so
this was left as-is. **Recommend a tracked issue** against the `snow_ice_column`
Flux BC if true energy closure is later required.

---

## B. Plan defects found and fixed during implementation (human-visible)

These are cases where the plan's literal code was wrong and the implementer
deviated, with the fix independently verified in review. Called out because the
shipped code intentionally differs from the written plan.

### B1. SeaIce/LandIce basal-flux sign inversion (ocean-ice)
The plan's `bottom_bc = Flux(heat_flux_into_sea_water_due_to_sea_ice)` inverts
the physics — `solve_column` treats a positive `Flux` as heat *entering* the
column, but that CF quantity is positive when heat *leaves* ice into the ocean.
Corrected to `Flux(-Q)` in both `SeaIce` and `LandIce`; verified that a warm
ocean thins the ice.

### B2. Singular mid-month SST solve (ocean-ice, `DataOcean`)
The plan's bidiagonal `0.5·(mmₘ + mmₘ₊₁) = meanₘ` system is singular for a
12-month cycle (det = 0). Replaced with the canonical Taylor–Williamson–Zwiers
cyclic-tridiagonal `(1/8, 3/4, 1/8)` system (the exact monthly average of the
piecewise-linear interpolant through mid-month values), with an independent
numerical-integration test. Residual simplification (unchanged): `mid_month`
assumes equal-length months while `interp_time` uses real calendar days —
acceptable AMIP-grade approximation.

### B3. IceSheet clobbered non-owned cells (ocean-ice, Critical — fixed)
During the `IceSheet`→shim refactor, `surface_temperature` on `"sea"` cells
(owned by neither `SeaIce` nor `LandIce`) was being zeroed. Fixed via a 3-way
masked merge that passes non-owned cells straight through; cache reseeded.

### B4. SecondBEST zeroed non-land columns (land-surface, Critical — fixed)
`SecondBEST.array_call` used `initialize_numpy_arrays_with_properties` (zero
fill) and only wrote land/land_ice columns, so `"sea"` columns came out at 0 K.
Fixed by pre-filling every output from the matching input before the per-column
loop (mirrors `IceSheet`'s convention); covered by a passthrough regression test.

---

## C. Deliberate deviations from the plans (by design)

- **Two-layer `BucketHydrology` uses scalar `deep_soil_moisture_content` /
  `deep_soil_temperature` fields**, not a new `soil_levels` grid dim
  (`get_land_grid` raises for 3-D land). Keeps the shallow fields at `["*"]` in
  both modes and preserves the `num_layers=1` bit-for-bit guarantee. Physics
  identical.
- **`BestSubsurfaceTransport` is self-contained** (its own scipy sparse
  tridiagonal) rather than sharing the ocean-ice `_core/snow_ice_column.py`, so
  the land-surface half runs independently of the ocean-ice half. Unifying the
  two solvers is a reasonable follow-up.
- **`SecondBEST` energy-conservation test** runs in a no-phase-change regime
  (zero soil water/ice, sub-freezing profile) so the freeze/melt clamp is a
  no-op and the base class's tight `atol=1e-3` closure check is meaningful,
  rather than relaxing `atol` to swallow the clamp.

---

## D. Deferred minor follow-ups (fix-in-follow-up / accept-as-is)

### Land-surface
- `SecondBEST` `num_soil_layers` kwarg is stored but unused (soil levels come
  from the grid dim); passing a non-default value silently no-ops. Consider
  raising for unsupported values, or wiring it up.
- `BestSurfaceLayer` `float()` coercion breaks on array inputs — an undocumented
  scalar-only contract (orchestrator only calls it with scalars today).
- `BestSurfaceFluxes` full `"B"`-present exfiltration path has no *standalone*
  unit test. Note: it *is* exercised in production (`BestSoilProperties` always
  returns `B`/`K_H0`, so every `SecondBEST()` step goes through it) and the
  integration test drives it; a targeted unit test would still be good insurance.
- `BucketHydrology` `num_layers` validation uses `in (1, 2)`, so `1.0`/`True`
  pass (both route to the correct branch; harmless). Tighten to `isinstance int`.
- `k_soil = 2.0` and the Bolton saturation-humidity coefficients are hardcoded
  in the bucket/`SecondBEST` code (plan-sanctioned; consider a kwarg/registry).
- `z0 = 0.01/0.001` surface-roughness lengths are hardcoded inline in
  `SecondBEST` (area-type-keyed, like the soil property tables).
- `SurfaceFluxes.momentum_flux` is computed but discarded by the orchestrator.
  (`surface_air_pressure` is now read — see A1.)
- Dead placeholder `B = 4.0` line in `fluxes.py` (always overwritten/unused).
- `docs/api/BucketHydrology.qmd` shows a stale signature — regenerate via
  `quartodoc build` (the `.qmd` pages are docstring-generated build artifacts;
  `SecondBEST` is now registered in `docs/api/_metadata.yml` + `index.qmd`).
- `SoilProperties.__call__` base has no docstring (only the class does).

### Ocean & ice
- Nearest-neighbour longitude matching uses planar (not circular) distance —
  latent wrong-neighbour bug near the 0/360 seam for *irregular user-supplied*
  masks/SST grids only; the bundled 2° grid is unaffected. Fix:
  `d = min(|Δ|, 360−|Δ|)`. Recurs in both `LandMask` and `DataOcean` — consider
  a shared fix.
- `LandMask._weights` / `DataOcean._weights` caches are built once per instance
  with no grid-shape/hash guard — stale if one instance is reused across
  differing grids (not exercised by tests).
- `DataOcean` rebuilds its `RegularGridInterpolator` + weights every
  `array_call` (efficiency; the `_weights` cache is currently a no-op).
- `DataOcean` `compute_fluxes=True` path is implemented but untested; its
  `isfinite` runtime check is an `assert` (disabled under `-O`).
- `_core/surface_fluxes.py` `bulk_fluxes` hardcodes `cp=1004`/`Lv=2.5e6` instead
  of `get_constant`; Ekman uses sea-water `cp≈3985` vs the mixed-layer `cp≈4185`
  (~5% internal inconsistency).
- `_core/horizontal_operators.py` uses a module-global `_RADIUS` (won't reflect
  a `planetary_radius` change; test-isolation smell).
- `SlabSurface` Ekman: single-column `lat2d.ndim==1` reshape branch is
  dead/untested; the equator `np.sign(0)` fix has no regression test (Gaussian
  grids avoid exact-0 latitude).
- `SlabSurface`'s `ocean_heat_transport_convergence` diagnostic reports the TOTAL
  (prescribed + Ekman) when `include_ekman=True`; the Ekman-only breakdown is the
  separate `ekman_heat_transport_convergence` diagnostic (see `HISTORY.rst`).
- `LandMask` generated `.nc`: lat/lon coords lack `units` attrs
  (`degrees_north`/`east`); smoke test checks only point samples, not
  `flag_values`/`flag_meanings`/dtype/ordering.
- `surface_downward_heat_flux_in_sea_ice` is now `SeaIce`-only, so it reads `0.0`
  on land/land_ice columns (the old monolith computed it there). Documented in
  the `IceSheet` deprecation note; a behavior change for callers using that
  diagnostic on land.

### Process / infrastructure notes
- `benchmarks/` is **not** in pytest collection — kernel signature changes there
  aren't caught by the suite. Check `benchmarks/` call sites when changing any
  jit kernel signature.
- Two `TestSlabSurface` cache-key failures involving
  `height_on_soil_interface_levels` are **pre-existing** (reproduce at the
  pre-work baseline), not caused by this feature.
