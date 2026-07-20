# Ocean & Ice Components Design

**Date:** 2026-07-20
**Branch:** feature/land-ocean-ice-components
**Goal:** Improve climt's surface ocean and cryosphere so it can support realistic and AMIP-style GCM runs: (1) a more realistic **slab ocean** that can carry prescribed arbitrary ocean heat-transport forcing and an optional wind-driven **Ekman transport**; (2) a **DataOcean** `DiagnosticComponent` that prescribes observed SSTs read from file (enabling AMIP-style integrations); and (3) splitting the monolithic `IceSheet` into dedicated **SeaIce** and **LandIce** components over a shared 1-D snow/ice column solver, with the identified defects fixed along the way.

Companion doc: `2026-07-20-land-surface-design.md` (BucketHydrology, SecondBEST).

---

## Problem

| Component | What it does today | Gaps |
|---|---|---|
| `SlabSurface` (`_components/slab_surface.py`) | Local surface energy balance; sea path uses `ocean_mixed_layer_thickness` as slab depth. `TendencyComponent`, numba `jit_compile` kernel with integer `AREA_MAP`. | No representation of ocean heat transport at all. Cannot proxy dynamical ocean heat convergence (needed to avoid unrealistic SST drift) and has no wind-driven transport. |
| — (no data ocean) | climt cannot ingest observed SSTs to force the atmosphere. | Blocks AMIP-style runs. |
| `IceSheet` (`_components/surface_ice.py`) | One `Stepper` handling `land_ice`, `sea_ice`, and snow-on-land in a single per-column Python loop with a tridiagonal implicit heat solve. | Several concrete defects (below); sea-ice basal boundary has no ocean heat flux; land and sea ice physics are entangled in one method. |

### `IceSheet` defects (assessment)

Read of `surface_ice.py`:

| # | Location | Issue |
|---|---|---|
| 1 | `surface_ice.py:384–389` | Two identical `elif area_type == "sea_ice" and sea_ice_thickness > 0` branches — the second is dead code. |
| 2 | `surface_ice.py:293` | `TODO Add ocean heat flux parameterization`. Sea-ice **basal** boundary is a bare melting-temperature Dirichlet (`rhs[0] = melting_temperature`); no ocean heat flux into the ice base, so basal growth/melt is unphysical when the mixed layer is not at freezing. |
| 3 | `surface_ice.py:382–393` | Albedos hardcoded (`0.8` snow, `0.5` ice, `0.2` melt). Not configurable, not state. |
| 4 | `surface_ice.py:352–353` | `print(...)` + `assert False` left in as debugging; crashes the run instead of handling the case. |
| 5 | `surface_ice.py:231` | Snow/ice layer split via `int((1 − snow_fraction)·num_layers) − 1` — coarse quantisation of the snow line onto the grid. |
| 6 | `surface_ice.py:359–364` | Melt can exceed available snow+ice; `sea_ice_thickness` can go negative (no clamp to zero). |
| 7 | `surface_ice.py:190` + `spsolve` | Per-column Python loop calling scipy `spsolve` — not vectorised, no JIT; slow for many columns. |
| 8 | `surface_ice.py:454–457` | Bottom BC branches on `soil_temperature is None` to distinguish land vs sea inside one method — the entanglement the split removes. |

---

## Design overview

Four workstreams:

1. **SlabSurface** gains an additive ocean heat-transport term over sea = prescribed q-flux (arbitrary forcing) + optional diagnosed Ekman convergence. Backward compatible (both zero by default).
2. **DataOcean** — new `DiagnosticComponent` that reads an observed SST dataset and prescribes sea-cell `surface_temperature` / `sea_surface_temperature` as diagnostics, plus surface fluxes.
3. **Ice split** — `SeaIce` + `LandIce` over a shared `_snow_ice_column` implicit solver, fixing defects 1–8.
4. **LandMask** — new `DiagnosticComponent` that sets `area_type` to the present-day Earth land/sea/ice configuration, at arbitrary model resolution (Part 4). Foundational: DataOcean and the land/ice components all key off `area_type`.

All follow climt/sympl conventions and ship with cached-result tests.

---

# Part 1 — Slab ocean with ocean heat transport

## Section 1.1 — Physical model

Over sea, the mixed-layer energy balance becomes:
```
C_ml dT/dt = Q_net + Q_transport
Q_transport = Q_flux_prescribed + Q_ekman
```
where `C_ml = ρ_sw c_sw h_ml` is the existing mixed-layer heat capacity and `Q_net` is today's local net surface flux. Land and ice paths are unchanged.

**Prescribed q-flux (`Q_flux_prescribed`)** — a new input field `ocean_heat_transport_convergence` (W m⁻², default 0). This is the standard slab-ocean "q-flux": arbitrary, user-prescribable forcing that proxies for dynamical ocean heat convergence. Zero default = exactly today's behaviour.

**Ekman energy transport (`Q_ekman`, optional)** — the convergence of wind-driven Ekman *energy* transport, computed from the **wind stress curl**. Define the vertically integrated Ekman **mass** transport per unit width (rotated from the stress by Coriolis):
```
M = (k̂ × τ) / f     →   M_x = τ_y / f ,   M_y = −τ_x / f          [kg m⁻¹ s⁻¹]
```
The Ekman **energy** transport per unit width is `F_Ek = c_sw θ_ml M`, and the heating of the mixed-layer column is its convergence:
```
Q_ekman = −∇·F_Ek = −c_sw ∇·(θ_ml M)
        = −c_sw [ θ_ml (∇·M)  +  M·∇θ_ml ]                          [W m⁻²]
```
- The first term is the **wind-stress-curl / Ekman-pumping** contribution: `∇·M = ρ_sw w_Ek = k̂·∇×(τ/f)` is exactly the curl of the stress over Coriolis (Ekman pumping moving water into/out of the mixed layer). This is the term the "compute transport from wind stress curl" requirement is about, and it is the dominant driver of the gyre-scale energy convergence pattern.
- The second term is horizontal Ekman advection of the SST gradient.

Both terms are retained. (Note: `M = τ/f` is already a *mass* transport per unit width, so the energy flux is `c_sw θ M` with **no** extra `ρ_sw` factor — this corrects a stray density in an earlier draft.)

Enabled by `include_ekman=True`. Requires surface wind stress (`surface_downward_eastward_stress`, `surface_downward_northward_stress`), `sea_water_density`, `heat_capacity_of_sea_water`, and latitude/`coriolis_parameter`. Near the equator `f → 0`; apply the conventional cap (bound `1/f`, or blend to zero within a few degrees of the equator) — documented. Computing the curl/divergence needs horizontal grid geometry (spherical `∂/∂x`, `∂/∂y` from climt's lat/lon spacing); in single-column / degenerate grids the horizontal derivatives are zero, so `Q_ekman` falls back to zero. Diagnostics expose `ekman_pumping` (the wind-stress-curl term `w_Ek`) and `ekman_heat_transport_convergence` separately.

## Section 1.2 — Component surface

Keep `SlabSurface` a `TendencyComponent` with the numba kernel. Changes:
- `__init__(include_ekman=False, equatorial_ekman_cap_latitude=5.0, **kwargs)`.
- New optional inputs: `ocean_heat_transport_convergence` (default 0 via state init), and — only when `include_ekman` — the two stress components, `coriolis_parameter`, `sea_water_density`, and `heat_capacity_of_sea_water` (plus the grid lat/lon already in state for the curl).
- New diagnostics: `ocean_heat_transport_convergence` (echo of total applied), `ekman_heat_transport_convergence` (the diagnosed Ekman energy convergence), and `ekman_pumping` (the wind-stress-curl term `w_Ek`), so runs can separate prescribed vs diagnosed contributions and inspect the curl-driven pumping directly.
- Kernel adds `Q_transport / C_ml` to the sea-cell tendency only.

## Section 1.3 — Tests

- q-flux: with a constant `ocean_heat_transport_convergence`, equilibrium SST shifts by exactly `Q_flux / (surface flux sensitivity)`; zero field reproduces current cached results.
- Ekman: an idealised zonal wind stress with nonzero curl (e.g. `τ_x = τ_0 cos(lat)` over a β/f-plane) reproduces the analytic wind-stress-curl Ekman pumping `w_Ek = k̂·∇×(τ/f)` (sign and magnitude), and the resulting `Q_ekman` heats/cools the mixed layer with the sign expected from Ekman convergence/divergence at the trade-wind and westerly latitudes. Zero-curl uniform stress gives zero pumping.
- Backward compat: `include_ekman=False`, no q-flux ⇒ bit-for-bit with current `SlabSurface`.

---

# Part 2 — DataOcean (AMIP-style prescribed SST)

## Section 2.1 — Type and role

`DataOcean` is a **`DiagnosticComponent`** — it prescribes, it does not prognose. Each step it reads the current model time from the state, interpolates an observed SST dataset to that time and the model grid, and **returns** `surface_temperature` and `sea_surface_temperature` over sea cells as **diagnostics** (plus surface fluxes). The set of sea cells comes from `area_type` — typically produced by `LandMask` (Part 4) for Earth runs. No timestep, no prognostic ownership of SST. This is the correct sympl idiom for AMIP forcing and composes cleanly: the atmosphere sees prescribed SSTs; sea-ice-covered cells are owned by `SeaIce` (Part 3) via `area_type`.

```python
DataOcean(
    sst_dataset,                       # path to netCDF, or opened xarray.Dataset
    sst_variable="tos",                # variable name in the dataset
    time_interpolation="mid_month",    # "mid_month" (AMIP, Taylor et al. 2000) | "linear"
    spatial_interpolation="bilinear",
    relaxation_timescale=None,         # None = hard prescription; else nudge toward obs
    compute_fluxes=True,               # emit surface SHF/LHF/albedo diagnostics
)
```

## Section 2.2 — Data ingestion

- **Reading**: open via xarray (netCDF). climt already depends on xarray/sympl, so no new hard dependency; SST files are optional at import.
- **Expected schema** (documented; HadISST / AMIP-II-like): coordinates `time`, `lat`, `lon`; SST variable in K or °C (units read from attrs, converted to K). Land/missing values masked.
- **Time interpolation**: default **mid-month** — monthly-mean SSTs are treated as mid-month values and interpolated with the AMIP-II correction (Taylor, Williamson & Zwiers, 2000) so monthly means are preserved. `"linear"` offered as a simpler option.
- **Spatial interpolation (resolution-agnostic, all sea cells)**: the interpolation target is **every model column with `area_type == "sea"`** — no other cell type is touched. Weights are built from the **model grid's own latitude/longitude coordinates read from the state**, not from any assumed resolution, so DataOcean works at arbitrary model resolution (coarser *or* finer than the source SST grid) and on regular or irregular grids. Bilinear interpolation from the observed source grid to each sea cell's (lat, lon).
  - **Coastline / masking guarantee**: because the model land–sea mask (`area_type`) generally does not align with the source dataset's mask, some sea cells can fall over source points that are land/missing. Every such cell is filled from the **nearest valid (unmasked) source SST** (a masked-array fill / nearest-neighbour extrapolation step) so **no `area_type=="sea"` cell is ever left masked or NaN**. This is asserted in the component.
  - Precomputed interpolation + fill weights are cached on first call (the model grid and `area_type` are fixed for a run) and reused every step.
- **Caching**: the dataset handle and the bracketing time slices are cached between steps; only re-read when the model time crosses into a new bracket.

## Section 2.3 — Outputs

Diagnostics: `sea_surface_temperature`, `surface_temperature` (sea cells set to obs; land/ice cells passed through unchanged), and when `compute_fluxes`: `surface_upward_sensible_heat_flux`, `surface_upward_latent_heat_flux`, `surface_albedo_for_direct_shortwave`, `surface_albedo_for_diffuse_shortwave` computed with the same bulk formulae used elsewhere (shared flux helper — see Part 4). With `relaxation_timescale` set, it emits a nudging tendency instead of a hard overwrite (then it is an `ImplicitTendencyComponent` variant; default hard-prescription stays a `DiagnosticComponent`).

## Section 2.4 — Tests

- Reads a small synthetic netCDF, prescribes SST at a known time, asserts sea cells match the (time-interpolated) obs and land/ice cells are untouched.
- **Arbitrary resolution**: the same source dataset is interpolated onto model grids both coarser and finer than the source; asserts every `area_type=="sea"` cell receives a finite SST and results are consistent (interpolated values bounded by local source neighbourhood).
- **Full sea coverage**: with a model land–sea mask deliberately offset from the source mask (sea cells over source-land points), assert no sea cell is NaN/masked — the nearest-valid fill covers all of them.
- Mid-month interpolation preserves monthly means to tolerance over a synthetic annual cycle.
- Spatial interpolation weight caching returns identical results across steps and is invariant to repeated calls.

---

# Part 3 — SeaIce + LandIce over a shared column solver

## Section 3.1 — Shared core: `_snow_ice_column`

Extract the 1-D implicit snow/ice heat-diffusion solver (`IceSheet.calculate_new_ice_temperature`) into a reusable core (`_core/snow_ice_column.py`), parameterised by its **boundary conditions** rather than branching on `soil_temperature is None` (defect 8). Contract:
```
solve_column(rho, c, kappa, T_profile, dt, dz,
             top_bc,      # Dirichlet(melting) when melting, else Neumann(net_flux)
             bottom_bc)   # BoundaryCondition object: Dirichlet(T) | Flux(Q)
    -> new T_profile
```
`top_bc` / `bottom_bc` are small value objects so `SeaIce` and `LandIce` (and SecondBEST's `SubsurfaceTransport`, see land-surface doc) supply their own boundaries without the core knowing which caller it serves. This is where the land-surface doc's recommendation to share one tridiagonal core lands. Vectorising the per-column loop / batching the solve (defect 7) is addressed here — the core takes column-stacked inputs.

## Section 3.2 — `SeaIce`

`Stepper` owning `sea_ice` cells:
- **Top BC**: surface energy balance (melting-Dirichlet when at freezing, else net-flux Neumann) — as today.
- **Bottom BC (fix defect 2)**: an **ocean heat flux** into the ice base. Basal growth/melt driven by `F_base = Q_ocean − k(T_ice_base − T_freeze)/dz`, where `Q_ocean` comes from the mixed layer (input `heat_flux_into_sea_water_due_to_sea_ice` / mixed-layer temperature) instead of a bare freezing-point Dirichlet. Couples to the slab ocean via the existing `heat_flux_into_sea_water_due_to_sea_ice` diagnostic.
- **Thermodynamics**: growth from basal freezing, surface + basal melt, snow-on-ice; **clamp thickness ≥ 0** (fix defect 6) with leftover melt energy passed to the ocean.
- **Albedo (fix defect 3)**: configurable via `__init__` (defaults preserving 0.8/0.5/0.2 as documented constants), optionally reading a snow-age / temperature-dependent albedo; exposed as the standard albedo diagnostics.
- Remove `print`/`assert False` (defect 4); replace with a clamped, logged handling path.

## Section 3.3 — `LandIce`

`Stepper` owning `land_ice` cells (and snow-on-`land`):
- **Bottom BC**: soil heat flux (Dirichlet at `soil_surface_temperature`), as the land branch does today.
- **Mass balance**: accumulation (snowfall) − surface melt/runoff, tracking `land_ice_thickness` and `surface_snow_thickness`; no basal ocean flux.
- Shares `_snow_ice_column`, the configurable albedo, and thickness clamping with `SeaIce`.

## Section 3.4 — Migration & compatibility

- `IceSheet` is **kept as a thin deprecated shim** that dispatches to `SeaIce` / `LandIce` internally (or is retained unchanged for one release with a deprecation note) so existing scripts don't break. New code uses the split components.
- Snow-layer indexing (defect 5) revisited in the shared core: keep the grid-quantised split initially but document it and expose `num_layers`; a finer sub-grid snow-line is a follow-up, not blocking.

## Section 3.5 — Tests

- Column solver: unit tests vs analytic steady-state conduction for Dirichlet and Flux bottom BCs; batched vs looped agreement.
- SeaIce: basal growth under a cold ocean flux, surface melt under strong insolation, thickness never negative, energy closure.
- LandIce: accumulation/ablation mass balance closes; land bottom flux matches the old `IceSheet` land path (regression).
- Defect regressions: assert no dead branch, no `assert False` path, albedo configurable.

---

# Part 4 — LandMask: default Earth land/sea/ice configuration

Realistic runs need every column tagged `land` / `sea` / `land_ice` — but climt's state default is `area_type = "sea"` everywhere (`_core/initialization.py:885`), fine for single columns and aquaplanets, useless for Earth geography. Rather than baking geography into state init (which would break the aquaplanet default and single-column use), provide a **separate `DiagnosticComponent`** that *sets* `area_type`.

### Section 4.1 — Type and role

`LandMask` is a `DiagnosticComponent` that reads the model grid's `latitude`/`longitude` from the state and **returns `area_type`** as a diagnostic, filled with the present-day Earth configuration. It is the natural companion to `DataOcean`: same pattern (bundled reference data → interpolate to the model grid → diagnostic), and `DataOcean` / the land / ice components all key off the `area_type` it produces.

```python
LandMask(
    mask_dataset=None,        # None = bundled default Earth mask; or path / xarray.Dataset
    include_land_ice=True,    # tag Greenland & Antarctica as land_ice
)
```

### Section 4.2 — Data and interpolation (arbitrary resolution)

- **Categories set**: `land`, `sea`, and (when `include_land_ice`) `land_ice` for the permanent ice sheets (Greenland, Antarctica). **`sea_ice` is *not* set here** — sea ice is dynamic and is owned by the `SeaIce` component (Part 3), which promotes `sea` → `sea_ice` as ice forms. `LandMask` sets the *static* geography only.
- **Arbitrary resolution**: `area_type` is categorical, so interpolation from the reference geography to the model grid is **nearest-neighbour** (never bilinear), driven by the model grid's own lat/lon read from the state. Works at any resolution, regular or Gaussian, coarser or finer than the source — mirroring DataOcean's resolution-agnostic design.
- **Bundled default**: ship a compact low-resolution land/sea/ice-sheet raster (kept small in the spirit of the recent `chore/wheel-slim` work); users can override with `mask_dataset`. Nearest-neighbour upsampling means even a coarse bundled mask yields a valid mask at high model resolution (with the expected blockiness — documented; a finer bundled source is a follow-up).

### Section 4.3 — Composition and ordering

`LandMask` runs **first** in a run's diagnostic set so downstream components (`DataOcean`, `SlabSurface`, `SeaIce`, `LandIce`, BucketHydrology, SecondBEST) see the correct `area_type`. Because it is a `DiagnosticComponent` returning `area_type`, it can also be evaluated once at setup to stamp the initial state. It composes with `DataOcean` (which then prescribes SST on exactly the `sea` cells this produces) and with the ice components (which further tag `sea_ice`).

### Section 4.4 — Tests

- Bundled mask on a coarse grid reproduces expected gross geography (e.g. a mid-Pacific cell is `sea`, a central-Asia cell is `land`, an Antarctic cell is `land_ice`).
- **Arbitrary resolution**: the same source mask stamped onto coarse and fine model grids yields only valid category strings, every cell assigned, and land fraction converging toward the source as resolution increases.
- Composition: `LandMask` → `DataOcean` prescribes SST on precisely the `sea` cells and leaves `land`/`land_ice` untouched.

---

# Part 5 — Shared surface-flux helper (cross-cutting)

`SlabSurface`, `DataOcean`, BucketHydrology and SecondBEST all compute bulk SHF/LHF. Factor the bulk aerodynamic formulae into one helper (`_core/surface_fluxes.py`) so the ocean/data-ocean paths and the land paths use a single, tested implementation. This is a light refactor, done when `DataOcean` needs it, not a big-bang change.

---

## Build order

1. `LandMask` (`DiagnosticComponent` setting `area_type`) + bundled default mask — foundational; everything else keys off `area_type`.
2. `_snow_ice_column` shared core (unblocks the ice split and the land-surface `SubsurfaceTransport`).
3. Ice split: `SeaIce`, `LandIce`, `IceSheet` shim + defect fixes.
4. `SlabSurface` q-flux + Ekman.
5. Shared surface-flux helper (Part 5), then `DataOcean` (depends on it and on `LandMask`).

## Open questions / decisions deferred to implementation

- Ekman divergence on climt's actual horizontal grid (spherical spacing, poles) — pick the discretisation during implementation; single-column reduces to zero.
- Whether `IceSheet` becomes a hard shim or is retained one release unchanged — decide with maintainers at implementation time.
- Exact observed-SST default dataset shipped (if any) vs user-supplied — likely user-supplied to keep the wheel slim (cf. recent `chore/wheel-slim` work); document the expected schema instead.
