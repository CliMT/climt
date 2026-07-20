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

Three workstreams:

1. **SlabSurface** gains an additive ocean heat-transport term over sea = prescribed q-flux (arbitrary forcing) + optional diagnosed Ekman convergence. Backward compatible (both zero by default).
2. **DataOcean** — new `DiagnosticComponent` that reads an observed SST dataset and prescribes sea-cell `surface_temperature` / `sea_surface_temperature` as diagnostics, plus surface fluxes.
3. **Ice split** — `SeaIce` + `LandIce` over a shared `_snow_ice_column` implicit solver, fixing defects 1–8.

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

**Ekman transport (`Q_ekman`, optional)** — wind-driven horizontal heat transport diagnosed from surface wind stress:
```
Ekman mass transport:  M = (τ_y, −τ_x) / f          (per unit length, f = Coriolis)
Q_ekman = −ρ_sw c_sw · ∇·( M · T_ml )                (convergence of Ekman heat transport)
```
Enabled by `include_ekman=True`. Requires surface wind stress (`surface_downward_eastward_stress`, `surface_downward_northward_stress`) and latitude/`coriolis_parameter`. Near the equator `f → 0`; apply the conventional cap (bound `1/f`, or blend to zero within a few degrees of the equator) — documented. The divergence needs horizontal grid geometry; computed with climt's existing lat/lon spacing, falling back to zero (single-column) when the grid is degenerate.

## Section 1.2 — Component surface

Keep `SlabSurface` a `TendencyComponent` with the numba kernel. Changes:
- `__init__(include_ekman=False, equatorial_ekman_cap_latitude=5.0, **kwargs)`.
- New optional inputs: `ocean_heat_transport_convergence` (default 0 via state init), and — only when `include_ekman` — the two stress components and `coriolis_parameter`.
- New diagnostics: `ocean_heat_transport_convergence` (echo of total applied), `ekman_heat_transport_convergence` (the diagnosed part), so runs can separate prescribed vs diagnosed contributions.
- Kernel adds `Q_transport / C_ml` to the sea-cell tendency only.

## Section 1.3 — Tests

- q-flux: with a constant `ocean_heat_transport_convergence`, equilibrium SST shifts by exactly `Q_flux / (surface flux sensitivity)`; zero field reproduces current cached results.
- Ekman: idealised zonal wind stress over an f-plane produces the analytic Ekman convergence sign (upwelling/downwelling) at the expected latitudes.
- Backward compat: `include_ekman=False`, no q-flux ⇒ bit-for-bit with current `SlabSurface`.

---

# Part 2 — DataOcean (AMIP-style prescribed SST)

## Section 2.1 — Type and role

`DataOcean` is a **`DiagnosticComponent`** — it prescribes, it does not prognose. Each step it reads the current model time from the state, interpolates an observed SST dataset to that time and the model grid, and **returns** `surface_temperature` and `sea_surface_temperature` over sea cells as **diagnostics** (plus surface fluxes). No timestep, no prognostic ownership of SST. This is the correct sympl idiom for AMIP forcing and composes cleanly: the atmosphere sees prescribed SSTs; sea-ice-covered cells are owned by `SeaIce` (Part 3) via `area_type`.

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
- **Spatial interpolation**: bilinear from the dataset grid to climt's grid; nearest-neighbour fallback at coastlines. Precomputed weights cached on first call (grid is fixed).
- **Caching**: the dataset handle and the bracketing time slices are cached between steps; only re-read when the model time crosses into a new bracket.

## Section 2.3 — Outputs

Diagnostics: `sea_surface_temperature`, `surface_temperature` (sea cells set to obs; land/ice cells passed through unchanged), and when `compute_fluxes`: `surface_upward_sensible_heat_flux`, `surface_upward_latent_heat_flux`, `surface_albedo_for_direct_shortwave`, `surface_albedo_for_diffuse_shortwave` computed with the same bulk formulae used elsewhere (shared flux helper — see Part 4). With `relaxation_timescale` set, it emits a nudging tendency instead of a hard overwrite (then it is an `ImplicitTendencyComponent` variant; default hard-prescription stays a `DiagnosticComponent`).

## Section 2.4 — Tests

- Reads a small synthetic netCDF, prescribes SST at a known time, asserts sea cells match the (time-interpolated) obs and land/ice cells are untouched.
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

## Build order

1. `_snow_ice_column` shared core (unblocks the ice split and the land-surface `SubsurfaceTransport`).
2. Ice split: `SeaIce`, `LandIce`, `IceSheet` shim + defect fixes.
3. `SlabSurface` q-flux + Ekman.
4. `DataOcean` (depends on a shared surface-flux helper — Part 4 below).

## Part 4 — Shared surface-flux helper (cross-cutting)

`SlabSurface`, `DataOcean`, BucketHydrology and SecondBEST all compute bulk SHF/LHF. Factor the bulk aerodynamic formulae into one helper (`_core/surface_fluxes.py`) so the ocean/data-ocean paths and the land paths use a single, tested implementation. This is a light refactor, done when `DataOcean` needs it, not a big-bang change.

## Open questions / decisions deferred to implementation

- Ekman divergence on climt's actual horizontal grid (spherical spacing, poles) — pick the discretisation during implementation; single-column reduces to zero.
- Whether `IceSheet` becomes a hard shim or is retained one release unchanged — decide with maintainers at implementation time.
- Exact observed-SST default dataset shipped (if any) vs user-supplied — likely user-supplied to keep the wheel slim (cf. recent `chore/wheel-slim` work); document the expected schema instead.
