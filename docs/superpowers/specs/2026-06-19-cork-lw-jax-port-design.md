# CORK Longwave — Differentiable JAX Port Design

**Date:** 2026-06-19
**Branch:** `jax-research` (permanent differentiable-port branch — never merged to `develop`/main)
**Status:** Design approved; pending implementation plan

## Goal

Provide a fully `jax.grad`-differentiable, `jax.jit`- and `jax.vmap`-compatible JAX
implementation of CORK's shared **correlated-k optics** layer and its **longwave (LW)**
component, validated numerically against the existing Numba kernels.

On this branch the JAX path is the **primary** path for researchers; the Numba/`@njit`
code is retained **only as a validation oracle** (parity reference), not for production
use here. This spec covers **optics + longwave**. CORK shortwave gets its own spec later
and will reuse the proven differentiable optics layer.

### Context: where this fits

This is the next component port in the `jax-research` roadmap. JAX versions of the GFS
dynamical core and EmanuelConvection already exist; the LW radiation port plus a later
shortwave port and a boundary-layer scheme are the remaining pieces needed for a full
**differentiable aquaplanet** simulation. The Emanuel differentiable port
(`docs/superpowers/plans/2026-04-04-emanuel-jax-differentiable.md`) is the template for
conventions: runtime `isinstance(jax.Array)` dispatch in `array_call`, JAX functional
kernels alongside the Numba path, `jax_enable_x64`, and phased parity-gated validation.

## Scope

**In scope**

- Differentiable JAX correlated-k optics: optical depth `tau` from the loaded k-table.
- Differentiable JAX longwave: Planck sources, radiative transport, heating rate.
- JAX-first dispatch in `CorkLongwaveRadiation` with the Numba path retained as oracle.
- Numerical parity, transform-safety (jit/vmap), and gradient-correctness tests.

**Out of scope** (named for later specs)

- CORK shortwave (own spec/plan).
- The boundary-layer scheme and the aquaplanet driver.
- Gradients w.r.t. scheme parameters: the k-tables, continuum, or Parmentier optical
  coefficients. Gradients are required **w.r.t. atmospheric state only**.

## Decisions (from brainstorming)

| Decision | Choice |
| --- | --- |
| Differentiability target | Full `jax.grad` (JIT + vmap + reverse-mode), like the Emanuel port |
| Spec scope | Optics + longwave first; shortwave later |
| Gradient inputs | Atmospheric state only (T, p, composition, T_surface); tables static |
| Code organization | **Approach B** — separate JAX kernel modules; leave Numba kernels untouched |
| Backend role on this branch | JAX = primary/default; Numba = validation oracle only |

## Architecture & module layout

```
climt/_components/cork/
  optics/
    correlated_k.py        # UNTOUCHED — Numba oracle + table loading (host I/O)
    correlated_k_jax.py    # NEW — differentiable tau via gather + multilinear blend
  lw/
    component.py           # MODIFIED — array_call dispatches JAX-first
    kernels.py             # UNTOUCHED — Numba oracle (planck + transport)
    kernels_jax.py         # NEW — differentiable planck sources + lax.scan transport
```

Key boundaries:

- **Host vs. traced.** Table *loading* (HDF5/netCDF → numpy) stays host-side in
  `correlated_k.py` and runs once at component construction. The loaded k-table, axis
  grids, g-point weights, continuum, and band data are passed into the JAX kernels as
  **static `jnp` constants** (closed over / non-differentiated). Only the atmospheric
  state (T, p, composition, T_surface) is traced.
- **Pure functional core.** The new `*_jax.py` modules expose side-effect-free functions
  `(state arrays, static table data) -> fluxes/heating` with no in-place writes, so
  `jit`, `vmap` (over columns), and `grad` compose cleanly.
- **Dispatch.** `CorkLongwaveRadiation.array_call` selects the JAX kernels when inputs are
  `jax.Array` (and `HAS_JAX`), else Numba — same convention as Emanuel V3, but JAX is the
  documented default on this branch.

## Differentiable optics layer (`correlated_k_jax.py`)

The Numba kernel `_ck_tau_additive_co2_kernel` computes `tau[nband, ngpt, nlev, ncol]`
via scalar loops over `(ncol, nlev, nband, ngpt, ngas)` with `_ck_bracket` (searchsorted)
on four axes (T, log-p, log-x, log-c) and an 8-corner multilinear blend (`_ck_txx7`),
plus optional continuum and a `co2_logk` log-interpolation path.

Differentiable reformulation — separate the discrete from the smooth:

1. **Bracketing (discrete, non-differentiated).** For each of the 4 axes, compute the
   integer bracket index `i` and fractional weight `f` for every (level, column).
   Index via `jnp.searchsorted` + `jnp.clip` to `[0, n-2]`; weight
   `f = clip((v − grid[i]) / (grid[i+1] − grid[i]), 0, 1)`. The **index** carries no
   gradient (integer gather); the **weight `f`** is a smooth function of the state and is
   where the gradient flows. This matches the Emanuel `searchsorted`-bracket treatment.
2. **Gather (vectorized).** Replace the 8 scalar corner reads with vectorized `jnp`
   fancy-indexing of the static k-table over the 8 (T,P,X) corners (and the 2 C-corners),
   producing corner arrays shaped to broadcast over `(nband, ngpt, nlev, ncol)`. The table
   is a static `jnp` constant, so this is a pure gather.
3. **Smooth blend (differentiated).** Multilinear combination with `(1−f, f)` weights —
   pure arithmetic, fully differentiable. The `co2_logk` path uses `exp(log · blend)` with
   a `jnp.maximum(k, FLOOR)` floor (smooth a.e.); continuum uses `exp(_ck_txx_cont(...))`.
   Final `tau = Σ_gas kv · gas_amount + continuum`.

Gradient-correctness notes:

- The blend is **piecewise-linear** in the state (weights clipped to `[0,1]`, indices
  stepwise). Gradients are correct almost everywhere; the non-differentiable points are
  exactly at grid nodes / clip boundaries — acceptable and identical in character to the
  accepted Emanuel behavior.
- `FLOOR` and weight clips use `jnp.maximum` / `jnp.clip`, whose subgradients JAX handles;
  no `nan`-producing ops (no bare `log(0)` or `0/0`).
- **vmap over columns** is native since everything is expressed as array ops over the
  `(…, nlev, ncol)` trailing axes.

Interface (pure function):

```
tau = compute_tau_jax(T, log_p, log_x, log_c, gas_amounts,   # traced state
                      table)                                  # static table data
# tau: (nband, ngpt, nlev, ncol)
```

`table` is a frozen container (NamedTuple of `jnp` arrays + python flags like `co2_logk`,
`has_cont`) built once at component construction from the existing host-side loader.

## Longwave component (`lw/kernels_jax.py` + `lw/component.py`)

After `tau`, the LW path has three stages; all are smooth except the two `searchsorted`
brackets, handled as in the optics layer.

### 1. Planck sources (differentiable)

The Numba `planck_sources_kernel` interpolates a band/g-point Planck *fraction* linearly
in T (a `searchsorted` bracket on `T_grid`) and multiplies by `σT⁴`.

JAX: bracket `T` (and `T_surface`) on `T_grid` → index (gather) + weight `f`
(differentiated); `planck_frac` is a static `jnp` constant gathered at the 2 T-corners;
blend with `(1−f, f)`; multiply by `σ·T⁴`. `σT⁴` and the linear blend are smooth, so
`∂(planck_source)/∂T` and `∂(surface_source)/∂T_surface` flow cleanly. Outputs
`planck_source[nband,ngpt,nlev,ncol]`, `surface_source[nband,ngpt,ncol]`.

### 2. Transport via `lax.scan` (differentiable)

The Numba `_lw_transport_kernel` is two sequential recurrences per (band, g-point):
`trans = exp(−1.66 · tau)`, `flux = prev · trans + planck · (1 − trans)` — smooth, no
scattering.

JAX: express each sweep as a `jax.lax.scan` over the `nlev` axis (up: surface→TOA with BC
`emissivity · surface_source`; down: TOA→surface with BC 0). Vectorize the
`(nband, ngpt, ncol)` axes inside the scan body so a single scan carries all
bands/g-points/columns. `scan` is reverse-mode differentiable, so gradients flow through
the full vertical recurrence. Accumulate weighted fluxes
(`up_broad = Σ_b Σ_g weights[b,g] · up`, likewise down) via `jnp.einsum`/sum; the
bit-for-bit accumulation *order* of the oracle is irrelevant for the JAX path — only
parity to tolerance is required.

### 3. Heating rate

Net flux divergence → temperature tendency, same algebra as the existing component, in
`jnp`.

### Component wiring (`lw/component.py`)

`array_call` keeps the existing `input_properties` / `tendency_properties` /
`diagnostic_properties`. It builds the static `table` container once (host load) and, when
inputs are `jax.Array`, calls a single jitted `cork_lw_jax(state, table) -> (tendencies,
diagnostics)` composed of stages 1–3. The Numba branch is untouched. Diagnostics
(`transmittance`, per-g-point fluxes) are **optional and forward-only** — computed only
when requested and kept out of the differentiated core so they never complicate gradients.

Differentiability scope: full `jax.grad` of the heating rate / fluxes w.r.t. atmospheric
state (T, p, composition, T_surface); tables static. `jit` + `vmap`-over-columns supported
end to end.

## Runtime dispatch & public API

**No public API change.** `CorkLongwaveRadiation` keeps its name, constructor,
`input_properties`, `tendency_properties`, and `diagnostic_properties`. The backend is
selected by the *data*, not by a new class.

Construction (once): `__init__` loads the k-table/grids/weights/continuum host-side via the
existing `correlated_k.py` loader, then builds two views — the numpy arrays the Numba
oracle already uses, and a frozen `table` NamedTuple of `jnp` arrays + python flags for the
JAX path (materialized lazily on first JAX call and cached).

Dispatch in `array_call` (mirrors Emanuel V3):

```
use_jax = HAS_JAX and isinstance(<a primary input array>, jax.Array)
if use_jax:
    out = cork_lw_jax(state_arrays, self._jax_table)   # jitted, differentiable
else:
    out = <existing numba path>                          # oracle / fallback
```

The jitted entry point is cached so JIT compiles once per shape (table closed over / marked
static). On this branch the JAX path is the documented default: docstrings and any
example/driver pass `jax.Array` state; Numba stays reachable only for parity tests.

State plumbing: a thin adapter converts the sympl `state` dict arrays to the kernel's
expected `(nlev, ncol)` / `(ncol,)` layouts (matching the orientations the Numba kernels
already use) and converts kernel outputs back into the tendency/diagnostic dict. This
adapter is backend-agnostic (works on `np` or `jnp`), so the only branch point is which
compute function runs.

`jax_enable_x64`: the kernels require float64 (radiation fluxes need the precision; the
oracle is float64). The module sets `jax.config.update("jax_enable_x64", True)` on import,
exactly as `pure_python_v3.py` does — a hard requirement for parity.

## Validation strategy

The Numba kernels are the oracle. Every JAX stage is validated against them on the *same*
loaded table and representative profiles before we rely on it. Three test layers:

1. **Forward parity (numerical, not bit-for-bit).** JAX vs Numba on a standard Earth-LW
   table + a representative column set (`ncol ∈ {1, 4}`, moist/dry profiles):
   - `compute_tau_jax` vs `_ck_tau_additive_co2_kernel`
   - JAX Planck sources vs `planck_sources_kernel`
   - JAX transport vs `lw_transport` (up/down band & broadband fluxes)
   - full-component heating rate vs the Numba component.

   Tolerance: `rtol/atol ≈ 1e-6` for tau and fluxes (float64 both sides; differences come
   only from reduction order and the float32 `planck_frac` promotion the oracle documents).
   Bit-for-bit is explicitly **not** required.
2. **Transform-safety.** `jax.jit` and `jax.vmap`-over-columns produce identical results to
   the un-transformed JAX path; recompilation happens once per shape.
3. **Gradient correctness.** `jax.grad` of a scalar reduction of the heating rate w.r.t.
   `T`, `T_surface`, and a composition field, checked against central finite differences
   (loose tol, e.g. `rtol ≈ 1e-3`), plus an assertion that gradients are finite (no
   `nan`/`inf`) across moist/dry/cold columns — i.e. the floors and clips behave.

## Phasing & milestones

Each phase is independently reviewable and parity-gated.

- **Phase 1 — Differentiable optics.** `correlated_k_jax.py`: brackets, gather,
  multilinear/log/continuum blend → `compute_tau_jax`. Tests: tau parity + `∂tau/∂T`
  finite-difference. *Exit:* tau matches oracle and is differentiable.
- **Phase 2 — LW kernels.** `kernels_jax.py`: differentiable Planck sources + `lax.scan`
  transport + flux accumulation. Tests: Planck/flux parity, grad finite. *Exit:* fluxes
  match oracle, differentiable.
- **Phase 3 — Component integration.** `lw/component.py` dispatch + state adapter + cached
  jitted entry point + x64. Tests: full heating-rate parity, jit/vmap-safety, end-to-end
  `jax.grad`. *Exit:* `CorkLongwaveRadiation` runs JAX-first, differentiable, parity-clean.

## Success criteria

- `CorkLongwaveRadiation` on JAX inputs returns heating rates matching the Numba oracle
  within tolerance.
- The whole LW path is `jit` / `vmap` / `grad`-compatible, with verified finite,
  finite-difference-correct gradients w.r.t. atmospheric state.
- The Numba kernels remain untouched and serve as the parity oracle.
