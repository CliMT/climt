# Design: `SimpleBoundaryLayer` component + canonical jit tridiagonal solver

**Date:** 2026-07-22
**Branch:** `feature/land-ocean-ice-components` (or a fresh branch off `develop`)

## Motivation

Two coupled goals:

1. **Add the Frierson (2006) boundary-layer scheme** as a first-class climt
   component (`SimpleBoundaryLayer`). Source is a working but off-style
   prototype (raw `numba`, `parallel=True`, nested helper functions, a
   hand-rolled Thomas solver).

2. **Establish one canonical, numba-jitable, JAX-portable tridiagonal solver
   in `_core`** and remove `scipy` from every diffusion *compute* path.
   climt's strategic direction is that *all* components are numba-jitable now
   and JAX-portable later; `scipy.sparse`/`spsolve` cannot be `njit`-ed or
   JAX-traced, and the Python-level per-column loops around it are a
   performance problem (no `prange`).

## Reference

The scheme is the boundary layer from the appendix of:

> Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A Gray-Radiation
> Aquaplanet Moist GCM. Part I: Static Stability and Eddy Scale.
> *J. Atmos. Sci.*, **63**(10), 2548–2566. doi:10.1175/JAS3753.1

It is the companion to climt's existing `Frierson06LongwaveOpticalDepth`:
simplified Monin–Obukhov bulk surface exchange, a K-profile eddy diffusivity
capped by a critical Richardson number, and a diagnosed boundary-layer height.

**Important physics modification (from the thesis, not canonical Frierson).**
The concrete implementation follows Abel's MSc thesis
(`Abel_Master_s_thesis_signed_JMM.pdf`, §2.1.1), which modifies the
surface-layer diffusion coefficient `Kb`. Canonical Frierson (Eqn 2.6) uses
the *local* Richardson number `Ri(z)` in the `Kb` multiplier, which is
discontinuous at `Ria = 0`. The thesis (Eqn 2.8) replaces it with the
*surface-layer* Richardson number `Ria`, giving a multiplier that is
continuous at `Ria = 0` and decays from 1 → 0 as `Ria → Ric`:

```
Kb ∝ [ 1 + (Ria/Ric) · ln(z/z0) / (1 - Ria/Ric) ]^-1     for Ria > 0
Kb ∝ 1                                                    for Ria <= 0
```

**The supplied prototype used the un-modified local-`Ri` form (canonical 2.6);
the climt component must implement the thesis `Ria` form (2.8).** So
`_richardson_diffusivity` takes `Ri_a` (not the local `Ri`) in the denominator.
The local Richardson profile is still computed to diagnose the boundary-layer
top, just not used in the `Kb` multiplier.

## Scope

### In scope
- `_core/tridiagonal.py`: canonical `solve_tridiagonal` (Thomas), `@jit_compile`.
- `SimpleBoundaryLayer` component (greenfield), using that solver.
- Migrate the two `scipy`-in-solver call sites onto it, with `prange` column
  loops:
  - `_core/snow_ice_column.py` (consumed by `SeaIce` and `LandIce`)
  - `_components/second_best/processes/subsurface.py`
- Tests: solver unit tests; component standard suite + I/O-property +
  conservation tests; regression re-verification of the migrated components.
- Docs: a detailed `SimpleBoundaryLayer` manual + user guide; `references.bib`.

### Explicitly out of scope
- **`data_ocean` cyclic solve** (`_sst_interpolation.mid_month_values`): uses
  `np.linalg.solve` on a 12×12 matrix once at setup — **no scipy, cold path,
  not a perf concern.** Left as-is; a jitable cyclic (Sherman–Morrison)
  variant is a possible future add, not needed now.
- Non-solver scipy uses (`initialization.CubicSpline`,
  `data_ocean` `RegularGridInterpolator`/`cKDTree`, `cork` `netcdf_file`).
  These are setup/IO, not diffusion compute loops.
- A working JAX backend. We only ensure the solver is *written* to be
  JAX-swappable (pure function, no input mutation), with the swap point
  documented.

## Architecture

### 1. Canonical solver — `climt/_core/tridiagonal.py`

```python
from .backend import jit_compile

@jit_compile
def solve_tridiagonal(lower, diag, upper, rhs):
    """Solve a tridiagonal system A x = rhs by the Thomas algorithm.

    lower : sub-diagonal, length n-1  (A[i, i-1])
    diag  : main diagonal, length n
    upper : super-diagonal, length n-1 (A[i, i+1])
    rhs   : right-hand side, length n
    Returns x, length n. Pure function: inputs are not mutated.
    """
```

- Straight port of the prototype's `TDMAsolver`, promoted to `_core`,
  documented, `@jit_compile`'d so it is callable from inside other
  `@jit_compile` kernels (numba njit → njit calls are fine).
- Re-exported from `climt/_core/__init__.py`.
- **JAX portability note (in the docstring/comment):** the forward/backward
  sweeps are loop-carried; a JAX implementation would express them via
  `lax.scan` (or parallel cyclic reduction) behind the same
  `jit_compile(backend=...)` hook. Keeping the function pure (allocating a new
  output, never writing into `rhs`/`diag`) is what makes that swap mechanical.
- Taxonomy comment documenting the three tridiagonal families in climt so no
  one hand-rolls a fifth: (a) this dense Thomas solver — default; (b) cyclic
  (data_ocean) — needs a periodic variant; (c) scipy-sparse — being removed.

### 2. `SimpleBoundaryLayer` — `climt/_components/simple_boundary_layer/`

`component.py` + `__init__.py`, class `SimpleBoundaryLayer(Stepper)`.
Rewrite of the prototype to climt conventions:

- Kernel `_boundary_layer_kernel` decorated `@jit_compile`; the column loop
  uses `prange` (both from `_core.backend`). No raw `numba`/`parallel=True`
  import, so the numba-optional fallback (pure Python when numba absent) works
  like every other component.
- Prototype nested helpers become module-level `@jit_compile` functions:
  - `_richardson_diffusivity` (the `K_b` profile function) — **modified to use
    the surface-layer `Ri_a` per thesis Eqn 2.8, not the prototype's local
    `Ri`** (see Reference section).
  - `_diffuse_profile` (assembles the implicit diffusion matrix, then calls
    `solve_tridiagonal`)
- Physics constants bundled in a `NamedTuple` (`BoundaryLayerParams`) passed
  into the kernel, matching `DryConvectiveAdjustment`, instead of ~8 loose
  scalar args.
- Clean the prototype's `'degK '` trailing-space units to `'degK'`.

**Properties** (unchanged from the prototype; all inputs verified to have
registered defaults in `initialization.py`):

- inputs: `air_temperature`, `specific_humidity`, `air_pressure`,
  `air_pressure_on_interface_levels`, `northward_wind`, `eastward_wind`,
  `surface_air_pressure`, `surface_temperature`, `surface_specific_humidity`
- outputs: `air_temperature`, `specific_humidity`, `northward_wind`,
  `eastward_wind`
- diagnostics: `northward_wind_stress`, `eastward_wind_stress`,
  `boundary_layer_height`

`__init__` kwargs (from the prototype): `von_karman_constant=0.4`,
`roughness_length=3.21e-5`, `specific_fraction=0.1`,
`reference_pressure=100000`, `critical_richardson_number=1`.

**Registration:** add `SimpleBoundaryLayer` to `climt/_components/__init__.py`
(import + `__all__`) and `climt/__init__.py` (import + `__all__`).

### 3. Migrate `snow_ice_column` off scipy (affects `SeaIce`, `LandIce`)

`solve_column` currently assembles two `scipy.sparse` matrices, does a sparse
matvec for the RHS, folds boundary conditions into the diagonals via
`isinstance(bc, Dirichlet/Flux)`, and calls `spsolve`.

- Replace the sparse LHS + `spsolve` with `solve_tridiagonal` on the
  fold-in diagonal arrays (`a_sub_lhs`, `dp_lhs`, `a_sup_lhs`). Replace the
  sparse RHS matvec (`mat_rhs * temperature`) with an explicit jitable
  tridiagonal matvec.
- `isinstance` on the `Dirichlet`/`Flux` dataclasses cannot run under `njit`.
  Split into: a thin Python wrapper `solve_column(...)` that translates the BC
  dataclasses to `(bc_type_flag: int, value: float)` scalars, and a
  `@jit_compile` core `_solve_column_kernel(rho, c, kappa, T, dt, dz,
  top_flag, top_val, bot_flag, bot_val)` that does the assembly + solve.
- **`prange`:** port each consumer's per-column body into a `@jit_compile`
  kernel with a `prange` column loop that calls `_solve_column_kernel`. For
  `SeaIce`/`LandIce` this means moving the melt-check / re-solve / clamp /
  flux-accounting arithmetic into the kernel and returning what is currently
  logged as an explicit array. Where a branch is genuinely non-numeric
  (debug logging), drop it or surface it as a returned diagnostic.

### 4. Migrate `second_best` subsurface off scipy

`BestSubsurfaceTransport.__call__` builds a `scipy.sparse` LIL matrix with
Neumann top/bottom rows, calls `spsolve`, then applies an explicit freeze/melt
source. Second_best's column loop (`component.py:93`) is Python-level.

- Replace matrix assembly + `spsolve` with diagonal arrays + `solve_tridiagonal`
  (Neumann rows folded in exactly as today).
- Keep the freeze/melt source arithmetic; make the whole per-column transport a
  `@jit_compile` kernel and drive the column loop with `prange`.
- Preserve `second_best`'s independence from the ice core (do **not** re-couple
  it to `snow_ice_column`); it simply calls the shared `solve_tridiagonal`.

## Testing

### Solver unit tests — `tests/test_tridiagonal.py`
- `solve_tridiagonal` vs `numpy.linalg.solve` on randomly generated
  diagonally-dominant tridiagonal systems (assert `np.allclose`, tight tol).
- Edge cases: n=1, n=2, a strongly diagonally-dominant system, a
  non-symmetric system.

### `SimpleBoundaryLayer` tests — `tests/test_components.py`
- `TestSimpleBoundaryLayer(ComponentBaseColumn, ComponentBase3D)` → inherits
  cached-output, no-NaN, transpose/reverse-invariance. Generates a new cache
  entry on first run (committed).
- **I/O property test:** assert input/output/diagnostic names, dims, and units
  match the spec above.
- **Conservation test:** the implicit no-flux scheme conserves the
  interface-pressure-weighted column integral. (Verified analytically: the
  off-diagonal flux terms are antisymmetric —
  `w_i * diag_p[i] == w_{i+1} * diag_m[i+1]` with `w_i = p_int_i - p_int_{i+1}`
  — and `diag_m[0] = diag_p[n-1] = 0`, so there is no flux through either
  boundary; the surface flux is assumed already applied by a separate
  component.) Test asserts mass-weighted ∫T dp, ∫q dp, ∫u dp, ∫v dp are each
  conserved across one step to tolerance.

### Re-verification for the migrated components
`SeaIce`, `LandIce`, and `SecondBEST` are **all new on this branch** — there
are no published/reference values that must be preserved, so a bit-level
regression is not a concern. After migrating each off scipy:
1. Run its existing `ComponentBase` suite; the Thomas solve and `spsolve`
   should agree to ~machine precision on these well-conditioned systems, so
   the cached-output tests will most likely pass unchanged.
2. If values drift at round-off level, simply regenerate and commit the cache
   after a quick physical sanity check (no NaNs, energy/thickness behaves
   sensibly).

## Docs
- **No API auto-listing entry** (auto-doc has had issues previously).
- `docs/user-guide/simple_boundary_layer.qmd`: a detailed manual — the physics
  (Monin–Obukhov surface exchange, K-profile, critical Ri, boundary-layer
  height), the assumption that a surface-flux component ran first, the inputs/
  outputs/diagnostics, the `__init__` parameters, and a worked usage snippet,
  citing Frierson et al. (2006).
- Register the manual page in `docs/_quarto.yml` **user-guide section only**
  (matching the existing `dry_convective_adjustment.qmd` / `bucket_hydrology.qmd`
  entries), not the API listing.
- Add the Frierson (2006) Part I entry to `references.bib`.

## Post-change housekeeping
Per `CLAUDE.md`: run `graphify update .` then `python scripts/augment_graph.py`
to refresh the knowledge graph after the code changes land.

## Risks
- These components are all new on this branch, so there is no reproducibility
  risk from shifting numerics — caches can be regenerated freely. The main
  cost is engineering effort, not regression exposure.
- Fully `njit`-ing the sea-ice per-column body (melt logic, logging) may
  require restructuring control flow that numba dislikes; fallback is to jit
  the numeric core and keep a thinner Python loop if a clean prange proves
  impractical — but the goal is a full prange column loop.
