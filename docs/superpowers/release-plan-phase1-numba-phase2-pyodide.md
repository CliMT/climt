# climt Release Plan

## Phase 1: Stabilize Numba Branch, Benchmark, Merge, Release

**Branch:** `numba-optimized-components` → merge to `develop`

### Context

The `numba-optimized-components` branch adds `@jit_compile`/`@njit` kernels to 8 pure Python components (HeldSuarez, GrayLongwaveRadiation, Frierson06, GridScaleCondensation, DryConvectiveAdjustment, BergerSolarInsolation, SlabSurface, Instellation). It also introduces the UnytBackend as a lightweight alternative to sympl's DataArrayBackend. Both need to be tested and benchmarked before merging.

### Step 1: Run all tests, fix failures

Run the full test suite in the `climt` conda env. Fix any failures.

```
conda activate climt
pytest tests/ -v
```

**Key test files:**
- `tests/test_*_optimization.py` (8 files) — parity tests for numba kernels
- `tests/test_components.py` — integration tests for all components
- `tests/test_initialization.py` — grid/state setup
- `tests/test_conservation.py` — conservation laws
- `tests/test_unyt_backend.py` — UnytBackend unit tests
- `tests/test_unyt_integration.py` — UnytBackend + component integration

### Step 2: Create combined benchmark (Numba x Backend)

**Goal:** A single benchmark script showing a 2x2 matrix: {Numba on/off} x {DataArrayBackend / UnytBackend}.

**New file:** `benchmarks/benchmark_numba_x_backend.py`

For each pure Python component, measure wall-clock time across 4 configurations:
1. **DataArrayBackend, no Numba** (`NUMBA_DISABLE_JIT=1`)
2. **DataArrayBackend, with Numba**
3. **UnytBackend, no Numba** (`NUMBA_DISABLE_JIT=1`)
4. **UnytBackend, with Numba**

This shows both the Numba speedup and the backend overhead independently. Use `NUMBA_DISABLE_JIT=1` env var to disable JIT without code changes. Run each config in a subprocess to ensure clean JIT state.

Output format:
```
Component               | DA+noJIT | DA+Numba | Unyt+noJIT | Unyt+Numba
HeldSuarez              |  2.34s   |  0.12s   |   1.85s    |   0.09s
GrayLongwaveRadiation   |  ...     |  ...     |   ...      |   ...
```

**Existing files to reference:**
- `benchmarks/benchmark_numba_speedups.py` — kernel-level Numba timing pattern
- `benchmarks/benchmark_all_components.py` — DataArray vs Unyt component-level pattern

### Step 3: Run benchmarks and record results

Run the new combined benchmark plus the existing kernel-level benchmark. Save output to `benchmarks/results/` for the release notes.

```
python benchmarks/benchmark_numba_x_backend.py > benchmarks/results/numba_x_backend.txt
python benchmarks/benchmark_numba_speedups.py > benchmarks/results/numba_kernel_speedups.txt
```

### Step 4: Fix issues found during testing/benchmarking

**Known issues to fix:**

#### 4a. Add pytest tests for Emanuel V3
V3 works correctly in isolation (max diff vs V1: `1.7e-4` temperature, `6e-7` humidity on moist profile; exact parity on dry). It has no pytest integration — only `tests/verify_v3_basic.py` (standalone script) which fails because it also imports V2.

- Add `tests/test_emanuel_v3_parity.py` mirroring the structure of `tests/test_emanuel_optimization.py` but comparing V3 (`EmanuelConvectionPythonV3`) against V1 (`EmanuelConvectionPythonComponent`) only. Parametrize over `ncol=[1,4]` and `moisture_type=['moist','dry','unstable']`. Set tolerances to `atol=1e-3, rtol=1e-3` (V3 has minor algorithmic differences from V1).
- Do NOT import V2 in this test (avoid triggering the V2 `set_item` Numba bug).

#### 4b. Fix flaky test_adding_constant
`tests/test_utils.py::test_adding_constant` fails when the full suite runs but passes in isolation — state leak from another test. Investigate and fix the ordering issue (likely a shared sympl constant being modified and not reset).

#### 4c. Re-run benchmark_numba_x_backend.py at 512 columns
The 4096-column run showed ~1x Unyt speedup for most components because kernel computation dominates at large grids. At 512 columns, framework overhead matters and Unyt shows ~2x. Add a 512-column run alongside the 4096 results:
```
python benchmarks/benchmark_numba_x_backend.py --ncol 512 --iters 50 > benchmarks/results/numba_x_backend_512col.txt
```

### Step 5: Merge to develop and release

1. Merge `numba-optimized-components` into `develop`
2. Tag a new version
3. Verify CI passes on develop

---

## Phase 2: Pyodide Pure Python Support (after Phase 1 release)

**Spec:** `docs/superpowers/specs/2026-03-29-pyodide-pure-python-support-design.md`
**Branch:** New branch off `develop` after Phase 1 merge.

### Step 1: Build system for pure Python wheel

**Goal:** `CLIMT_PURE_PYTHON=1 python -m build` produces a `py3-none-any` wheel.

**Files:**
- `setup.py`: Wrap `ext_modules` and `climt_build_ext` in `if not os.environ.get('CLIMT_PURE_PYTHON')` check.
- `pyproject.toml`: Make Cython a conditional build dependency.

### Step 2: Silent import degradation for Fortran components

Replace `logging.warning` + `print` with `_FORTRAN_AVAILABLE` flag and clear `ImportError` at instantiation.

**Files:** `simple_physics/component.py`, `rrtmg/lw/component.py`, `rrtmg/sw/component.py`, `emanuel/component.py`, `dcmip/component.py`

### Step 3: Clean up BergerSolarInsolation

Remove dead Cython import (class is already 100% pure Python).

### Step 4: Implement TDMA tridiagonal solver

**New file:** `climt/_core/tridiagonal.py` — Thomas algorithm with `@jit_compile` decorator.

**Tests:** `tests/test_tridiagonal.py` with:
- Known analytic solutions (finite-difference ODE discretization)
- Random diagonally-dominant systems at N = 5, 50, 500, 5000 vs `numpy.linalg.solve`
- Singular/near-singular detection
- N=1 edge case
- Comparison against `scipy.sparse.linalg.spsolve` for IceSheet's actual tridiagonal system
- Boundary condition patterns (Dirichlet and flux-style row modifications)

### Step 5: Refactor IceSheet to use TDMA

Replace `scipy.sparse` with `solve_tridiagonal` in `surface_ice.py`.

### Step 6: Replace CubicSpline with np.interp

In `climt/_core/initialization.py` for ozone profile interpolation.

### Step 7: Remove scipy from dependencies

Remove `scipy>=0.18.1` from `setup.py` `install_requires`.

### Step 8: Runtime discoverability

Add `climt.has_fortran_extensions()` and `climt.FORTRAN_AVAILABLE`.

### Step 9: CI pure Python wheel build

Add job to `.github/workflows/release_climt.yml`.

### Step 10: End-to-end validation

Build pure Python wheel, install without scipy/numba/Cython, verify silent import and component availability.
