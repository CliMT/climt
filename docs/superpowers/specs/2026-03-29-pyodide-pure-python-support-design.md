# Pyodide-Compatible Pure Python climt

**Date**: 2026-03-29
**Status**: Draft — **superseded in part by the 2026-07-19 revision below**

> **Revision 2026-07-19 (CORK-aware Phase 2).** This design predates the CORK
> correlated-k radiation scheme and the pure-Python Emanuel port, so its
> premise ("the pedagogical case only needs gray radiation") is now obsolete.
> The pure-Python **science** stack for a non-grey radiative-convective
> equilibrium (RCE) already exists and works. The remaining work is packaging
> and dependency plumbing. Read the [**CORK-aware Phase 2 revision**](#revision-2026-07-19--cork-aware-phase-2)
> at the end of this document; where it conflicts with the original text, the
> revision wins.

## Context

climt needs to run in Pyodide/JupyterLite for pedagogical purposes. The current package requires Fortran compilation for 5 of its 15+ components, making it incompatible with WebAssembly environments. The goal is to ship a single `climt` package (not a separate `climt_lite`) that works both as a full-featured package with Fortran extensions and as a pure Python package in Pyodide.

The pedagogical use case requires: gray radiation, surface components (SlabSurface), and Emanuel convection (pure Python variant) -- all of which are already implemented in pure Python.

## Phased Approach

### Phase 1: Stabilize and Merge Numba Branch

Merge `numba-optimized-components` into `develop` after testing. This establishes the Numba-optional JIT pattern (`@jit_compile` with graceful fallback) as the baseline.

### Phase 2: Pyodide Support (new branch off develop)

All changes below are Phase 2 work.

## Design

### 1. Two Wheels from One Source

Build two wheel types from the same source for each release:

- **Platform-specific wheels** (`climt-X.Y.Z-cpXY-cpXY-<platform>.whl`): Include compiled Fortran extensions. Built by CI on Linux/macOS as today.
- **Pure Python wheel** (`climt-X.Y.Z-py3-none-any.whl`): No compiled extensions. Built with `CLIMT_PURE_PYTHON=1`.

pip/micropip auto-selects the right one: platform-specific wheels are preferred on regular systems; Pyodide can only install `py3-none-any`.

**Files to modify:**
- `setup.py`: Check `CLIMT_PURE_PYTHON` env var. When set, skip `ext_modules` and Fortran compilation entirely. Ensure the wheel is tagged `py3-none-any`.
- `pyproject.toml`: Cython moves from build-requires to an optional/conditional dependency (only needed when building with extensions).
- `.github/workflows/release_climt.yml`: Add a job to build the pure Python wheel.

### 2. Silent Import Degradation

Replace noisy import-time warnings with silent flags. Clear errors only at instantiation.

**Current pattern (5 Fortran component files):**
```python
try:
    from . import _simple_physics as phys
except ImportError as error:
    logging.warning("Import failed...")
    print(error)
```

**New pattern:**
```python
_FORTRAN_AVAILABLE = False
try:
    from . import _simple_physics as phys
    _FORTRAN_AVAILABLE = True
except ImportError:
    pass

class SimplePhysics(Stepper):
    def __init__(self, ...):
        if not _FORTRAN_AVAILABLE:
            raise ImportError(
                "SimplePhysics requires compiled Fortran extensions. "
                "Install climt with Fortran support or use a pure Python alternative."
            )
        ...
```

**Files to modify:**
- `climt/_components/simple_physics/component.py`
- `climt/_components/rrtmg/lw/component.py`
- `climt/_components/rrtmg/sw/component.py`
- `climt/_components/emanuel/component.py`
- `climt/_components/dcmip/component.py`

### 3. BergerSolarInsolation: Remove Dead Cython Import

`berger_solar_insolation.py` already has a complete pure Python implementation. The Cython import (`from ._berger_solar_insolation import ...`) is unused -- the class uses `_get_orbital_parameters_functional` and `_get_solar_parameters_np`, both defined in the same file.

**Action:** Remove the try/except Cython import block. No functional change.

**File:** `climt/_components/berger_solar_insolation.py`

### 4. Remove scipy Dependency

scipy is used in exactly two places. Both are replaceable with pure numpy.

#### 4a. TDMA Tridiagonal Solver (shared utility)

Extract a reusable Thomas algorithm (TDMA) solver to replace `scipy.sparse.spdiags + spsolve` in IceSheet. This becomes a shared utility for IceSheet and upcoming diffusion components.

**New file:** `climt/_core/tridiagonal.py`

```python
@jit_compile(backend=np)
def solve_tridiagonal(a_lower, a_diag, a_upper, rhs):
    """Solve tridiagonal system Ax = rhs using Thomas algorithm.

    Args:
        a_lower: sub-diagonal (length n-1, first element unused)
        a_diag: main diagonal (length n)
        a_upper: super-diagonal (length n-1, last element unused)
        rhs: right-hand side vector (length n)

    Returns:
        x: solution vector (length n)
    """
    ...
```

Numba-JIT compatible for performance on regular systems; pure numpy fallback in Pyodide.

**Files to modify:**
- `climt/_core/tridiagonal.py` (new)
- `climt/_core/__init__.py` (export)
- `climt/_components/surface_ice.py` (replace scipy.sparse with TDMA calls)

#### 4b. Replace CubicSpline with numpy

`climt/_core/initialization.py` uses `scipy.interpolate.CubicSpline` for ozone profile interpolation onto model levels. Replace with `np.interp` (linear interpolation). The reference profile has 30 points -- linear interpolation is sufficient for this coarse initialization data. If higher fidelity is needed later, a pure-numpy cubic spline can be added to `climt/_core/`.

**File:** `climt/_core/initialization.py`

#### 4c. Update setup.py

Remove `scipy>=0.18.1` from `install_requires`.

**File:** `setup.py`

### 5. Runtime Discoverability

Add a simple mechanism for users/notebooks to check what's available:

**File:** `climt/__init__.py`

```python
def has_fortran_extensions():
    """Return True if compiled Fortran extensions are available."""
    try:
        from climt._components.simple_physics import _simple_physics  # noqa: F401
        return True
    except ImportError:
        return False
```

### 6. Dependency Adjustments for Pure Python Wheel

**Hard dependencies (work in Pyodide):**
- numpy
- xarray
- pint
- unyt
- sympl (pure Python)

**Removed:**
- scipy (replaced with numpy implementations)

**Build-only (not needed at runtime for pure Python wheel):**
- Cython (only for compiling extensions)

**Optional (graceful degradation):**
- numba (not available in Pyodide; JIT decorators already degrade to no-ops)

## Component Availability Summary

| Component | Pure Python Wheel | Full Wheel |
|-----------|------------------|------------|
| GrayLongwaveRadiation | Yes | Yes |
| Frierson06LongwaveOpticalDepth | Yes | Yes |
| SlabSurface | Yes | Yes |
| HeldSuarez | Yes | Yes |
| GridScaleCondensation | Yes | Yes |
| DryConvectiveAdjustment | Yes | Yes |
| Instellation | Yes | Yes |
| BucketHydrology | Yes | Yes |
| IceSheet | Yes | Yes |
| BergerSolarInsolation | Yes | Yes |
| EmanuelConvectionPython | Yes (WIP) | Yes |
| EmanuelConvectionPythonV3 | Yes (WIP) | Yes |
| EmanuelConvection (Fortran) | Error on init | Yes |
| SimplePhysics | Error on init | Yes |
| RRTMGLongwave | Error on init | Yes |
| RRTMGShortwave | Error on init | Yes |
| DcmipInitialConditions | Error on init | Yes |

## Verification

1. **Unit tests pass** with `CLIMT_PURE_PYTHON=1` (skip Fortran component tests)
2. **`import climt` produces no warnings** when Fortran extensions are absent
3. **Fortran components raise clear `ImportError`** on instantiation without extensions
4. **TDMA solver matches scipy results** for IceSheet test cases (validate before removing scipy)
5. **Pure Python wheel installs and runs** in a Pyodide environment (JupyterLite)
6. **Platform-specific wheel still works** identically to current behavior (regression test)
7. **Emanuel pure Python variants** tested against Fortran reference (ongoing, separate effort)

## Related: in-browser notebook target

The cork Experiments post
(`docs/experiments/2026-06-05-cork-co2-bands/`) is a target for the
website spec's deferred live-cells follow-up, which depends on this pure-Python
wheel. Its companion notebook already maps the boundary: the setup, the k(g) toy,
and all analysis/plotting of pre-baked arrays are pure-Python-live under Pyodide,
while the cells that *produce* the data (RRTMG Fortran, linepyline line-by-line)
read committed `_artifacts/*.npz` instead. When this wheel lands, that notebook is
ready to go live with no re-mapping.

---

## Revision 2026-07-19 — CORK-aware Phase 2

Since the original draft, two things changed the picture:

1. **CORK** (correlated-k non-grey radiation; `CorkLongwaveRadiation` /
   `CorkShortwaveRadiation`) landed. It is **100% pure Python** — no `.pyx`, no
   Fortran — with a numba-*optional* kernel layer (`climt/_components/cork/common.py`
   defines an `njit` no-op fallback when numba is absent).
2. The **pure-Python Emanuel** port (`EmanuelConvectionPython`, v3) is exported
   at the top level.

Together these mean the **science stack for a real non-grey RCE** —
`CorkLongwaveRadiation` + `CorkShortwaveRadiation` + `EmanuelConvectionPython`
+ `SlabSurface` + `get_default_state`/`get_grid` + sympl steppers — is already
pure Python and verified to instantiate and step (2026-07-19). The pedagogical
target is no longer "gray-only toy" but a genuine non-grey RCE in the browser.

### Performance is not the blocker

Benchmark (single column, 28 levels, native CPython, `climt` env, 2026-07-19),
comparing numba-JIT vs. `NUMBA_DISABLE_JIT=1` (the closest analog to Pyodide's
no-numba path):

| Call        | JIT (ms) | No-JIT (ms) |
|-------------|---------:|------------:|
| CORK-LW     |     2.7  |        1.9  |
| CORK-SW     |     1.9  |        2.7  |
| Emanuel(py) |     3.2  |        2.2  |
| **RCE step**|   **7.8**|      **6.8**|

At this problem size numba gives **no meaningful speedup** — Python/sympl
wrapper overhead dominates the njit inner loops — so losing numba in Pyodide
costs essentially nothing. A full equilibrium run (~2000–5000 steps) is ~15–35 s
native; apply a ~3–10× WASM-interpreter factor and it is tens of seconds to a
few minutes in-browser. Acceptable for a "run it and watch" notebook. Gray
radiation remains the fast fallback.

### Revised blocker list (packaging/plumbing, not physics)

The original Design sections 1–6 still apply. Additions/corrections:

- **[DONE 2026-07-19] `package_data` stale path.** `setup.py` still globbed
  `_data/picket_fence/correlated_k/*.nc|*.npz` after the picket-fence→CORK
  rename; the tables live at `_data/cork/correlated_k/`. Fixed to point at
  `_data/cork/...` (14 `.nc` + 2 `.npz` now matched). Latent bug: CORK tables
  shipped in *no* built wheel — masked because the dev install is editable.

- **`packages=["climt"]` excludes subpackages.** `setup.py` lists only the
  top-level package (no `find_packages`), so a built wheel would omit
  `climt._components`, `climt._core`, `climt._data`, etc. This works today only
  because the environment install is editable. **Must switch to
  `find_packages()`** (and confirm `_data` subpackages carry `__init__.py`, which
  they do) before any wheel — pure or platform — is reliable. Verify with an
  actual `pip wheel` + fresh-venv install, not just an editable checkout.

- **scipy: one usage the original spec missed.** Beyond §4a (IceSheet TDMA) and
  §4b (`initialization.py` `CubicSpline`→`np.interp`), CORK reads k-tables via
  `scipy.io.netcdf_file` (`cork/optics/correlated_k.py`). Options: (i) keep
  scipy — it *is* available in Pyodide, least work, but heavy and against the
  scipy-removal goal; or (ii) add a small numpy/`.npz` table reader and prefer
  shipping `.npz` tables (the loader already supports `.npz`), so the pure wheel
  drops scipy entirely. Recommend (ii) for the pure wheel.

- **`initialization.py` scipy import is top-level**, so it loads on `import climt`
  — §4b is therefore required for a clean scipy-free import, not just optional.

- **CORK not re-exported at top level.** `CorkLongwaveRadiation` /
  `CorkShortwaveRadiation` are only reachable via `climt._components`; the API
  reference lists them as public. Add them to `climt/__init__.py` `__all__`.

### Updated component availability (adds CORK)

| Component | Pure Python Wheel | Full Wheel |
|-----------|------------------|------------|
| CorkLongwaveRadiation | **Yes** | Yes |
| CorkShortwaveRadiation | **Yes** | Yes |
| EmanuelConvectionPython (v3) | **Yes** | Yes |
| GrayLongwaveRadiation / Frierson06 | Yes | Yes |
| SlabSurface / IceSheet / BucketHydrology | Yes | Yes |
| HeldSuarez / DryConvectiveAdjustment / Instellation | Yes | Yes |
| BergerSolarInsolation | Yes | Yes |
| EmanuelConvection (Fortran) | Error on init | Yes |
| RRTMGLongwave / RRTMGShortwave | Error on init | Yes |
| SimplePhysics / DcmipInitialConditions | Error on init | Yes |

### Suggested Phase 2 task order

1. `find_packages()` + `package_data` (done) — verify with a real wheel build/install.
2. Silent Fortran degradation (§2) + `has_fortran_extensions()` (§5).
3. scipy removal: §4b (`np.interp`, unblocks clean import) → §4a (TDMA) →
   CORK `.npz`/numpy netcdf reader → drop `scipy` from `install_requires` (§4c).
4. `CLIMT_PURE_PYTHON=1` build path + `py3-none-any` wheel + CI job (§1).
5. Ship CORK at top-level `__all__`.
6. Prove RCE end-to-end in JupyterLite; wire the deferred website live-cells.
