# Pyodide-Compatible Pure Python climt

**Date**: 2026-03-29
**Status**: Draft

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
