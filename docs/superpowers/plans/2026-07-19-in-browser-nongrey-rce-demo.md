# In-Browser Non-Grey Radiative Equilibrium Demo — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make climt install and run in the browser (Pyodide), and turn the website's Radiative Transfer theory track into editable, in-browser-runnable code — centred on a flagship page where gray vs non-grey radiative equilibrium is the *same* CORK component with a different k-table.

**Architecture:** Two halves on one branch. **Phase A** finishes the minimal pure-Python packaging so `climt-*-py3-none-any.whl` builds, imports without scipy/Fortran, and ships the k-tables. **Phase B** hosts that wheel as a GitHub release asset and adds `quarto-live` (Pyodide) editable cells to the docs — a flagship RE page plus Ch.1/Ch.6 retrofits — with the boot/install cell written as a reusable template.

**Tech Stack:** Python, numpy, sympl, setuptools, numba (optional), Cython (build-only), Quarto, `quarto-live` (Pyodide + micropip), matplotlib (in-browser via Pyodide).

## Global Constraints

- **Branch:** all work continues on `feature/pyodide-cork-prep` (open PR #217). Do **not** branch off `develop`.
- **Conda env:** run all Python/pytest commands in the `climt` conda env (`conda activate climt`).
- **Pure-Python trigger:** the pure wheel is built with env var `CLIMT_PURE_PYTHON=1`. In that mode setup.py must skip Cython import, compiler probing, and `ext_modules`.
- **No new hard dependencies.** Pure-wheel runtime deps: `numpy`, `pint`, `unyt`, `xarray`, `sympl`, `importlib_resources`. `scipy` is removed from `install_requires`. `numba` stays optional (already degrades to no-ops via `climt/_core/backend.py`).
- **Gray vs non-grey mechanism:** both flagship cases use `CorkLongwaveRadiation(optics="correlated_k", table=...)`; only the table string changes (`single_band_unit_lw` = gray, `earth_low_res_lw` = non-grey). No `GrayLongwaveRadiation` in the flagship.
- **Wheel hosting:** GitHub release asset, installed by pinned URL via `micropip`. The version/tag is defined once in a shared Quarto include.
- **Spec:** `docs/superpowers/specs/2026-07-19-in-browser-nongrey-rce-demo-design.md` is the source of truth; where this plan and the spec disagree, stop and reconcile.
- **Commit** after every task (frequent commits). End commit messages with the `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>` trailer.

---

# Phase A — Pure-Python packaging

## Task 1: `find_packages()` + `CLIMT_PURE_PYTHON` pure-wheel build path

**Files:**
- Modify: `setup.py`
- Test: `tests/test_pure_wheel_build.py` (create)

**Interfaces:**
- Produces: a buildable `climt-*-py3-none-any.whl` under `CLIMT_PURE_PYTHON=1` that contains all subpackages and the CORK data files; consumed by Tasks 5, 6, 8 and Phase B.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pure_wheel_build.py
import os
import subprocess
import sys
import glob
import zipfile
from pathlib import Path


def test_pure_wheel_builds_and_contains_subpackages_and_data(tmp_path):
    repo = Path(__file__).resolve().parent.parent
    env = dict(os.environ, CLIMT_PURE_PYTHON="1")
    out = tmp_path / "wheelhouse"
    subprocess.run(
        [sys.executable, "-m", "pip", "wheel", ".", "--no-deps", "-w", str(out)],
        cwd=repo, env=env, check=True,
    )
    wheels = glob.glob(str(out / "climt-*-py3-none-any.whl"))
    assert wheels, "expected a py3-none-any wheel"
    names = zipfile.ZipFile(wheels[0]).namelist()
    # subpackages present
    assert any(n.endswith("climt/_components/cork/lw/component.py") for n in names)
    assert any(n.endswith("climt/_core/initialization.py") for n in names)
    # data files present
    assert any("climt/_data/cork/correlated_k/" in n and n.endswith(".nc") for n in names)
    assert any(n.endswith("climt/_data/ozone_profile.npy") for n in names)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_pure_wheel_build.py -v`
Expected: FAIL — today `packages=["climt"]` omits subpackages, and there is no `CLIMT_PURE_PYTHON` path (Cython/compiler code runs).

- [ ] **Step 3: Add the pure-python guard and `find_packages` to setup.py**

Add the import and flag near the top of `setup.py` (after the stdlib imports, before the Cython block):

```python
from setuptools import Extension, find_packages, setup

PURE_PYTHON = os.environ.get("CLIMT_PURE_PYTHON") == "1"
```

Guard the Cython import block (currently ~lines 18–23) so it is skipped in pure mode:

```python
if not PURE_PYTHON:
    try:
        from Cython.Build.Distutils import build_ext as native_build_ext
    except ImportError:
        print("Suitable Cython unavailable, installing...")
        pip_main(["install", "cython"])
        from Cython.Build.Distutils import build_ext as native_build_ext
```

Guard the compiler-configuration and extension section. Replace the
`if os.environ.get("READTHEDOCS") == "True":` extension guard so both READTHEDOCS
and pure mode produce no extensions and no custom build classes:

```python
if PURE_PYTHON or os.environ.get("READTHEDOCS") == "True":
    ext_modules = []
    cmdclass = {}
else:
    ext_modules = [ ... existing Extension(...) list unchanged ... ]
    cmdclass = {
        "build_ext": climt_build_ext,
        "bdist_wheel": climt_bdist_wheel,
    }
```

Also wrap the compiler-probing block (everything from `operating_system = platform.system()`
through the `print("Lib Paths: ", lib_path_list)` line, and the `climt_build_ext` /
`climt_bdist_wheel` class definitions) in `if not PURE_PYTHON:` so pure builds do
not need gfortran. In the `setup(...)` call, replace:

```python
    packages=[
        "climt",
    ],
    package_dir={"climt": "climt"},
    ...
    cmdclass={
        "build_ext": climt_build_ext,
        "bdist_wheel": climt_bdist_wheel,
    },
    ext_modules=ext_modules,
```

with:

```python
    packages=find_packages(include=["climt", "climt.*"]),
    package_dir={"climt": "climt"},
    ...
    cmdclass=cmdclass,
    ext_modules=ext_modules,
```

Keep `package_data` unchanged (its CORK paths are already correct).

- [ ] **Step 4: Run test to verify it passes**

Run: `conda activate climt && pytest tests/test_pure_wheel_build.py -v`
Expected: PASS. Also confirm the non-pure build still lists extensions:
Run: `python -c "import os; os.environ.pop('CLIMT_PURE_PYTHON', None); print('ok')"`

- [ ] **Step 5: Fresh-venv smoke check (manual, recorded in commit body)**

```bash
python -m venv /tmp/climt_pure && . /tmp/climt_pure/bin/activate
pip install numpy scipy xarray pint unyt sympl importlib_resources
pip install $(ls /tmp/pytest-*/**/climt-*-py3-none-any.whl 2>/dev/null | head -1) || \
  { CLIMT_PURE_PYTHON=1 pip wheel . --no-deps -w /tmp/wh && pip install /tmp/wh/climt-*-py3-none-any.whl; }
python -c "import climt; from climt._components.cork import CorkLongwaveRadiation; print('import OK')"
deactivate
```
Expected: `import OK` (scipy still present here; scipy removal is Tasks 3–6).

- [ ] **Step 6: Commit**

```bash
git add setup.py tests/test_pure_wheel_build.py
git commit -m "build: find_packages + CLIMT_PURE_PYTHON pure-wheel path"
```

---

## Task 2: Silent Fortran degradation + `has_fortran_extensions()`

**Files:**
- Modify: `climt/_components/simple_physics/component.py`
- Modify: `climt/_components/emanuel/component.py`
- Modify: `climt/_components/rrtmg/lw/component.py`
- Modify: `climt/_components/rrtmg/sw/component.py`
- Modify: `climt/_components/dcmip/component.py`
- Modify: `climt/__init__.py`
- Test: `tests/test_fortran_degradation.py` (create)

**Interfaces:**
- Produces: `climt.has_fortran_extensions() -> bool`; each Fortran component raises `ImportError` at instantiation when its extension is absent, with no import-time warning.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_fortran_degradation.py
import climt


def test_has_fortran_extensions_is_bool():
    assert isinstance(climt.has_fortran_extensions(), bool)


def test_no_import_warning(recwarn):
    import importlib
    importlib.reload(climt)
    assert not [w for w in recwarn if "compiled" in str(w.message).lower()]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_fortran_degradation.py -v`
Expected: FAIL — `has_fortran_extensions` does not exist.

- [ ] **Step 3: Rewrite each Fortran component's import guard**

For `simple_physics/component.py`, replace the current try/except (which warns and prints) with:

```python
_FORTRAN_AVAILABLE = False
try:
    from . import _simple_physics as phys
    _FORTRAN_AVAILABLE = True
except ImportError:
    phys = None
```

and at the top of `SimplePhysics.__init__` (before any other work):

```python
        if not _FORTRAN_AVAILABLE:
            raise ImportError(
                "SimplePhysics requires compiled Fortran extensions, which are "
                "not available in this (pure-Python) install of climt."
            )
```

Apply the identical pattern to the other four, using each file's own import
target and class name, and remove the `logging.warning(...)`/`print(error)` lines:

- `emanuel/component.py`: `from . import _emanuel_convection` → flag `_FORTRAN_AVAILABLE`; guard in `EmanuelConvection.__init__` (message names "EmanuelConvection (Fortran)"; note the pure-Python `EmanuelConvectionPython` alternative).
- `rrtmg/lw/component.py`: `from . import _rrtmg_lw`; guard in `RRTMGLongwave.__init__`.
- `rrtmg/sw/component.py`: `from . import _rrtmg_sw`; guard in `RRTMGShortwave.__init__`.
- `dcmip/component.py`: `from . import _dcmip`; guard in `DcmipInitialConditions.__init__`.

- [ ] **Step 4: Add `has_fortran_extensions()` to `climt/__init__.py`**

Add (near the bottom, before `__version__`):

```python
def has_fortran_extensions():
    """Return True if compiled Fortran extensions are available."""
    try:
        from climt._components.simple_physics import _simple_physics  # noqa: F401
        return True
    except ImportError:
        return False
```

- [ ] **Step 5: Run test to verify it passes**

Run: `conda activate climt && pytest tests/test_fortran_degradation.py -v`
Expected: PASS (in the dev env with Fortran built, `has_fortran_extensions()` returns True; the warning assertion passes because warnings were removed).

- [ ] **Step 6: Commit**

```bash
git add climt/_components/*/component.py climt/__init__.py tests/test_fortran_degradation.py
git commit -m "feat: silent Fortran degradation + has_fortran_extensions()"
```

---

## Task 3: scipy removal — `initialization.py` CubicSpline → `np.interp`

**Files:**
- Modify: `climt/_core/initialization.py`
- Test: `tests/test_initialization.py` (add a case)

**Interfaces:**
- Produces: `climt/_core/initialization.py` no longer imports scipy at module load.

- [ ] **Step 1: Find the CubicSpline usage**

Run: `conda activate climt && grep -n "CubicSpline" climt/_core/initialization.py`
Expected: the top-level import (line 5) and one call site (ozone-profile interpolation onto model levels).

- [ ] **Step 2: Write the failing test**

```python
# tests/test_initialization.py  (append)
def test_initialization_does_not_import_scipy():
    import sys, importlib
    for m in list(sys.modules):
        if m == "scipy" or m.startswith("scipy."):
            del sys.modules[m]
    import climt._core.initialization as init
    importlib.reload(init)
    assert not any(m == "scipy" or m.startswith("scipy.") for m in sys.modules), \
        "initialization.py must not import scipy"
```

- [ ] **Step 3: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_initialization.py::test_initialization_does_not_import_scipy -v`
Expected: FAIL — scipy imported at module top.

- [ ] **Step 4: Replace CubicSpline with np.interp**

Remove `from scipy.interpolate import CubicSpline`. At the call site, replace the
spline evaluation with linear interpolation. The reference ozone profile is
coarse (30 points), so `np.interp` is sufficient:

```python
    # was: spline = CubicSpline(reference_pressure, reference_ozone)
    #      ozone_on_levels = spline(target_pressure)
    ozone_on_levels = np.interp(
        target_pressure, reference_pressure, reference_ozone
    )
```

`np.interp` requires the reference x-array to be increasing; if
`reference_pressure` is descending (TOA→surface), pass reversed views:
`np.interp(target_pressure, reference_pressure[::-1], reference_ozone[::-1])`.
Confirm ordering at the call site and pick the correct form.

- [ ] **Step 5: Run tests to verify they pass**

Run: `conda activate climt && pytest tests/test_initialization.py -v`
Expected: PASS (new no-scipy test and existing initialization tests).

- [ ] **Step 6: Commit**

```bash
git add climt/_core/initialization.py tests/test_initialization.py
git commit -m "refactor: ozone init uses np.interp, drops scipy CubicSpline"
```

---

## Task 4: TDMA solver + IceSheet rewrite (drop scipy.sparse)

**Files:**
- Create: `climt/_core/tridiagonal.py`
- Modify: `climt/_core/__init__.py` (export `solve_tridiagonal`)
- Modify: `climt/_components/surface_ice.py`
- Test: `tests/test_tridiagonal.py` (create)
- Test: `tests/test_surface_ice_tdma.py` (create — characterization)

**Interfaces:**
- Produces: `solve_tridiagonal(a_lower, a_diag, a_upper, rhs) -> x` where all band arrays have length `n` (`a_lower[0]` and `a_upper[n-1]` are unused), solving the tridiagonal system `A x = rhs`. Consumed by `IceSheet.calculate_new_ice_temperature`.

- [ ] **Step 1: Write the failing solver test**

```python
# tests/test_tridiagonal.py
import numpy as np
from climt._core.tridiagonal import solve_tridiagonal


def _dense(a_lower, a_diag, a_upper):
    n = a_diag.size
    A = np.diag(a_diag).astype(float)
    for i in range(1, n):
        A[i, i - 1] = a_lower[i]
    for i in range(n - 1):
        A[i, i + 1] = a_upper[i]
    return A


def test_matches_numpy_solve_on_random_systems():
    rng = np.random.default_rng(0)
    for n in (2, 5, 20):
        a_lower = np.zeros(n); a_lower[1:] = rng.normal(size=n - 1)
        a_upper = np.zeros(n); a_upper[:-1] = rng.normal(size=n - 1)
        a_diag = rng.normal(size=n) + 5.0  # diagonally dominant
        rhs = rng.normal(size=n)
        x = solve_tridiagonal(a_lower, a_diag, a_upper, rhs)
        np.testing.assert_allclose(x, np.linalg.solve(_dense(a_lower, a_diag, a_upper), rhs), rtol=1e-10)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_tridiagonal.py -v`
Expected: FAIL — module does not exist.

- [ ] **Step 3: Implement the Thomas algorithm**

```python
# climt/_core/tridiagonal.py
import numpy as np

from .backend import jit_compile


@jit_compile(backend=np)
def solve_tridiagonal(a_lower, a_diag, a_upper, rhs):
    """Solve the tridiagonal system A x = rhs via the Thomas algorithm.

    All band arrays have length n. a_lower[i] is A[i, i-1] (a_lower[0] unused);
    a_upper[i] is A[i, i+1] (a_upper[n-1] unused); a_diag[i] is A[i, i].
    """
    n = rhs.shape[0]
    cp = np.zeros(n)
    dp = np.zeros(n)
    x = np.zeros(n)
    cp[0] = a_upper[0] / a_diag[0]
    dp[0] = rhs[0] / a_diag[0]
    for i in range(1, n):
        m = a_diag[i] - a_lower[i] * cp[i - 1]
        cp[i] = a_upper[i] / m
        dp[i] = (rhs[i] - a_lower[i] * dp[i - 1]) / m
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x
```

Add to `climt/_core/__init__.py`: `from .tridiagonal import solve_tridiagonal` and include it in that module's `__all__` if one is defined.

- [ ] **Step 4: Run solver test to verify it passes**

Run: `conda activate climt && pytest tests/test_tridiagonal.py -v`
Expected: PASS.

- [ ] **Step 5: Write an IceSheet characterization test (golden, pre-refactor)**

```python
# tests/test_surface_ice_tdma.py
import numpy as np
from datetime import timedelta
from climt import IceSheet, get_default_state


def _run_once():
    ice = IceSheet()
    state = get_default_state([ice])
    # ensure a non-trivial ice column so calculate_new_ice_temperature runs
    state["area_type"].values[:] = "sea_ice"
    state["sea_ice_thickness"].values[:] = 2.0
    state["snow_and_ice_temperature"].values[:] = 260.0
    diag, new = ice(state, timedelta(hours=1))
    return new["snow_and_ice_temperature"].values.copy()


def test_icesheet_matches_golden():
    golden = np.load("tests/_golden/icesheet_snow_ice_temperature.npy")
    np.testing.assert_allclose(_run_once(), golden, rtol=1e-9, atol=1e-9)
```

Generate the golden array from the **current (scipy)** implementation, before
editing `surface_ice.py`:

```bash
conda activate climt
mkdir -p tests/_golden
python -c "
import numpy as np
from tests.test_surface_ice_tdma import _run_once
np.save('tests/_golden/icesheet_snow_ice_temperature.npy', _run_once())
"
```
(If the default `sea_ice` column produces an all-copy result with no solve,
adjust the seed state until `calculate_new_ice_temperature` is exercised, then
regenerate the golden.)

- [ ] **Step 6: Rewrite `calculate_new_ice_temperature` to use `solve_tridiagonal`**

Remove the scipy imports (`from scipy import sparse`, `from scipy.sparse.linalg import spsolve`).
Replace the matrix build/solve (current lines ~431–459) with band arrays. The
current `spdiags([a_sub, dp, a_sup], [-1,0,1], n, n)` maps to bands as
`A[i,i-1]=a_sub[i-1]`, `A[i,i]=dp[i]`, `A[i,i+1]=a_sup[i+1]`; the RHS matrix uses
`[-a_sub, dm, -a_sup]`:

```python
        from .._core.tridiagonal import solve_tridiagonal

        n = num_layers
        lower = np.zeros(n)      # lower[i] = A[i, i-1]
        upper = np.zeros(n)      # upper[i] = A[i, i+1]
        diag = dp.copy()         # diag[i]  = A[i, i]
        lower[1:] = a_sub[:-1]
        upper[:-1] = a_sup[1:]

        # rhs = mat_rhs @ temp_profile, mat_rhs diagonals [-a_sub, dm, -a_sup]
        rhs = dm * temp_profile
        rhs[1:] += -a_sub[:-1] * temp_profile[:-1]
        rhs[:-1] += -a_sup[1:] * temp_profile[1:]

        # boundary conditions (identical to the previous mat_lhs/rhs edits)
        if surf_temperature < self._melting_temperature - self._epsilon:
            diag[-1] = -1.0
            lower[-1] = 1.0
            rhs[-1] = -net_flux * dz / K_mid[-1]
        else:
            diag[-1] = 1.0
            lower[-1] = 0.0
            rhs[-1] = self._melting_temperature

        diag[0] = 1.0
        upper[0] = 0.0
        rhs[0] = self._melting_temperature if soil_temperature is None else soil_temperature

        return solve_tridiagonal(lower, diag, upper, rhs)
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `conda activate climt && pytest tests/test_tridiagonal.py tests/test_surface_ice_tdma.py -v`
Expected: PASS — the golden IceSheet output is reproduced bit-close by the TDMA path.

- [ ] **Step 8: Confirm scipy is gone from surface_ice**

Run: `conda activate climt && grep -n scipy climt/_components/surface_ice.py`
Expected: no matches (delete the commented `# from scipy...` line too).

- [ ] **Step 9: Commit**

```bash
git add climt/_core/tridiagonal.py climt/_core/__init__.py climt/_components/surface_ice.py tests/test_tridiagonal.py tests/test_surface_ice_tdma.py tests/_golden/
git commit -m "refactor: IceSheet uses pure-numpy TDMA solver, drops scipy.sparse"
```

---

## Task 5: CORK `.npz` table reader + ship demo tables as `.npz`

**Files:**
- Modify: `climt/_components/cork/optics/correlated_k.py`
- Create: `scripts/convert_ck_table_to_npz.py`
- Create (data): `climt/_data/cork/correlated_k/single_band_unit_lw.npz`, `earth_low_res_lw.npz`, `earth_low_res_sw.npz`
- Modify: `setup.py` `package_data` (already globs `*.npz`; confirm)
- Test: `tests/test_correlated_k_npz.py` (create)

**Interfaces:**
- Produces: `load_k_table(name)` resolves a `.npz` table without importing scipy; the three demo tables exist as `.npz`. Consumed by Phase B cells and Task 6.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_correlated_k_npz.py
import sys
import numpy as np
from climt._components.cork.optics.correlated_k import load_k_table


def test_npz_table_loads_without_scipy(monkeypatch):
    monkeypatch.setitem(sys.modules, "scipy", None)  # force ImportError on scipy use
    for name in ("single_band_unit_lw", "earth_low_res_lw", "earth_low_res_sw"):
        table = load_k_table(name)
        assert "k_coefficients" in table
        assert np.asarray(table["k_coefficients"]).ndim >= 5
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_correlated_k_npz.py -v`
Expected: FAIL — those `.npz` files do not exist yet; `load_k_table` prefers `.nc` (scipy).

- [ ] **Step 3: Write the `.nc`→`.npz` converter**

```python
# scripts/convert_ck_table_to_npz.py
"""Convert a CORK correlated-k .nc table to .npz (numpy-native, scipy-free reader)."""
import sys
from pathlib import Path
import numpy as np
from climt._components.cork.optics.correlated_k import _load_netcdf_table


def convert(nc_path):
    nc_path = Path(nc_path)
    table = _load_netcdf_table(str(nc_path))  # uses scipy at build time (dev env only)
    out = nc_path.with_suffix(".npz")
    np.savez_compressed(out, **{k: np.asarray(v) for k, v in table.items()})
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        convert(p)
```

Run it (dev env has scipy):

```bash
conda activate climt
python scripts/convert_ck_table_to_npz.py \
  climt/_data/cork/correlated_k/single_band_unit_lw.nc \
  climt/_data/cork/correlated_k/earth_low_res_lw.nc \
  climt/_data/cork/correlated_k/earth_low_res_sw.nc
```

- [ ] **Step 4: Make `load_k_table` resolve `.npz` first and give a clear scipy error**

In `load_k_table`, change the shipped-table resolution to prefer `.npz`, and in
`_load_netcdf_table` raise a clear error if scipy is unavailable. Replace the
name-resolution tail of `load_k_table`:

```python
    pkg = importlib_resources.files("climt._data.cork.correlated_k")
    npz_path = pkg.joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(npz_path) as f:
        if os.path.isfile(f):
            return np.load(f, allow_pickle=True)

    nc_path = pkg.joinpath(f"{name_or_path}.nc")
    with importlib_resources.as_file(nc_path) as f:
        if os.path.isfile(f):
            return _load_netcdf_table(str(f))

    raise FileNotFoundError(f"No k-table named {name_or_path!r} (.npz or .nc)")
```

And at the top of `_load_netcdf_table`:

```python
    try:
        from scipy.io import netcdf_file
    except ImportError as exc:
        raise ImportError(
            "Reading .nc k-tables requires scipy. Install scipy, or use a .npz "
            "table (the shipped demo tables are available as .npz)."
        ) from exc
```

- [ ] **Step 5: Run test to verify it passes**

Run: `conda activate climt && pytest tests/test_correlated_k_npz.py -v`
Expected: PASS.

- [ ] **Step 6: Verify the CORK components still run from the .npz tables**

Run: `conda activate climt && pytest tests/test_grey_limit.py tests/test_cork_lw.py -v`
Expected: PASS (grey-limit test still ties CORK single-band to GrayLongwave).

- [ ] **Step 7: Commit**

```bash
git add climt/_components/cork/optics/correlated_k.py scripts/convert_ck_table_to_npz.py climt/_data/cork/correlated_k/*.npz tests/test_correlated_k_npz.py
git commit -m "feat: scipy-free .npz CORK table reader; ship demo tables as .npz"
```

---

## Task 6: Drop scipy from `install_requires`; assert scipy-free import

**Files:**
- Modify: `setup.py`
- Test: `tests/test_no_scipy_import.py` (create)

**Interfaces:**
- Produces: `import climt` (+ constructing the demo components) imports no scipy.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_no_scipy_import.py
import subprocess
import sys


def test_import_climt_does_not_import_scipy():
    code = (
        "import sys, climt\n"
        "from climt import get_default_state, get_grid\n"
        "from climt._components.cork import CorkLongwaveRadiation\n"
        "lw = CorkLongwaveRadiation(optics='correlated_k', table='single_band_unit_lw')\n"
        "state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=20))\n"
        "lw(state)\n"
        "assert not any(m=='scipy' or m.startswith('scipy.') for m in sys.modules), "
        "sorted(m for m in sys.modules if m.startswith('scipy'))\n"
        "print('no-scipy OK')\n"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
```

- [ ] **Step 2: Run test to verify it fails or passes**

Run: `conda activate climt && pytest tests/test_no_scipy_import.py -v`
Expected: After Tasks 3–5 this may already PASS; if it FAILS, the failure message lists the offending scipy submodule — trace and remove that import.

- [ ] **Step 3: Remove scipy from requirements**

In `setup.py`, delete `"scipy>=0.18.1",` from the `requirements` list.

- [ ] **Step 4: Run test to verify it passes**

Run: `conda activate climt && pytest tests/test_no_scipy_import.py -v`
Expected: PASS — `no-scipy OK`.

- [ ] **Step 5: Commit**

```bash
git add setup.py tests/test_no_scipy_import.py
git commit -m "build: drop scipy from install_requires; guard scipy-free import"
```

---

## Task 7: Re-export CORK at top level

**Files:**
- Modify: `climt/__init__.py`
- Test: `tests/test_cork_toplevel_export.py` (create)

**Interfaces:**
- Produces: `climt.CorkLongwaveRadiation`, `climt.CorkShortwaveRadiation`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cork_toplevel_export.py
import climt


def test_cork_is_top_level():
    assert hasattr(climt, "CorkLongwaveRadiation")
    assert hasattr(climt, "CorkShortwaveRadiation")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_cork_toplevel_export.py -v`
Expected: FAIL — not exported at top level.

- [ ] **Step 3: Add the imports and `__all__` entries**

In `climt/__init__.py`, add `CorkLongwaveRadiation, CorkShortwaveRadiation` to the
`from ._components import (...)` block, and add both names to the `__all__` tuple.

- [ ] **Step 4: Run test to verify it passes**

Run: `conda activate climt && pytest tests/test_cork_toplevel_export.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/__init__.py tests/test_cork_toplevel_export.py
git commit -m "feat: re-export CorkLongwaveRadiation/CorkShortwaveRadiation at top level"
```

---

## Task 8: CI job — build & publish the pure wheel on release

**Files:**
- Modify: `.github/workflows/release_climt.yml`
- Test: manual (workflow lint + inspect)

**Interfaces:**
- Produces: on tagged release, a `climt-*-py3-none-any.whl` is built and attached to the GitHub release. This is the asset Phase B installs by URL.

- [ ] **Step 1: Read the current release workflow**

Run: `sed -n '1,200p' .github/workflows/release_climt.yml`
Identify the release trigger and how platform wheels are uploaded.

- [ ] **Step 2: Add a pure-wheel job**

Add a job that runs on `ubuntu-latest`, checks out, sets up Python, then:

```yaml
  build-pure-wheel:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: "3.11" }
      - name: Build pure-Python wheel
        env:
          CLIMT_PURE_PYTHON: "1"
        run: |
          python -m pip install --upgrade pip build
          python -m build --wheel
          ls dist/*-py3-none-any.whl
      - name: Attach to release
        uses: softprops/action-gh-release@v2
        with:
          files: dist/*-py3-none-any.whl
```

Match the auth/trigger conventions already used in the file (e.g. reuse the
existing `on: release` trigger and any `GITHUB_TOKEN` permissions block).

- [ ] **Step 3: Validate the workflow YAML**

Run: `python -c "import yaml,sys; yaml.safe_load(open('.github/workflows/release_climt.yml')); print('yaml OK')"`
Expected: `yaml OK`.

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/release_climt.yml
git commit -m "ci: build and attach py3-none-any pure wheel on release"
```

---

# Phase B — In-browser demo (quarto-live)

## Task 9: quarto-live spike (de-risk the runtime)

**Files:**
- Create: `docs/_spikes/quarto-live-smoke.qmd` (throwaway; deleted or kept as a smoke page)
- Modify: `docs/_quarto.yml` (add the `quarto-live` filter/extension)

**Interfaces:**
- Produces: a verified answer to "can quarto-live micropip-install a wheel from a URL, import climt, and render matplotlib?" — gate for Tasks 10–15. If it fails, switch to the hand-rolled Pyodide include fallback (documented in the spec) and note it here.

- [ ] **Step 1: Install the quarto-live extension**

Run: `cd docs && quarto add r-wasm/quarto-live` (accept prompts).
Expected: `_extensions/r-wasm/live/` created.

- [ ] **Step 2: Write a minimal live page that installs climt from a URL**

```markdown
---
title: quarto-live smoke test
engine: knitr
filters:
  - live
---

```{pyodide}
#| autorun: true
import micropip
await micropip.install("https://github.com/CliMT/climt/releases/download/PLACEHOLDER/climt-0.20.0-py3-none-any.whl")
import climt
print("climt", climt.__version__, "has_fortran:", climt.has_fortran_extensions())
```

```{pyodide}
import numpy as np, matplotlib.pyplot as plt
plt.plot(np.linspace(0, 1, 50), np.sqrt(np.linspace(0, 1, 50)))
plt.title("matplotlib in Pyodide")
```
```

Until Task 10 produces a real release URL, host the Task-1 wheel temporarily:
either `python -m http.server` in the wheelhouse and use `http://localhost:8000/...`,
or upload a pre-release asset. Record which was used.

- [ ] **Step 3: Render and open**

Run: `cd docs && quarto preview _spikes/quarto-live-smoke.qmd`
Expected: page loads; both cells run; version prints; plot renders. Manually confirm in a browser.

- [ ] **Step 4: Record the outcome and decide**

Write findings (load time, wheel size felt, any failures) into the commit body.
If quarto-live cannot install from a URL or render matplotlib, STOP and switch to
the hand-rolled Pyodide + CodeMirror include; note the decision here before proceeding.

- [ ] **Step 5: Commit**

```bash
git add docs/_extensions docs/_quarto.yml docs/_spikes/quarto-live-smoke.qmd
git commit -m "spike: quarto-live Pyodide runtime validated for climt wheel"
```

---

## Task 10: Release the wheel + shared boot include

**Files:**
- Create: `docs/_includes/climt-live-setup.qmd` (shared boot + helper; the reusable template)
- Modify: `docs/_quarto.yml` (register the include / filter once)
- Create: `docs/radiative-transfer/_live/rce_helpers.py` (Python source the setup cell reads, so it is testable natively — see Task 11)

**Interfaces:**
- Consumes: the pure wheel from Task 1/Task 8.
- Produces: an includable Quarto snippet that boots Pyodide, installs the pinned wheel, and defines `integrate_to_equilibrium(...)`. The wheel URL/tag is defined once here.

- [ ] **Step 1: Publish the wheel as a GitHub release asset**

Build (Task 1 path) and attach to a release (or pre-release) tag, e.g. `web-demo-v0.20.0`:

```bash
conda activate climt
CLIMT_PURE_PYTHON=1 python -m build --wheel
gh release create web-demo-v0.20.0 dist/climt-0.20.0-py3-none-any.whl \
  --title "climt web demo wheel 0.20.0" --notes "Pure-Python wheel for in-browser docs." --prerelease
```
Record the exact asset URL.

- [ ] **Step 2: Write the shared boot include**

```markdown
<!-- docs/_includes/climt-live-setup.qmd -->
```{pyodide}
#| autorun: true
#| edit: false
import micropip
await micropip.install(
    "https://github.com/CliMT/climt/releases/download/web-demo-v0.20.0/climt-0.20.0-py3-none-any.whl"
)
import numpy as np
import matplotlib.pyplot as plt
import climt
from sympl import AdamsBashforth
from datetime import timedelta


def integrate_to_equilibrium(tendency_components, stepper_components, state,
                             timestep, n_steps):
    """Step tendency + stepper components forward toward radiative equilibrium."""
    model = AdamsBashforth(list(tendency_components))
    for _ in range(n_steps):
        diagnostics, new_state = model(state, timestep)
        state.update(diagnostics)
        state.update(new_state)
        for stepper in stepper_components:
            s_diag, s_new = stepper(state, timestep)
            state.update(s_diag)
            state.update(s_new)
        state["time"] += timestep
    return state
```
```

- [ ] **Step 3: Keep the helper in sync with a native module**

Copy the same `integrate_to_equilibrium` body into
`docs/radiative-transfer/_live/rce_helpers.py` (plain importable module) so Task 11
can unit-test the exact loop. Add a comment in both files pointing at each other.

- [ ] **Step 4: Verify the include renders and installs the real wheel**

Point the Task 9 smoke page's URL at the real release asset; `quarto preview` and
confirm install + `integrate_to_equilibrium` is defined (`print(integrate_to_equilibrium)`).

- [ ] **Step 5: Commit**

```bash
git add docs/_includes/climt-live-setup.qmd docs/_quarto.yml docs/radiative-transfer/_live/rce_helpers.py
git commit -m "feat: shared quarto-live boot include + release-hosted climt wheel"
```

---

## Task 11: Flagship page `09-live-rce.qmd` + native RE-loop test

**Files:**
- Create: `docs/radiative-transfer/09-live-rce.qmd`
- Modify: `docs/_quarto.yml` (add to the Radiative Transfer sidebar)
- Test: `tests/test_live_rce_demo.py` (create — the science the page shows)

**Interfaces:**
- Consumes: `integrate_to_equilibrium` (Task 10), CORK `.npz` tables (Task 5).
- Produces: the flagship page and a native test asserting the non-grey column develops a stratosphere the gray column lacks.

- [ ] **Step 1: Write the failing native science test**

```python
# tests/test_live_rce_demo.py
import numpy as np
from datetime import timedelta
import sympl
from climt import (get_default_state, get_grid, SlabSurface,
                   CorkLongwaveRadiation, CorkShortwaveRadiation)
from docs.radiative_transfer._live.rce_helpers import integrate_to_equilibrium


def _equilibrium(table):
    sympl.set_backend(sympl.DataArrayBackend())
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
    sw = CorkShortwaveRadiation(optics="correlated_k", table="earth_low_res_sw")
    surface = SlabSurface()
    grid = get_grid(nx=1, ny=1, nz=28)
    state = get_default_state([lw, sw, surface], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 3
    state["surface_albedo_for_direct_shortwave"].values[:] = 0.3
    state = integrate_to_equilibrium([lw, sw], [surface], state,
                                     timedelta(hours=6), n_steps=400)
    return state["air_temperature"].values[:, 0, 0]  # index 0 = surface

def test_nongrey_has_stratosphere_gray_does_not():
    T_gray = _equilibrium("single_band_unit_lw")
    T_nongrey = _equilibrium("earth_low_res_lw")
    # Non-grey: upper atmosphere warmer than a pure-grey monotonic decline
    # => temperature gradient near TOA flattens/reverses (a stratosphere).
    dT_gray = np.diff(T_gray)      # surface -> TOA
    dT_nongrey = np.diff(T_nongrey)
    assert dT_nongrey[-3:].mean() > dT_gray[-3:].mean() + 1.0
```

Add `docs/radiative-transfer/_live/__init__.py` and
`docs/radiative-transfer/__init__.py` (empty) so the test import path resolves,
or adjust the import to load the module by file path. Confirm the assertion's
sign against a quick manual run and tune `n_steps`/`nz`/tolerance so it is a
robust, non-flaky check of the real physics.

- [ ] **Step 2: Run test to verify it fails**

Run: `conda activate climt && pytest tests/test_live_rce_demo.py -v`
Expected: FAIL — helper module import or the assertion (until parameters are tuned).

- [ ] **Step 3: Make it pass (tune parameters, confirm the science)**

Adjust `n_steps`, `timestep`, `nz`, and the stratosphere metric until the test
robustly passes and reflects real equilibrium (temperatures stable to <0.1 K
over the last 10% of steps). This is systematic verification of the physics the
page will show; do not weaken the assertion to force a pass.

- [ ] **Step 4: Write the flagship page**

```markdown
---
title: "Radiative equilibrium: gray vs non-grey — live"
engine: knitr
filters: [live]
---

{{< include ../_includes/climt-live-setup.qmd >}}

Gray and non-grey radiation are the *same code* here — one CORK component —
differing only in the absorption table. Change one string and watch a
stratosphere appear.

## Gray column

```{pyodide}
from climt import (get_default_state, get_grid, SlabSurface,
                   CorkLongwaveRadiation, CorkShortwaveRadiation)
import sympl; from datetime import timedelta
sympl.set_backend(sympl.DataArrayBackend())

lw = CorkLongwaveRadiation(optics="correlated_k", table="single_band_unit_lw")
sw = CorkShortwaveRadiation(optics="correlated_k", table="earth_low_res_sw")
surface = SlabSurface()
grid = get_grid(nx=1, ny=1, nz=28)
state = get_default_state([lw, sw, surface], grid_state=grid)
state["zenith_angle"].values[:] = np.pi/3
state = integrate_to_equilibrium([lw, sw], [surface], state, timedelta(hours=6), 400)

p = state["air_pressure"].values[:,0,0]; T_gray = state["air_temperature"].values[:,0,0]
plt.plot(T_gray, p/100); plt.gca().invert_yaxis(); plt.xlabel("T (K)"); plt.ylabel("p (hPa)")
plt.title("Gray radiative equilibrium")
```

## Non-grey column — change one string

```{pyodide}
lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")  # <- the only change
state = get_default_state([lw, sw, surface], grid_state=grid)
state["zenith_angle"].values[:] = np.pi/3
state = integrate_to_equilibrium([lw, sw], [surface], state, timedelta(hours=6), 400)
T_ng = state["air_temperature"].values[:,0,0]
plt.plot(T_gray, p/100, label="gray"); plt.plot(T_ng, p/100, label="non-grey")
plt.gca().invert_yaxis(); plt.legend(); plt.xlabel("T (K)"); plt.ylabel("p (hPa)")
plt.title("Non-grey grows a stratosphere")
```
```

Add a per-band flux / OLR cell (Cell C) and an editable CO₂ cell (Cell D) —
Cell D sets the CO₂ VMR in the state before integrating and re-plots. Keep the
default `n_steps` at the tuned value with a comment on increasing it.

- [ ] **Step 5: Add to the sidebar and render**

In `docs/_quarto.yml`, add `radiative-transfer/09-live-rce.qmd` after `08-multiplanet.qmd`.
Run: `cd docs && quarto render radiative-transfer/09-live-rce.qmd`
Expected: renders without error. Manually preview and confirm both plots draw and the non-grey curve shows a stratosphere.

- [ ] **Step 6: Commit**

```bash
git add docs/radiative-transfer/09-live-rce.qmd docs/radiative-transfer/_live/__init__.py docs/radiative-transfer/__init__.py docs/_quarto.yml tests/test_live_rce_demo.py
git commit -m "feat: flagship live gray-vs-non-grey radiative-equilibrium page"
```

---

## Task 12: Ch.1 retrofit — live cells

**Files:**
- Modify: `docs/radiative-transfer/01-why-nongrey.qmd`

**Interfaces:**
- Consumes: `quarto-live` filter; the shared setup include (only if the RE teaser is used).

- [ ] **Step 1: Make the "mean of exp" toy live**

Add the `live` filter to the page front matter, then convert the Fig 1.1 toy to a
`{pyodide}` cell (pure numpy — no climt needed):

```{pyodide}
import numpy as np, matplotlib.pyplot as plt
L = np.linspace(0, 5, 200)
sigma = np.where(np.random.default_rng(0).random(1000) < 0.1, 20.0, 0.05)  # 10% strong
true = np.array([np.mean(np.exp(-sigma*x)) for x in L])
grey = np.exp(-sigma.mean()*L)
plt.plot(L, true, label="band-averaged (true)"); plt.plot(L, grey, label="grey")
plt.xlabel("path length"); plt.ylabel("transmission"); plt.legend()
```

- [ ] **Step 2: Replace the "open the notebook" callout with an in-page cell**

Replace the `::: {.callout-tip} Try it yourself … examples/…ipynb` block with a
small live cell that runs `CorkLongwaveRadiation` per-band flux (mirror Task 13's
cell, kept short), including the shared setup include at the top of the page.

- [ ] **Step 3: Render and preview**

Run: `cd docs && quarto render radiative-transfer/01-why-nongrey.qmd`
Expected: renders; the toy cell runs client-side.

- [ ] **Step 4: Commit**

```bash
git add docs/radiative-transfer/01-why-nongrey.qmd
git commit -m "docs: live cells in Ch.1 (mean-of-exp toy + CORK per-band teaser)"
```

---

## Task 13: Ch.6 retrofit — live CORK per-band flux

**Files:**
- Modify: `docs/radiative-transfer/06-picket-fence.qmd`

**Interfaces:**
- Consumes: shared setup include; `earth_low_res_lw` `.npz` table.

- [ ] **Step 1: Add the live per-band flux cell**

Include the setup snippet, then:

```{pyodide}
from climt import get_default_state, get_grid, CorkLongwaveRadiation
import sympl; sympl.set_backend(sympl.DataArrayBackend())
lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw", diagnostics_level=1)
state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=28))
tend, diag = lw(state)
# plot per-band upward flux at TOA (window vs CO2 band contrast)
key = [k for k in diag if "per_band" in k and "flux" in k]
print("per-band diagnostics:", key)
```

Adjust the diagnostic key to the actual per-band flux name CORK emits at
`diagnostics_level=1` (inspect `diag.keys()` in the dev env first), then plot the
window-vs-CO₂-band contrast.

- [ ] **Step 2: Remove the stale static-notebook pointer**

Replace any "see the notebook" callout with this in-page cell.

- [ ] **Step 3: Render and preview**

Run: `cd docs && quarto render radiative-transfer/06-picket-fence.qmd`
Expected: renders; cell runs; per-band contrast visible.

- [ ] **Step 4: Commit**

```bash
git add docs/radiative-transfer/06-picket-fence.qmd
git commit -m "docs: live CORK per-band flux diagnostic in Ch.6"
```

---

## Task 14: Reusable-template appendix page

**Files:**
- Create: `docs/radiative-transfer/live-code-template.qmd`
- Modify: `docs/_quarto.yml` (add under Appendices)

**Interfaces:**
- Produces: a copy-paste recipe for embedding runnable climt in any Quarto site.

- [ ] **Step 1: Write the how-to page**

Document: (1) `quarto add r-wasm/quarto-live`; (2) the front-matter `filters: [live]`;
(3) the `micropip.install("<release-url>")` one-liner; (4) the
`integrate_to_equilibrium` helper; (5) the page-weight/first-load caveat. Point
readers at `model_tour_climate` as the motivating use case.

- [ ] **Step 2: Add to the sidebar (Appendices) and render**

In `docs/_quarto.yml`, add `radiative-transfer/live-code-template.qmd` under the
existing Appendices section.
Run: `cd docs && quarto render radiative-transfer/live-code-template.qmd`
Expected: renders.

- [ ] **Step 3: Commit**

```bash
git add docs/radiative-transfer/live-code-template.qmd docs/_quarto.yml
git commit -m "docs: reusable 'embed runnable climt' template appendix"
```

---

## Task 15: End-to-end browser verification + page-weight budget

**Files:**
- Create: `docs/_spikes/weight-budget.md` (recorded numbers; or fold into the PR description)

**Interfaces:**
- Consumes: everything above.
- Produces: recorded evidence that the flagship page works in a real browser and the download budget is acceptable.

- [ ] **Step 1: Full site render**

Run: `cd docs && quarto render`
Expected: whole site builds with no errors.

- [ ] **Step 2: Real-browser end-to-end check**

Run: `cd docs && quarto preview`, open the flagship page in a browser. Confirm:
Pyodide boots, the wheel installs from the release URL, the gray cell runs and
plots, changing the table string and re-running produces the non-grey stratosphere,
and the CO₂ cell responds. (Automate with Playwright if desired; otherwise record
a manual pass with screenshots.)

- [ ] **Step 3: Record the page-weight budget**

From the browser network panel on first load, record: Pyodide core, numpy,
the climt wheel, and each `.npz` table (esp. `earth_low_res_lw`), plus
time-to-first-plot. Note whether `earth_low_res_lw.npz` needs slimming. Write the
numbers into `docs/_spikes/weight-budget.md`.

- [ ] **Step 4: Full native test sweep**

Run: `conda activate climt && pytest tests/test_pure_wheel_build.py tests/test_fortran_degradation.py tests/test_initialization.py tests/test_tridiagonal.py tests/test_surface_ice_tdma.py tests/test_correlated_k_npz.py tests/test_no_scipy_import.py tests/test_cork_toplevel_export.py tests/test_live_rce_demo.py tests/test_grey_limit.py -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add docs/_spikes/weight-budget.md
git commit -m "docs: record in-browser RCE demo e2e verification + weight budget"
```

---

## Self-Review notes (for the executor)

- **Spec coverage:** Tasks 1–8 cover packaging §1–6 + CORK top-level export + CI; Tasks 9–15 cover the quarto-live runtime, hosting, flagship, Ch.1/Ch.6 retrofits, template, and end-to-end verification. Ch.8 multiplanet and RCE-with-convection are intentionally deferred (spec "Deferred").
- **Ordering risk:** Task 5's `.npz` tables must exist before Tasks 11/13 use them, and before Task 6's no-scipy assertion (which loads `single_band_unit_lw`). Keep the order.
- **Helper drift:** the `integrate_to_equilibrium` body is duplicated in the Quarto include (Task 10) and the native module (Task 10/11). If you change one, change both — the native test only guards the module copy.
- **Diagnostic names:** the exact CORK per-band flux diagnostic key (Tasks 11 Cell C, 13) must be read from `diag.keys()` in the dev env; do not guess.
