# Package Split (pure `climt` + `climt-fortran`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `climt` a pure-Python distribution and move the five Cython/Fortran extensions into a companion distribution, `climt-fortran`, in the same repository, without changing what `pip install climt` gives a user.

**Architecture:** The compiled modules and Fortran sources move to `fortran/` (import package `climt_fortran`), which builds either real extensions or, with `CLIMT_FORTRAN_STUB=1`, an empty `py3-none-any` stub wheel. `climt` depends on `climt-fortran` and reaches it only through `climt/_core/fortran.py`, which reports the companion's state and hands components their extension module or `None`. `climt.list_components()` tells users what is usable.

**Tech Stack:** setuptools (PEP 621 `pyproject.toml` for the core; `setup.py` + Cython + make for the companion), `importlib.metadata`, pytest, tox, GitHub Actions, cibuildwheel 3.3.1.

**Spec:** `docs/superpowers/specs/2026-11-14-package-split-design.md`

## Global Constraints

- Python 3.10 or newer. flake8 with the repo's `tox.ini` settings (max line length 150) must pass on `climt` and on `fortran/climt_fortran`.
- Distribution name `climt-fortran`, import name `climt_fortran`, first version `1.0.0`, `INTERFACE_VERSION = 1`.
- The five extension modules are exactly `_dcmip`, `_emanuel_convection`, `_rrtmg_lw`, `_rrtmg_sw`, `_simple_physics`.
- The core's requirement string is exactly `climt-fortran>=1.0,<2; sys_platform != "emscripten"`.
- `import climt` prints nothing and raises nothing, whatever the companion's state (`tests/test_fortran_degradation.py::test_no_import_warning` must keep passing).
- The only module in `climt/` that may contain the string `climt_fortran` in an import is `climt/_core/fortran.py`.
- The core wheel is `py3-none-any` and contains no file ending `.so`, `.pyd`, `.pyx`, `.c`, `.a`, `.f90`, `.F90`, `.f` or `.mod`, and nothing under `fortran/` or `climt_fortran/`.
- Files are moved with `git mv`. No Fortran or Cython source is edited except the `touch` paths in `fortran/_lib/Makefile`.
- Runtime dependencies of the core are unchanged apart from dropping `cython` and adding `climt-fortran`. Do not bump climt's version; add to the `Unreleased` section of `HISTORY.rst`.
- If `setup.py` on the branch you start from declares `extras_require`, carry every entry over verbatim into `[project.optional-dependencies]`.
- **Environment.** Work in a dedicated conda environment named `climt-split` (Python 3.11, `pip install -r requirements_dev.txt`). Never install into the environments `climt-dev`, `climt` or `gfs_dycore_develop`.
- **Toolchain on the maintainer's Mac.** Homebrew gcc 15 needs `export SDKROOT="$(xcrun --show-sdk-path)"` before any build of the extensions. If the extensions do not build even on the untouched starting commit (Task 2, Step 1), stop and report BLOCKED with the compiler output; do not work around a broken toolchain inside this plan.
- Nothing is published or pushed. Workflow files are written and syntax-checked only.

## Review Focus

1. **A binary that is present but will not load** (NumPy ABI mismatch, issue #174): `import climt` stays silent, the status is `"broken"`, and the component's error quotes the loader's message. Test in Task 1.
2. **`climt_fortran` not importable at all** (Pyodide, `--no-deps`, after `pip uninstall climt-fortran`): `import climt`, `has_fortran_extensions()` and `list_components()` all work. Tests in Tasks 1, 2 and 4.
3. **Core and companion out of step** (`INTERFACE_VERSION` differs): every extension is treated as unavailable and the error says to upgrade both. Test in Task 1.
4. **A developer checkout with stale build products** (old `.so` files beside the components, `.a` and `.c` files under `fortran/`): old `.so` files are never imported, and the core wheel built from that tree is still clean. Tests in Tasks 2 and 3.
5. **The with-binaries test job passing by skipping**: with `CLIMT_REQUIRE_FORTRAN=1` a `FortranUnavailableError` is a failure, not a skip. Test in Task 2.

## File Structure

| Path | Responsibility |
|---|---|
| `climt/_core/fortran.py` (new) | Sole contact with `climt_fortran`: status, module loading, error construction |
| `climt/_core/registry.py` (new) | `list_components()` and the `climt.components` entry-point group |
| `climt/_components/{dcmip,emanuel,simple_physics}/component.py`, `rrtmg/{lw,sw}/component.py` | Ask the loader for their extension |
| `climt/__init__.py` | Export `FortranUnavailableError`, `list_components`; `has_fortran_extensions()` via the loader |
| `pyproject.toml` | Core metadata, dependencies, extras, package data |
| `setup.py` (root) | Deleted in Task 3 |
| `fortran/setup.py`, `fortran/pyproject.toml`, `fortran/MANIFEST.in`, `fortran/README.rst` | Companion build |
| `fortran/climt_fortran/__init__.py` | Version, interface version, module table, `available()` |
| `fortran/climt_fortran/_*.pyx` | Moved from `climt/_components/...` |
| `fortran/_lib/` | Moved from `climt/_lib/` |
| `fortran/tests/test_climt_fortran.py` | Companion's own tests |
| `tests/conftest.py` | `FortranUnavailableError` → skip (or error under `CLIMT_REQUIRE_FORTRAN=1`) |
| `tests/test_fortran_loader.py`, `tests/test_fortran_skip.py`, `tests/fortran_skip_case.py`, `tests/test_list_components.py` (new); `tests/test_pure_wheel_build.py`, `tests/test_fortran_degradation.py` (extended) | Tests |
| `tox.ini`, `.github/workflows/python-package.yml`, `.github/workflows/release_climt.yml`, `.github/workflows/release_climt_fortran.yml` (new) | CI and release |
| `README.rst`, `docs/get-started/installation.qmd`, `CONTRIBUTING.rst`, `docs/user-guide/contributing.qmd`, `climt_release_checks.txt`, `HISTORY.rst`, `Makefile`, `setup.cfg`, `MANIFEST.in` | Docs and housekeeping |

---

### Task 1: The loader, `climt/_core/fortran.py`

Not wired into any component yet; nothing else changes behaviour in this task.

**Files:**
- Create: `climt/_core/fortran.py`
- Test: `tests/test_fortran_loader.py`

**Interfaces:**
- Consumes: nothing.
- Produces:

```python
EXPECTED_INTERFACE_VERSION = 1
MODULES = ("_dcmip", "_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw", "_simple_physics")

class FortranUnavailableError(ImportError): ...

def fortran_status() -> dict:
    """Keys: 'state', 'installed' (bool), 'version' (str | None),
    'interface_version' (int | None), 'modules' ({name: bool} for every name in
    MODULES), 'detail' (str, '' when available)."""

def load_extension(name: str):
    """Return the imported module climt_fortran.<name>, or None. Never raises
    for a missing, stub, mismatched or broken companion. Raises ValueError for
    a name not in MODULES."""

def unavailable_error(component: str, alternative: str | None = None) -> FortranUnavailableError:
    """Build (do not raise) the error for a component whose extension is missing."""

def _reset_cache() -> None:
    """Forget cached import results. For tests."""
```

`state` is one of:

| state | when |
|---|---|
| `"missing"` | `import climt_fortran` raises `ImportError` |
| `"interface_mismatch"` | `climt_fortran.INTERFACE_VERSION != EXPECTED_INTERFACE_VERSION` (or the attribute is absent); every entry of `modules` is `False` |
| `"stub"` | package imports, interface matches, and no module in `MODULES` can be found |
| `"broken"` | at least one module was found but importing it raised `ImportError`; `detail` holds `"<module>: <exception text>"` for each failure |
| `"available"` | all five modules import |

A companion in which some modules import and the rest are simply absent is `"broken"` with the absent ones named in `detail`.

Import results are cached per process (a failing binary must not be re-imported on every call). `version` comes from `climt_fortran.__version__`.

The error message must contain: the component name; a sentence for the state; and the fix command, exactly one of:

| state | command in the message |
|---|---|
| `missing` | `pip install climt-fortran` |
| `stub` | `pip install --no-binary climt-fortran climt-fortran` (and the words "needs gfortran and a C compiler") |
| `interface_mismatch` | `pip install -U climt climt-fortran` |
| `broken` | `pip install --force-reinstall climt-fortran`, preceded by the `detail` text |

When `alternative` is given the message ends with `Use <alternative> instead.`

- [ ] **Step 1: Write the failing tests**

`tests/test_fortran_loader.py`. A fake companion is a real package written to `tmp_path` and put first on `sys.path`, so the loader is exercised through the ordinary import system whatever its implementation. A module that must fail to load is a file that raises `ImportError`, as a binary built against another NumPy does. Every test starts and ends with `fortran._reset_cache()` and with no `climt_fortran*` entry in `sys.modules`.

```python
import sys

import pytest

from climt._core import fortran

ABI_MESSAGE = "numpy.core.multiarray failed to import"


@pytest.fixture(autouse=True)
def isolated(monkeypatch):
    for name in [n for n in sys.modules if n == "climt_fortran" or n.startswith("climt_fortran.")]:
        monkeypatch.delitem(sys.modules, name)
    fortran._reset_cache()
    yield
    for name in [n for n in sys.modules if n == "climt_fortran" or n.startswith("climt_fortran.")]:
        sys.modules.pop(name, None)
    fortran._reset_cache()


def install_fake(monkeypatch, tmp_path, interface=1, present=fortran.MODULES, failing=()):
    package = tmp_path / "climt_fortran"
    package.mkdir()
    lines = ['__version__ = "9.9.9"']
    if interface is not None:
        lines.append("INTERFACE_VERSION = %d" % interface)
    (package / "__init__.py").write_text("\n".join(lines) + "\n")
    for name in present:
        (package / (name + ".py")).write_text("LOADED = True\n")
    for name in failing:
        (package / (name + ".py")).write_text("raise ImportError(%r)\n" % ABI_MESSAGE)
    monkeypatch.syspath_prepend(str(tmp_path))


def test_missing(monkeypatch):
    monkeypatch.setitem(sys.modules, "climt_fortran", None)
    status = fortran.fortran_status()
    assert status["state"] == "missing" and status["installed"] is False
    assert status["version"] is None and status["interface_version"] is None
    assert status["modules"] == {name: False for name in fortran.MODULES}
    assert fortran.load_extension("_rrtmg_lw") is None
    message = str(fortran.unavailable_error("RRTMGLongwave", "CorkLongwaveRadiation"))
    assert "RRTMGLongwave" in message and "pip install climt-fortran" in message
    assert message.rstrip().endswith("Use CorkLongwaveRadiation instead.")


def test_available(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path)
    status = fortran.fortran_status()
    assert status == {
        "state": "available", "installed": True, "version": "9.9.9",
        "interface_version": 1, "detail": "",
        "modules": {name: True for name in fortran.MODULES},
    }
    assert fortran.load_extension("_dcmip").LOADED is True


def test_stub(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path, present=())
    assert fortran.fortran_status()["state"] == "stub"
    message = str(fortran.unavailable_error("SimplePhysics"))
    assert "pip install --no-binary climt-fortran climt-fortran" in message
    assert "gfortran" in message and "instead" not in message


@pytest.mark.parametrize("interface", [2, 0, None])
def test_interface_mismatch(monkeypatch, tmp_path, interface):
    install_fake(monkeypatch, tmp_path, interface=interface)
    status = fortran.fortran_status()
    assert status["state"] == "interface_mismatch"
    assert not any(status["modules"].values())
    assert fortran.load_extension("_rrtmg_sw") is None
    assert "pip install -U climt climt-fortran" in str(fortran.unavailable_error("RRTMGShortwave"))


def test_broken_binary_is_reported_not_raised(monkeypatch, tmp_path):
    present = tuple(n for n in fortran.MODULES if n != "_rrtmg_lw")
    install_fake(monkeypatch, tmp_path, present=present, failing=("_rrtmg_lw",))
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert status["modules"]["_rrtmg_lw"] is False and status["modules"]["_dcmip"] is True
    assert "_rrtmg_lw" in status["detail"] and ABI_MESSAGE in status["detail"]
    assert fortran.load_extension("_rrtmg_lw") is None
    assert fortran.load_extension("_dcmip") is not None
    message = str(fortran.unavailable_error("RRTMGLongwave"))
    assert ABI_MESSAGE in message and "pip install --force-reinstall climt-fortran" in message


def test_partly_built_companion_is_broken(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path, present=("_dcmip",))
    status = fortran.fortran_status()
    assert status["state"] == "broken"
    assert "_rrtmg_sw" in status["detail"]


def test_error_type_and_unknown_name(monkeypatch):
    monkeypatch.setitem(sys.modules, "climt_fortran", None)
    assert isinstance(fortran.unavailable_error("X"), ImportError)
    with pytest.raises(ValueError):
        fortran.load_extension("_not_a_module")


def test_results_are_cached(monkeypatch, tmp_path):
    install_fake(monkeypatch, tmp_path)
    assert fortran.fortran_status()["state"] == "available"
    monkeypatch.setitem(sys.modules, "climt_fortran", None)
    assert fortran.fortran_status()["state"] == "available"      # cached
    fortran._reset_cache()
    assert fortran.fortran_status()["state"] == "missing"
```

- [ ] **Step 2: Run and confirm failure**

Run: `python -m pytest tests/test_fortran_loader.py -q`
Expected: collection error, `cannot import name 'fortran' from 'climt._core'`.

- [ ] **Step 3: Implement `climt/_core/fortran.py`**

Standard library only (`importlib`). One private function does the imports once and stores the result; the three public functions read from it. Do not export anything from `climt/_core/__init__.py` in this task.

- [ ] **Step 4: Run the tests and flake8**

Run: `python -m pytest tests/test_fortran_loader.py -q && flake8 climt --config=tox.ini`
Expected: all pass; flake8 silent.

- [ ] **Step 5: Commit**

```bash
git add climt/_core/fortran.py tests/test_fortran_loader.py
git commit -m "feat: loader for the climt-fortran companion package"
```

---

### Task 2: Create `fortran/` by moving the compiled code, and wire the components to the loader

After this task the root `setup.py` builds nothing compiled, the companion builds from `fortran/`, and the full test suite passes with the companion installed from source.

**Files:**
- Move (`git mv`): `climt/_lib` → `fortran/_lib`; `climt/_components/dcmip/_dcmip.pyx`, `climt/_components/emanuel/_emanuel_convection.pyx`, `climt/_components/rrtmg/lw/_rrtmg_lw.pyx`, `climt/_components/rrtmg/sw/_rrtmg_sw.pyx`, `climt/_components/simple_physics/_simple_physics.pyx` → `fortran/climt_fortran/`
- Delete (`git rm`): `climt/_components/_berger_solar_insolation.pyx`
- Create: `fortran/setup.py`, `fortran/pyproject.toml`, `fortran/MANIFEST.in`, `fortran/README.rst`, `fortran/climt_fortran/__init__.py`, `fortran/tests/test_climt_fortran.py`, `tests/test_fortran_skip.py`, `tests/fortran_skip_case.py`
- Modify: `fortran/_lib/Makefile` (the five `touch` paths), the five `component.py` files, `climt/__init__.py`, `climt/_core/__init__.py`, root `setup.py`, root `MANIFEST.in`, `setup.cfg` (flake8 `exclude`), `Makefile`, `tox.ini`, `tests/conftest.py`, `tests/test_fortran_degradation.py`

**Interfaces:**
- Consumes: Task 1's `load_extension`, `fortran_status`, `unavailable_error`, `FortranUnavailableError`, `MODULES`.
- Produces:
  - `climt_fortran.__version__ == "1.0.0"`, `climt_fortran.INTERFACE_VERSION == 1`, `climt_fortran.MODULES` (same tuple as the loader's), `climt_fortran.available() -> dict[str, bool]` (uses `importlib.util.find_spec`; imports no extension).
  - `climt.FortranUnavailableError` (re-export), `climt.has_fortran_extensions()` returning `fortran_status()["state"] == "available"`.
  - Environment variables: `CLIMT_FORTRAN_STUB=1` (build the stub), `CLIMT_REQUIRE_FORTRAN=1` (tests: unavailability is an error).

- [ ] **Step 1: Baseline — prove the toolchain builds the extensions before moving anything**

```bash
export SDKROOT="$(xcrun --show-sdk-path 2>/dev/null)"   # macOS only; harmless elsewhere
pip install -e . && python setup.py build_ext --inplace
python -c "import climt; assert climt.has_fortran_extensions()"
python -m pytest tests -q -m "not slow" -x -q 2>&1 | tail -5
```

Expected: the assertion passes; note the pass/fail/skip counts in your report as the baseline. If the build fails here, stop: report BLOCKED with the compiler output.

Then remove the in-place products so later steps cannot pick them up by accident: `find climt -name '*.so' -delete; find climt/_components -name '*.c' -delete`.

- [ ] **Step 2: Write the companion's tests (failing)**

`fortran/tests/test_climt_fortran.py`:

```python
import importlib
import importlib.metadata
import os

import pytest

import climt_fortran

EXPECTED = ("_dcmip", "_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw", "_simple_physics")
REQUIRED = os.environ.get("CLIMT_REQUIRE_FORTRAN") == "1"
ENTRY_POINTS = {
    "_dcmip": ("get_baroclinic_wave_ics", "get_tropical_cyclone_ics"),
    "_emanuel_convection": ("init_emanuel_convection", "convect"),
    "_rrtmg_lw": ("set_constants", "initialise_rrtm_radiation", "rrtm_calculate_longwave_fluxes"),
    "_rrtmg_sw": ("set_constants", "initialise_rrtm_radiation", "rrtm_calculate_shortwave_fluxes"),
    "_simple_physics": ("init_simple_physics",),
}


def test_metadata():
    assert climt_fortran.INTERFACE_VERSION == 1
    assert climt_fortran.MODULES == EXPECTED
    assert climt_fortran.__version__ == importlib.metadata.version("climt-fortran") == "1.0.0"


def test_available_reports_every_module_without_importing():
    result = climt_fortran.available()
    assert tuple(sorted(result)) == EXPECTED
    assert all(isinstance(v, bool) for v in result.values())
    if REQUIRED:
        assert all(result.values()), result


@pytest.mark.parametrize("name", EXPECTED)
def test_extension_exposes_the_functions_climt_calls(name):
    if not climt_fortran.available()[name]:
        if REQUIRED:
            pytest.fail(name + " is not built")
        pytest.skip("stub install")
    module = importlib.import_module("climt_fortran." + name)
    for function in ENTRY_POINTS[name]:
        assert callable(getattr(module, function))
```

And one build test, in the same file, for the stub:

```python
def test_stub_wheel_has_no_binaries(tmp_path):
    import glob, subprocess, sys, zipfile
    from pathlib import Path
    project = Path(__file__).resolve().parent.parent
    env = dict(os.environ, CLIMT_FORTRAN_STUB="1")
    subprocess.run([sys.executable, "-m", "pip", "wheel", str(project), "--no-deps", "-w", str(tmp_path)],
                   env=env, check=True)
    wheels = glob.glob(str(tmp_path / "climt_fortran-1.0.0-py3-none-any.whl"))
    assert len(wheels) == 1
    names = zipfile.ZipFile(wheels[0]).namelist()
    assert "climt_fortran/__init__.py" in names
    assert not [n for n in names if n.endswith((".so", ".pyd", ".pyx", ".c", ".a", ".f90", ".F90"))]
```

- [ ] **Step 3: Move the files and write the companion**

```bash
mkdir -p fortran/climt_fortran fortran/tests
git mv climt/_lib fortran/_lib
git mv climt/_components/dcmip/_dcmip.pyx climt/_components/emanuel/_emanuel_convection.pyx \
       climt/_components/rrtmg/lw/_rrtmg_lw.pyx climt/_components/rrtmg/sw/_rrtmg_sw.pyx \
       climt/_components/simple_physics/_simple_physics.pyx fortran/climt_fortran/
git rm climt/_components/_berger_solar_insolation.pyx
```

Then:

- `fortran/_lib/Makefile`: each `touch ../../_components/<...>/<name>.pyx` becomes `touch ../../climt_fortran/<name>.pyx`. Nothing else in `_lib` changes.
- `fortran/climt_fortran/__init__.py`: the four names in Interfaces. No other imports than `importlib.util`.
- `fortran/pyproject.toml`: `[build-system] requires = ["setuptools", "wheel", "Cython", "numpy>=2.0.0"]`, `build-backend = "setuptools.build_meta"`.
- `fortran/setup.py`: take the build logic from the root `setup.py` as it stands (compiler guessing on Darwin, `FC`/`CC`/`CLIMT_OPT_FLAGS`/`FFLAGS`/`CFLAGS`, `LDSHARED` on Darwin, `make CLIMT_ARCH=<system>` run in `_lib`, the custom `build_ext` and `bdist_wheel` classes, and the five `Extension` definitions with their libraries, OpenMP flags for the two RRTMG extensions, and static-library link arguments). Changes: extension names are `climt_fortran.<name>` with sources `climt_fortran/<name>.pyx`; the library directory is `_lib/<system>/` and the include directory `_lib/include`; the Berger extension, the `READTHEDOCS` branches and the `pip_main(["install", ...])` fallbacks are dropped (build requirements now come from `pyproject.toml`). Metadata: `name="climt-fortran"`, `version="1.0.0"`, `packages=["climt_fortran"]`, `install_requires=["numpy>=1.16.0"]`, `python_requires=">=3.10"`, `license="BSD license"`, `zip_safe=False`, description "Compiled Fortran extensions for climt". When `os.environ.get("CLIMT_FORTRAN_STUB") == "1"`: do not import Cython or numpy, pass no `ext_modules` and no `cmdclass`, and `install_requires=[]`, so the wheel is `py3-none-any` and contains only `__init__.py`.
- `fortran/MANIFEST.in`: include `README.rst`, `recursive-include climt_fortran *.py *.pyx`, `recursive-include _lib *`, and the same prunes and excludes the root `MANIFEST.in` has for `_lib` build products (`*.o`, `*.mod`, `*.a`, `*.c`, the `Open*`, `fftw*`, `shtns*`, `lib`, `share`, `include`, `bin` directories), with paths relative to `fortran/`.
- `fortran/README.rst`: what the package is, that users normally get it through `pip install climt`, how to build from source, the stub, and that the Fortran sources keep their original licences.

- [ ] **Step 4: Strip the root build to pure Python**

Root `setup.py`: delete everything compiled (the `PURE_PYTHON` switch, the pip fallbacks, compiler guessing, `build_libraries`, both custom command classes, `ext_modules`, `cmdclass`, `include_dirs`, the numpy import) and remove `"cython>=0.25"` from `requirements`. It keeps metadata, `packages`, `package_data` and `install_requires` only. (Task 3 replaces the file with `pyproject.toml`.)

Root `pyproject.toml`: `requires = ["setuptools", "wheel"]`.

Root `MANIFEST.in`: remove `recursive-include climt/_lib *`, every `prune`/`exclude` line under `climt/_lib`, and `*.so *.dylib *.pyd` from the `recursive-include climt` line.

`setup.cfg`: flake8 `exclude = docs,build,fortran/_lib`. `Makefile`: `clean-cython` and `clean-libs` point at `fortran/climt_fortran` and `fortran/_lib`.

`tox.ini`, both `[testenv]` and `[testenv:cov]`: replace `pip install -e .` and `python setup.py build_ext --inplace` with `pip install -e ./fortran` then `pip install -e .`. (Task 5 restructures tox properly; this keeps it working meanwhile.)

- [ ] **Step 5: Wire the components**

In each of the five `component.py` files replace the `try: from . import _x ... except ImportError` block with the loader, keeping the existing local names so no call site changes:

```python
from ..._core.fortran import load_extension, unavailable_error

_emanuel_convection = load_extension("_emanuel_convection")
_FORTRAN_AVAILABLE = _emanuel_convection is not None
```

(`simple_physics` keeps its local name `phys`; the RRTMG files use `...._core` with the right number of dots.) Replace each `raise ImportError(...)` guarded by `_FORTRAN_AVAILABLE` with `raise unavailable_error(<component>, <alternative>)`:

| component | alternative |
|---|---|
| `EmanuelConvection` | `"EmanuelConvectionPython"` |
| `RRTMGLongwave` | `"CorkLongwaveRadiation"` |
| `RRTMGShortwave` | `"CorkShortwaveRadiation"` |
| `SimplePhysics` | `None` |
| `DcmipInitialConditions` | `None` |

`climt/_core/__init__.py`: export `FortranUnavailableError`, `fortran_status`. `climt/__init__.py`: import and list `FortranUnavailableError` in `__all__`; `has_fortran_extensions()` becomes `return fortran_status()["state"] == "available"`.

- [ ] **Step 6: Skip-or-fail hook and its tests**

`tests/conftest.py`: a `FortranUnavailableError` raised during a test's setup or call phase becomes a skip whose reason is the error's message, unless `os.environ.get("CLIMT_REQUIRE_FORTRAN") == "1"`, in which case it is left alone. Also register the hook's behaviour for fixtures (setup phase), since several test classes build components in fixtures.

`tests/fortran_skip_case.py` (not collected by default because of its name):

```python
import pytest
from climt import FortranUnavailableError


@pytest.fixture
def component():
    raise FortranUnavailableError("fixture needs the extensions")


def test_raises_in_body():
    raise FortranUnavailableError("body needs the extensions")


def test_raises_in_fixture(component):
    pass


def test_other_import_errors_still_fail():
    raise ImportError("an unrelated import problem")
```

`tests/test_fortran_skip.py`:

```python
import os, subprocess, sys
from pathlib import Path

CASE = str(Path(__file__).with_name("fortran_skip_case.py"))


def run(require):
    env = dict(os.environ)
    env.pop("CLIMT_REQUIRE_FORTRAN", None)
    if require:
        env["CLIMT_REQUIRE_FORTRAN"] = "1"
    return subprocess.run([sys.executable, "-m", "pytest", CASE, "-q", "-p", "no:cacheprovider"],
                          capture_output=True, text=True, env=env)


def test_unavailable_is_skipped_by_default():
    result = run(require=False)
    assert "2 skipped" in result.stdout and "1 failed" in result.stdout, result.stdout


def test_unavailable_is_an_error_when_required():
    result = run(require=True)
    assert "skipped" not in result.stdout, result.stdout
    assert "FortranUnavailableError" in result.stdout
    assert result.returncode != 0
```

Extend `tests/test_fortran_degradation.py`:

```python
def test_climt_works_without_the_companion():
    code = (
        "import sys; sys.modules['climt_fortran'] = None\n"
        "import climt\n"
        "assert climt.has_fortran_extensions() is False\n"
        "try:\n"
        "    climt.RRTMGLongwave()\n"
        "except climt.FortranUnavailableError as error:\n"
        "    assert 'pip install climt-fortran' in str(error)\n"
        "else:\n"
        "    raise SystemExit('RRTMGLongwave did not raise')\n"
        "climt.GrayLongwaveRadiation()\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == "" and result.stderr == ""


def test_components_never_import_extensions_from_the_old_location():
    import pathlib
    root = pathlib.Path(climt.__file__).parent
    offenders = [
        str(path) for path in root.rglob("*.py")
        if path.name != "fortran.py" and "climt_fortran" in path.read_text()
    ]
    offenders += [
        str(path) for path in (root / "_components").rglob("component.py")
        if "from . import _" in path.read_text()
    ]
    assert offenders == []
```

- [ ] **Step 7: Build from source and run everything**

```bash
export SDKROOT="$(xcrun --show-sdk-path 2>/dev/null)"
pip install -e ./fortran && pip install -e .
CLIMT_REQUIRE_FORTRAN=1 python -m pytest fortran/tests tests -q -m "not slow" 2>&1 | tail -5
```

Expected: `fortran/tests` all pass; `tests` gives the same pass count as the Step 1 baseline plus the new tests, with no new failures and no test skipped that passed in the baseline.

- [ ] **Step 8: Run without the binaries**

```bash
pip uninstall -y climt-fortran && find fortran/climt_fortran -name '*.so' -delete
CLIMT_FORTRAN_STUB=1 pip install ./fortran
python -c "import climt; s = climt._core.fortran.fortran_status(); assert s['state'] == 'stub', s"
python -m pytest fortran/tests tests -q -m "not slow" 2>&1 | tail -5
```

Expected: no failures; the Fortran-backed tests are reported as skipped with the `climt-fortran` message as the reason. A test module that fails at *collection* because it builds a Fortran component at import time gets, at the top of the module, `pytestmark = pytest.mark.skipif(not climt.has_fortran_extensions(), reason="needs climt-fortran binaries")` and its module-level construction moved into a fixture or guarded; list every module you had to change in your report.

Finish by restoring the source build: `pip uninstall -y climt-fortran && pip install -e ./fortran`.

- [ ] **Step 9: flake8 and commit**

```bash
flake8 climt fortran/climt_fortran --config=tox.ini
git add -A fortran climt tests setup.py pyproject.toml MANIFEST.in setup.cfg Makefile tox.ini
git commit -m "refactor: move compiled extensions into the climt-fortran companion package"
```

---

### Task 3: Declarative packaging for the core

**Files:**
- Modify: `pyproject.toml`, `setup.cfg`, `Makefile`, `tests/test_pure_wheel_build.py`
- Delete: `setup.py`

**Interfaces:**
- Consumes: the trimmed `setup.py` from Task 2 (its metadata is the source for `pyproject.toml`).
- Produces: a `climt` wheel whose metadata later tasks and CI rely on: `Requires-Dist: climt-fortran<2,>=1.0; sys_platform != "emscripten"`, `Requires-Python: >=3.10`, `Provides-Extra: numba`.

- [ ] **Step 1: Extend the wheel test (failing)**

Add to `tests/test_pure_wheel_build.py`; factor the existing wheel build into a module-scoped fixture `built_wheel(tmp_path_factory)` that builds once **without** `CLIMT_PURE_PYTHON` and returns the wheel path. Keep the existing test, pointed at a second build **with** `CLIMT_PURE_PYTHON=1` (it must still produce the same pure wheel, because the docs' boot script sets that variable).

```python
FORBIDDEN_SUFFIXES = (".so", ".pyd", ".pyx", ".c", ".a", ".f90", ".F90", ".f", ".mod")


def test_wheel_is_pure_and_clean(built_wheel):
    assert built_wheel.name.endswith("-py3-none-any.whl")
    names = zipfile.ZipFile(built_wheel).namelist()
    assert not [n for n in names if n.endswith(FORBIDDEN_SUFFIXES)]
    assert not [n for n in names if n.startswith(("fortran/", "climt_fortran/")) or "/_lib/" in n]


def test_wheel_metadata(built_wheel):
    archive = zipfile.ZipFile(built_wheel)
    metadata = archive.read(next(n for n in archive.namelist() if n.endswith(".dist-info/METADATA"))).decode()
    requires = [line.split(":", 1)[1].strip() for line in metadata.splitlines() if line.startswith("Requires-Dist:")]
    fortran = [r for r in requires if r.startswith("climt-fortran")]
    assert len(fortran) == 1
    assert ">=1.0" in fortran[0] and "<2" in fortran[0]
    assert "sys_platform" in fortran[0] and "emscripten" in fortran[0] and "extra" not in fortran[0]
    assert not [r for r in requires if r.lower().startswith("cython")]
    assert [r for r in requires if r.startswith("numba") and 'extra == "numba"' in r.replace("'", '"')]
    assert "Requires-Python: >=3.10" in metadata
    assert "Provides-Extra: numba" in metadata
    assert "Programming Language :: Python :: 2" not in metadata


def test_wheel_keeps_every_data_file(built_wheel):
    names = zipfile.ZipFile(built_wheel).namelist()
    for needle in (
        "climt/_data/ozone_profile.npy",
        "climt/_data/atmospheric_properties/earth.toml",
        "climt/_data/cork/correlated_k/earth_low_res_lw.nc",
        "climt/_data/cork/correlated_k/earth_low_res_lw.npz",
        "climt/_data/land_mask/earth_landmask_2deg.nc",
        "climt/_data/data_ocean/earth_sst_climatology_1deg.nc",
        "climt/_data/topography/earth_topography_2deg.nc",
        "climt/_core/fortran.py",
    ):
        assert any(n.endswith(needle) for n in names), needle
```

Before deleting `setup.py`, record the data-file list of the wheel it builds (`unzip -l <wheel> | grep _data | wc -l`) and state in your report that the new wheel has the same count.

- [ ] **Step 2: Run and confirm the metadata test fails**

Run: `python -m pytest tests/test_pure_wheel_build.py -q`
Expected: `test_wheel_metadata` fails (no `climt-fortran` requirement, no `Requires-Python`).

- [ ] **Step 3: Write `pyproject.toml` and delete `setup.py`**

`[build-system]`: `requires = ["setuptools>=68"]`, `build-backend = "setuptools.build_meta"`.

`[project]`: `name`, `version` (the current value in `setup.py`; unchanged), `description`, `authors`, `keywords`, `urls` carried over from `setup.py`; `readme = "README.rst"`; `requires-python = ">=3.10"`; `license = {text = "BSD license"}`; classifiers: Development Status, Intended Audience, the BSD licence classifier, Natural Language, and `Programming Language :: Python :: 3`, `3.10`, `3.11`, `3.12`.

`dependencies`: `numpy>=1.16.0`, `pint>=0.7.0`, `unyt`, `xarray>=0.8.0`, `sympl>=0.5.0`, `importlib_resources`, and `climt-fortran>=1.0,<2; sys_platform != "emscripten"`.

`[project.optional-dependencies]`: `numba = ["numba>=0.59.0"]`, plus anything carried over (Global Constraints).

`[tool.setuptools.packages.find]`: `include = ["climt", "climt.*"]`. `[tool.setuptools.package-data]`: `climt = [...]` with the same seven globs `setup.py` has.

`setup.cfg`: the `[bumpversion:file:setup.py]` section becomes `[bumpversion:file:pyproject.toml]` with `search = version = "{current_version}"` and `replace = version = "{new_version}"`.

`Makefile`: `release`, `dist` and `install` use `python -m build` and `pip install .`; the `upload` commands are removed (publishing is done by the release workflows).

- [ ] **Step 4: Run the tests**

Run: `pip install -e . && python -m pytest tests/test_pure_wheel_build.py tests/test_fortran_degradation.py -q`
Expected: all pass. Then `python -c "import climt, importlib.metadata as m; assert m.version('climt') == climt.__version__"`.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml setup.cfg Makefile tests/test_pure_wheel_build.py && git rm -q setup.py
git commit -m "build: declarative pyproject for a pure-Python climt"
```

---

### Task 4: `climt.list_components()`

**Files:**
- Create: `climt/_core/registry.py`
- Modify: `climt/_core/__init__.py`, `climt/__init__.py`
- Test: `tests/test_list_components.py`

**Interfaces:**
- Consumes: `fortran_status()` (Task 1); `climt._components.__all__` (a tuple of classes); `climt._core.backend.HAS_NUMBA`.
- Produces:

```python
ENTRY_POINT_GROUP = "climt.components"

class ComponentInfo(NamedTuple):
    name: str        # class name, or the entry point's name
    kind: str        # "Stepper" | "TendencyComponent" | "ImplicitTendencyComponent" | "DiagnosticComponent" | "other" | "external"
    available: bool
    requires: str    # "" when available; otherwise what to install, e.g. "climt-fortran"
    source: str      # "climt", or "<distribution> (<entry point value>)"

class ComponentList(list):
    """A list of ComponentInfo. repr() is an aligned text table followed by one
    line about numba."""

def list_components() -> ComponentList
```

Rules:
- Built-in records come from `climt._components.__all__`, sorted by name. `kind` is the first of the four sympl base classes found in the class's MRO, else `"other"`.
- A built-in is Fortran-backed if it is a key of the module-level table `FORTRAN_BACKED = {"DcmipInitialConditions": "_dcmip", "EmanuelConvection": "_emanuel_convection", "RRTMGLongwave": "_rrtmg_lw", "RRTMGShortwave": "_rrtmg_sw", "SimplePhysics": "_simple_physics"}`. Its `available` is `fortran_status()["modules"][<module>]`, and `requires` is `"climt-fortran"` when unavailable.
- External records come from `importlib.metadata.entry_points(group=ENTRY_POINT_GROUP)`, after the built-ins, sorted by name. Entry points are **not loaded**: `kind="external"`, `available=True`, `requires=""`, `source="<ep.dist.name> (<ep.value>)"` (just `ep.value` if the distribution is unknown). An entry point whose name equals a built-in's is listed as well, not merged.
- The numba line reads `numba: installed` or `numba: not installed (kernels run as plain Python; pip install "climt[numba]")`.

- [ ] **Step 1: Write the failing tests**

```python
import inspect
import sys

import climt
from climt._core import fortran, registry


def by_name(components):
    return {c.name: c for c in components}


def test_every_exported_component_is_listed():
    listed = by_name(climt.list_components())
    for cls in climt._components.__all__:
        assert listed[cls.__name__].source == "climt"
    assert listed["GrayLongwaveRadiation"] == registry.ComponentInfo(
        "GrayLongwaveRadiation", "TendencyComponent", True, "", "climt")


def test_fortran_backed_table_matches_the_code():
    """A component is in FORTRAN_BACKED exactly when its module asks the loader for an extension."""
    for cls in climt._components.__all__:
        uses_loader = "load_extension(" in inspect.getsource(sys.modules[cls.__module__])
        assert uses_loader == (cls.__name__ in registry.FORTRAN_BACKED), cls.__name__
    assert set(registry.FORTRAN_BACKED.values()) == set(fortran.MODULES)


def test_unavailable_components_say_what_they_need(monkeypatch):
    status = dict(fortran.fortran_status(), state="stub", modules={m: False for m in fortran.MODULES})
    monkeypatch.setattr(registry, "fortran_status", lambda: status)
    listed = by_name(climt.list_components())
    for name in registry.FORTRAN_BACKED:
        assert listed[name].available is False and listed[name].requires == "climt-fortran"
    assert listed["CorkLongwaveRadiation"].available is True


def test_entry_points_are_listed_without_being_loaded(monkeypatch):
    class FakeDist:
        name = "some-package"

    class FakeEntryPoint:
        name = "SomeComponent"
        value = "some_package.components:SomeComponent"
        dist = FakeDist()

        def load(self):
            raise AssertionError("entry points must not be loaded")

    monkeypatch.setattr(registry, "_entry_points", lambda: [FakeEntryPoint()])
    record = by_name(climt.list_components())["SomeComponent"]
    assert record == registry.ComponentInfo(
        "SomeComponent", "external", True, "", "some-package (some_package.components:SomeComponent)")


def test_repr_is_a_table_with_a_numba_line():
    text = repr(climt.list_components())
    lines = text.splitlines()
    assert lines[0].split()[:4] == ["name", "kind", "available", "requires"]
    assert any(line.startswith("RRTMGLongwave") for line in lines)
    assert lines[-1].startswith("numba: ")
```

(`registry._entry_points()` is the one private function that calls `importlib.metadata.entry_points`.)

- [ ] **Step 2: Run and confirm failure**

Run: `python -m pytest tests/test_list_components.py -q`
Expected: `AttributeError: module 'climt' has no attribute 'list_components'`.

- [ ] **Step 3: Implement `climt/_core/registry.py` and export `list_components` from `climt._core` and `climt`**

`registry` imports `climt._components` inside `list_components()` (not at module top) to avoid a circular import.

- [ ] **Step 4: Run the tests, the import-silence test and flake8**

Run: `python -m pytest tests/test_list_components.py tests/test_fortran_degradation.py -q && flake8 climt --config=tox.ini`
Expected: all pass.

- [ ] **Step 5: If `climt/_core/versions.py` exists on this branch** (it comes from the development-environment work), add `"climt-fortran"` to its `REPORTED_PACKAGES` and extend its test's expected names. If the file does not exist, skip this step and say so in your report.

- [ ] **Step 6: Commit**

```bash
git add climt/_core/registry.py climt/_core/__init__.py climt/__init__.py tests/test_list_components.py
git commit -m "feat: climt.list_components() reports what is usable and what it needs"
```

---

### Task 5: tox, CI and release workflows

Workflows cannot be run locally. This task is verified by tox locally, by YAML syntax checks, and by a test that pins the facts the workflows must agree on.

**Files:**
- Modify: `tox.ini`, `.github/workflows/python-package.yml`, `.github/workflows/release_climt.yml`
- Create: `.github/workflows/release_climt_fortran.yml`, `tests/test_workflows.py`

**Interfaces:**
- Consumes: `CLIMT_FORTRAN_STUB`, `CLIMT_REQUIRE_FORTRAN` (Task 2); `fortran/tests` (Task 2).
- Produces: tox environments `py310`, `py311`, `py312` (with binaries) and `py310-pure`, `py311-pure`, `py312-pure` (stub); `flake8`; `cov`.

Required behaviour:

**tox.** A non-`pure` environment runs `pip install -e ./fortran`, `pip install -e .`, then `py.test fortran/tests tests --basetemp={envtmpdir} -v -m "not slow"` with `CLIMT_REQUIRE_FORTRAN=1`. A `pure` environment runs `pip install ./fortran` with `CLIMT_FORTRAN_STUB=1`, `pip install -e .`, asserts `python -c "import climt; assert not climt.has_fortran_extensions()"`, then the same `py.test` line without `CLIMT_REQUIRE_FORTRAN`. `py39` leaves the `envlist`. `flake8` lints `climt` and `fortran/climt_fortran`. `cov` follows the non-`pure` recipe.

**`python-package.yml`.** The matrix gains `mode: [full, pure]`; the tox environment is `py<ver>` for `full` and `py<ver>-pure` for `pure`. The flake8 step lints both packages. Everything else is unchanged.

**`release_climt.yml`.** Runs on a published release whose tag starts with `v`. One build job: `python -m build` (sdist and wheel), check exactly one `*-py3-none-any.whl` and one `.tar.gz`, attach the wheel to the GitHub release (the in-browser docs demo uses it), upload both as an artifact. The publish job (PyPI trusted publishing, environment `pypi`) keeps its `if: ${{ !github.event.release.prerelease }}` guard and now needs only the build job. No `cibuildwheel`, no QEMU.

**`release_climt_fortran.yml`.** Runs on a published release whose tag starts with `fortran-v`. Jobs:
- `build_wheels`: the matrix, runner pin (`macos-15`, with the existing comment explaining why) and environment of today's `build_wheels` job, with `python -m cibuildwheel fortran --output-dir dist` and, added, `CIBW_TEST_REQUIRES: pytest` and `CIBW_TEST_COMMAND: "CLIMT_REQUIRE_FORTRAN=1 pytest {package}/tests -k 'not stub_wheel'"`.
- `build_stub_and_sdist`: `CLIMT_FORTRAN_STUB=1 python -m build --wheel fortran` and `python -m build --sdist fortran`; check for exactly one `climt_fortran-*-py3-none-any.whl`.
- `publish`: `needs: [build_wheels, build_stub_and_sdist]`, not for prereleases, trusted publishing with environment `pypi` and URL `https://pypi.org/p/climt-fortran`. Because it needs both jobs, a version is never published with a stub but without its binaries.

- [ ] **Step 1: Write the pinning test (failing)**

`tests/test_workflows.py` (uses PyYAML, already in `requirements_dev.txt`; `pytest.importorskip("yaml")`):

```python
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
ROOT = Path(__file__).resolve().parent.parent
WORKFLOWS = ROOT / ".github" / "workflows"


def load(name):
    return yaml.safe_load((WORKFLOWS / name).read_text())


def text(name):
    return (WORKFLOWS / name).read_text()


def test_ci_runs_both_modes():
    matrix = load("python-package.yml")["jobs"]["build"]["strategy"]["matrix"]
    assert matrix["mode"] == ["full", "pure"]
    assert matrix["python-version"] == ["3.10", "3.11", "3.12"]


def test_core_release_never_compiles():
    body = text("release_climt.yml")
    assert "cibuildwheel" not in body and "qemu" not in body.lower()
    assert "startsWith(github.event.release.tag_name, 'v')" in body
    assert "python -m build" in body


def test_fortran_release_publishes_only_complete_sets():
    workflow = load("release_climt_fortran.yml")
    assert set(workflow["jobs"]["publish"]["needs"]) == {"build_wheels", "build_stub_and_sdist"}
    body = text("release_climt_fortran.yml")
    assert "startsWith(github.event.release.tag_name, 'fortran-v')" in body
    assert "CLIMT_FORTRAN_STUB" in body and "CLIMT_REQUIRE_FORTRAN=1" in body
    assert "macos-15" in body and "cp310-* cp311-* cp312-*" in body
    assert "https://pypi.org/p/climt-fortran" in body


def test_tox_has_both_kinds_of_environment():
    tox = (ROOT / "tox.ini").read_text()
    for name in ("py310", "py311", "py312", "py310-pure", "py311-pure", "py312-pure"):
        assert name in tox
    assert "py39" not in tox
    assert "setup.py" not in tox
```

- [ ] **Step 2: Run and confirm failure**

Run: `python -m pytest tests/test_workflows.py -q`
Expected: failures (no `mode` in the matrix; no `release_climt_fortran.yml`).

- [ ] **Step 3: Write `tox.ini` and the three workflows as specified above**

- [ ] **Step 4: Verify**

```bash
python -m pytest tests/test_workflows.py -q
tox -e py311-pure 2>&1 | tail -5
export SDKROOT="$(xcrun --show-sdk-path 2>/dev/null)"; tox -e py311 2>&1 | tail -5
```

Expected: the pinning tests pass; both tox environments finish with no failures; the `py311` run reports no skipped Fortran tests. If `actionlint` is on the PATH, run it on the three workflow files and fix what it reports; if it is not, say so in your report.

- [ ] **Step 5: Commit**

```bash
git add tox.ini .github/workflows tests/test_workflows.py
git commit -m "ci: test with and without binaries; separate release workflows for climt and climt-fortran"
```

---

### Task 6: Documentation and housekeeping

**Files:**
- Modify: `docs/get-started/installation.qmd`, `README.rst`, `CONTRIBUTING.rst`, `docs/user-guide/contributing.qmd`, `climt_release_checks.txt`, `HISTORY.rst`, `docs/_includes/climt-live-boot.qmd`, `scripts/serve_wheel.py`
- Conditionally modify: `scripts/bootstrap_dev_env.sh`, `scripts/verify_dev_env.py`, `envs/climt-dev.in` (only if they exist on this branch)

**Interfaces:**
- Consumes: everything above.
- Produces: nothing other tasks use.

Content each file must carry:

- **`docs/get-started/installation.qmd`**: replace the section "Which wheel you get" with "What `pip install climt` installs": climt itself is pure Python; it brings in `climt-fortran`, which carries the Fortran-backed components (RRTMG, Fortran Emanuel, SimplePhysics, DCMIP initial conditions) as binaries on Linux x86-64 and macOS on Apple silicon (macOS 15 or newer) for CPython 3.10–3.12, and as an empty stub everywhere else. Show `climt.has_fortran_extensions()` and `climt.list_components()` (with a few lines of its output). Add: removing a faulty binary (`pip uninstall climt-fortran`; climt keeps working), building the binaries from source (`pip install --no-binary climt-fortran climt-fortran`, needs gfortran and a C compiler), and `pip install "climt[numba]"` with one sentence on why. Replace the "from source" instructions that use `python setup.py install` with `pip install ./fortran` (optional) followed by `pip install -e .`. Keep the Windows note, reworded for the stub.
- **`README.rst`**: the install section gains two sentences on `climt-fortran` and `list_components()`.
- **`CONTRIBUTING.rst`** and **`docs/user-guide/contributing.qmd`** (kept in step): development setup is `pip install -e ./fortran` (or `CLIMT_FORTRAN_STUB=1 pip install ./fortran` without a Fortran toolchain) then `pip install -e .`; tests are `pytest` or `tox -e py311` / `tox -e py311-pure`; on macOS with Homebrew gcc, `export SDKROOT="$(xcrun --show-sdk-path)"` first. Remove `python setup.py develop` and `python setup.py test`.
- **`climt_release_checks.txt`**: two checklists. Core: HISTORY, bumpversion, `python -m build`, `twine check dist/*`, push, publish a GitHub release tagged `vX.Y.Z`. Companion: only when Fortran, Cython, the supported Python versions or the NumPy ABI change; edit the version in `fortran/setup.py` and `fortran/climt_fortran/__init__.py`; bump `INTERFACE_VERSION` and the core's requirement range if a function the core calls changed; publish a GitHub release tagged `fortran-vX.Y.Z`; **the companion is released before any core release that needs it**. One-time setup: reserve `climt-fortran` on PyPI and add a trusted publisher for workflow `release_climt_fortran.yml`.
- **`HISTORY.rst`**, under `Unreleased`: climt is now a pure-Python package; the compiled extensions moved to `climt-fortran`, installed automatically; `climt.list_components()`; `FortranUnavailableError`; the `numba` extra; Python 3.10+ declared; `setup.py` removed.
- **`docs/_includes/climt-live-boot.qmd`** and **`scripts/serve_wheel.py`**: the documented build command drops `CLIMT_PURE_PYTHON=1` (it is no longer needed; the variable is still accepted). Do not change how the wheel is served or loaded.

- [ ] **Step 1: Write the documentation changes above**

- [ ] **Step 2: If `scripts/bootstrap_dev_env.sh` exists on this branch**, change its climt install step to install `"$CLIMT_DIR/fortran"` editable first (keeping its existing behaviour of reporting, not aborting, when the Fortran build fails; on failure it installs the stub with `CLIMT_FORTRAN_STUB=1` instead) and then `"$CLIMT_DIR"` editable; update `tests/test_dev_env_scripts.py` to match; `scripts/verify_dev_env.py`'s Fortran check keeps using `climt.has_fortran_extensions()`. If the scripts do not exist, skip this step and say so in your report.

- [ ] **Step 3: Check nothing still points at the old layout**

```bash
git grep -n "climt/_lib\|setup.py build_ext\|setup.py develop\|setup.py install\|setup.py test\|CLIMT_PURE_PYTHON" -- . ':!docs/superpowers' ':!HISTORY.rst' ':!fortran/_lib'
```

Expected: the only hits are `tests/test_pure_wheel_build.py` (which deliberately sets `CLIMT_PURE_PYTHON`) and `appveyor.yml` (a retired Windows configuration this plan does not touch). Fix anything else.

- [ ] **Step 4: Build the docs if Quarto is installed**

Run: `quarto render docs/get-started/installation.qmd` (skip, and say so, if `quarto` is not on the PATH).
Expected: renders without error.

- [ ] **Step 5: Full verification and commit**

```bash
flake8 climt fortran/climt_fortran --config=tox.ini
CLIMT_REQUIRE_FORTRAN=1 python -m pytest fortran/tests tests -q -m "not slow" 2>&1 | tail -3
git add -A docs README.rst CONTRIBUTING.rst climt_release_checks.txt HISTORY.rst scripts
git commit -m "docs: installation, contributing and release notes for the package split"
```

Expected: flake8 silent; no failures.
