# climt-fortran 2.0: cffi instead of Cython — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** climt-fortran 2.0.0 reaches RRTMG and Emanuel through cffi, with no Cython anywhere, the same numbers and no slower calls; climt moves to interface version 2.

**Architecture:** For each remaining library, a cffi out-of-line build declares the `bind(c)` routines (prototypes copied from the `.pyx` `cdef extern` blocks) and compiles a module `climt_fortran._<name>_cffi` linked against the static library `make` already builds. A plain-Python module `climt_fortran._<name>` keeps the functions climt calls today, doing what the `.pyx` did around each call (contiguity, level order, scratch arrays) and calling `lib.<routine>`. Reference outputs from the Cython build are recorded first and every new module must reproduce them exactly.

**Tech Stack:** cffi 2.x (API mode, setuptools `cffi_modules`), gfortran, NumPy, numba 0.68 (one compatibility test), pytest, cibuildwheel 3.3.1.

**Spec:** `docs/superpowers/specs/2026-11-19-cffi-fortran-design.md`. Runs after `docs/superpowers/plans/2026-11-19-numba-ports.md` is complete on this branch.

## Global Constraints

- No Cython at build or run time once Task 5 is done: `fortran/pyproject.toml` build requirements are `setuptools`, `wheel`, `cffi>=1.17`; no `.pyx` file remains.
- Same numbers as 1.0: each new module reproduces the Task 1 Cython reference outputs at `rtol = 0, atol = 0`.
- No slower calls: per call, at one column, each new module is at most **1.10×** its Task 1 Cython time (10% allows timing noise).
- No change to Fortran sources, except adding `bind(c)` to an entry point that lacks one.
- Fortran initialisation state (`rrtmg_set_constants`, `rrtmg_lw_ini_wrapper`, the shortwave equivalent, `init_emanuel_convection_fortran`) stays set once per process, as now.
- climt-fortran 2.0.0: `INTERFACE_VERSION = 2`; `MODULES = ("_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw")`; `cffi>=1.17` and `numpy>=1.16.0` in `install_requires`.
- climt: `EXPECTED_INTERFACE_VERSION = 2`; requirement `climt-fortran>=2,<3; sys_platform != "emscripten"`.
- **How to run.** `$PY` is shorthand: write it out as `/Users/joymonteiro/miniconda3/envs/climt-scm/bin/python` (the env from the numba-ports plan), from this worktree. After any change under `fortran/`, rebuild with `SDKROOT="$(xcrun --show-sdk-path)" $PY -m pip install --no-build-isolation --no-deps -e fortran` and check the output ends with `Successfully installed climt-fortran`. climt-fortran tests: `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q`. climt suite: `$PY -m pytest -m "not slow" -q`.
- flake8 (`tox.ini` settings) passes on `climt fortran/climt_fortran`.
- climt-fortran is not released and nothing is pushed: the release is the maintainer's step, after which the climt side (Task 6) can merge.

## Review Focus

1. **A stub install of 2.0** (`CLIMT_FORTRAN_STUB=1`: the Python wrappers present, no compiled `_cffi` modules): climt reports state `stub` with its existing message, not `broken`. Tests in Tasks 2 and 5.
2. **Input arrays that are not C-contiguous or not float64** (a transposed view; a float32 array): same numbers as contiguous float64, never a silent read of the wrong memory. Test in Task 2 for Emanuel; the shared helper covers RRTMG.
3. **Output arrays climt passes in to be filled** (the `cpdef void` functions write into caller arrays): they are filled in place as before, including when the caller's array is not contiguous. Tests in Tasks 3 and 4.
4. **climt 1.x with climt-fortran 2.0 installed** (a user upgrades only one): the existing `interface_mismatch` message, not an import error. Test in Task 6.
5. **The MCICA variants** of RRTMG longwave and shortwave (separate entry points, random cloud sampling with a fixed seed): exact parity like the plain ones. Reference cases in Task 1, checked in Tasks 3 and 4.

## File Structure

| Path | Responsibility |
|---|---|
| `fortran/_build_config.py` (new) | Compiler, flags, library paths and link arguments, moved out of `setup.py`; used by `setup.py` and the cffi builds |
| `fortran/climt_fortran/_cffi_build.py` (new) | One `cffi.FFI` builder per library: `emanuel_ffi`, `rrtmg_lw_ffi`, `rrtmg_sw_ffi` |
| `fortran/climt_fortran/_arrays.py` (new) | `as_double_pointer`, `as_int_pointer`: checked conversion of arrays and scalars for cffi calls |
| `fortran/climt_fortran/_emanuel_convection.py`, `_rrtmg_lw.py`, `_rrtmg_sw.py` (new, replacing the `.pyx`) | The functions climt calls, over `lib` |
| `fortran/setup.py`, `fortran/pyproject.toml`, `fortran/MANIFEST.in`, `fortran/README.rst`, `fortran/climt_fortran/__init__.py`, `fortran/_lib/Makefile` | Build, metadata, version 2.0.0 |
| `fortran/tests/generate_cython_reference.py` (new), `fortran/tests/fixtures/cython_reference.npz` (new) | Reference outputs and timings from the 1.0 build |
| `fortran/tests/benchmark_calls.py` (new), `fortran/tests/fixtures/cython_timings.json` (new) | Per-call timing, 1.0 and 2.0 |
| `fortran/tests/test_parity.py` (new), `fortran/tests/test_climt_fortran.py` | Tests |
| `climt/_core/fortran.py`, `pyproject.toml` (line 31 pins climt-fortran), `tests/test_fortran_loader.py`, `HISTORY.rst` | climt side |

---

### Task 1: Cython reference outputs, timings and a shared build config

**Files:**
- Create: `fortran/_build_config.py`, `fortran/tests/generate_cython_reference.py`, `fortran/tests/fixtures/cython_reference.npz`, `fortran/tests/benchmark_calls.py`, `fortran/tests/fixtures/cython_timings.json`
- Modify: `fortran/setup.py`

**Interfaces:**
- Produces:
  - `fortran/_build_config.py`: `def build_config() -> dict` with keys `libraries`, `library_dirs`, `extra_link_args`, `extra_compile_args`, `include_dirs` (no NumPy include), `lib_dir` (the `_lib/<system>/` path), and `def build_static_libraries() -> None` (runs `make` in `_lib`, raising `RuntimeError` on failure). Environment set-up (FC, CC, FFLAGS, CFLAGS, LDSHARED, Homebrew gcc discovery) moves here unchanged.
  - `cython_reference.npz`: keys `"<call>/in/<arg>"` and `"<call>/out/<arg>"` for calls `emanuel`, `rrtmg_lw`, `rrtmg_lw_mcica`, `rrtmg_sw`, `rrtmg_sw_mcica`, each on 4 columns × 30 levels with fixed inputs (`numpy.random.default_rng(20261119)`), clouds in levels 10–14 for the MCICA calls, and the arguments climt's components pass by default.
  - `benchmark_calls.py`: `def time_calls(n=2000, repeats=5) -> dict` mapping each call name to its best per-call time in microseconds at one column; run as a script it prints them, and with `--save <path>` writes them as JSON.

- [ ] **Step 1: Move the platform logic from `setup.py` into `_build_config.py`; `setup.py` imports it and still builds the five Cython extensions**

- [ ] **Step 2: Rebuild and run everything**

Run the rebuild command, then `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q` and `$PY -m pytest -m "not slow" -q 2>&1 | tail -2`
Expected: both pass with the same counts as before the change.

- [ ] **Step 3: Write and run the reference generator and the timer against the Cython build**

Run: `$PY fortran/tests/generate_cython_reference.py && NUMBA_NUM_THREADS=1 $PY fortran/tests/benchmark_calls.py --save fortran/tests/fixtures/cython_timings.json`
Expected: both files written; five calls in each.

- [ ] **Step 4: Commit**

```bash
git add fortran/_build_config.py fortran/setup.py fortran/tests
git commit -m "build(fortran): shared build config; Cython reference outputs and timings"
```

---

### Task 2: cffi infrastructure and Emanuel

**Files:**
- Create: `fortran/climt_fortran/_cffi_build.py`, `fortran/climt_fortran/_arrays.py`, `fortran/climt_fortran/_emanuel_convection.py`, `fortran/tests/test_parity.py`
- Delete: `fortran/climt_fortran/_emanuel_convection.pyx`
- Modify: `fortran/setup.py`, `fortran/pyproject.toml` (add `cffi>=1.17` to build requirements; Cython stays until Task 5), `fortran/climt_fortran/__init__.py`, `climt/_core/fortran.py`

**Interfaces:**
- Consumes: `build_config`, `build_static_libraries` (Task 1); the reference outputs.
- Produces:

```python
# fortran/climt_fortran/_cffi_build.py
emanuel_ffi: cffi.FFI      # set_source("climt_fortran._emanuel_convection_cffi", ...), linked to libemanuel.a
rrtmg_lw_ffi: cffi.FFI     # added in Task 3
rrtmg_sw_ffi: cffi.FFI     # added in Task 4

# fortran/climt_fortran/_arrays.py
def as_double_pointer(ffi, array):
    """cffi double* to a C-contiguous float64 array. Raises TypeError naming the
    dtype if not float64, ValueError if not C-contiguous (callers make arrays
    contiguous first, as the .pyx did)."""
def as_int_pointer(ffi, value) -> "int32_t *":
    """A new int32_t* holding value."""

# fortran/climt_fortran/_emanuel_convection.py
# Same functions, arguments, defaults, return values and in-place behaviour as the
# .pyx: init_emanuel_convection(...), convect(...).
```

- `setup.py` passes `cffi_modules=["climt_fortran/_cffi_build.py:emanuel_ffi"]` (more added in Tasks 3–4); its `build_ext` subclass still runs `build_static_libraries()` first.
- Presence test for interface 2: a library counts as installed when `climt_fortran._<name>_cffi` is found. `climt_fortran.available()` checks `_<name>_cffi`; climt's loader (`climt/_core/fortran.py`) checks `_<name>_cffi` with `PathFinder.find_spec` exactly as it checks `_<name>` today, then imports `_<name>`. Until Task 5 some modules are still Cython extensions, so a library also counts as installed when `_<name>` itself is found as an extension module (`spec.origin` ends with one of `importlib.machinery.EXTENSION_SUFFIXES`); a plain `_<name>.py` alone never counts, so a 2.0 stub reports `stub`.

- [ ] **Step 1: Write the failing tests** (`fortran/tests/test_parity.py`)

```python
@pytest.mark.parametrize("call", ["emanuel"])          # extended in Tasks 3 and 4
def test_matches_cython_reference(call):
    outputs = run_call(call, reference_inputs(call))   # calls the new module as climt does
    for name, expected in reference_outputs(call).items():
        np.testing.assert_array_equal(outputs[name], expected, err_msg=name)

def test_emanuel_rejects_float32():
    # convect with float32 temperature raises TypeError mentioning "float64"

def test_emanuel_transposed_input_matches():
    # the Task 1 inputs passed as np.asfortranarray copies give the same outputs

def test_emanuel_callable_from_numba():
    numba = pytest.importorskip("numba")
    # register climt_fortran._emanuel_convection_cffi with
    # numba.core.typing.cffi_utils.register_module and call
    # lib.emanuel_convection from an @numba.njit function on the reference inputs;
    # outputs equal the reference (keeps approach A of the build step open)

def test_emanuel_not_slower_than_cython():
    # time_calls()["emanuel"] <= 1.10 * cython_timings["emanuel"]
```

- [ ] **Step 2: Run them to see them fail**

Run: `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests/test_parity.py -q`
Expected: FAIL: `climt_fortran._emanuel_convection_cffi` does not exist.

- [ ] **Step 3: Implement `_cffi_build.py` (Emanuel), `_arrays.py`, `_emanuel_convection.py`; wire `cffi_modules`; delete the `.pyx`; update `available()` and climt's loader for the `_cffi` presence test**

Prototypes come from the `cdef extern` blocks of the deleted `.pyx`, with `cnp.int32_t` as `int32_t` and `double` unchanged. If the speed test fails because Python-level argument conversion dominates, move the call into an `@numba.njit` function using the registration from `test_emanuel_callable_from_numba`; numba then becomes a run-time dependency of climt-fortran, so add it to `install_requires`.

- [ ] **Step 4: Rebuild and run**

Run the rebuild command, then `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q` and `$PY -m pytest -m "not slow" -q -k "Emanuel or emanuel or fortran"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add fortran climt/_core/fortran.py
git commit -m "feat(fortran): Emanuel convection through cffi"
```

---

### Task 3: RRTMG longwave through cffi

**Files:**
- Create: `fortran/climt_fortran/_rrtmg_lw.py`
- Delete: `fortran/climt_fortran/_rrtmg_lw.pyx`
- Modify: `fortran/climt_fortran/_cffi_build.py` (`rrtmg_lw_ffi`), `fortran/setup.py`, `fortran/tests/test_parity.py`

**Interfaces:**
- Consumes: `as_double_pointer`, `as_int_pointer` (Task 2); `rrtmg_lw_ffi` linked to `librrtmg_lw.a` with OpenMP link flags as the Cython extension had.
- Produces: `climt_fortran._rrtmg_lw` with the same functions, arguments, defaults and in-place behaviour as the `.pyx`: `initialise_rrtm_radiation_mcica`, `initialise_rrtm_radiation`, `set_constants`, `rrtm_calculate_longwave_fluxes`, `rrtm_calculate_longwave_fluxes_mcica`.

- [ ] **Step 1: Add the failing tests**: extend the parity parametrisation with `"rrtmg_lw"` and `"rrtmg_lw_mcica"`; add `test_rrtmg_lw_fills_non_contiguous_output` (an output array passed as a strided view is filled with the reference values) and `test_rrtmg_lw_not_slower_than_cython` (both calls, 1.10×).

- [ ] **Step 2: Run them to see them fail**

Run: `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests/test_parity.py -q -k rrtmg_lw`
Expected: FAIL: `climt_fortran._rrtmg_lw_cffi` does not exist.

- [ ] **Step 3: Implement; delete the `.pyx`; add the builder to `cffi_modules`**

- [ ] **Step 4: Rebuild and run**

Run the rebuild command, then `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q` and `$PY -m pytest -m "not slow" -q -k "RRTMG or rrtmg or Longwave or fortran"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add fortran
git commit -m "feat(fortran): RRTMG longwave through cffi"
```

---

### Task 4: RRTMG shortwave through cffi

Same as Task 3 for the shortwave library.

**Files:**
- Create: `fortran/climt_fortran/_rrtmg_sw.py`
- Delete: `fortran/climt_fortran/_rrtmg_sw.pyx`
- Modify: `fortran/climt_fortran/_cffi_build.py` (`rrtmg_sw_ffi`), `fortran/setup.py`, `fortran/tests/test_parity.py`

**Interfaces:**
- Produces: `climt_fortran._rrtmg_sw` with `initialise_rrtm_radiation`, `initialise_rrtm_radiation_mcica`, `set_constants`, `rrtm_calculate_shortwave_fluxes`, `rrtm_calculate_shortwave_fluxes_mcica`, as the `.pyx`.

- [ ] **Step 1: Add the failing tests**: parity for `"rrtmg_sw"` and `"rrtmg_sw_mcica"`; `test_rrtmg_sw_fills_non_contiguous_output`; `test_rrtmg_sw_not_slower_than_cython`.

- [ ] **Step 2: Run them to see them fail**

Run: `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests/test_parity.py -q -k rrtmg_sw`
Expected: FAIL: `climt_fortran._rrtmg_sw_cffi` does not exist.

- [ ] **Step 3: Implement; delete the `.pyx`; add the builder to `cffi_modules`**

- [ ] **Step 4: Rebuild and run**

Run the rebuild command, then `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q` and `$PY -m pytest -m "not slow" -q -k "RRTMG or rrtmg or Shortwave or fortran"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add fortran
git commit -m "feat(fortran): RRTMG shortwave through cffi"
```

---

### Task 5: climt-fortran 2.0.0 without Cython

**Files:**
- Delete: `fortran/climt_fortran/_simple_physics.pyx`, `fortran/climt_fortran/_dcmip.pyx`, `fortran/_lib/simple_physics/`, `fortran/_lib/dcmip/`
- Modify: `fortran/_lib/Makefile` (drop their targets), `fortran/setup.py` (no Cython import, no `Extension` list, version `2.0.0`, `install_requires=["numpy>=1.16.0", "cffi>=1.17"]`), `fortran/pyproject.toml`, `fortran/MANIFEST.in` (no `*.pyx`), `fortran/climt_fortran/__init__.py` (`__version__ = "2.0.0"`, `INTERFACE_VERSION = 2`, `MODULES`), `fortran/README.rst`, `fortran/tests/test_climt_fortran.py` (`EXPECTED`, `ENTRY_POINTS`, interface 2)

**Interfaces:**
- Produces: climt-fortran 2.0.0 as stated in Global Constraints. The stub build (`CLIMT_FORTRAN_STUB=1`) still needs no compiler and ships the Python wrappers without `_cffi` modules.

- [ ] **Step 1: Update `fortran/tests/test_climt_fortran.py` to the 2.0 expectations, and add `test_no_cython_in_build_requirements` (reads `fortran/pyproject.toml`; `Cython` absent, `cffi` present) and `test_stub_has_wrappers_but_no_cffi_modules`; run them to see them fail**

Run: `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests/test_climt_fortran.py -q`
Expected: FAIL on interface version, `MODULES` and build requirements.

- [ ] **Step 2: Make the changes listed under Files**

- [ ] **Step 3: Rebuild from clean and run**

Run: `rm -rf fortran/build && find fortran -name "*.so" -delete`, then the rebuild command, then `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q`
Expected: PASS. climt itself now reports `interface_mismatch` until Task 6; that is expected here.

- [ ] **Step 4: Build both wheels**

Run: `$PY -m pip wheel --no-deps -w /private/tmp/claude-501/-Users-joymonteiro-github/a468073c-68b3-4f08-9d3c-dd919105d850/scratchpad/wheels ./fortran` and the same with `CLIMT_FORTRAN_STUB=1`
Expected: a platform wheel containing three `_*_cffi` binaries and no `.pyx`, and a `py3-none-any` stub wheel with no binaries.

- [ ] **Step 5: Commit**

```bash
git add -A fortran
git commit -m "feat(fortran)!: climt-fortran 2.0.0, cffi only, without simple physics and dcmip"
```

---

### Task 6: climt on interface version 2

**Files:**
- Modify: `climt/_core/fortran.py` (`EXPECTED_INTERFACE_VERSION = 2`), `pyproject.toml` (the climt-fortran pin, line 31), `tests/test_fortran_loader.py`, `HISTORY.rst`, `docs/` pages that describe climt-fortran 1.0 or Cython (`grep -rln "Cython\|climt-fortran" docs`)

**Interfaces:**
- Produces: climt requiring `climt-fortran>=2,<3; sys_platform != "emscripten"` and reporting `available` with the Task 5 build.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_fortran_loader.py`)

```python
def test_expects_interface_two():
    assert EXPECTED_INTERFACE_VERSION == 2

def test_interface_one_package_reports_mismatch(monkeypatch, tmp_path):
    # a fake climt_fortran package on sys.path with INTERFACE_VERSION = 1:
    # fortran_status()["state"] == "interface_mismatch" and unavailable_error's
    # message contains "pip install -U climt climt-fortran"

def test_stub_two_reports_stub(monkeypatch, tmp_path):
    # a fake package with INTERFACE_VERSION = 2 and the three wrapper .py files
    # but no _cffi modules: state == "stub"
```

- [ ] **Step 2: Run them to see them fail**

Run: `$PY -m pytest tests/test_fortran_loader.py -q`
Expected: FAIL on `test_expects_interface_two`.

- [ ] **Step 3: Make the changes; HISTORY under Unreleased: requires climt-fortran 2.0 (cffi, no Cython; simple physics and dcmip removed from it, now numba in climt)**

- [ ] **Step 4: Run the full suite, the climt-fortran tests and flake8**

Run: `$PY -m pytest -m "not slow" -q 2>&1 | tail -2`, `CLIMT_REQUIRE_FORTRAN=1 $PY -m pytest fortran/tests -q`, `$PY -m flake8 climt fortran/climt_fortran --config=tox.ini`
Expected: all pass; flake8 clean.

- [ ] **Step 5: Commit**

```bash
git add climt tests HISTORY.rst docs pyproject.toml
git commit -m "feat: climt requires climt-fortran 2 (interface version 2)"
```
