# JAX Merge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the Emanuel and CORK longwave JAX kernels from `jax-research` onto `develop` behind parity tests and a compile-time budget. (The build step that was planned here moved to `feature/scm-build` on 2026-11-19.)

**Architecture:** A real merge of `origin/jax-research` (minus `JaxBackend` and research notes), then clean-up: JAX code is isolated in lazily imported modules, and the Emanuel kernel's Python-unrolled level loops are regrouped to shrink the traced graph. numba stays the reference.

**Tech Stack:** JAX 0.10.1 (CPU, 64-bit), numba, pytest.

**Spec:** `docs/superpowers/specs/2026-11-14-jax-merge-build-step-design.md`

## Global Constraints

- `import climt` must not import `jax`, with or without JAX installed. No module that `import climt` loads may contain a module-level `import jax`.
- climt never calls `jax.config.update`. (Test modules may.)
- Not merged: `climt/_core/jax_backend.py`, `JAX_KERNEL_FUSION.md`, `PROPOSAL_JAX_BACKEND.md`, `docs/HYBRID_STRATEGY.md`, `docs/SENSITIVITY_ANALYSIS_STRATEGY.md`, `docs/predict_parameters.py`, `docs/predict_parameters_extremes.py`.
- Parity tolerances are fixed and may not be loosened: Emanuel JAX vs numba `rtol = atol = 1e-10`; CORK longwave JAX vs numba `rtol = 1e-6` (`atol = 1e-9` on tendencies, `1e-6` on fluxes).
- Budget, Emanuel kernel, 30 levels, 1 column, CPU, JAX 0.10.1, maintainer's machine: gradient first-call time at most **120 s**; forward first-call time at most **30 s**.
- The numba and NumPy code paths of every component must produce bit-identical results before and after this plan (the existing suite is the check). No numba kernel is edited.
- JAX extra: `jax>=0.10.1`, `jaxlib>=0.10.1`. Tests needing JAX start with `jax = pytest.importorskip("jax")` and carry `pytestmark = pytest.mark.jax`.
- flake8 with the repo's `tox.ini` settings (max line length 150) passes on `climt`.
- Add to the `Unreleased` section of `HISTORY.rst`; do not bump the version.
- **How to run.** Use the interpreter of the `climt-dev` environment with this worktree first on the path. In the commands below `$PY` is shorthand, not a real variable: write it out each time as `env PYTHONPATH="$PWD" /Users/joymonteiro/miniconda3/envs/climt-dev/bin/python` (shell variables do not survive between tool calls, and zsh does not split an unquoted variable into words). Do not `pip install` anything into `climt-dev`. This worktree has no compiled Fortran extensions, so tests of the Fortran-backed components fail with `ImportError` here; record that baseline in Task 1 and compare against it, do not try to fix it.
- **Running without JAX.** `tests/no_jax_shim/jax/__init__.py` (created in Task 2) raises `ImportError`. Prefixing `PYTHONPATH` with `$PWD/tests/no_jax_shim` gives a run in which JAX cannot be imported.
- Nothing is pushed.

## Review Focus

1. **JAX arrays reaching code with no JAX kernel** (a component without one, or CORK longwave in a configuration its JAX path does not cover): an error naming the component and the reason, never a quiet NumPy computation. Tests in Task 2.
2. **64-bit mode off** (`jax_enable_x64` false, arrays are `float32`): `JaxPrecisionError` with the line that fixes it, not a parity failure at the tenth digit. Test in Task 2.
3. **Convection switching on or off inside a differentiated run**: gradients stay finite on stable (non-convecting) as well as unstable profiles. Test in Task 4.

## File Structure

| Path | Responsibility |
|---|---|
| `climt/_core/jax_support.py` (new) | `is_jax_array`, `require_x64`, `JaxPrecisionError`, `count_jaxpr_equations` |
| `climt/_components/emanuel/kernel_jax.py` (new) | The Emanuel JAX kernel, moved out of `pure_python_v3.py` |
| `climt/_components/emanuel/pure_python_v3.py` | numba kernel and component; lazy dispatch to `kernel_jax` |
| `climt/_components/cork/lw/kernels_jax.py`, `climt/_components/cork/optics/correlated_k_jax.py` (from the merge) | CORK longwave JAX kernels |
| `climt/_components/cork/lw/component.py` | Lazy JAX table, explicit errors for uncovered configurations |
| `benchmarks/benchmark_jax_compile.py` (new), `benchmarks/results/emanuel_jax_compile.txt` (new) | Compile-time and graph-size measurements |
| `tests/jax_cases.py` (new) | Shared builders for JAX test inputs |
| `tests/no_jax_shim/jax/__init__.py` (new) | Makes `import jax` fail on demand |
| `tests/test_jax_isolation.py`, `tests/test_jax_dispatch.py`, `tests/test_emanuel_jax_budget.py` (new) | Tests |
| `tests/test_emanuel_jax_parity.py`, `tests/test_jax_differentiation.py`, `tests/test_cork_lw_jax_integration.py`, `tests/test_cork_lw_kernels_jax_parity.py`, `tests/test_cork_optics_jax_parity.py` (from the merge) | Parity and differentiation tests |
| `setup.py` or `pyproject.toml` (whichever declares extras), `tox.ini`, `.github/workflows/python-package.yml`, `HISTORY.rst` | Extra, marker, CI job, changelog |

---

### Task 1: Merge `jax-research`

**Files:**
- Merge: `origin/jax-research` into this branch
- Remove in the merge commit: the seven paths listed in Global Constraints
- Modify: `climt/_components/emanuel/pure_python_v3.py` (conflict), `tests/test_jax_differentiation.py`, `tests/test_emanuel_jax_parity.py`, `tox.ini`

**Interfaces:**
- Consumes: nothing.
- Produces: on this branch, `EmanuelConvectionPython.array_call` and `CorkLongwaveRadiation.array_call` dispatch to JAX kernels when given JAX arrays (still with module-level JAX imports; Task 2 removes those). `climt._components.cork.lw.kernels_jax.cork_lw_jax`, `build_cork_lw_table`, `CorkLWTable`; in `pure_python_v3.py`: `_jax_vectorized_convect`, `_convect_functional_jax`, `_tlift_functional_jax`, `compute_qs_jax`.

- [ ] **Step 1: Record the baseline before merging**

```bash
$PY -m pytest tests -q -m "not slow" -p no:cacheprovider 2>&1 | tail -3
```

Write the summary line (passed / failed / skipped / errors) into your report as BASELINE, with the list of failing test files. Failures from missing Fortran extensions are expected.

- [ ] **Step 2: Merge without committing**

```bash
git merge --no-ff --no-commit origin/jax-research
```

Expected: one conflict, in `climt/_components/emanuel/pure_python_v3.py`.

- [ ] **Step 3: Resolve the conflict in `array_call`**

Keep the branch's structure (`_use_jax` deciding between the JAX call and the NumPy call). In it:

- `delt = float(timestep.total_seconds())`, with `develop`'s comment explaining why it is a plain float (numba strips units, Pyodide has no numba), placed above that line. Both paths use this `delt`.
- The NumPy path is otherwise exactly `develop`'s.
- The JAX path is exactly the branch's.

Everything else in the file is the branch's side (it only adds code). The alias `EmanuelConvectionPythonV3 = EmanuelConvectionPython` must appear once, at the end of the file.

- [ ] **Step 4: Drop what is not being merged, and commit the merge**

```bash
git rm -q -f climt/_core/jax_backend.py JAX_KERNEL_FUSION.md PROPOSAL_JAX_BACKEND.md \
    docs/HYBRID_STRATEGY.md docs/SENSITIVITY_ANALYSIS_STRATEGY.md \
    docs/predict_parameters.py docs/predict_parameters_extremes.py
git grep -n "jax_backend\|JaxBackend\|JaxStateContainer" -- climt tests    # expected: no output
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "Merge jax-research: Emanuel and CORK longwave JAX kernels (without JaxBackend)"
```

- [ ] **Step 5: Make the merged tests collectable**

- `tests/test_jax_differentiation.py`: it imports `EmanuelConvectionPythonV3` from `climt`, which does not export it. Use `EmanuelConvectionPython`. Remove the `print` calls and the `if __name__ == "__main__"` block.
- `tests/test_emanuel_jax_parity.py`: use `EmanuelConvectionPython` from `climt`.
- All five merged JAX test modules: the first statements become `import pytest`, `jax = pytest.importorskip("jax")`, `pytestmark = pytest.mark.jax`, followed by the existing `jax.config.update(...)` lines.
- `tox.ini`, `[pytest] markers`: add `jax: needs JAX (pip install "climt[jax]"); skipped when JAX is absent.`

- [ ] **Step 6: Run**

```bash
$PY -m pytest tests/test_emanuel_jax_parity.py tests/test_jax_differentiation.py \
    tests/test_cork_lw_jax_integration.py tests/test_cork_lw_kernels_jax_parity.py \
    tests/test_cork_optics_jax_parity.py -q 2>&1 | tail -5
$PY -m pytest tests -q -m "not slow and not jax" -p no:cacheprovider 2>&1 | tail -3
flake8 climt --config=tox.ini
```

Expected: the five JAX modules pass. `test_jax_differentiation.py` compiles a gradient through the Emanuel kernel and takes about 14 minutes at this point; mark that one test `@pytest.mark.slow` now and run it once by itself to confirm it passes (`-m slow`). The second command matches BASELINE exactly. Fix any flake8 findings in merged files by reformatting only (no behaviour change); if a merged file has many, fix them in a separate commit.

- [ ] **Step 7: Commit**

```bash
git add tests tox.ini climt
git commit -m "test: make the merged JAX tests collectable and mark them"
```

---

### Task 2: Isolate JAX

**Files:**
- Create: `climt/_core/jax_support.py`, `climt/_components/emanuel/kernel_jax.py`, `tests/no_jax_shim/jax/__init__.py`, `tests/test_jax_isolation.py`, `tests/test_jax_dispatch.py`
- Modify: `climt/_components/emanuel/pure_python_v3.py`, `climt/_components/cork/lw/component.py`, `climt/_components/cork/lw/kernels_jax.py` and `climt/_components/cork/optics/correlated_k_jax.py` (remove their `jax.config.update` lines), `setup.py` or `pyproject.toml`, `tox.ini`, `.github/workflows/python-package.yml`

**Interfaces:**
- Consumes: Task 1's kernels.
- Produces:

```python
# climt/_core/jax_support.py — imports nothing from jax at module level
class JaxPrecisionError(RuntimeError): ...

def is_jax_array(value) -> bool:
    """True for a jax.Array (including tracers). Never imports jax: if 'jax'
    is not in sys.modules the answer is False."""

def require_x64() -> None:
    """Raise JaxPrecisionError unless jax.config.jax_enable_x64 is True. The
    message contains: jax.config.update("jax_enable_x64", True)"""

# climt/_components/emanuel/kernel_jax.py — imports jax at module level; imported lazily
def emanuel_convect_jax(t, q, qs, u, v, p, ph, nlev, nl, ntra, delt, cbmf, params, cond): ...   # was _jax_vectorized_convect
def compute_qs_jax(T, P, cond, RD): ...
```

  Class attribute `has_jax_kernel = True` on `EmanuelConvectionPython` and `CorkLongwaveRadiation` (absent on every other component; readers use `getattr(component, "has_jax_kernel", False)`).

Required behaviour:
- `pure_python_v3.py` has no JAX import. `array_call` does `if is_jax_array(t): require_x64(); from .kernel_jax import ...`. The JAX section of the file (from the comment `# JAX path:` to just before the `EmanuelConvectionPythonV3` alias) moves verbatim to `kernel_jax.py`, which imports what it needs from `pure_python_v3`'s siblings (`EmanuelParams` stays in `pure_python_v3.py`; `kernel_jax.py` must not import `pure_python_v3` at module level if that creates a cycle — take `params` and `cond` as arguments, as the functions already do).
- `cork/lw/component.py` has no module-level JAX import and no `HAS_JAX`. `__init__` sets `self._jax_table = None` and does nothing else for JAX. `array_call` begins:
  - if `is_jax_array(state["T"])`: call `require_x64()`; then raise `NotImplementedError` naming `CorkLongwaveRadiation` and the reason if `optics != "correlated_k"`, or the table's `k_coefficients` is not 7-dimensional or the table is not pre-mixed background (`self._premixed_bg`), or `self._diagnostics_level != 0`; otherwise build the table once (`build_cork_lw_table`) and return `self._array_call_jax(state)`.
  - The blanket `try/except Exception` around the table build is removed.
- The JAX kernel modules (`kernels_jax.py`, `correlated_k_jax.py`, `kernel_jax.py`) contain no `jax.config.update` call. The two as merged each have one directly under `import jax`; delete it. 64-bit mode is checked by `require_x64()` at the two component entry points instead.

- [ ] **Step 1: Write the failing tests**

`tests/no_jax_shim/jax/__init__.py`:

```python
raise ImportError("jax is blocked for this test run (tests/no_jax_shim)")
```

`tests/test_jax_isolation.py`:

```python
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SHIM = str(Path(__file__).resolve().parent / "no_jax_shim")


def run(code, block_jax=False):
    env = dict(os.environ)
    paths = [str(REPO)] + ([SHIM] if block_jax else [])
    env["PYTHONPATH"] = os.pathsep.join(paths[::-1] + [env.get("PYTHONPATH", "")])
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)


def test_import_climt_does_not_import_jax():
    result = run("import sys, climt; assert 'jax' not in sys.modules, 'jax was imported'")
    assert result.returncode == 0, result.stderr


def test_climt_works_when_jax_cannot_be_imported():
    code = (
        "import climt, numpy as np\n"
        "from climt._core.jax_support import is_jax_array\n"
        "assert is_jax_array(np.zeros(3)) is False\n"
        "climt.EmanuelConvectionPython(); climt.CorkLongwaveRadiation()\n"
    )
    result = run(code, block_jax=True)
    assert result.returncode == 0, result.stderr


def test_no_module_level_jax_import_in_loaded_modules():
    code = (
        "import sys, climt, re, pathlib\n"
        "bad = []\n"
        "for name, module in list(sys.modules.items()):\n"
        "    path = getattr(module, '__file__', None) or ''\n"
        "    if name.startswith('climt') and path.endswith('.py'):\n"
        "        if re.search(r'^(import jax|from jax)', pathlib.Path(path).read_text(), re.M):\n"
        "            bad.append(name)\n"
        "assert not bad, bad\n"
    )
    result = run(code)
    assert result.returncode == 0, result.stderr


def test_climt_never_changes_jax_config():
    offenders = [str(p) for p in (REPO / "climt").rglob("*.py") if "config.update" in p.read_text()]
    assert offenders == []
```

`tests/test_jax_dispatch.py` (a separate module, so that a run without JAX skips these and still runs the isolation tests above):

```python
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax

import jax.numpy as jnp       # noqa: E402
import numpy as np            # noqa: E402

import climt                  # noqa: E402
from climt._core.jax_support import is_jax_array   # noqa: E402
from tests.test_jax_isolation import run            # noqa: E402


def test_is_jax_array():
    assert is_jax_array(jnp.zeros(3)) is True
    assert is_jax_array(np.zeros(3)) is False
    assert jax.jit(lambda x: jnp.asarray(float(is_jax_array(x))))(jnp.zeros(3)) == 1.0   # tracers count


def test_x64_off_is_refused_with_the_fix():
    code = (
        "import jax, jax.numpy as jnp, climt\n"
        "from datetime import timedelta\n"
        "from climt._core.jax_support import JaxPrecisionError\n"
        "assert not jax.config.jax_enable_x64\n"
        "state = {k: jnp.ones((30, 1)) for k in ('air_temperature', 'specific_humidity', 'air_pressure',"
        " 'eastward_wind', 'northward_wind')}\n"
        "state['air_pressure_on_interface_levels'] = jnp.ones((31, 1))\n"
        "try:\n"
        "    climt.EmanuelConvectionPython().array_call(state, timedelta(minutes=15))\n"
        "except JaxPrecisionError as error:\n"
        "    assert 'jax.config.update(\"jax_enable_x64\", True)' in str(error)\n"
        "else:\n"
        "    raise SystemExit('no error')\n"
    )
    result = run(code)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("kwargs, word", [
    ({"optics": "correlated_k", "table": "earth_low_res_lw", "diagnostics_level": 1}, "diagnostics_level"),
    ({"optics": "parmentier"}, "optics"),
])
def test_cork_lw_refuses_uncovered_configurations(kwargs, word):
    jax.config.update("jax_enable_x64", True)
    component = climt.CorkLongwaveRadiation(**kwargs)
    with pytest.raises(NotImplementedError, match=word) as info:
        component.array_call({"T": jnp.ones((20, 2))})
    assert "CorkLongwaveRadiation" in str(info.value)


def test_has_jax_kernel_flags():
    flagged = {cls.__name__ for cls in climt._components.__all__ if getattr(cls, "has_jax_kernel", False)}
    assert flagged == {"EmanuelConvectionPython", "CorkLongwaveRadiation"}
```

(`"parmentier"` is `CorkLongwaveRadiation`'s default optics mode, which has no JAX kernel.)

- [ ] **Step 2: Run and confirm failure**

Run: `$PY -m pytest tests/test_jax_isolation.py tests/test_jax_dispatch.py -q`
Expected: failures (`jax was imported`; no module `climt._core.jax_support`).

- [ ] **Step 3: Implement** `jax_support.py`, the move to `kernel_jax.py`, and the CORK changes as specified above. `count_jaxpr_equations` is Task 3.

- [ ] **Step 4: Declare the extra and the CI job**

- Extras: add `jax = ["jax>=0.10.1", "jaxlib>=0.10.1"]` where this branch declares extras (`extras_require` in `setup.py`, creating it if absent; or `[project.optional-dependencies]` if `pyproject.toml` already has a `[project]` table).
- `tox.ini`: an environment `py312-jax` that installs as the default environment does plus `jax>=0.10.1 jaxlib>=0.10.1`, and runs `py.test --basetemp={envtmpdir} -v -m "jax and not slow"`.
- `.github/workflows/python-package.yml`: one additional job (or matrix `include` entry) on `ubuntu-latest`, Python 3.12, running `tox -e py312-jax`.

- [ ] **Step 5: Run everything, with and without JAX**

```bash
$PY -m pytest tests/test_jax_isolation.py tests/test_jax_dispatch.py tests/test_emanuel_jax_parity.py tests/test_cork_lw_jax_integration.py \
    tests/test_cork_lw_kernels_jax_parity.py tests/test_cork_optics_jax_parity.py -q 2>&1 | tail -3
env PYTHONPATH="$PWD/tests/no_jax_shim:$PWD" /Users/joymonteiro/miniconda3/envs/climt-dev/bin/python \
    -m pytest tests -q -m "not slow" -p no:cacheprovider 2>&1 | tail -3
flake8 climt --config=tox.ini
```

Expected: first command all pass. Second command: the same passed/failed counts as BASELINE plus the new non-JAX tests, with every `jax`-marked test skipped. flake8 silent.

- [ ] **Step 6: Commit**

```bash
git add climt tests tox.ini .github/workflows/python-package.yml setup.py pyproject.toml
git commit -m "refactor: import JAX lazily; explicit errors for precision and uncovered configurations"
```

---

### Task 3: Measure compile time and graph size

**Files:**
- Modify: `climt/_core/jax_support.py`
- Create: `tests/jax_cases.py`, `benchmarks/benchmark_jax_compile.py`, `benchmarks/results/emanuel_jax_compile.txt`
- Test: `tests/test_emanuel_jax_budget.py` (first part)

**Interfaces:**
- Consumes: `emanuel_convect_jax`, `EmanuelConvectionPython` (Task 2).
- Produces:

```python
# climt/_core/jax_support.py
def count_jaxpr_equations(function, *args) -> int:
    """Trace function(*args) with jax.make_jaxpr and return the number of
    equations, descending into every nested jaxpr (jit, scan, cond, while
    bodies) and counting each nested jaxpr object once."""

# tests/jax_cases.py
def emanuel_case(nlev=30, ncol=1, profile="unstable", as_jax=True):
    """Return (component, state, timestep): an EmanuelConvectionPython, the
    array_call state dictionary (pressures in hPa, arrays (nlev, ncol)), and
    timedelta(minutes=15). profile is 'moist', 'dry', 'unstable' (as in
    tests/test_emanuel_jax_parity.py) or 'stable' (isothermal 250 K,
    specific humidity 1e-6)."""

def emanuel_forward(component, state, timestep):
    """Function of air_temperature returning the temperature tendency array."""

def emanuel_loss(component, state, timestep):
    """Function of air_temperature returning sum(tendency**2); for jax.grad."""
```

`benchmarks/benchmark_jax_compile.py`, run as a script, prints and appends to `benchmarks/results/emanuel_jax_compile.txt` one block containing: date, JAX version, `platform.machine()`, `nlev`, `ncol`, and for the forward function and for `jax.grad` of the loss: equation count, trace time, first-call time, and the median of 20 steady calls. Options: `--nlev` (default 30), `--ncol` (default 1), `--skip-grad`, `--label TEXT`.

- [ ] **Step 1: Write the tests**

`tests/test_emanuel_jax_budget.py`:

```python
import pytest

jax = pytest.importorskip("jax")
pytestmark = pytest.mark.jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")

import jax.numpy as jnp                                  # noqa: E402
from climt._core.jax_support import count_jaxpr_equations   # noqa: E402


def test_count_is_exact_on_a_known_program():
    def body(carry, x):
        return carry + x * 2.0, carry

    def program(xs):
        return jax.lax.scan(body, 0.0, xs)[0] + 1.0

    # top level: scan, add  -> 2; scan body: mul, add -> 2
    assert count_jaxpr_equations(program, jnp.arange(5.0)) == 4


def test_count_descends_into_jitted_functions():
    inner = jax.jit(lambda x: jnp.sin(x) + jnp.cos(x))     # sin, cos, add -> 3
    assert count_jaxpr_equations(lambda x: inner(x) * 2.0, 1.0) == 5   # pjit, mul, + 3
```

If JAX's lowering of these two programs yields different counts on the installed version (for example an extra `convert_element_type`), keep the tests' structure, set the expected numbers to what `jax.make_jaxpr` prints, and paste that printout in a comment above each assertion.

- [ ] **Step 2: Run and confirm failure**

Run: `$PY -m pytest tests/test_emanuel_jax_budget.py -q`
Expected: `ImportError: cannot import name 'count_jaxpr_equations'`.

- [ ] **Step 3: Implement `count_jaxpr_equations`, `tests/jax_cases.py` and the benchmark script**

- [ ] **Step 4: Record the baseline**

```bash
$PY -m pytest tests/test_emanuel_jax_budget.py -q
$PY benchmarks/benchmark_jax_compile.py --label "baseline (as merged)"
```

The gradient measurement takes about 15 minutes. Expected: the tests pass; `benchmarks/results/emanuel_jax_compile.txt` now holds the baseline block. Put the four headline numbers (forward and gradient equation counts, forward and gradient first-call seconds) in your report.

- [ ] **Step 5: Commit**

```bash
git add climt/_core/jax_support.py tests/jax_cases.py tests/test_emanuel_jax_budget.py benchmarks
git commit -m "bench: measure JAX compile time and traced graph size for the Emanuel kernel"
```

---

### Task 4: Bring the Emanuel kernel inside the budget

**Files:**
- Modify: `climt/_components/emanuel/kernel_jax.py`, `tests/test_emanuel_jax_budget.py`, `tests/test_emanuel_jax_parity.py`, `tests/test_jax_differentiation.py`, `benchmarks/results/emanuel_jax_compile.txt`

**Interfaces:**
- Consumes: `count_jaxpr_equations`, `tests/jax_cases.py`, the benchmark (Task 3).
- Produces: the same `emanuel_convect_jax` signature and results (to parity tolerance), with a smaller traced graph. Constants `FORWARD_EQUATION_LIMIT`, `GRADIENT_EQUATION_LIMIT` in `tests/test_emanuel_jax_budget.py`.

The kernel as merged has 14 Python `for ... in range(...)` loops over levels and 64 `.at[...]` updates inside the jitted function; one block (the downdraft) already uses `lax.scan`. Regroup the others by these rules, one block at a time, in file order:

| Loop shape | Replace with |
|---|---|
| Body at level `i` reads only level `i` of arrays plus loop-invariant values, and writes level `i` (`X = X.at[i].set(jnp.where(mask_i, f_i, X[i]))`) | One whole-array expression: `idx = jnp.arange(NL + 1)`, `X = jnp.where(mask(idx), f(idx), X)` |
| Body at level `i` reads a value written at level `i - 1` (or `i + 1`) in the same loop | `jnp.cumsum` / `jnp.cumprod` / a cumulative maximum when the recurrence is one of those; otherwise `jax.lax.scan` (with `reverse=True` for downward sweeps) |
| Two nested level loops (`i`, `j`) with a masked body | Broadcasting over an `(NL + 1, NL + 1)` grid with the same masks, or `jax.vmap` over `i` of a whole-array body in `j` |
| A fixed, small iteration count unrelated to levels (the two saturation-adjustment iterations) | Leave unrolled |

Rules that do not bend: results agree with the numba kernel to `rtol = atol = 1e-10`; `jnp.where` branches must stay finite on both sides wherever a gradient flows (guard divisions and logs with safe operands, as the merged code does), so gradients have no NaN; `static_argnums` and the `vmap` over columns are unchanged.

- [ ] **Step 1: Extend the tests first**

`tests/test_emanuel_jax_parity.py`: build inputs with `tests/jax_cases.emanuel_case` and add `"stable"` to the profiles and `nlev` in `(30, 45)`.

`tests/test_jax_differentiation.py`: parametrize over `profile in ("unstable", "moist", "stable")`; assert the gradient is finite everywhere for all three, and non-zero for `"unstable"`. Keep it `slow` until Step 4.

`tests/test_emanuel_jax_budget.py`, added:

```python
import os
import time

from tests.jax_cases import emanuel_case, emanuel_forward, emanuel_loss

FORWARD_EQUATION_LIMIT = None    # set in Step 4
GRADIENT_EQUATION_LIMIT = None   # set in Step 4
SCALE = float(os.environ.get("CLIMT_JAX_BUDGET_SCALE", "1.0"))


def test_forward_graph_size():
    component, state, timestep = emanuel_case()
    count = count_jaxpr_equations(emanuel_forward(component, state, timestep), state["air_temperature"])
    assert count <= FORWARD_EQUATION_LIMIT, count


def test_gradient_graph_size():
    component, state, timestep = emanuel_case()
    count = count_jaxpr_equations(jax.grad(emanuel_loss(component, state, timestep)), state["air_temperature"])
    assert count <= GRADIENT_EQUATION_LIMIT, count


@pytest.mark.slow
def test_forward_first_call_within_budget():
    component, state, timestep = emanuel_case()
    start = time.perf_counter()
    jax.block_until_ready(jax.jit(emanuel_forward(component, state, timestep))(state["air_temperature"]))
    assert time.perf_counter() - start <= 30.0 * SCALE


@pytest.mark.slow
def test_gradient_first_call_within_budget():
    component, state, timestep = emanuel_case()
    start = time.perf_counter()
    jax.block_until_ready(jax.jit(jax.grad(emanuel_loss(component, state, timestep)))(state["air_temperature"]))
    assert time.perf_counter() - start <= 120.0 * SCALE
```

- [ ] **Step 2: Regroup, one block at a time**

After each block: `$PY -m pytest tests/test_emanuel_jax_parity.py -q` must pass, then commit (`git commit -m "perf(emanuel-jax): vectorise <block name>"`). Every few blocks run `$PY benchmarks/benchmark_jax_compile.py --skip-grad --label "<block>"` to see the forward equation count fall. Start with the blocks that dominate: find them by printing, for the baseline, the equation count after temporarily returning early from the kernel at each block boundary, or simply by the number of `.at[` updates per loop.

- [ ] **Step 3: Measure the gradient**

```bash
$PY benchmarks/benchmark_jax_compile.py --label "after regrouping"
```

Expected: gradient first call at most 120 s and forward first call at most 30 s. If the gradient is still over, continue Step 2 on the blocks that remain as Python loops, then on the largest contributors to the gradient equation count. If, with no per-level Python loop left, the gradient first call is still above 120 s: stop, do not change the limit, and report DONE_WITH_CONCERNS with the benchmark blocks and the per-block equation counts.

- [ ] **Step 4: Pin the achieved graph sizes and un-slow the gradient test**

Set `FORWARD_EQUATION_LIMIT` and `GRADIENT_EQUATION_LIMIT` to the achieved counts times 1.10, rounded up to the next hundred. In the module docstring of `tests/test_emanuel_jax_budget.py` record: baseline counts and first-call times (Task 3), achieved counts and times, JAX version, machine. If the gradient first call is now under 60 s, remove `@pytest.mark.slow` from the `"unstable"` case of `tests/test_jax_differentiation.py` (keep it on the other two profiles); otherwise leave all three slow and say so. If either graph-size test itself takes more than 60 s, mark it `slow` and say so.

- [ ] **Step 5: Full check**

```bash
$PY -m pytest tests/test_emanuel_jax_parity.py tests/test_emanuel_jax_budget.py tests/test_jax_differentiation.py -q -m "jax" 2>&1 | tail -3
$PY -m pytest tests/test_emanuel_jax_budget.py tests/test_jax_differentiation.py -q -m "slow" 2>&1 | tail -3
$PY -m pytest tests -q -m "not slow and not jax" -p no:cacheprovider 2>&1 | tail -3
flake8 climt --config=tox.ini
```

Expected: all JAX tests pass, including the two wall-time tests and all three gradient profiles; the third command matches BASELINE (plus tests added since); flake8 silent.

- [ ] **Step 6: Commit**

```bash
git add climt/_components/emanuel/kernel_jax.py tests benchmarks/results/emanuel_jax_compile.txt
git commit -m "perf(emanuel-jax): meet the compile budget; pin traced graph size"
```

---

### Tasks 5–7: moved

On 2026-11-19 the build step (`compile_component`, `compile_tendencies`) and
the worked example moved to `feature/scm-build`, where they are designed for
numba first (`docs/superpowers/specs/2026-11-19-scm-build-design.md` and its
plan on that branch). The JAX kernels merged here plug into that build step
afterwards. This plan ends at Task 4.
