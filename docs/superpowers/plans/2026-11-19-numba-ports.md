# simple_physics and dcmip in numba — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Fortran behind `SimplePhysics` and `DcmipInitialConditions` with numba kernels that match it, so climt no longer needs climt-fortran for either.

**Architecture:** Reference outputs are generated from the Fortran first and committed as fixtures. Each scheme is then ported line by line into one numba kernel module beside its component; the component keeps its properties and constructor and calls the kernel. SimplePhysics parameters become a `NamedTuple` whose fields the kernel reads per column, so the build step (third project) can vary them by scenario.

**Tech Stack:** numba 0.68, NumPy, sympl 0.5.x, pytest. climt-fortran 1.0 (Cython) only to generate fixtures and for the speed comparison.

**Spec:** `docs/superpowers/specs/2026-11-19-numba-ports-design.md` (context: `2026-11-19-scm-build-design.md`).

## Global Constraints

- No change to the physics: parity with the Fortran at `rtol = 1e-10` (SimplePhysics) and `rtol = 1e-12` (dcmip), with `atol = 1e-12` in output units for values that are zero in the Fortran.
- Speed: at one column and one thread, the numba SimplePhysics kernel call is no slower than **1.5×** the Fortran call (`climt_fortran._simple_physics.get_new_state` with constants set once beforehand).
- Constants are read with `get_constant` once, in `__init__`, never in `array_call`.
- Inputs, outputs, `input_properties`/`diagnostic_properties`/`output_properties` and constructor arguments of both components are unchanged.
- Emanuel and RRTMG are not touched.
- **How to run.** `$PY` is shorthand: write it out as `/Users/joymonteiro/miniconda3/envs/climt-scm/bin/python`, run from this worktree (`/Users/joymonteiro/github/climt/.worktrees/scm-build`). The env is created in Task 1; never `pip install` into `climt-dev`. Unit tests: `$PY -m pytest <path> -q`. The full suite: `$PY -m pytest -m "not slow" -q`.
- flake8 with the repo's `tox.ini` settings passes on `climt`.
- Add to the `Unreleased` section of `HISTORY.rst`; do not bump the version.
- Nothing is pushed.

## Review Focus

1. **Supersaturated layers** (specific humidity above saturation at several levels, large-scale condensation on): condensation, latent heating and precipitation match the Fortran. Fixture case `lsc_supersaturated` in Task 1, checked in Task 2.
2. **Surface winds above 20 m/s**, where the momentum drag coefficient switches to its maximum: matches the Fortran. Fixture case `strong_wind` (Task 1), checked in Task 2.
3. **Surface temperature below 271 K** with `use_external_surface_specific_humidity=False`: the surface humidity branch matches. Fixture case `cold_surface` (Task 1), checked in Task 2.
4. **Two SimplePhysics instances with different parameters called alternately** (the Fortran re-set its global constants before every call for this reason): each gives its own answer every time. Test in Task 3.
5. **Input arrays that are not C-contiguous** (a state whose `air_temperature` is a transposed view): same numbers as a contiguous state. Test in Task 3.

## File Structure

| Path | Responsibility |
|---|---|
| `scripts/generate_fortran_reference.py` (new) | Runs the Fortran-backed components on fixed cases and writes the fixtures |
| `tests/fixtures/fortran_reference/simple_physics.npz`, `dcmip.npz`, `cases.json` (new) | Fortran outputs, their inputs and constructor arguments |
| `climt/_core/parameters.py` (new) | `per_column(params, ncol)`; `PARALLEL_MIN_COLUMNS` |
| `climt/_components/simple_physics/kernel.py` (new) | `SimplePhysicsParams`, `SimplePhysicsFlags`, `simple_physics_kernel` |
| `climt/_components/simple_physics/component.py` | Calls the kernel; no Fortran |
| `climt/_components/dcmip/kernels.py` (new) | `baroclinic_wave_ics`, `tropical_cyclone_ics` |
| `climt/_components/dcmip/component.py` | Calls the kernels; no Fortran |
| `climt/_core/fortran.py` | `MODULES` drops `_simple_physics` and `_dcmip` |
| `benchmarks/benchmark_simple_physics_port.py` (new), `benchmarks/results/simple_physics_port.txt` (new) | Speed comparison with the Fortran |
| `tests/test_simple_physics_port.py`, `tests/test_dcmip_port.py`, `tests/test_parameters.py` (new) | Tests |

---

### Task 1: Environment and Fortran reference fixtures

**Files:**
- Create: `scripts/generate_fortran_reference.py`, `tests/fixtures/fortran_reference/simple_physics.npz`, `tests/fixtures/fortran_reference/dcmip.npz`, `tests/fixtures/fortran_reference/cases.json`

**Interfaces:**
- Produces: the fixture format every later task reads.
  - `cases.json`: `{"simple_physics": {case: {"kwargs": {...}, "timestep_seconds": 600.0}}, "dcmip": {case: {"kwargs": {...}}}}`.
  - `simple_physics.npz`: keys `"<case>/in/<quantity>"` (the raw input arrays given to `array_call`, in the component's dims) and `"<case>/out/<quantity>"` (every diagnostic and every new-state quantity).
  - `dcmip.npz`: keys `"<case>/in/<quantity>"` and `"<case>/out/<quantity>"` likewise.

- [ ] **Step 1: Create the environment**

Run: `CLIMT_DEV_ENV_NAME=climt-scm scripts/bootstrap_dev_env.sh`
Then: `$PY -c "from climt._core.fortran import fortran_status; print(fortran_status()['state'])"`
Expected: `available` (climt editable from this worktree, climt-fortran 1.0 built from `./fortran`).

- [ ] **Step 2: Record the baseline**

Run: `$PY -m pytest -m "not slow" -q 2>&1 | tail -3`
Expected: all pass; write the counts into the Task 1 commit message.

- [ ] **Step 3: Write `scripts/generate_fortran_reference.py`**

Runs `SimplePhysics(**kwargs).array_call(raw_state, timedelta(seconds=600))` and `DcmipInitialConditions(**kwargs).array_call(raw_state)` on the cases below and writes the three files. Raw input states come from `climt.get_default_state([component], grid_state=climt.get_grid(nx=4, ny=2, nz=28))` (8 columns) for SimplePhysics and `get_grid(nx=16, ny=8, nz=10)` for dcmip, converted to raw arrays with the component's input properties (`sympl.get_numpy_arrays_with_properties`), then modified per case with a fixed `numpy.random.default_rng(20261119)`. It refuses to run unless `fortran_status()['state'] == 'available'`.

SimplePhysics cases (every column of every case differs: temperature profiles 220–305 K at the surface, humidity 0–1.2 × saturation):
`default` (constructor defaults); `lsc_off` (`large_scale_condensation=False`); `pbl_off` (`boundary_layer=False`); `fluxes_off` (`surface_fluxes=False`); `internal_sst` (`use_external_surface_temperature=False`); `cyclone` (`simulate_cyclone=True, use_external_surface_temperature=False`); `external_qsurf` (`use_external_surface_specific_humidity=True`); `drag_params` (all six drag and boundary-layer arguments at 1.5× their defaults); `lsc_supersaturated` (humidity 1.1–1.3 × saturation in levels 5–15 of every column); `strong_wind` (surface `eastward_wind` 15–30 m/s across columns); `cold_surface` (`surface_temperature` 250–275 K).

dcmip cases: `bw_dry` (`moist=False`), `bw_moist` (`moist=True`), `bw_unperturbed` (`moist=False, add_perturbation=False`), `tc` (`condition_type="tropical_cyclone", moist=True`).

- [ ] **Step 4: Generate and check**

Run: `$PY scripts/generate_fortran_reference.py && $PY -c "import numpy as np; z=np.load('tests/fixtures/fortran_reference/simple_physics.npz'); print(len(z.files))"`
Expected: the script exits 0; the key count is printed; no array in the files contains NaN.

- [ ] **Step 5: Commit**

```bash
git add scripts/generate_fortran_reference.py tests/fixtures/fortran_reference
git commit -m "test: Fortran reference outputs for SimplePhysics and dcmip"
```

---

### Task 2: The SimplePhysics kernel

**Files:**
- Create: `climt/_core/parameters.py`, `climt/_components/simple_physics/kernel.py`, `tests/test_parameters.py`, `tests/test_simple_physics_port.py`, `benchmarks/benchmark_simple_physics_port.py`, `benchmarks/results/simple_physics_port.txt`

**Interfaces:**
- Consumes: the Task 1 fixtures.
- Produces:

```python
# climt/_core/parameters.py
PARALLEL_MIN_COLUMNS: int = 32     # kernels use their parallel variant at or above this many columns

def per_column(params: NamedTuple, ncol: int) -> NamedTuple:
    """Same type; every field as an array of shape (ncol,): int64 for fields
    annotated int in the NamedTuple, float64 otherwise. A scalar or
    shape-(1,) field is broadcast (np.broadcast_to, read-only); a shape-(ncol,)
    field is passed through; any other shape raises ValueError naming the field."""

# climt/_components/simple_physics/kernel.py
class SimplePhysicsParams(NamedTuple):
    gravitational_acceleration: float              # m s-2
    heat_capacity_of_dry_air_at_constant_pressure: float   # J kg-1 K-1
    gas_constant_of_dry_air: float                 # J kg-1 K-1
    gas_constant_of_vapor_phase: float             # J kg-1 K-1
    latent_heat_of_condensation: float             # J kg-1
    planetary_radius: float                        # m
    planetary_rotation_rate: float                 # s-1
    density_of_liquid_water: float                 # kg m-3
    top_of_boundary_layer: float                   # Pa
    boundary_layer_influence_height: float         # Pa
    drag_coefficient_heat_fluxes: float
    base_momentum_drag_coefficient: float
    wind_dependent_momentum_drag_coefficient: float
    maximum_momentum_drag_coefficient: float

class SimplePhysicsFlags(NamedTuple):
    simulate_cyclone: bool
    large_scale_condensation: bool
    boundary_layer: bool
    surface_fluxes: bool
    use_external_surface_temperature: bool
    use_external_surface_specific_humidity: bool

def simple_physics_kernel(u, v, t, q, p, p_int, ps, ts, qs, lat, dt, params, flags):
    """u, v, t, q, p: (nlev, ncol); p_int: (nlev+1, ncol); ps, ts, qs, lat: (ncol,),
    climt's level order (index 0 nearest the surface), lat in degrees as the
    component receives it; dt seconds; params a SimplePhysicsParams of (ncol,)
    arrays (see per_column); flags a SimplePhysicsFlags.
    Returns (t_out, u_out, v_out, q_out, precip, sensible_heat_flux, latent_heat_flux),
    new arrays, inputs unmodified. latent_heat_flux is not clipped (the component clips)."""
```

- [ ] **Step 1: Write the failing tests**

`tests/test_parameters.py`:
```python
def test_per_column_broadcasts_scalars_and_keeps_arrays():
    P = namedtuple("P", "a b")
    out = per_column(P(a=2.0, b=np.arange(3.0)), 3)
    assert type(out) is P and out.a.shape == (3,) and out.a.dtype == np.float64
    # (a NamedTuple class with a field annotated int gets an int64 array for it)
    np.testing.assert_array_equal(out.a, [2.0, 2.0, 2.0]); np.testing.assert_array_equal(out.b, [0.0, 1.0, 2.0])

def test_per_column_rejects_wrong_length():
    with pytest.raises(ValueError, match="b"):
        per_column(namedtuple("P", "a b")(a=1.0, b=np.ones(2)), 3)
```

`tests/test_simple_physics_port.py` (a helper `run_kernel(case)` builds params and flags from `cases.json` kwargs and climt's constants, and calls the kernel on the fixture inputs; the component's clip of negative latent heat flux is applied to the kernel output before comparing):
```python
@pytest.mark.parametrize("case", SIMPLE_PHYSICS_CASES)   # all 11 from cases.json
def test_kernel_matches_fortran(case):
    out = run_kernel(case)
    for name, expected in fixture_outputs(case).items():
        np.testing.assert_allclose(out[name], expected, rtol=1e-10, atol=1e-12, err_msg=name)

def test_per_column_parameters_match_columns_run_alone():
    # case "default", base_momentum_drag_coefficient = linspace(0.0005, 0.0012, ncol)
    # run all columns at once with the array; run each column alone with its scalar;
    # every output agrees at rtol=1e-14.

def test_parallel_and_serial_variants_agree():
    # 64 columns (fixture columns tiled 8x): output with ncol >= PARALLEL_MIN_COLUMNS
    # equals the serial variant's at rtol=0, atol=0.
```

- [ ] **Step 2: Run them to see them fail**

Run: `$PY -m pytest tests/test_parameters.py tests/test_simple_physics_port.py -q`
Expected: FAIL with `ModuleNotFoundError` (`climt._core.parameters`, `...simple_physics.kernel`).

- [ ] **Step 3: Implement `per_column` and `PARALLEL_MIN_COLUMNS` in `climt/_core/parameters.py`**

- [ ] **Step 4: Implement `simple_physics_kernel` in `climt/_components/simple_physics/kernel.py`**

A line-by-line port of `simple_physics_func` in `fortran/_lib/simple_physics/simple_physics_custom.f90` (the file the Makefile compiles), with the constants from `set_physical_constants_func` replaced by `params` fields indexed per column, the Fortran module's switches replaced by `flags`, and the layer thickness and inverse thickness computed as `fortran/climt_fortran/_simple_physics.pyx::get_new_state` does. The Fortran works top level first: reverse the level axis on entry and on exit as the `.pyx` does. Write the column loop with `prange` in one Python function and compile it twice, `njit(cache=True)` and `njit(parallel=True, cache=True)`; `simple_physics_kernel` calls the parallel one when `ncol >= PARALLEL_MIN_COLUMNS`. Without numba (`HAS_NUMBA` false in `climt/_core/backend.py`) both names are the plain Python function.

- [ ] **Step 5: Run the tests**

Run: `$PY -m pytest tests/test_parameters.py tests/test_simple_physics_port.py -q`
Expected: PASS. A parity failure is a porting bug; tolerances are not loosened.

- [ ] **Step 6: Write and run the speed comparison**

`benchmarks/benchmark_simple_physics_port.py`: one column from fixture case `default`; the Fortran through `climt_fortran._simple_physics` (`init_simple_physics` and `set_physical_constants` once, then `get_new_state` timed); the numba kernel timed after one warm-up call; 2,000 calls each, best of 5; prints both in microseconds and the ratio, and writes the same lines plus machine, numba version and commit to `benchmarks/results/simple_physics_port.txt`.
Run: `NUMBA_NUM_THREADS=1 $PY benchmarks/benchmark_simple_physics_port.py`
Expected: ratio (numba / Fortran) ≤ 1.5. If above, profile the kernel and fix it before continuing; do not change the bound.

- [ ] **Step 7: Commit**

```bash
git add climt/_core/parameters.py climt/_components/simple_physics/kernel.py tests/test_parameters.py tests/test_simple_physics_port.py benchmarks/benchmark_simple_physics_port.py benchmarks/results/simple_physics_port.txt
git commit -m "feat: numba port of the simple physics kernel"
```

---

### Task 3: SimplePhysics on the kernel

**Files:**
- Modify: `climt/_components/simple_physics/component.py`
- Test: `tests/test_simple_physics_port.py`

**Interfaces:**
- Consumes: `simple_physics_kernel`, `SimplePhysicsParams`, `SimplePhysicsFlags` (Task 2), `per_column` (Task 2).
- Produces, on `SimplePhysics`:
  - `_params: SimplePhysicsParams` of scalars, built in `__init__` from the constructor arguments and `get_constant` (the eight constants under the names of the fields).
  - `_flags: SimplePhysicsFlags` built in `__init__`.
  - class attribute `per_column_parameters = SimplePhysicsParams._fields`.
  - `array_call(state, timestep)` calls `simple_physics_kernel(..., per_column(self._params, ncol), self._flags)`; clips negative latent heat flux to 0 as today; returns `(diagnostics, new_state)` as today.
  - No `load_extension`, no `unavailable_error`, no `_set_fortran_constants`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_simple_physics_port.py`)

```python
@pytest.mark.parametrize("case", SIMPLE_PHYSICS_CASES)
def test_component_matches_fortran(case):
    component = SimplePhysics(**kwargs(case))
    diagnostics, new_state = component.array_call(fixture_inputs(case), timedelta(seconds=600))
    for name, expected in fixture_outputs(case).items():
        actual = diagnostics[name] if name in diagnostics else new_state[name]
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12, err_msg=name)

def test_alternating_instances_keep_their_parameters():
    a, b = SimplePhysics(), SimplePhysics(base_momentum_drag_coefficient=0.0014)
    raw, dt = fixture_inputs("default"), timedelta(seconds=600)
    first = a.array_call(raw, dt); b.array_call(raw, dt); again = a.array_call(raw, dt)
    # every array in first equals the one in again at rtol=0

def test_non_contiguous_input_matches():
    # raw = fixture_inputs("default"); a copy where every 2-D array is
    # np.asfortranarray (same values, different memory order) gives identical outputs.

def test_runs_without_climt_fortran():
    # subprocess: sys.modules["climt_fortran"] = None before importing climt;
    # SimplePhysics()(state, timedelta(minutes=10)) on a default state succeeds.

def test_constants_read_once(monkeypatch):
    component = SimplePhysics()
    monkeypatch.setattr("climt._components.simple_physics.component.get_constant",
                        lambda *a, **k: pytest.fail("get_constant called in array_call"))
    component.array_call(fixture_inputs("default"), timedelta(seconds=600))
```

- [ ] **Step 2: Run them to see them fail**

Run: `$PY -m pytest tests/test_simple_physics_port.py -q -k "component or alternating or contiguous or without or once"`
Expected: FAIL (`test_constants_read_once` fails through `pytest.fail`; `test_runs_without_climt_fortran` fails with the unavailable error).

- [ ] **Step 3: Rewrite `climt/_components/simple_physics/component.py` as specified under Interfaces**

Keep the docstrings and `@ensure_contiguous_state`.

- [ ] **Step 4: Run the new tests and the existing SimplePhysics tests**

Run: `$PY -m pytest tests/test_simple_physics_port.py tests/test_components.py -q -k "SimplePhysics or simple_physics or component or alternating or contiguous or without or once"`
Expected: PASS, including `TestSimplePhysics` and `TestSimplePhysicsImplicitPrognostic` against their cached outputs.

- [ ] **Step 5: Commit**

```bash
git add climt/_components/simple_physics/component.py tests/test_simple_physics_port.py
git commit -m "feat: SimplePhysics runs on the numba kernel, constants read once"
```

---

### Task 4: dcmip in numba

**Files:**
- Create: `climt/_components/dcmip/kernels.py`, `tests/test_dcmip_port.py`
- Modify: `climt/_components/dcmip/component.py`

**Interfaces:**
- Consumes: the Task 1 fixtures.
- Produces:

```python
# climt/_components/dcmip/kernels.py
def baroclinic_wave_ics(pressure, longitude, latitude, perturb, moist):
    """pressure (nlev, ncol) Pa; longitude, latitude (ncol,) radians; perturb, moist bool.
    Returns (u, v, t, q, surface_pressure, surface_geopotential): (nlev, ncol) x4, (ncol,) x2."""

def tropical_cyclone_ics(pressure, longitude, latitude):
    """Same shapes and return order; the tropical cyclone test is always moist."""
```

`DcmipInitialConditions.array_call` calls these in place of `_dcmip.get_baroclinic_wave_ics` / `get_tropical_cyclone_ics`; the rest of `array_call` is unchanged. No `load_extension` in the module.

- [ ] **Step 1: Write the failing tests** (`tests/test_dcmip_port.py`)

```python
@pytest.mark.parametrize("case", DCMIP_CASES)   # bw_dry, bw_moist, bw_unperturbed, tc
def test_component_matches_fortran(case):
    out = DcmipInitialConditions(**kwargs(case)).array_call(fixture_inputs(case))
    for name, expected in fixture_outputs(case).items():
        np.testing.assert_allclose(out[name], expected, rtol=1e-12, atol=1e-12, err_msg=name)

def test_runs_without_climt_fortran():
    # subprocess, as in Task 3, with DcmipInitialConditions on a default state.
```

- [ ] **Step 2: Run them to see them fail**

Run: `$PY -m pytest tests/test_dcmip_port.py -q`
Expected: FAIL (`ModuleNotFoundError` for the kernels module once the test imports it; the subprocess test fails with the unavailable error).

- [ ] **Step 3: Implement the kernels**

Ports of `test4_baroclinic_wave` (`fortran/_lib/dcmip/dcmip_initial_conditions_test_4_v3.f90`, with `X = 1` and `zcoords = 0` as the `.pyx` passes them) and `test5_tropical_cyclone` (`..._test_5_v1.f90`), each a per-point function called over levels and columns as `fortran/climt_fortran/_dcmip.pyx` does. numba `njit` (serial: dcmip runs once). Only tests 4 and 5 are ported; `..._test_1_2_3_v5.f90` is not built today and is not ported.

- [ ] **Step 4: Switch `climt/_components/dcmip/component.py` to the kernels**

- [ ] **Step 5: Run the tests and the existing dcmip tests**

Run: `$PY -m pytest tests/test_dcmip_port.py tests/test_components.py tests/test_initialization.py -q -k "dcmip or Dcmip or component or without"`
Expected: PASS, including `TestDcmip`, `test_dcmip_options` and the cached outputs.

- [ ] **Step 6: Commit**

```bash
git add climt/_components/dcmip tests/test_dcmip_port.py
git commit -m "feat: numba port of the DCMIP initial conditions"
```

---

### Task 5: Stop loading the two Fortran modules; changelog

**Files:**
- Modify: `climt/_core/fortran.py`, `tests/test_fortran_loader.py`, `HISTORY.rst`, and any page under `docs/` that says SimplePhysics or DcmipInitialConditions needs climt-fortran (find them with `grep -rln "SimplePhysics\|Dcmip" docs`)

**Interfaces:**
- Produces: `climt._core.fortran.MODULES == ("_emanuel_convection", "_rrtmg_lw", "_rrtmg_sw")`. climt then reports `available` with climt-fortran 1.0 (which still ships the two modules, unused) and with a future 2.0 without them.

- [ ] **Step 1: Write the failing test** (append to `tests/test_fortran_loader.py`)

```python
def test_ported_modules_are_not_loaded():
    assert "_simple_physics" not in MODULES and "_dcmip" not in MODULES
```

- [ ] **Step 2: Run it to see it fail**

Run: `$PY -m pytest tests/test_fortran_loader.py -q`
Expected: FAIL on the new test.

- [ ] **Step 3: Change `MODULES`; update any loader test that lists five modules; update the docs pages found above and add to `HISTORY.rst` under Unreleased: SimplePhysics and DcmipInitialConditions now run on numba and no longer need climt-fortran (so they also work in Pyodide); SimplePhysics no longer re-reads its constants on every call.**

- [ ] **Step 4: Run the full suite and flake8**

Run: `$PY -m pytest -m "not slow" -q 2>&1 | tail -3` and `$PY -m flake8 climt`
Expected: the Task 1 baseline count plus the new tests, all passing; flake8 clean.

- [ ] **Step 5: Commit**

```bash
git add climt/_core/fortran.py tests/test_fortran_loader.py HISTORY.rst docs
git commit -m "chore: climt no longer loads the simple physics and dcmip Fortran modules"
```
