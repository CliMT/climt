# Differentiable climt Branch Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a `differentiable-climt` branch where each component has a JAX kernel alongside its Numba kernel, enabling gradient computation through climt components for use in calibration, parameter estimation, and ML integration.

**Architecture:** Each component dispatches to either a Numba kernel (`_kernel_np`) or a JAX kernel (`_kernel_jax`) based on the array namespace of its inputs (detected via `get_array_namespace()`). The JAX kernels are pure functional rewrites using `jnp` — vectorized broadcasting instead of explicit loops, `jnp.maximum` instead of `max()`, no in-place mutation. A `JaxBackend` + `JaxStateContainer` provide the sympl-compatible interface for constructing JAX-array states. `jax.grad` / `jax.jacobian` then work through `array_call()`.

**Tech Stack:** Python, JAX (`jax`, `jax.numpy`), NumPy, Numba, sympl, pytest, conda env `climt`

---

## Prerequisites

Complete `2026-03-17-numba-branch-finalization.md` first (or at least have that branch with all tests passing). This branch starts from `numba-optimized-components` after it is clean.

---

## Background: How the JAX Dispatch Works

From the commit history (`a9e4486`), the pattern is:

```python
def array_call(self, raw_state):
    t = raw_state["air_temperature"]
    xp = get_array_namespace(t)   # returns np or jnp module
    u = raw_state["eastward_wind"]
    # ...extract other fields...
    if xp is np:
        result = _component_kernel_np(...)   # @njit decorated
    else:
        result = _component_kernel_jax(...)  # pure jnp, differentiable
    return result
```

`get_array_namespace()` lives in `climt/_core/backend.py`. It inspects the array type and returns the right module (`numpy` vs `jax.numpy`).

**JAX kernel rules:**
- No in-place mutation (`out[i, j] = ...` → return new array)
- No Python `max()`/`min()` → use `jnp.maximum()`/`jnp.minimum()`
- No explicit loops over columns → use vectorized broadcasting or `jax.vmap`
- NamedTuple params work unchanged (JAX traces through them)
- `jax.config.update("jax_enable_x64", True)` needed for float64 parity

---

## File Structure

**New files:**
- `climt/_core/jax_backend.py` — `JaxBackend`, `JaxStateContainer`, `JaxTimeDelta`, PyTree registration

**Modified files (each adds a `_kernel_jax` function and dispatch in `array_call`):**
- `climt/_components/held_suarez.py`
- `climt/_components/grid_scale_condensation.py`
- `climt/_components/radiation.py` (both `GrayLongwaveRadiation` and `Frierson06LongwaveOpticalDepth`)
- `climt/_components/dry_convection/component.py`
- `climt/_components/berger_solar_insolation.py`
- `climt/_components/slab_surface.py`
- `climt/_components/instellation/component.py`
- `climt/_components/emanuel/pure_python_v3.py`
- `climt/_core/backend.py` — add `get_array_namespace()`
- `climt/_core/__init__.py` — export `JaxBackend`, `JaxStateContainer`, `JaxTimeDelta`
- `climt/__init__.py` — export `JaxBackend`, `JaxStateContainer`, `JaxTimeDelta`

**New test files:**
- `tests/test_jax_differentiation.py` — gradient tests for each JAX-enabled component

---

## Task 1: Create the Branch

- [ ] **Step 1: Create `differentiable-climt` from the finalized numba branch**

```bash
git checkout numba-optimized-components
git checkout -b differentiable-climt
```

- [ ] **Step 2: Confirm base state**

```bash
conda run -n climt pytest tests/ -q 2>&1 | tail -5
```

Expected: 0 failures.

---

## Task 2: Restore `jax_backend.py`

**Files:**
- Create: `climt/_core/jax_backend.py`

The implementation already exists in the git history at commit `a9e4486`. It provides:
- `JaxTimeDelta` — timedelta subclass whose `total_seconds()` returns a JAX scalar
- `JaxStateContainer` — array wrapper with `.data`, `.dims`, `.units`, `.attrs`, arithmetic ops, and dimension alignment
- `JaxBackend` — implements sympl `StateBackend` interface; `get_array()` handles dim permutation, `create_quantity()` returns a `JaxStateContainer`
- PyTree registration so JAX can trace through `JaxStateContainer`

- [ ] **Step 1: Restore from git history**

```bash
git show a9e4486:climt/_core/jax_backend.py > climt/_core/jax_backend.py
```

- [ ] **Step 2: Read the restored file and verify it imports cleanly**

```bash
conda run -n climt python -c "from climt._core.jax_backend import JaxBackend, JaxStateContainer, JaxTimeDelta; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Restore exports**

In `climt/_core/__init__.py`, add:
```python
from .jax_backend import JaxBackend, JaxStateContainer, JaxTimeDelta
```

In `climt/__init__.py`, add to imports and `__all__`:
```python
from ._core import (
    ...
    JaxBackend,
    JaxStateContainer,
    JaxTimeDelta,
    ...
)
```

- [ ] **Step 4: Verify public import**

```bash
conda run -n climt python -c "from climt import JaxBackend, JaxStateContainer, JaxTimeDelta; print('OK')"
```

- [ ] **Step 5: Commit**

```bash
git add climt/_core/jax_backend.py climt/_core/__init__.py climt/__init__.py
git commit -m "feat: restore JaxBackend and JaxStateContainer infrastructure

Provides JAX-array-compatible state container and sympl backend
needed for gradient computation through component kernels.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 3: Add `get_array_namespace` to Backend

**Files:**
- Modify: `climt/_core/backend.py`

- [ ] **Step 1: Read the current `backend.py`**

```bash
cat climt/_core/backend.py
```

- [ ] **Step 2: Add `get_array_namespace()`**

The function accepts a single array-like value. It must handle:
- plain `np.ndarray` → return `np`
- `jax.Array` → return `jnp`
- `JaxStateContainer` (has `__jax_array__` and `.data` that is a `jax.Array`) → return `jnp`
- `xarray.DataArray` (has `.values` that is `np.ndarray`) → return `np`

```python
def get_array_namespace(arr):
    """Return numpy or jax.numpy based on the array type.

    Checks for __jax_array__ protocol first so JaxStateContainer is handled
    before inspecting .data (which might not exist on plain arrays).
    """
    try:
        import jax
        import jax.numpy as jnp
        # Catches JaxStateContainer and any JAX tracer/array
        if hasattr(arr, '__jax_array__'):
            return jnp
        # Plain jax.Array
        if isinstance(arr, jax.Array):
            return jnp
        # Container wrapping a JAX array (e.g. sympl DataArray .data is still np)
        data = getattr(arr, 'data', None)
        if data is not None and isinstance(data, jax.Array):
            return jnp
    except ImportError:
        pass
    return np
```

> **Note:** The historical `a9e4486` version is variadic (`*arrays`). This single-argument version is intentionally simpler — it is sufficient because each component's `array_call` inspects a single representative array (typically `air_temperature`).

- [ ] **Step 3: Verify — must test all four cases**

```bash
conda run -n climt python -c "
import numpy as np
import jax.numpy as jnp
from climt._core.backend import get_array_namespace
from climt import JaxStateContainer

arr_np   = np.array([1.0])
arr_jax  = jnp.array([1.0])
arr_jsc  = JaxStateContainer(jnp.array([1.0]), ('x',))

assert get_array_namespace(arr_np)  is np,  'numpy array failed'
assert get_array_namespace(arr_jax) is jnp, 'jax.Array failed'
assert get_array_namespace(arr_jsc) is jnp, 'JaxStateContainer failed'
print('OK')
"
```

- [ ] **Step 4: Commit**

```bash
git add climt/_core/backend.py
git commit -m "feat: add get_array_namespace for numpy/JAX dispatch

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 4: Add JAX Kernel to `HeldSuarez` (Reference Implementation)

`HeldSuarez` is the simplest component — a good first target. The JAX kernel already exists verbatim in the git history at `a9e4486`.

**Files:**
- Modify: `climt/_components/held_suarez.py`
- Test: `tests/test_jax_differentiation.py` (create)

- [ ] **Step 1: Write the failing test first**

Create `tests/test_jax_differentiation.py`:

```python
# -*- coding: utf-8 -*-
import pytest
import numpy as np

# IMPORTANT: jax_enable_x64 must be set before any JAX computation.
# This module-level call is safe as long as pytest collects this file
# before any other test file that runs JAX computations. Run this file
# in isolation if that becomes a problem: pytest tests/test_jax_differentiation.py
jax = pytest.importorskip("jax")
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from climt import get_grid, get_default_state, HeldSuarez, JaxBackend, JaxStateContainer


def _make_jax_state(component, grid):
    """Build a JAX state by wrapping all sympl DataArray values in JaxStateContainer.

    Uses canonical sympl names (from get_default_state), not component aliases.
    The component's input_properties may define aliases, but get_default_state
    always returns canonical names. array_call receives pre-aliased arrays from
    sympl's dispatch; when called directly we must use the aliased keys.
    See Task 4 note below.
    """
    sympl_state = get_default_state([component], grid_state=grid)
    return {
        k: JaxStateContainer(jnp.array(v.values), v.dims)
        for k, v in sympl_state.items()
        if k != 'time'
    }


def test_held_suarez_jax_parity():
    """JAX kernel produces same output as Numba kernel."""
    import sympl
    sympl.set_backend(sympl.DataArrayBackend())

    grid = get_grid(nx=4, ny=4, nz=28)
    comp = HeldSuarez()
    state_sympl = get_default_state([comp], grid_state=grid)
    tend_np, _ = comp(state_sympl)

    jax_state = _make_jax_state(comp, grid)
    tend_jax, _ = comp.array_call(jax_state)

    t_np = tend_np['air_temperature'].values
    t_jax = np.array(tend_jax['air_temperature'].data)
    assert np.allclose(t_np, t_jax, rtol=1e-5), "JAX and Numba outputs diverge"


def test_held_suarez_gradient():
    """jax.grad works through HeldSuarez.array_call."""
    grid = get_grid(nx=2, ny=2, nz=28)
    comp = HeldSuarez()
    jax_state = _make_jax_state(comp, grid)

    def loss(t):
        state = {**jax_state, 'air_temperature': JaxStateContainer(t, ('mid_levels', 'lat', 'lon'))}
        tend, _ = comp.array_call(state)
        return jnp.sum(tend['air_temperature'].data ** 2)

    t0 = jax_state['air_temperature'].data
    grad = jax.grad(loss)(t0)

    assert grad.shape == t0.shape
    assert not jnp.any(jnp.isnan(grad)), "Gradient contains NaN"
    assert jnp.any(grad != 0.0), "Gradient is all zeros"
```

- [ ] **Step 2: Run — confirm FAIL (JAX kernel not yet added)**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_held_suarez_jax_parity -v 2>&1 | tail -10
```

Expected: FAIL (dispatch goes to numpy path, not JAX).

- [ ] **Step 3: Add JAX kernel and dispatch to `held_suarez.py`**

Read `climt/_components/held_suarez.py`, then add at the end:

```python
def _held_suarez_kernel_jax(u, v, t, p, ps, lat, params):
    import jax.numpy as jnp
    ps_expanded = jnp.expand_dims(ps, -1)
    lat_expanded = jnp.expand_dims(lat, -1)
    sigma = p / ps_expanded
    lat_rad = jnp.deg2rad(lat_expanded)
    p_norm = p / params.p0
    Teq = jnp.maximum(
        200.0,
        (315.0 - params.delta_T_y * jnp.sin(lat_rad)**2
         - params.delta_theta_z * jnp.log(p_norm) * jnp.cos(lat_rad)**2)
        * p_norm**params.kappa
    )
    sigma_fac = jnp.maximum(
        0.0,
        (sigma - params.sigma_b) / (1.0 - params.sigma_b)
    )
    k_t = params.k_a + (params.k_s - params.k_a) * sigma_fac * jnp.cos(lat_rad)**4
    k_v = params.k_f * sigma_fac
    return -k_v * u, -k_v * v, -k_t * (t - Teq)
```

Update `array_call` to dispatch:

```python
def array_call(self, raw_state):
    from .._core.backend import get_array_namespace
    t = raw_state["air_temperature"]
    xp = get_array_namespace(t)
    u = raw_state["eastward_wind"]; v = raw_state["northward_wind"]
    p = raw_state["air_pressure"]; ps = raw_state["surface_air_pressure"]
    lat = raw_state["latitude"]
    # extract raw array from JaxStateContainer or DataArray
    def _d(x): return getattr(x, 'data', x)
    if xp is np:
        # NumPy/Numba path — return plain arrays so sympl can wrap them normally
        tend_u, tend_v, tend_t = _held_suarez_kernel_np(
            _d(u), _d(v), _d(t), _d(p), _d(ps), _d(lat), self._params)
        return (
            {"eastward_wind": tend_u,
             "northward_wind": tend_v,
             "air_temperature": tend_t},
            {}
        )
    else:
        # JAX path — wrap results in JaxStateContainer preserving dims
        from .._core.jax_backend import JaxStateContainer as JSC
        tend_u, tend_v, tend_t = _held_suarez_kernel_jax(
            _d(u), _d(v), _d(t), _d(p), _d(ps), _d(lat), self._params)
        def _wrap(arr, src):
            dims = getattr(src, 'dims', getattr(src, 'dims', None))
            return JSC(arr, dims) if dims is not None else arr
        return (
            {"eastward_wind": _wrap(tend_u, u),
             "northward_wind": _wrap(tend_v, v),
             "air_temperature": _wrap(tend_t, t)},
            {}
        )
```

> **Important:** The NumPy path returns plain `np.ndarray` values (not wrapped). This is required for sympl's normal result-handling code to work. Only the JAX path returns `JaxStateContainer`.

> **Note on `array_call` vs component `__call__`:** When calling `comp(state_sympl)`, sympl translates `input_properties` aliases before passing to `array_call`. When calling `comp.array_call(jax_state)` directly in tests, you must use the *aliased* key names that the component declares (e.g. `"u"` if the component aliases `"eastward_wind"` → `"u"`). Check `component.input_properties` for the aliases. `HeldSuarez` does not use aliases, so canonical names work for direct `array_call` testing.

- [ ] **Step 4: Run parity test — should PASS**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_held_suarez_jax_parity -v 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Run gradient test — should PASS**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_held_suarez_gradient -v 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 6: Verify Numba path still works**

```bash
conda run -n climt pytest tests/test_held_suarez_optimization.py tests/test_components.py -k "HeldSuarez" -v 2>&1 | tail -10
```

Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add climt/_components/held_suarez.py tests/test_jax_differentiation.py
git commit -m "feat: add JAX kernel to HeldSuarez with numpy/JAX dispatch

Dispatches to @njit kernel for NumPy arrays and jnp vectorized
kernel for JAX arrays. Enables jax.grad through array_call.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 5: Add JAX Kernel to `GrayLongwaveRadiation`

**Files:**
- Modify: `climt/_components/radiation.py`
- Modify: `tests/test_jax_differentiation.py`

The Numba kernel iterates over columns and levels computing fluxes from the top-of-atmosphere downward and surface upward. The JAX version uses cumulative products and sums.

- [ ] **Step 1: Write the failing test (add to `test_jax_differentiation.py`)**

```python
def test_gray_radiation_gradient():
    """jax.grad works through GrayLongwaveRadiation.array_call."""
    from climt import GrayLongwaveRadiation, Frierson06LongwaveOpticalDepth
    import sympl
    sympl.set_backend(sympl.DataArrayBackend())
    grid = get_grid(nx=2, ny=2, nz=28)
    tau_comp = Frierson06LongwaveOpticalDepth()
    lw_comp = GrayLongwaveRadiation()
    state_sympl = get_default_state([tau_comp, lw_comp], grid_state=grid)
    state_sympl.update(tau_comp(state_sympl))

    # Build JAX state using CANONICAL names (what get_default_state returns),
    # not the aliases declared in GrayLongwaveRadiation.input_properties.
    # GrayLongwaveRadiation uses aliases: tau, sl, T_surface, p, p_interface.
    # When calling array_call directly, pass the ALIASED keys; when using the
    # component via __call__, sympl handles the translation.
    # Here we call lw_comp.array_call() directly, so we must use aliased names.
    # Read lw_comp.input_properties to get alias → canonical mapping, then build:
    alias_map = {
        v.get('alias', k): k
        for k, v in lw_comp.input_properties.items()
    }
    # alias_map is e.g. {'tau': 'longwave_optical_depth_on_interface_levels', ...}
    # Build jax_state with aliased keys:
    jax_state = {}
    for alias, canonical in alias_map.items():
        da = state_sympl[canonical]
        jax_state[alias] = JaxStateContainer(jnp.array(da.values), da.dims)

    def loss(temp):
        s = {**jax_state, 'sl': JaxStateContainer(temp, jax_state['sl'].dims)}
        tend, _ = lw_comp.array_call(s)
        return jnp.sum(tend['air_temperature'].data ** 2)

    t0 = jax_state['sl'].data
    grad = jax.grad(loss)(t0)
    assert not jnp.any(jnp.isnan(grad)), "NaN in gradient"
    assert jnp.any(grad != 0.0), "Zero gradient"
```

> **Alias note:** `GrayLongwaveRadiation.input_properties` declares aliases (`tau`, `sl`, `T_surface`, `p`, `p_interface`). When calling `array_call` directly (bypassing sympl), these aliased keys are what `array_call` receives. The test must construct the state with aliased keys accordingly.

- [ ] **Step 2: Run — confirm FAIL**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_gray_radiation_gradient -v 2>&1 | tail -10
```

- [ ] **Step 3: Read the Numba radiation kernel before writing the JAX version**

```bash
grep -n "def _gray_lw_kernel_np\|def _frierson_tau_kernel_np\|def array_call" \
  climt/_components/radiation.py
```

Read the full kernel to understand the exact index conventions:
- Index 0 = **surface** (bottom), index `nlev` = **TOA** (top) — verify this carefully
- Note which constants are retrieved via `get_constant()` vs. hardcoded

**⚠️ Critical:** The JAX kernel must match the Numba kernel's array index convention exactly. Verify the convention by reading the kernel, not by assuming. The parity test will catch mismatches.

- [ ] **Step 4: Add JAX kernel to `radiation.py`**

The JAX kernel uses `jax.lax.scan` for sequential level accumulation. The sketch below follows the common convention where index 0 = surface. **Adjust if the Numba kernel uses 0 = TOA.**

Use `get_constant()` for physical constants, not hardcoded values:

```python
def _gray_lw_kernel_jax(tau, sl, T_surface, p_interface):
    """
    tau:         (nlev+1, ncol) optical depth at interfaces; index 0 = surface
    sl:          (nlev, ncol) air temperature at mid-levels; index 0 = surface
    T_surface:   (ncol,) surface temperature
    p_interface: (nlev+1, ncol) pressure at interfaces
    Returns: lw_down, lw_up (nlev+1, ncol), tend_T (nlev, ncol)
    All arrays indexed with 0 = surface unless Numba kernel uses opposite convention.
    """
    import jax.numpy as jnp
    import jax.lax as lax
    from sympl import get_constant
    sigma = get_constant('stefan_boltzmann_constant', 'W/m^2/K^4')
    Cpd   = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
    g     = get_constant('gravitational_acceleration', 'm/s^2')

    dtau = jnp.abs(jnp.diff(tau, axis=0))  # (nlev, ncol), positive definite

    B      = sigma * sl ** 4        # (nlev, ncol)
    B_surf = sigma * T_surface ** 4  # (ncol,)

    ncol = sl.shape[1]

    # --- Downwelling (TOA → surface) ---
    # Scan from TOA (index nlev-1) downward to surface (index 0)
    # dtau_toa_first: reverse so scan processes TOA first
    dtau_dn = dtau[::-1]   # (nlev, ncol), TOA first
    B_dn    = B[::-1]

    def scan_down(carry, x):
        dtau_k, B_k = x
        trans = jnp.exp(-dtau_k)
        new_flux = carry * trans + B_k * (1.0 - trans)
        return new_flux, new_flux

    _, lw_down_reversed = lax.scan(
        scan_down,
        jnp.zeros(ncol),    # zero incident flux at TOA
        (dtau_dn, B_dn)
    )
    # lw_down_reversed[0] = flux just below TOA, [nlev-1] = flux at surface
    # Restore surface-first order and prepend the TOA zero
    lw_down = jnp.concatenate([
        lw_down_reversed[::-1],
        jnp.zeros((1, ncol))   # TOA boundary: zero
    ], axis=0)  # (nlev+1, ncol), index 0 = surface

    # --- Upwelling (surface → TOA) ---
    def scan_up(carry, x):
        dtau_k, B_k = x
        trans = jnp.exp(-dtau_k)
        new_flux = carry * trans + B_k * (1.0 - trans)
        return new_flux, new_flux

    _, lw_up_levels = lax.scan(
        scan_up,
        B_surf,            # surface emission
        (dtau, B)          # surface first
    )
    lw_up = jnp.concatenate([
        B_surf[jnp.newaxis, :],  # surface interface
        lw_up_levels
    ], axis=0)  # (nlev+1, ncol), index 0 = surface

    # --- Temperature tendency ---
    dp    = jnp.abs(jnp.diff(p_interface, axis=0))  # (nlev, ncol)
    d_net = jnp.diff(lw_down - lw_up, axis=0)        # (nlev, ncol)
    tend_T = g / Cpd * d_net / dp

    return lw_down, lw_up, tend_T
```

> **⚠️ Warning:** The scan direction sketched above assumes index 0 = surface. After implementing, run the parity test (Step 5). If it fails, read the Numba kernel's scan direction carefully and flip the downwelling scan accordingly. The parity test is the ground truth.

Update `array_call` to dispatch using `get_array_namespace`, following the same NumPy-path / JAX-path pattern as `HeldSuarez`.

- [ ] **Step 5: Run parity test (JAX vs Numba) — should PASS**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_gray_radiation_parity -v 2>&1 | tail -10
```

(Add `test_gray_radiation_parity` to the test file following the same pattern as `test_held_suarez_jax_parity`.)

- [ ] **Step 6: Run gradient test — should PASS**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_gray_radiation_gradient -v 2>&1 | tail -10
```

- [ ] **Step 7: Verify Numba path unchanged**

```bash
conda run -n climt pytest tests/test_gray_radiation_optimization.py tests/test_components.py -k "GrayLongwave" -v 2>&1 | tail -10
```

- [ ] **Step 8: Commit**

```bash
git add climt/_components/radiation.py tests/test_jax_differentiation.py
git commit -m "feat: add JAX kernel to GrayLongwaveRadiation

Uses lax.scan for sequential level accumulation to remain
fully differentiable through jax.grad.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 6: Add JAX Kernels to Remaining Components

Repeat the same TDD pattern (failing test → JAX kernel → dispatch → pass → commit) for each remaining component. Each component follows the same architectural pattern.

### 6a: `GridScaleCondensation`

**Key JAX translation challenges:**
- Condensation loop iterates levels top-to-bottom — each level's q/T depends on condensate removed in the level above → sequential dependency → `jax.lax.scan`
- `q_sat` calculation: replace the Bolton formula with `jnp` equivalents (no Python math module)
- Scan carry state: `(q, T, precip)` — humidity, temperature, and accumulated precipitation at each step

**Scan carry structure:**
```python
# carry = (q_col, T_col, precip_col)  — shape (ncol,) each
# xs    = (q_sat_col, dtau_col, ...)  — stacked level slices
def scan_fn(carry, x):
    q, T, precip = carry
    # ... condensation logic for one level using jnp.where for conditionals ...
    return (new_q, new_T, new_precip), (new_q, new_T)
```

**Files:** `climt/_components/grid_scale_condensation.py`

- [ ] Read `_gsc_kernel_np` to confirm top-to-bottom iteration and carry state
- [ ] Add gradient test for GSC
- [ ] Run failing test
- [ ] Add `_gsc_kernel_jax()` + dispatch in `array_call`
- [ ] Run tests (parity + gradient)
- [ ] Commit

### 6b: `DryConvectiveAdjustment`

**Key JAX translation challenges:**
- Convective adjustment is iterative (convergence loop) → use `jax.lax.while_loop` or fixed iteration count
- Consider a fixed max-iteration version for differentiability

**Files:** `climt/_components/dry_convection/component.py`

- [ ] Add gradient test for DryConvection
- [ ] Run failing test
- [ ] Add `_dry_adj_kernel_jax()` with `jax.lax.while_loop` or fixed-iteration unrolling
- [ ] Run tests
- [ ] Commit

### 6c: `BergerSolarInsolation`

**Key JAX translation challenges:**
- Mostly trigonometric operations → straightforward `jnp` translation
- No sequential dependencies

**Files:** `climt/_components/berger_solar_insolation.py`

- [ ] Add gradient test
- [ ] Run failing test
- [ ] Add `_get_solar_parameters_jax()` + dispatch
- [ ] Run tests
- [ ] Commit

### 6d: `SlabSurface`

**Key JAX translation challenges:**
- Surface type indexing (land/sea/ice flags) → use `jnp.where` instead of integer indexing

**Files:** `climt/_components/slab_surface.py`

- [ ] Add gradient test
- [ ] Run failing test
- [ ] Add `_slab_surface_kernel_jax()` + dispatch
- [ ] Run tests
- [ ] Commit

### 6e: `Instellation`

**Key JAX translation challenges:**
- Multiple helper functions with trigonometric operations → all convert cleanly to `jnp`

**Files:** `climt/_components/instellation/component.py`

- [ ] Add gradient test
- [ ] Run failing test
- [ ] Add `_instellation_kernel_jax()` + dispatch
- [ ] Run tests
- [ ] Commit

### 6f: `EmanuelConvectionPythonV3`

**Key JAX translation challenges:**
- Most complex: nested loops with conditional logic
- Use `jax.lax.cond` for branches, `jax.lax.scan` for level loops
- May need a simplified/approximated JAX version

**Files:** `climt/_components/emanuel/pure_python_v3.py`

- [ ] Add gradient test (using the existing `test_jax_differentiation.py` pattern from commit `a9e4486`)
- [ ] Run failing test
- [ ] Add JAX kernel (consult prior implementation from `git show a9e4486:climt/_components/emanuel/pure_python_v3.py`)
- [ ] Run tests
- [ ] Commit

---

## Task 7: Integration Test — Full Pipeline Gradient

This verifies that a multi-component pipeline is end-to-end differentiable.

**Files:**
- Modify: `tests/test_jax_differentiation.py`

- [ ] **Step 1: Write the integration test**

```python
def test_full_pipeline_gradient():
    """
    jax.grad works through a radiation + held_suarez pipeline.
    Simulates parameter calibration: gradient of loss w.r.t. initial temperature.
    """
    from climt import GrayLongwaveRadiation, Frierson06LongwaveOpticalDepth, HeldSuarez
    import sympl
    sympl.set_backend(sympl.DataArrayBackend())

    grid = get_grid(nx=2, ny=2, nz=28)
    tau_comp = Frierson06LongwaveOpticalDepth()
    lw_comp = GrayLongwaveRadiation()
    hs_comp = HeldSuarez()

    state_sympl = get_default_state([tau_comp, lw_comp, hs_comp], grid_state=grid)
    state_sympl.update(tau_comp(state_sympl))

    all_keys = set(lw_comp.input_properties) | set(hs_comp.input_properties)
    base_jax_state = {
        k: JaxStateContainer(jnp.array(state_sympl[k].values), state_sympl[k].dims)
        for k in all_keys if k in state_sympl
    }

    def loss(temp):
        s = {**base_jax_state, 'air_temperature': JaxStateContainer(temp, ('mid_levels', 'lat', 'lon'))}
        tend_lw, _ = lw_comp.array_call(s)
        tend_hs, _ = hs_comp.array_call(s)
        combined = tend_lw['air_temperature'].data + tend_hs['air_temperature'].data
        return jnp.sum(combined ** 2)

    t0 = base_jax_state['air_temperature'].data
    grad = jax.grad(loss)(t0)
    assert not jnp.any(jnp.isnan(grad))
    assert jnp.any(grad != 0.0)
    print(f"Pipeline gradient norm: {jnp.linalg.norm(grad):.4f}")
```

- [ ] **Step 2: Run — should PASS (if Tasks 4 and 5 are complete)**

```bash
conda run -n climt pytest tests/test_jax_differentiation.py::test_full_pipeline_gradient -v 2>&1 | tail -10
```

- [ ] **Step 3: Run full test suite to check no regressions**

```bash
conda run -n climt pytest tests/ -q 2>&1 | tail -10
```

Expected: 0 failures.

- [ ] **Step 4: Commit**

```bash
git add tests/test_jax_differentiation.py
git commit -m "test: add full pipeline gradient integration test

Verifies end-to-end differentiability through radiation + held_suarez
pipeline, representative of calibration/ML use cases.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 8: Document the Differentiable Path

**Files:**
- Create: `docs/differentiable_climt.rst`

- [ ] **Step 1: Write documentation showing how to use JAX for calibration**

The document should cover:
1. Installing JAX (`pip install jax`)
2. Enabling float64 (`jax.config.update("jax_enable_x64", True)`)
3. Building a JAX state with `JaxStateContainer`
4. Using `jax.grad` for parameter sensitivity
5. Simple calibration example using `jax.scipy.optimize.minimize`

- [ ] **Step 2: Commit**

```bash
git add docs/differentiable_climt.rst
git commit -m "docs: add differentiable climt usage guide

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Notes

### JAX Kernel Writing Guide

**Replacing explicit loops:**
- Column loops (`for i in prange(ncol)`) → vectorize with broadcasting or `jax.vmap`
- Level loops with sequential dependence → `jax.lax.scan`
- Conditional branches (`if condition`) → `jax.lax.cond`

**Replacing in-place mutation:**
- `out[i, j] = value` → return `value` from scan carry or use `jnp.where`

**Replacing Python builtins:**
- `max(a, b)` → `jnp.maximum(a, b)`
- `min(a, b)` → `jnp.minimum(a, b)`
- `abs(a)` → `jnp.abs(a)`

**Numerical parity:**
- Always call `jax.config.update("jax_enable_x64", True)` in tests
- Results should match Numba to `rtol=1e-5` or better

### On `DryConvectiveAdjustment` Differentiability

Convective adjustment uses a convergence loop that can be tricky for autodiff. Options:
1. Use `jax.lax.while_loop` — not differentiable through the loop itself
2. Use `jax.lax.fori_loop` with fixed max iterations — differentiable
3. Implicit differentiation using the fixed-point theorem — most accurate but complex

Start with option 2 (fixed max iterations). Note the limitation in tests.

### On `EmanuelConvectionPythonV3`

The full Emanuel scheme has complex conditional logic. Consult:
```bash
git show a9e4486:climt/_components/emanuel/pure_python_v3.py
```
for the prior JAX implementation. It used `jax.lax.cond` for key branches and `jax.lax.scan` for level loops.
