# Boundary Layer + Canonical JIT Tridiagonal Solver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the Frierson-2006 `SimpleBoundaryLayer` component (with the thesis Eqn 2.8 `Ria` modification) built on a new canonical numba-jitable tridiagonal solver in `_core`, and remove `scipy` from the `snow_ice_column` and `second_best` diffusion paths in favour of `prange` column loops.

**Architecture:** One `@jit_compile` Thomas solver (`climt/_core/tridiagonal.py`) is the single dense tridiagonal solver for the whole library. The boundary-layer component and the migrated ice/soil solvers all call it. Column loops run under `@jit_compile(parallel=True)` + `prange`.

**Tech Stack:** Python, numpy, numba (via `climt/_core/backend.py`'s `jit_compile`/`prange`), sympl (`Stepper`, `get_constant`, `initialize_numpy_arrays_with_properties`), pytest, quarto docs.

## Global Constraints

- No new `scipy` in any diffusion/compute path. `scipy.sparse`/`spsolve` must be gone from `snow_ice_column.py` and `second_best/processes/subsurface.py`. (Non-solver scipy — `CubicSpline`, `RegularGridInterpolator`, `cKDTree`, `netcdf_file` — is out of scope and stays.)
- All new numerics use `from ..._core.backend import jit_compile, prange` — never `import numba` directly, never `@jit(nopython=True, parallel=True)`. Parallel kernels use `@jit_compile(parallel=True)` (the decorator forwards `parallel=True` to `njit`, and no-ops when numba is absent).
- Solver functions are pure: they allocate and return a new array; they never mutate their inputs (keeps a future JAX/`lax.scan` swap mechanical).
- `data_ocean/_sst_interpolation.py` is **not** touched (it uses `np.linalg.solve`, cold path, no scipy).
- The `Kb` diffusivity multiplier uses the **surface-layer `Ri_a`**, not the local `Ri` (thesis Eqn 2.8). See spec `docs/superpowers/specs/2026-07-22-boundary-layer-and-jit-tridiagonal-design.md`.
- Component/dir naming: package `climt/_components/simple_boundary_layer/`, class `SimpleBoundaryLayer`.
- Commit messages end with the `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>` trailer.

## File Structure

- Create `climt/_core/tridiagonal.py` — `solve_tridiagonal`, `tridiagonal_matvec`.
- Modify `climt/_core/__init__.py` — export the two functions.
- Create `tests/test_tridiagonal.py` — solver unit tests.
- Create `climt/_components/simple_boundary_layer/component.py` — kernel + `SimpleBoundaryLayer`.
- Create `climt/_components/simple_boundary_layer/__init__.py` — re-export.
- Modify `climt/_components/__init__.py` — import + `__all__`.
- Modify `climt/__init__.py` — import + `__all__`.
- Modify `tests/test_components.py` — `TestSimpleBoundaryLayer` + I/O-property + conservation tests.
- Modify `climt/_core/snow_ice_column.py` — jit rewrite, scipy removed.
- Modify `climt/_components/sea_ice/component.py` — `prange` column kernel.
- Modify `climt/_components/land_ice/component.py` — `prange` column kernel.
- Modify `climt/_components/second_best/processes/subsurface.py` — jit rewrite, scipy removed.
- Create `docs/user-guide/simple_boundary_layer.qmd` — manual.
- Modify `docs/_quarto.yml` — user-guide entry only.
- Modify `references.bib` — Frierson 2006 Part I.

---

### Task 1: Canonical tridiagonal solver + unit tests

**Files:**
- Create: `climt/_core/tridiagonal.py`
- Modify: `climt/_core/__init__.py`
- Test: `tests/test_tridiagonal.py`

**Interfaces:**
- Produces:
  - `solve_tridiagonal(lower, diag, upper, rhs) -> np.ndarray` — Thomas solve. `lower`/`upper` length `n-1`, `diag`/`rhs` length `n`. Returns length-`n` solution. Pure.
  - `tridiagonal_matvec(lower, diag, upper, x) -> np.ndarray` — returns `A @ x` for the same banded layout. Pure.

- [ ] **Step 1: Write the failing test** — `tests/test_tridiagonal.py`

```python
import numpy as np
import pytest

from climt._core.tridiagonal import solve_tridiagonal, tridiagonal_matvec


def _dense(lower, diag, upper):
    n = diag.shape[0]
    A = np.diag(diag).astype(float)
    for i in range(n - 1):
        A[i + 1, i] = lower[i]
        A[i, i + 1] = upper[i]
    return A


@pytest.mark.parametrize("n", [1, 2, 3, 5, 30])
def test_solve_matches_numpy(n):
    rng = np.random.RandomState(0)
    diag = rng.rand(n) + 3.0            # diagonally dominant => well-conditioned
    lower = rng.rand(max(n - 1, 0))
    upper = rng.rand(max(n - 1, 0))
    rhs = rng.rand(n)
    x = solve_tridiagonal(lower, diag, upper, rhs)
    if n == 1:
        expected = rhs / diag
    else:
        expected = np.linalg.solve(_dense(lower, diag, upper), rhs)
    assert np.allclose(x, expected, atol=1e-10)


def test_solve_does_not_mutate_inputs():
    lower = np.array([1.0, 1.0])
    diag = np.array([4.0, 4.0, 4.0])
    upper = np.array([1.0, 1.0])
    rhs = np.array([1.0, 2.0, 3.0])
    lower0, diag0, upper0, rhs0 = (a.copy() for a in (lower, diag, upper, rhs))
    solve_tridiagonal(lower, diag, upper, rhs)
    for a, a0 in ((lower, lower0), (diag, diag0), (upper, upper0), (rhs, rhs0)):
        assert np.array_equal(a, a0)


def test_matvec_matches_dense():
    rng = np.random.RandomState(1)
    n = 6
    diag = rng.rand(n)
    lower = rng.rand(n - 1)
    upper = rng.rand(n - 1)
    x = rng.rand(n)
    y = tridiagonal_matvec(lower, diag, upper, x)
    assert np.allclose(y, _dense(lower, diag, upper) @ x)


def test_solve_matvec_roundtrip():
    rng = np.random.RandomState(2)
    n = 10
    diag = rng.rand(n) + 3.0
    lower = rng.rand(n - 1)
    upper = rng.rand(n - 1)
    x = rng.rand(n)
    rhs = tridiagonal_matvec(lower, diag, upper, x)
    assert np.allclose(solve_tridiagonal(lower, diag, upper, rhs), x, atol=1e-10)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_tridiagonal.py -q`
Expected: FAIL / collection error — `climt._core.tridiagonal` does not exist.

- [ ] **Step 3: Write the implementation** — `climt/_core/tridiagonal.py`

```python
# -*- coding: utf-8 -*-
"""Canonical dense tridiagonal solver for climt.

This is THE dense tridiagonal solver for the library. Diffusion-type
components (boundary layer, snow/ice/soil columns, subsurface transport)
should call ``solve_tridiagonal`` rather than hand-rolling a Thomas sweep or
pulling in ``scipy.sparse``/``spsolve`` (which cannot be numba-jitted or, in
future, JAX-traced).

Tridiagonal families in climt, for reference:
  * dense, non-cyclic  -> ``solve_tridiagonal`` (this module) -- default.
  * cyclic/periodic    -> ``data_ocean._sst_interpolation`` (12x12
    ``np.linalg.solve``, cold path); needs a Sherman-Morrison variant if it
    ever moves onto a jit path.

JAX note: the forward/backward sweeps are loop-carried. A JAX backend would
express them via ``jax.lax.scan`` behind the ``jit_compile(backend=...)``
hook in ``backend.py``. Both functions here are pure (they allocate and
return; inputs are never mutated), which keeps that swap mechanical.
"""
import numpy as np

from .backend import jit_compile


@jit_compile
def solve_tridiagonal(lower, diag, upper, rhs):
    """Solve ``A x = rhs`` for a tridiagonal ``A`` via the Thomas algorithm.

    ``lower`` (sub-diagonal, ``A[i, i-1]``) and ``upper`` (super-diagonal,
    ``A[i, i+1]``) have length ``n-1``; ``diag`` and ``rhs`` have length ``n``.
    Returns the length-``n`` solution. Inputs are not mutated.
    """
    n = rhs.shape[0]
    x = np.zeros(n)
    if n == 1:
        x[0] = rhs[0] / diag[0]
        return x

    c_prime = np.zeros(n - 1)
    d_prime = np.zeros(n)
    c_prime[0] = upper[0] / diag[0]
    d_prime[0] = rhs[0] / diag[0]
    for i in range(1, n - 1):
        m = diag[i] - lower[i - 1] * c_prime[i - 1]
        c_prime[i] = upper[i] / m
        d_prime[i] = (rhs[i] - lower[i - 1] * d_prime[i - 1]) / m
    m = diag[n - 1] - lower[n - 2] * c_prime[n - 2]
    d_prime[n - 1] = (rhs[n - 1] - lower[n - 2] * d_prime[n - 2]) / m

    x[n - 1] = d_prime[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x


@jit_compile
def tridiagonal_matvec(lower, diag, upper, x):
    """Return ``A @ x`` for the banded layout used by ``solve_tridiagonal``.

    Same argument layout as ``solve_tridiagonal``. Pure.
    """
    n = x.shape[0]
    y = diag * x
    if n > 1:
        y[1:] = y[1:] + lower * x[:-1]
        y[:-1] = y[:-1] + upper * x[1:]
    return y
```

- [ ] **Step 4: Export from `climt/_core/__init__.py`**

Add near the other `_core` re-exports (append import line; keep alphabetical-ish grouping consistent with the file):

```python
from .tridiagonal import solve_tridiagonal, tridiagonal_matvec  # noqa: F401
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/test_tridiagonal.py -q`
Expected: PASS (all parametrizations).

- [ ] **Step 6: Commit**

```bash
git add climt/_core/tridiagonal.py climt/_core/__init__.py tests/test_tridiagonal.py
git commit -m "$(printf 'feat(core): add canonical jit tridiagonal solver\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 2: `SimpleBoundaryLayer` component + registration

**Files:**
- Create: `climt/_components/simple_boundary_layer/component.py`
- Create: `climt/_components/simple_boundary_layer/__init__.py`
- Modify: `climt/_components/__init__.py`
- Modify: `climt/__init__.py`

**Interfaces:**
- Consumes: `climt._core.tridiagonal.solve_tridiagonal` (Task 1); `climt._core.backend.jit_compile, prange`.
- Produces: `SimpleBoundaryLayer(Stepper)` with inputs `air_temperature`, `specific_humidity`, `air_pressure`, `air_pressure_on_interface_levels`, `northward_wind`, `eastward_wind`, `surface_air_pressure`, `surface_temperature`, `surface_specific_humidity`; outputs `air_temperature`, `specific_humidity`, `northward_wind`, `eastward_wind`; diagnostics `northward_wind_stress`, `eastward_wind_stress`, `boundary_layer_height`. Kwargs `von_karman_constant=0.4`, `roughness_length=3.21e-5`, `specific_fraction=0.1`, `reference_pressure=100000`, `critical_richardson_number=1`.

- [ ] **Step 1: Write the component** — `climt/_components/simple_boundary_layer/component.py`

```python
from typing import NamedTuple

import numpy as np
from sympl import (
    Stepper,
    get_constant,
    initialize_numpy_arrays_with_properties,
)

from ..._core.backend import jit_compile, prange
from ..._core.tridiagonal import solve_tridiagonal


class BoundaryLayerParams(NamedTuple):
    Rd: float
    Cp: float
    g: float
    k: float
    z0: float
    fb: float
    P0: float
    Ric: float


@jit_compile
def _richardson_diffusivity(Ri_a, u_fric, C, z, params):
    """Surface-layer diffusion coefficient K_b (thesis Eqn 2.8).

    Uses the *surface-layer* Richardson number ``Ri_a`` in the multiplier
    (NOT the local Ri of canonical Frierson 2006 Eqn 2.6), so K_b is
    continuous at ``Ri_a == 0`` and decays 1 -> 0 as ``Ri_a -> Ric``.
    """
    base = params.k * u_fric * np.sqrt(C) * z
    if Ri_a <= 0:
        return base
    return base / (
        1.0
        + Ri_a / params.Ric * np.log(z / params.z0) / (1.0 - Ri_a / params.Ric)
    )


@jit_compile
def _diffuse_profile(profile, p, p_int, rho, diff, dt, g):
    """Implicit vertical diffusion of ``profile`` with no-flux boundaries.

    Assembles the tridiagonal system and solves it with the shared Thomas
    solver. ``rho``/``diff`` have length ``num_layers - 1``.
    """
    num_layers = profile.shape[0]
    diag_m = np.zeros(num_layers)
    diag_p = np.zeros(num_layers)
    diag = np.zeros(num_layers)
    for i in range(num_layers):
        if i != 0:
            diag_m[i] = (
                g * g * rho[i - 1] * rho[i - 1] * diff[i - 1] * dt
                / (p[i - 1] - p[i])
                / (p_int[i] - p_int[i + 1])
            )
        if i != num_layers - 1:
            diag_p[i] = (
                g * g * rho[i] * rho[i] * diff[i] * dt
                / (p[i] - p[i + 1])
                / (p_int[i] - p_int[i + 1])
            )
        diag[i] = 1.0 + diag_m[i] + diag_p[i]
    return solve_tridiagonal(-diag_m[1:], diag, -diag_p[:-1], profile)


@jit_compile(parallel=True)
def _boundary_layer_kernel(
    air_temperature, surface_temperature, air_pressure, air_pressure_int,
    surface_pressure, specific_humidity, surface_humidity, northward_wind,
    eastward_wind, new_air_temperature, new_specific_humidity,
    new_northward_wind, new_eastward_wind, north_wind_stress, east_wind_stress,
    boundary_height, params, num_cols, timestep,
):
    Rd = params.Rd
    Cp = params.Cp
    g = params.g
    k = params.k
    z0 = params.z0
    fb = params.fb
    P0 = params.P0
    Ri_c = params.Ric

    for col in prange(num_cols):
        col_T = air_temperature[:, col]
        col_Ts = surface_temperature[col]
        col_p = air_pressure[:, col]
        col_p_int = air_pressure_int[:, col]
        col_ps = surface_pressure[col]
        col_q = specific_humidity[:, col]
        col_qs = surface_humidity[col]
        col_v = northward_wind[:, col]
        col_u = eastward_wind[:, col]

        col_v_int = 0.5 * (col_v[1:] + col_v[:-1])
        col_u_int = 0.5 * (col_u[1:] + col_u[:-1])
        col_T_int = 0.5 * (col_T[1:] + col_T[:-1])
        col_q_int = 0.5 * (col_q[1:] + col_q[:-1])
        col_rho = col_p_int[1:-1] / (
            Rd * (1.0 + 0.608 * col_q_int) * col_T_int
        )

        n = col_T_int.shape[0]

        pot_virt_temp = (
            col_T_int
            * np.power(P0 / col_p_int[1:-1], Rd / Cp)
            * (1.0 + 0.608 * col_q_int)
        )
        pot_virt_temp_surf = (
            col_Ts * np.power(P0 / col_ps, Rd / Cp) * (1.0 + 0.608 * col_qs)
        )

        z = np.zeros(n)
        z[0] = (
            Rd * (1.0 + 0.608 * col_q[0]) * col_T[0] / g
        ) * np.log(col_ps / col_p_int[1:-1][0])
        for i in range(1, n):
            z[i] = z[i - 1] + (
                Rd * (1.0 + 0.608 * col_q[i]) * col_T[i] / g
            ) * np.log(col_p_int[1:-1][i - 1] / col_p_int[1:-1][i])

        wind_int = np.sqrt(col_v_int * col_v_int + col_u_int * col_u_int)
        for i in range(wind_int.shape[0]):
            if wind_int[i] < 1.0:
                wind_int[i] = 1.0

        Ri_a = (
            g * z[0] * (pot_virt_temp[0] - pot_virt_temp_surf)
            / (pot_virt_temp_surf * wind_int[0] * wind_int[0])
        )
        if Ri_a < 0:
            C = k * k * np.power(np.log(z[0] / z0), -2)
        elif Ri_a < Ri_c:
            C = k * k * np.power(np.log(z[0] / z0), -2) * np.power(
                1.0 - Ri_a / Ri_c, 2
            )
        else:
            C = 0.0

        count = 0
        Rich = np.zeros(n)
        for i in range(n):
            Rich[i] = (
                g * z[i] * (pot_virt_temp[i] - pot_virt_temp[0])
                / (pot_virt_temp[0] * wind_int[i] * wind_int[i])
            )
            if Rich[i] > Ri_c:
                count = i + 1
                break
        h = z[count - 1]
        boundary_height[col] = h

        north_wind_stress[col] = col_rho[0] * C * wind_int[0] * col_v_int[0]
        east_wind_stress[col] = col_rho[0] * C * wind_int[0] * col_u_int[0]

        u_fric = wind_int[0]

        diff = np.zeros(n)
        for i in range(count):
            if z[i] < fb * h:
                diff[i] = _richardson_diffusivity(Ri_a, u_fric, C, z[i], params)
            else:
                diff[i] = (
                    _richardson_diffusivity(Ri_a, u_fric, C, fb * h, params)
                    * z[i] / (h * fb)
                    * np.power(1.0 - (z[i] - fb * h) / ((1.0 - fb) * h), 2)
                )

        new_air_temperature[:, col] = _diffuse_profile(
            col_T, col_p, col_p_int, col_rho, diff, timestep, g
        )
        new_specific_humidity[:, col] = _diffuse_profile(
            col_q, col_p, col_p_int, col_rho, diff, timestep, g
        )
        new_northward_wind[:, col] = _diffuse_profile(
            col_v, col_p, col_p_int, col_rho, diff, timestep, g
        )
        new_eastward_wind[:, col] = _diffuse_profile(
            col_u, col_p, col_p_int, col_rho, diff, timestep, g
        )


class SimpleBoundaryLayer(Stepper):
    """A simple boundary-layer scheme that diffuses heat, moisture and
    momentum upward from the lowest model level.

    This is the surface-flux / boundary-layer formulation of
    Frierson, Held & Zurita-Gotor (2006), with the surface-layer diffusion
    coefficient modified (thesis Eqn 2.8) to use the surface-layer Richardson
    number in its multiplier, making it continuous at ``Ri_a == 0``.

    The component assumes a surface-flux component has already applied the
    surface fluxes at the lowest model level; it then diffuses the resulting
    profiles using diffusion coefficients from a simplified Monin-Obukhov
    theory with a K-profile capped by a critical Richardson number.
    """

    input_properties = {
        'air_temperature': {'dims': ['mid_levels', '*'], 'units': 'degK'},
        'specific_humidity': {'dims': ['mid_levels', '*'], 'units': 'kg/kg'},
        'air_pressure': {'dims': ['mid_levels', '*'], 'units': 'Pa'},
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'], 'units': 'Pa',
        },
        'northward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
        'eastward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
        'surface_air_pressure': {'dims': ['*'], 'units': 'Pa'},
        'surface_temperature': {'dims': ['*'], 'units': 'degK'},
        'surface_specific_humidity': {'dims': ['*'], 'units': 'kg/kg'},
    }

    output_properties = {
        'air_temperature': {'dims': ['mid_levels', '*'], 'units': 'degK'},
        'specific_humidity': {'dims': ['mid_levels', '*'], 'units': 'kg/kg'},
        'northward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
        'eastward_wind': {'dims': ['mid_levels', '*'], 'units': 'm s^-1'},
    }

    diagnostic_properties = {
        'northward_wind_stress': {'dims': ['*'], 'units': 'Pa'},
        'eastward_wind_stress': {'dims': ['*'], 'units': 'Pa'},
        'boundary_layer_height': {'dims': ['*'], 'units': 'm'},
    }

    def __init__(self, von_karman_constant=0.4, roughness_length=0.0000321,
                 specific_fraction=0.1, reference_pressure=100000,
                 critical_richardson_number=1, **kwargs):
        """
        Args:
            von_karman_constant: von Karman constant k.
            roughness_length: surface roughness length z0 (m).
            specific_fraction: surface-layer fraction fb of the boundary
                layer depth.
            reference_pressure: reference pressure P0 (Pa) for potential
                temperature.
            critical_richardson_number: critical Richardson number Ric that
                caps the diffusion and sets the boundary-layer top.
        """
        self._k = von_karman_constant
        self._z0 = roughness_length
        self._fb = specific_fraction
        self._P0 = reference_pressure
        self._Ric = critical_richardson_number
        self._update_constants()
        super(SimpleBoundaryLayer, self).__init__(**kwargs)

    def _update_constants(self):
        self._Rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        self._Cp = get_constant(
            'heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1'
        )
        self._g = get_constant('gravitational_acceleration', 'm s^-2')

    def array_call(self, state, timestep):
        """Diffuse temperature, humidity and wind profiles for each column."""
        num_cols = state['air_temperature'].shape[1]

        new_state = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties
        )
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        params = BoundaryLayerParams(
            Rd=self._Rd, Cp=self._Cp, g=self._g, k=self._k, z0=self._z0,
            fb=self._fb, P0=self._P0, Ric=self._Ric,
        )

        _boundary_layer_kernel(
            state['air_temperature'],
            state['surface_temperature'],
            state['air_pressure'],
            state['air_pressure_on_interface_levels'],
            state['surface_air_pressure'],
            state['specific_humidity'],
            state['surface_specific_humidity'],
            state['northward_wind'],
            state['eastward_wind'],
            new_state['air_temperature'],
            new_state['specific_humidity'],
            new_state['northward_wind'],
            new_state['eastward_wind'],
            diagnostics['northward_wind_stress'],
            diagnostics['eastward_wind_stress'],
            diagnostics['boundary_layer_height'],
            params,
            num_cols,
            timestep.total_seconds(),
        )

        return diagnostics, new_state
```

- [ ] **Step 2: Write** `climt/_components/simple_boundary_layer/__init__.py`

```python
from .component import SimpleBoundaryLayer

__all__ = SimpleBoundaryLayer
```

- [ ] **Step 3: Register in `climt/_components/__init__.py`**

Add the import (alphabetical, after `.sea_ice`/before `.simple_physics` grouping is fine — match the existing style):

```python
from .simple_boundary_layer import SimpleBoundaryLayer
```

Add `SimpleBoundaryLayer,` into the `__all__` tuple.

- [ ] **Step 4: Register in `climt/__init__.py`**

Mirror the same two edits (import block + `__all__` tuple) in `climt/__init__.py`, matching how `DryConvectiveAdjustment`/`BucketHydrology` appear there.

- [ ] **Step 5: Smoke-test import and a single call**

Run:
```bash
python -c "
import numpy as np, climt
from climt import get_grid, get_default_state
c = climt.SimpleBoundaryLayer()
state = get_default_state([c], grid_state=get_grid(nx=None, ny=None, nz=30))
from datetime import timedelta
diag, new = c(state, timestep=timedelta(seconds=600))
print('diag keys:', sorted(diag.keys()))
print('out keys:', sorted(new.keys()))
print('any NaN T:', bool(np.any(np.isnan(np.asarray(new['air_temperature'])))))
"
```
Expected: prints the three diagnostic keys, four output keys, and `any NaN T: False`.

- [ ] **Step 6: Commit**

```bash
git add climt/_components/simple_boundary_layer climt/_components/__init__.py climt/__init__.py
git commit -m "$(printf 'feat(components): add SimpleBoundaryLayer (Frierson 2006, thesis Eqn 2.8)\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 3: `SimpleBoundaryLayer` tests (standard suite + I/O props + conservation)

**Files:**
- Modify: `tests/test_components.py`

**Interfaces:**
- Consumes: `climt.SimpleBoundaryLayer`; existing `ComponentBaseColumn`, `ComponentBase3D`, `call_with_timestep_if_needed`, `get_grid`, `climt.get_default_state`.

- [ ] **Step 1: Add the standard-suite test class** — near the other `Test*` classes in `tests/test_components.py`

```python
class TestSimpleBoundaryLayer(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.SimpleBoundaryLayer()
```

Also add `SimpleBoundaryLayer` to the `from climt import (...)` block at the top of the file.

- [ ] **Step 2: Add the I/O-property test** — same file

```python
def test_simple_boundary_layer_properties():
    c = climt.SimpleBoundaryLayer()
    assert set(c.input_properties) == {
        'air_temperature', 'specific_humidity', 'air_pressure',
        'air_pressure_on_interface_levels', 'northward_wind', 'eastward_wind',
        'surface_air_pressure', 'surface_temperature',
        'surface_specific_humidity',
    }
    assert set(c.output_properties) == {
        'air_temperature', 'specific_humidity', 'northward_wind',
        'eastward_wind',
    }
    assert set(c.diagnostic_properties) == {
        'northward_wind_stress', 'eastward_wind_stress', 'boundary_layer_height',
    }
    assert c.diagnostic_properties['boundary_layer_height']['units'] == 'm'
    assert c.output_properties['air_temperature']['units'] == 'degK'
```

- [ ] **Step 3: Add the conservation test** — same file

The implicit no-flux scheme conserves the interface-pressure-weighted column integral of each diffused field. Weight per layer is `p_int[k] - p_int[k+1]`.

```python
def test_simple_boundary_layer_conserves_column_integrals():
    import numpy as np
    from datetime import timedelta

    component = climt.SimpleBoundaryLayer()
    state = climt.get_default_state(
        [component], grid_state=get_grid(nx=None, ny=None, nz=30)
    )
    # perturb winds/humidity so there is something to diffuse
    np.asarray(state['northward_wind'])[:] = 5.0
    np.asarray(state['eastward_wind'])[:] = -3.0
    np.asarray(state['specific_humidity'])[:] = 0.005

    diagnostics, new_state = component(state, timestep=timedelta(seconds=600))

    p_int = np.asarray(state['air_pressure_on_interface_levels'])
    # mid-levels axis is axis 0 for a single column here
    dp = p_int[:-1] - p_int[1:]
    for name in ('air_temperature', 'specific_humidity',
                 'northward_wind', 'eastward_wind'):
        before = np.sum(np.asarray(state[name]) * dp)
        after = np.sum(np.asarray(new_state[name]) * dp)
        assert np.isclose(before, after, rtol=1e-8, atol=1e-6), name
```

- [ ] **Step 4: Run the boundary-layer tests**

Run: `python -m pytest tests/test_components.py -k SimpleBoundaryLayer -q`
Expected: the property and conservation tests PASS. The `ComponentBase` cached-output tests will FAIL on the first run with "Failed due to no cached output, cached current output." — this is by design; the run writes new cache files under `tests/cached_component_output/`.

- [ ] **Step 5: Re-run to confirm cached-output tests pass**

Run: `python -m pytest tests/test_components.py -k SimpleBoundaryLayer -q`
Expected: all PASS now that caches exist.

- [ ] **Step 6: Commit** (include the new cache files)

```bash
git add tests/test_components.py tests/cached_component_output/TestSimpleBoundaryLayer-*.cache
git commit -m "$(printf 'test(components): SimpleBoundaryLayer suite, I/O-property and conservation tests\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 4: Remove scipy from `snow_ice_column` (+ `prange` for `SeaIce`/`LandIce`)

**Files:**
- Modify: `climt/_core/snow_ice_column.py`
- Modify: `climt/_components/sea_ice/component.py`
- Modify: `climt/_components/land_ice/component.py`

**Interfaces:**
- Consumes: `solve_tridiagonal`, `tridiagonal_matvec` (Task 1); `jit_compile`, `prange`.
- Produces: unchanged public `solve_column(rho, c, kappa, temperature, dt, dz, top_bc, bottom_bc)` signature (still accepts `Dirichlet`/`Flux`); plus a `@jit_compile` core the components' `prange` loops call.

**Background (read before editing):** `solve_column` currently builds two `scipy.sparse` matrices, does `mat_rhs * temperature`, folds BC rows into the diagonals using `isinstance(bc, Dirichlet/Flux)`, and calls `spsolve`. BC dataclasses + `isinstance` cannot run under `njit`, so the numeric core must take **integer BC flags + scalar values**.

- [ ] **Step 1: Add a jitable core to `climt/_core/snow_ice_column.py`**

Keep the `Dirichlet`/`Flux` dataclasses. Add BC-flag constants and a `@jit_compile` kernel that reproduces the current assembly exactly, replacing the sparse RHS matvec with `tridiagonal_matvec` and `spsolve` with `solve_tridiagonal`:

```python
import numpy as np

from .backend import jit_compile
from .tridiagonal import solve_tridiagonal, tridiagonal_matvec

_BC_DIRICHLET = 0
_BC_FLUX = 1


@jit_compile
def _solve_column_kernel(rho, c, kappa, temperature, dt, dz,
                         top_flag, top_val, bot_flag, bot_val):
    num_layers = temperature.shape[0]
    K_mid = kappa
    K_interface = 0.5 * (kappa[:-1] + kappa[1:])
    heat_capacity = rho * c
    heat_capacity_int = 0.5 * (heat_capacity[:-1] + heat_capacity[1:])
    mu_inv_int = dt / (heat_capacity_int * 2.0 * dz * dz)

    r = np.zeros(num_layers)
    a_sub = np.zeros(num_layers)   # length-num_layers scratch, matches original
    a_sup = np.zeros(num_layers)
    r[1:-1] = K_interface * mu_inv_int
    dp = 1.0 + 2.0 * r
    dm = 1.0 - 2.0 * r
    a_sub[:-2] = -mu_inv_int * K_mid[:-1]
    a_sup[2:] = -mu_inv_int * K_mid[1:]

    # RHS = mat_rhs @ temperature, where mat_rhs has diagonal dm and
    # off-diagonals -a_sub / -a_sup (banded layout: sub len n-1, sup len n-1).
    rhs = tridiagonal_matvec(-a_sub[1:], dm, -a_sup[:-1], temperature)

    dp_lhs = dp.copy()
    a_sub_lhs = a_sub.copy()
    a_sup_lhs = a_sup.copy()

    # Top boundary (node n-1)
    if top_flag == _BC_DIRICHLET:
        dp_lhs[-1] = 1.0
        a_sub_lhs[-2] = 0.0
        rhs[-1] = top_val
    else:
        dp_lhs[-1] = -1.0
        a_sub_lhs[-2] = 1.0
        rhs[-1] = -top_val * dz / K_mid[-1]

    # Bottom boundary (node 0)
    if bot_flag == _BC_DIRICHLET:
        dp_lhs[0] = 1.0
        a_sup_lhs[1] = 0.0
        rhs[0] = bot_val
    else:
        dp_lhs[0] = -1.0
        a_sup_lhs[1] = 1.0
        rhs[0] = -bot_val * dz / K_mid[0]

    return solve_tridiagonal(a_sub_lhs[1:], dp_lhs, a_sup_lhs[:-1], rhs)
```

> Note the banded-index mapping: the original `a_sub`/`a_sup` are length-`num_layers` arrays aligned so that `spdiags([a_sub, dp, a_sup], [-1, 0, 1])` places `a_sub[1:]` on the sub-diagonal and `a_sup[:-1]` on the super-diagonal. `solve_tridiagonal`/`tridiagonal_matvec` take those `n-1`-length slices directly.

- [ ] **Step 2: Rewrite the public `solve_column` as a thin translator**

Replace the scipy body of `solve_column` with a BC-to-flag translation that calls the kernel (drop the `from scipy...` imports):

```python
def solve_column(rho, c, kappa, temperature, dt, dz, top_bc, bottom_bc):
    top_flag = _BC_DIRICHLET if isinstance(top_bc, Dirichlet) else _BC_FLUX
    bot_flag = _BC_DIRICHLET if isinstance(bottom_bc, Dirichlet) else _BC_FLUX
    return _solve_column_kernel(
        np.asarray(rho, dtype=float), np.asarray(c, dtype=float),
        np.asarray(kappa, dtype=float), np.asarray(temperature, dtype=float),
        float(dt), float(dz),
        top_flag, float(top_bc.value), bot_flag, float(bottom_bc.value),
    )
```

Delete `from scipy import sparse` and `from scipy.sparse.linalg import spsolve`.

- [ ] **Step 3: Verify the solver change in isolation**

Run: `python -m pytest tests/test_components.py -k "SeaIce or LandIce or IceSheet" -q`
Expected: PASS. (`solve_column`'s result is mathematically identical — Thomas vs `spsolve` on the same matrix — so cached outputs should match. If any test reports "no cached output" it means a cache is missing, not a failure; if values drifted at round-off, regenerate the cache per Step 5.)

- [ ] **Step 4: Move the `SeaIce` and `LandIce` column loops onto `prange`**

Read the full per-column loop in each `component.py` (`sea_ice`: around lines 240–360; `land_ice`: analogous). Extract the per-column numeric body — profile construction, the `solve_column` call(s), melt checks, clamps, and flux accounting — into a module-level `@jit_compile(parallel=True)` kernel that:
  - takes the stacked per-column input arrays + scalar params,
  - loops `for col in prange(ncol)`,
  - calls `_solve_column_kernel` directly (translate BCs to flags/values *outside* the kernel or compute them inline as flags/scalars),
  - writes results into preallocated output arrays,
  - returns any per-column quantities currently produced by logging as explicit output arrays instead of `logger` calls (numba cannot log). Keep any Python-level `logger.debug` summary outside the loop.

Preserve the existing physics and sign conventions exactly (see the extensive comments already in `sea_ice/component.py` about the basal-flux sign). This step is a mechanical extraction: the arithmetic per column is unchanged; only the loop is jitted and parallelised.

- [ ] **Step 5: Run the ice suites; regenerate caches only if needed**

Run: `python -m pytest tests/test_components.py -k "SeaIce or LandIce" -q`
Expected: PASS. If round-off drift makes a cached-output test fail, delete the affected `tests/cached_component_output/TestSeaIce-*.cache` / `TestLandIce-*.cache`, re-run to regenerate, eyeball that temperatures/thicknesses are physically sane (no NaNs, thickness ≥ 0), and commit the new caches. These components are new on this branch, so regenerating is acceptable.

- [ ] **Step 6: Confirm scipy is gone from this path**

Run: `grep -n "scipy" climt/_core/snow_ice_column.py climt/_components/sea_ice/component.py climt/_components/land_ice/component.py`
Expected: no matches.

- [ ] **Step 7: Commit**

```bash
git add climt/_core/snow_ice_column.py climt/_components/sea_ice/component.py climt/_components/land_ice/component.py tests/cached_component_output/TestSeaIce-*.cache tests/cached_component_output/TestLandIce-*.cache
git commit -m "$(printf 'refactor(ice): drop scipy from snow_ice_column, prange column loops\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 5: Remove scipy from `second_best` subsurface (+ `prange`)

**Files:**
- Modify: `climt/_components/second_best/processes/subsurface.py`
- Modify: `climt/_components/second_best/component.py` (column loop → `prange`)

**Interfaces:**
- Consumes: `solve_tridiagonal` (Task 1); `jit_compile`, `prange`.
- Produces: `BestSubsurfaceTransport.__call__` unchanged in behaviour; a `@jit_compile` numeric core replacing the `scipy.sparse` assembly.

**Background:** `BestSubsurfaceTransport.__call__` builds a `scipy.sparse` LIL matrix with Neumann top/bottom rows, calls `spsolve`, then applies an explicit freeze/melt source. Preserve `second_best`'s independence — call the shared `solve_tridiagonal`, do **not** import the ice core.

- [ ] **Step 1: Replace the sparse assembly + `spsolve` with `solve_tridiagonal`**

In `subsurface.py`, drop `from scipy import sparse` / `from scipy.sparse.linalg import spsolve`. Build the three diagonal arrays explicitly (same values that went into `sparse.diags`/the Neumann rows) and solve:

```python
import numpy as np

from ...._core.backend import jit_compile, prange
from ...._core.tridiagonal import solve_tridiagonal
```

Then, inside `__call__`, where the matrix was assembled (`off`, `main`, the `A[0,0]`, `A[-1,-1]` Neumann rows), construct:

```python
        # off has length n-1 (both sub and super are the same off array today)
        lower = -off[1:].copy() if False else None  # replaced below
```

(Read the current exact `off`/`main` construction; the LIL matrix uses `[off[1:], main, off[:-1]]` on diagonals `[-1, 0, 1]`, with Neumann overrides `A[0,0]=1+rr; A[0,1]=-rr` and `A[-1,-1]=1+rr; A[-1,-2]=-rr`.) Produce:

```python
        lower = off[1:].copy()          # sub-diagonal, length n-1
        upper = off[:-1].copy()         # super-diagonal, length n-1
        main = main.copy()
        # Neumann bottom (node 0)
        main[0] = 1 + rr; upper[0] = -rr
        # Neumann surface (node n-1)
        main[-1] = 1 + rr; lower[-1] = -rr
        rhs = T.copy()
        rhs[-1] = T[-1] + surface_flux_bc * dt / (self._cv * dz)
        T_diff = solve_tridiagonal(lower, main, upper, rhs)
```

Keep the freeze/melt `Gamma` source arithmetic that follows exactly as-is.

- [ ] **Step 2: Verify the subsurface change in isolation**

Run: `python -m pytest tests/test_components.py -k SecondBEST -q`
Expected: PASS (Thomas vs `spsolve` identical on the same system).

- [ ] **Step 3: Move the `second_best` column loop onto `prange`**

In `second_best/component.py` the loop is `for col in range(ncol):` (~line 93) calling `self._subsurface(...)`. Extract the per-column subsurface transport into a `@jit_compile(parallel=True)` kernel over `prange(ncol)` that calls `solve_tridiagonal`, writing into preallocated output arrays. If the surrounding per-column code mixes dict access and Python objects, first marshal the needed fields into plain 2-D arrays before the kernel, then scatter results back after. Preserve behaviour exactly.

- [ ] **Step 4: Run the SecondBEST suite; regenerate caches if needed**

Run: `python -m pytest tests/test_components.py -k SecondBEST -q`
Expected: PASS. Regenerate `TestSecondBEST-*.cache` only if round-off drift requires it (new component on this branch; acceptable).

- [ ] **Step 5: Confirm scipy is gone**

Run: `grep -n "scipy" climt/_components/second_best/processes/subsurface.py`
Expected: no matches.

- [ ] **Step 6: Commit**

```bash
git add climt/_components/second_best tests/cached_component_output/TestSecondBEST-*.cache
git commit -m "$(printf 'refactor(second_best): drop scipy from subsurface, prange column loop\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 6: Docs — manual, quarto entry, references.bib

**Files:**
- Create: `docs/user-guide/simple_boundary_layer.qmd`
- Modify: `docs/_quarto.yml`
- Modify: `references.bib`

**Interfaces:**
- Consumes: `climt.SimpleBoundaryLayer` (Task 2).

- [ ] **Step 1: Add the reference** — append to `references.bib`

```bibtex
@article{frierson2006gray,
  title   = {A Gray-Radiation Aquaplanet Moist {GCM}. {Part I}: Static Stability and Eddy Scale},
  author  = {Frierson, Dargan M. W. and Held, Isaac M. and Zurita-Gotor, Pablo},
  journal = {Journal of the Atmospheric Sciences},
  volume  = {63},
  number  = {10},
  pages   = {2548--2566},
  year    = {2006},
  doi     = {10.1175/JAS3753.1}
}
```

- [ ] **Step 2: Write the manual** — `docs/user-guide/simple_boundary_layer.qmd`

Model it on `docs/user-guide/dry_convective_adjustment.qmd` (open that first to match front-matter/format). Content must cover: what the scheme is (Frierson 2006 simplified Monin–Obukhov surface exchange + K-profile), the assumption that a surface-flux component ran first, the `Kb` thesis Eqn 2.8 modification (surface-layer `Ri_a`, continuity at `Ri_a = 0`), inputs/outputs/diagnostics, the five `__init__` parameters, a worked usage snippet, and the conserved quantity (mass-weighted column integrals). Skeleton:

````markdown
---
title: "Simple boundary layer"
---

`SimpleBoundaryLayer` diffuses heat, moisture and momentum upward from the
lowest model level using the surface-flux / boundary-layer formulation of
@frierson2006gray, with a modification to the surface-layer diffusion
coefficient described below.

## What it does

The component assumes a surface-flux component has already applied surface
fluxes at the lowest model level. It then computes a K-profile eddy
diffusivity from a simplified Monin–Obukhov theory, capped by a critical
Richardson number $Ri_c$, and applies an implicit vertical diffusion to
temperature, specific humidity and both horizontal wind components.

## Modification to the Frierson scheme

Canonical @frierson2006gray uses the *local* Richardson number in the
surface-layer diffusion coefficient $K_b$, which is discontinuous at
$Ri_a = 0$. Following the modification of the scheme, this component uses the
*surface-layer* Richardson number $Ri_a$ in the multiplier:

$$
K_b \propto \left[\,1 + \frac{Ri_a}{Ri_c}\,
\frac{\ln(z/z_0)}{1 - Ri_a/Ri_c}\right]^{-1}, \quad Ri_a > 0,
$$

so $K_b$ is continuous at $Ri_a = 0$ and decays from 1 to 0 as
$Ri_a \to Ri_c$.

## Inputs, outputs and diagnostics

| Kind | Quantities |
|------|------------|
| Inputs | air_temperature, specific_humidity, air_pressure, air_pressure_on_interface_levels, northward_wind, eastward_wind, surface_air_pressure, surface_temperature, surface_specific_humidity |
| Outputs | air_temperature, specific_humidity, northward_wind, eastward_wind |
| Diagnostics | northward_wind_stress, eastward_wind_stress, boundary_layer_height |

## Parameters

- `von_karman_constant` (0.4)
- `roughness_length` (3.21e-5 m)
- `specific_fraction` — surface-layer fraction $f_b$ (0.1)
- `reference_pressure` (100000 Pa)
- `critical_richardson_number` $Ri_c$ (1)

## Usage

```python
import climt
from datetime import timedelta

boundary_layer = climt.SimpleBoundaryLayer()
state = climt.get_default_state(
    [boundary_layer],
    grid_state=climt.get_grid(nx=None, ny=None, nz=30),
)
diagnostics, new_state = boundary_layer(state, timestep=timedelta(minutes=10))
```

## Conservation

The implicit no-flux diffusion conserves the interface-pressure-weighted
column integral of each diffused field (surface exchange is handled by the
separate surface-flux component).
````

- [ ] **Step 3: Register in `docs/_quarto.yml` (user-guide section only)**

Find the user-guide `Components` section (around the `dry_convective_adjustment.qmd` / `bucket_hydrology.qmd` entries near line 57–67) and add:

```yaml
            - user-guide/simple_boundary_layer.qmd
```

**Do NOT** add an entry to the API auto-listing section.

- [ ] **Step 4: Verify the manual renders / has no obvious errors**

Run (if quarto is available): `quarto render docs/user-guide/simple_boundary_layer.qmd --to html 2>&1 | tail -5`
If quarto is not installed, instead confirm the file is valid by checking the front-matter and that `@frierson2006gray` matches the bib key: `grep -n "frierson2006gray" references.bib docs/user-guide/simple_boundary_layer.qmd`
Expected: the bib key is present in both files.

- [ ] **Step 5: Commit**

```bash
git add docs/user-guide/simple_boundary_layer.qmd docs/_quarto.yml references.bib
git commit -m "$(printf 'docs: SimpleBoundaryLayer manual + user-guide entry\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 7: Full test sweep + graphify refresh

**Files:** none (verification + generated graph artifacts).

- [ ] **Step 1: Run the full component + solver test suite**

Run: `python -m pytest tests/test_tridiagonal.py tests/test_components.py -q`
Expected: PASS (all green). Investigate and fix any failure before continuing.

- [ ] **Step 2: Confirm no scipy remains in any diffusion/solver path**

Run: `grep -rn "spsolve\|scipy.sparse" climt/_core climt/_components/second_best climt/_components/sea_ice climt/_components/land_ice`
Expected: no matches.

- [ ] **Step 3: Refresh the knowledge graph** (per `CLAUDE.md`)

Run:
```bash
graphify update .
python scripts/augment_graph.py
```
Expected: both complete without error; `graphify-out/` updated.

- [ ] **Step 4: Commit the regenerated graph**

```bash
git add graphify-out
git commit -m "$(printf 'chore(graphify): regenerate graph for boundary layer + tridiagonal solver\n\nCo-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

## Self-Review notes

- **Spec coverage:** solver (T1), boundary layer + registration (T2), tests incl. conservation + I/O props (T3), snow_ice_column/sea_ice/land_ice scipy removal + prange (T4), second_best scipy removal + prange (T5), docs manual + quarto + bib, no API listing (T6), graphify refresh (T7), data_ocean excluded (Global Constraints). All spec sections map to a task.
- **Kb modification:** enforced in T2 `_richardson_diffusivity` (uses `Ri_a`) and documented in T6.
- **Type consistency:** `solve_tridiagonal(lower, diag, upper, rhs)` and `tridiagonal_matvec(lower, diag, upper, x)` used identically in T1, T2, T4, T5.
- **Known execution caveat:** T4/T5 Steps 4/3 (full `prange` extraction of the ice/second_best per-column bodies) require reading each component's current loop in full at execution time; the arithmetic is copied verbatim, only the loop is jitted/parallelised. If a clean full-column `njit` proves impractical for a branchy body, fall back to jitting the numeric core and keeping a thin Python loop (still scipy-free) — noted as a risk in the spec.
```
