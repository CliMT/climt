# SimpleBoundaryLayer surface fluxes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `SimpleBoundaryLayer` a surface boundary condition so heat, moisture and momentum actually enter the atmosphere — either from Frierson (2006) bulk formulae it computes itself, or from externally supplied surface heat/moisture fluxes.

**Architecture:** `_diffuse_profile` grows two scalar boundary-condition terms (`surface_exchange`, an implicit addition to `diag[0]`; `surface_source`, an explicit addition to `rhs[0]`). One code path serves all three modes; today's behaviour is both scalars zero. The kernel picks the scalars per column from an integer `flux_mode` (0 = none, 1 = bulk, 2 = external), reusing the `C`, `wind_int[0]` and `col_rho[0]` it already computes. The component exposes the mode as `surface_fluxes='bulk' | 'external' | None` and varies `input_properties` / `diagnostic_properties` per mode, following the instance-level-override pattern already used by `BucketHydrology(num_layers=2)`.

**Tech Stack:** Python, NumPy, numba (`@jit_compile` from `climt/_core/backend.py`), sympl (`Stepper`, `get_constant`, `initialize_numpy_arrays_with_properties`), pytest.

## Global Constraints

- Source spec: `docs/superpowers/specs/2026-08-07-sbl-surface-fluxes-design.md`.
- **Conda environment:** every Python/pytest command runs in the `climt` conda environment. Prefix with `conda run -n climt` or activate first.
- Everything inside `_boundary_layer_kernel` and `_diffuse_profile` must stay **numba-`njit`-compatible in `nopython` mode**: no dicts, no strings, no Python objects. The mode is an `int`, never a string. All branches must type-check for every call, so array arguments are always arrays (zeros when unused) — never `None`.
- Do not mutate input arrays. `profile` inside `_diffuse_profile` is a strided view into the caller's state array; copy before writing to it.
- `Δp₀ = p_int[0] − p_int[1]` (**not** `surface_air_pressure − p_int[1]`). `_diffuse_profile` already uses `p_int[i] − p_int[i+1]` as layer `i`'s mass, and exact budget closure requires the same value.
- Constants come from `sympl.get_constant`, never hard-coded numerals in the component:
  `gravitational_acceleration` (m s^-2), `heat_capacity_of_dry_air_at_constant_pressure` (J kg^-1 K^-1), `latent_heat_of_condensation` (J kg^-1, = 2.5e6).
- Surface flux sign convention: `surface_upward_sensible_heat_flux` and `surface_upward_latent_heat_flux` are W m^-2, **positive upward into the atmosphere**, so they enter `rhs[0]` with a **positive** sign.
- Default mode is `'bulk'`. This is a deliberate behaviour change, accepted in the spec.
- **Out of scope** (do not do these, even if tempting): switching heat diffusion from `T` to dry static energy; evaporative resistance / β for land; separate exchange coefficients for momentum vs heat vs moisture; changing `_richardson_diffusivity`'s use of `Ri_a`.
- Ignore `build/lib/climt/...` — it is a stale build copy, not source.

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `climt/_components/simple_boundary_layer/component.py` | `_diffuse_profile` BC terms, `BoundaryLayerParams.Lv`, kernel mode branch, `SimpleBoundaryLayer.__init__` API + per-mode properties, docstring | Modify |
| `tests/test_simple_boundary_layer.py` | New home for all surface-flux behaviour tests (`_diffuse_profile` unit tests, budget, sign, consistency, per-mode properties) | Create |
| `tests/test_components.py` | Existing conservation test pinned to `surface_fluxes=None`; existing properties test kept for the default mode | Modify (lines 781–834) |
| `docs/user-guide/simple_boundary_layer.qmd` | Mode-dependent framing, surface-flux physics section, updated tables | Modify |
| `examples/full_gcm_land_ocean_ice.py` | `surface_fluxes='external'` migration + the `Stepper`-in-dycore wrapping bug | Modify |

Tests for the surface-flux feature go in a new `tests/test_simple_boundary_layer.py` rather than growing `test_components.py`, matching the codebase's existing per-component test files (`test_bucket_hydrology.py`, `test_land_ice.py`, `test_sea_ice.py`).

---

## Task 1: `_diffuse_profile` surface boundary terms

Pure numerics, no component API change. After this task `_diffuse_profile` can apply a surface boundary condition, but nothing calls it with nonzero terms yet — the four call sites pass `0.0, 0.0`.

**Files:**
- Modify: `climt/_components/simple_boundary_layer/component.py:42-67` (`_diffuse_profile`) and `:173-184` (the four call sites)
- Create: `tests/test_simple_boundary_layer.py`

**Interfaces:**
- Consumes: `solve_tridiagonal(lower, diag, upper, rhs)` from `climt/_core/tridiagonal.py` — `lower`/`upper` length `n-1`, `diag`/`rhs` length `n`, returns length `n`, does not mutate inputs.
- Produces: `_diffuse_profile(profile, p, p_int, rho, diff, dt, g, surface_exchange, surface_source)` — nine **positional** arguments (no defaults; numba default-argument handling across a jitted call boundary is not worth the risk). Returns a new length-`num_layers` array.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_simple_boundary_layer.py`. The `pytest`, `climt`, `timedelta`
and `get_constant` imports are unused by this task's tests — Tasks 2 and 3 use
all of them, so write the import block once, now, rather than churning it.

```python
"""Surface-flux behaviour of climt.SimpleBoundaryLayer.

Covers the three ``surface_fluxes`` modes ('bulk', 'external', None), the
tridiagonal surface boundary condition that implements them, and the column
budgets they must close.
"""
from datetime import timedelta

import numpy as np
import pytest
from sympl import get_constant

import climt
from climt._components.simple_boundary_layer.component import _diffuse_profile


# ---------------------------------------------------------------- fixtures

def _test_column():
    """A small, well-conditioned 4-layer column for _diffuse_profile tests."""
    p_int = np.array([1.0e5, 9.0e4, 8.0e4, 7.0e4, 6.0e4])
    p = 0.5 * (p_int[:-1] + p_int[1:])
    rho = np.array([1.0, 0.9, 0.8])      # length num_layers - 1
    diff = np.array([10.0, 8.0, 6.0])    # length num_layers - 1
    return p, p_int, rho, diff


# -------------------------------------------- _diffuse_profile boundary terms

def test_diffuse_profile_zero_boundary_terms_reproduce_no_flux():
    """Both scalars zero must reproduce the pre-existing no-flux solve."""
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              0.0, 0.0)
    dp = p_int[:-1] - p_int[1:]
    assert np.isclose(np.sum(result * dp), np.sum(profile * dp), rtol=1e-12)


def test_diffuse_profile_does_not_mutate_profile():
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    before = profile.copy()
    _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665, 0.5, 140.0)
    assert np.array_equal(profile, before)


def test_diffuse_profile_bulk_term_matches_analytic_two_layer():
    """With diff == 0 the system is diagonal, so layer 0 has a closed form."""
    p, p_int, rho, _ = _test_column()
    diff = np.zeros(3)
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    beta = 0.25
    surface_value = 300.0
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              beta, beta * surface_value)
    expected0 = (profile[0] + beta * surface_value) / (1.0 + beta)
    assert np.isclose(result[0], expected0, rtol=1e-12)
    assert np.allclose(result[1:], profile[1:], rtol=1e-12)


def test_diffuse_profile_neumann_source_adds_exact_column_amount():
    """A pure Neumann source changes the mass-weighted integral by dp0 * S."""
    p, p_int, rho, diff = _test_column()
    profile = np.array([280.0, 275.0, 270.0, 265.0])
    source = 3.0
    result = _diffuse_profile(profile, p, p_int, rho, diff, 600.0, 9.80665,
                              0.0, source)
    dp = p_int[:-1] - p_int[1:]
    change = np.sum(result * dp) - np.sum(profile * dp)
    assert np.isclose(change, dp[0] * source, rtol=1e-10)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -v
```

Expected: all four FAIL. The three that pass nine arguments fail with a `TypeError` about `_diffuse_profile` taking 7 positional arguments; `test_diffuse_profile_zero_boundary_terms_reproduce_no_flux` fails the same way.

- [ ] **Step 3: Add the boundary terms to `_diffuse_profile`**

Replace `_diffuse_profile` in `climt/_components/simple_boundary_layer/component.py` with:

```python
@jit_compile
def _diffuse_profile(
    profile, p, p_int, rho, diff, dt, g, surface_exchange, surface_source
):
    """Implicit vertical diffusion of ``profile`` with a surface boundary term.

    Assembles the tridiagonal system and solves it with the shared Thomas
    solver. ``rho``/``diff`` have length ``num_layers - 1``.

    The lowest layer carries the surface boundary condition. With
    ``dp0 = p_int[0] - p_int[1]`` the layer-0 mass:

    * ``surface_exchange`` is the implicit bulk-exchange coefficient
      ``beta = g * rho_s * C * |v| * dt / dp0``, added to ``diag[0]``. It makes
      an unknown-dependent flux ``rho_s C |v| (X_s - X_0)`` backward-Euler.
    * ``surface_source`` is the explicit right-hand-side addition: ``beta * X_s``
      for that bulk flux, or ``g * dt * F / dp0`` for a prescribed flux ``F``
      (positive upward into the atmosphere).

    Both zero reproduces the no-flux solve.
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
    diag[0] = diag[0] + surface_exchange
    rhs = profile.copy()
    rhs[0] = rhs[0] + surface_source
    return solve_tridiagonal(-diag_m[1:], diag, -diag_p[:-1], rhs)
```

- [ ] **Step 4: Pass `0.0, 0.0` at the four existing call sites**

In `_boundary_layer_kernel`, replace the four `_diffuse_profile` calls (currently lines 173–184) with:

```python
        new_air_temperature[:, col] = _diffuse_profile(
            col_T, col_p, col_p_int, col_rho, diff, timestep, g, 0.0, 0.0
        )
        new_specific_humidity[:, col] = _diffuse_profile(
            col_q, col_p, col_p_int, col_rho, diff, timestep, g, 0.0, 0.0
        )
        new_northward_wind[:, col] = _diffuse_profile(
            col_v, col_p, col_p_int, col_rho, diff, timestep, g, 0.0, 0.0
        )
        new_eastward_wind[:, col] = _diffuse_profile(
            col_u, col_p, col_p_int, col_rho, diff, timestep, g, 0.0, 0.0
        )
```

- [ ] **Step 5: Run the tests to verify they pass**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -v
```

Expected: 4 passed.

- [ ] **Step 6: Verify the existing component tests still pass unchanged**

```bash
conda run -n climt python -m pytest tests/test_components.py -k "boundary_layer or BoundaryLayer" -v
```

Expected: all pass. Behaviour is unchanged at this point — the boundary terms are zero everywhere.

- [ ] **Step 7: Commit**

```bash
git add climt/_components/simple_boundary_layer/component.py tests/test_simple_boundary_layer.py
git commit -m "feat(sbl): add surface boundary terms to _diffuse_profile"
```

---

## Task 2: Bulk surface fluxes (`surface_fluxes='bulk'`, the new default)

**Files:**
- Modify: `climt/_components/simple_boundary_layer/component.py` — `BoundaryLayerParams`, `_boundary_layer_kernel`, `SimpleBoundaryLayer.__init__`, `_update_constants`, `array_call`
- Modify: `tests/test_simple_boundary_layer.py` (append)
- Modify: `tests/test_components.py:781-834`

**Interfaces:**
- Consumes: `_diffuse_profile(profile, p, p_int, rho, diff, dt, g, surface_exchange, surface_source)` from Task 1.
- Produces:
  - `BoundaryLayerParams(Rd, Cp, g, k, z0, fb, P0, Ric, Lv)` — a `NamedTuple` with `Lv: float` appended last.
  - `SimpleBoundaryLayer(surface_fluxes='bulk', von_karman_constant=0.4, roughness_length=0.0000321, specific_fraction=0.1, reference_pressure=100000, critical_richardson_number=1, **kwargs)`. `surface_fluxes` accepts `'bulk'`, `'external'` or `None`; anything else raises `ValueError`. Attributes set: `self._surface_fluxes` (the string/None) and `self._flux_mode` (int 0/1/2).
  - `_boundary_layer_kernel(air_temperature, surface_temperature, air_pressure, air_pressure_int, surface_pressure, specific_humidity, surface_humidity, northward_wind, eastward_wind, sensible_heat_flux, latent_heat_flux, new_air_temperature, new_specific_humidity, new_northward_wind, new_eastward_wind, north_wind_stress, east_wind_stress, boundary_height, applied_sensible_flux, applied_latent_flux, params, num_cols, timestep, flux_mode)`.
  - In `'bulk'` mode only, two extra diagnostics: `surface_upward_sensible_heat_flux` and `surface_upward_latent_heat_flux`, both `{'dims': ['*'], 'units': 'W m^-2'}`.
  - **Changed behaviour:** `northward_wind_stress` / `eastward_wind_stress` change value in every mode (see below).

**One rule for all four surface-flux diagnostics** (spec, "Surface-flux diagnostics"): each reports the flux the solve actually delivered, evaluated with the **post-solve** layer-0 value. The bulk flux is implicit, so that value is `X₀^new`, not `X₀^old`. With `κ = col_rho[0] * C * wind_int[0]` shared by all four:

```
surface_upward_sensible_heat_flux = cp  · κ · (T_s − T₀^new)
surface_upward_latent_heat_flux   = L_v · κ · (q_s − q₀^new)
eastward_wind_stress              = κ · u₀^new
northward_wind_stress             = κ · v₀^new
```

Using `X₀^old` would make the budget and consistency tests approximate instead of exact. Do not.

**The wind stresses are a deliberate value change, in every mode including `None`.** Today they are `κ · u_int[0]` — the pre-step *interface* wind `0.5·(u₀^old + u₁^old)` — which never matched the momentum the column actually lost. One formula now serves all three modes; in `None` mode nothing applies the stress, so it is advisory only. No in-tree code reads these diagnostics.

**Positional-argument note.** `surface_fluxes` goes **first** in the signature, which shifts the positional order of the existing tuning parameters. No caller in this repo passes them positionally (`examples/full_gcm_land_ocean_ice.py:51` and every test use `climt.SimpleBoundaryLayer()`), and the component is one commit old.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_simple_boundary_layer.py`:

```python
# ------------------------------------------------------------------ helpers

def _constants():
    return (
        get_constant('gravitational_acceleration', 'm s^-2'),
        get_constant(
            'heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1'
        ),
        get_constant('latent_heat_of_condensation', 'J kg^-1'),
    )


def _column_state(component, nz=30):
    return climt.get_default_state(
        [component], grid_state=climt.get_grid(nx=None, ny=None, nz=nz)
    )


def _layer_mass(state):
    """dp for each mid-level of a single column (Pa)."""
    p_int = np.asarray(state['air_pressure_on_interface_levels'])
    return p_int[:-1] - p_int[1:]


def _column_budgets(state, dp):
    """(column enthalpy J/m2, column water kg/m2) for a single column."""
    g, cp, _ = _constants()
    enthalpy = cp * np.sum(np.asarray(state['air_temperature']) * dp) / g
    water = np.sum(np.asarray(state['specific_humidity']) * dp) / g
    return enthalpy, water


def _column_momentum(state, dp, name):
    """Column-integrated momentum of one wind component, kg m^-1 s^-1."""
    g, _, _ = _constants()
    return np.sum(np.asarray(state[name]) * dp) / g


def _scalar(field):
    return float(np.asarray(field).ravel()[0])


def _warm_wet_surface_state(component, nz=30):
    """320 K saturated-ish surface under a 250 K, bone-dry, windy atmosphere.

    This is the empirical case from the design doc that demonstrated the bug.
    """
    state = _column_state(component, nz=nz)
    np.asarray(state['air_temperature'])[:] = 250.0
    np.asarray(state['specific_humidity'])[:] = 0.0
    np.asarray(state['eastward_wind'])[:] = 10.0
    np.asarray(state['northward_wind'])[:] = 0.0
    np.asarray(state['surface_temperature'])[:] = 320.0
    np.asarray(state['surface_specific_humidity'])[:] = 0.02
    return state


# --------------------------------------------------------------------- API

def test_invalid_surface_fluxes_mode_raises():
    with pytest.raises(ValueError, match='surface_fluxes'):
        climt.SimpleBoundaryLayer(surface_fluxes='bogus')


def test_bulk_is_the_default_mode():
    assert climt.SimpleBoundaryLayer()._surface_fluxes == 'bulk'


def test_bulk_mode_diagnoses_surface_fluxes():
    c = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    assert set(c.diagnostic_properties) == {
        'northward_wind_stress',
        'eastward_wind_stress',
        'boundary_layer_height',
        'surface_upward_sensible_heat_flux',
        'surface_upward_latent_heat_flux',
    }
    assert (
        c.diagnostic_properties['surface_upward_sensible_heat_flux']['units']
        == 'W m^-2'
    )
    # bulk mode must NOT require the fluxes as inputs
    assert 'surface_upward_sensible_heat_flux' not in c.input_properties


def test_none_mode_has_no_flux_diagnostics():
    c = climt.SimpleBoundaryLayer(surface_fluxes=None)
    assert set(c.diagnostic_properties) == {
        'northward_wind_stress',
        'eastward_wind_stress',
        'boundary_layer_height',
    }
    assert 'surface_upward_sensible_heat_flux' not in c.input_properties


# ------------------------------------------------------------ bulk physics

def test_bulk_column_budget_closes_against_diagnosed_fluxes():
    """Column enthalpy/water change == diagnosed surface flux x dt, exactly."""
    _, _, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)

    enthalpy_before, water_before = _column_budgets(state, dp)
    timestep = timedelta(seconds=1200)
    diagnostics, new_state = component(state, timestep=timestep)
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    sh = _scalar(diagnostics['surface_upward_sensible_heat_flux'])
    lh = _scalar(diagnostics['surface_upward_latent_heat_flux'])
    dt = timestep.total_seconds()

    assert np.isclose(enthalpy_after - enthalpy_before, sh * dt, rtol=1e-10)
    assert np.isclose(water_after - water_before, lh * dt / lv, rtol=1e-10)


def test_bulk_momentum_budget_closes_against_the_stress_diagnostic():
    """Column momentum change == -wind stress x dt, exactly.

    The pre-reconciliation diagnostic (built from the pre-step *interface*
    wind) could not have passed this.
    """
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    np.asarray(state['northward_wind'])[:] = -4.0
    dp = _layer_mass(state)

    timestep = timedelta(seconds=1200)
    dt = timestep.total_seconds()
    before = {
        name: _column_momentum(state, dp, name)
        for name in ('eastward_wind', 'northward_wind')
    }
    diagnostics, new_state = component(state, timestep=timestep)

    for name, stress in (('eastward_wind', 'eastward_wind_stress'),
                         ('northward_wind', 'northward_wind_stress')):
        change = _column_momentum(new_state, dp, name) - before[name]
        assert np.isclose(
            change, -_scalar(diagnostics[stress]) * dt, rtol=1e-10
        ), name


def test_bulk_wind_stress_uses_the_post_solve_layer_zero_wind():
    """Guard the reconciliation against a revert to the old formula.

    kappa is recovered from the momentum budget, which the solver satisfies
    regardless of how the diagnostic is written -- so this is not circular.
    The old formula used the pre-step interface wind
    ``0.5 * (u0_old + u1_old)``; under shear that is a different number.
    """
    g, _, _ = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    # sheared profile: mid-level 0 is the surface, so this is weak wind at the
    # bottom and strong aloft, making u0 and the 0/1 interface clearly differ
    u = np.asarray(state['eastward_wind'])
    u[:] = np.linspace(2.0, 30.0, u.shape[0])
    u_old_0, u_old_1 = u[0], u[1]
    u_interface_old = 0.5 * (u_old_0 + u_old_1)

    dp = _layer_mass(state)
    before = _column_momentum(state, dp, 'eastward_wind')
    timestep = timedelta(seconds=1200)
    dt = timestep.total_seconds()
    diagnostics, new_state = component(state, timestep=timestep)

    tau = _scalar(diagnostics['eastward_wind_stress'])
    u_new_0 = np.asarray(new_state['eastward_wind'])[0]
    momentum_lost = before - _column_momentum(new_state, dp, 'eastward_wind')
    kappa = momentum_lost / (u_new_0 * dt)

    assert tau > 0.0
    assert np.isclose(tau, kappa * u_new_0, rtol=1e-10)
    # the two candidate winds must be distinguishable, so that the assertion
    # below actually discriminates (which one is larger depends on how drag
    # and shear-driven diffusion compete -- do not assume an ordering)
    assert not np.isclose(u_interface_old, u_new_0, rtol=1e-2)
    # the old form would have reported a visibly different number
    assert not np.isclose(tau, kappa * u_interface_old, rtol=1e-3)


def test_none_mode_stress_is_diagnosed_but_not_applied():
    component = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)
    before = _column_momentum(state, dp, 'eastward_wind')

    diagnostics, new_state = component(state, timestep=timedelta(seconds=1200))

    # advisory: a real stress is reported ...
    assert _scalar(diagnostics['eastward_wind_stress']) > 0.0
    # ... but no momentum left the column
    after = _column_momentum(new_state, dp, 'eastward_wind')
    assert np.isclose(after, before, rtol=1e-12)


def test_bulk_mode_warms_moistens_and_decelerates():
    """The design doc's empirical case, with the assertions inverted."""
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)
    u_before = np.asarray(state['eastward_wind'])[0]

    for _ in range(20):
        _, new_state = component(state, timestep=timedelta(seconds=1200))
        state.update(new_state)

    enthalpy_after, water_after = _column_budgets(state, dp)
    assert enthalpy_after > enthalpy_before
    assert water_after > water_before
    assert np.asarray(state['air_temperature'])[0] > 250.0
    assert np.asarray(state['air_temperature'])[0] < 320.0
    assert abs(np.asarray(state['eastward_wind'])[0]) < abs(u_before)
    assert np.all(np.isfinite(np.asarray(state['air_temperature'])))


def test_bulk_mode_cools_column_over_a_cold_surface():
    component = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    state = _column_state(component)
    np.asarray(state['air_temperature'])[:] = 290.0
    np.asarray(state['surface_temperature'])[:] = 250.0
    np.asarray(state['eastward_wind'])[:] = 10.0
    dp = _layer_mass(state)
    enthalpy_before, _ = _column_budgets(state, dp)

    _, new_state = component(state, timestep=timedelta(seconds=1200))
    enthalpy_after, _ = _column_budgets(new_state, dp)
    assert enthalpy_after < enthalpy_before


def test_none_mode_conserves_over_a_warm_wet_surface():
    """The regression guard: opting out must restore exact conservation."""
    component = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = _warm_wet_surface_state(component)
    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)

    _, new_state = component(state, timestep=timedelta(seconds=1200))
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    assert np.isclose(enthalpy_after, enthalpy_before, rtol=1e-12)
    assert np.isclose(water_after, water_before, rtol=1e-12)
    assert np.asarray(new_state['air_temperature'])[0] == pytest.approx(250.0)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -v
```

Expected: the four Task 1 tests pass; the new ones FAIL — `test_invalid_surface_fluxes_mode_raises` with `TypeError: __init__() got an unexpected keyword argument 'surface_fluxes'`, and the rest likewise.

- [ ] **Step 3: Add `Lv` to `BoundaryLayerParams`**

```python
class BoundaryLayerParams(NamedTuple):
    Rd: float
    Cp: float
    g: float
    k: float
    z0: float
    fb: float
    P0: float
    Ric: float
    Lv: float
```

- [ ] **Step 4: Add the mode branch to the kernel**

Change the `_boundary_layer_kernel` signature to:

```python
@jit_compile(parallel=True)
def _boundary_layer_kernel(
    air_temperature, surface_temperature, air_pressure, air_pressure_int,
    surface_pressure, specific_humidity, surface_humidity, northward_wind,
    eastward_wind, sensible_heat_flux, latent_heat_flux,
    new_air_temperature, new_specific_humidity, new_northward_wind,
    new_eastward_wind, north_wind_stress, east_wind_stress, boundary_height,
    applied_sensible_flux, applied_latent_flux,
    params, num_cols, timestep, flux_mode,
):
```

Add `Lv = params.Lv` alongside the other `params` unpacking at the top of the function.

**Delete** the two existing stress-diagnostic lines (currently lines 157–158, immediately above `u_fric = wind_int[0]`):

```python
        north_wind_stress[col] = col_rho[0] * C * wind_int[0] * col_v_int[0]
        east_wind_stress[col] = col_rho[0] * C * wind_int[0] * col_u_int[0]
```

They move below the solve, where the post-solve winds exist (Step 4b). Note this makes `col_v_int` / `col_u_int` unused outside the `wind_int` calculation — leave them, `wind_int` still needs them.

Then, between the now-first `u_fric = wind_int[0]` line and the `diff = np.zeros(n)` block, insert the boundary-condition selection:

```python
        dp0 = col_p_int[0] - col_p_int[1]
        bulk_conductance = col_rho[0] * C * wind_int[0]
        beta = g * bulk_conductance * timestep / dp0

        if flux_mode == 1:
            scalar_exchange = beta
            source_T = beta * col_Ts
            source_q = beta * col_qs
        elif flux_mode == 2:
            scalar_exchange = 0.0
            source_T = g * timestep * sensible_heat_flux[col] / (Cp * dp0)
            source_q = g * timestep * latent_heat_flux[col] / (Lv * dp0)
        else:
            scalar_exchange = 0.0
            source_T = 0.0
            source_q = 0.0

        if flux_mode == 0:
            wind_exchange = 0.0
        else:
            wind_exchange = beta
```

Replace the four `_diffuse_profile` calls from Task 1 Step 4 with:

```python
        new_air_temperature[:, col] = _diffuse_profile(
            col_T, col_p, col_p_int, col_rho, diff, timestep, g,
            scalar_exchange, source_T,
        )
        new_specific_humidity[:, col] = _diffuse_profile(
            col_q, col_p, col_p_int, col_rho, diff, timestep, g,
            scalar_exchange, source_q,
        )
        new_northward_wind[:, col] = _diffuse_profile(
            col_v, col_p, col_p_int, col_rho, diff, timestep, g,
            wind_exchange, 0.0,
        )
        new_eastward_wind[:, col] = _diffuse_profile(
            col_u, col_p, col_p_int, col_rho, diff, timestep, g,
            wind_exchange, 0.0,
        )

        applied_sensible_flux[col] = Cp * bulk_conductance * (
            col_Ts - new_air_temperature[0, col]
        )
        applied_latent_flux[col] = Lv * bulk_conductance * (
            col_qs - new_specific_humidity[0, col]
        )
        # Post-solve layer-0 winds, so the stress equals the momentum the
        # column actually lost. The surface value for momentum is no-slip
        # (X_s == 0), hence no (0 - u) difference here -- just the drag.
        north_wind_stress[col] = (
            bulk_conductance * new_northward_wind[0, col]
        )
        east_wind_stress[col] = bulk_conductance * new_eastward_wind[0, col]
```

- [ ] **Step 5: Add the `surface_fluxes` API to the component**

Add a module-level mapping just above `class SimpleBoundaryLayer`:

```python
_FLUX_MODES = {None: 0, 'bulk': 1, 'external': 2}
```

Replace `__init__` and `_update_constants` with:

```python
    def __init__(self, surface_fluxes='bulk', von_karman_constant=0.4,
                 roughness_length=0.0000321, specific_fraction=0.1,
                 reference_pressure=100000, critical_richardson_number=1,
                 **kwargs):
        """
        Args:
            surface_fluxes: how surface fluxes enter the lowest model level.

                * ``'bulk'`` (default): the component computes Frierson (2006)
                  bulk fluxes of heat, moisture and momentum from its own
                  exchange coefficient and applies them implicitly. It then
                  reports the applied heat and moisture fluxes as diagnostics.
                * ``'external'``: ``surface_upward_sensible_heat_flux`` and
                  ``surface_upward_latent_heat_flux`` become required inputs and
                  are applied as prescribed fluxes. Momentum still uses the
                  internal bulk drag.
                * ``None``: no surface exchange at all. The diffusion has
                  no-flux boundaries and conserves every column integral. A
                  separate surface-flux component must supply the fluxes.

                Pairing ``'bulk'`` with a component that already applies
                surface fluxes -- ``SimplePhysics(surface_fluxes=True)``, which
                is its default -- applies them twice.
            von_karman_constant: von Karman constant k.
            roughness_length: surface roughness length z0 (m).
            specific_fraction: surface-layer fraction fb of the boundary
                layer depth.
            reference_pressure: reference pressure P0 (Pa) for potential
                temperature.
            critical_richardson_number: critical Richardson number Ric that
                caps the diffusion and sets the boundary-layer top.
        """
        if surface_fluxes not in _FLUX_MODES:
            raise ValueError(
                "surface_fluxes must be 'bulk', 'external' or None, got %r"
                % (surface_fluxes,)
            )
        self._surface_fluxes = surface_fluxes
        self._flux_mode = _FLUX_MODES[surface_fluxes]
        self._k = von_karman_constant
        self._z0 = roughness_length
        self._fb = specific_fraction
        self._P0 = reference_pressure
        self._Ric = critical_richardson_number

        if surface_fluxes == 'bulk':
            # Instance-level override: only bulk mode reports the fluxes it
            # applied. In 'external' mode these names are inputs, so they must
            # not also be diagnostics.
            self.diagnostic_properties = dict(self.diagnostic_properties)
            self.diagnostic_properties['surface_upward_sensible_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }
            self.diagnostic_properties['surface_upward_latent_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }

        self._update_constants()
        super(SimpleBoundaryLayer, self).__init__(**kwargs)

    def _update_constants(self):
        self._Rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        self._Cp = get_constant(
            'heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1'
        )
        self._g = get_constant('gravitational_acceleration', 'm s^-2')
        self._Lv = get_constant('latent_heat_of_condensation', 'J kg^-1')
```

(The `'external'` branch of `input_properties` lands in Task 3.)

- [ ] **Step 6: Wire `array_call` to the new kernel signature**

Replace `array_call`'s body from the `params = ...` line onward with:

```python
        params = BoundaryLayerParams(
            Rd=self._Rd, Cp=self._Cp, g=self._g, k=self._k, z0=self._z0,
            fb=self._fb, P0=self._P0, Ric=self._Ric, Lv=self._Lv,
        )

        # numba needs real arrays for every branch, so unused slots are zeros.
        zeros = np.zeros(num_cols)
        if self._flux_mode == 2:
            sensible_flux = state['surface_upward_sensible_heat_flux']
            latent_flux = state['surface_upward_latent_heat_flux']
        else:
            sensible_flux = zeros
            latent_flux = zeros

        if self._flux_mode == 1:
            applied_sensible = diagnostics['surface_upward_sensible_heat_flux']
            applied_latent = diagnostics['surface_upward_latent_heat_flux']
        else:
            applied_sensible = np.zeros(num_cols)
            applied_latent = np.zeros(num_cols)

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
            sensible_flux,
            latent_flux,
            new_state['air_temperature'],
            new_state['specific_humidity'],
            new_state['northward_wind'],
            new_state['eastward_wind'],
            diagnostics['northward_wind_stress'],
            diagnostics['eastward_wind_stress'],
            diagnostics['boundary_layer_height'],
            applied_sensible,
            applied_latent,
            params,
            num_cols,
            timestep.total_seconds(),
            self._flux_mode,
        )

        return diagnostics, new_state
```

- [ ] **Step 7: Run the tests to verify they pass**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -v
```

Expected: all pass. If numba raises a typing error about `sensible_flux`, confirm every branch receives a 1-D float64 array — never `None`, never a scalar.

- [ ] **Step 8: Update the two existing tests in `tests/test_components.py`**

Pin the conservation test to `None` mode and extend the properties test. Replace `tests/test_components.py:781-834` with:

```python
@pytest.mark.parametrize("mode", ["bulk", "external", None])
def test_simple_boundary_layer_properties(mode):
    c = climt.SimpleBoundaryLayer(surface_fluxes=mode)
    expected_inputs = {
        "air_temperature",
        "specific_humidity",
        "air_pressure",
        "air_pressure_on_interface_levels",
        "northward_wind",
        "eastward_wind",
        "surface_air_pressure",
        "surface_temperature",
        "surface_specific_humidity",
    }
    expected_diagnostics = {
        "northward_wind_stress",
        "eastward_wind_stress",
        "boundary_layer_height",
    }
    if mode == "external":
        expected_inputs |= {
            "surface_upward_sensible_heat_flux",
            "surface_upward_latent_heat_flux",
        }
    if mode == "bulk":
        expected_diagnostics |= {
            "surface_upward_sensible_heat_flux",
            "surface_upward_latent_heat_flux",
        }
    assert set(c.input_properties) == expected_inputs
    assert set(c.output_properties) == {
        "air_temperature",
        "specific_humidity",
        "northward_wind",
        "eastward_wind",
    }
    assert set(c.diagnostic_properties) == expected_diagnostics
    assert c.diagnostic_properties["boundary_layer_height"]["units"] == "m"
    assert c.output_properties["air_temperature"]["units"] == "degK"


def test_simple_boundary_layer_conserves_column_integrals():
    # Conservation is a property of the no-flux mode only. The 'bulk' and
    # 'external' modes deliberately add mass and energy at the surface; their
    # budgets are checked in tests/test_simple_boundary_layer.py.
    component = climt.SimpleBoundaryLayer(surface_fluxes=None)
    state = climt.get_default_state(
        [component], grid_state=get_grid(nx=None, ny=None, nz=30)
    )
    # perturb winds/humidity so there is something to diffuse
    np.asarray(state["northward_wind"])[:] = 5.0
    np.asarray(state["eastward_wind"])[:] = -3.0
    np.asarray(state["specific_humidity"])[:] = 0.005

    diagnostics, new_state = component(
        state, timestep=timedelta(seconds=600)
    )

    p_int = np.asarray(state["air_pressure_on_interface_levels"])
    # mid-levels axis is axis 0 for a single column here
    dp = p_int[:-1] - p_int[1:]
    for name in (
        "air_temperature",
        "specific_humidity",
        "northward_wind",
        "eastward_wind",
    ):
        before = np.sum(np.asarray(state[name]) * dp)
        after = np.sum(np.asarray(new_state[name]) * dp)
        assert np.isclose(before, after, rtol=1e-8, atol=1e-6), name
```

The `test_simple_boundary_layer_properties` parametrisation covers the `'external'` input set that Task 3 implements, so that case fails until then — that is intended and is the failing test Task 3 starts from.

- [ ] **Step 9: Run the component tests**

```bash
conda run -n climt python -m pytest tests/test_components.py -k "boundary_layer or BoundaryLayer" -v
```

Expected: all pass **except** `test_simple_boundary_layer_properties[external]`, which fails on the missing input properties. Task 3 fixes it. Everything else — including `TestSimpleBoundaryLayer`'s inherited `ComponentBaseColumn` / `ComponentBase3D` cases against the default (now bulk) component — must pass.

- [ ] **Step 10: Commit**

```bash
git add climt/_components/simple_boundary_layer/component.py tests/test_simple_boundary_layer.py tests/test_components.py
git commit -m "feat(sbl): apply Frierson bulk surface fluxes, default surface_fluxes='bulk'"
```

---

## Task 3: External surface fluxes (`surface_fluxes='external'`)

**Files:**
- Modify: `climt/_components/simple_boundary_layer/component.py` — the `__init__` mode branch only
- Modify: `tests/test_simple_boundary_layer.py` (append)

**Interfaces:**
- Consumes: everything from Task 2. The kernel's `flux_mode == 2` branch already exists and is already wired to `sensible_heat_flux` / `latent_heat_flux`; this task only makes the component declare the inputs and pass the real arrays.
- Produces: in `'external'` mode, `input_properties` additionally contains `surface_upward_sensible_heat_flux` and `surface_upward_latent_heat_flux`, both `{'dims': ['*'], 'units': 'W m^-2'}`. Both already have `get_default_state` defaults of `0.0` (`climt/_core/initialization.py:917-925`), so a default state constructs cleanly.

**Why bulk and external agree exactly.** Bulk assembles `(1 + Σ + β)·X₀ − … = X₀ᵒˡᵈ + β·X_s`. Rearranged, that is `(1 + Σ)·X₀ − … = X₀ᵒˡᵈ + β·(X_s − X₀ⁿᵉʷ)`, which is precisely the external system fed `F = ρ_s·C·|v|·(X_s − X₀ⁿᵉʷ)` — the flux Task 2 diagnoses. The consistency test below therefore holds to solver round-off, not to a tolerance.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_simple_boundary_layer.py`:

```python
# -------------------------------------------------------- external physics

def test_external_mode_requires_the_flux_inputs():
    c = climt.SimpleBoundaryLayer(surface_fluxes='external')
    assert 'surface_upward_sensible_heat_flux' in c.input_properties
    assert 'surface_upward_latent_heat_flux' in c.input_properties
    assert (
        c.input_properties['surface_upward_latent_heat_flux']['units']
        == 'W m^-2'
    )
    # they are inputs, not diagnostics -- the two must not collide
    assert 'surface_upward_sensible_heat_flux' not in c.diagnostic_properties


def test_external_mode_applies_prescribed_fluxes_exactly():
    _, _, lv = _constants()
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _column_state(component)
    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = 50.0
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = 100.0

    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)
    timestep = timedelta(seconds=1200)
    _, new_state = component(state, timestep=timestep)
    enthalpy_after, water_after = _column_budgets(new_state, dp)
    dt = timestep.total_seconds()

    assert np.isclose(enthalpy_after - enthalpy_before, 50.0 * dt, rtol=1e-10)
    assert np.isclose(
        water_after - water_before, 100.0 * dt / lv, rtol=1e-10
    )


def test_external_mode_zero_fluxes_conserve():
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _column_state(component)
    np.asarray(state['specific_humidity'])[:] = 0.005
    dp = _layer_mass(state)
    enthalpy_before, water_before = _column_budgets(state, dp)

    _, new_state = component(state, timestep=timedelta(seconds=1200))
    enthalpy_after, water_after = _column_budgets(new_state, dp)

    assert np.isclose(enthalpy_after, enthalpy_before, rtol=1e-12)
    assert np.isclose(water_after, water_before, rtol=1e-12)


def test_external_mode_still_applies_bulk_momentum_drag():
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _warm_wet_surface_state(component)
    # fluxes stay zero; only drag should act
    u_before = np.asarray(state['eastward_wind'])[0]
    _, new_state = component(state, timestep=timedelta(seconds=1200))
    assert abs(np.asarray(new_state['eastward_wind'])[0]) < abs(u_before)


def test_external_momentum_budget_closes_against_the_stress_diagnostic():
    """Momentum is bulk-internal in external mode, so its budget closes too."""
    component = climt.SimpleBoundaryLayer(surface_fluxes='external')
    state = _warm_wet_surface_state(component)
    np.asarray(state['northward_wind'])[:] = -4.0
    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = 50.0
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = 100.0
    dp = _layer_mass(state)

    timestep = timedelta(seconds=1200)
    dt = timestep.total_seconds()
    before = {
        name: _column_momentum(state, dp, name)
        for name in ('eastward_wind', 'northward_wind')
    }
    diagnostics, new_state = component(state, timestep=timestep)

    for name, stress in (('eastward_wind', 'eastward_wind_stress'),
                         ('northward_wind', 'northward_wind_stress')):
        change = _column_momentum(new_state, dp, name) - before[name]
        assert np.isclose(
            change, -_scalar(diagnostics[stress]) * dt, rtol=1e-10
        ), name


def test_bulk_and_external_agree_when_fed_the_bulk_fluxes():
    bulk = climt.SimpleBoundaryLayer(surface_fluxes='bulk')
    external = climt.SimpleBoundaryLayer(surface_fluxes='external')
    # build the state from the external component so it carries the flux
    # fields; bulk simply ignores the extras.
    state = _warm_wet_surface_state(external)

    diagnostics, bulk_new = bulk(state, timestep=timedelta(seconds=1200))

    np.asarray(state['surface_upward_sensible_heat_flux'])[:] = _scalar(
        diagnostics['surface_upward_sensible_heat_flux']
    )
    np.asarray(state['surface_upward_latent_heat_flux'])[:] = _scalar(
        diagnostics['surface_upward_latent_heat_flux']
    )
    external_diagnostics, external_new = external(
        state, timestep=timedelta(seconds=1200)
    )

    for name in ('air_temperature', 'specific_humidity',
                 'eastward_wind', 'northward_wind'):
        assert np.allclose(
            np.asarray(bulk_new[name]), np.asarray(external_new[name]),
            rtol=1e-10, atol=1e-12,
        ), name
    # momentum is bulk-internal in both modes, so the stresses match too
    for name in ('eastward_wind_stress', 'northward_wind_stress'):
        assert np.isclose(
            _scalar(diagnostics[name]),
            _scalar(external_diagnostics[name]),
            rtol=1e-10,
        ), name
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -k external -v
conda run -n climt python -m pytest "tests/test_components.py::test_simple_boundary_layer_properties[external]" -v
```

Expected: `test_external_mode_requires_the_flux_inputs` FAILs on the missing key; the budget tests FAIL with a `KeyError` on `surface_upward_sensible_heat_flux` (the default state does not carry a field the component never declared); the `test_components.py` case FAILs on the input set.

- [ ] **Step 3: Declare the external-mode inputs**

In `SimpleBoundaryLayer.__init__`, extend the mode branch added in Task 2 (`if surface_fluxes == 'bulk': ...`) with an `elif`:

```python
        elif surface_fluxes == 'external':
            # Instance-level override: only external mode consumes fluxes
            # computed by a separate surface component. The fluxes were
            # computed with that component's own drag coefficient (e.g.
            # BucketHydrology's bulk_coefficient=0.0011), which differs from
            # this component's Monin-Obukhov C. That inconsistency is inherent
            # to the modular split and is not reconciled here.
            self.input_properties = dict(self.input_properties)
            self.input_properties['surface_upward_sensible_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }
            self.input_properties['surface_upward_latent_heat_flux'] = {
                'dims': ['*'], 'units': 'W m^-2',
            }
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py -v
conda run -n climt python -m pytest tests/test_components.py -k "boundary_layer or BoundaryLayer" -v
```

Expected: everything passes, including `test_simple_boundary_layer_properties[external]`.

- [ ] **Step 5: Commit**

```bash
git add climt/_components/simple_boundary_layer/component.py tests/test_simple_boundary_layer.py
git commit -m "feat(sbl): accept externally supplied surface heat and moisture fluxes"
```

---

## Task 4: Documentation

The class docstring currently contradicts itself ("This is the surface-flux / boundary-layer formulation…" followed by "assumes a surface-flux component has already run"). The user guide makes an unconditional conservation claim that is now mode-dependent.

**Files:**
- Modify: `climt/_components/simple_boundary_layer/component.py` (class docstring, lines 188–200)
- Modify: `docs/user-guide/simple_boundary_layer.qmd` (lines 15–18, 73–103, and the usage section)

**Interfaces:**
- Consumes: the final `surface_fluxes` API from Tasks 2 and 3. No code changes.

- [ ] **Step 1: Rewrite the class docstring**

Replace the `SimpleBoundaryLayer` class docstring with:

```python
    """A simple boundary-layer scheme that exchanges heat, moisture and
    momentum with the surface and diffuses them through the boundary layer.

    This is the surface-flux / boundary-layer formulation of
    Frierson, Held & Zurita-Gotor (2006), with the surface-layer diffusion
    coefficient modified (thesis Eqn 2.8) to use the surface-layer Richardson
    number in its multiplier, making it continuous at ``Ri_a == 0``.

    Diffusivities come from a simplified Monin-Obukhov theory with a K-profile
    capped by a critical Richardson number. How the surface enters the lowest
    model level is set by ``surface_fluxes``:

    * ``'bulk'`` (default) -- the component computes the bulk fluxes itself
      from its own exchange coefficient ``C``, the interface wind speed and the
      surface density, and applies them implicitly. It reports the applied heat
      and moisture fluxes as ``surface_upward_sensible_heat_flux`` and
      ``surface_upward_latent_heat_flux`` diagnostics. Momentum uses a no-slip
      surface value, which turns the wind-stress calculation into real drag.
    * ``'external'`` -- the two flux quantities become required *inputs* and are
      applied as prescribed fluxes, so a surface component's fluxes reach the
      atmosphere. Momentum still uses the internal bulk drag.
    * ``None`` -- no surface exchange; no-flux boundaries, and every column
      integral is conserved exactly. This requires a separate component to have
      already applied the surface fluxes.

    Every surface-flux diagnostic -- the two above plus
    ``northward_wind_stress`` and ``eastward_wind_stress`` -- reports the flux
    the solve actually delivered, evaluated with the post-solve layer-0 value,
    so the column budgets close to round-off. In ``None`` mode the stresses are
    advisory: they are diagnosed but nothing applies them.

    Do not pair ``'bulk'`` with another component that also applies surface
    fluxes -- they will be applied twice.
    """
```

- [ ] **Step 2: Update the introduction in the user guide**

In `docs/user-guide/simple_boundary_layer.qmd`, replace the paragraph at lines 15–18 ("The component assumes that a **surface-flux component has already run**…") with:

```markdown
By default the component computes the **bulk surface fluxes itself** and applies
them at the lowest model level, then diffuses the result vertically using eddy
diffusivities from a simplified Monin–Obukhov theory, capped by a critical
Richardson number. It can alternatively apply surface fluxes supplied by another
component, or skip surface exchange entirely — see
[Surface fluxes](#surface-fluxes) below.
```

- [ ] **Step 3: Add a Surface fluxes section**

Insert this new section in `docs/user-guide/simple_boundary_layer.qmd` immediately before the "### Implicit diffusion and conservation" heading:

````markdown
### Surface fluxes {#surface-fluxes}

The `surface_fluxes` argument selects how the surface enters the lowest model
level. With $\Delta p_0 = p_{\text{int},0} - p_{\text{int},1}$ the mass of the
lowest layer, and $F_{\text{surf}}$ the upward flux of $X$ into it,

$$
\frac{\Delta p_0}{g}\,\frac{\mathrm{d}X_0}{\mathrm{d}t} = F_{\text{surf}} - F_{\text{top}}.
$$

**`'bulk'` (default).** $F_{\text{surf}} = \rho_s\,C\,|\mathbf{v}|\,(X_s - X_0)$
depends on the unknown, so it is taken implicitly: with
$\beta = g\,\rho_s\,C\,|\mathbf{v}|\,\Delta t / \Delta p_0$ the tridiagonal
system gains $\beta$ on its first diagonal entry and $\beta X_s$ on its first
right-hand-side entry. $X_s$ is $T_s$ for temperature, $q_s$ for humidity, and
$0$ for both wind components (no-slip) — which is what turns the wind-stress
calculation into actual drag on the flow. $\rho_s$ is the density at the first
interface and $|\mathbf{v}|$ the interface wind speed, subject to the same
1 m s$^{-1}$ gustiness floor used elsewhere in the scheme.

**All four surface-flux diagnostics follow one rule:** each reports the flux the
solve actually delivered, evaluated with the *post-solve* layer-0 value. With
$\kappa = \rho_s\,C\,|\mathbf{v}|$,

$$
\mathrm{SH} = c_p\,\kappa\,(T_s - T_0^{\text{new}}),
\quad
\mathrm{LH} = L_v\,\kappa\,(q_s - q_0^{\text{new}}),
\quad
\tau_u = \kappa\,u_0^{\text{new}},
\quad
\tau_v = \kappa\,v_0^{\text{new}}.
$$

This is what makes every column budget close to round-off rather than
approximately.

::: {.callout-important}
## Changed: the wind-stress diagnostics
`northward_wind_stress` and `eastward_wind_stress` previously used the
*pre-step interface* wind, $\kappa\,\cdot\,\tfrac{1}{2}(u_0 + u_1)$, which did
not match the momentum the column actually lost. They now use the post-solve
layer-0 wind, so the momentum budget closes like the heat and moisture ones.
**Their values change in every mode**, including `surface_fluxes=None` — where
no drag is applied, so the stress is advisory only: the bulk stress a separate
surface component could apply.
:::

**`'external'`.** The two flux quantities become required **inputs**, in
W m$^{-2}$ positive upward, and enter as a pure Neumann source:

$$
\text{rhs}_T[0] \mathrel{+}= \frac{g\,\Delta t\,\mathrm{SH}}{c_p\,\Delta p_0},
\qquad
\text{rhs}_q[0] \mathrel{+}= \frac{g\,\Delta t\,\mathrm{LH}}{L_v\,\Delta p_0},
$$

with $L_v = 2.5\times10^6$ J kg$^{-1}$. Momentum stays bulk-internal: surface
components do not compute wind stress, and `northward_wind_stress` /
`eastward_wind_stress` are already diagnostics of this component.

**`None`.** No surface exchange at all — the historical behaviour, requiring a
separate component to have already applied the surface fluxes.

::: {.callout-warning}
## Do not apply surface fluxes twice
`surface_fluxes='bulk'` alongside a component that already applies surface
fluxes (for example `SimplePhysics` with its surface fluxes enabled) applies
them twice. Use `surface_fluxes=None` in that case.
:::

::: {.callout-note}
## Two caveats in `'external'` mode
**Inconsistent exchange coefficients.** The surface component computed its
fluxes with its own drag coefficient — `BucketHydrology` uses
`bulk_coefficient=0.0011` — which differs from this component's Monin–Obukhov
$C$. That inconsistency is inherent to the modular split and is not reconciled.

**A one-step lag inside a dynamical core.** `sympl`'s component composite calls
every component with the *same* input state and merges diagnostics only into its
return value, never back into the state. So when `SimpleBoundaryLayer` is
wrapped in a `TimeDifferencingWrapper` and handed to `GFSDynamicalCore` next to
the surface component, it reads the surface component's *previous*-step fluxes
(and zero on the first step), provided the time loop does
`state.update(diagnostics)`. Calling `SimpleBoundaryLayer` directly in the time
loop, after the surface component and with a `state.update(...)` in between,
has no lag. Both patterns are legitimate — pick one knowingly.
:::
````

- [ ] **Step 4: Make the conservation claim mode-dependent**

In the "### Implicit diffusion and conservation" section, replace the sentence beginning "Because the boundary conditions are no-flux…" and the equation that follows with:

```markdown
Thomas solver. The interior discretisation is flux-form, so interior exchanges
telescope exactly and the only change to a column integral is what crosses the
surface:

$$
\sum_k \phi_k^{\text{(new)}}\,\Delta p_k - \sum_k \phi_k\,\Delta p_k
= g\,\Delta t\,F_{\text{surf}}(\phi),
\qquad \phi \in \{T,\ q,\ u,\ v\}.
$$

With `surface_fluxes=None` the right-hand side is zero and the scheme
**conserves the interface-pressure-weighted column integral** of every diffused
field. In the other two modes it closes the budget against the surface flux to
round-off.
```

- [ ] **Step 5: Update the inputs/outputs/diagnostics and parameters tables**

Replace the "Inputs, outputs and diagnostics" table with:

```markdown
| Kind | Quantities |
|------|------------|
| Inputs | `air_temperature`, `specific_humidity`, `air_pressure`, `air_pressure_on_interface_levels`, `northward_wind`, `eastward_wind`, `surface_air_pressure`, `surface_temperature`, `surface_specific_humidity` |
| Inputs (`'external'` only) | `surface_upward_sensible_heat_flux`, `surface_upward_latent_heat_flux` |
| Outputs | `air_temperature`, `specific_humidity`, `northward_wind`, `eastward_wind` |
| Diagnostics | `northward_wind_stress`, `eastward_wind_stress`, `boundary_layer_height` |
| Diagnostics (`'bulk'` only) | `surface_upward_sensible_heat_flux`, `surface_upward_latent_heat_flux` |
```

Add a first row to the Parameters table:

```markdown
| `surface_fluxes` | `'bulk'` | `'bulk'`, `'external'` or `None` — see [Surface fluxes](#surface-fluxes) |
```

- [ ] **Step 6: Update the usage example**

Replace the usage code block and the paragraph after it with:

````markdown
```python
import climt
from datetime import timedelta

# default: compute and apply bulk surface fluxes
boundary_layer = climt.SimpleBoundaryLayer()
state = climt.get_default_state(
    [boundary_layer],
    grid_state=climt.get_grid(nx=None, ny=None, nz=30),
)
diagnostics, new_state = boundary_layer(state, timestep=timedelta(minutes=10))
print(diagnostics["surface_upward_sensible_heat_flux"])

# or apply fluxes computed by a surface component
boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="external")

# or diffuse only, leaving surface exchange to someone else
boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes=None)
```

`SimpleBoundaryLayer` is a `Stepper`, so it returns the updated
temperature, humidity and wind profiles directly (not tendencies), alongside
the surface wind-stress and boundary-layer-height diagnostics — and, in the
default `'bulk'` mode, the surface heat and moisture fluxes it applied.
````

- [ ] **Step 7: Verify the docstring change did not break anything**

```bash
conda run -n climt python -m pytest tests/test_simple_boundary_layer.py tests/test_components.py -k "boundary_layer or BoundaryLayer" -q
```

Expected: all pass.

- [ ] **Step 8: Commit**

```bash
git add climt/_components/simple_boundary_layer/component.py docs/user-guide/simple_boundary_layer.qmd
git commit -m "docs(sbl): document the three surface_fluxes modes"
```

---

## Task 5: Migrate the full-GCM example (and fix the `Stepper`-in-dycore bug)

`examples/full_gcm_land_ocean_ice.py:59-63` hands `BucketHydrology` and `SimpleBoundaryLayer` — both `Stepper`s — to `GFSDynamicalCore`, whose composite accepts only `TendencyComponent` / `ImplicitTendencyComponent`. `build_model()` therefore raises `TypeError` and the example has never run, despite its docstring advertising a smoke test. `UpdateFrequencyWrapper` (used on the radiation components) preserves component type; `TimeDifferencingWrapper` is what converts a `Stepper`.

This task must both fix that and move the boundary layer to `surface_fluxes='external'` so `BucketHydrology`'s fluxes drive the atmosphere instead of the boundary layer computing a second, inconsistent set.

**Files:**
- Modify: `examples/full_gcm_land_ocean_ice.py:14-16` (docstring), `:36-37` (imports), `:44-67` (`build_model`)

**Interfaces:**
- Consumes: `SimpleBoundaryLayer(surface_fluxes='external')` from Task 3; `BucketHydrology`'s `surface_upward_sensible_heat_flux` / `surface_upward_latent_heat_flux` diagnostics (`climt/_components/bucket_hydrology/component.py:83-100`).

**This task is a real integration check, not a mechanical edit.** The verification step runs the example. If it fails, do not paper over it — see Step 4.

- [x] **Step 1: Wrap both `Stepper`s and switch to external fluxes**

In `examples/full_gcm_land_ocean_ice.py`, change the sympl import (line 37) to:

```python
from sympl import TimeDifferencingWrapper, UpdateFrequencyWrapper
```

Replace lines 50–63 of `build_model` with:

```python
    # BucketHydrology and SimpleBoundaryLayer are Steppers; GFSDynamicalCore's
    # component list accepts only TendencyComponent / ImplicitTendencyComponent,
    # so both need TimeDifferencingWrapper (UpdateFrequencyWrapper, used on the
    # radiation components below, preserves component type and would not do).
    land = TimeDifferencingWrapper(climt.BucketHydrology())
    # 'external' so the bucket's own surface fluxes drive the atmosphere rather
    # than the boundary layer computing a second, inconsistent set. Inside the
    # dycore this costs a one-step lag: sympl's composite hands every component
    # the same input state, so the boundary layer reads the previous step's
    # fluxes (zero on step 0).
    boundary_layer = TimeDifferencingWrapper(
        climt.SimpleBoundaryLayer(surface_fluxes="external")
    )
    convection = climt.EmanuelConvection()
    condensation = climt.GridScaleCondensation()
    radiation_lw = UpdateFrequencyWrapper(
        climt.RRTMGLongwave(), radiation_interval)
    radiation_sw = UpdateFrequencyWrapper(
        climt.RRTMGShortwave(), radiation_interval)

    dycore = GFSDynamicalCore(
        [land, radiation_sw, radiation_lw, convection, condensation,
         boundary_layer],
        number_of_damped_levels=5,
    )
```

- [x] **Step 2: Update the module docstring**

Replace lines 14–16 of the docstring with:

```
* ``SimpleBoundaryLayer`` -> boundary-layer diffusion (Frierson 2006), run with
                       ``surface_fluxes='external'`` so it applies the surface
                       heat and moisture fluxes ``BucketHydrology`` computes and
                       adds its own bulk momentum drag.
```

- [x] **Step 3: Run the example's smoke test**

```bash
conda run -n climt python examples/full_gcm_land_ocean_ice.py --steps 6 --nx 32 --ny 16
```

Expected: it runs to "Done. Everything ran and stayed finite." A reduced grid keeps the check quick; the point is that `build_model()` no longer raises and the coupled step is finite.

- [ ] **Step 4: If Step 3 fails, take the documented fallback** — ~~not needed~~

Step 3 passed on the primary path: `TimeDifferencingWrapper(BucketHydrology())`
inside the dycore was accepted and the run stayed finite. The fallback below was
**not** taken. Two extra fixes were needed that this task did not anticipate:

1. `GridScaleCondensation` is a `Stepper` too, so it also needed
   `TimeDifferencingWrapper`.
2. `get_default_state([dycore, mask, ocean, insolation], ...)` raised
   `InvalidPropertyDictError` on `latitude` — `GFSDynamicalCore` declares it
   `['lat', 'lon']` while `LandMask`/`Instellation` declare `['*']`. The state is
   now seeded from `[dycore]` alone, matching
   `examples/full_radiation_with_insolation_gcm.py`.

`TimeDifferencingWrapper(BucketHydrology())` turns the bucket's `surface_temperature` and `lwe_thickness_of_soil_moisture_content` outputs into tendencies that `GFSDynamicalCore` must be willing to integrate. If it rejects them, or the run goes non-finite because of them, move the bucket out of the dycore and into the time loop instead of forcing it — the bucket is a surface component and belongs there anyway. Note this does **not** remove the boundary layer's flux lag: the boundary layer stays inside the dycore, so it still reads the previous step's fluxes.

Fallback edit — drop `land` from the dycore list and return it separately:

```python
    land = climt.BucketHydrology()
    ...
    dycore = GFSDynamicalCore(
        [radiation_sw, radiation_lw, convection, condensation,
         boundary_layer],
        number_of_damped_levels=5,
    )
    grid = climt.get_grid(nx=nx, ny=ny)
    state = climt.get_default_state(
        [dycore, land, mask, ocean, insolation], grid_state=grid)
    return state, dycore, land, mask, ocean, insolation
```

and call it in `step()`, after the dycore so it sees this step's radiation:

```python
def step(state, dycore, land, ocean, insolation, timestep):
    """Advance the model one dynamics step."""
    state.update(insolation(state))     # zenith angle for this time
    state.update(ocean(state))          # prescribe looping SST on sea cells
    diagnostics, state = dycore(state, timestep)
    state.update(diagnostics)
    land_diagnostics, land_state = land(state, timestep)
    state.update(land_diagnostics)
    state.update(land_state)
    state["time"] += timestep
    return state
```

Update the two `build_model(...)` / `step(...)` call sites in `main()` (lines 154 and 163) to match the new signatures, then re-run Step 3.

**Report which path you took and paste the failure output if you took the fallback.** Do not silently change the physics.

- [x] **Step 5: Commit**

```bash
git add examples/full_gcm_land_ocean_ice.py
git commit -m "fix(example): wrap the Stepper physics and drive the atmosphere with bucket fluxes"
```

---

## Task 6: Full verification and knowledge-graph refresh

**Files:**
- Modify: `graphify-out/` (regenerated, not hand-edited)

**Interfaces:**
- Consumes: the complete implementation from Tasks 1–5.

- [x] **Step 1: Run the full test suite**

```bash
conda run -n climt python -m pytest tests/ -q
```

Result: **606 passed, 2 skipped, 1 failed** in 17m29s. The one failure is
`tests/test_hd209458b_reproduction.py::test_hd209458b_equilibrium_profile`
("FPI did not converge after 20 iterations in q & T"), which exercises
`SimplePhysics` + CORK and never imports `SimpleBoundaryLayer`. No commit in
this plan touches either. **The baseline comparison below was not run** — the
human partner asked for that test to be left alone — so it is untouched by this
work by inspection, not by measurement.

Expected: no failures. Compare against the pre-change baseline if anything looks pre-broken — `git stash` and re-run to confirm a failure is not yours before investigating.

- [x] **Step 2: Confirm no `scipy` import crept in**

```bash
conda run -n climt python -m pytest tests/test_no_scipy_import.py tests/test_pure_wheel_build.py -q
```

Expected: pass. The Pyodide two-wheel strategy depends on this.

- [x] **Step 3: Check the numba-disabled path still works**

The `jit_compile` decorator falls through to plain Python when numba is absent. Confirm the new kernel branch is not numba-only:

```bash
conda run -n climt env NUMBA_DISABLE_JIT=1 python -m pytest tests/test_simple_boundary_layer.py -q
```

Expected: pass.

- [x] **Step 4: Regenerate the knowledge graph** — ~~superseded~~

Dropped. graphify was removed from the repo during this task (`graphify-out/`
and `CLAUDE.md` deleted, `graphify-out/` added to `.gitignore`,
`scripts/augment_graph.py` kept for a possible revival). There is no graph left
to regenerate.

- [x] **Step 5: Commit** — done as the graphify-removal commit instead.

---

## Settled before implementation

Both items that were previously open are now decided in the spec — do not reopen them mid-task:

1. **The wind-stress diagnostics are reconciled** (spec, "Surface-flux diagnostics"). They now use the post-solve layer-0 wind, so the momentum budget closes like the heat and moisture budgets. This changes their value in every mode including `None`; Task 2 implements it, and Tasks 2–4 test and document it.

2. **`_richardson_diffusivity` keeps `Ri_a`** in the stability denominator (spec, "Resolved questions") — that is the choice that makes `K_b` continuous at `Ri_a == 0`. No code change, and no Abel thesis check needed.
