# Land Surface Components Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the land-surface half of the land-ocean-ice work: a two-layer (deep + shallow) mode for `BucketHydrology`, and a new modular, intermediate-complexity `SecondBEST` land surface model built from pluggable BEST process objects.

**Architecture:** `BucketHydrology` gains a backward-compatible `num_layers` option using two extra scalar state fields (no new grid dimension). `SecondBEST` is a thin `Stepper` orchestrator that composes five swappable process objects (`SoilProperties`, `SurfaceAlbedo`, `SurfaceLayer`, `SurfaceFluxes`, `SubsurfaceTransport`), each a plain class with a documented `__call__` and a `Best*` default. SecondBEST introduces a proper `soil` vertical grid, mirroring the existing ice grid. This plan is self-contained (no external data, no ocean dependency) and is a good autonomous-execution candidate.

**Tech Stack:** Python, numpy, scipy (sparse tridiagonal solve), sympl (`Stepper`, `get_constant`, `set_constant`, `initialize_numpy_arrays_with_properties`), pytest. Reference physics: `docs/Description_of_SecondBEST.ipynb` and `climt/_components/second_best/BEST-Description.pdf` (Pitman et al.).

## Global Constraints

- **Branch:** `feature/land-ocean-ice-components` (already created off `develop`).
- **No JAX / no accelerator path.** Match `BucketHydrology`'s existing style: plain-numpy `array_call` with `initialize_numpy_arrays_with_properties`. `SecondBEST` follows the same style.
- **Backward compatibility is non-negotiable for BucketHydrology:** the `num_layers=1` default MUST reproduce today's output bit-for-bit. The existing `tests/test_bucket_hydrology.py` and `tests/test_components.py::TestBucketHydrology` cached output are the guards.
- **`freezing_temperature_of_liquid_phase` is 273.0 K** in climt — always read via `get_constant`, never hardcode.
- **Verified constants** (read via `get_constant(name, units)`): `gravitational_acceleration`=9.80665 (m/s^2); `heat_capacity_of_dry_air_at_constant_pressure`=1004.64 (J/kg/degK); `latent_heat_of_condensation`=`latent_heat_of_vaporization`=2.5e6 (J/kg); `latent_heat_of_fusion`=333550.0 (J/kg); `gas_constant_of_dry_air`=287.0 (J/kg/degK); `density_of_liquid_water`=1000.0 (kg/m^3); `heat_capacity_of_liquid_water`=4185.5 (J/kg/degK). **`von_karman_constant` is NOT registered** — Task 4 registers it as 0.4.
- **`get_land_grid` raises `NotImplementedError` for 3-D land** today; SecondBEST adds a `soil` vertical grid (Task 3) rather than overloading the land grid.
- **Test env:** `~/miniconda3/envs/climt/bin/python -m pytest ...`.
- **Cached-output tests self-seed:** first run of a `ComponentBase*` subclass writes a `.cache` under `tests/cached_component_output/` and raises "Failed due to no cached output"; commit the cache, rerun → pass.
- **Registration required:** new components must be added to `climt/_components/__init__.py` and `climt/__init__.py`; new state quantities to `default_values` in `climt/_core/initialization.py`.

## File Structure

- Modify `climt/_components/bucket_hydrology/component.py` — add `num_layers`, deep water/thermal stores, fix the `0.15` clamp bug.
- Modify `climt/_core/initialization.py` — add `default_values` for `deep_soil_moisture_content`, `deep_soil_temperature`, `runoff_rate`, and the SecondBEST soil-profile quantities; add `get_soil_grid` + `soil` domain + `height_on_soil_interface_levels` grid quantity + `n_soil_interface_levels` param.
- Modify `climt/__init__.py` — register `von_karman_constant`; export `SecondBEST`.
- Create `climt/_components/second_best/component.py` — `SecondBEST` orchestrator (replaces the broken stub `__init__.py` import).
- Create `climt/_components/second_best/processes/__init__.py` — process base protocols + `Best*` defaults, split by responsibility:
  - `soil_properties.py`, `albedo.py`, `surface_layer.py`, `fluxes.py`, `subsurface.py`.
- Modify `climt/_components/second_best/__init__.py` — fix export (`from .component import SecondBEST`).
- Modify `climt/_components/__init__.py` — export `SecondBEST`.
- Create tests: `tests/test_bucket_two_layer.py`, `tests/test_best_processes.py`, `tests/test_second_best.py`; extend `tests/test_components.py`, `tests/test_conservation.py`.
- Modify `docs/Description_of_SecondBEST.ipynb` — add an "Extending SecondBEST" section (Task 10).

Build order: Bucket (Tasks 1–2) → SecondBEST soil grid (Task 3) → process objects (Tasks 4–8) → orchestrator + integration (Task 9) → docs (Task 10).

---

## Phase A — Two-layer BucketHydrology

### Task 1: Fix the `soil_moisture_max` clamp bug and add `num_layers` plumbing (still 1-layer, bit-for-bit)

**Files:**
- Modify: `climt/_components/bucket_hydrology/component.py`
- Test: `tests/test_bucket_two_layer.py`

**Interfaces:**
- Produces: `BucketHydrology(num_layers=1, soil_moisture_max=0.15, beta_parameter=0.75, specific_latent_heat_of_water=2260000, bulk_coefficient=0.0011, deep_soil_moisture_max=0.50, moisture_diffusion_timescale=None, deep_layer_thickness_ratio=10.0, deep_drainage_timescale=None)`. With `num_layers=1`, identical output to today; the soil-moisture clamp uses `self._smax` instead of the literal `0.15`.

- [ ] **Step 1: Write the failing test (clamp honours configured max)**

```python
# tests/test_bucket_two_layer.py
from datetime import timedelta
import numpy as np
from climt import BucketHydrology, get_default_state, get_grid


def _one_layer_state(smax):
    comp = BucketHydrology(soil_moisture_max=smax)
    state = get_default_state([comp], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["downwelling_shortwave_flux_in_air"].values[:] = 40.0
    state["downwelling_longwave_flux_in_air"].values[:] = 40.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 40.0
    state["upwelling_longwave_flux_in_air"].values[:] = 40.0
    state["soil_layer_thickness"].values[:] = 1.0
    state["surface_temperature"].values[:] = 300.0
    state["stratiform_precipitation_rate"].values[:] = 0.5
    state["convective_precipitation_rate"].values[:] = 0.5
    state["lwe_thickness_of_soil_moisture_content"].values[:] = 0.4
    return comp, state


def test_soil_moisture_clamped_to_configured_max():
    comp, state = _one_layer_state(smax=0.4)
    _, new = comp(state, timedelta(seconds=1))
    # heavy precip; moisture must be capped at the configured max (0.4), not 0.15
    assert np.all(new["lwe_thickness_of_soil_moisture_content"].values <= 0.4 + 1e-12)
    assert np.all(new["lwe_thickness_of_soil_moisture_content"].values > 0.15)
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_two_layer.py::test_soil_moisture_clamped_to_configured_max -v`
Expected: FAIL — moisture capped at 0.15 by the hardcoded clamp (`component.py:216`).

- [ ] **Step 3: Fix the clamp and add config plumbing**

In `bucket_hydrology/component.py` `__init__`, add the new kwargs and store them:

```python
    def __init__(
        self,
        num_layers=1,
        soil_moisture_max=0.15,
        beta_parameter=0.75,
        specific_latent_heat_of_water=2260000,
        bulk_coefficient=0.0011,
        deep_soil_moisture_max=0.50,
        moisture_diffusion_timescale=None,
        deep_layer_thickness_ratio=10.0,
        deep_drainage_timescale=None,
        **kwargs,
    ):
        if num_layers not in (1, 2):
            raise ValueError("num_layers must be 1 or 2")
        self._num_layers = num_layers
        self._smax = soil_moisture_max
        self._g = beta_parameter
        self._c = bulk_coefficient
        self._l = specific_latent_heat_of_water
        self._deep_smax = deep_soil_moisture_max
        self._tau_m = moisture_diffusion_timescale
        self._deep_ratio = deep_layer_thickness_ratio
        self._tau_drain = deep_drainage_timescale
        super(BucketHydrology, self).__init__(**kwargs)
```

Replace the clamp at the end of `array_call` (currently `np.minimum(..., 0.15)`):

```python
        new_state["lwe_thickness_of_soil_moisture_content"] = np.minimum(
            state["lwe_thickness_of_soil_moisture_content"]
            + (soil_moisture_tendency * timestep.total_seconds()),
            self._smax,
        )
```

- [ ] **Step 4: Run the test** → PASS.

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_two_layer.py::test_soil_moisture_clamped_to_configured_max -v`

- [ ] **Step 5: Backward-compat guard (bit-for-bit)**

Run the existing bucket tests — `num_layers` defaults to 1 and the default `soil_moisture_max` is still 0.15, so the fixed clamp is numerically identical for the default config.
Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_hydrology.py "tests/test_components.py::TestBucketHydrology" -v`
Expected: PASS unchanged. (If `TestBucketHydrology` doesn't exist, skip that one.)

- [ ] **Step 6: Commit**

```bash
git add climt/_components/bucket_hydrology/component.py tests/test_bucket_two_layer.py
git commit -m "fix(bucket): honour configured soil_moisture_max in clamp; add num_layers plumbing"
```

---

### Task 2: Two-layer water and thermal balance

**Files:**
- Modify: `climt/_components/bucket_hydrology/component.py`
- Modify: `climt/_core/initialization.py` (add `deep_soil_moisture_content`, `deep_soil_temperature`, `runoff_rate` defaults)
- Test: `tests/test_bucket_two_layer.py`, `tests/test_conservation.py` (`TestBucketTwoLayerWater`)

**Interfaces:**
- Consumes: Task 1 config.
- Produces: with `num_layers=2`, `BucketHydrology` reads/writes `deep_soil_moisture_content` (m) and `deep_soil_temperature` (degK) in addition to the shallow fields; emits `runoff_rate` (m s^-1) and `deep_soil_moisture_fraction` diagnostics. Shallow `lwe_thickness_of_soil_moisture_content` and `surface_temperature` keep their `["*"]` dims in both modes.

- [ ] **Step 1: Add state defaults**

In `climt/_core/initialization.py` `default_values`, near the land block (`~line 905`, next to `lwe_thickness_of_soil_moisture_content`), add:

```python
    "deep_soil_moisture_content": {
        "value": 0.25, "units": "m", "domain": "land_horizontal"},
    "deep_soil_temperature": {
        "value": 285.0, "units": "degK", "domain": "land_horizontal"},
    "runoff_rate": {
        "value": 0.0, "units": "m s^-1", "domain": "land_horizontal"},
```

- [ ] **Step 2: Write the failing two-layer conservation test (unit + conservation)**

```python
# tests/test_bucket_two_layer.py  (append)
def _two_layer_state():
    comp = BucketHydrology(num_layers=2, soil_moisture_max=0.15,
                           deep_soil_moisture_max=0.5,
                           moisture_diffusion_timescale=86400.0)
    state = get_default_state([comp], grid_state=get_grid(nx=1, ny=1, nz=10))
    state["downwelling_shortwave_flux_in_air"].values[:] = 200.0
    state["downwelling_longwave_flux_in_air"].values[:] = 300.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 50.0
    state["upwelling_longwave_flux_in_air"].values[:] = 300.0
    state["soil_layer_thickness"].values[:] = 1.0
    state["surface_temperature"].values[:] = 295.0
    state["deep_soil_temperature"].values[:] = 290.0
    state["stratiform_precipitation_rate"].values[:] = 0.001
    state["convective_precipitation_rate"].values[:] = 0.0
    state["lwe_thickness_of_soil_moisture_content"].values[:] = 0.08
    state["deep_soil_moisture_content"].values[:] = 0.20
    return comp, state


def test_two_layer_total_water_conserved():
    comp, state = _two_layer_state()
    dt = timedelta(seconds=100)
    diag, new = comp(state, dt)
    w0 = (state["lwe_thickness_of_soil_moisture_content"].values
          + state["deep_soil_moisture_content"].values)
    w1 = (new["lwe_thickness_of_soil_moisture_content"].values
          + new["deep_soil_moisture_content"].values)
    P = (state["stratiform_precipitation_rate"].values
         + state["convective_precipitation_rate"].values)
    E = diag["evaporation_rate"].values
    R = diag["runoff_rate"].values
    # no drainage set -> total water change = (P - E - R) * dt
    assert np.isclose((w1 - w0), (P - E - R) * dt.total_seconds(), rtol=0, atol=1e-9)
```

- [ ] **Step 3: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_two_layer.py::test_two_layer_total_water_conserved -v`
Expected: FAIL (`deep_soil_moisture_content` not consumed / `runoff_rate` missing / KeyError).

- [ ] **Step 4: Implement the two-layer path**

Add `deep_soil_moisture_content`, `deep_soil_temperature` to `input_properties` and `output_properties` (`{"dims": ["*"], "units": "m"}` / `"degK"`), and `runoff_rate` + `deep_soil_moisture_fraction` to `diagnostic_properties`. In `array_call`, branch on `self._num_layers`. Keep the existing single-layer code path untouched under `if self._num_layers == 1:`. Add the two-layer branch:

```python
        if self._num_layers == 2:
            dt = timestep.total_seconds()
            w_s = state["lwe_thickness_of_soil_moisture_content"]
            w_d = state["deep_soil_moisture_content"]
            S_s, S_d = self._smax, self._deep_smax
            tau_m = self._tau_m if self._tau_m is not None else 5 * 86400.0

            # beta keyed on shallow wetness (same limiter as 1-layer)
            beta = np.ones(w_s.shape)
            mask = w_s <= self._g * S_s
            beta[mask] = w_s[mask] / (self._g * S_s)
            evaporation_rate = beta * potential_evaporation  # from existing code above

            # shallow<->deep exchange toward equal relative saturation
            F_sd = (w_s / S_s - w_d / S_d) * (0.5 * (S_s + S_d)) / tau_m
            # optional deep drainage
            D = (w_d / self._tau_drain) if self._tau_drain is not None else 0.0

            w_s_new = w_s + (precipitation_rate - evaporation_rate - F_sd) * dt
            w_d_new = w_d + (F_sd - D) * dt

            # runoff from overflow of either layer
            runoff = np.zeros(w_s.shape)
            over_s = np.maximum(w_s_new - S_s, 0.0)
            over_d = np.maximum(w_d_new - S_d, 0.0)
            runoff += (over_s + over_d) / dt
            w_s_new = np.clip(w_s_new - over_s, 0.0, S_s)
            w_d_new = np.clip(w_d_new - over_d, 0.0, S_d)

            diagnostics["evaporation_rate"] = evaporation_rate
            diagnostics["runoff_rate"] = runoff
            diagnostics["deep_soil_moisture_fraction"] = w_d_new / S_d

            # thermal: shallow (skin) + deep store coupled by conduction
            k_soil = 2.0  # W/m/degK; documented constant
            dz_s = state["soil_layer_thickness"]
            dz_d = self._deep_ratio * dz_s
            C_s = (state["surface_material_density"] * dz_s
                   * state["heat_capacity_of_soil"])
            C_d = (state["surface_material_density"] * dz_d
                   * state["heat_capacity_of_soil"])
            T_s = state["surface_temperature"]
            T_d = state["deep_soil_temperature"]
            G_sd = k_soil * (T_s - T_d) / (0.5 * (dz_s + dz_d))
            # net_heat_flux computed above from radiation - SHF - LHF
            new_state["surface_temperature"] = T_s + (net_heat_flux - G_sd) / C_s * dt
            new_state["deep_soil_temperature"] = T_d + G_sd / C_d * dt
            new_state["lwe_thickness_of_soil_moisture_content"] = w_s_new
            new_state["deep_soil_moisture_content"] = w_d_new
            return diagnostics, new_state
```

Ensure `potential_evaporation`, `precipitation_rate`, `net_heat_flux`, the sensible/latent flux diagnostics, and `diagnostics["precipitation_rate"]` are computed **before** the `num_layers` branch (they are shared). Refactor so the shared prelude (wind speed, potential evaporation, precip, fluxes, net heat flux) runs first, then the 1-layer or 2-layer update.

- [ ] **Step 5: Run the conservation test** → PASS.

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_two_layer.py::test_two_layer_total_water_conserved -v`

- [ ] **Step 6: Property test (timescale separation)**

```python
# tests/test_bucket_two_layer.py  (append)
def test_two_layer_deep_store_is_slower_than_shallow():
    comp, state = _two_layer_state()
    dt = timedelta(seconds=3600)
    # drive with an oscillating precip; integrate and record variance
    shallow, deep = [], []
    import numpy as np
    for step in range(200):
        # subseasonal (fast) + seasonal (slow) precip signal
        fast = 0.001 * (1 + np.sin(2 * np.pi * step / 24))
        state["convective_precipitation_rate"].values[:] = max(fast, 0.0)
        diag, new = comp(state, dt)
        state.update(new); state.update(diag)
        shallow.append(float(new["lwe_thickness_of_soil_moisture_content"].values[0]))
        deep.append(float(new["deep_soil_moisture_content"].values[0]))
    # shallow store tracks the fast forcing more strongly than the deep store
    assert np.std(np.diff(shallow)) > np.std(np.diff(deep))
```

- [ ] **Step 7: Conservation-suite test class**

```python
# tests/test_conservation.py  (append)
class TestBucketTwoLayerWater(ConservationTestBase):
    def get_component_instance(self):
        return climt.BucketHydrology(num_layers=2,
                                     moisture_diffusion_timescale=86400.0)

    def modify_state(self, state):
        state["stratiform_precipitation_rate"].values[:] = 0.001
        state["lwe_thickness_of_soil_moisture_content"].values[:] = 0.08
        state["deep_soil_moisture_content"].values[:] = 0.2
        return state

    def get_quantity_amount(self, state):
        return (state["lwe_thickness_of_soil_moisture_content"].to_units("m").values
                + state["deep_soil_moisture_content"].to_units("m").values)

    def get_quantity_forcing(self, state):
        P = (state["convective_precipitation_rate"].to_units("m s^-1").values
             + state["stratiform_precipitation_rate"].to_units("m s^-1").values)
        E = state["evaporation_rate"].to_units("m s^-1").values
        R = state["runoff_rate"].to_units("m s^-1").values
        return P - E - R
```

- [ ] **Step 8: Add cached-output class for the 2-layer config (unit + invariance)**

```python
# tests/test_components.py  (append)
class TestBucketHydrologyTwoLayer(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.BucketHydrology(num_layers=2,
                                     moisture_diffusion_timescale=86400.0)
```

- [ ] **Step 9: Run all bucket tests, self-seed cache, commit**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_two_layer.py "tests/test_conservation.py::TestBucketTwoLayerWater" "tests/test_components.py -k BucketHydrology" -v`
Expected: cached seeds then passes; others PASS. Also confirm 1-layer tests still pass.

```bash
git add climt/_components/bucket_hydrology/component.py climt/_core/initialization.py \
        tests/test_bucket_two_layer.py tests/test_conservation.py tests/test_components.py \
        tests/cached_component_output/TestBucketHydrologyTwoLayer-*.cache
git commit -m "feat(bucket): two-layer (deep+shallow) water and thermal stores with num_layers=2"
```

---

## Phase B — SecondBEST soil grid

### Task 3: Register a `soil` vertical grid and soil-profile state quantities

**Files:**
- Modify: `climt/_core/initialization.py`
- Test: `tests/test_second_best.py` (grid smoke test)

**Interfaces:**
- Produces: `get_grid(..., n_soil_interface_levels=4)` adds `height_on_soil_interface_levels` to the grid state; `domain_shape_descriptor["soil"] = get_soil_grid`; the existing `soil_temperature` default is redomained to `soil_interface` (degK), and `soil_liquid_water_content` (m^3/m^3), `soil_ice_content` (m^3/m^3) are added on domain `soil_interface`; `surface_snow_thickness` already exists.

- [ ] **Step 1: Write the failing grid test**

```python
# tests/test_second_best.py
import numpy as np
from climt import get_grid


def test_soil_grid_present():
    grid = get_grid(nx=4, ny=2, nz=10, n_soil_interface_levels=4)
    assert "height_on_soil_interface_levels" in grid
    assert grid["height_on_soil_interface_levels"].shape[0] == 4
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_second_best.py::test_soil_grid_present -v`
Expected: FAIL (`n_soil_interface_levels` unknown / key absent).

- [ ] **Step 3: Add the soil grid (copy the ice-grid template)**

In `climt/_core/initialization.py`, next to `get_ice_grid` (~line 59), add:

```python
def get_soil_grid(grid_state, interface=False, horizontal=False):
    y, x = grid_state["latitude"].shape
    y_dim, x_dim = grid_state["latitude"].dims
    z = grid_state["height_on_soil_interface_levels"].shape[0]
    if horizontal:
        return (y, x), (y_dim, x_dim)
    if interface:
        return (z, y, x), ("soil_interface_levels", y_dim, x_dim)
    else:
        return (z - 1, y, x), ("soil_mid_levels", y_dim, x_dim)
```

Add to `domain_shape_descriptor` (~line 76): `"soil": get_soil_grid,`.

Add `n_soil_interface_levels=4` to the `get_grid(...)` signature (~line 436) and to `get_default_state(...)` (~line 1044, thread it through to `get_grid`). After the `n_ice_interface_levels` block (~line 538), add:

```python
    if n_soil_interface_levels is not None:
        return_state["height_on_soil_interface_levels"] = get_backend().create_quantity(
            np.linspace(0.0, 2.0, n_soil_interface_levels),
            "height_on_soil_interface_levels",
            "m",
            ("soil_interface_levels",),
        )
```

**Change** the existing `soil_temperature` default (currently at `initialization.py:907`, a scalar `land_horizontal` placeholder `{"value": 274.0, "units": "degK", "domain": "land_horizontal"}` that is defined but consumed nowhere — it was reserved for exactly this kind of land model) to a soil-profile quantity, and **add** the liquid/ice profile quantities alongside it:

```python
    "soil_temperature": {
        "value": 285.0, "units": "degK", "domain": "soil_interface"},
    "soil_liquid_water_content": {
        "value": 0.2, "units": "m^3/m^3", "domain": "soil_interface"},
    "soil_ice_content": {
        "value": 0.0, "units": "m^3/m^3", "domain": "soil_interface"},
```

`soil_temperature` is thus SecondBEST's prognostic profile temperature on the `soil_interface` domain. Verify nothing else references it as a scalar: `grep -rn '"soil_temperature"' climt tests` should show only `initialization.py` (which you are editing) — the two-layer bucket uses `surface_temperature`/`deep_soil_temperature`, not `soil_temperature`.

- [ ] **Step 4: Run the grid test** → PASS.

- [ ] **Step 5: Backward-compat guard**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_initialization.py -v`
Expected: PASS (soil grid is additive; `n_soil_interface_levels` defaults to a value, existing callers unaffected). If any test asserts an exact set of grid keys, update it to allow the new key.

- [ ] **Step 6: Commit**

```bash
git add climt/_core/initialization.py tests/test_second_best.py
git commit -m "feat(second-best): register soil vertical grid and soil-profile state quantities"
```

---

## Phase C — BEST process objects

Each process lives in `climt/_components/second_best/processes/<name>.py`. Define a tiny protocol per process (a base class with a documented `__call__`) plus the `Best*` default. All defaults implement the equations in `docs/Description_of_SecondBEST.ipynb`.

### Task 4: Process protocol + `BestSoilProperties`

**Files:**
- Create: `climt/_components/second_best/processes/__init__.py`
- Create: `climt/_components/second_best/processes/soil_properties.py`
- Modify: `climt/__init__.py` (register `von_karman_constant`)
- Test: `tests/test_best_processes.py`

**Interfaces:**
- Produces: `SoilProperties` base with `__call__(soil_type: str, area_type: str) -> dict`; `BestSoilProperties()` returning keys `colour`, `texture`, `porosity`, `field_capacity`, `wilting_point`, `B`, `K_H0`, `psi_0`.

- [ ] **Step 1: Register the von Karman constant**

In `climt/__init__.py`, next to `sympl.set_constant("top_of_model_pressure", 20.0, "Pa")` add:

```python
sympl.set_constant("von_karman_constant", 0.4, "dimensionless")
```

- [ ] **Step 2: Write the failing unit test (hand-computed notebook values)**

```python
# tests/test_best_processes.py
import numpy as np
from climt._components.second_best.processes.soil_properties import BestSoilProperties


def test_soil_properties_sand_and_clay():
    sp = BestSoilProperties()
    clay = sp("clay", "land")
    sand = sp("sand", "land")
    # clay: tex=0 -> porosity 0.6, X_FC 0.57, colour 0.2
    assert np.isclose(clay["texture"], 0.0)
    assert np.isclose(clay["porosity"], 0.6)
    assert np.isclose(clay["field_capacity"], 0.95 * 0.6)
    assert np.isclose(clay["colour"], 0.2)
    # sand: tex=9 -> porosity 0.33, X_FC (0.95-0.774)*0.33, colour 1.0
    assert np.isclose(sand["texture"], 9.0)
    assert np.isclose(sand["porosity"], 0.6 - 0.03 * 9)
    assert np.isclose(sand["field_capacity"], (0.95 - 0.086 * 9) * (0.6 - 0.03 * 9))
    assert np.isclose(sand["colour"], 1.0)
    # land_ice override on texture
    ice = sp("clay", "land_ice")
    assert np.isclose(ice["texture"], 0.07)
    assert np.isclose(ice["wilting_point"], 0.01)
```

- [ ] **Step 3: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_best_processes.py::test_soil_properties_sand_and_clay -v`
Expected: FAIL `ModuleNotFoundError`.

- [ ] **Step 4: Implement**

```python
# climt/_components/second_best/processes/__init__.py
"""Pluggable process objects for SecondBEST.

Each process is a plain class exposing a documented ``__call__``. The
``Best*`` defaults implement the equations in
``docs/Description_of_SecondBEST.ipynb`` (Pitman et al.). Swap a process by
passing an instance to ``SecondBEST(...)``.
"""
from .soil_properties import SoilProperties, BestSoilProperties
from .albedo import SurfaceAlbedo, BestSurfaceAlbedo
from .surface_layer import SurfaceLayer, BestSurfaceLayer
from .fluxes import SurfaceFluxes, BestSurfaceFluxes
from .subsurface import SubsurfaceTransport, BestSubsurfaceTransport

__all__ = [
    "SoilProperties", "BestSoilProperties",
    "SurfaceAlbedo", "BestSurfaceAlbedo",
    "SurfaceLayer", "BestSurfaceLayer",
    "SurfaceFluxes", "BestSurfaceFluxes",
    "SubsurfaceTransport", "BestSubsurfaceTransport",
]
```

```python
# climt/_components/second_best/processes/soil_properties.py
"""Soil property process (BEST Eqs 4.10-4.12)."""


class SoilProperties:
    """Contract: __call__(soil_type, area_type) -> dict of soil parameters."""
    def __call__(self, soil_type, area_type):
        raise NotImplementedError


class BestSoilProperties(SoilProperties):
    _COLOUR = {"clay": 0.2, "sand": 1.0}
    _TEXTURE = {"clay": 0.0, "sand": 9.0}
    _B = {"clay": 10.0, "sand": 4.0}
    _K_H0 = {"clay": 0.001, "sand": 0.1}

    def __call__(self, soil_type, area_type):
        colour = self._COLOUR[soil_type]
        texture = self._TEXTURE[soil_type]
        if area_type == "land_ice":
            texture = 0.07
        porosity = 0.6 - 0.03 * texture
        field_capacity = (0.95 - 0.086 * texture) * porosity
        wilting_point = 0.01 if area_type == "land_ice" else porosity - 0.03
        return {
            "colour": colour, "texture": texture, "porosity": porosity,
            "field_capacity": field_capacity, "wilting_point": wilting_point,
            "B": self._B[soil_type], "K_H0": self._K_H0[soil_type],
            "psi_0": -0.2,
        }
```

- [ ] **Step 5: Run the unit test** → PASS.

- [ ] **Step 6: Commit**

```bash
git add climt/_components/second_best/processes climt/__init__.py tests/test_best_processes.py
git commit -m "feat(second-best): process protocol + BestSoilProperties; register von_karman_constant"
```

---

### Task 5: `BestSurfaceAlbedo`

**Files:**
- Create: `climt/_components/second_best/processes/albedo.py`
- Test: `tests/test_best_processes.py`

**Interfaces:**
- Produces: `SurfaceAlbedo` base with `__call__(soil_props: dict, wetness: float, area_type: str) -> dict`; `BestSurfaceAlbedo()` returning `alpha_sw`, `alpha_lw` (BEST Eqs 5.5-5.8).

- [ ] **Step 1: Write the failing unit test (hand-computed)**

```python
# tests/test_best_processes.py  (append)
from climt._components.second_best.processes.albedo import BestSurfaceAlbedo
from climt._components.second_best.processes.soil_properties import BestSoilProperties


def test_albedo_bare_soil_and_land_ice():
    alb = BestSurfaceAlbedo()
    sp = BestSoilProperties()
    clay = sp("clay", "land")
    # wet clay (W=1): alpha_sw = 0.10 + 0.1*0.2 + 0.06*0 = 0.12; alpha_lw = 0.24
    out = alb(clay, wetness=1.0, area_type="land")
    assert np.isclose(out["alpha_sw"], 0.12)
    assert np.isclose(out["alpha_lw"], 0.24)
    # dry sand (W=0): 0.10 + 0.1*1.0 + 0.06*1 = 0.26
    sand = sp("sand", "land")
    out2 = alb(sand, wetness=0.0, area_type="land")
    assert np.isclose(out2["alpha_sw"], 0.26)
    # land ice (W=0): alpha_sw = 0.60 + 0.06 = 0.66; alpha_lw = alpha_sw/3
    out3 = alb(sp("clay", "land_ice"), wetness=0.0, area_type="land_ice")
    assert np.isclose(out3["alpha_sw"], 0.66)
    assert np.isclose(out3["alpha_lw"], 0.66 / 3.0)
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_best_processes.py::test_albedo_bare_soil_and_land_ice -v`
Expected: FAIL `ModuleNotFoundError`.

- [ ] **Step 3: Implement**

```python
# climt/_components/second_best/processes/albedo.py
"""Surface albedo process (BEST Eqs 5.5-5.8)."""


class SurfaceAlbedo:
    def __call__(self, soil_props, wetness, area_type):
        raise NotImplementedError


class BestSurfaceAlbedo(SurfaceAlbedo):
    def __call__(self, soil_props, wetness, area_type):
        if area_type == "land_ice":
            alpha_sw = 0.60 + 0.06 * (1.0 - wetness)
            alpha_lw = alpha_sw / 3.0
        else:
            alpha_sw = 0.10 + 0.1 * soil_props["colour"] + 0.06 * (1.0 - wetness)
            alpha_lw = 2.0 * alpha_sw
        return {"alpha_sw": alpha_sw, "alpha_lw": alpha_lw}
```

- [ ] **Step 4: Run the unit test** → PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_components/second_best/processes/albedo.py tests/test_best_processes.py
git commit -m "feat(second-best): BestSurfaceAlbedo (SW/LW bare soil and land ice)"
```

---

### Task 6: `BestSurfaceLayer` (drag / Richardson number)

**Files:**
- Create: `climt/_components/second_best/processes/surface_layer.py`
- Test: `tests/test_best_processes.py`

**Interfaces:**
- Produces: `SurfaceLayer` base with `__call__(z_mid, z0, wind_speed, T_surf, T_air, area_type) -> dict`; `BestSurfaceLayer()` returning `C_Dm`, `C_Dh`, `C_DN`, `Ri` (BEST Section 6).

- [ ] **Step 1: Write the failing unit test (neutral case is analytic)**

```python
# tests/test_best_processes.py  (append)
from climt._components.second_best.processes.surface_layer import BestSurfaceLayer


def test_surface_layer_neutral_matches_neutral_drag():
    sl = BestSurfaceLayer()
    z_mid, z0 = 10.0, 0.01
    out = sl(z_mid, z0, wind_speed=5.0, T_surf=300.0, T_air=300.0,
             area_type="land")   # neutral: T_surf == T_air -> Ri = 0
    kappa = 0.4
    C_DN = (kappa / (np.log(z_mid) - np.log(z0))) ** 2
    assert np.isclose(out["Ri"], 0.0)
    assert np.isclose(out["C_DN"], C_DN)
    # at Ri=0 the stable branch reduces to C_DN
    assert np.isclose(out["C_Dm"], C_DN)
    assert np.isclose(out["C_Dh"], C_DN)


def test_surface_layer_unstable_increases_drag():
    sl = BestSurfaceLayer()
    out = sl(10.0, 0.01, wind_speed=5.0, T_surf=290.0, T_air=300.0,
             area_type="land")   # T_surf < T_air is stable; flip for unstable
    unstable = sl(10.0, 0.01, wind_speed=5.0, T_surf=310.0, T_air=300.0,
                  area_type="land")
    assert unstable["Ri"] < 0.0
    assert unstable["C_Dh"] > unstable["C_DN"]   # unstable enhances exchange
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_best_processes.py -k surface_layer -v`
Expected: FAIL `ModuleNotFoundError`.

- [ ] **Step 3: Implement (BEST Section 6)**

```python
# climt/_components/second_best/processes/surface_layer.py
"""Surface-layer drag process (BEST Section 6)."""
import numpy as np
from sympl import get_constant


class SurfaceLayer:
    def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
        raise NotImplementedError


class BestSurfaceLayer(SurfaceLayer):
    def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
        kappa = get_constant("von_karman_constant", "dimensionless")
        g = get_constant("gravitational_acceleration", "m/s^2")
        U = np.maximum(wind_speed, 1e-3)
        C_DN = (kappa / (np.log(z_mid) - np.log(z0))) ** 2
        zeta = np.exp(-kappa / np.sqrt(C_DN))
        # Ri < 0 unstable (T_surf > T_air)
        Ri = -(g * z_mid / (T_surf * U * U)) * (T_surf - T_air)
        eps = 1.0 if area_type in ("sea", "sea_ice") else 0.01
        if np.ndim(Ri) == 0:
            Ri = float(Ri)
        unstable = Ri < 0.0
        # Momentum
        C_Dm = np.where(
            unstable,
            C_DN * (1 - 8 * Ri / (1 + 56.768 * C_DN * np.sqrt(np.abs(Ri) / zeta))),
            C_DN * ((1 - 4 * eps * Ri) ** 2) / (1 + 8 * (1 - eps) * Ri),
        )
        # Heat / scalars
        C_Dh = np.where(
            unstable,
            C_DN * (1 - 12 * Ri / (1 + 41.801 * C_DN * np.sqrt(np.abs(Ri) / zeta))),
            C_DN * ((1 - 4 * eps * Ri) / (1 + (6 - 4 * eps) * Ri)) ** 2,
        )
        return {"C_Dm": float(C_Dm), "C_Dh": float(C_Dh),
                "C_DN": float(C_DN), "Ri": float(Ri)}
```

- [ ] **Step 4: Run the unit tests** → PASS.

- [ ] **Step 5: Property test (drag positive & finite over Ri sweep)**

```python
# tests/test_best_processes.py  (append)
def test_surface_layer_drag_positive_over_ri_sweep():
    sl = BestSurfaceLayer()
    for dT in np.linspace(-20, 20, 41):
        out = sl(10.0, 0.01, 5.0, 300.0 + dT, 300.0, "land")
        assert out["C_Dm"] > 0.0 and np.isfinite(out["C_Dm"])
        assert out["C_Dh"] > 0.0 and np.isfinite(out["C_Dh"])
```

- [ ] **Step 6: Commit**

```bash
git add climt/_components/second_best/processes/surface_layer.py tests/test_best_processes.py
git commit -m "feat(second-best): BestSurfaceLayer stability-dependent drag (Richardson number)"
```

---

### Task 7: `BestSurfaceFluxes`

**Files:**
- Create: `climt/_components/second_best/processes/fluxes.py`
- Test: `tests/test_best_processes.py`

**Interfaces:**
- Produces: `SurfaceFluxes` base with `__call__(drag, atmos, soil, soil_props, timestep) -> dict`; `BestSurfaceFluxes()` returning `sensible_heat_flux`, `latent_heat_flux`, `momentum_flux`, `evaporation` (BEST Section 8). `drag` is the `BestSurfaceLayer` output dict; `atmos` is `{air_density, wind_speed, air_temperature, air_specific_humidity, u, v}`; `soil` is `{surface_temperature, saturation_specific_humidity, W_Lu, W_Fu}`.

- [ ] **Step 1: Write the failing unit test (hand-computed SHF, neutral)**

```python
# tests/test_best_processes.py  (append)
from climt._components.second_best.processes.fluxes import BestSurfaceFluxes


def test_sensible_heat_flux_bulk_formula():
    from sympl import get_constant
    fx = BestSurfaceFluxes()
    C_Dh = 0.003354
    drag = {"C_Dm": C_Dh, "C_Dh": C_Dh, "C_DN": C_Dh, "Ri": 0.0}
    atmos = {"air_density": 1.2, "wind_speed": 5.0, "air_temperature": 295.0,
             "air_specific_humidity": 0.01, "u": 5.0, "v": 0.0}
    soil = {"surface_temperature": 300.0, "saturation_specific_humidity": 0.02,
            "W_Lu": 1.0, "W_Fu": 0.0}
    out = fx(drag, atmos, soil, soil_props={"porosity": 0.6}, timestep=100.0)
    Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK")
    expected_shf = 1.2 * Cpd * 5.0 * C_Dh * (300.0 - 295.0)
    assert np.isclose(out["sensible_heat_flux"], expected_shf, rtol=1e-6)
    assert out["latent_heat_flux"] > 0.0   # surface wetter than air
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_best_processes.py -k surface_flux -v`
Expected: FAIL `ModuleNotFoundError`.

- [ ] **Step 3: Implement (BEST Section 8)**

```python
# climt/_components/second_best/processes/fluxes.py
"""Surface flux process (BEST Section 8): bulk SHF/LHF with beta wetness."""
import numpy as np
from sympl import get_constant


class SurfaceFluxes:
    def __call__(self, drag, atmos, soil, soil_props, timestep):
        raise NotImplementedError


class BestSurfaceFluxes(SurfaceFluxes):
    def __call__(self, drag, atmos, soil, soil_props, timestep):
        Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK")
        Lv = get_constant("latent_heat_of_vaporization", "J/kg")
        Lf = get_constant("latent_heat_of_fusion", "J/kg")
        Li = Lv + Lf  # sublimation
        rho = atmos["air_density"]
        U = atmos["wind_speed"]
        C_Dh = drag["C_Dh"]
        dT = soil["surface_temperature"] - atmos["air_temperature"]
        shf = rho * Cpd * U * C_Dh * dT

        # wetness factor beta_u (soil): frozen sublimation term + exfiltration ratio
        W_Lu = soil["W_Lu"]
        W_Fu = soil["W_Fu"]
        c_u = C_Dh * U                         # soil conductance
        dq = soil["saturation_specific_humidity"] - atmos["air_specific_humidity"]
        E_usP = rho * c_u * dq                 # potential evaporation rate
        # exfiltration limit E_usmax (Eq. in Section 8) — clamp Theta to [0.01, 1]
        B = 4.0                                # placeholder overwritten below if given
        beta_u = np.clip(W_Lu, 0.0, 1.0)       # simplified when soil_props lacks B/K_H0
        if "B" in soil_props:
            B = soil_props["B"]
            K_H0 = soil_props["K_H0"]
            Theta = np.clip((W_Lu - 0.01) / max(1.0 - W_Fu, 1e-6), 1e-3, 1.0)
            rho_w = get_constant("density_of_liquid_water", "kg/m^3")
            Xv = soil_props["porosity"]
            psi0 = soil_props["psi_0"]
            K_HD = (-4 * K_H0 * B * psi0 * rho_w * Xv * (1 - W_Fu)) / (np.pi * timestep)
            E_usmax = K_HD * Theta ** (0.5 * B + 2) - K_H0 * Theta ** (2 * B + 3)
            frozen_term = (W_Fu * Lv / Li) if Li > 0 else 0.0
            ratio = np.clip(E_usmax / E_usP, 0.0, 1.0) if abs(E_usP) > 1e-12 else 0.0
            beta_u = np.clip(frozen_term + ratio, 0.0, 1.0)

        evaporation = beta_u * E_usP / rho     # m/s equivalent (per unit density)
        lhf = Lv * rho * evaporation
        u, v = atmos["u"], atmos["v"]
        momentum_flux = -rho * drag["C_Dm"] * U * np.array([u, v])
        return {"sensible_heat_flux": float(shf), "latent_heat_flux": float(lhf),
                "momentum_flux": momentum_flux, "evaporation": float(evaporation)}
```

- [ ] **Step 4: Run the unit test** → PASS.

- [ ] **Step 5: Property test (LHF zero when surface as dry as air, beta bounded)**

```python
# tests/test_best_processes.py  (append)
def test_latent_flux_zero_when_no_humidity_gradient():
    fx = BestSurfaceFluxes()
    drag = {"C_Dm": 0.003, "C_Dh": 0.003, "C_DN": 0.003, "Ri": 0.0}
    atmos = {"air_density": 1.2, "wind_speed": 5.0, "air_temperature": 300.0,
             "air_specific_humidity": 0.02, "u": 5.0, "v": 0.0}
    soil = {"surface_temperature": 300.0, "saturation_specific_humidity": 0.02,
            "W_Lu": 1.0, "W_Fu": 0.0}
    out = fx(drag, atmos, soil, {"porosity": 0.6}, 100.0)
    assert np.isclose(out["latent_heat_flux"], 0.0)
```

- [ ] **Step 6: Commit**

```bash
git add climt/_components/second_best/processes/fluxes.py tests/test_best_processes.py
git commit -m "feat(second-best): BestSurfaceFluxes (bulk SHF/LHF with beta exfiltration limit)"
```

---

### Task 8: `BestSubsurfaceTransport` (coupled heat + liquid + ice)

**Files:**
- Create: `climt/_components/second_best/processes/subsurface.py`
- Test: `tests/test_best_processes.py`

**Interfaces:**
- Produces: `SubsurfaceTransport` base with `__call__(profiles, surface_flux_bc, timestep, dz) -> dict`; `BestSubsurfaceTransport()`; `profiles = {"T": array, "X_w": array, "X_i": array}`, `surface_flux_bc` a float (net surface heat flux, W/m^2); returns updated `{"T", "X_w", "X_i"}`. Implicit tridiagonal for diffusion, explicit Euler for freeze/melt source `Gamma`, Neumann bottom BC.

- [ ] **Step 1: Write the failing unit test (freezing produces ice, releases latent heat)**

```python
# tests/test_best_processes.py  (append)
from climt._components.second_best.processes.subsurface import BestSubsurfaceTransport


def test_subsurface_freezing_creates_ice_and_warms_toward_freezing():
    st = BestSubsurfaceTransport()
    n = 6
    profiles = {
        "T": np.full(n, 270.0),        # below freezing (273)
        "X_w": np.full(n, 0.2),        # liquid water present
        "X_i": np.zeros(n),
    }
    out = st(profiles, surface_flux_bc=-20.0, timestep=600.0, dz=0.1)
    assert np.all(out["X_i"] >= 0.0)
    assert out["X_i"].sum() > 0.0          # some water froze
    assert out["X_w"].sum() < profiles["X_w"].sum()  # liquid decreased
    assert np.all(out["T"] <= 273.0 + 1e-6)          # cannot exceed freezing here
    assert not np.any(np.isnan(out["T"]))
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_best_processes.py -k subsurface -v`
Expected: FAIL `ModuleNotFoundError`.

- [ ] **Step 3: Implement (self-contained tridiagonal solve)**

```python
# climt/_components/second_best/processes/subsurface.py
"""Subsurface heat + moisture + ice transport (BEST conduction section).

Implicit tridiagonal for the diffusive (density) terms; explicit Euler for
the freeze/melt source Gamma. Neumann bottom BC. Node 0 = bottom, node n-1
= surface. This solver is self-contained (its own scipy sparse tridiagonal)
so SecondBEST does not depend on the ocean/ice plan's snow_ice_column core.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from sympl import get_constant


class SubsurfaceTransport:
    def __call__(self, profiles, surface_flux_bc, timestep, dz):
        raise NotImplementedError


class BestSubsurfaceTransport(SubsurfaceTransport):
    def __init__(self, thermal_conductivity=2.0, volumetric_heat_capacity=2.0e6):
        self._kappa = thermal_conductivity
        self._cv = volumetric_heat_capacity

    def __call__(self, profiles, surface_flux_bc, timestep, dz):
        T = np.array(profiles["T"], dtype=float)
        X_w = np.array(profiles["X_w"], dtype=float)
        X_i = np.array(profiles["X_i"], dtype=float)
        n = T.shape[0]
        dt = timestep
        Tf = get_constant("freezing_temperature_of_liquid_phase", "degK")
        Lf = get_constant("latent_heat_of_fusion", "J/kg")
        rho_w = get_constant("density_of_liquid_water", "kg/m^3")

        # --- implicit heat diffusion ---
        rr = self._kappa * dt / (self._cv * dz * dz)
        main = np.full(n, 1 + 2 * rr)
        off = np.full(n, -rr)
        A = sparse.diags([off[1:], main, off[:-1]], [-1, 0, 1], format="csc")
        # Neumann bottom (node 0): dT/dz = 0
        A[0, 0] = 1 + rr; A[0, 1] = -rr
        rhs = T.copy()
        # surface flux Neumann (node n-1): kappa dT/dz = -flux
        A[-1, -1] = 1 + rr; A[-1, -2] = -rr
        rhs[-1] = T[-1] + surface_flux_bc * dt / (self._cv * dz)
        T_diff = spsolve(A, rhs)

        # --- freeze/melt source Gamma (explicit) ---
        # potential ice production rate: drives T toward Tf, converting phase
        Gamma = (self._cv / Lf) * (Tf - T_diff) / dt   # kg/m^3/s of ice
        # limit by available water (freezing) or ice (melting)
        max_freeze = rho_w * X_w / dt
        max_melt = rho_w * X_i / dt
        Gamma = np.clip(Gamma, -max_melt, max_freeze)

        X_i_new = X_i + Gamma * dt / rho_w
        X_w_new = X_w - Gamma * dt / rho_w
        # latent heat release/absorption adjusts temperature back
        T_new = T_diff + Lf * Gamma * dt / self._cv
        T_new = np.minimum(T_new, Tf) if surface_flux_bc <= 0 else T_new
        return {"T": T_new, "X_w": np.clip(X_w_new, 0.0, None),
                "X_i": np.clip(X_i_new, 0.0, None)}
```

- [ ] **Step 4: Run the unit test** → PASS.

- [ ] **Step 5: Property test (water+ice mass conserved by phase change)**

```python
# tests/test_best_processes.py  (append)
def test_subsurface_conserves_total_water_mass():
    st = BestSubsurfaceTransport()
    n = 8
    profiles = {"T": np.full(n, 268.0), "X_w": np.full(n, 0.15),
                "X_i": np.full(n, 0.05)}
    total0 = profiles["X_w"].sum() + profiles["X_i"].sum()
    out = st(profiles, surface_flux_bc=-10.0, timestep=300.0, dz=0.1)
    total1 = out["X_w"].sum() + out["X_i"].sum()
    assert np.isclose(total0, total1, atol=1e-9)   # phase change moves mass, conserves it
```

- [ ] **Step 6: Commit**

```bash
git add climt/_components/second_best/processes/subsurface.py tests/test_best_processes.py
git commit -m "feat(second-best): BestSubsurfaceTransport coupled heat/liquid/ice with freeze-melt"
```

---

## Phase D — SecondBEST orchestrator

### Task 9: `SecondBEST` component wiring, registration, integration

**Files:**
- Create: `climt/_components/second_best/component.py`
- Modify: `climt/_components/second_best/__init__.py` (fix broken export)
- Modify: `climt/_components/__init__.py`, `climt/__init__.py` (export `SecondBEST`)
- Test: `tests/test_second_best.py`, `tests/test_components.py` (`TestSecondBEST`), `tests/test_conservation.py` (`TestSecondBESTEnergy`)

**Interfaces:**
- Consumes: all five `Best*` processes (Tasks 4-8), the soil grid (Task 3).
- Produces: `SecondBEST(soil_type="clay", num_soil_layers=3, minimum_wind_speed=1.0, soil_properties=None, albedo=None, surface_layer=None, fluxes=None, subsurface=None, **kwargs)`, a `Stepper`. Inputs: lowest-level `air_temperature`, `specific_humidity`, `northward_wind`, `eastward_wind`, `air_pressure`, SW/LW down/up fluxes, `area_type`, and soil-profile state (`soil_temperature`, `soil_liquid_water_content`, `soil_ice_content`, `surface_snow_thickness`, `surface_temperature`). Diagnostics: `surface_upward_sensible_heat_flux`, `surface_upward_latent_heat_flux`, `evaporation_rate`, `surface_albedo_for_direct_shortwave`, `surface_albedo_for_diffuse_shortwave`, `surface_drag_coefficient_for_heat`, `surface_drag_coefficient_for_momentum`.

- [ ] **Step 1: Write the failing substitution test (proves the interface)**

```python
# tests/test_second_best.py  (append)
import numpy as np
from datetime import timedelta
from climt import SecondBEST, get_default_state, get_grid


class _StubFluxes:
    called = False
    def __call__(self, drag, atmos, soil, soil_props, timestep):
        _StubFluxes.called = True
        return {"sensible_heat_flux": 12.0, "latent_heat_flux": 34.0,
                "momentum_flux": np.array([0.0, 0.0]), "evaporation": 0.0}


def test_second_best_uses_injected_process():
    comp = SecondBEST(fluxes=_StubFluxes())
    state = get_default_state([comp], grid_state=get_grid(
        nx=1, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    diag, new = comp(state, timedelta(seconds=100))
    assert _StubFluxes.called
    assert np.isclose(diag["surface_upward_sensible_heat_flux"].values[0], 12.0)
```

- [ ] **Step 2: Run to verify it fails**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_second_best.py::test_second_best_uses_injected_process -v`
Expected: FAIL `ImportError: cannot import name 'SecondBEST'` (the current `__init__.py` imports a non-existent `component`).

- [ ] **Step 3: Implement the orchestrator**

```python
# climt/_components/second_best/component.py
import numpy as np
from sympl import Stepper, get_constant, initialize_numpy_arrays_with_properties
from .processes import (BestSoilProperties, BestSurfaceAlbedo, BestSurfaceLayer,
                        BestSurfaceFluxes, BestSubsurfaceTransport)


class SecondBEST(Stepper):
    """Modular intermediate-complexity BEST land surface model.

    A thin orchestrator: every physical calculation is delegated to a
    swappable process object (see ``climt._components.second_best.processes``).
    Pass an instance to override any process; ``None`` selects the BEST default.
    """

    input_properties = {
        "air_temperature": {"dims": ["mid_levels", "*"], "units": "degK"},
        "specific_humidity": {"dims": ["mid_levels", "*"], "units": "kg/kg"},
        "northward_wind": {"dims": ["mid_levels", "*"], "units": "m s^-1"},
        "eastward_wind": {"dims": ["mid_levels", "*"], "units": "m s^-1"},
        "air_pressure": {"dims": ["mid_levels", "*"], "units": "Pa"},
        "downwelling_shortwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "downwelling_longwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "upwelling_shortwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "upwelling_longwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "area_type": {"dims": ["*"], "units": "dimensionless"},
        "surface_temperature": {"dims": ["*"], "units": "degK"},
        "surface_air_pressure": {"dims": ["*"], "units": "Pa"},
        "soil_temperature": {"dims": ["soil_interface_levels", "*"], "units": "degK"},
        "soil_liquid_water_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "soil_ice_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "surface_snow_thickness": {"dims": ["*"], "units": "m"},
        "height_on_soil_interface_levels": {"dims": ["soil_interface_levels", "*"], "units": "m"},
    }
    output_properties = {
        "surface_temperature": {"dims": ["*"], "units": "degK"},
        "soil_temperature": {"dims": ["soil_interface_levels", "*"], "units": "degK"},
        "soil_liquid_water_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "soil_ice_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "surface_snow_thickness": {"dims": ["*"], "units": "m"},
    }
    diagnostic_properties = {
        "surface_upward_sensible_heat_flux": {"dims": ["*"], "units": "W m^-2"},
        "surface_upward_latent_heat_flux": {"dims": ["*"], "units": "W m^-2"},
        "evaporation_rate": {"dims": ["*"], "units": "m s^-1"},
        "surface_albedo_for_direct_shortwave": {"dims": ["*"], "units": "dimensionless"},
        "surface_albedo_for_diffuse_shortwave": {"dims": ["*"], "units": "dimensionless"},
        "surface_drag_coefficient_for_heat": {"dims": ["*"], "units": "dimensionless"},
        "surface_drag_coefficient_for_momentum": {"dims": ["*"], "units": "dimensionless"},
    }

    def __init__(self, soil_type="clay", num_soil_layers=3, minimum_wind_speed=1.0,
                 soil_properties=None, albedo=None, surface_layer=None,
                 fluxes=None, subsurface=None, **kwargs):
        self._soil_type = soil_type
        self._num_soil_layers = num_soil_layers
        self._min_wind = minimum_wind_speed
        self._soil_props = soil_properties or BestSoilProperties()
        self._albedo = albedo or BestSurfaceAlbedo()
        self._surface_layer = surface_layer or BestSurfaceLayer()
        self._fluxes = fluxes or BestSurfaceFluxes()
        self._subsurface = subsurface or BestSubsurfaceTransport()
        super(SecondBEST, self).__init__(**kwargs)

    def array_call(self, state, timestep):
        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties)
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties)

        Rd = get_constant("gas_constant_of_dry_air", "J/kg/degK")
        ncol = state["area_type"].shape[0]
        for col in range(ncol):
            area = state["area_type"][col].astype(str)
            if area not in ("land", "land_ice"):
                continue
            props = self._soil_props(self._soil_type, area)

            u = state["eastward_wind"][0, col]
            v = state["northward_wind"][0, col]
            wind = max(np.sqrt(u * u + v * v), self._min_wind)
            T_air = state["air_temperature"][0, col]
            p = state["air_pressure"][0, col]
            rho = p / (Rd * T_air)
            z_mid = Rd * T_air / get_constant("gravitational_acceleration", "m/s^2")
            z0 = 0.01 if area == "land" else 0.001

            T_surf = state["surface_temperature"][col]
            drag = self._surface_layer(z_mid, z0, wind, T_surf, T_air, area)

            X_w = state["soil_liquid_water_content"][:, col]
            W_Lu = X_w[-1] / props["porosity"]
            albedo = self._albedo(props, W_Lu, area)

            atmos = {"air_density": rho, "wind_speed": wind,
                     "air_temperature": T_air,
                     "air_specific_humidity": state["specific_humidity"][0, col],
                     "u": u, "v": v}
            q_sat = _saturation_specific_humidity(T_surf, p)
            soil = {"surface_temperature": T_surf,
                    "saturation_specific_humidity": q_sat,
                    "W_Lu": W_Lu, "W_Fu": state["soil_ice_content"][-1, col]
                    / props["porosity"]}
            flux = self._fluxes(drag, atmos, soil, props, timestep.total_seconds())

            net = (state["downwelling_shortwave_flux_in_air"][col, 0]
                   + state["downwelling_longwave_flux_in_air"][col, 0]
                   - state["upwelling_shortwave_flux_in_air"][col, 0]
                   - state["upwelling_longwave_flux_in_air"][col, 0]
                   - flux["sensible_heat_flux"] - flux["latent_heat_flux"])

            z = state["height_on_soil_interface_levels"][:, col]
            dz = float(abs(z[1] - z[0])) if z.shape[0] > 1 else 0.5
            new_prof = self._subsurface(
                {"T": state["soil_temperature"][:, col],
                 "X_w": X_w, "X_i": state["soil_ice_content"][:, col]},
                surface_flux_bc=net, timestep=timestep.total_seconds(), dz=dz)

            outputs["soil_temperature"][:, col] = new_prof["T"]
            outputs["soil_liquid_water_content"][:, col] = new_prof["X_w"]
            outputs["soil_ice_content"][:, col] = new_prof["X_i"]
            outputs["surface_temperature"][col] = new_prof["T"][-1]
            outputs["surface_snow_thickness"][col] = state["surface_snow_thickness"][col]

            diagnostics["surface_upward_sensible_heat_flux"][col] = flux["sensible_heat_flux"]
            diagnostics["surface_upward_latent_heat_flux"][col] = flux["latent_heat_flux"]
            diagnostics["evaporation_rate"][col] = flux["evaporation"]
            diagnostics["surface_albedo_for_direct_shortwave"][col] = albedo["alpha_sw"]
            diagnostics["surface_albedo_for_diffuse_shortwave"][col] = albedo["alpha_sw"]
            diagnostics["surface_drag_coefficient_for_heat"][col] = drag["C_Dh"]
            diagnostics["surface_drag_coefficient_for_momentum"][col] = drag["C_Dm"]
        return diagnostics, outputs


def _saturation_specific_humidity(T, p):
    # Bolton (1980) saturation vapour pressure over water, then specific humidity
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / (p - 0.378 * es)
```

Fix the stub package export:

```python
# climt/_components/second_best/__init__.py
from .component import SecondBEST

__all__ = (SecondBEST,)
```

- [ ] **Step 4: Register `SecondBEST`** in `_components/__init__.py` (`from .second_best import SecondBEST` + add to `__all__`) and `climt/__init__.py` (import block + `__all__`).

- [ ] **Step 5: Run the substitution test** → PASS.

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_second_best.py::test_second_best_uses_injected_process -v`

- [ ] **Step 6: Integration test (equilibrium under fixed forcing)**

```python
# tests/test_second_best.py  (append)
def test_second_best_reaches_reasonable_state():
    comp = SecondBEST()
    state = get_default_state([comp], grid_state=get_grid(
        nx=2, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    state["downwelling_shortwave_flux_in_air"].values[:] = 300.0
    state["downwelling_longwave_flux_in_air"].values[:] = 300.0
    state["upwelling_shortwave_flux_in_air"].values[:] = 60.0
    state["upwelling_longwave_flux_in_air"].values[:] = 350.0
    for _ in range(50):
        diag, new = comp(state, timedelta(seconds=600))
        state.update(new); state.update(diag)
    T = state["surface_temperature"].values
    assert np.all(np.isfinite(T))
    assert np.all((T > 200.0) & (T < 350.0))
```

- [ ] **Step 7: Cached-output class + conservation class**

```python
# tests/test_components.py  (append)
class TestSecondBEST(ComponentBaseColumn, ComponentBase3D):
    def get_component_instance(self):
        return climt.SecondBEST()

    def get_1d_input_state(self, component=None):
        state = super(TestSecondBEST, self).get_1d_input_state(component)
        state["area_type"].values[:] = "land"
        return state

    def get_3d_input_state(self, component=None):
        state = super(TestSecondBEST, self).get_3d_input_state(component)
        state["area_type"].values[:] = "land"
        return state
```

```python
# tests/test_conservation.py  (append)
class TestSecondBESTEnergy(SurfaceEnergyConservation):
    def get_component_instance(self):
        return climt.SecondBEST()

    def modify_state(self, state):
        state["area_type"].values[:] = "land"
        return state

    def get_quantity_amount(self, state):
        cv = 2.0e6  # matches BestSubsurfaceTransport default volumetric heat capacity
        z = state["height_on_soil_interface_levels"].to_units("m").values
        dz = abs(z[1] - z[0]) if z.shape[0] > 1 else 0.5
        T = state["soil_temperature"].to_units("degK").values
        return cv * dz * T.sum(axis=0)
```

If exact energy closure is defeated by the freeze/melt latent bookkeeping, assert closure with `atol` scaled to `cv*dz` and document that phase-change energy leaves the sensible store (same caveat as the ice components).

- [ ] **Step 8: Run everything, self-seed caches, commit**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_second_best.py tests/test_best_processes.py "tests/test_components.py::TestSecondBEST" "tests/test_conservation.py::TestSecondBESTEnergy" -v`
Expected: cached seeds then passes; others PASS.

```bash
git add climt/_components/second_best climt/_components/__init__.py climt/__init__.py \
        tests/test_second_best.py tests/test_components.py tests/test_conservation.py \
        tests/cached_component_output/TestSecondBEST-*.cache
git commit -m "feat(second-best): SecondBEST orchestrator composing pluggable BEST processes"
```

---

## Phase E — Documentation

### Task 10: "Extending SecondBEST" documentation

**Files:**
- Modify: `docs/Description_of_SecondBEST.ipynb`
- Test: manual (notebook executes) + a doctest-style example test in `tests/test_second_best.py`

**Interfaces:** none (docs only).

- [ ] **Step 1: Add a worked-example test that doubles as documentation**

```python
# tests/test_second_best.py  (append)
def test_extending_second_best_swap_example():
    """Documents how to swap a process — a constant-drag SurfaceLayer."""
    from climt._components.second_best.processes import SurfaceLayer

    class ConstantDrag(SurfaceLayer):
        def __call__(self, z_mid, z0, wind_speed, T_surf, T_air, area_type):
            return {"C_Dm": 0.002, "C_Dh": 0.002, "C_DN": 0.002, "Ri": 0.0}

    comp = SecondBEST(surface_layer=ConstantDrag())
    state = get_default_state([comp], grid_state=get_grid(
        nx=1, ny=1, nz=10, n_soil_interface_levels=4))
    state["area_type"].values[:] = "land"
    diag, _ = comp(state, timedelta(seconds=100))
    assert np.isclose(diag["surface_drag_coefficient_for_heat"].values[0], 0.002)
```

- [ ] **Step 2: Run it** → PASS.

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_second_best.py::test_extending_second_best_swap_example -v`

- [ ] **Step 3: Add the notebook section**

Append markdown + code cells to `docs/Description_of_SecondBEST.ipynb`:
- A markdown cell "## Extending SecondBEST" listing the five process contracts (`SoilProperties`, `SurfaceAlbedo`, `SurfaceLayer`, `SurfaceFluxes`, `SubsurfaceTransport`) with their `__call__` signatures.
- A code cell reproducing the `ConstantDrag` swap from Step 1.
- A markdown cell explaining how to add a *new* process (e.g. a vegetation/transpiration term): subclass the relevant contract, extend the orchestrator's `array_call` wiring, and pass the instance to `SecondBEST(...)`.

- [ ] **Step 4: Verify the notebook runs**

Run: `~/miniconda3/envs/climt/bin/jupyter nbconvert --to notebook --execute docs/Description_of_SecondBEST.ipynb --output /tmp/_sb_check.ipynb`
Expected: executes without error. (If `nbconvert`/`jupyter` is absent, `pip install nbconvert` in the env; this is a docs check only.)

- [ ] **Step 5: Commit**

```bash
git add docs/Description_of_SecondBEST.ipynb tests/test_second_best.py
git commit -m "docs(second-best): Extending SecondBEST — process contracts and swap/add examples"
```

---

## Final integration check

- [ ] **Step 1: Full suite**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/ -q`
Expected: all pass (pre-existing Cython-dependent skips for RRTMG/Emanuel/SimplePhysics if not compiled — unrelated; confirm the skip count matches a clean `develop`).

- [ ] **Step 2: Import + construction smoke test**

Run: `~/miniconda3/envs/climt/bin/python -c "from climt import BucketHydrology, SecondBEST; BucketHydrology(num_layers=2); SecondBEST(); print('ok')"`
Expected: `ok`.

- [ ] **Step 3: Confirm backward compatibility**

Run: `~/miniconda3/envs/climt/bin/python -m pytest tests/test_bucket_hydrology.py "tests/test_components.py::TestBucketHydrology" -v`
Expected: PASS (1-layer default unchanged).

---

## Self-review notes (coverage map)

| Spec section | Task(s) |
|---|---|
| Part 1 §1.1 config + `0.15` clamp bug | Task 1 |
| Part 1 §1.2 state (scalar deep fields — see deviation) | Tasks 1, 2 |
| Part 1 §1.3 water balance (2-layer) | Task 2 |
| Part 1 §1.4 thermal balance (2-layer) | Task 2 |
| Part 1 §1.5 diagnostics + conservation/timescale tests | Task 2 |
| Part 2 §2.1 orchestrator | Task 9 |
| Part 2 §2.2 process interfaces + Best defaults | Tasks 4–8 |
| Part 2 §2.3 component surface (inputs/outputs/diagnostics) | Tasks 3, 9 |
| Part 2 §2.4 "Extending SecondBEST" | Task 10 |
| Part 2 §2.5 unit/substitution/conservation/integration tests | Tasks 4–9 |

**Deliberate deviations from the spec (flagged):**
1. **Two-layer bucket uses scalar `deep_soil_moisture_content`/`deep_soil_temperature` fields, not a `soil_levels` grid dim** (spec §1.2). `get_land_grid` raises for 3-D land; scalar fields keep the existing `lwe_thickness_of_soil_moisture_content` at `["*"]` in both modes, preserving the bit-for-bit 1-layer guarantee with far less core-grid risk. Physics is identical.
2. **`BestSubsurfaceTransport` is self-contained** (its own sparse tridiagonal) rather than sharing the ocean/ice plan's `_core/snow_ice_column.py`, so this plan runs independently. Unifying the two solvers into one `_core` module is a reasonable follow-up.

**Deferred (spec open questions):** exact `moisture_diffusion_timescale` (default 5 days) and `k_soil` (2.0 W/m/degK) values are exposed as constants/config and can be tuned; the `soil` grid uses a simple `linspace(0, 2m)` default profile that a later task can refine.
