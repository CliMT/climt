# General Condensibles Design

**Date:** 2026-05-13  
**Branch:** feature/cork-radiation (will likely need its own branch)  
**Goal:** Make `EmanuelConvectionPython` (and future moist schemes) work for arbitrary condensibles — initially H₂O (Earth), CH₄ (Titan), CO₂ (Mars) — by pulling condensible physics out of hardcoded constants into a shared, profile-driven layer.

---

## Problem

`pure_python_v3.py` hardcodes H₂O properties throughout:

| Category | What's hardcoded | Where |
|---|---|---|
| Condensible physical constants | `RV`, `CPV`, `CL`, `LV0`, `ROWL` | `EmanuelParams` / `__init__` |
| Dry air constants | `CPD`, `RD`, `G` | `EmanuelParams` (already in TOMLs, not yet wired) |
| Saturation formula | Magnus/Bolton constants | inline in `_tlift_functional_np`, `bolton_q_sat` |
| Phase transition | `273.0` / `273.15` | `_tlift_functional_np`, `_convect_functional_np` |
| LCL formula | empirical constants `1669.0`, `122.0` | `_convect_functional_np` line 361 |
| Reference pressure | `1000.0` hPa | potential temperature calc |

The saturation formula is inside `@njit` kernels, so Python callables cannot be passed — dispatch must happen via an integer tag.

---

## Design

### Approach: `CondensibleParams` + shared `condensibles.py` core module

All condensible physics lives in `climt/_core/condensibles.py`, importable by any future scheme (grid-scale condensation, clouds, etc.). The atmospheric profile TOML is the single source of truth for physical values. `atmospheric_properties.py` is extended minimally to track the active condensible species string. Integer dispatch in `@njit` kernels selects the saturation formula per species.

---

## Section 1 — TOML changes

Each profile TOML gains two new sections.

**`[condensible]`** — string identifier, parsed separately (not a sympl constant):
```toml
[condensible]
species = "h2o"   # or "ch4", "co2"
```

**`[condensible_properties]`** — sympl-style `{value, units}` pairs, loaded via `set_constant()` / `get_constant()`:
```toml
[condensible_properties]
gas_constant_of_condensible_vapor         = { value = 461.5,   units = "J/kg/K" }
heat_capacity_of_condensible_vapor        = { value = 1870.0,  units = "J/kg/K" }
heat_capacity_of_condensible_liquid       = { value = 2500.0,  units = "J/kg/K" }
latent_heat_of_condensible_vaporization   = { value = 2.501e6, units = "J/kg" }
density_of_condensible_liquid             = { value = 1000.0,  units = "kg/m^3" }
freezing_point_of_condensible            = { value = 273.15,  units = "degK" }
```

Values shown are for H₂O (earth.toml). CH₄ (titan.toml) and CO₂ (mars.toml) get their respective values. The `[condensible_properties]` keys flow through the existing `_parse_toml` → `set_constant()` path unchanged.

---

## Section 2 — `atmospheric_properties.py` changes

Three additions, all following existing module patterns.

**Module-level state** (analogous to `_active_profile`):
```python
_active_condensible = {"species": "h2o"}  # default
```

**`load_atmospheric_properties()`** — after parsing TOML constants, also reads:
```python
_active_condensible["species"] = data.get("condensible", {}).get("species", "h2o")
```
The `[condensible_properties]` constants flow through `_parse_toml` automatically since they use the standard `{value, units}` format.

**Two new public functions:**
```python
def get_active_condensible() -> str:
    """Return the species string of the currently loaded condensible (e.g. 'ch4')."""
    return _active_condensible["species"]
```
`reset_atmospheric_properties()` is updated to also restore `_active_condensible["species"] = "h2o"` on reset.

`_parse_toml` and `_snapshot_constants` are unchanged — the string `species` field is parsed separately from the constants dict and never enters the sympl system.

---

## Section 3 — `climt/_core/condensibles.py` (new module)

Shared kernel-level library for condensible physics. Imported by any component needing saturation physics.

### Species registry
```python
SPECIES_ID = {"h2o": 0, "ch4": 1, "co2": 2}
```

### `CondensibleParams` NamedTuple
Numba-compatible (all scalar fields), passed into all kernels:
```python
class CondensibleParams(NamedTuple):
    species_id: int
    RV: float       # gas constant of condensible vapor (J/kg/K)
    CPV: float      # heat capacity of condensible vapor (J/kg/K)
    CL: float       # heat capacity of condensible liquid (J/kg/K)
    LV0: float      # latent heat of vaporization at reference T (J/kg)
    ROWL: float     # density of condensible liquid (kg/m^3)
    T_freeze: float # phase transition temperature (K)
```

### `get_condensible_params()` — factory
Called once per component `__init__`. Reads from `get_active_condensible()` and `get_constant()`:
```python
def get_condensible_params() -> CondensibleParams:
    species = get_active_condensible()
    return CondensibleParams(
        species_id=SPECIES_ID[species],
        RV=get_constant("gas_constant_of_condensible_vapor", "J/kg/K"),
        CPV=get_constant("heat_capacity_of_condensible_vapor", "J/kg/K"),
        CL=get_constant("heat_capacity_of_condensible_liquid", "J/kg/K"),
        LV0=get_constant("latent_heat_of_condensible_vaporization", "J/kg"),
        ROWL=get_constant("density_of_condensible_liquid", "kg/m^3"),
        T_freeze=get_constant("freezing_point_of_condensible", "degK"),
    )
```

### `@njit _sat_vap_pressure(T, species_id)` → `es` (Pa)
Scalar function, dispatches on `species_id`:
- `0` (H₂O): existing Magnus formula (above/below `T_freeze`), extracted from `_tlift_functional_np`
- `1` (CH₄): Clausius-Clapeyron with CH₄ reference constants
- `2` (CO₂): Clausius-Clapeyron with CO₂ reference constants

The `if/elif` branch on `species_id` is on a scalar constant for the entire run — CPU branch predictor achieves 100% accuracy; Numba may constant-fold with `numba.literally()` in future.

### `@njit _lv(T, cond)` → latent heat (J/kg)
Species-general given the right constants (no dispatch needed):
```python
return cond.LV0 - (cond.CL - cond.CPV) * (T - cond.T_freeze)
```

### `@njit _lcl_pressure(P_nk, RH, T_nk, cond)` → PLCL (hPa)
Receives `cond` (not just `species_id`) so Clausius-Clapeyron branches have access to `RV` and `LV0`:
- `0` (H₂O): existing empirical formula `CHI = T/(1669.0 - 122.0*RH - T)`, `PLCL = P * RH**CHI`
- `1`, `2`: generalized Clausius-Clapeyron form `PLCL = P * exp((cond.RV/cond.LV0) * T * log(RH))`

### `compute_qs(T, P, cond, RD)` → saturation specific humidity array
Replaces `bolton_q_sat` in `array_call`. Parallelized across columns, same pattern as `_numpy_vectorized_convect`:
```python
@jit_compile(backend=np, parallel=True)
def compute_qs(T, P, cond, RD):
    nlev, ncol = T.shape
    qs = np.zeros_like(T)
    EPS = RD / cond.RV
    for i in prange(ncol):
        for k in range(nlev):
            es = _sat_vap_pressure(T[k, i], cond.species_id)
            qs[k, i] = EPS * es / (P[k, i] - (1 - EPS) * es)
    return qs
```
`RD` is passed explicitly (not taken from either params struct) so `compute_qs` remains self-contained and reusable.

---

## Section 4 — `pure_python_v3.py` changes

### `EmanuelParams`
Drops: `RV`, `CPV`, `CL`, `LV0`, `ROWL` (moved to `CondensibleParams`).  
Keeps: `CPD`, `RD`, `G` (dry air / planetary properties, independent of condensible).

### `__init__`
- Adds: `self._cond = get_condensible_params()`
- Drops hardcoded values for the five removed fields
- `self.RD` and `self.RV` instance attributes removed (now on `self._params.RD` and `self._cond.RV`)

### `array_call`
Replaces:
```python
qs = bolton_q_sat(t, p * 100, self.RD, self.RV)
```
with:
```python
qs = compute_qs(t, p * 100, self._cond, self._params.RD)
```

### `_tlift_functional_np`
- Gains `cond` argument
- Inline Magnus block (lines 258–263) → `_sat_vap_pressure(TG, cond.species_id)`
- `273.15` reference → `cond.T_freeze`
- `CPVMCL`, `EPS` now computed from `cond` fields: `cond.CL - cond.CPV`, `params.RD / cond.RV`

### `_convect_functional_np`
- Gains `cond` argument, passes it to `_tlift_functional_np`
- Rain/snow phase check `T[i] > 273.0` (lines 602–603) → `T[i] > cond.T_freeze`
- LCL lines 361–362 → `_lcl_pressure(P[NK], RH, T[NK], cond)`
- `TPRIME` (line 659): `params.LV0 * QPRIME / params.CPD` → `cond.LV0 * QPRIME / params.CPD`

### `_numpy_vectorized_convect`
- Gains `cond` argument, passes to `_convect_functional_np`

### Sympl interface
`input_properties`, `tendency_properties`, `diagnostic_properties` are unchanged. Units like `kg/kg` for specific humidity remain valid for any condensible mass mixing ratio.

---

## What is NOT changed

- `bolton_q_sat` in `util.py` — left in place; used by other components that haven't been generalized yet
- `35.0` K minimum parcel temperature guard — becomes condensible-specific in a future pass

## Dead code

`TH` (potential temperature array, lines 314–319 of `_convect_functional_np`) is computed but never read — a carry-over from the original Fortran translation. The block (`TH`, its loop, and the `RDCP` local variable) is **commented out** rather than deleted. This also eliminates the Earth-specific `1000.0` hPa reference pressure concern for now. Other components that need `reference_pressure` will add it to the TOML as `reference_pressure = {value = ..., units = "hPa"}` in a future pass.

---

## Adding a new condensible in future

1. Add entry to `SPECIES_ID` in `condensibles.py`
2. Add `elif` branch to `_sat_vap_pressure` and `_lcl_pressure`
3. Add a TOML with `[condensible]` + `[condensible_properties]` sections
4. No changes needed to any component
