# General Condensibles Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generalise `EmanuelConvectionPython` (and the climt moist-physics infrastructure) to support arbitrary condensibles — initially H₂O (Earth), CH₄ (Titan), CO₂ (Mars) — driven by atmospheric profile TOMLs, with a shared `condensibles.py` core module that any future scheme can import.

**Architecture:** A new `climt/_core/condensibles.py` registers H₂O constants as sympl defaults at import time, exposes a `CondensibleParams` NamedTuple and `get_condensible_params()` factory, and provides `@njit` kernel helpers (`_sat_vap_pressure`, `_lv`, `_lcl_pressure`) with integer dispatch on `species_id`. Each atmospheric profile TOML gains `[condensible]` (species string) and `[condensible_properties]` (sympl constants) sections. `EmanuelParams` loses the five condensible fields; `_tlift_functional_np` and `_convect_functional_np` gain a `cond: CondensibleParams` argument. The dead `TH` potential-temperature block (lines 314–319 in `_convect_functional_np`) is commented out, eliminating the Earth-specific 1000 hPa reference pressure.

**Tech Stack:** Python 3.10+, NumPy, Numba (optional, falls back to Python via `@njit` shim), sympl constants, TOML (tomllib / tomli).

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `climt/_core/condensibles.py` | `CondensibleParams`, `SPECIES_ID`, H₂O defaults, `get_condensible_params()`, `_sat_vap_pressure`, `_lv`, `_lcl_pressure`, `compute_qs` |
| Modify | `climt/_core/atmospheric_properties.py` | `_active_condensible` state, `get_active_condensible()`, update `load_`/`reset_atmospheric_properties` |
| Modify | `climt/_core/__init__.py` | Export `get_active_condensible`, `get_condensible_params`, `CondensibleParams` |
| Modify | `climt/_data/atmospheric_properties/earth.toml` | Add `[condensible]` + `[condensible_properties]` for H₂O |
| Modify | `climt/_data/atmospheric_properties/titan.toml` | Add `[condensible]` + `[condensible_properties]` for CH₄ |
| Modify | `climt/_data/atmospheric_properties/mars.toml` | Add `[condensible]` + `[condensible_properties]` for CO₂ |
| Modify | `climt/_components/emanuel/pure_python_v3.py` | Remove condensible fields from `EmanuelParams`, add `self._cond`, thread `cond` into all kernels, comment out TH |
| Create | `tests/test_condensibles.py` | Unit tests for `condensibles.py` |
| Modify | `tests/test_atmospheric_properties.py` | Tests for `get_active_condensible()` and condensible constants from TOML |

---

## Task 1: Create `condensibles.py` — CondensibleParams, defaults, factory

**Files:**
- Create: `climt/_core/condensibles.py`
- Create: `tests/test_condensibles.py`

- [x] **Step 1: Write the failing test**

```python
# tests/test_condensibles.py
import pytest


def test_condensible_params_h2o_defaults():
    """get_condensible_params() returns H2O values without loading any profile."""
    from climt._core.condensibles import get_condensible_params, SPECIES_ID

    cond = get_condensible_params()
    assert cond.species_id == SPECIES_ID["h2o"]  # 0
    assert abs(cond.RV - 461.5) < 0.1
    assert abs(cond.CPV - 1870.0) < 1.0
    assert abs(cond.CL - 2500.0) < 1.0
    assert abs(cond.LV0 - 2.501e6) < 1e3
    assert abs(cond.ROWL - 1000.0) < 0.1
    assert abs(cond.T_freeze - 273.15) < 0.01


def test_species_id_mapping():
    from climt._core.condensibles import SPECIES_ID

    assert SPECIES_ID["h2o"] == 0
    assert SPECIES_ID["ch4"] == 1
    assert SPECIES_ID["co2"] == 2
```

- [x] **Step 2: Run test to confirm it fails**

```bash
conda run -n climt pytest tests/test_condensibles.py -v
```
Expected: `ModuleNotFoundError` or `ImportError` — `condensibles` does not exist yet.

- [x] **Step 3: Create `climt/_core/condensibles.py` with CondensibleParams, defaults, factory**

```python
# climt/_core/condensibles.py
from typing import NamedTuple

import numpy as np
from sympl import get_constant, set_constant

from .atmospheric_properties import get_active_condensible

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    def njit(x, **kwargs):
        return x

SPECIES_ID = {"h2o": 0, "ch4": 1, "co2": 2}

# Register H2O condensible constants as sympl defaults.
# load_atmospheric_properties() overrides these when a profile is loaded.
_H2O_DEFAULTS = {
    "gas_constant_of_condensible_vapor":       (461.5,   "J/kg/K"),
    "heat_capacity_of_condensible_vapor":      (1870.0,  "J/kg/K"),
    "heat_capacity_of_condensible_liquid":     (2500.0,  "J/kg/K"),
    "latent_heat_of_condensible_vaporization": (2.501e6, "J/kg"),
    "density_of_condensible_liquid":           (1000.0,  "kg/m^3"),
    "freezing_point_of_condensible":           (273.15,  "degK"),
}
for _name, (_val, _units) in _H2O_DEFAULTS.items():
    set_constant(_name, _val, _units)


class CondensibleParams(NamedTuple):
    species_id: int
    RV: float        # gas constant of condensible vapor (J/kg/K)
    CPV: float       # heat capacity of condensible vapor (J/kg/K)
    CL: float        # heat capacity of condensible liquid/solid (J/kg/K)
    LV0: float       # latent heat of vaporization at T_freeze (J/kg)
    ROWL: float      # density of condensible liquid/solid (kg/m^3)
    T_freeze: float  # phase transition temperature (K)


def get_condensible_params() -> CondensibleParams:
    """Read active condensible constants from sympl and return a CondensibleParams.

    Call once per component __init__. Requires either H2O defaults (registered
    at module import) or a profile loaded via load_atmospheric_properties().
    """
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

- [x] **Step 4: Run tests to confirm they pass**

```bash
conda run -n climt pytest tests/test_condensibles.py::test_condensible_params_h2o_defaults tests/test_condensibles.py::test_species_id_mapping -v
```
Expected: both PASS.

- [x] **Step 5: Commit**

```bash
git add climt/_core/condensibles.py tests/test_condensibles.py
git commit -m "feat(condensibles): CondensibleParams, H2O defaults, get_condensible_params"
```

---

## Task 2: Add saturation kernel functions to `condensibles.py`

**Files:**
- Modify: `climt/_core/condensibles.py`
- Modify: `tests/test_condensibles.py`

- [x] **Step 1: Write the failing tests**

Add to `tests/test_condensibles.py`:

```python
def test_sat_vap_pressure_h2o_above_freezing():
    """H2O saturation pressure at 273.15 K (0°C) is ~6.112 hPa."""
    from climt._core.condensibles import _sat_vap_pressure

    es = _sat_vap_pressure(273.15, 0)
    assert abs(es - 6.112) < 0.01


def test_sat_vap_pressure_h2o_below_freezing():
    """H2O ice saturation formula is active below 273.15 K."""
    from climt._core.condensibles import _sat_vap_pressure

    # At -20°C (253.15 K), ice saturation should be ~1.03 hPa
    es = _sat_vap_pressure(253.15, 0)
    assert 0.9 < es < 1.2


def test_sat_vap_pressure_ch4_at_boiling():
    """CH4 saturation pressure at 111.65 K equals 1 atm (1013.25 hPa) by construction."""
    from climt._core.condensibles import _sat_vap_pressure

    es = _sat_vap_pressure(111.65, 1)
    assert abs(es - 1013.25) < 1.0


def test_sat_vap_pressure_co2_at_sublimation():
    """CO2 saturation pressure at 194.7 K equals 1 atm by construction."""
    from climt._core.condensibles import _sat_vap_pressure

    es = _sat_vap_pressure(194.7, 2)
    assert abs(es - 1013.25) < 1.0


def test_lv_at_freeze():
    """_lv at T_freeze returns LV0."""
    from climt._core.condensibles import _lv, get_condensible_params

    cond = get_condensible_params()
    lv = _lv(cond.T_freeze, cond)
    assert abs(lv - cond.LV0) < 1.0


def test_lcl_pressure_saturated():
    """_lcl_pressure with RH=1.0 returns P unchanged (both H2O and CC formula)."""
    from climt._core.condensibles import _lcl_pressure, get_condensible_params

    cond = get_condensible_params()
    P, T = 1000.0, 300.0
    plcl = _lcl_pressure(P, 1.0, T, cond)
    assert abs(plcl - P) < 0.1


def test_lcl_pressure_subsaturated():
    """_lcl_pressure with RH < 1.0 returns PLCL < P."""
    from climt._core.condensibles import _lcl_pressure, get_condensible_params

    cond = get_condensible_params()
    plcl = _lcl_pressure(1000.0, 0.7, 300.0, cond)
    assert plcl < 1000.0
    assert plcl > 500.0  # sanity bound
```

- [x] **Step 2: Run tests to confirm they fail**

```bash
conda run -n climt pytest tests/test_condensibles.py -v
```
Expected: all new tests FAIL with `ImportError` — functions not defined yet.

- [x] **Step 3: Add `_sat_vap_pressure`, `_lv`, `_lcl_pressure` to `condensibles.py`**

Append after `get_condensible_params()`:

```python
@njit
def _sat_vap_pressure(T, species_id):
    """Saturation vapor pressure in hPa for the given species.

    H2O: Magnus formula (above/below freezing).
    CH4, CO2: Clausius-Clapeyron anchored at 1 atm reference point.
    """
    if species_id == 0:  # H2O
        TC = T - 273.15
        if TC >= 0.0:
            return 6.112 * np.exp(17.67 * TC / (243.5 + TC))
        else:
            return 0.01 * np.exp(23.33086 - 6111.72784 / T + 0.15215 * np.log(T))
    elif species_id == 1:  # CH4 — anchored at normal boiling point 111.65 K, 1 atm
        return 1013.25 * np.exp((5.1e5 / 518.3) * (1.0 / 111.65 - 1.0 / T))
    else:  # CO2 — anchored at sublimation point 194.7 K, 1 atm
        return 1013.25 * np.exp((5.71e5 / 188.9) * (1.0 / 194.7 - 1.0 / T))


@njit
def _lv(T, cond):
    """Latent heat of vaporization/sublimation (J/kg) as a function of temperature."""
    return cond.LV0 - (cond.CL - cond.CPV) * (T - cond.T_freeze)


@njit
def _lcl_pressure(P_nk, RH, T_nk, cond):
    """Lifting Condensation Level pressure (hPa).

    H2O: Bolton/Bohren empirical formula.
    CH4, CO2: Clausius-Clapeyron approximation PLCL = P * RH^(Rv*T/Lv).
    """
    if cond.species_id == 0:  # H2O empirical
        CHI = T_nk / (1669.0 - 122.0 * RH - T_nk)
        return P_nk * (RH ** CHI)
    else:  # Clausius-Clapeyron
        chi = cond.RV * T_nk / cond.LV0
        return P_nk * (RH ** chi)
```

- [x] **Step 4: Run tests**

```bash
conda run -n climt pytest tests/test_condensibles.py -v
```
Expected: all tests PASS.

- [x] **Step 5: Commit**

```bash
git add climt/_core/condensibles.py tests/test_condensibles.py
git commit -m "feat(condensibles): _sat_vap_pressure, _lv, _lcl_pressure kernel helpers"
```

---

## Task 3: Add `compute_qs` to `condensibles.py`

**Files:**
- Modify: `climt/_core/condensibles.py`
- Modify: `tests/test_condensibles.py`

- [x] **Step 1: Write the failing test**

Add to `tests/test_condensibles.py`:

```python
def test_compute_qs_matches_bolton_h2o():
    """compute_qs must match bolton_q_sat for H2O above freezing to within rtol=1e-5."""
    import numpy as np
    from climt._core.condensibles import compute_qs, get_condensible_params
    from climt._core.util import bolton_q_sat

    cond = get_condensible_params()
    # All temperatures above 273.15 K so both formulas use the same (liquid) branch
    T = np.array([[280.0], [290.0], [300.0], [310.0]])   # (4, 1)
    P_hpa = np.array([[950.0], [900.0], [850.0], [800.0]])

    qs_new = compute_qs(T, P_hpa, cond, 287.04)
    qs_ref = bolton_q_sat(T, P_hpa * 100.0, 287.04, 461.5)
    np.testing.assert_allclose(qs_new, qs_ref, rtol=1e-5)
```

- [x] **Step 2: Run test to confirm it fails**

```bash
conda run -n climt pytest tests/test_condensibles.py::test_compute_qs_matches_bolton_h2o -v
```
Expected: FAIL — `compute_qs` not defined.

- [x] **Step 3: Add `compute_qs` to `condensibles.py`**

Add the following imports at the top of `condensibles.py` (after the existing imports):

```python
from ..._core.backend import jit_compile, prange  # noqa: F401  (prange used inside jit_compile)
```

Wait — `condensibles.py` is *in* `climt/_core/`, so the import is:

```python
from .backend import jit_compile, prange
```

Append `compute_qs` after `_lcl_pressure`:

```python
@jit_compile(backend=np, parallel=True)
def compute_qs(T, P, cond, RD):
    """Saturation specific humidity (kg/kg) for a 2-D (nlev, ncol) array.

    Args:
        T:   air temperature (K), shape (nlev, ncol)
        P:   air pressure (hPa), shape (nlev, ncol)
        cond: CondensibleParams for the active condensible
        RD:  specific gas constant of dry air (J/kg/K)

    Returns:
        qs array with same shape as T.
    """
    nlev, ncol = T.shape
    qs = np.zeros_like(T)
    EPS = RD / cond.RV
    for i in prange(ncol):
        for k in range(nlev):
            es = _sat_vap_pressure(T[k, i], cond.species_id)
            qs[k, i] = EPS * es / (P[k, i] - (1.0 - EPS) * es)
    return qs
```

- [x] **Step 4: Run tests**

```bash
conda run -n climt pytest tests/test_condensibles.py -v
```
Expected: all PASS.

- [x] **Step 5: Commit**

```bash
git add climt/_core/condensibles.py tests/test_condensibles.py
git commit -m "feat(condensibles): compute_qs — parallelised saturation specific humidity"
```

---

## Task 4: Extend `atmospheric_properties.py` with condensible tracking

**Files:**
- Modify: `climt/_core/atmospheric_properties.py`
- Modify: `tests/test_atmospheric_properties.py`

- [x] **Step 1: Write the failing tests**

Add to `tests/test_atmospheric_properties.py`:

```python
def test_get_active_condensible_default():
    """get_active_condensible() returns 'h2o' before any profile is loaded."""
    from climt._core.atmospheric_properties import get_active_condensible

    assert get_active_condensible() == "h2o"


def test_reset_restores_condensible():
    """reset_atmospheric_properties() restores condensible to 'h2o'."""
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from climt._core.atmospheric_properties import get_active_condensible

    load_atmospheric_properties("earth")
    try:
        load_atmospheric_properties("titan")
        assert get_active_condensible() == "ch4"
        reset_atmospheric_properties()
        assert get_active_condensible() == "h2o"
    finally:
        reset_atmospheric_properties()
```

Note: the `titan` test will only pass after Task 5 adds `[condensible]` to titan.toml. For now, run only the default test.

- [x] **Step 2: Run the default test to confirm it fails**

```bash
conda run -n climt pytest tests/test_atmospheric_properties.py::test_get_active_condensible_default -v
```
Expected: FAIL — `get_active_condensible` not defined.

- [x] **Step 3: Add `_active_condensible` state and `get_active_condensible()` to `atmospheric_properties.py`**

After the `_active_profile` dict (around line 31), add:

```python
_active_condensible = {"species": "h2o"}
```

After the `_active_profile` dict definition, add the public function (can be placed near the bottom with the other public functions):

```python
def get_active_condensible() -> str:
    """Return the species string of the currently loaded condensible (e.g. 'ch4')."""
    return _active_condensible["species"]
```

In `load_atmospheric_properties()`, after `_active_profile["path"] = path`, add:

```python
    _active_condensible["species"] = data.get("condensible", {}).get("species", "h2o")
```

Note: `data` is the raw TOML dict from `tomllib.load()`. Currently `_parse_toml` only returns the constants dict — you need access to the raw TOML dict in `load_atmospheric_properties`. Modify `load_atmospheric_properties` to call `tomllib.load()` directly (before `_parse_toml`) or pass the raw data through. The cleanest approach: open the file once and pass raw data:

```python
def load_atmospheric_properties(name_or_path):
    path = _resolve_profile_path(name_or_path)
    with open(path, "rb") as f:
        raw = tomllib.load(f)
    constants = _parse_toml_from_dict(raw)   # rename _parse_toml to accept dict

    _snapshot_stack.append(_snapshot_constants(constants))
    _active_profile["name"] = os.path.splitext(os.path.basename(path))[0]
    _active_profile["path"] = path
    _active_condensible["species"] = raw.get("condensible", {}).get("species", "h2o")

    for key, (value, units) in constants.items():
        set_constant(key, value, units)
```

Rename the existing `_parse_toml` to `_parse_toml_from_dict(data)` that accepts a dict instead of a path, and update the existing `_parse_toml(path)` to call it:

```python
def _parse_toml_from_dict(data):
    """Extract sympl constants from a parsed TOML dict."""
    constants = {}
    for section in data.values():
        if isinstance(section, dict):
            for key, val in section.items():
                if isinstance(val, dict) and "value" in val and "units" in val:
                    constants[key] = (val["value"], val["units"])
    return constants


def _parse_toml(path):
    """Parse a TOML profile and return a flat dict of {name: (value, units)}."""
    with open(path, "rb") as f:
        data = tomllib.load(f)
    return _parse_toml_from_dict(data)
```

In `reset_atmospheric_properties()`, after restoring constants, add:

```python
    _active_condensible["species"] = "h2o"
    _active_profile["name"] = None
    _active_profile["path"] = None
```

- [x] **Step 4: Run the default test**

```bash
conda run -n climt pytest tests/test_atmospheric_properties.py::test_get_active_condensible_default -v
```
Expected: PASS.

- [x] **Step 5: Run full existing atmospheric properties test suite to check for regressions**

```bash
conda run -n climt pytest tests/test_atmospheric_properties.py -v
```
Expected: all PASS.

- [x] **Step 6: Commit**

```bash
git add climt/_core/atmospheric_properties.py tests/test_atmospheric_properties.py
git commit -m "feat(atmospheric_properties): condensible species tracking via get_active_condensible"
```

---

## Task 5: Add condensible sections to TOML files

**Files:**
- Modify: `climt/_data/atmospheric_properties/earth.toml`
- Modify: `climt/_data/atmospheric_properties/titan.toml`
- Modify: `climt/_data/atmospheric_properties/mars.toml`

- [x] **Step 1: Add `[condensible]` and `[condensible_properties]` to `earth.toml`**

Append to `climt/_data/atmospheric_properties/earth.toml`:

```toml
[condensible]
species = "h2o"

[condensible_properties]
gas_constant_of_condensible_vapor         = { value = 461.5,   units = "J/kg/K" }
heat_capacity_of_condensible_vapor        = { value = 1870.0,  units = "J/kg/K" }
heat_capacity_of_condensible_liquid       = { value = 2500.0,  units = "J/kg/K" }
latent_heat_of_condensible_vaporization   = { value = 2.501e6, units = "J/kg" }
density_of_condensible_liquid             = { value = 1000.0,  units = "kg/m^3" }
freezing_point_of_condensible            = { value = 273.15,  units = "degK" }
```

Note: `CL = 2500.0` follows the original Emanuel scheme convention (not the true liquid water value of ~4186 J/kg/K). Do not change this — it preserves numerical parity with existing tests.

- [x] **Step 2: Add to `titan.toml`**

Append to `climt/_data/atmospheric_properties/titan.toml`:

```toml
[condensible]
species = "ch4"

[condensible_properties]
# Methane condensible properties
# Latent heat anchored at normal boiling point (111.65 K)
gas_constant_of_condensible_vapor         = { value = 518.3,  units = "J/kg/K" }
heat_capacity_of_condensible_vapor        = { value = 2232.0, units = "J/kg/K" }
heat_capacity_of_condensible_liquid       = { value = 3469.0, units = "J/kg/K" }
latent_heat_of_condensible_vaporization   = { value = 5.1e5,  units = "J/kg" }
density_of_condensible_liquid             = { value = 422.6,  units = "kg/m^3" }
freezing_point_of_condensible            = { value = 90.69,  units = "degK" }
```

- [x] **Step 3: Add to `mars.toml`**

Append to `climt/_data/atmospheric_properties/mars.toml`:

```toml
[condensible]
species = "co2"

[condensible_properties]
# CO2 condensible properties — condensate is solid (dry ice) at Mars surface conditions
# Latent heat is sublimation heat anchored at 194.7 K (1 atm sublimation point)
# T_freeze set to CO2 triple point; at Mars pressures CO2 always precipitates as dry ice
gas_constant_of_condensible_vapor         = { value = 188.9,  units = "J/kg/K" }
heat_capacity_of_condensible_vapor        = { value = 846.0,  units = "J/kg/K" }
heat_capacity_of_condensible_liquid       = { value = 709.0,  units = "J/kg/K" }
latent_heat_of_condensible_vaporization   = { value = 5.71e5, units = "J/kg" }
density_of_condensible_liquid             = { value = 1562.0, units = "kg/m^3" }
freezing_point_of_condensible            = { value = 216.6,  units = "degK" }
```

- [x] **Step 4: Run TOML-loading tests**

```bash
conda run -n climt pytest tests/test_atmospheric_properties.py -v
```
Expected: all PASS, including the new `test_reset_restores_condensible` (titan now has `[condensible]`).

- [x] **Step 5: Run a quick smoke test that constants load correctly after loading titan**

```bash
conda run -n climt python -c "
from climt import load_atmospheric_properties, reset_atmospheric_properties
from sympl import get_constant
from climt._core.atmospheric_properties import get_active_condensible
load_atmospheric_properties('titan')
print('species:', get_active_condensible())
print('RV:', get_constant('gas_constant_of_condensible_vapor', 'J/kg/K'))
reset_atmospheric_properties()
print('after reset:', get_active_condensible())
"
```
Expected output:
```
species: ch4
RV: 518.3
after reset: h2o
```

- [x] **Step 6: Commit**

```bash
git add climt/_data/atmospheric_properties/earth.toml \
        climt/_data/atmospheric_properties/titan.toml \
        climt/_data/atmospheric_properties/mars.toml
git commit -m "feat(toml): add condensible species + properties to earth/titan/mars profiles"
```

---

## Task 6: Export new public API from `climt/_core/__init__.py`

**Files:**
- Modify: `climt/_core/__init__.py`

- [x] **Step 1: Add imports and exports**

In `climt/_core/__init__.py`, add the following import after the `atmospheric_properties` imports:

```python
from .condensibles import (
    CondensibleParams,
    get_condensible_params,
)
from .atmospheric_properties import get_active_condensible
```

Add to the `__all__` tuple:

```python
    CondensibleParams,
    get_condensible_params,
    get_active_condensible,
```

- [x] **Step 2: Confirm the imports work**

```bash
conda run -n climt python -c "
from climt._core import CondensibleParams, get_condensible_params, get_active_condensible
print('imports ok')
cond = get_condensible_params()
print('species_id:', cond.species_id, 'RV:', cond.RV)
"
```
Expected:
```
imports ok
species_id: 0 RV: 461.5
```

- [x] **Step 3: Commit**

```bash
git add climt/_core/__init__.py
git commit -m "feat(core): export CondensibleParams, get_condensible_params, get_active_condensible"
```

---

## Task 7: Update `EmanuelParams` and `__init__`

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

- [x] **Step 1: Run baseline parity test to establish a reference**

```bash
conda run -n climt pytest tests/test_emanuel_v3_parity.py -v
```
Expected: all PASS. Record this as the baseline.

- [x] **Step 2: Remove condensible fields from `EmanuelParams`**

In `pure_python_v3.py`, remove the following five fields from `EmanuelParams`:

```python
# Remove these lines from class EmanuelParams(NamedTuple):
    CPV: float
    CL: float
    RV: float
    LV0: float
    ROWL: float
```

`EmanuelParams` should now end at `DELT0: float`.

- [x] **Step 3: Add condensibles import and update `__init__`**

At the top of `pure_python_v3.py`, add:

```python
from ..._core.condensibles import CondensibleParams, get_condensible_params
```

In `__init__`, remove the five hardcoded lines:
```python
# Remove:
self.CPV = 1870.0
self.CL = 2500.0
self.RV = 461.5
self.LV0 = 2.501e6
self.ROWL = 1000.0
```

And remove their corresponding entries in `EmanuelParams(...)`:
```python
# Remove from self._params = EmanuelParams(...):
            CPV=float(self.CPV),
            CL=float(self.CL),
            RV=float(self.RV),
            LV0=float(self.LV0),
            ROWL=float(self.ROWL),
```

After the `self._params = EmanuelParams(...)` assignment, add:

```python
        self._cond = get_condensible_params()
```

- [x] **Step 4: Run parity test — expect failure** (kernels still reference `params.RV` etc.)

```bash
conda run -n climt pytest tests/test_emanuel_v3_parity.py -v 2>&1 | head -40
```
Expected: errors about `params` missing `CPV` / `RV` etc. This is correct — kernel updates follow in Tasks 8–10.

- [x] **Step 5: Commit the params change as a work-in-progress**

```bash
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "refactor(emanuel): remove condensible fields from EmanuelParams, add self._cond"
```

---

## Task 8: Update `_tlift_functional_np`

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

- [x] **Step 1: Add `cond` parameter and update all condensible references**

Replace the function signature:

```python
# Old:
def _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params):
# New:
def _tlift_functional_np(P, T, Q, QS, GZ, ICB, NK, ND, NL, KK, TVP, TPK, CLW, params, cond):
```

Make the following replacements inside the function body (use line numbers as guidance; locate each by the expression shown):

| Find | Replace |
|---|---|
| `CPVMCL = params.CL - params.CPV` | `CPVMCL = cond.CL - cond.CPV` |
| `EPS = params.RD / params.RV` | `EPS = params.RD / cond.RV` |
| `params.CL * Q[NK]` (in `AH0`) | `cond.CL * Q[NK]` |
| `params.LV0 - CPVMCL * (T[NK] - 273.15)` (in `AH0`) | `cond.LV0 - CPVMCL * (T[NK] - cond.T_freeze)` |
| `params.CPD * (1.0 - Q[NK]) + Q[NK] * params.CPV` (CPP) | `params.CPD * (1.0 - Q[NK]) + Q[NK] * cond.CPV` |
| `ALV = params.LV0 - CPVMCL * (T[i] - 273.15)` | `ALV = cond.LV0 - CPVMCL * (T[i] - cond.T_freeze)` |
| `params.RV * T[i] * T[i]` (in `S`) | `cond.RV * T[i] * T[i]` |
| `(params.CL - params.CPD) * Q[NK] * T[i]` (in `AHG`) | `(cond.CL - params.CPD) * Q[NK] * T[i]` |

Replace the Magnus formula block (find by `TC = TG - 273.15`):

```python
# Remove these lines:
        TC = TG - 273.15
        DENOM = 243.5 + TC
        ES = (
            6.112 * np.exp(17.67 * TC / DENOM)
            if TC >= 0.0
            else np.exp(23.33086 - 6111.72784 / TG + 0.15215 * np.log(TG))
        )
        QG = EPS * ES / (P[i] - ES * (1.0 - EPS))

# Replace with:
        ES = _sat_vap_pressure(TG, cond.species_id)
        QG = EPS * ES / (P[i] - ES * (1.0 - EPS))
```

Add import of `_sat_vap_pressure` at top of the `@njit` function (it's already in the same file, so no import needed — just confirm `_sat_vap_pressure` is defined before `_tlift_functional_np`).

Replace the `TPK` line that references `params.CL`:
```python
# Old:
    TPK[i] = (
        AH0 - (params.CL - params.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG
    ) / params.CPD
# New:
    TPK[i] = (
        AH0 - (cond.CL - params.CPD) * Q[NK] * T[i] - GZ[i] - ALV * QG
    ) / params.CPD
```

- [x] **Step 2: Update the two call sites of `_tlift_functional_np` in `_convect_functional_np`**

There are two calls (search for `_tlift_functional_np(`). Both need `cond` appended:

```python
# Old:
    TVP, TP, CLW = _tlift_functional_np(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params
    )
# New:
    TVP, TP, CLW = _tlift_functional_np(
        P, T, Q, QS, GZ, ICB, NK, ND, NL, 1, TVP, TP, CLW, params, cond
    )
```

Do the same for the second call (KK=2). Note: `_convect_functional_np` doesn't have `cond` yet — that comes in Task 9. For now add `cond` as an argument to `_convect_functional_np` signature only (to unblock compilation), and pass `None` temporarily in `_numpy_vectorized_convect`. Actually, it's cleaner to update `_convect_functional_np` signature now even if its body changes come in Task 9. Add `cond` to its signature and pass it through immediately.

- [x] **Step 3: Run parity test**

```bash
conda run -n climt pytest tests/test_emanuel_v3_parity.py -v
```
Expected: still failing (more references to `params.RV` etc. remain in `_convect_functional_np`). That's fine — continue.

- [x] **Step 4: Commit**

```bash
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "refactor(emanuel): _tlift_functional_np — use cond for condensible properties"
```

---

## Task 9: Update `_convect_functional_np` — body changes + TH comment-out

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

This task replaces all remaining `params.CPV`, `params.CL`, `params.RV`, `params.LV0`, `params.ROWL` references in `_convect_functional_np` with `cond.*` equivalents, and comments out the dead `TH` block.

- [x] **Step 1: Comment out the `TH` block**

Locate the block starting with `TH = np.zeros(NL + 1)` (around line 314). Comment it out:

```python
    # TH (potential temperature) is computed but never used in this function.
    # It is a carry-over from the original Fortran code. Commented out rather
    # than deleted to preserve the original structure. The 1000 hPa reference
    # pressure is Earth-specific; generalisation deferred to a future pass.
    # TH = np.zeros(NL + 1)
    # for i in range(NL + 1):
    #     RDCP = (params.RD * (1.0 - Q[i]) + Q[i] * params.RV) / (
    #         params.CPD * (1.0 - Q[i]) + Q[i] * params.CPV
    #     )
    #     TH[i] = T[i] * (1000.0 / P[i]) ** RDCP
```

- [x] **Step 2: Replace condensible references in `_convect_functional_np` body**

Make these replacements (each `find` is unique in the function):

| Find | Replace |
|---|---|
| `CPVMCL = params.CL - params.CPV` | `CPVMCL = cond.CL - cond.CPV` |
| `EPS = params.RD / params.RV` | `EPS = params.RD / cond.RV` |
| `CPN[0] = params.CPD * (1.0 - Q[0]) + Q[0] * params.CPV` | `CPN[0] = params.CPD * (1.0 - Q[0]) + Q[0] * cond.CPV` |
| `LV[0] = params.LV0 - CPVMCL * (T[0] - 273.15)` | `LV[0] = cond.LV0 - CPVMCL * (T[0] - cond.T_freeze)` |
| `CPN[i] = params.CPD * (1.0 - Q[i]) + params.CPV * Q[i]` | `CPN[i] = params.CPD * (1.0 - Q[i]) + cond.CPV * Q[i]` |
| `LV[i] = params.LV0 - CPVMCL * (T[i] - 273.15)` | `LV[i] = cond.LV0 - CPVMCL * (T[i] - cond.T_freeze)` |
| `(params.CPD * (1.0 - Q[i]) + params.CL * Q[i]) * (T[i] - T[0])` (in HM) | `(params.CPD * (1.0 - Q[i]) + cond.CL * Q[i]) * (T[i] - T[0])` |
| `CHI = T[NK] / (1669.0 - 122.0 * RH - T[NK])` + `PLCL = P[NK] * (RH**CHI)` | `PLCL = _lcl_pressure(P[NK], RH, T[NK], cond)` |
| `T[i] > 273.0` (rain/snow check, line ~602) | `T[i] > cond.T_freeze` |
| `(LV[i] + (params.CPD - params.CPV) * T[i]) * EP[i]` | `(LV[i] + (params.CPD - cond.CPV) * T[i]) * EP[i]` |
| `LV[j] * LV[j] * QS[j] / (params.RV * T[j]` | `LV[j] * LV[j] * QS[j] / (cond.RV * T[j]` |
| `(params.CPV - params.CPD) * T[j] * (QTI - Q[j])` (in ANUM) | `(cond.CPV - params.CPD) * T[j] * (QTI - Q[j])` |
| `(params.CPD - params.CPV) * (Q[i] - QTI) * T[j]` (in DENOM) | `(params.CPD - cond.CPV) * (Q[i] - QTI) * T[j]` |
| `QP[i + 1] * (LV[i + 1] + T[i + 1] * (params.CL - params.CPD))` | `QP[i + 1] * (LV[i + 1] + T[i + 1] * (cond.CL - params.CPD))` |
| `/ (LV[i] + T[i] * (params.CL - params.CPD))` | `/ (LV[i] + T[i] * (cond.CL - params.CPD))` |
| `/ (params.ROWL * params.G)` (PRECIP line) | `/ (cond.ROWL * params.G)` |
| `TPRIME = params.LV0 * QPRIME / params.CPD` | `TPRIME = cond.LV0 * QPRIME / params.CPD` |
| `T[i] * (params.CPV - params.CPD) * (Q[i] - QENT[i, i])` (in FT[i]) | `T[i] * (cond.CPV - params.CPD) * (Q[i] - QENT[i, i])` |
| `(params.CL - params.CPD) * WATER[1]` (in FT[0]) | `(cond.CL - params.CPD) * WATER[1]` |
| `(params.CL - params.CPD) * WATER[i + 1]` (in FT[i] loop) | `(cond.CL - params.CPD) * WATER[i + 1]` |

Add `_lcl_pressure` to the imports from `condensibles` at the top of the file:

```python
from ..._core.condensibles import CondensibleParams, get_condensible_params, _lcl_pressure, _sat_vap_pressure
```

- [x] **Step 3: Run parity test**

```bash
conda run -n climt pytest tests/test_emanuel_v3_parity.py -v
```
Expected: still failing — `_numpy_vectorized_convect` and `array_call` haven't been updated yet.

- [x] **Step 4: Commit**

```bash
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "refactor(emanuel): _convect_functional_np — cond param, comment TH, condensible props"
```

---

## Task 10: Update `_numpy_vectorized_convect` and `array_call`

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

- [x] **Step 1: Thread `cond` through `_numpy_vectorized_convect`**

Update the function signature:

```python
# Old:
def _numpy_vectorized_convect(
    t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, tra, params
):
# New:
def _numpy_vectorized_convect(
    t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, tra, params, cond
):
```

Inside the function, add `cond` to the `_convect_functional_np` call:

```python
        res = _convect_functional_np(
            t[:, i],
            q[:, i],
            qs[:, i],
            u[:, i],
            v[:, i],
            p[:, i],
            ph[:, i],
            ND,
            NL,
            NTRA,
            DELT,
            cbmf[i],
            tra[:, :, i],
            params,
            cond,       # ← add this
        )
```

- [x] **Step 2: Update `array_call`**

Add `compute_qs` to the condensibles import at the top of `pure_python_v3.py`:

```python
from ..._core.condensibles import (
    CondensibleParams,
    get_condensible_params,
    _lcl_pressure,
    _sat_vap_pressure,
    compute_qs,
)
```

In `array_call`, replace:

```python
        from climt._core import bolton_q_sat
        qs = bolton_q_sat(t, p * 100, self.RD, self.RV)
```

with:

```python
        qs = compute_qs(t, p, self._cond, self._params.RD)
```

Note: `compute_qs` receives pressure in hPa (same units as `p` in the Emanuel kernels) — do **not** multiply by 100.

Update the call to `_numpy_vectorized_convect` to include `cond`:

```python
        results = _numpy_vectorized_convect(
            t, q, qs, u, v, p, ph, nlev, nlev - 3, ntra, delt, cbmf, tra_vector,
            self._params,
            self._cond,   # ← add this
        )
```

- [x] **Step 3: Run parity test**

```bash
conda run -n climt pytest tests/test_emanuel_v3_parity.py -v
```
Expected: all PASS. If any fail, check for remaining `params.RV` / `params.CL` / `params.LV0` references with:

```bash
grep -n "params\.RV\|params\.CPV\|params\.CL\b\|params\.LV0\|params\.ROWL" \
    climt/_components/emanuel/pure_python_v3.py
```

Expected: no output.

- [x] **Step 4: Run full existing test suite for Emanuel**

```bash
conda run -n climt pytest tests/test_emanuel_v3_parity.py tests/test_emanuel_optimization.py tests/test_components.py -v
```
Expected: all PASS.

- [x] **Step 5: Commit**

```bash
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "feat(emanuel): wire cond through vectorized kernel and array_call; use compute_qs"
```

---

## Task 11: Titan (CH₄) smoke test

**Files:**
- Modify: `tests/test_condensibles.py`

- [x] **Step 1: Write Titan integration smoke test**

Add to `tests/test_condensibles.py`:

```python
def test_emanuel_titan_ch4_smoke():
    """EmanuelConvectionPython runs without error on a Titan-like CH4 atmosphere."""
    import numpy as np
    from datetime import timedelta
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPython

    load_atmospheric_properties("titan")
    try:
        nlev, ncol = 30, 1
        # Titan surface ~94 K, tropopause ~70 K, lapse rate ~1 K/km
        T = np.linspace(94.0, 70.0, nlev).reshape(nlev, ncol)
        # CH4 specific humidity: ~5% near surface, decreasing with altitude
        q = np.linspace(0.05, 0.001, nlev).reshape(nlev, ncol)
        u = np.zeros((nlev, ncol))
        v = np.zeros((nlev, ncol))
        # Titan surface pressure ~1467 hPa
        p = np.linspace(1467.0, 10.0, nlev).reshape(nlev, ncol)
        ph = np.linspace(1500.0, 5.0, nlev + 1).reshape(nlev + 1, ncol)
        cbmf = np.zeros(ncol)

        conv = EmanuelConvectionPython()
        state = {
            "air_temperature": T,
            "specific_humidity": q,
            "eastward_wind": u,
            "northward_wind": v,
            "air_pressure": p,
            "air_pressure_on_interface_levels": ph,
            "cloud_base_mass_flux": cbmf,
        }
        # Should complete without exception; output values not checked here
        tendencies, diagnostics = conv.array_call(state, timedelta(minutes=15))
        assert "air_temperature" in tendencies
        assert tendencies["air_temperature"].shape == (nlev, ncol)
    finally:
        reset_atmospheric_properties()
```

- [x] **Step 2: Run the smoke test**

```bash
conda run -n climt pytest tests/test_condensibles.py::test_emanuel_titan_ch4_smoke -v
```
Expected: PASS (no exception; tendencies have correct shape).

- [x] **Step 3: Run the full condensibles test suite**

```bash
conda run -n climt pytest tests/test_condensibles.py -v
```
Expected: all PASS.

- [x] **Step 4: Run the complete test suite to check for regressions**

```bash
conda run -n climt pytest tests/ -v --ignore=tests/test_rrtmg_comparison.py -x
```
Expected: all PASS. (`test_rrtmg_comparison.py` is excluded as it requires RRTMG data files that may not be present.)

- [x] **Step 5: Commit**

```bash
git add tests/test_condensibles.py
git commit -m "test(condensibles): Titan CH4 smoke test for EmanuelConvectionPython"
```

---

## Self-Review Checklist (completed inline)

**Spec coverage:**
- ✅ TOML `[condensible]` + `[condensible_properties]` — Tasks 5
- ✅ `atmospheric_properties.py` `_active_condensible` + `get_active_condensible()` — Task 4
- ✅ `condensibles.py` — Tasks 1–3
- ✅ `CondensibleParams` NamedTuple — Task 1
- ✅ `get_condensible_params()` factory — Task 1
- ✅ `_sat_vap_pressure` with H2O/CH4/CO2 dispatch — Task 2
- ✅ `_lv` — Task 2
- ✅ `_lcl_pressure` receives `cond` not just `species_id` — Task 2
- ✅ `compute_qs` parallelised with `prange` — Task 3
- ✅ `EmanuelParams` loses condensible fields — Task 7
- ✅ `__init__` adds `self._cond = get_condensible_params()` — Task 7
- ✅ `array_call` uses `compute_qs(t, p, ...)` not `bolton_q_sat` — Task 10
- ✅ `_tlift_functional_np` gains `cond`, uses `_sat_vap_pressure`, `cond.T_freeze` — Task 8
- ✅ `_convect_functional_np` gains `cond`, TH commented out, LCL uses `_lcl_pressure(cond)`, rain/snow uses `cond.T_freeze` — Task 9
- ✅ `_numpy_vectorized_convect` threads `cond` — Task 10
- ✅ H2O defaults registered at import — Task 1 (backward compatible with all existing tests)
- ✅ Titan smoke test — Task 11
- ✅ Exports from `_core/__init__.py` — Task 6

**Placeholder scan:** None found.

**Type consistency:** `_lcl_pressure(P_nk, RH, T_nk, cond)` defined in Task 2, called with same signature in Task 9. `_sat_vap_pressure(T, species_id)` defined in Task 2, called in Task 8. `compute_qs(T, P, cond, RD)` defined in Task 3, called in Task 10. `CondensibleParams` defined in Task 1, used throughout. ✓
