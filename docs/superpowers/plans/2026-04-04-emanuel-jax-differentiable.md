# Emanuel Convection JAX Differentiable Port — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring `EmanuelConvectionPythonV3` to full parity with the develop branch, then make the JAX path fully differentiable so `jax.grad` works through `array_call`.

**Architecture:** Three phases. Phase 1 copies the complete numpy/numba kernel from develop (mechanical). Phase 2 writes a full JAX version using Python loops — correct but not JIT-able, serving as a reference. Phase 3 converts each algorithmic block to use `jax.lax.scan`, `jnp.where`, and masked array ops so the whole thing is `jax.jit`-able and differentiable. Each phase has parity tests against the previous.

**Tech Stack:** Python, JAX, NumPy, Numba, sympl, pytest, conda env `climt`

**Key simplification:** `NTRA=0` is hardcoded in `array_call`, so all tracer-related code (FTRA, TRAENT, TRAP) can be dropped from the JAX path.

---

## Algorithm Map

The full `_convect_functional_np` on develop has these blocks (line refs are develop `pure_python_v3.py`):

| # | Block | Lines | JAX difficulty | Status on jax-research |
|---|-------|-------|----------------|----------------------|
| 1 | Init + constants | 289–313 | trivial | done (stub) |
| 2 | Thermodynamic profiles (GZ, CPN, H, LV, HM, TV) | 314–349 | easy (scan for GZ, vectorize rest) | done (stub) |
| 3 | Launch level NK, IHMIN | 350–356 | easy (argmax with mask) | done (stub) |
| 4 | Early exits (T[NK], PLCL, ICB) | 357–374 | easy (masks instead of returns) | done (stub) |
| 5 | TLIFT (KK=1, KK=2) | 375–389 | medium (already ported) | done (stub) |
| 6 | EP/SIGP | 390–406 | easy (vectorized) | done (stub) |
| 7 | Setup matrices (QENT, SIJ, MENT, etc.) | 407–445 | easy (broadcast init) | **missing** |
| 8 | CAPE loop → INB | 446–464 | medium (scan) | done (stub) |
| 9 | HP update, CBMF | 465–483 | easy | done (stub) |
| 10 | Updraft mass flux M | 484–493 | easy (vectorized) | **missing** |
| 11 | Entrainment/detrainment (SIJ/MENT double loop) | 494–590 | **hard** (nested loop, conditionals) | **missing** |
| 12 | Downdraft (WDTRAIN, WATER, EVAP, MP, QP) | 591–656 | **hard** (reverse sequential scan) | **missing** |
| 13 | Tendency accumulation (FT, FQ, FU, FV) | 657–791 | medium (matrix sums with masks) | **missing** |
| 14 | FRAC smoothing at INB | 792–815 | easy (dynamic index) | **missing** |
| 15 | Conservation correction | 822–844 | easy (masked sums) | **missing** |

---

## Phase 1: Numpy Kernel Parity with Develop

### Task 1: Copy complete `_convect_functional_np` from develop

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

The current `_convect_functional_np` (lines 139–187) is a stub that returns zeros for FT/FQ. Replace it with the complete version from develop.

- [x] **Step 1: Replace the stub `_convect_functional_np`**

Delete lines 139–187 of `pure_python_v3.py` (the stub `_convect_functional_np`) and replace with the complete function from `develop:climt/_components/emanuel/pure_python_v3.py` lines 272–845. The function signature is identical. Keep the `@njit` decorator.

To get the code:
```bash
git show develop:climt/_components/emanuel/pure_python_v3.py | sed -n '272,845p'
```

Paste that in place of the stub. Keep the formatting from develop (it's more readable than the compressed style on jax-research).

- [x] **Step 2: Verify it runs**

```bash
conda activate climt && python -c "
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPythonV3
from climt import get_grid, get_default_state
from datetime import timedelta
import numpy as np
conv = EmanuelConvectionPythonV3()
grid = get_grid(nx=1, ny=1, nz=30)
state = get_default_state([conv], grid_state=grid)
state['air_temperature'].values[:] = 300.0
state['specific_humidity'].values[0:5,:,:] = 0.02
tend, diag = conv.array_call(
    {k: v.values if hasattr(v,'values') else v for k,v in state.items()},
    timedelta(minutes=15)
)
print('FT max:', np.max(np.abs(tend['air_temperature'])))
print('Precip:', diag['convective_precipitation_rate'])
"
```

Expected: Non-zero FT and possibly non-zero precipitation (depends on profile).

- [x] **Step 3: Commit**

```bash
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "feat(emanuel): port complete _convect_functional_np from develop"
```

---

### Task 2: Parity test — V3 numpy vs V1 reference

**Files:**
- Create: `tests/test_emanuel_v3_parity.py`

Port the parity test from develop. This compares V3 against the V1 pure-Python reference (`EmanuelConvectionPython`) to catch any transcription errors.

- [x] **Step 1: Write the parity test**

```python
# tests/test_emanuel_v3_parity.py
from datetime import timedelta
import numpy as np
import pytest
from climt import EmanuelConvection, get_default_state, get_grid
from climt._components.emanuel.pure_python import EmanuelConvectionPython
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPythonV3
from climt._core.util import numpy_version_of


def create_test_state(nlev, ncol, moisture_type="moist"):
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state["air_temperature"].values[:] = 300.0
    q = np.zeros_like(state["specific_humidity"].values)
    if moisture_type == "moist":
        q[0:5, :, :] = 0.020
    elif moisture_type == "dry":
        q[:] = 1e-6
    elif moisture_type == "unstable":
        for i in range(nlev):
            state["air_temperature"].values[i, :, :] = 310.0 - (i * 2.0)
        q[0:10, :, :] = 0.015
    state["specific_humidity"].values[:] = q
    return state


@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("moisture_type", ["moist", "dry", "unstable"])
def test_emanuel_v3_vs_v1_parity(ncol, moisture_type):
    nlev = 30
    state = create_test_state(nlev, ncol, moisture_type)
    timestep = timedelta(minutes=15)
    raw_state = numpy_version_of(state)
    python_state = {
        "air_temperature": raw_state["air_temperature"].reshape(nlev, ncol),
        "specific_humidity": raw_state["specific_humidity"].reshape(nlev, ncol),
        "air_pressure": (raw_state["air_pressure"] / 100.0).reshape(nlev, ncol),
        "air_pressure_on_interface_levels": (
            raw_state["air_pressure_on_interface_levels"] / 100.0
        ).reshape(nlev + 1, ncol),
        "eastward_wind": raw_state["eastward_wind"].reshape(nlev, ncol),
        "northward_wind": raw_state["northward_wind"].reshape(nlev, ncol),
    }
    conv_v1 = EmanuelConvectionPython()
    conv_v3 = EmanuelConvectionPythonV3()
    tend_v1, diag_v1 = conv_v1.array_call(python_state, timestep)
    tend_v3, diag_v3 = conv_v3.array_call(python_state, timestep)
    for key in tend_v1:
        np.testing.assert_allclose(
            tend_v1[key], tend_v3[key], atol=1e-3, rtol=1e-3,
            err_msg=f"Tendency '{key}' (ncol={ncol}, {moisture_type})",
        )
    for key in diag_v1:
        np.testing.assert_allclose(
            diag_v1[key], diag_v3[key], atol=1e-3, rtol=1e-3,
            err_msg=f"Diagnostic '{key}' (ncol={ncol}, {moisture_type})",
        )
```

- [x] **Step 2: Run test**

```bash
conda activate climt && python -m pytest tests/test_emanuel_v3_parity.py -v
```

Expected: All 6 parametrizations pass.

- [x] **Step 3: Commit**

```bash
git add tests/test_emanuel_v3_parity.py
git commit -m "test(emanuel): add V3 vs V1 parity tests"
```

---

## Phase 2: Full JAX Port (Non-Differentiable, Python Loops)

The goal here is a `_convect_functional_jax` that uses Python loops and produces the same numerical results as `_convect_functional_np`. It won't be JIT-able (Python control flow with traced values), but it establishes correctness.

### Task 3: Write `_convect_functional_jax` with Python loops

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

Replace the current stub `_convect_functional_jax` (lines 189–202) with a full implementation. This is a mechanical translation of `_convect_functional_np`: replace `np` → `jnp`, `max()` → `jnp.maximum()`, in-place assignment → `.at[].set()`, keep all the Python `for` loops and `if` statements as-is.

- [x] **Step 1: Write the full JAX-with-loops kernel**

Replace `_convect_functional_jax` with a new function. Key translation rules:
- `np.zeros(N)` → `jnp.zeros(N)`
- `arr[i] = val` → `arr = arr.at[i].set(val)`
- `arr[i] += val` → `arr = arr.at[i].add(val)`
- `arr[i, j] = val` → `arr = arr.at[i, j].set(val)`
- `max(a, b)` → `float(jnp.maximum(a, b))` (use `float()` to extract scalar for Python control flow)
- `min(a, b)` → `float(jnp.minimum(a, b))`
- `np.exp` → `jnp.exp`, `np.log` → `jnp.log`, `np.sqrt` → `jnp.sqrt`, `np.abs` → `jnp.abs`
- `arr.copy()` → `jnp.array(arr)` (JAX arrays are immutable, so copy is free)
- Early `return` statements remain as-is (Python control flow, not traced)
- Comparisons like `if T[NK] < 250.0` use `float(T[NK])` to get a concrete Python value

The function should call `_tlift_functional_jax` (already implemented in the stub) instead of `_tlift_functional_np`.

Important: since NTRA=0 always, skip all `for k in range(NTRA)` loops (they're zero-iteration). Still allocate FTRA as `jnp.zeros((ND, 1))` for return-signature compatibility, but don't compute it.

The function body follows `_convect_functional_np` from develop exactly — same variable names, same order — just with the jnp translations above. It should be ~300 lines (shorter than the numpy version because we skip NTRA loops).

- [x] **Step 2: Update `array_call` dispatch**

The existing `array_call` already dispatches based on `xp`. Ensure the JAX path calls this new function. The stub already does this via `jax.vmap` — keep that pattern but call the new complete function.

- [x] **Step 3: Smoke test**

```bash
conda activate climt && python -c "
import jax; jax.config.update('jax_enable_x64', True); jax.config.update('jax_platform_name', 'cpu')
import jax.numpy as jnp
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPythonV3
from climt import get_grid, get_default_state
from climt._core.util import numpy_version_of
from datetime import timedelta
import numpy as np
conv = EmanuelConvectionPythonV3()
grid = get_grid(nx=1, ny=1, nz=30)
state = get_default_state([conv], grid_state=grid)
state['air_temperature'].values[:] = 300.0
state['specific_humidity'].values[0:5,:,:] = 0.02
raw = numpy_version_of(state)
nlev, ncol = 30, 1
jax_state = {
    'air_temperature': jnp.array(raw['air_temperature'].reshape(nlev, ncol)),
    'specific_humidity': jnp.array(raw['specific_humidity'].reshape(nlev, ncol)),
    'air_pressure': jnp.array(raw['air_pressure'].reshape(nlev, ncol) / 100.0),
    'air_pressure_on_interface_levels': jnp.array(raw['air_pressure_on_interface_levels'].reshape(nlev+1, ncol) / 100.0),
    'eastward_wind': jnp.zeros((nlev, ncol)),
    'northward_wind': jnp.zeros((nlev, ncol)),
}
tend, diag = conv.array_call(jax_state, timedelta(minutes=15))
print('FT max:', jnp.max(jnp.abs(tend['air_temperature'])))
print('Precip:', diag['convective_precipitation_rate'])
"
```

Expected: Non-zero results matching numpy path.

- [ ] **Step 4: Commit**

```bash
git add climt/_components/emanuel/pure_python_v3.py
git commit -m "feat(emanuel): full JAX port with Python loops (non-JIT)"
```

---

### Task 4: Parity test — JAX-loops vs numpy

**Files:**
- Create: `tests/test_emanuel_jax_parity.py`

- [ ] **Step 1: Write the JAX vs numpy parity test**

```python
# tests/test_emanuel_jax_parity.py
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest
from datetime import timedelta
from climt import EmanuelConvection, get_default_state, get_grid
from climt._components.emanuel.pure_python_v3 import EmanuelConvectionPythonV3
from climt._core.util import numpy_version_of


def create_test_state(nlev, ncol, moisture_type="moist"):
    grid = get_grid(nx=ncol, ny=1, nz=nlev, p_surf_in_Pa=101325.0)
    state = get_default_state([EmanuelConvection()], grid_state=grid)
    state["air_temperature"].values[:] = 300.0
    q = np.zeros_like(state["specific_humidity"].values)
    if moisture_type == "moist":
        q[0:5, :, :] = 0.020
    elif moisture_type == "dry":
        q[:] = 1e-6
    elif moisture_type == "unstable":
        for i in range(nlev):
            state["air_temperature"].values[i, :, :] = 310.0 - (i * 2.0)
        q[0:10, :, :] = 0.015
    state["specific_humidity"].values[:] = q
    return state


def make_column_state(raw_state, nlev, ncol, as_jax=False):
    xp = jnp if as_jax else np
    return {
        "air_temperature": xp.array(raw_state["air_temperature"].reshape(nlev, ncol)),
        "specific_humidity": xp.array(raw_state["specific_humidity"].reshape(nlev, ncol)),
        "air_pressure": xp.array(raw_state["air_pressure"].reshape(nlev, ncol) / 100.0),
        "air_pressure_on_interface_levels": xp.array(
            raw_state["air_pressure_on_interface_levels"].reshape(nlev + 1, ncol) / 100.0
        ),
        "eastward_wind": xp.zeros((nlev, ncol)),
        "northward_wind": xp.zeros((nlev, ncol)),
    }


@pytest.mark.parametrize("ncol", [1, 4])
@pytest.mark.parametrize("moisture_type", ["moist", "dry", "unstable"])
def test_jax_vs_numpy_parity(ncol, moisture_type):
    nlev = 30
    state = create_test_state(nlev, ncol, moisture_type)
    timestep = timedelta(minutes=15)
    raw_state = numpy_version_of(state)
    conv = EmanuelConvectionPythonV3()

    np_state = make_column_state(raw_state, nlev, ncol, as_jax=False)
    jax_state = make_column_state(raw_state, nlev, ncol, as_jax=True)

    tend_np, diag_np = conv.array_call(np_state, timestep)
    tend_jax, diag_jax = conv.array_call(jax_state, timestep)

    for key in tend_np:
        np.testing.assert_allclose(
            np.asarray(tend_jax[key]), tend_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Tendency '{key}' (ncol={ncol}, {moisture_type})",
        )
    for key in ["convective_precipitation_rate", "cloud_base_mass_flux",
                 "atmosphere_convective_available_potential_energy"]:
        np.testing.assert_allclose(
            np.asarray(diag_jax[key]), diag_np[key], atol=1e-10, rtol=1e-10,
            err_msg=f"Diagnostic '{key}' (ncol={ncol}, {moisture_type})",
        )
```

- [ ] **Step 2: Run tests**

```bash
conda activate climt && python -m pytest tests/test_emanuel_jax_parity.py -v
```

Expected: All 6 pass at tight tolerance (1e-10) — the algorithm is identical, only the array library differs.

- [ ] **Step 3: Commit**

```bash
git add tests/test_emanuel_jax_parity.py
git commit -m "test(emanuel): JAX vs numpy parity tests"
```

---

## Phase 3: Differentiable JAX (JIT-Compatible)

The strategy: replace `_convect_functional_jax` with a version that uses no Python control flow on traced values. Convert block by block, testing parity after each.

Key JAX patterns used throughout:
- **Data-dependent branch** → `jnp.where(condition, true_val, false_val)` or mask arrays
- **Early return** → compute everything, multiply final result by `active` mask
- **Loop with dynamic bounds** `for i in range(ICB, INB)` → `jax.lax.scan` over all levels with `(i >= ICB) & (i < INB)` mask
- **Sequential accumulation** → `jax.lax.scan`
- **Double loop** (i over levels, j over levels) → 2D array operations with mask matrices
- **Dynamic indexing** `arr[NK]` → `arr[NK]` (works in JAX, just not JIT-able if NK is traced; use `jax.lax.dynamic_slice` if needed)

### Task 5: Differentiable updraft mass flux M (Block 10)

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

This is the easiest missing block. The numpy version:
```python
M[ICB] = 0.0
DBOSUM = 0.0
for i in range(ICB + 1, INB + 1):
    k = min(i, INB1)
    DBO = abs(TV[k] - TVP[k]) + params.ENTP * 0.02 * (PH[k] - PH[k + 1])
    DBOSUM += DBO
    M[i] = CBMF * DBO
if DBOSUM > 0:
    for i in range(ICB + 1, INB + 1):
        M[i] /= DBOSUM
```

- [x] **Step 1: Write the vectorized JAX equivalent**

```python
# Inside _convect_functional_jax, after CBMF computation:
lvl = jnp.arange(ND + 1)
k_idx = jnp.minimum(lvl, INB1)
DBO_all = jnp.abs(TV[k_idx] - TVP[k_idx]) + params.ENTP * 0.02 * (PH[k_idx] - PH[k_idx + 1])
mass_mask = (lvl >= ICB + 1) & (lvl <= INB)
DBO_masked = jnp.where(mass_mask, DBO_all, 0.0)
DBOSUM = jnp.sum(DBO_masked)
M = jnp.where(mass_mask, CBMF * DBO_masked / jnp.maximum(DBOSUM, 1e-30), 0.0)
```

- [x] **Step 2: Verify parity**

Run `tests/test_emanuel_jax_parity.py`. Expected: still passes. ✓ 6/6 passed.

- [ ] **Step 3: Commit**

```bash
git commit -am "feat(emanuel): differentiable updraft mass flux M"
```

---

### Task 6: Differentiable entrainment/detrainment — SIJ/MENT (Block 11)

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

This is the hardest block. The numpy version is a double loop `i in [ICB+1, INB+1]`, `j in [ICB, INB+1]` computing mixing fractions SIJ, entrainment MENT, and mixing properties QENT/UENT/VENT. It has data-dependent branching (recalculation when SIJ out of range).

The JAX approach: compute everything as 2D arrays indexed by `(i, j)`, using masks.

- [x] **Step 1: Write the 2D vectorized SIJ computation**

```python
# Create 2D index arrays
ii = jnp.arange(ND + 1)[:, None]  # (ND+1, 1) — updraft level
jj = jnp.arange(ND + 1)[None, :]  # (1, ND+1) — environment level

# Active mask: i in [ICB+1, INB], j in [ICB, INB]
i_active = (ii >= ICB + 1) & (ii <= INB)
j_active = (jj >= ICB) & (jj <= INB)
ij_active = i_active & j_active

# QTI for each updraft level i
QTI = Q[NK] - EP * CLW  # (ND,) — QTI[i]

# BF2 for each environment level j
BF2 = 1.0 + LV * LV * QS / (params.RV * T * T * params.CPD)  # (ND,)

# First-pass SIJ
# ANUM[i,j] = H[j] - HP[i] + (CPV - CPD) * T[j] * (QTI[i] - Q[j])
ANUM = H[None, :ND+1] - HP[:, None] + (params.CPV - params.CPD) * T[None, :ND+1] * (QTI[:, None] - Q[None, :ND+1])
DENOM_ij = H[:ND+1, None] - HP[:, None] + (params.CPD - params.CPV) * (Q[:ND+1, None] - QTI[:, None]) * T[None, :ND+1]
# Note: the index domains need care; H, Q etc. are (ND+1,) but HP, QTI are (ND,)
# This needs to be worked out carefully for the actual index ranges.
```

The full implementation of this block is complex. The key insight is that the conditional recalculation (when `SIJ < 0 or SIJ > 1 or ALTEM > CWAT` and `j > i`) can be computed for all (i,j) pairs and selected with `jnp.where`. Similarly, the SCRIT/SJMAX/SJMIN normalization loop can be expressed as array operations on the SIJ matrix.

This step requires careful index arithmetic. Write it to match the numpy kernel exactly, then verify with the parity test.

- [x] **Step 2: Write the SCRIT normalization (second half of block 11)**

The second half (develop lines 541–590) normalizes MENT using SCRIT, SJMAX, SJMIN. This involves:
1. Computing SCRIT for each level i
2. Computing ASIJ (weighted sum of SIJ differences) for each i
3. Normalizing MENT by ASIJ
4. Fallback: if sum of MENT for level i is tiny, set diagonal entry

This can be done as:
- SCRIT: vectorized over i
- SJMAX/SJMIN: use `jnp.roll` on SIJ to get SIJ[i, j+1] and SIJ[i, j-1]
- ASIJ: sum over j for each i
- Normalization: divide MENT by ASIJ with mask

- [x] **Step 3: Verify parity** ✓ 6/6 passed.

- [x] **Step 4: Commit** ✓ 8462894

---

### Task 7: Differentiable downdraft (Block 12)

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

The downdraft loop runs from INB down to 0 (reversed). Each level depends on the level above: `WATER[i]` depends on `WATER[i+1]`, `MP[i]` depends on `MP[i+1]`, `QP[i]` depends on `QP[i+1]`. This is a natural fit for `jax.lax.scan` running in reverse.

- [x] **Step 1: Write the reverse scan for downdraft** ✓ `jax.lax.scan` with 9-element carry (WATER, WT, MP, QP, UP, VP, JTT, MP_JTT, P_JTT). WDTRAIN pre-computed as vectorized 2D operation. BL correction carries MP_JTT/P_JTT through scan.

- [x] **Step 2: Compute PRECIP from WATER[0]** ✓ Integrated into scan output (WT_sc[INB], WATER_sc[INB] = level 0 values).

- [x] **Step 3: Verify parity** ✓ 6/6 passed.

- [x] **Step 4: Commit** ✓ da84ad1

---

### Task 8: Differentiable tendency accumulation (Blocks 13–15)

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`

Block 13 (FT, FQ, FU, FV computation) has the structure:
1. Surface level (i=0): direct formulas using AM, MP[1], EVAP[0], MENT sums
2. Interior levels (i=1 to INB): AMP1 (sum of M and MENT above), AD (sum of MENT below), then tendency formulas

Block 14: FRAC smoothing between INB and INB-1
Block 15: Conservation correction (ENTS, UAV, VAV averages subtracted)

- [x] **Step 1: Compute AM, AMP1, AD as prefix sums**

```python
# AM = sum of M[1:INB+1] (only when NK==0)
AM = jnp.where(NK == 0, jnp.sum(jnp.where((lvl >= 1) & (lvl <= INB), M, 0.0)), 0.0)

# AMP1[i] = sum of M[ii] for ii in [i+1, INB+1] (when ii >= NK)
#         + sum of MENT[k, jj] for k in [0, i], jj in [i+1, INB+1]
# This is a 1D array computed for each level i.

# Use cumulative sums:
M_cumsum_rev = jnp.cumsum(M[::-1])[::-1]  # M_cumsum_rev[i] = sum of M[i:]
# AMP1_mass[i] = sum of M[ii] for ii >= max(i+1, NK) and ii <= INB+1
# For the MENT part, sum MENT[k, jj] where k <= i and jj > i
# This requires a 2D prefix sum on MENT.

MENT_upper = jnp.triu(MENT, k=1)  # entries where j > i
# For each level i, sum over k <= i, j > i:
# This is cumsum along axis 0 of the column sums of the upper triangle... 
# Simpler: compute it as a matrix operation
# AMP1_ment[i] = sum_{k=0}^{i} sum_{jj=i+1}^{INB} MENT[k, jj]
```

The exact formulation needs careful linear-algebra equivalents of the nested sums. The pattern is: each sum over a triangular region of MENT can be expressed as a cumulative sum of row/column sums.

- [x] **Step 2: Compute FT, FQ, FU, FV for all levels**

Once AMP1 and AD are computed as arrays, the tendency formulas are element-wise:
```python
DPINV = 0.01 / (PH[:ND] - PH[1:ND+1])
CPINV = 1.0 / CPN[:ND]

# Interior levels (i >= 1):
FT_interior = (
    params.G * DPINV * (
        AMP1 * (T_shifted_up - T) + ... # using jnp.roll or shifted arrays
    )
    - params.SIGD * LVCP * EVAP
    + ...  # MENT[i,i] self-entrainment term
    + ...  # water loading term
)
# Plus MENT contributions: sum over k of MENT[k,i] * (QENT[k,i] - Q[i])
# These are column sums of MENT * (QENT - Q) matrices
```

- [x] **Step 3: FRAC smoothing (Block 14)**

```python
# Smooth between INB and INB-1
ph_ratio = (PH[INB] - PH[INB+1]) / (PH[INB-1] - PH[INB])
FQ = FQ.at[INB-1].add(FRAC * FQ[INB] * ph_ratio * LV[INB] / LV[INB-1])
FQ = FQ.at[INB].multiply(1.0 - FRAC)
# Same for FT, FU, FV
```

Note: Order matters — save old values before modifying:
```python
FQOLD = FQ[INB]
FQ = FQ.at[INB].set(FQOLD * (1.0 - FRAC))
FQ = FQ.at[INB-1].add(FRAC * FQOLD * ph_ratio * LV[INB] / LV[INB-1])
```

- [x] **Step 4: Conservation correction (Block 15)**

```python
active_levels = lvl[:ND] <= INB
dp = PH[:ND] - PH[1:ND+1]
DENOM2 = jnp.sum(jnp.where(active_levels, dp, 0.0))
ENTS = jnp.sum(jnp.where(active_levels, (CPN[:ND] * FT + LV[:ND] * FQ) * dp, 0.0)) / DENOM2
UAV = jnp.sum(jnp.where(active_levels, FU * dp, 0.0)) / DENOM2
VAV = jnp.sum(jnp.where(active_levels, FV * dp, 0.0)) / DENOM2
FT = jnp.where(active_levels, FT - ENTS / CPN[:ND], FT)
FU = jnp.where(active_levels, (1.0 - params.CU) * (FU - UAV), FU)
FV = jnp.where(active_levels, (1.0 - params.CU) * (FV - VAV), FV)
```

- [x] **Step 5: Apply the `active` mask to final outputs** (N/A — early returns handle this)

All tendency outputs should be zeroed when `active` is False (no convection):
```python
FT = jnp.where(active & (CBMF > 0), FT, 0.0)
FQ = jnp.where(active & (CBMF > 0), FQ, 0.0)
FU = jnp.where(active & (CBMF > 0), FU, 0.0)
FV = jnp.where(active & (CBMF > 0), FV, 0.0)
PRECIP = jnp.where(active & (CBMF > 0), PRECIP, 0.0)
```

- [x] **Step 6: Verify parity**

```bash
conda activate climt && python -m pytest tests/test_emanuel_jax_parity.py -v
```

- [x] **Step 7: Commit**

```bash
git commit -am "feat(emanuel): differentiable tendency accumulation and conservation"
```

---

### Task 9: JIT + vmap integration

**Files:**
- Modify: `climt/_components/emanuel/pure_python_v3.py`
- Modify: `tests/test_emanuel_jax_parity.py`

**Overview:** The current `_convect_functional_jax` (lines 967–1664) and `_tlift_functional_jax` (lines 914–964) use Python `float()` extraction, Python `if/else` on traced values, and data-dependent loop bounds (`ICB`, `NK`, `INB`). These all break `jax.jit` tracing. This task rewrites blocks 2–9 and `_tlift_functional_jax` to be trace-safe, then wraps with `@jax.jit` and replaces the column loop with `jax.vmap`.

**Strategy:**
- `ND`, `NL`, `NTRA` are passed as `static_argnums` — loops over `range(NL)` are fine.
- `ICB`, `NK`, `INB`, `IHMIN` are data-dependent integers used as indices → keep as JAX integers, use `jnp.where` masks instead of `range(ICB, INB)`.
- Early returns → accumulate an `active` flag, compute everything unconditionally, multiply results by `active` at the end.
- Python `float()` on traced values → remove (JAX scalar ops work directly).
- Python `max()`/`min()` on traced values → `jnp.maximum()`/`jnp.minimum()`.
- Python `if cond:` on traced values → `jnp.where(cond, ...)`.

---

#### Step 9a: Remove all `float()` calls from `_convect_functional_jax`

~50 `float()` calls extract concrete values from traced arrays. Remove them all — JAX scalar indexing (`T[i]`, `Q[NK]`) returns JAX scalars that work in arithmetic.

Also replace Python `max()`/`min()` with `jnp.maximum()`/`jnp.minimum()` where arguments may be traced.

- [x] Remove `float()` calls in blocks 2–9 (lines 978–1205)
- [x] Replace `max()`/`min()` with `jnp.maximum()`/`jnp.minimum()` on traced values
- [x] Remove `float()` calls in blocks 13–15 (lines 1504, 1524, 1572)

---

#### Step 9b: Make `_tlift_functional_jax` JIT-compatible

Current issues (lines 914–964):
- `float(Q[NK])`, `float(T[NK])`, `float(GZ[NK])` — remove `float()`
- `if KK == 1:` — `KK` is always a Python literal (called with `1` or `2`), so this is fine if KK is static
- `for i in range(ICB):` and `for i in range(NK, ICB):` — `ICB`/`NK` are traced → iterate over `range(NL+1)` with `(i < ICB)` mask
- `for i in range(NSB, NST + 1):` — `NSB`/`NST` depend on `ICB` → full range + mask
- `if TC >= 0.0:` on traced value → `jnp.where`
- `max(TG + ..., 35.0)` → `jnp.maximum`
- `max(0.0, qNK - QG)` → `jnp.maximum`

- [x] Rewrite `_tlift_functional_jax` with full-range loops and masks
- [x] Replace Python `if TC >= 0.0` with `jnp.where` for ES computation
- [x] Make `KK` a static argument (it's always literal 1 or 2)

---

#### Step 9c: Rewrite Block 2 (thermodynamic profiles) for JIT

Current: Python loop `for i in range(NL+1)` with `float()` — loop bound is static so only need to remove `float()`.

- [x] Remove `float()` calls (already done in 9a, verify correctness) ✓ verified

---

#### Step 9d: Rewrite Block 3 (IHMIN/NK) for JIT

Current (lines 1024–1060): Two Python loops that compute `IHMIN` (argmin of HM with conditions) and `NK` (argmax of HM in a subrange), using Python `if` on traced values.

Approach:
- Vectorize IHMIN: compute a masked array where invalid entries are `+inf`, then `jnp.argmin`
- Vectorize NK: compute a masked array where invalid entries are `-inf`, then `jnp.argmax`

- [x] Replace IHMIN loop with `jnp.argmin` + mask
- [x] Replace NK loop with `jnp.argmax` + mask

---

#### Step 9e: Rewrite Block 4 (early exits) for JIT

Current (lines 1062–1084): Three `if`/`return` blocks based on traced conditions (`T[NK] < 250`, `PLCL` range, `ICB >= NL-1`).

Approach: Accumulate a boolean `active` flag. Compute everything unconditionally. At the end (before return), select between zeros and computed values based on `active`.

- `ICB` computation (lines 1076–1079): Python loop with `if float(P[i]) < PLCL` → vectorize with `jnp.argmax(P < PLCL)` or similar.

- [x] Convert ICB search to vectorized `jnp.where` + index computation
- [x] Replace early returns with `active` flag accumulation
- [x] Track `IFLAG` values for each exit condition

---

#### Step 9f: Rewrite Block 5 (TLIFT calls) + Block 6 (EP/SIGP) for JIT

Block 5 (lines 1086–1104): Calls `_tlift_functional_jax` (fixed in 9b), then loops `for i in range(NK, ICB+1)` and `for i in range(ICB+1, NL+1)` with traced bounds.

Block 6 (lines 1106–1125): Loops with traced bounds (`NK`, `ICB`, `NL`), Python `if` on TCA.

- [x] Convert variable-bound loops to full-range loops with `jnp.where` masks
- [x] Replace `if TCA >= 0.0` / `max`/`min` with `jnp.where`/`jnp.maximum`/`jnp.minimum`
- [x] Replace `if CBMF == 0.0 and TVP[ICB] <= ...` early return with `active` flag

---

#### Step 9g: Rewrite Block 7 (setup matrices) for JIT

Current (lines 1126–1159): Loops with traced bounds to fill QENT, UENT, VENT, QP, UP, VP.

Approach: Use broadcasting. `QENT[i,j] = Q[j]` for all i,j → `QENT = jnp.broadcast_to(Q[None,:], ...)` etc.

- [x] Replace matrix-fill loops with broadcasting
- [x] Replace QP/UP/VP init loops with array ops

---

#### Step 9h: Rewrite Block 8 (CAPE → INB) for JIT

Current (lines 1161–1180): Loop `for i in range(ICB+1, NL)` computing CAPE, INB, BYP with Python `if` on traced values.

Approach: Iterate over full range with masks. Use `jnp.where` for conditional updates to INB, CAPE, etc.

- [x] Convert to full-range loop (0 to NL) with `(i > ICB) & (i < NL)` mask
- [x] Replace Python `if BY >= 0` / `if CAPE > 0` with `jnp.where`

---

#### Step 9i: Rewrite Block 9 (HP update, CBMF) for JIT

Current (lines 1182–1205): Loop with traced bounds, final early return `if CBMF == 0.0 and CBMFOLD == 0.0`.

- [x] Convert HP update loop to full-range + mask
- [x] Replace TVPPLCL/TVAPLCL/DTPBL computations (use dynamic indexing, already works)
- [x] Replace final early return with `active` flag update

---

#### Step 9j: Fix remaining JIT issues in Blocks 10–15

Blocks 10–15 are already differentiable but still have a few issues:
- Line 1304–1313: Python `for` loop with `if` on `NENT_count[i]` (traced) in Block 11
- Line 1319–1378: Python `for` loop with `if not has_ent` (traced) — Block 11 SCRIT normalization
- Line 1335–1358: Nested Python `for` with `if j > i` — Block 11 inner loop
- Line 1475–1476: `if ep_active` where `EP[INB]` is traced — Block 12 entry
- Line 1493–1501: Python `for` scattering scan results — Block 12
- Line 1524, 1572: `if` on traced value for IFLAG=4 — Block 13

- [x] Replace Block 11 fallback loop (lines 1304–1313) with vectorized `jnp.where`
- [x] Replace Block 11 SCRIT normalization loop with vectorized or `lax.scan` approach
- [x] Replace Block 12 `if ep_active` with unconditional compute + `jnp.where`
- [x] Replace Block 12 scatter loop with array indexing
- [ ] Replace Block 13 IFLAG checks with `jnp.where`

---

#### Step 9k: Add `@jax.jit` and test

- [x] Add `@functools.partial(jax.jit, static_argnums=(7, 8, 9))` to `_convect_functional_jax`
- [x] Run single-column JIT test to verify no tracing errors
- [x] Run parity test: `conda activate climt && python -m pytest tests/test_emanuel_jax_parity.py -v`

---

#### Step 9l: Replace column loop with `jax.vmap`

Replace `_jax_vectorized_convect` (lines 1667–1702) which loops over columns with a `jax.vmap` call:

```python
def _jax_vectorized_convect(t, q, qs, u, v, p, ph, ND, NL, NTRA, DELT, cbmf, params):
    def column_call(T, Q, QS, U, V, P, PH, CBMF):
        return _convect_functional_jax(T, Q, QS, U, V, P, PH, ND, NL, NTRA, DELT, CBMF, params)
    results = jax.vmap(column_call, in_axes=(1, 1, 1, 1, 1, 1, 1, 0))(
        t, q, qs, u, v, p, ph, cbmf)
    # vmap returns (ncol, nlev) for 1D outputs → transpose to (nlev, ncol)
    ft, fq, fu, fv = results[0].T, results[1].T, results[2].T, results[3].T
    precip, wd, tprime, qprime, cbmf_new, outcape, iflag = results[4:]
    return ft, fq, fu, fv, precip, wd, tprime, qprime, cbmf_new, outcape, iflag
```

- [x] Replace `_jax_vectorized_convect` with `vmap` version
- [x] Run parity test
- [x] Commit: `feat(emanuel): JIT + vmap integration for JAX kernel`

---

### Task 10: Gradient test

**Files:**
- Modify: `tests/test_jax_differentiation.py`

- [x] **Step 1: Update the existing gradient test**

The existing test in `tests/test_jax_differentiation.py` computes `jax.grad` of a loss over temperature tendencies. Now that the kernel is complete, it should produce non-zero gradients for unstable profiles.

```python
def test_emanuel_jax_gradient():
    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_platform_name", "cpu")
    nlev, ncol = 30, 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    conv = EmanuelConvectionPythonV3()
    state = get_default_state([conv], grid_state=grid)

    # Unstable profile
    temp = jnp.linspace(310, 250, nlev).reshape(nlev, ncol)
    q = jnp.zeros((nlev, ncol)).at[0:10].set(0.015)
    jax_state = {
        "air_temperature": temp,
        "specific_humidity": q,
        "air_pressure": jnp.array(state["air_pressure"].values / 100.0).reshape(nlev, ncol),
        "air_pressure_on_interface_levels": jnp.array(
            state["air_pressure_on_interface_levels"].values / 100.0
        ).reshape(nlev + 1, ncol),
        "eastward_wind": jnp.zeros((nlev, ncol)),
        "northward_wind": jnp.zeros((nlev, ncol)),
    }

    def compute_loss(temp_values):
        s = {**jax_state, "air_temperature": temp_values}
        tend, _ = conv.array_call(s, timedelta(minutes=15))
        return jnp.sum(jnp.square(tend["air_temperature"]))

    grad_fn = jax.grad(compute_loss)
    gradient = grad_fn(temp)
    assert jnp.any(gradient != 0.0), "Gradient is all zeros"
    assert not jnp.any(jnp.isnan(gradient)), "Gradient contains NaNs"
```

- [x] **Step 2: Run gradient test**

```bash
conda activate climt && python -m pytest tests/test_jax_differentiation.py -v
```

- [x] **Step 3: Commit**

```bash
git commit -am "test(emanuel): gradient test with complete differentiable kernel"
```

---

## Notes for Implementation

1. **NTRA simplification**: Since `NTRA=0`, all `for k in range(NTRA)` loops are zero-iteration. Don't port FTRA, TRAENT, TRAP to JAX. Allocate dummy zeros for return compatibility.

2. **Index domains**: The numpy kernel uses arrays sized `ND+1` (e.g., GZ, PH, H) and `ND` (e.g., T, Q, P, FT). Keep this convention. `ND = nlev`, `NL = nlev - 3`.

3. **The `_tlift_functional_jax` already exists** in the stub and handles KK=1 and KK=2 correctly. Reuse it.

4. **Debugging tip**: When a block fails parity, test it in isolation. Extract the intermediate arrays (GZ, TV, TVP, M, SIJ, etc.) from both numpy and JAX paths and compare them element-by-element to find where they diverge.

5. **`lax.scan` indexing**: When scanning in reverse (downdraft), the scan index `i_rev` runs 0, 1, 2, ... and the actual level is `INB - i_rev`. The carry propagates downward.

6. **Dynamic indexing**: Operations like `FT[INB]` where INB is a traced integer work in JAX via `FT[INB]` for reads and `FT.at[INB].set(val)` for writes. These are differentiable.
