# Cork Longwave Ozone (O₃) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add longwave ozone absorption to `CorkLongwaveRadiation` as a dedicated additive gas with a precomputed `k_O3(T,p)` cross-section and a runtime-variable O₃ profile, so CORK can reproduce a realistic (warm) stratosphere.

**Architecture:** O₃ becomes a second gas in the Earth LW correlated-k table alongside the existing premixed-air gas. Its cross-section is stored per-kg-O₃ and broadcast-constant across the table's H₂O VMR axis, so `interpolate_k`'s existing 6-D `ngas`-loop returns it unchanged. Optical depth is `τ = k_air·m_air + k_O3·m_O3` via the existing additive-overlap path. The O₃ amount comes from `mole_fraction_of_ozone_in_air`, which climt already ships via `init_ozone` — the same field RRTMG consumes, so validation is apples-to-apples.

**Tech Stack:** numpy, scipy.io.netcdf, climt (`climt` conda env), linepyline (`radiation`/`linepyline` env, table generation only), pytest, matplotlib.

**Design reference:** `docs/superpowers/specs/2026-06-01-cork-ozone-lw-design.md`

**Conda envs:** Use `climt` (`/Users/joymonteiro/miniconda3/envs/climt/bin/python`) for all code/tests. Use the `radiation`/`linepyline` env ONLY for the table-generation task (Task 8). See memory `feedback_conda_env.md`.

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `scripts/cork_table_builder/kappa_sampling.py` | line-by-line κ sampling | **add** `sample_secondary_gas_kappa()` — raw per-kg-species κ(T,p,ν) |
| `scripts/cork_table_builder/k_distribution.py` | κ → k-distribution | **add** `broadcast_and_stack_gases()` — pure helper, no linepyline |
| `scripts/cork_table_builder/netcdf_writer.py` | write `.nc` schema | **add** `background_is_premixed` attribute to `write_lw_table` |
| `scripts/generate_cork_tables_linepyline.py` | scenario registry + `build_table` | **add** `lw_secondary_gases` to `earth`; stack O₃ gas in `build_table` |
| `climt/_components/cork/lw/component.py` | LW component | **modify** `input_properties` + `array_call` premixed-bg branch for secondary gases |
| `tests/cork_table_builder/test_secondary_gas.py` | builder unit tests | **create** |
| `tests/test_cork_optics.py` | optics unit tests | **add** 6-D × ngas=2 additive regression test |
| `tests/test_cork_lw.py` | component tests | **add** O₃ on/off behaviour tests |
| `climt/_data/cork/correlated_k/earth_low_res_lw_co2refined_o3.nc` | shipped table | **create** (Task 8, linepyline env) |
| `scripts/experiments/radiative_equilibrium_cork_o3.py` | RRTMG validation | **create** (Task 9) |

---

## Task 1: Writer — `background_is_premixed` attribute

The component decides premixed-background mode from `gas_names == ["effective"]` **or** a `background_is_premixed="true"` netCDF attribute (`climt/_components/cork/lw/component.py:47-51`). A 2-gas table (`["effective","o3"]`) fails the first test, so the writer must emit the attribute. The writer (`scripts/cork_table_builder/netcdf_writer.py:14`) does not currently write it.

**Files:**
- Modify: `scripts/cork_table_builder/netcdf_writer.py:14-84`
- Test: `tests/cork_table_builder/test_netcdf_writer.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/cork_table_builder/test_netcdf_writer.py`:

```python
def test_write_lw_table_background_is_premixed_attr(tmp_path):
    import numpy as np
    from scipy.io import netcdf_file
    from cork_table_builder.netcdf_writer import write_lw_table

    out = str(tmp_path / "bg.nc")
    nband, ngpt, nT, nP, nX = 2, 2, 3, 3, 4
    write_lw_table(
        out,
        k_coefficients=np.ones((1, nband, ngpt, nT, nP, nX)) * 0.01,
        gpoint_weights=np.ones((nband, ngpt)) * 0.5,
        T_grid=np.linspace(200, 300, nT),
        log_p_grid=np.linspace(0, 11, nP),
        band_edges=np.array([10.0, 800.0, 3250.0]),
        planck_fraction=np.ones((nband, ngpt, nT)) * 0.5,
        h2o_vmr_grid=np.array([1e-6, 1e-4, 1e-2, 1.0]),
        gas_names=("effective", "o3"),
        background_is_premixed=True,
    )
    with netcdf_file(out, "r", mmap=False) as nc:
        assert nc.gas_names.decode() == "effective,o3"
        assert nc.background_is_premixed.decode().lower() == "true"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_netcdf_writer.py::test_write_lw_table_background_is_premixed_attr -v`
Expected: FAIL — `write_lw_table() got an unexpected keyword argument 'background_is_premixed'`.

- [ ] **Step 3: Add the parameter and attribute**

In `scripts/cork_table_builder/netcdf_writer.py`, add `background_is_premixed=False,` to the `write_lw_table` signature (after `gas_names=("effective",),`). Then, immediately after the line `nc.gas_names = ",".join(gas_names)`:

```python
        nc.gas_names = ",".join(gas_names)
        nc.background_is_premixed = "true" if background_is_premixed else "false"
        nc.overlap_method = overlap_method
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_netcdf_writer.py::test_write_lw_table_background_is_premixed_attr -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/cork_table_builder/netcdf_writer.py tests/cork_table_builder/test_netcdf_writer.py
git commit -m "feat(cork-tables): write background_is_premixed attr in LW table"
```

---

## Task 2: Builder — `broadcast_and_stack_gases` pure helper

Combine a primary 6-D `k` array `(1, nband, ngpt, nT, nP, nX)` with one or more secondary 4-D gas arrays `(nband, ngpt, nT, nP)` by broadcasting each secondary gas constant across the `nX` H₂O axis and concatenating along the gas axis. Pure numpy — fully unit-testable without linepyline.

**Files:**
- Modify: `scripts/cork_table_builder/k_distribution.py` (append function)
- Test: `tests/cork_table_builder/test_secondary_gas.py` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/cork_table_builder/test_secondary_gas.py`:

```python
import numpy as np
from cork_table_builder.k_distribution import broadcast_and_stack_gases


def test_broadcast_and_stack_gases_shape_and_xconstant():
    nband, ngpt, nT, nP, nX = 2, 3, 4, 5, 7
    primary = np.random.rand(1, nband, ngpt, nT, nP, nX)
    secondary = np.random.rand(nband, ngpt, nT, nP)  # one O3 gas, no X axis

    stacked = broadcast_and_stack_gases(primary, [secondary])

    # gas axis grew from 1 to 2
    assert stacked.shape == (2, nband, ngpt, nT, nP, nX)
    # gas 0 unchanged
    np.testing.assert_array_equal(stacked[0], primary[0])
    # gas 1 is constant along the X axis and equals the secondary at every X
    for iX in range(nX):
        np.testing.assert_array_equal(stacked[1, ..., iX], secondary)


def test_broadcast_and_stack_gases_no_x_axis():
    nband, ngpt, nT, nP = 2, 3, 4, 5
    primary = np.random.rand(1, nband, ngpt, nT, nP)  # 5-D, no X
    secondary = np.random.rand(nband, ngpt, nT, nP)
    stacked = broadcast_and_stack_gases(primary, [secondary])
    assert stacked.shape == (2, nband, ngpt, nT, nP)
    np.testing.assert_array_equal(stacked[1], secondary)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_secondary_gas.py -v`
Expected: FAIL — `ImportError: cannot import name 'broadcast_and_stack_gases'`.

- [ ] **Step 3: Implement the helper**

Append to `scripts/cork_table_builder/k_distribution.py`:

```python
def broadcast_and_stack_gases(primary_k, secondary_k_list):
    """Stack secondary additive gases onto a primary k-coefficient array.

    Args:
        primary_k: (1, nband, ngpt, nT, nP[, nX]) — the premixed primary gas,
            with the leading ngas=1 axis already present.
        secondary_k_list: list of (nband, ngpt, nT, nP) arrays. Each is a
            per-kg-species cross-section with NO H2O X-axis dependence; it is
            broadcast constant across nX when the primary has an X axis.

    Returns:
        (ngas, nband, ngpt, nT, nP[, nX]) with ngas = 1 + len(secondary_k_list).
    """
    has_x = primary_k.ndim == 6
    gases = [primary_k[0]]  # strip the leading axis; re-add via stack below
    for sec in secondary_k_list:
        if has_x:
            nX = primary_k.shape[-1]
            sec = np.broadcast_to(sec[..., None], sec.shape + (nX,))
        gases.append(np.asarray(sec))
    return np.stack(gases, axis=0)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_secondary_gas.py -v`
Expected: PASS (both tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/cork_table_builder/k_distribution.py tests/cork_table_builder/test_secondary_gas.py
git commit -m "feat(cork-tables): add broadcast_and_stack_gases helper"
```

---

## Task 3: Builder — `sample_secondary_gas_kappa`

Sample a single gas's air-broadened cross-section **per kg of that gas** (no mass-fraction weighting, no H₂O axis). This is the minor-gas treatment: foreign (air) broadening dominates for a trace gas, so the self partial pressure is set to ~0. Output units: m²/kg-species.

**Files:**
- Modify: `scripts/cork_table_builder/kappa_sampling.py` (append function)
- Test: `tests/cork_table_builder/test_secondary_gas.py`

- [ ] **Step 1: Write the failing test (linepyline-gated)**

Append to `tests/cork_table_builder/test_secondary_gas.py`:

```python
import pytest


def test_sample_secondary_gas_kappa_o3_shape_and_positive():
    pytest.importorskip("linepyline")
    from cork_table_builder.kappa_sampling import sample_secondary_gas_kappa

    T_grid = np.array([220.0, 260.0])
    p_grid = np.array([1e3, 1e4])
    # 9.6 um O3 band ~ 980-1080 cm^-1
    nu_grid = np.arange(950.0, 1100.0 + 0.05, 0.1)

    kappa = sample_secondary_gas_kappa(
        species="O3",
        background_gas="air",
        T_grid=T_grid,
        p_grid=p_grid,
        nu_grid=nu_grid,
        line_shape="pseudovoigt",
        binning=False,
    )
    assert kappa.shape == (len(T_grid), len(p_grid), len(nu_grid))
    assert np.all(kappa >= 0.0)
    # O3 absorbs somewhere in its 9.6 um band
    assert kappa.max() > 0.0
```

- [ ] **Step 2: Run test to verify it fails (or skips without linepyline)**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_secondary_gas.py::test_sample_secondary_gas_kappa_o3_shape_and_positive -v`
Expected: FAIL — `ImportError: cannot import name 'sample_secondary_gas_kappa'` (in `climt` env where linepyline IS installed for radiation work it will reach the call; if linepyline is absent it SKIPS). Either way, not PASS yet.

- [ ] **Step 3: Implement the sampler**

Append to `scripts/cork_table_builder/kappa_sampling.py` (it already imports numpy as np at module top and imports `linepyline` lazily inside `sample_kappa_grid`):

```python
def sample_secondary_gas_kappa(
    *,
    species,
    background_gas,
    T_grid,
    p_grid,
    nu_grid,
    line_shape="pseudovoigt",
    binning=False,
    surface_gravity=9.81,
):
    """Air-broadened cross-section of a minor gas, per kg of that gas.

    Unlike ``sample_kappa_grid`` (which premixes absorbers into a per-kg-air
    kappa via mass fractions), this returns the RAW per-kg-species kappa so the
    gas can be carried as a separate additive gas whose runtime amount is its
    own column mass. Self partial pressure is ~0 (foreign/air broadening),
    appropriate for a trace gas. No H2O axis.

    Returns:
        kappa: (nT, nP, nNu), units m^2/kg of `species`.
    """
    import linepyline

    rtm = linepyline.rtm(background_gas=background_gas,
                         surface_gravity=surface_gravity)
    nu_min, nu_max = float(nu_grid[0]), float(nu_grid[-1])
    dnu = float(np.diff(nu_grid).mean())
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)

    kappa = np.zeros((nT, nP, nNu))
    for iT, Tval in enumerate(T_grid):
        T_arr = np.full(nP, Tval)
        k_da = rtm.get_kappa_hitran(
            species, nu_min, nu_max, dnu,
            p_grid, T_arr, p_self=np.zeros(nP),
            line_shape=line_shape, binning=binning,
            remove_plinth=False,
        )
        kappa[iT] = np.maximum(
            k_da.transpose("p", "nu").interp(nu=nu_grid).values, 0.0
        )
    return kappa
```

- [ ] **Step 4: Run test to verify it passes (or skips)**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_secondary_gas.py::test_sample_secondary_gas_kappa_o3_shape_and_positive -v`
Expected: PASS if linepyline is importable in the `climt` env; SKIP otherwise. (Full numeric validation happens in Task 8 under the `radiation` env.)

- [ ] **Step 5: Commit**

```bash
git add scripts/cork_table_builder/kappa_sampling.py tests/cork_table_builder/test_secondary_gas.py
git commit -m "feat(cork-tables): add sample_secondary_gas_kappa (per-kg minor gas)"
```

---

## Task 4: Builder — wire O₃ into the `earth` scenario and `build_table`

Add a `lw_secondary_gases` config key and stack secondary-gas k-distributions during LW table builds. SW is out of scope (O₃ UV bands), so only the LW branch consumes it.

**Files:**
- Modify: `scripts/generate_cork_tables_linepyline.py:85-100` (earth scenario) and `:123-194` (`build_table`)

- [ ] **Step 1: Add the config key to the `earth` scenario**

In `scripts/generate_cork_tables_linepyline.py`, inside the `"earth": dict(...)` block, add after the `cia_files=[],` line:

```python
        lw_secondary_gases=["O3"],   # additive minor gas; runtime profile from init_ozone
```

- [ ] **Step 2: Stack secondary gases in `build_table` (LW only)**

In `build_table`, the imports block at the top of the file already pulls from `cork_table_builder.k_distribution import kappa_to_k_coeffs` and `cork_table_builder.kappa_sampling import sample_kappa_grid`. Extend those imports to also bring in the new helpers:

```python
from cork_table_builder.kappa_sampling import sample_kappa_grid, sample_secondary_gas_kappa
from cork_table_builder.k_distribution import kappa_to_k_coeffs, broadcast_and_stack_gases
```

Then, in `build_table`, immediately AFTER the existing block that produces `k_coeffs` with the gas axis and rearranges it (the lines ending at `k_coeffs = np.moveaxis(...)`, i.e. just before `log_p_grid = np.log(cfg["p_grid"])`), insert:

```python
    # Secondary additive gases (LW only): each carried per-kg-species, broadcast
    # constant across the H2O X-axis. Runtime amount comes from its own profile.
    gas_names = ["effective"]
    secondary = cfg.get("lw_secondary_gases", []) if kind == "lw" else []
    if secondary:
        sec_k_list = []
        for sp in secondary:
            sec_kappa = sample_secondary_gas_kappa(
                species=sp,
                background_gas=cfg["background_gas"],
                T_grid=cfg["T_grid"],
                p_grid=cfg["p_grid"],
                nu_grid=nu_grid,
                line_shape=line_shape,
                binning=binning,
            )
            sec_k, _ = kappa_to_k_coeffs(sec_kappa, nu_grid, band_edges, ngpt)
            # (nT, nP, nband, ngpt) -> (nband, ngpt, nT, nP)
            sec_k = np.moveaxis(sec_k, (0, 1, 2, 3), (2, 3, 0, 1))
            sec_k_list.append(sec_k)
            gas_names.append(sp.lower())
        k_coeffs = broadcast_and_stack_gases(k_coeffs, sec_k_list)
```

- [ ] **Step 3: Pass `gas_names` and `background_is_premixed` to the writer**

In the `if kind == "lw":` branch of `build_table`, update the `write_lw_table(...)` call to pass the new arguments (add these two keyword args; keep the existing ones):

```python
        write_lw_table(
            output_path,
            k_coefficients=k_coeffs,
            gpoint_weights=gpt_weights,
            T_grid=cfg["T_grid"],
            log_p_grid=log_p_grid,
            band_edges=band_edges,
            planck_fraction=planck_fraction,
            h2o_vmr_grid=cfg["h2o_vmr_grid"],
            gas_names=tuple(gas_names),
            background_is_premixed=bool(secondary),
            source=f"linepyline:{scenario_name}",
        )
```

- [ ] **Step 4: Syntax/import check (no linepyline call)**

Run: `conda run -n climt python -c "import ast; ast.parse(open('scripts/generate_cork_tables_linepyline.py').read()); print('ok')"`
Expected: `ok`.

- [ ] **Step 5: Commit**

```bash
git add scripts/generate_cork_tables_linepyline.py
git commit -m "feat(cork-tables): build O3 as additive secondary LW gas in earth scenario"
```

---

## Task 5: Optics — regression test for the 6-D × ngas=2 additive path

No shipped table has ever combined the H₂O X-axis (6-D `k`) with `ngas>1`. The code paths (`interpolate_k`'s 6-D `ngas`-loop and `_compute_ck_optical_depth_additive`) already support it, but it is untested. Pin it. **No production code change is expected; if the test fails, fix the optics code.**

**Files:**
- Test: `tests/test_cork_optics.py` (add)

- [ ] **Step 1: Write the test**

Add to `tests/test_cork_optics.py`:

```python
def test_interpolate_k_6d_two_gas_additive():
    """6-D H2O-axis table with a second X-constant gas: gas1 is X-independent
    and additive overlap sums both gases."""
    import numpy as np
    from climt._components.cork.optics.correlated_k import (
        interpolate_k, compute_ck_optical_depth,
    )

    ngas, nband, ngpt, nT, nP, nX = 2, 1, 1, 3, 3, 4
    k = np.zeros((ngas, nband, ngpt, nT, nP, nX))
    # gas 0 (air): varies along X (so X interpolation is exercised)
    k[0] = np.linspace(1.0, 4.0, nX).reshape(1, 1, 1, 1, nX)
    # gas 1 (o3): constant 7.0 across ALL axes incl. X
    k[1] = 7.0
    table = {
        "k_coefficients": k,
        "temperature_grid": np.linspace(200.0, 300.0, nT),
        "pressure_grid_log": np.linspace(0.0, 11.0, nP),
        "h2o_vmr_grid": np.array([1e-6, 1e-4, 1e-2, 1.0]),
        "gpoint_weights": np.ones((nband, ngpt)),
        "overlap_method": np.array("additive"),
    }

    T = np.array([250.0, 250.0])
    p = np.array([1e3, 1e4])
    # two columns with very different H2O -> gas0 differs, gas1 must not
    h2o = np.array([1e-6, 1.0])
    ki = interpolate_k(table, T, p, h2o_vmr=h2o)  # (ngas, nband, ngpt, ncol)
    assert ki[1, 0, 0, 0] == 7.0 and ki[1, 0, 0, 1] == 7.0  # gas1 X-independent
    assert ki[0, 0, 0, 0] != ki[0, 0, 0, 1]                 # gas0 X-dependent

    # additive overlap: tau = k0*amount0 + k1*amount1 on the shared g-point
    T2 = np.array([[250.0, 250.0]])           # (nlev=1, ncol=2)
    p2 = np.array([[1e3, 1e4]])
    h2o2 = np.array([[1e-6, 1.0]])
    amounts = np.stack([                       # (ngas, nlev, ncol)
        np.full((1, 2), 10.0),                 # air column mass
        np.full((1, 2), 2.0),                  # o3 column mass
    ])
    tau = compute_ck_optical_depth(table, T2, p2, amounts, h2o_vmr=h2o2)
    # column 0: k0=1.0*10 + 7*2 = 24 ; column 1: k0=4.0*10 + 7*2 = 54
    np.testing.assert_allclose(tau[0, 0, 0, 0], 24.0)
    np.testing.assert_allclose(tau[0, 0, 0, 1], 54.0)
```

- [ ] **Step 2: Run the test**

Run: `conda run -n climt python -m pytest tests/test_cork_optics.py::test_interpolate_k_6d_two_gas_additive -v`
Expected: PASS. If it FAILS, debug `interpolate_k` / `_compute_ck_optical_depth_additive` until the documented arithmetic holds, then re-run.

- [ ] **Step 3: Commit**

```bash
git add tests/test_cork_optics.py
git commit -m "test(cork): pin 6-D H2O-axis + ngas=2 additive overlap path"
```

---

## Task 6: Component — declare `mole_fraction_of_ozone_in_air` in `input_properties`

When the table is premixed-background AND carries secondary gases (e.g. `o3`), the component must request their runtime profiles. The premixed-bg branch (`climt/_components/cork/lw/component.py:101-107`) currently only adds `specific_humidity`.

**Files:**
- Modify: `climt/_components/cork/lw/component.py:100-107`
- Test: `tests/test_cork_lw.py`

- [ ] **Step 1: Write the failing test (with a synthetic 2-gas table fixture)**

Add to `tests/test_cork_lw.py` (imports `numpy as np`, `os`, and the writer from `scripts`):

```python
def _write_o3_test_table(path):
    """Tiny 6-D, 2-gas (effective+o3) premixed-bg LW table for component tests."""
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    from cork_table_builder.netcdf_writer import write_lw_table

    nband, ngpt, nT, nP, nX = 4, 2, 4, 4, 4
    band_edges = np.array([10.0, 500.0, 800.0, 1800.0, 3250.0])  # O3 -> band idx 2
    k = np.zeros((2, nband, ngpt, nT, nP, nX))
    k[0] = 1e-3                       # air: weak grey background per kg-air
    k[1, 2] = 5.0                     # o3: absorbs only in 800-1800 band, per kg-O3
    write_lw_table(
        path,
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(180.0, 320.0, nT),
        log_p_grid=np.linspace(np.log(1.0), np.log(1e5), nP),
        band_edges=band_edges,
        planck_fraction=np.full((nband, ngpt, nT), 0.25),
        h2o_vmr_grid=np.array([1e-6, 1e-4, 1e-2, 1.0]),
        gas_names=("effective", "o3"),
        background_is_premixed=True,
    )


def test_cork_lw_o3_table_requests_ozone_input(tmp_path):
    from climt import CorkLongwaveRadiation
    tbl = str(tmp_path / "o3.nc")
    _write_o3_test_table(tbl)
    comp = CorkLongwaveRadiation(optics="correlated_k", table=tbl)
    assert "mole_fraction_of_ozone_in_air" in comp.input_properties
    assert "specific_humidity" in comp.input_properties
    assert comp.input_properties["mole_fraction_of_ozone_in_air"]["alias"] == "o3"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py::test_cork_lw_o3_table_requests_ozone_input -v`
Expected: FAIL — `mole_fraction_of_ozone_in_air` not in `input_properties`.

- [ ] **Step 3: Implement — add secondary-gas inputs in the premixed-bg branch**

In `climt/_components/cork/lw/component.py`, replace the premixed-bg block (the `if self._premixed_bg:` body that adds `specific_humidity`) with:

```python
            if self._premixed_bg:
                # Background is pre-mixed; H2O VMR is the live X-axis.
                props["specific_humidity"] = {
                    "dims": ["mid_levels", "*"],
                    "units": "kg/kg",
                    "alias": "h2o",
                }
                # Secondary additive gases (gas index >= 1) carry their own
                # runtime mole-fraction profile (e.g. O3 from init_ozone).
                _SECONDARY_CF_NAME = {"o3": "mole_fraction_of_ozone_in_air"}
                for gas in self._gas_names[1:]:
                    cf_name = _SECONDARY_CF_NAME.get(
                        gas, f"mole_fraction_of_{gas}_in_air"
                    )
                    props[cf_name] = {
                        "dims": ["mid_levels", "*"],
                        "units": "mole/mole",
                        "alias": gas,
                    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py::test_cork_lw_o3_table_requests_ozone_input -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_components/cork/lw/component.py tests/test_cork_lw.py
git commit -m "feat(cork-lw): request mole_fraction_of_ozone_in_air for O3 tables"
```

---

## Task 7: Component — fill `gas_amounts` for secondary gases in `array_call`

The premixed-bg branch in `array_call` (`climt/_components/cork/lw/component.py:236-249`) sets only `gas_amounts[0]` (air mass) and `h2o_vmr`. Extend it to convert each secondary gas's mole fraction to column mass.

**Files:**
- Modify: `climt/_components/cork/lw/component.py:236-249`
- Test: `tests/test_cork_lw.py`

- [ ] **Step 1: Write the failing test (O₃ on/off behaviour)**

Add to `tests/test_cork_lw.py` (reuses `_write_o3_test_table`):

```python
def _o3_state(ncol=1, nlev=10):
    """Minimal state dict for CorkLongwaveRadiation correlated_k call."""
    import numpy as np
    p_int = np.linspace(1e5, 1.0, nlev + 1)[:, None] * np.ones((1, ncol))
    p = 0.5 * (p_int[1:] + p_int[:-1])
    T = np.full((nlev, ncol), 250.0)
    return {
        "T": T, "p": p, "p_int": p_int,
        "T_surf": np.full((ncol,), 288.0),
        "emissivity": np.ones((4, ncol)),
        "h2o": np.full((nlev, ncol), 1e-6),
        "tau_cloud_lw": np.zeros((nlev, ncol, 4)),
    }


def test_cork_lw_o3_increases_window_band_optical_depth(tmp_path):
    import numpy as np
    from climt import CorkLongwaveRadiation
    tbl = str(tmp_path / "o3.nc")
    _write_o3_test_table(tbl)
    comp = CorkLongwaveRadiation(optics="correlated_k", table=tbl)

    state = _o3_state()
    nlev, ncol = state["T"].shape
    # O3 zero everywhere
    state_off = dict(state, o3=np.zeros((nlev, ncol)))
    _, diag_off = comp.array_call(state_off)
    # O3 present (1 ppm)
    state_on = dict(state, o3=np.full((nlev, ncol), 1e-6))
    _, diag_on = comp.array_call(state_on)

    tau_off = diag_off["longwave_optical_depth_per_band"]
    tau_on = diag_on["longwave_optical_depth_per_band"]
    # band index 2 is the 800-1800 window where the synthetic O3 absorbs
    assert np.all(tau_on[..., 2] > tau_off[..., 2])
    # other bands unchanged (no O3 absorption there)
    np.testing.assert_allclose(tau_on[..., 0], tau_off[..., 0])
    np.testing.assert_allclose(tau_on[..., 3], tau_off[..., 3])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py::test_cork_lw_o3_increases_window_band_optical_depth -v`
Expected: FAIL — either a `KeyError`/missing `o3` handling or `tau_on == tau_off` because `gas_amounts[1]` is never populated.

- [ ] **Step 3: Implement — fill secondary-gas amounts**

In `climt/_components/cork/lw/component.py`, in the `elif self._premixed_bg:` block of `array_call`, after the existing line that computes `h2o_vmr = q_h2o / np.maximum(...)`, append:

```python
                # Secondary additive gases: mole fraction -> mass mixing ratio
                # -> column mass. k for these gases is per-kg-species.
                for ig, gas in enumerate(self._gas_names):
                    if ig == 0:
                        continue
                    x_gas = state[gas].reshape(nlev, -1)  # mole fraction
                    M_gas = MOLAR_MASS.get(gas, MOLAR_MASS_DRY_AIR)
                    q_gas = x_gas * (M_gas / MOLAR_MASS_DRY_AIR)
                    gas_amounts[ig, :, :] = compute_column_amount(
                        q_gas, p_int_flat, g
                    )
```

(`MOLAR_MASS`, `MOLAR_MASS_DRY_AIR`, and `compute_column_amount` are already imported at the top of the file.)

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py::test_cork_lw_o3_increases_window_band_optical_depth -v`
Expected: PASS.

- [ ] **Step 5: Run the full CORK LW suite for regressions**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py tests/test_cork_optics.py -q`
Expected: all PASS (existing single-gas tables still use the unchanged single-gas path).

- [ ] **Step 6: Commit**

```bash
git add climt/_components/cork/lw/component.py tests/test_cork_lw.py
git commit -m "feat(cork-lw): fill column amounts for secondary additive gases (O3)"
```

---

## Task 8: Generate the shipped Earth+O₃ table (radiation/linepyline env)

This is the long-running line-by-line generation. It runs in the `radiation`/`linepyline` env, NOT `climt`. Uses the validated linepyline recipe (`binning=False`, `dnu=0.1`, `pseudovoigt`, `ngpt=8`, Gauss-Legendre quadrature — the linepyline path's `kappa_to_k_coeffs`).

**Files:**
- Create: `climt/_data/cork/correlated_k/earth_low_res_lw_co2refined_o3.nc`

- [ ] **Step 1: Generate the table**

Run from repo root:

```bash
conda run -n radiation python scripts/generate_cork_tables_linepyline.py \
  --scenario earth --kind lw \
  --output climt/_data/cork/correlated_k/earth_low_res_lw_co2refined_o3.nc \
  --ngpt 8 --dnu 0.1 --line-shape pseudovoigt
```

Expected final line: `[earth/lw] wrote .../earth_low_res_lw_co2refined_o3.nc  k_coefficients.shape=(2, 6, 8, 12, 8, 7)` — note `shape[0] == 2` (two gases).

- [ ] **Step 2: Verify dimensions, gas names, and X-constant O₃**

Run (in `climt` env):

```bash
conda run -n climt python -c "
import numpy as np
from scipy.io import netcdf_file
nc = netcdf_file('climt/_data/cork/correlated_k/earth_low_res_lw_co2refined_o3.nc','r',mmap=False)
k = nc.variables['k_coefficients'][:]
print('shape', k.shape, 'gas_names', nc.gas_names, 'bg_premixed', nc.background_is_premixed)
# O3 (gas 1) must be constant along the X (last) axis
print('o3 X-constant:', np.allclose(k[1], k[1][...,:1]))
# O3 absorbs only in the 9.6um window band (edges 10,500,630,700,800,1800,3250 -> band idx 4 = 800-1800)
band_max = k[1].reshape(k.shape[1], -1).max(axis=1)
print('o3 per-band max:', band_max)
"
```

Expected: `shape (2, 6, 8, 12, 8, 7)`, `gas_names b'effective,o3'`, `bg_premixed b'true'`, `o3 X-constant: True`, and the O₃ per-band max concentrated in the 800–1800 cm⁻¹ band (index 4), near-zero elsewhere.

- [ ] **Step 3: Smoke-test loading through the component**

Run (in `climt` env):

```bash
conda run -n climt python -c "
from climt import CorkLongwaveRadiation
c = CorkLongwaveRadiation(optics='correlated_k',
    table='earth_low_res_lw_co2refined_o3')
print('bands', c.num_longwave_bands, 'gases', c._gas_names, 'premixed_bg', c._premixed_bg)
assert 'mole_fraction_of_ozone_in_air' in c.input_properties
print('OK')
"
```

Expected: `bands 6 gases ['effective', 'o3'] premixed_bg True` then `OK`.

- [ ] **Step 4: Commit the table**

```bash
git add climt/_data/cork/correlated_k/earth_low_res_lw_co2refined_o3.nc
git commit -m "data(cork): add Earth LW table with O3 9.6um additive gas (6-band, GL-8)"
```

---

## Task 9: RRTMG validation — recover the ~30 K stratospheric warming

End-to-end check that CORK-with-O₃ recovers the stratospheric warming O₃ contributes in RRTMG, using climt's `init_ozone` profile on BOTH sides (apples-to-apples). Output figures to `debug_data/` per memory `feedback_experiment_scripts_location.md`.

**Files:**
- Create: `scripts/experiments/radiative_equilibrium_cork_o3.py`

- [ ] **Step 1: Write the validation script**

Create `scripts/experiments/radiative_equilibrium_cork_o3.py`. Base it on the existing `scripts/experiments/radiative_equilibrium_cork_convergence.py` (read it first for the column setup, surface/SW configuration, and integration loop). Adapt it to a 3-way comparison at radiative equilibrium (dry, CO₂=376 ppm, SW=240 W/m², ~1000 days):

1. **RRTMG** with `init_ozone` ozone profile (default `get_default_state`).
2. **CORK (6-band CO₂-only)** using `earth_low_res_lw_co2refined_linepyline` — no O₃.
3. **CORK (6-band CO₂+O₃)** using `earth_low_res_lw_co2refined_o3` — feed the SAME `mole_fraction_of_ozone_in_air` profile RRTMG uses.

Plot equilibrium T(p) for all three and save to `debug_data/cork_o3_rce_comparison.png`. Print T at the exp #10 diagnostic nodes (90, 40, 10, 2.4 hPa).

- [ ] **Step 2: Run the validation**

Run: `conda run -n climt python scripts/experiments/radiative_equilibrium_cork_o3.py`
Expected: writes `debug_data/cork_o3_rce_comparison.png` and prints a node table.

- [ ] **Step 3: Check acceptance**

Acceptance criteria (record the actual numbers in the experiment log):
- CORK-with-O₃ is **warmer than CORK-without-O₃** through the stratosphere (the O₃ warming has the right sign and is order ~10–30 K at 40–90 hPa).
- Under matched O₃, the CORK-vs-RRTMG stratospheric gap at 40–90 hPa is materially smaller than the ~26 K no-O₃ artifact gap reported in exp #10 (target: within roughly ±10 K, acknowledging CORK's coarse bands and the residual upper-strat integrator bias exp #10 flagged separately).

If the sign is wrong or O₃ has no effect, debug with `superpowers:systematic-debugging` (likely suspects: O₃ amount unit conversion in Task 7, or the O₃ band placement in Task 8).

- [ ] **Step 4: Log the experiment**

Append an experiment entry (tries / findings / numbers) to `docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md` per memory `feedback_log_experiments.md`, cross-referencing this plan and the design spec.

- [ ] **Step 5: Commit**

```bash
git add scripts/experiments/radiative_equilibrium_cork_o3.py docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md
git commit -m "test(cork): RCE validation of O3 stratospheric warming vs RRTMG"
```

---

## Task 10: Refresh the knowledge graph

Per `CLAUDE.md`, after modifying code update graphify.

- [ ] **Step 1: Update and augment**

```bash
graphify update .
conda run -n climt python scripts/augment_graph.py
```

Expected: both complete without error; new `sample_secondary_gas_kappa` / `broadcast_and_stack_gases` nodes and the O₃ component edges appear.

- [ ] **Step 2: Commit graph changes (if any tracked files changed)**

```bash
git add graphify-out/
git commit -m "chore(graph): refresh after CORK O3 absorber"
```

---

## Self-Review Notes

- **Spec coverage:** builder (Tasks 1–4, 8) ↔ spec "What changes" #1/#3; component (Tasks 6–7) ↔ #2; untested-combination risk (Task 5) ↔ spec risk table; validation (Task 9) ↔ spec "Validation"; profile reuse via `init_ozone` ↔ spec de-risk. SW O₃, η-binning, table retirement remain explicitly out of scope.
- **Type/name consistency:** `broadcast_and_stack_gases(primary_k, secondary_k_list)`, `sample_secondary_gas_kappa(species=..., ...)`, `lw_secondary_gases`, `background_is_premixed`, gas alias `o3`, CF name `mole_fraction_of_ozone_in_air` — used consistently across tasks.
- **Quadrature note:** the linepyline `kappa_to_k_coeffs` is Gauss-Legendre; the spec's "two-stretch g_split=0.97" referred to the separate exo_k `generate_cork_tables.py` path and does not apply to this linepyline-built table.
```

