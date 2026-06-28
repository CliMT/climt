# CO₂-adjustable, RRTMG-fidelity Earth-LW cork table — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a production Earth-LW correlated-k scheme that decouples the H₂O continuum, uses a refined 14-band partition, and is CO₂-adjustable at runtime (10–10000 ppm) from a single table, becoming the default `earth_low_res_lw`.

**Architecture:** Extend the existing premixed-background H₂O-axis machinery with a second runtime absorber axis for CO₂. The line k-distribution gains a `(…, nX_H₂O, nX_CO₂)` shape; the H₂O continuum stays decoupled and CO₂-independent. CO₂ is read from state, interpolated geometrically (log-k vs log-X_CO₂) and quadrilinearly with T, log p, log X_H₂O. The LW transport kernel is consolidated into a single `@njit(parallel=True)` routine to absorb the ~3.5× band-count growth, reproducing current fluxes bit-for-bit.

**Tech Stack:** Python, NumPy, numba (njit/prange), scipy.io.netcdf, xarray (generation only), linepyline + MT_CKD (offline table generation only, in the `linepyline` conda env). All pytest runs use the `climt` conda env.

---

## Background the engineer needs

The Earth LW table is a **premixed-background** correlated-k table: `gas_names == ["effective"]`, k is in m²/kg-of-air, and H₂O is the only *live* (runtime) absorber axis. The component branch that handles this is `_premixed_bg` in `climt/_components/cork/lw/component.py`. This plan adds CO₂ as a **second** live axis to that same path. CO₂ today is baked into the table at 376 ppm and not read from state; after this work it is read from `mole_fraction_of_carbon_dioxide_in_air` (already in the default state at 330 ppm, see `climt/_core/initialization.py:730`).

**Data layout.** A premixed-bg table's `k_coefficients` is currently 6-D: `(ngas=1, nband, ngpt, nT, nP, nX_H₂O)`. After this work, when a CO₂ axis is present it is 7-D: `(ngas=1, nband, ngpt, nT, nP, nX_H₂O, nX_CO₂)`. The decoupled `continuum_kappa` stays 4-D `(nband, nT, nP, nX_H₂O)` — H₂O-only, no CO₂ axis. All new axes are **optional and detected by presence**, so existing tables behave exactly as today.

**Interpolation scheme (resolving a spec ambiguity).** The design spec is internally inconsistent: the Design section (line 46) says the CO₂ axis defaults to **log-k vs log-X_CO₂ (geometric)** and calls this "the starting design," while the code-changes summary (line 58) says "linear-in-k vs log-X_CO₂." We follow the **Design section: geometric (log-k) in CO₂**, because that section is the detailed/authoritative one and explicitly justifies log-k (CO₂ band-mean k is convex/saturating in amount — linear-in-value over log-spaced nodes over-estimates it, the exact failure mode that caused the continuum bug, exp #25). The H₂O line axis stays **linear-in-k** (switching it is a deferred non-goal). We implement CO₂ interpolation so the k-vs-X mode is a single named constant, because the spec says the final log-k/linear-k choice is "confirmed by the A5 CO₂-accuracy check."

**What this plan does NOT do in code:** the actual multi-hour line-by-line table *generation* run (needs linepyline + an adequate machine). Tasks 1–8 are fully test-driven against tiny synthetic tables and never call linepyline. Task 9 documents and runs the real generation. Tasks 10–11 validate against LBL/RRTMG and flip the default.

---

## File structure

| File | Responsibility | Change |
|------|----------------|--------|
| `scripts/cork_table_builder/netcdf_writer.py` | netCDF schema writer | add `co2_vmr_grid` var + CO₂ dim on `k_coefficients` (Task 1) |
| `climt/_components/cork/optics/correlated_k.py` | loader, `interpolate_k`, optical-depth assembly | read `co2_vmr_grid`; quadrilinear interp; thread `co2_vmr` (Tasks 2–4) |
| `climt/_components/cork/lw/component.py` | LW component | declare CO₂ input, pass VMR into optics (Task 5) |
| `scripts/cork_table_builder/kappa_sampling.py` | line-by-line κ sampler | add CO₂ sweep axis (Task 6) |
| `scripts/generate_cork_tables_linepyline.py` | table generation driver | `earth` 14-band + CO₂-grid config; 4-D moveaxis (Task 7) |
| `climt/_components/cork/lw/kernels.py` | LW transport | consolidate into single `@njit(parallel=True)` (Task 8) |
| `scripts/experiments/eval_band_structure.py` | LBL-vs-CORK evaluator | CO₂×H₂O grid (Task 10) |
| `tests/`, `tests/cork_table_builder/` | unit tests | added throughout |

---

## Task 1: Writer — add CO₂ axis to `write_lw_table`

**Files:**
- Modify: `scripts/cork_table_builder/netcdf_writer.py:14-104`
- Test: `tests/cork_table_builder/test_netcdf_writer.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/cork_table_builder/test_netcdf_writer.py`:

```python
def test_write_lw_table_with_co2_axis_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.netcdf_writer import write_lw_table
    from climt._components.cork.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP, nX, nC = 1, 4, 2, 3, 5, 4, 6
    k = np.random.RandomState(3).uniform(
        1e-7, 1e-2, size=(ngas, nband, ngpt, nT, nP, nX, nC))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    h2o_vmr = np.array([1e-6, 1e-4, 1e-2, 1e-1])
    co2_vmr = np.logspace(-5, -2, nC)
    planck_fraction = np.full((nband, ngpt, nT), 1.0 / ngpt)
    # continuum stays H2O-only (no CO2 axis): (nband, nT, nP, nX)
    cont = np.random.RandomState(4).uniform(1e-6, 1e-3, size=(nband, nT, nP, nX))

    path = tmp_path / "test_lw_co2.nc"
    write_lw_table(
        str(path),
        k_coefficients=k,
        gpoint_weights=weights,
        T_grid=T_grid,
        log_p_grid=log_p_grid,
        band_edges=band_edges,
        planck_fraction=planck_fraction,
        h2o_vmr_grid=h2o_vmr,
        co2_vmr_grid=co2_vmr,
        continuum_kappa=cont,
        gas_names=("effective",),
        source="linepyline:test_co2",
    )
    loaded = load_k_table(str(path))
    assert loaded["k_coefficients"].shape == (ngas, nband, ngpt, nT, nP, nX, nC)
    np.testing.assert_allclose(loaded["co2_vmr_grid"], co2_vmr, rtol=1e-5)
    np.testing.assert_allclose(loaded["h2o_vmr_grid"], h2o_vmr, rtol=1e-5)
    assert loaded["continuum_kappa"].shape == (nband, nT, nP, nX)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_netcdf_writer.py::test_write_lw_table_with_co2_axis_roundtrip -v`
Expected: FAIL — `write_lw_table() got an unexpected keyword argument 'co2_vmr_grid'`.

- [ ] **Step 3: Implement the CO₂ axis in `write_lw_table`**

In `scripts/cork_table_builder/netcdf_writer.py`, change the signature and body of `write_lw_table`. Replace the signature line `h2o_vmr_grid=None,` to add `co2_vmr_grid=None,` right after it, and replace the shape-unpacking + dimension/variable block. The new function reads:

```python
def write_lw_table(
    out_path,
    *,
    k_coefficients,
    gpoint_weights,
    T_grid,
    log_p_grid,
    band_edges,
    planck_fraction,
    h2o_vmr_grid=None,
    co2_vmr_grid=None,
    gas_names=("effective",),
    overlap_method="additive",
    resolution="low",
    source="linepyline",
    continuum_kappa=None,
):
    """Write LW k-table.

    k_coefficients shape:
        (ngas, nband, ngpt, nT, nP)                  no live-X axis
        (ngas, nband, ngpt, nT, nP, nX_H2O)          H2O axis
        (ngas, nband, ngpt, nT, nP, nX_H2O, nX_CO2)  H2O + CO2 axes

    co2_vmr_grid (optional): (nX_CO2,) log-spaced CO2 volume mixing ratios.
        Requires h2o_vmr_grid. Adds a trailing CO2 dimension to k_coefficients.

    continuum_kappa (optional): band-grey H2O continuum mass-absorption
        (m^2/kg-air), shape (nband, nT, nP, nX_H2O). H2O-only and CO2-independent
        (no CO2 axis) — decouples its quadratic-in-H2O scaling from the line CDF
        (cf. RRTMG's separate tauself/taufor).
    """
    _ensure_parent(out_path)
    has_x = (h2o_vmr_grid is not None)
    has_co2 = (co2_vmr_grid is not None)
    if has_co2 and not has_x:
        raise ValueError("co2_vmr_grid requires an h2o_vmr axis")
    if has_co2:
        ngas, nband, ngpt, nT, nP, nX, nC = k_coefficients.shape
    elif has_x:
        ngas, nband, ngpt, nT, nP, nX = k_coefficients.shape
    else:
        ngas, nband, ngpt, nT, nP = k_coefficients.shape
    edges = np.asarray(band_edges)
    limits = np.column_stack([edges[:-1], edges[1:]])

    with netcdf_file(out_path, "w") as nc:
        nc.createDimension("gas", ngas)
        nc.createDimension("band", nband)
        nc.createDimension("gpoint", ngpt)
        nc.createDimension("temperature", nT)
        nc.createDimension("pressure", nP)
        nc.createDimension("bounds", 2)
        if has_x:
            nc.createDimension("h2o_vmr", nX)
        if has_co2:
            nc.createDimension("co2_vmr", nC)
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure", "h2o_vmr", "co2_vmr"))
        elif has_x:
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure", "h2o_vmr"))
        else:
            v = nc.createVariable("k_coefficients", "f4",
                                  ("gas", "band", "gpoint", "temperature",
                                   "pressure"))
        v[:] = k_coefficients.astype("f4")

        w = nc.createVariable("gpoint_weights", "f4", ("band", "gpoint"))
        w[:] = gpoint_weights.astype("f4")

        tg = nc.createVariable("temperature_grid", "f4", ("temperature",))
        tg[:] = np.asarray(T_grid, dtype="f4")

        pg = nc.createVariable("pressure_grid_log", "f4", ("pressure",))
        pg[:] = np.asarray(log_p_grid, dtype="f4")

        bl = nc.createVariable("band_wavenumber_limits", "f4", ("band", "bounds"))
        bl[:] = limits.astype("f4")

        cork = nc.createVariable("planck_fraction", "f4",
                               ("band", "gpoint", "temperature"))
        cork[:] = planck_fraction.astype("f4")

        if has_x:
            xv = nc.createVariable("h2o_vmr_grid", "f4", ("h2o_vmr",))
            xv[:] = np.asarray(h2o_vmr_grid, dtype="f4")

        if has_co2:
            cvg = nc.createVariable("co2_vmr_grid", "f4", ("co2_vmr",))
            cvg[:] = np.asarray(co2_vmr_grid, dtype="f4")

        if continuum_kappa is not None:
            if not has_x:
                raise ValueError("continuum_kappa requires an h2o_vmr axis")
            ck = np.asarray(continuum_kappa)
            if ck.shape != (nband, nT, nP, nX):
                raise ValueError(
                    f"continuum_kappa shape {ck.shape} != "
                    f"(nband,nT,nP,nX)={(nband, nT, nP, nX)}")
            cv = nc.createVariable(
                "continuum_kappa", "f4",
                ("band", "temperature", "pressure", "h2o_vmr"))
            cv[:] = ck.astype("f4")

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = overlap_method
        nc.resolution = resolution
        nc.source = source
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_netcdf_writer.py -v`
Expected: PASS (all existing + new test). The new test also exercises the loader (`co2_vmr_grid` must already be read — see Task 2; if you run Task 1 alone the `loaded["co2_vmr_grid"]` assertion will KeyError, so do Task 2 in the same commit).

- [ ] **Step 5: Commit** (after Task 2)

---

## Task 2: Loader — read `co2_vmr_grid`

**Files:**
- Modify: `climt/_components/cork/optics/correlated_k.py:10-21`
- Test: covered by Task 1's roundtrip test.

- [ ] **Step 1: Add `co2_vmr_grid` to the loaded variables**

In `climt/_components/cork/optics/correlated_k.py`, edit the `_NETCDF_VARS` tuple to add `"co2_vmr_grid"`:

```python
_NETCDF_VARS = (
    "k_coefficients",
    "gpoint_weights",
    "temperature_grid",
    "pressure_grid_log",
    "h2o_vmr_grid",
    "co2_vmr_grid",
    "band_wavenumber_limits",
    "planck_fraction",
    "solar_source_per_gpoint",
    "rayleigh_coefficient",
    "continuum_kappa",
)
```

- [ ] **Step 2: Run the Task 1 test to verify it passes**

Run: `conda run -n climt python -m pytest tests/cork_table_builder/test_netcdf_writer.py::test_write_lw_table_with_co2_axis_roundtrip -v`
Expected: PASS.

- [ ] **Step 3: Commit Tasks 1+2 together**

```bash
git add scripts/cork_table_builder/netcdf_writer.py \
        climt/_components/cork/optics/correlated_k.py \
        tests/cork_table_builder/test_netcdf_writer.py
git commit -m "feat(cork): add optional CO2 axis to LW k-table format

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Interpolation — quadrilinear CO₂ axis in `interpolate_k`

The 7-D case reuses the existing trilinear `(T, log p, log X_H₂O)` interpolation **twice** (once on each adjacent CO₂ slice, linear-in-k as today) and then blends the two results **geometrically** in `(log-k vs log X_CO₂)`. A module-level constant `_CO2_INTERP_LOGK = True` selects geometric vs linear-k so the A5 accuracy check can flip it.

**Files:**
- Modify: `climt/_components/cork/optics/correlated_k.py:97-195`
- Test: `tests/test_cork_co2_interp.py` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/test_cork_co2_interp.py`:

```python
import numpy as np


def _build_co2_table(tmp_path, k):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from scripts.cork_table_builder.netcdf_writer import write_lw_table
    from climt._components.cork.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP, nX, nC = k.shape
    write_lw_table(
        str(tmp_path / "t.nc"),
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(200.0, 300.0, nT),
        log_p_grid=np.linspace(np.log(1e2), np.log(1e5), nP),
        band_edges=np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0]),
        planck_fraction=np.full((nband, ngpt, nT), 1.0 / ngpt),
        h2o_vmr_grid=np.array([1e-6, 1e-3, 1e-1]),
        co2_vmr_grid=np.array([1e-4, 1e-3, 1e-2]),
        gas_names=("effective",),
    )
    return load_k_table(str(tmp_path / "t.nc"))


def test_interpolate_k_recovers_node_value(tmp_path):
    """At an exact grid node, quadrilinear interp returns the stored k."""
    from climt._components.cork.optics.correlated_k import interpolate_k

    rng = np.random.RandomState(0)
    k = rng.uniform(1e-6, 1e-2, size=(1, 4, 2, 3, 4, 3, 3))
    table = _build_co2_table(tmp_path, k)

    # Pick node iT=1, iP=2, iX=1 (1e-3), iC=2 (1e-2).
    T = np.array([250.0])
    p = np.array([np.exp(table["pressure_grid_log"][2])])
    h2o = np.array([1e-3])
    co2 = np.array([1e-2])
    out = interpolate_k(table, T, p, h2o_vmr=h2o, co2_vmr=co2)
    # out shape (ngas, nband, ngpt, ncol)
    np.testing.assert_allclose(
        out[:, :, :, 0], k[:, :, :, 1, 2, 1, 2].astype(np.float32), rtol=2e-6
    )


def test_interpolate_k_co2_geometric_midpoint(tmp_path):
    """Halfway (in log-X_CO2) between two CO2 nodes, log-k interp returns the
    geometric mean of the two node values (T,P,X_H2O all on-node)."""
    from climt._components.cork.optics.correlated_k import interpolate_k

    k = np.ones((1, 1, 1, 3, 4, 3, 3))
    # CO2 nodes 1e-4,1e-3,1e-2; set k=2 at iC=1 and k=8 at iC=2 for one node.
    k[0, 0, 0, 1, 2, 1, 1] = 2.0
    k[0, 0, 0, 1, 2, 1, 2] = 8.0
    table = _build_co2_table(tmp_path, k)

    T = np.array([250.0])
    p = np.array([np.exp(table["pressure_grid_log"][2])])
    h2o = np.array([1e-3])
    co2 = np.array([np.sqrt(1e-3 * 1e-2)])  # geometric midpoint of the two nodes
    out = interpolate_k(table, T, p, h2o_vmr=h2o, co2_vmr=co2)
    # geometric mean of 2 and 8 = 4
    np.testing.assert_allclose(out[0, 0, 0, 0], 4.0, rtol=1e-4)
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_cork_co2_interp.py -v`
Expected: FAIL — `interpolate_k() got an unexpected keyword argument 'co2_vmr'`.

- [ ] **Step 3: Implement the CO₂ axis in `interpolate_k`**

In `climt/_components/cork/optics/correlated_k.py`, add a module-level constant just below `_NETCDF_VARS`:

```python
# CO2 runtime axis interpolation: True = geometric (log-k vs log-X_CO2),
# the design default (CO2 band-mean k is convex/saturating in amount, so
# linear-in-value over log-spaced nodes over-estimates it). Flip to False to
# evaluate linear-in-k during the A5 CO2-accuracy check.
_CO2_INTERP_LOGK = True
```

Replace the whole body of `interpolate_k` (currently lines 97-195) with the version below. It adds the `co2_vmr` parameter, handles `k.ndim == 7`, factors the existing `(T, P, X_H₂O)` trilinear blend into an inline helper applied per CO₂ slice, and geometrically blends in CO₂:

```python
def interpolate_k(table, T, p, h2o_vmr=None, co2_vmr=None):
    """Interpolate k-coefficients in (T, log p[, log X_H2O[, log X_CO2]]) space.

    Args:
        table: loaded k-table. ``k_coefficients`` may be shaped
            ``(ngas, nband, ngpt, nT, nP)`` (bilinear),
            ``(..., nP, nX)`` (trilinear, H2O axis), or
            ``(..., nP, nX, nC)`` (quadrilinear, H2O + CO2 axes).
        T: (ncol,) temperature values, K
        p: (ncol,) pressure values, Pa
        h2o_vmr: (ncol,) H2O mole fraction per column. Required if the table
            has an H2O axis.
        co2_vmr: (ncol,) CO2 mole fraction per column. Required if the table
            has a CO2 axis.

    Returns:
        k_interp: (ngas, nband, ngpt, ncol)
    """
    k = table["k_coefficients"]
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]
    has_co2_axis = k.ndim == 7
    has_x_axis = k.ndim >= 6

    ncol = len(T)
    log_p = np.log(np.maximum(p, 1.0))

    nX = nC = None
    if has_co2_axis:
        ngas, nband, ngpt, nT, nP, nX, nC = k.shape
    elif has_x_axis:
        ngas, nband, ngpt, nT, nP, nX = k.shape
    else:
        ngas, nband, ngpt, nT, nP = k.shape

    if has_x_axis:
        x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
        log_x_grid = np.log(np.maximum(x_grid, 1e-30))
        if h2o_vmr is None:
            raise ValueError(
                "k-table has an h2o_vmr_grid axis but h2o_vmr was not provided"
            )
        x_clamped = np.clip(h2o_vmr, float(x_grid[0]), float(x_grid[-1]))
        log_x = np.log(np.maximum(x_clamped, 1e-30))

    if has_co2_axis:
        c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
        log_c_grid = np.log(np.maximum(c_grid, 1e-30))
        if co2_vmr is None:
            raise ValueError(
                "k-table has a co2_vmr_grid axis but co2_vmr was not provided"
            )
        c_clamped = np.clip(co2_vmr, float(c_grid[0]), float(c_grid[-1]))
        log_c = np.log(np.maximum(c_clamped, 1e-30))

    result = np.zeros((ngas, nband, ngpt, ncol))

    for col in range(ncol):
        iT = max(0, min(np.searchsorted(T_grid, T[col]) - 1, nT - 2))
        fT = (T[col] - T_grid[iT]) / (T_grid[iT + 1] - T_grid[iT])
        fT = max(0.0, min(1.0, fT))

        iP = max(0, min(np.searchsorted(p_grid_log, log_p[col]) - 1, nP - 2))
        fP = (log_p[col] - p_grid_log[iP]) / (p_grid_log[iP + 1] - p_grid_log[iP])
        fP = max(0.0, min(1.0, fP))

        if has_x_axis:
            iX = max(0, min(np.searchsorted(log_x_grid, log_x[col]) - 1, nX - 2))
            fX = (log_x[col] - log_x_grid[iX]) / (
                log_x_grid[iX + 1] - log_x_grid[iX])
            fX = max(0.0, min(1.0, fX))

        if has_co2_axis:
            iC = max(0, min(np.searchsorted(log_c_grid, log_c[col]) - 1, nC - 2))
            fC = (log_c[col] - log_c_grid[iC]) / (
                log_c_grid[iC + 1] - log_c_grid[iC])
            fC = max(0.0, min(1.0, fC))

        def _txx(base):
            # Trilinear (T, P, X_H2O), linear-in-k. base shape (nT, nP, nX).
            v000 = base[iT, iP, iX]
            v100 = base[iT + 1, iP, iX]
            v010 = base[iT, iP + 1, iX]
            v110 = base[iT + 1, iP + 1, iX]
            v001 = base[iT, iP, iX + 1]
            v101 = base[iT + 1, iP, iX + 1]
            v011 = base[iT, iP + 1, iX + 1]
            v111 = base[iT + 1, iP + 1, iX + 1]
            x0 = (v000 * (1 - fT) * (1 - fP) + v100 * fT * (1 - fP)
                  + v010 * (1 - fT) * fP + v110 * fT * fP)
            x1 = (v001 * (1 - fT) * (1 - fP) + v101 * fT * (1 - fP)
                  + v011 * (1 - fT) * fP + v111 * fT * fP)
            return x0 * (1 - fX) + x1 * fX

        def _tp(base):
            # Bilinear (T, P), linear-in-k. base shape (nT, nP).
            return (base[iT, iP] * (1 - fT) * (1 - fP)
                    + base[iT + 1, iP] * fT * (1 - fP)
                    + base[iT, iP + 1] * (1 - fT) * fP
                    + base[iT + 1, iP + 1] * fT * fP)

        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    if has_co2_axis:
                        c0 = _txx(k[ig, ib, igp, :, :, :, iC])
                        c1 = _txx(k[ig, ib, igp, :, :, :, iC + 1])
                        if _CO2_INTERP_LOGK:
                            lc0 = np.log(max(c0, 1e-40))
                            lc1 = np.log(max(c1, 1e-40))
                            result[ig, ib, igp, col] = np.exp(
                                lc0 * (1 - fC) + lc1 * fC)
                        else:
                            result[ig, ib, igp, col] = c0 * (1 - fC) + c1 * fC
                    elif has_x_axis:
                        result[ig, ib, igp, col] = _txx(k[ig, ib, igp])
                    else:
                        result[ig, ib, igp, col] = _tp(k[ig, ib, igp])

    return result
```

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n climt python -m pytest tests/test_cork_co2_interp.py -v`
Expected: PASS (both tests).

- [ ] **Step 5: Run the existing optics tests to confirm no regression**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py -v`
Expected: PASS (the 6-D/5-D paths are unchanged in behavior).

- [ ] **Step 6: Commit**

```bash
git add climt/_components/cork/optics/correlated_k.py \
        tests/test_cork_co2_interp.py
git commit -m "feat(cork): quadrilinear CO2 axis in interpolate_k (geometric log-k)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Thread `co2_vmr` through optical-depth assembly

**Files:**
- Modify: `climt/_components/cork/optics/correlated_k.py:230-257, 323-388`
- Test: `tests/test_cork_co2_interp.py` (add)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_cork_co2_interp.py`:

```python
def test_compute_ck_optical_depth_threads_co2(tmp_path):
    """Higher CO2 VMR yields larger optical depth where the table's CO2 axis
    increases k (additive premixed-bg path)."""
    from climt._components.cork.optics.correlated_k import (
        compute_ck_optical_depth,
    )

    # k grows with CO2 index in band 0; flat elsewhere.
    k = np.full((1, 2, 1, 3, 4, 3, 3), 1e-4)
    k[0, 0, 0, :, :, :, :] = np.array([1e-4, 1e-3, 1e-2])[None, None, None, :]
    table = _build_co2_table(tmp_path, k)

    nlev, ncol = 5, 1
    T = np.full((nlev, ncol), 250.0)
    p = np.full((nlev, ncol), 1e4)
    air = np.ones((1, nlev, ncol)) * 100.0  # kg/m^2 air per layer
    h2o = np.full((nlev, ncol), 1e-3)

    tau_lo = compute_ck_optical_depth(
        table, T, p, air, h2o_vmr=h2o, co2_vmr=np.full((nlev, ncol), 1e-4))
    tau_hi = compute_ck_optical_depth(
        table, T, p, air, h2o_vmr=h2o, co2_vmr=np.full((nlev, ncol), 1e-2))
    # Band 0 optical depth must rise with CO2; band 1 unchanged.
    assert tau_hi[0].sum() > tau_lo[0].sum() * 5
    np.testing.assert_allclose(tau_hi[1], tau_lo[1], rtol=1e-5)
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n climt python -m pytest tests/test_cork_co2_interp.py::test_compute_ck_optical_depth_threads_co2 -v`
Expected: FAIL — `compute_ck_optical_depth() got an unexpected keyword argument 'co2_vmr'`.

- [ ] **Step 3: Add `co2_vmr` to the optical-depth functions**

In `climt/_components/cork/optics/correlated_k.py`:

(a) `compute_ck_optical_depth` — add `co2_vmr=None` and forward it. Replace the signature and dispatch:

```python
def compute_ck_optical_depth(table, T, p, gas_amounts, h2o_vmr=None,
                             co2_vmr=None):
    """Compute optical depths from correlated-k table.

    ... (existing docstring) ...
        co2_vmr: (nlev, ncol) CO2 mole fraction per level/column, required when
            the table has a co2_vmr_grid axis.
    """
    overlap = str(table.get("overlap_method", np.array("additive")))

    if overlap == "esft":
        return _compute_ck_optical_depth_esft(
            table, T, p, gas_amounts, h2o_vmr, co2_vmr)
    else:
        return _compute_ck_optical_depth_additive(
            table, T, p, gas_amounts, h2o_vmr, co2_vmr)
```

(b) `_compute_ck_optical_depth_additive` — accept `co2_vmr=None` and pass the per-level slice into `interpolate_k`. Change the signature and the `interpolate_k` call:

```python
def _compute_ck_optical_depth_additive(table, T, p, gas_amounts, h2o_vmr=None,
                                       co2_vmr=None):
    """Additive overlap: sum optical depths from all gases on the same g-point grid."""
    k_data = table["k_coefficients"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    tau = np.zeros((nband, ngpt, nlev, ncol))
    has_continuum = (
        "continuum_kappa" in table
        and np.asarray(table["continuum_kappa"]).ndim == 4
        and h2o_vmr is not None
    )

    for k_lev in range(nlev):
        x_lev = h2o_vmr[k_lev, :] if h2o_vmr is not None else None
        c_lev = co2_vmr[k_lev, :] if co2_vmr is not None else None
        k_interp = interpolate_k(
            table, T[k_lev, :], p[k_lev, :], h2o_vmr=x_lev, co2_vmr=c_lev)
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    for icol in range(ncol):
                        tau[ib, igp, k_lev, icol] += (
                            k_interp[ig, ib, igp, icol] * gas_amounts[ig, k_lev, icol]
                        )
        if has_continuum:
            cont = interpolate_continuum(table, T[k_lev, :], p[k_lev, :], x_lev)
            if cont is not None:
                for ib in range(nband):
                    for igp in range(ngpt):
                        for icol in range(ncol):
                            tau[ib, igp, k_lev, icol] += (
                                cont[ib, icol] * gas_amounts[0, k_lev, icol]
                            )

    return tau
```

(c) `_compute_ck_optical_depth_esft` — accept and forward `co2_vmr=None` (Earth uses additive, but keep the signature consistent). Change its signature line and the `interpolate_k` call:

```python
def _compute_ck_optical_depth_esft(table, T, p, gas_amounts, h2o_vmr=None,
                                   co2_vmr=None):
    """ESFT overlap: outer product of g-points across gases."""
    k_data = table["k_coefficients"]
    gpoint_weights = table["gpoint_weights"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    ngpt_combined = ngpt ** ngas
    tau = np.zeros((nband, ngpt_combined, nlev, ncol))

    combined_weights = compute_esft_weights(gpoint_weights, ngas)

    for k_lev in range(nlev):
        x_lev = h2o_vmr[k_lev, :] if h2o_vmr is not None else None
        c_lev = co2_vmr[k_lev, :] if co2_vmr is not None else None
        k_interp = interpolate_k(
            table, T[k_lev, :], p[k_lev, :], h2o_vmr=x_lev, co2_vmr=c_lev)

        for ib in range(nband):
            for idx in range(ngpt_combined):
                remainder = idx
                for ig in range(ngas):
                    g_idx = remainder % ngpt
                    remainder //= ngpt
                    for icol in range(ncol):
                        tau[ib, idx, k_lev, icol] += (
                            k_interp[ig, ib, g_idx, icol] * gas_amounts[ig, k_lev, icol]
                        )

    return tau, combined_weights
```

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n climt python -m pytest tests/test_cork_co2_interp.py -v`
Expected: PASS (all three tests).

- [ ] **Step 5: Commit**

```bash
git add climt/_components/cork/optics/correlated_k.py \
        tests/test_cork_co2_interp.py
git commit -m "feat(cork): thread co2_vmr through CK optical-depth assembly

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Component — read `mole_fraction_of_carbon_dioxide_in_air`

The premixed-bg branch gains CO₂ detection (`_has_co2_axis`), declares CO₂ as input when the axis is present, and passes a `(nlev, ncol)` CO₂ VMR (mole fraction is already a VMR — no conversion) into `_correlated_k_optics`.

**Files:**
- Modify: `climt/_components/cork/lw/component.py:38-51, 100-108, 228-262, 404-411`
- Test: `tests/test_cork_lw.py` (add)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_cork_lw.py`:

```python
def test_cork_lw_co2_axis_changes_olr(tmp_path):
    """A premixed-bg table with a CO2 axis reads CO2 from state: more CO2 ->
    less OLR (more trapping)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from scripts.cork_table_builder.netcdf_writer import write_lw_table
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkLongwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    # Tiny 2-band CO2-axis table; band 0 k rises steeply with CO2 index.
    nband, ngpt, nT, nP, nX, nC = 2, 2, 4, 4, 3, 3
    k = np.full((1, nband, ngpt, nT, nP, nX, nC), 1e-5)
    k[0, 0] = np.array([1e-5, 1e-3, 1e-1])[None, None, None, None, :]
    write_lw_table(
        str(tmp_path / "co2tab.nc"),
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(180.0, 320.0, nT),
        log_p_grid=np.linspace(np.log(1e2), np.log(1e5), nP),
        band_edges=np.array([10.0, 1000.0, 3250.0]),
        planck_fraction=np.full((nband, ngpt, nT), 1.0 / ngpt),
        h2o_vmr_grid=np.array([1e-5, 1e-3, 1e-1]),
        co2_vmr_grid=np.array([1e-4, 1e-3, 1e-2]),
        gas_names=("effective",),
    )

    lw = CorkLongwaveRadiation(optics="correlated_k", table=str(tmp_path / "co2tab.nc"))
    grid = get_grid(nx=1, ny=1, nz=20)

    state_lo = get_default_state([lw], grid_state=grid)
    state_lo["specific_humidity"].values[:] = 1e-3
    state_lo["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1e-4
    _, diag_lo = lw(state_lo)
    olr_lo = diag_lo["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    state_hi = get_default_state([lw], grid_state=grid)
    state_hi["specific_humidity"].values[:] = 1e-3
    state_hi["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1e-2
    _, diag_hi = lw(state_hi)
    olr_hi = diag_hi["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    assert olr_hi < olr_lo, f"more CO2 should trap more: lo={olr_lo}, hi={olr_hi}"
```

- [ ] **Step 2: Run to verify it fails**

Run: `conda run -n climt python -m pytest "tests/test_cork_lw.py::test_cork_lw_co2_axis_changes_olr" -v`
Expected: FAIL — CO₂ axis is ignored (OLR identical, assertion fails), or a `co2_vmr` KeyError.

- [ ] **Step 3: Detect the CO₂ axis at init**

In `climt/_components/cork/lw/component.py`, in the `elif optics == "correlated_k":` block, after the `_has_h2o_axis` line, add CO₂ detection:

```python
            _has_h2o_axis = "h2o_vmr_grid" in self._table
            self._has_co2_axis = "co2_vmr_grid" in self._table
            self._fully_premixed = (self._gas_names == ["effective"]) and not _has_h2o_axis
            self._premixed_bg = (
                (self._gas_names == ["effective"] and _has_h2o_axis)
                or str(self._table.get("background_is_premixed", np.array("")))
                .lower() == "true"
            )
```

For the `parmentier` branch (which never sets `_has_co2_axis`), add a default at the top of `__init__` after `self._num_bands = ...`. Simplest: initialize `self._has_co2_axis = False` right after `self._optics_mode = optics`:

```python
        self._optics_mode = optics
        self._has_co2_axis = False
        self._num_bands = 2 if optics == "parmentier" else None
```

- [ ] **Step 4: Declare CO₂ as an input in the premixed-bg branch**

In `input_properties`, inside the `if self._premixed_bg:` block, after the `specific_humidity` entry, add:

```python
            if self._premixed_bg:
                # Background is pre-mixed; H2O VMR is the only live axis.
                props["specific_humidity"] = {
                    "dims": ["mid_levels", "*"],
                    "units": "kg/kg",
                    "alias": "h2o",
                }
                if self._has_co2_axis:
                    props["mole_fraction_of_carbon_dioxide_in_air"] = {
                        "dims": ["mid_levels", "*"],
                        "units": "mole/mole",
                        "alias": "co2",
                    }
```

- [ ] **Step 5: Build the CO₂ VMR array and pass it into the optics**

In `array_call`, the `elif self._premixed_bg:` block builds `h2o_vmr`. Add CO₂ right after it. Find the block ending with the `h2o_vmr = q_h2o / np.maximum(...)` assignment and add:

```python
            elif self._premixed_bg:
                # ... existing h2o_vmr computation ...
                h2o_vmr = q_h2o / np.maximum(
                    q_h2o + (1.0 - q_h2o) * M_H2O, 1e-30
                )
                if self._has_co2_axis:
                    # mole_fraction_of_carbon_dioxide_in_air IS already a VMR.
                    co2_vmr = state["co2"].reshape(nlev, -1)
```

Initialize `co2_vmr = None` next to `h2o_vmr = None` near the top of the `elif self._optics_mode == "correlated_k":` block:

```python
        elif self._optics_mode == "correlated_k":
            ngas = len(self._gas_names)
            gas_amounts = np.zeros((ngas, nlev, ncol))
            h2o_vmr = None
            co2_vmr = None
```

Then change the `_correlated_k_optics` call to forward `co2_vmr`:

```python
            tau, planck_src, surf_src, weights = self._correlated_k_optics(
                T_flat, p_flat, gas_amounts, T_surf_flat, sigma,
                h2o_vmr=h2o_vmr, co2_vmr=co2_vmr,
            )
```

- [ ] **Step 6: Forward `co2_vmr` from `_correlated_k_optics`**

Change `_correlated_k_optics`'s signature and its `compute_ck_optical_depth` call:

```python
    def _correlated_k_optics(self, T, p, gas_amounts, T_surf, sigma,
                             h2o_vmr=None, co2_vmr=None):
        from ..optics.correlated_k import compute_ck_optical_depth

        nlev, ncol = T.shape

        result = compute_ck_optical_depth(
            self._table, T, p, gas_amounts, h2o_vmr=h2o_vmr, co2_vmr=co2_vmr
        )
```

(Leave the rest of `_correlated_k_optics` unchanged.)

- [ ] **Step 7: Run to verify it passes**

Run: `conda run -n climt python -m pytest "tests/test_cork_lw.py::test_cork_lw_co2_axis_changes_olr" -v`
Expected: PASS.

- [ ] **Step 8: Run the full LW + optics suite for regressions**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py tests/test_cork_co2_interp.py -v`
Expected: PASS (all). The existing premixed-bg-without-CO₂ and parmentier paths are unaffected (`_has_co2_axis` False → `co2_vmr` stays None).

- [ ] **Step 9: Commit**

```bash
git add climt/_components/cork/lw/component.py tests/test_cork_lw.py
git commit -m "feat(cork): read CO2 from state as runtime axis in LW component

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Generation — CO₂ sweep axis in `sample_kappa_grid`

Add an optional `co2_vmr_grid` to the sampler. When given, CO₂ is treated as a swept premixed absorber: for each CO₂ node the full κ(T, p, X_H₂O, ν) cube is sampled with that CO₂ amount, producing `(nT, nP, nX_H₂O, nX_CO₂, nNu)`. This requires `h2o_vmr_grid` to also be set (the Earth use-case). H₂O continuum is sampled once (CO₂-independent) by the caller (Task 7), so this function still respects `include_mtckd_continuum`.

**Files:**
- Modify: `scripts/cork_table_builder/kappa_sampling.py:14-105`
- Test: `tests/cork_table_builder/test_kappa_sampling.py` (add; skipped unless linepyline present)

- [ ] **Step 1: Write the failing test**

Append to `tests/cork_table_builder/test_kappa_sampling.py`:

```python
def test_sample_kappa_grid_with_co2_axis():
    """With both H2O and CO2 grids, output is (nT, nP, nX, nC, nNu) and opacity
    rises with CO2 in the CO2 nu2 band."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.cork_table_builder.kappa_sampling import sample_kappa_grid

    T_grid = np.array([250.0, 300.0])
    p_grid = np.array([1e3, 1e5])
    X_h2o = np.array([1e-6, 1e-2])
    X_co2 = np.array([1e-4, 1e-2])
    nu_grid = np.arange(600.0, 700.0, 1.0)  # CO2 nu2 wing

    kappa = sample_kappa_grid(
        background_gas="air",
        absorbers={},  # CO2 supplied via co2_vmr_grid sweep
        h2o_vmr_grid=X_h2o,
        co2_vmr_grid=X_co2,
        T_grid=T_grid,
        p_grid=p_grid,
        nu_grid=nu_grid,
        line_shape="lorentz",
        binning=True,
    )
    assert kappa.shape == (2, 2, 2, 2, 100)
    # More CO2 -> more opacity in this band.
    assert kappa[:, :, :, -1, :].mean() > kappa[:, :, :, 0, :].mean()
```

- [ ] **Step 2: Run to verify it fails (or skips if no linepyline)**

Run: `conda run -n linepyline python -m pytest tests/cork_table_builder/test_kappa_sampling.py::test_sample_kappa_grid_with_co2_axis -v`
Expected: FAIL — `sample_kappa_grid() got an unexpected keyword argument 'co2_vmr_grid'`. (If linepyline is unavailable the module-level `importorskip` skips it; verify the logic via Task 7's dry-run instead.)

- [ ] **Step 3: Add the CO₂ sweep to `sample_kappa_grid`**

In `scripts/cork_table_builder/kappa_sampling.py`, add a `co2_vmr_grid` keyword and a sweep that wraps the existing base sampling. Edit the signature to add `co2_vmr_grid: np.ndarray | None = None,` after `h2o_vmr_grid`, extend the docstring's Returns to mention the 5-D shape, and insert the CO₂ branch near the top of the body — right after `rtm = linepyline.rtm(...)` and the `nu_min/nu_max/dnu/nT/nP/nNu` setup, **before** the `non_h2o = {...}` line:

```python
    nu_min, nu_max = float(nu_grid[0]), float(nu_grid[-1])
    dnu = float(np.diff(nu_grid).mean())
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)

    # CO2 runtime axis: sweep CO2 as a premixed absorber. For each CO2 node,
    # recurse with CO2 folded into `absorbers` (no co2_vmr_grid) and stack along
    # a new axis just before nu. Requires an H2O axis (the Earth use-case).
    if co2_vmr_grid is not None:
        if h2o_vmr_grid is None:
            raise ValueError("co2_vmr_grid requires h2o_vmr_grid")
        nX, nC = len(h2o_vmr_grid), len(co2_vmr_grid)
        out = np.zeros((nT, nP, nX, nC, nNu))
        for iC, co2 in enumerate(co2_vmr_grid):
            sub = sample_kappa_grid(
                background_gas=background_gas,
                absorbers={**absorbers, "CO2": float(co2)},
                T_grid=T_grid,
                p_grid=p_grid,
                nu_grid=nu_grid,
                h2o_vmr_grid=h2o_vmr_grid,
                line_shape=line_shape,
                binning=binning,
                include_mtckd_continuum=include_mtckd_continuum,
                surface_gravity=surface_gravity,
            )  # (nT, nP, nX, nNu)
            out[:, :, :, iC, :] = sub
        return out

    # Pre-compute kappa(T, p, nu) for each non-H2O absorber.
    non_h2o = {s: v for s, v in absorbers.items() if s != "H2O"}
```

(The rest of the function is unchanged.)

- [ ] **Step 4: Run to verify it passes**

Run: `conda run -n linepyline python -m pytest tests/cork_table_builder/test_kappa_sampling.py -v`
Expected: PASS (or SKIP only if linepyline truly unavailable — run this in the `linepyline` env).

- [ ] **Step 5: Commit**

```bash
git add scripts/cork_table_builder/kappa_sampling.py \
        tests/cork_table_builder/test_kappa_sampling.py
git commit -m "feat(cork-gen): CO2 sweep axis in line-by-line kappa sampler

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: Generation driver — 14-band Earth scenario with CO₂ grid

Wire the CO₂ axis through `build_table` (sampling, the `decouple_continuum` continuum sampled once without CO₂, the 5-D moveaxis, and the writer call), and add the production `earth_hifi` scenario with the 14-band partition and 10-node CO₂ grid.

**Files:**
- Modify: `scripts/generate_cork_tables_linepyline.py:36-101, 123-244`
- Test: manual dry-run (full run is Task 9).

- [ ] **Step 1: Add the `earth_hifi` scenario**

In `scripts/generate_cork_tables_linepyline.py`, add to the `SCENARIOS` dict (after the existing `"earth"` entry):

```python
    # Production CO2-adjustable, RRTMG-fidelity Earth LW (sub-project A).
    # 14-band partition isolating the window + far-IR/nu2; decoupled+log-X
    # H2O continuum; CO2 swept on a 10-node 1/3-decade log grid (10-10000 ppm).
    "earth_hifi": dict(
        background_gas="air",
        absorbers={},  # CO2 supplied via co2_vmr_grid sweep
        h2o_vmr_grid=np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]),
        co2_vmr_grid=np.logspace(-5, -2, 10),  # 10,21.5,...,10000 ppm
        T_grid=np.array([50.0, 110.0, 170.0, 230.0, 290.0, 350.0,
                         410.0, 470.0, 530.0, 590.0, 650.0, 1000.0]),
        p_grid=np.array([1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]),
        stellar_spectrum="sun",
        toa_irradiance_W_m2=1361.0,
        mean_molar_mass_g=28.964,
        rayleigh_refractivity=2.77e-4,
        cia_files=[],
        lw_band_edges=np.array([10.0, 250.0, 350.0, 500.0, 630.0, 700.0,
                                800.0, 980.0, 1080.0, 1180.0, 1250.0, 1400.0,
                                1600.0, 1800.0, 3250.0]),
        sw_band_edges=np.array([3250.0, 8000.0, 14000.0, 30000.0]),
    ),
```

- [ ] **Step 2: Thread `co2_vmr_grid` through `build_table`**

In `build_table`, the `_sample` closure must pass `co2_vmr_grid`, the continuum must be sampled **without** the CO₂ axis (it is CO₂-independent), the band-grey continuum reshape must stay 4-D `(nband, nT, nP, nX_H₂O)`, the moveaxis must handle the 5-D κ case, and the writer call must pass `co2_vmr_grid`. Apply these edits:

(a) Replace the `_sample` closure so it forwards the CO₂ grid:

```python
    co2_vmr_grid = cfg.get("co2_vmr_grid")

    def _sample(continuum, with_co2):
        return sample_kappa_grid(
            background_gas=cfg["background_gas"],
            absorbers=cfg["absorbers"],
            h2o_vmr_grid=cfg["h2o_vmr_grid"],
            co2_vmr_grid=(co2_vmr_grid if with_co2 else None),
            T_grid=cfg["T_grid"],
            p_grid=cfg["p_grid"],
            nu_grid=nu_grid,
            line_shape=line_shape,
            binning=binning,
            include_mtckd_continuum=continuum,
        )
```

(b) In the `if decouple_continuum:` branch, sample lines **with** CO₂ but the continuum **without** CO₂ (it has no CO₂ axis). The continuum is the difference of two no-CO₂ samples (lines+cont) − (lines), each at the table's reference CO₂ of 0 — i.e. CO₂-independent by construction:

```python
    continuum_band = None
    if decouple_continuum:
        # Line-only k-distribution, CO2-swept if requested.
        kappa = _sample(continuum=False, with_co2=(co2_vmr_grid is not None))
        # MT_CKD H2O continuum is CO2-INDEPENDENT: sample once with NO CO2 axis.
        cont_lines = _sample(continuum=False, with_co2=False)   # (nT,nP,nX,nNu)
        cont_total = _sample(continuum=True, with_co2=False)
        kappa_cont = cont_total - cont_lines
        np.clip(kappa_cont, 0.0, None, out=kappa_cont)
        nband = len(band_edges) - 1
        nT, nP, nX = (len(cfg["T_grid"]), len(cfg["p_grid"]),
                      len(cfg["h2o_vmr_grid"]))
        continuum_band = np.zeros((nband, nT, nP, nX))
        for b in range(nband):
            m = (nu_grid >= band_edges[b]) & (nu_grid < band_edges[b + 1])
            continuum_band[b] = kappa_cont[:, :, :, m].mean(axis=-1)
    else:
        kappa = _sample(continuum=include_mtckd_continuum,
                        with_co2=(co2_vmr_grid is not None))
```

(c) After `kappa_to_k_coeffs` and `k_coeffs = k_coeffs[np.newaxis]`, handle the 5-D κ case in the moveaxis block:

```python
    print(f"[{scenario_name}/{kind}] building k-distribution (ngpt={ngpt})")
    k_coeffs, gpt_weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt)
    k_coeffs = k_coeffs[np.newaxis]

    if kappa.ndim == 5:
        # k_coeffs: (1, nT, nP, nX, nC, nband, ngpt)
        #        -> (1, nband, ngpt, nT, nP, nX, nC)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4, 5, 6), (3, 4, 5, 6, 1, 2))
    elif kappa.ndim == 4:
        # k_coeffs: (1, nT, nP, nX, nband, ngpt) -> (1, nband, ngpt, nT, nP, nX)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4, 5), (3, 4, 5, 1, 2))
    else:
        # k_coeffs: (1, nT, nP, nband, ngpt) -> (1, nband, ngpt, nT, nP)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4), (3, 4, 1, 2))
```

(d) In the `if kind == "lw":` branch, pass `co2_vmr_grid` to `write_lw_table`:

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
            co2_vmr_grid=co2_vmr_grid,
            source=f"linepyline:{scenario_name}",
            continuum_kappa=continuum_band,
        )
```

(e) Note: the CIA block (`if cfg["cia_files"]:`) only handles κ of ndim 3 or 4. `earth_hifi` has `cia_files=[]`, so it is skipped — leave that block as-is.

- [ ] **Step 3: Dry-run the moveaxis logic without linepyline**

Create `scripts/experiments/dryrun_earth_hifi_shapes.py` to validate axis bookkeeping with a fake κ (no linepyline):

```python
"""Validate earth_hifi k_coefficients axis ordering with synthetic kappa."""
import os, sys
import numpy as np

_HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "scripts"))
from cork_table_builder.k_distribution import kappa_to_k_coeffs

nT, nP, nX, nC, ngpt = 12, 8, 7, 10, 8
band_edges = np.array([10.0, 250.0, 350.0, 500.0, 630.0, 700.0, 800.0,
                       980.0, 1080.0, 1180.0, 1250.0, 1400.0, 1600.0,
                       1800.0, 3250.0])
nband = len(band_edges) - 1
nu = np.arange(band_edges[0], band_edges[-1], 1.0)
kappa = np.abs(np.random.RandomState(0).standard_normal(
    (nT, nP, nX, nC, len(nu)))) * 1e-4

k, w = kappa_to_k_coeffs(kappa, nu, band_edges, ngpt)
k = k[np.newaxis]
k = np.moveaxis(k, (1, 2, 3, 4, 5, 6), (3, 4, 5, 6, 1, 2))
assert k.shape == (1, nband, ngpt, nT, nP, nX, nC), k.shape
print("OK", k.shape)
```

Run: `conda run -n climt python scripts/experiments/dryrun_earth_hifi_shapes.py`
Expected: `OK (1, 14, 8, 12, 8, 7, 10)`.

- [ ] **Step 4: Commit**

```bash
git add scripts/generate_cork_tables_linepyline.py \
        scripts/experiments/dryrun_earth_hifi_shapes.py
git commit -m "feat(cork-gen): earth_hifi 14-band CO2-grid scenario + build wiring

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8: Performance — consolidate LW transport into a single `@njit(parallel=True)` kernel

Move the per-band/per-g-point accumulation (currently a Python loop over bands × g-points calling a single-g-point njit, plus Python broadband summation) into one compiled kernel that loops over a parallel column axis and accumulates weighted fluxes inside. Must reproduce current fluxes **bit-for-bit** on existing tables. Accumulation order is preserved (g ascending into `up_band`/`down_band`; b ascending into broadband) to guarantee identical floating-point results.

**Files:**
- Modify: `climt/_components/cork/lw/kernels.py`
- Test: `tests/test_lw_kernel_consolidation.py` (create)

- [ ] **Step 1: Capture a bit-for-bit golden snapshot of current output**

Create `tests/test_lw_kernel_consolidation.py`. This test computes the current `lw_transport` output on a fixed synthetic input, then (after the refactor) asserts identical results. To make it self-contained and refactor-proof, the golden values are generated inside the test from the SAME inputs via a frozen reference implementation embedded in the test:

```python
import numpy as np


def _reference_lw_transport(tau, planck_source, surface_source, emissivity, weights):
    """Frozen copy of the ORIGINAL python-loop transport (pre-consolidation),
    used as the bit-for-bit oracle."""
    D = 1.66
    nband, ngpt, nlev, ncol = tau.shape
    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))
    for b in range(nband):
        for g in range(ngpt):
            up = np.zeros((nlev + 1, ncol))
            dn = np.zeros((nlev + 1, ncol))
            for i in range(ncol):
                up[0, i] = emissivity[b, i] * surface_source[b, g, i]
                for k in range(nlev):
                    t = np.exp(-D * tau[b, g, k, i])
                    up[k + 1, i] = up[k, i] * t + planck_source[b, g, k, i] * (1.0 - t)
                dn[nlev, i] = 0.0
                for k in range(nlev - 1, -1, -1):
                    t = np.exp(-D * tau[b, g, k, i])
                    dn[k, i] = dn[k + 1, i] * t + planck_source[b, g, k, i] * (1.0 - t)
            w = weights[b, g]
            for k in range(nlev + 1):
                for i in range(ncol):
                    up_band[b, k, i] += w * up[k, i]
                    down_band[b, k, i] += w * dn[k, i]
    up_broad = up_band.sum(axis=0)
    down_broad = down_band.sum(axis=0)
    return up_band, down_band, up_broad, down_broad


def _make_inputs(seed=0):
    rng = np.random.RandomState(seed)
    nband, ngpt, nlev, ncol = 14, 8, 20, 3
    tau = rng.uniform(0.0, 2.0, size=(nband, ngpt, nlev, ncol))
    planck = rng.uniform(50.0, 400.0, size=(nband, ngpt, nlev, ncol))
    surf = rng.uniform(50.0, 400.0, size=(nband, ngpt, ncol))
    emis = rng.uniform(0.8, 1.0, size=(nband, ncol))
    weights = rng.uniform(0.1, 0.5, size=(nband, ngpt))
    return tau, planck, surf, emis, weights


def test_consolidated_kernel_matches_reference_bitwise():
    from climt._components.cork.lw.kernels import lw_transport

    tau, planck, surf, emis, weights = _make_inputs()
    T = np.zeros((20, 3)); T_surf = np.zeros(3); sigma = 5.67e-8
    ub_r, db_r, ubr_r, dbr_r = _reference_lw_transport(
        tau, planck, surf, emis, weights)
    ub, db, ubr, dbr = lw_transport(
        T, T_surf, tau, planck, surf, emis, weights, sigma,
        diagnostics_level=0)
    np.testing.assert_array_equal(ub, ub_r)
    np.testing.assert_array_equal(db, db_r)
    np.testing.assert_array_equal(ubr, ubr_r)
    np.testing.assert_array_equal(dbr, dbr_r)


def test_diagnostics_level_one_still_returns_diag():
    from climt._components.cork.lw.kernels import lw_transport

    tau, planck, surf, emis, weights = _make_inputs(seed=1)
    T = np.zeros((20, 3)); T_surf = np.zeros(3); sigma = 5.67e-8
    out = lw_transport(T, T_surf, tau, planck, surf, emis, weights, sigma,
                       diagnostics_level=1)
    assert len(out) == 5
    diag = out[4]
    assert diag["transmittance"].shape == (14, 8, 20, 3)
    assert diag["up_per_gpoint"].shape == (14, 8, 21, 3)
    # Diagnostic weighted up-flux sums to up_band.
    np.testing.assert_allclose(
        diag["up_per_gpoint"].sum(axis=1), out[0], rtol=1e-12)
```

- [ ] **Step 2: Run to verify the bitwise test PASSES against the current code**

Run: `conda run -n climt python -m pytest tests/test_lw_kernel_consolidation.py -v`
Expected: PASS — the current implementation already matches the frozen reference. This locks in behavior **before** refactoring (a guard, not a red test).

- [ ] **Step 3: Replace `kernels.py` with the consolidated kernel**

Rewrite `climt/_components/cork/lw/kernels.py` as:

```python
# climt/_components/cork/lw/kernels.py
import numpy as np

from ..common import njit, prange

DIFFUSIVITY_FACTOR = 1.66


@njit(parallel=True)
def _lw_transport_kernel(
    tau, planck_source, surface_source, emissivity, weights,
    up_band, down_band, up_broad, down_broad,
    diag_trans, diag_up_gpt, diag_dn_gpt, want_diag,
):
    """Consolidated multi-band, multi-g-point LW transport.

    Loops over columns in parallel; for each (band, g-point) runs the up/down
    diffusivity sweeps and accumulates weighted fluxes into up_band/down_band
    inside the compiled kernel. Accumulation order (g ascending, then b
    ascending for broadband) matches the original python loops bit-for-bit.

    When want_diag != 0, also fills per-layer transmittance and weighted
    per-g-point fluxes. All output arrays are preallocated by the caller.
    """
    nband, ngpt, nlev, ncol = tau.shape
    for i in prange(ncol):
        for k in range(nlev + 1):
            up_broad[k, i] = 0.0
            down_broad[k, i] = 0.0
        for b in range(nband):
            for k in range(nlev + 1):
                up_band[b, k, i] = 0.0
                down_band[b, k, i] = 0.0
            for g in range(ngpt):
                w = weights[b, g]
                # Upward sweep: surface -> TOA
                up_prev = emissivity[b, i] * surface_source[b, g, i]
                up_band[b, 0, i] += w * up_prev
                if want_diag != 0:
                    diag_up_gpt[b, g, 0, i] = w * up_prev
                for k in range(nlev):
                    trans = np.exp(-DIFFUSIVITY_FACTOR * tau[b, g, k, i])
                    up_cur = up_prev * trans + planck_source[b, g, k, i] * (1.0 - trans)
                    up_band[b, k + 1, i] += w * up_cur
                    if want_diag != 0:
                        diag_trans[b, g, k, i] = trans
                        diag_up_gpt[b, g, k + 1, i] = w * up_cur
                    up_prev = up_cur
                # Downward sweep: TOA -> surface
                dn_prev = 0.0  # down_band[b, nlev, i] stays 0
                if want_diag != 0:
                    diag_dn_gpt[b, g, nlev, i] = 0.0
                for k in range(nlev - 1, -1, -1):
                    trans = np.exp(-DIFFUSIVITY_FACTOR * tau[b, g, k, i])
                    dn_cur = dn_prev * trans + planck_source[b, g, k, i] * (1.0 - trans)
                    down_band[b, k, i] += w * dn_cur
                    if want_diag != 0:
                        diag_dn_gpt[b, g, k, i] = w * dn_cur
                    dn_prev = dn_cur
            for k in range(nlev + 1):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]


def lw_transport(
    T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma,
    diagnostics_level=0,
):
    """Multi-band, multi-g-point LW radiative transfer (consolidated kernel).

    Args/Returns: unchanged from the previous implementation. See
    _lw_transport_kernel for the compiled core.
    """
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))
    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))

    want_diag = 1 if diagnostics_level >= 1 else 0
    if want_diag:
        diag_trans = np.zeros((nband, ngpt, nlev, ncol))
        diag_up_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))
        diag_dn_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))
    else:
        # Dummy minimal arrays (numba needs concrete dtypes/shapes).
        diag_trans = np.zeros((1, 1, 1, 1))
        diag_up_gpt = np.zeros((1, 1, 1, 1))
        diag_dn_gpt = np.zeros((1, 1, 1, 1))

    _lw_transport_kernel(
        tau, planck_source, surface_source, emissivity, weights,
        up_band, down_band, up_broad, down_broad,
        diag_trans, diag_up_gpt, diag_dn_gpt, want_diag,
    )

    if want_diag:
        diag = {
            "transmittance": diag_trans,
            "up_per_gpoint": diag_up_gpt,
            "down_per_gpoint": diag_dn_gpt,
        }
        return up_band, down_band, up_broad, down_broad, diag

    return up_band, down_band, up_broad, down_broad
```

- [ ] **Step 4: Run the consolidation tests to verify bit-for-bit reproduction**

Run: `conda run -n climt python -m pytest tests/test_lw_kernel_consolidation.py -v`
Expected: PASS — `assert_array_equal` confirms the consolidated kernel matches the frozen oracle exactly, and the diagnostics path returns correctly shaped arrays.

- [ ] **Step 5: Run the full LW suite for regressions**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py tests/test_cork_sw.py -v`
Expected: PASS (the per-band-sum and isothermal tests rely on the kernel; they must still pass).

- [ ] **Step 6: Commit**

```bash
git add climt/_components/cork/lw/kernels.py \
        tests/test_lw_kernel_consolidation.py
git commit -m "perf(cork): consolidate LW transport into single njit(parallel) kernel

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 9: Generate the production table (documented offline run)

This is the one multi-hour step. It runs line-by-line in the `linepyline` conda env on the 4-D `(T, p, X_H₂O, X_CO₂)` grid; offline cost is ~10× the single-CO₂ Earth table. Run on an adequate machine (overnight). **Not test-driven** — gated downstream by Tasks 10–11.

**Files:**
- Create: `climt/_data/cork/correlated_k/earth_low_res_lw.nc` (regenerated in place)
- Preserve: rename current shipped files as "before" references for sub-project B.

- [ ] **Step 1: Preserve the current shipped table as a "before" reference**

```bash
cd /Users/joymonteiro/github/climt
cp climt/_data/cork/correlated_k/earth_low_res_lw.nc \
   climt/_data/cork/correlated_k/earth_low_res_lw_4band_ngpt2_before.nc
```

(`earth_low_res_lw_8gpt.nc` already exists as the 8-g-point "before" variant — leave it.)

- [ ] **Step 2: Run the production generation (linepyline env, overnight)**

```bash
conda run -n linepyline python scripts/generate_cork_tables_linepyline.py \
    --scenario earth_hifi --kind lw \
    --output climt/_data/cork/correlated_k/earth_hifi_lw.nc \
    --ngpt 8 --dnu 0.1 --line-shape pseudovoigt --decouple-continuum
```

Expected stdout (final line): `k_coefficients.shape=(1, 14, 8, 12, 8, 7, 10)`.

- [ ] **Step 3: Sanity-load the generated table**

Run:
```bash
conda run -n climt python -c "
from climt._components.cork.optics.correlated_k import load_k_table
t = load_k_table('climt/_data/cork/correlated_k/earth_hifi_lw.nc')
print('k', t['k_coefficients'].shape)
print('co2_vmr_grid', t['co2_vmr_grid'])
print('h2o_vmr_grid', t['h2o_vmr_grid'])
print('continuum_kappa', t['continuum_kappa'].shape)
assert t['k_coefficients'].ndim == 7
assert t['continuum_kappa'].ndim == 4
print('OK')
"
```
Expected: shapes `(1, 14, 8, 12, 8, 7, 10)` and `(14, 12, 8, 7)`, `co2_vmr_grid` of 10 values from 1e-5 to 1e-2, then `OK`.

- [ ] **Step 4: Promote to the default name once Tasks 10–11 pass** (deferred to Task 11)

Do NOT overwrite `earth_low_res_lw.nc` yet — the default switch is gated by validation. Keep the new table as `earth_hifi_lw.nc` until Task 11.

- [ ] **Step 5: Commit the generated artifact + the preserved reference**

```bash
git add climt/_data/cork/correlated_k/earth_hifi_lw.nc \
        climt/_data/cork/correlated_k/earth_low_res_lw_4band_ngpt2_before.nc
git commit -m "data(cork): generated CO2-adjustable 14-band Earth-LW table (earth_hifi_lw)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 10: Validation — LBL fidelity across a CO₂×H₂O grid + CO₂-interp accuracy

Extend the single-column evaluator to sweep CO₂ and confirm `|LBL − CORK|` holds away from 376 ppm, and verify off-node CO₂ interpolation fidelity (the A5 check that confirms log-k vs linear-k).

**Files:**
- Modify: `scripts/experiments/eval_band_structure.py`
- Create: `scripts/experiments/eval_co2_interp_accuracy.py`

- [ ] **Step 1: Add a CO₂ argument to `eval_band_structure.py`**

In `scripts/experiments/eval_band_structure.py`, make `cork_olr` already accepts `co2`; add an optional CLI 3rd-positional CO₂-override and print it. Change `main()` to read an optional `--co2` (parse from argv) and pass it into `cork_olr` for both moist/dry cases instead of the profile's baked CO₂:

```python
def main():
    table = sys.argv[1]
    lbl_suffix = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("co2=") else ""
    co2_override = None
    for a in sys.argv[2:]:
        if a.startswith("co2="):
            co2_override = float(a.split("=", 1)[1])
    d = np.load(os.path.join(_DBG, "forward_profile.npz"), allow_pickle=True)
    p_mid, T, q_moist, Ts, co2 = (d["p_mid_Pa"], d["T_ref"],
                                  d["q_ref_moist"], float(d["T_surf"]),
                                  float(d["CO2"]))
    if co2_override is not None:
        co2 = co2_override
    grid = get_grid(nx=1, ny=1, nz=NZ)
    # ... rest unchanged (uses `co2`) ...
```

Run (after Task 9): `conda run -n climt python scripts/experiments/eval_band_structure.py earth_hifi_lw co2=0.000376`
Expected: per-band table prints; total `LBL-CORK` within a few W/m² (target < 5, matching RRTMG's own LBL agreement of 4.1 moist / 7.3 dry). Repeat with `co2=0.00112` (1120 ppm) and confirm the residual stays bounded.

- [ ] **Step 2: Write the CO₂-interpolation accuracy probe**

Create `scripts/experiments/eval_co2_interp_accuracy.py`. It compares the table interpolated at OFF-node CO₂ (e.g. 400, 2000 ppm) to a freshly LBL-sampled κ→k at the same CO₂, isolating interpolation error. This requires linepyline. Skeleton:

```python
"""A5 CO2-interpolation accuracy check.

For each off-node CO2 (400, 2000 ppm), compare the shipped table's
quadrilinearly-interpolated band-mean optical depth on the fixed exp #16
profile against a fresh LBL sample at that exact CO2. Reports per-band and
total OLR error and the geometric-vs-linear-k choice. Run in linepyline env.
"""
import os, sys
import numpy as np

# Reuse forward_profile.npz + eval_band_structure.cork_olr for the CORK side.
# For the LBL side, call scripts/cork_table_builder/kappa_sampling.sample_kappa_grid
# at the single (profile) column with the exact CO2, k-distribute on the table's
# band edges, run lw_transport, and compare OLR.
# Print: CO2(ppm) | CORK_interp_OLR | LBL_OLR | diff, for log-k and linear-k
# (toggle correlated_k._CO2_INTERP_LOGK) to confirm log-k is the better choice.
...
```

Implement it concretely by mirroring `eval_band_structure.cork_olr` for the CORK side and `cork_table_builder` for the LBL side; run:

```bash
conda run -n linepyline python scripts/experiments/eval_co2_interp_accuracy.py earth_hifi_lw
```
Expected: at 400 and 2000 ppm, log-k interpolation error < linear-k error and within a few W/m² of LBL. Record numbers in the investigation log (Step 4).

- [ ] **Step 3: RCE at ≥2 CO₂ values vs RRTMG (moist + dry)**

Use the existing RCE drivers `scripts/experiments/rce_moist_cork_vs_rrtmg.py` and `rce_dry_cork_vs_rrtmg.py`, pointing them at `earth_hifi_lw` and setting `mole_fraction_of_carbon_dioxide_in_air` to 280 ppm and 1120 ppm. Confirm converged surface-temperature bias vs RRTMG is a few K (design target: moist +2.7 K, dry +1.0 K at the reference CO₂). Also run the **true baseline** with the actual current default (`earth_low_res_lw`, ngpt=2) for the honest "before" number.

```bash
conda run -n climt python scripts/experiments/rce_moist_cork_vs_rrtmg.py --table earth_hifi_lw --co2 280e-6
conda run -n climt python scripts/experiments/rce_moist_cork_vs_rrtmg.py --table earth_hifi_lw --co2 1120e-6
conda run -n climt python scripts/experiments/rce_dry_cork_vs_rrtmg.py   --table earth_hifi_lw --co2 280e-6
# True baseline:
conda run -n climt python scripts/experiments/rce_moist_cork_vs_rrtmg.py --table earth_low_res_lw
```

(If those drivers don't yet accept `--table`/`--co2`, add the argparse flags — small, mechanical — and pass through to `CorkLongwaveRadiation(table=...)` and the state's CO₂ field.)

- [ ] **Step 4: Log results to the investigation plan**

Append a dated results block (CO₂×H₂O LBL-CORK table, CO₂-interp errors, RCE biases, true baseline number) to `docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md` (the active investigation log, per project convention). Include the final log-k vs linear-k decision.

- [ ] **Step 5: Commit the validation tooling + logged results**

```bash
git add scripts/experiments/eval_band_structure.py \
        scripts/experiments/eval_co2_interp_accuracy.py \
        scripts/experiments/rce_moist_cork_vs_rrtmg.py \
        scripts/experiments/rce_dry_cork_vs_rrtmg.py \
        docs/superpowers/plans/2026-05-16-cork-co2-band-refinement.md
git commit -m "test(cork): CO2xH2O LBL fidelity + CO2-interp accuracy + RCE validation

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 11: Flip the default to the new table

Gated by Task 10 passing (RCE bias a few K; CO₂-interp fidelity confirmed). Promotes `earth_hifi_lw.nc` to the canonical `earth_low_res_lw.nc`, adds the loader/interp/component coverage as a shipped-table integration test, and verifies docs/examples still work with CO₂ now settable.

**Files:**
- Modify: `climt/_data/cork/correlated_k/earth_low_res_lw.nc` (replace in place)
- Modify: `climt/_data/cork/correlated_k/MANIFEST.md`
- Test: `tests/test_cork_lw.py` (add shipped-table test)

- [ ] **Step 1: Replace the shipped default in place**

```bash
cp climt/_data/cork/correlated_k/earth_hifi_lw.nc \
   climt/_data/cork/correlated_k/earth_low_res_lw.nc
```

- [ ] **Step 2: Write the shipped-table integration test**

Append to `tests/test_cork_lw.py`:

```python
def test_shipped_earth_lw_is_co2_adjustable():
    """The default earth_low_res_lw exposes a CO2 axis and responds to CO2."""
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkLongwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    assert lw._has_co2_axis is True
    assert lw.num_longwave_bands == 14

    grid = get_grid(nx=1, ny=1, nz=30)
    s_lo = get_default_state([lw], grid_state=grid)
    s_lo["specific_humidity"].values[:] = 5e-3
    s_lo["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 280e-6
    _, d_lo = lw(s_lo)
    olr_lo = d_lo["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    s_hi = get_default_state([lw], grid_state=grid)
    s_hi["specific_humidity"].values[:] = 5e-3
    s_hi["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1120e-6
    _, d_hi = lw(s_hi)
    olr_hi = d_hi["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    assert olr_hi < olr_lo  # CO2 quadrupling traps more LW
    assert not np.isnan(olr_lo) and not np.isnan(olr_hi)
```

- [ ] **Step 3: Run the integration test + full CORK suite**

Run: `conda run -n climt python -m pytest tests/test_cork_lw.py tests/test_cork_co2_interp.py tests/test_lw_kernel_consolidation.py tests/cork_table_builder/ -v`
Expected: PASS (all). The 14-band count and CO₂ response confirm the default is live.

- [ ] **Step 4: Update the MANIFEST and check docs/examples**

In `climt/_data/cork/correlated_k/MANIFEST.md`, update the `earth_low_res_lw` entry: 14 bands, ngpt=8, decoupled+log-X H₂O continuum, runtime CO₂ axis 10–10000 ppm (10 nodes), source `linepyline:earth_hifi`. Note `earth_low_res_lw_4band_ngpt2_before.nc` and `earth_low_res_lw_8gpt.nc` as preserved "before" references.

Grep for examples/docs referencing the table and confirm they still run (CO₂ is now settable but defaults to 330 ppm from the state, so existing examples need no change):

```bash
grep -rn "earth_low_res_lw" examples/ docs/ 2>/dev/null
```

- [ ] **Step 5: Refresh the knowledge graph (per CLAUDE.md)**

```bash
graphify update . && conda run -n climt python scripts/augment_graph.py
```

- [ ] **Step 6: Commit**

```bash
git add climt/_data/cork/correlated_k/earth_low_res_lw.nc \
        climt/_data/cork/correlated_k/MANIFEST.md \
        tests/test_cork_lw.py graphify-out/
git commit -m "feat(cork): make CO2-adjustable 14-band table the default earth_low_res_lw

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-review against the spec

**Spec coverage:**
- Goal 1 (decoupled + log-X continuum) — already implemented in the codebase (`--decouple-continuum`, `interpolate_continuum`); Task 7 wires it for `earth_hifi` and keeps the continuum CO₂-independent (4-D). ✓
- Goal 2 (14-band partition) — Task 7 `earth_hifi` `lw_band_edges` (exact spec edges). ✓
- Goal 3 (CO₂-adjustable 10–10000 ppm, single table) — Tasks 1–7 (format, loader, interp, optics, component, sampler, scenario); 10-node `logspace(-5,-2,10)`. ✓
- Goal 4 (becomes default, matches RRTMG to a few K) — Tasks 9–11 (generate, validate, flip). ✓
- Design "true premixed mixture per (X_H₂O,X_CO₂) node" — Task 6 samples CO₂ as a premixed absorber per node (no gas-overlap approx). ✓
- Design "quadrilinear; log-X_CO₂; default log-k" — Task 3, with `_CO2_INTERP_LOGK` switch for the A5 confirmation. ✓ (Spec internal inconsistency on log-k vs linear-k flagged and resolved in Background.)
- Design "continuum has no CO₂ axis" — Tasks 1, 7 keep `continuum_kappa` 4-D. ✓
- Code change 6 (backward compatibility) — new axes optional/by-presence; Tasks 2–5 guard on `ndim`/`"co2_vmr_grid" in table`. ✓
- Performance pass (single njit(parallel), bit-for-bit) — Task 8 with `assert_array_equal` oracle. ✓
- Validation gates (true baseline, LBL CO₂×H₂O, CO₂-interp accuracy, RCE ≥2 CO₂, test suite) — Task 10; default switch gated in Task 11. ✓
- Naming (regenerate canonical name, preserve "before" refs) — Tasks 9, 11. ✓

**Non-goals respected:** no SW table/axis, no line-side log-X, no tropopause diagnostic, no further band optimization, no Titan/TRAPPIST regen. ✓

**Type consistency:** `co2_vmr` is the parameter name in `interpolate_k`, `compute_ck_optical_depth`, `_compute_ck_optical_depth_additive/_esft`, and `_correlated_k_optics`; `co2_vmr_grid` is the netCDF var / writer kwarg / scenario key; `_has_co2_axis` is the component flag; `_CO2_INTERP_LOGK` the interp-mode constant. Consistent across tasks. ✓

**Placeholder scan:** Tasks 1–8 contain complete code. Task 9 is a documented offline run (no code). Task 10's `eval_co2_interp_accuracy.py` is given as a skeleton with a precise spec rather than full code — this is the one intentional non-complete artifact, because its LBL side depends on the exact `forward_profile.npz` contents and linepyline call signatures that must be confirmed at implementation time; the engineer should complete it by mirroring the cited existing scripts.
