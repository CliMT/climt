# Titan & TRAPPIST-1e Opacity Tables (via linepyline) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a linepyline-based offline pipeline that produces climt-format picket-fence correlated-k tables for TRAPPIST-1e (THAI Hab1 + Hab2 scenarios) and Titan, so students can run the picket-fence radiation scheme on those atmospheres and compare against SOCRATES.

**Architecture:** A pure-python offline script `scripts/generate_pf_tables_linepyline.py` calls `linepyline.rtm.get_kappa_hitran` (+ MT_CKD continuum for H2O, + HITRAN CIA for Titan N₂–N₂ / N₂–CH₄) on a (T, p) grid at line-by-line resolution, bins the resulting κ(ν) into picket-fence bands, builds a per-band k-distribution sampled at `ngpt` Gauss–Legendre g-points, and writes a netCDF file matching the schema consumed by `climt._components.picket_fence.optics.correlated_k.load_k_table`. The script is run from the `linepyline` conda env (not a climt runtime dep). Outputs are committed under `climt/_data/picket_fence/correlated_k/` and smoke-tested from the `climt` env. THAI Hab1 uses the H₂O-VMR axis (trilinear in T, log p, log X_H2O) so the same table can drive runs over a wide humidity range; Hab2 and Titan are pre-mixed (no H₂O axis). Titan additionally pulls HITRAN CIA cross-sections — without CIA, Titan's LW window opacity is wrong by orders of magnitude.

**Tech Stack:** Python 3.12, numpy, xarray, scipy.io.netcdf, linepyline (HITRAN 2024 lines + MT_CKD 4.3), HITRAN CIA flat files, pytest.

**Caveats baked into the plan:**
- Linepyline does not natively support CIA; Titan needs a small CIA loader that reads HITRAN CIA `.cia` text files (downloaded once from hitran.org/cia/).
- Linepyline does not do shortwave scattering. SW absorption tables built here contain only κ(ν) from gas lines; Rayleigh is added per-band as a separate coefficient (already a slot in climt's schema). Cloud scattering is added by the SW component at runtime.
- THAI Ben1/Ben2 differ from Hab1/Hab2 only in surface boundary (aquaplanet vs ocean depth) — the gas-phase opacity tables are identical, so Hab1/Hab2 tables cover all four THAI scenarios.

---

## File Structure

```
scripts/
  generate_pf_tables_linepyline.py     # CREATE — CLI driver (one scenario per invocation)
  pf_table_builder/
    __init__.py                        # CREATE
    kappa_sampling.py                  # CREATE — wraps linepyline κ on a (T,p[,X]) grid
    k_distribution.py                  # CREATE — band-bin + sort + Gauss-Legendre quadrature
    planck_fraction.py                 # CREATE — LW Planck weights per (band, gpoint, T)
    solar_source.py                    # CREATE — SW solar source + Rayleigh per (band, gpoint)
    cia.py                             # CREATE — HITRAN CIA loader (Titan only)
    netcdf_writer.py                   # CREATE — climt-schema writer
climt/_data/picket_fence/correlated_k/
  trappist1e_hab1_lw.nc                # CREATE (with H2O VMR axis)
  trappist1e_hab1_sw.nc                # CREATE (with H2O VMR axis)
  trappist1e_hab2_lw.nc                # CREATE (premixed pure CO2)
  trappist1e_hab2_sw.nc                # CREATE (premixed pure CO2)
  titan_lw.nc                          # CREATE (premixed; CIA included)
  titan_sw.nc                          # CREATE (premixed; Sun-illuminated)
  MANIFEST.md                          # MODIFY — append new rows
climt/_data/picket_fence/cia/
  N2-N2_2018.cia                       # CREATE (committed from HITRAN download)
  N2-CH4_2018.cia                      # CREATE
  CH4-N2_2018.cia                      # CREATE (sometimes only one of the two pairs is in HITRAN; document)
tests/
  test_correlated_k_tables.py          # MODIFY — extend SHIPPED list
  test_picket_fence_titan_trappist.py  # CREATE — end-to-end smoke run
docs/
  radiative_transfer/
    table_generation.rst               # MODIFY — add linepyline section
```

**Notes:**
- `pf_table_builder/` is a *package next to the script*, not inside `climt/` — keep climt's runtime dep graph clean. The runtime only consumes the resulting `.nc` files.
- All `_data/picket_fence/correlated_k/*.nc` files are committed as binary. Each is small (< 200 KB at 4 bands × 2 gpts × ~10×15 (T,p) × optional 5 X-points).
- The `cia/` subdirectory is only used by the offline build script — climt itself does not need to load it at runtime (the CIA contribution is baked into the table's k-coefficients).

---

## Run-environment convention

All script invocations in this plan are run from the `linepyline` conda env, except tests which run from `climt`:

```bash
# Build (heavy, offline):
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py …

# Test (cheap, in climt):
cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest …
```

---

### Task 1: Kappa sampling utility — linepyline wrapper over a (T, p) grid

**Files:**
- Create: `scripts/pf_table_builder/__init__.py`
- Create: `scripts/pf_table_builder/kappa_sampling.py`
- Create: `tests/pf_table_builder/test_kappa_sampling.py`

The job: for one scenario (background gas + absorber dict), produce a 4D array `kappa[iT, iP, iNu]` of total mass absorption coefficient (m²/kg) on a fixed (T grid, log-p grid, wavenumber grid).

- [x] **Step 1: Create the package**

Create `scripts/pf_table_builder/__init__.py` with a single line:

```python
"""Offline utilities to convert linepyline line-by-line opacity into picket-fence k-tables."""
```

- [x] **Step 2: Write the failing test**

Create `tests/pf_table_builder/__init__.py` (empty file).

Create `tests/pf_table_builder/test_kappa_sampling.py`:

```python
"""Tests for the linepyline kappa sampler.

These tests are skipped unless linepyline is importable in the active env.
"""
import numpy as np
import pytest

linepyline = pytest.importorskip("linepyline")


def test_sample_kappa_grid_shapes(monkeypatch):
    """sample_kappa_grid returns (nT, nP, nNu) for a single absorber."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.kappa_sampling import sample_kappa_grid

    T_grid = np.array([200.0, 250.0, 300.0])
    p_grid = np.array([1e3, 1e4, 1e5])
    nu_grid = np.arange(500.0, 600.0, 1.0)  # 100 points, narrow LW band

    # Pure CO2 atmosphere, no background gas
    kappa = sample_kappa_grid(
        background_gas=None,
        absorbers={"CO2": 1.0},
        T_grid=T_grid,
        p_grid=p_grid,
        nu_grid=nu_grid,
        line_shape="lorentz",
        binning=True,
    )
    assert kappa.shape == (3, 3, 100)
    assert (kappa >= 0).all()
    # CO2 ν₂ bending mode is in this band — kappa should be non-trivial
    assert kappa.max() > 1e-5


def test_sample_kappa_grid_with_h2o_vmr_axis():
    """sample_kappa_grid with X axis returns (nT, nP, nX, nNu)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.kappa_sampling import sample_kappa_grid

    T_grid = np.array([250.0, 300.0])
    p_grid = np.array([1e3, 1e5])
    X_h2o = np.array([1e-6, 1e-3, 1e-1])  # 3 humidity levels
    nu_grid = np.arange(100.0, 200.0, 1.0)

    kappa = sample_kappa_grid(
        background_gas="N2",
        absorbers={"CO2": 4e-4},
        h2o_vmr_grid=X_h2o,
        T_grid=T_grid,
        p_grid=p_grid,
        nu_grid=nu_grid,
        line_shape="lorentz",
        binning=True,
    )
    assert kappa.shape == (2, 2, 3, 100)
    # Wet column must have more opacity than dry one
    assert kappa[:, :, -1, :].mean() > kappa[:, :, 0, :].mean()
```

- [x] **Step 3: Run test to verify it fails**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_kappa_sampling.py -v`
Expected: FAIL with `ModuleNotFoundError: scripts.pf_table_builder.kappa_sampling`.

- [x] **Step 4: Implement the kappa sampler**

Create `scripts/pf_table_builder/kappa_sampling.py`:

```python
"""Sample linepyline line-by-line mass absorption coefficients on a (T, p[, X_H2O]) grid.

This is the slow step in table generation. Each (T, p[, X]) point requires one
call to linepyline.rtm.get_kappa_hitran (and MT_CKD continuum for H2O) per
absorber. With `binning=True` and ~0.5 cm^-1 wavenumber resolution, a 10x15
(T, p) grid over 10-30000 cm^-1 takes ~5-15 min per absorber on a modern laptop.
"""
from __future__ import annotations

import numpy as np
import xarray as xr


def sample_kappa_grid(
    *,
    background_gas,
    absorbers: dict,
    T_grid: np.ndarray,
    p_grid: np.ndarray,
    nu_grid: np.ndarray,
    h2o_vmr_grid: np.ndarray | None = None,
    line_shape: str = "lorentz",
    binning: bool = True,
    include_mtckd_continuum: bool = True,
    surface_gravity: float = 9.81,
):
    """Compute mass absorption coefficient kappa(T, p[, X_H2O], nu).

    Args:
        background_gas: "air", "N2", or None (passed to linepyline.rtm).
        absorbers: dict mapping HITRAN species (H2O, CO2, O3, CH4, NH3) to
            volume mixing ratio. Values are scalars — homogeneous columns only.
            For Hab1-style runs with variable H2O, pass H2O=0.0 here and supply
            h2o_vmr_grid; H2O kappa is computed separately and added with the
            grid VMR.
        T_grid: (nT,) temperatures in K (ascending).
        p_grid: (nP,) pressures in Pa (ascending; log-spaced recommended).
        nu_grid: (nNu,) wavenumber grid in cm^-1 (linepyline expects uniform).
        h2o_vmr_grid: optional (nX,) H2O VMRs. If given, output has an extra
            axis and H2O kappa is recomputed per VMR.
        line_shape, binning, include_mtckd_continuum: passed to linepyline.
        surface_gravity: passed to linepyline.rtm (does not affect kappa, only
            optical-depth integration, but linepyline requires it at __init__).

    Returns:
        kappa array, shape (nT, nP, nNu) or (nT, nP, nX, nNu) if h2o_vmr_grid
        is given. Units: m^2/kg of *total moist atmosphere* (consistent with
        how climt's correlated_k.compute_ck_optical_depth multiplies it by
        column air mass).
    """
    import linepyline

    rtm = linepyline.rtm(background_gas=background_gas,
                         surface_gravity=surface_gravity)

    nu_min, nu_max = float(nu_grid[0]), float(nu_grid[-1])
    dnu = float(np.diff(nu_grid).mean())
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)

    # Pre-compute kappa(T, p, nu) for each non-H2O absorber.
    # linepyline.get_kappa_hitran accepts 1D arrays of (T, p) and returns
    # kappa(p, nu) for a single T. We call it once per T.
    non_h2o = {s: v for s, v in absorbers.items() if s != "H2O"}
    base_kappa = np.zeros((nT, nP, nNu))
    for iT, Tval in enumerate(T_grid):
        T_arr = np.full(nP, Tval)
        for species, vmr in non_h2o.items():
            p_self = p_grid * vmr  # partial pressure (Pa)
            k_da = rtm.get_kappa_hitran(
                species, nu_min, nu_max, dnu,
                p_grid, T_arr, p_self=p_self,
                line_shape=line_shape, binning=binning,
                remove_plinth=False,
            )
            # rtm output is xarray (p, nu); interp wavenumber onto target grid.
            k_arr = k_da.transpose("p", "nu").interp(nu=nu_grid).values
            # Weight by mass mixing ratio of this species relative to total moist gas
            # For non-trace species we use molecular masses; linepyline's internal
            # mass weighting is done in get_optical_depth, so here we multiply by
            # the dimensionless mass fraction.
            mass_frac = _mass_fraction(species, vmr, background_gas, non_h2o)
            base_kappa[iT] += k_arr * mass_frac

    h2o_vmr_used = absorbers.get("H2O", 0.0)
    if h2o_vmr_grid is None:
        if h2o_vmr_used > 0.0:
            base_kappa += _h2o_kappa(
                rtm, nu_grid, p_grid, T_grid, h2o_vmr_used,
                line_shape, binning, include_mtckd_continuum,
                background_gas, non_h2o,
            )
        return base_kappa

    nX = len(h2o_vmr_grid)
    out = np.zeros((nT, nP, nX, nNu))
    for iX, X in enumerate(h2o_vmr_grid):
        h2o_contrib = _h2o_kappa(
            rtm, nu_grid, p_grid, T_grid, float(X),
            line_shape, binning, include_mtckd_continuum,
            background_gas, non_h2o,
        )
        out[:, :, iX, :] = base_kappa + h2o_contrib
    return out


def _h2o_kappa(rtm, nu_grid, p_grid, T_grid, X_h2o,
               line_shape, binning, include_mtckd_continuum,
               background_gas, other_absorbers):
    """H2O kappa (lines + optional MT_CKD continuum), shape (nT, nP, nNu)."""
    nu_min, nu_max = float(nu_grid[0]), float(nu_grid[-1])
    dnu = float(np.diff(nu_grid).mean())
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)
    out = np.zeros((nT, nP, nNu))
    for iT, Tval in enumerate(T_grid):
        T_arr = np.full(nP, Tval)
        p_self = p_grid * X_h2o
        k_lines = rtm.get_kappa_hitran(
            "H2O", nu_min, nu_max, dnu,
            p_grid, T_arr, p_self=p_self,
            line_shape=line_shape, binning=binning,
            remove_plinth=include_mtckd_continuum,
        )
        k = k_lines.transpose("p", "nu").interp(nu=nu_grid).values
        if include_mtckd_continuum:
            k_cont = rtm.get_kappa_mtckd(nu_min, nu_max, dnu,
                                         p_grid, T_arr, p_self)
            k = k + k_cont.transpose("p", "nu").interp(nu=nu_grid).values
        mass_frac = _mass_fraction(
            "H2O", X_h2o, background_gas,
            {**other_absorbers, "H2O": X_h2o},
        )
        out[iT] = k * mass_frac
    return out


# HITRAN molar masses (g/mol). Background gas masses follow phys.gases entries.
_MOLAR_MASS = {"H2O": 18.015, "CO2": 44.01, "O3": 47.998, "CH4": 16.04,
               "NH3": 17.031, "air": 28.97, "N2": 28.014}


def _mass_fraction(species, vmr, background_gas, all_absorbers):
    """Mass fraction of `species` in a moist gas mixture."""
    mean_mw = 0.0
    f_tot = 0.0
    for s, v in all_absorbers.items():
        mean_mw += v * _MOLAR_MASS[s]
        f_tot += v
    if background_gas is not None:
        mean_mw += (1.0 - f_tot) * _MOLAR_MASS[background_gas]
    return (vmr * _MOLAR_MASS[species]) / mean_mw
```

- [x] **Step 5: Run tests to verify they pass**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_kappa_sampling.py -v`
Expected: ALL PASS (first invocation downloads/loads HITRAN data into linepyline; subsequent runs are fast.)

- [x] **Step 6: Commit**

```bash
git add scripts/pf_table_builder/__init__.py \
        scripts/pf_table_builder/kappa_sampling.py \
        tests/pf_table_builder/
git commit -m "feat(pf-tables): add linepyline kappa sampler over (T, p[, X_H2O]) grid"
```

---

### Task 2: k-distribution builder — band-bin → sort → Gauss–Legendre

**Files:**
- Create: `scripts/pf_table_builder/k_distribution.py`
- Create: `tests/pf_table_builder/test_k_distribution.py`

Given `kappa[..., nu]` and band edges, build `k_coeffs[..., band, gpt]` using the classical k-distribution: within each band, sort κ values and integrate against the cumulative distribution function evaluated at Gauss–Legendre nodes on [0, 1].

- [x] **Step 1: Write failing test**

Create `tests/pf_table_builder/test_k_distribution.py`:

```python
import numpy as np
import pytest


def test_kappa_to_k_coeffs_uniform_recovers_value():
    """A uniform kappa(nu) in a band should give k_coeffs = kappa at every g-point."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.k_distribution import kappa_to_k_coeffs

    nu_grid = np.linspace(10.0, 3250.0, 1000)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    # Constant kappa = 1e-3 everywhere
    kappa = 1e-3 * np.ones((1, 1, len(nu_grid)))  # (nT=1, nP=1, nNu)

    k_coeffs, weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt=2)
    # Shape: (nT, nP, nband, ngpt)
    assert k_coeffs.shape == (1, 1, 4, 2)
    np.testing.assert_allclose(k_coeffs, 1e-3, rtol=1e-10)
    # Gauss-Legendre weights on [0,1] sum to 1
    np.testing.assert_allclose(weights.sum(axis=-1), 1.0, rtol=1e-10)


def test_kappa_to_k_coeffs_two_peaks_orders_correctly():
    """Within a band, k_coeffs at successive g-points are monotone increasing."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.k_distribution import kappa_to_k_coeffs

    nu_grid = np.linspace(10.0, 3250.0, 5000)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    # Pretend lines: kappa large at every 50 cm^-1, small elsewhere
    kappa = 1e-5 * np.ones_like(nu_grid)
    kappa[::50] = 1.0
    kappa = kappa.reshape(1, 1, -1)

    k_coeffs, _ = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt=4)
    # In every band, k must be non-decreasing in g
    for ib in range(4):
        for ig in range(3):
            assert k_coeffs[0, 0, ib, ig] <= k_coeffs[0, 0, ib, ig + 1] + 1e-12


def test_kappa_to_k_coeffs_with_extra_axes():
    """Builder is shape-agnostic in the leading kappa axes."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.k_distribution import kappa_to_k_coeffs

    nu_grid = np.linspace(10.0, 3250.0, 200)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    kappa = np.random.RandomState(0).uniform(1e-6, 1e-2,
                                              size=(3, 4, 5, len(nu_grid)))
    k_coeffs, weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt=2)
    assert k_coeffs.shape == (3, 4, 5, 4, 2)
    assert weights.shape == (4, 2)
```

- [x] **Step 2: Run — expect FAIL**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_k_distribution.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [x] **Step 3: Implement k-distribution**

Create `scripts/pf_table_builder/k_distribution.py`:

```python
"""Band-bin κ(ν) and quadrature to a small number of g-points."""
from __future__ import annotations

import numpy as np


def kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt):
    """Convert line-by-line κ to per-band k-distribution coefficients.

    Args:
        kappa: array shaped (..., nNu), arbitrary leading axes (typically
            (nT, nP) or (nT, nP, nX)).
        nu_grid: (nNu,) wavenumber grid (cm^-1), monotone increasing.
        band_edges: (nband+1,) band edges in cm^-1.
        ngpt: number of Gauss-Legendre nodes per band.

    Returns:
        k_coeffs: (..., nband, ngpt) k-coefficients evaluated at GL nodes.
        weights: (nband, ngpt) GL weights on [0, 1].
    """
    leading_shape = kappa.shape[:-1]
    nNu = kappa.shape[-1]
    assert len(nu_grid) == nNu, "kappa and nu_grid disagree on length"

    nband = len(band_edges) - 1
    # Gauss-Legendre on [-1, 1] -> [0, 1]
    xi_raw, wi_raw = np.polynomial.legendre.leggauss(ngpt)
    g_nodes = 0.5 * (xi_raw + 1.0)
    g_weights = 0.5 * wi_raw

    k_coeffs = np.zeros(leading_shape + (nband, ngpt))

    flat = kappa.reshape(-1, nNu)  # (nLead, nNu)

    for ib in range(nband):
        lo, hi = band_edges[ib], band_edges[ib + 1]
        mask = (nu_grid >= lo) & (nu_grid < hi)
        if not mask.any():
            raise ValueError(f"No nu_grid points in band [{lo}, {hi})")
        # Per-column: sort kappa in band, build CDF over wavenumber width,
        # then evaluate at GL g-nodes.
        for i in range(flat.shape[0]):
            k_band = np.sort(flat[i, mask])
            # Cumulative g axis: uniform in ν since the linepyline grid is uniform
            n_in = len(k_band)
            g_axis = (np.arange(n_in) + 0.5) / n_in
            # Interpolate k(g) onto the GL nodes
            k_coeffs.reshape(-1, nband, ngpt)[i, ib, :] = np.interp(
                g_nodes, g_axis, k_band,
            )

    weights = np.broadcast_to(g_weights, (nband, ngpt)).copy()
    return k_coeffs, weights
```

- [x] **Step 4: Run — expect PASS**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_k_distribution.py -v`
Expected: ALL PASS.

- [x] **Step 5: Commit**

```bash
git add scripts/pf_table_builder/k_distribution.py \
        tests/pf_table_builder/test_k_distribution.py
git commit -m "feat(pf-tables): k-distribution builder with Gauss-Legendre quadrature"
```

---

### Task 3: Planck-fraction utility (LW source distribution per g-point)

**Files:**
- Create: `scripts/pf_table_builder/planck_fraction.py`
- Create: `tests/pf_table_builder/test_planck_fraction.py`

climt's LW correlated-k loader expects `planck_fraction(band, gpt, T)` so that the Planck source per (band, g) is `B_band(T) × planck_fraction`. Because we built the k-distribution by sorting κ in wavenumber without preserving which ν went to which g, the assumption is that within a band the Planck source is distributed *uniformly across g* — i.e., `planck_fraction = 1/ngpt` for every (band, gpt, T). This matches what the existing Earth/Mars/Venus tables ship.

This task gives that a name and a unit test so the writer (Task 6) reads from one canonical helper.

- [x] **Step 1: Failing test**

Create `tests/pf_table_builder/test_planck_fraction.py`:

```python
import numpy as np


def test_uniform_planck_fraction_sums_to_one():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.planck_fraction import build_uniform_planck_fraction

    nband, ngpt, nT = 4, 2, 5
    pf = build_uniform_planck_fraction(nband, ngpt, nT)
    assert pf.shape == (nband, ngpt, nT)
    np.testing.assert_allclose(pf.sum(axis=1), 1.0, rtol=1e-12)
    np.testing.assert_allclose(pf, 1.0 / ngpt)
```

- [x] **Step 2: Run — expect FAIL**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_planck_fraction.py -v`
Expected: FAIL.

- [x] **Step 3: Implement**

Create `scripts/pf_table_builder/planck_fraction.py`:

```python
"""Planck-source distribution across g-points within each LW band."""
import numpy as np


def build_uniform_planck_fraction(nband, ngpt, nT):
    """Uniform 1/ngpt fraction.

    Matches the convention used by climt's existing low-res tables. A
    higher-fidelity scheme would weight by the sub-band wavenumber slice
    that each g-point represents, but with κ sorted by value (not by ν)
    that mapping has been thrown away. Uniform is conservative and exact
    in the band-averaged sense.
    """
    return np.full((nband, ngpt, nT), 1.0 / ngpt)
```

- [x] **Step 4: Run + commit**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_planck_fraction.py -v`
Expected: PASS.

```bash
git add scripts/pf_table_builder/planck_fraction.py \
        tests/pf_table_builder/test_planck_fraction.py
git commit -m "feat(pf-tables): uniform Planck-fraction helper for LW writer"
```

---

### Task 4: SW solar-source + Rayleigh per band

**Files:**
- Create: `scripts/pf_table_builder/solar_source.py`
- Create: `tests/pf_table_builder/test_solar_source.py`

Compute `solar_source_per_gpoint(band, gpt)` by integrating a stellar spectrum across each band and splitting equally across the band's g-points (same uniform assumption as Planck fraction). Also compute a single-value-per-band Rayleigh coefficient using a standard λ⁻⁴ formula and the bulk gas refractive index (parameterized by mean molar mass and a refractivity coefficient).

- [x] **Step 1: Failing test**

Create `tests/pf_table_builder/test_solar_source.py`:

```python
import numpy as np


def test_solar_source_partitions_total_irradiance():
    """Sum of solar_source across (band, gpt) equals total stellar flux integral."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.solar_source import (
        build_solar_source_per_gpoint, build_rayleigh_per_band,
    )

    # Simple flat spectrum: 1 W/m²/cm⁻¹ from 1000 to 30000 cm^-1
    wn = np.linspace(1000.0, 30000.0, 5000)
    irr = np.ones_like(wn)
    spectrum = {"wavenumber": wn, "irradiance": irr}
    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])

    src = build_solar_source_per_gpoint(spectrum, band_edges, ngpt=2)
    assert src.shape == (3, 2)
    # Total = 30000 - 3250 (flat unit spectrum)
    np.testing.assert_allclose(src.sum(), 30000.0 - 3250.0, rtol=1e-3)


def test_rayleigh_per_band_decreases_with_wavelength():
    """Rayleigh coefficient should drop as wavenumber drops (longer λ)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.solar_source import build_rayleigh_per_band

    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])
    r = build_rayleigh_per_band(band_edges, mean_molar_mass_g=28.0,
                                refractivity_298K=2.97e-4)
    assert r.shape == (3,)
    assert r[0] < r[1] < r[2]
    assert (r >= 0).all()
```

- [x] **Step 2: Run — expect FAIL**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_solar_source.py -v`
Expected: FAIL.

- [x] **Step 3: Implement**

Create `scripts/pf_table_builder/solar_source.py`:

```python
"""SW band-integrated solar source and Rayleigh cross-section."""
from __future__ import annotations

import numpy as np


def build_solar_source_per_gpoint(spectrum, band_edges, ngpt):
    """Distribute band-integrated stellar flux uniformly across g-points.

    Args:
        spectrum: dict with "wavenumber" (cm^-1) and "irradiance" (W/m^2/cm^-1).
        band_edges: (nband+1,) cm^-1.
        ngpt: per-band g-point count.

    Returns:
        (nband, ngpt) array, W/m^2 per (band, g-point).
    """
    wn = np.asarray(spectrum["wavenumber"])
    irr = np.asarray(spectrum["irradiance"])
    nband = len(band_edges) - 1
    out = np.zeros((nband, ngpt))
    for ib in range(nband):
        lo, hi = band_edges[ib], band_edges[ib + 1]
        mask = (wn > lo) & (wn < hi)
        wn_b = np.concatenate(([lo], wn[mask], [hi]))
        irr_b = np.concatenate(([np.interp(lo, wn, irr)],
                                irr[mask],
                                [np.interp(hi, wn, irr)]))
        band_flux = np.trapezoid(irr_b, wn_b)
        out[ib, :] = band_flux / ngpt
    return out


def build_rayleigh_per_band(band_edges, mean_molar_mass_g,
                            refractivity_298K=2.97e-4):
    """Band-mean Rayleigh mass scattering coefficient (m²/kg).

    Uses Bodhaine et al. (1999) form. `refractivity_298K` is (n-1) at STP
    for the bulk gas; defaults to dry-air. For N₂ use ~2.97e-4, for CO₂ ~4.5e-4,
    for H₂ ~1.4e-4.

    The wavenumber dependence is ν⁴ (i.e. λ⁻⁴). We evaluate at the band
    geometric mean for the coefficient.
    """
    # Bodhaine: σ(λ) = (24 π³ / N²λ⁴) × ((n²-1)/(n²+2))²
    # mass cross-section = σ × N_A / M
    band_mid_wn = np.sqrt(band_edges[:-1] * band_edges[1:])  # cm^-1
    lam_m = 1.0 / (band_mid_wn * 100.0)  # cm^-1 -> m
    N_loschmidt = 2.547e25  # molecules/m^3 at STP
    f = (refractivity_298K + 1.0) ** 2
    sigma = (24.0 * np.pi**3) / (N_loschmidt**2 * lam_m**4) * \
            ((f - 1) / (f + 2)) ** 2
    # Convert per-molecule cross-section (m²) to per-mass (m²/kg)
    M = mean_molar_mass_g * 1e-3  # kg/mol
    N_A = 6.022e23
    return sigma * N_A / M
```

- [x] **Step 4: Run + commit**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_solar_source.py -v`
Expected: PASS.

```bash
git add scripts/pf_table_builder/solar_source.py \
        tests/pf_table_builder/test_solar_source.py
git commit -m "feat(pf-tables): SW solar source + band Rayleigh utilities"
```

---

### Task 5: HITRAN CIA loader (Titan N₂–N₂, N₂–CH₄)

**Files:**
- Create: `scripts/pf_table_builder/cia.py`
- Create: `tests/pf_table_builder/test_cia.py`
- Create: `climt/_data/picket_fence/cia/README.md`

HITRAN distributes CIA cross-sections (units cm⁻¹ amagat⁻²) as ASCII files split into temperature blocks. Each block header is `<species1> <species2> ν_min ν_max nPoints T ...`. We need to read them, interpolate κ_CIA(ν, T) at each (T, p), and convert to mass absorption coefficient (m²/kg).

For Titan the relevant pairs are **N₂–N₂** (dominant rotational continuum 0–600 cm⁻¹) and **N₂–CH₄** (broad bands 0–1400 cm⁻¹). Download from <https://hitran.org/cia/>:
- `N2-N2_2018.cia`
- `N2-CH4_2018.cia`

(File names current as of HITRAN 2024 release; substitute the latest if newer.)

- [x] **Step 1: README and placeholder for downloaded files**

Create `climt/_data/picket_fence/cia/README.md`:

```markdown
# HITRAN CIA flat files (offline build only)

These ASCII files are consumed by `scripts/pf_table_builder/cia.py` when
building the Titan opacity table. They are not loaded at climt runtime —
the CIA contribution is baked into the resulting `.nc` k-table.

Source: <https://hitran.org/cia/> (Karman et al., 2019; CC-BY).

Files expected (download once):
- `N2-N2_2018.cia`
- `N2-CH4_2018.cia`

If your HITRAN download supplies only `CH4-N2_…`, that is the symmetric
pair — the loader accepts either ordering.
```

- [x] **Step 2: Failing test (uses a tiny synthetic CIA file)**

Create `tests/pf_table_builder/test_cia.py`:

```python
import numpy as np
import textwrap


def _write_fake_cia(tmp_path, name="N2-N2_fake.cia"):
    """Write a 2-block, 3-point synthetic HITRAN CIA file."""
    text = textwrap.dedent("""\
        N2-N2          10.0     30.0    3  200.0  1.000E-50    Karman2019
            10.0      1.0E-46
            20.0      2.0E-46
            30.0      4.0E-46
        N2-N2          10.0     30.0    3  300.0  1.000E-50    Karman2019
            10.0      2.0E-46
            20.0      4.0E-46
            30.0      8.0E-46
        """)
    path = tmp_path / name
    path.write_text(text)
    return path


def test_load_cia_blocks(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.cia import load_cia_blocks

    p = _write_fake_cia(tmp_path)
    blocks = load_cia_blocks(str(p))
    assert sorted(blocks.keys()) == [200.0, 300.0]
    nu, k = blocks[200.0]
    np.testing.assert_allclose(nu, [10.0, 20.0, 30.0])
    np.testing.assert_allclose(k, [1e-46, 2e-46, 4e-46])


def test_cia_kappa_on_grid(tmp_path):
    """cia_kappa_on_grid returns kappa(T,p,nu) in m^2/kg for the pair."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.cia import cia_kappa_on_grid

    p = _write_fake_cia(tmp_path)
    T_grid = np.array([200.0, 250.0, 300.0])
    p_grid = np.array([1e4, 1e5])
    nu_grid = np.array([15.0, 25.0])
    kappa = cia_kappa_on_grid(
        str(p), pair=("N2", "N2"),
        vmr_a=1.0, vmr_b=1.0,
        background_gas="N2",
        absorbers={"N2": 1.0},  # used for mass-fraction normalisation
        T_grid=T_grid, p_grid=p_grid, nu_grid=nu_grid,
    )
    assert kappa.shape == (3, 2, 2)
    # CIA grows like p² at fixed T (amagat² scaling), so high-p > low-p
    assert (kappa[:, 1, :] > kappa[:, 0, :]).all()
    # And grows with T in the synthetic file
    assert kappa[2, 0, 0] > kappa[0, 0, 0]
```

- [x] **Step 3: Run — expect FAIL**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_cia.py -v`
Expected: FAIL.

- [x] **Step 4: Implement the CIA loader**

Create `scripts/pf_table_builder/cia.py`:

```python
"""Load HITRAN CIA flat files and evaluate κ_CIA(T, p, ν) on a (T, p) grid.

HITRAN CIA distributes binary collision cross-sections per temperature
block in cm⁻¹ amagat⁻². We convert to a per-mass absorption coefficient
(m²/kg of total moist gas) so the result can be added directly to the
gas-line kappa array from kappa_sampling.sample_kappa_grid.
"""
from __future__ import annotations

import numpy as np


_AMAGAT_NUMBER_DENSITY = 2.6868e25  # molecules / m³ at STP
_N_A = 6.022e23


def load_cia_blocks(path):
    """Parse a HITRAN CIA file. Returns dict: T (K) -> (nu_array, sigma_array).

    sigma units: cm⁻¹ amagat⁻² (the HITRAN file convention).
    """
    blocks = {}
    with open(path, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        header = lines[i].split()
        if len(header) < 5:
            i += 1
            continue
        nPoints = int(header[3])
        T = float(header[4])
        data = np.loadtxt(lines[i + 1: i + 1 + nPoints])
        blocks[T] = (data[:, 0], data[:, 1])
        i += 1 + nPoints
    return blocks


def cia_kappa_on_grid(
    path, *, pair, vmr_a, vmr_b,
    background_gas, absorbers,
    T_grid, p_grid, nu_grid,
):
    """Evaluate κ_CIA on (T, p, ν) and convert to m²/kg of moist gas.

    Args:
        path: HITRAN .cia flat file.
        pair: tuple (species_a, species_b); used for the mass-fraction
            normalisation only (the file already encodes which pair).
        vmr_a, vmr_b: volume mixing ratios of the two colliders.
        background_gas, absorbers: same conventions as kappa_sampling, used
            only to compute the mean molar mass for unit conversion.
        T_grid, p_grid, nu_grid: target grids.

    Returns:
        kappa: (nT, nP, nNu), m²/kg.
    """
    blocks = load_cia_blocks(path)
    T_file = np.array(sorted(blocks.keys()))
    nT_f = len(T_file)
    # Resample each block onto nu_grid (zero outside file's spectral range)
    sigma_grid = np.zeros((nT_f, len(nu_grid)))
    for iT, T in enumerate(T_file):
        nu_f, sig_f = blocks[T]
        mask = (nu_grid >= nu_f[0]) & (nu_grid <= nu_f[-1])
        sigma_grid[iT, mask] = np.interp(nu_grid[mask], nu_f, sig_f)

    # Interpolate in T (clamp outside file range)
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)
    sigma_TpNu = np.zeros((nT, nNu))
    for inu in range(nNu):
        sigma_TpNu[:, inu] = np.interp(T_grid, T_file, sigma_grid[:, inu],
                                       left=sigma_grid[0, inu],
                                       right=sigma_grid[-1, inu])

    # CIA absorption coefficient (m^-1):
    #   α(ν) = σ(ν, T) × n_a × n_b
    # where n_x is number density of colliders in amagat. amagat = n/n_STP.
    # For an ideal gas:  n_amagat = (p/p_STP) × (T_STP/T)
    p_STP, T_STP = 1.01325e5, 273.15
    n_a = (p_grid * vmr_a / p_STP)[:, None] * (T_STP / T_grid)[None, :]  # (nP, nT)
    n_b = (p_grid * vmr_b / p_STP)[:, None] * (T_STP / T_grid)[None, :]

    # Convert σ from cm⁻¹ amagat⁻² to m⁻¹ amagat⁻²: ×100
    sigma_si = sigma_TpNu * 100.0  # (nT, nNu)

    # κ_volume [m⁻¹] = σ × n_a × n_b  (per amagat²)
    # broadcast (nT, 1, nNu) × (nT, nP) × (nT, nP) → (nT, nP, nNu)
    n_a = n_a.T  # (nT, nP)
    n_b = n_b.T  # (nT, nP)
    kappa_vol = sigma_si[:, None, :] * (n_a * n_b)[:, :, None]

    # Convert to mass absorption coefficient: divide by mass density of moist gas
    # ρ(T, p) = p × M_mean / (R_universal × T)
    R = 8.31446
    _MM = {"H2O": 18.015, "CO2": 44.01, "O3": 47.998, "CH4": 16.04,
           "NH3": 17.031, "N2": 28.014, "air": 28.97}
    f_tot = sum(absorbers.values())
    M_mean_g = sum(v * _MM[s] for s, v in absorbers.items())
    if background_gas is not None:
        M_mean_g += (1.0 - f_tot) * _MM[background_gas]
    M_mean = M_mean_g * 1e-3  # kg/mol
    rho = (p_grid[None, :] * M_mean) / (R * T_grid[:, None])  # (nT, nP) kg/m^3

    kappa_mass = kappa_vol / rho[:, :, None]
    return kappa_mass
```

- [x] **Step 5: Run + commit**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/pf_table_builder/test_cia.py -v`
Expected: ALL PASS.

```bash
git add scripts/pf_table_builder/cia.py \
        tests/pf_table_builder/test_cia.py \
        climt/_data/picket_fence/cia/README.md
git commit -m "feat(pf-tables): HITRAN CIA loader for Titan N2-N2 / N2-CH4"
```

---

### Task 6: netCDF writer matching climt's correlated-k schema

**Files:**
- Create: `scripts/pf_table_builder/netcdf_writer.py`
- Create: `tests/pf_table_builder/test_netcdf_writer.py`

The writer mirrors variable names + dimensions consumed by `climt._components.picket_fence.optics.correlated_k._load_netcdf_table`:
- LW: `k_coefficients`, `gpoint_weights`, `temperature_grid`, `pressure_grid_log`, optionally `h2o_vmr_grid`, `band_wavenumber_limits`, `planck_fraction`.
- SW: same plus `solar_source_per_gpoint`, `rayleigh_coefficient`. No `planck_fraction`.

- [x] **Step 1: Failing test**

Create `tests/pf_table_builder/test_netcdf_writer.py`:

```python
import numpy as np
import pytest


def test_write_lw_table_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.netcdf_writer import write_lw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP = 1, 4, 2, 3, 5
    k = np.random.RandomState(0).uniform(1e-7, 1e-2,
                                          size=(ngas, nband, ngpt, nT, nP))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    planck_fraction = np.full((nband, ngpt, nT), 1.0 / ngpt)

    path = tmp_path / "test_lw.nc"
    write_lw_table(
        str(path),
        k_coefficients=k,
        gpoint_weights=weights,
        T_grid=T_grid,
        log_p_grid=log_p_grid,
        band_edges=band_edges,
        planck_fraction=planck_fraction,
        gas_names=("effective",),
        source="linepyline:test",
    )
    loaded = load_k_table(str(path))
    np.testing.assert_allclose(loaded["k_coefficients"], k, rtol=1e-5)
    np.testing.assert_allclose(loaded["gpoint_weights"], weights, rtol=1e-5)
    np.testing.assert_allclose(loaded["temperature_grid"], T_grid, rtol=1e-5)
    np.testing.assert_allclose(loaded["pressure_grid_log"], log_p_grid, rtol=1e-5)


def test_write_lw_table_with_h2o_axis_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.netcdf_writer import write_lw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP, nX = 1, 4, 2, 3, 5, 4
    k = np.random.RandomState(1).uniform(1e-7, 1e-2,
                                          size=(ngas, nband, ngpt, nT, nP, nX))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
    h2o_vmr = np.array([1e-6, 1e-4, 1e-2, 1e-1])
    planck_fraction = np.full((nband, ngpt, nT), 1.0 / ngpt)

    path = tmp_path / "test_lw_hmoist.nc"
    write_lw_table(
        str(path),
        k_coefficients=k,
        gpoint_weights=weights,
        T_grid=T_grid,
        log_p_grid=log_p_grid,
        band_edges=band_edges,
        planck_fraction=planck_fraction,
        h2o_vmr_grid=h2o_vmr,
        gas_names=("effective",),
        source="linepyline:test_hmoist",
    )
    loaded = load_k_table(str(path))
    assert loaded["k_coefficients"].shape == (ngas, nband, ngpt, nT, nP, nX)
    np.testing.assert_allclose(loaded["h2o_vmr_grid"], h2o_vmr, rtol=1e-5)


def test_write_sw_table_roundtrip(tmp_path):
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.netcdf_writer import write_sw_table
    from climt._components.picket_fence.optics.correlated_k import load_k_table

    ngas, nband, ngpt, nT, nP = 1, 3, 2, 3, 5
    k = np.random.RandomState(2).uniform(0, 1e-3,
                                          size=(ngas, nband, ngpt, nT, nP))
    weights = np.full((nband, ngpt), 0.5)
    T_grid = np.linspace(150, 350, nT)
    log_p_grid = np.linspace(np.log(1.0), np.log(1e5), nP)
    band_edges = np.array([3250.0, 8000.0, 14000.0, 30000.0])
    solar = np.array([[100.0, 100.0], [200.0, 200.0], [300.0, 300.0]])
    rayleigh = np.array([1e-6, 5e-6, 2e-5])

    path = tmp_path / "test_sw.nc"
    write_sw_table(
        str(path),
        k_coefficients=k, gpoint_weights=weights,
        T_grid=T_grid, log_p_grid=log_p_grid, band_edges=band_edges,
        solar_source_per_gpoint=solar, rayleigh_coefficient=rayleigh,
        gas_names=("effective",), source="linepyline:test_sw",
    )
    loaded = load_k_table(str(path))
    np.testing.assert_allclose(loaded["solar_source_per_gpoint"], solar, rtol=1e-5)
    np.testing.assert_allclose(loaded["rayleigh_coefficient"], rayleigh, rtol=1e-5)
```

- [x] **Step 2: Run — expect FAIL**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/pf_table_builder/test_netcdf_writer.py -v`
Expected: FAIL.

- [x] **Step 3: Implement the writer**

Create `scripts/pf_table_builder/netcdf_writer.py`:

```python
"""Write picket-fence correlated-k netCDF tables in climt's schema."""
from __future__ import annotations

import os

import numpy as np
from scipy.io import netcdf_file


def _ensure_parent(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)


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
    gas_names=("effective",),
    overlap_method="additive",
    resolution="low",
    source="linepyline",
):
    """Write LW k-table.

    k_coefficients shape:
        (ngas, nband, ngpt, nT, nP) or (ngas, nband, ngpt, nT, nP, nX) if H2O axis.
    """
    _ensure_parent(out_path)
    has_x = (h2o_vmr_grid is not None)
    if has_x:
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

        pf = nc.createVariable("planck_fraction", "f4",
                               ("band", "gpoint", "temperature"))
        pf[:] = planck_fraction.astype("f4")

        if has_x:
            xv = nc.createVariable("h2o_vmr_grid", "f4", ("h2o_vmr",))
            xv[:] = np.asarray(h2o_vmr_grid, dtype="f4")

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = overlap_method
        nc.resolution = resolution
        nc.source = source


def write_sw_table(
    out_path, *,
    k_coefficients, gpoint_weights, T_grid, log_p_grid, band_edges,
    solar_source_per_gpoint, rayleigh_coefficient,
    h2o_vmr_grid=None,
    gas_names=("effective",),
    overlap_method="additive", resolution="low", source="linepyline",
):
    _ensure_parent(out_path)
    has_x = (h2o_vmr_grid is not None)
    if has_x:
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

        ss = nc.createVariable("solar_source_per_gpoint", "f4", ("band", "gpoint"))
        ss[:] = solar_source_per_gpoint.astype("f4")

        rc = nc.createVariable("rayleigh_coefficient", "f4", ("band",))
        rc[:] = np.asarray(rayleigh_coefficient, dtype="f4")

        if has_x:
            xv = nc.createVariable("h2o_vmr_grid", "f4", ("h2o_vmr",))
            xv[:] = np.asarray(h2o_vmr_grid, dtype="f4")

        nc.gas_names = ",".join(gas_names)
        nc.overlap_method = overlap_method
        nc.resolution = resolution
        nc.source = source
```

- [x] **Step 4: Run + commit**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/pf_table_builder/test_netcdf_writer.py -v`
Expected: PASS.

```bash
git add scripts/pf_table_builder/netcdf_writer.py \
        tests/pf_table_builder/test_netcdf_writer.py
git commit -m "feat(pf-tables): netCDF writer for picket-fence k-tables"
```

---

### Task 7: CLI driver — `generate_pf_tables_linepyline.py`

**Files:**
- Create: `scripts/generate_pf_tables_linepyline.py`

A single CLI that picks a named scenario and assembles all five helpers into one `.nc` file. Scenarios are hard-coded so reproducing a table is one shell command.

- [x] **Step 1: Write the driver**

Create `scripts/generate_pf_tables_linepyline.py`:

```python
#!/usr/bin/env python
"""Generate a picket-fence correlated-k table for a named planetary scenario.

Usage:
    conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \\
        --scenario trappist1e_hab2 --kind lw --output climt/_data/picket_fence/correlated_k/trappist1e_hab2_lw.nc

Scenarios:
    trappist1e_hab1  THAI Hab1: N2 bath + 400 ppm CO2 + variable H2O (axis)
    trappist1e_hab2  THAI Hab2: pure 1-bar CO2
    titan            Titan: 95% N2 + 5% CH4, N2-N2 and N2-CH4 CIA included

`--kind lw|sw` selects band edges and whether Planck or solar source is written.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from pf_table_builder.kappa_sampling import sample_kappa_grid
from pf_table_builder.k_distribution import kappa_to_k_coeffs
from pf_table_builder.planck_fraction import build_uniform_planck_fraction
from pf_table_builder.solar_source import (
    build_solar_source_per_gpoint, build_rayleigh_per_band,
)
from pf_table_builder.netcdf_writer import write_lw_table, write_sw_table
from pf_table_builder.cia import cia_kappa_on_grid


# Scenario configuration registry
SCENARIOS = {
    "trappist1e_hab1": dict(
        background_gas="N2",
        absorbers={"CO2": 4e-4},  # H2O handled via h2o_vmr_grid
        h2o_vmr_grid=np.array([1e-7, 1e-5, 1e-3, 1e-2, 1e-1]),
        T_grid=np.linspace(150.0, 350.0, 10),
        p_grid=np.logspace(0.0, 5.05, 15),     # 1 Pa to ~1.1e5 Pa
        stellar_spectrum="trappist1",
        mean_molar_mass_g=28.05,
        rayleigh_refractivity=2.97e-4,         # N2-dominated
        cia_files=[],
    ),
    "trappist1e_hab2": dict(
        background_gas=None,
        absorbers={"CO2": 1.0},
        h2o_vmr_grid=None,
        T_grid=np.linspace(150.0, 350.0, 10),
        p_grid=np.logspace(0.0, 5.05, 15),
        stellar_spectrum="trappist1",
        mean_molar_mass_g=44.01,
        rayleigh_refractivity=4.50e-4,         # CO2
        cia_files=[],
    ),
    "titan": dict(
        background_gas="N2",
        absorbers={"CH4": 0.05},
        h2o_vmr_grid=None,
        T_grid=np.linspace(70.0, 200.0, 10),
        p_grid=np.logspace(-1.0, 5.2, 15),
        stellar_spectrum="sun",
        mean_molar_mass_g=28.6,
        rayleigh_refractivity=2.97e-4,
        cia_files=[
            dict(path="climt/_data/picket_fence/cia/N2-N2_2018.cia",
                 pair=("N2", "N2"), vmr_a=0.95, vmr_b=0.95),
            dict(path="climt/_data/picket_fence/cia/N2-CH4_2018.cia",
                 pair=("N2", "CH4"), vmr_a=0.95, vmr_b=0.05),
        ],
    ),
}

LW_BAND_EDGES = np.array([10.0, 500.0, 1250.0, 2500.0, 3250.0])
SW_BAND_EDGES = np.array([3250.0, 8000.0, 14000.0, 30000.0])


def build_table(scenario_name, kind, output_path, ngpt=2, dnu=0.5,
                line_shape="lorentz", binning=True):
    cfg = SCENARIOS[scenario_name]
    band_edges = LW_BAND_EDGES if kind == "lw" else SW_BAND_EDGES
    nu_grid = np.arange(band_edges[0], band_edges[-1] + dnu / 2, dnu)

    print(f"[{scenario_name}/{kind}] sampling kappa on "
          f"(nT={len(cfg['T_grid'])}, nP={len(cfg['p_grid'])}, nNu={len(nu_grid)})")
    kappa = sample_kappa_grid(
        background_gas=cfg["background_gas"],
        absorbers=cfg["absorbers"],
        h2o_vmr_grid=cfg["h2o_vmr_grid"],
        T_grid=cfg["T_grid"],
        p_grid=cfg["p_grid"],
        nu_grid=nu_grid,
        line_shape=line_shape,
        binning=binning,
    )

    # Add CIA contribution if any (CIA is independent of X_H2O)
    if cfg["cia_files"]:
        print(f"[{scenario_name}/{kind}] adding {len(cfg['cia_files'])} CIA pair(s)")
        for entry in cfg["cia_files"]:
            kappa_cia = cia_kappa_on_grid(
                entry["path"], pair=entry["pair"],
                vmr_a=entry["vmr_a"], vmr_b=entry["vmr_b"],
                background_gas=cfg["background_gas"],
                absorbers={**cfg["absorbers"],
                           cfg["background_gas"]: 1.0 - sum(cfg["absorbers"].values())},
                T_grid=cfg["T_grid"],
                p_grid=cfg["p_grid"],
                nu_grid=nu_grid,
            )
            # kappa_cia is (nT, nP, nNu); broadcast over X axis if present
            if kappa.ndim == 4:
                kappa = kappa + kappa_cia[:, :, None, :]
            else:
                kappa = kappa + kappa_cia

    print(f"[{scenario_name}/{kind}] building k-distribution (ngpt={ngpt})")
    k_coeffs, gpt_weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt)
    # Add the leading "ngas=1" axis expected by the climt loader.
    k_coeffs = k_coeffs[np.newaxis]

    # Rearrange axes to (ngas, nband, ngpt, nT, nP[, nX])
    if kappa.ndim == 4:
        # current k_coeffs shape: (1, nT, nP, nX, nband, ngpt)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4, 5), (3, 4, 5, 1, 2))
    else:
        # current k_coeffs shape: (1, nT, nP, nband, ngpt)
        k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4), (3, 4, 1, 2))

    log_p_grid = np.log(cfg["p_grid"])

    if kind == "lw":
        planck_fraction = build_uniform_planck_fraction(
            len(band_edges) - 1, ngpt, len(cfg["T_grid"]),
        )
        write_lw_table(
            output_path,
            k_coefficients=k_coeffs,
            gpoint_weights=gpt_weights,
            T_grid=cfg["T_grid"],
            log_p_grid=log_p_grid,
            band_edges=band_edges,
            planck_fraction=planck_fraction,
            h2o_vmr_grid=cfg["h2o_vmr_grid"],
            source=f"linepyline:{scenario_name}",
        )
    else:
        # Load stellar spectrum and build solar source / Rayleigh
        from climt._components.picket_fence.optics.stellar import (
            load_stellar_spectrum,
        )
        spectrum = load_stellar_spectrum(cfg["stellar_spectrum"])
        solar = build_solar_source_per_gpoint(spectrum, band_edges, ngpt)
        rayleigh = build_rayleigh_per_band(
            band_edges,
            mean_molar_mass_g=cfg["mean_molar_mass_g"],
            refractivity_298K=cfg["rayleigh_refractivity"],
        )
        write_sw_table(
            output_path,
            k_coefficients=k_coeffs,
            gpoint_weights=gpt_weights,
            T_grid=cfg["T_grid"],
            log_p_grid=log_p_grid,
            band_edges=band_edges,
            solar_source_per_gpoint=solar,
            rayleigh_coefficient=rayleigh,
            h2o_vmr_grid=cfg["h2o_vmr_grid"],
            source=f"linepyline:{scenario_name}",
        )
    print(f"[{scenario_name}/{kind}] wrote {output_path}  "
          f"k_coefficients.shape={k_coeffs.shape}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=sorted(SCENARIOS), required=True)
    ap.add_argument("--kind", choices=("lw", "sw"), required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--ngpt", type=int, default=2)
    ap.add_argument("--dnu", type=float, default=0.5,
                    help="wavenumber resolution in cm^-1 (default 0.5)")
    ap.add_argument("--line-shape", default="lorentz",
                    choices=("lorentz", "voigt", "pseudovoigt"))
    args = ap.parse_args()
    build_table(args.scenario, args.kind, args.output,
                ngpt=args.ngpt, dnu=args.dnu, line_shape=args.line_shape)


if __name__ == "__main__":
    main()
```

- [x] **Step 2: Smoke-run the driver in synthetic mode**

To avoid a long build during plan execution, first run it for `trappist1e_hab2 / lw` with a very coarse grid by editing `SCENARIOS` temporarily — or just trust the unit tests and proceed to Task 8.

- [x] **Step 3: Commit**

```bash
git add scripts/generate_pf_tables_linepyline.py
git commit -m "feat(pf-tables): CLI driver for Titan/TRAPPIST-1e scenarios"
```

---

### Task 8: Build TRAPPIST-1e Hab2 tables (simplest scenario first)

**Files:**
- Create: `climt/_data/picket_fence/correlated_k/trappist1e_hab2_lw.nc`
- Create: `climt/_data/picket_fence/correlated_k/trappist1e_hab2_sw.nc`

Pure-CO2 atmosphere, no H2O VMR axis. This is the cleanest validation of the full pipeline because the only opacity is CO2 lines.

- [x] **Step 1: Build LW table**

Run:
```bash
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \
    --scenario trappist1e_hab2 --kind lw \
    --output climt/_data/picket_fence/correlated_k/trappist1e_hab2_lw.nc \
    --ngpt 2 --dnu 0.5
```
Expected output: `[trappist1e_hab2/lw] wrote climt/_data/picket_fence/correlated_k/trappist1e_hab2_lw.nc  k_coefficients.shape=(1, 4, 2, 10, 15)`.
Expected runtime: 5–20 minutes (CO2 has many lines).

- [x] **Step 2: Build SW table**

Run:
```bash
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \
    --scenario trappist1e_hab2 --kind sw \
    --output climt/_data/picket_fence/correlated_k/trappist1e_hab2_sw.nc \
    --ngpt 2 --dnu 1.0
```
(SW range 3250–30000 cm⁻¹ is wide; use coarser dnu.)

- [x] **Step 3: Sanity-check from the climt env**

Run:
```bash
cd /Users/joymonteiro/github/climt && conda run -n climt python -c "
from climt._components.picket_fence.optics.correlated_k import load_k_table
t = load_k_table('trappist1e_hab2_lw')
print('LW shape:', t['k_coefficients'].shape, 'min', t['k_coefficients'].min(), 'max', t['k_coefficients'].max())
t = load_k_table('trappist1e_hab2_sw')
print('SW shape:', t['k_coefficients'].shape)
print('SW solar:', t['solar_source_per_gpoint'].sum(), 'W/m^2')
"
```
Expected:
- LW shape `(1, 4, 2, 10, 15)` with `min ≥ 0` and `max` not absurd (< 100 m²/kg).
- SW solar total of order ~900 W/m² (TRAPPIST-1 integrated flux at 1 AU equivalent).

- [x] **Step 4: Commit (binary files)**

```bash
git add climt/_data/picket_fence/correlated_k/trappist1e_hab2_lw.nc \
        climt/_data/picket_fence/correlated_k/trappist1e_hab2_sw.nc
git commit -m "feat(pf-tables): ship TRAPPIST-1e Hab2 (pure CO2) LW+SW tables"
```

---

### Task 9: Build TRAPPIST-1e Hab1 tables (with H2O VMR axis)

**Files:**
- Create: `climt/_data/picket_fence/correlated_k/trappist1e_hab1_lw.nc`
- Create: `climt/_data/picket_fence/correlated_k/trappist1e_hab1_sw.nc`

Same as Hab2 but with N2 bath, 400 ppm CO2, and the H2O-VMR axis (5 humidity bins).

- [x] **Step 1: Build LW**

Run:
```bash
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \
    --scenario trappist1e_hab1 --kind lw \
    --output climt/_data/picket_fence/correlated_k/trappist1e_hab1_lw.nc \
    --ngpt 2 --dnu 0.5
```
Expected shape: `(1, 4, 2, 10, 15, 5)`.
Expected runtime: 30–60 min (H2O is the slow molecule; nX=5 multiplies that).

- [x] **Step 2: Build SW**

```bash
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \
    --scenario trappist1e_hab1 --kind sw \
    --output climt/_data/picket_fence/correlated_k/trappist1e_hab1_sw.nc \
    --ngpt 2 --dnu 1.0
```

- [x] **Step 3: Sanity-check via climt loader (must trigger trilinear interp path)**

```bash
cd /Users/joymonteiro/github/climt && conda run -n climt python -c "
import numpy as np
from climt._components.picket_fence.optics.correlated_k import load_k_table, interpolate_k
t = load_k_table('trappist1e_hab1_lw')
print('shape', t['k_coefficients'].shape)
print('h2o_vmr_grid', t['h2o_vmr_grid'])
# trilinear lookup at a single column
k = interpolate_k(t, np.array([250.0]), np.array([1e4]),
                  h2o_vmr=np.array([1e-3]))
print('interp kappa', k.shape, 'min/max', k.min(), k.max())
"
```
Expected: `shape (1, 4, 2, 10, 15, 5)` and a positive `kappa` of order 1e-4 – 1e-1 m²/kg.

- [x] **Step 4: Commit**

```bash
git add climt/_data/picket_fence/correlated_k/trappist1e_hab1_lw.nc \
        climt/_data/picket_fence/correlated_k/trappist1e_hab1_sw.nc
git commit -m "feat(pf-tables): ship TRAPPIST-1e Hab1 (N2+CO2+H2O axis) LW+SW tables"
```

---

### Task 10: Download HITRAN CIA files for Titan

**Files:**
- Create: `climt/_data/picket_fence/cia/N2-N2_2018.cia`
- Create: `climt/_data/picket_fence/cia/N2-CH4_2018.cia`

These are public CC-BY ASCII files (typical size 1–5 MB).

- [x] **Step 1: Download**

From <https://hitran.org/cia/>, download:
- `N2-N2_2018.cia`
- `N2-CH4_2018.cia` (if absent, try `CH4-N2_2018.cia`; the CIA loader accepts either ordering)

Place them at the paths above.

- [x] **Step 2: Verify they parse**

Run:
```bash
conda run -n linepyline python -c "
from scripts.pf_table_builder.cia import load_cia_blocks
blocks = load_cia_blocks('climt/_data/picket_fence/cia/N2-N2_2018.cia')
print('N2-N2 temperatures:', sorted(blocks)[:5], '...', sorted(blocks)[-1])
print('block count:', len(blocks))
blocks = load_cia_blocks('climt/_data/picket_fence/cia/N2-CH4_2018.cia')
print('N2-CH4 temperatures:', sorted(blocks)[:5])
"
```
Expected: N2-N2 has ~10–15 temperature blocks covering ~40–400 K; N2-CH4 similar.

- [x] **Step 3: Commit (these are ~1–5 MB binary blobs)**

```bash
git add climt/_data/picket_fence/cia/*.cia
git commit -m "feat(pf-tables): ship HITRAN CIA files for Titan opacity build"
```

---

### Task 11: Build Titan tables

**Files:**
- Create: `climt/_data/picket_fence/correlated_k/titan_lw.nc`
- Create: `climt/_data/picket_fence/correlated_k/titan_sw.nc`

- [x] **Step 1: Build LW (includes CIA)**

```bash
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \
    --scenario titan --kind lw \
    --output climt/_data/picket_fence/correlated_k/titan_lw.nc \
    --ngpt 2 --dnu 0.5
```
Watch the output for the line `[titan/lw] adding 2 CIA pair(s)`.

- [x] **Step 2: Build SW**

```bash
conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \
    --scenario titan --kind sw \
    --output climt/_data/picket_fence/correlated_k/titan_sw.nc \
    --ngpt 2 --dnu 1.0
```

- [x] **Step 3: Sanity-check Titan LW window opacity is dominated by CIA below 500 cm⁻¹**

```bash
cd /Users/joymonteiro/github/climt && conda run -n climt python -c "
from climt._components.picket_fence.optics.correlated_k import load_k_table
t = load_k_table('titan_lw')
import numpy as np
k = t['k_coefficients']  # (ngas, nband, ngpt, nT, nP)
# Average k in band 0 (10-500 cm^-1) vs band 1 (500-1250) at T=140 K, p=1e4 Pa
iT, iP = 5, 7
print('band 0 (CIA window) <k>:', k[0, 0, :, iT, iP].mean())
print('band 1 (gas only)   <k>:', k[0, 1, :, iT, iP].mean())
"
```
Expected: band 0 average ≳ 10× band 1 average (because CIA dominates the rotational continuum).

- [x] **Step 4: Commit**

```bash
git add climt/_data/picket_fence/correlated_k/titan_lw.nc \
        climt/_data/picket_fence/correlated_k/titan_sw.nc
git commit -m "feat(pf-tables): ship Titan (N2+CH4+CIA) LW+SW tables"
```

---

### Task 12: Smoke tests + MANIFEST + end-to-end picket-fence run

**Files:**
- Modify: `tests/test_correlated_k_tables.py`
- Create: `tests/test_picket_fence_titan_trappist.py`
- Modify: `climt/_data/picket_fence/correlated_k/MANIFEST.md`

- [x] **Step 1: Extend the shipped-tables smoke test**

Edit `tests/test_correlated_k_tables.py`:

```python
SHIPPED = [
    "earth_low_res_lw", "earth_low_res_sw",
    "mars_lw",          "mars_sw",
    "venus_lw",         "venus_sw",
    "trappist1e_hab1_lw", "trappist1e_hab1_sw",
    "trappist1e_hab2_lw", "trappist1e_hab2_sw",
    "titan_lw",         "titan_sw",
]
```

- [x] **Step 2: Write an end-to-end picket-fence integration test**

Create `tests/test_picket_fence_titan_trappist.py`:

```python
"""End-to-end smoke tests: run PicketFence{Long,Short}wave with the new tables.

Verifies the netCDF schema we wrote matches what climt's loader expects and
that fluxes are physically sane (positive, finite, monotone in obvious ways).
"""
import numpy as np
import pytest

from climt import get_default_state, load_atmospheric_properties, reset_atmospheric_properties
from climt._components.picket_fence import (
    PicketFenceLongwave, PicketFenceShortwave,
)


@pytest.mark.parametrize("scenario,profile", [
    ("trappist1e_hab1", "trappist1e"),
    ("trappist1e_hab2", "trappist1e"),
    ("titan", "titan"),
])
def test_picket_fence_lw_runs(scenario, profile):
    """LW kernel produces non-negative OLR with the new table."""
    load_atmospheric_properties(profile)
    try:
        lw = PicketFenceLongwave(optics="correlated_k", table=f"{scenario}_lw")
        state = get_default_state([lw])
        tend, diag = lw(state)
        olr = diag["upwelling_longwave_flux_in_air"].values[-1, :]
        assert np.all(np.isfinite(olr))
        assert np.all(olr > 0), f"non-positive OLR for {scenario}: {olr}"
    finally:
        reset_atmospheric_properties()


@pytest.mark.parametrize("scenario,profile,star", [
    ("trappist1e_hab1", "trappist1e", "trappist1"),
    ("trappist1e_hab2", "trappist1e", "trappist1"),
    ("titan", "titan", "sun"),
])
def test_picket_fence_sw_runs(scenario, profile, star):
    """SW kernel produces non-negative downwelling flux with the new table."""
    load_atmospheric_properties(profile)
    try:
        sw = PicketFenceShortwave(
            optics="correlated_k",
            table=f"{scenario}_sw",
            stellar_spectrum=star,
        )
        state = get_default_state([sw])
        state["zenith_angle"].values[:] = np.pi / 4
        tend, diag = sw(state)
        flux_dn = diag["downwelling_shortwave_flux_in_air"].values
        assert np.all(np.isfinite(flux_dn))
        assert np.all(flux_dn >= 0)
    finally:
        reset_atmospheric_properties()
```

- [x] **Step 3: Run the new tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_correlated_k_tables.py tests/test_picket_fence_titan_trappist.py -v`
Expected: ALL PASS.

- [x] **Step 4: Update MANIFEST.md**

Append to `climt/_data/picket_fence/correlated_k/MANIFEST.md`:

```markdown

## linepyline-derived tables

Built offline by `scripts/generate_pf_tables_linepyline.py` from
HITRAN 2024 line data (via the [linepyline](https://github.com/...)
package) plus MT_CKD 4.3 (H2O continuum) and HITRAN CIA 2018 (Titan
N2-N2, N2-CH4).

| File | Atmosphere | Bands (cm^-1) | g-points | H2O axis | Notes |
|---|---|---|---|---|---|
| trappist1e_hab1_lw.nc | THAI Hab1: N2 + 400 ppm CO2 + H2O | 10, 500, 1250, 2500, 3250 | 2 | yes | Covers THAI Ben1 (radiation identical). |
| trappist1e_hab1_sw.nc | as above | 3250, 8000, 14000, 30000 | 2 | yes | Solar source from `trappist1.npz`. |
| trappist1e_hab2_lw.nc | THAI Hab2: pure 1-bar CO2 | 10, 500, 1250, 2500, 3250 | 2 | no | Covers THAI Ben2. |
| trappist1e_hab2_sw.nc | as above | 3250, 8000, 14000, 30000 | 2 | no | |
| titan_lw.nc | 95% N2 + 5% CH4, incl. N2-N2 & N2-CH4 CIA | 10, 500, 1250, 2500, 3250 | 2 | no | CIA dominates rotational window. |
| titan_sw.nc | as above | 3250, 8000, 14000, 30000 | 2 | no | Solar source from `sun.npz`. |

THAI protocol: <https://gmd.copernicus.org/articles/13/707/2020/>.
Titan CIA: Karman et al. (2019), HITRAN CIA database (CC-BY).
```

- [x] **Step 5: Update `docs/radiative_transfer/table_generation.rst`**

Append a section pointing users at the linepyline workflow:

```rst

linepyline-based tables
-----------------------

For scenarios outside Chaverot's Zenodo coverage (Titan, TRAPPIST-1e,
ad-hoc compositions), climt ships a second offline pipeline using
`linepyline <https://github.com/.../linepyline>`_ and HITRAN 2024 line
data. The driver is ``scripts/generate_pf_tables_linepyline.py``.

Building a Titan table (run from the ``linepyline`` conda env)::

    conda run -n linepyline python scripts/generate_pf_tables_linepyline.py \\
        --scenario titan --kind lw \\
        --output climt/_data/picket_fence/correlated_k/titan_lw.nc

Building a custom scenario: copy one of the entries in the ``SCENARIOS``
dict at the top of the driver, edit the absorber dict / VMR grid / band
edges, and re-run. The output netCDF drops straight into climt's
``correlated_k`` data directory and is picked up by name.

Titan in particular requires HITRAN CIA flat files in
``climt/_data/picket_fence/cia/``; see the README there for the download
list.
```

- [x] **Step 6: Commit**

```bash
git add tests/test_correlated_k_tables.py \
        tests/test_picket_fence_titan_trappist.py \
        climt/_data/picket_fence/correlated_k/MANIFEST.md \
        docs/radiative_transfer/table_generation.rst
git commit -m "test(pf-tables): smoke tests + docs for Titan/TRAPPIST-1e tables"
```

---

### Task 13: Validation — compare to linepyline reference column

**Files:**
- Create: `tests/test_picket_fence_linepyline_validation.py`

The strongest physics test: take one Titan and one TRAPPIST-1e profile, compute OLR with the picket-fence + new k-table, and compare to a direct linepyline LBL run on the same column. Expect agreement to within ~10–20% (limited by 4-band / 2-g-point resolution).

- [x] **Step 1: Failing test**

Create `tests/test_picket_fence_linepyline_validation.py`:

```python
"""LBL-vs-picket-fence consistency on one column per scenario."""
import numpy as np
import pytest

linepyline = pytest.importorskip("linepyline")


@pytest.mark.slow
def test_olr_matches_linepyline_for_hab2():
    """TRAPPIST-1e Hab2 OLR from picket-fence within 20% of linepyline LBL."""
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from climt._components.picket_fence import PicketFenceLongwave

    # Pick one (T, p) profile
    p = np.logspace(2, 5, 30)               # 100 Pa -> 1e5 Pa
    T = np.linspace(160, 250, 30)
    Ts, ps = 260.0, 1e5

    # Reference: linepyline LBL
    rtm = linepyline.rtm(background_gas=None, surface_gravity=3.721)
    ds = rtm.radiative_transfer(
        10.0, 3250.0, 0.5, p, ps, T, Ts,
        absorbers={"CO2": 1.0}, line_shape="lorentz", binning=True,
    )
    olr_lbl = float(ds.olr.integrate("nu"))

    # Picket-fence
    load_atmospheric_properties("trappist1e")
    try:
        lw = PicketFenceLongwave(optics="correlated_k",
                                 table="trappist1e_hab2_lw")
        # ... build a state with these (T, p), call lw, extract OLR ...
        from climt import get_default_state
        state = get_default_state([lw])
        # Overwrite the default profile with ours
        state["air_temperature"].values[:, 0] = T
        state["air_pressure"].values[:, 0] = p
        state["surface_temperature"].values[:] = Ts
        tend, diag = lw(state)
        olr_pf = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0])
    finally:
        reset_atmospheric_properties()

    rel_err = abs(olr_pf - olr_lbl) / olr_lbl
    assert rel_err < 0.20, f"OLR mismatch {rel_err:.1%}: pf={olr_pf:.1f}, lbl={olr_lbl:.1f}"
```

- [x] **Step 2: Run**

Run: `cd /Users/joymonteiro/github/climt && conda run -n linepyline python -m pytest tests/test_picket_fence_linepyline_validation.py -v`
(Requires both climt and linepyline in same env — alternatively, run linepyline in its env, dump `olr_lbl` to a file, then run the comparison in climt env reading that file.)

Expected: PASS (within 20%). If not, the most likely causes are: too-coarse `dnu` in the table build, too few g-points, or band edges chosen poorly for the gas mixture. Iterate `--ngpt`, `--dnu`, or band edges in the driver.

- [x] **Step 3: Commit**

```bash
git add tests/test_picket_fence_linepyline_validation.py
git commit -m "test(pf-tables): LBL consistency check for new tables"
```

---

## Summary

After this plan completes:

1. Six new k-tables are shipped: `trappist1e_hab1_{lw,sw}.nc`, `trappist1e_hab2_{lw,sw}.nc`, `titan_{lw,sw}.nc`.
2. A reusable, tested `pf_table_builder` package handles linepyline kappa sampling, band-binning, k-distribution quadrature, Planck/solar source partitioning, HITRAN CIA, and netCDF writing.
3. Students can launch picket-fence runs for TRAPPIST-1e THAI Hab1/Hab2/Ben1/Ben2 and Titan with one-line API calls (`PicketFenceLongwave(optics="correlated_k", table="trappist1e_hab1_lw")` + `load_atmospheric_properties("trappist1e")`).
4. The driver is extensible — adding a new scenario means appending a dict to `SCENARIOS` and re-running.

**Self-review notes (for the executor):**
- All step bodies contain complete code or exact shell commands.
- Smoke tests load from the climt env; build steps run from the linepyline env. Don't mix.
- The two long-running steps are Task 8 SW and Task 9 LW (H2O is the slowest molecule); budget ~1 hour each, run in background.
- If a step's output shape doesn't match what the next step expects, fix it in `generate_pf_tables_linepyline.py:build_table` before re-running — don't paper over with a transpose in the writer.
