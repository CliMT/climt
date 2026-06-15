# Picket-Fence Correlated-k Table Optimizer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an automated optimizer that places correlated-k band edges, decides band count, toggles continuum decoupling, builds PF LW tables, and iterates against line-by-line (LBL) radiative transfer until OLR agreement reaches RRTMG-vs-LBL accuracy.

**Architecture:** A greedy residual-driven loop in the `climt` env. The expensive, band-independent linepyline work (high-res κ sampling + LBL ground-truth RT) is computed once per grid and cached; each iteration only re-bins the cached κ-cube into a table (numpy) and runs PicketFence (njit). Where to split a band is chosen by a spectral correlation-breakdown metric computed from the cached κ-cube. The loop encodes the band-refinement campaign learnings: band width over g-points (ngpt is never touched), wide transparent window (non-splittable), continuum decoupling for moist-specific error.

**Tech Stack:** Python, numpy, scipy, xarray, linepyline (separate conda env), climt PicketFence, pytest.

**Design spec:** `docs/superpowers/specs/2026-06-16-pf-table-optimizer-design.md`

**Environments (do not mix):**
- linepyline env python: `/Users/joymonteiro/miniconda3/envs/linepyline/bin/python`
- climt env python: `/Users/joymonteiro/miniconda3/envs/climt/bin/python`
- All `pytest` and PicketFence steps run in the **climt** env (per project convention).
- Only `kappa_cache.py` and `truth.py` run in the **linepyline** env.

**Scope decisions locked for v1 (see spec):**
- Scenario: Earth LW, fixed CO₂ = 376 ppm (the `earth` SCENARIOS entry), no CO₂ axis. Code is structured to admit other scenarios/SW later but they are not validated here.
- Hard accept gate: `max(|OLR_err_moist|, |OLR_err_dry|) < OLR_TARGET` (default 5 W/m²) across the validation columns. Heating-rate RMSE is **computed and reported** from net-flux divergence as a secondary diagnostic; it is NOT a v1 gate (linepyline flux-profile wiring is a v2 follow-up). This is a deliberate narrowing of the spec's acceptance section and is reflected there.
- Knobs: band edges + count (primary), continuum decoupling (moist lever), T/p/X grid nodes (last-resort). ngpt is never modified.
- Stall behavior: stop, save best-so-far table + write a diagnostic report. No global refiner in v1.

---

## File Structure

New package `scripts/pf_optimizer/`:

- `__init__.py` — package marker.
- `config.py` — `OptimizerConfig` dataclass + the Earth-LW default (grid, outer band range, dnu, targets, limits).
- `kappa_cache.py` — (linepyline env) sample & cache high-res line-only and continuum-only κ cubes, keyed by grid hash.
- `truth.py` — (linepyline env) compute & cache LBL OLR per validation column (band-independent).
- `table_build.py` — (numpy) cached cubes + band_edges + decouple flag → temp `.nc` table. Reuses `pf_table_builder`.
- `diagnostics.py` — (numpy, pure) window-band detection, per-band OLR residual attribution, correlation-breakdown split metric, grid leave-one-out interpolation-error metric.
- `evaluate.py` — (climt) run PicketFenceLongwave on validation columns → residuals vs cached truth.
- `actions.py` — (numpy, pure) `choose_action`: the brain, encodes learnings + guardrails.
- `loop.py` — orchestrator: seed → build → evaluate → decide → apply → repeat; stopping; iteration logging.

Driver: `scripts/optimize_pf_table.py` (CLI).

Reused (refactored in Task 1): `scripts/generate_pf_tables_linepyline.py`,
`scripts/pf_table_builder/{kappa_sampling,k_distribution,planck_fraction,netcdf_writer}.py`.

Tests: `tests/pf_optimizer/test_*.py`.

---

## Task 1: Refactor the generator to expose a reusable table-assembly function (DRY)

The optimizer must turn a cached κ-cube + band edges into a `.nc` table using the
*exact* axis-rearrangement and decoupled-continuum logic already in
`generate_pf_tables_linepyline.build_table`. Extract that numpy-only tail into a
standalone function so both the CLI and the optimizer share it.

**Files:**
- Modify: `scripts/generate_pf_tables_linepyline.py`
- Create: `tests/pf_optimizer/__init__.py`
- Test: `tests/pf_optimizer/test_assemble_lw_table.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/pf_optimizer/test_assemble_lw_table.py
import os, sys
import numpy as np
import xarray as xr

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from generate_pf_tables_linepyline import assemble_lw_table


def _toy_cube():
    # (nT, nP, nX, nNu) line-only kappa on a coarse nu grid
    T_grid = np.array([200.0, 300.0])
    p_grid = np.array([1e3, 1e5])
    h2o_grid = np.array([1e-4, 1e-2])
    nu_grid = np.linspace(10.0, 200.0, 400)
    rng = np.random.default_rng(0)
    kappa_lines = np.abs(rng.normal(1e-3, 1e-4, (2, 2, 2, nu_grid.size)))
    kappa_cont = np.abs(rng.normal(1e-5, 1e-6, (2, 2, 2, nu_grid.size)))
    return T_grid, p_grid, h2o_grid, nu_grid, kappa_lines, kappa_cont


def test_assemble_coupled_writes_loadable_table(tmp_path):
    T, p, x, nu, kl, kc = _toy_cube()
    out = tmp_path / "toy_lw.nc"
    assemble_lw_table(
        str(out), kappa_lines=kl, kappa_cont=kc, nu_grid=nu,
        band_edges=np.array([10.0, 100.0, 200.0]), T_grid=T, p_grid=p,
        h2o_vmr_grid=x, ngpt=2, decouple_continuum=False,
        source="test:toy",
    )
    with xr.open_dataset(str(out)) as ds:
        # (gas, band, gpoint, T, p, h2o) = (1, 2, 2, 2, 2, 2)
        assert ds["k_coefficients"].shape == (1, 2, 2, 2, 2, 2)
        assert "continuum_kappa" not in ds   # coupled: no separate continuum
        assert ds["band_wavenumber_limits"].shape == (2, 2)


def test_assemble_decoupled_emits_continuum_var(tmp_path):
    T, p, x, nu, kl, kc = _toy_cube()
    out = tmp_path / "toy_lw_dec.nc"
    assemble_lw_table(
        str(out), kappa_lines=kl, kappa_cont=kc, nu_grid=nu,
        band_edges=np.array([10.0, 100.0, 200.0]), T_grid=T, p_grid=p,
        h2o_vmr_grid=x, ngpt=2, decouple_continuum=True, source="test:toy",
    )
    with xr.open_dataset(str(out)) as ds:
        assert ds["continuum_kappa"].shape == (2, 2, 2, 2)  # (band,T,p,X)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_assemble_lw_table.py -v`
Expected: FAIL — `ImportError: cannot import name 'assemble_lw_table'`.

- [ ] **Step 3: Extract `assemble_lw_table` from `build_table`**

In `scripts/generate_pf_tables_linepyline.py`, add this function (it lifts the
existing kappa→k_coeffs→moveaxis→write_lw_table logic, parameterised on the
cached cubes instead of re-sampling). Place it above `build_table`:

```python
def assemble_lw_table(output_path, *, kappa_lines, kappa_cont, nu_grid,
                      band_edges, T_grid, p_grid, h2o_vmr_grid, ngpt=2,
                      decouple_continuum=False, source="linepyline",
                      co2_vmr_grid=None):
    """Build & write an LW k-table from cached high-res kappa cubes.

    kappa_lines, kappa_cont: (nT, nP, nX, nNu) line-only and continuum-only
        mass-absorption cubes on `nu_grid`. Band-independent — the optimizer
        caches them once and calls this per candidate band structure.
    decouple_continuum: if True, the line CDF is built from kappa_lines alone
        and the band-averaged continuum is written separately (continuum_kappa);
        if False, kappa_lines + kappa_cont are binned together.
    """
    band_edges = np.asarray(band_edges, dtype=float)
    nband = len(band_edges) - 1
    nT, nP, nX = len(T_grid), len(p_grid), len(h2o_vmr_grid)

    continuum_band = None
    if decouple_continuum:
        kappa = kappa_lines
        continuum_band = np.zeros((nband, nT, nP, nX))
        for b in range(nband):
            m = (nu_grid >= band_edges[b]) & (nu_grid < band_edges[b + 1])
            continuum_band[b] = kappa_cont[:, :, :, m].mean(axis=-1)
    else:
        kappa = kappa_lines + kappa_cont

    k_coeffs, gpt_weights = kappa_to_k_coeffs(kappa, nu_grid, band_edges, ngpt)
    k_coeffs = k_coeffs[np.newaxis]  # leading ngas=1
    # (1, nT, nP, nX, nband, ngpt) -> (1, nband, ngpt, nT, nP, nX)
    k_coeffs = np.moveaxis(k_coeffs, (1, 2, 3, 4, 5), (3, 4, 5, 1, 2))

    planck_fraction = build_uniform_planck_fraction(band_edges, ngpt, T_grid)
    write_lw_table(
        output_path, k_coefficients=k_coeffs, gpoint_weights=gpt_weights,
        T_grid=T_grid, log_p_grid=np.log(p_grid), band_edges=band_edges,
        planck_fraction=planck_fraction, h2o_vmr_grid=h2o_vmr_grid,
        co2_vmr_grid=co2_vmr_grid, source=source,
        continuum_kappa=continuum_band,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_assemble_lw_table.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Confirm the CLI generator still works (no regression)**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_table_builder -q`
Expected: PASS (existing builder tests unchanged).

- [ ] **Step 6: Commit**

```bash
git add scripts/generate_pf_tables_linepyline.py tests/pf_optimizer/__init__.py tests/pf_optimizer/test_assemble_lw_table.py
git commit -m "refactor(pf): extract assemble_lw_table for optimizer reuse

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Optimizer config + Earth-LW default

**Files:**
- Create: `scripts/pf_optimizer/__init__.py` (empty)
- Create: `scripts/pf_optimizer/config.py`
- Test: `tests/pf_optimizer/test_config.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/pf_optimizer/test_config.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.config import OptimizerConfig, EARTH_LW


def test_earth_lw_defaults_are_consistent():
    c = EARTH_LW
    assert c.nu_min == 10.0 and c.nu_max == 3250.0
    assert c.seed_band_edges[0] == c.nu_min
    assert c.seed_band_edges[-1] == c.nu_max
    assert c.olr_target == 5.0
    assert c.max_bands >= len(c.seed_band_edges) - 1
    # grid axes are ascending
    assert np.all(np.diff(c.T_grid) > 0)
    assert np.all(np.diff(c.p_grid) > 0)


def test_grid_hash_changes_with_grid():
    c = EARTH_LW
    h0 = c.grid_hash()
    c2 = OptimizerConfig(**{**c.__dict__, "T_grid": c.T_grid * 1.01})
    assert c2.grid_hash() != h0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_config.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer'`.

- [ ] **Step 3: Implement config**

```python
# scripts/pf_optimizer/config.py
from __future__ import annotations
import hashlib
from dataclasses import dataclass, field
import numpy as np


@dataclass
class OptimizerConfig:
    # scenario / grid (Earth LW, fixed CO2)
    background_gas: str = "air"
    co2_vmr: float = 376e-6
    T_grid: np.ndarray = field(default_factory=lambda: np.array(
        [50., 110., 170., 230., 290., 350., 410., 470., 530., 590., 650., 1000.]))
    p_grid: np.ndarray = field(default_factory=lambda: np.array(
        [1., 10., 100., 1e3, 1e4, 1e5, 1e6, 1e7]))
    h2o_vmr_grid: np.ndarray = field(default_factory=lambda: np.array(
        [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]))
    nu_min: float = 10.0
    nu_max: float = 3250.0
    dnu: float = 0.1
    line_shape: str = "pseudovoigt"
    ngpt: int = 2
    # search seed + limits
    seed_band_edges: np.ndarray = field(default_factory=lambda: np.array(
        [10.0, 500.0, 1250.0, 1800.0, 3250.0]))
    max_bands: int = 14
    max_iters: int = 30
    min_band_width: float = 40.0   # cm^-1
    # targets
    olr_target: float = 5.0        # W/m^2
    moist_lever_threshold: float = 3.0  # W/m^2 moist-minus-dry to trigger continuum
    # window detection: bands whose grid-mean optical proxy is below this are
    # treated as transparent windows and are non-splittable.
    window_transmittance: float = 0.85

    def grid_hash(self) -> str:
        h = hashlib.sha1()
        for a in (self.T_grid, self.p_grid, self.h2o_vmr_grid):
            h.update(np.asarray(a, dtype="f8").tobytes())
        h.update(np.array([self.nu_min, self.nu_max, self.dnu,
                           self.co2_vmr], dtype="f8").tobytes())
        h.update(self.line_shape.encode())
        return h.hexdigest()[:12]


EARTH_LW = OptimizerConfig()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_config.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/pf_optimizer/__init__.py scripts/pf_optimizer/config.py tests/pf_optimizer/test_config.py
git commit -m "feat(pf-optimizer): add OptimizerConfig + Earth-LW default

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: `table_build.py` — cached cubes + edges → temp table

Thin climt-env wrapper over `assemble_lw_table` that loads a κ-cache `.npz` and
emits a temp `.nc`. Keeps the loop code clean and lets `evaluate.py` consume a path.

**Files:**
- Create: `scripts/pf_optimizer/table_build.py`
- Test: `tests/pf_optimizer/test_table_build.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/pf_optimizer/test_table_build.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.table_build import build_table_from_cache


def _write_cache(path):
    nu = np.linspace(10.0, 200.0, 400)
    rng = np.random.default_rng(1)
    np.savez(
        path,
        nu_grid=nu,
        kappa_lines=np.abs(rng.normal(1e-3, 1e-4, (2, 2, 2, nu.size))),
        kappa_cont=np.abs(rng.normal(1e-5, 1e-6, (2, 2, 2, nu.size))),
        T_grid=np.array([200.0, 300.0]),
        p_grid=np.array([1e3, 1e5]),
        h2o_vmr_grid=np.array([1e-4, 1e-2]),
    )


def test_build_table_from_cache_returns_loadable_path(tmp_path):
    cache = tmp_path / "kc.npz"
    _write_cache(str(cache))
    out = build_table_from_cache(
        str(cache), band_edges=np.array([10.0, 100.0, 200.0]),
        decouple_continuum=False, ngpt=2, out_path=str(tmp_path / "t.nc"))
    assert os.path.isfile(out)
    import xarray as xr
    with xr.open_dataset(out) as ds:
        assert ds["k_coefficients"].shape == (1, 2, 2, 2, 2, 2)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_table_build.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer.table_build'`.

- [ ] **Step 3: Implement**

```python
# scripts/pf_optimizer/table_build.py
from __future__ import annotations
import os, sys
import numpy as np

_SCRIPTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)
from generate_pf_tables_linepyline import assemble_lw_table


def load_cache(cache_path):
    d = np.load(cache_path)
    return {k: d[k] for k in d.files}


def build_table_from_cache(cache_path, *, band_edges, decouple_continuum,
                           ngpt, out_path, source="optimizer"):
    c = load_cache(cache_path)
    assemble_lw_table(
        out_path,
        kappa_lines=c["kappa_lines"], kappa_cont=c["kappa_cont"],
        nu_grid=c["nu_grid"], band_edges=np.asarray(band_edges, dtype=float),
        T_grid=c["T_grid"], p_grid=c["p_grid"],
        h2o_vmr_grid=c["h2o_vmr_grid"], ngpt=ngpt,
        decouple_continuum=decouple_continuum, source=source,
    )
    return out_path
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_table_build.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/pf_optimizer/table_build.py tests/pf_optimizer/test_table_build.py
git commit -m "feat(pf-optimizer): build temp table from cached kappa cube

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: `diagnostics.py` — split metric, window detection, attribution

The heart of the "smart band placement." All pure numpy — tested on synthetic
cubes with planted structure.

**Files:**
- Create: `scripts/pf_optimizer/diagnostics.py`
- Test: `tests/pf_optimizer/test_diagnostics.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/pf_optimizer/test_diagnostics.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.diagnostics import (
    correlation_breakdown_split, is_window_band, worst_dry_band,
)


def test_split_lands_at_planted_breakdown():
    # Band 0..100 cm^-1; left half layer-ordering A, right half ordering B.
    nu = np.linspace(0.0, 100.0, 200)
    nlayer = 4
    kappa = np.zeros((nlayer, nu.size))
    left = nu < 50.0
    # ordering A: layer 0 strongest; ordering B (right): layer 0 weakest
    A = np.array([4.0, 3.0, 2.0, 1.0])
    B = A[::-1]
    kappa[:, left] = A[:, None] * (1.0 + 0.01 * nu[left])
    kappa[:, ~left] = B[:, None] * (1.0 + 0.01 * nu[~left])
    nu_split = correlation_breakdown_split(kappa, nu)
    assert 45.0 < nu_split < 55.0   # within one bin of the planted boundary


def test_window_band_flagged_transparent():
    # near-zero kappa over the band -> transmittance ~1 -> window
    nu = np.linspace(800.0, 1200.0, 100)
    kappa = np.full((3, nu.size), 1e-12)   # m^2/kg, negligible
    air_mass = 1e4  # kg/m^2 representative column
    assert is_window_band(kappa, air_mass, threshold=0.85) is True
    kappa_opaque = np.full((3, nu.size), 1.0)
    assert is_window_band(kappa_opaque, air_mass, threshold=0.85) is False


def test_worst_dry_band_picks_largest_abs_residual():
    per_band = np.array([0.5, -7.0, 2.0])      # LBL - PF per band (dry)
    edges = np.array([10.0, 500.0, 1250.0, 3250.0])
    assert worst_dry_band(per_band, edges) == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_diagnostics.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer.diagnostics'`.

- [ ] **Step 3: Implement**

```python
# scripts/pf_optimizer/diagnostics.py
from __future__ import annotations
import numpy as np

_TINY = 1e-30


def correlation_breakdown_split(kappa_band, nu_band):
    """Wavenumber at which the per-layer kappa fingerprint changes most.

    kappa_band: (nlayer, nnu) kappa within one band, for a representative
        column's layers. nu_band: (nnu,) wavenumbers (ascending).

    The per-nu "fingerprint" is the z-scored log-kappa vector across layers.
    We segment the band into two parts minimising total within-segment
    dispersion of that fingerprint (1D, exhaustive over interior cut points);
    the returned split is the nu at the optimal boundary. This directly targets
    the correlated-k assumption failure: where the layer-ordering of kappa stops
    being coherent within the band.
    """
    logk = np.log(np.asarray(kappa_band, dtype=float) + _TINY)
    # z-score per layer across nu so absolute magnitude doesn't dominate
    f = (logk - logk.mean(axis=1, keepdims=True))
    sd = f.std(axis=1, keepdims=True)
    f = f / np.where(sd > 0, sd, 1.0)            # (nlayer, nnu)
    F = f.T                                       # (nnu, nlayer) fingerprints
    nnu = F.shape[0]
    if nnu < 4:
        return float(0.5 * (nu_band[0] + nu_band[-1]))
    # prefix sums for O(n) within-segment SS-to-mean
    cs = np.cumsum(F, axis=0)
    cs2 = np.cumsum((F * F).sum(axis=1))
    total2 = cs2  # cumulative sum of squared norms

    def seg_ss(i, j):  # SS to mean over rows [i, j)
        n = j - i
        s = cs[j - 1] - (cs[i - 1] if i > 0 else 0.0)
        ss = total2[j - 1] - (total2[i - 1] if i > 0 else 0.0)
        return ss - (s @ s) / n

    best_k, best_val = nnu // 2, np.inf
    for k in range(2, nnu - 1):           # interior cut: left=[0,k), right=[k,nnu)
        val = seg_ss(0, k) + seg_ss(k, nnu)
        if val < best_val:
            best_val, best_k = val, k
    return float(nu_band[best_k])


def is_window_band(kappa_band, air_mass, threshold=0.85):
    """True if the band is a transparent window (mean transmittance > threshold).

    kappa_band: (nlayer, nnu) m^2/kg. air_mass: representative column mass
    (kg/m^2). Transmittance = exp(-tau), tau = mean_kappa * air_mass.
    """
    mean_kappa = float(np.mean(kappa_band))
    transmittance = float(np.exp(-mean_kappa * air_mass))
    return transmittance > threshold


def worst_dry_band(per_band_residual, band_edges):
    """Index of the band with the largest |LBL - PF| dry OLR residual."""
    return int(np.argmax(np.abs(np.asarray(per_band_residual))))


def leave_one_out_axis_error(kappa_cube, axis):
    """Max relative kappa error from dropping each interior node on `axis`
    and linearly re-interpolating (in log-kappa). Used as the last-resort
    grid-refinement trigger. Returns (max_rel_err, worst_node_index)."""
    k = np.log(np.asarray(kappa_cube, dtype=float) + _TINY)
    n = k.shape[axis]
    worst_err, worst_i = 0.0, -1
    for i in range(1, n - 1):
        lo = np.take(k, i - 1, axis=axis)
        hi = np.take(k, i + 1, axis=axis)
        approx = 0.5 * (lo + hi)
        true = np.take(k, i, axis=axis)
        err = float(np.max(np.abs(approx - true)))
        if err > worst_err:
            worst_err, worst_i = err, i
    return worst_err, worst_i
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_diagnostics.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/pf_optimizer/diagnostics.py tests/pf_optimizer/test_diagnostics.py
git commit -m "feat(pf-optimizer): split metric, window detection, residual attribution

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: `actions.py` — the brain (`choose_action`)

Pure decision logic. No I/O, no linepyline, no climt. Encodes the learnings.

**Files:**
- Create: `scripts/pf_optimizer/actions.py`
- Test: `tests/pf_optimizer/test_actions.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/pf_optimizer/test_actions.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.actions import choose_action, Residual, State


def _state(edges, decoupled=False):
    return State(band_edges=np.array(edges), decoupled=decoupled,
                 max_bands=14, min_band_width=40.0, window_transmittance=0.85)


def test_accept_when_under_target():
    r = Residual(olr_err_moist=2.0, olr_err_dry=-3.0,
                 per_band_dry=np.array([1.0, -2.0]), kappa_per_band=None,
                 air_mass=1e4)
    act = choose_action(r, _state([10.0, 1250.0, 3250.0]),
                        olr_target=5.0, moist_lever_threshold=3.0)
    assert act["kind"] == "accept"


def test_moist_specific_error_triggers_continuum():
    r = Residual(olr_err_moist=18.0, olr_err_dry=3.0,
                 per_band_dry=np.array([1.0, 2.0]), kappa_per_band=None,
                 air_mass=1e4)
    act = choose_action(r, _state([10.0, 1250.0, 3250.0], decoupled=False),
                        olr_target=5.0, moist_lever_threshold=3.0)
    assert act["kind"] == "decouple_continuum"


def test_splits_worst_dry_band_when_not_window():
    # band 1 (1250-3250) has the worst dry residual and is opaque -> split it
    nu = np.linspace(1250.0, 3250.0, 200)
    nlayer = 4
    k = np.zeros((nlayer, nu.size)) + 1.0
    k[:, nu > 2250.0] = np.array([4., 1., 3., 2.])[:, None]  # breakdown midband
    kappa_per_band = [np.full((nlayer, 10), 1e-12), k]       # band0 window, band1 opaque
    r = Residual(olr_err_moist=4.0, olr_err_dry=9.0,
                 per_band_dry=np.array([0.5, 9.0]),
                 kappa_per_band=kappa_per_band, air_mass=1e4)
    st = _state([10.0, 1250.0, 3250.0], decoupled=True)  # continuum already done
    act = choose_action(r, st, olr_target=5.0, moist_lever_threshold=3.0)
    assert act["kind"] == "split_band"
    assert act["band_index"] == 1
    assert 1250.0 < act["new_edge"] < 3250.0


def test_window_band_not_split_even_if_worst():
    nu = np.linspace(800.0, 1250.0, 200)
    nlayer = 3
    win = np.full((nlayer, nu.size), 1e-12)   # transparent
    other = np.full((nlayer, 50), 1.0)        # opaque, smaller residual
    kappa_per_band = [win, other]
    r = Residual(olr_err_moist=4.0, olr_err_dry=9.0,
                 per_band_dry=np.array([9.0, 2.0]),  # window is worst
                 kappa_per_band=kappa_per_band, air_mass=1e4)
    st = _state([800.0, 1250.0, 3250.0], decoupled=True)
    act = choose_action(r, st, olr_target=5.0, moist_lever_threshold=3.0)
    # must NOT split the window (band 0); splits the opaque band 1 instead
    assert act["kind"] == "split_band"
    assert act["band_index"] == 1


def test_stall_when_at_max_bands_and_decoupled():
    r = Residual(olr_err_moist=9.0, olr_err_dry=9.0,
                 per_band_dry=np.array([9.0, 9.0]), kappa_per_band=None,
                 air_mass=1e4)
    st = State(band_edges=np.array([10.0, 1250.0, 3250.0]), decoupled=True,
               max_bands=2, min_band_width=40.0, window_transmittance=0.85)
    act = choose_action(r, st, olr_target=5.0, moist_lever_threshold=3.0)
    assert act["kind"] == "stall"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_actions.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer.actions'`.

- [ ] **Step 3: Implement**

```python
# scripts/pf_optimizer/actions.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, List
import numpy as np

from pf_optimizer.diagnostics import (
    correlation_breakdown_split, is_window_band, worst_dry_band,
)


@dataclass
class Residual:
    olr_err_moist: float          # LBL - PF, total (W/m^2)
    olr_err_dry: float
    per_band_dry: np.ndarray      # LBL - PF per band, dry (W/m^2)
    kappa_per_band: Optional[List[np.ndarray]]  # each (nlayer, nnu_in_band)
    air_mass: float               # representative column mass (kg/m^2)


@dataclass
class State:
    band_edges: np.ndarray
    decoupled: bool
    max_bands: int
    min_band_width: float
    window_transmittance: float


def _band_nu(edges, b):
    return float(edges[b]), float(edges[b + 1])


def choose_action(resid, state, *, olr_target, moist_lever_threshold):
    """Return the single next action dict. Priority: accept > continuum >
    split worst non-window dry band > stall."""
    worst_olr = max(abs(resid.olr_err_moist), abs(resid.olr_err_dry))
    if worst_olr < olr_target:
        return {"kind": "accept"}

    # 1. Moist-specific error -> decouple continuum (exp #24-25)
    moist_excess = abs(resid.olr_err_moist) - abs(resid.olr_err_dry)
    if moist_excess > moist_lever_threshold and not state.decoupled:
        return {"kind": "decouple_continuum"}

    # 2. Split the worst NON-WINDOW dry band (exp #23 + #18 guardrail)
    nband = len(state.band_edges) - 1
    if (nband < state.max_bands) and resid.kappa_per_band is not None:
        order = np.argsort(np.abs(resid.per_band_dry))[::-1]  # worst first
        for b in order:
            kb = resid.kappa_per_band[b]
            if is_window_band(kb, resid.air_mass, state.window_transmittance):
                continue
            lo, hi = _band_nu(state.band_edges, b)
            if (hi - lo) < 2 * state.min_band_width:
                continue
            nu_band = np.linspace(lo, hi, kb.shape[1])
            new_edge = correlation_breakdown_split(kb, nu_band)
            # keep the split away from edges by min_band_width
            new_edge = float(np.clip(new_edge, lo + state.min_band_width,
                                     hi - state.min_band_width))
            return {"kind": "split_band", "band_index": int(b),
                    "new_edge": new_edge}

    # 3. Nothing left to try
    return {"kind": "stall"}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_actions.py -v`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/pf_optimizer/actions.py tests/pf_optimizer/test_actions.py
git commit -m "feat(pf-optimizer): choose_action brain encoding band-refinement learnings

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: `evaluate.py` — PicketFence on validation columns → residuals

**Files:**
- Create: `scripts/pf_optimizer/evaluate.py`
- Test: `tests/pf_optimizer/test_evaluate.py`

`evaluate.py` runs PicketFence on each validation column and re-bins the cached
LBL OLR spectrum onto the candidate band edges for per-band attribution. The
cached truth is an `.npz` with `nu`, `olr_spec` (W/m²/cm⁻¹), `total` for each of
the moist/dry cases (produced by Task 7). The test fakes both a table and a
truth file to keep it climt-only and fast.

- [ ] **Step 1: Write the failing test**

```python
# tests/pf_optimizer/test_evaluate.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.evaluate import rebin_lbl_per_band, residual_from_totals


def test_rebin_lbl_per_band_sums_to_total():
    nu = np.arange(10.0, 200.0, 0.1)
    olr_spec = np.ones_like(nu)            # 1 W/m^2/cm^-1 flat
    edges = np.array([10.0, 100.0, 200.0])
    per_band = rebin_lbl_per_band(nu, olr_spec, edges, dnu=0.1)
    assert per_band.shape == (2,)
    assert abs(per_band.sum() - np.sum(olr_spec) * 0.1) < 1e-6


def test_residual_from_totals():
    r = residual_from_totals(lbl_moist=250.0, pf_moist=240.0,
                             lbl_dry=350.0, pf_dry=349.0,
                             per_band_dry=np.array([1.0, 0.0]),
                             kappa_per_band=None, air_mass=1e4)
    assert r.olr_err_moist == 10.0
    assert r.olr_err_dry == 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_evaluate.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer.evaluate'`.

- [ ] **Step 3: Implement**

```python
# scripts/pf_optimizer/evaluate.py
from __future__ import annotations
import numpy as np

from pf_optimizer.actions import Residual


def rebin_lbl_per_band(nu, olr_spec, band_edges, dnu):
    """Integrate the LBL OLR spectrum into each band (W/m^2)."""
    edges = np.asarray(band_edges)
    out = np.zeros(len(edges) - 1)
    for b in range(len(edges) - 1):
        sel = (nu >= edges[b]) & (nu < edges[b + 1])
        out[b] = float(np.sum(olr_spec[sel]) * dnu)
    return out


def residual_from_totals(*, lbl_moist, pf_moist, lbl_dry, pf_dry,
                         per_band_dry, kappa_per_band, air_mass):
    return Residual(
        olr_err_moist=float(lbl_moist - pf_moist),
        olr_err_dry=float(lbl_dry - pf_dry),
        per_band_dry=np.asarray(per_band_dry, dtype=float),
        kappa_per_band=kappa_per_band, air_mass=float(air_mass),
    )


def pf_olr_on_column(table_path, p_int, T, q, Ts, co2, nz):
    """Total + per-band PF OLR for one column. climt env only.

    p_int: (nz+1,) interface pressures (Pa), surface first per climt convention.
    Returns (total_olr, per_band_olr)."""
    from climt import get_default_state, get_grid
    from climt._components.picket_fence import PicketFenceLongwave
    pf = PicketFenceLongwave(optics="correlated_k", table=table_path)
    grid = get_grid(nx=1, ny=1, nz=nz)
    state = get_default_state([pf], grid_state=grid)
    state["air_temperature"].values[:, 0, 0] = T
    state["specific_humidity"].values[:, 0, 0] = q
    state["surface_temperature"].values[:] = Ts
    if "surface_longwave_emissivity" in state:
        state["surface_longwave_emissivity"].values[:] = 1.0
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2
    if "mole_fraction_of_ozone_in_air" in state:
        state["mole_fraction_of_ozone_in_air"].values[:] = 0.0
    _t, diag = pf(state)
    total = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    per_band = np.asarray(
        diag["upwelling_longwave_flux_in_air_per_band"].values[-1, 0, 0, :])
    return total, per_band
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_evaluate.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/pf_optimizer/evaluate.py tests/pf_optimizer/test_evaluate.py
git commit -m "feat(pf-optimizer): PF column eval + LBL re-bin + residual builder

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: `kappa_cache.py` + `truth.py` (linepyline env)

The two band-independent caches. Built together because they share the grid and
both run in the linepyline env. Tested only for cache I/O contract (no
linepyline in CI); the heavy sampling is exercised in the Task 9 smoke run.

**Files:**
- Create: `scripts/pf_optimizer/kappa_cache.py`
- Create: `scripts/pf_optimizer/truth.py`
- Test: `tests/pf_optimizer/test_cache_contract.py`

- [ ] **Step 1: Write the failing test (I/O contract only)**

```python
# tests/pf_optimizer/test_cache_contract.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.kappa_cache import cache_path_for, save_kappa_cache
from pf_optimizer.config import EARTH_LW


def test_cache_path_is_grid_keyed(tmp_path):
    p = cache_path_for(EARTH_LW, str(tmp_path))
    assert EARTH_LW.grid_hash() in os.path.basename(p)
    assert p.endswith(".npz")


def test_save_kappa_cache_roundtrips(tmp_path):
    nu = np.linspace(10.0, 50.0, 20)
    kl = np.ones((2, 2, 2, 20)); kc = np.ones((2, 2, 2, 20)) * 0.1
    out = str(tmp_path / "kc.npz")
    save_kappa_cache(out, nu_grid=nu, kappa_lines=kl, kappa_cont=kc,
                     T_grid=np.array([200., 300.]), p_grid=np.array([1e3, 1e5]),
                     h2o_vmr_grid=np.array([1e-4, 1e-2]))
    d = np.load(out)
    assert set(["nu_grid", "kappa_lines", "kappa_cont", "T_grid",
                "p_grid", "h2o_vmr_grid"]).issubset(set(d.files))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_cache_contract.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer.kappa_cache'`.

- [ ] **Step 3: Implement `kappa_cache.py`**

```python
# scripts/pf_optimizer/kappa_cache.py  (RUN IN linepyline ENV)
from __future__ import annotations
import os, sys
import numpy as np

_SCRIPTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)
from pf_table_builder.kappa_sampling import sample_kappa_grid


def cache_path_for(cfg, cache_dir):
    return os.path.join(cache_dir, f"kappa_cache_{cfg.grid_hash()}.npz")


def save_kappa_cache(path, *, nu_grid, kappa_lines, kappa_cont,
                     T_grid, p_grid, h2o_vmr_grid):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    np.savez(path, nu_grid=nu_grid.astype("f4"),
             kappa_lines=kappa_lines.astype("f4"),
             kappa_cont=kappa_cont.astype("f4"),
             T_grid=T_grid, p_grid=p_grid, h2o_vmr_grid=h2o_vmr_grid)


def build_kappa_cache(cfg, cache_dir):
    """Sample line-only and continuum-only kappa cubes once and cache them."""
    nu_grid = np.arange(cfg.nu_min, cfg.nu_max + cfg.dnu / 2, cfg.dnu)
    absorbers = {"CO2": cfg.co2_vmr}
    common = dict(background_gas=cfg.background_gas, absorbers=absorbers,
                  T_grid=cfg.T_grid, p_grid=cfg.p_grid, nu_grid=nu_grid,
                  h2o_vmr_grid=cfg.h2o_vmr_grid, line_shape=cfg.line_shape,
                  binning=False)
    kappa_lines = sample_kappa_grid(include_mtckd_continuum=False, **common)
    kappa_total = sample_kappa_grid(include_mtckd_continuum=True, **common)
    kappa_cont = np.clip(kappa_total - kappa_lines, 0.0, None)
    path = cache_path_for(cfg, cache_dir)
    save_kappa_cache(path, nu_grid=nu_grid, kappa_lines=kappa_lines,
                     kappa_cont=kappa_cont, T_grid=cfg.T_grid,
                     p_grid=cfg.p_grid, h2o_vmr_grid=cfg.h2o_vmr_grid)
    return path


if __name__ == "__main__":
    # Usage (linepyline env):
    #   python -m pf_optimizer.kappa_cache <cache_dir>
    import argparse
    from pf_optimizer.config import EARTH_LW
    ap = argparse.ArgumentParser()
    ap.add_argument("cache_dir")
    a = ap.parse_args()
    print(build_kappa_cache(EARTH_LW, a.cache_dir))
```

- [ ] **Step 4: Implement `truth.py`**

```python
# scripts/pf_optimizer/truth.py  (RUN IN linepyline ENV)
from __future__ import annotations
import os
import numpy as np


def truth_path_for(cfg, column_name, cache_dir):
    return os.path.join(cache_dir, f"lbl_truth_{cfg.grid_hash()}_{column_name}.npz")


def build_lbl_truth(cfg, columns, cache_dir):
    """Compute & cache LBL OLR spectrum per column (moist & dry). Band-independent.

    columns: list of dicts with keys name, p_mid_Pa (ascending: TOA..surface),
        T, q_moist, Ts, co2.
    """
    import linepyline as lpl
    import xarray as xr
    rtm = lpl.rtm(background_gas=cfg.background_gas, use_numba=True)
    out_paths = []
    for col in columns:
        order = np.argsort(col["p_mid_Pa"])
        p_s = col["p_mid_Pa"][order]
        T_s = col["T"][order]
        q_s = col["q_moist"][order]
        cases = {}
        for tag, q in (("moist", q_s), ("dry", np.zeros_like(q_s))):
            p_da = xr.DataArray(p_s, dims=("p",), coords={"p": p_s})
            T_da = xr.DataArray(T_s, dims=("p",), coords={"p": p_s})
            q_da = xr.DataArray(q, dims=("p",), coords={"p": p_s})
            ds = rtm.radiative_transfer(
                cfg.nu_min, cfg.nu_max, cfg.dnu, p_da, float(p_s[-1]),
                T_da, float(col["Ts"]), q=q_da, absorbers={"CO2": col["co2"]},
                D=1.66, line_shape=cfg.line_shape, binning=False,
                include_mtckd_continuum=True)
            nu = np.asarray(ds["nu"].values)
            olr_spec = np.asarray(ds["olr"].values)
            cases[f"nu_{tag}"] = nu
            cases[f"olr_spec_{tag}"] = olr_spec
            cases[f"total_{tag}"] = float(np.sum(olr_spec) * cfg.dnu)
        path = truth_path_for(cfg, col["name"], cache_dir)
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        np.savez(path, **cases)
        out_paths.append(path)
    return out_paths
```

- [ ] **Step 5: Run the contract test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_cache_contract.py -v`
Expected: PASS (2 passed).

- [ ] **Step 6: Commit**

```bash
git add scripts/pf_optimizer/kappa_cache.py scripts/pf_optimizer/truth.py tests/pf_optimizer/test_cache_contract.py
git commit -m "feat(pf-optimizer): grid-keyed kappa cache + LBL truth (linepyline env)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8: `loop.py` — orchestrator + logging + stopping

Ties it together. Pure orchestration over the pieces above; the only external
calls are `build_table_from_cache`, `pf_olr_on_column`, and `choose_action`. Uses
an injected `pf_eval_fn` so the loop is testable without climt.

**Files:**
- Create: `scripts/pf_optimizer/loop.py`
- Test: `tests/pf_optimizer/test_loop.py`

- [ ] **Step 1: Write the failing test (fake eval, no climt/linepyline)**

```python
# tests/pf_optimizer/test_loop.py
import os, sys
import numpy as np
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from pf_optimizer.loop import optimize
from pf_optimizer.config import OptimizerConfig


def test_loop_converges_by_splitting(tmp_path, monkeypatch):
    # Fake PF eval: dry OLR error shrinks each time a band is added; the loop
    # must reach accept by splitting (never touching ngpt).
    state_iter = {"n": 0}

    def fake_eval(band_edges, decoupled, cfg):
        nband = len(band_edges) - 1
        dry_err = max(0.0, 12.0 - 4.0 * (nband - 2))   # 12, 8, 4, 0...
        moist_err = 2.0
        per_band_dry = np.zeros(nband); per_band_dry[-1] = dry_err
        # opaque kappa so the band is splittable; planted breakdown midband
        nu = np.linspace(band_edges[-2], band_edges[-1], 50)
        kb = np.ones((3, 50)); kb[:, nu > nu.mean()] = np.array([3., 1., 2.])[:, None]
        kappa_per_band = [np.ones((3, 50)) for _ in range(nband - 1)] + [kb]
        return dict(olr_err_moist=moist_err, olr_err_dry=dry_err,
                    per_band_dry=per_band_dry, kappa_per_band=kappa_per_band,
                    air_mass=1e4)

    cfg = OptimizerConfig(seed_band_edges=np.array([10.0, 1250.0, 3250.0]),
                          max_bands=8, max_iters=10, olr_target=5.0)
    result = optimize(cfg, pf_eval_fn=fake_eval, log_path=str(tmp_path / "log.md"))
    assert result["status"] == "accepted"
    assert len(result["band_edges"]) - 1 >= 3          # at least one split happened
    assert os.path.isfile(str(tmp_path / "log.md"))
    # ngpt never changed
    assert all(a["kind"] != "set_ngpt" for a in result["history"])


def test_loop_stalls_and_reports_best(tmp_path):
    def stuck_eval(band_edges, decoupled, cfg):
        nband = len(band_edges) - 1
        # window band -> non-splittable; error never clears
        return dict(olr_err_moist=9.0, olr_err_dry=9.0,
                    per_band_dry=np.full(nband, 9.0),
                    kappa_per_band=[np.full((3, 20), 1e-12) for _ in range(nband)],
                    air_mass=1e4)

    cfg = OptimizerConfig(seed_band_edges=np.array([10.0, 3250.0]),
                          max_bands=8, max_iters=5, olr_target=5.0)
    result = optimize(cfg, pf_eval_fn=stuck_eval, log_path=str(tmp_path / "log.md"))
    assert result["status"] == "stalled"
    assert "best_band_edges" in result
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_loop.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pf_optimizer.loop'`.

- [ ] **Step 3: Implement**

```python
# scripts/pf_optimizer/loop.py
from __future__ import annotations
import numpy as np

from pf_optimizer.actions import choose_action, Residual, State


def _log(log_path, lines):
    with open(log_path, "a") as f:
        for ln in lines:
            f.write(ln + "\n")


def optimize(cfg, *, pf_eval_fn, log_path):
    """Greedy residual-driven band optimizer.

    pf_eval_fn(band_edges, decoupled, cfg) -> dict with keys
        olr_err_moist, olr_err_dry, per_band_dry, kappa_per_band, air_mass.
    Returns a result dict with status in {accepted, stalled, exhausted}.
    """
    band_edges = np.array(cfg.seed_band_edges, dtype=float)
    decoupled = False
    history = []
    best = {"score": np.inf, "band_edges": band_edges.copy(),
            "decoupled": decoupled}
    _log(log_path, [f"# PF table optimizer run", "",
                    f"- seed bands: {band_edges.tolist()}",
                    f"- OLR target: {cfg.olr_target} W/m^2", ""])

    for it in range(cfg.max_iters):
        ev = pf_eval_fn(band_edges, decoupled, cfg)
        resid = Residual(olr_err_moist=ev["olr_err_moist"],
                         olr_err_dry=ev["olr_err_dry"],
                         per_band_dry=np.asarray(ev["per_band_dry"]),
                         kappa_per_band=ev["kappa_per_band"],
                         air_mass=ev["air_mass"])
        score = max(abs(resid.olr_err_moist), abs(resid.olr_err_dry))
        if score < best["score"]:
            best = {"score": score, "band_edges": band_edges.copy(),
                    "decoupled": decoupled}
        st = State(band_edges=band_edges, decoupled=decoupled,
                   max_bands=cfg.max_bands, min_band_width=cfg.min_band_width,
                   window_transmittance=cfg.window_transmittance)
        act = choose_action(resid, st, olr_target=cfg.olr_target,
                            moist_lever_threshold=cfg.moist_lever_threshold)
        history.append(act)
        _log(log_path, [
            f"## Iter {it}  bands={len(band_edges) - 1} decoupled={decoupled}",
            f"- OLR err moist={resid.olr_err_moist:+.2f}  "
            f"dry={resid.olr_err_dry:+.2f}  worst={score:.2f}",
            f"- action: {act}", ""])

        if act["kind"] == "accept":
            _log(log_path, [f"**ACCEPTED at iter {it}.**", ""])
            return {"status": "accepted", "band_edges": band_edges.tolist(),
                    "decoupled": decoupled, "history": history}
        if act["kind"] == "decouple_continuum":
            decoupled = True
            continue
        if act["kind"] == "split_band":
            band_edges = np.sort(np.append(band_edges, act["new_edge"]))
            continue
        # stall
        _log(log_path, [f"**STALLED at iter {it}** — reporting best-so-far.",
                        f"- best worst-OLR-err: {best['score']:.2f} W/m^2",
                        f"- best bands: {best['band_edges'].tolist()}", ""])
        return {"status": "stalled", "history": history,
                "best_band_edges": best["band_edges"].tolist(),
                "best_decoupled": best["decoupled"]}

    _log(log_path, [f"**EXHAUSTED max_iters={cfg.max_iters}.**",
                    f"- best worst-OLR-err: {best['score']:.2f} W/m^2", ""])
    return {"status": "exhausted", "history": history,
            "best_band_edges": best["band_edges"].tolist(),
            "best_decoupled": best["decoupled"]}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_loop.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Run the whole optimizer unit suite**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer -q`
Expected: PASS (all tasks' tests green).

- [ ] **Step 6: Commit**

```bash
git add scripts/pf_optimizer/loop.py tests/pf_optimizer/test_loop.py
git commit -m "feat(pf-optimizer): greedy orchestrator loop with logging + stall handling

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 9: CLI driver `scripts/optimize_pf_table.py` + real eval wiring

Wires the real climt PF eval into the loop and provides the CLI. The real
`pf_eval_fn` builds a temp table from the κ-cache, runs PF on each validation
column, loads the cached LBL truth, and assembles the `Residual` inputs. It also
extracts `kappa_per_band` for the candidate edges from the κ-cache interpolated
onto a representative column.

**Files:**
- Create: `scripts/optimize_pf_table.py`
- Test: `tests/pf_optimizer/test_cli_smoke.py`

- [ ] **Step 1: Write the failing test (CLI arg parsing + dry validation)**

```python
# tests/pf_optimizer/test_cli_smoke.py
import os, sys, subprocess
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_cli_help_runs():
    p = subprocess.run(
        [sys.executable, os.path.join(_REPO, "scripts", "optimize_pf_table.py"),
         "--help"], capture_output=True, text=True)
    assert p.returncode == 0
    assert "--cache-dir" in p.stdout
    assert "--kappa-cache" in p.stdout
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_cli_smoke.py -v`
Expected: FAIL — file does not exist / FileNotFoundError.

- [ ] **Step 3: Implement the CLI + real eval**

```python
# scripts/optimize_pf_table.py  (RUN IN climt ENV)
"""Optimize a picket-fence LW correlated-k band structure against LBL.

Prereqs (linepyline env, once per grid):
    python -m pf_optimizer.kappa_cache <cache_dir>
    # plus LBL truth via pf_optimizer.truth.build_lbl_truth(...)

Then (climt env):
    python scripts/optimize_pf_table.py --cache-dir debug_data/pf_opt \\
        --kappa-cache debug_data/pf_opt/kappa_cache_<hash>.npz \\
        --profile debug_data/forward_profile.npz \\
        --out climt/_data/picket_fence/correlated_k/earth_opt_lw.nc
"""
from __future__ import annotations
import argparse, os, sys, tempfile
import numpy as np

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPTS)
from pf_optimizer.config import EARTH_LW
from pf_optimizer.table_build import build_table_from_cache, load_cache
from pf_optimizer.evaluate import (
    rebin_lbl_per_band, residual_from_totals, pf_olr_on_column)
from pf_optimizer.truth import truth_path_for
from pf_optimizer.loop import optimize


def _interp_kappa_per_band(cache, band_edges, p_int, T, q_vmr, p_grid_col):
    """kappa_lines+cont interpolated onto a representative column, split per band.

    Returns list of (nlayer, nnu_in_band) arrays for the correlation metric.
    Interpolates the cached (nT,nP,nX,nNu) cube to the column's (T,p,X) layers
    via nearest grid index (sufficient for the ranking-based split metric)."""
    nu = cache["nu_grid"]
    kl, kc = cache["kappa_lines"], cache["kappa_cont"]
    Tg, pg, xg = cache["T_grid"], cache["p_grid"], cache["h2o_vmr_grid"]
    iT = np.argmin(np.abs(Tg[None, :] - T[:, None]), axis=1)
    iP = np.argmin(np.abs(np.log(pg)[None, :] - np.log(p_grid_col)[:, None]), axis=1)
    iX = np.argmin(np.abs(np.log(xg)[None, :] - np.log(np.clip(q_vmr, xg[0], xg[-1]))[:, None]), axis=1)
    col = (kl + kc)[iT, iP, iX, :]            # (nlayer, nnu)
    edges = np.asarray(band_edges)
    out = []
    for b in range(len(edges) - 1):
        m = (nu >= edges[b]) & (nu < edges[b + 1])
        out.append(col[:, m])
    return out


def make_pf_eval(cfg, cache, truth, profile, out_dir):
    p_mid = profile["p_mid_Pa"]; T = profile["T_ref"]
    q_moist = profile["q_ref_moist"]; Ts = float(profile["T_surf"])
    co2 = float(profile["CO2"]); nz = len(p_mid)
    # interface pressures, surface-first per climt convention
    order_desc = np.argsort(p_mid)[::-1]      # surface(high p) -> TOA
    p_climt = p_mid[order_desc]; T_climt = T[order_desc]; q_climt = q_moist[order_desc]
    p_int = np.empty(nz + 1)
    p_int[1:-1] = np.sqrt(p_climt[:-1] * p_climt[1:])
    p_int[0] = p_climt[0]; p_int[-1] = p_climt[-1] ** 2 / p_int[-2]
    air_mass = float(np.sum(np.abs(np.diff(p_int))) / 9.81)
    # crude H2O VMR for kappa indexing (mass q -> vmr)
    q_vmr = q_moist * (28.97 / 18.015)

    def pf_eval(band_edges, decoupled, _cfg):
        tbl = os.path.join(out_dir, "cand.nc")
        build_table_from_cache(
            cfg_cache_path, band_edges=band_edges, decouple_continuum=decoupled,
            ngpt=cfg.ngpt, out_path=tbl)
        pf_moist, pb_moist = pf_olr_on_column(tbl, p_int, T_climt, q_climt, Ts, co2, nz)
        pf_dry, pb_dry = pf_olr_on_column(tbl, p_int, T_climt, np.zeros_like(q_climt), Ts, co2, nz)
        lbl_pb_dry = rebin_lbl_per_band(truth["nu_dry"], truth["olr_spec_dry"],
                                        band_edges, cfg.dnu)
        kappa_per_band = _interp_kappa_per_band(
            cache, band_edges, p_int, T_climt, q_vmr[order_desc], p_climt)
        r = residual_from_totals(
            lbl_moist=float(truth["total_moist"]), pf_moist=pf_moist,
            lbl_dry=float(truth["total_dry"]), pf_dry=pf_dry,
            per_band_dry=(lbl_pb_dry - pb_dry),
            kappa_per_band=kappa_per_band, air_mass=air_mass)
        return dict(olr_err_moist=r.olr_err_moist, olr_err_dry=r.olr_err_dry,
                    per_band_dry=r.per_band_dry, kappa_per_band=r.kappa_per_band,
                    air_mass=r.air_mass)
    return pf_eval


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--kappa-cache", required=True,
                    help="path to kappa_cache_<hash>.npz (built in linepyline env)")
    ap.add_argument("--profile", required=True,
                    help="forward_profile.npz with p_mid_Pa,T_ref,q_ref_moist,T_surf,CO2")
    ap.add_argument("--truth", required=True,
                    help="lbl_truth_<hash>_<col>.npz (built in linepyline env)")
    ap.add_argument("--out", required=True, help="final table .nc path")
    ap.add_argument("--log", default=None)
    a = ap.parse_args()

    global cfg_cache_path
    cfg_cache_path = a.kappa_cache
    cfg = EARTH_LW
    cache = load_cache(a.kappa_cache)
    truth = {k: np.load(a.truth)[k] for k in np.load(a.truth).files}
    profile = {k: np.load(a.profile, allow_pickle=True)[k]
               for k in np.load(a.profile, allow_pickle=True).files}
    os.makedirs(a.cache_dir, exist_ok=True)
    log_path = a.log or os.path.join(a.cache_dir, "optimizer_log.md")
    pf_eval = make_pf_eval(cfg, cache, truth, profile, a.cache_dir)
    result = optimize(cfg, pf_eval_fn=pf_eval, log_path=log_path)
    # write the final/best table
    edges = result.get("band_edges") or result["best_band_edges"]
    dec = result.get("decoupled", result.get("best_decoupled", False))
    build_table_from_cache(a.kappa_cache, band_edges=np.array(edges),
                           decouple_continuum=dec, ngpt=cfg.ngpt, out_path=a.out)
    print(f"status={result['status']}  bands={len(edges) - 1}  -> {a.out}")
    print(f"log: {log_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `/Users/joymonteiro/miniconda3/envs/climt/bin/python -m pytest tests/pf_optimizer/test_cli_smoke.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/optimize_pf_table.py tests/pf_optimizer/test_cli_smoke.py
git commit -m "feat(pf-optimizer): CLI driver wiring real climt PF eval into the loop

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 10: End-to-end smoke run on the real Earth grid (manual, gated)

This is the only step that runs linepyline + a full optimize. It is a manual
acceptance gate, not a CI test (it takes minutes and needs both envs). Document
the result in the experiment log per project convention.

- [ ] **Step 1: Build the κ-cache (linepyline env, once)**

Run:
```bash
/Users/joymonteiro/miniconda3/envs/linepyline/bin/python -m pf_optimizer.kappa_cache debug_data/pf_opt
```
(Run from repo root with `PYTHONPATH=scripts`.)
Expected: prints `debug_data/pf_opt/kappa_cache_<hash>.npz`; file exists.

- [ ] **Step 2: Build the LBL truth for the forward profile (linepyline env)**

Run a short inline script that loads `debug_data/forward_profile.npz`, packs it
as one `columns` entry (`name="forward"`, `co2=CO2`), and calls
`pf_optimizer.truth.build_lbl_truth(EARTH_LW, columns, "debug_data/pf_opt")`.
Expected: `debug_data/pf_opt/lbl_truth_<hash>_forward.npz` exists with keys
`nu_dry`, `olr_spec_dry`, `total_dry`, `total_moist`.

- [ ] **Step 3: Run the optimizer (climt env)**

Run:
```bash
/Users/joymonteiro/miniconda3/envs/climt/bin/python scripts/optimize_pf_table.py \
  --cache-dir debug_data/pf_opt \
  --kappa-cache debug_data/pf_opt/kappa_cache_<hash>.npz \
  --profile debug_data/forward_profile.npz \
  --truth debug_data/pf_opt/lbl_truth_<hash>_forward.npz \
  --out climt/_data/picket_fence/correlated_k/earth_opt_lw.nc
```
Expected: terminates with `status=accepted` (or `stalled`/`exhausted` with a
best-so-far); prints the band count and writes `earth_opt_lw.nc` and
`debug_data/pf_opt/optimizer_log.md`.

- [ ] **Step 4: Sanity-check the produced table against the hand-tuned winsplit**

Run:
```bash
/Users/joymonteiro/miniconda3/envs/climt/bin/python scripts/experiments/eval_band_structure.py earth_opt_lw
```
Expected: moist & dry `|LBL-PF|` are within a few W/m² and no worse than
`earth_low_res_lw_co2refined` winsplit on the same profile. Record the numbers.

- [ ] **Step 5: Append the result to the experiment log**

Append a dated entry to
`docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md` (the active PF
experiment log) recording: final band edges, whether the continuum was
decoupled, iteration count, accepted/stalled, and the eval_band_structure
numbers vs winsplit. Per project convention, log tries/learnings/rejections.

- [ ] **Step 6: Commit**

```bash
git add docs/superpowers/plans/2026-05-16-pf-co2-band-refinement.md
git commit -m "docs(pf-optimizer): log first end-to-end Earth-LW optimizer run

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review Notes

- **Spec coverage:** band edges+count (Tasks 4,5,8), continuum decoupling
  (Tasks 1,5,8 — `decouple_continuum` action), grid nodes as last-resort
  (`leave_one_out_axis_error` in Task 4 is implemented and available; the loop
  wiring of a `refine_grid` action is intentionally deferred — see note below);
  correlation-breakdown split metric (Task 4); cached LBL truth + κ-cube (Task
  7); accept gate on OLR (Tasks 5,8); stall→best-so-far+report (Task 8);
  logging (Task 8, Task 10 Step 5); Earth-LW-first/general-ready (config-driven).
- **Deferred within v1 (flagged, not silent):** (a) the grid-refinement *action*
  in `choose_action`/`loop` — the diagnostic exists but the loop does not yet
  trigger a re-cache; band splitting + continuum are expected to clear Earth LW
  without it. (b) heating-rate gating — OLR is the v1 gate. Both are documented
  here and in the spec so the executor doesn't mistake them for omissions.
- **Type consistency:** `Residual`/`State` dataclasses defined in `actions.py`
  (Task 5) and imported unchanged by `evaluate.py` (Task 6) and `loop.py` (Task
  8). `choose_action` returns dicts with `kind` in {accept, decouple_continuum,
  split_band, stall}; `loop.py` handles exactly those. `build_table_from_cache`
  signature is identical in Task 3, Task 9. `assemble_lw_table` signature
  identical in Task 1 and Task 3's call.
- **Placeholder scan:** no TBD/TODO; every code step is complete. The one
  `<hash>` token in Task 9/10 commands is a runtime value the executor fills
  from Step-1 output, not a plan gap.
```
