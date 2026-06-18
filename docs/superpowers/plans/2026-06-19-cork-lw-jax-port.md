# CORK Longwave Differentiable JAX Port — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fully `jax.grad`-differentiable, `jit`/`vmap`-compatible JAX implementation of CORK's correlated-k optics and longwave radiative transfer, validated against the existing Numba kernels, dispatched automatically when `CorkLongwaveRadiation` receives `jax.Array` inputs.

**Architecture:** Separate JAX kernel modules sit beside the untouched Numba kernels (Approach B). `optics/correlated_k_jax.py` computes optical depth `tau` via integer-bracket gather + differentiable multilinear blend; `lw/kernels_jax.py` computes differentiable Planck sources and a `lax.scan` transport solver; `lw/component.py` gains a JAX dispatch path that reuses these. The Numba code is kept solely as the parity oracle.

**Tech Stack:** Python, JAX (`jax.numpy`, `jax.lax.scan`, `jax.grad`, `jax.vmap`), NumPy, Numba (oracle only), sympl, pytest.

## Global Constraints

- Branch `jax-research` — never merged to `develop`/main. JAX path is primary; Numba is the validation oracle and must remain **untouched**.
- Every new JAX module sets float64 at import: `jax.config.update("jax_enable_x64", True)` (radiation fluxes require it; oracle is float64). Tests also set `jax.config.update("jax_platform_name", "cpu")`.
- Scope is the **correlated-k** optics mode with a **7-D additive-CO2, premixed-background** k-table (the `earth_low_res_lw` family: `k.shape == (ngas, nband, ngpt, nT, nP, nX, nC)`, `gas_names == ["effective"]`, H2O+CO2 axes, additive overlap, optional 4-D continuum). The JAX path activates only for this case; all other cases (parmentier optics, ESFT overlap, non-premixed multi-gas, `diagnostics_level > 0`) fall back to the existing Numba path.
- Gradients are required **w.r.t. atmospheric state only** (T, surface T, composition, pressure). The k-table and all grids are static `jnp` constants — never differentiated.
- Forward parity tolerance vs the oracle: `rtol=1e-6, atol` per-quantity as specified (not bit-for-bit). Gradient finite-difference tolerance: `rtol≈1e-3`.
- Multilinear blend conventions copied verbatim from the oracle: bracket index `i = clip(searchsorted(grid, v) - 1, 0, n-2)`, weight `f = clip((v-grid[i])/(grid[i+1]-grid[i]), 0, 1)`; CO2 axis uses log-interpolation with `FLOOR = 1e-40`; continuum uses `exp(trilinear(log_cont))` with floor `1e-40`; diffusivity factor `D = 1.66`.

**Reference files (read, do not modify):**
- Optics oracle: `climt/_components/cork/optics/correlated_k.py` (`_ck_bracket`, `_ck_txx7`, `_ck_txx_cont`, `_ck_tau_additive_co2_kernel`, `_additive_co2_fast`, `_CO2_INTERP_LOGK=True`).
- LW oracle: `climt/_components/cork/lw/kernels.py` (`planck_sources_kernel`, `_lw_transport_kernel`, `lw_transport`).
- Component: `climt/_components/cork/lw/component.py` (`array_call`, `_correlated_k_optics`).
- Helpers: `climt/_components/cork/common.py` (`compute_heating_rate`, `compute_column_amount`, `MOLAR_MASS`, `MOLAR_MASS_DRY_AIR`).
- Test idiom: `tests/test_cork_ck_njit.py` (`_build`, `_inputs`, `_reference_tau`).

---

## File Structure

- Create `climt/_components/cork/optics/correlated_k_jax.py` — differentiable `compute_tau_jax` + bracket/gather helpers.
- Create `climt/_components/cork/lw/kernels_jax.py` — `planck_sources_jax`, `lw_transport_jax`, `heating_rate_jax`, and the `cork_lw_jax` orchestrator + `CorkLWTable` container.
- Modify `climt/_components/cork/lw/component.py` — build a `CorkLWTable` at construction (correlated-k mode only); dispatch to the JAX path in `array_call` when inputs are `jax.Array`.
- Create `tests/test_cork_optics_jax_parity.py` — optics τ parity + gradient.
- Create `tests/test_cork_lw_kernels_jax_parity.py` — Planck + transport parity + gradient.
- Create `tests/test_cork_lw_jax_integration.py` — full-component parity, jit/vmap-safety, end-to-end grad.

---

## Phase 1 — Differentiable optics (`correlated_k_jax.py`)

### Task 1: `compute_tau_jax` — bracket + gather + multilinear blend

**Files:**
- Create: `climt/_components/cork/optics/correlated_k_jax.py`
- Test: `tests/test_cork_optics_jax_parity.py`

**Interfaces:**
- Produces:
  - `_bracket(grid, v) -> (idx, frac)` where `grid: (n,)`, `v: (...)`, `idx: int (...)`, `frac: float64 (...)`.
  - `compute_tau_jax(T, log_p, log_x, log_c, gas_amounts, k, T_grid, p_grid_log, log_x_grid, log_c_grid, has_cont, log_cont, co2_logk) -> tau` with shapes: `T,log_p,log_x,log_c: (nlev,ncol)`; `gas_amounts: (ngas,nlev,ncol)`; `k: (ngas,nband,ngpt,nT,nP,nX,nC)`; grids 1-D; `log_cont: (nband,nT,nP,nX)`; `has_cont,co2_logk: bool`; returns `tau: (nband,ngpt,nlev,ncol)`.
  - The arguments are exactly the prepared arrays `_ck_tau_additive_co2_kernel` consumes, so parity is on identical inputs.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cork_optics_jax_parity.py
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest


def _build(tmp_path, k, cont=None):
    from scripts.cork_table_builder.netcdf_writer import write_lw_table
    from climt._components.cork.optics.correlated_k import load_k_table
    ngas, nband, ngpt, nT, nP, nX, nC = k.shape
    write_lw_table(
        str(tmp_path / "t.nc"),
        k_coefficients=k,
        gpoint_weights=np.full((nband, ngpt), 0.5),
        T_grid=np.linspace(200.0, 320.0, nT),
        log_p_grid=np.linspace(np.log(1e2), np.log(1e5), nP),
        band_edges=np.linspace(10.0, 3250.0, nband + 1),
        planck_fraction=np.full((nband, ngpt, nT), 1.0 / ngpt),
        h2o_vmr_grid=np.array([1e-6, 1e-3, 1e-1])[:nX],
        co2_vmr_grid=np.logspace(-5, -2, nC),
        continuum_kappa=cont,
        gas_names=("effective",),
    )
    return load_k_table(str(tmp_path / "t.nc"))


def _inputs(nlev=6, ncol=4, seed=0):
    rng = np.random.RandomState(seed)
    T = rng.uniform(210.0, 310.0, size=(nlev, ncol))
    p = np.exp(rng.uniform(np.log(2e2), np.log(8e4), size=(nlev, ncol)))
    h2o = np.exp(rng.uniform(np.log(1e-6), np.log(1e-1), size=(nlev, ncol)))
    co2 = np.exp(rng.uniform(np.log(1e-5), np.log(1e-2), size=(nlev, ncol)))
    return T, p, h2o, co2


def _prep(table, T, p, h2o, co2):
    """Replicate _additive_co2_fast host prep so JAX gets identical inputs."""
    x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
    c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
    log_p = np.log(np.maximum(p, 1.0))
    x_clamped = np.clip(h2o, float(x_grid[0]), float(x_grid[-1]))
    log_x = np.log(np.maximum(x_clamped, 1e-30))
    c_clamped = np.clip(co2, float(c_grid[0]), float(c_grid[-1]))
    log_c = np.log(np.maximum(c_clamped, 1e-30))
    return log_p, log_x, log_c


def _jax_args(table, T, log_p, log_x, log_c, gas):
    has_cont = ("continuum_kappa" in table
                and np.asarray(table["continuum_kappa"]).ndim == 4)
    log_cont = (np.log(np.maximum(np.asarray(table["continuum_kappa"], np.float64), 1e-40))
                if has_cont else np.zeros((table["k_coefficients"].shape[1], 1, 1, 1)))
    return dict(
        T=jnp.asarray(T), log_p=jnp.asarray(log_p), log_x=jnp.asarray(log_x),
        log_c=jnp.asarray(log_c), gas_amounts=jnp.asarray(gas),
        k=jnp.asarray(table["k_coefficients"], dtype=jnp.float64),
        T_grid=jnp.asarray(table["temperature_grid"], dtype=jnp.float64),
        p_grid_log=jnp.asarray(table["pressure_grid_log"], dtype=jnp.float64),
        log_x_grid=jnp.log(jnp.maximum(jnp.asarray(table["h2o_vmr_grid"], jnp.float64), 1e-30)),
        log_c_grid=jnp.log(jnp.maximum(jnp.asarray(table["co2_vmr_grid"], jnp.float64), 1e-30)),
        has_cont=has_cont, log_cont=jnp.asarray(log_cont), co2_logk=True,
    )


@pytest.mark.parametrize("with_cont", [False, True])
def test_tau_jax_matches_oracle(tmp_path, with_cont):
    from climt._components.cork.optics.correlated_k import compute_ck_optical_depth
    from climt._components.cork.optics.correlated_k_jax import compute_tau_jax
    rng = np.random.RandomState(1)
    k = rng.uniform(1e-6, 1e-2, size=(1, 5, 3, 4, 5, 3, 6))
    cont = rng.uniform(1e-7, 1e-3, size=(5, 4, 5, 3)) if with_cont else None
    table = _build(tmp_path, k, cont=cont)
    T, p, h2o, co2 = _inputs()
    gas = np.full((1,) + T.shape, 90.0)
    tau_oracle = compute_ck_optical_depth(table, T, p, gas, h2o_vmr=h2o, co2_vmr=co2)
    log_p, log_x, log_c = _prep(table, T, p, h2o, co2)
    tau_jax = compute_tau_jax(**_jax_args(table, T, log_p, log_x, log_c, gas))
    np.testing.assert_allclose(np.asarray(tau_jax), tau_oracle, rtol=1e-6, atol=1e-12)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cork_optics_jax_parity.py -v`
Expected: FAIL — `ModuleNotFoundError: climt._components.cork.optics.correlated_k_jax`.

- [ ] **Step 3: Write minimal implementation**

```python
# climt/_components/cork/optics/correlated_k_jax.py
"""Differentiable JAX correlated-k optics (7-D additive-CO2 path).

Mirrors the Numba oracle ``_ck_tau_additive_co2_kernel`` in
``correlated_k.py``. Integer bracket indices are gathered from the static
k-table; gradients flow through the fractional interpolation weights, which
are smooth functions of the atmospheric state.
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

FLOOR = 1e-40


def _bracket(grid, v):
    """Index + fraction for multilinear interpolation. Matches _ck_bracket."""
    n = grid.shape[0]
    idx = jnp.clip(jnp.searchsorted(grid, v) - 1, 0, n - 2)
    lo = grid[idx]
    frac = jnp.clip((v - lo) / (grid[idx + 1] - lo), 0.0, 1.0)
    return idx, frac


def _txx7(k, iT, fT, iP, fP, iX, fX, iC):
    """8-corner trilinear over (T,P,X) at fixed C index. Returns
    (ngas,nband,ngpt,nlev,ncol)."""
    def g(a, b, c):
        return k[:, :, :, iT + a, iP + b, iX + c, iC]
    x0 = (g(0, 0, 0) * (1 - fT) * (1 - fP) + g(1, 0, 0) * fT * (1 - fP)
          + g(0, 1, 0) * (1 - fT) * fP + g(1, 1, 0) * fT * fP)
    x1 = (g(0, 0, 1) * (1 - fT) * (1 - fP) + g(1, 0, 1) * fT * (1 - fP)
          + g(0, 1, 1) * (1 - fT) * fP + g(1, 1, 1) * fT * fP)
    return x0 * (1 - fX) + x1 * fX


def _txx_cont(log_cont, iT, fT, iP, fP, iX, fX):
    """8-corner trilinear over (T,P,X) of log continuum. Returns (nband,nlev,ncol)."""
    def g(a, b, c):
        return log_cont[:, iT + a, iP + b, iX + c]
    x0 = (g(0, 0, 0) * (1 - fT) * (1 - fP) + g(1, 0, 0) * fT * (1 - fP)
          + g(0, 1, 0) * (1 - fT) * fP + g(1, 1, 0) * fT * fP)
    x1 = (g(0, 0, 1) * (1 - fT) * (1 - fP) + g(1, 0, 1) * fT * (1 - fP)
          + g(0, 1, 1) * (1 - fT) * fP + g(1, 1, 1) * fT * fP)
    return x0 * (1 - fX) + x1 * fX


def compute_tau_jax(T, log_p, log_x, log_c, gas_amounts,
                    k, T_grid, p_grid_log, log_x_grid, log_c_grid,
                    has_cont, log_cont, co2_logk):
    iT, fT = _bracket(T_grid, T)
    iP, fP = _bracket(p_grid_log, log_p)
    iX, fX = _bracket(log_x_grid, log_x)
    iC, fC = _bracket(log_c_grid, log_c)

    c0 = _txx7(k, iT, fT, iP, fP, iX, fX, iC)       # (ngas,nband,ngpt,nlev,ncol)
    c1 = _txx7(k, iT, fT, iP, fP, iX, fX, iC + 1)
    if co2_logk:
        l0 = jnp.log(jnp.maximum(c0, FLOOR))
        l1 = jnp.log(jnp.maximum(c1, FLOOR))
        kv = jnp.exp(l0 * (1 - fC) + l1 * fC)
    else:
        kv = c0 * (1 - fC) + c1 * fC

    # acc[ib,igp,nlev,ncol] = sum over gas of kv * gas_amount
    tau = jnp.sum(kv * gas_amounts[:, None, None, :, :], axis=0)

    if has_cont:
        cont_val = jnp.exp(_txx_cont(log_cont, iT, fT, iP, fP, iX, fX))  # (nband,nlev,ncol)
        tau = tau + cont_val[:, None, :, :] * gas_amounts[0][None, None, :, :]
    return tau
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_optics_jax_parity.py -v`
Expected: PASS (both `with_cont` cases).

- [ ] **Step 5: Commit**

```bash
git add climt/_components/cork/optics/correlated_k_jax.py tests/test_cork_optics_jax_parity.py
git commit -m "feat(cork): differentiable JAX correlated-k optical depth (Phase 1)"
```

### Task 2: Gradient check for `compute_tau_jax`

**Files:**
- Test: `tests/test_cork_optics_jax_parity.py` (add to existing file)

**Interfaces:**
- Consumes: `compute_tau_jax`, `_jax_args`, `_prep`, `_build`, `_inputs` from Task 1.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_cork_optics_jax_parity.py
def test_tau_jax_grad_matches_finite_difference(tmp_path):
    from climt._components.cork.optics.correlated_k_jax import compute_tau_jax
    rng = np.random.RandomState(7)
    k = rng.uniform(1e-6, 1e-2, size=(1, 3, 2, 4, 5, 3, 6))
    table = _build(tmp_path, k)
    T, p, h2o, co2 = _inputs(nlev=5, ncol=2, seed=4)
    gas = np.full((1,) + T.shape, 100.0)
    log_p, log_x, log_c = _prep(table, T, p, h2o, co2)
    args = _jax_args(table, T, log_p, log_x, log_c, gas)
    static = {key: args[key] for key in args if key != "T"}

    def scalar(Tj):
        return jnp.sum(compute_tau_jax(T=Tj, **static))

    g = np.asarray(jax.grad(scalar)(args["T"]))
    assert np.all(np.isfinite(g))
    # central finite difference at one interior point well inside grid brackets
    eps = 1e-3
    i, j = 2, 1
    Tp = args["T"].at[i, j].add(eps)
    Tm = args["T"].at[i, j].add(-eps)
    fd = (float(scalar(Tp)) - float(scalar(Tm))) / (2 * eps)
    np.testing.assert_allclose(g[i, j], fd, rtol=1e-3, atol=1e-6)
```

- [ ] **Step 2: Run test to verify it fails (or passes immediately)**

Run: `pytest tests/test_cork_optics_jax_parity.py::test_tau_jax_grad_matches_finite_difference -v`
Expected: PASS if Task 1 is correct (the function is already differentiable). If it FAILS with `nan`, inspect for a `log(0)`/`0/0` — fix by confirming the `FLOOR`/`clip` guards from Task 1 are present.

- [ ] **Step 3: (only if failing) fix guards**

No new code expected. If `nan` gradients appear, verify `_bracket` uses `jnp.clip` and the `co2_logk` path uses `jnp.maximum(c, FLOOR)` before `log` (as written in Task 1).

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_optics_jax_parity.py -v`
Expected: PASS (all tests).

- [ ] **Step 5: Commit**

```bash
git add tests/test_cork_optics_jax_parity.py
git commit -m "test(cork): gradient finite-difference check for JAX optical depth"
```

---

## Phase 2 — LW kernels (`kernels_jax.py`)

### Task 3: `planck_sources_jax` — differentiable Planck sources

**Files:**
- Create: `climt/_components/cork/lw/kernels_jax.py`
- Test: `tests/test_cork_lw_kernels_jax_parity.py`

**Interfaces:**
- Produces: `planck_sources_jax(planck_frac, T_grid, T, T_surf, sigma) -> (planck_src, surf_src)` with `planck_frac: (nband,ngpt,nT)`, `T_grid: (nT,)`, `T: (nlev,ncol)`, `T_surf: (ncol,)`, `sigma: float`; returns `planck_src: (nband,ngpt,nlev,ncol)`, `surf_src: (nband,ngpt,ncol)`. Valid only for the additive case where `nband==nband_orig`, `ngpt==ngpt_orig`, `is_esft=False`.
- Consumes: `_bracket` from `correlated_k_jax`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cork_lw_kernels_jax_parity.py
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest


def test_planck_sources_jax_matches_oracle():
    from climt._components.cork.lw.kernels import planck_sources_kernel
    from climt._components.cork.lw.kernels_jax import planck_sources_jax
    rng = np.random.RandomState(0)
    nband, ngpt, nT, nlev, ncol = 4, 3, 6, 7, 5
    planck_frac = rng.uniform(0.0, 1.0, size=(nband, ngpt, nT))
    T_grid = np.linspace(180.0, 330.0, nT)
    T = rng.uniform(200.0, 315.0, size=(nlev, ncol))
    T_surf = rng.uniform(260.0, 305.0, size=(ncol,))
    sigma = 5.670374419e-8

    planck_src = np.zeros((nband, ngpt, nlev, ncol))
    surf_src = np.zeros((nband, ngpt, ncol))
    planck_sources_kernel(
        np.ascontiguousarray(planck_frac),
        np.ascontiguousarray(T_grid, dtype=np.float64),
        np.ascontiguousarray(T, dtype=np.float64),
        np.ascontiguousarray(T_surf, dtype=np.float64),
        float(sigma), nband, ngpt, False, ngpt, nband,
        planck_src, surf_src,
    )
    ps_jax, ss_jax = planck_sources_jax(
        jnp.asarray(planck_frac), jnp.asarray(T_grid),
        jnp.asarray(T), jnp.asarray(T_surf), sigma)
    np.testing.assert_allclose(np.asarray(ps_jax), planck_src, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(np.asarray(ss_jax), surf_src, rtol=1e-6, atol=1e-9)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cork_lw_kernels_jax_parity.py -v`
Expected: FAIL — `ModuleNotFoundError: climt._components.cork.lw.kernels_jax`.

- [ ] **Step 3: Write minimal implementation**

```python
# climt/_components/cork/lw/kernels_jax.py
"""Differentiable JAX longwave kernels (additive correlated-k path).

Mirrors the Numba oracles ``planck_sources_kernel`` and
``_lw_transport_kernel``. Valid for the non-ESFT additive case where the
g-point grid is not expanded (nband==nband_orig, ngpt==ngpt_orig).
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from jax import lax

from ..optics.correlated_k_jax import _bracket

DIFFUSIVITY_FACTOR = 1.66


def planck_sources_jax(planck_frac, T_grid, T, T_surf, sigma):
    iT, fT = _bracket(T_grid, T)            # (nlev,ncol)
    pf0 = planck_frac[:, :, iT]             # (nband,ngpt,nlev,ncol)
    pf1 = planck_frac[:, :, iT + 1]
    frac_l = pf0 * (1.0 - fT) + pf1 * fT
    planck_src = frac_l * (sigma * T ** 4)

    iTs, fTs = _bracket(T_grid, T_surf)     # (ncol,)
    pfs0 = planck_frac[:, :, iTs]           # (nband,ngpt,ncol)
    pfs1 = planck_frac[:, :, iTs + 1]
    frac_s = pfs0 * (1.0 - fTs) + pfs1 * fTs
    surf_src = frac_s * (sigma * T_surf ** 4)
    return planck_src, surf_src
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_lw_kernels_jax_parity.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_components/cork/lw/kernels_jax.py tests/test_cork_lw_kernels_jax_parity.py
git commit -m "feat(cork): differentiable JAX Planck sources (Phase 2)"
```

### Task 4: `lw_transport_jax` — `lax.scan` transport solver

**Files:**
- Modify: `climt/_components/cork/lw/kernels_jax.py`
- Test: `tests/test_cork_lw_kernels_jax_parity.py` (add)

**Interfaces:**
- Produces: `lw_transport_jax(tau, planck_source, surface_source, emissivity, weights) -> (up_band, down_band, up_broad, down_broad)` with `tau,planck_source: (nband,ngpt,nlev,ncol)`, `surface_source: (nband,ngpt,ncol)`, `emissivity: (nband,ncol)`, `weights: (nband,ngpt)`; returns `up_band,down_band: (nband,nlev+1,ncol)`, `up_broad,down_broad: (nlev+1,ncol)`. Forward-only (no diagnostics).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_cork_lw_kernels_jax_parity.py
def test_lw_transport_jax_matches_oracle():
    from climt._components.cork.lw.kernels import lw_transport
    from climt._components.cork.lw.kernels_jax import lw_transport_jax
    rng = np.random.RandomState(2)
    nband, ngpt, nlev, ncol = 4, 3, 8, 5
    tau = rng.uniform(0.0, 2.0, size=(nband, ngpt, nlev, ncol))
    planck = rng.uniform(50.0, 400.0, size=(nband, ngpt, nlev, ncol))
    surf = rng.uniform(200.0, 450.0, size=(nband, ngpt, ncol))
    emis = rng.uniform(0.8, 1.0, size=(nband, ncol))
    weights = np.full((nband, ngpt), 0.5)

    up_b, dn_b, up_br, dn_br = lw_transport(
        None, None, tau, planck, surf, emis, weights, 5.67e-8, diagnostics_level=0)
    upj, dnj, ubrj, dbrj = lw_transport_jax(
        jnp.asarray(tau), jnp.asarray(planck), jnp.asarray(surf),
        jnp.asarray(emis), jnp.asarray(weights))
    np.testing.assert_allclose(np.asarray(upj), up_b, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(np.asarray(dnj), dn_b, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(np.asarray(ubrj), up_br, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(np.asarray(dbrj), dn_br, rtol=1e-7, atol=1e-9)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cork_lw_kernels_jax_parity.py::test_lw_transport_jax_matches_oracle -v`
Expected: FAIL — `ImportError: cannot import name 'lw_transport_jax'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to climt/_components/cork/lw/kernels_jax.py
def lw_transport_jax(tau, planck_source, surface_source, emissivity, weights):
    nband, ngpt, nlev, ncol = tau.shape
    D = DIFFUSIVITY_FACTOR
    trans = jnp.exp(-D * tau)                       # (nband,ngpt,nlev,ncol)
    trans_lev = jnp.moveaxis(trans, 2, 0)           # (nlev,nband,ngpt,ncol)
    planck_lev = jnp.moveaxis(planck_source, 2, 0)  # (nlev,nband,ngpt,ncol)

    def step(prev, layer):
        tr, pl = layer
        cur = prev * tr + pl * (1.0 - tr)
        return cur, cur

    # Upward sweep: surface -> TOA. BC = emissivity * surface_source.
    up0 = emissivity[:, None, :] * surface_source   # (nband,ngpt,ncol)
    _, up_stack = lax.scan(step, up0, (trans_lev, planck_lev))   # (nlev,nband,ngpt,ncol)
    up_iface = jnp.concatenate([up0[None], up_stack], axis=0)    # (nlev+1,...)

    # Downward sweep: TOA -> surface. BC = 0 at TOA.
    dn0 = jnp.zeros((nband, ngpt, ncol))
    _, dn_stack = lax.scan(step, dn0, (trans_lev, planck_lev), reverse=True)
    down_iface = jnp.concatenate([dn_stack, dn0[None]], axis=0)  # (nlev+1,...)

    up_g = jnp.moveaxis(up_iface, 0, 2)             # (nband,ngpt,nlev+1,ncol)
    down_g = jnp.moveaxis(down_iface, 0, 2)
    w = weights[:, :, None, None]
    up_band = jnp.sum(w * up_g, axis=1)             # (nband,nlev+1,ncol)
    down_band = jnp.sum(w * down_g, axis=1)
    up_broad = jnp.sum(up_band, axis=0)             # (nlev+1,ncol)
    down_broad = jnp.sum(down_band, axis=0)
    return up_band, down_band, up_broad, down_broad
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_lw_kernels_jax_parity.py -v`
Expected: PASS (all tests).

- [ ] **Step 5: Commit**

```bash
git add climt/_components/cork/lw/kernels_jax.py tests/test_cork_lw_kernels_jax_parity.py
git commit -m "feat(cork): differentiable lax.scan LW transport solver (Phase 2)"
```

### Task 5: Gradient check through Planck + transport

**Files:**
- Test: `tests/test_cork_lw_kernels_jax_parity.py` (add)

**Interfaces:**
- Consumes: `planck_sources_jax`, `lw_transport_jax`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_cork_lw_kernels_jax_parity.py
def test_transport_grad_through_temperature_is_finite_and_correct():
    from climt._components.cork.lw.kernels_jax import (
        planck_sources_jax, lw_transport_jax)
    rng = np.random.RandomState(5)
    nband, ngpt, nT, nlev, ncol = 3, 2, 6, 6, 2
    planck_frac = jnp.asarray(rng.uniform(0.0, 1.0, size=(nband, ngpt, nT)))
    T_grid = jnp.asarray(np.linspace(180.0, 330.0, nT))
    tau = jnp.asarray(rng.uniform(0.0, 1.5, size=(nband, ngpt, nlev, ncol)))
    emis = jnp.asarray(rng.uniform(0.85, 1.0, size=(nband, ncol)))
    weights = jnp.asarray(np.full((nband, ngpt), 0.5))
    T_surf = jnp.asarray(rng.uniform(270.0, 300.0, size=(ncol,)))
    sigma = 5.670374419e-8
    T0 = jnp.asarray(rng.uniform(210.0, 310.0, size=(nlev, ncol)))

    def olr(T):
        ps, ss = planck_sources_jax(planck_frac, T_grid, T, T_surf, sigma)
        _, _, up_broad, _ = lw_transport_jax(tau, ps, ss, emis, weights)
        return up_broad[-1].sum()   # outgoing longwave at TOA

    g = np.asarray(jax.grad(olr)(T0))
    assert np.all(np.isfinite(g))
    eps = 1e-2
    i, j = 3, 1
    fd = (float(olr(T0.at[i, j].add(eps))) - float(olr(T0.at[i, j].add(-eps)))) / (2 * eps)
    np.testing.assert_allclose(g[i, j], fd, rtol=1e-3, atol=1e-6)
```

- [ ] **Step 2: Run test to verify it fails or passes**

Run: `pytest tests/test_cork_lw_kernels_jax_parity.py::test_transport_grad_through_temperature_is_finite_and_correct -v`
Expected: PASS (kernels from Tasks 3–4 are already differentiable). If `nan`, check that `_bracket` clips and no zero-division enters.

- [ ] **Step 3: (only if failing) fix**

No new code expected; resolve any `nan` by confirming guard clips from Tasks 1 & 3.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_lw_kernels_jax_parity.py -v`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add tests/test_cork_lw_kernels_jax_parity.py
git commit -m "test(cork): gradient check through JAX Planck + transport"
```

---

## Phase 3 — Component integration (`lw/component.py`)

### Task 6: `CorkLWTable` container + `cork_lw_jax` orchestrator

**Files:**
- Modify: `climt/_components/cork/lw/kernels_jax.py`
- Test: `tests/test_cork_lw_jax_integration.py`

**Interfaces:**
- Produces:
  - `CorkLWTable` — a `typing.NamedTuple` of static `jnp` arrays + python flags: `k, T_grid, p_grid_log, log_x_grid, log_c_grid, planck_frac, gpoint_weights, has_cont, log_cont, co2_logk`.
  - `build_cork_lw_table(table) -> CorkLWTable` — builds the container from a loaded numpy k-table dict (host-side, once).
  - `cork_lw_jax(T, p, p_int, T_surf, h2o, co2, emissivity, tau_cloud, jtable, sigma, g, cpd) -> (tendency, up_broad, down_broad, up_band, down_band, tau)`. Shapes: `T,p: (nlev,ncol)`; `p_int: (nlev+1,ncol)`; `T_surf,h2o*?`: see body; `h2o,co2: (nlev,ncol)`; `emissivity: (nband,ncol)`; `tau_cloud: (nlev,ncol,nband)`; returns `tendency: (nlev,ncol)` (K/s), fluxes as in `lw_transport_jax`, `tau: (nband,ngpt,nlev,ncol)`.
- Consumes: `compute_tau_jax`, `planck_sources_jax`, `lw_transport_jax`.

This orchestrator replicates the `correlated_k` + premixed-bg branch of `component.array_call` (lines 234–305) in `jnp`: build `gas_amounts[0] = column air mass`, derive `h2o_vmr` from specific humidity and `co2_vmr` from state, prepare `log_p/log_x/log_c`, compute τ, add cloud τ, Planck sources, transport, heating rate.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cork_lw_jax_integration.py
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_platform_name", "cpu")
import jax.numpy as jnp
import numpy as np
import pytest
from climt._components.cork.optics.correlated_k import load_k_table


def test_cork_lw_table_roundtrip_shapes():
    from climt._components.cork.lw.kernels_jax import build_cork_lw_table
    table = load_k_table("earth_low_res_lw")
    jt = build_cork_lw_table(table)
    assert jt.k.shape == table["k_coefficients"].shape
    assert jt.planck_frac.shape == np.asarray(table["planck_fraction"]).shape
    assert jt.gpoint_weights.shape == np.asarray(table["gpoint_weights"]).shape
    assert jt.has_cont is True and jt.co2_logk is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cork_lw_jax_integration.py::test_cork_lw_table_roundtrip_shapes -v`
Expected: FAIL — `ImportError: cannot import name 'build_cork_lw_table'`.

- [ ] **Step 3: Write minimal implementation**

```python
# append to climt/_components/cork/lw/kernels_jax.py
from typing import NamedTuple
import numpy as np

from ..optics.correlated_k_jax import compute_tau_jax


class CorkLWTable(NamedTuple):
    k: jnp.ndarray
    T_grid: jnp.ndarray
    p_grid_log: jnp.ndarray
    log_x_grid: jnp.ndarray
    log_c_grid: jnp.ndarray
    planck_frac: jnp.ndarray
    gpoint_weights: jnp.ndarray
    has_cont: bool
    log_cont: jnp.ndarray
    co2_logk: bool


def build_cork_lw_table(table):
    """Host-side: pack a loaded numpy k-table dict into static jnp arrays.

    Only valid for the 7-D additive-CO2 premixed-background case.
    """
    k = np.asarray(table["k_coefficients"])
    assert k.ndim == 7, "JAX LW path requires a 7-D additive-CO2 k-table"
    x_grid = np.asarray(table["h2o_vmr_grid"], dtype=np.float64)
    c_grid = np.asarray(table["co2_vmr_grid"], dtype=np.float64)
    has_cont = ("continuum_kappa" in table
                and np.asarray(table["continuum_kappa"]).ndim == 4)
    if has_cont:
        log_cont = np.log(np.maximum(
            np.asarray(table["continuum_kappa"], dtype=np.float64), 1e-40))
    else:
        log_cont = np.zeros((k.shape[1], 1, 1, 1))
    return CorkLWTable(
        k=jnp.asarray(k, dtype=jnp.float64),
        T_grid=jnp.asarray(table["temperature_grid"], dtype=jnp.float64),
        p_grid_log=jnp.asarray(table["pressure_grid_log"], dtype=jnp.float64),
        log_x_grid=jnp.log(jnp.maximum(jnp.asarray(x_grid), 1e-30)),
        log_c_grid=jnp.log(jnp.maximum(jnp.asarray(c_grid), 1e-30)),
        planck_frac=jnp.asarray(table["planck_fraction"], dtype=jnp.float64),
        gpoint_weights=jnp.asarray(table["gpoint_weights"], dtype=jnp.float64),
        has_cont=bool(has_cont),
        log_cont=jnp.asarray(log_cont),
        co2_logk=True,
    )


# Mass ratio M_H2O / M_dry_air, for specific humidity -> H2O VMR conversion.
from ..common import MOLAR_MASS, MOLAR_MASS_DRY_AIR
_M_H2O_RATIO = MOLAR_MASS["h2o"] / MOLAR_MASS_DRY_AIR


def _column_amount_jax(q, p_int, g):
    dp = jnp.abs(p_int[1:] - p_int[:-1])
    return q * dp / g


def _heating_rate_jax(net_flux, p_int, g, cpd):
    dp = p_int[1:] - p_int[:-1]
    dflux = net_flux[1:] - net_flux[:-1]
    return g / cpd * dflux / dp


def cork_lw_jax(T, p, p_int, T_surf, h2o, co2, emissivity, tau_cloud,
                jtable, sigma, g, cpd):
    nlev, ncol = T.shape
    # gas amount is column air mass (premixed-bg); H2O enters via the X axis.
    air = _column_amount_jax(jnp.ones((nlev, ncol)), p_int, g)
    gas_amounts = air[None, :, :]                       # (1,nlev,ncol)

    h2o_vmr = h2o / jnp.maximum(h2o + (1.0 - h2o) * _M_H2O_RATIO, 1e-30)
    x_lo = jnp.exp(jtable.log_x_grid[0]); x_hi = jnp.exp(jtable.log_x_grid[-1])
    c_lo = jnp.exp(jtable.log_c_grid[0]); c_hi = jnp.exp(jtable.log_c_grid[-1])
    log_p = jnp.log(jnp.maximum(p, 1.0))
    log_x = jnp.log(jnp.maximum(jnp.clip(h2o_vmr, x_lo, x_hi), 1e-30))
    log_c = jnp.log(jnp.maximum(jnp.clip(co2, c_lo, c_hi), 1e-30))

    tau = compute_tau_jax(
        T, log_p, log_x, log_c, gas_amounts,
        jtable.k, jtable.T_grid, jtable.p_grid_log,
        jtable.log_x_grid, jtable.log_c_grid,
        jtable.has_cont, jtable.log_cont, jtable.co2_logk)

    # cloud optical depth: (nlev,ncol,nband) -> (nband,1,nlev,ncol)
    tau = tau + jnp.transpose(tau_cloud, (2, 0, 1))[:, None, :, :]

    planck_src, surf_src = planck_sources_jax(
        jtable.planck_frac, jtable.T_grid, T, T_surf, sigma)

    up_band, down_band, up_broad, down_broad = lw_transport_jax(
        tau, planck_src, surf_src, emissivity, jtable.gpoint_weights)

    net_flux = up_broad - down_broad
    tendency = _heating_rate_jax(net_flux, p_int, g, cpd)
    return tendency, up_broad, down_broad, up_band, down_band, tau
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_lw_jax_integration.py::test_cork_lw_table_roundtrip_shapes -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_components/cork/lw/kernels_jax.py tests/test_cork_lw_jax_integration.py
git commit -m "feat(cork): cork_lw_jax orchestrator + static table container (Phase 3)"
```

### Task 7: Dispatch JAX-first in `array_call` + full-component parity

**Files:**
- Modify: `climt/_components/cork/lw/component.py`
- Test: `tests/test_cork_lw_jax_integration.py` (add)

**Interfaces:**
- Consumes: `build_cork_lw_table`, `cork_lw_jax`.
- Produces: `CorkLongwaveRadiation.array_call` returns identical `({"T": tendency}, diagnostics)` structure whether inputs are numpy or `jax.Array`; JAX path active only for correlated-k + 7-D premixed-bg + `diagnostics_level == 0`.

The dispatch keeps the numpy oracle path byte-for-byte. It adds, at the **top** of `array_call`, a guard that routes to a new `_array_call_jax` method when `HAS_JAX and isinstance(state["T"], jax.Array)`.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_cork_lw_jax_integration.py
def _make_state(nlev, ncol, as_jax):
    from climt import get_grid, get_default_state
    from climt._components.cork import CorkLongwaveRadiation
    from climt._core.util import numpy_version_of
    comp = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    state = get_default_state([comp], grid_state=grid)
    raw = numpy_version_of(state)
    nband = comp.num_longwave_bands
    xp = jnp if as_jax else np

    def arr(name, shape):
        return xp.asarray(raw[name].reshape(shape))

    s = {
        "T": arr("air_temperature", (nlev, ncol)),
        "p": arr("air_pressure", (nlev, ncol)),
        "p_int": arr("air_pressure_on_interface_levels", (nlev + 1, ncol)),
        "T_surf": arr("surface_temperature", (ncol,)),
        "emissivity": arr("surface_longwave_emissivity", (nband, ncol)),
        "h2o": arr("specific_humidity", (nlev, ncol)),
        "co2": arr("mole_fraction_of_carbon_dioxide_in_air", (nlev, ncol)),
        "tau_cloud_lw": arr("longwave_optical_thickness_due_to_cloud", (nlev, ncol, nband)),
    }
    return comp, s


def test_array_call_jax_matches_numpy():
    nlev, ncol = 20, 4
    comp, s_np = _make_state(nlev, ncol, as_jax=False)
    _, s_jx = _make_state(nlev, ncol, as_jax=True)
    tend_np, diag_np = comp.array_call(s_np)
    tend_jx, diag_jx = comp.array_call(s_jx)
    np.testing.assert_allclose(np.asarray(tend_jx["T"]), tend_np["T"], rtol=1e-6, atol=1e-9)
    for key in ("upwelling_longwave_flux_in_air", "downwelling_longwave_flux_in_air"):
        np.testing.assert_allclose(
            np.asarray(diag_jx[key]), diag_np[key], rtol=1e-6, atol=1e-6)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cork_lw_jax_integration.py::test_array_call_jax_matches_numpy -v`
Expected: FAIL — JAX inputs hit the numpy code path and raise (numpy ops on `jax.Array`) or produce a type error; no JAX dispatch exists yet.

- [ ] **Step 3: Write minimal implementation**

Add imports near the top of `climt/_components/cork/lw/component.py`:

```python
try:
    import jax
    import jax.numpy as jnp
    HAS_JAX = True
except ImportError:
    HAS_JAX = False
```

In `__init__`, after the `correlated_k` table is loaded (after line 51), build and cache the JAX table when the case is supported:

```python
        self._jax_table = None
        if optics == "correlated_k" and HAS_JAX:
            try:
                if self._table["k_coefficients"].ndim == 7 and self._premixed_bg:
                    from .kernels_jax import build_cork_lw_table
                    self._jax_table = build_cork_lw_table(self._table)
            except Exception:
                self._jax_table = None
```

At the very top of `array_call` (before `T = state["T"]`), add the dispatch:

```python
        if (HAS_JAX and self._jax_table is not None
                and self._diagnostics_level == 0
                and isinstance(state["T"], jax.Array)):
            return self._array_call_jax(state)
```

Add the new method (place after `array_call`):

```python
    def _array_call_jax(self, state):
        from .kernels_jax import cork_lw_jax
        sigma = get_constant("stefan_boltzmann_constant", "W/m^2/K^4")
        g = get_constant("gravitational_acceleration", "m/s^2")
        cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")

        T = state["T"]; orig_shape_T = T.shape
        p_int = state["p_int"]; orig_shape_pint = p_int.shape
        nlev = T.shape[0]
        T_flat = T.reshape(nlev, -1)
        ncol = T_flat.shape[1]
        p_flat = state["p"].reshape(nlev, -1)
        p_int_flat = p_int.reshape(nlev + 1, -1)
        T_surf_flat = state["T_surf"].reshape(-1)
        nband = self._num_bands
        emissivity = state["emissivity"].reshape(nband, ncol)
        h2o = state["h2o"].reshape(nlev, ncol)
        co2 = state["co2"].reshape(nlev, ncol)
        tau_cloud = state["tau_cloud_lw"].reshape(nlev, ncol, nband)

        tendency, up_broad, down_broad, up_band, down_band, tau = cork_lw_jax(
            T_flat, p_flat, p_int_flat, T_surf_flat, h2o, co2,
            emissivity, tau_cloud, self._jax_table, sigma, g, cpd)

        ngpt = tau.shape[1]
        weights = self._jax_table.gpoint_weights
        D = 1.66
        tau_band = jnp.einsum("bgnc,bg->bnc", tau, weights)          # (nband,nlev,ncol)
        trans_band = jnp.exp(-D * tau_band)
        hr_band = jnp.stack([
            self._hr_band_jax(up_band[b] - down_band[b], p_int_flat, g, cpd)
            for b in range(nband)], axis=0) * 86400.0

        def lay(a):
            return jnp.moveaxis(a, 0, -1).reshape(orig_shape_T + (nband,))

        def ifc(a):
            return jnp.moveaxis(a, 0, -1).reshape(orig_shape_pint + (nband,))

        tendency_out = tendency.reshape(orig_shape_T)
        diagnostics = {
            "upwelling_longwave_flux_in_air": up_broad.reshape(orig_shape_pint),
            "downwelling_longwave_flux_in_air": down_broad.reshape(orig_shape_pint),
            "upwelling_longwave_flux_in_air_per_band": ifc(up_band),
            "downwelling_longwave_flux_in_air_per_band": ifc(down_band),
            "air_temperature_tendency_from_longwave": tendency_out * 86400.0,
            "longwave_optical_depth_per_band": lay(tau_band),
            "longwave_transmittance_per_band": lay(trans_band),
            "air_temperature_tendency_from_longwave_per_band": lay(hr_band),
        }
        return ({"T": tendency_out}, diagnostics)

    @staticmethod
    def _hr_band_jax(net_flux, p_int, g, cpd):
        dp = p_int[1:] - p_int[:-1]
        dflux = net_flux[1:] - net_flux[:-1]
        return g / cpd * dflux / dp
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cork_lw_jax_integration.py::test_array_call_jax_matches_numpy -v`
Expected: PASS.

- [ ] **Step 5: Run the full CORK LW suite to confirm the numpy oracle is untouched**

Run: `pytest tests/test_cork_lw.py tests/test_cork_ck_njit.py -q`
Expected: PASS (no regressions in the numpy path).

- [ ] **Step 6: Commit**

```bash
git add climt/_components/cork/lw/component.py tests/test_cork_lw_jax_integration.py
git commit -m "feat(cork): JAX-first dispatch in CorkLongwaveRadiation.array_call (Phase 3)"
```

### Task 8: jit/vmap-safety + end-to-end `jax.grad`

**Files:**
- Test: `tests/test_cork_lw_jax_integration.py` (add)

**Interfaces:**
- Consumes: `cork_lw_jax`, `build_cork_lw_table`, `_make_state` (from Task 7).

- [ ] **Step 1: Write the failing test**

```python
# append to tests/test_cork_lw_jax_integration.py
def _core_args(nlev, ncol):
    from climt._components.cork.lw.kernels_jax import build_cork_lw_table
    comp, s = _make_state(nlev, ncol, as_jax=True)
    jt = comp._jax_table
    args = dict(
        T=s["T"], p=s["p"], p_int=s["p_int"], T_surf=s["T_surf"],
        h2o=s["h2o"], co2=s["co2"], emissivity=s["emissivity"],
        tau_cloud=s["tau_cloud_lw"], sigma=5.670374419e-8, g=9.80665, cpd=1004.64)
    return jt, args


def test_cork_lw_jax_jit_matches_eager():
    import jax
    from climt._components.cork.lw.kernels_jax import cork_lw_jax
    jt, args = _core_args(20, 4)
    eager = cork_lw_jax(jtable=jt, **args)[0]
    jitted = jax.jit(lambda **kw: cork_lw_jax(jtable=jt, **kw))(**args)[0]
    np.testing.assert_allclose(np.asarray(jitted), np.asarray(eager), rtol=1e-9, atol=1e-12)


def test_cork_lw_jax_vmap_over_columns_matches_loop():
    import jax
    from climt._components.cork.lw.kernels_jax import cork_lw_jax
    jt, args = _core_args(20, 3)

    def one_col(T, p, p_int, T_surf, h2o, co2, emis, tau_cloud):
        # promote a single column (…,1) and take tendency
        t, *_ = cork_lw_jax(T[:, None], p[:, None], p_int[:, None], T_surf[None],
                            h2o[:, None], co2[:, None], emis[:, None],
                            tau_cloud[:, None, :], jt, args["sigma"], args["g"], args["cpd"])
        return t[:, 0]

    vmapped = jax.vmap(one_col, in_axes=(1, 1, 1, 0, 1, 1, 1, 1))(
        args["T"], args["p"], args["p_int"], args["T_surf"],
        args["h2o"], args["co2"], args["emissivity"], args["tau_cloud"])
    full = cork_lw_jax(jtable=jt, **args)[0]
    np.testing.assert_allclose(np.asarray(vmapped.T), np.asarray(full), rtol=1e-7, atol=1e-9)


def test_cork_lw_jax_grad_wrt_state_is_finite():
    import jax
    from climt._components.cork.lw.kernels_jax import cork_lw_jax
    jt, args = _core_args(20, 2)

    def olr(T, T_surf):
        kw = dict(args); kw.pop("T"); kw.pop("T_surf")
        _, up_broad, *_ = cork_lw_jax(T=T, T_surf=T_surf, jtable=jt, **kw)
        return up_broad[-1].sum()

    gT, gTs = jax.grad(olr, argnums=(0, 1))(args["T"], args["T_surf"])
    assert np.all(np.isfinite(np.asarray(gT)))
    assert np.all(np.isfinite(np.asarray(gTs)))
    # OLR increases with surface temperature
    assert float(np.asarray(gTs).sum()) > 0.0
```

- [ ] **Step 2: Run tests to verify they fail or pass**

Run: `pytest tests/test_cork_lw_jax_integration.py -k "jit or vmap or grad" -v`
Expected: PASS if the kernels are transform-clean. A FAIL here points to a concrete defect — e.g. a python-`if` on a traced value (replace with `jnp.where`), or a shape that assumed a fixed `ncol` (must be generic). Fix in the relevant `*_jax.py` file and re-run.

- [ ] **Step 3: (only if failing) fix transform issues**

Typical fixes: remove data-dependent python branching inside the kernels; ensure all internal shapes derive from input shapes, not captured constants. No fix expected if Phases 1–2 followed the plan.

- [ ] **Step 4: Run the whole new JAX suite**

Run: `pytest tests/test_cork_optics_jax_parity.py tests/test_cork_lw_kernels_jax_parity.py tests/test_cork_lw_jax_integration.py -q`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add tests/test_cork_lw_jax_integration.py
git commit -m "test(cork): jit/vmap-safety and end-to-end jax.grad for CORK LW"
```

---

## Definition of Done

- `CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")` dispatches to the JAX path on `jax.Array` input and matches the numpy oracle within `rtol=1e-6` for tendency and broadband fluxes.
- The full LW path is `jit`/`vmap`(over columns)/`grad`-compatible; state gradients are finite and finite-difference-correct.
- Numba kernels (`correlated_k.py`, `lw/kernels.py`) are unmodified; `tests/test_cork_lw.py` and `tests/test_cork_ck_njit.py` still pass.
- New modules set `jax_enable_x64` at import.

## Out of Scope (future plans)

- CORK shortwave JAX port (own spec/plan).
- Parmentier-optics JAX path; ESFT overlap; non-premixed multi-gas tables; `diagnostics_level > 0` in the JAX path.
- Gradients w.r.t. k-table / continuum / Parmentier parameters.
- Boundary-layer scheme and the aquaplanet driver.
