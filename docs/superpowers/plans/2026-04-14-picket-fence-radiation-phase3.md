# Picket-Fence Radiation Phase 3: Pedagogical Features, Clouds, Multi-Planet, and Validation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the picket-fence radiation scheme by adding pedagogical introspection (diagnostics levels, component-level per-band diagnostics), SW cloud inputs, multi-planet atmospheric properties (TOML profiles), and validation against analytical solutions and RRTMG.

**Architecture:** Phase 1 delivered Parmentier + correlated-k LW/SW with direct beam. Phase 2 added the full Meador & Weaver two-stream SW kernel, ESFT overlap, correlated-k SW, and stellar spectrum loading. Phase 3 adds the remaining design spec features: optional intermediate diagnostics inside the solver kernels (`diagnostics_level`), component-level per-band optical depth/transmittance/heating rate diagnostics for pedagogical use, SW cloud scattering inputs, a `load_atmospheric_properties` API backed by TOML profiles, and validation notebooks/tests comparing against analytical solutions and RRTMG.

**Tech Stack:** Python 3.12, NumPy, numba, sympl, tomli (for TOML parsing), pytest

---

## File Structure

```
climt/
  _components/
    picket_fence/
      sw/
        kernels.py             # MODIFY: add diagnostics_level to sw_two_stream
        component.py           # MODIFY: add SW cloud inputs, pass diagnostics_level, add per-band diagnostics
      lw/
        kernels.py             # MODIFY: add diagnostics_level to lw_transport
        component.py           # MODIFY: pass diagnostics_level, add per-band diagnostics
      common.py                # no changes
  _core/
    atmospheric_properties.py  # CREATE: load_atmospheric_properties, reset_atmospheric_properties
    constants.py               # no changes (uses set_constants_from_dict)
  _data/
    atmospheric_properties/
      earth.toml               # CREATE: Earth constants
      mars.toml                # CREATE: Mars constants
      hot_jupiter.toml         # CREATE: hot Jupiter constants
  __init__.py                  # MODIFY: export atmospheric_properties API
tests/
  test_picket_fence_kernels.py # MODIFY: add diagnostics_level tests
  test_picket_fence_sw.py      # MODIFY: add SW cloud tests
  test_atmospheric_properties.py # CREATE: TOML loading tests
  test_picket_fence_validation.py # CREATE: validation against analytical solutions
```

**Notes:**
- The `diagnostics_level` parameter is added to the existing kernel functions but defaults to 0 (no overhead). When non-zero, the kernels return an additional dictionary of intermediate arrays.
- Component-level per-band diagnostics (optical depth, transmittance, heating rate) are always computed and returned — they add negligible cost since the tau arrays are already available inside `array_call`.
- SW cloud inputs follow the same pattern as the LW cloud optical depth already in the LW component.
- The atmospheric properties API wraps sympl's `set_constant`/`get_constant` with TOML file loading and snapshot/restore.
- Validation tests are separated into their own file since they test physics accuracy, not code correctness.

---

### Task 1: SW kernel diagnostics — `diagnostics_level` parameter

**Files:**
- Modify: `climt/_components/picket_fence/sw/kernels.py`
- Modify: `tests/test_picket_fence_kernels.py`

Add an optional `diagnostics_level` parameter to `sw_two_stream` that returns intermediate solver quantities. When `diagnostics_level=0` (default), behavior is unchanged. When non-zero, the function returns a 5th element: a dictionary of arrays.

Because numba `@njit` functions cannot return Python dicts, we implement this by having `sw_two_stream` remain a plain Python function (not jitted) that calls the jitted inner functions. The current `sw_two_stream` is already structured this way — the `@njit` decorator on it works but the function orchestrates calls to `_delta_scale`, `_sw_dif_and_source`, and `_adding`. We'll split: the core loop stays jitted, and a thin wrapper handles the diagnostics dict.

- [x] **Step 1: Write failing tests for diagnostics_level**

Add to `tests/test_picket_fence_kernels.py`:

```python
def test_sw_two_stream_diagnostics_level0_unchanged():
    """diagnostics_level=0 returns the same 4-tuple as before."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev, ncol, nband, ngpt = 5, 1, 1, 1
    tau = 0.3 * np.ones((nband, ngpt, nlev, ncol))
    ssa = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    asymmetry = 0.1 * np.ones((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 4])
    albedo = np.array([0.3])
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    result = sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights,
                           diagnostics_level=0)
    assert len(result) == 4, "diagnostics_level=0 should return 4-tuple"


def test_sw_two_stream_diagnostics_level1():
    """diagnostics_level=1 returns per-layer diffuse R/T and direct beam profile."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev, ncol, nband, ngpt = 5, 1, 1, 1
    tau = 0.3 * np.ones((nband, ngpt, nlev, ncol))
    ssa = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    asymmetry = 0.1 * np.ones((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 4])
    albedo = np.array([0.3])
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    result = sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights,
                           diagnostics_level=1)
    assert len(result) == 5, "diagnostics_level=1 should return 5-tuple"
    up_band, down_band, up_broad, down_broad, diag = result

    # Level 1 diagnostics keys
    assert "Rdif" in diag
    assert "Tdif" in diag
    assert "Tnoscat" in diag
    assert "direct_beam_flux" in diag

    # Shapes: (nband, ngpt, nlev, ncol) for per-layer, (nband, ngpt, nlev+1, ncol) for interfaces
    assert diag["Rdif"].shape == (nband, ngpt, nlev, ncol)
    assert diag["Tdif"].shape == (nband, ngpt, nlev, ncol)
    assert diag["Tnoscat"].shape == (nband, ngpt, nlev, ncol)
    assert diag["direct_beam_flux"].shape == (nband, ngpt, nlev + 1, ncol)

    # Physical sanity: diffuse reflectance in [0, 1]
    assert np.all(diag["Rdif"] >= 0)
    assert np.all(diag["Rdif"] <= 1)


def test_sw_two_stream_diagnostics_level2():
    """diagnostics_level=2 additionally returns Rdir, Tdir, combined albedo, delta-scaled properties."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev, ncol, nband, ngpt = 5, 1, 1, 1
    tau = 0.3 * np.ones((nband, ngpt, nlev, ncol))
    ssa = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    asymmetry = 0.1 * np.ones((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 4])
    albedo = np.array([0.3])
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    result = sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights,
                           diagnostics_level=2)
    assert len(result) == 5
    _, _, _, _, diag = result

    # Level 2 adds these keys on top of Level 1
    assert "Rdir" in diag
    assert "Tdir" in diag
    assert "combined_albedo" in diag
    assert "tau_delta" in diag
    assert "ssa_delta" in diag
    assert "g_delta" in diag

    assert diag["Rdir"].shape == (nband, ngpt, nlev, ncol)
    assert diag["combined_albedo"].shape == (nband, ngpt, nlev + 1, ncol)
    assert diag["tau_delta"].shape == (nband, ngpt, nlev, ncol)
```

- [x] **Step 2: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_kernels.py::test_sw_two_stream_diagnostics_level0_unchanged tests/test_picket_fence_kernels.py::test_sw_two_stream_diagnostics_level1 tests/test_picket_fence_kernels.py::test_sw_two_stream_diagnostics_level2 -v`
Expected: FAIL (sw_two_stream does not accept diagnostics_level)

- [x] **Step 3: Implement diagnostics_level in sw_two_stream**

The strategy: rename the current jitted `sw_two_stream` to `_sw_two_stream_core` (keep `@njit`), and create a new plain-Python `sw_two_stream` wrapper that calls the core and optionally collects diagnostics.

In `climt/_components/picket_fence/sw/kernels.py`, rename the existing `sw_two_stream` function to `_sw_two_stream_core` (keep the `@njit` decorator and all existing code unchanged). Then add a new `sw_two_stream` wrapper:

```python
def sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights,
                  diagnostics_level=0):
    """Multi-band, multi-g-point SW radiative transfer.

    Meador & Weaver (1980) two-stream with delta-Eddington scaling and
    adding method for inter-layer coupling. Follows RTE/RRTMGP.

    Args:
        tau: (nband, ngpt, nlev, ncol) extinction optical depth
        ssa: (nband, ngpt, nlev, ncol) single scattering albedo
        asymmetry: (nband, ngpt, nlev, ncol) asymmetry parameter
        zenith: (ncol,) solar zenith angle, radians
        albedo: (ncol,) surface albedo
        solar_flux: (nband, ngpt) TOA solar flux per g-point, W/m^2
        weights: (nband, ngpt) g-point weights
        diagnostics_level: 0 (fluxes only), 1 (per-layer R/T + direct beam),
            2 (additionally delta-scaled properties + combined albedo)

    Returns:
        If diagnostics_level == 0:
            (up_band, down_band, up_broad, down_broad)
        If diagnostics_level > 0:
            (up_band, down_band, up_broad, down_broad, diagnostics_dict)
    """
    if diagnostics_level == 0:
        return _sw_two_stream_core(tau, ssa, asymmetry, zenith, albedo,
                                   solar_flux, weights)

    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    # Diagnostic arrays
    diag_Rdif = np.zeros((nband, ngpt, nlev, ncol))
    diag_Tdif = np.zeros((nband, ngpt, nlev, ncol))
    diag_Tnoscat = np.zeros((nband, ngpt, nlev, ncol))
    diag_direct = np.zeros((nband, ngpt, nlev + 1, ncol))

    if diagnostics_level >= 2:
        diag_Rdir = np.zeros((nband, ngpt, nlev, ncol))
        diag_Tdir = np.zeros((nband, ngpt, nlev, ncol))
        diag_combined_albedo = np.zeros((nband, ngpt, nlev + 1, ncol))
        diag_tau_d = np.zeros((nband, ngpt, nlev, ncol))
        diag_ssa_d = np.zeros((nband, ngpt, nlev, ncol))
        diag_g_d = np.zeros((nband, ngpt, nlev, ncol))

    for b in range(nband):
        for g in range(ngpt):
            w = weights[b, g]
            for i in range(ncol):
                mu0 = np.cos(zenith[i])
                if mu0 <= 1e-4:
                    continue

                Rdif = np.zeros(nlev)
                Tdif = np.zeros(nlev)
                Rdir = np.zeros(nlev)
                Tdir = np.zeros(nlev)
                Tnoscat = np.zeros(nlev)
                flux_dn_dir = np.zeros(nlev + 1)
                flux_dn_dir[nlev] = 1.0
                src_up = np.zeros(nlev)
                src_dn = np.zeros(nlev)

                for k in range(nlev - 1, -1, -1):
                    tau_s, ssa_s, g_s = _delta_scale(
                        tau[b, g, k, i], ssa[b, g, k, i], asymmetry[b, g, k, i]
                    )
                    Rdif[k], Tdif[k], Rdir[k], Tdir[k], Tnoscat[k] = (
                        _sw_dif_and_source(tau_s, ssa_s, g_s, mu0)
                    )
                    flux_dn_dir[k] = Tnoscat[k] * flux_dn_dir[k + 1]
                    src_up[k] = Rdir[k] * flux_dn_dir[k + 1]
                    src_dn[k] = Tdir[k] * flux_dn_dir[k + 1]

                    if diagnostics_level >= 2:
                        diag_tau_d[b, g, k, i] = tau_s
                        diag_ssa_d[b, g, k, i] = ssa_s
                        diag_g_d[b, g, k, i] = g_s

                src_sfc = flux_dn_dir[0] * albedo[i]
                flux_up_dif, flux_dn_dif = _adding(
                    nlev, albedo[i], Rdif, Tdif, src_up, src_dn, src_sfc,
                    flux_dn_dir
                )

                # Store diagnostics (last g-point/column wins — these are per-gpt)
                diag_Rdif[b, g, :, i] = Rdif
                diag_Tdif[b, g, :, i] = Tdif
                diag_Tnoscat[b, g, :, i] = Tnoscat
                diag_direct[b, g, :, i] = flux_dn_dir * solar_flux[b, g] * mu0

                if diagnostics_level >= 2:
                    diag_Rdir[b, g, :, i] = Rdir
                    diag_Tdir[b, g, :, i] = Tdir
                    # Reconstruct combined albedo from the adding method
                    comb_alb = np.zeros(nlev + 1)
                    comb_alb[0] = albedo[i]
                    for k in range(nlev):
                        d = 1.0 - Rdif[k] * comb_alb[k]
                        if abs(d) < 1e-30:
                            d = 1e-30
                        comb_alb[k + 1] = Rdif[k] + Tdif[k] * Tdif[k] * comb_alb[k] / d
                    diag_combined_albedo[b, g, :, i] = comb_alb

                scale = solar_flux[b, g] * mu0 * w
                for k in range(nlev + 1):
                    up_band[b, k, i] += flux_up_dif[k] * scale
                    down_band[b, k, i] += (flux_dn_dir[k] + flux_dn_dif[k]) * scale

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    diag = {
        "Rdif": diag_Rdif,
        "Tdif": diag_Tdif,
        "Tnoscat": diag_Tnoscat,
        "direct_beam_flux": diag_direct,
    }
    if diagnostics_level >= 2:
        diag["Rdir"] = diag_Rdir
        diag["Tdir"] = diag_Tdir
        diag["combined_albedo"] = diag_combined_albedo
        diag["tau_delta"] = diag_tau_d
        diag["ssa_delta"] = diag_ssa_d
        diag["g_delta"] = diag_g_d

    return up_band, down_band, up_broad, down_broad, diag
```

- [x] **Step 4: Run tests to verify they pass**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_kernels.py -v`
Expected: ALL PASS (including existing tests — the `diagnostics_level=0` default preserves backward compatibility)

- [x] **Step 5: Commit**

```bash
git add climt/_components/picket_fence/sw/kernels.py tests/test_picket_fence_kernels.py
git commit -m "feat(picket-fence): add diagnostics_level to SW two-stream kernel"
```

---

### Task 2: LW kernel diagnostics — `diagnostics_level` parameter

**Files:**
- Modify: `climt/_components/picket_fence/lw/kernels.py`
- Modify: `tests/test_picket_fence_kernels.py`

Same pattern as Task 1 but for the LW kernel. Intermediate diagnostics for LW are: per-layer transmittance, per-g-point upward/downward flux profiles.

- [x] **Step 1: Write failing tests for LW diagnostics**

Add to `tests/test_picket_fence_kernels.py`:

```python
def test_lw_transport_diagnostics_level0_unchanged():
    """diagnostics_level=0 returns the same 4-tuple as before."""
    from climt._components.picket_fence.lw.kernels import lw_transport

    nlev, ncol, nband, ngpt = 5, 1, 2, 1
    tau = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    sigma = 5.670374419e-8
    T = 250.0 * np.ones((nlev, ncol))
    T_surf = np.array([280.0])
    planck_src = sigma * T[np.newaxis, np.newaxis, :, :] ** 4 * np.array([[0.6], [0.4]])[:, :, np.newaxis, np.newaxis]
    surf_src = sigma * T_surf[np.newaxis, np.newaxis, :] ** 4 * np.array([[0.6], [0.4]])[:, :, np.newaxis]
    emissivity = np.ones((nband, ncol))
    weights = np.ones((nband, ngpt))

    result = lw_transport(T, T_surf, tau, planck_src, surf_src, emissivity, weights, sigma,
                          diagnostics_level=0)
    assert len(result) == 4, "diagnostics_level=0 should return 4-tuple"


def test_lw_transport_diagnostics_level1():
    """diagnostics_level=1 returns per-layer transmittance and per-gpoint fluxes."""
    from climt._components.picket_fence.lw.kernels import lw_transport

    nlev, ncol, nband, ngpt = 5, 1, 2, 1
    tau = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    sigma = 5.670374419e-8
    T = 250.0 * np.ones((nlev, ncol))
    T_surf = np.array([280.0])
    planck_src = sigma * T[np.newaxis, np.newaxis, :, :] ** 4 * np.array([[0.6], [0.4]])[:, :, np.newaxis, np.newaxis]
    surf_src = sigma * T_surf[np.newaxis, np.newaxis, :] ** 4 * np.array([[0.6], [0.4]])[:, :, np.newaxis]
    emissivity = np.ones((nband, ncol))
    weights = np.ones((nband, ngpt))

    result = lw_transport(T, T_surf, tau, planck_src, surf_src, emissivity, weights, sigma,
                          diagnostics_level=1)
    assert len(result) == 5
    _, _, _, _, diag = result

    assert "transmittance" in diag
    assert "up_per_gpoint" in diag
    assert "down_per_gpoint" in diag

    assert diag["transmittance"].shape == (nband, ngpt, nlev, ncol)
    assert diag["up_per_gpoint"].shape == (nband, ngpt, nlev + 1, ncol)
    assert diag["down_per_gpoint"].shape == (nband, ngpt, nlev + 1, ncol)

    # Transmittance must be in [0, 1]
    assert np.all(diag["transmittance"] >= 0)
    assert np.all(diag["transmittance"] <= 1)
```

- [x] **Step 2: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_kernels.py::test_lw_transport_diagnostics_level0_unchanged tests/test_picket_fence_kernels.py::test_lw_transport_diagnostics_level1 -v`
Expected: FAIL

- [x] **Step 3: Implement diagnostics_level in lw_transport**

Same pattern as SW: keep `_lw_transport_single_gpt` jitted, make `lw_transport` a plain-Python function that optionally collects intermediate arrays. In `climt/_components/picket_fence/lw/kernels.py`:

```python
def lw_transport(
    T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma,
    diagnostics_level=0,
):
    """Multi-band, multi-g-point LW radiative transfer.

    Args:
        ... (existing args unchanged) ...
        diagnostics_level: 0 (fluxes only), 1 (per-layer transmittance + per-gpoint fluxes)

    Returns:
        If diagnostics_level == 0:
            (up_band, down_band, up_broad, down_broad)
        If diagnostics_level > 0:
            (up_band, down_band, up_broad, down_broad, diagnostics_dict)
    """
    DIFFUSIVITY_FACTOR = 1.66
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    if diagnostics_level >= 1:
        diag_trans = np.zeros((nband, ngpt, nlev, ncol))
        diag_up_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))
        diag_dn_gpt = np.zeros((nband, ngpt, nlev + 1, ncol))

    for b in range(nband):
        for g in range(ngpt):
            up_gpt = np.zeros((nlev + 1, ncol))
            down_gpt = np.zeros((nlev + 1, ncol))
            _lw_transport_single_gpt(
                tau[b, g, :, :],
                planck_source[b, g, :, :],
                surface_source[b, g, :],
                emissivity[b, :],
                up_gpt,
                down_gpt,
                nlev,
                ncol,
            )

            if diagnostics_level >= 1:
                for k in range(nlev):
                    for i in range(ncol):
                        diag_trans[b, g, k, i] = np.exp(
                            -DIFFUSIVITY_FACTOR * tau[b, g, k, i]
                        )
                diag_up_gpt[b, g, :, :] = up_gpt * weights[b, g]
                diag_dn_gpt[b, g, :, :] = down_gpt * weights[b, g]

            w = weights[b, g]
            for k in range(nlev + 1):
                for i in range(ncol):
                    up_band[b, k, i] += w * up_gpt[k, i]
                    down_band[b, k, i] += w * down_gpt[k, i]

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    if diagnostics_level >= 1:
        diag = {
            "transmittance": diag_trans,
            "up_per_gpoint": diag_up_gpt,
            "down_per_gpoint": diag_dn_gpt,
        }
        return up_band, down_band, up_broad, down_broad, diag

    return up_band, down_band, up_broad, down_broad
```

- [x] **Step 4: Run tests to verify they pass**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_kernels.py -v`
Expected: ALL PASS

- [x] **Step 5: Commit**

```bash
git add climt/_components/picket_fence/lw/kernels.py tests/test_picket_fence_kernels.py
git commit -m "feat(picket-fence): add diagnostics_level to LW transport kernel"
```

---

### Task 3: Component-level per-band pedagogical diagnostics

**Files:**
- Modify: `climt/_components/picket_fence/lw/component.py`
- Modify: `climt/_components/picket_fence/sw/component.py`
- Modify: `tests/test_picket_fence_lw.py`
- Modify: `tests/test_picket_fence_sw.py`

The primary pedagogical value of the picket-fence scheme is exposing per-band optical depths, transmittances, and heating rates — letting students see the spectral physics that compiled codes hide. These diagnostics are always returned (no special flag needed) since the tau arrays are already computed in `array_call`.

- [x] **Step 1: Write failing tests for LW per-band diagnostics**

Add to `tests/test_picket_fence_lw.py`:

```python
def test_lw_per_band_optical_depth_diagnostics(get_default_state_lw):
    """LW component returns per-band optical depth and transmittance."""
    from climt._components.picket_fence import PicketFenceLongwave

    lw = PicketFenceLongwave(optics="parmentier")
    state = get_default_state_lw(lw)

    tend, diag = lw(state)
    assert "longwave_optical_depth_per_band" in diag
    assert "longwave_transmittance_per_band" in diag
    assert "longwave_heating_rate_per_band" in diag

    tau = diag["longwave_optical_depth_per_band"].values
    trans = diag["longwave_transmittance_per_band"].values

    # Optical depth must be non-negative
    assert np.all(tau >= 0)
    # Transmittance must be in [0, 1]
    assert np.all(trans >= 0)
    assert np.all(trans <= 1)
    # Transmittance should be consistent with optical depth (exp(-D*tau))
    D = 1.66
    expected_trans = np.exp(-D * tau)
    np.testing.assert_allclose(trans, expected_trans, rtol=1e-10)


def test_lw_per_band_heating_rate_sums_to_broadband(get_default_state_lw):
    """Per-band heating rates should sum to the broadband heating rate."""
    from climt._components.picket_fence import PicketFenceLongwave

    lw = PicketFenceLongwave(optics="parmentier")
    state = get_default_state_lw(lw)

    tend, diag = lw(state)
    hr_broad = diag["longwave_heating_rate"].values
    hr_band = diag["longwave_heating_rate_per_band"].values
    # Sum over band dimension (last axis)
    hr_band_sum = hr_band.sum(axis=-1)
    np.testing.assert_allclose(hr_band_sum, hr_broad, rtol=1e-10)
```

- [x] **Step 2: Write failing tests for SW per-band diagnostics**

Add to `tests/test_picket_fence_sw.py`:

```python
def test_sw_per_band_optical_depth_diagnostics(get_default_state_sw):
    """SW component returns per-band optical depth."""
    from climt._components.picket_fence import PicketFenceShortwave

    sw = PicketFenceShortwave(optics="parmentier")
    state = get_default_state_sw(sw)
    state["zenith_angle"].values[:] = np.pi / 4

    tend, diag = sw(state)
    assert "shortwave_optical_depth_per_band" in diag
    assert "shortwave_heating_rate_per_band" in diag

    tau = diag["shortwave_optical_depth_per_band"].values
    assert np.all(tau >= 0)


def test_sw_per_band_heating_rate_sums_to_broadband(get_default_state_sw):
    """Per-band SW heating rates should sum to the broadband heating rate."""
    from climt._components.picket_fence import PicketFenceShortwave

    sw = PicketFenceShortwave(optics="parmentier")
    state = get_default_state_sw(sw)
    state["zenith_angle"].values[:] = np.pi / 4

    tend, diag = sw(state)
    hr_broad = diag["shortwave_heating_rate"].values
    hr_band = diag["shortwave_heating_rate_per_band"].values
    hr_band_sum = hr_band.sum(axis=-1)
    np.testing.assert_allclose(hr_band_sum, hr_broad, rtol=1e-10)
```

- [x] **Step 3: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_lw.py::test_lw_per_band_optical_depth_diagnostics tests/test_picket_fence_sw.py::test_sw_per_band_optical_depth_diagnostics -v`
Expected: FAIL (diagnostic keys not present)

- [x] **Step 4: Implement LW per-band diagnostics**

In `climt/_components/picket_fence/lw/component.py`:

Add to `diagnostic_properties`:
```python
            "longwave_optical_depth_per_band": {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "dimensionless",
            },
            "longwave_transmittance_per_band": {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "dimensionless",
            },
            "longwave_heating_rate_per_band": {
                "dims": ["mid_levels", "*", "num_longwave_bands"],
                "units": "degK day^-1",
            },
```

In `array_call`, after the kernel call, compute and add:
```python
        # Per-band optical depth: sum over g-points (weighted)
        # tau shape: (nband, ngpt, nlev, ncol)
        D = 1.66  # diffusivity factor
        tau_band = np.zeros((nband, nlev, ncol))
        for b in range(nband):
            for g in range(ngpt):
                tau_band[b] += weights[b, g] * tau[b, g]
        tau_band_out = np.moveaxis(tau_band, 0, -1).reshape(orig_shape_T + (nband,))

        trans_band = np.exp(-D * tau_band)
        trans_band_out = np.moveaxis(trans_band, 0, -1).reshape(orig_shape_T + (nband,))

        # Per-band heating rate from per-band net flux divergence
        hr_band = np.zeros((nband, nlev, ncol))
        for b in range(nband):
            net_band = up_band[b] - down_band[b]  # (nlev+1, ncol)
            hr_band[b] = compute_heating_rate(net_band, p_int_flat, g, cpd) * 86400.0
        hr_band_out = np.moveaxis(hr_band, 0, -1).reshape(orig_shape_T + (nband,))
```

Then add to the diagnostics dict:
```python
                "longwave_optical_depth_per_band": tau_band_out,
                "longwave_transmittance_per_band": trans_band_out,
                "longwave_heating_rate_per_band": hr_band_out,
```

- [x] **Step 5: Implement SW per-band diagnostics**

Same pattern in `climt/_components/picket_fence/sw/component.py`:

Add to `diagnostic_properties`:
```python
            "shortwave_optical_depth_per_band": {
                "dims": ["mid_levels", "*", "num_shortwave_bands"],
                "units": "dimensionless",
            },
            "shortwave_heating_rate_per_band": {
                "dims": ["mid_levels", "*", "num_shortwave_bands"],
                "units": "degK day^-1",
            },
```

In `array_call`, compute per-band optical depth (weighted sum over g-points) and per-band heating rate from per-band net flux, following the same pattern as LW.

- [x] **Step 6: Run tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_lw.py tests/test_picket_fence_sw.py -v`
Expected: ALL PASS

- [x] **Step 7: Commit**

```bash
git add climt/_components/picket_fence/lw/component.py climt/_components/picket_fence/sw/component.py tests/test_picket_fence_lw.py tests/test_picket_fence_sw.py
git commit -m "feat(picket-fence): add per-band optical depth, transmittance, and heating rate diagnostics"
```

---

### Task 4: SW cloud scattering inputs

**Files:**
- Modify: `climt/_components/picket_fence/sw/component.py`
- Modify: `tests/test_picket_fence_sw.py`

The LW component already accepts `longwave_optical_thickness_due_to_cloud`. The SW component needs three cloud inputs: optical depth, single-scattering albedo, and asymmetry parameter. These default to zero/neutral values and are combined with the gas optical properties before the solver.

- [x] **Step 1: Write failing tests for SW cloud inputs**

Add to `tests/test_picket_fence_sw.py`:

```python
def test_sw_parmentier_cloud_increases_reflection(get_default_state_sw):
    """Adding cloud optical depth in SW increases upwelling flux."""
    from climt._components.picket_fence import PicketFenceShortwave

    sw = PicketFenceShortwave(optics="parmentier")
    state = get_default_state_sw(sw)
    state["zenith_angle"].values[:] = np.pi / 4

    # Clear sky
    tend_clear, diag_clear = sw(state)

    # Cloudy: add a scattering cloud in the middle of the atmosphere
    nlev = state["air_temperature"].shape[0]
    cloud_tau = state["shortwave_optical_thickness_due_to_cloud"].values.copy()
    cloud_tau[nlev // 2, :, :] = 5.0  # thick cloud in one layer
    state["shortwave_optical_thickness_due_to_cloud"].values[:] = cloud_tau

    cloud_ssa = state["single_scattering_albedo_due_to_cloud"].values.copy()
    cloud_ssa[nlev // 2, :, :] = 0.99  # highly scattering
    state["single_scattering_albedo_due_to_cloud"].values[:] = cloud_ssa

    tend_cloud, diag_cloud = sw(state)

    # Cloud should increase TOA upwelling (reflection)
    up_clear = diag_clear["upwelling_shortwave_flux_in_air"].values[-1, :]
    up_cloud = diag_cloud["upwelling_shortwave_flux_in_air"].values[-1, :]
    assert np.all(up_cloud > up_clear), (
        f"Cloud should increase TOA upwelling: clear={up_clear}, cloudy={up_cloud}"
    )
```

- [x] **Step 2: Run test to verify it fails**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_sw.py::test_sw_parmentier_cloud_increases_reflection -v`
Expected: FAIL (KeyError: `shortwave_optical_thickness_due_to_cloud` not in input_properties)

- [x] **Step 3: Add SW cloud inputs to component**

In `climt/_components/picket_fence/sw/component.py`, add to the `input_properties` property (after the existing properties):

```python
        # Cloud optical properties for SW (default zero = clear sky)
        props["shortwave_optical_thickness_due_to_cloud"] = {
            "dims": ["mid_levels", "*", "num_shortwave_bands"],
            "units": "dimensionless",
            "alias": "tau_cloud",
        }
        props["single_scattering_albedo_due_to_cloud"] = {
            "dims": ["mid_levels", "*", "num_shortwave_bands"],
            "units": "dimensionless",
            "alias": "ssa_cloud",
        }
        props["cloud_asymmetry_parameter"] = {
            "dims": ["mid_levels", "*", "num_shortwave_bands"],
            "units": "dimensionless",
            "alias": "g_cloud",
        }
```

Then in `array_call`, just before the call to `sw_two_stream`, add cloud combination:

```python
        # Combine gas and cloud optical properties
        nband_cur = tau.shape[0]
        tau_cloud_flat = state["tau_cloud"].reshape(nlev, ncol, nband_cur)
        ssa_cloud_flat = state["ssa_cloud"].reshape(nlev, ncol, nband_cur)
        g_cloud_flat = state["g_cloud"].reshape(nlev, ncol, nband_cur)

        # Rearrange cloud arrays to (nband, nlev, ncol) then broadcast over ngpt
        tau_c = tau_cloud_flat.transpose(2, 0, 1)[:, np.newaxis, :, :]  # (nband, 1, nlev, ncol)
        ssa_c = ssa_cloud_flat.transpose(2, 0, 1)[:, np.newaxis, :, :]
        g_c = g_cloud_flat.transpose(2, 0, 1)[:, np.newaxis, :, :]

        # Combined optical properties (gas + cloud):
        # tau_total = tau_gas + tau_cloud
        # ssa_total = (tau_gas * ssa_gas + tau_cloud * ssa_cloud) / tau_total
        # g_total = (tau_gas * ssa_gas * g_gas + tau_cloud * ssa_cloud * g_cloud) / (tau_total * ssa_total)
        tau_total = tau + tau_c
        scat_gas = tau * ssa
        scat_cloud = tau_c * ssa_c
        scat_total = scat_gas + scat_cloud
        ssa_total = np.where(tau_total > 0, scat_total / tau_total, 0.0)
        g_total = np.where(scat_total > 0,
                           (scat_gas * asym + scat_cloud * g_c) / scat_total,
                           0.0)
        tau = tau_total
        ssa = ssa_total
        asym = g_total
```

Also add defaults for the new cloud properties in `climt/_core/initialization.py` (they should default to zero). Check if the existing `default_values` dict needs entries; if it uses sympl's default of 0 for unspecified fields this may already work.

- [x] **Step 4: Run tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_sw.py -v`
Expected: ALL PASS

- [x] **Step 5: Run full picket-fence test suite for regressions**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_*.py -v`
Expected: ALL PASS

- [x] **Step 6: Commit**

```bash
git add climt/_components/picket_fence/sw/component.py tests/test_picket_fence_sw.py
git commit -m "feat(picket-fence): add SW cloud optical property inputs (tau, ssa, g)"
```

---

### Task 5: Multi-planet atmospheric properties — TOML profiles and API

**Files:**
- Create: `climt/_core/atmospheric_properties.py`
- Create: `climt/_data/atmospheric_properties/earth.toml`
- Create: `climt/_data/atmospheric_properties/mars.toml`
- Create: `climt/_data/atmospheric_properties/hot_jupiter.toml`
- Modify: `climt/__init__.py`
- Create: `tests/test_atmospheric_properties.py`

- [x] **Step 1: Write failing tests for atmospheric properties**

Create `tests/test_atmospheric_properties.py`:

```python
import pytest
import numpy as np


def test_load_earth_profile():
    """load_atmospheric_properties('earth') sets Earth constants."""
    from climt import load_atmospheric_properties
    from sympl import get_constant

    load_atmospheric_properties("earth")
    g = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g - 9.80665) < 1e-3


def test_load_mars_profile():
    """load_atmospheric_properties('mars') sets Mars constants."""
    from climt import load_atmospheric_properties
    from sympl import get_constant

    load_atmospheric_properties("mars")
    g = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g - 3.721) < 0.01
    cp = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")
    assert abs(cp - 735.0) < 10.0


def test_reset_restores_previous():
    """reset_atmospheric_properties restores prior constants."""
    from climt import load_atmospheric_properties, reset_atmospheric_properties
    from sympl import get_constant

    load_atmospheric_properties("earth")
    g_earth = get_constant("gravitational_acceleration", "m/s^2")

    load_atmospheric_properties("mars")
    g_mars = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g_mars - 3.721) < 0.01

    reset_atmospheric_properties()
    g_restored = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g_restored - g_earth) < 1e-6


def test_load_custom_toml(tmp_path):
    """load_atmospheric_properties with a file path loads custom TOML."""
    from climt import load_atmospheric_properties
    from sympl import get_constant

    toml_content = """
[planetary]
gravitational_acceleration = { value = 24.79, units = "m/s^2" }

[bulk_atmosphere]
molar_mass_of_dry_air = { value = 2.2, units = "g/mol" }
"""
    p = tmp_path / "jupiter.toml"
    p.write_text(toml_content)

    load_atmospheric_properties(str(p))
    g = get_constant("gravitational_acceleration", "m/s^2")
    assert abs(g - 24.79) < 0.01


def test_load_nonexistent_profile_raises():
    """Unknown profile name raises FileNotFoundError."""
    from climt import load_atmospheric_properties

    with pytest.raises(FileNotFoundError):
        load_atmospheric_properties("nonexistent_planet_xyz")
```

- [x] **Step 2: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_atmospheric_properties.py -v`
Expected: FAIL (ImportError — load_atmospheric_properties not defined)

- [x] **Step 3: Create Earth TOML profile**

Create `climt/_data/atmospheric_properties/earth.toml`:

```toml
# Earth standard atmosphere constants

[planetary]
gravitational_acceleration = { value = 9.80665, units = "m/s^2" }
planetary_radius = { value = 6371000.0, units = "m" }
planetary_rotation_rate = { value = 7.292e-05, units = "s^-1" }

[bulk_atmosphere]
molar_mass_of_dry_air = { value = 28.970, units = "g/mol" }
gas_constant_of_dry_air = { value = 287.0, units = "J/kg/K" }
heat_capacity_of_dry_air_at_constant_pressure = { value = 1004.64, units = "J/kg/K" }

[gas_species]
molar_mass_of_water_vapor = { value = 18.015, units = "g/mol" }
molar_mass_of_carbon_dioxide = { value = 44.010, units = "g/mol" }
molar_mass_of_ozone = { value = 47.998, units = "g/mol" }
molar_mass_of_methane = { value = 16.043, units = "g/mol" }
molar_mass_of_nitrous_oxide = { value = 44.013, units = "g/mol" }
```

- [x] **Step 4: Create Mars TOML profile**

Create `climt/_data/atmospheric_properties/mars.toml`:

```toml
# Mars CO2-dominated atmosphere constants

[planetary]
gravitational_acceleration = { value = 3.721, units = "m/s^2" }
planetary_radius = { value = 3389500.0, units = "m" }
planetary_rotation_rate = { value = 7.088e-05, units = "s^-1" }

[bulk_atmosphere]
molar_mass_of_dry_air = { value = 43.34, units = "g/mol" }
gas_constant_of_dry_air = { value = 191.8, units = "J/kg/K" }
heat_capacity_of_dry_air_at_constant_pressure = { value = 735.0, units = "J/kg/K" }

[gas_species]
molar_mass_of_carbon_dioxide = { value = 44.010, units = "g/mol" }
molar_mass_of_water_vapor = { value = 18.015, units = "g/mol" }
```

- [x] **Step 5: Create hot Jupiter TOML profile**

Create `climt/_data/atmospheric_properties/hot_jupiter.toml`:

```toml
# Generic H2/He-dominated hot Jupiter atmosphere

[planetary]
gravitational_acceleration = { value = 9.42, units = "m/s^2" }
planetary_radius = { value = 9.44e7, units = "m" }
planetary_rotation_rate = { value = 2.06e-05, units = "s^-1" }

[bulk_atmosphere]
molar_mass_of_dry_air = { value = 2.3, units = "g/mol" }
gas_constant_of_dry_air = { value = 3614.0, units = "J/kg/K" }
heat_capacity_of_dry_air_at_constant_pressure = { value = 13000.0, units = "J/kg/K" }

[gas_species]
molar_mass_of_water_vapor = { value = 18.015, units = "g/mol" }
molar_mass_of_carbon_dioxide = { value = 44.010, units = "g/mol" }
```

- [x] **Step 6: Implement atmospheric_properties.py**

Create `climt/_core/atmospheric_properties.py`:

```python
"""Multi-planet atmospheric properties via TOML profiles.

Loads named or custom atmospheric profiles that set sympl constants
for gravitational acceleration, molar masses, heat capacities, etc.
Supports snapshot/restore so users can switch between planets in a
single session.
"""
import copy
import os

from sympl import get_constant, set_constant

try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:
    import tomli as tomllib  # fallback for 3.10

# Path to built-in profile directory
_BUILTIN_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    "_data",
    "atmospheric_properties",
)

# Snapshot stack (list of dicts). Each entry is a dict of
# {constant_name: (value, units)} captured before a profile was loaded.
_snapshot_stack = []


def _resolve_profile_path(name_or_path):
    """Return the TOML file path for a profile name or file path."""
    if os.path.isfile(name_or_path):
        return name_or_path
    builtin = os.path.join(_BUILTIN_DIR, f"{name_or_path}.toml")
    if os.path.isfile(builtin):
        return builtin
    raise FileNotFoundError(
        f"Atmospheric profile '{name_or_path}' not found. "
        f"Looked for file '{name_or_path}' and built-in '{builtin}'."
    )


def _parse_toml(path):
    """Parse a TOML profile and return a flat dict of {name: (value, units)}."""
    with open(path, "rb") as f:
        data = tomllib.load(f)
    constants = {}
    for section in data.values():
        if isinstance(section, dict):
            for key, val in section.items():
                if isinstance(val, dict) and "value" in val and "units" in val:
                    constants[key] = (val["value"], val["units"])
    return constants


def _snapshot_constants(keys):
    """Capture current values of the given constants from sympl."""
    snap = {}
    for key, (_, units) in keys.items():
        try:
            snap[key] = (get_constant(key, units), units)
        except (KeyError, ValueError):
            # Constant didn't exist before — record as absent
            snap[key] = None
    return snap


def load_atmospheric_properties(name_or_path):
    """Load an atmospheric profile and set sympl constants accordingly.

    Takes a snapshot of the current constants before applying changes,
    so that reset_atmospheric_properties() can restore them.

    Args:
        name_or_path: Built-in profile name (e.g., "earth", "mars",
            "hot_jupiter") or path to a custom .toml file.
    """
    path = _resolve_profile_path(name_or_path)
    constants = _parse_toml(path)

    # Snapshot current state of the constants we're about to overwrite
    _snapshot_stack.append(_snapshot_constants(constants))

    # Apply new constants
    for key, (value, units) in constants.items():
        set_constant(key, value, units)


def reset_atmospheric_properties():
    """Restore constants to the state before the last load_atmospheric_properties call.

    Raises:
        RuntimeError: If no profile has been loaded (nothing to reset).
    """
    if not _snapshot_stack:
        raise RuntimeError(
            "No atmospheric profile snapshot to restore. "
            "Call load_atmospheric_properties() first."
        )
    snap = _snapshot_stack.pop()
    for key, val in snap.items():
        if val is not None:
            value, units = val
            set_constant(key, value, units)
```

- [x] **Step 7: Export the API from climt/__init__.py**

Add to `climt/__init__.py` imports:

```python
from ._core.atmospheric_properties import (
    load_atmospheric_properties,
    reset_atmospheric_properties,
)
```

- [x] **Step 8: Run tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_atmospheric_properties.py -v`
Expected: ALL PASS

- [x] **Step 9: Commit**

```bash
git add climt/_core/atmospheric_properties.py climt/_data/atmospheric_properties/ climt/__init__.py tests/test_atmospheric_properties.py
git commit -m "feat(climt): add multi-planet atmospheric properties API with TOML profiles"
```

---

### Task 6: Validation — analytical solutions and RRTMG comparison

**Files:**
- Create: `tests/test_picket_fence_validation.py`

These tests verify physics accuracy rather than code correctness. They compare picket-fence results against known analytical solutions and (where RRTMG is available) against the full correlated-k scheme.

- [x] **Step 1: Write validation tests**

Create `tests/test_picket_fence_validation.py`:

```python
"""Validation tests for picket-fence radiation scheme.

These tests verify physics accuracy against analytical solutions.
They are slower than unit tests and focus on scientific correctness.

Tests call the kernel functions directly (lw_transport, sw_two_stream)
with single-column arrays reshaped to the (nband, ngpt, nlev, ncol) format.
"""
import numpy as np
import pytest

from sympl import get_constant


def _run_lw_single_col(tau_1d, T_layer_1d, T_surface, emissivity=1.0):
    """Helper: run LW transport on a single grey column.

    Reshapes 1D arrays to (nband=1, ngpt=1, nlev, ncol=1) and calls
    lw_transport directly.
    """
    from climt._components.picket_fence.lw.kernels import lw_transport

    sigma = 5.670374419e-8
    nlev = len(tau_1d)

    tau = np.asarray(tau_1d).reshape(1, 1, nlev, 1)
    planck_src = sigma * np.asarray(T_layer_1d).reshape(1, 1, nlev, 1) ** 4
    T_s = float(T_surface)
    surf_src = np.array([[[sigma * T_s ** 4]]])
    emis = np.array([[float(emissivity)]])
    weights = np.ones((1, 1))
    T_2d = np.asarray(T_layer_1d).reshape(nlev, 1)
    T_surf_1d = np.array([T_s])

    up_band, down_band, up_broad, down_broad = lw_transport(
        T_2d, T_surf_1d, tau, planck_src, surf_src, emis, weights, sigma
    )
    return {
        "flux_up": up_broad[:, 0],
        "flux_down": down_broad[:, 0],
        "heating_rate": np.diff(up_broad[:, 0] - down_broad[:, 0]),
    }


def _run_sw_single_col(tau_1d, ssa_1d, g_1d, zenith, albedo, solar_flux):
    """Helper: run SW two-stream on a single grey column."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev = len(tau_1d)
    tau = np.asarray(tau_1d).reshape(1, 1, nlev, 1)
    ssa = np.asarray(ssa_1d).reshape(1, 1, nlev, 1)
    g = np.asarray(g_1d).reshape(1, 1, nlev, 1)
    zenith_1d = np.array([float(zenith)])
    albedo_1d = np.array([float(albedo)])
    solar_2d = np.array([[float(solar_flux)]])
    weights = np.ones((1, 1))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, g, zenith_1d, albedo_1d, solar_2d, weights
    )
    return {
        "flux_up": up_broad[:, 0],
        "flux_down": down_broad[:, 0],
    }


class TestLWAnalytical:
    """Compare LW solver against analytical solutions."""

    def test_grey_isothermal_equilibrium(self):
        """In an isothermal grey atmosphere, net flux should be zero at all levels.

        An isothermal atmosphere with emissivity=1 surface at the same temperature
        should have zero net flux everywhere (thermal equilibrium).
        """
        nlev = 30
        T_iso = 250.0
        tau_per_layer = 2.0 / nlev * np.ones(nlev)

        result = _run_lw_single_col(
            tau_per_layer, T_iso * np.ones(nlev), T_iso, emissivity=1.0
        )

        hr = result["heating_rate"]
        assert np.all(np.abs(hr) < 0.05 * np.max(np.abs(result["flux_up"]))), (
            f"Isothermal atmosphere should have near-zero heating rate, "
            f"max |dF_net| = {np.max(np.abs(hr)):.6f}"
        )

    def test_optically_thick_limit_olr(self):
        """In the optically thick limit, OLR approaches sigma * T_skin^4.

        With very large optical depth, the outgoing LW radiation at TOA should
        approach the Planck emission from the topmost layer.
        """
        nlev = 20
        T_layer = np.linspace(300, 200, nlev)
        tau_thick = 50.0 / nlev * np.ones(nlev)

        result = _run_lw_single_col(tau_thick, T_layer, 300.0, emissivity=1.0)

        sigma = 5.670374419e-8
        olr = result["flux_up"][-1]
        T_top = T_layer[-1]
        expected_olr = sigma * T_top ** 4
        assert abs(olr - expected_olr) / expected_olr < 0.20, (
            f"OLR={olr:.2f} W/m^2, expected ~{expected_olr:.2f} W/m^2 "
            f"(sigma * T_top^4 for T_top={T_top:.1f} K)"
        )

    def test_optically_thin_limit_olr(self):
        """In the optically thin limit, OLR approaches sigma * T_surface^4.

        With negligible optical depth, the atmosphere is transparent and OLR
        equals the surface Planck emission.
        """
        nlev = 10
        T_layer = np.linspace(300, 200, nlev)
        tau_thin = 1e-6 / nlev * np.ones(nlev)

        result = _run_lw_single_col(tau_thin, T_layer, 300.0, emissivity=1.0)

        sigma = 5.670374419e-8
        olr = result["flux_up"][-1]
        expected = sigma * 300.0 ** 4
        np.testing.assert_allclose(olr, expected, rtol=1e-3)


class TestSWAnalytical:
    """Compare SW solver against analytical solutions."""

    def test_transparent_atmosphere_flux_conservation(self):
        """With zero optical depth, downward flux equals mu0 * S at all levels."""
        nlev = 10
        S0 = 1361.0
        zenith = np.pi / 4
        mu0 = np.cos(zenith)

        result = _run_sw_single_col(
            np.zeros(nlev), np.zeros(nlev), np.zeros(nlev),
            zenith, 0.0, S0,
        )

        expected_down = S0 * mu0
        np.testing.assert_allclose(result["flux_down"], expected_down, rtol=1e-10)
        np.testing.assert_allclose(result["flux_up"], 0.0, atol=1e-10)

    def test_conservative_scattering_energy_conservation(self):
        """With ssa=1, total absorbed energy is zero, net flux constant."""
        nlev = 15
        result = _run_sw_single_col(
            0.5 * np.ones(nlev), np.ones(nlev), np.zeros(nlev),
            np.pi / 3, 0.3, 1000.0,
        )

        net_flux = result["flux_down"] - result["flux_up"]
        np.testing.assert_allclose(net_flux, net_flux[0], rtol=1e-3)

    def test_beer_law_direct_beam(self):
        """With ssa=0, downward flux follows Beer's law exponential decay."""
        nlev = 10
        tau_per_layer = 0.2
        zenith = np.pi / 6
        mu0 = np.cos(zenith)
        S0 = 1000.0

        result = _run_sw_single_col(
            tau_per_layer * np.ones(nlev), np.zeros(nlev), np.zeros(nlev),
            zenith, 0.0, S0,
        )

        total_tau = nlev * tau_per_layer
        expected_surface = S0 * mu0 * np.exp(-total_tau / mu0)
        actual_surface = result["flux_down"][0]
        # With ssa=0, all flux is direct: flux_down = direct beam
        np.testing.assert_allclose(actual_surface, expected_surface, rtol=1e-3)


@pytest.mark.skipif(
    not _rrtmg_available(),
    reason="RRTMG not available for comparison"
)
class TestRRTMGComparison:
    """Compare picket-fence against RRTMG for Earth standard atmosphere."""

    pass  # Placeholder: populated when RRTMG comparison is implemented


def _rrtmg_available():
    """Check if RRTMGLongwave/Shortwave can be imported."""
    try:
        from climt import RRTMGLongwave, RRTMGShortwave
        return True
    except (ImportError, OSError):
        return False
```

- [x] **Step 2: Run validation tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_validation.py -v`
Expected: ALL PASS (analytical tests should pass; RRTMG tests skipped if not available)

- [x] **Step 3: Commit**

```bash
git add tests/test_picket_fence_validation.py
git commit -m "test(picket-fence): add validation tests against analytical LW/SW solutions"
```

---

### Task 7: Wire diagnostics_level through the components

**Files:**
- Modify: `climt/_components/picket_fence/sw/component.py`
- Modify: `climt/_components/picket_fence/lw/component.py`

Pass `diagnostics_level` from the component constructor through to the kernel calls, and include the diagnostic arrays in the component's output when non-zero.

- [x] **Step 1: Write failing test for component-level diagnostics**

Add to `tests/test_picket_fence_sw.py`:

```python
def test_sw_component_diagnostics_level(get_default_state_sw):
    """SW component with diagnostics_level=1 includes extra diagnostics."""
    from climt._components.picket_fence import PicketFenceShortwave

    sw = PicketFenceShortwave(optics="parmentier", diagnostics_level=1)
    state = get_default_state_sw(sw)
    state["zenith_angle"].values[:] = np.pi / 4

    tend, diag = sw(state)
    assert "sw_layer_diffuse_reflectance" in diag
    assert "sw_direct_beam_profile" in diag
```

- [x] **Step 2: Run test to verify it fails**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_sw.py::test_sw_component_diagnostics_level -v`
Expected: FAIL (PicketFenceShortwave does not accept diagnostics_level)

- [x] **Step 3: Wire diagnostics_level through SW component**

In `PicketFenceShortwave.__init__`, add:

```python
        self._diagnostics_level = kwargs.pop("diagnostics_level", 0)
```

(Add this before the `super().__init__(**kwargs)` call.)

In `array_call`, change the `sw_two_stream` call to pass `diagnostics_level=self._diagnostics_level`, and handle the extra return value:

```python
        result = sw_two_stream(
            tau, ssa, asym, zenith_flat, albedo_flat, solar_flux, weights,
            diagnostics_level=self._diagnostics_level,
        )
        if self._diagnostics_level > 0:
            up_band, down_band, up_broad, down_broad, kernel_diag = result
        else:
            up_band, down_band, up_broad, down_broad = result
```

Then after building the diagnostics dict, add the kernel diagnostics with component-level names:

```python
        diagnostics = {
            "upwelling_shortwave_flux_in_air": up_broad_out,
            "downwelling_shortwave_flux_in_air": down_broad_out,
            "upwelling_shortwave_flux_in_air_per_band": up_band_out,
            "downwelling_shortwave_flux_in_air_per_band": down_band_out,
            "shortwave_heating_rate": heating_kday,
        }
        if self._diagnostics_level > 0:
            _DIAG_NAMES = {
                "Rdif": "sw_layer_diffuse_reflectance",
                "Tdif": "sw_layer_diffuse_transmittance",
                "Tnoscat": "sw_layer_direct_transmittance",
                "direct_beam_flux": "sw_direct_beam_profile",
                "Rdir": "sw_layer_direct_reflectance",
                "Tdir": "sw_layer_direct_source_transmittance",
                "combined_albedo": "sw_combined_albedo",
                "tau_delta": "sw_delta_scaled_optical_depth",
                "ssa_delta": "sw_delta_scaled_ssa",
                "g_delta": "sw_delta_scaled_asymmetry",
            }
            for kkey, cname in _DIAG_NAMES.items():
                if kkey in kernel_diag:
                    diagnostics[cname] = kernel_diag[kkey]
```

Also update `diagnostic_properties` to include these when `diagnostics_level > 0`. Add a dynamic property:

```python
    @property
    def diagnostic_properties(self):
        props = {
            # ... existing props ...
        }
        if self._diagnostics_level >= 1:
            props["sw_layer_diffuse_reflectance"] = {
                "dims": ["num_shortwave_bands", "*", "mid_levels", "*"],
                "units": "dimensionless",
            }
            props["sw_direct_beam_profile"] = {
                "dims": ["num_shortwave_bands", "*", "interface_levels", "*"],
                "units": "W m^-2",
            }
            # Add remaining Level 1 diagnostics similarly
        return props
```

Note: The exact dim specification for diagnostic arrays needs to match the kernel output shape. Since kernel diagnostic arrays have shape `(nband, ngpt, nlev, ncol)` which doesn't fit neatly into sympl's dim conventions, the simplest approach is to return them as raw numpy arrays in the diagnostics dict rather than trying to describe their dims in `diagnostic_properties`. sympl will pass through any extra keys in the diagnostics dict.

- [x] **Step 4: Wire diagnostics_level through LW component**

Apply the same pattern to `PicketFenceLongwave`:
- Add `self._diagnostics_level = kwargs.pop("diagnostics_level", 0)` in `__init__`
- Pass `diagnostics_level=self._diagnostics_level` to `lw_transport`
- Include kernel diagnostics with names like `lw_layer_transmittance`, `lw_up_per_gpoint`, `lw_down_per_gpoint`

- [x] **Step 5: Run tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_sw.py tests/test_picket_fence_lw.py -v`
Expected: ALL PASS

- [x] **Step 6: Commit**

```bash
git add climt/_components/picket_fence/sw/component.py climt/_components/picket_fence/lw/component.py tests/test_picket_fence_sw.py
git commit -m "feat(picket-fence): wire diagnostics_level through LW and SW components"
```

---

### Task 8: Run full test suite and verify no regressions

**Files:** None (test-only task)

- [x] **Step 1: Run all picket-fence tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_picket_fence_*.py tests/test_atmospheric_properties.py -v`
Expected: ALL PASS

- [x] **Step 2: Run the broader climt test suite for regressions**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/ -v --timeout=120 -x`
Expected: ALL PASS (or only pre-existing failures unrelated to picket-fence)
