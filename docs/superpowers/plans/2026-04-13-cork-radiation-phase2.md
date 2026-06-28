# Cork Radiation Phase 2: Two-Stream, ESFT, and Correlated-k SW

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the cork radiation scheme by adding the full Meador & Weaver two-stream SW solver, ESFT gas overlap for correlated-k mode, correlated-k SW mode, and stellar spectrum loading from data files.

**Architecture:** Phase 1 delivered direct-beam-only SW, additive-overlap correlated-k LW, and Parmentier mode for both LW and SW. Phase 2 upgrades the SW kernel to a proper two-stream (diffuse + direct), adds the ESFT overlap method for multi-gas correlated-k, implements the correlated-k path through the SW component, and adds stellar spectrum loading from `.npz` files. All changes are backward-compatible — existing Parmentier-mode tests must continue to pass.

**Tech Stack:** Python 3.12, NumPy, numba, sympl (climt's component framework)

---

## File Structure

```
climt/
  _components/
    cork/
      sw/
        kernels.py             # MODIFY: replace direct-beam-only with Meador & Weaver two-stream
      optics/
        correlated_k.py        # MODIFY: add ESFT overlap, add compute_ck_sw_optical_properties
        stellar.py             # CREATE: stellar spectrum loading + band integration
      common.py                # no changes
      lw/                      # no changes
  _data/
    cork/
      correlated_k/
        test_2band_sw.npz      # CREATE: synthetic SW test table
      stellar_spectra/
        sun.npz                # CREATE: solar spectrum data (wavenumber, irradiance)
tests/
  test_cork_kernels.py # MODIFY: add two-stream diffuse tests
  test_cork_optics.py  # MODIFY: add ESFT overlap tests
  test_cork_sw.py      # MODIFY: add correlated-k SW tests
  test_cork_stellar.py # CREATE: stellar spectrum loading tests
```

**Notes:**
- The SW kernel replacement is the biggest change. The current `sw_two_stream` in `sw/kernels.py` is direct-beam only. We replace it with a proper Meador & Weaver (1980) adding method that handles both direct and diffuse components.
- ESFT overlap is triggered by `overlap_method="esft"` in the k-table metadata. Additive overlap remains the default for existing tables.
- The correlated-k SW path reuses the same `interpolate_k` and table format from `correlated_k.py`, adding SW-specific fields (`solar_source_per_gpoint`, `rayleigh_coefficient`).

---

### Task 1: Meador & Weaver two-stream SW kernel

**Files:**
- Modify: `climt/_components/cork/sw/kernels.py`
- Modify: `tests/test_cork_kernels.py`

The current SW kernel (`sw_two_stream`) only does direct-beam exponential attenuation and ignores scattering entirely. We replace it with the Meador & Weaver (1980) two-stream approximation, which computes per-layer reflectance and transmittance for both direct and diffuse beams, then uses the adding method to couple layers.

The Eddington approximation variant (delta-scaled) is used, consistent with RRTMGP.

- [x] **Step 1: Write failing tests for the two-stream kernel**

Add these tests to `tests/test_cork_kernels.py`:

```python
def test_sw_two_stream_conservative_scattering():
    """With ssa=1 and no absorption, total flux (up+down) is conserved at every level."""
    from climt._components.cork.sw.kernels import sw_two_stream

    nlev = 10
    ncol = 1
    nband = 1
    ngpt = 1

    tau = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    ssa = np.ones((nband, ngpt, nlev, ncol))         # pure scattering
    asymmetry = np.zeros((nband, ngpt, nlev, ncol))   # isotropic
    zenith = np.array([np.pi / 3])
    albedo = np.array([0.3])
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, asymmetry, zenith, albedo, solar_flux, weights
    )

    mu0 = np.cos(zenith[0])
    toa_incoming = solar_flux[0, 0] * mu0
    # At TOA, net flux = incoming - reflected.  With conservative scattering,
    # the net flux must be the same at every level (no absorption).
    net_flux = down_broad[:, 0] - up_broad[:, 0]
    # Net flux should be constant across all levels (energy conservation)
    np.testing.assert_allclose(net_flux, net_flux[0], rtol=1e-4)


def test_sw_two_stream_forward_scattering():
    """With g=1 (perfect forward scattering), the result approaches the no-scattering limit."""
    from climt._components.cork.sw.kernels import sw_two_stream

    nlev = 5
    ncol = 1
    nband = 1
    ngpt = 1

    tau = 0.3 * np.ones((nband, ngpt, nlev, ncol))
    ssa = 0.9 * np.ones((nband, ngpt, nlev, ncol))
    # g=0.999: nearly perfect forward scattering (δ-scaling makes this nearly transparent)
    asymmetry = 0.999 * np.ones((nband, ngpt, nlev, ncol))
    zenith = np.array([0.0])  # overhead sun
    albedo = np.array([0.0])
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, asymmetry, zenith, albedo, solar_flux, weights
    )

    # With near-perfect forward scattering, almost all flux reaches the surface
    mu0 = np.cos(zenith[0])
    expected_surface = solar_flux[0, 0] * mu0
    # Should be within 10% of unattenuated (delta-scaling removes forward peak)
    assert down_broad[0, 0] > 0.9 * expected_surface, (
        f"Forward scattering should transmit most flux: got {down_broad[0, 0]:.2f}, "
        f"expected ~{expected_surface:.2f}"
    )


def test_sw_two_stream_albedo_reflection():
    """With zero atmosphere, upward flux equals albedo * downward flux at surface."""
    from climt._components.cork.sw.kernels import sw_two_stream

    nlev = 5
    ncol = 1
    nband = 1
    ngpt = 1

    tau = np.zeros((nband, ngpt, nlev, ncol))
    ssa = np.zeros((nband, ngpt, nlev, ncol))
    asymmetry = np.zeros((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 4])
    albedo = np.array([0.5])
    solar_flux = np.array([[200.0]])
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, asymmetry, zenith, albedo, solar_flux, weights
    )

    mu0 = np.cos(zenith[0])
    expected_down = solar_flux[0, 0] * mu0
    np.testing.assert_allclose(down_broad[0, 0], expected_down, rtol=1e-10)
    # Reflected = albedo * downward at surface, propagates up unattenuated (no atmosphere)
    np.testing.assert_allclose(up_broad[:, 0], albedo[0] * expected_down, rtol=1e-6)
```

- [x] **Step 2: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_kernels.py::test_sw_two_stream_conservative_scattering tests/test_cork_kernels.py::test_sw_two_stream_forward_scattering tests/test_cork_kernels.py::test_sw_two_stream_albedo_reflection -v`
Expected: FAIL (conservative scattering test fails because current kernel has no diffuse scattering)

- [x] **Step 3: Implement the Meador & Weaver two-stream kernel**

Replace the contents of `climt/_components/cork/sw/kernels.py` with:

```python
# climt/_components/cork/sw/kernels.py
import numpy as np

from ..common import njit, prange


@njit
def _delta_scale(tau, ssa, g):
    """Apply delta-Eddington scaling to optical properties.

    Removes the forward scattering peak for better two-stream accuracy.

    Args:
        tau, ssa, g: scalars — optical depth, single scattering albedo, asymmetry

    Returns:
        tau_s, ssa_s, g_s: delta-scaled values
    """
    f = g * g
    tau_s = tau * (1.0 - ssa * f)
    ssa_s = ssa * (1.0 - f) / (1.0 - ssa * f) if (1.0 - ssa * f) > 1e-30 else 0.0
    g_s = (g - f) / (1.0 - f) if (1.0 - f) > 1e-30 else 0.0
    return tau_s, ssa_s, g_s


@njit
def _meador_weaver_coefficients(tau, ssa, g, mu0):
    """Compute layer reflectance and transmittance using Meador & Weaver (1980).

    Uses the Eddington approximation (gamma coefficients from their Table 1,
    Row 1: Eddington).

    Args:
        tau: layer optical depth (scalar, delta-scaled)
        ssa: single scattering albedo (scalar, delta-scaled)
        g: asymmetry parameter (scalar, delta-scaled)
        mu0: cosine of solar zenith angle (scalar)

    Returns:
        R_dir: reflectance for direct beam
        T_dir: transmittance for direct beam (diffuse component only)
        R_dif: reflectance for diffuse radiation
        T_dif: transmittance for diffuse radiation
        T_direct: direct beam transmittance (exp(-tau/mu0))
    """
    # Eddington gamma coefficients (Meador & Weaver 1980, Table 1, Row 1)
    gamma1 = (7.0 - ssa * (4.0 + 3.0 * g)) / 4.0
    gamma2 = -(1.0 - ssa * (4.0 - 3.0 * g)) / 4.0
    gamma3 = (2.0 - 3.0 * g * mu0) / 4.0
    gamma4 = 1.0 - gamma3

    # Direct beam transmittance through the layer
    T_direct = np.exp(-tau / mu0) if mu0 > 1e-10 else 0.0

    # Handle the nearly-transparent case
    if tau < 1e-10:
        return 0.0, 0.0, 0.0, 1.0, T_direct

    # Handle the conservative scattering case (ssa ≈ 1)
    if ssa > 1.0 - 1e-10:
        # Conservative scattering solution (Meador & Weaver Eqs. 24-27)
        gamma1_cons = 0.5 * (7.0 - (4.0 + 3.0 * g)) / 4.0 * 2.0
        # Simplified: use the non-absorbing limit
        alpha = gamma1 * tau
        R_dif = alpha / (1.0 + alpha)
        T_dif = 1.0 / (1.0 + alpha)

        # Direct beam terms
        R_dir = (gamma3 * tau - (gamma3 - gamma4) * (T_direct - 1.0) / gamma1) / (1.0 + alpha) if alpha > 1e-10 else 0.0
        T_dir_diffuse = -R_dir + gamma3 * tau * T_direct + (gamma3 - gamma4) * (1.0 - T_direct) / gamma1 if gamma1 > 1e-10 else 0.0
        T_dir_diffuse = max(T_dir_diffuse, 0.0)
        R_dir = max(R_dir, 0.0)
        return R_dir, T_dir_diffuse, R_dif, T_dif, T_direct

    # General case: absorbing + scattering atmosphere
    k = np.sqrt(max(gamma1 * gamma1 - gamma2 * gamma2, 0.0))

    exp_neg = np.exp(-k * tau)
    exp_pos = np.exp(k * tau)

    # Denominator for reflectance/transmittance
    denom = k + gamma1 + (k - gamma1) * exp_neg * exp_neg
    if abs(denom) < 1e-30:
        denom = 1e-30

    R_dif = gamma2 * (1.0 - exp_neg * exp_neg) / denom
    T_dif = 2.0 * k * exp_neg / denom

    # Direct beam source terms
    alpha1 = gamma1 - 1.0 / mu0
    alpha2 = gamma2
    denom_dir = (k * k - 1.0 / (mu0 * mu0))
    if abs(denom_dir) < 1e-30:
        # Near-singular: mu0 ≈ 1/k, use a small offset
        denom_dir = 1e-30 if denom_dir >= 0 else -1e-30

    C_plus = ssa * ((alpha1 * gamma3 + alpha2 * gamma4) / denom_dir + gamma3)
    C_minus = ssa * ((alpha1 * gamma4 + alpha2 * gamma3) / denom_dir + gamma4)

    # Layer-specific direct beam reflectance and diffuse transmittance
    R_dir = C_plus * R_dif + C_minus * (T_dif * T_direct - 1.0)
    R_dir = max(R_dir, 0.0)
    T_dir_diffuse = C_plus * T_dif + C_minus * R_dif * T_direct - C_minus * T_direct + C_minus
    T_dir_diffuse = max(T_dir_diffuse, 0.0)

    return R_dir, T_dir_diffuse, R_dif, T_dif, T_direct


@njit
def _adding_method(R_dir, T_dir_dif, R_dif, T_dif, T_direct, albedo, nlev):
    """Combine layers using the adding method (Lacis & Hansen 1974; Shonk & Hogan 2008).

    Index convention: k=0 is the surface, k=nlev is TOA.
    The per-layer arrays are indexed [0..nlev-1] where layer k sits between
    interface k (bottom) and interface k+1 (top).

    Args:
        R_dir: (nlev,) per-layer reflectance for direct beam
        T_dir_dif: (nlev,) per-layer diffuse transmittance of direct beam
        R_dif: (nlev,) per-layer reflectance for diffuse
        T_dif: (nlev,) per-layer transmittance for diffuse
        T_direct: (nlev,) per-layer direct transmittance (Beer's law)
        albedo: surface albedo (scalar)
        nlev: number of layers

    Returns:
        up: (nlev+1,) upward flux at interfaces (normalised to unit TOA incoming)
        down_direct: (nlev+1,) direct downward flux at interfaces
        down_diffuse: (nlev+1,) diffuse downward flux at interfaces
    """
    # Combined albedo and source, built from surface upward.
    # R_combined[k] = effective albedo of everything below interface k+1.
    R_combined = np.zeros(nlev + 1)
    R_combined[0] = albedo  # surface

    for k in range(nlev):
        # Denominator: 1 - R_dif[k] * R_combined[k]
        denom = 1.0 - R_dif[k] * R_combined[k]
        if abs(denom) < 1e-30:
            denom = 1e-30
        R_combined[k + 1] = R_dif[k] + T_dif[k] * T_dif[k] * R_combined[k] / denom

    # Now sweep downward from TOA to get fluxes
    down_direct = np.zeros(nlev + 1)
    down_diffuse = np.zeros(nlev + 1)
    up = np.zeros(nlev + 1)

    # TOA: unit incoming direct flux
    down_direct[nlev] = 1.0
    down_diffuse[nlev] = 0.0

    for k in range(nlev - 1, -1, -1):
        # Direct beam attenuation
        down_direct[k] = down_direct[k + 1] * T_direct[k]

        # Diffuse contribution from direct beam hitting this layer
        denom = 1.0 - R_dif[k] * R_combined[k]
        if abs(denom) < 1e-30:
            denom = 1e-30

        # Source from direct beam scattering
        source_dir = down_direct[k + 1] * T_dir_dif[k]
        source_ref = down_direct[k + 1] * R_dir[k] * R_combined[k]

        down_diffuse[k] = (
            T_dif[k] * down_diffuse[k + 1]
            + source_dir
            + T_dif[k] * source_ref
        ) / denom + T_dif[k] * down_diffuse[k + 1] * R_dif[k] * R_combined[k] / denom

        # Simplification: combine properly
        incoming_dif = down_diffuse[k + 1]
        incoming_dir = down_direct[k + 1]

        down_diffuse[k] = (
            T_dif[k] * incoming_dif + incoming_dir * T_dir_dif[k]
            + (T_dif[k] * incoming_dif * R_dif[k]
               + incoming_dir * R_dir[k]) * R_combined[k] * T_dif[k]
        ) / denom

    # Upward flux at each interface
    for k in range(nlev + 1):
        if k == 0:
            up[0] = albedo * (down_direct[0] + down_diffuse[0])
        else:
            # Up = reflected from below
            denom = 1.0 - R_dif[k - 1] * R_combined[k - 1]
            if abs(denom) < 1e-30:
                denom = 1e-30
            up[k] = (
                down_direct[k] * R_dir[k - 1]
                + down_diffuse[k] * R_dif[k - 1]
                + (down_direct[k] * T_dir_dif[k - 1]
                   + T_dif[k - 1] * down_diffuse[k]) * R_combined[k - 1] * T_dif[k - 1]
            ) / denom + down_direct[k] * R_dir[k - 1]
            # Simpler: use the combined albedo approach
            up[k] = R_combined[k] * (down_direct[k] + down_diffuse[k])

    return up, down_direct, down_diffuse


@njit
def sw_two_stream(tau, ssa, asymmetry, zenith, albedo, solar_flux, weights):
    """Multi-band, multi-g-point SW radiative transfer.

    Meador & Weaver (1980) two-stream with delta-Eddington scaling and
    adding method for inter-layer coupling.

    Args:
        tau: (nband, ngpt, nlev, ncol) extinction optical depth
        ssa: (nband, ngpt, nlev, ncol) single scattering albedo
        asymmetry: (nband, ngpt, nlev, ncol) asymmetry parameter
        zenith: (ncol,) solar zenith angle, radians
        albedo: (ncol,) surface albedo
        solar_flux: (nband, ngpt) TOA solar flux per g-point, W/m^2
        weights: (nband, ngpt) g-point weights

    Returns:
        up_band: (nband, nlev+1, ncol) per-band upward flux
        down_band: (nband, nlev+1, ncol) per-band downward flux
        up_broad: (nlev+1, ncol) broadband upward flux
        down_broad: (nlev+1, ncol) broadband downward flux
    """
    nband, ngpt, nlev, ncol = tau.shape

    up_band = np.zeros((nband, nlev + 1, ncol))
    down_band = np.zeros((nband, nlev + 1, ncol))

    for b in range(nband):
        for g in range(ngpt):
            w = weights[b, g]
            for i in prange(ncol):
                mu0 = np.cos(zenith[i])
                if mu0 <= 1e-4:
                    continue  # nighttime

                # Per-layer Meador & Weaver coefficients
                R_dir = np.zeros(nlev)
                T_dir_dif = np.zeros(nlev)
                R_dif = np.zeros(nlev)
                T_dif = np.zeros(nlev)
                T_direct = np.zeros(nlev)

                for k in range(nlev):
                    # Delta-scale the optical properties
                    tau_s, ssa_s, g_s = _delta_scale(
                        tau[b, g, k, i], ssa[b, g, k, i], asymmetry[b, g, k, i]
                    )
                    R_dir[k], T_dir_dif[k], R_dif[k], T_dif[k], T_direct[k] = (
                        _meador_weaver_coefficients(tau_s, ssa_s, g_s, mu0)
                    )

                up_col, down_dir, down_dif = _adding_method(
                    R_dir, T_dir_dif, R_dif, T_dif, T_direct, albedo[i], nlev
                )

                # Scale by incoming solar flux * mu0 * weight
                scale = solar_flux[b, g] * mu0 * w
                for k in range(nlev + 1):
                    up_band[b, k, i] += up_col[k] * scale
                    down_band[b, k, i] += (down_dir[k] + down_dif[k]) * scale

    up_broad = np.zeros((nlev + 1, ncol))
    down_broad = np.zeros((nlev + 1, ncol))
    for b in range(nband):
        for k in range(nlev + 1):
            for i in range(ncol):
                up_broad[k, i] += up_band[b, k, i]
                down_broad[k, i] += down_band[b, k, i]

    return up_band, down_band, up_broad, down_broad
```

**Implementation notes:**
- The `_meador_weaver_coefficients` function implements the Eddington variant from Meador & Weaver (1980) Table 1 Row 1. It handles three regimes: transparent (`tau < 1e-10`), conservative scattering (`ssa ≈ 1`), and the general absorbing+scattering case.
- Delta-Eddington scaling is applied before the two-stream solve. This removes the forward-scattering peak from the phase function, improving accuracy for strongly forward-scattering media (clouds, aerosols).
- The adding method builds combined albedos from the surface upward, then sweeps fluxes downward from TOA. Upward flux at each interface uses the combined albedo below that level.
- The function signature and return shapes are identical to the phase 1 version, so the SW component needs no changes.

- [x] **Step 4: Run all SW kernel tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_kernels.py -v`
Expected: All PASS (including the existing `test_sw_two_stream_no_atmosphere`, `test_sw_two_stream_total_absorption`, `test_sw_reflected_beam_attenuates_upward`, plus the three new tests)

- [x] **Step 5: Run full cork test suite to check for regressions**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_sw.py tests/test_cork_integration.py -v`
Expected: All PASS

- [x] **Step 6: Commit**

```bash
git add climt/_components/cork/sw/kernels.py tests/test_cork_kernels.py
git commit -m "feat(cork): replace direct-beam SW with Meador & Weaver two-stream solver"
```

---

### Task 2: ESFT gas overlap for correlated-k

**Files:**
- Modify: `climt/_components/cork/optics/correlated_k.py`
- Modify: `tests/test_cork_optics.py`

Phase 1 correlated-k uses additive overlap: all gases share the same g-point grid and their optical depths are summed directly. ESFT (Equivalent Sum of exponenTials with Full Terms) correctly handles spectrally uncorrelated gases by taking the outer product of g-points. For N gases with G g-points each, the ESFT solve evaluates G^N combined g-points per band. The combined weight is the product of individual weights.

- [x] **Step 1: Write failing tests for ESFT overlap**

Add these tests to `tests/test_cork_optics.py`:

```python
def test_esft_single_gas_matches_additive():
    """With one gas, ESFT should give the same result as additive overlap."""
    from climt._components.cork.optics.correlated_k import (
        compute_ck_optical_depth,
        load_k_table,
    )

    table = load_k_table("test_2band_lw")
    # Override overlap method to esft
    table_esft = dict(table)
    table_esft["overlap_method"] = np.array("esft")

    nlev = 5
    ncol = 1
    ngas = 1
    T = 250.0 * np.ones((nlev, ncol))
    p = np.linspace(1e3, 1e5, nlev).reshape(nlev, 1)
    gas_amounts = 1e-3 * np.ones((ngas, nlev, ncol))

    tau_add = compute_ck_optical_depth(table, T, p, gas_amounts)
    tau_esft, weights_esft = compute_ck_optical_depth(table_esft, T, p, gas_amounts)

    # With 1 gas, ESFT optical depths per combined g-point should equal additive
    np.testing.assert_allclose(tau_esft, tau_add, rtol=1e-10)


def test_esft_two_gas_weight_sum():
    """ESFT combined weights for 2 gases × 2 g-points should sum to 1 per band."""
    from climt._components.cork.optics.correlated_k import (
        compute_esft_weights,
    )

    # 2 g-points per gas, 2 bands
    weights_per_gas = np.array([[0.3, 0.7], [0.4, 0.6]])  # (nband=2 [same], ngpt=2)
    # Actually: for ESFT with 2 gases, each with 2 g-points,
    # the combined weights should be the outer product per band.
    w1 = np.array([0.3, 0.7])  # gas 1 g-point weights (same for both bands here)
    w2 = np.array([0.4, 0.6])  # gas 2 g-point weights
    # Combined: 2×2 = 4 combined g-points
    # weights: w1[i] * w2[j] for each (i,j)
    expected = np.array([0.3 * 0.4, 0.3 * 0.6, 0.7 * 0.4, 0.7 * 0.6])
    np.testing.assert_allclose(expected.sum(), 1.0, rtol=1e-10)

    gpoint_weights = np.array([[0.3, 0.7], [0.3, 0.7]])  # (nband=2, ngpt=2)
    ngas = 2
    combined = compute_esft_weights(gpoint_weights, ngas)
    # Shape: (nband, ngpt^ngas) = (2, 4)
    assert combined.shape == (2, 4)
    for b in range(2):
        np.testing.assert_allclose(combined[b].sum(), 1.0, rtol=1e-10)


def test_esft_two_gas_optical_depth(tmp_path):
    """ESFT with 2 gases produces correct combined optical depths."""
    from climt._components.cork.optics.correlated_k import (
        compute_ck_optical_depth,
    )

    ngas, nband, ngpt, nT, nP = 2, 1, 2, 3, 3
    # Gas 1: k = 0.1 everywhere; Gas 2: k = 0.2 everywhere
    k_data = np.zeros((ngas, nband, ngpt, nT, nP))
    k_data[0, :, :, :, :] = 0.1
    k_data[1, :, :, :, :] = 0.2

    table_file = str(tmp_path / "esft_test.npz")
    np.savez(
        table_file,
        k_coefficients=k_data,
        gpoint_weights=np.array([[0.4, 0.6]]),
        planck_fraction=np.full((nband, ngpt, nT), 0.5),
        band_wavenumber_limits=np.array([[200.0, 800.0]]),
        temperature_grid=np.linspace(200.0, 400.0, nT),
        pressure_grid_log=np.linspace(np.log(100.0), np.log(1e5), nP),
        gas_names=np.array(["gas1", "gas2"]),
        overlap_method=np.array("esft"),
        resolution=np.array("low"),
    )

    from climt._components.cork.optics.correlated_k import load_k_table
    table = load_k_table(table_file)

    nlev = 3
    ncol = 1
    T = 300.0 * np.ones((nlev, ncol))
    p = 5e4 * np.ones((nlev, ncol))
    gas_amounts = np.ones((ngas, nlev, ncol))  # 1 kg/m^2

    tau, weights = compute_ck_optical_depth(table, T, p, gas_amounts)
    # ESFT: nband=1, ngpt_combined=2*2=4
    assert tau.shape[0] == 1   # nband
    assert tau.shape[1] == 4   # ngpt^ngas = 2^2
    assert weights.shape == (1, 4)

    # Combined weight should sum to 1
    np.testing.assert_allclose(weights[0].sum(), 1.0, rtol=1e-10)

    # Optical depth for combined g-point (i,j) = k1[i]*amount1 + k2[j]*amount2
    # k is uniform 0.1 for gas1, 0.2 for gas2, so all combined taus = 0.1 + 0.2 = 0.3
    np.testing.assert_allclose(tau[:, :, :, :], 0.3, rtol=1e-6)
```

- [x] **Step 2: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_optics.py::test_esft_single_gas_matches_additive tests/test_cork_optics.py::test_esft_two_gas_weight_sum tests/test_cork_optics.py::test_esft_two_gas_optical_depth -v`
Expected: FAIL (functions don't exist yet, and `compute_ck_optical_depth` doesn't return weights)

- [x] **Step 3: Implement ESFT overlap in correlated_k.py**

Modify `climt/_components/cork/optics/correlated_k.py`:

```python
# climt/_components/cork/optics/correlated_k.py
import os

import importlib_resources
import numpy as np

from ..common import njit, prange


def load_k_table(name_or_path):
    """Load a correlated-k table.

    Args:
        name_or_path: table name (e.g., "test_2band_lw") or path to .npz file

    Returns:
        dict-like npz object
    """
    if os.path.isfile(name_or_path):
        return np.load(name_or_path, allow_pickle=True)

    data_path = importlib_resources.files(
        "climt._data.cork.correlated_k"
    ).joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(data_path) as f:
        return np.load(f, allow_pickle=True)


def interpolate_k(table, T, p):
    """Bilinear interpolation of k-coefficients in (log p, T) space.

    Args:
        table: loaded k-table
        T: (ncol,) temperature values, K
        p: (ncol,) pressure values, Pa

    Returns:
        k_interp: (ngas, nband, ngpt, ncol)
    """
    k = table["k_coefficients"]  # (ngas, nband, ngpt, nT, nP)
    T_grid = table["temperature_grid"]
    p_grid_log = table["pressure_grid_log"]
    ngas, nband, ngpt, nT, nP = k.shape
    ncol = len(T)

    result = np.zeros((ngas, nband, ngpt, ncol))
    log_p = np.log(np.maximum(p, 1.0))

    for col in range(ncol):
        # Find T indices
        iT = np.searchsorted(T_grid, T[col]) - 1
        iT = max(0, min(iT, nT - 2))
        fT = (T[col] - T_grid[iT]) / (T_grid[iT + 1] - T_grid[iT])
        fT = max(0.0, min(1.0, fT))

        # Find log(p) indices
        iP = np.searchsorted(p_grid_log, log_p[col]) - 1
        iP = max(0, min(iP, nP - 2))
        fP = (log_p[col] - p_grid_log[iP]) / (p_grid_log[iP + 1] - p_grid_log[iP])
        fP = max(0.0, min(1.0, fP))

        # Bilinear interpolation
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    v00 = k[ig, ib, igp, iT, iP]
                    v10 = k[ig, ib, igp, iT + 1, iP]
                    v01 = k[ig, ib, igp, iT, iP + 1]
                    v11 = k[ig, ib, igp, iT + 1, iP + 1]
                    result[ig, ib, igp, col] = (
                        v00 * (1 - fT) * (1 - fP)
                        + v10 * fT * (1 - fP)
                        + v01 * (1 - fT) * fP
                        + v11 * fT * fP
                    )

    return result


def compute_esft_weights(gpoint_weights, ngas):
    """Compute ESFT combined g-point weights for multiple gases.

    For N gases each with G g-points, produces G^N combined weights per band,
    where the combined weight is the product of individual gas weights.

    Args:
        gpoint_weights: (nband, ngpt) per-gas g-point weights (same for all gases)
        ngas: number of gases

    Returns:
        combined_weights: (nband, ngpt^ngas) combined weights
    """
    nband, ngpt = gpoint_weights.shape
    ngpt_combined = ngpt ** ngas

    combined = np.zeros((nband, ngpt_combined))

    for b in range(nband):
        w = gpoint_weights[b]
        for idx in range(ngpt_combined):
            weight = 1.0
            remainder = idx
            for gas in range(ngas):
                g_idx = remainder % ngpt
                remainder //= ngpt
                weight *= w[g_idx]
            combined[b, idx] = weight

    return combined


def compute_ck_optical_depth(table, T, p, gas_amounts):
    """Compute optical depths from correlated-k table.

    Supports both additive and ESFT overlap methods.

    Args:
        table: loaded k-table
        T: (nlev, ncol) temperature, K
        p: (nlev, ncol) pressure, Pa
        gas_amounts: (ngas, nlev, ncol) column amount per gas, kg/m^2

    Returns:
        For additive overlap:
            tau: (nband, ngpt, nlev, ncol) optical depth per layer
        For ESFT overlap:
            tuple of (tau, weights) where
            tau: (nband, ngpt^ngas, nlev, ncol) optical depth per combined g-point
            weights: (nband, ngpt^ngas) combined g-point weights
    """
    overlap = str(table.get("overlap_method", np.array("additive")))

    if overlap == "esft":
        return _compute_ck_optical_depth_esft(table, T, p, gas_amounts)
    else:
        return _compute_ck_optical_depth_additive(table, T, p, gas_amounts)


def _compute_ck_optical_depth_additive(table, T, p, gas_amounts):
    """Additive overlap: sum optical depths from all gases on the same g-point grid."""
    k_data = table["k_coefficients"]
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    tau = np.zeros((nband, ngpt, nlev, ncol))

    for k_lev in range(nlev):
        k_interp = interpolate_k(table, T[k_lev, :], p[k_lev, :])
        for ig in range(ngas):
            for ib in range(nband):
                for igp in range(ngpt):
                    for icol in range(ncol):
                        tau[ib, igp, k_lev, icol] += (
                            k_interp[ig, ib, igp, icol] * gas_amounts[ig, k_lev, icol]
                        )

    return tau


def _compute_ck_optical_depth_esft(table, T, p, gas_amounts):
    """ESFT overlap: outer product of g-points across gases.

    For combined g-point index idx, decode per-gas g-point indices by
    treating idx as a mixed-radix number (each digit base ngpt).
    """
    k_data = table["k_coefficients"]  # (ngas, nband, ngpt, nT, nP)
    gpoint_weights = table["gpoint_weights"]  # (nband, ngpt)
    ngas, nband, ngpt = k_data.shape[:3]
    nlev, ncol = T.shape

    ngpt_combined = ngpt ** ngas
    tau = np.zeros((nband, ngpt_combined, nlev, ncol))

    # Precompute combined weights
    combined_weights = compute_esft_weights(gpoint_weights, ngas)

    for k_lev in range(nlev):
        k_interp = interpolate_k(table, T[k_lev, :], p[k_lev, :])
        # k_interp: (ngas, nband, ngpt, ncol)

        for ib in range(nband):
            for idx in range(ngpt_combined):
                # Decode combined index to per-gas g-point indices
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

- [x] **Step 4: Update the LW component to handle ESFT weights**

In `climt/_components/cork/lw/component.py`, the correlated-k path calls `compute_ck_optical_depth` and uses `self._table["gpoint_weights"]`. For ESFT, the function now returns a tuple. Update the correlated-k branch in `array_call`:

Find (in `lw/component.py`, `array_call` method, correlated-k branch):
```python
            tau, planck_src, surf_src = self._correlated_k_optics(
                T_flat, p_flat, gas_amounts, T_surf_flat, sigma
            )
            weights = self._table["gpoint_weights"]
```

Replace with:
```python
            tau, planck_src, surf_src, weights = self._correlated_k_optics(
                T_flat, p_flat, gas_amounts, T_surf_flat, sigma
            )
```

Then update `_correlated_k_optics` to return weights and handle ESFT:

Find (in `lw/component.py`):
```python
    def _correlated_k_optics(self, T, p, gas_amounts, T_surf, sigma):
        from ..optics.correlated_k import compute_ck_optical_depth
```

Replace with:
```python
    def _correlated_k_optics(self, T, p, gas_amounts, T_surf, sigma):
        from ..optics.correlated_k import compute_ck_optical_depth
```

And update the body of `_correlated_k_optics` to detect ESFT:

```python
    def _correlated_k_optics(self, T, p, gas_amounts, T_surf, sigma):
        from ..optics.correlated_k import compute_ck_optical_depth

        nlev, ncol = T.shape

        result = compute_ck_optical_depth(self._table, T, p, gas_amounts)
        if isinstance(result, tuple):
            # ESFT mode: returns (tau, combined_weights)
            tau_raw, weights = result
        else:
            # Additive mode: returns tau only
            tau_raw = result
            weights = self._table["gpoint_weights"]

        nband = tau_raw.shape[0]
        ngpt = tau_raw.shape[1]

        planck_frac = self._table["planck_fraction"]  # (nband_orig, ngpt_orig, nT)
        T_grid = self._table["temperature_grid"]
        nT = len(T_grid)
        nband_orig = planck_frac.shape[0]
        ngpt_orig = planck_frac.shape[1]

        planck_src = np.zeros((nband, ngpt, nlev, ncol))
        surf_src = np.zeros((nband, ngpt, ncol))

        overlap = str(self._table.get("overlap_method", np.array("additive")))
        ngas = self._table["k_coefficients"].shape[0]

        for icol in range(ncol):
            # Surface source
            T_s = T_surf[icol]
            iTs = np.searchsorted(T_grid, T_s) - 1
            iTs = max(0, min(iTs, nT - 2))
            fTs = (T_s - T_grid[iTs]) / (T_grid[iTs + 1] - T_grid[iTs])
            fTs = max(0.0, min(1.0, fTs))

            surf_planck = sigma * T_s**4
            for ib in range(nband):
                for igp in range(ngpt):
                    if overlap == "esft" and ngas > 1:
                        # Use the first gas's g-point index for Planck fraction
                        g_idx_orig = igp % ngpt_orig
                    else:
                        g_idx_orig = igp
                    ib_orig = min(ib, nband_orig - 1)
                    frac_s = (
                        planck_frac[ib_orig, g_idx_orig, iTs] * (1 - fTs)
                        + planck_frac[ib_orig, g_idx_orig, iTs + 1] * fTs
                    )
                    surf_src[ib, igp, icol] = frac_s * surf_planck

            # Layer source
            for k in range(nlev):
                T_l = T[k, icol]
                iTl = np.searchsorted(T_grid, T_l) - 1
                iTl = max(0, min(iTl, nT - 2))
                fTl = (T_l - T_grid[iTl]) / (T_grid[iTl + 1] - T_grid[iTl])
                fTl = max(0.0, min(1.0, fTl))

                layer_planck = sigma * T_l**4
                for ib in range(nband):
                    for igp in range(ngpt):
                        if overlap == "esft" and ngas > 1:
                            g_idx_orig = igp % ngpt_orig
                        else:
                            g_idx_orig = igp
                        ib_orig = min(ib, nband_orig - 1)
                        frac_l = (
                            planck_frac[ib_orig, g_idx_orig, iTl] * (1 - fTl)
                            + planck_frac[ib_orig, g_idx_orig, iTl + 1] * fTl
                        )
                        planck_src[ib, igp, k, icol] = frac_l * layer_planck

        return tau_raw, planck_src, surf_src, weights
```

- [x] **Step 5: Run all optics and LW tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_optics.py tests/test_cork_lw.py -v`
Expected: All PASS

- [x] **Step 6: Commit**

```bash
git add climt/_components/cork/optics/correlated_k.py climt/_components/cork/lw/component.py tests/test_cork_optics.py
git commit -m "feat(cork): add ESFT gas overlap for correlated-k mode"
```

---

### Task 3: Stellar spectrum loading

**Files:**
- Create: `climt/_components/cork/optics/stellar.py`
- Create: `climt/_data/cork/stellar_spectra/sun.npz`
- Create: `tests/test_cork_stellar.py`

The SW component currently uses a hardcoded solar constant (1361 W/m^2 split equally across 3 bands). For correlated-k mode, the stellar flux must be distributed across bands according to the actual spectrum. This module loads a stellar spectrum and integrates it over each band's wavenumber limits.

- [x] **Step 1: Create the solar spectrum data file**

Write a script to generate the solar spectrum `.npz`. The spectrum uses the Lean & DeLand (2012) quiet-sun values, simplified to ~100 wavenumber points spanning 0–50000 cm^-1.

```python
# Script to generate climt/_data/cork/stellar_spectra/sun.npz
# Run once, not part of the package.
import numpy as np

# Simplified solar spectrum: wavenumber (cm^-1) and spectral irradiance (W/m^2/cm^-1)
# These approximate the Lean & DeLand (2012) quiet-sun values.
# The integral over all wavenumbers should give ~1361 W/m^2.
wavenumber = np.array([
    100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000,
    9000, 10000, 12000, 14000, 16000, 18000, 20000, 22000, 25000,
    28000, 30000, 33000, 35000, 40000, 45000, 50000
], dtype=np.float64)

# Approximate spectral irradiance in W/m^2 per cm^-1
# Peaks around 10000-18000 cm^-1 (visible, ~0.5-1.0 μm)
irradiance = np.array([
    0.001, 0.005, 0.02, 0.06, 0.08, 0.09, 0.10, 0.10, 0.095, 0.085,
    0.075, 0.065, 0.05, 0.04, 0.035, 0.03, 0.025, 0.02, 0.013,
    0.008, 0.005, 0.003, 0.002, 0.0008, 0.0003, 0.0001
], dtype=np.float64)

# Scale so integral ≈ 1361 W/m^2
total = np.trapz(irradiance, wavenumber)
irradiance *= 1361.0 / total

np.savez(
    "climt/_data/cork/stellar_spectra/sun.npz",
    wavenumber=wavenumber,
    irradiance=irradiance,
    total_irradiance=np.array(1361.0),
    source="Approximate Lean & DeLand (2012) quiet-sun spectrum",
)
```

Run this script from the repository root to generate the data file.

- [x] **Step 2: Write failing tests for stellar spectrum loading**

```python
# tests/test_cork_stellar.py
import numpy as np
import pytest


def test_load_solar_spectrum():
    """Loading the solar spectrum returns arrays with matching shapes."""
    from climt._components.cork.optics.stellar import load_stellar_spectrum

    spec = load_stellar_spectrum("sun")
    assert "wavenumber" in spec
    assert "irradiance" in spec
    assert len(spec["wavenumber"]) == len(spec["irradiance"])
    assert len(spec["wavenumber"]) > 5


def test_solar_spectrum_integrates_to_tsi():
    """The solar spectrum should integrate to approximately 1361 W/m^2."""
    from climt._components.cork.optics.stellar import load_stellar_spectrum

    spec = load_stellar_spectrum("sun")
    total = np.trapz(spec["irradiance"], spec["wavenumber"])
    np.testing.assert_allclose(total, 1361.0, rtol=0.02)


def test_integrate_spectrum_over_bands():
    """Band integration distributes total flux across bands."""
    from climt._components.cork.optics.stellar import (
        integrate_spectrum_over_bands,
        load_stellar_spectrum,
    )

    spec = load_stellar_spectrum("sun")
    # 3 equal-width bands spanning the spectrum
    band_limits = np.array([
        [100.0, 10000.0],
        [10000.0, 25000.0],
        [25000.0, 50000.0],
    ])
    flux_per_band = integrate_spectrum_over_bands(spec, band_limits)
    assert flux_per_band.shape == (3,)
    assert np.all(flux_per_band >= 0)
    # Sum should be close to TSI (bands cover entire range)
    np.testing.assert_allclose(flux_per_band.sum(), 1361.0, rtol=0.05)


def test_integrate_spectrum_narrow_band():
    """A very narrow band should get a small fraction of total flux."""
    from climt._components.cork.optics.stellar import (
        integrate_spectrum_over_bands,
        load_stellar_spectrum,
    )

    spec = load_stellar_spectrum("sun")
    band_limits = np.array([[14000.0, 14100.0]])  # 100 cm^-1 narrow band
    flux = integrate_spectrum_over_bands(spec, band_limits)
    assert flux[0] > 0
    assert flux[0] < 100.0  # much less than TSI
```

- [x] **Step 3: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_stellar.py -v`
Expected: FAIL (`ImportError: cannot import name 'load_stellar_spectrum'`)

- [x] **Step 4: Implement the stellar spectrum module**

```python
# climt/_components/cork/optics/stellar.py
import os

import importlib_resources
import numpy as np


def load_stellar_spectrum(name_or_path):
    """Load a stellar spectrum from a built-in name or file path.

    Args:
        name_or_path: "sun", "trappist1", or path to .npz file

    Returns:
        dict with keys "wavenumber" (cm^-1) and "irradiance" (W/m^2/cm^-1)
    """
    if os.path.isfile(name_or_path):
        data = np.load(name_or_path)
    else:
        data_path = importlib_resources.files(
            "climt._data.cork.stellar_spectra"
        ).joinpath(f"{name_or_path}.npz")
        with importlib_resources.as_file(data_path) as f:
            data = np.load(f)

    return {
        "wavenumber": np.array(data["wavenumber"]),
        "irradiance": np.array(data["irradiance"]),
    }


def integrate_spectrum_over_bands(spectrum, band_wavenumber_limits):
    """Integrate a stellar spectrum over spectral bands.

    Args:
        spectrum: dict with "wavenumber" (cm^-1, ascending) and
                  "irradiance" (W/m^2/cm^-1)
        band_wavenumber_limits: (nband, 2) lower and upper wavenumber limits, cm^-1

    Returns:
        flux_per_band: (nband,) integrated flux per band, W/m^2
    """
    wn = spectrum["wavenumber"]
    irr = spectrum["irradiance"]
    nband = band_wavenumber_limits.shape[0]
    flux = np.zeros(nband)

    for b in range(nband):
        wn_lo, wn_hi = band_wavenumber_limits[b]
        mask = (wn >= wn_lo) & (wn <= wn_hi)
        if mask.sum() >= 2:
            flux[b] = np.trapz(irr[mask], wn[mask])
        elif mask.sum() == 1:
            # Single point in band: approximate with rectangle
            idx = np.where(mask)[0][0]
            flux[b] = irr[idx] * (wn_hi - wn_lo)
        else:
            # No spectral points in band: interpolate at band edges
            irr_lo = np.interp(wn_lo, wn, irr)
            irr_hi = np.interp(wn_hi, wn, irr)
            flux[b] = 0.5 * (irr_lo + irr_hi) * (wn_hi - wn_lo)

    return flux
```

- [x] **Step 5: Run stellar spectrum tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_stellar.py -v`
Expected: All PASS

- [x] **Step 6: Commit**

```bash
git add climt/_components/cork/optics/stellar.py climt/_data/cork/stellar_spectra/sun.npz tests/test_cork_stellar.py
git commit -m "feat(cork): add stellar spectrum loading and band integration"
```

---

### Task 4: Create SW correlated-k test table

**Files:**
- Create: `climt/_data/cork/correlated_k/test_2band_sw.npz`
- Modify: `tests/test_cork_optics.py`

The correlated-k SW path needs a test table that includes SW-specific fields: `solar_source_per_gpoint` and `rayleigh_coefficient`.

- [x] **Step 1: Write a test that loads the SW table**

Add to `tests/test_cork_optics.py`:

```python
def test_correlated_k_sw_table_loads():
    """Loading the SW test table returns expected dimensions and SW-specific fields."""
    from climt._components.cork.optics.correlated_k import load_k_table

    table = load_k_table("test_2band_sw")
    assert table["k_coefficients"].shape[0] == 1  # 1 gas
    assert table["k_coefficients"].shape[1] == 2  # 2 bands
    assert table["k_coefficients"].shape[2] == 2  # 2 g-points
    assert "solar_source_per_gpoint" in table
    assert table["solar_source_per_gpoint"].shape == (2, 2)  # (nband, ngpt)
    assert "rayleigh_coefficient" in table
    assert table["rayleigh_coefficient"].shape == (2,)  # (nband,)
    # Solar source should sum to something reasonable per band
    for b in range(2):
        assert table["solar_source_per_gpoint"][b].sum() > 0
```

- [x] **Step 2: Run test to verify it fails**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_optics.py::test_correlated_k_sw_table_loads -v`
Expected: FAIL (file not found)

- [x] **Step 3: Create the SW test table**

```python
# Script to generate test_2band_sw.npz — run once from repo root
import numpy as np

ngas = 1
nband = 2
ngpt = 2
nT = 5
nP = 5

# Absorption coefficients: moderate values for a single gas
k_data = np.ones((ngas, nband, ngpt, nT, nP)) * 0.005  # m^2/kg
# Band 1 (visible): slightly more absorbing at high g-point
k_data[0, 0, 1, :, :] = 0.01
# Band 2 (near-IR): more absorbing
k_data[0, 1, :, :, :] = 0.02
k_data[0, 1, 1, :, :] = 0.05

np.savez(
    "climt/_data/cork/correlated_k/test_2band_sw.npz",
    k_coefficients=k_data,
    gpoint_weights=np.array([[0.6, 0.4], [0.5, 0.5]]),
    band_wavenumber_limits=np.array([[10000.0, 25000.0], [4000.0, 10000.0]]),
    temperature_grid=np.linspace(200.0, 400.0, nT),
    pressure_grid_log=np.linspace(np.log(100.0), np.log(1e5), nP),
    gas_names=np.array(["h2o"]),
    overlap_method=np.array("additive"),
    resolution=np.array("low"),
    # SW-specific fields:
    solar_source_per_gpoint=np.array([[300.0, 200.0], [250.0, 250.0]]),  # W/m^2
    rayleigh_coefficient=np.array([1e-6, 5e-7]),  # per band
)
```

- [x] **Step 4: Run the test**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_optics.py::test_correlated_k_sw_table_loads -v`
Expected: PASS

- [x] **Step 5: Commit**

```bash
git add climt/_data/cork/correlated_k/test_2band_sw.npz tests/test_cork_optics.py
git commit -m "test(cork): add SW correlated-k test table"
```

---

### Task 5: CorkShortwaveRadiation correlated-k mode

**Files:**
- Modify: `climt/_components/cork/sw/component.py`
- Modify: `tests/test_cork_sw.py`

Enable `CorkShortwaveRadiation(optics="correlated_k", table="test_2band_sw")`. This path loads the SW k-table, computes absorption optical depths per band/g-point, adds Rayleigh scattering, and passes everything to the two-stream solver. Stellar flux per g-point comes from the table's `solar_source_per_gpoint` field.

- [x] **Step 1: Write failing tests for correlated-k SW**

Add to `tests/test_cork_sw.py`:

```python
def test_cork_sw_correlated_k_runs():
    """CorkShortwaveRadiation in correlated-k mode produces non-zero fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkShortwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    sw = CorkShortwaveRadiation(optics="correlated_k", table="test_2band_sw")
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)

    state["zenith_angle"].values[:] = np.pi / 4

    tendencies, diagnostics = sw(state)

    assert np.any(diagnostics["downwelling_shortwave_flux_in_air"].values != 0.0)
    assert not np.any(np.isnan(tendencies["air_temperature"].values))


def test_cork_sw_correlated_k_band_sum():
    """Per-band fluxes sum to broadband in correlated-k SW mode."""
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkShortwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    sw = CorkShortwaveRadiation(optics="correlated_k", table="test_2band_sw")
    grid = get_grid(nx=1, ny=1, nz=20)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 3

    _, diagnostics = sw(state)

    down_band = diagnostics["downwelling_shortwave_flux_in_air_per_band"].values
    down_broad = diagnostics["downwelling_shortwave_flux_in_air"].values
    np.testing.assert_allclose(down_band.sum(axis=-1), down_broad, rtol=1e-10)


def test_cork_sw_correlated_k_nighttime():
    """Correlated-k SW produces zero flux at nighttime."""
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkShortwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    sw = CorkShortwaveRadiation(optics="correlated_k", table="test_2band_sw")
    grid = get_grid(nx=1, ny=1, nz=15)
    state = get_default_state([sw], grid_state=grid)
    state["zenith_angle"].values[:] = np.pi / 2 + 0.1

    tendencies, diagnostics = sw(state)

    np.testing.assert_allclose(
        diagnostics["downwelling_shortwave_flux_in_air"].values, 0.0, atol=1e-10
    )
```

- [x] **Step 2: Run tests to verify they fail**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_sw.py::test_cork_sw_correlated_k_runs -v`
Expected: FAIL (`NotImplementedError: Correlated-k SW mode not yet implemented`)

- [x] **Step 3: Implement correlated-k mode in the SW component**

Update `climt/_components/cork/sw/component.py`:

```python
# climt/_components/cork/sw/component.py
import numpy as np
from sympl import TendencyComponent, get_constant

from ..common import compute_column_amount, compute_heating_rate, njit, prange
from ..optics.parmentier import (
    compute_rosseland_mean_opacity,
    load_freedman2014_coefficients,
    load_parmentier_coefficients,
    lookup_ratio_coefficients,
)
from .kernels import sw_two_stream


class CorkShortwaveRadiation(TendencyComponent):
    def __init__(
        self,
        optics="parmentier",
        table=None,
        coefficients="solar_composition",
        stellar_spectrum="sun",
        rosseland_mean_fit="freedman2014",
        **kwargs,
    ):
        self._optics_mode = optics

        if optics == "parmentier":
            self._coefficients = load_parmentier_coefficients(coefficients)
            self._freedman_coeffs = load_freedman2014_coefficients()
            self._num_bands = 3
            self._default_solar_flux_per_band = np.array([1361.0 / 3.0] * 3)
        elif optics == "correlated_k":
            from ..optics.correlated_k import load_k_table

            self._table = load_k_table(table)
            self._num_bands = self._table["k_coefficients"].shape[1]
            self._num_gpts = self._table["k_coefficients"].shape[2]
            self._gas_names = list(self._table["gas_names"])
            self._solar_source = self._table["solar_source_per_gpoint"]
            self._rayleigh = self._table.get("rayleigh_coefficient", None)
        else:
            raise ValueError(f"Unknown optics mode: {optics}")

        from climt._core.initialization import set_num_shortwave_bands
        set_num_shortwave_bands(self._num_bands)
        super(CorkShortwaveRadiation, self).__init__(**kwargs)

    @property
    def input_properties(self):
        props = {
            "air_temperature": {
                "dims": ["mid_levels", "*"],
                "units": "degK",
                "alias": "T",
            },
            "air_pressure": {"dims": ["mid_levels", "*"], "units": "Pa", "alias": "p"},
            "air_pressure_on_interface_levels": {
                "dims": ["interface_levels", "*"],
                "units": "Pa",
                "alias": "p_int",
            },
            "surface_temperature": {"dims": ["*"], "units": "degK", "alias": "T_surf"},
            "zenith_angle": {"dims": ["*"], "units": "radians", "alias": "zenith"},
            "surface_albedo_for_direct_shortwave": {
                "dims": ["*"],
                "units": "dimensionless",
                "alias": "albedo",
            },
            "flux_adjustment_for_earth_sun_distance": {
                "dims": ["*"],
                "units": "dimensionless",
                "alias": "earth_sun_factor",
            },
        }
        if self._optics_mode == "parmentier":
            props["irradiation_temperature"] = {
                "dims": ["*"],
                "units": "degK",
                "alias": "T_irr",
            }
            props["internal_temperature"] = {
                "dims": ["*"],
                "units": "degK",
                "alias": "T_int",
            }
        elif self._optics_mode == "correlated_k":
            _GAS_CF_NAME = {
                "h2o": "specific_humidity",
                "co2": "mole_fraction_of_carbon_dioxide_in_air",
            }
            _GAS_UNITS = {"h2o": "kg/kg"}
            for gas in self._gas_names:
                cf_name = _GAS_CF_NAME.get(gas, f"mole_fraction_of_{gas}_in_air")
                units = _GAS_UNITS.get(gas, "mole/mole")
                props[cf_name] = {
                    "dims": ["mid_levels", "*"],
                    "units": units,
                    "alias": gas,
                }
        return props

    @property
    def tendency_properties(self):
        return {"air_temperature": {"units": "degK s^-1"}}

    @property
    def diagnostic_properties(self):
        return {
            "upwelling_shortwave_flux_in_air": {
                "dims": ["interface_levels", "*"],
                "units": "W m^-2",
            },
            "downwelling_shortwave_flux_in_air": {
                "dims": ["interface_levels", "*"],
                "units": "W m^-2",
            },
            "upwelling_shortwave_flux_in_air_per_band": {
                "dims": ["interface_levels", "*", "num_shortwave_bands"],
                "units": "W m^-2",
            },
            "downwelling_shortwave_flux_in_air_per_band": {
                "dims": ["interface_levels", "*", "num_shortwave_bands"],
                "units": "W m^-2",
            },
            "shortwave_heating_rate": {
                "dims": ["mid_levels", "*"],
                "units": "degK day^-1",
            },
        }

    @property
    def num_shortwave_bands(self):
        return self._num_bands

    def array_call(self, state):
        T = state["T"]
        p = state["p"]
        p_int = state["p_int"]

        orig_shape_T = T.shape
        orig_shape_pint = p_int.shape
        nlev = T.shape[0]

        T_flat = T.reshape(nlev, -1)
        p_flat = p.reshape(nlev, -1)
        p_int_flat = p_int.reshape(nlev + 1, -1)
        zenith_flat = state["zenith"].reshape(-1)
        albedo_flat = state["albedo"].reshape(-1)
        ncol = T_flat.shape[1]

        sigma = get_constant("stefan_boltzmann_constant", "W/m^2/K^4")
        g_const = get_constant("gravitational_acceleration", "m/s^2")
        cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")

        if self._optics_mode == "parmentier":
            nband = self._num_bands
            ngpt = 1
            T_irr_flat = state["T_irr"].reshape(-1)
            T_int_flat = state["T_int"].reshape(-1)
            tau, ssa, asym = self._parmentier_sw_optics(
                T_flat, p_flat, p_int_flat, T_irr_flat, T_int_flat, g_const
            )
            T_irr_max = T_irr_flat.max()
            if T_irr_max > 0:
                F0 = sigma * T_irr_max**4
                solar_flux_per_band = np.array([F0 / 3.0] * 3)
            else:
                solar_flux_per_band = self._default_solar_flux_per_band

            earth_sun_factor = float(state["earth_sun_factor"].reshape(-1)[0])
            solar_flux = solar_flux_per_band.reshape(nband, 1) * np.ones((nband, ngpt)) * earth_sun_factor
            weights = np.ones((nband, ngpt))

        elif self._optics_mode == "correlated_k":
            tau, ssa, asym, weights = self._correlated_k_sw_optics(
                T_flat, p_flat, p_int_flat, g_const
            )
            nband = tau.shape[0]
            ngpt = tau.shape[1]

            earth_sun_factor = float(state["earth_sun_factor"].reshape(-1)[0])
            # Solar flux per g-point from table, scaled by earth-sun factor
            # For ESFT, solar_source may need to be expanded to combined g-points
            if ngpt == self._solar_source.shape[1]:
                solar_flux = self._solar_source * earth_sun_factor
            else:
                # ESFT: replicate the solar source across combined g-points
                ngpt_orig = self._solar_source.shape[1]
                ngas = len(self._gas_names)
                solar_flux = np.zeros((nband, ngpt))
                for b in range(nband):
                    for idx in range(ngpt):
                        g_idx = idx % ngpt_orig
                        solar_flux[b, idx] = self._solar_source[b, g_idx] * earth_sun_factor

        up_band, down_band, up_broad, down_broad = sw_two_stream(
            tau, ssa, asym, zenith_flat, albedo_flat, solar_flux, weights
        )

        net_flux = up_broad - down_broad
        heating_rate = compute_heating_rate(net_flux, p_int_flat, g_const, cpd)

        tendency = heating_rate.reshape(orig_shape_T)
        up_broad_out = up_broad.reshape(orig_shape_pint)
        down_broad_out = down_broad.reshape(orig_shape_pint)
        heating_kday = heating_rate.reshape(orig_shape_T) * 86400.0

        up_band_out = np.moveaxis(up_band, 0, -1).reshape(orig_shape_pint + (nband,))
        down_band_out = np.moveaxis(down_band, 0, -1).reshape(
            orig_shape_pint + (nband,)
        )

        return (
            {"T": tendency},
            {
                "upwelling_shortwave_flux_in_air": up_broad_out,
                "downwelling_shortwave_flux_in_air": down_broad_out,
                "upwelling_shortwave_flux_in_air_per_band": up_band_out,
                "downwelling_shortwave_flux_in_air_per_band": down_band_out,
                "shortwave_heating_rate": heating_kday,
            },
        )

    def _parmentier_sw_optics(self, T, p, p_int, T_irr, T_int, g):
        """Compute SW optical depths for Parmentier mode (3 visible bands)."""
        nlev, ncol = T.shape
        nband = 3
        ngpt = 1

        tau = np.zeros((nband, ngpt, nlev, ncol))
        ssa = np.zeros((nband, ngpt, nlev, ncol))  # pure absorption for Parmentier
        asym = np.zeros((nband, ngpt, nlev, ncol))

        for i in range(ncol):
            A_B = 0.0
            mu_star = 0.25
            T_eff = (T_int[i] ** 4 + (1.0 - A_B) * mu_star * T_irr[i] ** 4) ** 0.25
            T_eff = max(T_eff, 100.0)
            gv1, gv2, gv3, beta, gamma_P, R = lookup_ratio_coefficients(
                self._coefficients, T_eff
            )
            gamma_vs = [gv1, gv2, gv3]

            for k in range(nlev):
                kappa_R = compute_rosseland_mean_opacity(
                    T[k, i], p[k, i], self._freedman_coeffs
                )
                dp = abs(p_int[k + 1, i] - p_int[k, i])
                mass = dp / g

                for b in range(3):
                    kappa_v = gamma_vs[b] * kappa_R
                    tau[b, 0, k, i] = kappa_v * mass

        return tau, ssa, asym

    def _correlated_k_sw_optics(self, T, p, p_int, g):
        """Compute SW optical properties for correlated-k mode."""
        from ..optics.correlated_k import compute_ck_optical_depth

        nlev, ncol = T.shape

        _MOLAR_MASS = {
            "h2o": 18.015, "co2": 44.010, "o3": 47.998,
            "ch4": 16.043, "n2o": 44.013, "o2": 31.998,
        }
        _M_AIR = 28.970

        ngas = len(self._gas_names)
        gas_amounts = np.zeros((ngas, nlev, ncol))
        for ig, gas_name in enumerate(self._gas_names):
            from sympl import get_constant as _gc
            q_gas_flat = np.zeros((nlev, ncol))
            # Gas amounts are already extracted by sympl into state[alias]
            # We need the raw state — but array_call only gets the aliased dict.
            # The gas amount computation uses compute_column_amount.
            # For now, we replicate the LW pattern.

        # Actually, we need the state dict. Let's refactor array_call to pass gas amounts.
        # This is handled inside array_call — see the updated version above.
        raise NotImplementedError("Called internal method directly — use array_call")
```

Wait — the correlated-k SW path needs gas amounts which come from the state dict. The cleanest approach is to compute them inside `array_call` (same pattern as the LW component). Let me revise `array_call` to handle this:

In `array_call`, the correlated-k branch should be:

```python
        elif self._optics_mode == "correlated_k":
            _MOLAR_MASS = {
                "h2o": 18.015, "co2": 44.010, "o3": 47.998,
                "ch4": 16.043, "n2o": 44.013, "o2": 31.998,
            }
            _M_AIR = 28.970

            ngas = len(self._gas_names)
            gas_amounts = np.zeros((ngas, nlev, ncol))
            for ig, gas in enumerate(self._gas_names):
                q_gas_flat = state[gas].reshape(nlev, -1)
                if gas != "h2o":
                    M_gas = _MOLAR_MASS.get(gas, _M_AIR)
                    q_gas_flat = q_gas_flat * (M_gas / _M_AIR)
                gas_amounts[ig, :, :] = compute_column_amount(q_gas_flat, p_int_flat, g_const)

            from ..optics.correlated_k import compute_ck_optical_depth
            result = compute_ck_optical_depth(self._table, T_flat, p_flat, gas_amounts)
            if isinstance(result, tuple):
                tau_abs, weights = result
            else:
                tau_abs = result
                weights = self._table["gpoint_weights"]

            nband = tau_abs.shape[0]
            ngpt = tau_abs.shape[1]

            # Add Rayleigh scattering: increases tau and sets ssa > 0
            ssa = np.zeros((nband, ngpt, nlev, ncol))
            asym = np.zeros((nband, ngpt, nlev, ncol))
            tau = tau_abs.copy()

            if self._rayleigh is not None:
                for b in range(nband):
                    for k in range(nlev):
                        dp = abs(p_int_flat[k + 1, :] - p_int_flat[k, :])
                        # Rayleigh optical depth per layer = coefficient * column air mass
                        tau_ray = self._rayleigh[b] * dp / g_const
                        for gp in range(ngpt):
                            tau_total = tau_abs[b, gp, k, :] + tau_ray
                            ssa[b, gp, k, :] = np.where(
                                tau_total > 0,
                                tau_ray / tau_total,
                                0.0,
                            )
                            tau[b, gp, k, :] = tau_total

            earth_sun_factor = float(state["earth_sun_factor"].reshape(-1)[0])
            if ngpt == self._solar_source.shape[1]:
                solar_flux = self._solar_source * earth_sun_factor
            else:
                ngpt_orig = self._solar_source.shape[1]
                solar_flux = np.zeros((nband, ngpt))
                for b in range(nband):
                    for idx in range(ngpt):
                        g_idx = idx % ngpt_orig
                        solar_flux[b, idx] = self._solar_source[b, g_idx] * earth_sun_factor
```

So the full updated `sw/component.py` integrates both paths. The `_parmentier_sw_optics` method stays as-is. The `_correlated_k_sw_optics` method is eliminated — all logic lives in `array_call`.

- [x] **Step 4: Run all SW tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_sw.py -v`
Expected: All PASS (Parmentier tests still pass, new correlated-k tests pass)

- [x] **Step 5: Run full test suite**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_kernels.py tests/test_cork_optics.py tests/test_cork_lw.py tests/test_cork_sw.py tests/test_cork_integration.py tests/test_cork_stellar.py -v`
Expected: All PASS

- [x] **Step 6: Commit**

```bash
git add climt/_components/cork/sw/component.py tests/test_cork_sw.py
git commit -m "feat(cork): add correlated-k mode to CorkShortwaveRadiation"
```

---

### Task 6: Wire stellar spectrum into SW component

**Files:**
- Modify: `climt/_components/cork/sw/component.py`
- Modify: `tests/test_cork_sw.py`

Currently the Parmentier mode uses a hardcoded 1361/3 W/m^2 per band, and the correlated-k mode uses `solar_source_per_gpoint` from the table. Add an option to compute the Parmentier-mode solar flux from the loaded stellar spectrum integrated over equal-width visible bands, and use `integrate_spectrum_over_bands` for correlated-k tables that don't include `solar_source_per_gpoint`.

- [x] **Step 1: Write a test for stellar spectrum integration in the component**

Add to `tests/test_cork_sw.py`:

```python
def test_cork_sw_stellar_spectrum_loads():
    """CorkShortwaveRadiation with stellar_spectrum='sun' uses non-default solar flux."""
    from climt._components.cork import CorkShortwaveRadiation

    sw = CorkShortwaveRadiation(optics="parmentier", stellar_spectrum="sun")
    # The component should have loaded the spectrum and have non-trivial flux per band
    # (not exactly 1361/3 per band, since the spectrum is non-uniform)
    flux = sw._solar_flux_per_band
    assert flux.shape == (3,)
    assert np.all(flux > 0)
    # Total should be approximately 1361 W/m^2
    np.testing.assert_allclose(flux.sum(), 1361.0, rtol=0.05)
    # Bands should NOT be exactly equal (spectrum is non-uniform)
    assert not np.allclose(flux, flux[0], rtol=0.01)
```

- [x] **Step 2: Run test to verify it fails**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_sw.py::test_cork_sw_stellar_spectrum_loads -v`
Expected: FAIL (currently the stored flux is exactly 1361/3 per band)

- [x] **Step 3: Wire stellar spectrum loading into the SW component**

In `sw/component.py`, modify `__init__` to load the spectrum when `stellar_spectrum` is provided:

In the Parmentier `__init__` branch, replace:
```python
            self._default_solar_flux_per_band = np.array([1361.0 / 3.0] * 3)
```

With:
```python
            # Compute per-band solar flux from stellar spectrum.
            # Parmentier mode: 3 equal-width visible bands spanning the full spectrum.
            from ..optics.stellar import load_stellar_spectrum, integrate_spectrum_over_bands
            try:
                spec = load_stellar_spectrum(stellar_spectrum)
                wn = spec["wavenumber"]
                wn_lo, wn_hi = wn.min(), wn.max()
                band_width = (wn_hi - wn_lo) / 3.0
                band_limits = np.array([
                    [wn_lo, wn_lo + band_width],
                    [wn_lo + band_width, wn_lo + 2 * band_width],
                    [wn_lo + 2 * band_width, wn_hi],
                ])
                self._solar_flux_per_band = integrate_spectrum_over_bands(spec, band_limits)
            except (FileNotFoundError, KeyError):
                self._solar_flux_per_band = np.array([1361.0 / 3.0] * 3)
```

And update `array_call` to use `self._solar_flux_per_band` instead of `self._default_solar_flux_per_band` for the Earth fallback:

Replace:
```python
                solar_flux_per_band = self._default_solar_flux_per_band
```
With:
```python
                solar_flux_per_band = self._solar_flux_per_band
```

- [x] **Step 4: Run all SW tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_sw.py -v`
Expected: All PASS

- [x] **Step 5: Commit**

```bash
git add climt/_components/cork/sw/component.py tests/test_cork_sw.py
git commit -m "feat(cork): wire stellar spectrum loading into SW component"
```

---

### Task 7: Integration tests for phase 2 features

**Files:**
- Modify: `tests/test_cork_integration.py`

- [x] **Step 1: Write integration tests**

Add to `tests/test_cork_integration.py`:

```python
def test_correlated_k_lw_sw_combined():
    """Combined LW+SW in correlated-k mode produces physically sensible results."""
    import os
    import tempfile
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkLongwaveRadiation, CorkShortwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    # Build a minimal LW+SW table pair in a temp directory
    ngas, nband, ngpt, nT, nP = 1, 2, 2, 5, 5
    with tempfile.TemporaryDirectory() as tmpdir:
        lw_file = os.path.join(tmpdir, "lw.npz")
        sw_file = os.path.join(tmpdir, "sw.npz")

        k_data = np.ones((ngas, nband, ngpt, nT, nP)) * 0.01
        np.savez(
            lw_file,
            k_coefficients=k_data,
            gpoint_weights=np.full((nband, ngpt), 0.5),
            planck_fraction=np.full((nband, ngpt, nT), 0.5),
            band_wavenumber_limits=np.array([[200.0, 800.0], [800.0, 1400.0]]),
            temperature_grid=np.linspace(200.0, 400.0, nT),
            pressure_grid_log=np.linspace(np.log(100.0), np.log(1e5), nP),
            gas_names=np.array(["h2o"]),
            overlap_method=np.array("additive"),
            resolution=np.array("low"),
        )
        np.savez(
            sw_file,
            k_coefficients=k_data,
            gpoint_weights=np.full((nband, ngpt), 0.5),
            band_wavenumber_limits=np.array([[10000.0, 25000.0], [4000.0, 10000.0]]),
            temperature_grid=np.linspace(200.0, 400.0, nT),
            pressure_grid_log=np.linspace(np.log(100.0), np.log(1e5), nP),
            gas_names=np.array(["h2o"]),
            overlap_method=np.array("additive"),
            resolution=np.array("low"),
            solar_source_per_gpoint=np.array([[300.0, 200.0], [250.0, 250.0]]),
            rayleigh_coefficient=np.array([1e-6, 5e-7]),
        )

        lw = CorkLongwaveRadiation(optics="correlated_k", table=lw_file)
        sw = CorkShortwaveRadiation(optics="correlated_k", table=sw_file)
        grid = get_grid(nx=1, ny=1, nz=20)
        state = get_default_state([lw, sw], grid_state=grid)

        state["zenith_angle"].values[:] = np.pi / 4

        tend_lw, diag_lw = lw(state)
        tend_sw, diag_sw = sw(state)

        # LW OLR positive
        olr = diag_lw["upwelling_longwave_flux_in_air"].values[-1, 0, 0]
        assert olr > 0

        # SW surface flux positive
        sw_sfc = diag_sw["downwelling_shortwave_flux_in_air"].values[0, 0, 0]
        assert sw_sfc > 0

        # No NaNs
        assert not np.any(np.isnan(tend_lw["air_temperature"].values))
        assert not np.any(np.isnan(tend_sw["air_temperature"].values))


def test_sw_scattering_increases_upward_flux():
    """A scattering atmosphere should produce more upward SW flux than pure absorption."""
    from climt._components.cork.sw.kernels import sw_two_stream

    nlev = 10
    ncol = 1
    nband = 1
    ngpt = 1

    tau = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 4])
    albedo = np.array([0.0])  # black surface
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    # Pure absorption
    ssa_abs = np.zeros((nband, ngpt, nlev, ncol))
    asym_abs = np.zeros((nband, ngpt, nlev, ncol))
    up_abs, _, _, _ = sw_two_stream(
        tau, ssa_abs, asym_abs, zenith, albedo, solar_flux, weights
    )

    # With scattering (ssa=0.5)
    ssa_scat = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    asym_scat = np.zeros((nband, ngpt, nlev, ncol))
    up_scat, _, _, _ = sw_two_stream(
        tau, ssa_scat, asym_scat, zenith, albedo, solar_flux, weights
    )

    # Scattering atmosphere should reflect more radiation upward at TOA
    assert up_scat[nlev, 0] > up_abs[nlev, 0], (
        f"Scattering should increase upward flux: "
        f"absorbing={up_abs[nlev, 0]:.4f}, scattering={up_scat[nlev, 0]:.4f}"
    )
```

- [x] **Step 2: Run integration tests**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_integration.py -v`
Expected: All PASS

- [x] **Step 3: Run full test suite**

Run: `cd /Users/joymonteiro/github/climt && conda run -n climt python -m pytest tests/test_cork_kernels.py tests/test_cork_optics.py tests/test_cork_lw.py tests/test_cork_sw.py tests/test_cork_integration.py tests/test_cork_stellar.py -v`
Expected: All PASS

- [x] **Step 4: Commit**

```bash
git add tests/test_cork_integration.py
git commit -m "test(cork): add phase 2 integration tests for two-stream and correlated-k SW"
```

---

## Summary

| Task | Component | Tests |
|---|---|---|
| 1 | Meador & Weaver two-stream SW kernel | 3 kernel tests |
| 2 | ESFT gas overlap for correlated-k | 3 optics tests |
| 3 | Stellar spectrum loading | 4 stellar tests |
| 4 | SW correlated-k test table | 1 table test |
| 5 | CorkShortwaveRadiation correlated-k mode | 3 component tests |
| 6 | Stellar spectrum wiring into SW | 1 component test |
| 7 | Integration tests | 2 integration tests |

**Total: 7 tasks, ~17 tests**

### What's deferred (not in this plan):

- **Pre-built correlated-k tables for Earth/Mars/Venus/TRAPPIST**: Require external `linepyline` table generation. The infrastructure to load and use them is complete.
- **Cloud SW optical properties**: The `shortwave_optical_thickness_due_to_cloud`, `single_scattering_albedo_due_to_cloud`, and `cloud_asymmetry_parameter` inputs could be added to the SW component following the same pattern as LW cloud tau. Deferred until the two-stream solver is validated.
- **Bond albedo feedback for Parmentier mode**: The `A_B` parameter in T_eff calculation is currently hardcoded to 0. Could be computed self-consistently from SW results.
- **Additional stellar spectra**: Only `sun.npz` is provided. `trappist1.npz` and others can be added as data files.
- **Documentation**: Detailed pedagogical documentation as specified in Section 11 of the design spec.
