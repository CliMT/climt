# tests/test_picket_fence_kernels.py
import numpy as np
import pytest


def test_lw_transport_isothermal_no_absorption():
    """An isothermal atmosphere with zero optical depth: fluxes equal sigma*T^4 everywhere."""
    from climt._components.picket_fence.lw.kernels import lw_transport

    nlev = 10
    ncol = 1
    nband = 1
    ngpt = 1
    sigma = 5.670374419e-8

    T = 250.0 * np.ones((nlev, ncol))
    T_surface = 250.0 * np.ones((ncol,))
    tau = np.zeros((nband, ngpt, nlev, ncol))
    planck_source = sigma * T[np.newaxis, np.newaxis, :, :] ** 4
    surface_source = sigma * T_surface[np.newaxis, np.newaxis, :] ** 4
    emissivity = np.ones((nband, ncol))
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = lw_transport(
        T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma
    )

    expected = sigma * 250.0**4
    # With zero optical depth, upward flux = surface emission everywhere,
    # downward flux = 0 everywhere (no absorption/emission by atmosphere)
    np.testing.assert_allclose(up_broad[:, 0], expected, rtol=1e-10)
    np.testing.assert_allclose(down_broad[:, 0], 0.0, atol=1e-10)


def test_lw_transport_opaque_atmosphere():
    """With very large optical depth, each layer emits as a blackbody at local T."""
    from climt._components.picket_fence.lw.kernels import lw_transport

    nlev = 5
    ncol = 1
    nband = 1
    ngpt = 1
    sigma = 5.670374419e-8

    T = np.linspace(200.0, 300.0, nlev).reshape(nlev, ncol)
    T_surface = np.array([310.0])
    # Very large optical depth per layer — each layer is opaque
    tau = 100.0 * np.ones((nband, ngpt, nlev, ncol))
    planck_source = sigma * T[np.newaxis, np.newaxis, :, :] ** 4
    surface_source = sigma * T_surface[np.newaxis, np.newaxis, :] ** 4
    emissivity = np.ones((nband, ncol))
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = lw_transport(
        T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma
    )

    # At TOA (index nlev), upward flux should be ~sigma * T_top^4
    expected_toa = sigma * T[-1, 0] ** 4
    np.testing.assert_allclose(up_broad[nlev, 0], expected_toa, rtol=0.01)

    # At surface (index 0), downward flux should be ~sigma * T_bottom^4
    expected_sfc = sigma * T[0, 0] ** 4
    np.testing.assert_allclose(down_broad[0, 0], expected_sfc, rtol=0.01)


def test_lw_transport_two_bands():
    """Two bands with different optical depths produce per-band fluxes that sum to broadband."""
    from climt._components.picket_fence.lw.kernels import lw_transport

    nlev = 10
    ncol = 2
    nband = 2
    ngpt = 1
    sigma = 5.670374419e-8

    T = 260.0 * np.ones((nlev, ncol))
    T_surface = 280.0 * np.ones((ncol,))
    tau = np.zeros((nband, ngpt, nlev, ncol))
    tau[0, 0, :, :] = 0.5  # band 0: moderate optical depth
    tau[1, 0, :, :] = 0.01  # band 1: nearly transparent
    # Equal Planck weighting: each band gets half the Planck function
    planck_source = 0.5 * sigma * T[np.newaxis, np.newaxis, :, :] ** 4
    planck_source = np.broadcast_to(planck_source, (nband, ngpt, nlev, ncol)).copy()
    surface_source = 0.5 * sigma * T_surface[np.newaxis, np.newaxis, :] ** 4
    surface_source = np.broadcast_to(surface_source, (nband, ngpt, ncol)).copy()
    emissivity = np.ones((nband, ncol))
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = lw_transport(
        T, T_surface, tau, planck_source, surface_source, emissivity, weights, sigma
    )

    # Per-band fluxes should sum to broadband
    np.testing.assert_allclose(
        up_band[0, :, :] + up_band[1, :, :], up_broad, rtol=1e-12
    )
    np.testing.assert_allclose(
        down_band[0, :, :] + down_band[1, :, :], down_broad, rtol=1e-12
    )
