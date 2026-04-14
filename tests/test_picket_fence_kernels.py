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


def test_sw_two_stream_no_atmosphere():
    """With zero optical depth, direct beam reaches surface unattenuated."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev = 10
    ncol = 1
    nband = 1
    ngpt = 1

    tau = np.zeros((nband, ngpt, nlev, ncol))
    ssa = np.zeros((nband, ngpt, nlev, ncol))
    asymmetry = np.zeros((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 3])  # 60 degrees
    albedo = np.array([0.0])
    solar_flux = np.array([[100.0]])  # (nband, ngpt) W/m^2
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, asymmetry, zenith, albedo, solar_flux, weights
    )

    mu0 = np.cos(zenith[0])
    # With no atmosphere, downward flux at surface = solar_flux * mu0
    np.testing.assert_allclose(down_broad[0, 0], 100.0 * mu0, rtol=1e-10)
    # Upward flux everywhere = 0 (zero albedo)
    np.testing.assert_allclose(up_broad[:, 0], 0.0, atol=1e-10)


def test_sw_two_stream_total_absorption():
    """With very large optical depth and zero scattering, no flux reaches surface."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev = 5
    ncol = 1
    nband = 1
    ngpt = 1

    tau = 100.0 * np.ones((nband, ngpt, nlev, ncol))
    ssa = np.zeros((nband, ngpt, nlev, ncol))  # pure absorption
    asymmetry = np.zeros((nband, ngpt, nlev, ncol))
    zenith = np.array([np.pi / 4])
    albedo = np.array([0.3])
    solar_flux = np.array([[1000.0]])
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, asymmetry, zenith, albedo, solar_flux, weights
    )

    # Nearly all radiation absorbed before reaching surface
    np.testing.assert_allclose(down_broad[0, 0], 0.0, atol=1e-6)


def test_sw_reflected_beam_attenuates_upward():
    """In a pure-absorption atmosphere, reflected photons are absorbed going back up."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev = 5
    ncol = 1
    nband = 1
    ngpt = 1

    # Moderate absorption (tau_total = 2.5 per layer path) with overhead sun
    tau = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    ssa = np.zeros((nband, ngpt, nlev, ncol))  # pure absorption
    asymmetry = np.zeros((nband, ngpt, nlev, ncol))
    zenith = np.array([0.0])  # overhead sun (mu0=1)
    albedo = np.array([1.0])  # perfect reflector: makes the effect maximal
    solar_flux = np.array([[100.0]])
    weights = np.ones((nband, ngpt))

    up_band, down_band, up_broad, down_broad = sw_two_stream(
        tau, ssa, asymmetry, zenith, albedo, solar_flux, weights
    )

    # Reflected beam is strongest at the surface and must attenuate going upward
    # because the atmosphere absorbs it (ssa=0).
    assert up_broad[0, 0] > up_broad[nlev, 0], (
        f"Reflected beam must attenuate upward: "
        f"up[surface]={up_broad[0, 0]:.4f}, up[TOA]={up_broad[nlev, 0]:.4f}"
    )
    # In the two-stream framework, the reflected flux at the surface is diffuse
    # and propagates upward via Tdif (which uses a ~2x diffusivity factor, not mu0).
    # Check that TOA upward flux is much less than surface upward flux.
    assert up_broad[nlev, 0] < 0.1 * up_broad[0, 0], (
        f"Reflected beam should be heavily attenuated: "
        f"up[TOA]={up_broad[nlev, 0]:.4f}, up[surface]={up_broad[0, 0]:.4f}"
    )


def test_sw_two_stream_conservative_scattering():
    """With ssa=1 and no absorption, total flux (up+down) is conserved at every level."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

    nlev = 10
    ncol = 1
    nband = 1
    ngpt = 1

    tau = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    ssa = np.ones((nband, ngpt, nlev, ncol))  # pure scattering
    asymmetry = np.zeros((nband, ngpt, nlev, ncol))  # isotropic
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
    from climt._components.picket_fence.sw.kernels import sw_two_stream

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
    # Delta-scaling with g=0.999 greatly reduces tau but also reduces ssa,
    # so some absorption remains.  Check that most flux reaches the surface (>80%).
    assert down_broad[0, 0] > 0.8 * expected_surface, (
        f"Forward scattering should transmit most flux: got {down_broad[0, 0]:.2f}, "
        f"expected ~{expected_surface:.2f}"
    )


def test_sw_two_stream_albedo_reflection():
    """With zero atmosphere, upward flux equals albedo * downward flux at surface."""
    from climt._components.picket_fence.sw.kernels import sw_two_stream

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

    assert "Rdif" in diag
    assert "Tdif" in diag
    assert "Tnoscat" in diag
    assert "direct_beam_flux" in diag

    assert diag["Rdif"].shape == (nband, ngpt, nlev, ncol)
    assert diag["Tdif"].shape == (nband, ngpt, nlev, ncol)
    assert diag["Tnoscat"].shape == (nband, ngpt, nlev, ncol)
    assert diag["direct_beam_flux"].shape == (nband, ngpt, nlev + 1, ncol)

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

    assert "Rdir" in diag
    assert "Tdir" in diag
    assert "combined_albedo" in diag
    assert "tau_delta" in diag
    assert "ssa_delta" in diag
    assert "g_delta" in diag

    assert diag["Rdir"].shape == (nband, ngpt, nlev, ncol)
    assert diag["combined_albedo"].shape == (nband, ngpt, nlev + 1, ncol)
    assert diag["tau_delta"].shape == (nband, ngpt, nlev, ncol)


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

    assert np.all(diag["transmittance"] >= 0)
    assert np.all(diag["transmittance"] <= 1)


def test_sw_single_column_basic():
    """sw_single_column runs on 1D inputs and returns a dict."""
    from climt._components.picket_fence.sw.kernels import sw_single_column

    nlev = 10
    tau = 0.3 * np.ones(nlev)
    ssa = 0.5 * np.ones(nlev)
    g = 0.1 * np.ones(nlev)

    result = sw_single_column(
        tau, ssa, g,
        zenith=np.pi / 4,
        albedo=0.3,
        solar_flux=1361.0,
    )
    assert isinstance(result, dict)
    assert "flux_up" in result
    assert "flux_down" in result
    assert "flux_down_direct" in result
    assert "heating_rate" in result

    assert result["flux_up"].shape == (nlev + 1,)
    assert result["flux_down"].shape == (nlev + 1,)
    assert result["heating_rate"].shape == (nlev,)

    # With scattering atmosphere: some upward flux at TOA
    assert result["flux_up"][-1] > 0


def test_sw_single_column_diagnostics():
    """sw_single_column passes diagnostics_level through."""
    from climt._components.picket_fence.sw.kernels import sw_single_column

    nlev = 5
    result = sw_single_column(
        tau=0.3 * np.ones(nlev),
        ssa=0.5 * np.ones(nlev),
        g=0.1 * np.ones(nlev),
        zenith=np.pi / 4,
        albedo=0.3,
        solar_flux=1361.0,
        diagnostics_level=1,
    )
    assert "Rdif" in result
    assert result["Rdif"].shape == (nlev,)


def test_lw_single_column_basic():
    """lw_single_column runs on 1D inputs and returns a dict."""
    from climt._components.picket_fence.lw.kernels import lw_single_column

    nlev = 10
    T_layer = np.linspace(300, 200, nlev)
    tau = 0.5 * np.ones(nlev)

    result = lw_single_column(
        tau=tau,
        T_layer=T_layer,
        T_surface=300.0,
        emissivity=1.0,
    )
    assert isinstance(result, dict)
    assert "flux_up" in result
    assert "flux_down" in result
    assert "heating_rate" in result

    assert result["flux_up"].shape == (nlev + 1,)
    assert result["flux_down"].shape == (nlev + 1,)

    # Surface emits upward: flux_up at surface should be positive
    assert result["flux_up"][0] > 0


def test_lw_single_column_diagnostics():
    """lw_single_column passes diagnostics_level through."""
    from climt._components.picket_fence.lw.kernels import lw_single_column

    nlev = 5
    result = lw_single_column(
        tau=0.5 * np.ones(nlev),
        T_layer=np.linspace(300, 200, nlev),
        T_surface=300.0,
        diagnostics_level=1,
    )
    assert "transmittance" in result
    assert result["transmittance"].shape == (nlev,)
    assert np.all(result["transmittance"] >= 0)
    assert np.all(result["transmittance"] <= 1)
