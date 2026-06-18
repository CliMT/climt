"""Validation tests for cork radiation scheme.

These tests verify physics accuracy against analytical solutions.
They are slower than unit tests and focus on scientific correctness.

Tests call the kernel functions directly (lw_transport, sw_two_stream)
with single-column arrays reshaped to the (nband, ngpt, nlev, ncol) format.
"""
import numpy as np
import pytest


def _rrtmg_available():
    """Check if RRTMGLongwave/Shortwave can be imported."""
    try:
        from climt import RRTMGLongwave, RRTMGShortwave  # noqa: F401
        return True
    except (ImportError, OSError):
        return False


def _run_lw_single_col(tau_1d, T_layer_1d, T_surface, emissivity=1.0):
    """Helper: run LW transport on a single grey column.

    Reshapes 1D arrays to (nband=1, ngpt=1, nlev, ncol=1) and calls
    lw_transport directly.
    """
    from climt._components.cork.lw.kernels import lw_transport

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
    from climt._components.cork.sw.kernels import sw_two_stream

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

    def test_grey_isothermal_upward_flux_constant(self):
        """In an isothermal grey atmosphere, upward flux is constant at sigma*T^4.

        Each layer absorbs and re-emits the same amount it transmits, so the
        upward flux leaving each interface equals the surface emission sigma*T^4.
        The downward flux (bounded above by 0 from cold space) is NOT symmetric,
        but the upward flux should be exactly sigma*T^4 at every level.
        """
        nlev = 30
        T_iso = 250.0
        tau_per_layer = 2.0 / nlev * np.ones(nlev)

        result = _run_lw_single_col(
            tau_per_layer, T_iso * np.ones(nlev), T_iso, emissivity=1.0
        )

        sigma = 5.670374419e-8
        expected_up = sigma * T_iso ** 4
        np.testing.assert_allclose(result["flux_up"], expected_up, rtol=1e-10,
            err_msg="Upward flux must equal sigma*T^4 everywhere in isothermal atmosphere"
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
    """Compare cork against RRTMG for Earth standard atmosphere."""

    pass  # Placeholder: populated when RRTMG comparison is implemented
