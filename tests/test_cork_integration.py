# tests/test_cork_integration.py
import numpy as np
import pytest
import sympl


def test_lw_sw_combined_energy_balance():
    """Combined LW+SW fluxes produce physically sensible results."""
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkLongwaveRadiation, CorkShortwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    lw = CorkLongwaveRadiation(optics="parmentier")
    sw = CorkShortwaveRadiation(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=30)
    state = get_default_state([lw, sw], grid_state=grid)

    state["zenith_angle"].values[:] = np.pi / 3  # daytime
    state["surface_albedo_for_direct_shortwave"].values[:] = 0.3
    state["irradiation_temperature"].values[:] = (
        300.0  # SW flux to ensure net heating at TOA
    )

    tend_lw, diag_lw = lw(state)
    tend_sw, diag_sw = sw(state)

    # OLR should be positive (upward LW at TOA)
    olr = diag_lw["upwelling_longwave_flux_in_air"].values[-1, 0, 0]
    assert olr > 0

    # Downward SW at surface should be positive
    sw_sfc = diag_sw["downwelling_shortwave_flux_in_air"].values[0, 0, 0]
    assert sw_sfc > 0

    # Net tendency should have cooling at some levels and heating at others
    total_tend = tend_lw["air_temperature"].values + tend_sw["air_temperature"].values
    assert np.any(total_tend > 0)  # some heating
    assert np.any(total_tend < 0)  # some cooling


def test_cloud_optical_depth_modifies_fluxes():
    """Adding cloud optical depth changes the LW fluxes."""
    from climt import get_default_state, get_grid
    from climt._components.cork import CorkLongwaveRadiation

    sympl.set_backend(sympl.DataArrayBackend())

    lw = CorkLongwaveRadiation(optics="parmentier")
    grid = get_grid(nx=1, ny=1, nz=20)
    state_clear = get_default_state([lw], grid_state=grid)

    _, diag_clear = lw(state_clear)
    olr_clear = diag_clear["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    # Now add cloud in the middle of the atmosphere
    state_cloudy = get_default_state([lw], grid_state=grid)
    assert "longwave_optical_thickness_due_to_cloud" in state_cloudy, (
        "Cloud tau not in default state — input_properties missing it"
    )
    cloud_tau = state_cloudy["longwave_optical_thickness_due_to_cloud"].values
    cloud_tau[10, ...] = 5.0  # thick cloud at level 10, all columns, all bands

    _, diag_cloudy = lw(state_cloudy)
    olr_cloudy = diag_cloudy["upwelling_longwave_flux_in_air"].values[-1, 0, 0]

    # Cloud should reduce OLR (more LW trapped)
    assert olr_cloudy < olr_clear, (
        f"Cloud should reduce OLR: clear={olr_clear:.2f}, cloudy={olr_cloudy:.2f}"
    )


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
    _, _, up_abs, _ = sw_two_stream(
        tau, ssa_abs, asym_abs, zenith, albedo, solar_flux, weights
    )

    # With scattering (ssa=0.5)
    ssa_scat = 0.5 * np.ones((nband, ngpt, nlev, ncol))
    asym_scat = np.zeros((nband, ngpt, nlev, ncol))
    _, _, up_scat, _ = sw_two_stream(
        tau, ssa_scat, asym_scat, zenith, albedo, solar_flux, weights
    )

    # Scattering atmosphere should reflect more radiation upward at TOA
    assert up_scat[nlev, 0] > up_abs[nlev, 0], (
        f"Scattering should increase upward flux: "
        f"absorbing={up_abs[nlev, 0]:.4f}, scattering={up_scat[nlev, 0]:.4f}"
    )
