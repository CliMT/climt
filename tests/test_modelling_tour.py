"""The science the Modelling Tour pages claim, checked natively.

Pyodide cells cannot run in CI, so each page's computational core lives in
``docs/modelling-tour/_tour/`` as importable Python and is exercised here. This
mirrors the pattern of ``tests/test_live_rce_demo.py``.

``docs/modelling-tour`` contains a hyphen and is not a valid package path, so
the helpers are loaded by file path.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import sympl

import climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid

REPO_ROOT = Path(__file__).resolve().parent.parent
TOUR = REPO_ROOT / "docs/modelling-tour/_tour"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, TOUR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def _unyt_backend():
    sympl.set_backend(climt.UnytBackend())


@pytest.fixture
def soundings():
    return _load("soundings")


def test_lapse_rate_sounding_shape_and_bounds(soundings):
    p = np.linspace(1e5, 1e3, 40)
    T, q = soundings.lapse_rate_sounding(p, 1e5, T_surf=288.0)
    assert T.shape == p.shape and q.shape == p.shape
    assert T[0] == pytest.approx(288.0, abs=0.5)     # bottom ~ surface
    assert T.min() == pytest.approx(200.0, abs=1e-6)  # stratosphere floor
    assert np.all(np.diff(T) <= 1e-9)                 # monotone with height
    assert np.all(q > 0) and q.max() < 0.05           # sane humidity


def test_saturation_vapour_pressure_matches_known_values(soundings):
    # Bolton (1980): ~611.2 Pa at 0 C, ~2339 Pa at 20 C
    assert soundings.saturation_vapour_pressure(273.15) == pytest.approx(611.2, rel=1e-3)
    assert soundings.saturation_vapour_pressure(293.15) == pytest.approx(2339.0, rel=2e-2)


def test_analytic_gray_equilibrium_reproduces_the_notes(soundings):
    """T_skin = 2^-0.25 T_e and T_g = T_e (1 + tau_inf/2)^0.25."""
    p = np.linspace(1e5, 1.0, 60)
    T, T_surf, tau_star = soundings.analytic_gray_equilibrium(p, 1e5, tau_inf=4.0,
                                                              T_e=255.0)
    assert T[-1] == pytest.approx(255.0 * 0.5 ** 0.25, rel=1e-3)
    assert T_surf == pytest.approx(255.0 * (1 + 4.0 / 2) ** 0.25, rel=1e-9)
    assert tau_star[0] == pytest.approx(4.0, rel=1e-9)   # deepest at the surface
    assert tau_star[-1] == pytest.approx(0.0, abs=1e-3)  # zero at TOA


@pytest.fixture
def spectra():
    return _load("spectra")


def test_brightness_temperature_inverts_planck(spectra):
    nu = np.array([300.0, 667.0, 900.0, 1400.0])
    for T in (200.0, 255.0, 288.0, 320.0):
        flux = spectra.planck_flux(nu, T)
        np.testing.assert_allclose(
            spectra.brightness_temperature(flux, nu), T, rtol=1e-8)


def test_per_band_olr_sums_to_the_broadband_diagnostic(spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=30))
    _, diag = lw(state)
    limits = spectra.band_limits_of(lw)
    olr_band, _ = spectra.spectral_olr(diag, limits)
    broadband = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert olr_band.sum() == pytest.approx(broadband, rel=1e-9)


def test_band_limits_cover_the_expected_features(spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    limits = spectra.band_limits_of(lw)
    assert limits.shape == (14, 2)
    # the 15 um CO2 core and the window edges must be exact band boundaries
    edges = set(np.round(limits.ravel(), 1))
    for edge in (630.0, 700.0, 800.0, 980.0, 1080.0, 1180.0):
        assert edge in edges


def _exponential_absorber(tau_inf=10.0, scale_height=8000.0, nz=200, n_scale=8):
    """Per-layer optical depth for a well-mixed absorber on a uniform-z grid.

    The notes' W ∝ tau* exp(-tau*) assumes layers of equal THICKNESS over an
    atmosphere whose density falls off exponentially, so each layer's optical
    depth grows towards the surface. Layers of equal optical depth are a
    different problem, and their weighting function is monotone.
    """
    z = np.linspace(0.0, n_scale * scale_height, nz + 1)  # edges, surface first
    above = tau_inf * np.exp(-z / scale_height)           # tau above each edge
    return above[:-1] - above[1:]


@pytest.mark.parametrize("D", [1.0, 1.66, 2.0])
def test_emission_weight_peaks_where_tau_star_is_one(spectra, D):
    """The notes' weighting function: emission * transmission peaks at tau*=1."""
    tau_layer = _exponential_absorber(tau_inf=10.0)
    w = spectra.emission_weight(tau_layer, diffusivity_factor=D)
    tau_star = spectra.tau_star_cumulative(tau_layer, diffusivity_factor=D)
    assert tau_star[-1] == pytest.approx(0.0, abs=0.06)     # ~0 at TOA
    assert tau_star[0] == pytest.approx(10.0 * D, rel=0.03)  # deepest at surface
    np.testing.assert_allclose(tau_star[np.argmax(w)], 1.0, atol=0.1)


def test_emission_weight_is_monotone_for_equal_optical_depth_layers(spectra):
    """Guards the trap above: equal-tau layers have no interior peak.

    Every layer emits the same amount and the ones higher up are less
    attenuated, so the weight rises monotonically to the top of the column.
    A "radiating level" needs a grid that resolves height, not optical depth.
    """
    w = spectra.emission_weight(np.full(200, 0.05), diffusivity_factor=1.0)
    assert np.all(np.diff(w) > 0)
    assert np.argmax(w) == len(w) - 1


@pytest.mark.slow
def test_page1_window_is_warmer_than_the_co2_core(soundings, spectra):
    """The page's central claim: brightness temperature swings ~200 K to ~285 K."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, q, T_surf=288.0)

    _, diag = lw(state)
    limits = spectra.band_limits_of(lw)
    _, flux_density = spectra.spectral_olr(diag, limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))

    centres = spectra.band_centres(limits)
    co2_core = tb[(centres > 630) & (centres < 700)]
    window = tb[(centres > 800) & (centres < 1180)]

    assert co2_core.max() < 230.0          # radiating from the stratosphere
    assert window.min() > 270.0            # radiating from near the surface
    assert window.min() < 288.0            # but not from the surface itself
    assert window.mean() - co2_core.mean() > 60.0


@pytest.mark.slow
def test_page2_window_is_transparent_and_co2_core_is_opaque(soundings, spectra):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, q, T_surf=288.0)
    _, diag = lw(state)

    tau = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]  # (nz, nband)
    column_tau = tau.sum(axis=0)
    centres = spectra.band_centres(spectra.band_limits_of(lw))

    co2_core = column_tau[(centres > 630) & (centres < 700)]
    window = column_tau[(centres > 800) & (centres < 1180)]
    assert co2_core.min() > 1.0                 # opaque
    assert window.max() < co2_core.min()        # window is the transparent part


@pytest.mark.slow
def test_page2_removing_absorbers_opens_the_window(soundings):
    """A near-N2 atmosphere: OLR should approach the surface blackbody flux."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)

    soundings.apply_sounding(state, T, q, T_surf=288.0)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 280e-6
    _, diag_earth = lw(state)
    olr_earth = float(diag_earth["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    soundings.apply_sounding(state, T, np.full_like(q, 1e-6), T_surf=288.0)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 10e-6
    _, diag_thin = lw(state)
    olr_thin = float(diag_thin["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    sigma_T4 = 5.670374419e-8 * 288.0 ** 4
    assert olr_thin > olr_earth + 100.0        # the window opened
    assert olr_thin > 0.85 * sigma_T4          # approaching a bare surface


PAGE3_D = 1.66


def _emission_pressure(spectra, tau_layer, p, D):
    """Pressure where tau* = 1, or nan if the column never reaches it.

    Interpolated in log tau* against log p, which is what page 3 does inline.
    """
    tau_star = spectra.tau_star_cumulative(tau_layer, D)
    if tau_star[0] < 1.0:
        return np.nan
    return float(np.exp(np.interp(0.0, np.log(tau_star[::-1]), np.log(p[::-1]))))


def _page3_column(soundings, spectra, humidity_scale=1.0, nz=60):
    """Page 3's column: the standard sounding, with a humidity knob."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=nz))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    soundings.apply_sounding(state, T, np.maximum(q * humidity_scale, 1e-7),
                             T_surf=288.0)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 280e-6
    _, diag = lw(state)

    limits = spectra.band_limits_of(lw)
    tau = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]
    p_emit = np.array([_emission_pressure(spectra, tau[:, b], p, PAGE3_D)
                       for b in range(tau.shape[1])])
    return dict(lw=lw, p=p, T=T, diag=diag, limits=limits, tau=tau, p_emit=p_emit)


@pytest.mark.slow
def test_page3_radiating_levels_are_a_distribution_not_a_height(soundings, spectra):
    """The page's central claim: there is no single radiating level."""
    col = _page3_column(soundings, spectra)
    limits, p_emit = col["limits"], col["p_emit"]
    centres = spectra.band_centres(limits)

    # The CO2 core emits from the stratosphere.
    core = p_emit[(centres > 630) & (centres < 700)]
    assert np.all(core < 2000.0)                    # Pa, i.e. above ~20 hPa

    # The clearest window bands never reach tau* = 1 at all: no radiating level
    # exists for them, the surface is simply visible.
    window = p_emit[(centres > 800) & (centres < 1080)]
    assert np.all(np.isnan(window))

    # Across the bands that do have one, the levels span the whole atmosphere.
    have = p_emit[np.isfinite(p_emit)]
    assert have.max() / have.min() > 100.0


@pytest.mark.slow
def test_page3_brightness_temperature_tracks_the_emission_level(soundings, spectra):
    """T_b(band) follows the emission-weighted temperature -- with a cold bias.

    The band-mean optical depth cork reports is the g-weighted mean over eight
    g-points whose k values span three to eight orders of magnitude. It is
    therefore dominated by the band's *strongest* absorption, while the OLR
    escapes preferentially through the band's weakest. So the band-mean
    weighting function always places the emission too high, and T_b comes out
    warmer than the weighted temperature -- never colder. That one-sidedness is
    the page's punchline, so it is what gets asserted.
    """
    col = _page3_column(soundings, spectra)
    T, tau, limits = col["T"], col["tau"], col["limits"]
    _, flux_density = spectra.spectral_olr(col["diag"], limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))

    t_weighted = np.array([
        float((spectra.emission_weight(tau[:, b], PAGE3_D) * T).sum()
              / spectra.emission_weight(tau[:, b], PAGE3_D).sum())
        for b in range(tau.shape[1])])

    sane = tb <= 288.0            # see the band-centre test below
    bias = (tb - t_weighted)[sane]
    assert np.all(bias > 0.0)     # always warmer, never colder
    assert bias.max() < 45.0
    assert np.corrcoef(t_weighted[sane], tb[sane])[0, 1] > 0.9

    # The page plots T at the tau* = 1 level rather than the weighted mean, so
    # guard that comparison too. Same story, with one documented exception:
    # 1080-1180 cm^-1 has column tau = 0.68, so its level lands inside the
    # bottom model layer and the construction degenerates into "the ground".
    p_emit = col["p_emit"]
    has_level = np.isfinite(p_emit) & sane
    t_level = np.exp(np.interp(np.log(p_emit[has_level]),
                               np.log(col["p"][::-1]), np.log(T[::-1])))
    level_bias = tb[has_level] - t_level
    assert (level_bias > 0.0).sum() == level_bias.size - 1
    assert level_bias.min() > -5.0
    assert level_bias.max() < 45.0
    assert np.corrcoef(t_level, tb[has_level])[0, 1] > 0.9


@pytest.mark.slow
def test_page3_band_centre_brightness_temperature_fails_on_the_widest_band(
        soundings, spectra):
    """1800-3250 cm^-1 returns T_b above the surface -- an inversion artifact.

    Inverting a band-mean flux density at the band *centre* only works if the
    Planck function is near-linear across the band. Over 1450 cm^-1 of shortwave
    tail it is not, and the answer comes back hotter than the ground, which is
    impossible. The page says so rather than plotting it silently.
    """
    col = _page3_column(soundings, spectra)
    limits = col["limits"]
    _, flux_density = spectra.spectral_olr(col["diag"], limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))

    bogus = np.flatnonzero(tb > 288.0)
    assert bogus.tolist() == [len(tb) - 1]
    np.testing.assert_allclose(limits[-1], [1800.0, 3250.0])


@pytest.mark.slow
def test_page3_moistening_raises_the_emission_levels(soundings, spectra):
    """The knob: a wetter atmosphere emits from higher, colder air."""
    dry = _page3_column(soundings, spectra, humidity_scale=1.0)
    wet = _page3_column(soundings, spectra, humidity_scale=4.0)
    centres = spectra.band_centres(dry["limits"])

    # Every band that had a radiating level keeps it, and moves up.
    had = np.isfinite(dry["p_emit"])
    assert np.all(np.isfinite(wet["p_emit"][had]))
    assert np.all(wet["p_emit"][had] <= dry["p_emit"][had] + 1.0)

    # The water vapour bands move a long way; the CO2 core barely moves.
    rotational = (centres > 10) & (centres < 250)
    assert (dry["p_emit"][rotational] / wet["p_emit"][rotational]).min() > 1.5
    core = (centres > 630) & (centres < 700)
    assert (dry["p_emit"][core] / wet["p_emit"][core]).max() < 1.2

    # And the window closes: bands with no radiating level acquire one.
    assert np.any(np.isfinite(wet["p_emit"][~had]))

    olr_dry = float(dry["diag"]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    olr_wet = float(wet["diag"]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert olr_wet < olr_dry - 20.0
