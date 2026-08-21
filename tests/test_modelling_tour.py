"""The science the Modelling Tour pages claim, checked natively.

Pyodide cells cannot run in CI, so each page's computational core lives in
``docs/modelling-tour/_tour/`` as importable Python and is exercised here.

``docs/modelling-tour`` contains a hyphen and is not a valid package path, so
the helpers are loaded by file path.
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pytest
import sympl

import climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid

REPO_ROOT = Path(__file__).resolve().parent.parent
TOUR = REPO_ROOT / "docs/modelling-tour/_tour"

SOLAR = 240.0     # prescribed absorbed shortwave at the surface, W m^-2, on
                  # every tranche 2 page. There is no shortwave component.


# The pages do `sys.path.insert(0, "_tour")` and then import helpers by bare
# name; _tour modules import each other the same way. Match that here so a
# module loaded by path can still `import assets`.
if str(TOUR) not in sys.path:
    sys.path.insert(0, str(TOUR))


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


def test_lapse_rate_sounding_stratosphere_gradient(soundings):
    """``gamma_strat`` warms the stratosphere with height, page 5's knob.

    Default is 0.0, which must leave the profile bit-identical to the isothermal
    one -- page 5 quotes numbers from both, and only the difference between them
    means anything.
    """
    p = np.linspace(1e5, 1e2, 60)
    T_flat, q_flat = soundings.lapse_rate_sounding(p, 1e5, T_surf=288.0)
    T_zero, _ = soundings.lapse_rate_sounding(p, 1e5, T_surf=288.0,
                                              gamma_strat=0.0)
    np.testing.assert_array_equal(T_flat, T_zero)

    EPSILON = soundings.EPSILON
    T_warm, q_warm = soundings.lapse_rate_sounding(p, 1e5, T_surf=288.0,
                                                   gamma_strat=2.5e-3)
    strat = T_flat <= 200.0 + 1e-9        # where the isothermal cap had bitten
    assert np.all(np.diff(T_warm[strat]) > 0)      # warms upward, monotonically
    assert T_warm[-1] > T_flat[-1] + 20.0          # and by a lot at the lid

    # The troposphere is untouched, humidity included.
    np.testing.assert_allclose(T_warm[~strat], T_flat[~strat])
    np.testing.assert_allclose(q_warm[~strat], q_flat[~strat])

    # Fixed RH over warm, thin air would otherwise return q of order 100 kg/kg:
    # at 1 hPa the saturation vapour pressure alone exceeds the total pressure.
    assert np.all(q_warm > 0) and q_warm.max() <= EPSILON

    # A cooling stratosphere works too, and must not cross the sign convention.
    T_cold, _ = soundings.lapse_rate_sounding(p, 1e5, T_surf=288.0,
                                              gamma_strat=-1e-3)
    assert np.all(np.diff(T_cold[strat]) < 0)


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


# Pages 1-3 call ``tables.spectrum_table()`` and run on the 56-band table; pages
# 4-6 name the shipped 14-band one directly. Testing pages 1-3 on both is the
# point: the prose quotes numbers from the 56-band run, and the 14-band run is
# the documented fallback a reader gets if the asset is missing, so neither path
# may go unguarded. Resolved at import, because parametrize needs it at
# collection time.
LOW_RES = "earth_low_res_lw"
_TABLES = _load("tables")
HI_RES = _TABLES.spectrum_table()


@pytest.fixture
def tables():
    return _TABLES


PAGE123_TABLES = [
    pytest.param(LOW_RES, id="14band"),
    pytest.param(HI_RES, id="56band",
                 marks=pytest.mark.skipif(
                     HI_RES == _TABLES.FALLBACK,
                     reason="the 56-band spectrum table asset is not present")),
]


def test_spectrum_table_finds_the_asset_from_any_working_directory(
        tables, tmp_path, monkeypatch):
    """The switch that decides whether pages 1-3 get 56 bands or silently 14.

    It resolves by ``os.path.isfile`` on two candidates, the second of which is
    anchored to this file rather than the caller's cwd -- that is what makes it
    work in a test run and a native render, where the working directory is not
    the page's. A regression here degrades all three pages' resolution and does
    so *silently*, by design, so it is worth pinning.
    """
    monkeypatch.chdir(tmp_path)          # nowhere near the docs tree
    table = tables.spectrum_table()
    assert table != tables.FALLBACK
    assert Path(table).is_file()
    assert Path(table).name == tables.ASSET


def test_spectrum_table_honours_prefer_hires(tables):
    assert tables.spectrum_table(prefer_hires=False) == tables.FALLBACK


def test_spectrum_table_falls_back_when_the_asset_is_missing(tables, tmp_path,
                                                             monkeypatch):
    """Page 1's exercise 3 and the browser-without-the-asset case."""
    monkeypatch.chdir(tmp_path)
    assert tables.spectrum_table(base_url="no-such-directory") == tables.FALLBACK


def test_spectrum_table_result_is_usable_as_a_table_argument(tables):
    """Whatever it returns must construct a component -- name or path alike."""
    table = tables.spectrum_table()
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
    assert lw.num_longwave_bands == (14 if table == tables.FALLBACK else 56)


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


@pytest.mark.parametrize("table, nband, required_edges", [
    pytest.param(LOW_RES, 14, (630.0, 700.0, 800.0, 980.0, 1080.0, 1180.0),
                 id="14band"),
    # The 56-band table resolves the CO2 core into four bands and the window
    # into eleven; pages 2 and 3 quote both, so their boundaries are pinned.
    pytest.param(HI_RES, 56,
                 (630.0, 647.5, 665.0, 682.5, 700.0, 800.0, 980.0, 1005.0,
                  1155.0, 2525.0, 2887.5),
                 id="56band",
                 marks=pytest.mark.skipif(
                     HI_RES == _TABLES.FALLBACK,
                     reason="the 56-band spectrum table asset is not present")),
])
def test_band_limits_cover_the_expected_features(spectra, table, nband,
                                                 required_edges):
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
    limits = spectra.band_limits_of(lw)
    assert limits.shape == (nband, 2)
    # the 15 um CO2 core and the window edges must be exact band boundaries
    edges = set(np.round(limits.ravel(), 1))
    for edge in required_edges:
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
@pytest.mark.parametrize("table", PAGE123_TABLES)
def test_page1_window_is_warmer_than_the_co2_core(soundings, spectra, table):
    """The page's central claim: brightness temperature swings ~200 K to ~285 K."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
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
@pytest.mark.parametrize("table", PAGE123_TABLES)
def test_page2_window_is_transparent_and_co2_core_is_opaque(soundings, spectra,
                                                            table):
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
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
@pytest.mark.parametrize("table", PAGE123_TABLES)
def test_page2_removing_absorbers_opens_the_window(soundings, table):
    """A near-N2 atmosphere: OLR should approach the surface blackbody flux."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
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


def _page3_column(soundings, spectra, table=LOW_RES, humidity_scale=1.0, nz=60):
    """Page 3's column: the standard sounding, with a humidity knob."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
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
@pytest.mark.parametrize("table, n_levels, min_span, max_core_pa", [
    pytest.param(LOW_RES, 12, 100.0, 2000.0, id="14band"),
    # The numbers page 3 quotes in prose: forty-six of the fifty-six bands have
    # a tau* = 1 level, and those levels run 1.3 hPa to 1005 hPa -- a factor of
    # 770. The core reaches 46 hPa because the 56-band table resolves the band's
    # outer wing at 630-648, which the 14-band table averages away.
    pytest.param(HI_RES, 46, 700.0, 5000.0, id="56band",
                 marks=pytest.mark.skipif(
                     HI_RES == _TABLES.FALLBACK,
                     reason="the 56-band spectrum table asset is not present")),
])
def test_page3_radiating_levels_are_a_distribution_not_a_height(
        soundings, spectra, table, n_levels, min_span, max_core_pa):
    """The page's central claim: there is no single radiating level."""
    col = _page3_column(soundings, spectra, table=table)
    limits, p_emit = col["limits"], col["p_emit"]
    centres = spectra.band_centres(limits)

    # The CO2 core emits from the stratosphere and the top of the troposphere.
    core = p_emit[(centres > 630) & (centres < 700)]
    assert np.all(core < max_core_pa)

    # The clearest window bands never reach tau* = 1 at all: no radiating level
    # exists for them, the surface is simply visible.
    window = p_emit[(centres > 980) & (centres < 1080)]
    assert np.all(np.isnan(window))

    # Across the bands that do have one, the levels span the whole atmosphere.
    have = p_emit[np.isfinite(p_emit)]
    assert have.size == n_levels
    assert have.max() / have.min() > min_span


@pytest.mark.slow
@pytest.mark.skipif(HI_RES == _TABLES.FALLBACK,
                    reason="the 56-band spectrum table asset is not present")
def test_page3_nine_of_the_eleven_window_bands_have_no_radiating_level(
        soundings, spectra):
    """Page 3 states the count outright, so the count is what gets asserted.

    The two exceptions are 800-845 and 1130-1155 cm^-1 at the window's edges,
    whose column tau only just scrapes past 1/D, so their "radiating level"
    lands at or below the bottom model layer and has degenerated into "the
    ground". Page 2 measured the same eleven bands as the window.
    """
    col = _page3_column(soundings, spectra, table=HI_RES)
    limits, p_emit = col["limits"], col["p_emit"]
    window = (limits[:, 0] >= 800.0) & (limits[:, 1] <= 1155.0)

    assert window.sum() == 11
    assert np.isnan(p_emit[window]).sum() == 9
    # The two that do have one are at the edges, and both sit near the ground.
    edge = limits[window][np.isfinite(p_emit[window])]
    np.testing.assert_allclose(edge, [[800.0, 845.0], [1130.0, 1155.0]])
    assert np.nanmin(p_emit[window]) > 9e4        # Pa: at or below the surface layer


@pytest.mark.slow
@pytest.mark.parametrize("table, n_below_weighted, n_below_level, max_bias, min_corr", [
    pytest.param(LOW_RES, 0, 1, 45.0, 0.9, id="14band"),
    # Page 3 names all six: three window-edge bands whose level degenerates into
    # the ground, and the three saturated core bands emitting from an isothermal
    # 200 K stratosphere, where "the temperature at tau* = 1" barely constrains.
    # The miss is larger here (the page quotes +13, +25 and +80 K) because the
    # narrower bands separate the strongly and weakly absorbing g-points that
    # the 14-band table averages together.
    pytest.param(HI_RES, 3, 6, 90.0, 0.85, id="56band",
                 marks=pytest.mark.skipif(
                     HI_RES == _TABLES.FALLBACK,
                     reason="the 56-band spectrum table asset is not present")),
])
def test_page3_brightness_temperature_tracks_the_emission_level(
        soundings, spectra, table, n_below_weighted, n_below_level, max_bias,
        min_corr):
    """T_b(band) follows the emission-weighted temperature -- with a cold bias.

    The band-mean optical depth cork reports is the g-weighted mean over eight
    g-points whose k values span three to eight orders of magnitude. It is
    therefore dominated by the band's *strongest* absorption, while the OLR
    escapes preferentially through the band's weakest. So the band-mean
    weighting function places the emission too high, and T_b comes out warmer
    than the weighted temperature. That one-sidedness is the page's punchline,
    so it is what gets asserted.

    It has exactly one class of exception, and only on the 56-band table: the
    saturated CO2 core bands emit from the isothermal 200 K stratosphere, where
    the weighted temperature is 200.0 K whatever the weights do, so the
    comparison has no leverage and the residual band-centre inversion error
    (about 2 K) sets the sign instead. Page 3 names the same three bands in its
    tau* = 1 discussion, at 0.3 to 2.2 K.
    """
    col = _page3_column(soundings, spectra, table=table)
    T, tau, limits = col["T"], col["tau"], col["limits"]
    _, flux_density = spectra.spectral_olr(col["diag"], limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))

    t_weighted = np.array([
        float((spectra.emission_weight(tau[:, b], PAGE3_D) * T).sum()
              / spectra.emission_weight(tau[:, b], PAGE3_D).sum())
        for b in range(tau.shape[1])])

    sane = tb <= 288.0            # see the band-centre test below
    bias = (tb - t_weighted)[sane]
    assert (bias <= 0.0).sum() == n_below_weighted
    assert bias.min() > -2.5      # and the exceptions miss by very little
    assert bias.max() < max_bias
    assert np.corrcoef(t_weighted[sane], tb[sane])[0, 1] > min_corr

    # The page plots T at the tau* = 1 level rather than the weighted mean, so
    # guard that comparison too. Same story, with a handful of bands going the
    # other way -- the ones whose level has degenerated into the ground, or that
    # sit in the isothermal stratosphere. Page 3 names each of them; the count is
    # the parameter, and the page quotes the size of the miss as 0.3 to 6.6 K.
    p_emit = col["p_emit"]
    has_level = np.isfinite(p_emit) & sane
    t_level = np.exp(np.interp(np.log(p_emit[has_level]),
                               np.log(col["p"][::-1]), np.log(T[::-1])))
    level_bias = tb[has_level] - t_level
    assert (level_bias <= 0.0).sum() == n_below_level
    assert level_bias.min() > -7.0
    assert level_bias.max() < max_bias
    assert np.corrcoef(t_level, tb[has_level])[0, 1] > min_corr


@pytest.mark.slow
@pytest.mark.parametrize("table, band, tb_approx", [
    pytest.param(LOW_RES, (1800.0, 3250.0), 300.0, id="14band"),
    # The 56-band table splits that 1450 cm^-1 monster into three; the artifact
    # survives in the one that carries the flux, and shrinks from 300 K to 290 K
    # because the band is narrower. Page 1 shows both side by side; page 3's
    # callout is written about this one.
    pytest.param(HI_RES, (2525.0, 2887.5), 290.0, id="56band",
                 marks=pytest.mark.skipif(
                     HI_RES == _TABLES.FALLBACK,
                     reason="the 56-band spectrum table asset is not present")),
])
def test_page3_band_centre_brightness_temperature_fails_on_one_wide_band(
        soundings, spectra, table, band, tb_approx):
    """One band returns T_b above the surface -- an inversion artifact.

    Inverting a band-mean flux density at the band *centre* only works if the
    Planck function is near-linear across the band. Over hundreds of cm^-1 of
    shortwave tail it is not, and the answer comes back hotter than the ground,
    which is impossible. The pages say so rather than plotting it silently.
    """
    col = _page3_column(soundings, spectra, table=table)
    limits = col["limits"]
    _, flux_density = spectra.spectral_olr(col["diag"], limits)
    tb = spectra.brightness_temperature(flux_density, spectra.band_centres(limits))

    bogus = np.flatnonzero(tb > 288.0)
    assert bogus.size == 1                    # exactly one, and we know which
    np.testing.assert_allclose(limits[bogus[0]], band)
    assert tb[bogus[0]] == pytest.approx(tb_approx, abs=1.0)


@pytest.mark.slow
@pytest.mark.parametrize("table", PAGE123_TABLES)
def test_page3_moistening_raises_the_emission_levels(soundings, spectra, table):
    """The knob: a wetter atmosphere emits from higher, colder air."""
    dry = _page3_column(soundings, spectra, table=table, humidity_scale=1.0)
    wet = _page3_column(soundings, spectra, table=table, humidity_scale=4.0)
    centres = spectra.band_centres(dry["limits"])

    # Every band that had a radiating level keeps it, and moves up.
    had = np.isfinite(dry["p_emit"])
    assert np.all(np.isfinite(wet["p_emit"][had]))
    assert np.all(wet["p_emit"][had] <= dry["p_emit"][had] + 1.0)

    # The water vapour bands move a long way; the CO2 core barely moves. The
    # strongest mover is what gets asserted, not every band: on the 56-band
    # table the 70-130 cm^-1 band is already emitting from the model lid when
    # dry, so it has nowhere left to rise.
    rotational = (centres > 10) & (centres < 250)
    assert (dry["p_emit"][rotational] / wet["p_emit"][rotational]).max() > 2.0
    core = (centres > 630) & (centres < 700)
    assert (dry["p_emit"][core] / wet["p_emit"][core]).max() < 1.2

    # And the window closes: bands with no radiating level acquire one.
    assert np.any(np.isfinite(wet["p_emit"][~had]))

    olr_dry = float(dry["diag"]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    olr_wet = float(wet["diag"]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert olr_wet < olr_dry - 20.0


TAU_INF, D_NOTES, T_E = 4.0, 2.0, 255.0


def _gray_equilibrium_state(soundings, nz=60):
    lw = CorkLongwaveRadiation(optics="correlated_k", table="tour_gray_lw",
                               diffusivity_factor=D_NOTES)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=nz))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, T_surf, _ = soundings.analytic_gray_equilibrium(p, ps, TAU_INF, T_E)
    soundings.apply_sounding(state, T, T_surf=T_surf)
    return lw, state, T, T_surf


@pytest.mark.slow
def test_page4_analytic_profile_is_an_equilibrium(soundings):
    """The notes' closed form must give ~zero heating and OLR = sigma T_e^4."""
    lw, state, _, _ = _gray_equilibrium_state(soundings)
    tendencies, diag = lw(state)
    H = tendencies["air_temperature"].values[:, 0, 0] * 86400.0    # K/day
    olr = float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    assert np.abs(H).max() < 0.05
    assert olr == pytest.approx(5.670374419e-8 * T_E ** 4, rel=2e-3)


@pytest.mark.slow
def test_page4_same_profile_is_not_an_equilibrium_for_a_real_spectrum(soundings):
    """The whole motivation for non-grey radiation, in one comparison."""
    lw_gray, state_gray, T, T_surf = _gray_equilibrium_state(soundings)
    H_gray = lw_gray(state_gray)[0]["air_temperature"].values[:, 0, 0] * 86400.0

    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=60))
    soundings.apply_sounding(state, T, np.full(T.shape, 1e-6), T_surf=T_surf)
    H_real = lw(state)[0]["air_temperature"].values[:, 0, 0] * 86400.0

    rms = lambda x: float(np.sqrt((x ** 2).mean()))
    assert rms(H_real) > 100 * rms(H_gray)
    assert np.abs(H_real).max() > 5.0


def _page5_column(soundings, nz=40):
    """Page 5's column: the standard sounding, ready for a CO2 sweep."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=nz))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    T, q = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)
    return lw, state, p, T, q


def _olr_at_co2(lw, state, soundings, T, q, T_surf, co2_ppm):
    soundings.apply_sounding(state, T, q, T_surf=T_surf)
    state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2_ppm * 1e-6
    _, diag = lw(state)
    return float(diag["upwelling_longwave_flux_in_air"].values[-1, 0, 0])


@pytest.mark.slow
def test_page5_co2_doubling_forcing_is_canonical(soundings):
    """Students measure ~3.7 W/m2 per doubling themselves."""
    lw, state, _, T, q = _page5_column(soundings)

    olr = {c: _olr_at_co2(lw, state, soundings, T, q, 288.0, c)
           for c in (280, 560, 1120)}
    first, second = olr[280] - olr[560], olr[560] - olr[1120]

    assert 3.0 < first < 4.5           # canonical ~3.7 W/m2
    assert 3.0 < second < 4.5
    assert abs(first - second) < 1.0   # logarithmic, not linear


@pytest.mark.slow
def test_page5_core_optical_depth_is_linear_in_co2(soundings, spectra):
    """The other half of the page's beat: tau_core is *proportional* to CO2.

    Forcing is logarithmic in concentration while the absorber it comes from is
    linear in it. The page reconciles the two with page 3's weighting function,
    so both halves are guarded: this one fits log(tau_core) against log(ppm) and
    demands an exponent of 1. A dry column isolates CO2 from water vapour, which
    overlaps the wings of the 15 um band.
    """
    lw, state, _, T, _ = _page5_column(soundings)
    dry = np.full(T.shape, 1e-7)
    centres = spectra.band_centres(spectra.band_limits_of(lw))
    core = (centres > 630) & (centres < 700)

    co2 = np.array([30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0])
    tau_core = []
    for ppm in co2:
        soundings.apply_sounding(state, T, dry, T_surf=288.0)
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = ppm * 1e-6
        _, diag = lw(state)
        tau = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]
        tau_core.append(float(tau.sum(axis=0)[core].sum()))

    exponent, _ = np.polyfit(np.log(co2), np.log(np.array(tau_core)), 1)
    assert exponent == pytest.approx(1.0, abs=0.03)


@pytest.mark.slow
def test_page5_wing_brightness_temperature_falls_by_a_fixed_step(soundings, spectra):
    """Why the forcing is logarithmic, in the one quantity that is unbiased.

    Page 3 showed that a band-mean tau* = 1 level is biased high, so the page
    cannot argue from where the emission level *is*. Brightness temperature has
    no such bias: it is inverted straight from the band's escaping flux. In the
    700-800 cm^-1 wing it drops by a near-constant step for every doubling --
    a fixed temperature step per FACTOR of two in concentration is exactly what
    "logarithmic" means. The saturated core does the opposite and flattens out.
    """
    lw, state, _, T, q = _page5_column(soundings)
    limits = spectra.band_limits_of(lw)
    centres = spectra.band_centres(limits)
    wing = int(np.flatnonzero((centres > 700) & (centres < 800))[0])
    core = int(np.flatnonzero((centres > 630) & (centres < 700))[0])

    ladder = [140, 280, 560, 1120, 2240, 4480]
    tb = {}
    for ppm in ladder:
        soundings.apply_sounding(state, T, q, T_surf=288.0)
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = ppm * 1e-6
        _, fd = spectra.spectral_olr(lw(state)[1], limits)
        tb[ppm] = spectra.brightness_temperature(fd, centres)

    steps = np.array([tb[b][wing] - tb[a][wing] for a, b in zip(ladder, ladder[1:])])
    assert np.all(steps < 0)                       # every doubling cools it
    assert steps.std() < 0.25                      # by very nearly the same step
    assert -5.5 < steps.mean() < -4.0

    core_steps = np.array([tb[b][core] - tb[a][core]
                           for a, b in zip(ladder, ladder[1:])])
    # The core saturates instead: each doubling buys less than the last.
    assert np.all(np.diff(np.abs(core_steps)) < 0)
    assert abs(core_steps[-1]) < 0.2 * abs(core_steps[0])


@pytest.mark.slow
def test_page5_forcing_comes_from_the_wings_not_the_core(soundings, spectra):
    lw, state, _, T, q = _page5_column(soundings)
    limits = spectra.band_limits_of(lw)
    centres = spectra.band_centres(limits)

    per_band = {}
    for co2 in (280, 1120):
        soundings.apply_sounding(state, T, q, T_surf=288.0)
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = co2 * 1e-6
        per_band[co2], _ = spectra.spectral_olr(lw(state)[1], limits)

    delta = per_band[280] - per_band[1120]
    core = delta[(centres > 630) & (centres < 700)].sum()
    wings = delta[((centres > 500) & (centres < 630))
                  | ((centres > 700) & (centres < 800))].sum()
    assert wings > core          # the saturated core cannot contribute much


@pytest.mark.slow
def test_page6_fixed_rh_suppresses_olr_relative_to_fixed_q(soundings):
    """The water vapour feedback: OLR rises more slowly when moisture responds."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])
    _, q_ref = soundings.lapse_rate_sounding(p, ps, T_surf=288.0, rh=0.8)

    def olr(T_surf, fixed_q):
        T, q = soundings.lapse_rate_sounding(p, ps, T_surf=T_surf, rh=0.8)
        soundings.apply_sounding(state, T, q_ref if fixed_q else q, T_surf=T_surf)
        return float(lw(state)[1]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    d_fixed_q = olr(298.0, True) - olr(288.0, True)
    d_fixed_rh = olr(298.0, False) - olr(288.0, False)
    assert d_fixed_q > d_fixed_rh > 0        # moisture damps the OLR response


@pytest.mark.slow
def test_page6_olr_saturates_at_high_surface_temperature(soundings):
    """Approach to the Simpson-Nakajima limit on saturated soundings."""
    lw = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=40))
    p = state["air_pressure"].values[:, 0, 0]
    ps = float(state["surface_air_pressure"].values.ravel()[0])

    def olr(T_surf):
        T, q = soundings.lapse_rate_sounding(p, ps, T_surf=T_surf, rh=1.0)
        soundings.apply_sounding(state, T, q, T_surf=T_surf)
        return float(lw(state)[1]["upwelling_longwave_flux_in_air"].values[-1, 0, 0])

    warm = olr(310.0) - olr(300.0)
    hot = olr(340.0) - olr(330.0)
    assert hot < warm            # the OLR response flattens
    assert olr(340.0) < 400.0    # nowhere near sigma T^4 = 757 W/m2


@pytest.fixture
def assets():
    return _load("assets")


def test_assets_resolve_finds_a_staged_file_from_any_working_directory(
        assets, monkeypatch, tmp_path):
    """Resolution must not depend on the caller's working directory.

    In the browser the page's working directory *is* the page directory, so
    `_data/x.npz` resolves directly. Natively -- a test run, a static-figure
    render, someone poking at it from the repo root -- it does not, and the
    fallback is the path relative to this module's own location.
    """
    monkeypatch.chdir(tmp_path)
    found = assets.resolve("earth_spectrum_lw.npz")
    assert found is not None and os.path.isfile(found)


def test_assets_resolve_returns_none_for_a_missing_asset(assets, tmp_path,
                                                        monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert assets.resolve("no_such_asset.npz") is None


def test_assets_locate_falls_back_to_resolve_outside_pyodide(assets, tmp_path,
                                                             monkeypatch):
    """`locate` is `resolve` plus a browser-only fetch. Natively they agree."""
    monkeypatch.chdir(tmp_path)
    assert assets.locate("earth_spectrum_lw.npz") == \
        assets.resolve("earth_spectrum_lw.npz")
    assert assets.locate("no_such_asset.npz") is None


def test_assets_honours_a_non_default_base_url(assets, tmp_path, monkeypatch):
    """Page 11 and 12 pass the same `_data`, but the argument must work."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "elsewhere").mkdir()
    (tmp_path / "elsewhere" / "probe.npz").write_bytes(b"x")
    assert assets.resolve("probe.npz", base_url="elsewhere") == \
        os.path.join("elsewhere", "probe.npz")


# ------------------------------------------------------- _tour/stepping.py


@pytest.fixture
def stepping():
    return _load("stepping")


def _gray_column(nz=28, slab_depth=2.0, table="single_band_gray_lw"):
    """Page 07's column: gray longwave over a slab, prescribed absorbed solar."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k", table=table)
    surface = climt.SlabSurface()
    state = climt.get_default_state([longwave, surface],
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, surface], state


def test_integrate_advances_the_clock_and_the_state(stepping):
    components, state = _gray_column()
    start = state["time"]
    T0 = state["air_temperature"].values.copy()

    timestep = climt.UnytTimeDelta(hours=12)
    out = stepping.integrate(components, [], state, timestep, 10)

    assert out is state                       # stepped in place, and returned
    assert (out["time"] - start).total_seconds() == pytest.approx(10 * 12 * 3600)
    assert not np.allclose(out["air_temperature"].values, T0)


def test_integrate_keeps_the_flux_diagnostics(stepping):
    """The update-order gotcha that broke the old demo once already.

    The prognostic state is applied *before* the diagnostics. Reverse it and
    the stepper carries the longwave fluxes forward at their pre-step value,
    SlabSurface consumes stale fluxes, and the surface heats without bound.
    """
    components, state = _gray_column()
    stepping.integrate(components, [], state,
                       climt.UnytTimeDelta(hours=12), 100)

    for flux in ("upwelling_longwave_flux_in_air",
                 "downwelling_longwave_flux_in_air"):
        values = state[flux].values
        assert np.all(np.isfinite(values)) and np.all(values >= 0.0)
    surface = float(state["surface_temperature"].values.ravel()[0])
    assert 150.0 < surface < 400.0


def test_integrate_runs_stepper_components_too(stepping):
    """Page 08 onward puts a Stepper in the same loop as the tendencies."""
    components, state = _gray_column()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    state = climt.get_default_state(components + [boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = 2.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    stepping.integrate(components, [boundary_layer], state,
                       climt.UnytTimeDelta(hours=1), 10)

    # The stepper's own diagnostics survive into the state.
    assert "surface_upward_sensible_heat_flux" in state
    assert np.all(np.isfinite(
        state["surface_upward_sensible_heat_flux"].values))
    assert np.all(np.isfinite(state["air_temperature"].values))


def test_integrate_rejects_a_plain_timedelta(stepping):
    """The UnytTimeDelta gotcha, asserted rather than left as prose.

    Page 07 teaches this; if sympl ever stops raising, the page's craft note is
    wrong and this test is where we find out.
    """
    from datetime import timedelta

    components, state = _gray_column()
    with pytest.raises(Exception) as caught:
        stepping.integrate(components, [], state, timedelta(hours=12), 1)
    assert "degK" in str(caught.value) or "unit" in str(caught.value).lower()


def test_history_records_every_step_and_n_snapshots(stepping):
    components, state = _gray_column()
    state, history = stepping.integrate_with_history(
        components, [], state, climt.UnytTimeDelta(hours=12), 12,
        n_snapshots=4)

    assert history["days"].shape == (12,)
    assert history["olr"].shape == (12,)
    assert history["surface_temperature"].shape == (12,)
    assert history["days"][0] == pytest.approx(0.5)     # 12 h, in days
    assert history["days"][-1] == pytest.approx(6.0)
    assert len(history["snapshots"]) == 4
    first, last = history["snapshots"][0], history["snapshots"][-1]
    assert first["day"] < last["day"]
    for key in ("T", "H", "U", "D"):
        assert np.all(np.isfinite(first[key]))
    assert first["T"].shape == (28,)          # mid levels
    assert first["U"].shape == (29,)          # interface levels


def test_history_olr_matches_the_final_state(stepping):
    """The recorded series is the state's own diagnostic, not a recomputation."""
    components, state = _gray_column()
    state, history = stepping.integrate_with_history(
        components, [], state, climt.UnytTimeDelta(hours=12), 8)

    assert history["olr"][-1] == pytest.approx(
        float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0]))
    assert history["surface_temperature"][-1] == pytest.approx(
        float(state["surface_temperature"].values.ravel()[0]))


def test_wind_relaxation_holds_the_wind_up(stepping):
    """Without a momentum source the column spins down; with one it does not.

    Measured without relaxation: the lowest-level wind reaches 0.000 m/s and
    the equilibrium is independent of how fast it started.
    """
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table=PAGE7_GRAY_TABLE,
                                           diffusivity_factor=PAGE7_DIFFUSIVITY)
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                               roughness_length=1e-3)
    state = climt.get_default_state([longwave, surface, boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = 2.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    relaxation = stepping.wind_relaxation(state, 8.0, timescale_hours=24.0)
    assert "equilibrium_eastward_wind" in state
    assert state["eastward_wind_relaxation_timescale"].attrs["units"] == "s"

    stepping.integrate([longwave, surface, relaxation], [boundary_layer],
                       state, climt.UnytTimeDelta(hours=1), 500)

    lowest = float(state["eastward_wind"].values[0, 0, 0])
    assert lowest > 1.0, (
        f"lowest-level wind {lowest:.3f} m/s — the relaxation is not holding "
        "the wind up against the surface drag")
    assert lowest < 8.0, (
        "the lowest level should sit below the target: drag is still acting, "
        "which is the point")


def test_wind_relaxation_leaves_a_loaded_wind_alone(stepping):
    """``initialise=False`` is how pages 11/12 attach to a state from disk.

    That state carries a wind profile already sheared by the surface drag;
    overwriting it with a uniform ``speed`` would throw away part of the
    equilibrium, and nothing would error.
    """
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    state = climt.get_default_state([boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    profile = np.linspace(2.0, 9.0, 28).reshape(28, 1, 1)
    state["eastward_wind"].values[:] = profile

    stepping.wind_relaxation(state, 5.0, initialise=False)

    np.testing.assert_allclose(state["eastward_wind"].values, profile)
    assert np.all(state["equilibrium_eastward_wind"].values == 5.0)


def test_unyt_relaxation_units_are_parseable(stepping):
    """Pins the reason the subclass exists.

    The stock component's tendency units come from pint and read
    '1.0 meter / second ** 2', which unyt refuses. If sympl ever emits
    something unyt can parse, this shim can go -- and this is where we notice.
    """
    import sympl as _sympl

    stock = _sympl.RelaxationTendencyComponent("eastward_wind", "m s^-1")
    assert "**" in stock.tendency_properties["eastward_wind"]["units"], (
        "sympl's tendency units are no longer pint-formatted — re-check "
        "whether stepping.UnytRelaxation is still needed")

    ours = stepping.UnytRelaxation("eastward_wind", "m s^-1", "m s^-2")
    assert ours.tendency_properties["eastward_wind"]["units"] == "m s^-2"


def test_draw_evolution_builds_the_four_panel_figure(stepping):
    """The figure the pages end on: it must draw without a display attached."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    components, state = _gray_column()
    state, history = stepping.integrate_with_history(
        components, [], state, climt.UnytTimeDelta(hours=12), 8, n_snapshots=3)

    plt.close("all")
    try:
        with pytest.warns(UserWarning, match="non-interactive"):
            stepping.draw_evolution(history, state, SOLAR, title="test")
        figure = plt.gcf()
        # three profile panels, the OLR panel with its twin, and the colorbar
        assert len(figure.axes) == 6
        assert figure.axes[0].get_yscale() == "log"
    finally:
        plt.close("all")


# ------------------------------------------------------------------ page 7
#
# Inherited from the deleted tests/test_live_rce_demo.py, on the deleted
# radiative-transfer/09-live-rce.qmd. Same two claims, re-measured on page 7's
# configuration: UnytBackend rather than DataArrayBackend, nz=28 rather than 18,
# dt=12 h rather than 2 h, a 2 m slab rather than 5 m.
#
# READ THIS BEFORE WRITING PAGE 07.
#
# 1. These tests hard-code page 07's five configuration choices, and the page
#    must declare exactly those: `tour_gray_lw` at `diffusivity_factor=2.0`,
#    `earth_low_res_lw` at the default 1.66, nz=28, dt=12 h, and a 2 m slab.
#
# 2. THE OLR SEPARATION IS A TRANSIENT OF THE 300-STEP RUN, NOT AN EQUILIBRIUM
#    RESULT. At PAGE7_STEPS the gray column is essentially converged (OLR
#    234.4 against SOLAR = 240) while the non-grey one is not (OLR 265.8, i.e.
#    26 W/m^2 above the forcing and still cooling hard). At the true
#    equilibrium page 07 is heading for, *both* columns radiate 240 W/m^2 and
#    the OLR gap closes to zero by construction. So page 07 must either run at
#    PAGE7_STEPS, or add a second test at its own step count, or simply not
#    quote the OLR gap in its text. The surface-temperature and top-gradient
#    separations do not have this problem: both strengthen as the non-grey
#    column converges.
#
# 3. Page 07 must NOT reuse the deleted page's "change one string and re-run"
#    framing. Its two columns differ in table *and* in diffusivity factor
#    *and* in total optical depth -- it is no longer one string.

PAGE7_STEPS = 300      # far enough to separate the two columns unambiguously;
                       # the page itself runs longer. See the log below.

# The tables page 7 actually declares, and the diffusivity that goes with the
# gray one. NOT `single_band_gray_lw`, which is what the deleted page used:
# page 7's gray half must reproduce page 4's analytic profile, and that only
# works against the table page 4 calibrated (D = 2, D * sum(tau) = 4).
# Per the global constraint, a test asserts on the table its page declares.
PAGE7_GRAY_TABLE = "tour_gray_lw"
PAGE7_DIFFUSIVITY = 2.0
PAGE7_NONGREY_TABLE = "earth_low_res_lw"


def _page7_gray_column(nz=28, slab_depth=2.0):
    """Page 7's gray column: page 4's table and diffusivity, over a slab."""
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table=PAGE7_GRAY_TABLE,
        diffusivity_factor=PAGE7_DIFFUSIVITY)
    surface = climt.SlabSurface()
    state = climt.get_default_state([longwave, surface],
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, surface], state


@pytest.fixture(scope="module")
def page7_columns():
    """Both of page 7's integrations, run once and shared."""
    # Load-bearing: this module-scoped fixture is instantiated before the
    # function-scoped autouse `_unyt_backend` one.
    sympl.set_backend(climt.UnytBackend())
    stepping_module = _load("stepping")
    timestep = climt.UnytTimeDelta(hours=12)

    gray_components, gray_state = _page7_gray_column()
    # `_gray_column` is the shared stepping helper, named for its default
    # table; handed the non-grey table it builds the non-grey column.
    nongrey_components, nongrey_state = _gray_column(
        table=PAGE7_NONGREY_TABLE)
    return {
        PAGE7_GRAY_TABLE: stepping_module.integrate(
            gray_components, [], gray_state, timestep, PAGE7_STEPS),
        PAGE7_NONGREY_TABLE: stepping_module.integrate(
            nongrey_components, [], nongrey_state, timestep, PAGE7_STEPS),
    }


def _olr(state):
    return float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0])


def _surface_temperature(state):
    return float(state["surface_temperature"].values.ravel()[0])


@pytest.mark.slow
def test_page7_nongrey_column_is_the_more_efficient_radiator(page7_columns):
    """Page 7's comparison: its two configurations separate, and which way.

    These are page 7's two whole configurations, not one knob: they differ in
    table, in diffusivity factor (2.0 against the default 1.66) and in total
    optical depth (the gray table is calibrated to tau_inf = 4). The
    separation asserted here is their combined effect -- the gray column's
    opacity dominates the surface gap; spectral windows are one contributor,
    not the whole of it.

    Measured separations at PAGE7_STEPS (see the log): OLR 265.8 vs 234.4
    W/m^2 and surface 271.9 vs 334.1 K. Thresholds are half of each. NOTE the
    OLR half of this is transient -- see point 2 of the page-7 comment block.
    """
    gray = page7_columns[PAGE7_GRAY_TABLE]
    nongrey = page7_columns[PAGE7_NONGREY_TABLE]

    assert _olr(nongrey) > _olr(gray) + 15.0, (
        f"non-grey OLR {_olr(nongrey):.1f} should exceed gray "
        f"{_olr(gray):.1f} W/m^2 — windows radiate to space more efficiently")
    assert _surface_temperature(nongrey) < _surface_temperature(gray) - 30.0, (
        f"non-grey surface {_surface_temperature(nongrey):.1f} K should be "
        f"cooler than gray {_surface_temperature(gray):.1f} K")


@pytest.mark.slow
def test_page7_nongrey_column_cools_faster_aloft(page7_columns):
    """Page 7's second comparison: a steeper temperature drop-off aloft.

    Band-resolved absorption concentrates cooling in the strongly absorbing
    bands high in the column, which the gray column smears out. Measured:
    -6.84 K/level non-grey against -2.82 K/level gray; threshold is half.
    """
    top_gradient = {}
    for name, state in page7_columns.items():
        temperature = state["air_temperature"].values[:, 0, 0]
        # Level index increases upward; mean of the top three level-to-level
        # differences, in K per level.
        top_gradient[name] = float(np.mean(np.diff(temperature)[-3:]))

    assert top_gradient[PAGE7_NONGREY_TABLE] < \
        top_gradient[PAGE7_GRAY_TABLE] - 2.0, (
        f"non-grey top-of-column gradient {top_gradient[PAGE7_NONGREY_TABLE]:.2f}"
        f" K/level should be markedly steeper than gray "
        f"{top_gradient[PAGE7_GRAY_TABLE]:.2f} K/level")


@pytest.mark.slow
@pytest.mark.parametrize("table", [PAGE7_GRAY_TABLE, PAGE7_NONGREY_TABLE])
def test_page7_flux_diagnostics_survive_the_update_order(page7_columns, table):
    """The update-order gotcha, on both of page 7's tables.

    ``stepping.integrate`` applies the prognostic state *before* the
    diagnostics. Reversing that clobbers the longwave fluxes SlabSurface
    consumes, and the surface heats without bound.
    """
    state = page7_columns[table]

    for flux in ("upwelling_longwave_flux_in_air",
                 "downwelling_longwave_flux_in_air"):
        assert flux in state, f"{flux} missing — diagnostics were clobbered"
        values = state[flux].values
        assert np.all(np.isfinite(values)), f"{flux} has non-finite values"
        assert np.all(values >= 0.0), f"{flux} must be non-negative"

    temperature = state["air_temperature"].values
    assert np.all(np.isfinite(temperature))
    assert temperature.min() > 150.0 and temperature.max() < 400.0, (
        f"air temperature ran to [{temperature.min():.1f}, "
        f"{temperature.max():.1f}] K — suspect a flux-coupling regression")
    assert 150.0 < _surface_temperature(state) < 400.0


@pytest.fixture
def budgets():
    return _load("budgets")


def test_toa_imbalance_is_absorbed_solar_minus_olr(budgets):
    """With no atmospheric shortwave absorption, TOA net = surface SW - OLR.

    Every page in this tranche prescribes the absorbed shortwave at the
    surface and runs no shortwave component, so the planet's absorbed
    shortwave *is* the surface's. The function must read that from the state
    rather than take it as an argument, or a page can quote a budget that
    disagrees with its own configuration.
    """
    components, state = _gray_column()
    stepping_module = _load("stepping")
    stepping_module.integrate(components, [], state,
                              climt.UnytTimeDelta(hours=12), 5)

    olr = float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert budgets.toa_imbalance(state) == pytest.approx(SOLAR - olr, abs=1e-9)


def test_toa_imbalance_closes_on_a_converged_column(budgets):
    """The whole point: near equilibrium, the budget is near zero."""
    components, state = _gray_column()
    _load("stepping").integrate(components, [], state,
                                climt.UnytTimeDelta(hours=12), 700)
    assert abs(budgets.toa_imbalance(state)) < 1.0


def test_surface_imbalance_counts_every_term(budgets):
    """Radiation in, radiation out, and the two turbulent fluxes out."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="single_band_gray_lw")
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    state = climt.get_default_state([longwave, surface, boundary_layer],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["ocean_mixed_layer_thickness"].values[:] = 2.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    _load("stepping").integrate([longwave, surface], [boundary_layer], state,
                                climt.UnytTimeDelta(hours=1), 10)

    expected = (
        SOLAR
        + float(state["downwelling_longwave_flux_in_air"].values[0, 0, 0])
        - float(state["upwelling_longwave_flux_in_air"].values[0, 0, 0])
        - float(state["surface_upward_sensible_heat_flux"].values.ravel()[0])
        - float(state["surface_upward_latent_heat_flux"].values.ravel()[0])
    )
    assert budgets.surface_imbalance(state) == pytest.approx(expected, abs=1e-9)


def test_surface_imbalance_works_without_turbulent_fluxes(budgets):
    """Page 7 has no boundary layer; the two flux terms are then zero."""
    components, state = _gray_column()
    _load("stepping").integrate(components, [], state,
                                climt.UnytTimeDelta(hours=12), 5)
    expected = (
        SOLAR
        + float(state["downwelling_longwave_flux_in_air"].values[0, 0, 0])
        - float(state["upwelling_longwave_flux_in_air"].values[0, 0, 0])
    )
    assert budgets.surface_imbalance(state) == pytest.approx(expected, abs=1e-9)


def test_column_enthalpy_matches_the_conservation_suite_formula(budgets):
    """Same integral tests/test_conservation.py uses, so page 10 can quote it.

    Cp_dry * T + Lv * q, mass-weighted by dp/g. If these two ever disagree,
    page 10's enthalpy-conservation claim is guarded by a test measuring
    something else.
    """
    from sympl import get_constant

    components, state = _gray_column()
    # A gray-longwave column carries no humidity of its own; give it one, so
    # the latent half of the integral is actually exercised.
    state["specific_humidity"] = state["air_temperature"].copy(deep=True)
    state["specific_humidity"].values[:] = 4e-3

    Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                       "J/kg/degK")
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    g = get_constant("gravitational_acceleration", "m/s^2")
    p_int = state["air_pressure_on_interface_levels"].values[:, 0, 0]
    dp = p_int[:-1] - p_int[1:]
    T = state["air_temperature"].values[:, 0, 0]
    q = state["specific_humidity"].values[:, 0, 0]
    expected = float(np.sum((Cpd * T + Lv * q) * dp / g))

    assert budgets.column_enthalpy(state) == pytest.approx(expected, rel=1e-12)


def test_pressure_thickness_is_positive_and_sums_to_the_column(budgets):
    """Bottom-first layer masses, in Pa, summing to the full column depth."""
    _components, state = _gray_column()
    p_int = state["air_pressure_on_interface_levels"].values[:, 0, 0]
    dp = budgets.pressure_thickness(state)

    assert dp.shape == (28,)
    assert np.all(dp > 0.0)
    assert float(np.sum(dp)) == pytest.approx(
        float(p_int[0] - p_int[-1]), rel=1e-12)


def test_evaporation_rate_inverts_the_latent_heat_flux(budgets):
    """Page 12 closes the moisture budget against this."""
    from sympl import get_constant

    components, state = _gray_column()
    state["surface_upward_latent_heat_flux"].values[:] = 80.0
    Lv = get_constant("latent_heat_of_condensation", "J/kg")
    assert budgets.evaporation_rate(state) == pytest.approx(
        80.0 / float(Lv) * 86400.0, rel=1e-12)


def test_evaporation_rate_is_zero_without_a_boundary_layer(budgets):
    """A page with no surface flux component has no latent heat flux at all.

    Missing quantities read as zero rather than raising, so page 10 can print
    the same budget table as page 12 without branching on its component list.
    """
    _components, state = _gray_column()
    del state["surface_upward_latent_heat_flux"]
    assert budgets.evaporation_rate(state) == 0.0


def test_precipitation_rate_adds_convective_and_grid_scale(budgets):
    """mm/day from Emanuel plus kg/m^2-per-step from GridScaleCondensation."""
    _components, state = _gray_column()
    assert budgets.precipitation_rate(
        state, climt.UnytTimeDelta(hours=1)) == 0.0

    state["convective_precipitation_rate"] = \
        state["surface_temperature"].copy(deep=True)
    state["convective_precipitation_rate"].values[:] = 3.0
    state["precipitation_amount"] = state["surface_temperature"].copy(deep=True)
    state["precipitation_amount"].values[:] = 0.5    # kg/m^2 per 1-hour step

    # 0.5 mm per hour is 12 mm/day, on top of the convective 3 mm/day.
    assert budgets.precipitation_rate(
        state, climt.UnytTimeDelta(hours=1)) == pytest.approx(15.0, rel=1e-12)
    assert budgets.precipitation_rate(
        state, climt.UnytTimeDelta(hours=2)) == pytest.approx(9.0, rel=1e-12)


def test_summary_reports_the_five_numbers_the_pages_quote(budgets):
    components, state = _gray_column()
    _load("stepping").integrate(components, [], state,
                                climt.UnytTimeDelta(hours=12), 5)
    report = budgets.summary(state)
    assert set(report) == {"surface_temperature", "olr", "absorbed_shortwave",
                           "toa_imbalance", "surface_imbalance"}
    assert all(np.isfinite(value) for value in report.values())
    assert report["absorbed_shortwave"] == pytest.approx(SOLAR)


def test_precipitation_rate_from_a_real_grid_scale_condensation(budgets):
    """The loop the hand-set-key test leaves open: run the real component.

    ``GridScaleCondensation`` declares ``precipitation_amount`` in kg m^-2,
    and ``precipitation_rate`` trusts that label. This steps the real
    component on a supersaturated column and checks the mm/day it produces
    against the condensed water mass computed by hand.
    """
    condensation = climt.GridScaleCondensation()
    state = climt.get_default_state([condensation],
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["air_temperature"].values[:] = 280.0
    state["specific_humidity"].values[:] = 0.0
    state["specific_humidity"].values[:5] = 0.03

    q_before = state["specific_humidity"].values.copy()
    p_int = state["air_pressure_on_interface_levels"].values.copy()

    timestep = climt.UnytTimeDelta(hours=1)
    diagnostics, outputs = condensation(state, timestep)
    state.update(diagnostics)

    # Hand-computed: condensed specific humidity times layer mass dp/g.
    dp = np.asarray(p_int[:-1, ...] - p_int[1:, ...])
    condensed = np.asarray(q_before - outputs["specific_humidity"].values)
    expected_kg_per_m2 = float(np.sum(condensed * dp / 9.80665))
    expected_mm_per_day = expected_kg_per_m2 * 86400.0 / float(
        timestep.total_seconds())

    rate = budgets.precipitation_rate(state, timestep)

    assert expected_kg_per_m2 > 0
    assert rate == pytest.approx(expected_mm_per_day, rel=1e-10)
    # An hour of condensing a supersaturated boundary layer is a heavy but
    # physical rain rate; a factor-of-1000 units error cannot land in here.
    assert 1e2 < rate < 1e5


@pytest.fixture
def states():
    return _load("states")


def _page11_components():
    """Page 11's exact component list, used by the round-trip tests."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    surface = climt.SlabSurface()
    boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
    adjustment = climt.DryConvectiveAdjustment()
    return [longwave, surface, boundary_layer, adjustment]


def test_state_round_trips_values_dims_and_units(states, tmp_path):
    """Save a state, load it back, and get the same state.

    Values, dimension names and units all have to survive. Values alone is not
    enough: a page that reloads a state with the right numbers under the wrong
    dims produces a plausible, wrong figure.
    """
    components = _page11_components()
    state = climt.get_default_state(components,
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    state["air_temperature"].values[:] = np.linspace(288.0, 200.0, 28
                                                     ).reshape(28, 1, 1)
    state["surface_temperature"].values[:] = 291.5
    state["specific_humidity"].values[:] = 3e-3

    path = states.save(str(tmp_path / "probe.npz"), state,
                       dict(note="round-trip probe"))
    reloaded, provenance = states.load(
        path, _page11_components(), grid_state=get_grid(nx=1, ny=1, nz=28))

    assert provenance["note"] == "round-trip probe"
    for name in ("air_temperature", "surface_temperature", "specific_humidity",
                 "air_pressure", "ocean_mixed_layer_thickness"):
        np.testing.assert_allclose(reloaded[name].values, state[name].values)
        assert tuple(reloaded[name].dims) == tuple(state[name].dims)
        assert reloaded[name].attrs["units"] == state[name].attrs["units"]


def test_load_records_the_climt_version_and_the_configuration(states, tmp_path):
    """Provenance travels with the state or the state means nothing."""
    components = _page11_components()
    state = climt.get_default_state(components,
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    path = states.save(str(tmp_path / "probe.npz"), state,
                       dict(table="earth_low_res_lw", nz=28, dt_hours=1.0,
                            n_steps=1700, slab_depth_m=2.0, solar=SOLAR,
                            co2_ppm=330.0,
                            components=["CorkLongwaveRadiation", "SlabSurface"],
                            toa_imbalance=0.02, surface_imbalance=-0.01))
    _, provenance = states.load(path, _page11_components(),
                                grid_state=get_grid(nx=1, ny=1, nz=28))

    assert provenance["climt_version"] == climt.__version__
    assert provenance["saved_at"]                      # ISO timestamp, non-empty
    assert provenance["table"] == "earth_low_res_lw"
    assert provenance["nz"] == 28
    assert provenance["components"][0] == "CorkLongwaveRadiation"


def test_load_fails_loudly_on_a_dimension_mismatch(states, tmp_path):
    """The whole reason load rebuilds instead of unpickling.

    A state saved at nz=28 must not quietly load into an nz=20 grid.
    """
    components = _page11_components()
    state = climt.get_default_state(components,
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    path = states.save(str(tmp_path / "probe.npz"), state, {})

    with pytest.raises(ValueError, match="air_temperature"):
        states.load(path, _page11_components(),
                    grid_state=get_grid(nx=1, ny=1, nz=20))


def test_load_fails_loudly_on_a_missing_quantity(states, tmp_path):
    """A component list that needs something the file does not carry."""
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    thin_state = climt.get_default_state(
        [longwave], grid_state=get_grid(nx=1, ny=1, nz=28))
    path = states.save(str(tmp_path / "thin.npz"), thin_state, {})

    with pytest.raises(ValueError, match="not in the saved state"):
        states.load(path, _page11_components(),
                    grid_state=get_grid(nx=1, ny=1, nz=28))


def test_load_returns_none_for_a_missing_asset(states):
    """A page whose data asset did not stage degrades, it does not crash."""
    assert states.load("no_such_equilibrium.npz", _page11_components(),
                       grid_state=get_grid(nx=1, ny=1, nz=28)) == (None, None)


def test_describe_names_the_configuration(states, tmp_path):
    text = states.describe(dict(
        climt_version="0.31.0", saved_at="2026-08-20T10:00:00",
        table="earth_low_res_lw", nz=28, dt_hours=1.0, n_steps=1700,
        slab_depth_m=2.0, solar=240.0, co2_ppm=330.0,
        wind_m_s=5.0, wind_timescale_hours=24.0, roughness_length_m=1e-3,
        components=["CorkLongwaveRadiation", "SlabSurface"],
        toa_imbalance=0.02, surface_imbalance=-0.01))
    for token in ("0.31.0", "earth_low_res_lw", "28", "1700", "240", "330",
                  "5.0 m/s"):
        assert token in text
