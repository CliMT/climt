"""The science the Modelling Tour pages claim, checked natively.

Pyodide cells cannot run in CI, so each page's computational core lives in
``docs/modelling-tour/_tour/`` as importable Python and is exercised here.

``docs/modelling-tour`` contains a hyphen and is not a valid package path, so
the helpers are loaded by file path.
"""
import contextlib
import copy
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
# name; _tour modules import each other the same way. Make that work here too,
# so a module loaded by path can still `import assets`.
#
# Appended, not inserted at 0: this stays on sys.path for the whole pytest
# session, and `_tour/tables.py` would otherwise shadow PyTables (import name
# `tables`, a common transitive dependency of the pandas/xarray HDF5 stack)
# for every test that runs after this file is collected. Appending still
# resolves the bare `import assets`, while a real installed package wins.
if str(TOUR) not in sys.path:
    sys.path.append(str(TOUR))


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
    """Both of page 7's integrations, run once and shared.

    Restores the backend itself. Being module-scoped, this fixture is set up
    before the function-scoped autouse ones -- including
    ``conftest.reset_sympl_backend``, which would otherwise capture the
    ``UnytBackend`` set here as its "saved" value and restore *to* it for the
    remainder of the session, leaking into every later file that sets no
    backend of its own.
    """
    saved_backend = sympl.get_backend()
    sympl.set_backend(climt.UnytBackend())
    try:
        stepping_module = _load("stepping")
        timestep = climt.UnytTimeDelta(hours=12)

        gray_components, gray_state = _page7_gray_column()
        # `_gray_column` is the shared stepping helper, named for its default
        # table; handed the non-grey table it builds the non-grey column.
        nongrey_components, nongrey_state = _gray_column(
            table=PAGE7_NONGREY_TABLE)
        yield {
            PAGE7_GRAY_TABLE: stepping_module.integrate(
                gray_components, [], gray_state, timestep, PAGE7_STEPS),
            PAGE7_NONGREY_TABLE: stepping_module.integrate(
                nongrey_components, [], nongrey_state, timestep, PAGE7_STEPS),
        }
    finally:
        sympl.set_backend(saved_backend)


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


# Page 7's reveal. `PAGE7_GRAY_TABLE`, `PAGE7_DIFFUSIVITY` and
# `_page7_gray_column` above are its configuration; these are chapter 8's two
# analytic constants, the ones page 4 prescribed and page 7 must recover.
PAGE7_TAU_INF = 4.0
PAGE7_TE = 255.0


@pytest.fixture(scope="module")
def page7_gray_equilibrium():
    """The gray column, stepped to equilibrium once and shared.

    900 steps at 12 h, 450 simulated days: page 7's own run. Restores the
    backend itself, for the reason ``page7_columns`` gives.
    """
    saved_backend = sympl.get_backend()
    sympl.set_backend(climt.UnytBackend())
    try:
        components, state = _page7_gray_column()
        yield _load("stepping").integrate(
            components, [], state, climt.UnytTimeDelta(hours=12), 900)
    finally:
        sympl.set_backend(saved_backend)


@pytest.mark.slow
def test_page7_gray_column_finds_page4s_analytic_profile(
        page7_gray_equilibrium, soundings):
    """Page 7's reveal, and the strongest claim in the tranche.

    Page 4 checked that the heating rate vanishes on a profile it was handed.
    Page 7 hands the model nothing and gets the same profile back. Measured
    max |dT| is 0.103 K; the threshold is 0.5 K, comfortable but not vacuous.
    """
    state = page7_gray_equilibrium
    p = state["air_pressure"].values[:, 0, 0]
    surface_pressure = float(state["surface_air_pressure"].values.ravel()[0])

    T_analytic, T_ground_analytic, _ = soundings.analytic_gray_equilibrium(
        p, surface_pressure, tau_inf=PAGE7_TAU_INF, T_e=PAGE7_TE)
    T = state["air_temperature"].values[:, 0, 0]

    assert np.max(np.abs(T - T_analytic)) < 0.5, (
        f"max |T - T(tau)| = {np.max(np.abs(T - T_analytic)):.3f} K; the "
        "column did not find chapter 8's profile")
    surface = _surface_temperature(state)
    assert abs(surface - T_ground_analytic) < 0.5, (
        f"surface {surface:.2f} K vs analytic {T_ground_analytic:.2f} K")


@pytest.mark.slow
def test_page7_top_level_sits_near_the_skin_temperature(
        page7_gray_equilibrium):
    """The isothermal top the figure labels is the skin temperature.

    Measured +0.63 K above Te/2^(1/4); the threshold is 1.0 K. The residual is
    real and physical -- the top level is a finite layer, not the tau -> 0
    limit -- so this is not a tolerance to tighten toward zero.
    """
    T_top = float(page7_gray_equilibrium["air_temperature"].values[-1, 0, 0])
    skin = PAGE7_TE / 2.0 ** 0.25
    assert abs(T_top - skin) < 1.0, (
        f"top level {T_top:.2f} K vs skin temperature {skin:.2f} K")


@pytest.mark.slow
def test_page7_gray_column_reaches_energy_balance(page7_gray_equilibrium,
                                                  budgets):
    """The convergence claim, read off the budget rather than a curve."""
    assert abs(budgets.toa_imbalance(page7_gray_equilibrium)) < 0.2


@pytest.mark.slow
def test_page7_mixed_layer_depth_changes_speed_not_equilibrium():
    """Page 7's knob, and the cleanest lesson in the tranche.

    1 m and 5 m slabs reach the *same* equilibrium at very different speeds.
    Heat capacity sets response time; it does not set where you end up.
    """
    saved_backend = sympl.get_backend()
    sympl.set_backend(climt.UnytBackend())
    try:
        stepping_module = _load("stepping")
        finals = {}
        for depth in (1.0, 5.0):
            components, state = _page7_gray_column(slab_depth=depth)
            stepping_module.integrate(components, [], state,
                                      climt.UnytTimeDelta(hours=12), 1400)
            finals[depth] = _surface_temperature(state)
    finally:
        sympl.set_backend(saved_backend)

    assert abs(finals[1.0] - finals[5.0]) < 0.5, (
        f"1 m settled at {finals[1.0]:.2f} K and 5 m at {finals[5.0]:.2f} K — "
        "slab depth must not change the equilibrium")


# ------------------------------------------------------------------ page 8
#
# Page 7's gray column, plus SimpleBoundaryLayer(surface_fluxes='bulk') and a
# wind held up by stepping.wind_relaxation. dt = 1 h, NOT page 7's 12 h: the
# turbulent column's equilibrium depends on the timestep (sensible heat flux
# 34.7 W/m^2 at 1 h, 34.2 at 30 min, 39.6 at 12 h, 30-day means), so the page steps hourly.
# The radiative-only column's does not (15.15 K at both), so the page runs
# that one at page 7's 12 h.
#
# Numbers are 30-day means, as on the page: at equilibrium the boundary layer
# deepens for a step every six hours or so, and the jump and the flux flicker
# with it by about +-0.3 K and +-2 W/m^2. Measured with
# scripts/experiments/tour_page8_measurements.py.

PAGE8_Z0 = 1e-3
PAGE8_WIND = 5.0
PAGE8_STEPS = 12000          # dt = 1 h, 500 days; the cold start converges
PAGE8_RESTART_STEPS = 1000   # the knob cells restart from the headline
                             # equilibrium and re-settle in this many steps
PAGE8_DT = dict(hours=1)

# (target wind, minimum sensible heat flux, maximum jump). Measured 30-day
# means: 20.0, 34.7, 51.2 W/m^2 and 9.33, 7.47, 6.16 K. The thresholds allow
# ~25% -- they check the monotone response, not a fit.
PAGE8_WIND_CASES = [
    (2.0, 13.5, 9.5),
    (5.0, 23.0, 7.5),
    (10.0, 31.0, 5.5),
]


@contextlib.contextmanager
def _unyt_backend_restored():
    """UnytBackend for the duration, and whatever was set before, after."""
    saved_backend = sympl.get_backend()
    sympl.set_backend(climt.UnytBackend())
    try:
        yield
    finally:
        sympl.set_backend(saved_backend)


def _page8_column(z0=PAGE8_Z0, wind=PAGE8_WIND, slab_depth=2.0, nz=28,
                  boundary_layer=True):
    """Page 8's column: page 7's gray column, a boundary layer, and a wind.

    ``wind=None`` omits the momentum source, which is how the page shows what
    happens without one -- the column spins itself down to dead calm.
    """
    stepping_module = _load("stepping")
    longwave = climt.CorkLongwaveRadiation(
        optics="correlated_k", table=PAGE7_GRAY_TABLE,
        diffusivity_factor=PAGE7_DIFFUSIVITY)
    surface = climt.SlabSurface()
    steppers = ([climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                           roughness_length=z0)]
                if boundary_layer else [])
    state = climt.get_default_state([longwave, surface] + steppers,
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    tendencies = [longwave, surface]
    if boundary_layer and wind is not None:
        tendencies.append(stepping_module.wind_relaxation(state, wind))
    return tendencies, steppers, state


def _lowest_wind(state):
    return float(state["eastward_wind"].values[0, 0, 0])


def _page8_recorder(stepping_module, budgets_module):
    return stepping_module.Recorder(
        climt.UnytTimeDelta(**PAGE8_DT),
        jump=budgets_module.surface_air_jump,
        flux=budgets_module.sensible_heat_flux,
        wind=_lowest_wind)


@pytest.fixture(scope="module")
def page8_equilibrium():
    """The headline turbulent column at 5 m/s, cold-started to equilibrium.

    Shared by every page-8 test that, like the page's own knob cells, starts
    from it. Returns ``(tendencies, steppers, state, record)``; tests must
    deep-copy the state before stepping it.
    """
    with _unyt_backend_restored():
        stepping_module = _load("stepping")
        record = _page8_recorder(stepping_module, _load("budgets"))
        tendencies, steppers, state = _page8_column()
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(**PAGE8_DT),
                                  PAGE8_STEPS, after_step=record)
        yield tendencies, steppers, state, record


def _page8_restart(page8_equilibrium, wind=PAGE8_WIND, z0=PAGE8_Z0,
                   reset=None):
    """What the page's knob cells do: copy the equilibrium, change one thing,
    and step PAGE8_RESTART_STEPS. ``reset`` assigns the wind back every step
    instead of relaxing it. Call inside ``_unyt_backend_restored``."""
    stepping_module = _load("stepping")
    tendencies, steppers, base, _ = page8_equilibrium
    state = copy.deepcopy(base)
    longwave, surface = tendencies[:2]
    if z0 != PAGE8_Z0:
        steppers = [climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                              roughness_length=z0)]
    record = _page8_recorder(stepping_module, _load("budgets"))
    if reset is None:
        here = [longwave, surface,
                stepping_module.wind_relaxation(state, wind,
                                                initialise=False)]
        after_step = record
    else:
        here = [longwave, surface]

        def after_step(state):
            state["eastward_wind"].values[:] = reset
            record(state)
    stepping_module.integrate(here, steppers, state,
                              climt.UnytTimeDelta(**PAGE8_DT),
                              PAGE8_RESTART_STEPS, after_step=after_step)
    return state, record


@pytest.mark.slow
def test_page8_turbulence_shrinks_the_surface_air_discontinuity(
        page8_equilibrium):
    """Page 8's reveal. Measured 15.15 K -> 5.95 K (30-day mean) at z0=1e-3,
    5 m/s.

    Thresholds sit either side of the measurement with room, so this checks
    that turbulence erodes the discontinuity, not by exactly how much. The
    radiative column runs at 12 h, as on the page -- its equilibrium does not
    depend on the timestep.
    """
    with _unyt_backend_restored():
        budgets_module = _load("budgets")
        tendencies, steppers, state = _page8_column(boundary_layer=False)
        _load("stepping").integrate(tendencies, steppers, state,
                                    climt.UnytTimeDelta(hours=12), 900)
        radiative = budgets_module.surface_air_jump(state)
    turbulent = page8_equilibrium[3].mean("jump")

    assert radiative > 13.0, (
        f"radiative-equilibrium jump {radiative:.2f} K — page 4's "
        "discontinuity should be there before anything erodes it")
    assert turbulent < radiative - 4.0, (
        f"turbulent jump {turbulent:.2f} K vs radiative {radiative:.2f} K — "
        "the boundary layer should erode it")


@pytest.mark.slow
def test_page8_turbulent_column_is_converged(page8_equilibrium):
    """The headline is an equilibrium, read off the budget."""
    with _unyt_backend_restored():
        imbalance = _load("budgets").toa_imbalance(page8_equilibrium[2])
    assert abs(imbalance) < 0.1, f"TOA imbalance {imbalance:+.3f} W/m^2"


@pytest.mark.slow
def test_page8_a_column_with_no_momentum_source_spins_itself_down():
    """Why the page needs a wind forcing at all.

    Nothing drives a wind in a column with no dynamics, and the boundary
    layer's own surface drag removes any wind prescribed as an initial
    condition. Measured: three initialisations spanning 0-10 m/s all end at a
    lowest-level wind of 0.000 and the same equilibrium.
    """
    with _unyt_backend_restored():
        stepping_module = _load("stepping")
        budgets_module = _load("budgets")
        finals = {}
        for wind in (0.0, 5.0, 10.0):
            tendencies, steppers, state = _page8_column(wind=None)
            state["eastward_wind"].values[:] = wind
            record = _page8_recorder(stepping_module, budgets_module)
            stepping_module.integrate(tendencies, steppers, state,
                                      climt.UnytTimeDelta(**PAGE8_DT),
                                      PAGE8_STEPS, after_step=record)
            finals[wind] = (record.mean("jump"), _lowest_wind(state))

    jumps = [jump for jump, _ in finals.values()]
    assert max(jumps) - min(jumps) < 0.5, (
        f"initial wind changed the equilibrium jump: {finals} — without a "
        "momentum source it must not, and page 8's argument for adding one "
        "depends on that")
    for _, lowest_wind in finals.values():
        assert abs(lowest_wind) < 0.1, (
            f"lowest-level wind {lowest_wind:.3f} m/s — the drag should have "
            "removed it entirely")


@pytest.mark.slow
@pytest.mark.parametrize("wind, min_flux, max_jump", PAGE8_WIND_CASES)
def test_page8_wind_speed_is_a_knob_once_it_is_held_up(
        page8_equilibrium, wind, min_flux, max_jump):
    """Page 8's knob, over the range the page says was tested, by the page's
    own route: restart from the 5 m/s equilibrium at the new wind.

    Thresholds allow ~25% on the measured 30-day means -- this checks the
    monotone response, not a fit.
    """
    with _unyt_backend_restored():
        if wind == PAGE8_WIND:
            record = page8_equilibrium[3]
        else:
            _, record = _page8_restart(page8_equilibrium, wind=wind)
    flux, jump = record.mean("flux"), record.mean("jump")
    assert flux > min_flux, f"wind={wind}: sensible heat flux {flux:.2f} W/m^2"
    assert jump < max_jump, f"wind={wind}: jump {jump:.2f} K"


@pytest.mark.slow
def test_page8_restart_lands_where_a_cold_start_does(page8_equilibrium):
    """The knob cells' shortcut is honest.

    They restart from the headline equilibrium and step PAGE8_RESTART_STEPS
    rather than cold-starting 12 000. At 10 m/s, the largest change the page
    makes, the two must agree on the 30-day-mean jump and flux.
    """
    with _unyt_backend_restored():
        _, restarted = _page8_restart(page8_equilibrium, wind=10.0)
        stepping_module = _load("stepping")
        record = _page8_recorder(stepping_module, _load("budgets"))
        tendencies, steppers, state = _page8_column(wind=10.0)
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(**PAGE8_DT),
                                  PAGE8_STEPS, after_step=record)
    assert restarted.mean("jump") == pytest.approx(record.mean("jump"),
                                                   abs=0.15)
    assert restarted.mean("flux") == pytest.approx(record.mean("flux"),
                                                   abs=1.5)


@pytest.mark.slow
def test_page8_rougher_surface_moves_more_heat(page8_equilibrium):
    """The second knob, at the page's own wind, over the range the page says
    was tested. Asserts the monotone ordering of the 30-day-mean flux."""
    with _unyt_backend_restored():
        fluxes = {PAGE8_Z0: page8_equilibrium[3].mean("flux")}
        for z0 in (3.21e-5, 1e-1):
            fluxes[z0] = _page8_restart(page8_equilibrium, z0=z0)[1].mean(
                "flux")
    assert fluxes[3.21e-5] < fluxes[1e-3] < fluxes[1e-1], (
        f"sensible heat flux by z0: {fluxes} — a rougher surface must move "
        "more heat")


@pytest.mark.slow
def test_page8_relaxation_lets_the_drag_win_near_the_ground(
        page8_equilibrium):
    """The comparison that distinguishes a forcing from an assignment.

    Relaxing toward 5 m/s leaves the lowest level near 3.4, because the
    surface drag is still acting and the relaxation only pulls. Assigning the
    wind back every step -- the idiom in examples/column_code_with_slab.py --
    pins it at exactly 5.0 and so overrides the drag at the one level where
    the exchange happens, delivering more flux from the same nominal wind.
    """
    relaxed = page8_equilibrium[3]
    with _unyt_backend_restored():
        reset_state, reset = _page8_restart(page8_equilibrium,
                                            reset=PAGE8_WIND)

    assert 1.0 < relaxed.mean("wind") < 4.5, (
        f"relaxed lowest-level wind {relaxed.mean('wind'):.2f} m/s — it "
        "should sit well below the 5 m/s target, because the drag is still "
        "acting")
    assert _lowest_wind(reset_state) == pytest.approx(PAGE8_WIND)
    assert reset.mean("flux") > relaxed.mean("flux") + 3.0, (
        f"hard reset {reset.mean('flux'):.1f} vs relaxation "
        f"{relaxed.mean('flux'):.1f} W/m^2 — pinning the lowest level "
        "overrides the drag where the exchange happens, so it must move more "
        "heat")


def test_page8_no_flux_mode_conserves_the_column():
    """The third mode: with surface_fluxes=None the diffusion conserves.

    Cheap, so unmarked: ten steps is enough, because conservation is exact
    rather than asymptotic.
    """
    with _unyt_backend_restored():
        budgets_module = _load("budgets")
        longwave = climt.CorkLongwaveRadiation(
            optics="correlated_k", table=PAGE7_GRAY_TABLE,
            diffusivity_factor=PAGE7_DIFFUSIVITY)
        boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes=None)
        state = climt.get_default_state([longwave, boundary_layer],
                                        grid_state=get_grid(nx=1, ny=1, nz=28))
        state["specific_humidity"].values[:] = 4e-3
        state["eastward_wind"].values[:] = 5.0
        before = budgets_module.column_enthalpy(state)

        timestep = climt.UnytTimeDelta(**PAGE8_DT)
        for _ in range(10):
            diagnostics, new_state = boundary_layer(state, timestep)
            state.update(new_state)
            state.update(diagnostics)

        after = budgets_module.column_enthalpy(state)
    assert after == pytest.approx(before, rel=1e-9), (
        "surface_fluxes=None must conserve every column integral exactly")


def test_page8_helpers_record_and_reset_every_step():
    """``integrate(after_step=...)`` runs once per step, after the clock, and
    ``Recorder`` keeps what it is given. Cheap and unmarked: the slow tests
    above all lean on these two."""
    with _unyt_backend_restored():
        stepping_module = _load("stepping")
        budgets_module = _load("budgets")
        tendencies, steppers, state = _page8_column(wind=None)
        record = stepping_module.Recorder(
            climt.UnytTimeDelta(**PAGE8_DT), profiles=["eastward_wind"],
            profile_every=2, jump=budgets_module.surface_air_jump,
            flux=budgets_module.sensible_heat_flux, wind=_lowest_wind)

        def after_step(state):
            state["eastward_wind"].values[:] = 7.0
            record(state)

        start = state["time"]
        stepping_module.integrate(tendencies, steppers, state,
                                  climt.UnytTimeDelta(**PAGE8_DT), 5,
                                  after_step=after_step)
        elapsed = float((state["time"] - start).total_seconds())

    assert elapsed == pytest.approx(5 * 3600.0)
    np.testing.assert_allclose(record["days"], np.arange(1, 6) / 24.0)
    np.testing.assert_allclose(record["wind"], 7.0)
    assert record["flux"][-1] == pytest.approx(
        float(state["surface_upward_sensible_heat_flux"].values.ravel()[0]))
    assert record["eastward_wind"].shape == (3, 28)
    np.testing.assert_allclose(record["profile_days"], [1 / 24, 3 / 24,
                                                        5 / 24])
    assert record.mean("wind", days=1.0) == pytest.approx(7.0)


def test_page8_headline_figure_draws(tmp_path):
    """``draw_discontinuity`` builds its three panels from short runs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with _unyt_backend_restored():
        stepping_module = _load("stepping")
        budgets_module = _load("budgets")
        timestep = climt.UnytTimeDelta(**PAGE8_DT)
        tendencies, steppers, radiative = _page8_column(boundary_layer=False)
        stepping_module.integrate(tendencies, steppers, radiative, timestep, 3)
        tendencies, steppers, turbulent = _page8_column()
        record = _page8_recorder(stepping_module, budgets_module)
        stepping_module.integrate(tendencies, steppers, turbulent, timestep,
                                  50, after_step=record)
        plt.close("all")
        try:
            stepping_module.draw_discontinuity(radiative, turbulent, record,
                                               title="test")
            figure = plt.gcf()
            assert len(figure.axes) == 3
            assert "jump" in "".join(
                text.get_text() for text in figure.axes[0].texts)
        finally:
            plt.close("all")


# ------------------------------------------------------------------ page 9
#
# Page 8's column with the 14-band earth_low_res_lw table, a wet surface and
# GridScaleCondensation. The surface humidity is held at a fixed relative
# humidity *of the surface's own temperature* by stepping.SurfaceHumidity --
# not at a fixed number, which the column would leave behind as it cooled
# (0.015 kg/kg is 131 % of saturation at the 289 K this column settles at).
#
# dt = 1 h, as on page 8: converged sensible heat flux 49.3 W/m^2 at 6 h,
# 46.8 at 3 h, 44.7 at 1 h, 44.1 at 30 min. The column takes ~900 days to
# converge, which is over an hour in a browser, so the page's cells run 30-day
# cold starts and quote the converged numbers from
# scripts/experiments/tour_page9_measurements.py. These tests guard both: the
# converged claims on converged runs, and the cells' claims on the cells'
# own route. Numbers are 30-day means for converged runs and 10-day means for
# the page's 30-day runs, as the page quotes them.

PAGE9_DT = dict(hours=1)
PAGE9_STEPS = 21600          # 900 days; the surface is within 0.03 K of its
                             # day-1000 value, and |TOA| < 0.25 W/m^2
PAGE9_CELL_STEPS = 720       # the page's cold-start cells: 30 days
PAGE9_RESTART_STEPS = 240    # its condensation cell: 10 more days


def _page9_column(surface_relative_humidity=1.0, nz=28, slab_depth=2.0,
                  condensation=True, wind=PAGE8_WIND):
    """Page 9's column: 14-band radiation, boundary layer, wet surface.

    Carries page 8's wind relaxation. Every configuration from page 8 on does:
    a column with no momentum source spins itself down to dead calm, and its
    surface fluxes -- including the latent flux this page is about -- are then
    those of a windless planet.
    """
    stepping_module = _load("stepping")
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                           table="earth_low_res_lw")
    surface = climt.SlabSurface()
    steppers = [stepping_module.SurfaceHumidity(surface_relative_humidity),
                climt.SimpleBoundaryLayer(surface_fluxes="bulk",
                                          roughness_length=PAGE8_Z0)]
    if condensation:
        steppers.append(climt.GridScaleCondensation())
    state = climt.get_default_state([longwave, surface] + steppers,
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    state["ocean_mixed_layer_thickness"].values[:] = slab_depth
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    relaxation = stepping_module.wind_relaxation(state, wind)
    return [longwave, surface, relaxation], steppers, state


def _virtual_potential_temperature(state):
    """theta_v = theta * (1 + 0.61 q), from state quantities."""
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    Cp = float(get_constant("heat_capacity_of_dry_air_at_constant_pressure",
                            "J/kg/degK"))
    p = state["air_pressure"].values[:, 0, 0]
    T = state["air_temperature"].values[:, 0, 0]
    q = state["specific_humidity"].values[:, 0, 0]
    theta = T * (1.0e5 / p) ** (Rd / Cp)
    return theta, theta * (1.0 + 0.61 * q)


def _peak_relative_humidity(state):
    """Largest relative humidity in the column, by the saturation formula
    GridScaleCondensation uses (``soundings.saturation_specific_humidity``),
    so the page, this test and the measurement script all agree."""
    return float(np.max(_load("soundings").relative_humidity(state)))


def _page9_recorder(stepping_module, budgets_module):
    return stepping_module.Recorder(
        climt.UnytTimeDelta(**PAGE9_DT),
        sh=budgets_module.sensible_heat_flux,
        lh=budgets_module.latent_heat_flux,
        ts=_surface_temperature,
        blh=lambda state: float(
            state["boundary_layer_height"].values.ravel()[0]),
        peak_rh=_peak_relative_humidity)


def _page9_run(n_steps, **column_kwargs):
    """A cold start of ``n_steps``. Call inside ``_unyt_backend_restored``.
    Returns ``(tendencies, steppers, state, record)``."""
    stepping_module = _load("stepping")
    record = _page9_recorder(stepping_module, _load("budgets"))
    tendencies, steppers, state = _page9_column(**column_kwargs)
    stepping_module.integrate(tendencies, steppers, state,
                              climt.UnytTimeDelta(**PAGE9_DT), n_steps,
                              after_step=record)
    return tendencies, steppers, state, record


def _bowen(record, days):
    return record.mean("sh", days) / record.mean("lh", days)


@pytest.fixture(scope="module")
def page9_equilibrium():
    """Page 9's column at surface relative humidity 1.0, converged."""
    with _unyt_backend_restored():
        yield _page9_run(PAGE9_STEPS)


@pytest.fixture(scope="module")
def page9_cell():
    """What the page's headline cell leaves behind: a 30-day cold start."""
    with _unyt_backend_restored():
        yield _page9_run(PAGE9_CELL_STEPS)


def test_page9_virtual_potential_temperature_identity():
    """theta_v - theta = 0.61 q theta, to machine precision.

    Cheap and exact: this is the definition, and the page derives it in a
    cell. If a reader's arithmetic disagrees with the page, it is the page
    that is wrong, so pin it.
    """
    with _unyt_backend_restored():
        tendencies, steppers, state = _page9_column()
        state["specific_humidity"].values[:] = 8e-3
        theta, theta_v = _virtual_potential_temperature(state)
    np.testing.assert_allclose(theta_v - theta, 0.61 * 8e-3 * theta,
                               rtol=1e-12)


def test_page9_surface_humidity_follows_the_surface_temperature():
    """``SurfaceHumidity`` writes RH * q_sat(Ts, ps) -- and q_sat is the
    saturation GridScaleCondensation condenses to, so "relative humidity 1"
    means one thing at the surface and in the air."""
    with _unyt_backend_restored():
        soundings_module = _load("soundings")
        timestep = climt.UnytTimeDelta(**PAGE9_DT)
        tendencies, steppers, state = _page9_column(
            surface_relative_humidity=0.7)
        surface_humidity, condensation = steppers[0], steppers[-1]
        q_surface = {}
        for surface_temperature in (280.0, 300.0):
            state["surface_temperature"].values[:] = surface_temperature
            _, new_state = surface_humidity(state, timestep)
            q_surface[surface_temperature] = float(np.asarray(
                new_state["surface_specific_humidity"]).ravel()[0])
            expected = 0.7 * float(soundings_module.saturation_specific_humidity(
                surface_temperature,
                float(state["surface_air_pressure"].values.ravel()[0])))
            assert q_surface[surface_temperature] == pytest.approx(
                expected, rel=1e-12)
        assert q_surface[300.0] > 2.0 * q_surface[280.0]

        # A column 2 % supersaturated everywhere below 500 hPa comes out of
        # GridScaleCondensation at relative humidity 1 by the same formula.
        # (Its adjustment is linearised in temperature, so it is exact only
        # for small excesses -- which is what one hourly step leaves it.)
        p = state["air_pressure"].values[:, 0, 0]
        T = state["air_temperature"].values[:, 0, 0]
        low = p > 5.0e4
        q = np.zeros_like(T)
        q[low] = 1.02 * soundings_module.saturation_specific_humidity(T, p)[low]
        state["specific_humidity"].values[:, 0, 0] = q
        _, condensed = condensation(state, timestep)
        state.update(condensed)
        relative_humidity = soundings_module.relative_humidity(state)
    np.testing.assert_allclose(relative_humidity[low], 1.0, atol=1e-3)


@pytest.mark.slow
def test_page9_latent_flux_exceeds_sensible_over_a_saturated_surface(
        page9_equilibrium):
    """Page 9's headline claim, at this column's own converged temperature.

    The spec phrases it "over a saturated ~288 K surface". Measured: the
    column converges at 286.21 K with a Bowen ratio of 0.20 (30-day means,
    1000 days, dt 1 h), so the claim holds as the spec states it. The page
    quotes both.
    """
    _, _, state, record = page9_equilibrium
    with _unyt_backend_restored():
        imbalance = _load("budgets").toa_imbalance(state)
    surface = record.mean("ts")
    bowen = _bowen(record, 30.0)

    assert abs(imbalance) < 0.5, f"TOA imbalance {imbalance:+.3f} W/m^2"
    assert 286.0 < surface < 292.0, (
        f"converged surface {surface:.2f} K -- the page quotes 286.21 K")
    assert bowen < 1.0, (
        f"Bowen ratio {bowen:.2f} at surface {surface:.1f} K -- over a "
        "saturated surface the latent flux should dominate")
    assert 0.15 < bowen < 0.25, f"Bowen ratio {bowen:.3f}; the page quotes 0.20"


@pytest.mark.slow
def test_page9_headline_cell_is_already_latent_dominated(page9_cell):
    """The headline cell's 30 days: far from equilibrium, and the page says
    so, but the partition it shows is already the converged one's sign.
    Measured 10-day means: SH 22.0, LH 102.4 W/m^2, Bowen 0.21, TOA -82."""
    _, _, state, record = page9_cell
    with _unyt_backend_restored():
        imbalance = _load("budgets").toa_imbalance(state)
    assert _bowen(record, 10.0) < 0.8
    assert imbalance < -50.0, (
        f"TOA {imbalance:+.1f} W/m^2 -- the page tells readers the 30-day "
        "column is far from equilibrium")


@pytest.mark.slow
def test_page9_condensation_removes_the_supersaturation(page9_cell):
    """The reason GridScaleCondensation is in the stack from this page on,
    by the page's own route: from the headline cell's state, ten more days
    with the sink and without it. Measured: peak relative humidity 100 %
    with it, 2876 % without, and the precipitation 3.57 mm/day."""
    tendencies, steppers, base, _ = page9_cell
    peaks, rain = {}, {}
    with _unyt_backend_restored():
        stepping_module = _load("stepping")
        budgets_module = _load("budgets")
        timestep = climt.UnytTimeDelta(**PAGE9_DT)
        for label, keep in (("with", True), ("without", False)):
            state = copy.deepcopy(base)
            here = steppers if keep else steppers[:2]
            if not keep:
                del state["precipitation_amount"]
            record = stepping_module.Recorder(
                timestep, peak_rh=_peak_relative_humidity,
                precip=lambda s: budgets_module.precipitation_rate(s,
                                                                   timestep))
            stepping_module.integrate(tendencies, here, state, timestep,
                                      PAGE9_RESTART_STEPS, after_step=record)
            peaks[label] = float(np.max(record["peak_rh"]))
            rain[label] = record.mean("precip", 10.0)

    # GridScaleCondensation's adjustment is linearised in temperature, so it
    # leaves ~1e-5 of the excess behind. That is saturated, for any purpose
    # this page has.
    assert peaks["with"] < 1.0 + 1e-4, (
        f"peak RH {peaks['with']:.4f} with condensation last in the stepper "
        "list -- nothing the step leaves behind should be supersaturated")
    assert peaks["without"] > 4.0, (
        f"peak RH without condensation {peaks['without']:.2f} -- the page "
        "quotes 2876 % after ten days")
    assert rain["with"] > 1.0 and rain["without"] == 0.0, rain


@pytest.mark.slow
def test_page9_a_column_with_no_moisture_sink_has_no_equilibrium():
    """What the page says happens to a sink-free column left alone: the
    vapour keeps accumulating, and the surface keeps warming. Measured at
    day 60 of a cold start: 31.3 g/kg at the lowest level, surface 310.5 K,
    both still rising; the run fails outright at day 464."""
    with _unyt_backend_restored():
        _, _, state, record = _page9_run(24 * 60, condensation=False)
    q_lowest = float(state["specific_humidity"].values[0, 0, 0])
    assert q_lowest > 0.018, f"lowest-level q {q_lowest * 1e3:.1f} g/kg"
    assert record.mean("ts", 5.0) > record["ts"][24 * 30] + 2.0, (
        "the surface should still be warming at day 60")


@pytest.mark.slow
def test_page9_surface_relative_humidity_moves_the_bowen_ratio(page9_cell):
    """Page 9's knob, by the knob cell's route (30-day cold starts), over the
    range the page says was tested. Measured 10-day means: Bowen 0.89 at RH
    0.4 and 0.21 at 1.0 -- the partition moves toward sensible heat (it flips
    only once the column has converged)."""
    with _unyt_backend_restored():
        _, _, _, dry = _page9_run(PAGE9_CELL_STEPS,
                                  surface_relative_humidity=0.4)
    wet = page9_cell[3]
    assert _bowen(dry, 10.0) > 3.0 * _bowen(wet, 10.0), (
        f"Bowen ratio {_bowen(dry, 10.0):.2f} at RH 0.4, "
        f"{_bowen(wet, 10.0):.2f} at RH 1.0 -- a drier surface partitions "
        "more into sensible heat")


@pytest.mark.slow
def test_page9_a_drier_surface_ends_colder(page9_equilibrium):
    """The knob's converged half, which the page quotes because it reverses
    the 30-day answer: at RH 0.4 the column settles at 283.45 K, 2.76 K
    colder than at RH 1.0, at Bowen 1.48. Less vapour, less greenhouse."""
    with _unyt_backend_restored():
        _, _, _, dry = _page9_run(PAGE9_STEPS, surface_relative_humidity=0.4)
    wet = page9_equilibrium[3]
    assert _bowen(dry, 30.0) > 1.0
    assert dry.mean("ts") < wet.mean("ts") - 1.5, (
        f"RH 0.4 at {dry.mean('ts'):.2f} K vs RH 1.0 at "
        f"{wet.mean('ts'):.2f} K -- the page says the drier column ends "
        "colder")


@pytest.mark.slow
def test_page9_the_surface_supplies_a_greenhouse_and_a_deeper_mixed_layer(
        page9_equilibrium):
    """Two converged comparisons against the same column over a dry surface.

    Measured: dry 266.69 K against moist 286.21 K, the water-vapour
    greenhouse the surface supplied; and a boundary layer about 8.1 km deep
    dry and 12.0 km moist (medians, a habit from when the depth jumped for
    single steps; it now holds steady).
    """
    with _unyt_backend_restored():
        _, _, _, dry = _page9_run(PAGE9_STEPS, surface_relative_humidity=0.0)
    moist = page9_equilibrium[3]
    warming = moist.mean("ts") - dry.mean("ts")
    assert 18.0 < warming < 28.0, (
        f"moist minus dry surface temperature {warming:.2f} K -- the page "
        "quotes about 20 K")

    def median_depth(record):
        days = record["days"]
        return float(np.median(record["blh"][days > days[-1] - 30.0]))

    assert median_depth(moist) > median_depth(dry) + 2000.0, (
        f"median boundary-layer depth {median_depth(moist):.0f} m moist vs "
        f"{median_depth(dry):.0f} m dry")


def test_page9_headline_figure_draws():
    """``draw_moisture`` builds its three panels from a short run."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with _unyt_backend_restored():
        _, _, state, record = _page9_run(50)
        plt.close("all")
        try:
            _load("stepping").draw_moisture(state, record, title="test")
            figure = plt.gcf()
            assert len(figure.axes) == 3
            assert "θ" in "".join(t.get_text()
                                  for t in figure.axes[0].texts)
        finally:
            plt.close("all")


# --- Page 10: dry convection ------------------------------------------------
#
# One call, no time loop. The profile is the page's: a constant lapse rate
# written exactly in pressure, T = Ts (p/ps)^(Rd Gamma/g), superadiabatic in the
# lowest layers and continuous into 6.5 K/km above them, floored at 200 K.
# Numbers are those the page quotes; its cells print them.

PAGE10_DT = dict(hours=1)


def _page10_column(unstable_levels=8, gamma_unstable=14e-3, gamma=6.5e-3,
                   T_surface=300.0, T_min=200.0, specific_humidity=0.0,
                   nz=28):
    """Page 10's prescribed profile: superadiabatic near the ground."""
    from sympl import get_constant

    adjustment = climt.DryConvectiveAdjustment()
    state = climt.get_default_state([adjustment],
                                    grid_state=get_grid(nx=1, ny=1, nz=nz))
    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    g = float(get_constant("gravitational_acceleration", "m/s^2"))
    p = state["air_pressure"].to_units("Pa").values[:, 0, 0]
    ps = float(state["surface_air_pressure"].to_units("Pa").values.ravel()[0])

    T = T_surface * (p / ps) ** (Rd * gamma_unstable / g)
    top = unstable_levels - 1
    T[unstable_levels:] = (T[top] * (p[unstable_levels:] / p[top])
                           ** (Rd * gamma / g))
    state["air_temperature"].values[:, 0, 0] = np.maximum(T, T_min)
    state["specific_humidity"].values[:] = specific_humidity
    return adjustment, state


def _page10_theta(state):
    """Dry potential temperature, K: page 9's helper, without the humidity."""
    return _virtual_potential_temperature(state)[0]


def _page10_adjust(**column_kwargs):
    """The page's one call: (state before, state after, levels adjusted)."""
    adjustment, before = _page10_column(**column_kwargs)
    diagnostics, new_state = adjustment(before,
                                        climt.UnytTimeDelta(**PAGE10_DT))
    after = copy.deepcopy(before)
    after.update(new_state)
    adjusted = np.abs(after["air_temperature"].values[:, 0, 0]
                      - before["air_temperature"].values[:, 0, 0]) > 1e-6
    return before, after, adjusted


def test_page10_the_prescribed_profile_is_unstable_where_the_page_says():
    """theta falls with height through the eight superadiabatic layers only."""
    with _unyt_backend_restored():
        _, state = _page10_column()
        theta = _page10_theta(state)
    assert np.all(np.diff(theta[:8]) < 0)
    assert np.all(np.diff(theta[7:]) > 0)
    assert theta[0] == pytest.approx(298.78, abs=0.005)
    assert theta[7] == pytest.approx(291.83, abs=0.005)


def test_page10_adjustment_makes_theta_uniform():
    """Page 10's first claim: ten layers share one theta, 294.80 K, and the
    column above them is left stable."""
    with _unyt_backend_restored():
        before, after, adjusted = _page10_adjust()
        theta = _page10_theta(after)
        p = after["air_pressure"].to_units("Pa").values[:, 0, 0]
        dT = (after["air_temperature"].values[:, 0, 0]
              - before["air_temperature"].values[:, 0, 0])

    assert adjusted.sum() == 10, (
        f"{adjusted.sum()} levels adjusted; the page says ten -- the eight "
        "unstable layers and two of the stable air above them")
    assert np.array_equal(np.where(adjusted)[0], np.arange(10))
    assert p[adjusted][-1] / 100 == pytest.approx(745.0, abs=0.5)
    spread = float(theta[adjusted].max() - theta[adjusted].min())
    assert spread < 1e-9, (
        f"theta spread {spread:.2e} K across the adjusted layers -- an "
        "adjusted layer is by definition neutrally stratified")
    assert theta[0] == pytest.approx(294.80, abs=0.005)
    assert np.all(np.diff(theta[9:]) > 0), "the air above is left stable"
    assert -dT[0] == pytest.approx(3.99, abs=0.005)
    assert dT.max() == pytest.approx(2.82, abs=0.005)
    assert np.argmax(dT) == 7


def test_page10_adjustment_conserves_column_enthalpy(budgets):
    """Page 10's second claim, and the test of whether you understood it.

    Measured relative change: exactly 0.0 on the page's profile, JIT on or
    off. tests/test_conservation.py::TestDryConvectionConservation asserts
    the same thing through a different route; this one asserts it on page
    10's own profile, with budgets.column_enthalpy -- the function the page
    calls.
    """
    with _unyt_backend_restored():
        before, after, _ = _page10_adjust()
        h_before = budgets.column_enthalpy(before)
        h_after = budgets.column_enthalpy(after)
    assert h_after == pytest.approx(h_before, rel=1e-14), (
        f"enthalpy changed by {abs(h_after - h_before) / h_before:.2e} "
        "relative -- dry adjustment mixes, it does not heat")


def test_page10_a_stable_column_is_left_alone():
    """The control. Nothing to adjust means nothing adjusted."""
    with _unyt_backend_restored():
        adjustment, state = _page10_column(gamma_unstable=6.5e-3)
        before = state["air_temperature"].values.copy()
        diagnostics, new_state = adjustment(state,
                                            climt.UnytTimeDelta(**PAGE10_DT))
        np.testing.assert_allclose(new_state["air_temperature"].values,
                                   before)


def test_page10_adjustment_returns_a_state_not_tendencies():
    """The craft claim: a Stepper's signature is different, and visibly so."""
    with _unyt_backend_restored():
        adjustment, state = _page10_column()
        before = state["air_temperature"].values.copy()
        result = adjustment(state, climt.UnytTimeDelta(**PAGE10_DT))
        after_call = state["air_temperature"].values.copy()
    assert isinstance(result, tuple) and len(result) == 2
    diagnostics, new_state = result
    assert dict(diagnostics) == {}
    assert sorted(new_state) == ["air_temperature", "specific_humidity"]
    assert not any("tendency" in key for key in new_state)
    np.testing.assert_array_equal(after_call, before)   # input left alone


def test_page10_the_timestep_is_ignored():
    """The craft callout: adjustment is instantaneous, whatever DT is."""
    with _unyt_backend_restored():
        adjustment, state = _page10_column()
        _, one_hour = adjustment(state, climt.UnytTimeDelta(hours=1))
        _, one_day = adjustment(state, climt.UnytTimeDelta(days=1))
        np.testing.assert_array_equal(one_hour["air_temperature"].values,
                                      one_day["air_temperature"].values)


@pytest.mark.parametrize("levels, gamma, n_adjusted, top_hPa, max_dT", [
    (4, 10e-3, 4, 968, 0.05), (4, 14e-3, 5, 942, 0.92),
    (4, 20e-3, 5, 942, 2.58),
    (8, 10e-3, 8, 836, 0.21), (8, 14e-3, 10, 745, 3.99),
    (8, 20e-3, 11, 695, 10.69),
    (12, 10e-3, 12, 643, 0.46), (12, 14e-3, 14, 534, 8.93),
    (12, 20e-3, 17, 370, 23.12),
])
def test_page10_knob_reach(levels, gamma, n_adjusted, top_hPa, max_dT):
    """The knob table: how far the adjusted layer reaches."""
    with _unyt_backend_restored():
        before, after, adjusted = _page10_adjust(unstable_levels=levels,
                                                 gamma_unstable=gamma)
        p = after["air_pressure"].to_units("Pa").values[:, 0, 0]
        dT = np.abs(after["air_temperature"].values[:, 0, 0]
                    - before["air_temperature"].values[:, 0, 0])
    assert adjusted.sum() == n_adjusted
    assert p[adjusted][-1] / 100 == pytest.approx(top_hPa, abs=0.5)
    assert dT.max() == pytest.approx(max_dT, abs=0.005)


def test_page10_knob_extremes():
    """The prose around the knob table: the mixed theta at 20 K/km over twelve
    layers, and the physics exercise's floor-limited deep instability."""
    with _unyt_backend_restored():
        _, after, _ = _page10_adjust(unstable_levels=12, gamma_unstable=20e-3)
        assert _page10_theta(after)[0] == pytest.approx(275.59, abs=0.005)
        for levels in (20, 24, 28):
            _, after, adjusted = _page10_adjust(unstable_levels=levels)
            p = after["air_pressure"].to_units("Pa").values[:, 0, 0]
            assert p[adjusted][-1] / 100 == pytest.approx(318, abs=0.5)


def test_page10_mixed_theta_by_hand():
    """Physics exercise 1: the enthalpy-conserving mixed theta over n layers.

    Eight is still unstable against level 8; nine would already be stable;
    the scheme mixes ten.
    """
    from sympl import get_constant

    with _unyt_backend_restored():
        _, state = _page10_column()
        kappa = (float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
                 / float(get_constant(
                     "heat_capacity_of_dry_air_at_constant_pressure",
                     "J/kg/degK")))
        p = state["air_pressure"].to_units("Pa").values[:, 0, 0]
        p_int = state["air_pressure_on_interface_levels"].to_units(
            "Pa").values[:, 0, 0]
        T = state["air_temperature"].values[:, 0, 0]
        theta = _page10_theta(state)
    dp = p_int[:-1] - p_int[1:]

    def mixed(n):
        return ((T[:n] * dp[:n]).sum()
                / ((p[:n] / 1.0e5) ** kappa * dp[:n]).sum())

    assert mixed(8) == pytest.approx(295.06, abs=0.005)
    assert theta[8] == pytest.approx(293.33, abs=0.005)
    assert mixed(9) == pytest.approx(294.75, abs=0.005)
    assert theta[9] == pytest.approx(295.06, abs=0.005)
    assert mixed(10) == pytest.approx(294.80, abs=0.005)


def test_page10_the_extent_is_decided_by_a_plain_theta_average():
    """Why ten layers, not nine: the page's account of the scheme.

    Replays ``_dry_adj_kernel_np`` top-down. Every mix but the last stops at
    level 8; at k = 0, levels 1-8 already sit at 294.65 K, and the thin, warm
    level 0 lifts the *unweighted* mean over 0-9 above theta at level 9, while
    the mass-weighted mixed value is 294.80 K.
    """
    from sympl import get_constant

    with _unyt_backend_restored():
        adjustment, state = _page10_column()
        kappa = (float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
                 / float(get_constant(
                     "heat_capacity_of_dry_air_at_constant_pressure",
                     "J/kg/degK")))
        p = state["air_pressure"].to_units("Pa").values[:, 0, 0]
        p_int = state["air_pressure_on_interface_levels"].to_units(
            "Pa").values[:, 0, 0]
        T = state["air_temperature"].values[:, 0, 0].copy()
        _, new_state = adjustment(state, climt.UnytTimeDelta(**PAGE10_DT))
        scheme = new_state["air_temperature"].values[:, 0, 0]
    dp = p_int[:-1] - p_int[1:]
    exner = (p / 1.0e5) ** kappa

    mixes = {}
    for k in range(len(T) - 1, -1, -1):
        theta = T / exner
        top = -1
        for m in range(k + 1, len(T)):
            if theta[k:m + 1].mean() > theta[m]:
                top = m
        if top == -1:
            continue
        if k == 0:
            assert theta[1:9] == pytest.approx(294.65, abs=0.005)
            assert theta[0] == pytest.approx(298.78, abs=0.005)
            assert theta[:10].mean() == pytest.approx(295.11, abs=0.005)
            assert theta[9] == pytest.approx(295.06, abs=0.005)
            assert theta[:10].mean() > theta[9]
        layer = slice(k, top + 1)
        mixed = (T[layer] * dp[layer]).sum() / (exner[layer] * dp[layer]).sum()
        T[layer] = mixed * exner[layer]
        mixes[k] = (top, mixed)

    assert {k: top for k, (top, _) in mixes.items()} == {
        6: 7, 5: 8, 4: 8, 3: 8, 2: 8, 1: 8, 0: 9}
    assert mixes[0][1] == pytest.approx(294.80, abs=0.005)
    np.testing.assert_allclose(T, scheme, rtol=1e-12)
    assert dp[0] / 100 == pytest.approx(5.5, abs=0.05)
    assert dp[9] / 100 == pytest.approx(48.7, abs=0.05)


def test_page10_moisture_makes_a_moist_theta_uniform(budgets):
    """Code exercise 2: with uniform humidity, enthalpy is still conserved but
    the dry theta is no longer uniform."""
    with _unyt_backend_restored():
        _, dry, _ = _page10_adjust()
        results = {}
        for q in (0.01, 0.02):
            before, after, adjusted = _page10_adjust(specific_humidity=q)
            theta = _page10_theta(after)
            results[q] = (theta[adjusted].max() - theta[adjusted].min(),
                          budgets.column_enthalpy(before),
                          budgets.column_enthalpy(after),
                          after["air_temperature"].values[0, 0, 0])
        dry_T0 = dry["air_temperature"].values[0, 0, 0]
    assert results[0.01][0] == pytest.approx(0.06, abs=0.005)
    assert results[0.02][0] == pytest.approx(0.12, abs=0.005)
    for spread, h_before, h_after, _ in results.values():
        assert h_after == pytest.approx(h_before, rel=1e-14)
    assert dry_T0 - results[0.01][3] == pytest.approx(0.03, abs=0.005)


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

    # Hand-computed: condensed specific humidity times layer mass dp/g. Read
    # the constant the component reads, so a constants change cannot fail this
    # test for the wrong reason.
    g = float(sympl.get_constant("gravitational_acceleration", "m/s^2"))
    dp = np.asarray(p_int[:-1, ...] - p_int[1:, ...])
    condensed = np.asarray(q_before - outputs["specific_humidity"].values)
    expected_kg_per_m2 = float(np.sum(condensed * dp / g))
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


# --- The browser constraint, checked ----------------------------------------

_NO_NUMBA_TOUR = '''
import sys
sys.path.append({tour!r})

import numpy as np
import sympl
import climt

import budgets
import states
import stepping

sympl.set_backend(climt.UnytBackend())

longwave = climt.CorkLongwaveRadiation(
    optics="correlated_k", table="tour_gray_lw", diffusivity_factor=2.0)
surface = climt.SlabSurface()
boundary_layer = climt.SimpleBoundaryLayer(surface_fluxes="bulk")
adjustment = climt.DryConvectiveAdjustment()
condensation = climt.GridScaleCondensation()


def build():
    return [climt.CorkLongwaveRadiation(optics="correlated_k",
                                        table="tour_gray_lw",
                                        diffusivity_factor=2.0),
            climt.SlabSurface(),
            climt.SimpleBoundaryLayer(surface_fluxes="bulk"),
            climt.DryConvectiveAdjustment(),
            climt.GridScaleCondensation()]


grid = climt.get_grid(nx=1, ny=1, nz=28)
state = climt.get_default_state(
    [longwave, surface, boundary_layer, adjustment, condensation],
    grid_state=grid)
state["ocean_mixed_layer_thickness"].values[:] = 2.0
state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
state["downwelling_shortwave_flux_in_air"].values[0, ...] = {solar!r}
state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
state["air_temperature"].values[:] = np.linspace(
    288.0, 200.0, 28).reshape(28, 1, 1)
state["specific_humidity"].values[:] = 3e-3

stepping.integrate([longwave, surface],
                   [stepping.SurfaceHumidity(1.0), boundary_layer, adjustment,
                    condensation],
                   state, climt.UnytTimeDelta(hours=1), 3)

summary = budgets.summary(state)
path = states.save({path!r}, state, dict(table="tour_gray_lw", nz=28))
reloaded, provenance = states.load(
    path, build(), grid_state=climt.get_grid(nx=1, ny=1, nz=28))
assert provenance["nz"] == 28
np.testing.assert_allclose(reloaded["air_temperature"].values,
                           state["air_temperature"].values)
print(summary["toa_imbalance"], summary["surface_temperature"])
'''


def test_the_tour_helpers_run_with_the_jit_disabled(tmp_path):
    """The tranche's own code, on the path Pyodide will take: no numba.

    Every page runs this code in a browser where numba does not exist, so the
    kernels execute as plain Python and nothing strips units off the timestep
    on the way in. numba reads NUMBA_DISABLE_JIT at import time, so this has
    to be a fresh interpreter. It steps a real column through
    ``stepping.integrate`` with all five browser-safe component kinds and
    page 9's ``SurfaceHumidity``, reads ``budgets.summary`` off the result,
    and round-trips it through
    ``states.save``/``states.load``.
    """
    import os
    import subprocess

    script = _NO_NUMBA_TOUR.format(
        tour=str(TOUR), solar=SOLAR, path=str(tmp_path / "nonumba.npz"))
    result = subprocess.run(
        [sys.executable, "-c", script],
        env={**os.environ, "NUMBA_DISABLE_JIT": "1"},
        capture_output=True, text=True)

    assert result.returncode == 0, (
        "the _tour helpers failed with the JIT disabled -- they will fail the "
        "same way in Pyodide, which has no numba at all:\n" + result.stderr)
    last_line = result.stdout.strip().splitlines()[-1]
    imbalance, surface_temperature = (float(v) for v in last_line.split())
    assert np.isfinite(imbalance)
    assert 150.0 < surface_temperature < 400.0


def test_saved_at_is_naive_utc_and_not_a_deprecated_call(states, tmp_path):
    """`saved_at` stays "YYYY-MM-DDTHH:MM:SS", with no deprecation warning.

    ``datetime.utcnow()`` is deprecated from Python 3.12, which CI builds, so
    ``save`` asks for UTC explicitly -- but the aware datetime it gets back
    would print a "+00:00" offset that the shipped equilibria and
    ``describe`` do not carry. The tz is dropped again, and this pins that.
    """
    import datetime
    import warnings

    state = climt.get_default_state(_page11_components(),
                                    grid_state=get_grid(nx=1, ny=1, nz=28))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        path = states.save(str(tmp_path / "stamp.npz"), state, {})

    _, provenance = states.load(path, _page11_components(),
                                grid_state=get_grid(nx=1, ny=1, nz=28))
    stamp = provenance["saved_at"]

    parsed = datetime.datetime.fromisoformat(stamp)
    assert parsed.tzinfo is None, stamp
    assert stamp == parsed.isoformat(timespec="seconds")
    now = datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None)
    assert abs((now - parsed).total_seconds()) < 600, (
        "saved_at should be UTC, not local time: " + stamp)


# --------------------------------------------------- the shipped equilibria
#
# These states are not regenerated by build_experiments.py's dependency
# hashing. A content hash over cork/**/*.py cannot tell a no-op from a physics
# change, and tranche 1 paid for that twice. The guard is instead: load the
# shipped state, step it a little under the page's own components, and check
# it does not move. Seconds to run, and it fails exactly when the physics did.

DATA = REPO_ROOT / "docs/modelling-tour/_data"


def _generator():
    """The generator script, imported for its component-list definitions."""
    spec = importlib.util.spec_from_file_location(
        "generate_tour_equilibria",
        REPO_ROOT / "scripts/generate_tour_equilibria.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("asset, build", [
    ("rce_dry_equilibrium.npz", "dry_components"),
    ("rce_moist_equilibrium.npz", "moist_components"),
    ("rce_moist_2xco2_equilibrium.npz", "moist_components"),
])
def test_shipped_equilibrium_is_still_an_equilibrium(states, asset, build):
    """The residual test. Ten steps must not move the shipped state.

    If this fails, the physics changed under a state that was computed before
    the change. Do not relax the thresholds: re-run
    `scripts/generate_tour_equilibria.py` (its defaults are what shipped) and
    commit the new file, then re-check every number the page quotes.
    """
    generator = _generator()
    components = getattr(generator, build)()
    tendencies, steppers = generator.split(components)

    state, provenance = states.load(
        str(DATA / asset), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))
    assert state is not None, f"{asset} did not load"

    # The wind relaxation is rebuilt on the loaded state, exactly as the pages
    # do it and as the generator did. Omitting it would let the drag start
    # spinning the column down, and the "shipped state does not move" check
    # would be measuring that instead of the physics it is guarding.
    stepping_module = _load("stepping")
    tendencies = tendencies + [stepping_module.wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]

    before = float(state["surface_temperature"].values.ravel()[0])
    timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
    stepping_module.integrate(tendencies, steppers, state, timestep, 10)
    after = float(state["surface_temperature"].values.ravel()[0])

    budgets_module = _load("budgets")
    imbalance = budgets_module.toa_imbalance(state)

    # The moist column never reaches TOA = 0: DryConvectiveAdjustment
    # conserves a moist-cp enthalpy the rest of the stack does not count
    # (page 12), so it settles with a steady residual (about
    # -1.1 W/m^2 over its saturated surface). Its file records the residual
    # it settled at, as the mean over the gate's last window, and the check
    # is against that. The dry files carry no such key: they settle at 0.
    settled = provenance.get("window_mean_toa_w_m2", 0.0)
    assert abs(imbalance - settled) < 1.0, (
        f"{asset}: TOA imbalance {imbalance:+.3f} W/m^2 after 10 steps, "
        f"against the {settled:+.3f} it settled at — the shipped equilibrium "
        "is stale. Re-run scripts/generate_tour_equilibria.py.")
    assert abs(after - before) < 0.05, (
        f"{asset}: surface temperature drifted {after - before:+.4f} K in 10 "
        "steps — the shipped equilibrium is stale. Re-run "
        "scripts/generate_tour_equilibria.py.")


@pytest.mark.parametrize("asset", ["rce_dry_equilibrium.npz",
                                   "rce_moist_equilibrium.npz",
                                   "rce_moist_2xco2_equilibrium.npz"])
def test_shipped_equilibrium_provenance_is_complete(states, asset):
    """Every field `states.describe` prints must actually be there.

    A page prints this block above its first figure. A `None` in it is a page
    telling a reader nothing while looking like it told them something.
    """
    generator = _generator()
    build = ("dry_components" if "dry" in asset else "moist_components")
    state, provenance = states.load(
        str(DATA / asset), getattr(generator, build)(),
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    for field in ("climt_version", "saved_at", "table", "nz", "dt_hours",
                  "n_steps", "slab_depth_m", "solar", "co2_ppm", "components",
                  "wind_m_s", "wind_timescale_hours", "roughness_length_m",
                  "toa_imbalance", "surface_imbalance"):
        assert provenance.get(field) is not None, f"{asset} lacks {field}"
    assert "None" not in states.describe(provenance)


@pytest.mark.parametrize("asset", ["rce_moist_equilibrium.npz",
                                   "rce_moist_2xco2_equilibrium.npz"])
def test_shipped_moist_states_have_a_saturated_surface(states, soundings,
                                                       asset):
    """The moist states were made over a surface saturated at its own
    temperature, the way page 9 teaches -- not over a fixed specific humidity.

    Until 2026-09-27 they were made with surface_specific_humidity fixed at
    0.015 kg/kg, which at the ~280 K they settled at was 246 % of saturation
    (Bowen ratio 0.02). ``SurfaceHumidity`` in ``moist_components()`` is what
    prevents that now; this pins that the shipped files were made with it.
    """
    generator = _generator()
    state, provenance = states.load(
        str(DATA / asset), generator.moist_components(),
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))

    assert provenance.get("surface_relative_humidity") == 1.0, (
        f"{asset} was not made with stepping.SurfaceHumidity(1.0)")
    assert "SurfaceHumidity" in provenance["components"]
    surface_temperature = float(state["surface_temperature"].values.ravel()[0])
    surface_pressure = float(state["surface_air_pressure"].values.ravel()[0])
    q_sat = float(soundings.saturation_specific_humidity(
        np.array(surface_temperature), np.array(surface_pressure)))
    q_surface = float(state["surface_specific_humidity"].values.ravel()[0])
    # Written from the surface temperature one slab update earlier, so equal
    # to within that step's change, not exactly.
    assert q_surface / q_sat == pytest.approx(1.0, abs=0.02), (
        f"{asset}: surface at {100 * q_surface / q_sat:.0f} % RH")


def test_shipped_2xco2_state_is_the_shipped_moist_state_perturbed(states):
    """Page 12 quotes the warming from this file, so it must be *this* base.

    The 2xCO2 state is made from rce_moist_equilibrium.npz by
    ``generate_tour_equilibria.py --moist-2xco2``. If the base is regenerated
    and the 2xCO2 state is not, the warming page 12 quotes is measured from an
    equilibrium the page no longer loads. The base's own ``saved_at`` stamp,
    recorded in the perturbed file, is what ties the two together.
    """
    generator = _generator()
    grid = get_grid(nx=1, ny=1, nz=generator.NZ)
    base_state, base = states.load(str(DATA / "rce_moist_equilibrium.npz"),
                                   generator.moist_components(),
                                   grid_state=grid)
    state, doubled = states.load(
        str(DATA / "rce_moist_2xco2_equilibrium.npz"),
        generator.moist_components(), grid_state=grid)

    assert doubled["perturbed_from"] == "rce_moist_equilibrium.npz"
    assert doubled["perturbed_from_saved_at"] == base["saved_at"], (
        "rce_moist_2xco2_equilibrium.npz was made from a different moist "
        "equilibrium than the one shipped. Re-run "
        "scripts/generate_tour_equilibria.py --moist-2xco2.")
    assert doubled["co2_ppm"] == generator.CO2_DOUBLING * base["co2_ppm"]
    for key in ("table", "nz", "dt_hours", "slab_depth_m", "solar",
                "wind_m_s", "wind_timescale_hours", "roughness_length_m",
                "components"):
        assert doubled[key] == base[key], key

    co2 = state["mole_fraction_of_carbon_dioxide_in_air"].values
    np.testing.assert_allclose(co2, doubled["co2_ppm"] * 1e-6)
    before = float(base_state["surface_temperature"].values.ravel()[0])
    after = float(state["surface_temperature"].values.ravel()[0])
    assert doubled["base_surface_temperature_k"] == pytest.approx(before)
    assert doubled["warming_k"] == pytest.approx(after - before)
    assert after > before, "doubling CO2 cooled the moist column"


# ------------------------------------------------------------------ page 11
#
# Page 11 loads rce_dry_equilibrium.npz and perturbs it. Its component list is
# the generator's dry_components(), and the wind relaxation is rebuilt on the
# loaded state with initialise=False, as the page does it. Numbers are those
# the page's cells print; the ones it quotes from longer runs come from
# scripts/experiments/tour_page11_measurements.py.

PAGE11_PAGE = REPO_ROOT / "docs/modelling-tour/11-dry-rce.qmd"
PAGE11_2XCO2_STEPS = 1000     # 500 days at 12 h; ~2 min in a browser
PAGE11_MONTH_STEPS = 60       # 30 days at 12 h: the averaging window
DRY_ADIABAT_K_PER_KM = 9.76   # g / c_p, to the precision the page quotes
TRANCHE_1_LAPSE_K_PER_KM = 6.5


def _lapse_rate_profile(state):
    """-dT/dz in K/km between mid levels, from the hypsometric thickness."""
    from sympl import get_constant

    Rd = float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
    g = float(get_constant("gravitational_acceleration", "m/s^2"))
    T = state["air_temperature"].values[:, 0, 0]
    p = state["air_pressure"].values[:, 0, 0]
    T_mean = 0.5 * (T[:-1] + T[1:])
    dz = (Rd * T_mean / g) * np.log(p[:-1] / p[1:])
    return -np.diff(T) / dz * 1000.0


def _page11_equilibrium():
    """(tendencies, steppers, state, provenance), as the page builds them.

    Call inside ``_unyt_backend_restored``.
    """
    generator = _generator()
    components = generator.dry_components()
    tendencies, steppers = generator.split(components)
    state, provenance = _load("states").load(
        str(DATA / "rce_dry_equilibrium.npz"), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))
    tendencies = tendencies + [_load("stepping").wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]
    return tendencies, steppers, state, provenance


def _page11_cells(upto, monkeypatch):
    """Exec the page's own cells 0..upto, as a reader would; return the
    namespace. The cells find _tour/ and _data/ relative to the page."""
    import re

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cells = re.findall(r"```\{pyodide\}\n(.*?)```", PAGE11_PAGE.read_text(),
                       flags=re.DOTALL)
    monkeypatch.chdir(PAGE11_PAGE.parent)
    # Cell 0 does sys.path.insert(0, "_tour"); keep that out of the session.
    monkeypatch.setattr(sys, "path", list(sys.path))
    namespace = {"__name__": "__main__"}
    try:
        for index in range(upto + 1):
            exec(compile(cells[index], f"11-dry-rce.qmd[cell {index}]",
                         "exec"), namespace)
    finally:
        plt.close("all")
    return namespace


def _page11_theta(state):
    from sympl import get_constant

    kappa = (float(get_constant("gas_constant_of_dry_air", "J/kg/degK"))
             / float(get_constant(
                 "heat_capacity_of_dry_air_at_constant_pressure",
                 "J/kg/degK")))
    p = state["air_pressure"].to_units("Pa").values[:, 0, 0]
    return state["air_temperature"].values[:, 0, 0] * (1.0e5 / p) ** kappa


def test_page11_the_shipped_column_is_dry():
    """The page's second section, and the plan's first suspect if the lapse
    rate ever comes out near 6.5: moisture would pull it toward moist
    adiabatic, page 12's result arriving a page early."""
    with _unyt_backend_restored():
        _, _, state, provenance = _page11_equilibrium()
    assert np.all(state["specific_humidity"].values == 0.0)
    assert np.all(state["surface_specific_humidity"].values == 0.0)
    assert provenance["co2_ppm"] == 330.0
    assert float(state["surface_temperature"].values.ravel()[0]) == \
        pytest.approx(266.62, abs=0.005)


@pytest.mark.slow     # one step, but ~11 s run alone: numba compiling the stack
def test_page11_convecting_layer_is_dry_adiabatic():
    """Page 11's reveal, and the number tranche 1 assumed away.

    The convecting layer is the levels dry adjustment actually touches. The
    shipped state was saved straight after an adjustment call, so calling the
    scheme on it again moves nothing by more than rounding (~1e-13 K). So:
    take one step of everything *but* the adjustment -- which re-creates the
    instability radiation and the boundary layer make each step -- and then
    see which levels the adjustment moves.
    """
    with _unyt_backend_restored():
        tendencies, steppers, state, provenance = _page11_equilibrium()
        adjustment = steppers[-1]
        assert isinstance(adjustment, climt.DryConvectiveAdjustment)

        _, again = adjustment(state, climt.UnytTimeDelta(hours=12))
        assert np.abs(again["air_temperature"].values
                      - state["air_temperature"].values).max() < 1e-8

        timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
        _load("stepping").integrate(tendencies, steppers[:-1], state,
                                    timestep, 1)
        before = state["air_temperature"].values[:, 0, 0].copy()
        _, adjusted = adjustment(state, timestep)
        state.update(adjusted)
        touched = np.abs(state["air_temperature"].values[:, 0, 0]
                         - before) > 1e-8
        lapse = _lapse_rate_profile(state)

    convecting = touched[:-1] & touched[1:]
    assert convecting.sum() >= 3, (
        "dry adjustment touched fewer than four levels after one step from "
        "the shipped equilibrium -- there is no convecting layer to make a "
        "claim about")
    mean_lapse = float(np.mean(lapse[convecting]))
    assert mean_lapse == pytest.approx(DRY_ADIABAT_K_PER_KM, rel=0.01), (
        f"convecting-layer lapse rate {mean_lapse:.2f} K/km is not within "
        f"1% of dry adiabatic ({DRY_ADIABAT_K_PER_KM})")
    assert mean_lapse > 8.0, (
        f"lapse rate {mean_lapse:.2f} K/km -- page 11's whole point is that "
        "it is NOT the 6.5 K/km tranche 1 prescribed")


def test_page11_tropopause_emerges_rather_than_being_prescribed():
    """Above the convecting layer the profile is stable, and nothing capped it.

    Tranche 1 prescribed a 200 K isothermal stratosphere. Here the top of
    convection is wherever the scheme stops, and the air above it is whatever
    radiation makes it: the lapse rate falls off gradually and the
    temperature keeps falling to the lid, with no minimum.
    """
    with _unyt_backend_restored():
        _, _, state, _ = _page11_equilibrium()
        T = state["air_temperature"].values[:, 0, 0]
        p_hPa = state["air_pressure"].values[:, 0, 0] / 100.0
        lapse = _lapse_rate_profile(state)

    assert not np.any(np.isclose(T, 200.0, atol=0.05)), (
        "a level sitting at exactly 200.0 K suggests a prescribed cap leaked "
        "in from tranche 1's soundings")
    assert lapse[-1] < lapse[3], (
        "the upper column should be less steeply lapsing than the convecting "
        "layer below it")
    assert np.all(np.diff(T[3:]) < 0), "no temperature minimum above 970 hPa"
    assert T[-1] == pytest.approx(111.13, abs=0.005)
    upper = np.sqrt(p_hPa[:-1] * p_hPa[1:]) < 340.0
    assert np.all(lapse[upper] < 9.0)
    assert lapse[-1] == pytest.approx(1.3, abs=0.05)


def test_page11_headline_cells_print_what_the_page_says(monkeypatch, capsys):
    """Cells 0-2, run as the page runs them, print the numbers its prose
    quotes: the convecting layer, its lapse rate, the boundary layer's
    adiabatic bottom, and the comparison with page 7's column."""
    with _unyt_backend_restored():
        namespace = _page11_cells(2, monkeypatch)
    out = capsys.readouterr().out

    for text in ("TOA -0.013, surface -0.000 W/m^2",
                 "convecting layer   1010 to 370 hPa, 17 levels",
                 "lapse rate there   9.76 K/km",
                 "dry adiabat g/Cp   9.76 K/km",
                 "jump 2.88 K",
                 "9.8, 9.8, 9.8 K/km",
                 "top of the model   111.13 K",
                 "surface           -1.16 K",
                 "air, largest      +23.70 K at 695 hPa",
                 "air, top level    +9.93 K"):
        assert text in out, f"{text!r} not printed:\n{out}"
    assert namespace["top"] == 16 and namespace["bottom"] == 0


@pytest.mark.slow
def test_page11_page7_profile_is_page7s_column_at_equilibrium(monkeypatch):
    """The dashed profile the headline cell writes in, recomputed.

    Page 7's 14-band column -- LW and slab only, from climt's default state
    -- stepped 10 000 steps at 12 h. The cell rounds to 0.1 K. ~30 s native.
    """
    with _unyt_backend_restored():
        namespace = _page11_cells(1, monkeypatch)
        longwave = climt.CorkLongwaveRadiation(
            optics="correlated_k", table="earth_low_res_lw",
            diffusivity_factor=1.66)
        slab = climt.SlabSurface()
        state = climt.get_default_state(
            [longwave, slab], grid_state=get_grid(nx=1, ny=1, nz=28))
        state["ocean_mixed_layer_thickness"].values[:] = 2.0
        state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
        state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
        state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
        stepping_module = _load("stepping")
        stepping_module.integrate([longwave, slab], [], state,
                                  climt.UnytTimeDelta(hours=12), 10000)
        imbalance = _load("budgets").toa_imbalance(state)

    assert abs(imbalance) < 0.005
    np.testing.assert_allclose(namespace["PAGE7_T"],
                               state["air_temperature"].values[:, 0, 0],
                               atol=0.06)
    assert namespace["PAGE7_SURFACE"] == pytest.approx(
        float(state["surface_temperature"].values.ravel()[0]), abs=0.005)
    assert _lapse_rate_profile(state)[0] == pytest.approx(81.5, abs=0.05)


def test_page11_the_surface_shines_through_the_window(monkeypatch):
    """Physics exercise 1 and the "Against page 7" prose: one longwave call
    on hybrids of the two states.

    The split holds the air -- and so its temperature-dependent absorption --
    fixed and chills the surface. Chilling the air instead makes the column
    more transparent, and "the surface alone" then exceeds the whole OLR.
    """
    with _unyt_backend_restored():
        namespace = _page11_cells(1, monkeypatch)
        state, longwave = namespace["state"], namespace["components"][0]
        Ts = float(state["surface_temperature"].values.ravel()[0])

        def olr(column, air=None, surface=None):
            column = copy.deepcopy(column)
            if air is not None:
                column["air_temperature"].values[:, 0, 0] = air
            if surface is not None:
                column["surface_temperature"].values[:] = surface
            return float(longwave(column)[1][
                "upwelling_longwave_flux_in_air"].values[-1, 0, 0])

        base = olr(state)
        page7_air = olr(state, air=namespace["PAGE7_T"])
        air_only = olr(state, surface=1.0)
        direct = olr(state, surface=Ts + 1.0) - base
        wrong_split = olr(state, air=1.0)
    sigma = 5.670374419e-8
    assert base == pytest.approx(240.01, abs=0.005)
    assert page7_air == pytest.approx(236.02, abs=0.005)
    assert base - page7_air == pytest.approx(3.99, abs=0.005)
    assert air_only == pytest.approx(26.95, abs=0.005)
    assert base - air_only == pytest.approx(213.06, abs=0.005)
    assert (base - air_only) / (sigma * Ts ** 4) == pytest.approx(0.74,
                                                                  abs=0.005)
    assert direct == pytest.approx(3.38, abs=0.005)
    assert direct / (4 * sigma * Ts ** 3) == pytest.approx(0.79, abs=0.005)
    assert 1.16 * direct == pytest.approx(3.92, abs=0.005)
    assert wrong_split == pytest.approx(244.61, abs=0.005)
    assert wrong_split > base, "the artefact the exercise warns about"


def _page11_unperturbed_mean(tendencies, steppers, state, provenance):
    """The knob cell's baseline: the 30-day mean surface temperature of a
    copy of ``state`` stepped 60 steps unperturbed. ``state`` is untouched."""
    stepping_module = _load("stepping")
    timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
    record = stepping_module.Recorder(
        timestep,
        surface=lambda s: float(s["surface_temperature"].values.ravel()[0]))
    stepping_module.integrate(tendencies, steppers, copy.deepcopy(state),
                              timestep, PAGE11_MONTH_STEPS, after_step=record)
    return record.mean("surface", days=30)


@pytest.mark.slow
def test_page11_co2_doubling_warms_the_surface_by_a_measured_amount():
    """Page 11's knob: a *measured* dry climate sensitivity.

    Tranche 1's page 5 could only compute the forcing. This integrates to the
    new equilibrium, as the page's knob cell does -- one 1000-step
    ``integrate`` -- and reads the warming off: the 30-day mean surface
    temperature minus the unperturbed column's 30-day mean. Measured +1.16 K.
    """
    with _unyt_backend_restored():
        tendencies, steppers, state, provenance = _page11_equilibrium()
        stepping_module = _load("stepping")
        budgets_module = _load("budgets")
        longwave = tendencies[0]

        def olr(column):
            return float(longwave(column)[1][
                "upwelling_longwave_flux_in_air"].values[-1, 0, 0])

        before = _page11_unperturbed_mean(tendencies, steppers, state,
                                          provenance)
        perturbed = copy.deepcopy(state)
        perturbed["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= 2.0
        forcing = olr(state) - olr(perturbed)

        timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
        record = stepping_module.Recorder(
            timestep,
            surface=lambda s: float(
                s["surface_temperature"].values.ravel()[0]),
            toa=budgets_module.toa_imbalance,
            surface_imbalance=budgets_module.surface_imbalance)
        stepping_module.integrate(tendencies, steppers, perturbed, timestep,
                                  PAGE11_2XCO2_STEPS, after_step=record)
        imbalance = budgets_module.toa_imbalance(perturbed)

    final = record.mean("surface", days=30)
    warming = final - before
    assert 0.69 < warming < 1.61, (
        f"2xCO2 dry surface warming {warming:+.2f} K is outside +-40% of the "
        "measured +1.16 K -- check that the run reached equilibrium")
    assert warming == pytest.approx(1.16, abs=0.005), "the page quotes +1.16 K"
    assert forcing == pytest.approx(4.09, abs=0.005)
    assert forcing / warming == pytest.approx(3.52, abs=0.005)
    assert forcing / 3.38 == pytest.approx(1.21, abs=0.005)
    assert abs(imbalance) < 0.5, (
        f"TOA imbalance {imbalance:+.3f} W/m^2 -- the perturbed run has not "
        f"equilibrated in {PAGE11_2XCO2_STEPS} steps, so the warming is a "
        "lower bound, not a sensitivity")
    assert abs(record.mean("toa", days=30)) < 0.05
    assert abs(record.mean("surface_imbalance", days=30)) < 0.05
    off = np.where(np.abs(record["surface"] - final) > 0.1)[0].max() + 1
    assert record["days"][off] == pytest.approx(120, abs=0.5)
    assert np.all(np.abs(record["toa"][299::100]) < 0.3)


@pytest.mark.slow
def test_page11_halving_co2_is_roughly_symmetric():
    """Code exercise 1: -1.12 K for halving against +1.16 K for doubling."""
    with _unyt_backend_restored():
        tendencies, steppers, state, provenance = _page11_equilibrium()
        before = _page11_unperturbed_mean(tendencies, steppers, state,
                                          provenance)
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] *= 0.5
        timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])
        record = _load("stepping").Recorder(
            timestep,
            surface=lambda s: float(
                s["surface_temperature"].values.ravel()[0]))
        _load("stepping").integrate(tendencies, steppers, state, timestep,
                                    PAGE11_2XCO2_STEPS, after_step=record)
    assert record.mean("surface", days=30) - before == pytest.approx(
        -1.12, abs=0.01)


@pytest.mark.slow
def test_page11_adjustment_last_leaves_every_step_stable():
    """The ordering callout. In the page's order no step ends with theta
    falling anywhere; with the adjustment before the boundary layer, every
    step does, somewhere between 1010 and 424 hPa."""
    with _unyt_backend_restored():
        stepping_module = _load("stepping")
        results = {}
        for name, reorder in (("shipped", lambda s: s),
                              ("swapped", lambda s: s[::-1])):
            tendencies, steppers, state, provenance = _page11_equilibrium()
            unstable = []
            stepping_module.integrate(
                tendencies, reorder(steppers), state,
                climt.UnytTimeDelta(hours=provenance["dt_hours"]),
                PAGE11_MONTH_STEPS,
                after_step=lambda column: unstable.append(
                    np.where(np.diff(_page11_theta(column)) < -1e-6)[0]))
            results[name] = unstable
            p_hPa = state["air_pressure"].values[:, 0, 0] / 100.0
    assert not any(len(levels) for levels in results["shipped"])
    assert all(len(levels) for levels in results["swapped"])
    where = p_hPa[np.concatenate(results["swapped"])]
    assert np.all((where > 420) & (where < 1015)), where


def test_page11_saved_state_round_trips(monkeypatch, tmp_path):
    """Cell 4's craft: save a state with its provenance and load it back.

    Uses the loaded equilibrium itself rather than cell 3's 1000-step run;
    the round trip is what is being checked.
    """
    with _unyt_backend_restored():
        states_module = _load("states")
        _, _, state, provenance = _page11_equilibrium()
        path = states_module.save(
            str(tmp_path / "rce_dry_2xco2.npz"), state,
            dict(provenance, co2_ppm=660.0, n_steps=1000,
                 perturbed_from="rce_dry_equilibrium.npz"))
        reloaded, meta = states_module.load(
            path, _generator().dry_components(),
            grid_state=get_grid(nx=1, ny=1, nz=28))
    assert meta["perturbed_from"] == "rce_dry_equilibrium.npz"
    assert "660.0 ppm" in states_module.describe(meta)
    np.testing.assert_array_equal(reloaded["air_temperature"].values,
                                  state["air_temperature"].values)


# ------------------------------------------------------------------ page 12
#
# Page 12 loads rce_moist_equilibrium.npz, its doubled-CO2 twin, and page 11's
# dry state for comparison. Its component list is the generator's
# moist_components(), built in the open in cell 0. Numbers the cells print are
# checked by exec'ing the cells; numbers the page quotes from longer runs come
# from scripts/experiments/tour_page12_measurements.py, and the slow tests at
# the end re-measure the ones that can be afforded.
#
# The page's reveal is not the one its spec expected. The spec says latent
# heating relaxes the lapse rate to near 6.5 K/km. In this column it relaxes it
# by ~0.3 K/km: the grid-scale condensation does ~99 % of the raining, and the
# dry adjustment keeps the lowest ~180 hPa on the dry adiabat. The tests pin
# what the column does, and what the page says about it.

PAGE12_PAGE = REPO_ROOT / "docs/modelling-tour/12-moist-rce.qmd"
PAGE12_BAND_HPA = (500.0, 792.2)    # above the adjusted layer, to 500 hPa


def _page_cells(page, indices, monkeypatch):
    """Exec a page's cells ``indices``, in order, as a reader would run them;
    return the namespace. The cells find _tour/ and _data/ relative to the
    page."""
    import re

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cells = re.findall(r"```\{pyodide\}\n(.*?)```", page.read_text(),
                       flags=re.DOTALL)
    monkeypatch.chdir(page.parent)
    monkeypatch.setattr(sys, "path", list(sys.path))
    namespace = {"__name__": "__main__"}
    try:
        for index in indices:
            exec(compile(cells[index], f"{page.name}[cell {index}]", "exec"),
                 namespace)
    finally:
        plt.close("all")
    return namespace


def _page12_equilibrium(filename="rce_moist_equilibrium.npz"):
    """(tendencies, steppers, state, provenance), as the page builds them.

    Call inside ``_unyt_backend_restored``.
    """
    generator = _generator()
    components = generator.moist_components()
    tendencies, steppers = generator.split(components)
    state, provenance = _load("states").load(
        str(DATA / filename), components,
        grid_state=get_grid(nx=1, ny=1, nz=generator.NZ))
    tendencies = tendencies + [_load("stepping").wind_relaxation(
        state, provenance["wind_m_s"], provenance["wind_timescale_hours"],
        initialise=False)]
    return tendencies, steppers, state, provenance


def _band(p):
    """Layers whose midpoint lies above the adjusted layer, below 500 hPa."""
    p_mid = np.sqrt(p[:-1] * p[1:]) / 100.0
    return (p_mid > PAGE12_BAND_HPA[0]) & (p_mid < PAGE12_BAND_HPA[1])


def test_moist_adiabatic_lapse_rate(soundings):
    """The rate is small in warm air, climbs back to g/cp in cold air, and
    matches the numbers page 12 quotes."""
    with _unyt_backend_restored():
        warm = soundings.moist_adiabatic_lapse_rate(300.0, 1.0e5)
        cold = soundings.moist_adiabatic_lapse_rate(200.0, 2.0e4)
        lowest = soundings.moist_adiabatic_lapse_rate(280.69, 101040.0)
    assert warm < 4.0
    assert cold == pytest.approx(DRY_ADIABAT_K_PER_KM, abs=0.05)
    assert lowest == pytest.approx(5.57, abs=0.005)


@pytest.mark.parametrize("T_start, mean_lapse", [(284.0, 6.52),
                                                 (286.0, 6.19)])
def test_moist_adiabat_averages_what_page12_says(soundings, T_start,
                                                 mean_lapse):
    """'A parcel saturated at 284 K at 1000 hPa averages 6.52 K/km to 500'."""
    with _unyt_backend_restored():
        p = np.geomspace(1.0e5, 5.0e4, 400)
        T = soundings.moist_adiabat(T_start, 1.0e5, p)
    dz = np.sum(287.0 * 0.5 * (T[1:] + T[:-1]) / 9.80665
                * np.log(p[:-1] / p[1:]))
    assert 1e3 * (T[0] - T[-1]) / dz == pytest.approx(mean_lapse, abs=0.005)


def test_moist_adiabat_starts_where_it_starts(soundings):
    with _unyt_backend_restored():
        T = soundings.moist_adiabat(285.0, 9.0e4, np.array([1.0e5, 9.0e4,
                                                            5.0e4]))
    assert np.isnan(T[0])
    assert T[1] == pytest.approx(285.0)
    assert T[2] < T[1]


def test_describe_adds_the_moist_lines_only_for_moist_states(states):
    """Page 12's describe block shows the surface humidity and the settled
    TOA; page 11's is unchanged."""
    with _unyt_backend_restored():
        generator = _generator()
        grid = get_grid(nx=1, ny=1, nz=generator.NZ)
        _, moist = states.load(str(DATA / "rce_moist_equilibrium.npz"),
                               generator.moist_components(), grid_state=grid)
        _, dry = states.load(str(DATA / "rce_dry_equilibrium.npz"),
                             generator.dry_components(), grid_state=grid)
    moist_text, dry_text = states.describe(moist), states.describe(dry)
    assert "295800 steps at dt = 5 min" in moist_text
    assert "surface RH     100 %" in moist_text
    assert "settled TOA    -1.095 W/m^2" in moist_text
    assert "3500 steps at dt = 12.0 h" in dry_text
    for token in ("surface RH", "settled TOA", "perturbed from"):
        assert token not in dry_text


def test_page12_headline_cells_print_what_the_page_says(monkeypatch, capsys):
    """Cells 0-2, run as the page runs them, print the numbers its prose
    quotes."""
    with _unyt_backend_restored():
        namespace = _page_cells(PAGE12_PAGE, range(3), monkeypatch)
    out = capsys.readouterr().out

    for text in ("'CorkLongwaveRadiation', 'SlabSurface', "
                 "'EmanuelConvectionPython', 'UnytRelaxation'",
                 "295800 steps at dt = 5 min",
                 "settled TOA    -1.095 W/m^2",
                 "surface_temperature     285.99",
                 "toa_imbalance            -1.27",
                 "dry surface (p. 11)     266.48",
                 "adjusted layer      968 to 792 hPa, lapse rate 9.76 K/km, "
                 "relative humidity 32 to 83 %",
                 "above it, to 500    lapse rate 9.48 K/km (page 11's dry "
                 "column: 9.75); the moist adiabat there 7.5 to 9.3",
                 "(280.7 K): 5.57 K/km; it passes 6.5 at 876 hPa, where "
                 "T = 270.6 K",
                 "(286.0 K): 6.24 K/km on average to 500 hPa, where it is "
                 "251.3 K and this column 231.4 K",
                 "+19.51 K"):
        assert text in out, f"{text!r} not printed:\n{out}"
    assert (namespace["bottom"], namespace["top"]) == (3, 8)


def test_page12_lapse_rate_lies_between_dry_and_moist_adiabatic(soundings):
    """The spec's claim, as far as it holds.

    Above the adjusted layer the moist column lapses less steeply than dry
    adiabatic and more steeply than its own moist adiabat. It is *not* near
    the moist adiabat -- the spec expected that, and the page explains why it
    does not happen in this stack. If this ever comes out near 6.5, the page's
    'Why this column is not on its moist adiabat' section is wrong.
    """
    with _unyt_backend_restored():
        _, _, state, _ = _page12_equilibrium()
        T = state["air_temperature"].values[:, 0, 0]
        p = state["air_pressure"].values[:, 0, 0]
        lapse = _lapse_rate_profile(state)
        moist = soundings.moist_adiabatic_lapse_rate(
            0.5 * (T[:-1] + T[1:]), np.sqrt(p[:-1] * p[1:]))
    band = _band(p)
    assert band.sum() == 6
    assert np.all(np.abs(lapse[3:8] - DRY_ADIABAT_K_PER_KM) < 0.01), (
        "the adjusted layer should be on the dry adiabat")
    assert moist[band].mean() < lapse[band].mean() < DRY_ADIABAT_K_PER_KM - 0.2
    assert lapse[band].mean() > 8.5, (
        f"{lapse[band].mean():.2f} K/km above the adjusted layer: the column "
        "has moved toward its moist adiabat -- page 12's explanation of why "
        "it does not needs re-measuring")


def test_page12_moist_column_is_warmer_than_the_dry_one(states):
    """The water vapour greenhouse plus its feedback, as a difference."""
    generator = _generator()
    with _unyt_backend_restored():
        grid = get_grid(nx=1, ny=1, nz=generator.NZ)
        dry, _ = states.load(str(DATA / "rce_dry_equilibrium.npz"),
                             generator.dry_components(), grid_state=grid)
        moist, _ = states.load(str(DATA / "rce_moist_equilibrium.npz"),
                               generator.moist_components(), grid_state=grid)
    difference = (float(moist["surface_temperature"].values.ravel()[0])
                  - float(dry["surface_temperature"].values.ravel()[0]))
    assert difference == pytest.approx(19.51, abs=0.005)


def test_page12_knob_cell_takes_the_warming_apart(monkeypatch, capsys):
    """Cell 4: the offline 2xCO2 state, and the radiation-only decomposition
    that says the water vapour feedback is the whole difference from page
    11's +1.15 K."""
    with _unyt_backend_restored():
        _page_cells(PAGE12_PAGE, (0, 1, 4), monkeypatch)
    out = capsys.readouterr().out
    for text in ("CO2            660.0 ppm",
                 "180550 steps at dt = 5 min",
                 "perturbed from rce_moist_equilibrium.npz "
                 "(saved 2026-09-27T02:22:56)",
                 "+2.26 K file to file; +2.22 K between",
                 "forcing, CO2 doubled and nothing else   +4.69 W/m^2",
                 "(+3.95 per K)", "(-1.89 per K)",
                 "warming if the vapour had not changed   +1.19 K"):
        assert text in out, f"{text!r} not printed:\n{out}"


def test_page12_emanuel_inside_adamsbashforth_warns_and_still_works():
    """Pins the warning the page explains -- and that it is the only one.

    If sympl ever stops emitting it, page 12's warning callout describes
    something the reader will not see. The second warning sympl can emit
    here, about being handed a list, is silenced by stepping.py splatting
    its components; the page says nothing about it, so it must stay silent.
    """
    import warnings

    from sympl import AdamsBashforth

    with _unyt_backend_restored():
        convection = climt.EmanuelConvectionPython()
        longwave = climt.CorkLongwaveRadiation(optics="correlated_k",
                                               table="earth_low_res_lw")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            AdamsBashforth(longwave, climt.SlabSurface(), convection)

    messages = " ".join(str(w.message) for w in caught)
    assert "ImplicitTendencyComponent" in messages, (
        "sympl no longer warns about an ImplicitTendencyComponent inside a "
        "TendencyStepper — page 12's warning callout needs rewriting")
    assert "individual Prognostics" not in messages


@pytest.mark.slow
def test_page12_budget_cell_closes_the_moisture_budget(monkeypatch, capsys):
    """Cell 3: two days at 5 min. P balances E, the rain is grid-scale, the
    condensation heats above the adjusted layer, and the dry adjustment adds
    the ~1 W/m^2 the TOA residual is made of."""
    with _unyt_backend_restored():
        namespace = _page_cells(PAGE12_PAGE, range(4), monkeypatch)
    out = capsys.readouterr().out
    for text in ("precipitation    2.62 mm/day  (Emanuel 0.03, grid-scale "
                 "2.59)",
                 "evaporation      2.59 mm/day",
                 "P - E           +0.03 mm/day",
                 "sensible 27.8, latent 74.9 W/m^2: Bowen ratio 0.37",
                 "  745 hPa   2.57",
                 "  370 hPa   0.79",
                 "column enthalpy added: adjustment +0.95, condensation "
                 "+0.00 W/m^2"):
        assert text in out, f"{text!r} not printed:\n{out}"
    P, E = namespace["P"], namespace["E"]
    assert P > 0.0, "an equilibrium moist column must precipitate"
    assert abs(P - E) < 0.05 * E, "the moisture budget does not close"
    assert namespace["P_conv"] < 0.05 * P, "Emanuel is doing the raining"


@pytest.mark.slow
def test_page12_timestep_cell_moves_emanuel_not_the_column(monkeypatch,
                                                           capsys):
    """Cell 5: at 10 min Emanuel's rain collapses; the column hardly moves."""
    with _unyt_backend_restored():
        _page_cells(PAGE12_PAGE, (0, 1, 3, 5), monkeypatch)
    out = capsys.readouterr().out
    for text in ("surface temperature (K)    285.99   286.02",
                 "precipitation (mm/day)       2.62     2.75",
                 "  from Emanuel              0.029    0.002",
                 "evaporation (mm/day)         2.59     2.54",
                 "TOA imbalance (W/m^2)       -1.03    -0.95",
                 "largest difference: -0.26 K at 268 hPa"):
        assert text in out, f"{text!r} not printed:\n{out}"


@pytest.mark.slow
def test_page12_thirty_day_means_the_page_quotes():
    """The 30-day numbers page 12 quotes for the base state
    (``tour_page12_measurements.py month``): P 2.65 and E 2.64 mm/day, Emanuel's
    0.03 of it, the Bowen ratio, and the lapse rates."""
    with _unyt_backend_restored():
        tendencies, steppers, state, provenance = _page12_equilibrium()
        stepping_module = _load("stepping")
        budgets_module = _load("budgets")
        timestep = climt.UnytTimeDelta(hours=provenance["dt_hours"])

        def surface(name):
            return lambda s: float(s[name].values.ravel()[0])

        record = stepping_module.Recorder(
            timestep, profiles=("air_temperature",), profile_every=1,
            precipitation=lambda s: budgets_module.precipitation_rate(
                s, timestep),
            convective=surface("convective_precipitation_rate"),
            evaporation=budgets_module.evaporation_rate,
            sensible=budgets_module.sensible_heat_flux,
            latent=budgets_module.latent_heat_flux)
        stepping_module.integrate(tendencies, steppers, state, timestep,
                                  30 * 288, after_step=record)
        p = state["air_pressure"].values[:, 0, 0]
    T = record["air_temperature"].mean(axis=0)
    dz = (287.0 * 0.5 * (T[:-1] + T[1:]) / 9.80665) * np.log(p[:-1] / p[1:])
    lapse = -np.diff(T) / dz * 1000.0
    assert record.mean("precipitation") == pytest.approx(2.65, abs=0.005)
    assert record.mean("evaporation") == pytest.approx(2.64, abs=0.005)
    assert record.mean("convective") == pytest.approx(0.03, abs=0.005)
    assert (record.mean("sensible") / record.mean("latent")
            == pytest.approx(0.37, abs=0.005))
    assert lapse[_band(p)].mean() == pytest.approx(9.48, abs=0.005)
    assert lapse[p[:-1] > 5.0e4].mean() == pytest.approx(9.07, abs=0.005)
