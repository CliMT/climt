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
