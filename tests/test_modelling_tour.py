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
