# tests/test_tropopause.py
"""Unit tests for the robust tropopause locator (scripts/experiments/tropopause.py).

Profiles are ordered surface -> top (index 0 = surface), matching the RCE
scripts' state arrays. The old constant-theta+10K marker mislocated the
tropopause for runaway-warm columns (the 583 hPa artifact); these tests pin the
replacement's behaviour.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "experiments"))
from tropopause import find_tropopause, potential_temperature  # noqa: E402


def _synthetic_column(p_tropopause_hpa, nz=40):
    """RCE-like column: constant-theta troposphere, then a sharp stratospheric
    theta ramp above the prescribed tropopause pressure. Returns (T, p_hPa)."""
    p = np.linspace(1000.0, 50.0, nz)              # surface -> top, decreasing
    theta_trop = 300.0
    # theta constant below the tropopause, ramping up steeply above it.
    theta = np.where(
        p >= p_tropopause_hpa,
        theta_trop,
        theta_trop + 2.0 * (p_tropopause_hpa - p),  # steep rise into stratosphere
    )
    T = theta * (p / 1000.0) ** (287.0 / 1004.0)
    return T, p


def test_potential_temperature_surface_identity():
    T = np.array([300.0])
    p = np.array([1000.0])
    assert potential_temperature(T, p)[0] == pytest.approx(300.0)


def test_curvature_finds_prescribed_kink():
    T, p = _synthetic_column(p_tropopause_hpa=200.0, nz=40)
    out = find_tropopause(T, p)
    # The kink is at 200 hPa; the located level must be within one grid step.
    grid_step = abs(p[1] - p[0])
    assert abs(out["p_curvature"] - 200.0) <= grid_step


def test_runaway_warm_column_is_not_mislocated_low():
    """A deep, warm, moist-adiabat-ish troposphere whose theta drifts upward
    through the column -- the case where +theta0+10K fired too low (~583 hPa).
    The curvature locator must land high (low pressure), not in the mid-troposphere."""
    nz = 40
    p = np.linspace(1010.0, 40.0, nz)
    # theta drifts +30 K across the troposphere (moist adiabat), then a hard
    # stratospheric cap above 150 hPa.
    theta_trop = 290.0 + 30.0 * (1010.0 - np.clip(p, 150.0, 1010.0)) / (1010.0 - 150.0)
    theta_strat_extra = np.where(p < 150.0, 3.0 * (150.0 - p), 0.0)
    theta = theta_trop + theta_strat_extra
    T = theta * (p / 1000.0) ** (287.0 / 1004.0)
    out = find_tropopause(T, p)
    assert out["p_curvature"] < 250.0      # high up, near the 150 hPa cap
    assert out["p_curvature"] > 80.0       # but not pinned at the model top


def test_degenerate_short_column_returns_sentinel():
    T = np.array([280.0, 270.0, 260.0])    # < 4 levels
    p = np.array([1000.0, 700.0, 400.0])
    out = find_tropopause(T, p)
    assert out["p_curvature"] == pytest.approx(p[-1])   # model-top sentinel
    assert out["k_coldpoint"] == 2                       # coldest = top here


def test_coldpoint_reported_alongside_curvature():
    T, p = _synthetic_column(p_tropopause_hpa=200.0, nz=40)
    out = find_tropopause(T, p)
    assert set(out) == {"p_curvature", "p_coldpoint", "k_curvature", "k_coldpoint"}
    assert out["p_coldpoint"] == pytest.approx(p[int(np.argmin(T))])
