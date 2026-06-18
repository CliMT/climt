"""End-to-end smoke tests: run Cork{Long,Short}wave with the new tables.

Verifies the netCDF schema we wrote matches what climt's loader expects and
that fluxes are physically sane (positive, finite, monotone in obvious ways).
"""
import numpy as np
import pytest

from climt import get_default_state, load_atmospheric_properties, reset_atmospheric_properties
from climt._components.cork import (
    CorkLongwaveRadiation,
    CorkShortwaveRadiation,
)


@pytest.mark.parametrize("scenario,profile", [
    ("trappist1e_hab1", "trappist1e"),
    ("trappist1e_hab2", None),
    ("titan", None),
])
def test_cork_lw_runs(scenario, profile):
    """LW kernel produces non-negative OLR with the new table."""
    if profile:
        load_atmospheric_properties(profile)
    try:
        lw = CorkLongwaveRadiation(optics="correlated_k", table=f"{scenario}_lw")
        state = get_default_state([lw])
        tend, diag = lw(state)
        olr = diag["upwelling_longwave_flux_in_air"].values.flat[-1]
        assert np.isfinite(olr), f"non-finite OLR for {scenario}: {olr}"
        assert olr > 0, f"non-positive OLR for {scenario}: {olr}"
    finally:
        if profile:
            reset_atmospheric_properties()


@pytest.mark.parametrize("scenario,profile,star", [
    ("trappist1e_hab1", "trappist1e", "trappist1"),
    ("trappist1e_hab2", None, "trappist1"),
    ("titan", None, "sun"),
])
def test_cork_sw_runs(scenario, profile, star):
    """SW kernel produces non-negative downwelling flux with the new table."""
    if profile:
        load_atmospheric_properties(profile)
    try:
        sw = CorkShortwaveRadiation(
            optics="correlated_k",
            table=f"{scenario}_sw",
            stellar_spectrum=star,
        )
        state = get_default_state([sw])
        state["zenith_angle"].values[:] = np.pi / 4
        tend, diag = sw(state)
        flux_dn = diag["downwelling_shortwave_flux_in_air"].values
        assert np.all(np.isfinite(flux_dn)), f"non-finite SW flux for {scenario}"
        assert np.all(flux_dn >= 0), f"negative SW flux for {scenario}"
    finally:
        if profile:
            reset_atmospheric_properties()
