"""Compare Cork (correlated-k) against RRTMG on an Earth standard atmosphere.

Cork ships low-resolution (4 LW / 3 SW bands) correlated-k tables
derived from Chaverot's high-resolution R=500 tables (Zenodo 16795590),
with the H2O VMR axis preserved so specific_humidity drives runtime
optical-depth lookup. Given the drastic spectral coarsening, we allow up
to 25% relative error at any level against RRTMG's 16-band LW / 14-band
SW computation.

Note on state building: ``CorkLongwaveRadiation.__init__`` calls
``set_num_longwave_bands(n)`` as a side-effect to size ``num_longwave_bands``
dimensions (e.g. surface emissivity). RRTMG needs 16; CORK needs 4. We run
RRTMG *before* instantiating CORK, then build CORK's state independently.
"""
import numpy as np
import pytest


def _rrtmg_available():
    try:
        from climt import RRTMGLongwave, RRTMGShortwave  # noqa: F401
        return True
    except (ImportError, OSError):
        return False


@pytest.mark.skipif(not _rrtmg_available(), reason="RRTMG not available")
def test_cork_lw_vs_rrtmg_earth_standard():
    """Broadband LW fluxes agree with RRTMG within 25% at all levels."""
    from climt import RRTMGLongwave, get_default_state
    from climt._components.cork import CorkLongwaveRadiation

    rrtmg = RRTMGLongwave()
    state_rr = get_default_state([rrtmg])
    _, d_rr = rrtmg(state_rr)
    up_rr = d_rr["upwelling_longwave_flux_in_air"].values

    cork = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")
    state_cork = get_default_state([cork])
    _, d_cork = cork(state_cork)
    up_cork = d_cork["upwelling_longwave_flux_in_air"].values

    rel_err = np.abs(up_cork - up_rr) / np.maximum(np.abs(up_rr), 1e-3)
    assert np.nanmax(rel_err) < 0.25, (
        f"max relative LW error {np.nanmax(rel_err):.2%} exceeds 25%"
    )


@pytest.mark.skipif(not _rrtmg_available(), reason="RRTMG not available")
def test_cork_sw_vs_rrtmg_earth_standard():
    """Broadband SW down-flux agrees with RRTMG within 25% at all levels."""
    from climt import RRTMGShortwave, get_default_state
    from climt._components.cork import CorkShortwaveRadiation

    rrtmg = RRTMGShortwave()
    state_rr = get_default_state([rrtmg])
    state_rr["zenith_angle"].values[:] = np.pi / 4
    _, d_rr = rrtmg(state_rr)
    dn_rr = d_rr["downwelling_shortwave_flux_in_air"].values

    cork = CorkShortwaveRadiation(optics="correlated_k", table="earth_low_res_sw")
    state_cork = get_default_state([cork])
    state_cork["zenith_angle"].values[:] = np.pi / 4
    _, d_cork = cork(state_cork)
    dn_cork = d_cork["downwelling_shortwave_flux_in_air"].values

    rel_err = np.abs(dn_cork - dn_rr) / np.maximum(np.abs(dn_rr), 1e-3)
    assert np.nanmax(rel_err) < 0.25, (
        f"max relative SW error {np.nanmax(rel_err):.2%} exceeds 25%"
    )


@pytest.mark.skipif(not _rrtmg_available(), reason="RRTMG not available")
def test_cork_lw_responds_to_humidity():
    """OLR should drop as specific humidity rises — proof that the H2O axis
    in the correlated-k table is live at runtime. If the table had been
    collapsed to a fixed VMR (as the first iteration did), OLR would be
    invariant to specific_humidity and this test would fail.
    """
    from climt import get_default_state
    from climt._components.cork import CorkLongwaveRadiation

    cork = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw")

    state_dry = get_default_state([cork])
    state_dry["specific_humidity"].values[:] = 1e-6
    _, d_dry = cork(state_dry)
    olr_dry = float(
        d_dry["upwelling_longwave_flux_in_air"].values.flat[-1]
    )

    state_moist = get_default_state([cork])
    state_moist["specific_humidity"].values[:] = 1e-2
    _, d_moist = cork(state_moist)
    olr_moist = float(
        d_moist["upwelling_longwave_flux_in_air"].values.flat[-1]
    )

    assert olr_moist < olr_dry, (
        f"OLR should fall with more H2O but got dry={olr_dry:.2f}, "
        f"moist={olr_moist:.2f}"
    )
    # Require at least a 5 W/m^2 drop — small enough to accommodate the
    # coarse 4-band spectral resolution, large enough to prove the table
    # is genuinely responding to humidity (not a rounding artefact).
    assert (olr_dry - olr_moist) > 5.0, (
        f"OLR only dropped {olr_dry - olr_moist:.2f} W/m^2 from dry→moist"
    )
