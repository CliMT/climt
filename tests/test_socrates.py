import numpy as np
from climt import SocratesLongwave, SocratesShortwave, get_default_state


def test_socrates_lw_co2_sensitivity():
    """Verify SocratesLongwave runs and OLR decreases as CO2 increases."""
    lw = SocratesLongwave()

    # Low CO2 state
    state_low = get_default_state([lw])
    state_low["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 0.0
    _, diag_low = lw(state_low)
    olr_low = diag_low["upwelling_longwave_flux_in_air"].values.flat[-1]

    # High CO2 state
    state_high = get_default_state([lw])
    state_high["mole_fraction_of_carbon_dioxide_in_air"].values[:] = 1000e-6
    _, diag_high = lw(state_high)
    olr_high = diag_high["upwelling_longwave_flux_in_air"].values.flat[-1]

    # OLR should drop with higher CO2 (greenhouse effect)
    assert olr_high < olr_low
    assert not np.isnan(olr_high)
    assert not np.isnan(olr_low)


def test_socrates_sw_zenith():
    """Verify SocratesShortwave runs and SW flux decreases with zenith."""
    sw = SocratesShortwave()

    # Overhead sun
    state_overhead = get_default_state([sw])
    state_overhead["zenith_angle"].values[:] = 0.0
    _, diag_overhead = sw(state_overhead)
    sw_dn_overhead = (
        diag_overhead["downwelling_shortwave_flux_in_air"].values.flat[0]
    )

    # Slanted sun
    state_slanted = get_default_state([sw])
    state_slanted["zenith_angle"].values[:] = np.pi / 3.0
    _, diag_slanted = sw(state_slanted)
    sw_dn_slanted = (
        diag_slanted["downwelling_shortwave_flux_in_air"].values.flat[0]
    )

    # Overhead sun should yield more downwelling shortwave radiation
    assert sw_dn_slanted < sw_dn_overhead
    assert not np.isnan(sw_dn_overhead)
    assert not np.isnan(sw_dn_slanted)
