"""``single_band_gray_lw`` reproduces climt's *default* gray longwave radiation.

``get_default_state`` seeds ``longwave_optical_depth_on_interface_levels`` with
the linear profile ``tau0 * (1 - p/ps)`` (tau0 = 1). The ``single_band_gray_lw``
cork table is a constant-k, single-band fixture calibrated (see
``scripts/generate_gray_default_table.py``) so that cork's diffusivity-scaled
column optical depth equals that default gray tau. Running the two components on
the default state must then give the same heating rates and OLR.
"""
import numpy as np

from climt import GrayLongwaveRadiation, get_default_state, get_grid
from climt._components.cork import CorkLongwaveRadiation

DIFFUSIVITY_FACTOR = 1.66  # matches climt/_components/cork/lw/kernels.py


def test_gray_default_table_matches_default_gray_radiation():
    gray = GrayLongwaveRadiation()
    cork = CorkLongwaveRadiation(optics="correlated_k", table="single_band_gray_lw")
    state = get_default_state([gray, cork], grid_state=get_grid(nx=1, ny=1, nz=28))

    tend_gray, diag_gray = gray(state)   # consumes default-state tau (tau0=1)
    tend_cork, diag_cork = cork(state)

    # cork's diffusivity-scaled cumulative optical depth reproduces default tau.
    tau_layer = diag_cork["longwave_optical_depth_per_band"].values[:, 0, 0, 0]
    tau_cork_int = np.concatenate([[0.0], np.cumsum(DIFFUSIVITY_FACTOR * tau_layer)])
    tau_default = state["longwave_optical_depth_on_interface_levels"].values[:, 0, 0]
    np.testing.assert_allclose(tau_cork_int, tau_default, atol=2e-3)

    # Same longwave heating rates and OLR.
    np.testing.assert_allclose(
        tend_cork["air_temperature"].values,
        tend_gray["air_temperature"].values,
        rtol=1e-2, atol=1e-7,
    )
    olr_gray = diag_gray["upwelling_longwave_flux_in_air"].values.flat[-1]
    olr_cork = diag_cork["upwelling_longwave_flux_in_air"].values.flat[-1]
    np.testing.assert_allclose(olr_cork, olr_gray, rtol=1e-3)
