"""The longwave diffusivity factor is settable per component.

The EC2213 course notes derive the two-stream equations with D = 2 (mu = 0.5);
cork's default is the Elsasser value 1.66. The modelling-tour pages need to
reproduce the notes' algebra exactly, so D must be a constructor argument.
"""
import numpy as np
import pytest
import sympl

import climt
from climt import CorkLongwaveRadiation, get_default_state, get_grid


@pytest.fixture(autouse=True)
def _unyt_backend():
    sympl.set_backend(climt.UnytBackend())


def _run(diffusivity_factor=None):
    kwargs = {"optics": "correlated_k", "table": "earth_low_res_lw"}
    if diffusivity_factor is not None:
        kwargs["diffusivity_factor"] = diffusivity_factor
    lw = CorkLongwaveRadiation(**kwargs)
    state = get_default_state([lw], grid_state=get_grid(nx=1, ny=1, nz=18))
    return lw(state)


def test_default_is_elsasser_166():
    _, diag = _run()
    tau = diag["longwave_optical_depth_per_band"].values
    trans = diag["longwave_transmittance_per_band"].values
    np.testing.assert_allclose(trans, np.exp(-1.66 * tau), rtol=1e-6)


def test_transmittance_follows_the_requested_factor():
    _, diag = _run(diffusivity_factor=2.0)
    tau = diag["longwave_optical_depth_per_band"].values
    trans = diag["longwave_transmittance_per_band"].values
    np.testing.assert_allclose(trans, np.exp(-2.0 * tau), rtol=1e-6)


def test_larger_diffusivity_lowers_olr():
    """A longer effective path means a higher, colder emission level."""
    _, diag_166 = _run(diffusivity_factor=1.66)
    _, diag_2 = _run(diffusivity_factor=2.0)
    olr_166 = float(diag_166["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    olr_2 = float(diag_2["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    assert olr_2 < olr_166
