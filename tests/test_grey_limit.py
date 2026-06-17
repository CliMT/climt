"""Grey-limit consistency between CorkLongwaveRadiation and GrayLongwaveRadiation.

A cork scheme with a single band, a single g-point, constant mass-
absorption k, and planck_fraction=1 reduces analytically to grey radiation.
The two climt components then solve the *same* two-stream equations, except
for two conventions:

1. Optical-depth profile: ``GrayLongwaveRadiation`` consumes a prescribed
   ``longwave_optical_depth_on_interface_levels``, while Cork computes
   tau from ``k * column_mass(gas)``. So we build CORK's per-layer tau and feed
   its cumulative surface-to-TOA integral to Gray.

2. Diffusivity factor: Cork multiplies vertical tau by D=1.66 inside
   its two-stream kernel; GrayLongwaveRadiation applies D=1 (it expects tau
   already scaled). So we pre-multiply by D when building Gray's tau input.

With these two matching steps applied, the OLR and layer heating-rate
tendencies from the two schemes agree to machine precision - which is what
"the grey limit of the cork scheme reproduces gray radiation" means
rigorously.
"""
import os

import numpy as np

from climt import GrayLongwaveRadiation, get_default_state
from climt._components.cork import CorkLongwaveRadiation


DIFFUSIVITY_FACTOR = 1.66  # matches climt/_components/cork/lw/kernels.py


def _ensure_unit_table_exists():
    path = os.path.join(
        os.path.dirname(__file__),
        "..", "climt", "_data", "cork", "correlated_k",
        "single_band_unit_lw.nc",
    )
    if not os.path.isfile(path):
        from scripts.generate_single_band_unit_table import (
            write_single_band_unit_lw,
        )
        write_single_band_unit_lw(path, k_value=0.1)


def test_grey_limit_matches_gray_longwave_bit_exact():
    """Single-band, single-gpt, constant-k cork == GrayLongwave
    when tau profile and diffusivity factor are matched."""
    _ensure_unit_table_exists()

    gray = GrayLongwaveRadiation()
    cork = CorkLongwaveRadiation(optics="correlated_k", table="single_band_unit_lw")

    state = get_default_state([gray, cork])

    # Call CORK to compute its vertical optical depth per layer.
    tend_cork, diag_cork = cork(state)
    tau_layer = diag_cork["longwave_optical_depth_per_band"].values[:, 0, 0, 0]

    # Build a cumulative interface-level tau for Gray, pre-multiplied by D
    # to absorb CORK's diffusivity-factor convention.
    nlev = tau_layer.size
    tau_int = np.zeros(nlev + 1)
    tau_int[1:] = np.cumsum(DIFFUSIVITY_FACTOR * tau_layer)

    state["longwave_optical_depth_on_interface_levels"].values[:, 0, 0] = tau_int

    tend_gray, diag_gray = gray(state)

    g = tend_gray["air_temperature"].values
    p = tend_cork["air_temperature"].values
    np.testing.assert_allclose(p, g, rtol=1e-10, atol=1e-18)

    olr_gray = diag_gray["upwelling_longwave_flux_in_air"].values.flat[-1]
    olr_cork = diag_cork["upwelling_longwave_flux_in_air"].values.flat[-1]
    np.testing.assert_allclose(olr_cork, olr_gray, rtol=1e-10)
