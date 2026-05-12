"""Planck-source distribution across g-points within each LW band."""
import numpy as np


def build_uniform_planck_fraction(nband, ngpt, nT):
    """Uniform 1/ngpt fraction.

    Matches the convention used by climt's existing low-res tables. A
    higher-fidelity scheme would weight by the sub-band wavenumber slice
    that each g-point represents, but with κ sorted by value (not by ν)
    that mapping has been thrown away. Uniform is conservative and exact
    in the band-averaged sense.
    """
    return np.full((nband, ngpt, nT), 1.0 / ngpt)
