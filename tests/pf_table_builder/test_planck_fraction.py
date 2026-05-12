import numpy as np


def test_uniform_planck_fraction_sums_to_one():
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from scripts.pf_table_builder.planck_fraction import build_uniform_planck_fraction

    nband, ngpt, nT = 4, 2, 5
    pf = build_uniform_planck_fraction(nband, ngpt, nT)
    assert pf.shape == (nband, ngpt, nT)
    np.testing.assert_allclose(pf.sum(axis=1), 1.0, rtol=1e-12)
    np.testing.assert_allclose(pf, 1.0 / ngpt)
