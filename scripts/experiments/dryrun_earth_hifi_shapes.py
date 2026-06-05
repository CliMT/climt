"""Validate earth_hifi k_coefficients axis ordering with synthetic kappa.

Run: conda run -n climt python scripts/experiments/dryrun_earth_hifi_shapes.py
"""
import os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_HERE, ".."))
from pf_table_builder.k_distribution import kappa_to_k_coeffs

nT, nP, nX, nC, ngpt = 12, 8, 7, 10, 8
band_edges = np.array([10.0, 250.0, 350.0, 500.0, 630.0, 700.0, 800.0,
                       980.0, 1080.0, 1180.0, 1250.0, 1400.0, 1600.0,
                       1800.0, 3250.0])
nband = len(band_edges) - 1
nu = np.arange(band_edges[0], band_edges[-1], 1.0)
kappa = np.abs(np.random.RandomState(0).standard_normal(
    (nT, nP, nX, nC, len(nu)))) * 1e-4

k, w = kappa_to_k_coeffs(kappa, nu, band_edges, ngpt)
k = k[np.newaxis]
k = np.moveaxis(k, (1, 2, 3, 4, 5, 6), (3, 4, 5, 6, 1, 2))
assert k.shape == (1, nband, ngpt, nT, nP, nX, nC), f"Got {k.shape}"
print(f"OK shape = {k.shape}")
