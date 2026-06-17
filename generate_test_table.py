# Script to generate test_2band_lw.npz — run once
import os

import numpy as np

ngas = 1
nband = 2
ngpt = 2
nT = 5
nP = 5

T_grid = np.linspace(200.0, 320.0, nT)
p_grid_log = np.linspace(np.log(100.0), np.log(1e5), nP)

# Synthetic k-coefficients: band 0 is opaque, band 1 is transparent
# k increases with T (like H2O absorption)
k = np.zeros((ngas, nband, ngpt, nT, nP))
for iT in range(nT):
    for iP in range(nP):
        # Band 0, g-point 0 (weak): 1e-4
        k[0, 0, 0, iT, iP] = 1e-4 * (1.0 + 0.01 * (T_grid[iT] - 260.0))
        # Band 0, g-point 1 (strong): 1e-1
        k[0, 0, 1, iT, iP] = 1e-1 * (1.0 + 0.01 * (T_grid[iT] - 260.0))
        # Band 1, g-point 0 (weak): 1e-5
        k[0, 1, 0, iT, iP] = 1e-5
        # Band 1, g-point 1 (strong): 1e-3
        k[0, 1, 1, iT, iP] = 1e-3

gpoint_weights = np.array(
    [
        [0.6, 0.4],  # band 0
        [0.7, 0.3],  # band 1
    ]
)

band_limits = np.array(
    [
        [10.0, 500.0],  # band 0: 10-500 cm^-1
        [500.0, 1500.0],  # band 1: 500-1500 cm^-1
    ]
)

planck_fraction = np.zeros((nband, ngpt, nT))
for iT in range(nT):
    planck_fraction[0, 0, iT] = 0.6 * 0.5  # band 0 gets 50% of Planck, split by g-point
    planck_fraction[0, 1, iT] = 0.4 * 0.5
    planck_fraction[1, 0, iT] = 0.7 * 0.5
    planck_fraction[1, 1, iT] = 0.3 * 0.5

data_dir = os.path.join("climt", "_data", "cork", "correlated_k")
os.makedirs(data_dir, exist_ok=True)

np.savez(
    os.path.join(data_dir, "test_2band_lw.npz"),
    k_coefficients=k,
    gpoint_weights=gpoint_weights,
    band_wavenumber_limits=band_limits,
    temperature_grid=T_grid,
    pressure_grid_log=p_grid_log,
    planck_fraction=planck_fraction,
    gas_names=np.array(["h2o"]),
    overlap_method="additive",
    resolution="test",
)
