import os

import numpy as np

# T_eff boundaries for piecewise fits (from Parmentier et al. 2015, Table 1)
T_eff_boundaries = np.array([0.0, 200.0, 300.0, 600.0, 1400.0, 2000.0, np.inf])

# Coefficients: each is shape (n_regions, 2) for (a, b) in log10(coeff) = a + b*X
# where X = log10(T_eff). Regions are: <200, 200-300, 300-600, 600-1400, 1400-2000, >2000

log10_gamma_v3_ab = np.array(
    [
        [-3.03, -0.2],  # T_eff < 200 K
        [-13.87, 4.51],  # 200-300 K
        [-11.95, 3.74],  # 300-600 K
        [-6.97, 1.94],  # 600-1400 K
        [-3.65, 0.89],  # 1400-2000 K
        [-6.02, 1.61],  # > 2000 K
    ]
)

log10_gamma_v2_ab = np.array(
    [
        [-7.37, 2.53],
        [13.99, -6.75],
        [-15.18, 5.02],
        [-10.41, 3.31],
        [-19.95, 6.34],
        [13.56, -3.81],
    ]
)

log10_gamma_v1_ab = np.array(
    [
        [-5.51, 1.23],
        [8.65, -0.45],
        [-3.45, 4.33],  # Note: this was corrected from paper
        [-12.96, 4.33],
        [-23.75, 7.76],
        [12.65, -3.27],
    ]
)

beta_ab = np.array(
    [
        [0.84, 0.0],
        [0.84, 0.0],
        [0.84, 0.0],
        [0.84, 0.0],
        [0.84, 0.0],
        [6.21, -1.63],
    ]
)

log10_gamma_P_ab = np.array(
    [
        [-2.36, 0.0],  # T_eff < 200 K: constant
        [-2.36, 0.0],  # 200-300 K: use same as <200 (gamma_P not well defined here)
        [13.92, 0.0],  # 300-600 K: a + bX + cX^2, but we store as quadratic below
        [13.92, 0.0],  # placeholder — handled by quadratic
        [13.92, 0.0],  # placeholder
        [13.92, 0.0],  # placeholder
    ]
)

# gamma_P is a full quadratic: log10(gamma_P) = a + b*X + c*X^2
# From Parmentier et al. 2015, Table 1 footnote: a=-19.38, b=13.92, c=-2.36
log10_gamma_P_quad = np.array([-19.38, 13.92, -2.36])  # a + b*X + c*X^2

data_dir = os.path.join("climt", "_data", "picket_fence", "parmentier")
os.makedirs(data_dir, exist_ok=True)

np.savez(
    os.path.join(data_dir, "solar_composition.npz"),
    T_eff_boundaries=T_eff_boundaries,
    log10_gamma_v3_ab=log10_gamma_v3_ab,
    log10_gamma_v2_ab=log10_gamma_v2_ab,
    log10_gamma_v1_ab=log10_gamma_v1_ab,
    beta_ab=beta_ab,
    log10_gamma_P_quad=log10_gamma_P_quad,
)
print(f"Saved solar_composition.npz in {data_dir}")
