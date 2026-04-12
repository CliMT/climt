# climt/_components/picket_fence/optics/parmentier.py
import os

import importlib_resources
import numpy as np

from ..common import njit


def compute_thermal_opacities(kappa_R, gamma_P, beta, R):
    """Compute the two thermal band opacities from Parmentier parameters.

    From Parmentier & Guillot (2014) Eqs. 87-96:
        kappa_R = kappa_1 * kappa_2 / (beta * kappa_2 + (1-beta) * kappa_1)
        R = kappa_1 / kappa_2

    Solving:
        kappa_2 = kappa_R * (beta/R + 1.0 - beta)
        kappa_1 = R * kappa_2

    Args:
        kappa_R: Rosseland mean opacity, m^2/kg
        gamma_P: Planck-to-Rosseland mean ratio (unused directly, but validates)
        beta: fraction of Planck function in high-opacity band
        R: ratio kappa_1/kappa_2

    Returns:
        kappa_1, kappa_2: high and low thermal opacities, m^2/kg
    """
    kappa_2 = kappa_R * (beta / R + 1.0 - beta)
    kappa_1 = R * kappa_2
    return kappa_1, kappa_2


def load_parmentier_coefficients(name_or_path):
    """Load Parmentier coefficient table.

    Args:
        name_or_path: "solar_composition" for built-in, or path to .npz file

    Returns:
        dict-like npz object with coefficient arrays
    """
    if os.path.isfile(name_or_path):
        return np.load(name_or_path)

    data_path = importlib_resources.files(
        "climt._data.picket_fence.parmentier"
    ).joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(data_path) as f:
        return np.load(f)


def lookup_ratio_coefficients(coeffs, T_eff):
    """Look up Parmentier ratio coefficients for a given T_eff.

    Args:
        coeffs: loaded coefficient data (from load_parmentier_coefficients)
        T_eff: effective temperature, K (scalar)

    Returns:
        gamma_v1, gamma_v2, gamma_v3, beta, gamma_P, R
    """
    X = np.log10(max(T_eff, 10.0))  # avoid log10(0)
    boundaries = coeffs["T_eff_boundaries"]

    # Find region index
    region = 0
    for i in range(len(boundaries) - 1):
        if T_eff >= boundaries[i] and T_eff < boundaries[i + 1]:
            region = i
            break

    def _eval_linear(ab_array, region_idx, X_val):
        a, b = ab_array[region_idx]
        return a + b * X_val

    log10_gv3 = _eval_linear(coeffs["log10_gamma_v3_ab"], region, X)
    log10_gv2 = _eval_linear(coeffs["log10_gamma_v2_ab"], region, X)
    log10_gv1 = _eval_linear(coeffs["log10_gamma_v1_ab"], region, X)
    beta_val = _eval_linear(coeffs["beta_ab"], region, X)
    beta_val = np.clip(beta_val, 0.01, 0.99)

    quad = coeffs["log10_gamma_P_quad"]
    log10_gP = quad[0] + quad[1] * X + quad[2] * X**2

    gamma_v1 = 10.0**log10_gv1
    gamma_v2 = 10.0**log10_gv2
    gamma_v3 = 10.0**log10_gv3
    gamma_P = 10.0**log10_gP
    gamma_P = max(gamma_P, 1.0)

    # R from gamma_P and beta using Parmentier & Guillot (2014) Eq. 96
    gp_minus_1 = gamma_P - 1.0
    discriminant = gp_minus_1**2 + 4.0 * beta_val * (1.0 - beta_val) * gp_minus_1
    if discriminant < 0:
        R = 1.0
    else:
        R = (
            1.0
            + gp_minus_1 / (2.0 * beta_val * (1.0 - beta_val))
            + np.sqrt(discriminant) / (2.0 * beta_val * (1.0 - beta_val))
        )
        R = max(R, 1.0)

    return gamma_v1, gamma_v2, gamma_v3, beta_val, gamma_P, R


@njit
def parmentier_lw_optics(T, p, T_eff_arr, kappa_R_params, ratio_coeffs_arr):
    """Compute LW optical depths for Parmentier mode.

    This is a simplified version that takes pre-computed per-column coefficients.

    Args:
        T: (nlev, ncol) temperature, K
        p: (nlev, ncol) pressure, Pa
        T_eff_arr: (ncol,) effective temperature per column
        kappa_R_params: tuple of Freedman 2014 fit parameters
        ratio_coeffs_arr: (ncol, 6) — gamma_v1, gamma_v2, gamma_v3, beta, gamma_P, R per column

    Returns:
        tau: (2, 1, nlev, ncol) optical depth per layer for 2 LW bands, 1 g-point each
        planck_source: (2, 1, nlev, ncol) Planck source per band
        surface_source: (2, 1, ncol) surface Planck source per band
    """
    nlev, ncol = T.shape
    sigma = 5.670374419e-8
    tau = np.zeros((2, 1, nlev, ncol))
    planck_source = np.zeros((2, 1, nlev, ncol))
    surface_source = np.zeros((2, 1, ncol))

    for i in prange(ncol):
        beta_val = ratio_coeffs_arr[i, 3]
        R_val = ratio_coeffs_arr[i, 5]

        for k in range(nlev):
            # Simplified Rosseland mean opacity: constant for now
            # Full Freedman 2014 fit would go here
            kappa_R = 1e-2  # placeholder, m^2/kg
            kappa_2 = kappa_R * (beta_val / R_val + 1.0 - beta_val)
            kappa_1 = R_val * kappa_2

            tau[0, 0, k, i] = kappa_1  # will be multiplied by dp/g in component
            tau[1, 0, k, i] = kappa_2

            planck_val = sigma * T[k, i] ** 4
            planck_source[0, 0, k, i] = beta_val * planck_val
            planck_source[1, 0, k, i] = (1.0 - beta_val) * planck_val

    return tau, planck_source, surface_source
