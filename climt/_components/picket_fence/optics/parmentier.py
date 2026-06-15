# climt/_components/picket_fence/optics/parmentier.py
import os

import importlib_resources
import numpy as np


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


def load_freedman2014_coefficients():
    """Load the Freedman et al. (2014) Rosseland mean opacity fit coefficients.

    Returns a dict with keys: T_boundary, a_hi, b_hi, c_hi, a_lo, b_lo, c_lo.
    The formula is:
        log10(kappa_R [cm^2/g]) = a * log10(T [K]) + b * log10(P [dyne/cm^2]) + c
    where the (a, b, c) set is chosen based on whether T < T_boundary or not.
    Convert to m^2/kg by multiplying by 0.1.

    Reference: Freedman et al. (2014), ApJS 214, 25. The coefficients here are
    approximate values calibrated to reproduce the qualitative behaviour of their
    Table 3 (solar metallicity, grain-free Rosseland mean); exact published values
    from Table 1 should be substituted once verified.
    """
    data_path = importlib_resources.files(
        "climt._data.picket_fence.parmentier"
    ).joinpath("freedman2014.npz")
    with importlib_resources.as_file(data_path) as f:
        with np.load(f) as npz:
            return {k: npz[k] for k in npz.files}


def compute_rosseland_mean_opacity(T, p, coeffs):
    """Compute the Freedman et al. (2014) Rosseland mean opacity.

    Args:
        T: temperature, K (scalar)
        p: pressure, Pa (scalar)
        coeffs: coefficient dict from load_freedman2014_coefficients()

    Returns:
        kappa_R: Rosseland mean opacity, m^2/kg (scalar)
    """
    log_T = np.log10(max(float(T), 10.0))
    log_P = np.log10(max(float(p) * 10.0, 1.0))  # Pa → dyne/cm^2
    if T < float(coeffs["T_boundary"]):
        log_k = float(coeffs["a_lo"]) * log_T + float(coeffs["b_lo"]) * log_P + float(coeffs["c_lo"])
    else:
        log_k = float(coeffs["a_hi"]) * log_T + float(coeffs["b_hi"]) * log_P + float(coeffs["c_hi"])
    return 10.0**log_k * 0.1  # cm^2/g → m^2/kg


def load_parmentier_coefficients(name_or_path):
    """Load Parmentier coefficient table.

    Args:
        name_or_path: "solar_composition" for built-in, or path to .npz file

    Returns:
        dict-like npz object with coefficient arrays
    """
    if os.path.isfile(name_or_path):
        with np.load(name_or_path) as npz:
            return {k: npz[k] for k in npz.files}

    data_path = importlib_resources.files(
        "climt._data.picket_fence.parmentier"
    ).joinpath(f"{name_or_path}.npz")
    with importlib_resources.as_file(data_path) as f:
        with np.load(f) as npz:
            return {k: npz[k] for k in npz.files}


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


def bond_albedo_from_fluxes(up_toa, down_toa):
    """Estimate the Bond albedo from top-of-atmosphere SW fluxes.

    Args:
        up_toa: upwelling SW flux at TOA, W/m^2 (ncol,)
        down_toa: downwelling SW flux at TOA, W/m^2 (ncol,)

    Returns:
        A_B: bond albedo per column, dimensionless, clipped to [0, 1].
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        A_B = np.where(down_toa > 0, up_toa / down_toa, 0.0)
    return np.clip(A_B, 0.0, 1.0)
