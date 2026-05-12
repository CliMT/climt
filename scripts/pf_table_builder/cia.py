"""Load HITRAN CIA flat files and evaluate κ_CIA(T, p, ν) on a (T, p) grid.

HITRAN CIA distributes binary collision cross-sections per temperature
block in cm⁻¹ amagat⁻². We convert to a per-mass absorption coefficient
(m²/kg of total moist gas) so the result can be added directly to the
gas-line kappa array from kappa_sampling.sample_kappa_grid.
"""
from __future__ import annotations

import numpy as np


_AMAGAT_NUMBER_DENSITY = 2.6868e25  # molecules / m³ at STP
_N_A = 6.022e23


def load_cia_blocks(path):
    """Parse a HITRAN CIA file. Returns dict: T (K) -> (nu_array, sigma_array).

    sigma units: cm⁻¹ amagat⁻² (the HITRAN file convention).
    """
    blocks = {}
    with open(path, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        header = lines[i].split()
        if len(header) < 5:
            i += 1
            continue
        nPoints = int(header[3])
        T = float(header[4])
        data = np.loadtxt(lines[i + 1: i + 1 + nPoints])
        blocks[T] = (data[:, 0], data[:, 1])
        i += 1 + nPoints
    return blocks


def cia_kappa_on_grid(
    path, *, pair, vmr_a, vmr_b,
    background_gas, absorbers,
    T_grid, p_grid, nu_grid,
):
    """Evaluate κ_CIA on (T, p, ν) and convert to m²/kg of moist gas.

    Args:
        path: HITRAN .cia flat file.
        pair: tuple (species_a, species_b); used for the mass-fraction
            normalisation only (the file already encodes which pair).
        vmr_a, vmr_b: volume mixing ratios of the two colliders.
        background_gas, absorbers: same conventions as kappa_sampling, used
            only to compute the mean molar mass for unit conversion.
        T_grid, p_grid, nu_grid: target grids.

    Returns:
        kappa: (nT, nP, nNu), m²/kg.
    """
    blocks = load_cia_blocks(path)
    T_file = np.array(sorted(blocks.keys()))
    nT_f = len(T_file)
    # Resample each block onto nu_grid (zero outside file's spectral range)
    sigma_grid = np.zeros((nT_f, len(nu_grid)))
    for iT, T in enumerate(T_file):
        nu_f, sig_f = blocks[T]
        mask = (nu_grid >= nu_f[0]) & (nu_grid <= nu_f[-1])
        sigma_grid[iT, mask] = np.interp(nu_grid[mask], nu_f, sig_f)

    # Interpolate in T (clamp outside file range)
    nT, nP, nNu = len(T_grid), len(p_grid), len(nu_grid)
    sigma_TpNu = np.zeros((nT, nNu))
    for inu in range(nNu):
        sigma_TpNu[:, inu] = np.interp(T_grid, T_file, sigma_grid[:, inu],
                                       left=sigma_grid[0, inu],
                                       right=sigma_grid[-1, inu])

    # CIA absorption coefficient (m^-1):
    #   α(ν) = σ(ν, T) × n_a × n_b
    # where n_x is number density of colliders in amagat. amagat = n/n_STP.
    # For an ideal gas:  n_amagat = (p/p_STP) × (T_STP/T)
    p_STP, T_STP = 1.01325e5, 273.15
    n_a = (p_grid * vmr_a / p_STP)[:, None] * (T_STP / T_grid)[None, :]  # (nP, nT)
    n_b = (p_grid * vmr_b / p_STP)[:, None] * (T_STP / T_grid)[None, :]

    # Convert σ from cm⁻¹ amagat⁻² to m⁻¹ amagat⁻²: ×100
    sigma_si = sigma_TpNu * 100.0  # (nT, nNu)

    # κ_volume [m⁻¹] = σ × n_a × n_b  (per amagat²)
    # broadcast (nT, 1, nNu) × (nT, nP) × (nT, nP) → (nT, nP, nNu)
    n_a = n_a.T  # (nT, nP)
    n_b = n_b.T  # (nT, nP)
    kappa_vol = sigma_si[:, None, :] * (n_a * n_b)[:, :, None]

    # Convert to mass absorption coefficient: divide by mass density of moist gas
    # ρ(T, p) = p × M_mean / (R_universal × T)
    R = 8.31446
    _MM = {"H2O": 18.015, "CO2": 44.01, "O3": 47.998, "CH4": 16.04,
           "NH3": 17.031, "N2": 28.014, "air": 28.97}
    f_tot = sum(absorbers.values())
    M_mean_g = sum(v * _MM[s] for s, v in absorbers.items())
    if background_gas is not None:
        M_mean_g += (1.0 - f_tot) * _MM[background_gas]
    M_mean = M_mean_g * 1e-3  # kg/mol
    rho = (p_grid[None, :] * M_mean) / (R * T_grid[:, None])  # (nT, nP) kg/m^3

    kappa_mass = kappa_vol / rho[:, :, None]
    return kappa_mass
