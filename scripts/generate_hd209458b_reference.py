"""Compute the Parmentier & Guillot (2014) analytic T-p profile for HD 209458b.

Uses the same solar_composition coefficients and Freedman (2014) opacity fit
that the climt Parmentier-mode picket-fence uses, so the converged RCE and this
analytic profile should agree to the level of the time-stepping discretisation.

Output: docs/radiative_transfer/figures/hd209458b_reference.csv
Columns: pressure_bar, temperature_K

Reference: Parmentier & Guillot (2014), A&A 562, A133
           Guillot (2010), A&A 520, A27  (analytic T(tau) formula, Eq. 27)
           Freedman et al. (2014), ApJS 214, 25  (Rosseland mean opacity fit)

HD 209458b parameters from Torres et al. (2008, ApJ 677, 1324):
  g = 9.42 m/s^2,  T_eq = 1459 K,  T_int ~ 500 K (evolutionary model estimate)

Run from the climt conda env:
  conda run -n climt python scripts/generate_hd209458b_reference.py
"""

import os
import numpy as np
from scipy.special import expn          # expn(2, x) = E_2(x)
import importlib_resources

# ---------------------------------------------------------------------------
# HD 209458b parameters (Torres et al. 2008)
# ---------------------------------------------------------------------------
G_PLANET   = 9.42      # m/s^2
T_IRR      = 1450.0    # K  (irradiation temperature, as used in the climt test)
T_INT      = 500.0     # K  (internal temperature, as used in the climt test)
MU_STAR    = 1.0       # cos(zenith_angle) = cos(0)
A_BOND     = 0.0       # bond albedo (zero for reference; test uses bond_albedo_feedback)
F_REDIST   = 1.0       # redistribution weight (1.0 for mu_star=1 direct beam, consistent with test)


# ---------------------------------------------------------------------------
# Load opacity tables from climt data
# ---------------------------------------------------------------------------
def _load(pkg, fname):
    path = importlib_resources.files(pkg).joinpath(fname)
    with importlib_resources.as_file(path) as f:
        return np.load(f)


coeffs   = _load("climt._data.picket_fence.parmentier", "solar_composition.npz")
freeco   = _load("climt._data.picket_fence.parmentier", "freedman2014.npz")


# ---------------------------------------------------------------------------
# Freedman (2014) Rosseland mean opacity  [m^2/kg]
# ---------------------------------------------------------------------------
def kappa_R(T, p_Pa):
    """Rosseland mean opacity in m^2/kg.  T in K, p in Pa."""
    log_T = np.log10(max(float(T), 10.0))
    log_P = np.log10(max(float(p_Pa) * 10.0, 1.0))   # Pa -> dyne/cm^2
    if T < float(freeco["T_boundary"]):
        log_k = float(freeco["a_lo"]) * log_T + float(freeco["b_lo"]) * log_P + float(freeco["c_lo"])
    else:
        log_k = float(freeco["a_hi"]) * log_T + float(freeco["b_hi"]) * log_P + float(freeco["c_hi"])
    return 10.0 ** log_k * 0.1   # cm^2/g -> m^2/kg


# ---------------------------------------------------------------------------
# Parmentier coefficients at T_eff
# ---------------------------------------------------------------------------
T_eff = (T_INT**4 + (1.0 - A_BOND) * MU_STAR * T_IRR**4) ** 0.25
print(f"T_eff = {T_eff:.1f} K")

X = np.log10(max(T_eff, 10.0))
boundaries = coeffs["T_eff_boundaries"]
region = 0
for i in range(len(boundaries) - 1):
    if boundaries[i] <= T_eff < boundaries[i + 1]:
        region = i
        break

def _lin(ab, X):
    return ab[region, 0] + ab[region, 1] * X

gamma_v1 = 10.0 ** _lin(coeffs["log10_gamma_v1_ab"], X)
gamma_v2 = 10.0 ** _lin(coeffs["log10_gamma_v2_ab"], X)
gamma_v3 = 10.0 ** _lin(coeffs["log10_gamma_v3_ab"], X)
gammas   = [gamma_v1, gamma_v2, gamma_v3]
w        = 1.0 / 3.0   # equal weight per visible band (P&G 2014 Eq. 27)
print(f"gamma_v1={gamma_v1:.4f}  gamma_v2={gamma_v2:.4f}  gamma_v3={gamma_v3:.4f}")


# ---------------------------------------------------------------------------
# Guillot (2010) Eq. 27: G(gamma, tau_th)
# See Guillot (2010, A&A 520, A27), also P&G 2014 Eq. 27
# ---------------------------------------------------------------------------
def G_func(gamma, tau):
    """Contribution of one visible band to the irradiation term."""
    tau   = np.asarray(tau, dtype=float)
    E2    = expn(2, gamma * tau)
    term1 = 2.0 / 3.0
    term2 = 2.0 / (3.0 * gamma) * (1.0 + (gamma * tau / 2.0 - 1.0) * np.exp(-gamma * tau))
    term3 = 2.0 * gamma / 3.0 * (1.0 - tau**2 / 2.0) * E2
    return term1 + term2 + term3


def T_analytic(tau_th):
    """Analytic T(tau_th) from P&G 2014 / Guillot 2010 [K]."""
    tau_th = np.asarray(tau_th, dtype=float)
    irr    = (T_IRR**4 / 4.0) * F_REDIST * (1.0 - A_BOND) * MU_STAR
    vis    = sum(w * G_func(gv, tau_th) for gv in gammas)
    T4     = 0.75 * T_INT**4 * (tau_th + 2.0 / 3.0) + irr * vis
    return T4 ** 0.25


# ---------------------------------------------------------------------------
# Convert tau_th grid -> pressure grid
# tau_th(p) = integral_0^p kappa_R(T(p'), p') dp'/g
# Invert numerically: given tau_th, find p such that integral = tau_th
# ---------------------------------------------------------------------------
# Build a fine tau grid then map to pressure via cumulative trapezoidal integration
N_TAU  = 300
tau_grid = np.concatenate([
    np.logspace(-6, -1, N_TAU // 2),
    np.logspace(-1, 2.0, N_TAU // 2),
])
T_grid   = T_analytic(tau_grid)   # K

# Integrate dp/dtau = g / kappa_R(T, p)  with an initial guess p_0 for tau=0 being 0
# We march forward: p[i] = p[i-1] + (tau[i]-tau[i-1]) * g / kappa_R(T[i], p[i-1])
p_grid = np.zeros_like(tau_grid)
p_grid[0] = tau_grid[0] * G_PLANET / kappa_R(T_grid[0], 1.0)   # bootstrap at low tau

for i in range(1, len(tau_grid)):
    dtau   = tau_grid[i] - tau_grid[i - 1]
    p_prev = p_grid[i - 1]
    T_mid  = 0.5 * (T_grid[i] + T_grid[i - 1])
    kR     = kappa_R(T_mid, max(p_prev, 1.0))
    p_grid[i] = p_prev + dtau * G_PLANET / kR

p_bar = p_grid / 1e5   # Pa -> bar

# ---------------------------------------------------------------------------
# Save CSV
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(__file__), "..",
                       "docs", "radiative_transfer", "figures")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "hd209458b_reference.csv")

header = (
    "# Parmentier & Guillot (2014) analytic T-p profile for HD 209458b\n"
    "# Parameters: T_irr=1450 K, T_int=500 K, mu_star=1.0, A_bond=0.0\n"
    "# Opacity: Freedman et al. (2014) Rosseland mean, solar composition\n"
    "# Gamma coefficients: Parmentier et al. (2015) solar_composition.npz\n"
    "# References:\n"
    "#   Torres et al. (2008), ApJ 677, 1324  (g = 9.42 m/s^2)\n"
    "#   Parmentier & Guillot (2014), A&A 562, A133\n"
    "#   Guillot (2010), A&A 520, A27\n"
    "#   Freedman et al. (2014), ApJS 214, 25\n"
    "pressure_bar,temperature_K"
)

rows = "\n".join(f"{p:.8e},{T:.4f}" for p, T in zip(p_bar, T_grid))
with open(out_path, "w") as f:
    f.write(header + "\n" + rows + "\n")

print(f"\nWrote {len(tau_grid)} points to {out_path}")
print(f"Pressure range: {p_bar.min():.4e} – {p_bar.max():.4e} bar")
print(f"Temperature range: {T_grid.min():.1f} – {T_grid.max():.1f} K")
print("\nFirst 5 rows:")
for p, T in zip(p_bar[:5], T_grid[:5]):
    print(f"  p={p:.4e} bar   T={T:.1f} K")
print("\nLast 5 rows:")
for p, T in zip(p_bar[-5:], T_grid[-5:]):
    print(f"  p={p:.4e} bar   T={T:.1f} K")
