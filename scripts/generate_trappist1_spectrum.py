"""Generate climt/_data/picket_fence/stellar_spectra/trappist1.npz.

Approach: Approximate TRAPPIST-1 (M8V, T_eff ~ 2566 K) as a blackbody,
scaled so the integrated flux at 1 AU equals 55.0 W/m^2 (= L_star / 4 pi AU^2
with L_star = 0.000553 L_sun). For a higher-fidelity SED, replace this
blackbody scaffold with the MUSCLES / HST observed spectrum.

Run once from repo root; commit the resulting .npz.
"""
import os
import numpy as np

OUT = os.path.join(
    "climt", "_data", "picket_fence", "stellar_spectra", "trappist1.npz"
)

T_EFF = 2566.0          # effective temperature, K
TARGET_AT_1AU = 55.0    # W/m^2, integrated irradiance at 1 AU

h = 6.62607015e-34
c = 2.99792458e8
k = 1.380649e-23

# Wavenumber grid (cm^-1), log-spaced 10 -> 30000
wavenumber = np.geomspace(10.0, 30000.0, 4000)

# Planck function B(nu, T) in W/m^2/sr/(cm^-1):
freq = wavenumber * 100.0 * c    # Hz
B_Hz = (2.0 * h * freq**3 / c**2) / (np.exp(h * freq / (k * T_EFF)) - 1.0)
B_cm = B_Hz * (100.0 * c)        # W/m^2/sr/(cm^-1)

# Hemispheric flux at the stellar surface = pi * B; then scale to 1 AU.
irradiance = np.pi * B_cm        # W/m^2/(cm^-1) at the stellar photosphere

# Normalise so integral irradiance dwn = TARGET_AT_1AU
total = np.trapezoid(irradiance, wavenumber)
irradiance *= TARGET_AT_1AU / total

np.savez(OUT, wavenumber=wavenumber, irradiance=irradiance, source="blackbody_T2566K")
print(f"wrote {OUT} ({wavenumber.size} points, total {TARGET_AT_1AU:.1f} W/m^2)")
