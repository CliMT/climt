"""Sweep g_split for the two-stretch quadrature on isothermal 250 K column.

Tests whether pushing GL nodes further into the strong-line tail (g→1)
recovers the RRTMG-like cool-to-space at low pressure.

g_split = 0.90, 0.95, 0.97, 0.99. All use 8 g-points total, split 4+4.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from climt import RRTMGLongwave, get_default_state, get_grid
from climt._components.picket_fence import PicketFenceLongwave

CO2_MOLE_FRAC = 376e-6
T_ISO = 250.0
NZ = 30

LABELS = [
    "RRTMG (16 bands)",
    "PF 4-band GL (shipped)",
    "PF 4-band 2s g=0.97",
    "PF 6-band 2s g=0.97 (refined)",
    "PF 6-band GL (refined)",
]
TABLES = {
    "PF 4-band GL (shipped)":        "earth_low_res_lw_8gpt",
    "PF 4-band 2s g=0.97":           "earth_low_res_lw_8gpt_2s097",
    "PF 6-band 2s g=0.97 (refined)": "earth_low_res_lw_co2refined",
    "PF 6-band GL (refined)":        "earth_low_res_lw_co2refined_gl",
}
HR_KEY_RRTMG = "air_temperature_tendency_from_longwave"
HR_KEY_PF = "air_temperature_tendency_from_longwave"
COLORS = ["black", "#999999", "#ff7f0e", "#d62728", "#1f77b4"]
STYLES = ["-", "--", "-", "-", "-"]

grid = get_grid(nx=1, ny=1, nz=NZ)

rad_lw = RRTMGLongwave()
rrtmg_state = get_default_state([rad_lw], grid_state=grid)

pf_components = {}
pf_states = {}
for label, tname in TABLES.items():
    c = PicketFenceLongwave(optics="correlated_k", table=tname)
    pf_components[label] = c
    pf_states[label] = get_default_state([c], grid_state=grid)

all_states = {LABELS[0]: rrtmg_state, **pf_states}
all_components = {LABELS[0]: rad_lw, **pf_components}

# Force isothermal, dry, fixed CO2.
for label, st in all_states.items():
    st["air_temperature"].values[:] = T_ISO
    st["surface_temperature"].values[:] = T_ISO
    st["specific_humidity"].values[:] = 0.0
    if "mole_fraction_of_carbon_dioxide_in_air" in st:
        st["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_MOLE_FRAC

results = {}
for label in LABELS:
    tendencies, diagnostics = all_components[label](all_states[label])
    merged = {**tendencies, **diagnostics}
    key = HR_KEY_RRTMG if label == LABELS[0] else HR_KEY_PF
    hr_day = merged[key].values[:, 0, 0].copy()
    p_hpa = all_states[label]["air_pressure"].values[:, 0, 0].copy() / 100.0
    olr = float(merged["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    results[label] = {"p_hpa": p_hpa, "hr_day": hr_day, "olr": olr}

print(f"Cooling-to-space at T_iso={T_ISO:.0f} K, 376 ppm CO2, dry:")
sigma_T4 = 5.670374419e-8 * T_ISO**4
print(f"  σT⁴ = {sigma_T4:.2f} W/m²")
for label in LABELS:
    r = results[label]
    print(f"  {label:22s}  OLR = {r['olr']:7.2f}  max|HR| = {np.max(np.abs(r['hr_day'])):.3f} K/day")

fig, ax = plt.subplots(figsize=(9, 7))
for label, color, ls in zip(LABELS, COLORS, STYLES):
    r = results[label]
    ax.plot(r["hr_day"], r["p_hpa"], color=color, ls=ls, lw=2, label=label)
ax.set_yscale("log")
ax.invert_yaxis()
ax.axvline(0.0, color="grey", lw=0.5)
ax.set_xlabel("LW heating rate (K day⁻¹)")
ax.set_ylabel("Pressure (hPa)")
ax.set_title(
    f"Cool-to-space: 4-band vs CO2-refined 6-band PF — T_iso = {T_ISO:.0f} K"
)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=9)
plt.tight_layout()
out_path = os.path.join(os.path.dirname(__file__), "cooling_to_space_band_refinement.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nSaved to {out_path}")
