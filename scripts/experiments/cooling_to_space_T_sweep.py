"""Cooling-to-space sweep over T_iso: does PF's oscillatory low-p cooling track temperature?

If PF's stratospheric oscillations reflect *physics* (Planck-weighted opacity),
the cooling profile should scale and shift with T_iso. If they reflect
*table-interpolation artifacts*, the shape will be largely T-independent.

Three isothermal columns at T = 200 / 250 / 300 K, dry, 376 ppm CO2.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from climt import RRTMGLongwave, get_default_state, get_grid
from climt._components.picket_fence import PicketFenceLongwave

CO2_MOLE_FRAC = 376e-6
T_ISOS = [200.0, 250.0, 300.0]
NZ = 30

LABELS = ["RRTMG (16 bands)", "PF 8 g-pts (GL)", "PF 8 g-pts (2-stretch)"]
HR_KEY = {
    "RRTMG (16 bands)": "air_temperature_tendency_from_longwave",
    "PF 8 g-pts (GL)": "air_temperature_tendency_from_longwave",
    "PF 8 g-pts (2-stretch)": "air_temperature_tendency_from_longwave",
}
T_COLORS = {200.0: "#1f77b4", 250.0: "#2ca02c", 300.0: "#d62728"}

grid = get_grid(nx=1, ny=1, nz=NZ)

rad_lw = RRTMGLongwave()
rrtmg_state = get_default_state([rad_lw], grid_state=grid)

pf_gl = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_8gpt")
pf_2s = PicketFenceLongwave(optics="correlated_k", table="earth_low_res_lw_8gpt_twostretch")

pf_state_gl = get_default_state([pf_gl], grid_state=grid)
pf_state_2s = get_default_state([pf_2s], grid_state=grid)

components = {LABELS[0]: rad_lw, LABELS[1]: pf_gl, LABELS[2]: pf_2s}
states = {LABELS[0]: rrtmg_state, LABELS[1]: pf_state_gl, LABELS[2]: pf_state_2s}

# results[label][T_iso] = {"p_hpa", "hr_day", "olr"}
results = {label: {} for label in LABELS}

for T_iso in T_ISOS:
    for label, st in states.items():
        st["air_temperature"].values[:] = T_iso
        st["surface_temperature"].values[:] = T_iso
        st["specific_humidity"].values[:] = 0.0
        if "mole_fraction_of_carbon_dioxide_in_air" in st:
            st["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_MOLE_FRAC

    for label in LABELS:
        tendencies, diagnostics = components[label](states[label])
        merged = {**tendencies, **diagnostics}
        hr_day = merged[HR_KEY[label]].values[:, 0, 0].copy()
        p_hpa = states[label]["air_pressure"].values[:, 0, 0].copy() / 100.0
        olr = float(merged["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
        results[label][T_iso] = {"p_hpa": p_hpa, "hr_day": hr_day, "olr": olr}

print("Cooling-to-space sweep — OLR vs σT⁴:")
for T_iso in T_ISOS:
    sigma_T4 = 5.670374419e-8 * T_iso**4
    print(f"  T_iso = {T_iso:.0f} K   σT⁴ = {sigma_T4:7.2f} W m⁻²")
    for label in LABELS:
        r = results[label][T_iso]
        print(f"    {label:18s}  OLR = {r['olr']:7.2f}  max|HR| = {np.max(np.abs(r['hr_day'])):.3f} K/day")

# One panel per scheme; overlay the three T_iso curves.
fig, axes = plt.subplots(1, 3, figsize=(14, 6), sharey=True)
for ax, label in zip(axes, LABELS):
    for T_iso in T_ISOS:
        r = results[label][T_iso]
        ax.plot(
            r["hr_day"], r["p_hpa"],
            color=T_COLORS[T_iso], lw=2,
            label=f"T_iso = {T_iso:.0f} K",
        )
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.axvline(0.0, color="grey", lw=0.5)
    ax.set_xlabel("LW heating rate (K day⁻¹)")
    ax.set_title(label)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
axes[0].set_ylabel("Pressure (hPa)")

fig.suptitle(
    "Cooling-to-space vs T_iso — if shape is T-invariant, "
    "it's table-interpolation, not Planck-weighted physics",
    fontsize=11,
)
plt.tight_layout(rect=[0, 0, 1, 0.95])
out_path = os.path.join(os.path.dirname(__file__), "cooling_to_space_T_sweep.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nSaved to {out_path}")
