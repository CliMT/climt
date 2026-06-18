"""Experiment: boost CO₂-wing opacity in the 6-band-GL CORK table and check
whether T_strat moves toward RRTMG's ~190 K.

Hypothesis (from 2026-05-18 per-band diagnostic): CORK's stratospheric T is
pinned at ~155 K by the CO₂ core (630-700 cm⁻¹). The wing bands (500-630
and 700-800 cm⁻¹) at 100 hPa look like windows (BT > 240 K), implying κ
≈ 0 aloft. If we artificially boost wing-band k by a factor B, the wings
should become "active" at 100 hPa, their BT should drop toward 190 K, and
T_strat should warm.

Method:
  1. Load earth_low_res_lw_co2refined_gl.nc
  2. Write three boosted variants: B = 3, 10, 30 on bands 1 and 3
  3. Run 1000-day pure-LW radiative equilibrium for each + the baseline
  4. Report final T_strat and per-band BT at 100 hPa

Outputs: examples/experiment_co2_wing_boost.png and console log.
"""
import os
import shutil
from datetime import timedelta

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from sympl import AdamsBashforth
from tqdm import trange

from climt import SlabSurface, get_default_state, get_grid
from climt._components.cork import CorkLongwaveRadiation

# ── Physical constants for band-integrated Planck brightness temperature ──────
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
K_BOLTZ = 1.380649e-23
SIGMA_SB = 5.670374419e-8


def planck_band_flux(T, nu1_cm, nu2_cm, n_sub=40):
    nu = np.linspace(nu1_cm, nu2_cm, n_sub) * 100.0
    T = np.atleast_1d(T).astype(float)
    x = (H_PLANCK * C_LIGHT * nu[:, None]) / (K_BOLTZ * np.maximum(T[None, :], 1.0))
    B = (2.0 * H_PLANCK * C_LIGHT ** 2 * nu[:, None] ** 3) / np.expm1(x)
    flux = np.pi * np.trapz(B, nu, axis=0)
    return flux if flux.size > 1 else flux[0]


def band_brightness_T(F_band, nu1_cm, nu2_cm):
    if F_band <= 0:
        return np.nan
    try:
        return brentq(lambda T: planck_band_flux(T, nu1_cm, nu2_cm) - F_band, 50.0, 400.0)
    except ValueError:
        return np.nan


# ── Configuration ────────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "climt", "_data",
                        "cork", "correlated_k")
BASE_TABLE = "earth_low_res_lw_co2refined_gl"
WING_BAND_IDX = [1, 3]  # 500-630 and 700-800 cm⁻¹ (CO₂ wings)
BOOST_FACTORS = [1.0, 3.0, 10.0, 30.0]
LABELS = [f"baseline" if b == 1.0 else f"wing×{int(b)}" for b in BOOST_FACTORS]

SW_FLUX = 240.0
CO2_MOLE_FRAC = 376e-6
MIXED_LAYER_DEPTH = 1.0
TIMESTEP = timedelta(hours=3)
N_STEPS = 8000   # 1000 days
NZ = 30
T_EFF = (SW_FLUX / SIGMA_SB) ** 0.25
T_SKIN = T_EFF / (2.0 ** 0.25)
DIAG_PRESSURE_HPA = 100.0


def write_boosted_table(boost):
    """Write a copy of BASE_TABLE with k_coefficients × boost in WING_BAND_IDX."""
    src = os.path.join(DATA_DIR, f"{BASE_TABLE}.nc")
    name = f"{BASE_TABLE}_wingboost{int(boost)}x"
    dst = os.path.join(DATA_DIR, f"{name}.nc")
    if boost == 1.0:
        return BASE_TABLE
    # Load fully into memory before writing back to the same path.
    ds = xr.open_dataset(src).load()
    k = ds["k_coefficients"].values.copy()
    for b in WING_BAND_IDX:
        k[:, b, :, :, :, :] *= boost
    ds["k_coefficients"].values = k
    if os.path.exists(dst):
        os.remove(dst)
    ds.to_netcdf(dst, mode="w", format="NETCDF3_CLASSIC")
    ds.close()
    return name


# ── Build all models ─────────────────────────────────────────────────────────
grid = get_grid(nx=1, ny=1, nz=NZ)

table_names = [write_boosted_table(b) for b in BOOST_FACTORS]
print("Tables prepared:")
for lbl, tn in zip(LABELS, table_names):
    print(f"  {lbl:12s} → {tn}")

components, surfaces, states = {}, {}, {}
for lbl, tn in zip(LABELS, table_names):
    lw = CorkLongwaveRadiation(optics="correlated_k", table=tn)
    sf = SlabSurface()
    st = get_default_state([lw, sf], grid_state=grid)
    st["specific_humidity"].values[:] = 0.0
    st["downwelling_shortwave_flux_in_air"].values[:] = SW_FLUX
    st["ocean_mixed_layer_thickness"].values[:] = MIXED_LAYER_DEPTH
    st["area_type"].values[:] = "sea"
    st["air_temperature"].values[:] = T_SKIN
    st["surface_temperature"].values[:] = T_SKIN
    components[lbl] = lw
    surfaces[lbl] = sf
    states[lbl] = st

steppers = {lbl: AdamsBashforth([components[lbl], surfaces[lbl]]) for lbl in LABELS}
history = {lbl: [] for lbl in LABELS}
STORE_EVERY = 100

# ── Integrate ────────────────────────────────────────────────────────────────
for step in trange(N_STEPS, desc="Integrating"):
    for lbl in LABELS:
        diag, new_state = steppers[lbl](states[lbl], TIMESTEP)
        new_state.update(diag)
        states[lbl] = new_state
        if step % STORE_EVERY == 0:
            up_lw = diag["upwelling_longwave_flux_in_air"].values[:, 0, 0].copy()
            history[lbl].append({
                "step": step,
                "p_hpa": new_state["air_pressure"].values[:, 0, 0].copy() / 100.0,
                "p_int_hpa": new_state["air_pressure_on_interface_levels"].values[:, 0, 0].copy() / 100.0,
                "T_profile": new_state["air_temperature"].values[:, 0, 0].copy(),
                "up_lw": up_lw,
                "up_lw_per_band": diag["upwelling_longwave_flux_in_air_per_band"]
                                  .values[:, 0, 0, :].copy(),
                "olr": float(up_lw[-1]),
                "T_surface": float(new_state["surface_temperature"].values[0, 0]),
            })

# ── Diagnostic report ────────────────────────────────────────────────────────
band_edges = xr.open_dataset(os.path.join(DATA_DIR, f"{BASE_TABLE}.nc"))[
    "band_wavenumber_limits"].values

print(f"\nFinal state (1000 days, p ≈ {DIAG_PRESSURE_HPA:.0f} hPa):")
print(f"{'label':14s} {'T_sfc':>8s} {'OLR':>8s} {'T_strat':>8s}   per-band BT (K)")
for lbl in LABELS:
    snap = history[lbl][-1]
    p_int = snap["p_int_hpa"]
    p_mid = snap["p_hpa"]
    k_diag = int(np.argmin(np.abs(p_int - DIAG_PRESSURE_HPA)))
    k_mid = int(np.argmin(np.abs(p_mid - DIAG_PRESSURE_HPA)))
    F_diag = snap["up_lw_per_band"][k_diag, :]
    BT = [band_brightness_T(F_diag[b], band_edges[b, 0], band_edges[b, 1])
          for b in range(band_edges.shape[0])]
    BT_str = "  ".join(f"b{b}:{BT[b]:6.1f}" for b in range(len(BT)))
    print(f"{lbl:14s} {snap['T_surface']:8.2f} {snap['olr']:8.2f} "
          f"{snap['T_profile'][k_mid]:8.2f}   {BT_str}")

# ── Plot final T profiles ────────────────────────────────────────────────────
fig, (ax_T, ax_BT) = plt.subplots(1, 2, figsize=(13, 7))
colors = plt.get_cmap("viridis")(np.linspace(0, 0.9, len(LABELS)))

for lbl, c in zip(LABELS, colors):
    snap = history[lbl][-1]
    ax_T.plot(snap["T_profile"], snap["p_hpa"], color=c, lw=2, label=lbl)
ax_T.set_yscale("log"); ax_T.invert_yaxis()
ax_T.set_xlabel("Temperature (K)"); ax_T.set_ylabel("Pressure (hPa)")
ax_T.set_title("Final T profile vs CO₂-wing boost factor")
ax_T.axvline(190, color="k", ls=":", lw=1, alpha=0.5, label="RRTMG T_strat ≈ 190 K")
ax_T.grid(True, alpha=0.3); ax_T.legend(fontsize=9)

# BT bar chart per band, grouped by boost factor.
nband = band_edges.shape[0]
width = 0.8 / len(LABELS)
x = np.arange(nband)
for i, (lbl, c) in enumerate(zip(LABELS, colors)):
    snap = history[lbl][-1]
    k_diag = int(np.argmin(np.abs(snap["p_int_hpa"] - DIAG_PRESSURE_HPA)))
    F_diag = snap["up_lw_per_band"][k_diag, :]
    BT = [band_brightness_T(F_diag[b], band_edges[b, 0], band_edges[b, 1])
          for b in range(nband)]
    ax_BT.bar(x + i * width - 0.4, BT, width, color=c, label=lbl)
ax_BT.set_xticks(x)
ax_BT.set_xticklabels([f"b{b}\n{band_edges[b,0]:.0f}–{band_edges[b,1]:.0f}"
                       for b in range(nband)], fontsize=8)
ax_BT.set_ylabel(f"Brightness T at {DIAG_PRESSURE_HPA:.0f} hPa (K)")
ax_BT.set_title("Per-band BT — wings (b1, b3) should drop if boost works")
ax_BT.legend(fontsize=9)
ax_BT.grid(True, alpha=0.3, axis="y")

fig.suptitle("CO₂-wing opacity boost: does T_strat warm toward RRTMG?", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(os.path.dirname(__file__), "experiment_co2_wing_boost.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved {out}")
