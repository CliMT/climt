"""Radiative equilibrium: RRTMG LW vs. Cork LW (2 / 4 / 8 g-points).

Each column is coupled to a SlabSurface and forced by a constant absorbed
shortwave flux.  At true equilibrium, OLR balances the absorbed SW, so the
TOA imbalance is the natural convergence diagnostic.

Shows g-point convergence: as g-points increase from 2→4→8, the Cork
correlated-k solution should converge toward RRTMG longwave.

CO2 is set to 376 ppm for RRTMG — matching the value baked into the Earth
correlated-k tables (Chaverot dataset).  Specific humidity is zero in all
models (dry, clear-sky comparison).

Run with:
    conda run -n climt python examples/radiative_equilibrium_cork_convergence.py
"""
import os
from datetime import timedelta

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from sympl import AdamsBashforth
from tqdm import trange

from climt import (
    RRTMGLongwave,
    SlabSurface,
    get_default_state,
    get_grid,
)
from climt._components.cork import CorkLongwaveRadiation

# Physical constants for band-integrated Planck brightness temperature.
H_PLANCK = 6.62607015e-34   # J s
C_LIGHT  = 2.99792458e8     # m s^-1
K_BOLTZ  = 1.380649e-23     # J K^-1
SIGMA_SB = 5.670374419e-8   # W m^-2 K^-4

def _planck_band_flux(T, nu1_cm, nu2_cm, n_sub=40):
    """π ∫_{ν1}^{ν2} B_ν(T) dν, with ν in cm^-1 → flux in W/m². Vectorised in T."""
    nu = np.linspace(nu1_cm, nu2_cm, n_sub) * 100.0    # cm^-1 → m^-1
    T = np.atleast_1d(T).astype(float)
    # B_ν in [W m^-2 sr^-1 / m^-1]; integrand × π gives flux density per ν.
    x = (H_PLANCK * C_LIGHT * nu[:, None]) / (K_BOLTZ * np.maximum(T[None, :], 1.0))
    B = (2.0 * H_PLANCK * C_LIGHT**2 * nu[:, None]**3) / (np.expm1(x))
    # Integrate over ν (in m^-1) then × π for hemispheric flux:
    flux = np.pi * np.trapz(B, nu, axis=0)
    return flux if flux.size > 1 else flux[0]

def band_brightness_T(F_band, nu1_cm, nu2_cm):
    """Invert F_band = π ∫ B_ν(T) dν over [ν1, ν2] (cm^-1) for T (K)."""
    if F_band <= 0:
        return np.nan
    def f(T):
        return _planck_band_flux(T, nu1_cm, nu2_cm) - F_band
    try:
        return brentq(f, 50.0, 400.0)
    except ValueError:
        return np.nan

def load_band_edges(table_name):
    """Return (nband, 2) array of [ν_lower, ν_upper] in cm^-1 for a CORK k-table."""
    path = os.path.join(
        os.path.dirname(__file__), "..", "..", "climt", "_data",
        "cork", "correlated_k", f"{table_name}.nc",
    )
    with xr.open_dataset(path) as ds:
        return ds["band_wavenumber_limits"].values  # (nband, 2), cm^-1

# ── Parameters ────────────────────────────────────────────────────────────────
CO2_MOLE_FRAC = 376e-6   # Match CO2 baked into earth_low_res_lw tables (Chaverot)
SW_FLUX = 240.0          # Absorbed SW at the surface (W m⁻²)
MIXED_LAYER_DEPTH = 1.0  # m — short radiative timescale so equilibrium is fast
TIMESTEP = timedelta(hours=3)
import argparse as _argparse  # local-only; keeps top imports stable
_ap = _argparse.ArgumentParser(add_help=False)
_ap.add_argument("--days", type=float, default=1000.0,
                 help="simulated days (default 1000)")
_args, _ = _ap.parse_known_args()
N_STEPS = max(1, int(round(_args.days * 86400.0 / TIMESTEP.total_seconds())))
STORE_EVERY = 100        # finer sampling so we can see TOA balance converge
DIAG_PRESSURE_HPA = 100.0  # stratospheric reference level for per-band BT diag

NZ = 30

# Skin temperature: T_eff / 2^(1/4) where σ T_eff^4 = SW_FLUX.
# Initialize the whole column isothermally at T_skin so the upper levels start
# near their expected radiative-equilibrium value and convergence is bounded
# by the dynamic range each model's true equilibrium differs from skin temp.
SIGMA_SB = 5.670374419e-8
T_EFF = (SW_FLUX / SIGMA_SB) ** 0.25
T_SKIN = T_EFF / (2.0 ** 0.25)

LABELS = [
    "RRTMG (16 bands)",
    "CORK 4-band GL (shipped)",
    "CORK 6-band Chaverot (defective)",
    "CORK 6-band linepyline (fix)",
]
# Table name for each CORK label (None for RRTMG — no per-band diagnostic).
CORK_TABLE = {
    "RRTMG (16 bands)":              None,
    "CORK 4-band GL (shipped)":        "earth_low_res_lw_8gpt",
    "CORK 6-band Chaverot (defective)": "earth_low_res_lw_co2refined_gl",
    "CORK 6-band linepyline (fix)":    "earth_low_res_lw_co2refined_linepyline",
}
COLORS = ["black", "#999999", "#ff7f0e", "#d62728"]
STYLES = ["-", "--", "-", "-"]
# Diagnostic key for LW heating rate (K day⁻¹) differs between RRTMG and CORK
HR_KEY = {label: "air_temperature_tendency_from_longwave" for label in LABELS}

# ── Grid ──────────────────────────────────────────────────────────────────────
grid = get_grid(nx=1, ny=1, nz=NZ)

# ── RRTMG LW ──────────────────────────────────────────────────────────────────
# Must be created before any CorkLongwaveRadiation so that get_default_state uses
# 16 LW bands for this state; CORK constructors reset the global to 4 bands.
rad_lw = RRTMGLongwave()
surface_rrtmg = SlabSurface()
rrtmg_state = get_default_state([rad_lw, surface_rrtmg], grid_state=grid)
rrtmg_state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_MOLE_FRAC
rrtmg_state["specific_humidity"].values[:] = 0.0
# RRTMG's get_default_state populates a climatological ozone profile (max ~3.9
# ppm in the stratosphere). The Cork CO2(+H2O) tables have NO ozone, so
# leaving RRTMG's O3 on makes the "head-to-head" compare RRTMG-with-O3 vs
# CORK-without — O3's 9.6 um LW band is a major stratospheric absorber. Zero it so
# both sides see the same absorbers (CO2 only, dry).
rrtmg_state["mole_fraction_of_ozone_in_air"].values[:] = 0.0

# ── Cork LW — three table variants ───────────────────────────────────
# Each CORK constructor resets the global num_longwave_bands, which get_default_state
# reads when sizing the emissivity field. Build each state immediately after its
# constructor so the band count matches.
cork_shipped = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw_8gpt")
surface_cork_shipped = SlabSurface()
cork_state_shipped = get_default_state([cork_shipped, surface_cork_shipped], grid_state=grid)

cork_6b_chav = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw_co2refined_gl")
surface_cork_6b_chav = SlabSurface()
cork_state_6b_chav = get_default_state([cork_6b_chav, surface_cork_6b_chav], grid_state=grid)

cork_6b_lpl  = CorkLongwaveRadiation(optics="correlated_k", table="earth_low_res_lw_co2refined_linepyline")
surface_cork_6b_lpl  = SlabSurface()
cork_state_6b_lpl  = get_default_state([cork_6b_lpl,  surface_cork_6b_lpl],  grid_state=grid)

all_cork_states = [cork_state_shipped, cork_state_6b_chav, cork_state_6b_lpl]
for st in all_cork_states:
    st["specific_humidity"].values[:] = 0.0

# Constant SW forcing and shallow slab on every column.
# Initialize the column isothermally at the skin temperature so the upper
# levels don't spend the run drifting from an arbitrary default initial T.
for st in [rrtmg_state] + all_cork_states:
    st["downwelling_shortwave_flux_in_air"].values[:] = SW_FLUX
    st["ocean_mixed_layer_thickness"].values[:] = MIXED_LAYER_DEPTH
    st["area_type"].values[:] = "sea"
    st["air_temperature"].values[:] = T_SKIN
    st["surface_temperature"].values[:] = T_SKIN

print(f"Isothermal init: T_eff = {T_EFF:.2f} K, T_skin = {T_SKIN:.2f} K")

# ── Sympl time steppers (LW radiation + slab surface) ─────────────────────────
steppers = {
    "RRTMG (16 bands)":               AdamsBashforth([rad_lw,     surface_rrtmg]),
    "CORK 4-band GL (shipped)":         AdamsBashforth([cork_shipped, surface_cork_shipped]),
    "CORK 6-band Chaverot (defective)": AdamsBashforth([cork_6b_chav, surface_cork_6b_chav]),
    "CORK 6-band linepyline (fix)":     AdamsBashforth([cork_6b_lpl,  surface_cork_6b_lpl]),
}
states = {
    "RRTMG (16 bands)":               rrtmg_state,
    "CORK 4-band GL (shipped)":         cork_state_shipped,
    "CORK 6-band Chaverot (defective)": cork_state_6b_chav,
    "CORK 6-band linepyline (fix)":     cork_state_6b_lpl,
}

# ── Time-stepping loop ────────────────────────────────────────────────────────
# history[label] = list of {step, T_profile, hr_day, olr, T_surface, imbalance, ...}
history = {label: [] for label in LABELS}

for step in trange(N_STEPS, desc="Integrating"):
    for label in LABELS:
        diag, new_state = steppers[label](states[label], TIMESTEP)
        new_state.update(diag)
        states[label] = new_state

        if step % STORE_EVERY == 0:
            up_lw = diag["upwelling_longwave_flux_in_air"].values[:, 0, 0].copy()
            olr = float(up_lw[-1])
            entry = {
                "step": step,
                "p_hpa": new_state["air_pressure"].values[:, 0, 0].copy() / 100.0,
                "p_int_hpa": new_state["air_pressure_on_interface_levels"].values[:, 0, 0].copy() / 100.0,
                "T_profile": new_state["air_temperature"].values[:, 0, 0].copy(),
                "hr_day": diag[HR_KEY[label]].values[:, 0, 0].copy(),
                "up_lw": up_lw,
                "olr": olr,
                "toa_imbalance": SW_FLUX - olr,
                "T_surface": float(new_state["surface_temperature"].values[0, 0]),
            }
            # CORK variants expose per-band upwelling fluxes; capture them so we
            # can compute κ-weighted brightness temperatures aloft.
            if "upwelling_longwave_flux_in_air_per_band" in diag:
                # shape: (interface_levels, x, y, num_longwave_bands)
                entry["up_lw_per_band"] = (
                    diag["upwelling_longwave_flux_in_air_per_band"]
                    .values[:, 0, 0, :].copy()
                )
            history[label].append(entry)

# ── Final report: TOA balance ─────────────────────────────────────────────────
print()
print(f"TOA balance (target: OLR = {SW_FLUX:.1f} W m⁻²):")
for label in LABELS:
    snap = history[label][-1]
    print(
        f"  {label:18s}  OLR = {snap['olr']:7.2f} W m⁻²"
        f"  imbalance = {snap['toa_imbalance']:+6.3f} W m⁻²"
        f"  T_sfc = {snap['T_surface']:6.2f} K"
    )

# ── Per-level RRTMG vs CORK diagnostic ──────────────────────────────────────────
# What we want: for each CORK variant, the per-level temperature and LW heating
# rate compared to RRTMG at the same level, so we can localize the stratospheric
# cooling-rate gap that produces the cold T bias.
RRTMG_LABEL = "RRTMG (16 bands)"
rr_snap = history[RRTMG_LABEL][-1]
T_rr = np.asarray(rr_snap["T_profile"])
HR_rr = np.asarray(rr_snap["hr_day"])
p_hpa = np.asarray(rr_snap["p_hpa"])

print()
print(f"Per-level RRTMG vs CORK (final snapshot, sim day "
      f"{N_STEPS * TIMESTEP.total_seconds() / 86400.0:.0f}):")
print(f"  RRTMG OLR = {rr_snap['olr']:.2f} W m^-2, "
      f"T_sfc = {rr_snap['T_surface']:.2f} K, "
      f"T_strat_min = {T_rr.min():.2f} K @ p = {p_hpa[np.argmin(T_rr)]:.2f} hPa")

for label in LABELS:
    if label == RRTMG_LABEL:
        continue
    cork_snap = history[label][-1]
    T_cork = np.asarray(cork_snap["T_profile"])
    HR_cork = np.asarray(cork_snap["hr_day"])
    dT = T_cork - T_rr
    dHR = HR_cork - HR_rr  # CORK − RRTMG; negative = CORK cools more
    print()
    print(f"--- {label} ---")
    print(f"  OLR = {cork_snap['olr']:.2f} W m^-2, "
          f"T_sfc = {cork_snap['T_surface']:.2f} K, "
          f"T_strat_min = {T_cork.min():.2f} K @ p = {p_hpa[np.argmin(T_cork)]:.2f} hPa")
    print(f"  max |ΔT| = {np.abs(dT).max():.2f} K @ p = "
          f"{p_hpa[np.argmax(np.abs(dT))]:.2f} hPa")
    print(f"  max |ΔHR| = {np.abs(dHR).max()*1000:.3f} mK/day @ p = "
          f"{p_hpa[np.argmax(np.abs(dHR))]:.2f} hPa  "
          f"(CORK − RRTMG; negative = CORK cools more)")
    print(f"  {'p [hPa]':>9s}  {'T_RR':>7s}  {'T_CORK':>7s}  {'ΔT':>7s}  "
          f"{'HR_RR':>9s}  {'HR_CORK':>9s}  {'ΔHR':>9s}    units: K, K/day")
    for k in range(len(p_hpa)):
        print(f"  {p_hpa[k]:9.3f}  "
              f"{T_rr[k]:7.2f}  {T_cork[k]:7.2f}  {dT[k]:+7.2f}  "
              f"{HR_rr[k]:9.4f}  {HR_cork[k]:9.4f}  {dHR[k]:+9.4f}")

# Save the history for later analysis without re-running.
_npz_path = os.path.join(os.path.dirname(__file__),
                         "radiative_equilibrium_cork_convergence_history.npz")
_to_save = {}
for label in LABELS:
    snap = history[label][-1]
    safe = label.replace(" ", "_").replace("(", "").replace(")", "")
    for k in ("T_profile", "hr_day", "up_lw", "p_hpa", "p_int_hpa"):
        if k in snap:
            _to_save[f"{safe}__{k}"] = np.asarray(snap[k])
np.savez(_npz_path, **_to_save)
print(f"\nSaved final snapshots to {_npz_path}")

# ── Final plot ────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
ax_T, ax_HR, ax_UP = axes[0]
ax_OLR, ax_IMB, ax_TTOP = axes[1]

# Index of the top model level (TOA-most layer mid-point). Used for the
# stratospheric-T-vs-time convergence diagnostic.
N_TOP_LEVELS = 3  # plot the top 3 levels

for label, color, ls in zip(LABELS, COLORS, STYLES):
    snap = history[label][-1]
    all_snaps = history[label]
    times_days = [s["step"] * TIMESTEP.total_seconds() / 86400.0 for s in all_snaps]

    ax_T.plot(snap["T_profile"], snap["p_hpa"], color=color, ls=ls, lw=2, label=label)
    ax_HR.plot(snap["hr_day"], snap["p_hpa"], color=color, ls=ls, lw=2, label=label)
    ax_UP.plot(snap["up_lw"], snap["p_int_hpa"], color=color, ls=ls, lw=2, label=label)
    ax_OLR.plot(times_days, [s["olr"] for s in all_snaps],
                color=color, ls=ls, lw=2, label=label)
    ax_IMB.plot(times_days, [s["toa_imbalance"] for s in all_snaps],
                color=color, ls=ls, lw=2, label=label)

    # Top-of-atmosphere temperature time series — convergence diagnostic.
    # If these are still drifting at end of run, equilibrium not yet reached.
    p_top = all_snaps[-1]["p_hpa"][-N_TOP_LEVELS:]
    for k in range(N_TOP_LEVELS):
        T_top_t = [s["T_profile"][-(k + 1)] for s in all_snaps]
        # Solid for top level, fading transparency for levels just below
        alpha = 1.0 - 0.3 * k
        ax_TTOP.plot(times_days, T_top_t, color=color, ls=ls, lw=1.5,
                     alpha=alpha,
                     label=f"{label} (p={p_top[-(k+1)]:.2f} hPa)" if k == 0 else None)

for ax in (ax_T, ax_HR, ax_UP):
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel("Pressure (hPa)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

ax_UP.set_xlabel("Upwelling LW flux (W m⁻²)")
ax_UP.set_title("Upwelling LW (radiating height)")
ax_UP.axvline(SW_FLUX, color="grey", lw=1, ls=":")

ax_T.set_xlabel("Temperature (K)")
ax_T.set_title("Temperature (final)")

ax_HR.set_xlabel("K day⁻¹")
ax_HR.set_title("LW heating rate (final)")

ax_OLR.axhline(SW_FLUX, color="grey", lw=1, ls=":", label=f"Absorbed SW = {SW_FLUX:.0f} W m⁻²")
ax_OLR.set_xlabel("Time (days)")
ax_OLR.set_ylabel("OLR (W m⁻²)")
ax_OLR.set_title("OLR convergence")
ax_OLR.legend(fontsize=8)
ax_OLR.grid(True, alpha=0.3)

ax_IMB.axhline(0.0, color="grey", lw=1, ls=":")
ax_IMB.set_xlabel("Time (days)")
ax_IMB.set_ylabel("SW − OLR (W m⁻²)")
ax_IMB.set_title("TOA energy imbalance")
ax_IMB.legend(fontsize=8)
ax_IMB.grid(True, alpha=0.3)

ax_TTOP.axhline(T_SKIN, color="grey", lw=1, ls=":",
                label=f"T_skin = {T_SKIN:.1f} K")
ax_TTOP.set_xlabel("Time (days)")
ax_TTOP.set_ylabel("Temperature (K)")
ax_TTOP.set_title(f"Top-{N_TOP_LEVELS}-level T(t)  (drift = not converged)")
ax_TTOP.legend(fontsize=7, loc="best")
ax_TTOP.grid(True, alpha=0.3)

fig.suptitle(
    f"Radiative equilibrium (dry, clear-sky, slab surface) — "
    f"RRTMG vs. Cork LW\n"
    f"Earth, {int(CO2_MOLE_FRAC * 1e6)} ppm CO₂, "
    f"SW = {SW_FLUX:.0f} W m⁻², "
    f"{N_STEPS} × {TIMESTEP.total_seconds() / 3600:.0f} h steps",
    fontsize=11,
)
plt.tight_layout(rect=[0, 0, 1, 0.95])

out_path = os.path.join(os.path.dirname(__file__), "radiative_equilibrium_cork_convergence.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"\nSaved to {out_path}")

# ── Per-band diagnostic: where does each band radiate from? ────────────────────
# Pure-LW rad-eq stratospheric T is set by the κ-weighted brightness temperature
# of the upwelling field. For each CORK variant, convert the per-band upwelling
# flux at DIAG_PRESSURE_HPA to a brightness temperature using the band's
# wavenumber range. Bands carrying little flux at this level are the ones that
# saturate low in the atmosphere — those control T_strat in equilibrium.
cork_labels = [lbl for lbl in LABELS if CORK_TABLE[lbl] is not None]
n_cork = len(cork_labels)

fig2, axes2 = plt.subplots(2, n_cork, figsize=(4 * n_cork, 9))
if n_cork == 1:
    axes2 = axes2[:, None]

print(f"\nPer-band upwelling at p ≈ {DIAG_PRESSURE_HPA:.0f} hPa "
      f"(brightness temperature inverted from band-integrated Planck):")

for j, label in enumerate(cork_labels):
    snap = history[label][-1]
    band_edges = load_band_edges(CORK_TABLE[label])  # (nband, 2) cm^-1
    nband = band_edges.shape[0]

    p_int = snap["p_int_hpa"]
    up_band = snap["up_lw_per_band"]          # (nlev+1, nband)
    p_mid = snap["p_hpa"]
    T_prof = snap["T_profile"]

    # Find the interface index closest to DIAG_PRESSURE_HPA.
    k_diag = int(np.argmin(np.abs(p_int - DIAG_PRESSURE_HPA)))
    p_actual = p_int[k_diag]
    F_diag = up_band[k_diag, :]               # (nband,)

    # Local layer temperature at the closest mid-level — for comparison with
    # κ-weighted brightness temperature target.
    k_mid = int(np.argmin(np.abs(p_mid - DIAG_PRESSURE_HPA)))
    T_local = T_prof[k_mid]

    # Per-band brightness temperatures at the diagnostic level.
    BT_band = np.array([
        band_brightness_T(F_diag[b], band_edges[b, 0], band_edges[b, 1])
        for b in range(nband)
    ])

    # Top row: per-band upwelling flux profile.
    ax_up = axes2[0, j]
    cmap = plt.get_cmap("turbo", nband)
    for b in range(nband):
        ax_up.plot(up_band[:, b], p_int,
                   color=cmap(b), lw=1.5,
                   label=f"b{b}: {band_edges[b,0]:.0f}-{band_edges[b,1]:.0f}")
    ax_up.set_yscale("log")
    ax_up.invert_yaxis()
    ax_up.axhline(p_actual, color="k", ls=":", lw=0.8)
    ax_up.set_xlabel("Upwelling LW (W m⁻²)")
    ax_up.set_ylabel("Pressure (hPa)")
    ax_up.set_title(label)
    ax_up.grid(True, alpha=0.3)
    ax_up.legend(fontsize=6, loc="best")

    # Bottom row: per-band brightness T at diagnostic level (bar chart).
    ax_bt = axes2[1, j]
    nu_mid = 0.5 * (band_edges[:, 0] + band_edges[:, 1])
    bars = ax_bt.bar(np.arange(nband), BT_band, color=[cmap(b) for b in range(nband)])
    ax_bt.axhline(T_local, color="red", ls="--", lw=1.5,
                  label=f"Local T = {T_local:.1f} K")
    ax_bt.set_xticks(np.arange(nband))
    ax_bt.set_xticklabels([f"{nu_mid[b]:.0f}" for b in range(nband)],
                          fontsize=7, rotation=45)
    ax_bt.set_xlabel("Band ν center (cm⁻¹)")
    ax_bt.set_ylabel(f"Brightness T at {p_actual:.1f} hPa (K)")
    ax_bt.set_title(f"Flux total = {F_diag.sum():.1f} W m⁻²")
    ax_bt.legend(fontsize=7)
    ax_bt.grid(True, alpha=0.3, axis="y")

    # Console summary.
    print(f"  {label}  (p_actual = {p_actual:.2f} hPa, T_local = {T_local:.2f} K):")
    for b in range(nband):
        print(f"    band {b} [{band_edges[b,0]:6.0f}-{band_edges[b,1]:6.0f} cm⁻¹]  "
              f"F = {F_diag[b]:7.3f} W m⁻²   BT = {BT_band[b]:6.2f} K")

fig2.suptitle(
    f"Per-band upwelling LW & brightness temperature at p ≈ {DIAG_PRESSURE_HPA:.0f} hPa\n"
    f"Low-BT bands are the ones constraining stratospheric T (κ-weighted Planck balance)",
    fontsize=11,
)
plt.tight_layout(rect=[0, 0, 1, 0.94])
out_path2 = os.path.join(os.path.dirname(__file__), "radiative_equilibrium_cork_per_band.png")
fig2.savefig(out_path2, dpi=150, bbox_inches="tight")
print(f"\nSaved per-band diagnostic to {out_path2}")
