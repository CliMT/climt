"""Moist radiative-convective equilibrium: RRTMG-LW vs Cork-LW tropopause.

Third rung of the model hierarchy (above pure radiative equilibrium and the dry
RCE in `rce_dry_cork_vs_rrtmg.py`). Each single column is driven by a prescribed
absorbed SW flux and relaxed to equilibrium under:

    LW radiation (RRTMG  or  CorkLongwaveRadiation)
  + SlabSurface
  + SimplePhysics(boundary_layer=True)    -- surface fluxes + BL diffusion
  + EmanuelConvection                     -- moist convection

The column is MOIST: specific humidity is allowed to spin up from surface
evaporation (SimplePhysics latent flux). Moist convection sets the tropospheric
lapse rate toward a moist adiabat; the tropopause emerges where convection stops
and radiative cooling balances. Diagnostic of interest: does CORK-LW place the
tropopause at a similar height/temperature as RRTMG-LW under moist convection?

Apples-to-apples LW test (per user choices):
  * Shortwave is a fixed prescribed flux, identical on every column.
  * Ozone is zeroed on RRTMG (CORK has none) -> both see CO2(+H2O) only.

The DRY counterpart (SimplePhysics no-BL + DryConvectiveAdjustment) is the
separate `rce_dry_cork_vs_rrtmg.py`; the two give different tropopauses and are
meant to be compared.

Run:  conda run --no-capture-output -n climt python scripts/experiments/rce_moist_cork_vs_rrtmg.py [--days N]
"""
import argparse
import os
from datetime import timedelta

import numpy as np
import matplotlib.pyplot as plt
from sympl import AdamsBashforth
from tqdm import trange

from climt import (
    RRTMGLongwave, SlabSurface, SimplePhysics, EmanuelConvection,
    get_default_state, get_grid,
)
from climt._components.cork import CorkLongwaveRadiation

# ── Parameters ─────────────────────────────────────────────────────────────
CO2_MOLE_FRAC = 376e-6
SW_FLUX = 240.0
MIXED_LAYER_DEPTH = 5.0      # m
TIMESTEP = timedelta(minutes=5)   # Emanuel is stiff -> shorter step
NZ = 30
T_INIT = 270.0
T_SFC_INIT = 285.0

ap = argparse.ArgumentParser()
ap.add_argument("--days", type=float, default=200.0, help="simulated days")
ap.add_argument("--co2", type=float, default=None,
                help="CO2 in ppm; overrides the default 376 ppm")
ap.add_argument("--table", default=None,
                help="single CORK table name to compare vs RRTMG; if given, "
                     "replaces the COLUMNS dict with {RRTMG, that table}")
ap.add_argument("--save", default=None,
                help="if given, write final-state arrays to this .npz "
                     "(keys: <label>__{p_hpa,T,q,hr_lw,T_sfc,olr})")
args = ap.parse_args()
if args.co2 is not None:
    CO2_MOLE_FRAC = args.co2 * 1e-6
N_STEPS = max(1, int(round(args.days * 86400.0 / TIMESTEP.total_seconds())))
STORE_EVERY = max(1, N_STEPS // 400)

COLUMNS = {
    "RRTMG-LW":               None,
    "CORK-LW shipped default":  "earth_low_res_lw_8gpt",
    "CORK-LW candD decoupcont": "earth_lw_candD_decoupcont",
}
COLORS = {"RRTMG-LW": "black", "CORK-LW shipped default": "#999999",
          "CORK-LW candD decoupcont": "#d62728"}
if args.table is not None:
    COLUMNS = {"RRTMG-LW": None, f"CORK-LW {args.table}": args.table}

grid = get_grid(nx=1, ny=1, nz=NZ)


def make_lw(table):
    return RRTMGLongwave() if table is None else \
        CorkLongwaveRadiation(optics="correlated_k", table=table)


def build_column(table):
    """One moist-RCE column. LW comp is built here so the CORK band-count global is
    set immediately before get_default_state (RRTMG column must be built first)."""
    lw = make_lw(table)
    slab = SlabSurface()
    simple_physics = SimplePhysics(boundary_layer=True)
    convection = EmanuelConvection()
    convection.current_time_step = TIMESTEP
    state = get_default_state([lw, slab, simple_physics, convection], grid_state=grid)

    state["downwelling_shortwave_flux_in_air"].values[:] = SW_FLUX
    state["ocean_mixed_layer_thickness"].values[:] = MIXED_LAYER_DEPTH
    state["area_type"].values[:] = "sea"
    state["air_temperature"].values[:] = T_INIT
    state["surface_temperature"].values[:] = T_SFC_INIT
    # moist: leave specific_humidity at its default; it spins up via evaporation
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_MOLE_FRAC
    if "mole_fraction_of_ozone_in_air" in state:
        state["mole_fraction_of_ozone_in_air"].values[:] = 0.0

    stepper = AdamsBashforth([convection, lw, slab])
    return {"stepper": stepper, "sp": simple_physics, "conv": convection,
            "state": state}


models = {}
for label, table in COLUMNS.items():
    models[label] = build_column(table)

# ── Time-stepping (tendencies incl. Emanuel -> simple physics w/ BL) ─────────
history = {label: [] for label in COLUMNS}
for step in trange(N_STEPS, desc="moist RCE"):
    for label in COLUMNS:
        m = models[label]
        st = m["state"]
        m["conv"].current_time_step = TIMESTEP

        diag, st = m["stepper"](st, TIMESTEP)    # Emanuel + LW + slab
        st.update(diag)

        diag, new = m["sp"](st, TIMESTEP)         # SimplePhysics (with BL)
        st.update(diag)
        st.update(new)

        st["eastward_wind"].values[:] = 3.0
        m["state"] = st

        if step % STORE_EVERY == 0 or step == N_STEPS - 1:
            up_lw = st["upwelling_longwave_flux_in_air"].values[:, 0, 0]
            entry = {
                "step": step,
                "p_hpa": st["air_pressure"].values[:, 0, 0].copy() / 100.0,
                "T": st["air_temperature"].values[:, 0, 0].copy(),
                "q": st["specific_humidity"].values[:, 0, 0].copy(),
                "hr_lw": st["air_temperature_tendency_from_longwave"].values[:, 0, 0].copy(),
                "olr": float(up_lw[-1]),
                "T_sfc": float(st["surface_temperature"].values[0, 0]),
            }
            if "air_temperature_tendency_from_convection" in st:
                # this diagnostic has dims ('*', 'mid_levels') -> flatten for the
                # single column (nx=ny=1) to align with the (mid_levels,) profile
                entry["hr_conv"] = st["air_temperature_tendency_from_convection"].values.reshape(-1).copy()
            history[label].append(entry)


def tropopause(p_hpa, T, dtheta=10.0):
    """Top of the convective (near-constant-theta) layer.

    A cold-point definition is meaningless here: with O3=0 and no stratospheric
    heating, T decreases monotonically to the model top, so argmin(T) just
    returns the top level. Instead mark the tropopause by potential temperature:
    scan up from the surface to the first level where theta exceeds the surface
    value by `dtheta` K — the base of the sharp stratospheric theta rise.

    Moist convection follows a moist adiabat, so theta drifts a few K through the
    troposphere (unlike the dry, well-mixed +2 K case in rce_dry_cork_vs_rrtmg.py);
    +10 K cleanly clears that drift and lands at the stratospheric cap.
    """
    theta = T * (1000.0 / p_hpa) ** (287.0 / 1004.0)
    above = np.where(theta > theta[0] + dtheta)[0]   # index 0 = surface
    k = int(above[0]) if above.size else len(p_hpa) - 1
    return p_hpa[k], T[k]


# ── Report ─────────────────────────────────────────────────────────────────
print(f"\nMoist RCE — {args.days:.0f} days, SW={SW_FLUX:.0f} W/m^2, O3=0\n")
print(f"{'column':22s} {'OLR':>7s} {'T_sfc':>7s} {'p_tp':>8s} {'T_tp':>7s} {'q_sfc[g/kg]':>11s}")
for label in COLUMNS:
    s = history[label][-1]
    p_tp, T_tp = tropopause(s["p_hpa"], s["T"])
    print(f"{label:22s} {s['olr']:7.2f} {s['T_sfc']:7.2f} {p_tp:8.2f} {T_tp:7.2f} {s['q'][0]*1e3:11.3f}")

# ── Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 7))
axT, axHR, axQ = axes
for label in COLUMNS:
    s = history[label][-1]
    c = COLORS.get(label, "#1f77b4")
    axT.plot(s["T"], s["p_hpa"], "-o", ms=3, color=c, label=label)
    axHR.plot(s["hr_lw"], s["p_hpa"], "-", color=c, label=f"{label} LW")
    if "hr_conv" in s:
        axHR.plot(s["hr_conv"], s["p_hpa"], "--", color=c, alpha=0.6)
    axQ.plot(s["q"] * 1e3, s["p_hpa"], "-", color=c, label=label)
    p_tp, T_tp = tropopause(s["p_hpa"], s["T"])
    axT.plot(T_tp, p_tp, "*", color=c, ms=16, mec="k")
for ax in axes:
    ax.set_yscale("log"); ax.invert_yaxis(); ax.grid(alpha=0.3)
    ax.set_ylabel("Pressure (hPa)"); ax.legend(fontsize=8)
axT.set_xlabel("Temperature (K)")
axT.set_title("Moist RCE temperature (★ = convective θ-top tropopause)")
axHR.set_xlabel("Heating rate (K/day)  [solid=LW, dashed=conv]")
axHR.set_title("LW & convective heating")
axQ.set_xlabel("Specific humidity (g/kg)")
axQ.set_title("Moisture")
fig.suptitle(f"Moist RCE: RRTMG-LW vs CORK-LW — SimplePhysics(BL) + EmanuelConvection\n"
             f"CO2={int(CO2_MOLE_FRAC*1e6)} ppm, O3=0, SW={SW_FLUX:.0f} W/m^2, {args.days:.0f} days")
plt.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "debug_data", "rce_moist_cork_vs_rrtmg.png"))
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved {out}")

npz = os.path.join(os.path.dirname(out), "rce_moist_cork_vs_rrtmg.npz")
np.savez(npz, **{f"{l.replace(' ','_')}__{k}": np.asarray(history[l][-1][k])
                 for l in COLUMNS for k in ("p_hpa", "T", "q", "hr_lw")})
print(f"Saved {npz}")

if args.save:
    payload = {}
    for label in COLUMNS:
        s = history[label][-1]
        safe = label.replace(" ", "_")
        payload[f"{safe}__p_hpa"] = np.asarray(s["p_hpa"])
        payload[f"{safe}__T"] = np.asarray(s["T"])
        payload[f"{safe}__q"] = np.asarray(s["q"])
        payload[f"{safe}__hr_lw"] = np.asarray(s["hr_lw"])
        payload[f"{safe}__T_sfc"] = np.asarray(s["T_sfc"])
        payload[f"{safe}__olr"] = np.asarray(s["olr"])
    np.savez(args.save, **payload)
    print(f"saved final state -> {args.save}")
