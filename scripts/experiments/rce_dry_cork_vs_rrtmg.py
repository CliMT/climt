"""Dry radiative-convective equilibrium: RRTMG-LW vs Cork-LW tropopause.

Second rung of the model hierarchy (above pure radiative equilibrium in
`radiative_equilibrium_cork_convergence.py`). Each single column is driven by a
prescribed absorbed SW flux at the surface and relaxed to equilibrium under:

    LW radiation (RRTMG  or  CorkLongwaveRadiation)
  + SlabSurface
  + SimplePhysics(boundary_layer=False)   -- surface fluxes, no BL diffusion
  + DryConvectiveAdjustment               -- keeps the troposphere dry-adiabatic

The column is DRY (q = 0). With dry convective adjustment the troposphere
follows the dry adiabat and a tropopause emerges where radiative cooling takes
over. The diagnostic of interest is the tropopause: does CORK-LW place it at a
similar height/temperature as RRTMG-LW?

Apples-to-apples LW test (per user choices):
  * Shortwave is a fixed prescribed flux, identical on every column -> only the
    LW scheme differs.
  * Ozone is zeroed on RRTMG (CORK tables have none) -> both see CO2(+H2O) only.

Run:  conda run --no-capture-output -n climt python scripts/experiments/rce_dry_cork_vs_rrtmg.py [--days N]
"""
import argparse
import os
from datetime import timedelta

import numpy as np
import matplotlib.pyplot as plt
from sympl import AdamsBashforth
from tqdm import trange

from climt import (
    RRTMGLongwave, SlabSurface, SimplePhysics, DryConvectiveAdjustment,
    get_default_state, get_grid,
)
from climt._components.cork import CorkLongwaveRadiation

# ── Parameters ─────────────────────────────────────────────────────────────
CO2_MOLE_FRAC = 376e-6
SW_FLUX = 240.0
MIXED_LAYER_DEPTH = 1.0      # m — thin slab so the surface equilibrates fast
TIMESTEP = timedelta(minutes=10)
NZ = 30
SIGMA_SB = 5.670374419e-8
T_INIT = 260.0               # isothermal start

ap = argparse.ArgumentParser()
ap.add_argument("--days", type=float, default=400.0, help="simulated days")
ap.add_argument("--co2", type=float, default=None,
                help="CO2 in ppm; overrides the default 376 ppm")
ap.add_argument("--table", default=None,
                help="single CORK table name to compare vs RRTMG; if given, "
                     "replaces the COLUMNS dict with {RRTMG, that table}")
ap.add_argument("--save", default=None,
                help="if given, write final-state arrays to this .npz "
                     "(keys: <label>__{p_hpa,T,q_max,hr_lw,T_sfc,olr})")
args = ap.parse_args()
if args.co2 is not None:
    CO2_MOLE_FRAC = args.co2 * 1e-6
N_STEPS = max(1, int(round(args.days * 86400.0 / TIMESTEP.total_seconds())))
STORE_EVERY = max(1, N_STEPS // 400)

# label -> CORK table name (None = RRTMG)
COLUMNS = {
    "RRTMG-LW":            None,
    "CORK-LW candB 12band":  "earth_lw_candB",
}
COLORS = {"RRTMG-LW": "black", "CORK-LW candB 12band": "#d62728"}
if args.table is not None:
    COLUMNS = {"RRTMG-LW": None, f"CORK-LW {args.table}": args.table}

grid = get_grid(nx=1, ny=1, nz=NZ)


def make_lw(table):
    return RRTMGLongwave() if table is None else \
        CorkLongwaveRadiation(optics="correlated_k", table=table)


def build_column(table):
    """One dry-RCE column. LW comp is built here so the CORK band-count global is
    set immediately before get_default_state (RRTMG column must be built first)."""
    lw = make_lw(table)
    slab = SlabSurface()
    # use_external_surface_specific_humidity=True -> SimplePhysics reads the
    # surface humidity from state instead of computing saturation. Held at 0 so
    # there is NO evaporation (truly dry); sensible-heat flux + drag still act.
    simple_physics = SimplePhysics(boundary_layer=False,
                                   use_external_surface_specific_humidity=True)
    dry_conv = DryConvectiveAdjustment()
    state = get_default_state([lw, slab, simple_physics, dry_conv], grid_state=grid)

    # forcing + initial condition (guard absorber fields that only RRTMG has)
    state["downwelling_shortwave_flux_in_air"].values[:] = SW_FLUX
    state["ocean_mixed_layer_thickness"].values[:] = MIXED_LAYER_DEPTH
    state["area_type"].values[:] = "sea"
    state["air_temperature"].values[:] = T_INIT
    state["surface_temperature"].values[:] = T_INIT
    state["specific_humidity"].values[:] = 0.0   # DRY
    if "surface_specific_humidity" in state:
        state["surface_specific_humidity"].values[:] = 0.0   # no evaporation
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2_MOLE_FRAC
    if "mole_fraction_of_ozone_in_air" in state:
        state["mole_fraction_of_ozone_in_air"].values[:] = 0.0   # apples-to-apples

    stepper = AdamsBashforth([lw, slab])
    return {"stepper": stepper, "sp": simple_physics, "dry": dry_conv,
            "state": state}


# Build RRTMG first (sets 16 LW bands for its own state), then CORK columns.
models = {}
for label, table in COLUMNS.items():
    models[label] = build_column(table)

# ── Time-stepping (operator splitting: tendencies -> simple physics -> dry adj)
history = {label: [] for label in COLUMNS}
for step in trange(N_STEPS, desc="dry RCE"):
    for label in COLUMNS:
        m = models[label]
        st = m["state"]

        diag, st = m["stepper"](st, TIMESTEP)   # LW + slab tendencies
        st.update(diag)

        diag, new = m["sp"](st, TIMESTEP)        # SimplePhysics (no BL)
        st.update(diag)
        st.update(new)
        if "surface_specific_humidity" in st:    # keep it dry (no evaporation)
            st["surface_specific_humidity"].values[:] = 0.0

        diag, new = m["dry"](st, TIMESTEP)       # dry convective adjustment
        st.update(diag)
        st.update(new)

        st["eastward_wind"].values[:] = 3.0      # sustain wind-dependent fluxes
        m["state"] = st

        if step % STORE_EVERY == 0 or step == N_STEPS - 1:
            up_lw = st["upwelling_longwave_flux_in_air"].values[:, 0, 0]
            m_hist = history[label]
            m_hist.append({
                "step": step,
                "p_hpa": st["air_pressure"].values[:, 0, 0].copy() / 100.0,
                "T": st["air_temperature"].values[:, 0, 0].copy(),
                "hr_lw": st["air_temperature_tendency_from_longwave"].values[:, 0, 0].copy(),
                "olr": float(up_lw[-1]),
                "T_sfc": float(st["surface_temperature"].values[0, 0]),
                "q_max": float(st["specific_humidity"].values.max()),
            })


def tropopause(p_hpa, T):
    """Convective-layer top via potential temperature.

    With O3=0 and no stratospheric heating the temperature decreases
    monotonically to the model top, so a cold-point tropopause is meaningless
    (it just returns the top level). The physically meaningful tropopause in
    dry RCE is the top of the well-mixed (constant-theta) convecting layer:
    scan up from the surface to the first level where theta rises >2 K above
    the surface value (the stable cap).
    """
    theta = T * (1000.0 / p_hpa) ** (287.0 / 1004.0)
    above = np.where(theta > theta[0] + 2.0)[0]
    k = int(above[0]) if above.size else len(p_hpa) - 1
    return p_hpa[k], T[k]


# ── Report ─────────────────────────────────────────────────────────────────
print(f"\nDry RCE — {args.days:.0f} days, SW={SW_FLUX:.0f} W/m^2, dry, O3=0\n")
print(f"{'column':22s} {'OLR':>7s} {'imbal':>7s} {'T_sfc':>7s} {'p_cnvtop':>9s} {'T_tp':>7s} {'q_max':>7s}")
for label in COLUMNS:
    s = history[label][-1]
    p_tp, T_tp = tropopause(s["p_hpa"], s["T"])
    imbal = SW_FLUX - s["olr"]   # >0 => still warming (not converged)
    qmax = s.get("q_max", float("nan")) * 1e3
    print(f"{label:22s} {s['olr']:7.2f} {imbal:+7.2f} {s['T_sfc']:7.2f} "
          f"{p_tp:9.2f} {T_tp:7.2f} {qmax:7.4f}")
print("  (imbal = SW-OLR; |imbal| >~1 W/m^2 => not yet converged.  q_max g/kg should be ~0)")

# ── Plot ───────────────────────────────────────────────────────────────────
fig, (axT, axHR) = plt.subplots(1, 2, figsize=(12, 7))
for label in COLUMNS:
    s = history[label][-1]
    c = COLORS.get(label, "#1f77b4")
    axT.plot(s["T"], s["p_hpa"], "-o", ms=3, color=c, label=label)
    axHR.plot(s["hr_lw"], s["p_hpa"], "-", color=c, label=label)
    p_tp, T_tp = tropopause(s["p_hpa"], s["T"])
    axT.plot(T_tp, p_tp, "*", color=c, ms=16, mec="k")
for ax in (axT, axHR):
    ax.set_yscale("log"); ax.invert_yaxis(); ax.grid(alpha=0.3)
    ax.set_ylabel("Pressure (hPa)"); ax.legend(fontsize=9)
axT.set_xlabel("Temperature (K)")
axT.set_title("Dry RCE temperature (★ = convective-layer top)")
axHR.set_xlabel("LW heating rate (K/day)")
axHR.set_title("LW heating rate")
fig.suptitle(f"Dry RCE: RRTMG-LW vs CORK-LW — SimplePhysics(no BL) + DryConvectiveAdjustment\n"
             f"dry, CO2={int(CO2_MOLE_FRAC*1e6)} ppm, O3=0, SW={SW_FLUX:.0f} W/m^2, {args.days:.0f} days")
plt.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(os.path.dirname(__file__), "..", "..", "debug_data", "rce_dry_cork_vs_rrtmg.png")
out = os.path.abspath(out)
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved {out}")

npz = os.path.join(os.path.dirname(out), "rce_dry_cork_vs_rrtmg.npz")
np.savez(npz, **{f"{l.replace(' ','_')}__{k}": np.asarray(history[l][-1][k])
                 for l in COLUMNS for k in ("p_hpa", "T", "hr_lw")})
print(f"Saved {npz}")

if args.save:
    payload = {}
    for label in COLUMNS:
        s = history[label][-1]
        safe = label.replace(" ", "_")
        payload[f"{safe}__p_hpa"] = np.asarray(s["p_hpa"])
        payload[f"{safe}__T"] = np.asarray(s["T"])
        # Dry columns track only a scalar q_max (not a profile); keep the key
        # name distinct from the moist script's profile-valued __q.
        payload[f"{safe}__q_max"] = np.asarray(s["q_max"])
        payload[f"{safe}__hr_lw"] = np.asarray(s["hr_lw"])
        payload[f"{safe}__T_sfc"] = np.asarray(s["T_sfc"])
        payload[f"{safe}__olr"] = np.asarray(s["olr"])
    np.savez(args.save, **payload)
    print(f"saved final state -> {args.save}")
