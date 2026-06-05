"""Experiment #16 — fixed-profile LW forward comparison: RRTMG-LW vs PF-LW.

Both dry (exp #14) and moist (exp #15) RCE show PF over-traps relative to RRTMG
(warmer surface, lower OLR), table-independent -> the cause is the radiative
scheme (thin bands / two-stream integrator), not the CO2/H2O k-tables (validated
exps #7/#13). This isolates that: ONE identical (T, p, q) profile, each scheme
called ONCE (no time-stepping, no convection, no feedback). Surface emissivity=1
and identical surface T -> the surface upward flux is identical, so any flux
divergence ALOFT is purely atmospheric trapping.

Compares broadband up/down LW flux and per-level heating rate at every level,
plus OLR and surface back-radiation. PF per-band TOA flux localizes spectrally
which PF band lets least escape (e.g. is the 800-1800 window leaking enough?).

Run: conda run --no-capture-output -n climt python scripts/experiments/lw_forward_pf_vs_rrtmg.py
"""
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from climt import RRTMGLongwave, get_default_state, get_grid
from climt._components.picket_fence import PicketFenceLongwave

_ap = argparse.ArgumentParser()
_ap.add_argument("--dry", action="store_true",
                 help="zero specific humidity (CO2-only) to isolate the solver "
                      "from H2O/window handling")
_ap.add_argument("--table", default="earth_low_res_lw_co2refined_linepyline",
                 help="PF correlated-k table name (without .nc)")
_ap.add_argument("--nz", type=int, default=40,
                 help="number of vertical levels (resolution test, exp #19)")
_args = _ap.parse_args()
TAG = "dry" if _args.dry else "moist"

CO2 = 376e-6
NZ = _args.nz
PF_TABLE = _args.table

grid = get_grid(nx=1, ny=1, nz=NZ)

# ── One shared, realistic reference profile (288 K moist troposphere + ~200 K
#    isothermal stratosphere; RH-like moisture decaying upward). ──────────────
def reference_profile(p_hpa):
    p_s = p_hpa.max()
    T = np.maximum(200.0, 288.0 * (p_hpa / p_s) ** 0.18)     # moist-adiabat-ish
    q = 0.015 * (p_hpa / p_s) ** 3.0                          # kg/kg, ~15 g/kg sfc
    return T, q


def run_scheme(make_comp, T_ref, q_ref, T_surf):
    """Build a column for one LW scheme, impose the shared profile, call once."""
    comp = make_comp()
    state = get_default_state([comp], grid_state=grid)
    state["air_temperature"].values[:, 0, 0] = T_ref
    state["specific_humidity"].values[:, 0, 0] = q_ref
    state["surface_temperature"].values[:] = T_surf
    if "surface_longwave_emissivity" in state:
        state["surface_longwave_emissivity"].values[:] = 1.0
    if "mole_fraction_of_carbon_dioxide_in_air" in state:
        state["mole_fraction_of_carbon_dioxide_in_air"].values[:] = CO2
    if "mole_fraction_of_ozone_in_air" in state:
        state["mole_fraction_of_ozone_in_air"].values[:] = 0.0   # PF has no O3
    _tendencies, diag = comp(state)   # TendencyComponent: (tendencies, diagnostics)
    return state, diag


# Build RRTMG first (sets 16 LW bands), grab the grid's pressure, build profile.
rr = RRTMGLongwave()
rr_state = get_default_state([rr], grid_state=grid)
p_hpa = rr_state["air_pressure"].values[:, 0, 0] / 100.0
p_int = rr_state["air_pressure_on_interface_levels"].values[:, 0, 0] / 100.0
T_ref, q_ref = reference_profile(p_hpa)
if _args.dry:
    q_ref = np.zeros_like(q_ref)   # CO2-only — isolate the solver from H2O/window
T_SURF = 288.0

_, dr = run_scheme(lambda: rr, T_ref, q_ref, T_SURF)          # reuse rr instance
_, dp = run_scheme(lambda: PicketFenceLongwave(optics="correlated_k", table=PF_TABLE),
                   T_ref, q_ref, T_SURF)

up_r = dr["upwelling_longwave_flux_in_air"].values[:, 0, 0]
dn_r = dr["downwelling_longwave_flux_in_air"].values[:, 0, 0]
hr_r = dr["air_temperature_tendency_from_longwave"].values[:, 0, 0]
up_p = dp["upwelling_longwave_flux_in_air"].values[:, 0, 0]
dn_p = dp["downwelling_longwave_flux_in_air"].values[:, 0, 0]
hr_p = dp["air_temperature_tendency_from_longwave"].values[:, 0, 0]

# ── Scalars ─────────────────────────────────────────────────────────────────
print(f"Fixed-profile LW forward comparison [{TAG.upper()}] (NZ={NZ}, "
      f"T_surf={T_SURF} K, CO2={int(CO2*1e6)} ppm, O3=0, emissivity=1, "
      f"q_sfc={q_ref[0]*1e3:.1f} g/kg)\n")
print(f"{'quantity':32s} {'RRTMG':>9s} {'PF':>9s} {'PF-RRTMG':>9s}")
print(f"{'OLR (TOA up)':32s} {up_r[-1]:9.2f} {up_p[-1]:9.2f} {up_p[-1]-up_r[-1]:+9.2f}")
print(f"{'surface back-radiation (sfc down)':32s} {dn_r[0]:9.2f} {dn_p[0]:9.2f} {dn_p[0]-dn_r[0]:+9.2f}")
print(f"{'surface up (sfc, should match)':32s} {up_r[0]:9.2f} {up_p[0]:9.2f} {up_p[0]-up_r[0]:+9.2f}")
print("  (PF OLR < RRTMG  AND/OR  PF back-radiation > RRTMG  => PF traps more)\n")

# ── Per-level up-flux divergence (where does PF emit less to space?) ──────────
print(f"{'p_hPa':>8s} {'up_RR':>8s} {'up_PF':>8s} {'d_up':>7s} | "
      f"{'dn_RR':>8s} {'dn_PF':>8s} {'d_dn':>7s}")
for k in range(len(p_int) - 1, -1, -1):
    print(f"{p_int[k]:8.2f} {up_r[k]:8.2f} {up_p[k]:8.2f} {up_p[k]-up_r[k]:+7.2f} | "
          f"{dn_r[k]:8.2f} {dn_p[k]:8.2f} {dn_p[k]-dn_r[k]:+7.2f}")

# ── PF per-band TOA flux: which PF band lets least escape? ────────────────────
if "upwelling_longwave_flux_in_air_per_band" in dp:
    ckdir = os.path.join("climt", "_data", "picket_fence", "correlated_k", PF_TABLE + ".nc")
    with xr.open_dataset(ckdir) as ds:
        be = ds["band_wavenumber_limits"].values
    up_band_toa = dp["upwelling_longwave_flux_in_air_per_band"].values[-1, 0, 0, :]
    print(f"\nPF per-band OLR (TOA upward), total={up_band_toa.sum():.2f} W/m^2:")
    for b in range(len(be)):
        print(f"  band {b} [{be[b,0]:6.0f}-{be[b,1]:6.0f} cm^-1]  {up_band_toa[b]:7.2f} W/m^2")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(18, 6))
axT, axUp, axDn, axHR = axes
axT.plot(T_ref, p_hpa, "k-o", ms=3, label="T")
axT2 = axT.twiny(); axT2.plot(q_ref * 1e3, p_hpa, "b--", label="q"); axT2.set_xlabel("q (g/kg)", color="b")
axT.set_xlabel("T (K)"); axT.set_title("Shared reference profile")
axUp.plot(up_r, p_int, "k-", label="RRTMG"); axUp.plot(up_p, p_int, "r-", label="PF")
axUp.set_xlabel("Upwelling LW (W/m^2)"); axUp.set_title("Upward flux (lower PF = more trapping)")
axDn.plot(dn_r, p_int, "k-", label="RRTMG"); axDn.plot(dn_p, p_int, "r-", label="PF")
axDn.set_xlabel("Downwelling LW (W/m^2)"); axDn.set_title("Downward flux")
axHR.plot(hr_r, p_hpa, "k-", label="RRTMG"); axHR.plot(hr_p, p_hpa, "r-", label="PF")
axHR.set_xlabel("LW heating rate (K/day)"); axHR.set_title("Heating rate")
for ax in (axT, axUp, axDn, axHR):
    ax.set_yscale("log"); ax.invert_yaxis(); ax.grid(alpha=0.3); ax.set_ylabel("p (hPa)")
    if ax is not axT:
        ax.legend()
fig.suptitle("Fixed-profile LW forward: RRTMG vs PF (identical T,p,q) — isolates the radiative scheme")
plt.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                   "debug_data", f"lw_forward_pf_vs_rrtmg_{TAG}_{PF_TABLE}.png"))
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved {out}")
