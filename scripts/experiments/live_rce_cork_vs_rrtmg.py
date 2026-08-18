"""Is the live demo's cold, non-isothermal top a CORK artifact, or the setup?

The committed cork-vs-RRTMG experiment (docs/experiments/2026-06-05-cork-co2-bands)
compares MOIST RCE with convection. The live demo is a different animal: dry,
no convection, pure radiative equilibrium with prescribed surface shortwave. So
run that exact configuration with RRTMG-LW in place of CORK-LW and compare the
equilibrium profiles. If RRTMG also cools monotonically to a very cold lid, the
behaviour belongs to the setup (an atmosphere still optically thick in the CO2
15 um core at the model top), not to CORK.

    conda run -n climt python scripts/experiments/live_rce_cork_vs_rrtmg.py
"""
import os
from datetime import timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympl

from climt import (CorkLongwaveRadiation, RRTMGLongwave, SlabSurface,
                   get_default_state, get_grid)

SOLAR = 240.0
NZ = 18
DT = timedelta(hours=6)
N_STEPS = 4000            # 1000 days; the top is slow to relax


def make_column(lw):
    surface = SlabSurface()
    state = get_default_state([lw, surface], grid_state=get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = 5.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [lw, surface], state


def run(lw, label):
    components, state = make_column(lw)
    model = sympl.AdamsBashforth(list(components))
    for _ in range(N_STEPS):
        diagnostics, new_state = model(state, DT)
        state.update(new_state)
        state.update(diagnostics)
        state["time"] += DT
    T = state["air_temperature"].values[:, 0, 0].copy()
    p = state["air_pressure"].values[:, 0, 0] / 100.0
    olr = float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0])
    ts = float(state["surface_temperature"].values.ravel()[0])
    print(f"{label:10s} Ts={ts:7.2f} K  OLR={olr:7.2f} W/m2  "
          f"T_top={T[-1]:7.2f} K  dT(top2)={abs(T[-1] - T[-2]):5.2f} K")
    return T, p


def main():
    sympl.set_backend(sympl.DataArrayBackend())
    print(f"dry pure radiative equilibrium, nz={NZ}, {N_STEPS} steps @ {DT} "
          f"= {N_STEPS * DT.total_seconds() / 86400:.0f} days\n")

    T_cork, p = run(CorkLongwaveRadiation(optics="correlated_k",
                                          table="earth_low_res_lw"), "CORK-14b")
    T_rrtmg, _ = run(RRTMGLongwave(), "RRTMG")
    T_gray, _ = run(CorkLongwaveRadiation(optics="correlated_k",
                                          table="single_band_gray_lw"), "CORK-gray")

    print("\n lev   p(hPa)     CORK-14b      RRTMG   diff      CORK-gray")
    for k in range(NZ - 1, -1, -1):
        print(f" {k:3d} {p[k]:9.2f} {T_cork[k]:11.2f} {T_rrtmg[k]:10.2f} "
              f"{T_cork[k] - T_rrtmg[k]:+7.2f} {T_gray[k]:13.2f}")

    fig, ax = plt.subplots(figsize=(6.5, 7))
    ax.plot(T_cork, p, "o-", label="CORK-LW (14 band)", color="#1c7ed6")
    ax.plot(T_rrtmg, p, "s--", label="RRTMG-LW", color="#2f9e44")
    ax.plot(T_gray, p, "^-", label="CORK-LW (gray)", color="#c92a2a", alpha=0.7)
    ax.set_yscale("log")
    ax.set_ylim(p.max(), p.min())
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("pressure (hPa)")
    ax.set_title("Dry radiative equilibrium: the live demo's configuration")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    os.makedirs("debug_data", exist_ok=True)
    out = "debug_data/live_rce_cork_vs_rrtmg.png"
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
