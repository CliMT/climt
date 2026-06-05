"""Diagnose the vertical opacity structure of the PF k-tables, using the
component's own `longwave_optical_depth_per_band` diagnostic (no offline
interp — guaranteed consistent with the production calc).

For each band, accumulate τ_layer top-down to get τ(TOA → p). The level
where τ_cum = 1 is the band's emission-to-space height.

Outputs: examples/diagnose_pf_opacity_structure.png + console table.
"""
import os
from datetime import timedelta

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from sympl import AdamsBashforth
from tqdm import trange

from climt import SlabSurface, get_default_state, get_grid
from climt._components.picket_fence import PicketFenceLongwave

SW_FLUX = 240.0
SIGMA_SB = 5.670374419e-8
T_EFF = (SW_FLUX / SIGMA_SB) ** 0.25
T_SKIN = T_EFF / (2.0 ** 0.25)
TIMESTEP = timedelta(hours=3)
N_STEPS = 8000
NZ = 30
TABLES = {
    "PF 6-band GL":       "earth_low_res_lw_co2refined_gl",
    "PF 10-band GL":      "earth_low_res_lw_10band_gl",
    "PF 4-band 16gpt 2s": "earth_low_res_lw_4band_16gpt",
    "PF 4-band 32gpt 2s": "earth_low_res_lw_4band_32gpt",
    "PF 6-band 16gpt GL": "earth_low_res_lw_6band_16gpt",
}

data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "climt", "_data",
                        "picket_fence", "correlated_k")


def equilibrate_and_diagnose(table_name):
    grid = get_grid(nx=1, ny=1, nz=NZ)
    lw = PicketFenceLongwave(optics="correlated_k", table=table_name)
    sf = SlabSurface()
    st = get_default_state([lw, sf], grid_state=grid)
    st["specific_humidity"].values[:] = 0.0
    st["downwelling_shortwave_flux_in_air"].values[:] = SW_FLUX
    st["ocean_mixed_layer_thickness"].values[:] = 1.0
    st["area_type"].values[:] = "sea"
    st["air_temperature"].values[:] = T_SKIN
    st["surface_temperature"].values[:] = T_SKIN
    stepper = AdamsBashforth([lw, sf])
    for _ in trange(N_STEPS, desc=f"Equilibrate {table_name}"):
        diag, st = stepper(st, TIMESTEP)
        st.update(diag)
    # One last forward call to grab the per-band optical depth at the converged state.
    diag, _ = stepper(st, TIMESTEP)
    tau_band = diag["longwave_optical_depth_per_band"].values[:, 0, 0, :]  # (nlev, nband)
    p_mid = st["air_pressure"].values[:, 0, 0]                              # Pa
    T_prof = st["air_temperature"].values[:, 0, 0]
    band_edges = xr.open_dataset(os.path.join(data_dir, f"{table_name}.nc"))[
        "band_wavenumber_limits"].values
    return p_mid, T_prof, tau_band, band_edges


fig, axes = plt.subplots(len(TABLES), 2, figsize=(13, 5 * len(TABLES)))
if len(TABLES) == 1:
    axes = axes[None, :]

for i, (label, tn) in enumerate(TABLES.items()):
    p_mid, T_prof, tau_band, band_edges = equilibrate_and_diagnose(tn)
    nlev, nband = tau_band.shape
    # Determine vertical orientation (level 0 = TOA?)
    descending = p_mid[-1] > p_mid[0]
    if descending:
        tau_cum = np.cumsum(tau_band, axis=0)            # along level (TOA → sfc)
    else:
        tau_cum = np.cumsum(tau_band[::-1, :], axis=0)[::-1, :]

    print(f"\n=== {label} ({tn}) ===")
    print(f"T_strat = {T_prof.min():.1f} K   T_sfc = {T_prof.max():.1f} K   "
          f"p_top = {p_mid.min()/100:.2f} hPa")
    print(f"{'band':>4s}  {'[ν1-ν2 cm⁻¹]':>18s}  {'τ(sfc)':>10s}  {'p(τ=1) hPa':>12s}")
    for b in range(nband):
        nu1, nu2 = band_edges[b]
        tc = tau_cum[:, b]
        tau_sfc = tc[-1] if descending else tc[0]
        # Find p where τ_cum = 1
        # tau_cum is monotonic from TOA (small) → surface (large)
        order = np.argsort(p_mid)  # ascending p (TOA → sfc)
        tc_ord = tc[order]
        p_ord = p_mid[order]
        if tau_sfc < 1.0:
            p_tau1_s = "  > sfc"
        else:
            idx = int(np.searchsorted(tc_ord, 1.0))
            if idx == 0:
                p_tau1 = p_ord[0]
            else:
                denom = max(tc_ord[idx] - tc_ord[idx - 1], 1e-30)
                f = (1.0 - tc_ord[idx - 1]) / denom
                p_tau1 = p_ord[idx - 1] + f * (p_ord[idx] - p_ord[idx - 1])
            p_tau1_s = f"{p_tau1/100:7.2f}"
        print(f"  b{b}  [{nu1:6.0f}-{nu2:6.0f}]  {tau_sfc:10.3e}  {p_tau1_s}")

    cmap = plt.get_cmap("turbo", nband)
    ax_tau, ax_layer = axes[i, 0], axes[i, 1]
    for b in range(nband):
        nu1, nu2 = band_edges[b]
        lbl = f"b{b}: {nu1:.0f}-{nu2:.0f}"
        ax_tau.plot(tau_cum[:, b], p_mid / 100.0, color=cmap(b), lw=1.8, label=lbl)
        ax_layer.plot(tau_band[:, b], p_mid / 100.0, color=cmap(b), lw=1.8, label=lbl)
    ax_tau.set_xscale("log"); ax_tau.set_yscale("log"); ax_tau.invert_yaxis()
    ax_tau.axvline(1.0, color="k", ls=":", lw=1)
    ax_tau.set_xlabel("τ (TOA → p)")
    ax_tau.set_ylabel("Pressure (hPa)")
    ax_tau.set_title(f"{label}: cumulative τ\n(where τ=1 → emission to space)")
    ax_tau.grid(True, alpha=0.3); ax_tau.legend(fontsize=7)

    ax_layer.set_xscale("log"); ax_layer.set_yscale("log"); ax_layer.invert_yaxis()
    ax_layer.set_xlabel("Layer τ per band")
    ax_layer.set_title(f"{label}: layer-by-layer τ")
    ax_layer.grid(True, alpha=0.3); ax_layer.legend(fontsize=7)

fig.suptitle("PF opacity vertical structure (from component diagnostic) — "
             "does the wing collapse aloft?", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(__file__), "diagnose_pf_opacity_structure.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved {out}")
