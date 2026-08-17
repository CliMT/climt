"""Is the isothermal stratosphere invisible because of the linear-p y axis?

Reproduces the two columns of docs/radiative-transfer/09-live-rce.qmd natively
(gray `single_band_gray_lw` and non-grey `earth_low_res_lw`) and draws each
final temperature profile twice: once on the page's current linear-pressure
axis, once on a log-pressure axis. Prints where the isothermal top actually
sits and what fraction of each axis it gets.

    conda run -n climt python scripts/experiments/live_rce_axis_check.py

Writes debug_data/live_rce_axis_check.png.
"""
import os
from datetime import timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympl

from climt import (CorkLongwaveRadiation, SlabSurface, get_default_state,
                   get_grid)

SOLAR = 240.0
NZ = 18

# Near-equilibrium budgets (see the RCE convergence notes): the docs page's
# default 400 steps at 2 h is only ~33 days and nowhere near settled, so run
# each column long enough that the isothermal top genuinely exists.
RUNS = [
    ("single_band_gray_lw", "Gray", 12, 300),
    ("earth_low_res_lw", "Non-grey (14 band)", 6, 1950),
]


def make_column(table):
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
    surface = SlabSurface()
    state = get_default_state([lw, surface], grid_state=get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = 5.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [lw, surface], state


def integrate(components, state, timestep, n_steps):
    # Matches integrate_with_snapshots in docs/_includes/climt-live-setup.qmd:
    # both cork LW and SlabSurface are tendency components, so both go into the
    # single AdamsBashforth stepper.
    model = sympl.AdamsBashforth(list(components))
    for _ in range(n_steps):
        diagnostics, new_state = model(state, timestep)
        state.update(new_state)
        state.update(diagnostics)
        state["time"] += timestep
    return state


def isothermal_top_span(T, p, tol=1.0):
    """Index of the lowest level of the near-isothermal top, by T within tol of T[-1]."""
    k = len(T) - 1
    while k > 0 and abs(T[k - 1] - T[-1]) < tol:
        k -= 1
    return k


def main():
    sympl.set_backend(sympl.DataArrayBackend())
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    for row, (table, label, dt_hours, n_steps) in enumerate(RUNS):
        components, state = make_column(table)
        state = integrate(components, state, timedelta(hours=dt_hours), n_steps)
        T = state["air_temperature"].values[:, 0, 0]
        p = state["air_pressure"].values[:, 0, 0] / 100.0

        k = isothermal_top_span(T, p)
        p_lo, p_hi = p[-1], p[k]
        frac_linear = (p_hi - p_lo) / (p.max() - p.min())
        frac_log = ((np.log10(p_hi) - np.log10(p_lo))
                    / (np.log10(p.max()) - np.log10(p.min())))
        print(f"\n{label} ({table}), {n_steps} steps @ {dt_hours} h "
              f"= {n_steps * dt_hours / 24:.0f} days")
        print(f"  surface T      = {float(state['surface_temperature'].values.ravel()[0]):.1f} K")
        print(f"  T range        = {T.min():.1f} - {T.max():.1f} K")
        print(f"  isothermal top = levels {k}..{len(T)-1} "
              f"({p_hi:.1f} -> {p_lo:.1f} hPa), T = {T[k]:.1f} - {T[-1]:.1f} K")
        print(f"  share of y axis: linear-p {100*frac_linear:5.1f}%   "
              f"log-p {100*frac_log:5.1f}%")

        for col, logscale in enumerate((False, True)):
            ax = axes[row, col]
            ax.plot(T, p, "o-", color="#1c7ed6", lw=2, ms=4)
            ax.plot(T[k:], p[k:], "o-", color="#c92a2a", lw=2.6, ms=5,
                    label="isothermal top")
            if logscale:
                ax.set_yscale("log")
            ax.invert_yaxis()
            ax.grid(alpha=0.3, which="both")
            ax.set_xlabel("temperature (K)")
            ax.set_ylabel("pressure (hPa)")
            ax.set_title(f"{label} — {'log' if logscale else 'linear'} p "
                         f"({100 * (frac_log if logscale else frac_linear):.0f}% of axis)")
            ax.legend(loc="lower left", fontsize=8)

    fig.tight_layout()
    os.makedirs("debug_data", exist_ok=True)
    out = "debug_data/live_rce_axis_check.png"
    fig.savefig(out, dpi=120)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
