"""Render page 7's evolution figure natively, from the shipped code itself.

The plotting used to live inside a ```{pyodide} cell in a Quarto include and
had to be extracted and exec'd. It is now `_tour/stepping.py`, so this just
imports it -- what lands in debug_data/ is what a reader's browser draws.

    conda run -n climt python scripts/experiments/render_page7_figure.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402  (must follow the backend choice)
import sympl                       # noqa: E402

import climt                       # noqa: E402

sys.path.insert(0, os.path.join("docs", "modelling-tour", "_tour"))
import stepping                    # noqa: E402

SOLAR = 240.0
NZ, DT_HOURS, N_STEPS, SLAB_DEPTH = 28, 12, 675, 2.0
OUT_DIR = "debug_data"


def column(table, diffusivity):
    longwave = climt.CorkLongwaveRadiation(optics="correlated_k", table=table,
                                           diffusivity_factor=diffusivity)
    surface = climt.SlabSurface()
    state = climt.get_default_state(
        [longwave, surface], grid_state=climt.get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = SLAB_DEPTH
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [longwave, surface], state


def main():
    sympl.set_backend(climt.UnytBackend())
    os.makedirs(OUT_DIR, exist_ok=True)
    # Page 7's own tables and diffusivities: `tour_gray_lw` with D = 2 is what
    # reproduces page 4's analytic profile, and is what the page runs.
    for table, diffusivity, slug, title in (
            ("tour_gray_lw", 2.0, "gray", "Gray radiative equilibrium"),
            ("earth_low_res_lw", 1.66, "nongrey",
             "Non-grey radiative equilibrium (14 bands)")):
        components, state = column(table, diffusivity)
        state, history = stepping.integrate_with_history(
            components, [], state, climt.UnytTimeDelta(hours=DT_HOURS),
            N_STEPS)
        stepping.draw_evolution(history, state, SOLAR, title=title)
        path = os.path.join(OUT_DIR, f"page7_{slug}.png")
        plt.savefig(path, dpi=110, bbox_inches="tight")
        plt.close("all")
        print(f"wrote {path}  Tsurf="
              f"{float(state['surface_temperature'].values.ravel()[0]):.2f} K")


if __name__ == "__main__":
    main()
