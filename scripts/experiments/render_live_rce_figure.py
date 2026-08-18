"""Render the live-RCE page's figure natively, from the shipped code itself.

The page's plotting lives inside a ```{pyodide} cell in
docs/_includes/climt-live-setup.qmd, so it cannot be imported. This script
extracts that cell and execs it, then drives the real
``integrate_with_snapshots`` for both columns at the page's own defaults --
so what lands in debug_data/ is what a reader's browser would draw, not a
paraphrase of it.

    conda run -n climt python scripts/experiments/render_live_rce_figure.py
"""
import os
import re
from datetime import timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympl

from climt import (CorkLongwaveRadiation, SlabSurface, get_default_state,
                   get_grid)

INCLUDE = "docs/_includes/climt-live-setup.qmd"
SOLAR = 240.0
NZ, DT_HOURS, N_STEPS = 18, 2, 400        # the page's shipped defaults


def load_page_code():
    """Exec every ```{pyodide} cell of the include; return its namespace."""
    text = open(INCLUDE).read()
    cells = re.findall(r"```\{pyodide\}\n(.*?)```", text, flags=re.DOTALL)
    if not cells:
        raise SystemExit(f"no pyodide cells found in {INCLUDE}")
    ns = {}
    for cell in cells:
        body = "\n".join(ln for ln in cell.splitlines()
                         if not ln.startswith("#|"))
        exec(compile(body, INCLUDE, "exec"), ns)
    return ns


def make_column(table):
    lw = CorkLongwaveRadiation(optics="correlated_k", table=table)
    surface = SlabSurface()
    state = get_default_state([lw, surface], grid_state=get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = 5.0
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0
    return [lw, surface], state


def main():
    sympl.set_backend(sympl.DataArrayBackend())
    ns = load_page_code()
    integrate_with_snapshots = ns["integrate_with_snapshots"]

    os.makedirs("debug_data", exist_ok=True)
    for table, title, slug in [
        ("single_band_gray_lw", "Gray radiative equilibrium", "gray"),
        ("earth_low_res_lw", "Non-grey radiative equilibrium (14 bands)", "nongrey"),
    ]:
        components, state = make_column(table)
        integrate_with_snapshots(components, state, timedelta(hours=DT_HOURS),
                                 N_STEPS, SOLAR, title=title)
        out = f"debug_data/live_rce_{slug}.png"
        plt.gcf().savefig(out, dpi=110, bbox_inches="tight")
        plt.close("all")
        Ts = float(state["surface_temperature"].values.ravel()[0])
        print(f"wrote {out}   (surface temperature {Ts:.1f} K)")


if __name__ == "__main__":
    main()
