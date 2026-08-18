"""Static fallback figures for the Modelling Tour pages.

Every figure on these pages is drawn by a ``{pyodide}`` cell, which needs
WebAssembly. Where that is unavailable -- a locked-down browser, a PDF print, a
reader skimming on a phone -- the page still has to show what it is about, so
each headline cell has a collapsed callout holding the same figure as a PNG.

Those PNGs are rendered here **from the pages' own cells**, not from a
reimplementation of them: this script extracts the ``{pyodide}`` blocks out of
the ``.qmd`` files and execs them in order, exactly as a reader would run them,
then saves whatever figure the target cell drew. So a fallback cannot quietly
drift from the cell it stands in for -- edit the cell and re-run this, and the
PNG follows. The same trick as ``scripts/experiments/render_live_rce_figure.py``.

Cells run natively here (numba on, real climt) rather than under Pyodide, which
is the point: this is a fast check that the pages' code still runs at all, and
it fails loudly if a page's cell is broken.

    conda run -n climt python docs/modelling-tour/_artifacts/generate.py
    conda run -n climt python docs/modelling-tour/_artifacts/generate.py 04 06

Naming: ``NN-slug.png``, one per headline cell. Page 6 has two, because its
title makes two promises -- the feedback and the limit -- and they are drawn by
different cells.
"""
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")           # no display, and no figure windows to manage
import matplotlib.pyplot as plt  # noqa: E402  (must follow the backend choice)

PAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(PAGE_DIR, "_artifacts")

# page file -> (cell index to capture, output name). Cells are executed from 0
# up to the largest index requested, in order, sharing one namespace: later
# cells on these pages build on the state the earlier ones left behind.
FIGURES = {
    "01-emissivity-spectrum.qmd": [(1, "01-spectrum.png")],
    "02-window-measured.qmd": [(1, "02-window.png")],
    "03-radiating-level.qmd": [(1, "03-weighting.png")],
    "04-gray-equilibrium-tested.qmd": [(2, "04-gray-tested.png")],
    "05-co2-knob.qmd": [(1, "05-co2-knob.png")],
    "06-water-vapour-limit.qmd": [(1, "06-feedback.png"),
                                  (2, "06-runaway.png")],
}


def page_cells(qmd):
    """The page's ``{pyodide}`` cells, in document order, minus cell options."""
    text = open(os.path.join(PAGE_DIR, qmd)).read()
    cells = re.findall(r"```\{pyodide\}\n(.*?)```", text, flags=re.DOTALL)
    if not cells:
        raise SystemExit(f"no pyodide cells found in {qmd}")
    return ["\n".join(line for line in cell.splitlines()
                      if not line.startswith("#|")) for cell in cells]


def render(qmd, targets):
    """Exec ``qmd``'s cells in order, saving a PNG after each target cell."""
    cells = page_cells(qmd)
    wanted = dict(targets)
    namespace = {"__name__": "__main__"}
    for index in range(max(wanted) + 1):
        plt.close("all")
        exec(compile(cells[index], f"{qmd}[cell {index}]", "exec"), namespace)
        if index not in wanted:
            continue
        if not plt.get_fignums():
            raise SystemExit(f"{qmd} cell {index} drew no figure")
        # The last figure the cell opened is the one it ends by showing.
        figure = plt.figure(plt.get_fignums()[-1])
        path = os.path.join(OUT_DIR, wanted[index])
        figure.savefig(path, dpi=110, bbox_inches="tight")
        print(f"  wrote {wanted[index]} ({os.path.getsize(path) / 1024:.0f} kB)")
    plt.close("all")


def main():
    # The cells say `sys.path.insert(0, "_tour")` and reach the spectrum table
    # at `_data/...`, both relative to the page -- which is the working
    # directory in the browser. Match it.
    os.chdir(PAGE_DIR)
    selected = sys.argv[1:]
    for qmd, targets in FIGURES.items():
        if selected and not any(qmd.startswith(arg) for arg in selected):
            continue
        print(qmd)
        render(qmd, targets)


if __name__ == "__main__":
    main()
