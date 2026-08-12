"""The science shown by the in-browser RCE demo page, checked natively.

``docs/radiative-transfer/09-live-rce.qmd`` runs two single-column radiative
equilibrium integrations in the browser under Pyodide: a *gray* column
(``single_band_gray_lw``, one band) and a *non-grey* column
(``earth_low_res_lw``, 14 bands of real H2O/CO2 line absorption). Nothing else
differs — same component, same state, same stepping — so the page's whole
teaching claim rests on what the spectral structure alone does.

Two things here are untestable from the browser and so are pinned down natively:

1. **The twin stays a twin.** ``integrate_to_equilibrium`` is defined twice: in
   the Pyodide boot cell of ``docs/_includes/climt-live-setup.qmd`` (what the
   page actually runs — Pyodide cannot import the ``.py``, since
   ``docs/radiative-transfer`` has a hyphen and is not a valid package path) and
   in ``docs/radiative-transfer/_live/rce_helpers.py`` (what is importable and
   testable). Both docstrings instruct that the bodies be kept identical; only
   this test enforces it. It is deliberately cheap, so it runs in routine CI.

2. **The physics the page claims.** Marked ``slow`` (see tox.ini) — the two
   integrations take ~20 s together.

Note what is deliberately *not* asserted: that either column reaches
equilibrium, or that the non-grey column has an isothermal stratosphere. At the
page's step count neither column is equilibrated (the page says so), and the
isothermal-stratosphere discussion describes the limiting state, not what 200
steps reach.
"""
import ast
import importlib.util
from datetime import timedelta
from pathlib import Path

import numpy as np
import pytest
import sympl

from climt import CorkLongwaveRadiation, SlabSurface, get_default_state, get_grid

REPO_ROOT = Path(__file__).resolve().parent.parent
INCLUDE_QMD = REPO_ROOT / "docs/_includes/climt-live-setup.qmd"
HELPERS_PY = REPO_ROOT / "docs/radiative-transfer/_live/rce_helpers.py"

# Mirrors the flagship page's own cell constants, so this test exercises the
# configuration readers actually see. Steps are halved from the page's 400: the
# gray/non-grey separation below is already unambiguous at 200 and it keeps the
# test near 20 s.
SOLAR = 240.0     # prescribed absorbed shortwave at the surface, W m^-2
NZ = 18
DT_HOURS = 2
N_STEPS = 200


def _load_helpers():
    """Import rce_helpers.py by path (its directory is not a package)."""
    spec = importlib.util.spec_from_file_location("rce_helpers", HELPERS_PY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pyodide_cell_source(qmd_path):
    """Concatenate the Python inside the ```{pyodide} fences of a .qmd."""
    lines, inside = [], False
    for line in qmd_path.read_text().splitlines():
        if line.startswith("```{pyodide}"):
            inside = True
        elif inside and line.startswith("```"):
            inside = False
        elif inside:
            # Skip quarto cell options (`#| autorun: true`), which are not Python.
            if not line.startswith("#|"):
                lines.append(line)
    return "\n".join(lines)


def _function_fingerprint(source, name):
    """Signature + docstring-stripped body of `name`, as a comparable string.

    Compares the parsed AST rather than text, so the two copies may differ in
    indentation, comments and docstring (they cite each other's paths) while
    any difference in actual logic still fails.
    """
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            body = node.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                body = body[1:]
            return (ast.dump(node.args),
                    ast.dump(ast.Module(body=body, type_ignores=[])))
    raise AssertionError(f"no function named {name!r} found")


def test_helper_twin_matches_the_page_include():
    """The .py twin and the .qmd copy of integrate_to_equilibrium must agree.

    If this fails, the page and the tested helper have drifted: whichever was
    edited, port the change to the other rather than relaxing this test.
    """
    from_qmd = _function_fingerprint(
        _pyodide_cell_source(INCLUDE_QMD), "integrate_to_equilibrium")
    from_py = _function_fingerprint(
        HELPERS_PY.read_text(), "integrate_to_equilibrium")

    assert from_qmd[0] == from_py[0], (
        "integrate_to_equilibrium signatures differ between "
        f"{INCLUDE_QMD.name} and {HELPERS_PY.name}")
    assert from_qmd[1] == from_py[1], (
        "integrate_to_equilibrium bodies differ between "
        f"{INCLUDE_QMD.name} and {HELPERS_PY.name} — keep them identical")


def _run_column(table):
    """Integrate one column exactly as the page's `make_column` cell does."""
    sympl.set_backend(sympl.DataArrayBackend())
    longwave = CorkLongwaveRadiation(optics="correlated_k", table=table)
    surface = SlabSurface()
    state = get_default_state([longwave, surface],
                             grid_state=get_grid(nx=1, ny=1, nz=NZ))
    state["ocean_mixed_layer_thickness"].values[:] = 5.0
    # Fixed shortwave heating of the slab; there is no shortwave component.
    state["downwelling_shortwave_flux_in_air"].values[:] = 0.0
    state["downwelling_shortwave_flux_in_air"].values[0, ...] = SOLAR
    state["upwelling_shortwave_flux_in_air"].values[:] = 0.0

    # The RCE demo has no stepper-kind components: SlabSurface is a
    # TendencyComponent, so both go in the tendency list.
    return _load_helpers().integrate_to_equilibrium(
        [longwave, surface], [], state, timedelta(hours=DT_HOURS), N_STEPS)


@pytest.fixture(scope="module")
def columns():
    """Both integrations, run once and shared (each is seconds, not ms)."""
    return {table: _run_column(table)
            for table in ("single_band_gray_lw", "earth_low_res_lw")}


def _surface_temperature(state):
    return float(state["surface_temperature"].values.ravel()[0])


def _olr(state):
    """Upwelling longwave flux at the top of the column."""
    return float(state["upwelling_longwave_flux_in_air"].values[-1, 0, 0])


@pytest.mark.slow
def test_nongrey_column_is_the_more_efficient_radiator(columns):
    """The page's core claim: spectral windows let the column radiate better.

    The non-grey table resolves atmospheric windows through which surface
    emission escapes more or less directly to space. The single-band gray
    column has no such windows, so at matched conditions the non-grey column
    emits more to space and runs cooler at the surface.
    """
    gray, nongrey = columns["single_band_gray_lw"], columns["earth_low_res_lw"]

    # Measured separation at these settings is ~38 W m^-2 and ~5.9 K; the
    # thresholds keep a wide margin so this checks the physics, not a fit.
    assert _olr(nongrey) > _olr(gray) + 20.0, (
        f"non-grey OLR {_olr(nongrey):.1f} should exceed gray "
        f"{_olr(gray):.1f} W/m^2 — windows radiate to space more efficiently")
    assert _surface_temperature(nongrey) < _surface_temperature(gray) - 3.0, (
        f"non-grey surface {_surface_temperature(nongrey):.1f} K should be "
        f"cooler than gray {_surface_temperature(gray):.1f} K")


@pytest.mark.slow
def test_nongrey_column_cools_faster_aloft(columns):
    """The page's second claim: a steeper temperature drop-off aloft.

    Band-resolved absorption concentrates cooling in the strongly absorbing
    bands high in the column, which the gray column smears out.
    """
    top_gradient = {}
    for name, state in columns.items():
        temperature = state["air_temperature"].values[:, 0, 0]
        # Level index increases upward; mean of the top three level-to-level
        # differences, in K per level.
        top_gradient[name] = float(np.mean(np.diff(temperature)[-3:]))

    assert top_gradient["earth_low_res_lw"] < top_gradient["single_band_gray_lw"] - 4.0, (
        f"non-grey top-of-column gradient {top_gradient['earth_low_res_lw']:.2f} "
        f"K/level should be markedly steeper than gray "
        f"{top_gradient['single_band_gray_lw']:.2f} K/level")


@pytest.mark.slow
@pytest.mark.parametrize("table", ["single_band_gray_lw", "earth_low_res_lw"])
def test_flux_diagnostics_survive_the_update_order(columns, table):
    """Guards the update-order gotcha that broke this demo once already.

    ``integrate_to_equilibrium`` applies the prognostic state *before* the
    diagnostics. Reversing that clobbers the longwave fluxes SlabSurface
    consumes, and the surface heats without bound. So: the flux diagnostics
    must still be present and finite at the end, and the temperatures must stay
    physical rather than running away.
    """
    state = columns[table]

    for flux in ("upwelling_longwave_flux_in_air",
                 "downwelling_longwave_flux_in_air"):
        assert flux in state, f"{flux} missing — diagnostics were clobbered"
        values = state[flux].values
        assert np.all(np.isfinite(values)), f"{flux} has non-finite values"
        assert np.all(values >= 0.0), f"{flux} must be non-negative"

    temperature = state["air_temperature"].values
    assert np.all(np.isfinite(temperature))
    assert temperature.min() > 150.0 and temperature.max() < 400.0, (
        f"air temperature ran to [{temperature.min():.1f}, "
        f"{temperature.max():.1f}] K — suspect a flux-coupling regression")
    surface = _surface_temperature(state)
    assert 150.0 < surface < 400.0, f"surface temperature {surface:.1f} K"
