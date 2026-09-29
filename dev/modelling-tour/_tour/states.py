"""Save a climt state to a small .npz, ship it, and load it back.

Pages 11 and 12 perturb an equilibrium rather than spending several thousand
in-browser steps finding one. The equilibrium therefore has to travel: as a
file in ``docs/modelling-tour/_data/``, published and staged by the same
mechanism as the spectrum table (see ``_tour/assets.py``).

**An equilibrium is only meaningful attached to its configuration**, so every
file carries provenance -- climt version, table, grid, timestep, step count,
slab depth, prescribed shortwave, CO2, the component list, and the TOA and
surface imbalances it finished at. :func:`describe` turns that into the block
the pages print above their first figure, so a reader can never be looking at
an equilibrium without also looking at what produced it.

**Loading rebuilds; it does not unpickle.** :func:`load` calls
``get_default_state`` for the caller's own component list and overwrites the
arrays, checking dims and units per quantity. A sympl or climt change that
renames a quantity, reorders its dimensions or changes its units then fails at
a named quantity with a message -- instead of yielding a state that looks fine
and is not an equilibrium.
"""
import datetime
import json

import numpy as np
from climt import get_default_state

import assets

# Rebuilt by get_default_state, or not an array at all. Saving them would only
# create a second place for the grid to be defined.
SKIPPED = (
    "time",
    "atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels",
    "atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels",
    "longitude",
    "latitude",
    "height_on_ice_interface_levels",
    "height_on_soil_interface_levels",
)


def save(path, state, provenance):
    """Write ``state`` and its provenance to ``path`` as an .npz.

    Args:
        path: output file path.
        state: a climt state dict.
        provenance: dict of configuration facts. ``climt_version`` and
            ``saved_at`` are added here; anything else is the caller's.

    Returns:
        ``path``.
    """
    import climt

    meta = dict(provenance)
    meta["climt_version"] = climt.__version__
    # `datetime.utcnow()` is deprecated from Python 3.12 (which CI builds), so
    # ask for UTC explicitly. The tzinfo is then dropped again so the stored
    # string keeps its naive-UTC shape, "2026-08-20T10:00:00": already-saved
    # states and `describe` both print it verbatim, and an appended "+00:00"
    # would make the shipped equilibria disagree with the ones a reader saves.
    meta["saved_at"] = (
        datetime.datetime.now(datetime.timezone.utc)
        .replace(tzinfo=None)
        .isoformat(timespec="seconds"))

    arrays = {}
    layout = {}
    for name, quantity in state.items():
        if name in SKIPPED:
            continue
        values = np.asarray(quantity.values)
        arrays["value__" + name] = values
        layout[name] = dict(dims=list(quantity.dims),
                            units=quantity.attrs.get("units"))

    arrays["__meta__"] = np.array(json.dumps(meta))
    arrays["__layout__"] = np.array(json.dumps(layout))
    np.savez_compressed(path, **arrays)
    return path


def load(name, components, grid_state=None, base_url=assets.DEFAULT_BASE_URL):
    """Rebuild a saved state for ``components``.

    Args:
        name: the asset's file name (resolved through ``_tour/assets.py``) or
            a direct path.
        components: the component list the page will run. The state is built
            for *these*, not for whatever produced the file.
        grid_state: passed to ``get_default_state``; must match the saved grid.
        base_url: directory holding the asset, relative to the page.

    Returns:
        ``(state, provenance)``, or ``(None, None)`` if the asset is not
        available -- so a page whose resource failed to stage can say so
        rather than traceback.

    Raises:
        ValueError: if the file is there but does not fit ``components`` --
            a missing quantity, or a dims/units mismatch.
    """
    path = name if "/" in name or "\\" in name else None
    if path is None:
        path = assets.locate(name, base_url)
    elif not _exists(path):
        path = assets.locate(name, base_url)
    if path is None:
        return None, None

    with np.load(path, allow_pickle=False) as handle:
        meta = json.loads(str(handle["__meta__"]))
        layout = json.loads(str(handle["__layout__"]))
        saved = {key[len("value__"):]: handle[key][...]
                 for key in handle.files if key.startswith("value__")}

    state = get_default_state(components, grid_state=grid_state)

    problems = []
    for quantity_name, quantity in state.items():
        if quantity_name in SKIPPED:
            continue
        if quantity_name not in saved:
            raise ValueError(
                f"{quantity_name!r} is required by these components but is "
                f"not in the saved state {path!r}. The file was written for a "
                f"different component list ({meta.get('components')}); "
                "regenerate it with scripts/generate_tour_equilibria.py.")

        expected = layout[quantity_name]
        if tuple(expected["dims"]) != tuple(quantity.dims):
            problems.append(
                f"{quantity_name!r} was saved with dims "
                f"{tuple(expected['dims'])} but this state wants "
                f"{tuple(quantity.dims)}")
        elif saved[quantity_name].shape != quantity.values.shape:
            problems.append(
                f"{quantity_name!r} was saved with shape "
                f"{saved[quantity_name].shape} but this state wants "
                f"{quantity.values.shape}")
        elif expected["units"] != quantity.attrs.get("units"):
            problems.append(
                f"{quantity_name!r} was saved in {expected['units']!r} but "
                f"this state wants {quantity.attrs.get('units')!r}")

    # Reported together, and by name: one wrong grid shows up in every
    # quantity on it, and a reader needs to see that it is the grid rather
    # than whichever quantity happened to be first in the dict.
    if problems:
        raise ValueError(
            "the saved state does not fit these components (the file says "
            f"nz={meta.get('nz')}; check grid_state, or a units convention "
            "changed under this file and it needs regenerating):\n  "
            + "\n  ".join(problems))

    for quantity_name, quantity in state.items():
        if quantity_name in SKIPPED:
            continue
        quantity.values[...] = saved[quantity_name]

    return state, meta


def describe(provenance):
    """A one-block human summary of what an equilibrium state is.

    Pages print this above their first figure. It is deliberately plain text
    rather than markdown: it goes through ``print()`` in a ``{pyodide}`` cell.
    """
    components = ", ".join(provenance.get("components", []))
    dt_hours = provenance.get("dt_hours")
    # A sub-hour step reads better in minutes: "5 min", not "0.0833... h".
    dt = ("{:g} min".format(round(60 * dt_hours, 6))
          if isinstance(dt_hours, float) and dt_hours < 1
          else "{} h".format(dt_hours))
    lines = [
        "equilibrium state, as shipped",
        "  climt          {}".format(provenance.get("climt_version")),
        "  saved          {}".format(provenance.get("saved_at")),
        "  LW table       {}".format(provenance.get("table")),
        "  grid           nz = {}, single column".format(provenance.get("nz")),
        "  spin-up        {} steps at dt = {}".format(
            provenance.get("n_steps"), dt),
        "  slab depth     {} m".format(provenance.get("slab_depth_m")),
        "  absorbed SW    {} W/m^2 (prescribed at the surface)".format(
            provenance.get("solar")),
        "  CO2            {} ppm".format(provenance.get("co2_ppm")),
        "  wind           relaxed to {} m/s over {} h, z0 = {} m".format(
            provenance.get("wind_m_s"), provenance.get("wind_timescale_hours"),
            provenance.get("roughness_length_m")),
        "  components     {}".format(components),
        "  residual       TOA {:+.3f}, surface {:+.3f} W/m^2".format(
            provenance.get("toa_imbalance", float("nan")),
            provenance.get("surface_imbalance", float("nan"))),
    ]
    # Moist states only; a dry file has none of these, so its block is
    # unchanged. The surface humidity is part of the configuration (page 9's
    # stepper), and the settled TOA is the residual the moist column keeps.
    if "surface_relative_humidity" in provenance:
        lines.insert(-2, "  surface RH     {:.0f} %, at the surface's own "
                     "temperature".format(
                         100 * provenance["surface_relative_humidity"]))
    if "window_mean_toa_w_m2" in provenance:
        lines.append("  settled TOA    {:+.3f} W/m^2, mean over the last "
                     "{} steps".format(provenance["window_mean_toa_w_m2"],
                                       provenance.get("steady_window_steps")))
    if "perturbed_from" in provenance:
        saved = provenance.get("perturbed_from_saved_at")
        lines.append("  perturbed from {}{}".format(
            provenance["perturbed_from"],
            " (saved {})".format(saved) if saved else ""))
    return "\n".join(lines)


def _exists(path):
    import os

    return os.path.isfile(path)
