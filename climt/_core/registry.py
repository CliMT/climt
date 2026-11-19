"""List climt's components, what is usable now, and what to install for the rest."""
from importlib.metadata import entry_points
from typing import NamedTuple

import sympl

from .backend import HAS_NUMBA
from .fortran import fortran_status

ENTRY_POINT_GROUP = "climt.components"

FORTRAN_BACKED = {
    "DcmipInitialConditions": "_dcmip",
    "EmanuelConvection": "_emanuel_convection",
    "RRTMGLongwave": "_rrtmg_lw",
    "RRTMGShortwave": "_rrtmg_sw",
    "SimplePhysics": "_simple_physics",
}

_KINDS = ("Stepper", "TendencyComponent", "ImplicitTendencyComponent", "DiagnosticComponent")
_COLUMNS = ("name", "kind", "available", "requires", "source")


class ComponentInfo(NamedTuple):
    name: str
    kind: str
    available: bool
    requires: str
    source: str


class ComponentList(list):
    """A list of ComponentInfo. repr() is an aligned text table followed by one
    line about numba."""

    def __repr__(self):
        rows = [_COLUMNS] + [tuple(str(field) for field in info) for info in self]
        widths = [max(len(row[i]) for row in rows) for i in range(len(_COLUMNS))]
        lines = ["  ".join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip() for row in rows]
        if HAS_NUMBA:
            lines.append("numba: installed")
        else:
            lines.append('numba: not installed (kernels run as plain Python; pip install "climt[numba]")')
        return "\n".join(lines)


def _entry_points():
    return entry_points(group=ENTRY_POINT_GROUP)


def _kind(cls):
    for base in cls.__mro__:
        for kind in _KINDS:
            if base is getattr(sympl, kind, None):
                return kind
    return "other"


def list_components():
    """Return a ComponentList of climt's components and any registered through entry points.

    Fortran-backed components are marked unavailable, with requires="climt-fortran",
    when their extension cannot be loaded. Entry points are listed, not loaded.
    """
    from .. import _components

    modules = fortran_status()["modules"]
    result = ComponentList()
    for cls in sorted(_components.__all__, key=lambda c: c.__name__):
        name = cls.__name__
        available = modules.get(FORTRAN_BACKED[name], False) if name in FORTRAN_BACKED else True
        requires = "" if available else "climt-fortran"
        result.append(ComponentInfo(name, _kind(cls), bool(available), requires, "climt"))
    for ep in sorted(_entry_points(), key=lambda e: e.name):
        dist = getattr(ep, "dist", None)
        source = "{} ({})".format(dist.name, ep.value) if dist is not None else ep.value
        result.append(ComponentInfo(ep.name, "external", True, "", source))
    return result
