"""Pluggable process objects for SecondBEST.

Each process is a plain class exposing a documented ``__call__``. The
``Best*`` defaults implement the equations in
``docs/Description_of_SecondBEST.ipynb`` (Pitman et al.). Swap a process by
passing an instance to ``SecondBEST(...)``.

Note: ``soil_properties`` (Task 4), ``albedo`` (Task 5), and
``surface_layer`` (Task 6) exist so far. ``fluxes`` and ``subsurface`` land
in Tasks 7-8 and should be imported/exported here as each module is added.
"""
from .soil_properties import SoilProperties, BestSoilProperties
from .albedo import SurfaceAlbedo, BestSurfaceAlbedo
from .surface_layer import SurfaceLayer, BestSurfaceLayer

__all__ = [
    "SoilProperties", "BestSoilProperties",
    "SurfaceAlbedo", "BestSurfaceAlbedo",
    "SurfaceLayer", "BestSurfaceLayer",
]
