"""Pluggable process objects for SecondBEST.

Each process is a plain class exposing a documented ``__call__``. The
``Best*`` defaults implement the equations in
``docs/Description_of_SecondBEST.ipynb`` (Pitman et al.). Swap a process by
passing an instance to ``SecondBEST(...)``.

Note: only ``soil_properties`` exists so far (Task 4). ``albedo``,
``surface_layer``, ``fluxes``, and ``subsurface`` land in Tasks 5-8 and
should be imported/exported here as each module is added.
"""
from .soil_properties import SoilProperties, BestSoilProperties

__all__ = [
    "SoilProperties", "BestSoilProperties",
]
