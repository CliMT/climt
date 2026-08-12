# -*- coding: utf-8 -*-
import sympl

from ._components import (
    BergerSolarInsolation,
    BucketHydrology,
    CorkLongwaveRadiation,
    CorkShortwaveRadiation,
    DataOcean,
    DcmipInitialConditions,
    DryConvectiveAdjustment,
    EmanuelConvection,
    EmanuelConvectionPython,
    Frierson06LongwaveOpticalDepth,
    GrayLongwaveRadiation,
    GridScaleCondensation,
    HeldSuarez,
    IceSheet,
    Instellation,
    LandIce,
    LandMask,
    RRTMGLongwave,
    RRTMGShortwave,
    SeaIce,
    SecondBEST,
    SimpleBoundaryLayer,
    SimplePhysics,
    SlabSurface,
)
from ._core import (
    ConstantNotFoundError,
    UnytBackend,
    UnytStateContainer,
    UnytTimeDelta,
    bolton_dqsat_dT,
    bolton_q_sat,
    calculate_q_sat,
    get_constant_checked,
    get_default_state,
    get_grid,
    get_interface_values,
    list_available_constants,
    load_atmospheric_properties,
    reset_atmospheric_properties,
    mass_to_volume_mixing_ratio,
    numpy_version_of,
    set_constants_from_dict,
)


sympl.set_constant("top_of_model_pressure", 20.0, "Pa")
sympl.set_constant("von_karman_constant", 0.4, "dimensionless")
sympl.set_constant("heat_capacity_of_sea_water", 3985.0, "J/kg/degK")

__all__ = (
    get_default_state,
    get_grid,
    mass_to_volume_mixing_ratio,
    numpy_version_of,
    get_interface_values,
    list_available_constants,
    load_atmospheric_properties,
    reset_atmospheric_properties,
    set_constants_from_dict,
    get_constant_checked,
    ConstantNotFoundError,

    bolton_q_sat,
    bolton_dqsat_dT,
    calculate_q_sat,
    UnytBackend,
    UnytStateContainer,
    UnytTimeDelta,
    Frierson06LongwaveOpticalDepth,
    GrayLongwaveRadiation,
    HeldSuarez,
    GridScaleCondensation,
    BergerSolarInsolation,
    SimplePhysics,
    RRTMGLongwave,
    RRTMGShortwave,
    CorkLongwaveRadiation,
    CorkShortwaveRadiation,
    EmanuelConvection,
    EmanuelConvectionPython,
    SlabSurface,
    DcmipInitialConditions,
    IceSheet,
    Instellation,
    LandMask,
    DryConvectiveAdjustment,
    BucketHydrology,
    SecondBEST,
    SeaIce,
    LandIce,
    DataOcean,
    SimpleBoundaryLayer,
)


def has_fortran_extensions():
    """Return True if compiled Fortran extensions are available."""
    try:
        from climt._components.simple_physics import _simple_physics  # noqa: F401
        return True
    except ImportError:
        return False


__version__ = "0.30.0"
