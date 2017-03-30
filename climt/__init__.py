# -*- coding: utf-8 -*-

try:
    import sympl

    from ._core import (
        get_default_state, mass_to_volume_mixing_ratio,
        get_input_arrays_from_state)

    from ._components import (
        Frierson06LongwaveOpticalDepth, GrayLongwaveRadiation, HeldSuarez,
        GridScaleCondensation, BergerSolarInsolation, SimplePhysics, RRTMLongwave)

    sympl.set_dimension_names(z=['interface_levels', 'mid_levels'])
    __all__ = (
        sympl, get_default_state, mass_to_volume_mixing_ratio, get_input_arrays_from_state,
        Frierson06LongwaveOpticalDepth, GrayLongwaveRadiation, HeldSuarez,
        GridScaleCondensation, BergerSolarInsolation,
        SimplePhysics, RRTMLongwave
    )
except:
    print('Could not import climt components and core!')

__version__ = '1.0.0'
