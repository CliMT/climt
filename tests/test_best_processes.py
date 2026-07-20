import numpy as np
from climt._components.second_best.processes.soil_properties import BestSoilProperties


def test_soil_properties_sand_and_clay():
    sp = BestSoilProperties()
    clay = sp("clay", "land")
    sand = sp("sand", "land")
    # clay: tex=0 -> porosity 0.6, X_FC 0.57, colour 0.2
    assert np.isclose(clay["texture"], 0.0)
    assert np.isclose(clay["porosity"], 0.6)
    assert np.isclose(clay["field_capacity"], 0.95 * 0.6)
    assert np.isclose(clay["colour"], 0.2)
    # sand: tex=9 -> porosity 0.33, X_FC (0.95-0.774)*0.33, colour 1.0
    assert np.isclose(sand["texture"], 9.0)
    assert np.isclose(sand["porosity"], 0.6 - 0.03 * 9)
    assert np.isclose(sand["field_capacity"], (0.95 - 0.086 * 9) * (0.6 - 0.03 * 9))
    assert np.isclose(sand["colour"], 1.0)
    # land_ice override on texture
    ice = sp("clay", "land_ice")
    assert np.isclose(ice["texture"], 0.07)
    assert np.isclose(ice["wilting_point"], 0.01)


from climt._components.second_best.processes.albedo import BestSurfaceAlbedo


def test_albedo_bare_soil_and_land_ice():
    alb = BestSurfaceAlbedo()
    sp = BestSoilProperties()
    clay = sp("clay", "land")
    # wet clay (W=1): alpha_sw = 0.10 + 0.1*0.2 + 0.06*0 = 0.12; alpha_lw = 0.24
    out = alb(clay, wetness=1.0, area_type="land")
    assert np.isclose(out["alpha_sw"], 0.12)
    assert np.isclose(out["alpha_lw"], 0.24)
    # dry sand (W=0): 0.10 + 0.1*1.0 + 0.06*1 = 0.26
    sand = sp("sand", "land")
    out2 = alb(sand, wetness=0.0, area_type="land")
    assert np.isclose(out2["alpha_sw"], 0.26)
    # land ice (W=0): alpha_sw = 0.60 + 0.06 = 0.66; alpha_lw = alpha_sw/3
    out3 = alb(sp("clay", "land_ice"), wetness=0.0, area_type="land_ice")
    assert np.isclose(out3["alpha_sw"], 0.66)
    assert np.isclose(out3["alpha_lw"], 0.66 / 3.0)
