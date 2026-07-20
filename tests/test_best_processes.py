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


from climt._components.second_best.processes.surface_layer import BestSurfaceLayer


def test_surface_layer_neutral_matches_neutral_drag():
    sl = BestSurfaceLayer()
    z_mid, z0 = 10.0, 0.01
    out = sl(z_mid, z0, wind_speed=5.0, T_surf=300.0, T_air=300.0,
             area_type="land")   # neutral: T_surf == T_air -> Ri = 0
    kappa = 0.4
    C_DN = (kappa / (np.log(z_mid) - np.log(z0))) ** 2
    assert np.isclose(out["Ri"], 0.0)
    assert np.isclose(out["C_DN"], C_DN)
    # at Ri=0 the stable branch reduces to C_DN
    assert np.isclose(out["C_Dm"], C_DN)
    assert np.isclose(out["C_Dh"], C_DN)


def test_surface_layer_unstable_increases_drag():
    sl = BestSurfaceLayer()
    out = sl(10.0, 0.01, wind_speed=5.0, T_surf=290.0, T_air=300.0,
             area_type="land")   # T_surf < T_air is stable; flip for unstable
    unstable = sl(10.0, 0.01, wind_speed=5.0, T_surf=310.0, T_air=300.0,
                  area_type="land")
    assert unstable["Ri"] < 0.0
    assert unstable["C_Dh"] > unstable["C_DN"]   # unstable enhances exchange


def test_surface_layer_drag_positive_over_ri_sweep():
    sl = BestSurfaceLayer()
    for dT in np.linspace(-20, 20, 41):
        out = sl(10.0, 0.01, 5.0, 300.0 + dT, 300.0, "land")
        assert out["C_Dm"] > 0.0 and np.isfinite(out["C_Dm"])
        assert out["C_Dh"] > 0.0 and np.isfinite(out["C_Dh"])


from climt._components.second_best.processes.fluxes import BestSurfaceFluxes


def test_sensible_heat_flux_bulk_formula():
    from sympl import get_constant
    fx = BestSurfaceFluxes()
    C_Dh = 0.003354
    drag = {"C_Dm": C_Dh, "C_Dh": C_Dh, "C_DN": C_Dh, "Ri": 0.0}
    atmos = {"air_density": 1.2, "wind_speed": 5.0, "air_temperature": 295.0,
             "air_specific_humidity": 0.01, "u": 5.0, "v": 0.0}
    soil = {"surface_temperature": 300.0, "saturation_specific_humidity": 0.02,
            "W_Lu": 1.0, "W_Fu": 0.0}
    out = fx(drag, atmos, soil, soil_props={"porosity": 0.6}, timestep=100.0)
    Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK")
    expected_shf = 1.2 * Cpd * 5.0 * C_Dh * (300.0 - 295.0)
    assert np.isclose(out["sensible_heat_flux"], expected_shf, rtol=1e-6)
    assert out["latent_heat_flux"] > 0.0   # surface wetter than air


def test_latent_flux_zero_when_no_humidity_gradient():
    fx = BestSurfaceFluxes()
    drag = {"C_Dm": 0.003, "C_Dh": 0.003, "C_DN": 0.003, "Ri": 0.0}
    atmos = {"air_density": 1.2, "wind_speed": 5.0, "air_temperature": 300.0,
             "air_specific_humidity": 0.02, "u": 5.0, "v": 0.0}
    soil = {"surface_temperature": 300.0, "saturation_specific_humidity": 0.02,
            "W_Lu": 1.0, "W_Fu": 0.0}
    out = fx(drag, atmos, soil, {"porosity": 0.6}, 100.0)
    assert np.isclose(out["latent_heat_flux"], 0.0)
