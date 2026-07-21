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


from climt._components.second_best.processes.subsurface import BestSubsurfaceTransport


def test_subsurface_freezing_creates_ice_and_warms_toward_freezing():
    st = BestSubsurfaceTransport()
    n = 6
    profiles = {
        "T": np.full(n, 270.0),        # below freezing (273)
        "X_w": np.full(n, 0.2),        # liquid water present
        "X_i": np.zeros(n),
    }
    out = st(profiles, surface_flux_bc=-20.0, timestep=600.0, dz=0.1)
    assert np.all(out["X_i"] >= 0.0)
    assert out["X_i"].sum() > 0.0          # some water froze
    assert out["X_w"].sum() < profiles["X_w"].sum()  # liquid decreased
    assert np.all(out["T"] <= 273.0 + 1e-6)          # cannot exceed freezing here
    assert not np.any(np.isnan(out["T"]))


def test_subsurface_conserves_total_water_mass():
    st = BestSubsurfaceTransport()
    n = 8
    profiles = {"T": np.full(n, 268.0), "X_w": np.full(n, 0.15),
                "X_i": np.full(n, 0.05)}
    total0 = profiles["X_w"].sum() + profiles["X_i"].sum()
    out = st(profiles, surface_flux_bc=-10.0, timestep=300.0, dz=0.1)
    total1 = out["X_w"].sum() + out["X_i"].sum()
    assert np.isclose(total0, total1, atol=1e-9)   # phase change moves mass, conserves it


def test_interpolate_neutral_reduces_to_log_law():
    sl = BestSurfaceLayer()
    z0, z_mid = 0.01, 69.0
    cdn = (0.4 / np.log(z_mid / z0)) ** 2
    neutral = {"C_Dm": cdn, "C_Dh": cdn, "C_DN": cdn, "Ri": 0.0}
    # scalar: weight is the neutral log ratio
    w = np.log(2.0 / z0) / np.log(z_mid / z0)
    t2m = sl.interpolate_to_height(neutral, z0, z_mid, 2.0,
                                   surface_value=300.0, level_value=290.0,
                                   kind="scalar")
    assert np.isclose(t2m, 300.0 + (290.0 - 300.0) * w)
    # wind: neutral log ratio of the lowest-level speed
    u10 = sl.interpolate_to_height(neutral, z0, z_mid, 10.0,
                                   surface_value=0.0, level_value=8.0,
                                   kind="wind")
    assert np.isclose(u10, 8.0 * np.log(10.0 / z0) / np.log(z_mid / z0))


def test_interpolate_bounded_and_stability_direction():
    sl = BestSurfaceLayer()
    z0, z_mid = 0.01, 69.0
    cdn = (0.4 / np.log(z_mid / z0)) ** 2
    neutral = {"C_Dm": cdn, "C_Dh": cdn, "C_DN": cdn, "Ri": 0.0}
    unstable = {"C_Dm": cdn * 1.5, "C_Dh": cdn * 1.8, "C_DN": cdn, "Ri": -0.3}
    stable = {"C_Dm": cdn * 0.7, "C_Dh": cdn * 0.6, "C_DN": cdn, "Ri": 0.15}
    Ts, Ta = 300.0, 290.0
    t_neu = sl.interpolate_to_height(neutral, z0, z_mid, 2.0, Ts, Ta, "scalar")
    t_uns = sl.interpolate_to_height(unstable, z0, z_mid, 2.0, Ts, Ta, "scalar")
    t_sta = sl.interpolate_to_height(stable, z0, z_mid, 2.0, Ts, Ta, "scalar")
    # all bounded between surface and air
    for t in (t_neu, t_uns, t_sta):
        assert Ta <= t <= Ts
    # unstable = stronger mixing -> 2m closer to the (well-mixed) air value;
    # stable = strong near-surface gradient -> 2m closer to the surface value
    assert t_uns < t_neu < t_sta


def test_interpolate_wind_not_exceeding_level_under_instability():
    sl = BestSurfaceLayer()
    z0, z_mid = 0.01, 30.0
    # mild instability produced by the scheme itself (light wind, warm surface)
    drag = sl(z_mid, z0, wind_speed=1.0, T_surf=295.0, T_air=290.0,
              area_type="land")
    assert drag["Ri"] < 0.0
    spd10 = sl.interpolate_to_height(drag, z0, z_mid, 10.0,
                                     surface_value=0.0, level_value=1.0,
                                     kind="wind")
    # the 10 m wind cannot exceed the lowest-level wind it is interpolated from
    assert 0.0 < spd10 <= 1.0
