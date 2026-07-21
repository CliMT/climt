import numpy as np
from sympl import Stepper, get_constant, initialize_numpy_arrays_with_properties
from .processes import (BestSoilProperties, BestSurfaceAlbedo, BestSurfaceLayer,
                        BestSurfaceFluxes, BestSubsurfaceTransport)


class SecondBEST(Stepper):
    """Modular intermediate-complexity BEST land surface model.

    A thin orchestrator: every physical calculation is delegated to a
    swappable process object (see ``climt._components.second_best.processes``).
    Pass an instance to override any process; ``None`` selects the BEST default.
    """

    input_properties = {
        "air_temperature": {"dims": ["mid_levels", "*"], "units": "degK"},
        "specific_humidity": {"dims": ["mid_levels", "*"], "units": "kg/kg"},
        "northward_wind": {"dims": ["mid_levels", "*"], "units": "m s^-1"},
        "eastward_wind": {"dims": ["mid_levels", "*"], "units": "m s^-1"},
        "air_pressure": {"dims": ["mid_levels", "*"], "units": "Pa"},
        "downwelling_shortwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "downwelling_longwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "upwelling_shortwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "upwelling_longwave_flux_in_air": {"dims": ["*", "interface_levels"], "units": "W m^-2"},
        "area_type": {"dims": ["*"], "units": "dimensionless"},
        "surface_temperature": {"dims": ["*"], "units": "degK"},
        "surface_air_pressure": {"dims": ["*"], "units": "Pa"},
        "soil_temperature": {"dims": ["soil_interface_levels", "*"], "units": "degK"},
        "soil_liquid_water_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "soil_ice_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "surface_snow_thickness": {"dims": ["*"], "units": "m"},
        "height_on_soil_interface_levels": {"dims": ["soil_interface_levels", "*"], "units": "m"},
    }
    output_properties = {
        "surface_temperature": {"dims": ["*"], "units": "degK"},
        "soil_temperature": {"dims": ["soil_interface_levels", "*"], "units": "degK"},
        "soil_liquid_water_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "soil_ice_content": {"dims": ["soil_interface_levels", "*"], "units": "m^3/m^3"},
        "surface_snow_thickness": {"dims": ["*"], "units": "m"},
    }
    diagnostic_properties = {
        "surface_upward_sensible_heat_flux": {"dims": ["*"], "units": "W m^-2"},
        "surface_upward_latent_heat_flux": {"dims": ["*"], "units": "W m^-2"},
        "evaporation_rate": {"dims": ["*"], "units": "m s^-1"},
        "surface_albedo_for_direct_shortwave": {"dims": ["*"], "units": "dimensionless"},
        "surface_albedo_for_diffuse_shortwave": {"dims": ["*"], "units": "dimensionless"},
        "surface_drag_coefficient_for_heat": {"dims": ["*"], "units": "dimensionless"},
        "surface_drag_coefficient_for_momentum": {"dims": ["*"], "units": "dimensionless"},
        "air_temperature_at_2m": {"dims": ["*"], "units": "degK"},
        "specific_humidity_at_2m": {"dims": ["*"], "units": "kg/kg"},
        "eastward_wind_at_10m": {"dims": ["*"], "units": "m s^-1"},
        "northward_wind_at_10m": {"dims": ["*"], "units": "m s^-1"},
    }

    def __init__(self, soil_type="clay", num_soil_layers=3, minimum_wind_speed=1.0,
                 soil_properties=None, albedo=None, surface_layer=None,
                 fluxes=None, subsurface=None, **kwargs):
        self._soil_type = soil_type
        self._num_soil_layers = num_soil_layers
        self._min_wind = minimum_wind_speed
        self._soil_props = soil_properties or BestSoilProperties()
        self._albedo = albedo or BestSurfaceAlbedo()
        self._surface_layer = surface_layer or BestSurfaceLayer()
        self._fluxes = fluxes or BestSurfaceFluxes()
        self._subsurface = subsurface or BestSubsurfaceTransport()
        super(SecondBEST, self).__init__(**kwargs)

    def array_call(self, state, timestep):
        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties)
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties)

        # Pre-fill every output from the matching input so columns whose
        # area_type is not "land"/"land_ice" (e.g. "sea") pass through
        # unchanged instead of being left at the zero-fill that
        # initialize_numpy_arrays_with_properties applies above. Mirrors the
        # convention in IceSheet.array_call (climt/_components/surface_ice.py).
        for name in self.output_properties:
            outputs[name][:] = state[name]

        Rd = get_constant("gas_constant_of_dry_air", "J/kg/degK")
        g = get_constant("gravitational_acceleration", "m/s^2")
        ncol = state["area_type"].shape[0]
        # timestep.total_seconds() returns a plain float for datetime.timedelta,
        # but a unyt_quantity (units "s") for sympl's UnytTimeDelta when the
        # UnytBackend is active. The Best* process objects do raw scalar
        # arithmetic assuming a dimensionless/plain-float timestep (e.g.
        # BestSurfaceFluxes mixes it with dimensionless floats pulled straight
        # out of the soil_properties dict), so it must be coerced here at the
        # orchestrator boundary rather than inside each process.
        dt = float(timestep.total_seconds())
        for col in range(ncol):
            area = state["area_type"][col].astype(str)
            if area not in ("land", "land_ice"):
                continue
            props = self._soil_props(self._soil_type, area)

            u = state["eastward_wind"][0, col]
            v = state["northward_wind"][0, col]
            wind = max(np.sqrt(u * u + v * v), self._min_wind)
            T_air = state["air_temperature"][0, col]
            p = state["air_pressure"][0, col]
            rho = p / (Rd * T_air)
            # Surface-layer reference height = geometric height of the lowest
            # model level above the surface, from the hypsometric equation
            # z = (Rd*T/g) * ln(p_surface / p_lowest). This is the height the
            # bulk drag profile is anchored to (tens of metres for a typical
            # column), not the atmospheric scale height. Dry temperature is
            # used (consistent with the dry Rd used elsewhere here); the floor
            # guards the degenerate p_surface <= p_lowest case so ln(z_mid) in
            # the drag law stays finite.
            p_surf = state["surface_air_pressure"][col]
            z_mid = (Rd * T_air / g) * np.log(p_surf / p)
            z_mid = max(float(z_mid), 2.0)
            z0 = 0.01 if area == "land" else 0.001

            T_surf = state["surface_temperature"][col]
            drag = self._surface_layer(z_mid, z0, wind, T_surf, T_air, area)

            X_w = state["soil_liquid_water_content"][:, col]
            W_Lu = X_w[-1] / props["porosity"]
            albedo = self._albedo(props, W_Lu, area)

            atmos = {"air_density": rho, "wind_speed": wind,
                     "air_temperature": T_air,
                     "air_specific_humidity": state["specific_humidity"][0, col],
                     "u": u, "v": v}
            q_sat = _saturation_specific_humidity(T_surf, p)
            soil = {"surface_temperature": T_surf,
                    "saturation_specific_humidity": q_sat,
                    "W_Lu": W_Lu, "W_Fu": state["soil_ice_content"][-1, col]
                    / props["porosity"]}
            flux = self._fluxes(drag, atmos, soil, props, dt)

            net = (state["downwelling_shortwave_flux_in_air"][col, 0]
                   + state["downwelling_longwave_flux_in_air"][col, 0]
                   - state["upwelling_shortwave_flux_in_air"][col, 0]
                   - state["upwelling_longwave_flux_in_air"][col, 0]
                   - flux["sensible_heat_flux"] - flux["latent_heat_flux"])

            z = state["height_on_soil_interface_levels"][:, col]
            dz = float(abs(z[1] - z[0])) if z.shape[0] > 1 else 0.5
            new_prof = self._subsurface(
                {"T": state["soil_temperature"][:, col],
                 "X_w": X_w, "X_i": state["soil_ice_content"][:, col]},
                surface_flux_bc=net, timestep=dt, dz=dz)

            outputs["soil_temperature"][:, col] = new_prof["T"]
            outputs["soil_liquid_water_content"][:, col] = new_prof["X_w"]
            outputs["soil_ice_content"][:, col] = new_prof["X_i"]
            outputs["surface_temperature"][col] = new_prof["T"][-1]
            outputs["surface_snow_thickness"][col] = state["surface_snow_thickness"][col]

            diagnostics["surface_upward_sensible_heat_flux"][col] = flux["sensible_heat_flux"]
            diagnostics["surface_upward_latent_heat_flux"][col] = flux["latent_heat_flux"]
            diagnostics["evaporation_rate"][col] = flux["evaporation"]
            diagnostics["surface_albedo_for_direct_shortwave"][col] = albedo["alpha_sw"]
            diagnostics["surface_albedo_for_diffuse_shortwave"][col] = albedo["alpha_sw"]
            diagnostics["surface_drag_coefficient_for_heat"][col] = drag["C_Dh"]
            diagnostics["surface_drag_coefficient_for_momentum"][col] = drag["C_Dm"]

            # Screen/anemometer-level diagnostics, interpolated between the
            # surface and the lowest model level with the surface-layer's own
            # stability profile (see SurfaceLayer.interpolate_to_height). q2m
            # interpolates from the *effective* surface humidity implied by the
            # evaporative beta, since dry soil is not saturated.
            q_air = atmos["air_specific_humidity"]
            beta = flux.get("beta", 1.0)
            q_surf_eff = beta * q_sat + (1.0 - beta) * q_air
            diagnostics["air_temperature_at_2m"][col] = \
                self._surface_layer.interpolate_to_height(
                    drag, z0, z_mid, 2.0, T_surf, T_air, "scalar")
            diagnostics["specific_humidity_at_2m"][col] = \
                self._surface_layer.interpolate_to_height(
                    drag, z0, z_mid, 2.0, q_surf_eff, q_air, "scalar")
            spd10 = self._surface_layer.interpolate_to_height(
                drag, z0, z_mid, 10.0, 0.0, wind, "wind")
            spd = np.sqrt(u * u + v * v)
            if spd > 0.0:
                diagnostics["eastward_wind_at_10m"][col] = spd10 * u / spd
                diagnostics["northward_wind_at_10m"][col] = spd10 * v / spd
        return diagnostics, outputs


def _saturation_specific_humidity(T, p):
    # Bolton (1980) saturation vapour pressure over water, then specific humidity
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    return 0.622 * es / (p - 0.378 * es)
