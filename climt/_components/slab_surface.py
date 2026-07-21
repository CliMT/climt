import numpy as np
from sympl import TendencyComponent, get_constant

from .._core.backend import jit_compile, prange
from .._core.horizontal_operators import curl_z, divergence

# Area type mapping for JIT compatibility
# land: 0, land_ice: 1, sea: 2, sea_ice: 3
AREA_MAP = {"land": 0, "land_ice": 1, "sea": 2, "sea_ice": 3}


class SlabSurface(TendencyComponent):
    """
    Calculate the surface energy balance.

    This component assumes the surface is a slab of possibly
    varying heat capacity, and calculates the surface temperature.
    """

    input_properties = {
        "downwelling_longwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "downwelling_shortwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "upwelling_longwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "upwelling_shortwave_flux_in_air": {
            "dims": ["*", "interface_levels"],
            "units": "W m^-2",
        },
        "surface_upward_latent_heat_flux": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
        "surface_upward_sensible_heat_flux": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "surface_thermal_capacity": {
            "dims": ["*"],
            "units": "J kg^-1 degK^-1",
        },
        "surface_material_density": {
            "dims": ["*"],
            "units": "kg m^-3",
        },
        "upward_heat_flux_at_ground_level_in_soil": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "heat_flux_into_sea_water_due_to_sea_ice": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "area_type": {
            "dims": ["*"],
            "units": "dimensionless",
        },
        "soil_layer_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "ocean_mixed_layer_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "heat_capacity_of_soil": {
            "dims": ["*"],
            "units": "J kg^-1 degK^-1",
        },
        "sea_water_density": {
            "dims": ["*"],
            "units": "kg m^-3",
        },
        "ocean_heat_transport_convergence": {
            "dims": ["*"],
            "units": "W m^-2",
        },
    }

    tendency_properties = {
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK s^-1",
        }
    }

    diagnostic_properties = {
        "depth_of_slab_surface": {
            "dims": ["*"],
            "units": "m",
        },
        "ocean_heat_transport_convergence": {
            "dims": ["*"],
            "units": "W m^-2",
        },
    }

    def __init__(self, include_ekman=False, equatorial_ekman_cap_latitude=5.0, **kwargs):
        """
        Args:
            include_ekman (bool, optional): If True, additionally compute an
                Ekman heat-transport convergence from wind-stress curl and
                fold it into the ocean heat-transport q-flux. Requires
                ``surface_downward_eastward_stress``,
                ``surface_downward_northward_stress``, ``latitude`` and
                ``longitude`` as 2-D (lat, lon) inputs. Off by default for
                backward compatibility.
            equatorial_ekman_cap_latitude (float, optional): Latitude (degrees)
                below which the Coriolis parameter is capped in magnitude, to
                avoid a singularity in the Ekman transport at the equator.
        """
        self._include_ekman = include_ekman
        self._eq_cap = equatorial_ekman_cap_latitude

        if include_ekman:
            # Instance-level copies so that include_ekman=False (the default)
            # leaves the class-level dictionaries, and therefore this
            # component's declared inputs/outputs, completely unchanged.
            self.input_properties = dict(self.input_properties)
            self.input_properties.update(
                {
                    "surface_downward_eastward_stress": {
                        "dims": ["lat", "lon"],
                        "units": "N m^-2",
                    },
                    "surface_downward_northward_stress": {
                        "dims": ["lat", "lon"],
                        "units": "N m^-2",
                    },
                    "latitude": {
                        "dims": ["lat", "lon"],
                        "units": "degrees_north",
                    },
                    "longitude": {
                        "dims": ["lat", "lon"],
                        "units": "degrees_east",
                    },
                }
            )
            self.diagnostic_properties = dict(self.diagnostic_properties)
            self.diagnostic_properties.update(
                {
                    "ekman_heat_transport_convergence": {
                        "dims": ["*"],
                        "units": "W m^-2",
                    },
                    "ekman_pumping": {
                        "dims": ["*"],
                        "units": "m s^-1",
                    },
                }
            )

        super(SlabSurface, self).__init__(**kwargs)

    def array_call(self, state):
        sw_down_raw = state["downwelling_shortwave_flux_in_air"]

        area_type_raw = state["area_type"]

        # Convert area_type to integer codes
        # land: 0, land_ice: 1, sea: 2, sea_ice: 3
        # Standard climt state uses strings/bytes for area_type.
        if hasattr(area_type_raw, "astype"):
            area_type_str = area_type_raw.astype(str)
        else:
            area_type_str = np.asarray(area_type_raw).astype(str)

        area_type_code = np.zeros(area_type_str.shape, dtype=np.int32)
        for k, v in AREA_MAP.items():
            area_type_code[area_type_str == k] = v

        # Flatten inputs for kernel
        def flat(x):
            return np.reshape(x, (-1,))

        lw_down_raw = state["downwelling_longwave_flux_in_air"]
        sw_up_raw = state["upwelling_shortwave_flux_in_air"]
        lw_up_raw = state["upwelling_longwave_flux_in_air"]

        def get_surface(arr):
            data = getattr(arr, "data", arr)
            data = np.asarray(data)

            if data.ndim > 1:
                surf = data[..., 0]
            else:
                surf = data

            return flat(surf)

        sw_down = get_surface(sw_down_raw)
        lw_down = get_surface(lw_down_raw)
        sw_up = get_surface(sw_up_raw)
        lw_up = get_surface(lw_up_raw)

        lh = flat(
            getattr(
                state["surface_upward_latent_heat_flux"],
                "data",
                state["surface_upward_latent_heat_flux"],
            )
        )
        sh = flat(
            getattr(
                state["surface_upward_sensible_heat_flux"],
                "data",
                state["surface_upward_sensible_heat_flux"],
            )
        )

        up_heat_soil = flat(
            getattr(
                state["upward_heat_flux_at_ground_level_in_soil"],
                "data",
                state["upward_heat_flux_at_ground_level_in_soil"],
            )
        )
        heat_flux_sea_ice = flat(
            getattr(
                state["heat_flux_into_sea_water_due_to_sea_ice"],
                "data",
                state["heat_flux_into_sea_water_due_to_sea_ice"],
            )
        )

        sea_water_dens = flat(
            getattr(state["sea_water_density"], "data", state["sea_water_density"])
        )
        surf_dens = flat(
            getattr(
                state["surface_material_density"],
                "data",
                state["surface_material_density"],
            )
        )

        heat_cap_soil = flat(
            getattr(
                state["heat_capacity_of_soil"], "data", state["heat_capacity_of_soil"]
            )
        )
        surf_therm_cap = flat(
            getattr(
                state["surface_thermal_capacity"],
                "data",
                state["surface_thermal_capacity"],
            )
        )

        ocean_mix_thick = flat(
            getattr(
                state["ocean_mixed_layer_thickness"],
                "data",
                state["ocean_mixed_layer_thickness"],
            )
        )
        soil_layer_thick = flat(
            getattr(
                state["soil_layer_thickness"], "data", state["soil_layer_thickness"]
            )
        )

        ocean_heat_transport_raw = state["ocean_heat_transport_convergence"]
        ocean_heat_transport = flat(
            getattr(
                ocean_heat_transport_raw, "data", ocean_heat_transport_raw
            )
        )

        # Ekman transport terms are computed here, on the 2-D (lat, lon)
        # grid, before everything gets flattened to "*" for the kernel.
        # They are folded additively into the Task-7 ocean_heat_transport
        # q-flux that is handed to the (unchanged) kernel, so the kernel
        # signature does not need to change.
        total_ocean_heat_transport = ocean_heat_transport
        ekman_heat_transport_flat = None
        ekman_pumping_flat = None

        if self._include_ekman:
            lat2d = np.asarray(
                getattr(state["latitude"], "data", state["latitude"]),
                dtype=float,
            )
            lon2d = np.asarray(
                getattr(state["longitude"], "data", state["longitude"]),
                dtype=float,
            )
            if lat2d.ndim == 1:
                # Defensive fallback for a flattened single-column view:
                # both horizontal_operators return zeros when a dimension
                # has < 3 points, so Q_ekman/w_ek collapse to 0 here.
                lat2d = lat2d.reshape(-1, 1)
                lon2d = lon2d.reshape(-1, 1)

            tau_x = np.asarray(
                getattr(
                    state["surface_downward_eastward_stress"],
                    "data",
                    state["surface_downward_eastward_stress"],
                ),
                dtype=float,
            ).reshape(lat2d.shape)
            tau_y = np.asarray(
                getattr(
                    state["surface_downward_northward_stress"],
                    "data",
                    state["surface_downward_northward_stress"],
                ),
                dtype=float,
            ).reshape(lat2d.shape)

            surf_temp_flat = flat(
                getattr(
                    state["surface_temperature"],
                    "data",
                    state["surface_temperature"],
                )
            )
            theta_2d = np.asarray(surf_temp_flat, dtype=float).reshape(lat2d.shape)

            omega = get_constant("planetary_rotation_rate", "s^-1")
            c_sw = get_constant("heat_capacity_of_sea_water", "J/kg/degK")

            f = 2.0 * omega * np.sin(np.deg2rad(lat2d))
            f_floor = 2.0 * omega * np.sin(np.deg2rad(self._eq_cap))
            # np.sign(0.0) == 0.0, which would zero out f_capped exactly at
            # the equator; use a sign that never vanishes instead.
            f_sign = np.where(f >= 0.0, 1.0, -1.0)
            f_capped = f_sign * np.maximum(np.abs(f), f_floor)

            # Ekman mass transport per unit width (kg m^-1 s^-1).
            Mx = tau_y / f_capped
            My = -tau_x / f_capped

            # Ekman pumping w_ek = curl_z(tau) / f: take the curl of the
            # raw wind stress first, then divide by f. Dividing by f before
            # differentiating (i.e. curl_z(tau/f, ...)) would spuriously
            # pick up the meridional variation of f itself even for a
            # spatially uniform wind stress, which has zero physical curl.
            w_ek = curl_z(tau_x, tau_y, lat2d, lon2d) / f_capped
            q_ekman_2d = -c_sw * divergence(
                theta_2d * Mx, theta_2d * My, lat2d, lon2d
            )

            # Only open-ocean (sea, not sea-ice) cells get the Ekman
            # convergence added to their q-flux, matching the Task-7
            # `sea_mask and not sea_ice_mask` convention.
            open_ocean_mask_2d = (
                area_type_code.reshape(lat2d.shape) == AREA_MAP["sea"]
            )
            q_ekman_2d = np.where(open_ocean_mask_2d, q_ekman_2d, 0.0)

            ekman_heat_transport_flat = flat(q_ekman_2d)
            ekman_pumping_flat = flat(w_ek)

            total_ocean_heat_transport = ocean_heat_transport + ekman_heat_transport_flat

        tend_ts, depth = _slab_surface_kernel_np(
            sw_down,
            lw_down,
            sw_up,
            lw_up,
            lh,
            sh,
            flat(area_type_code),
            up_heat_soil,
            heat_flux_sea_ice,
            sea_water_dens,
            surf_dens,
            heat_cap_soil,
            surf_therm_cap,
            ocean_mix_thick,
            soil_layer_thick,
            total_ocean_heat_transport,
        )

        diagnostics = {
            "depth_of_slab_surface": np.reshape(depth, area_type_raw.shape),
            "ocean_heat_transport_convergence": np.reshape(
                ocean_heat_transport, area_type_raw.shape
            ),
        }
        if self._include_ekman:
            diagnostics["ekman_heat_transport_convergence"] = np.reshape(
                ekman_heat_transport_flat, area_type_raw.shape
            )
            diagnostics["ekman_pumping"] = np.reshape(
                ekman_pumping_flat, area_type_raw.shape
            )

        return {
            "surface_temperature": np.reshape(tend_ts, area_type_raw.shape)
        }, diagnostics


@jit_compile(backend=np)
def _slab_surface_kernel_np(
    sw_down,
    lw_down,
    sw_up,
    lw_up,
    lh,
    sh,
    area_type,
    up_heat_soil,
    heat_flux_sea_ice,
    sea_water_dens,
    surf_dens,
    heat_cap_soil,
    surf_therm_cap,
    ocean_mix_thick,
    soil_layer_thick,
    ocean_heat_transport,
):

    ncol = area_type.size
    tend_ts = np.zeros(ncol)
    depth = np.zeros(ncol)

    for i in prange(ncol):
        net_heat_flux = sw_down[i] + lw_down[i] - sw_up[i] - lw_up[i] - sh[i] - lh[i]

        at = area_type[i]
        land_mask = (at == 0) or (at == 1)
        sea_mask = (at == 2) or (at == 3)
        land_ice_mask = at == 1
        sea_ice_mask = at == 3

        if land_ice_mask:
            net_heat_flux = -up_heat_soil[i]
        elif sea_ice_mask:
            net_heat_flux = heat_flux_sea_ice[i]

        if sea_mask and not sea_ice_mask:
            net_heat_flux = net_heat_flux + ocean_heat_transport[i]

        if sea_mask:
            final_dens = sea_water_dens[i]
            d = ocean_mix_thick[i]
        else:
            final_dens = surf_dens[i]
            d = soil_layer_thick[i] if land_mask else 0.0

        if land_mask:
            final_therm_cap = heat_cap_soil[i]
        else:
            final_therm_cap = surf_therm_cap[i]

        depth[i] = d
        mass_slab = final_dens * d
        heat_cap_slab = mass_slab * final_therm_cap

        if heat_cap_slab != 0:
            val = net_heat_flux / heat_cap_slab
        else:
            val = 0.0

        if land_ice_mask or sea_ice_mask:
            val = 0.0

        tend_ts[i] = val

    return tend_ts, depth
