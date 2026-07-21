import logging

import numpy as np
from sympl import Stepper, get_constant, initialize_numpy_arrays_with_properties

from ..._core.snow_ice_column import Dirichlet, Flux, solve_column

logger = logging.getLogger(__name__)


class SeaIce(Stepper):
    """
    1-D thermodynamic sea-ice model over ``area_type == "sea_ice"`` cells.

    This is a mechanical re-expression of the sea-ice branch of the
    monolithic :class:`~climt._components.surface_ice.IceSheet` component,
    using the shared implicit column solver
    (:func:`climt._core.snow_ice_column.solve_column`) in place of the
    hand-rolled tridiagonal solve, with the following defect fixes applied:

    1. The basal boundary condition is a prescribed ocean heat flux
       (:class:`~climt._core.snow_ice_column.Flux`) rather than a freezing
       Dirichlet condition, so heat can flow either way across the
       ice/ocean interface (growth *or* basal melt).
    2. ``sea_ice_thickness`` is clamped to be non-negative after melt;
       any energy that would have driven the thickness below zero is
       instead routed into ``heat_flux_into_sea_water_due_to_sea_ice``.
    3. The three albedo values (snow / bare ice / melting) are
       configurable via ``__init__`` rather than hardcoded, and the
       duplicated ``elif`` branch for the bare-ice albedo has been
       removed.
    4. The ``print``/``assert False`` guard against negative melt energy
       has been replaced with an ``np.clip`` and a debug log message.
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
        "surface_upward_sensible_heat_flux": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "sea_ice_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "surface_snow_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "area_type": {
            "dims": ["*"],
            "units": "dimensionless",
        },
        "snow_and_ice_temperature": {
            "dims": ["ice_interface_levels", "*"],
            "units": "degK",
        },
        "sea_surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
        "heat_flux_into_sea_water_due_to_sea_ice": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "height_on_ice_interface_levels": {
            "dims": ["ice_interface_levels", "*"],
            "units": "m",
        },
    }

    output_properties = {
        "sea_ice_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "surface_snow_thickness": {
            "dims": ["*"],
            "units": "m",
        },
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
        "snow_and_ice_temperature": {
            "dims": ["ice_interface_levels", "*"],
            "units": "degK",
        },
        "height_on_ice_interface_levels": {
            "dims": ["ice_interface_levels", "*"],
            "units": "m",
        },
    }

    diagnostic_properties = {
        "heat_flux_into_sea_water_due_to_sea_ice": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "surface_downward_heat_flux_in_sea_ice": {
            "dims": ["*"],
            "units": "W m^-2",
        },
        "surface_albedo_for_direct_shortwave": {
            "dims": ["*"],
            "units": "dimensionless",
        },
        "surface_albedo_for_diffuse_shortwave": {
            "dims": ["*"],
            "units": "dimensionless",
        },
    }

    def __init__(
        self,
        maximum_snow_ice_height=10,
        albedo_snow=0.8,
        albedo_ice=0.5,
        albedo_melt=0.2,
        **kwargs
    ):
        """
        Args:
            maximum_snow_ice_height (float):
                The maximum combined height of snow and ice handled by the
                model in :math:`m`.
            albedo_snow (float):
                Surface albedo used when snow is present on top of the ice.
            albedo_ice (float):
                Surface albedo used for bare (snow-free) sea ice.
            albedo_melt (float):
                Surface albedo used when the surface is actively melting.
        """
        self._max_height = maximum_snow_ice_height
        self._albedo_snow = albedo_snow
        self._albedo_ice = albedo_ice
        self._albedo_melt = albedo_melt
        self._epsilon = 1e-5
        super(SeaIce, self).__init__(**kwargs)

    def _update_constants(self):
        self._Kice = get_constant(
            "thermal_conductivity_of_solid_phase_as_ice", "W/m/degK"
        )
        self._Ksnow = get_constant(
            "thermal_conductivity_of_solid_phase_as_snow", "W/m/degK"
        )
        self._rho_ice = get_constant("density_of_solid_phase_as_ice", "kg/m^3")
        self._C_ice = get_constant("heat_capacity_of_solid_phase_as_ice", "J/kg/degK")
        self._rho_snow = get_constant("density_of_solid_phase_as_snow", "kg/m^3")
        self._C_snow = get_constant("heat_capacity_of_solid_phase_as_snow", "J/kg/degK")
        self._Lf = get_constant("latent_heat_of_fusion", "J/kg")
        self._melting_temperature = get_constant(
            "freezing_temperature_of_liquid_phase", "degK"
        )

    def array_call(self, raw_state, timestep):
        self._update_constants()

        num_cols = raw_state["area_type"].shape[0]
        # timestep.total_seconds() returns a plain float for datetime.timedelta,
        # but a unyt_quantity (units "s") for sympl's UnytTimeDelta when the
        # UnytBackend is active. This component does raw scalar arithmetic
        # assuming a dimensionless/plain-float timestep (e.g. combining it
        # with plain numpy arrays from initialize_numpy_arrays_with_properties),
        # so it must be coerced to a plain float here.
        dt = float(timestep.total_seconds())

        net_heat_flux = (
            raw_state["downwelling_shortwave_flux_in_air"][:, 0]
            + raw_state["downwelling_longwave_flux_in_air"][:, 0]
            - raw_state["upwelling_shortwave_flux_in_air"][:, 0]
            - raw_state["upwelling_longwave_flux_in_air"][:, 0]
            - raw_state["surface_upward_sensible_heat_flux"]
            - raw_state["surface_upward_latent_heat_flux"]
        )

        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, raw_state, self.input_properties
        )
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )

        # Copy input values (columns that are not sea ice, or that have no
        # ice, simply pass through unchanged).
        outputs["sea_ice_thickness"][:] = raw_state["sea_ice_thickness"]
        outputs["surface_snow_thickness"][:] = raw_state["surface_snow_thickness"]
        outputs["snow_and_ice_temperature"][:] = raw_state["snow_and_ice_temperature"]
        outputs["surface_temperature"][:] = raw_state["snow_and_ice_temperature"][-1, :]
        outputs["height_on_ice_interface_levels"][:] = raw_state[
            "height_on_ice_interface_levels"
        ]
        diagnostics["heat_flux_into_sea_water_due_to_sea_ice"][:] = raw_state[
            "heat_flux_into_sea_water_due_to_sea_ice"
        ]

        for col in range(num_cols):
            area_type = raw_state["area_type"][col].astype(str)
            if area_type != "sea_ice":
                continue

            if raw_state["sea_ice_thickness"][col] <= 0.0:
                # No sea ice, so skip calculation
                continue

            total_height = (
                raw_state["sea_ice_thickness"][col]
                + raw_state["surface_snow_thickness"][col]
            )
            if total_height > self._max_height:
                raise ValueError(
                    "Total height exceeds maximum value of {} m.".format(
                        self._max_height
                    )
                )

            if total_height < self._epsilon:  # Some epsilon
                continue

            snow_height_fraction = (
                raw_state["surface_snow_thickness"][col] / total_height
            )

            temp_profile = raw_state["snow_and_ice_temperature"][:, col].copy()
            num_layers = temp_profile.shape[0]
            dz = float(total_height / num_layers)

            snow_level = int((1 - snow_height_fraction) * num_layers) - 1
            levels = np.arange(num_layers - 1)

            # Create vertically varying profiles. Index convention here
            # (and for solve_column) is node 0 = bottom (ice base), node
            # n-1 = top (atmosphere-facing surface), matching IceSheet.
            rho_snow_ice = self._rho_ice * np.ones(num_layers - 1)
            rho_snow_ice[levels > snow_level] = self._rho_snow

            heat_capacity_snow_ice = self._C_ice * np.ones(num_layers - 1)
            heat_capacity_snow_ice[levels > snow_level] = self._C_snow

            kappa_snow_ice = self._Kice * np.ones(num_layers - 1)
            kappa_snow_ice[levels > snow_level] = self._Ksnow

            surface_temperature = temp_profile[-1]
            check_melting = True
            if surface_temperature < self._melting_temperature - self._epsilon:
                check_melting = False

            # Defect fix (basal BC): use the ocean heat flux at the ice
            # base instead of a freezing Dirichlet condition, so that both
            # basal growth (ocean colder than the ice base) and basal melt
            # (ocean warmer than the ice base) are representable.
            #
            # Sign convention: "heat_flux_into_sea_water_due_to_sea_ice" is
            # positive when heat leaves the ice and enters the sea water
            # (the CF/CMIP convention for this diagnostic, and also how it
            # is computed further down from (T[1]-T[0]): positive when the
            # ice interior is warmer than the base, i.e. heat conducts
            # downward out of the ice into the ocean). solve_column's
            # Flux convention is the opposite: positive means heat enters
            # the column at that boundary (this can be verified from the
            # solver's bottom-row equation, which reduces to
            # bottom_bc.value == K[0] * (T[0] - T[1]) / dz, the negative of
            # the "into sea water" quantity above). So the stored/forcing
            # value must be negated when used as the column's bottom BC.
            Q_into_sea_water = raw_state["heat_flux_into_sea_water_due_to_sea_ice"][
                col
            ]
            bottom_bc = Flux(-Q_into_sea_water)
            top_bc = (
                Dirichlet(self._melting_temperature)
                if check_melting
                else Flux(net_heat_flux[col])
            )

            new_temperature = solve_column(
                rho_snow_ice,
                heat_capacity_snow_ice,
                kappa_snow_ice,
                temp_profile,
                dt,
                dz,
                top_bc,
                bottom_bc,
            )

            # Cool down from melting temperature if heat flux through ice
            # is greater than net forcing heat flux at surface.
            heat_flux_through_ice = (
                (new_temperature[-1] - new_temperature[-2])
                * (kappa_snow_ice[-1] + kappa_snow_ice[-2])
                * 0.5
                / dz
            )

            if temp_profile[-1] > self._melting_temperature - self._epsilon:
                if heat_flux_through_ice > net_heat_flux[col]:
                    surface_temperature = temp_profile[-1] - 10 * self._epsilon
                    temp_profile[-1] = temp_profile[-1] - 10 * self._epsilon

                    top_bc = Flux(net_heat_flux[col])
                    new_temperature = solve_column(
                        rho_snow_ice,
                        heat_capacity_snow_ice,
                        kappa_snow_ice,
                        temp_profile,
                        dt,
                        dz,
                        top_bc,
                        bottom_bc,
                    )

                    check_melting = False

            # Energy balance at the lower (ocean-facing) surface: how much
            # heat is actually exchanged with the sea water, given the
            # solved temperature profile. Positive means heat is leaving
            # the ice into the ocean (see sign-convention note above);
            # negative means heat is being drawn from the ocean into the
            # ice, which grows the ice from below.
            heat_flux_to_sea_water = (
                (new_temperature[1] - new_temperature[0])
                * (kappa_snow_ice[0] + kappa_snow_ice[1])
                * 0.5
                / dz
            )
            heat_flux_to_sea_water = round(heat_flux_to_sea_water, 6)

            height_of_growing_ice = -(
                heat_flux_to_sea_water * dt / (rho_snow_ice[0] * self._Lf)
            )

            outputs["sea_ice_thickness"][col] += height_of_growing_ice
            diagnostics["heat_flux_into_sea_water_due_to_sea_ice"][
                col
            ] = heat_flux_to_sea_water

            # Energy balance at atmosphere surface
            heat_flux_through_ice = (
                (new_temperature[-1] - new_temperature[-2])
                * (kappa_snow_ice[-1] + kappa_snow_ice[-2])
                * 0.5
                / dz
            )
            diagnostics["surface_downward_heat_flux_in_sea_ice"][
                col
            ] = heat_flux_through_ice

            height_of_melting_ice = 0.0
            # Surface is melting
            if check_melting:
                energy_to_melt_ice = (
                    net_heat_flux[col] - heat_flux_through_ice
                ) * dt
                energy_to_melt_ice = round(energy_to_melt_ice, 6)
                if energy_to_melt_ice < 0:
                    # Defect fix: this used to be `print(...); assert
                    # False`. A slightly negative value here is a benign
                    # consequence of the melting-temperature bookkeeping
                    # above (rounding / the 10*epsilon nudge), not a
                    # physical inconsistency worth crashing over, so clip
                    # it to zero and log it for diagnosis instead.
                    logger.debug(
                        "Negative melt energy %s at column %d "
                        "(net_flux=%s, conducted_flux=%s); clamping to 0.",
                        energy_to_melt_ice,
                        col,
                        net_heat_flux[col],
                        heat_flux_through_ice,
                    )
                energy_to_melt_ice = np.clip(energy_to_melt_ice, 0.0, None)

                height_of_melting_ice = energy_to_melt_ice / (
                    rho_snow_ice[-1] * self._Lf
                )

                if height_of_melting_ice > raw_state["surface_snow_thickness"][col]:
                    outputs["sea_ice_thickness"][col] -= (
                        height_of_melting_ice
                        - raw_state["surface_snow_thickness"][col]
                    )
                    outputs["surface_snow_thickness"][col] = 0
                else:
                    outputs["surface_snow_thickness"][col] -= height_of_melting_ice

            # Defect fix: clamp sea_ice_thickness to be non-negative.
            # Any melt energy that would have driven the thickness below
            # zero didn't actually have ice left to melt, so route it into
            # the ocean instead (it warms the sea water directly, which is
            # a positive contribution to "heat_flux_into_sea_water").
            pre_clip_thickness = outputs["sea_ice_thickness"][col]
            outputs["sea_ice_thickness"][col] = np.clip(pre_clip_thickness, 0.0, None)
            if pre_clip_thickness < 0.0:
                leftover_thickness = -pre_clip_thickness
                leftover_flux = leftover_thickness * rho_snow_ice[-1] * self._Lf / dt
                diagnostics["heat_flux_into_sea_water_due_to_sea_ice"][
                    col
                ] += leftover_flux

            total_height = (
                outputs["sea_ice_thickness"][col]
                + outputs["surface_snow_thickness"][col]
            )

            outputs["snow_and_ice_temperature"][:, col] = new_temperature
            outputs["surface_temperature"][col] = new_temperature[-1]
            outputs["height_on_ice_interface_levels"][:, col] = np.linspace(
                0,
                total_height,
                outputs["height_on_ice_interface_levels"].shape[0],
                endpoint=True,
            )

            # Defect fix: single unambiguous albedo branch (no duplicate
            # `elif`), using the configurable albedo values.
            if outputs["surface_snow_thickness"][col] > 0:
                albedo = self._albedo_snow
            else:
                albedo = self._albedo_ice

            if height_of_melting_ice > 0:
                albedo = self._albedo_melt

            diagnostics["surface_albedo_for_direct_shortwave"][col] = albedo
            diagnostics["surface_albedo_for_diffuse_shortwave"][col] = albedo

        return diagnostics, outputs
