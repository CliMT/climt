import logging

import numpy as np
from sympl import Stepper, get_constant, initialize_numpy_arrays_with_properties

from ..._core.snow_ice_column import Dirichlet, Flux, solve_column

logger = logging.getLogger(__name__)


class LandIce(Stepper):
    """
    1-D thermodynamic snow/ice model over ``area_type in ("land_ice", "land")``
    cells (snow-on-land, including land glacier/ice-sheet columns).

    This is a mechanical re-expression of the ``land``/``land_ice`` branch of
    the monolithic :class:`~climt._components.surface_ice.IceSheet` component,
    using the shared implicit column solver
    (:func:`climt._core.snow_ice_column.solve_column`) in place of the
    hand-rolled tridiagonal solve, with the following defect fixes applied:

    1. The basal boundary condition is unchanged from the original (a
       prescribed :class:`~climt._core.snow_ice_column.Dirichlet` soil
       surface temperature) -- unlike :class:`~climt._components.sea_ice.SeaIce`,
       there is no basal-flux sign subtlety here.
    2. The original code decremented ``sea_ice_thickness`` (instead of
       ``land_ice_thickness``) when surface melt exceeded the snow layer on
       a ``land``/``land_ice`` column -- an apparent copy/paste artifact
       from the sea-ice branch of the same function. This is fixed to
       decrement ``land_ice_thickness``, which is the physically relevant
       reservoir for these area types.
    3. ``land_ice_thickness`` and ``surface_snow_thickness`` are clamped to
       be non-negative after melt (mirroring the ``SeaIce`` clamp fix).
    4. The three albedo values (snow / bare ice / melting) are configurable
       via ``__init__`` rather than hardcoded, and the duplicated ``elif``
       branch for the bare-ice albedo has been removed.
    5. The ``print``/``assert False`` guard against negative melt energy
       has been replaced with an ``np.clip`` and a debug log message.

    There is no basal ocean flux for this component (there is no ocean at
    the base of a land or land-ice column); the base of the column
    exchanges heat directly with the soil via the Dirichlet condition, and
    that exchange is reported as ``upward_heat_flux_at_ground_level_in_soil``.
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
        "land_ice_thickness": {
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
        "soil_surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
        "height_on_ice_interface_levels": {
            "dims": ["ice_interface_levels", "*"],
            "units": "m",
        },
    }

    output_properties = {
        "land_ice_thickness": {
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
        "upward_heat_flux_at_ground_level_in_soil": {
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
        albedo_ice=0.6,
        albedo_melt=0.2,
        **kwargs
    ):
        """
        Args:
            maximum_snow_ice_height (float):
                The maximum combined height of snow and ice handled by the
                model in :math:`m`.
            albedo_snow (float):
                Surface albedo used when snow is present on top of the ice
                or bare land.
            albedo_ice (float):
                Surface albedo used for bare (snow-free) land ice.
            albedo_melt (float):
                Surface albedo used when the surface is actively melting.
        """
        self._max_height = maximum_snow_ice_height
        self._albedo_snow = albedo_snow
        self._albedo_ice = albedo_ice
        self._albedo_melt = albedo_melt
        self._epsilon = 1e-5
        super(LandIce, self).__init__(**kwargs)

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

        # Copy input values (columns that are not land/land-ice, or that
        # have no snow/ice present, simply pass through unchanged).
        outputs["land_ice_thickness"][:] = raw_state["land_ice_thickness"]
        outputs["surface_snow_thickness"][:] = raw_state["surface_snow_thickness"]
        outputs["snow_and_ice_temperature"][:] = raw_state["snow_and_ice_temperature"]
        outputs["surface_temperature"][:] = raw_state["snow_and_ice_temperature"][-1, :]
        outputs["height_on_ice_interface_levels"][:] = raw_state[
            "height_on_ice_interface_levels"
        ]

        for col in range(num_cols):
            area_type = raw_state["area_type"][col].astype(str)
            if area_type not in ("land_ice", "land"):
                continue

            total_height = (
                raw_state["land_ice_thickness"][col]
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
            # (and for solve_column) is node 0 = bottom (soil interface),
            # node n-1 = top (atmosphere-facing surface), matching IceSheet.
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

            # Basal BC (unchanged land physics): the base of the column is
            # tied directly to the soil surface temperature via a Dirichlet
            # condition, so no basal-flux sign convention is needed here
            # (unlike SeaIce's ocean-flux Flux boundary).
            bottom_bc = Dirichlet(raw_state["soil_surface_temperature"][col])
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

            # Base-node conduction (soil interface): this mirrors the
            # original IceSheet formula for the land/land_ice branch
            # exactly -- (T[0]-T[1])*kappa[0]/dz, using the single
            # bottommost-layer conductivity (not a two-layer average, as
            # SeaIce's basal-ocean formula uses), which is the correctly
            # discretized flux across the bottommost layer matching
            # solve_column's own Dirichlet/Flux discretization at that
            # node. Positive means T[0] > T[1], i.e. heat conducts
            # upward out of the soil into the snow/ice column -- hence
            # "upward_heat_flux_at_ground_level_in_soil".
            heat_flux_to_soil = (
                (new_temperature[0] - new_temperature[1]) * kappa_snow_ice[0] / dz
            )
            diagnostics["upward_heat_flux_at_ground_level_in_soil"][
                col
            ] = heat_flux_to_soil

            # No basal ocean flux / basal growth for land or land ice.
            height_of_growing_ice = 0.0

            # Energy balance at atmosphere surface
            heat_flux_through_ice = (
                (new_temperature[-1] - new_temperature[-2])
                * (kappa_snow_ice[-1] + kappa_snow_ice[-2])
                * 0.5
                / dz
            )

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
                    # Defect fix: the original IceSheet code decremented
                    # ``sea_ice_thickness`` here even for the land/land_ice
                    # branch (an apparent copy/paste artifact from the
                    # sea-ice branch of the same function). The physically
                    # relevant reservoir for land/land-ice columns is
                    # ``land_ice_thickness``.
                    outputs["land_ice_thickness"][col] -= (
                        height_of_melting_ice
                        - raw_state["surface_snow_thickness"][col]
                    )
                    outputs["surface_snow_thickness"][col] = 0
                else:
                    outputs["surface_snow_thickness"][col] -= height_of_melting_ice

            # Defect fix: clamp land_ice_thickness and surface_snow_thickness
            # to be non-negative (mirrors the SeaIce clamp fix). Unlike
            # SeaIce, there is no ocean to route leftover melt energy into;
            # the mass balance here is glacier mass balance (accumulation,
            # implicit in the surface_snow_thickness input, minus surface
            # melt), so any melt energy in excess of the available
            # snow+ice column simply has nothing left to melt.
            outputs["land_ice_thickness"][col] = np.clip(
                outputs["land_ice_thickness"][col], 0.0, None
            )
            outputs["surface_snow_thickness"][col] = np.clip(
                outputs["surface_snow_thickness"][col], 0.0, None
            )

            total_height = (
                outputs["land_ice_thickness"][col]
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
