import warnings

from sympl import Stepper, initialize_numpy_arrays_with_properties

from .land_ice.component import LandIce
from .sea_ice.component import SeaIce


class IceSheet(Stepper):
    """
    Deprecated monolithic 1-d snow-ice energy balance model.

    .. deprecated::
        ``IceSheet`` used to handle ``sea_ice``, ``land_ice`` and ``land``
        columns in a single component with a hand-rolled tridiagonal solve.
        It has been split into :class:`~climt._components.sea_ice.SeaIce`
        (owns ``area_type == "sea_ice"``) and
        :class:`~climt._components.land_ice.LandIce` (owns
        ``area_type in ("land", "land_ice")``), each built on the shared
        implicit column solver in :mod:`climt._core.snow_ice_column`, with
        defect fixes relative to this monolith documented in their
        docstrings (basal boundary condition, thickness clamping, albedo
        configuration, melt-energy handling).

        ``IceSheet`` is now a thin dispatching shim that constructs a
        ``SeaIce`` and a ``LandIce`` instance and merges their per-column
        results, kept only so that existing callers keep working. New code
        should use ``SeaIce``/``LandIce`` directly.
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
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
        "snow_and_ice_temperature": {
            "dims": ["ice_interface_levels", "*"],
            "units": "degK",
        },
        "sea_surface_temperature": {
            "dims": ["*"],
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
        # Added relative to the pre-split monolith: SeaIce's basal boundary
        # condition is a prescribed ocean heat flux (a genuine defect fix,
        # see SeaIce's docstring) rather than the old hardcoded freezing
        # Dirichlet condition, so it needs this as a real input rather than
        # just a diagnostic. It already has a registered default value of
        # 0.0 W/m^2 in climt's default-state table, so get_default_state
        # is unaffected.
        "heat_flux_into_sea_water_due_to_sea_ice": {
            "dims": ["*"],
            "units": "W m^-2",
        },
    }

    output_properties = {
        "land_ice_thickness": {
            "dims": ["*"],
            "units": "m",
        },
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
        "sea_surface_temperature": {
            "dims": ["*"],
            "units": "degK",
        },
    }

    diagnostic_properties = {
        "upward_heat_flux_at_ground_level_in_soil": {
            "dims": ["*"],
            "units": "W m^-2",
        },
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

    def __init__(self, maximum_snow_ice_height=10, **kwargs):
        """
        Args:
            maximum_snow_ice_height (float):
                The maximum combined height of snow and ice handled by the
                model in :math:`m`. Forwarded to both the internal
                ``SeaIce`` and ``LandIce`` instances.
            levels (int):
                The number of levels on which temperature must be output.

        .. deprecated::
            Use :class:`~climt._components.sea_ice.SeaIce` and
            :class:`~climt._components.land_ice.LandIce` instead.
        """
        warnings.warn(
            "IceSheet is deprecated; use SeaIce and LandIce.",
            DeprecationWarning,
        )
        self._sea = SeaIce(maximum_snow_ice_height=maximum_snow_ice_height, **kwargs)
        self._land = LandIce(maximum_snow_ice_height=maximum_snow_ice_height, **kwargs)
        super(IceSheet, self).__init__(**kwargs)

    def array_call(self, raw_state, timestep):
        # Dispatch/merge design: run BOTH sub-components on the full raw
        # state (not a pre-sliced subset) and merge their results. Each
        # sub-component's per-column loop `continue`s past cells outside
        # its own ownership (`area_type == "sea_ice"` for SeaIce,
        # `area_type in ("land", "land_ice")` for LandIce), and both seed
        # their output/diagnostic arrays from `raw_state` before that loop
        # runs -- so for most shared fields, a non-owned cell's value in
        # each component's result is simply the input passed straight
        # through (confirmed by reading both `array_call`s' "Copy input
        # values" blocks).
        #
        # The one field that is NOT a simple input pass-through on
        # non-owned cells is `surface_temperature`: neither SeaIce nor
        # LandIce takes it as an input at all -- both instead *derive* it,
        # for every column, from `snow_and_ice_temperature[-1, col]` (the
        # top of the shared temperature profile), only overwriting it with
        # the newly solved value on their own owned cells. Since both
        # components compute this the same way from the same shared input,
        # they agree with each other on cells neither owns, but a plain
        # merge-by-union would silently prefer whichever component's dict
        # is copied last rather than reflecting real per-cell ownership.
        # So instead of assuming disjoint writes, every shared
        # output/diagnostic field below is combined with an explicit
        # `area_type` mask: default to the LandIce result (correct for
        # `land`/`land_ice` columns, and equal to an input pass-through for
        # any other area type, e.g. plain "sea"/"ocean" with no ice), then
        # overwrite `sea_ice`-owned columns with the SeaIce result.
        sea_diagnostics, sea_outputs = self._sea.array_call(raw_state, timestep)
        land_diagnostics, land_outputs = self._land.array_call(raw_state, timestep)

        area_type = raw_state["area_type"].astype(str)
        sea_mask = area_type == "sea_ice"

        outputs = initialize_numpy_arrays_with_properties(
            self.output_properties, raw_state, self.input_properties
        )
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )

        # Shared 1-D output fields: mask by area_type.
        for key in ("surface_snow_thickness", "surface_temperature"):
            outputs[key][:] = land_outputs[key]
            outputs[key][sea_mask] = sea_outputs[key][sea_mask]

        # Shared 2-D (level, column) output fields: mask on the column axis.
        for key in ("snow_and_ice_temperature", "height_on_ice_interface_levels"):
            outputs[key][:] = land_outputs[key]
            outputs[key][:, sea_mask] = sea_outputs[key][:, sea_mask]

        # sea_ice_thickness / land_ice_thickness are each declared by
        # exactly one sub-component (the other doesn't have the key at
        # all), so there is a single owner and nothing to mask.
        outputs["sea_ice_thickness"][:] = sea_outputs["sea_ice_thickness"]
        outputs["land_ice_thickness"][:] = land_outputs["land_ice_thickness"]

        # sea_surface_temperature was a pure pass-through in the original
        # monolith (never written by either area-type branch); neither
        # SeaIce nor LandIce declares it as an output, so carry it through
        # from the input unchanged, as IceSheet always did.
        outputs["sea_surface_temperature"][:] = raw_state["sea_surface_temperature"]

        # Diagnostics with a single owner: heat_flux_into_sea_water_due_to_
        # sea_ice and surface_downward_heat_flux_in_sea_ice are SeaIce-only
        # in the split components (LandIce doesn't declare them);
        # upward_heat_flux_at_ground_level_in_soil is LandIce-only.
        diagnostics["heat_flux_into_sea_water_due_to_sea_ice"][:] = sea_diagnostics[
            "heat_flux_into_sea_water_due_to_sea_ice"
        ]
        diagnostics["surface_downward_heat_flux_in_sea_ice"][:] = sea_diagnostics[
            "surface_downward_heat_flux_in_sea_ice"
        ]
        diagnostics["upward_heat_flux_at_ground_level_in_soil"][:] = land_diagnostics[
            "upward_heat_flux_at_ground_level_in_soil"
        ]

        # Albedo diagnostics are computed by both sub-components (each
        # over its own owned cells); mask the same way as the shared
        # outputs above.
        for key in (
            "surface_albedo_for_direct_shortwave",
            "surface_albedo_for_diffuse_shortwave",
        ):
            diagnostics[key][:] = land_diagnostics[key]
            diagnostics[key][sea_mask] = sea_diagnostics[key][sea_mask]

        return diagnostics, outputs
