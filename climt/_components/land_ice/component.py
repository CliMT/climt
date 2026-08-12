import logging

import numpy as np
from sympl import Stepper, get_constant, initialize_numpy_arrays_with_properties

from ..._core.backend import jit_compile, prange
from ..._core.snow_ice_column import (
    _BC_DIRICHLET,
    _BC_FLUX,
    _solve_column_kernel,
)

logger = logging.getLogger(__name__)


@jit_compile
def _round6(x):
    # Match Python's round(x, 6) (round-half-to-even) inside njit.
    return np.round(x * 1e6) / 1e6


@jit_compile(parallel=True)
def _land_ice_columns(
    active, temperature_in, land_ice_thickness_in, surface_snow_thickness_in,
    soil_surface_temperature, net_flux, dt,
    rho_ice, rho_snow, C_ice, C_snow, Kice, Ksnow, Lf, T_melt, eps,
    albedo_snow, albedo_ice, albedo_melt,
    land_ice_thickness_out, surface_snow_thickness_out, temperature_out,
    surface_temperature_out, height_out, ground_flux_out, albedo_out,
    neg_energy_flag, neg_energy_value,
):
    """Per-column land-ice / snow-on-land thermodynamics, parallel over columns.

    Pass-through columns must be pre-filled by the caller; this kernel only
    writes columns flagged ``active[col] != 0``. Columns are independent.
    """
    num_cols = active.shape[0]
    num_layers = temperature_in.shape[0]
    n_iface = height_out.shape[0]

    for col in prange(num_cols):
        if active[col] == 0:
            continue

        total_height = (
            land_ice_thickness_in[col] + surface_snow_thickness_in[col]
        )
        snow_height_fraction = surface_snow_thickness_in[col] / total_height

        temp_profile = temperature_in[:, col].copy()
        dz = total_height / num_layers
        snow_level = int((1.0 - snow_height_fraction) * num_layers) - 1

        rho_snow_ice = np.empty(num_layers - 1)
        heat_capacity_snow_ice = np.empty(num_layers - 1)
        kappa_snow_ice = np.empty(num_layers - 1)
        for i in range(num_layers - 1):
            if i > snow_level:
                rho_snow_ice[i] = rho_snow
                heat_capacity_snow_ice[i] = C_snow
                kappa_snow_ice[i] = Ksnow
            else:
                rho_snow_ice[i] = rho_ice
                heat_capacity_snow_ice[i] = C_ice
                kappa_snow_ice[i] = Kice

        check_melting = temp_profile[-1] >= T_melt - eps

        bot_flag = _BC_DIRICHLET
        bot_val = soil_surface_temperature[col]
        if check_melting:
            top_flag = _BC_DIRICHLET
            top_val = T_melt
        else:
            top_flag = _BC_FLUX
            top_val = net_flux[col]

        new_temperature = _solve_column_kernel(
            rho_snow_ice, heat_capacity_snow_ice, kappa_snow_ice,
            temp_profile, dt, dz, top_flag, top_val, bot_flag, bot_val,
        )

        heat_flux_through_ice = (
            (new_temperature[-1] - new_temperature[-2])
            * (kappa_snow_ice[-1] + kappa_snow_ice[-2]) * 0.5 / dz
        )
        if temp_profile[-1] > T_melt - eps:
            if heat_flux_through_ice > net_flux[col]:
                temp_profile[-1] = temp_profile[-1] - 10.0 * eps
                top_flag = _BC_FLUX
                top_val = net_flux[col]
                new_temperature = _solve_column_kernel(
                    rho_snow_ice, heat_capacity_snow_ice, kappa_snow_ice,
                    temp_profile, dt, dz, top_flag, top_val, bot_flag, bot_val,
                )
                check_melting = False

        ground_flux_out[col] = (
            (new_temperature[0] - new_temperature[1]) * kappa_snow_ice[0] / dz
        )

        heat_flux_through_ice = (
            (new_temperature[-1] - new_temperature[-2])
            * (kappa_snow_ice[-1] + kappa_snow_ice[-2]) * 0.5 / dz
        )

        height_of_melting_ice = 0.0
        if check_melting:
            energy_to_melt_ice = _round6(
                (net_flux[col] - heat_flux_through_ice) * dt
            )
            if energy_to_melt_ice < 0:
                neg_energy_flag[col] = 1
                neg_energy_value[col] = energy_to_melt_ice
            if energy_to_melt_ice < 0.0:
                energy_to_melt_ice = 0.0

            height_of_melting_ice = energy_to_melt_ice / (
                rho_snow_ice[-1] * Lf
            )
            if height_of_melting_ice > surface_snow_thickness_in[col]:
                land_ice_thickness_out[col] -= (
                    height_of_melting_ice - surface_snow_thickness_in[col]
                )
                surface_snow_thickness_out[col] = 0.0
            else:
                surface_snow_thickness_out[col] = (
                    surface_snow_thickness_in[col] - height_of_melting_ice
                )

        if land_ice_thickness_out[col] < 0.0:
            land_ice_thickness_out[col] = 0.0
        if surface_snow_thickness_out[col] < 0.0:
            surface_snow_thickness_out[col] = 0.0

        total_height = (
            land_ice_thickness_out[col] + surface_snow_thickness_out[col]
        )
        for j in range(num_layers):
            temperature_out[j, col] = new_temperature[j]
        surface_temperature_out[col] = new_temperature[-1]
        for j in range(n_iface):
            height_out[j, col] = total_height * j / (n_iface - 1)

        if surface_snow_thickness_out[col] > 0:
            albedo = albedo_snow
        else:
            albedo = albedo_ice
        if height_of_melting_ice > 0:
            albedo = albedo_melt
        albedo_out[col] = albedo


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

        # Active columns: land / land-ice cells with a non-negligible
        # snow+ice column. area_type is a string array, so the mask is built
        # here (in Python), not inside the jit kernel.
        area_type = raw_state["area_type"].astype(str)
        land_ice_thickness = np.asarray(
            raw_state["land_ice_thickness"], dtype=float
        )
        snow = np.asarray(raw_state["surface_snow_thickness"], dtype=float)
        total_height_in = land_ice_thickness + snow
        is_land = (area_type == "land_ice") | (area_type == "land")
        active = (is_land & (total_height_in >= self._epsilon)).astype(np.int8)

        if np.any(is_land & (total_height_in > self._max_height)):
            raise ValueError(
                "Total height exceeds maximum value of {} m.".format(
                    self._max_height
                )
            )

        neg_energy_flag = np.zeros(num_cols, dtype=np.int8)
        neg_energy_value = np.zeros(num_cols)

        _land_ice_columns(
            active,
            np.ascontiguousarray(
                raw_state["snow_and_ice_temperature"], dtype=float
            ),
            land_ice_thickness,
            snow,
            np.ascontiguousarray(
                raw_state["soil_surface_temperature"], dtype=float
            ),
            np.ascontiguousarray(net_heat_flux, dtype=float),
            dt,
            self._rho_ice, self._rho_snow, self._C_ice, self._C_snow,
            self._Kice, self._Ksnow, self._Lf, self._melting_temperature,
            self._epsilon, self._albedo_snow, self._albedo_ice,
            self._albedo_melt,
            outputs["land_ice_thickness"],
            outputs["surface_snow_thickness"],
            outputs["snow_and_ice_temperature"],
            outputs["surface_temperature"],
            outputs["height_on_ice_interface_levels"],
            diagnostics["upward_heat_flux_at_ground_level_in_soil"],
            diagnostics["surface_albedo_for_direct_shortwave"],
            neg_energy_flag,
            neg_energy_value,
        )
        diagnostics["surface_albedo_for_diffuse_shortwave"][:] = diagnostics[
            "surface_albedo_for_direct_shortwave"
        ]

        for col in np.nonzero(neg_energy_flag)[0]:
            logger.debug(
                "Negative melt energy %s at column %d (net_flux=%s, "
                "conducted_flux=%s); clamping to 0.",
                neg_energy_value[col], int(col), net_heat_flux[col],
                diagnostics["upward_heat_flux_at_ground_level_in_soil"][col],
            )

        return diagnostics, outputs
