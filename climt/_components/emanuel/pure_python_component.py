# -*- coding: utf-8 -*-
import numpy as np
from sympl import (
    ImplicitTendencyComponent,
    get_constant,
    initialize_numpy_arrays_with_properties,
)

from ..._core import bolton_q_sat, ensure_contiguous_state
from .pure_python_v2 import EmanuelConvectionPythonV2


class EmanuelConvectionPythonComponent(ImplicitTendencyComponent):
    """
    The Emanuel convection scheme from `[Emanuel and Zivkovic-Rothman]`_
    implemented in pure Python with support for various array backends.

    .. _[Emanuel and Zivkovic-Rothman]:
            http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1999)056%3C1766%3ADAEOAC%3E2.0.CO%3B2
    """

    input_properties = {
        "air_temperature": {
            "dims": ["*", "mid_levels"],
            "units": "degK",
        },
        "specific_humidity": {
            "dims": ["*", "mid_levels"],
            "units": "kg/kg",
        },
        "eastward_wind": {
            "dims": ["*", "mid_levels"],
            "units": "m s^-1",
        },
        "northward_wind": {
            "dims": ["*", "mid_levels"],
            "units": "m s^-1",
        },
        "air_pressure": {
            "dims": ["*", "mid_levels"],
            "units": "mbar",
        },
        "air_pressure_on_interface_levels": {
            "dims": ["*", "interface_levels"],
            "units": "mbar",
        },
        "cloud_base_mass_flux": {
            "dims": ["*"],
            "units": "kg m^-2 s^-1",
        },
    }

    diagnostic_properties = {
        "convective_state": {
            "dims": ["*"],
            "units": "dimensionless",
            "dtype": np.int32,
        },
        "convective_precipitation_rate": {
            "dims": ["*"],
            "units": "mm day^-1",
        },
        "convective_downdraft_velocity_scale": {
            "dims": ["*"],
            "units": "m s^-1",
        },
        "convective_downdraft_temperature_scale": {
            "dims": ["*"],
            "units": "degK",
        },
        "convective_downdraft_specific_humidity_scale": {
            "dims": ["*"],
            "units": "kg/kg",
        },
        "cloud_base_mass_flux": {
            "dims": ["*"],
            "units": "kg m^-2 s^-1",
        },
        "atmosphere_convective_available_potential_energy": {
            "dims": ["*"],
            "units": "J kg^-1",
        },
        "air_temperature_tendency_from_convection": {
            "dims": ["*", "mid_levels"],
            "units": "degK day^-1",
        },
    }

    tendency_properties = {
        "air_temperature": {"units": "degK s^-1"},
        "specific_humidity": {"units": "kg/kg s^-1"},
        "eastward_wind": {"units": "m s^-2"},
        "northward_wind": {"units": "m s^-2"},
    }

    def __init__(
        self,
        minimum_convecting_layer=1,
        autoconversion_water_content_threshold=0.0011,
        autoconversion_temperature_threshold=-55,
        entrainment_mixing_coefficient=1.5,
        downdraft_area_fraction=0.05,
        precipitation_fraction_outside_cloud=0.12,
        speed_water_droplets=50.0,
        speed_snow=5.5,
        rain_evaporation_coefficient=1.0,
        snow_evaporation_coefficient=0.8,
        convective_momentum_transfer_coefficient=0.7,
        downdraft_surface_velocity_coefficient=10.0,
        convection_bouyancy_threshold=0.9,
        mass_flux_relaxation_rate=0.1,
        mass_flux_damping_rate=0.1,
        reference_mass_flux_timescale=300.0,
        **kwargs,
    ):
        self._con_mom_txfr = convective_momentum_transfer_coefficient
        self._downdraft_area_frac = downdraft_area_fraction
        self._precip_frac_outside_cloud = precipitation_fraction_outside_cloud
        self._min_conv_layer = minimum_convecting_layer
        self._crit_humidity = autoconversion_water_content_threshold
        self._crit_temp = autoconversion_temperature_threshold
        self._entrain_coeff = entrainment_mixing_coefficient
        self._droplet_speed = speed_water_droplets
        self._snow_speed = speed_snow
        self._rain_evap = rain_evaporation_coefficient
        self._snow_evap = snow_evaporation_coefficient
        self._beta = downdraft_surface_velocity_coefficient
        self._dtmax = convection_bouyancy_threshold
        self._mf_damp = mass_flux_damping_rate
        self._alpha = mass_flux_relaxation_rate
        self._mf_timescale = reference_mass_flux_timescale
        self._ntracers = 0

        # Get constants from sympl
        self._g = get_constant("gravitational_acceleration", "m/s^2")
        self._Cpd = get_constant(
            "heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"
        )
        self._Cpv = get_constant("heat_capacity_of_vapor_phase", "J/kg/degK")
        self._Rdair = get_constant("gas_constant_of_dry_air", "J/kg/degK")
        self._Rcond = get_constant("gas_constant_of_vapor_phase", "J/kg/degK")
        self._Lv = get_constant("latent_heat_of_condensation", "J/kg")
        self._rho_condensible = get_constant("density_of_liquid_phase", "kg/m^3")
        self._Cl = get_constant("specific_enthalpy_of_vapor_phase", "J/kg")

        # Create the internal logic object
        self._python_scheme = EmanuelConvectionPythonV2(
            MINORIG=self._min_conv_layer,
            ELCRIT=self._crit_humidity,
            TLCRIT=self._crit_temp,
            ENTP=self._entrain_coeff,
            SIGD=self._downdraft_area_frac,
            SIGS=self._precip_frac_outside_cloud,
            OMTRAIN=self._droplet_speed,
            OMTSNOW=self._snow_speed,
            COEFFR=self._rain_evap,
            COEFFS=self._snow_evap,
            CU=self._con_mom_txfr,
            DTMAX=self._dtmax,
            BETA=self._beta,
            ALPHA=self._alpha,
            DAMP=self._mf_damp,
            CPD=self._Cpd,
            CPV=self._Cpv,
            CL=self._Cl,
            RV=self._Rcond,
            RD=self._Rdair,
            LV0=self._Lv,
            G=self._g,
            ROWL=self._rho_condensible,
            DELT0=self._mf_timescale,
            **kwargs,
        )

        super(EmanuelConvectionPythonComponent, self).__init__(**kwargs)

    @ensure_contiguous_state
    def array_call(self, raw_state, timestep):
        num_cols, num_levs = raw_state["air_temperature"].shape

        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, raw_state, self.input_properties
        )
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )

        # EmanuelConvectionPython expects (nlev, ncol) while raw_state is (ncol, nlev)
        python_state = {
            "air_temperature": raw_state["air_temperature"].T,
            "specific_humidity": raw_state["specific_humidity"].T,
            "eastward_wind": raw_state["eastward_wind"].T,
            "northward_wind": raw_state["northward_wind"].T,
            "air_pressure": raw_state["air_pressure"].T,
            "air_pressure_on_interface_levels": raw_state[
                "air_pressure_on_interface_levels"
            ].T,
            "cloud_base_mass_flux": raw_state["cloud_base_mass_flux"],
        }

        python_tendencies, python_diagnostics = self._python_scheme.array_call(
            python_state, timestep
        )

        # Map back to sympl tendencies and diagnostics (transpose back)
        for key in tendencies:
            tendencies[key][:] = python_tendencies[key].T

        for key in diagnostics:
            if key != "air_temperature_tendency_from_convection":
                if diagnostics[key].ndim == 2:
                    diagnostics[key][:] = python_diagnostics[key].T
                else:
                    diagnostics[key][:] = python_diagnostics[key]

        diagnostics["air_temperature_tendency_from_convection"][:] = (
            tendencies["air_temperature"] * 86400.0
        )

        return tendencies, diagnostics
