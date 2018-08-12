from sympl import (
    TendencyComponent, get_constant, initialize_numpy_arrays_with_properties
)
from ...._core import (
    mass_to_volume_mixing_ratio, get_interface_values, ensure_contiguous_state
)
import numpy as np
from numpy import pi as numpy_pi
from ..rrtmg_common import (
    rrtmg_cloud_overlap_method_dict, rrtmg_cloud_props_dict,
    rrtmg_cloud_ice_props_dict, rrtmg_cloud_liquid_props_dict,
    rrtmg_aerosol_input_dict)
import logging
try:
    from . import _rrtmg_sw
except ImportError:
    logging.warning(
        'Import failed. RRTMG Shortwave is likely not compiled and will not '
        'be available.'
    )


class RRTMGShortwave(TendencyComponent):
    """
    The Rapid Radiative Transfer Model (RRTMG).

    This module wraps RRTMG for shortwave radiation
    (i.e, emission from the sun).

    """

    num_shortwave_bands = 14
    num_ecmwf_aerosols = 6

    input_properties = {
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'mbar'
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'mbar'
        },
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK'
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless'
        },
        'mole_fraction_of_ozone_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless'
        },
        'mole_fraction_of_carbon_dioxide_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless'
        },
        'mole_fraction_of_methane_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless'
        },
        'mole_fraction_of_nitrous_oxide_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless'
        },
        'mole_fraction_of_oxygen_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless'
        },
        'mass_content_of_cloud_ice_in_atmosphere_layer': {
            'dims': ['mid_levels', '*'],
            'units': 'g m^-2'
        },
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {
            'dims': ['mid_levels', '*'],
            'units': 'g m^-2'
        },
        'cloud_ice_particle_size': {
            'dims': ['mid_levels', '*'],
            'units': 'micrometer'
        },
        'cloud_water_droplet_radius': {
            'dims': ['mid_levels', '*'],
            'units': 'micrometer'
        },
        'cloud_area_fraction_in_atmosphere_layer': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK'
        },
        'zenith_angle': {
            'dims': ['*'],
            'units': 'radians'
        },
        'surface_albedo_for_direct_shortwave': {
            'dims': ['*'],
            'units': 'dimensionless'
        },
        'surface_albedo_for_direct_near_infrared': {
            'dims': ['*'],
            'units': 'dimensionless'
        },
        'surface_albedo_for_diffuse_near_infrared': {
            'dims': ['*'],
            'units': 'dimensionless'},
        'surface_albedo_for_diffuse_shortwave': {
            'dims': ['*'],
            'units': 'dimensionless'
        },
        'shortwave_optical_thickness_due_to_cloud': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless'
        },
        'shortwave_optical_thickness_due_to_aerosol': {
            'dims': ['num_shortwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless'
        },
        'single_scattering_albedo_due_to_cloud': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless'
        },
        'single_scattering_albedo_due_to_aerosol': {
            'dims': ['num_shortwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless'
        },
        'cloud_asymmetry_parameter': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless'
        },
        'aerosol_asymmetry_parameter': {
            'dims': ['num_shortwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless'
        },
        'cloud_forward_scattering_fraction': {
            'dims': ['mid_levels', '*', 'num_shortwave_bands'],
            'units': 'dimensionless'
        },
        'aerosol_optical_depth_at_55_micron': {
            'dims': ['num_ecmwf_aerosols', 'mid_levels', '*'],
            'units': 'dimensionless'
        },
        'solar_cycle_fraction': {
            'dims': [],
            'units': 'dimensionless'
        },
        'flux_adjustment_for_earth_sun_distance': {
            'dims': [],
            'units': 'dimensionless'
        },
    }

    tendency_properties = {
        'air_temperature': {'units': 'degK day^-1'},
    }

    diagnostic_properties = {
        'upwelling_shortwave_flux_in_air': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'downwelling_shortwave_flux_in_air': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'upwelling_shortwave_flux_in_air_assuming_clear_sky': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'downwelling_shortwave_flux_in_air_assuming_clear_sky': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'air_temperature_tendency_from_shortwave_assuming_clear_sky': {
            'dims': ['mid_levels', '*'],
            'units': 'degK day^-1',
        },
        'air_temperature_tendency_from_shortwave': {
            'dims': ['mid_levels', '*'],
            'units': 'degK day^-1',
        },
    }

    def __init__(
            self,
            cloud_overlap_method='random',
            cloud_optical_properties='liquid_and_ice_clouds',
            cloud_ice_properties='ebert_curry_two',
            cloud_liquid_water_properties='radius_dependent_absorption',
            solar_variability_method=0,
            use_solar_constant_from_fortran=False,
            ignore_day_of_year=False,
            facular_sunspot_amplitude=None,
            solar_variability_by_band=None,
            aerosol_type='no_aerosol', **kwargs):
        """
        Args:

            cloud_overlap_method (int):
                Choose the method to do overlap with:

                * 'clear_only' = Clear only (no clouds)
                * 'random' = Random
                * 'maximum_random' = Maximum/Random
                * 'maximum' = Maximum.

            cloud_optical_properties (string):
                Choose how cloud optical properties are calculated:

                * :code:`direct_input` = Both cloud fraction and cloud optical depth are input directly
                * :code:`single_cloud_type` = Cloud fraction and cloud physical properties are input,
                  ice and liquid clouds are treated together, cloud absorptivity is a constant value (0.060241)
                * :code:`liquid_and_ice_clouds` = Cloud fraction and cloud physical properties are input, ice and liquid
                  clouds are treated separately

            cloud_ice_properties (string):
                set bounds on ice particle size.

                * :code:`ebert_curry_one` = ice particle has effective radius >= 10.0 micron `[Ebert and Curry 1992]`_
                * :code:`ebert_curry_two` = ice particle has effective radius between 13.0 and
                  130.0 micron `[Ebert and Curry 1992]`_
                * :code:`key_streamer_manual` = ice particle has effective radius between 5.0 and 131.0 micron
                  `[Key, Streamer Ref. Manual, 1996]`_
                * :code:`fu` = ice particle has generalised effective size (dge) between 5.0 and 140.0 micron
                  `[Fu, 1996]`_. (dge = 1.0315 * r_ec)

                Default value is 0.

            cloud_liquid_water_properties (string):
                set treatment of cloud liquid water.

                * :code:`radius_independent_absorption` = use radius independent absorption coefficient
                * :code:`radius_dependent_absorption` = use radius dependent absorption coefficient (radius between 2.5 and 60 micron)


            solar_variability_method (int):
                set the solar variability model used by RRTMG.

                * solar_variability_method = -1:

                    * If :code:`use_solar_constant_from_fortran = True`: No solar variability and no solar cycle
                      with a solar constant of 1368.22 :math:`W m^{-2}`.
                    * If :code:`use_solar_constant_from_fortran = False`: Solar variability defined by setting
                      non-zero scale factors in :code:`solar_variability_by_band`.

                * solar_variability_method = 0:

                    * If :code:`use_solar_constant_from_fortran = True`: No solar variability and no solar cycle
                      with a solar constant of 1360.85 :math:`W m^{-2}`, with facular and
                      sunspot effects fixed to the mean of solar cycles 13-24.
                    * If :code:`use_solar_constant_from_fortran = False`: No solar variability and no solar cycle.

                * solar_variability_method = 1: Solar variability using the NRLSSI2 solar model
                  with solar cycle contribution determined by :code:`solar_cycle_fraction` in
                  the model state, and facular and sunspot adjustment scale factors specified
                  in :code:`facular_sunspot_amplitude`.

                * solar_variability_method = 2: Solar variability using the NRLSSI2 solar model
                  using solar cycle determined by direct specification of **Mg** (facular)
                  and **SB** (sunspot) indices provided in :code:`facular_sunspot_amplitude`.
                  :code:`solar_constant` is ignored.

                * solar_variability_method = 3:
                 * If :code:`use_internal_solar_constant = True`: No solar variability and no solar cycle
                   with a solar constant of 1360.85 :math:`W m^{-2}`.
                 * If :code:`use_internal_solar_constant = False`: scale factors in :code:`solar_variability_by_band`.

            use_solar_constant_from_fortran (bool):
                If :code:`False`, the solar constant is taken from the constants library. The default
                value is :code:`False`.

            ignore_day_of_year (bool):
                If :code:`True`, the solar output does not vary by day of year (i.e, higher close to the
                solstices and lesser close to the equinoxes). Default value is :code:`False`.


            facular_sunspot_amplitude (array of dimension 2):
                Facular and Sunspot amplitude variability parameters, described previously.

            solar_variability_by_band (array of dimension 14 = number of spectral bands):
                scale factors for solar variability in all spectral bands.

            aerosol_type (string):
                Type of aerosol inputs to RRTMG.

                * :code:`no_aerosol`: No Aerosol.
                * :code:`ecmwf`: ECMWF method. Requires aerosol optical depth at 55 micron as the
                  state quantity :code:`aerosol_optical_depth_at_55_micron`.
                * :code:`all_aerosol_properties`: Input all aerosol optical properties.

        .. _[Ebert and Curry 1992]:
            http://onlinelibrary.wiley.com/doi/10.1029/91JD02472/abstract

        .. _[Key, Streamer Ref. Manual, 1996]:
            https://stratus.ssec.wisc.edu/streamer/userman.pdf

        .. _[Fu, 1996]:
             http://journals.ametsoc.org/doi/abs/10.1175/1520-0442(1996)009%3C2058%3AAAPOTS%3E2.0.CO%3B2
            """

        self._cloud_overlap = rrtmg_cloud_overlap_method_dict[cloud_overlap_method.lower()]
        self._cloud_optics = rrtmg_cloud_props_dict[cloud_optical_properties.lower()]
        self._ice_props = rrtmg_cloud_ice_props_dict[cloud_ice_properties.lower()]
        self._liq_props = rrtmg_cloud_liquid_props_dict[cloud_liquid_water_properties.lower()]
        self._solar_var_flag = solar_variability_method
        self._ignore_day_of_year = ignore_day_of_year

        if facular_sunspot_amplitude is None:
            self._fac_sunspot_coeff = np.ones(2)
        else:
            self._fac_sunspot_coeff = facular_sunspot_amplitude

        if solar_variability_by_band is None:
            self._solar_var_by_band = np.ones(16)
        else:
            self._solar_var_by_band = solar_variability_by_band

        self._aerosol_type = rrtmg_aerosol_input_dict[aerosol_type.lower()]

        if use_solar_constant_from_fortran:
            self._solar_const = 0
        else:
            self._solar_const = get_constant('stellar_irradiance', 'W/m^2')

        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._planck = get_constant('planck_constant', 'erg s')
        self._boltzmann = get_constant('boltzmann_constant', 'erg K^-1')
        self._c = get_constant('speed_of_light', 'cm s^-1')
        self._Na = get_constant('avogadro_constant', 'mole^-1')
        self._loschmidt = get_constant('loschmidt_constant', 'cm^-3')
        self._R = get_constant('universal_gas_constant', 'erg mol^-1 K^-1')
        self._stef_boltz = get_constant('stefan_boltzmann_constant', 'W cm^-2 K^-4')
        self._secs_per_day = get_constant('seconds_per_day', 'dimensionless')
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K')

        _rrtmg_sw.set_constants(
            numpy_pi, self._g,
            self._planck,
            self._boltzmann,
            self._c,
            self._Na,
            self._loschmidt,
            self._R,
            self._stef_boltz,
            self._secs_per_day)

        _rrtmg_sw.initialise_rrtm_radiation(
            self._Cpd,
            self._solar_const,
            self._fac_sunspot_coeff,
            self._solar_var_by_band,
            self._cloud_overlap,
            self._cloud_optics,
            self._ice_props,
            self._liq_props,
            self._aerosol_type,
            self._solar_var_flag)
        super(RRTMGShortwave, self).__init__(**kwargs)

    @ensure_contiguous_state
    def array_call(self, state):
        """
        Get heating tendencies and shortwave fluxes.

        Args:

            state (dict):
                The model state dictionary.

        Returns:

            tendencies (dict), diagnostics (dict):

                * The shortwave heating tendency.
                * The upward/downward shortwave fluxes for cloudy and clear
                  sky conditions.

        """
        Q = mass_to_volume_mixing_ratio(state['specific_humidity'], 18.02)
        assert state['air_pressure'].shape[0] + 1 == state['air_pressure_on_interface_levels'].shape[0]

        Tint = get_interface_values(
            state['air_temperature'],
            state['surface_temperature'],
            state['air_pressure'],
            state['air_pressure_on_interface_levels']
        )

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )
        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, state, self.input_properties
        )
        model_time = state['time']
        if self._ignore_day_of_year:
            day_of_year = 0
        else:
            day_of_year = model_time.timetuple().tm_yday
        cos_zenith_angle = np.cos(state['zenith_angle'])

        _rrtmg_sw.rrtm_calculate_shortwave_fluxes(
            state['air_temperature'].shape[1],
            state['air_temperature'].shape[0],
            day_of_year,
            state['solar_cycle_fraction'],
            state['flux_adjustment_for_earth_sun_distance'],
            state['air_pressure'],
            state['air_pressure_on_interface_levels'],
            state['air_temperature'],
            Tint,
            state['surface_temperature'],
            Q,
            state['mole_fraction_of_ozone_in_air'],
            state['mole_fraction_of_carbon_dioxide_in_air'],
            state['mole_fraction_of_methane_in_air'],
            state['mole_fraction_of_nitrous_oxide_in_air'],
            state['mole_fraction_of_oxygen_in_air'],
            state['surface_albedo_for_direct_shortwave'],
            state['surface_albedo_for_direct_near_infrared'],
            state['surface_albedo_for_diffuse_shortwave'],
            state['surface_albedo_for_diffuse_near_infrared'],
            cos_zenith_angle,
            state['cloud_area_fraction_in_atmosphere_layer'],
            diagnostics['upwelling_shortwave_flux_in_air'],
            diagnostics['downwelling_shortwave_flux_in_air'],
            tendencies['air_temperature'],
            diagnostics['upwelling_shortwave_flux_in_air_assuming_clear_sky'],
            diagnostics['downwelling_shortwave_flux_in_air_assuming_clear_sky'],
            diagnostics['air_temperature_tendency_from_shortwave_assuming_clear_sky'],
            state['shortwave_optical_thickness_due_to_aerosol'],
            state['single_scattering_albedo_due_to_aerosol'],
            state['aerosol_asymmetry_parameter'],
            state['aerosol_optical_depth_at_55_micron'],
            state['shortwave_optical_thickness_due_to_cloud'],
            state['single_scattering_albedo_due_to_cloud'],
            state['cloud_asymmetry_parameter'],
            state['cloud_forward_scattering_fraction'],
            state['mass_content_of_cloud_ice_in_atmosphere_layer'],
            state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'],
            state['cloud_ice_particle_size'],
            state['cloud_water_droplet_radius'])

        diagnostics['air_temperature_tendency_from_shortwave'][:] = tendencies['air_temperature']

        return tendencies, diagnostics
