from sympl import get_numpy_array
from ...._core import (
    mass_to_volume_mixing_ratio, get_interface_values,
    ClimtPrognostic, numpy_version_of, get_constant)
import numpy as np
from numpy import pi as numpy_pi
from ..rrtmg_common import (
    rrtmg_cloud_overlap_method_dict, rrtmg_cloud_props_dict,
    rrtmg_cloud_ice_props_dict, rrtmg_cloud_liquid_props_dict,
    rrtmg_aerosol_input_dict)
try:
    from . import _rrtmg_sw
except ImportError:
    print('Import failed. RRTMG Shortwave will not be available!')


class RRTMGShortwave(ClimtPrognostic):
    """
    The Rapid Radiative Transfer Model (RRTMG).

    This module wraps RRTMG for shortwave radiation
    (i.e, emission from the sun).

    """

    _climt_inputs = {
        'air_pressure': 'mbar',
        'air_pressure_on_interface_levels': 'mbar',
        'air_temperature': 'degK',
        'surface_temperature': 'degK',
        'specific_humidity': 'g/g',
        'mole_fraction_of_ozone_in_air': 'dimensionless',
        'mole_fraction_of_carbon_dioxide_in_air': 'dimensionless',
        'mole_fraction_of_methane_in_air': 'dimensionless',
        'mole_fraction_of_nitrous_oxide_in_air': 'dimensionless',
        'mole_fraction_of_oxygen_in_air': 'dimensionless',
        'cloud_area_fraction_in_atmosphere_layer': 'dimensionless',
        'shortwave_optical_thickness_due_to_cloud': 'dimensionless',
        'mass_content_of_cloud_ice_in_atmosphere_layer': 'g m^-2',
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': 'g m^-2',
        'cloud_ice_particle_size': 'micrometer',
        'cloud_water_droplet_radius': 'micrometer',
        'shortwave_optical_thickness_due_to_aerosol': 'dimensionless',
        'zenith_angle': 'radians',
        'flux_adjustment_for_earth_sun_distance': 'dimensionless',
        'surface_albedo_for_direct_shortwave': 'dimensionless',
        'surface_albedo_for_direct_near_infrared': 'dimensionless',
        'surface_albedo_for_diffuse_near_infrared': 'dimensionless',
        'surface_albedo_for_diffuse_shortwave': 'dimensionless',
        'single_scattering_albedo_due_to_cloud': 'dimensionless',
        'single_scattering_albedo_due_to_aerosol': 'dimensionless',
        'cloud_asymmetry_parameter': 'dimensionless',
        'aerosol_asymmetry_parameter': 'dimensionless',
        'cloud_forward_scattering_fraction': 'dimensionless',
        'aerosol_optical_depth_at_55_micron': 'dimensionless',
        'solar_cycle_fraction': 'dimensionless',
    }

    _climt_tendencies = {
        'air_temperature': 'degK day^-1'
    }

    _climt_diagnostics = {
        'upwelling_shortwave_flux_in_air': 'W m^-2',
        'downwelling_shortwave_flux_in_air': 'W m^-2',
        'upwelling_shortwave_flux_in_air_assuming_clear_sky': 'W m^-2',
        'downwelling_shortwave_flux_in_air_assuming_clear_sky': 'W m^-2',
        'shortwave_heating_rate_assuming_clear_sky': 'degK day^-1',
        'shortwave_heating_rate': 'degK day^-1',
    }

    extra_dimensions = {'num_shortwave_bands': np.arange(14),
                        'num_ecmwf_aerosols': np.arange(6),
                        'scalar': np.arange(1)}
    '''
    RRTM without MCICA requires certain arrays on spectral bands
    '''

    quantity_descriptions = {
        'shortwave_optical_thickness_due_to_cloud': {
            'dims': ['x', 'num_shortwave_bands', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'default_value': 0.
        },

        'cloud_asymmetry_parameter': {
            'dims': ['x', 'num_shortwave_bands', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'default_value': 0.85
        },

        'cloud_forward_scattering_fraction': {
            'dims': ['x', 'num_shortwave_bands', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'default_value': 0.8
        },

        'single_scattering_albedo_due_to_cloud': {
            'dims': ['x', 'num_shortwave_bands', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'default_value': 0.9
        },

        'shortwave_optical_thickness_due_to_aerosol': {
            'dims': ['x', 'y', 'mid_levels', 'num_shortwave_bands'],
            'units': 'dimensionless',
            'default_value': 0.
        },

        'aerosol_asymmetry_parameter': {
            'dims': ['x', 'y', 'mid_levels', 'num_shortwave_bands'],
            'units': 'dimensionless',
            'default_value': 0.
        },

        'single_scattering_albedo_due_to_aerosol': {
            'dims': ['x', 'y', 'mid_levels', 'num_shortwave_bands'],
            'units': 'dimensionless',
            'default_value': 0.5
        },

        'aerosol_optical_depth_at_55_micron': {
            'dims': ['x', 'y', 'mid_levels', 'num_ecmwf_aerosols'],
            'units': 'dimensionless',
            'default_value': 0.0
        },

        'solar_cycle_fraction': {
            'dims': ['scalar'],
            'units': 'dimensionless',
            'default_value': 0.0
        },

        'flux_adjustment_for_earth_sun_distance': {
            'dims': ['scalar'],
            'units': 'dimensionless',
            'default_value': 1.0
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
            aerosol_type='no_aerosol'):

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

    def __call__(self, state):
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

        raw_arrays = self.get_numpy_arrays_from_state('_climt_inputs', state)

        Q = get_numpy_array(state['specific_humidity'].to_units('g/g'), ['x', 'y', 'z'])
        Q = mass_to_volume_mixing_ratio(Q, 18.02)

        mid_level_shape = raw_arrays['air_temperature'].shape

        Tint = get_interface_values(raw_arrays['air_temperature'],
                                    raw_arrays['surface_temperature'],
                                    raw_arrays['air_pressure'],
                                    raw_arrays['air_pressure_on_interface_levels'])

        diag_dict = self.create_state_dict_for('_climt_diagnostics', state)

        diag_arrays = self.get_numpy_arrays_from_state('_climt_diagnostics', diag_dict)

        tend_dict = self.create_state_dict_for('_climt_tendencies', state)

        tend_arrays = numpy_version_of(tend_dict)

        model_time = state['time']

        if self._ignore_day_of_year:
            day_of_year = 0
        else:
            day_of_year = model_time.timetuple().tm_yday

        cos_zenith_angle = np.asfortranarray(np.cos(raw_arrays['zenith_angle']))

        raw_f_arrays = {}
        diag_f_arrays = {}
        tend_f_arrays = {}
        for quantity in raw_arrays.keys():
            if quantity not in ['flux_adjustment_for_earth_sun_distance',
                                'solar_cycle_fraction']:
                raw_f_arrays[quantity] = np.asfortranarray(raw_arrays[quantity].transpose())
            else:
                raw_f_arrays[quantity] = raw_arrays[quantity]

        for quantity in diag_arrays.keys():
            diag_f_arrays[quantity] = np.asfortranarray(diag_arrays[quantity].transpose())

        for quantity in tend_arrays.keys():
            tend_f_arrays[quantity] = np.asfortranarray(tend_arrays[quantity].transpose())

        Tint_f = np.asfortranarray(Tint.transpose())
        Q_f = np.asfortranarray(Q.transpose())
        zenith_f = np.asfortranarray(cos_zenith_angle.transpose())

        _rrtmg_sw.rrtm_calculate_shortwave_fluxes(
            mid_level_shape[0],
            mid_level_shape[1],
            mid_level_shape[2],
            day_of_year,
            raw_f_arrays['solar_cycle_fraction'],
            raw_f_arrays['flux_adjustment_for_earth_sun_distance'],
            raw_f_arrays['air_pressure'],
            raw_f_arrays['air_pressure_on_interface_levels'],
            raw_f_arrays['air_temperature'],
            Tint_f,
            raw_f_arrays['surface_temperature'],
            Q_f,
            raw_f_arrays['mole_fraction_of_ozone_in_air'],
            raw_f_arrays['mole_fraction_of_carbon_dioxide_in_air'],
            raw_f_arrays['mole_fraction_of_methane_in_air'],
            raw_f_arrays['mole_fraction_of_nitrous_oxide_in_air'],
            raw_f_arrays['mole_fraction_of_oxygen_in_air'],
            raw_f_arrays['surface_albedo_for_direct_shortwave'],
            raw_f_arrays['surface_albedo_for_direct_near_infrared'],
            raw_f_arrays['surface_albedo_for_diffuse_shortwave'],
            raw_f_arrays['surface_albedo_for_diffuse_near_infrared'],
            zenith_f,
            raw_f_arrays['cloud_area_fraction_in_atmosphere_layer'],
            diag_f_arrays['upwelling_shortwave_flux_in_air'],
            diag_f_arrays['downwelling_shortwave_flux_in_air'],
            tend_f_arrays['air_temperature'],
            diag_f_arrays['upwelling_shortwave_flux_in_air_assuming_clear_sky'],
            diag_f_arrays['downwelling_shortwave_flux_in_air_assuming_clear_sky'],
            diag_f_arrays['shortwave_heating_rate_assuming_clear_sky'],
            raw_f_arrays['shortwave_optical_thickness_due_to_aerosol'],
            raw_f_arrays['single_scattering_albedo_due_to_aerosol'],
            raw_f_arrays['aerosol_asymmetry_parameter'],
            raw_f_arrays['aerosol_optical_depth_at_55_micron'],
            raw_f_arrays['shortwave_optical_thickness_due_to_cloud'],
            raw_f_arrays['single_scattering_albedo_due_to_cloud'],
            raw_f_arrays['cloud_asymmetry_parameter'],
            raw_f_arrays['cloud_forward_scattering_fraction'],
            raw_f_arrays['mass_content_of_cloud_ice_in_atmosphere_layer'],
            raw_f_arrays['mass_content_of_cloud_liquid_water_in_atmosphere_layer'],
            raw_f_arrays['cloud_ice_particle_size'],
            raw_f_arrays['cloud_water_droplet_radius'])

        for quantity in diag_arrays.keys():
            diag_dict[quantity].values = diag_f_arrays[quantity].transpose()
        for quantity in tend_arrays.keys():
            tend_dict[quantity].values = tend_f_arrays[quantity].transpose()

        diag_dict['shortwave_heating_rate'].values[:] = tend_dict['air_temperature'].values

        return tend_dict, diag_dict
