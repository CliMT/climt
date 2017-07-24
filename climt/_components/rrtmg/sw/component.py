from sympl import (replace_none_with_default,
                   get_numpy_array)
from ...._core import (
    mass_to_volume_mixing_ratio, get_interface_values,
    ClimtPrognostic, numpy_version_of)
import numpy as np
from numpy import pi as numpy_pi
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
                        'num_ecmwf_aerosols': np.arange(6)}
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
            'dims': [],
            'units': 'dimensionless',
            'default_value': 0.0
        },

        'flux_adjustment_for_earth_sun_distance': {
            'dims': [],
            'units': 'dimensionless',
            'default_value': 1.0
        },
    }

    def __init__(
            self,
            cloud_overlap_method=1,
            cloud_optical_properties=2,
            cloud_ice_properties=1,
            cloud_liquid_water_properties=1,
            solar_variability_method=0,
            solar_constant=0.,
            facular_sunspot_amplitude=np.ones(2),
            solar_variability_by_band=np.ones(16),
            aerosol_type=0,
            acceleration_gravity=None,
            planck_constant=None,
            boltzmann_constant=None,
            speed_of_light=None,
            avogadro_constant=None,
            loschmidt_constant=None,
            universal_gas_constant=None,
            stefan_boltzmann_constant=None,
            seconds_per_day=None,
            specific_heat_dry_air=None):

        """

        Args:

            cloud_overlap_method (int):
                Choose the method to do overlap with:

                * 0 = Clear only (no clouds)
                * 1 = Random
                * 2 = Maximum/Random
                * 3 = Maximum.

                Default value is 1, corresponding to random cloud overlap.

            cloud_optical_properties (int):
                Choose how cloud optical properties are calculated:

                * 0 = Both cloud fraction and cloud optical depth are input directly
                * 1 = Cloud fraction and cloud physical properties are input, ice and liquid clouds are
                  treated together, cloud absorptivity is a constant value (0.060241)
                * 2 = Cloud fraction and cloud physical properties are input, ice and liquid
                  clouds are treated separately

                Default value is 0.

            cloud_ice_properties (int):
                set bounds on ice particle size.

                * 0 = ice particle has effective radius >= 10.0 micron `[Ebert and Curry 1992]`_
                * 1 = ice particle has effective radius between 13.0 and 130.0 micron `[Ebert and Curry 1992]`_
                * 2 = ice particle has effective radius between 5.0 and 131.0 micron
                  `[Key, Streamer Ref. Manual, 1996]`_
                * 3 = ice particle has generalised effective size (dge) between 5.0 and 140.0 micron
                  `[Fu, 1996]`_. (dge = 1.0315 * r_ec)

                Default value is 0.

            cloud_liquid_water_properties (int):
                set treatment of cloud liquid water.

                * 0 = use radius independent absorption coefficient
                * 1 = use radius dependent absorption coefficient (radius between 2.5 and 60 micron)

                Default value is 0.

            solar_variability_method (int):
                set the solar variability model used by RRTMG.

                * solar_variability_method = -1:

                    * If :code:`solar_constant = 0`: No solar variability and no solar cycle
                      with a solar constant of 1368.22 :math:`W m^{-2}`.
                    * If :code:`solar_constant != 0`: Solar variability defined by setting
                      non-zero scale factors in :code:`solar_variability_by_band`.

                * solar_variability_method = 0:

                    * If :code:`solar_constant = 0`: No solar variability and no solar cycle
                      with a solar constant of 1360.85 :math:`W m^{-2}`, with facular and
                      sunspot effects fixed to the mean of solar cycles 13-24.
                    * If :code:`solar_constant != 0`: No solar variability and no solar cycle.

                * solar_variability_method = 1: Solar variability using the NRLSSI2 solar model
                  with solar cycle contribution determined by :code:`solar_cycle_fraction` in
                  the model state, and facular and sunspot adjustment scale factors specified
                  in :code:`facular_sunspot_amplitude`.

                * solar_variability_method = 2: Solar variability using the NRLSSI2 solar model
                  using solar cycle determined by direct specification of **Mg** (facular)
                  and **SB** (sunspot) indices provided in :code:`facular_sunspot_amplitude`.
                  :code:`solar_constant` is ignored.

                * solar_variability_method = 3:
                 * If :code:`solar_constant = 0`: No solar variability and no solar cycle
                   with a solar constant of 1360.85 :math:`W m^{-2}`.
                 * If :code:`solar_constant != 0`: scale factors in :code:`solar_variability_by_band`.

            solar_constant (float):
                Solar constant -- solar irradiance averaged over a solar cycle -- in units of
                :code:`W m^{-2}`.

            facular_sunspot_amplitude (array of dimension 2):
                Facular and Sunspot amplitude variability parameters, described previously.

            solar_variability_by_band (array of dimension 14 = number of spectral bands):
                scale factors for solar variability in all spectral bands.

            aerosol_type (int):
                Type of aerosol inputs to RRTMG.

                * 0: No Aerosol.
                * 6 -- ECMWF method. Requires aerosol optical depth at 55 micron as the
                  state quantity :code:`aerosol_optical_depth_at_55_micron`.
                * 10: Input all aerosol optical properties.

            acceleration_gravity (float):
                value of acceleration due to gravity in
                :math:`m s^{-1}`. Default value from :code:`sympl.default.constants` is used if None.

            planck_constant (float):
                value of the planck constant in :math:`J s`.
                Default value from :code:`sympl.default.constants` is used if None.

            boltzmann_constant (float):
                value of the Boltzmann constant in :math:`J K^{-1}`.
                Default value from :code:`sympl.default.constants` is used if None.

            speed_of_light (float):
                value of the speed of light in :math:`m s^{-1}`.
                Default value from :code:`sympl.default.constants` is used if None.

            avogadro_constant (float):
                value of the Avogadro constant.
                Default value from :code:`sympl.default.constants` is used if None.

            loschmidt_constant (float):
                value of the Loschmidt constant.
                Default value from :code:`sympl.default.constants` is used if None.

            universal_gas_constant (float):
                value of the gas constant in :math:`J K^{-1} mol^{-1}`.
                Default value from :code:`sympl.default.constants` is used if None.

            stefan_boltzmann_constant (float):
                value of the Stefan-Boltzmann constant
                in :math:`W m^{-2} K^{-4}`. Default value from :code:`sympl.default.constants` is
                used if None.

            seconds_per_day (float):
                number of seconds per day.
                Default value from :code:`sympl.default.constants` (for earth) is used if None.

            specific_heat_dry_air (float):
                The specific heat of dry air in :math:`J K^{-1} kg^{-1}`.
                Default value from :code:`sympl.default.constants` is used if None.

        .. _[Ebert and Curry 1992]:
            http://onlinelibrary.wiley.com/doi/10.1029/91JD02472/abstract

        .. _[Key, Streamer Ref. Manual, 1996]:
            https://stratus.ssec.wisc.edu/streamer/userman.pdf

        .. _[Fu, 1996]:
             http://journals.ametsoc.org/doi/abs/10.1175/1520-0442(1996)009%3C2058%3AAAPOTS%3E2.0.CO%3B2

            """

        self._cloud_overlap = cloud_overlap_method

        self._cloud_optics = cloud_optical_properties

        self._ice_props = cloud_ice_properties

        self._liq_props = cloud_liquid_water_properties

        self._solar_var_flag = solar_variability_method

        self._solar_const = solar_constant

        self._fac_sunspot_coeff = facular_sunspot_amplitude

        self._solar_var_by_band = solar_variability_by_band

        self._aerosol_type = aerosol_type

        self._g = replace_none_with_default(
            'gravitational_acceleration', acceleration_gravity)

        self._planck = replace_none_with_default(
            'planck_constant', planck_constant).to_units('erg s')

        self._boltzmann = replace_none_with_default(
            'boltzmann_constant', boltzmann_constant).to_units('erg K^-1')

        self._c = replace_none_with_default(
            'speed_of_light', speed_of_light).to_units('cm s^-1')

        self._Na = replace_none_with_default(
            'avogadro_constant', avogadro_constant)

        self._loschmidt = replace_none_with_default(
            'loschmidt_constant', loschmidt_constant).to_units('cm^-3')

        self._R = replace_none_with_default(
            'universal_gas_constant', universal_gas_constant).to_units('erg mol^-1 K^-1')

        self._stef_boltz = replace_none_with_default(
            'stefan_boltzmann_constant', stefan_boltzmann_constant).to_units('W cm^-2 K^-4')

        self._secs_per_day = replace_none_with_default(
            'seconds_per_day', seconds_per_day)

        self._Cpd = replace_none_with_default(
            'heat_capacity_of_dry_air_at_constant_pressure', specific_heat_dry_air)

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

        day_of_year = model_time.timetuple().tm_yday

        cos_zenith_angle = np.asfortranarray(np.cos(raw_arrays['zenith_angle']))

        for lon in range(mid_level_shape[0]):

            raw_f_arrays = {}
            diag_f_arrays = {}
            tend_f_arrays = {}
            for quantity in raw_arrays.keys():
                if quantity not in ['flux_adjustment_for_earth_sun_distance',
                                    'solar_cycle_fraction']:
                    raw_f_arrays[quantity] = np.asfortranarray(raw_arrays[quantity][lon, :])
                else:
                    raw_f_arrays[quantity] = raw_arrays[quantity]
            for quantity in diag_arrays.keys():
                diag_f_arrays[quantity] = np.asfortranarray(diag_arrays[quantity][lon, :])
            for quantity in tend_arrays.keys():
                tend_f_arrays[quantity] = np.asfortranarray(tend_arrays[quantity][lon, :])

            Tint_f = np.asfortranarray(Tint[lon, :])
            Q_f = np.asfortranarray(Q[lon, :])
            zenith_f = np.asfortranarray(cos_zenith_angle[lon, :])

            _rrtmg_sw.rrtm_calculate_shortwave_fluxes(
                mid_level_shape[1],
                mid_level_shape[2],
                day_of_year,
                raw_f_arrays['solar_cycle_fraction'].item(),
                raw_f_arrays['flux_adjustment_for_earth_sun_distance'].item(),
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
                diag_arrays[quantity][lon, :] = diag_f_arrays[quantity]
            for quantity in tend_arrays.keys():
                tend_arrays[quantity][lon, :] = tend_f_arrays[quantity]

            tend_arrays['air_temperature'][lon, :, -1] = 0

        diag_dict['shortwave_heating_rate'].values[:] = tend_dict['air_temperature'].values

        return tend_dict, diag_dict
