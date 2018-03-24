from sympl import get_numpy_array
from ...._core import (
    mass_to_volume_mixing_ratio, get_interface_values,
    ClimtPrognostic, numpy_version_of, get_constant)
import numpy as np
from numpy import pi as PI
from ..rrtmg_common import (
    rrtmg_cloud_overlap_method_dict, rrtmg_cloud_props_dict,
    rrtmg_cloud_ice_props_dict, rrtmg_cloud_liquid_props_dict)
try:
    from . import _rrtmg_lw
except ImportError:
    print('Import failed. RRTMG Longwave will not be available!')


class RRTMGLongwave(ClimtPrognostic):
    """
    The Rapid Radiative Transfer Model (RRTMG).

    This module wraps RRTMG for longwave radiation
    (i.e, emission from the earth's surface).
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
        'mole_fraction_of_cfc11_in_air': 'dimensionless',
        'mole_fraction_of_cfc12_in_air': 'dimensionless',
        'mole_fraction_of_cfc22_in_air': 'dimensionless',
        'mole_fraction_of_carbon_tetrachloride_in_air': 'dimensionless',
        'surface_longwave_emissivity': 'dimensionless',
        'cloud_area_fraction_in_atmosphere_layer': 'dimensionless',
        'longwave_optical_thickness_due_to_cloud': 'dimensionless',
        'mass_content_of_cloud_ice_in_atmosphere_layer': 'g m^-2',
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': 'g m^-2',
        'cloud_ice_particle_size': 'micrometer',
        'cloud_water_droplet_radius': 'micrometer',
        'longwave_optical_thickness_due_to_aerosol': 'dimensionless'
    }

    _climt_tendencies = {
        'air_temperature': 'degK day^-1'
    }

    _climt_diagnostics = {
        'upwelling_longwave_flux_in_air': 'W m^-2',
        'downwelling_longwave_flux_in_air': 'W m^-2',
        'upwelling_longwave_flux_in_air_assuming_clear_sky': 'W m^-2',
        'downwelling_longwave_flux_in_air_assuming_clear_sky': 'W m^-2',
        'longwave_heating_rate_assuming_clear_sky': 'degK day^-1',
        'longwave_heating_rate': 'degK day^-1',
        # TODO Need to add those final two quantities from the code
    }

    extra_dimensions = {'num_longwave_bands': np.arange(16)}
    '''
    RRTM without MCICA requires certain arrays on spectral bands
    '''

    quantity_descriptions = {
        'surface_longwave_emissivity': {
            'dims': ['x', 'y', 'num_longwave_bands'],
            'units': 'dimensionless',
            'default_value': 1.
        },

        'longwave_optical_thickness_due_to_cloud': {
            'dims': ['x', 'num_longwave_bands', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'default_value': 0.
        },

        'longwave_optical_thickness_due_to_aerosol': {
            'dims': ['x', 'y', 'mid_levels', 'num_longwave_bands'],
            'units': 'dimensionless',
            'default_value': 0.
        }
    }

    def __init__(
            self,
            calculate_change_up_flux=False,
            cloud_overlap_method='random',
            cloud_optical_properties='liquid_and_ice_clouds',
            cloud_ice_properties='ebert_curry_two',
            cloud_liquid_water_properties='radius_dependent_absorption',
            calculate_interface_temperature=True):

        """

        Args:

            calculate_change_up_flux (bool):
                calculate derivative of flux change with respect to
                surface temperature alone. Can be used to adjust fluxes in between radiation calls
                only due to change of surface temperature. Default value is :code:`False`, meaning this quantity
                is not calculated.

            cloud_overlap_method (string):
                Choose the method to do overlap with:

                * :code:`clear_only` = Clear only (no clouds)
                * :code:`random` = Random
                * :code:`maximum_random` = Maximum/Random
                * :code:`maximum` = Maximum.

            cloud_optical_properties (string):
                Choose how cloud optical properties are calculated:

                * :code:`direct_input` = Both cloud fraction and cloud optical depth are input directly
                * :code:`single_cloud_type` = Cloud fraction and cloud physical properties are input, ice
                  and liquid clouds are treated together, cloud absorptivity is a constant value (0.060241)
                * :code:`liquid_and_ice_clouds` = Cloud fraction and cloud physical properties are input, ice and liquid
                  clouds are treated separately

            cloud_ice_properties (string):
                set bounds on ice particle size.

                * :code:`ebert_curry_one` = ice particle has effective radius >= 10.0 micron `[Ebert and Curry 1992]`_
                * :code:`ebert_curry_two` = ice particle has effective radius between 13.0 and 130.0 micron `[Ebert and Curry 1992]`_
                * :code:`key_streamer_manual` = ice particle has effective radius between 5.0 and 131.0 micron
                  `[Key, Streamer Ref. Manual, 1996]`_
                * :code:`fu` = ice particle has generalised effective size (dge) between 5.0 and 140.0 micron
                  `[Fu, 1996]`_. (dge = 1.0315 * r_ec)

                Default value is 0.

            cloud_liquid_water_properties (string):
                set treatment of cloud liquid water.

                * :code:`radius_independent_absorption` = use radius independent absorption coefficient
                * :code:`radius_dependent_absorption` = use radius dependent absorption coefficient (radius between 2.5 and 60 micron)

            calculate_interface_temperature (bool):
                if :code:`True`, the interface temperature is calculated internally using a weighted
                interpolation routine. If :code:`False`, the quantity called
                :code:`air_temperature_on_interface_levels` in the input state needs to be manually
                updated by user code.

        .. _[Ebert and Curry 1992]:
            http://onlinelibrary.wiley.com/doi/10.1029/91JD02472/abstract

        .. _[Key, Streamer Ref. Manual, 1996]:
            https://stratus.ssec.wisc.edu/streamer/userman.pdf

        .. _[Fu, 1996]:
             http://journals.ametsoc.org/doi/abs/10.1175/1520-0442(1996)009%3C2058%3AAAPOTS%3E2.0.CO%3B2

            """

        if calculate_change_up_flux:
            self._calc_dflxdt = 1
        else:
            self._calc_dflxdt = 0

        self._cloud_overlap = rrtmg_cloud_overlap_method_dict[cloud_overlap_method.lower()]

        self._cloud_optics = rrtmg_cloud_props_dict[cloud_optical_properties.lower()]

        self._ice_props = rrtmg_cloud_ice_props_dict[cloud_ice_properties.lower()]

        self._liq_props = rrtmg_cloud_liquid_props_dict[cloud_liquid_water_properties.lower()]

        self._calc_Tint = calculate_interface_temperature

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

        # if self._cloud_optics == 0:  # Cloud optical depth directly input
        #     for input_quantity in ['mass_content_of_cloud_ice_in_atmosphere_layer',
        #                            'mass_content_of_cloud_liquid_water_in_atmosphere_layer',
        #                            'cloud_ice_particle_size',
        #                            'cloud_water_droplet_radius']:
        #         copy_inputs = self._climt_inputs.copy()
        #         copy_inputs.pop(input_quantity)
        #         self._climt_inputs = copy_inputs
        # else:
        #     copy_inputs = self._climt_inputs.copy()
        #     copy_inputs.pop('longwave_optical_thickness_due_to_cloud')
        #     self._climt_inputs = copy_inputs

        if not self._calc_Tint:
            self._climt_inputs['air_temperature_on_interface_levels'] = 'degK'

        _rrtmg_lw.set_constants(
            PI, self._g,
            self._planck,
            self._boltzmann,
            self._c,
            self._Na,
            self._loschmidt,
            self._R,
            self._stef_boltz,
            self._secs_per_day)

        # TODO Add all other flags as well
        _rrtmg_lw.initialise_rrtm_radiation(
            self._Cpd,
            self._cloud_overlap,
            self._calc_dflxdt,
            self._cloud_optics,
            self._ice_props,
            self._liq_props)

    def __call__(self, state):
        """
        Get heating tendencies and longwave fluxes.

        Args:

            state (dict):
                The model state dictionary.

        Returns:

            tendencies (dict), diagnostics (dict):

                * The longwave heating tendency.
                * The upward/downward longwave fluxes for cloudy and clear
                  sky conditions.

        """

        raw_arrays = self.get_numpy_arrays_from_state('_climt_inputs', state)
        Q = get_numpy_array(state['specific_humidity'].to_units('g/g'), ['x', 'y', 'z'])
        Q = mass_to_volume_mixing_ratio(Q, 18.02)

        mid_level_shape = raw_arrays['air_temperature'].shape

        if self._calc_Tint:
            Tint = get_interface_values(raw_arrays['air_temperature'],
                                        raw_arrays['surface_temperature'],
                                        raw_arrays['air_pressure'],
                                        raw_arrays['air_pressure_on_interface_levels'])
        else:
            Tint = raw_arrays['air_temperature_on_interface_levels']

        diag_dict = self.create_state_dict_for('_climt_diagnostics', state)

        diag_arrays = numpy_version_of(diag_dict)

        tend_dict = self.create_state_dict_for('_climt_tendencies', state)

        tend_arrays = numpy_version_of(tend_dict)

        # TODO add dflx_dt as well

        raw_f_arrays = {}
        diag_f_arrays = {}
        tend_f_arrays = {}
        for quantity in raw_arrays.keys():
            raw_f_arrays[quantity] = np.asfortranarray(raw_arrays[quantity].transpose())

        for quantity in diag_arrays.keys():
            diag_f_arrays[quantity] = np.asfortranarray(diag_arrays[quantity].transpose())

        for quantity in tend_arrays.keys():
            tend_f_arrays[quantity] = np.asfortranarray(tend_arrays[quantity].transpose())

        Tint_f = np.asfortranarray(Tint.transpose())
        Q_f = np.asfortranarray(Q.transpose())

        _rrtmg_lw.rrtm_calculate_longwave_fluxes(
            mid_level_shape[0],
            mid_level_shape[1],
            mid_level_shape[2],
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
            raw_f_arrays['mole_fraction_of_cfc11_in_air'],
            raw_f_arrays['mole_fraction_of_cfc12_in_air'],
            raw_f_arrays['mole_fraction_of_cfc22_in_air'],
            raw_f_arrays['mole_fraction_of_carbon_tetrachloride_in_air'],
            raw_f_arrays['surface_longwave_emissivity'],
            raw_f_arrays['cloud_area_fraction_in_atmosphere_layer'],
            raw_f_arrays['longwave_optical_thickness_due_to_aerosol'],
            diag_f_arrays['upwelling_longwave_flux_in_air'],
            diag_f_arrays['downwelling_longwave_flux_in_air'],
            tend_f_arrays['air_temperature'],
            diag_f_arrays['upwelling_longwave_flux_in_air_assuming_clear_sky'],
            diag_f_arrays['downwelling_longwave_flux_in_air_assuming_clear_sky'],
            diag_f_arrays['longwave_heating_rate_assuming_clear_sky'],
            raw_f_arrays['longwave_optical_thickness_due_to_cloud'],
            raw_f_arrays['mass_content_of_cloud_ice_in_atmosphere_layer'],
            raw_f_arrays['mass_content_of_cloud_liquid_water_in_atmosphere_layer'],
            raw_f_arrays['cloud_ice_particle_size'],
            raw_f_arrays['cloud_water_droplet_radius'])

        for quantity in diag_dict.keys():
            diag_dict[quantity].values = diag_f_arrays[quantity].transpose()

        for quantity in tend_dict.keys():
            tend_dict[quantity].values = tend_f_arrays[quantity].transpose()

        diag_dict['longwave_heating_rate'].values = tend_dict['air_temperature'].values

        return tend_dict, diag_dict
