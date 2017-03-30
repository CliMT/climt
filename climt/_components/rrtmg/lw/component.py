from sympl import (Prognostic,
                   DataArray,
                   replace_none_with_default,
                   get_numpy_array,
                   combine_dimensions)
import climt
from . import _rrtm_lw
import numpy as np
from numpy import pi as PI


class RRTMLongwave(Prognostic):
    """
    Interface to the Rapid Radiative Transfer Model (RRTM) for
    longwave radiation (i.e, emission from the earth's surface).
    """

    inputs = {
        'air_pressure': 'mbar',
        'air_pressure_on_interface_levels': 'mbar',
        'air_temperature': 'degK',
        'air_temperature_on_interface_levels': 'degK',
        'surface_temperature': 'degK',
        'specific_humidity': 'g/g',
        'ozone_mixing_ratio': 'dimensionless',
        'carbon_dioxide_mixing_ratio': 'dimensionless',
        'methane_mixing_ratio': 'dimensionless',
        'nitrous_oxide_mixing_ratio': 'dimensionless',
        'oxygen_mixing_ratio': 'dimensionless',
        'cfc11_mixing_ratio': 'dimensionless',
        'cfc12_mixing_ratio': 'dimensionless',
        'cfc22_mixing_ratio': 'dimensionless',
        'ccl4_mixing_ratio': 'dimensionless',
        'surface_emissivity': 'dimensionless',
        'cloud_fraction': 'dimensionless',
        'cloud_optical_depth': 'dimensionless',
        'cloud_ice_water_path': 'g m^-2',
        'cloud_liquid_water_path': 'g m^-2',
        'cloud_ice_particle_size': 'micrometer',
        'cloud_water_droplet_radius': 'micrometer',
        'aerosol_optical_depth': 'dimensionless'
    }

    tendencies = {
        'longwave_heating_rate': 'K day^-1'
    }

    diagnostics = {
        'upward_longwave_flux': 'W m^-2',
        'downward_longwave_flux': 'W m^-2',
        'upward_longwave_flux_clearsky': 'W m^-2',
        'downward_longwave_flux_clearsky': 'W m^-2',
        'longwave_heating_rate_clearsky': 'K day^-1',
        # TODO Need to add those final two quantities from the code
    }

    '''
    RRTM without MCICA requires certain arrays on spectral bands
    '''
    extra_dimensions = {'num_longwave_bands': np.arange(16)}

    quantity_descriptions = {
        'surface_emissivity': {
            'dims': ['x', 'y', 'num_longwave_bands'],
            'units': 'dimensionless',
            'init_value': 1.
        },

        'cloud_optical_depth': {
            'dims': ['x', 'num_longwave_bands', 'y', 'mid_levels'],
            'units': 'dimensionless',
            'init_value': 0.
        },

        'aerosol_optical_depth': {
            'dims': ['x', 'y', 'mid_levels', 'num_longwave_bands'],
            'units': 'dimensionless',
            'init_value': 0.
        }
    }

    def __init__(
            self,
            calculate_change_up_flux=0,
            cloud_overlap_method=1,
            cloud_optical_properties=0,
            cloud_ice_properties=0,
            cloud_liquid_water_properties=0,
            gravitational_acceleration=None,
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

        calculate_change_up_flux (int): calculate derivative of flux change with respect to
        surface temperature alone. Can be used to adjust fluxes in between radiation calls
        only due to change of surface temperature. Default value is 0, meaning this quantity
        is not calculated.

        cloud_overlap_method (int): Choose the method to do overlap with:
        0 = Clear only (no clouds)
        1 = Random
        2 = Maximum/Random
        3 = Maximum.

        Default value is 1, corresponding to random cloud overlap.

        cloud_optical_properties (int): Choose how cloud optical properties are
        calculated.
        0 = Both cloud fraction and cloud optical depth are input directly
        1 = Cloud fraction and cloud physical properties are input, ice and liquid
            clouds are treated together, cloud absorptivity is a constant value (0.060241)
        2 = Cloud fraction and cloud physical properties are input, ice and liquid
            clouds are treated separately

        Default value is 0.

        cloud_ice_properties (int): set bounds on ice particle size.
        0 = ice particle has effective radius >= 10.0 micron (Ebert and Curry 1992)
        1 = ice particle has effective radius between 13.0 and 130.0 micron (Ebert and Curry 1992)
        2 = ice particle has effective radius between 5.0 and 131.0 micron (Key, Streamer Ref.
            Manual, 1996)
        3 = ice particle has generalised effective size (dge) between 5.0 and 140.0 micron (Fu, 1996)
            dge = 1.0315* effective radius

        Default value is 0.

        cloud_liquid_water_properties (int): set treatment of cloud liquid water.
        0 = use radius independent absorption coefficient
        1 = use radius dependent absorption coefficient (radius between 2.5 and 60 micron)

        Default value is 0.

        gravitational_acceleration (float): value of acceleration due to gravity in
                $m s^{-1}$. Default value from climt.default_constants is used if None.

        planck_constant (float): value of the planck constant in $J s$.
                Default value from climt.default_constants is used if None.

        boltzmann_constant (float): value of the Boltzmann constant in $J {K^-1}$.
                Default value from climt.default_constants is used if None.

        speed_of_light (float): value of the speed of light in $m {s^-1}$.
                Default value from climt.default_constants is used if None.

        avogadro_constant (float): value of the Avogadro number.
                Default value from climt.default_constants is used if None.

        loschmidt_constant (float): value of the Loschmidt constant.
                Default value from climt.default_constants is used if None.

        universal_gas_constant (float): value of the gas constant in $J {K^-1} mol^{-1}$.
                Default value from climt.default_constants is used if None.

        stefan_boltzmann_constant (float): value of the Stefan-Boltzmann constant
                in $W m^{-2} K^{-4}$. Default value from climt.default_constants is
                used if None.

        seconds_per_day (float): number of seconds per day.
                Default value from climt.default_constants (for earth) is used if None.

        specific_heat_dry_air (float): The specific heat of dry air in $J {K^-1} kg^{-1}$.
                Default value from climt.default_constants is used if None.
        """


        self._calc_dflxdt = calculate_change_up_flux

        self._cloud_overlap = cloud_overlap_method

        self._cloud_optics = cloud_optical_properties

        self._ice_props = cloud_ice_properties

        self._liq_props = cloud_liquid_water_properties

        self._g = replace_none_with_default(
            'gravitational_acceleration', gravitational_acceleration)

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

        if self._cloud_optics == 0:  # Cloud optical depth directly input
            for input_quantity in ['cloud_ice_water_path', 'cloud_liquid_water_path',
                                   'cloud_ice_particle_size', 'cloud_water_droplet_radius']:
                copy_inputs = self.inputs.copy()
                copy_inputs.pop(input_quantity)
                self.inputs = copy_inputs
        else:
            copy_inputs = self.inputs.copy()
            copy_inputs.pop('cloud_optical_depth')
            self.inputs = copy_inputs

        _rrtm_lw.set_constants(
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
        _rrtm_lw.initialise_rrtm_radiation(
            self._Cpd,
            self._cloud_overlap,
            self._calc_dflxdt,
            self._cloud_optics,
            self._ice_props,
            self._liq_props)

    def __call__(self, state):
        """
        Returns radiative heating tendencies and longwave flux diagnostics.

        Args:
            state (dict) : The model state dictionary

        Returns:
            tendencies (dict) : The longwave heating tendency

            diagnostics (dict) : The upward/downward longwave fluxes for cloudy and clear
                                sky conditions.

        """

        P = get_numpy_array(state['air_pressure'].to_units('mbar'), ['x', 'y', 'z'])

        Pint = get_numpy_array(state['air_pressure_on_interface_levels'].to_units('mbar'),
                               ['x', 'y', 'z'])

        T = get_numpy_array(state['air_temperature'].to_units('degK'), ['x', 'y', 'z'])

        Tint = get_numpy_array(state['air_temperature_on_interface_levels'].to_units('degK'),
                               ['x', 'y', 'z'])

        Ts = get_numpy_array(state['surface_temperature'].to_units('degK'), ['x', 'y'])

        Q = climt.mass_to_volume_mixing_ratio(state['specific_humidity'].to_units('g/g'), 18.02)
        Q = get_numpy_array(Q, ['x', 'y', 'z'])

        O3 = get_numpy_array(state['ozone_mixing_ratio'].to_units('dimensionless'),
                             ['x', 'y', 'z'])

        CO2 = get_numpy_array(state['carbon_dioxide_mixing_ratio'].to_units('dimensionless'),
                              ['x', 'y', 'z'])

        CH4 = get_numpy_array(state['methane_mixing_ratio'].to_units('dimensionless'),
                              ['x', 'y', 'z'])

        N2O = get_numpy_array(state['nitrous_oxide_mixing_ratio'].to_units('dimensionless'),
                              ['x', 'y', 'z'])

        O2 = get_numpy_array(state['oxygen_mixing_ratio'].to_units('dimensionless'),
                             ['x', 'y', 'z'])

        CFC11 = get_numpy_array(state['cfc11_mixing_ratio'].to_units('dimensionless'),
                                ['x', 'y', 'z'])

        CFC12 = get_numpy_array(state['cfc12_mixing_ratio'].to_units('dimensionless'),
                                ['x', 'y', 'z'])

        CFC22 = get_numpy_array(state['cfc22_mixing_ratio'].to_units('dimensionless'),
                                ['x', 'y', 'z'])

        CCL4 = get_numpy_array(state['ccl4_mixing_ratio'].to_units('dimensionless'),
                               ['x', 'y', 'z'])

        Emiss = get_numpy_array(state['surface_emissivity'].to_units('dimensionless'),
                                ['x', 'y', 'num_longwave_bands'])

        ClFrac = get_numpy_array(state['cloud_fraction'].to_units('dimensionless'),
                                 ['x', 'y', 'z'])

        if self._cloud_optics == 0:  # Optical depth is part of input
            ClTau = get_numpy_array(state['cloud_optical_depth'].to_units('dimensionless'),
                                    ['x', 'num_longwave_bands', 'y', 'z'])
        else:
            ClIWP = get_numpy_array(state['cloud_ice_water_path'].to_units('g m^-2'),
                                    ['x', 'y', 'z'])

            ClLWP = get_numpy_array(state['cloud_liquid_water_path'].to_units('g m^-2'),
                                    ['x', 'y', 'z'])

            ClIceSize = get_numpy_array(state['ice_particle_size'].to_units('micrometer'),
                                        ['x', 'y', 'z'])

            ClDropSize = get_numpy_array(state['cloud_droplet_radius'].to_units('micrometer'),
                                         ['x', 'y', 'z'])

        AerTau = get_numpy_array(state['aerosol_optical_depth'].to_units('dimensionless'),
                                 ['x', 'y', 'z', 'num_longwave_bands'])

        up_flux = np.zeros(Tint.shape, order='F')
        down_flux = np.zeros(Tint.shape, order='F')
        heating_rate = np.zeros(T.shape, order='F')
        up_flux_clear = np.zeros(Tint.shape, order='F')
        down_flux_clear = np.zeros(Tint.shape, order='F')
        heating_rate_clear = np.zeros(T.shape, order='F')

        # TODO add dflx_dt as well
        for lon in range(T.shape[0]):
            if self._cloud_optics == 0:
                _rrtm_lw.rrtm_calculate_longwave_fluxes(T.shape[1],
                                                        T.shape[2],
                                                        P[lon, :],
                                                        Pint[lon, :],
                                                        T[lon, :],
                                                        Tint[lon, :],
                                                        Ts[lon, :],
                                                        Q[lon, :],
                                                        O3[lon, :],
                                                        CO2[lon, :],
                                                        CH4[lon, :],
                                                        N2O[lon, :],
                                                        O2[lon, :],
                                                        CFC11[lon, :],
                                                        CFC12[lon, :],
                                                        CFC22[lon, :],
                                                        CCL4[lon, :],
                                                        Emiss[lon, :],
                                                        ClFrac[lon, :],
                                                        AerTau[lon, :],
                                                        up_flux[lon, :],
                                                        down_flux[lon, :],
                                                        heating_rate[lon, :],
                                                        up_flux_clear[lon, :],
                                                        down_flux_clear[lon, :],
                                                        heating_rate_clear[lon, :],
                                                        ClTau[lon, :])
            else:
                _rrtm_lw.rrtm_calculate_longwave_fluxes(T.shape[1],
                                                        T.shape[2],
                                                        P[lon, :],
                                                        Pint[lon, :],
                                                        T[lon, :],
                                                        Tint[lon, :],
                                                        Ts[lon, :],
                                                        Q[lon, :],
                                                        O3[lon, :],
                                                        CO2[lon, :],
                                                        CH4[lon, :],
                                                        N2O[lon, :],
                                                        O2[lon, :],
                                                        CFC11[lon, :],
                                                        CFC12[lon, :],
                                                        CFC22[lon, :],
                                                        CCL4[lon, :],
                                                        Emiss[lon, :],
                                                        ClFrac[lon, :],
                                                        AerTau[lon, :],
                                                        up_flux[lon, :],
                                                        down_flux[lon, :],
                                                        heating_rate[lon, :],
                                                        up_flux_clear[lon, :],
                                                        down_flux_clear[lon, :],
                                                        heating_rate_clear[lon, :],
                                                        0,
                                                        ClIWP[lon, :],
                                                        ClLWP[lon, :],
                                                        ClIceSize,
                                                        ClDropSize)

        dims_mid = combine_dimensions([state['air_temperature']], ['x', 'y', 'z'])
        # dims_int = combine_dimensions([state['air_temperature_on_interface_levels']], ['x', 'y', 'z'])

        tendencies = {
            'air_temperature': DataArray(
                heating_rate, dims=dims_mid, attrs={'units': 'K day^-1'}).to_units('K s^-1')
        }

        diagnostics = {}

        return tendencies, diagnostics
