import numpy as np
from sympl import (
    get_constant, initialize_numpy_arrays_with_properties,
    TendencyComponent
)
from ...._core import (
    mass_to_volume_mixing_ratio, get_interface_values, ensure_contiguous_state)
from numpy import pi as PI
from ..rrtmg_common import (
    rrtmg_cloud_overlap_method_dict, rrtmg_cloud_props_dict,
    rrtmg_cloud_ice_props_dict, rrtmg_cloud_liquid_props_dict,
    rrtmg_random_number_dict)
import logging
try:
    from . import _rrtmg_lw
except ImportError as error:
    logging.warning(
        'Import failed. RRTMG Longwave is likely not compiled and '
        'will not be available.'
    )
    print(error)


class RRTMGLongwave(TendencyComponent):
    """
    The Rapid Radiative Transfer Model (RRTMG).

    This module wraps RRTMG for longwave radiation
    (i.e, emission from the earth's surface).
    """

    num_longwave_bands = 16
    num_reduced_g_intervals = 140
    rrtm_iplon = 1

    input_properties = {
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'mbar',
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'mbar',
        },
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'g/g',
        },
        'mole_fraction_of_ozone_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_carbon_dioxide_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_methane_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_nitrous_oxide_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_oxygen_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_cfc11_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_cfc12_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_cfc22_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'mole_fraction_of_carbon_tetrachloride_in_air': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'surface_longwave_emissivity': {
            'dims': ['num_longwave_bands', '*'],
            'units': 'dimensionless',
        },
        'cloud_area_fraction_in_atmosphere_layer': {
            'dims': ['mid_levels', '*'],
            'units': 'dimensionless',
        },
        'longwave_optical_thickness_due_to_cloud': {
            'dims': ['mid_levels', '*', 'num_longwave_bands'],
            'units': 'dimensionless',
        },
        'mass_content_of_cloud_ice_in_atmosphere_layer': {
            'dims': ['mid_levels', '*'],
            'units': 'g m^-2',
        },
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {
            'dims': ['mid_levels', '*'],
            'units': 'g m^-2',
        },
        'cloud_ice_particle_size': {
            'dims': ['mid_levels', '*'],
            'units': 'micrometer',
        },
        'cloud_water_droplet_radius': {
            'dims': ['mid_levels', '*'],
            'units': 'micrometer',
        },
        'longwave_optical_thickness_due_to_aerosol': {
            'dims': ['num_longwave_bands', 'mid_levels', '*'],
            'units': 'dimensionless',
        },
    }

    tendency_properties = {
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK day^-1',
        },
    }

    diagnostic_properties = {
        'upwelling_longwave_flux_in_air': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'downwelling_longwave_flux_in_air': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'upwelling_longwave_flux_in_air_assuming_clear_sky': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'downwelling_longwave_flux_in_air_assuming_clear_sky': {
            'dims': ['interface_levels', '*'],
            'units': 'W m^-2',
        },
        'air_temperature_tendency_from_longwave_assuming_clear_sky': {
            'dims': ['mid_levels', '*'],
            'units': 'degK day^-1',
        },
        'air_temperature_tendency_from_longwave': {
            'dims': ['mid_levels', '*'],
            'units': 'degK day^-1',
        },
    }

    def __init__(
            self,
            calculate_change_up_flux=False,
            cloud_overlap_method=None,
            cloud_optical_properties='liquid_and_ice_clouds',
            cloud_ice_properties='ebert_curry_two',
            cloud_liquid_water_properties='radius_dependent_absorption',
            calculate_interface_temperature=True,
            mcica=False,
            random_number_generator='mersenne_twister',
            **kwargs):
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

                * :code:`direct_input` = Both cloud fraction and cloud optical depth are input directly.
                  Other cloud properties (ie cloud particle size) are irrelevant.
                * :code:`single_cloud_type` = Cloud fraction and cloud physical properties are input, ice
                  and liquid clouds are treated together, cloud absorptivity is a constant value (0.060241).
                  Not available with McICA.
                * :code:`liquid_and_ice_clouds` = Cloud fraction and cloud physical properties are input, ice and liquid
                  clouds are treated separately. Cloud optical depth is calculated from the cloud ice and water particle
                  sizes and the mass content of cloud ice and cloud water.


            cloud_ice_properties (string):
                set bounds on ice particle size. This is not used if 'cloud_optical_properties' == 'direct_input'

                * :code:`ebert_curry_one` = ice particle has effective radius >= 10.0 micron `[Ebert and Curry 1992]`_
                * :code:`ebert_curry_two` = ice particle has effective radius between 13.0 and 130.0 micron `[Ebert and Curry 1992]`_
                * :code:`key_streamer_manual` = ice particle has effective radius between 5.0 and 131.0 micron
                  `[Key, Streamer Ref. Manual, 1996]`_
                * :code:`fu` = ice particle has generalised effective size (dge) between 5.0 and 140.0 micron
                  `[Fu, 1996]`_. (dge = 1.0315 * r_ec)

                Default value is 0.

            cloud_liquid_water_properties (string):
                set treatment of cloud liquid water. This is not used if 'cloud_optical_properties' == 'direct_input'.

                * :code:`radius_independent_absorption` = use radius independent absorption coefficient
                * :code:`radius_dependent_absorption` = use radius dependent absorption coefficient (radius between 2.5 and 60 micron)

            calculate_interface_temperature (bool):
                if :code:`True`, the interface temperature is calculated internally using a weighted
                interpolation routine. If :code:`False`, the quantity called
                :code:`air_temperature_on_interface_levels` in the input state needs to be manually
                updated by user code.

            mcica (bool):
                * mcica = True: use the McICA version of the longwave component of RRTMG
                * mcica = False: use the nomcica version of the longwave component of RRTMG

            random_number_generator (string):
                Different methods of generating random numbers for McICA.
                * :code:`kissvec`
                * :code:`mersenne_twister`

        .. _[Ebert and Curry 1992]:
            http://onlinelibrary.wiley.com/doi/10.1029/91JD02472/abstract

        .. _[Key, Streamer Ref. Manual, 1996]:
            https://stratus.ssec.wisc.edu/streamer/userman.pdf

        .. _[Fu, 1996]:
             http://journals.ametsoc.org/doi/abs/10.1175/1520-0442(1996)009%3C2058%3AAAPOTS%3E2.0.CO%3B2

            """
        self.input_properties = RRTMGLongwave.input_properties.copy()
        if calculate_change_up_flux:
            self._calc_dflxdt = 1
        else:
            self._calc_dflxdt = 0

        self._mcica = mcica
        if mcica:
            self._permute_seed = None
            self._random_number_generator = rrtmg_random_number_dict[
                random_number_generator.lower()]
            if type(cloud_overlap_method) is str:
                if cloud_overlap_method.lower() == 'clear_only':
                    logging.info(
                        "cloud_overlap_method == 'clear_only'."
                        " This overrides all other properties. "
                        "There are no clouds."
                    )
            if cloud_optical_properties.lower() == 'single_cloud_type':
                logging.warning(
                    "cloud_optical_properties must be 'direct_input' or "
                    "'liquid_and_ice_clouds' for radiative calculations with "
                    "clouds using McICA."
                )

        if cloud_overlap_method is None:
            cloud_overlap_method = 'random'
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

        if not self._calc_Tint:
            self.input_properties['air_temperature_on_interface_levels'] = {
                'dims': ['interface_levels', '*'],
                'units': 'degK',
            }

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

        if not mcica:
            # TODO Add all other flags as well
            _rrtmg_lw.initialise_rrtm_radiation(
                self._Cpd,
                self._cloud_overlap,
                self._calc_dflxdt,
                self._cloud_optics,
                self._ice_props,
                self._liq_props)
        super(RRTMGLongwave, self).__init__(**kwargs)

    @ensure_contiguous_state
    def array_call(self, state):
        Q = mass_to_volume_mixing_ratio(state['specific_humidity'], 18.02)
        n_layers, n_columns = state['air_temperature'].shape

        if self._calc_Tint:
            T_interface = get_interface_values(
                state['air_temperature'],
                state['surface_temperature'],
                state['air_pressure'],
                state['air_pressure_on_interface_levels']
            )
        else:
            T_interface = state['air_temperature_on_interface_levels']

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )
        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, state, self.input_properties
        )

        if self._mcica:

            # First, define extra arrays needed for mcica.
            # The values for these arrays are calculated from state in the
            # first part of _rrtmg_sw.rrtm_calculate_longwave_fluxes_mcica.
            # Specifically they are calculated by mcica_subcol_gen_lw.f90
            # and are input to rrtmg_lw_rad.f90
            num_reduced_g_intervals = self.num_reduced_g_intervals
            mid_levels = state['air_pressure'].shape[0]

            try:
                num_cols = state['air_pressure'].shape[1]
            except IndexError:
                num_cols = 1

            mcica_properties = {
                'cloud_area_fraction_in_atmosphere_layer': np.zeros(
                    (mid_levels, num_cols, num_reduced_g_intervals)),
                'mass_content_of_cloud_ice_in_atmosphere_layer': np.zeros(
                    (mid_levels, num_cols, num_reduced_g_intervals)),
                'mass_content_of_cloud_liquid_water_in_atmosphere_layer': np.zeros(
                    (mid_levels, num_cols, num_reduced_g_intervals)),
                'cloud_ice_particle_size': np.zeros((mid_levels, num_cols)),
                'cloud_water_droplet_radius': np.zeros((mid_levels, num_cols)),
                'longwave_optical_thickness_due_to_cloud': np.zeros(
                    (mid_levels, num_cols, num_reduced_g_intervals))
            }

            # Change parameter for random number generator - each time the
            # radiation is called, with the same state / input properties,
            # a different result is obtained, because the wavelengths which
            # see cloud differ between each call.
            if self._random_number_generator == 0:
                # KISS algorithm: The seed determines the number of times
                # the random number generator is called iteratively to create a
                # new random number. The value range of the seed is limited to
                # avoid a performance decrease.
                self._permute_seed = np.random.randint(0, 1024)
            elif self._random_number_generator == 1:
                # Mersenne Twister: Use random seed from the full 32bit range.
                self._permute_seed = np.random.randint(0, 2**31 - 1)

            _rrtmg_lw.initialise_rrtm_radiation_mcica(
                self._Cpd,
                self._cloud_overlap,
                self._calc_dflxdt,
                self._cloud_optics,
                self._ice_props,
                self._liq_props,
                self._permute_seed,
                self._random_number_generator)

            _rrtmg_lw.rrtm_calculate_longwave_fluxes_mcica(
                self.rrtm_iplon,
                n_columns,
                n_layers,
                state['air_pressure'],
                state['air_pressure_on_interface_levels'],
                state['air_temperature'],
                T_interface,
                state['surface_temperature'],
                Q,
                state['mole_fraction_of_ozone_in_air'],
                state['mole_fraction_of_carbon_dioxide_in_air'],
                state['mole_fraction_of_methane_in_air'],
                state['mole_fraction_of_nitrous_oxide_in_air'],
                state['mole_fraction_of_oxygen_in_air'],
                state['mole_fraction_of_cfc11_in_air'],
                state['mole_fraction_of_cfc12_in_air'],
                state['mole_fraction_of_cfc22_in_air'],
                state['mole_fraction_of_carbon_tetrachloride_in_air'],
                state['surface_longwave_emissivity'],
                state['cloud_area_fraction_in_atmosphere_layer'],
                state['longwave_optical_thickness_due_to_aerosol'],
                diagnostics['upwelling_longwave_flux_in_air'],
                diagnostics['downwelling_longwave_flux_in_air'],
                tendencies['air_temperature'],
                diagnostics['upwelling_longwave_flux_in_air_assuming_clear_sky'],
                diagnostics['downwelling_longwave_flux_in_air_assuming_clear_sky'],
                diagnostics['air_temperature_tendency_from_longwave_assuming_clear_sky'],
                state['longwave_optical_thickness_due_to_cloud'],
                state['mass_content_of_cloud_ice_in_atmosphere_layer'],
                state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'],
                state['cloud_ice_particle_size'],
                state['cloud_water_droplet_radius'],
                mcica_properties['cloud_area_fraction_in_atmosphere_layer'],
                mcica_properties['longwave_optical_thickness_due_to_cloud'],
                mcica_properties['mass_content_of_cloud_ice_in_atmosphere_layer'],
                mcica_properties['mass_content_of_cloud_liquid_water_in_atmosphere_layer'],
                mcica_properties['cloud_ice_particle_size'],
                mcica_properties['cloud_water_droplet_radius']
            )

        else:
            _rrtmg_lw.rrtm_calculate_longwave_fluxes(
                n_columns,
                n_layers,
                state['air_pressure'],
                state['air_pressure_on_interface_levels'],
                state['air_temperature'],
                T_interface,
                state['surface_temperature'],
                Q,
                state['mole_fraction_of_ozone_in_air'],
                state['mole_fraction_of_carbon_dioxide_in_air'],
                state['mole_fraction_of_methane_in_air'],
                state['mole_fraction_of_nitrous_oxide_in_air'],
                state['mole_fraction_of_oxygen_in_air'],
                state['mole_fraction_of_cfc11_in_air'],
                state['mole_fraction_of_cfc12_in_air'],
                state['mole_fraction_of_cfc22_in_air'],
                state['mole_fraction_of_carbon_tetrachloride_in_air'],
                state['surface_longwave_emissivity'],
                state['cloud_area_fraction_in_atmosphere_layer'],
                state['longwave_optical_thickness_due_to_aerosol'],
                diagnostics['upwelling_longwave_flux_in_air'],
                diagnostics['downwelling_longwave_flux_in_air'],
                tendencies['air_temperature'],
                diagnostics['upwelling_longwave_flux_in_air_assuming_clear_sky'],
                diagnostics['downwelling_longwave_flux_in_air_assuming_clear_sky'],
                diagnostics['air_temperature_tendency_from_longwave_assuming_clear_sky'],
                state['longwave_optical_thickness_due_to_cloud'],
                state['mass_content_of_cloud_ice_in_atmosphere_layer'],
                state['mass_content_of_cloud_liquid_water_in_atmosphere_layer'],
                state['cloud_ice_particle_size'],
                state['cloud_water_droplet_radius'])

        diagnostics['air_temperature_tendency_from_longwave'] = tendencies['air_temperature']

        return tendencies, diagnostics
