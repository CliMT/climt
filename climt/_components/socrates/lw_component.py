import os
import numpy as np
from sympl import (
    get_constant, set_constant, initialize_numpy_arrays_with_properties,
    TendencyComponent
)
from ..._core import get_interface_values
from . import _runesfile as rf
from .utils import (
    resolve_spectral_file,
    get_num_bands,
    gas_name_mapping,
    aerosol_name_mapping,
    base_input_properties
)

class SocratesLongwave(TendencyComponent):
    """
    Socrates longwave radiation model wrapper.
    """
    input_properties = base_input_properties.copy()

    tendency_properties = {
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK s^-1',
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
        'air_temperature_tendency_from_longwave': {
            'dims': ['mid_levels', '*'],
            'units': 'degK s^-1',
        },
    }

    def __init__(
        self,
        spectral_file='sp_lw_ga9',
        n_band=None,
        albedo_lw=0.0,
        albedo_sw=0.5,
        **kwargs):
        
        self.spectral_file = spectral_file
        
        if n_band is None:
            resolved_path = resolve_spectral_file(spectral_file)
            n_band = get_num_bands(resolved_path)
            if n_band is None:
                n_band = 9
                
        self.n_band = n_band
        self.albedo_lw = albedo_lw
        self.albedo_sw = albedo_sw
        self.spectrum_setter = None

        self.input_properties = SocratesLongwave.input_properties.copy()

        try:
            set_constant('zenith_angle', np.pi*0.25, 'dimensionless')
        except ValueError:
            pass # Already set

        self._g = np.array(get_constant('gravitational_acceleration', 'm/s^2'))
        self._R = np.array(get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1'))
        self._secs_per_day = np.array(get_constant('seconds_per_day', 'dimensionless'))
        self._Cpd = np.array(get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K'))
        self._solar_constant = np.array(get_constant('solar_constant', 'W m^-2'))
        self._zenith_angle = np.array(get_constant('zenith_angle', 'dimensionless'))

        super(SocratesLongwave, self).__init__(**kwargs)

    def array_call(self, state):
        runesobject = rf._runesfile()

        temp_np = np.asarray(state['air_temperature'])
        n_layer = temp_np.shape[0]
        n_profile = temp_np.size // n_layer
        n_band = self.n_band

        za = state.get('zenith_angle')
        if za is None:
            za = self._zenith_angle

        flux_ground = np.zeros((n_profile, n_band))

        grey_albedo_sw = np.array(self.albedo_sw)
        grey_albedo_lw = np.array(self.albedo_lw)

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, state, self.input_properties
        )

        if self.spectrum_setter is None:
            resolved_path = resolve_spectral_file(self.spectral_file)
            runesobject._input_spectrum_lw(resolved_path)
            self.spectrum_setter = 1

        runesobject._set_planetary_parameters(
            self._g, 
            self._secs_per_day, 
            self._R, 
            self._Cpd,
            np.asarray(za), 
            self._solar_constant,
            grey_albedo_lw,
            grey_albedo_sw
        )

        T_interface = get_interface_values(
            state['air_temperature'],
            state['surface_temperature'],
            state['air_pressure'],
            state['air_pressure_on_interface_levels']
        )

        state['air_temperature_on_interface_levels'][:] = T_interface

        cf_shape = np.shape(state['air_temperature'])
        cld_frac = np.asarray(state.get('cloud_area_fraction_in_atmosphere_layer', np.zeros(cf_shape)))
        liquid_water_path = np.asarray(state.get('mass_content_of_cloud_liquid_water_in_atmosphere_layer', np.zeros(cf_shape)))
        ice_water_path = np.asarray(state.get('mass_content_of_cloud_ice_in_atmosphere_layer', np.zeros(cf_shape)))
        radius_water_droplet = np.asarray(state.get('cloud_water_droplet_radius', np.zeros(cf_shape)))
        radius_ice_particle = np.asarray(state.get('cloud_ice_particle_size', np.zeros(cf_shape)))

        # Gather gases in correct positional order
        gases = [
            'h2o', 'co2', 'o3', 'n2o', 'co', 'ch4', 'o2', 'no', 'so2', 'no2', 'nh3', 'hno3', 'n2',
            'cfc11', 'cfc12', 'cfc113', 'hcfc22', 'hfc125', 'hfc134a', 'cfc114', 'tio', 'vo', 'h2',
            'he', 'ocs', 'na', 'k', 'feh', 'crh', 'li', 'rb', 'cs', 'ph3', 'c2h2', 'hcn', 'h2s', 'ar', 'o', 'n', 'no3',
            'n2o5', 'hono', 'ho2no2', 'h2o2', 'c2h6', 'ch3', 'h2co', 'ho2',
            'hdo', 'hcl', 'hf', 'cosso', 'tosso', 'yosos',
            'ch3cho', 'ch3ooh', 'ch3coch3', 'ch3cocho', 'chocho', 'c2h5cho',
            'hoch2cho', 'c2h5coch3', 'mvk', 'macr', 'pan', 'ch3ono2'
        ]
        gas_arrays = []
        for gas in gases:
            gas_arrays.append(np.asarray(state.get(gas_name_mapping[gas], np.zeros(cf_shape))))

        # Gather aerosols in correct positional order
        aerosols = [
            'water_soluble', 'dust_like', 'oceanic', 'soot', 'ash', 'sulphuric',
            'ammonium_sulphate', 'saharan_dust',
            'accum_sulphate', 'aitken_sulphate',
            'fresh_soot', 'aged_soot',
            'sodium_chloride', 'seasalt_film', 'seasalt_jet',
            'dust_div1', 'dust_div2', 'dust_div3',
            'dust_div4', 'dust_div5', 'dust_div6',
            'biomass_1', 'biomass_2',
            'biogenic',
            'ocff_fresh', 'ocff_aged',
            'delta', 'murk',
            'nitrate',
            'twobindust_1', 'twobindust_2'
        ]
        aerosol_arrays = []
        for aero in aerosols:
            aerosol_arrays.append(np.asarray(state.get(aerosol_name_mapping[aero], np.zeros(cf_shape))))

        runesobject._set_fields_lw(
            n_profile, n_layer, n_band,
            *(gas_arrays + aerosol_arrays),
            np.asarray(state['air_pressure']), 
            np.asarray(state['air_temperature']), 
            np.asarray(state['air_pressure_on_interface_levels']), 
            np.asarray(state['air_temperature_on_interface_levels']),
            np.asarray(state['surface_temperature']),
            np.asarray(flux_ground), 
            np.asarray(diagnostics['downwelling_longwave_flux_in_air']), 
            np.asarray(diagnostics['upwelling_longwave_flux_in_air']), 
            np.asarray(tendencies['air_temperature']),
            cld_frac,
            liquid_water_path,
            ice_water_path,
            radius_water_droplet,
            radius_ice_particle,
            True 
        )        

        diagnostics['air_temperature_tendency_from_longwave'] = tendencies['air_temperature']

        return tendencies, diagnostics
