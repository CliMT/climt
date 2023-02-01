from sympl import Stepper
from sympl import initialize_numpy_arrays_with_properties
import numpy as np


class BucketHydrology(Stepper):
    """
    Manages surface energy and moisture balance
    This component assumes that the surface is a slab with some heat capacity and moisture holding capacity.
    Calculates the sensible and latent heat flux, takes precipitation values as input.
    """

    input_properties = {
        'downwelling_longwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'downwelling_shortwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'upwelling_longwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'upwelling_shortwave_flux_in_air': {
            'dims': ['*', 'interface_levels'],
            'units': 'W m^-2',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'surface_material_density': {
            'dims': ['*'],
            'units': 'kg m^-3',
        },
        'soil_layer_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'heat_capacity_of_soil': {
            'dims': ['*'],
            'units': 'J kg^-1 degK^-1',
        },
        'lwe_thickness_of_soil_moisture_content': {
            'dims': ['*'],
            'units': 'm',
        },
        'convective_precipitation_rate': {
            'dims': ['*'],
            'units': 'm s^-1 ',
        },
        'stratiform_precipitation_rate': {
            'dims': ['*'],
            'units': 'm s^-1 ',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'kg/kg ',
        },
        'surface_specific_humidity': {
            'dims': ['*'],
            'units': 'kg/kg ',
        },
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK ',
        },
        'northward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'eastward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'area_type': {
            'dims': ['*'],
            'units': 'dimensionless',
        },
    }

    diagnostic_properties = {
        'precipitation_rate': {
            'dims': ['*'],
            'units': 'm s^-1 ',
        },
        'surface_upward_latent_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_upward_sensible_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'evaporation_rate': {
            'dims': ['*'],
            'units': 'm s^-1',
        },
    }

    output_properties = {
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'lwe_thickness_of_soil_moisture_content': {
            'dims': ['*'],
            'units': 'm'},
    }

    def __init__(self, soil_moisture_max=0.15, beta_parameter=0.75,
                 specific_latent_heat_of_water=2260000, bulk_coefficient=0.0011, **kwargs):
        """
        Args:
        soil_moisture_max:
            The maximum moisture that can be held by the surface_temperature
        beta_parameter:
            A constant value that is used in the beta_factor calculation.
        bulk_coefficient:
            The bulk transfer coefficient that is used to calculate
            maximum evaporation rate and sensible heat flux
        """
        self._smax = soil_moisture_max
        self._g = beta_parameter
        self._c = bulk_coefficient
        self._l = specific_latent_heat_of_water
        super(BucketHydrology, self).__init__(**kwargs)

    def array_call(self, state, timestep):
        '''
        Calculates sensible and latent heat flux and returns
        surface temperature and soil moisture after timestep.
        '''

        beta_factor = 0

        new_state = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties
        )

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        north_wind_speed = state['northward_wind'][0]
        east_wind_speed = state['eastward_wind'][0]

        wind_speed = np.sqrt(np.power(north_wind_speed, 2) +
                             np.power(east_wind_speed, 2))

        potential_evaporation = self._c * wind_speed * (state['surface_specific_humidity'] - state['specific_humidity'][0])

        precipitation_rate = state['convective_precipitation_rate'] + state['stratiform_precipitation_rate']
        diagnostics['precipitation_rate'] = precipitation_rate

        soil_moisture = state['lwe_thickness_of_soil_moisture_content']

        beta_factor = np.ones(soil_moisture.shape)

        mask = soil_moisture <= self._g * self._smax
        beta_factor[mask] = soil_moisture[mask]/(self._g*self._smax)

        evaporation_rate = beta_factor * potential_evaporation
        diagnostics['evaporation_rate'] = evaporation_rate

        soil_moisture_tendency = np.zeros(soil_moisture.shape)
        mask = np.logical_or(soil_moisture < self._smax, precipitation_rate <= evaporation_rate)
        soil_moisture_tendency[mask] = precipitation_rate[mask] - evaporation_rate[mask]

        surface_upward_latent_heat_flux = self._l * evaporation_rate
        surface_upward_sensible_heat_flux = self._c * wind_speed * (state['surface_temperature'] - state['air_temperature'][0])
        diagnostics['surface_upward_sensible_heat_flux'] = surface_upward_sensible_heat_flux
        diagnostics['surface_upward_latent_heat_flux'] = surface_upward_latent_heat_flux

        net_heat_flux = (
            state['downwelling_shortwave_flux_in_air'][:, 0] +
            state['downwelling_longwave_flux_in_air'][:, 0] -
            state['upwelling_shortwave_flux_in_air'][:, 0] -
            state['upwelling_longwave_flux_in_air'][:, 0] -
            surface_upward_sensible_heat_flux -
            surface_upward_latent_heat_flux)

        mass_surface_slab = state['surface_material_density'] * \
            state['soil_layer_thickness']
        heat_capacity_surface = mass_surface_slab * state['heat_capacity_of_soil']

        new_state['surface_temperature'] = state['surface_temperature'] + \
            (net_heat_flux/heat_capacity_surface * timestep.total_seconds())
        new_state['lwe_thickness_of_soil_moisture_content'] = np.minimum(state['lwe_thickness_of_soil_moisture_content'] +
                                                                         (soil_moisture_tendency * timestep.total_seconds()), 0.15)

        return diagnostics, new_state
