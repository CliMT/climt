from math import sqrt
from sympl import TendencyComponent, initialize_numpy_arrays_with_properties


class BucketSurface(TendencyComponent):
    """
    Calculates the surface energy and the hydrology balance.

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
        'lwe_thickness_of_soil_moisture_content':{
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
    }

    tendency_properties = {
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK s^-1',
        },
        'lwe_thickness_of_soil_moisture_content':{
            'dims': ['*'],
            'units': 'm s^-1',
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
        'beta_factor': {
            'dims': ['*'],
            'units': 'dimensionless',
        },
    }


    def __init__(self, soil_moisture_max=0.15, g=0.75, specific_latent_heat_of_water=2260000,
                 bulk_coefficient=0.0011, **kwargs):

        self._smax = soil_moisture_max
        self._g = g
        self._c = bulk_coefficient
        self._l = specific_latent_heat_of_water
        super(BucketSurface, self).__init__(**kwargs)


    def array_call(self, raw_state):
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )
        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, raw_state, self.input_properties
        )

        diagnostics['beta_factor'] = 0

        wind_speed = sqrt(pow(raw_state['northward_wind'][0], 2) + \
           pow(raw_state['eastward_wind'][0], 2))
        evaporation_rate_max = self._c * wind_speed * \
           (raw_state['surface_specific_humidity'] - raw_state['specific_humidity'][0])


        precipitation_rate = raw_state['convective_precipitation_rate'] + \
                             raw_state['stratiform_precipitation_rate']


        soil_moisture = raw_state['lwe_thickness_of_soil_moisture_content']

        soil_moisture_tendency = 0

        if soil_moisture >= self._g * self._smax:
            diagnostics['beta_factor'] = 1
        else:
            diagnostics['beta_factor'] = soil_moisture/(self._g*self._smax)

        evaporation_rate = diagnostics['beta_factor'] * evaporation_rate_max

        if soil_moisture < self._smax or precipitation_rate <= evaporation_rate:
            soil_moisture_tendency = precipitation_rate - evaporation_rate
        else:
            soil_moisture_tendency = 0


        diagnostics['surface_upward_latent_heat_flux'] = self._l * evaporation_rate
        diagnostics['surface_upward_sensible_heat_flux'] = self._c * wind_speed * \
        (raw_state['surface_temperature'] - raw_state['air_temperature'][0])

        net_heat_flux = (
            raw_state['downwelling_shortwave_flux_in_air'][:, 0] +
            raw_state['downwelling_longwave_flux_in_air'][:, 0] -
            raw_state['upwelling_shortwave_flux_in_air'][:, 0] -
            raw_state['upwelling_longwave_flux_in_air'][:, 0] -
            diagnostics['surface_upward_sensible_heat_flux'] -
            diagnostics['surface_upward_latent_heat_flux']
        )


        mass_surface_slab = raw_state['surface_material_density'] * \
            raw_state['soil_layer_thickness']
        heat_capacity_surface = mass_surface_slab * raw_state['heat_capacity_of_soil']


        tendencies = {'surface_temperature': net_heat_flux/heat_capacity_surface,
                      'lwe_thickness_of_soil_moisture_content' : soil_moisture_tendency}


        return tendencies, diagnostics
