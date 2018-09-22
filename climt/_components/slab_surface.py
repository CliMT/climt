from sympl import TendencyComponent, initialize_numpy_arrays_with_properties
import numpy as np


class SlabSurface(TendencyComponent):
    """
    Calculate the surface energy balance.

    This component assumes the surface is a slab of possibly
    varying heat capacity, and calculates the surface temperature.
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
        'surface_upward_latent_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK',
        },
        'surface_upward_sensible_heat_flux': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'surface_thermal_capacity': {
            'dims': ['*'],
            'units': 'J kg^-1 degK^-1',
        },
        'surface_material_density': {
            'dims': ['*'],
            'units': 'kg m^-3',
        },
        'upward_heat_flux_at_ground_level_in_soil': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'heat_flux_into_sea_water_due_to_sea_ice': {
            'dims': ['*'],
            'units': 'W m^-2',
        },
        'area_type': {
            'dims': ['*'],
            'units': 'dimensionless',
        },
        'soil_layer_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'ocean_mixed_layer_thickness': {
            'dims': ['*'],
            'units': 'm',
        },
        'heat_capacity_of_soil': {
            'dims': ['*'],
            'units': 'J kg^-1 degK^-1',
        },
        'sea_water_density': {
            'dims': ['*'],
            'units': 'kg m^-3',
        },
    }

    tendency_properties = {
        'surface_temperature': {
            'dims': ['*'],
            'units': 'degK s^-1',
        }
    }

    diagnostic_properties = {
        'depth_of_slab_surface': {
            'dims': ['*'],
            'units': 'm',
        },
    }

    def array_call(self, raw_state):
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )
        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, raw_state, self.input_properties
        )

        net_heat_flux = (
            raw_state['downwelling_shortwave_flux_in_air'][:, 0] +
            raw_state['downwelling_longwave_flux_in_air'][:, 0] -
            raw_state['upwelling_shortwave_flux_in_air'][:, 0] -
            raw_state['upwelling_longwave_flux_in_air'][:, 0] -
            raw_state['surface_upward_sensible_heat_flux'] -
            raw_state['surface_upward_latent_heat_flux']
        )

        area_type = raw_state['area_type'].astype(str)

        land_mask = np.logical_or(area_type == 'land', area_type == 'land_ice')
        sea_mask = np.logical_or(area_type == 'sea', area_type == 'sea_ice')
        land_ice_mask = area_type == 'land_ice'
        sea_ice_mask = area_type == 'sea_ice'

        net_heat_flux[land_ice_mask] = -raw_state['upward_heat_flux_at_ground_level_in_soil'][land_ice_mask]
        net_heat_flux[sea_ice_mask] = raw_state['heat_flux_into_sea_water_due_to_sea_ice'][sea_ice_mask]
        raw_state['surface_material_density'][sea_mask] = raw_state['sea_water_density'][sea_mask]
        raw_state['surface_thermal_capacity'][land_mask] = raw_state['heat_capacity_of_soil'][land_mask]
        diagnostics['depth_of_slab_surface'][sea_mask] =\
            raw_state['ocean_mixed_layer_thickness'][sea_mask]
        diagnostics['depth_of_slab_surface'][land_mask] =\
            raw_state['soil_layer_thickness'][land_mask]

        mass_surface_slab = raw_state['surface_material_density'] * \
            diagnostics['depth_of_slab_surface']
        heat_capacity_surface = mass_surface_slab * raw_state['surface_thermal_capacity']

        tendencies['surface_temperature'][:] = net_heat_flux/heat_capacity_surface
        tendencies['surface_temperature'][land_ice_mask] = 0
        tendencies['surface_temperature'][sea_ice_mask] = 0

        return tendencies, diagnostics
