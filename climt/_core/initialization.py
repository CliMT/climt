from sympl import (DataArray, add_direction_names,
                   ensure_no_shared_keys)
import numpy as np
import copy
from datetime import datetime
from scipy.interpolate import CubicSpline


def init_mid_level_pressures(array_dims, quantity_description, initial_state):
    # We assume the vertical level is the last dimension
    vert_levels = array_dims[-1]

    p_surf = quantity_description['surface_air_pressure']['default_value']
    spacing = np.linspace(0.998, 0.001, vert_levels)
    single_column = p_surf*spacing
    single_column = single_column[np.newaxis, np.newaxis, :]

    return single_column*np.ones(array_dims, order='F')


def init_interface_level_pressures(array_dims, quantity_description, initial_state):

    vert_levels = array_dims[-1]
    p_surf = quantity_description['surface_air_pressure']['default_value']
    spacing = np.linspace(0.998, 0.001, vert_levels-1)
    midlevel = p_surf*spacing

    interface = np.zeros(vert_levels)
    interface[1:-1] = 0.5*(midlevel[:-1] + midlevel[1:])
    interface[0] = p_surf
    interface[-1] = 0.0005*p_surf

    return interface*np.ones(array_dims, order='F')


def init_interface_level_sigma(array_dims, quantity_description, initial_state):

    vert_levels = array_dims[-1]
    spacing = np.linspace(0.998, 0.001, vert_levels-1)
    midlevel = spacing

    interface = np.zeros(vert_levels)
    interface[1:-1] = 0.5*(midlevel[:-1] + midlevel[1:])
    interface[0] = 1.
    interface[-1] = 0.0005

    return interface*np.ones(array_dims, order='F')


def init_mid_level_sigma(array_dims, quantity_description, initial_state):

    vert_levels = array_dims[-1]
    spacing = np.linspace(0.998, 0.001, vert_levels)
    midlevel = spacing

    return midlevel*np.ones(array_dims, order='F')


def init_interface_level_tau_longwave(array_dims, quantity_description, initial_state):

    vert_levels = array_dims[-1]
    spacing = np.linspace(0.998, 0.001, vert_levels-1)
    midlevel = spacing

    interface = np.zeros(vert_levels)
    interface[1:-1] = 0.5*(midlevel[:-1] + midlevel[1:])
    interface[0] = 1.
    interface[-1] = 0.0005
    tau_longwave = 1.*(1 - interface)

    return tau_longwave*np.ones(array_dims, order='F')


def init_ozone(array_dims, quantity_description, initial_state):

    import pkg_resources

    init_array = np.ones(array_dims, order='F')
    current_levels = np.linspace(0.998, 0.001, 30)[::-1]

    target_levels = np.linspace(0.998, 0.001, array_dims[-1])[::-1]

    file_name = 'ozone_profile.npy'
    file_path = 'climt._data'

    resource_path = pkg_resources.resource_filename(file_path, file_name)

    profile = np.load(resource_path)

    target_profile = CubicSpline(current_levels, profile[::-1])(target_levels)[::-1]

    if array_dims[-1] == 30:
        target_profile = profile

    # Ensure ozone concentration at top of model is not excessive
    # This is more in line with observations.
    target_profile[-1] /= 10

    init_array[:] = target_profile[np.newaxis, np.newaxis, :]

    return init_array


climt_quantity_descriptions = {
    'air_pressure': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'Pa',
        'init_func': init_mid_level_pressures
    },
    'air_pressure_on_interface_levels': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'Pa',
        'init_func': init_interface_level_pressures
    },
    'sigma_levels': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_func': init_mid_level_sigma
    },
    'sigma_on_interface_levels': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'dimensionless',
        'init_func': init_interface_level_sigma
    },
    'longwave_optical_depth_on_interface_levels': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'dimensionless',
        'init_func': init_interface_level_tau_longwave
    },
    'surface_air_pressure': {
        'dims': ['x', 'y'],
        'units': 'Pa',
        'default_value': 1.0132e5
    },
    'air_temperature': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'degK',
        'default_value': 290.
    },
    'air_temperature_on_interface_levels': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'degK',
        'default_value': 290.
    },
    'surface_temperature': {
        'dims': ['x', 'y'],
        'units': 'degK',
        'default_value': 300.
    },
    'sea_surface_temperature': {
        'dims': ['x', 'y'],
        'units': 'degK',
        'default_value': 300.
    },
    'soil_surface_temperature': {
        'dims': ['x', 'y'],
        'units': 'degK',
        'default_value': 300.
    },
    'northward_wind': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'm s^-1',
        'default_value': 0.
    },
    'eastward_wind': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'm s^-1',
        'default_value': 0.
    },
    'divergence_of_wind': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 's^-1',
        'default_value': 0.
    },
    'atmosphere_relative_vorticity': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 's^-1',
        'default_value': 0.
    },
    'surface_geopotential': {
        'dims': ['x', 'y'],
        'units': 'm^2 s^-2',
        'default_value': 0.
    },
    'surface_longwave_emissivity': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 1.
    },
    'specific_humidity': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'g/g',
        'default_value': 0.
    },
    'surface_specific_humidity': {
        'dims': ['x', 'y'],
        'units': 'g/g',
        'default_value': 0.
    },
    'mole_fraction_of_ozone_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_func': init_ozone
    },
    'mole_fraction_of_carbon_dioxide_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 330e-6
    },
    'mole_fraction_of_methane_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.
    },
    'mole_fraction_of_nitrous_oxide_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.
    },
    'mole_fraction_of_oxygen_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.21
    },
    'mole_fraction_of_nitrogen_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.78
    },
    'mole_fraction_of_hydrogen_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 500e-9
    },
    'mole_fraction_of_cfc11_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.
    },
    'mole_fraction_of_cfc12_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.
    },
    'mole_fraction_of_cfc22_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.
    },
    'mole_fraction_of_carbon_tetrachloride_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'default_value': 0.
    },
    'cloud_area_fraction_in_atmosphere_layer': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'default_value': 0.
    },
    'shortwave_optical_thickness_due_to_aerosol': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'default_value': 0.
    },
    'longwave_optical_thickness_due_to_aerosol': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'default_value': 0.
    },
    'mass_content_of_cloud_ice_in_atmosphere_layer': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'kg m^-2',
        'default_value': 0.
    },
    'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'kg m^-2',
        'default_value': 0.
    },
    'cloud_ice_particle_size': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'micrometer',
        'default_value': 20.
    },
    'cloud_water_droplet_radius': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'micrometer',
        'default_value': 10.
    },
    'longwave_optical_thickness_due_to_cloud': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'default_value': 0.
    },
    'shortwave_optical_thickness_due_to_cloud': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'default_value': 0.
    },
    'longwave_heating_rate': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'default_value': 0.
    },
    'longwave_heating_rate_assuming_clear_sky': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'default_value': 0.
    },
    'upwelling_longwave_flux_in_air': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'upwelling_longwave_flux_in_air_assuming_clear_sky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'downwelling_longwave_flux_in_air': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'downwelling_longwave_flux_in_air_assuming_clear_sky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'shortwave_heating_rate': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'default_value': 0.
    },
    'shortwave_heating_rate_assuming_clear_sky': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'default_value': 0.
    },
    'upwelling_shortwave_flux_in_air': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'upwelling_shortwave_flux_in_air_assuming_clear_sky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'downwelling_shortwave_flux_in_air': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'downwelling_shortwave_flux_in_air_assuming_clear_sky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'surface_upward_sensible_heat_flux': {
        'dims': ['x', 'y'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'surface_upward_latent_heat_flux': {
        'dims': ['x', 'y'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'upward_heat_flux_at_ground_level_in_soil': {
        'dims': ['x', 'y'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'heat_flux_into_sea_water_due_to_sea_ice': {
        'dims': ['x', 'y'],
        'units': 'W m^-2',
        'default_value': 0.
    },
    'precipitation_amount': {
        'dims': ['x', 'y'],
        'units': 'kg m^-2',
        'default_value': 0.
    },
    'convective_precipitation_amount': {
        'dims': ['x', 'y'],
        'units': 'kg m^-2',
        'default_value': 0.
    },
    'stratiform_precipitation_amount': {
        'dims': ['x', 'y'],
        'units': 'kg m^-2',
        'default_value': 0.
    },
    'precipitation_rate': {
        'dims': ['x', 'y'],
        'units': 'm s^-1',
        'default_value': 0.
    },
    'convective_precipitation_rate': {
        'dims': ['x', 'y'],
        'units': 'm s^-1',
        'default_value': 0.
    },
    'stratiform_precipitation_rate': {
        'dims': ['x', 'y'],
        'units': 'm s^-1',
        'default_value': 0.
    },
    'atmosphere_convective_mass_flux': {
        'dims': ['x', 'y'],
        'units': 'kg m^-2 s^-1',
        'default_value': 0.
    },
    'atmosphere_convective_available_potential_energy': {
        'dims': ['x', 'y'],
        'units': 'J kg^-1',
        'default_value': 0.
    },
    'convective_heating_rate': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'default_value': 0.
    },
    'zenith_angle': {
        'dims': ['x', 'y'],
        'units': 'radians',
        'default_value': 0.
    },
    'land_ice_thickness': {
        'dims': ['x', 'y'],
        'units': 'm',
        'default_value': 0.
    },
    'sea_ice_thickness': {
        'dims': ['x', 'y'],
        'units': 'm',
        'default_value': 0.
    },
    'surface_snow_thickness': {
        'dims': ['x', 'y'],
        'units': 'm',
        'default_value': 0.
    },
    'ocean_mixed_layer_thickness': {
        'dims': ['x', 'y'],
        'units': 'm',
        'default_value': 50.
    },
    'soil_layer_thickness': {
        'dims': ['x', 'y'],
        'units': 'm',
        'default_value': 50.
    },
    'soil_thermal_capacity': {
        'dims': ['x', 'y'],
        'units': 'J kg^-1 degK^-1',
        'default_value': 2000.
    },
    'area_type': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 'sea',
        'dtype': 'a100'
    },
    'snow_and_ice_temperature_spline': {
        'dims': ['x', 'y'],
        'units': 'degK',
        'default_value': CubicSpline(
            np.linspace(0, 50, 50), 260.*np.ones(50)),
        'dtype': object
    },
    'sea_water_density': {
        'dims': ['x', 'y'],
        'units': 'kg m^-3',
        'default_value': 1.029e3
    },
    'surface_albedo_for_direct_shortwave': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.06
    },

    'surface_albedo_for_diffuse_shortwave': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.06
    },

    'surface_albedo_for_direct_near_infrared': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.06
    },

    'surface_albedo_for_diffuse_near_infrared': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.06
    },
    'soil_type': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 'clay',
        'dtype': 'a100'
    },
    'surface_roughness_length': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.0002,
    },
    'surface_drag_coefficient_for_heat_in_air': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.0012,
    },
    'surface_drag_coefficient_for_momentum_in_air': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'default_value': 0.0012,
    },
    'soil_temperature': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'degK',
        'default_value': 274.,
    },
}


def get_default_state(component_list,
                      x={}, y={}, mid_levels={}, interface_levels={},
                      initial_state={}):
    """
    Return a state dictionary required to run the model.

    The model comprises of components in :code:`component_list`. If coordinate
    values in :code:`x`, :code:`y`, :code:`mid_levels` and :code:`interface_levels`
    are not provided, a single column
    centered at 0 degrees south, 0 degrees east, with 30 vertical levels is returned.

    Args:
        component_list (iterable): The _components for which a default
            state is required, in the order that they are called.
            It is assumed diagnostic outputs are updated into
            the state and passed into the next component, and that the state
            from the previous component is passed into the next component. The
            function will attempt to determine required inputs from the series
            of _components.

        x (dict,optional): A dictionary containing keys :code:`label`, :code:`values`,
            and :code:`units`. :code:`label` refers to the name the coordinate axis will assume.
            :code:`values` refers to the array of coordinate values. :code:`units` refers to the
            units of the coordinate values.
            If :code:`x` is an empty dictionary, a single default value of 0 degrees longitude
            is used.

        y (dict,optional): A dictionary containing keys :code:`label`, :code:`values`,
            and :code:`units`. :code:`label` refers to the name the coordinate axis will assume.
            :code:`values` refers to the array of coordinate values. :code:`units` refers to the
            units of the coordinate values.
            If :code:`y` is an empty dictionary, a single default value of 0 degrees latitude
            is used.

        mid_levels (dict,optional): A dictionary containing keys :code:`label`, :code:`values`,
            and :code:`units`. :code:`label` refers to the name the coordinate axis will assume.
            :code:`values` refers to the array of coordinate values. :code:`units` refers to the
            units of the coordinate values.
            If :code:`mid_levels` is an empty dictionary, 30 levels of arbitrary units are used.

        interface_levels (dict,optional): A dictionary containing keys :code:`label`, :code:`values`,
            and :code:`units`. :code:`label` refers to the name the coordinate axis will assume.
            :code:`values` refers to the array of coordinate values. :code:`units` refers to the
            units of the coordinate values.
            If :code:`interface_levels` is an empty dictionary, 31 levels of arbitrary units are used.

       initial_state (dict, optional): A dictionary containing some quantities that will
            also be added to the final output state. Must not contain any quantities that
            the output state will overwrite.


    Returns:
        default_state (dict):
            A state dictionary containing the requested
            quantities using the provided coordinate state.

    Raises:
        ValueError:
            if any of the following conditions are satisfied:

                * if :code:`component_list` is empty
                * if the shape of :code:`x['values']` and :code:`y['values']` is not the same
                * if only one of :code:`mid_levels` or :code:`interface_levels` is specified
                * if vertical coordinates are not one dimensional
                * if length of :code:`mid_levels['values']` is not one less than length
                  of :code:`interface['values']`
    """

    if len(component_list) == 0:
        raise ValueError('Component list must contain at least one component')

    if (len(mid_levels.keys()) == 0) != (len(interface_levels.keys()) == 0):
        raise ValueError('Both mid and interface levels must be specified')

    output_state = {}

    # Create 2D coordinate arrays
    if len(x.keys()) == 0:
        x_coordinate_values = np.zeros((1,))
        x_coordinate_label = 'longitude'
        x_coordinate_units = 'degrees_east'
    else:
        x_coordinate_values = x['values']
        x_coordinate_label = x['label']
        x_coordinate_units = x['units']

    if len(y.keys()) == 0:
        y_coordinate_values = np.zeros((1,))
        y_coordinate_label = 'latitude'
        y_coordinate_units = 'degrees_north'
    else:
        y_coordinate_values = y['values']
        y_coordinate_label = y['label']
        y_coordinate_units = y['units']

    if len(mid_levels.keys()) == 0:
        mid_levels_coordinate_values = np.arange(30)
        mid_levels_coordinate_label = 'mid_levels'
        mid_levels_coordinate_units = ''
    else:
        mid_levels_coordinate_values = mid_levels['values']
        mid_levels_coordinate_label = mid_levels['label']
        mid_levels_coordinate_units = mid_levels['units']

    if len(interface_levels.keys()) == 0:
        interface_levels_coordinate_values = np.arange(31)
        interface_levels_coordinate_label = 'interface_levels'
        interface_levels_coordinate_units = ''
    else:
        interface_levels_coordinate_values = interface_levels['values']
        interface_levels_coordinate_label = interface_levels['label']
        interface_levels_coordinate_units = interface_levels['units']

    if not x_coordinate_values.ndim == y_coordinate_values.ndim:
        raise ValueError('x and y coordinates must have the same shape')

    if mid_levels_coordinate_values.ndim > 1:
        raise ValueError(
            'vertical coordinate mid_levels must be one dimensional.')

    if interface_levels_coordinate_values.ndim > 1:
        raise ValueError(
            'vertical coordinate interface_levels must be one dimensional.')

    if len(mid_levels_coordinate_values) != len(interface_levels_coordinate_values)-1:
        raise ValueError('Interface levels must have one value more than mid levels')

    use_2d_coordinate = False
    two_dim_coord_dict = {}

    if x_coordinate_values.ndim == 2:
        if not x_coordinate_values.shape == y_coordinate_values.shape:
            raise ValueError(
                'If x and y are 2d coordinates, they must have the same shape')

        two_dim_coord_dict[x_coordinate_label] = {}
        two_dim_coord_dict[x_coordinate_label]['values'] = x_coordinate_values
        two_dim_coord_dict[x_coordinate_label]['logical_dims'] = (
            'logical_x_coordinate', 'logical_y_coordinate')
        two_dim_coord_dict[x_coordinate_label]['units'] = x_coordinate_units

        two_dim_coord_dict[y_coordinate_label] = {}
        two_dim_coord_dict[y_coordinate_label]['values'] = y_coordinate_values
        two_dim_coord_dict[y_coordinate_label]['logical_dims'] = (
            'logical_x_coordinate', 'logical_y_coordinate')
        two_dim_coord_dict[y_coordinate_label]['units'] = y_coordinate_units

        output_state['x'] = DataArray(
            x_coordinate_values,
            dims=two_dim_coord_dict[x_coordinate_label]['logical_dims'],
            attrs={
                'units': two_dim_coord_dict[x_coordinate_label]['units'],
                'label': x_coordinate_label
            }
        )

        output_state['y'] = DataArray(
            y_coordinate_values,
            dims=two_dim_coord_dict[y_coordinate_label]['logical_dims'],
            attrs={
                'units': two_dim_coord_dict[y_coordinate_label]['units'],
                'label': y_coordinate_label
            }
        )

        x_coordinate_values = np.arange(
            two_dim_coord_dict[x_coordinate_label]['values'].shape[0])

        y_coordinate_values = np.arange(
            two_dim_coord_dict[x_coordinate_label]['values'].shape[1])

        x_coordinate_label = 'logical_x_coordinate'
        y_coordinate_label = 'logical_y_coordinate'

        x_coordinate_units = ''
        y_coordinate_units = ''

        use_2d_coordinate = True

    x_coordinate = DataArray(
        x_coordinate_values,
        dims=(x_coordinate_label),
        attrs={'units': x_coordinate_units, 'label': x_coordinate_label})

    output_state[x_coordinate_label] = x_coordinate
    add_direction_names(x=x_coordinate_label)

    y_coordinate = DataArray(
        y_coordinate_values,
        dims=(y_coordinate_label),
        attrs={'units': y_coordinate_units, 'label': y_coordinate_label})

    output_state[y_coordinate_label] = y_coordinate
    add_direction_names(y=y_coordinate_label)

    mid_levels_coordinate = DataArray(
        mid_levels_coordinate_values,
        dims=(mid_levels_coordinate_label,),
        attrs={'units': mid_levels_coordinate_units,
               'label': mid_levels_coordinate_label})

    output_state[mid_levels_coordinate_label] = mid_levels_coordinate
    output_state['mid_levels'] = mid_levels_coordinate
    add_direction_names(
        z=mid_levels_coordinate_label)

    interface_levels_coordinate = DataArray(
        interface_levels_coordinate_values,
        dims=(interface_levels_coordinate_label,),
        attrs={'units': interface_levels_coordinate_units,
               'label': interface_levels_coordinate_label})

    output_state[interface_levels_coordinate_label] = interface_levels_coordinate
    output_state['interface_levels'] = interface_levels_coordinate
    add_direction_names(
        z=interface_levels_coordinate_label)

    if not use_2d_coordinate:
        output_state['x'] = x_coordinate
        output_state['y'] = y_coordinate

    quantity_list = set()
    temporary_description = copy.deepcopy(climt_quantity_descriptions)
    additional_dimensions = {}
    additional_descriptions = {}

    for component in component_list:
        quantity_list = quantity_list.union(set(component.inputs))

        if hasattr(component, 'extra_dimensions'):
            ensure_no_shared_keys(additional_dimensions,
                                  component.extra_dimensions)

            for dimension in component.extra_dimensions.keys():
                if component.extra_dimensions[dimension].ndim > 1:
                    raise NotImplementedError(
                        'Two dimensional coordinates in extra_dimensions not yet supported')

                output_state[dimension] = DataArray(
                    component.extra_dimensions[dimension], dims=(dimension,))
                additional_dimensions[dimension] = output_state[dimension]

        if hasattr(component, 'quantity_descriptions'):
            ensure_no_shared_keys(additional_descriptions,
                                  component.quantity_descriptions)
            additional_descriptions.update(component.quantity_descriptions)

    temporary_description.update(additional_descriptions)
    for name in quantity_list:

        if name in output_state.keys():
            continue

        output_state[name] = get_default_values(
            name, x_coordinate, y_coordinate,
            mid_levels_coordinate,
            interface_levels_coordinate,
            initial_state,
            temporary_description,
            additional_dimensions)

        if use_2d_coordinate:

            for physical_dimension in two_dim_coord_dict.keys():

                output_state[name].coords[physical_dimension] = (
                    two_dim_coord_dict[physical_dimension]['logical_dims'],
                    two_dim_coord_dict[physical_dimension]['values'])

                output_state[name][physical_dimension].attrs['units'] = \
                    two_dim_coord_dict[physical_dimension]['units']

    ensure_no_shared_keys(initial_state, output_state)
    output_state.update(initial_state)
    if 'time' not in initial_state:
        output_state['time'] = datetime(1, 1, 1)
    return output_state


def get_default_values(quantity_name, x, y, mid_levels, interface_levels,
                       initial_state, quantity_description, additional_dimensions={}):
    """
    Returns default values for individual quantities.

    Args:
    quantity_name (string): string corresponding to a quantity in the pre-defined
        dictionary of names in :code:`quantity_description`.

    x (DataArray): DataArray containing the definition of the x coordinates.

    y (DataArray): DataArray containing the definition of the y coordinates.

    z (DataArray): DataArray containing the definition of the mid-level z coordinates.
        Quantities on interface levels will have a vertical dimension of size :code:`len(z)+1`.

    quantity_description (dict): Description of dimensions, units and default values of
        quantities to be used to create DataArrays.

    additional_dimensions (dict): Dictionary containing names and coordinates of dimensions
        that are not x,y,z (used for spectral bands in radiative codes, for example).

    """

    description = quantity_description[quantity_name]

    dimension_length = {}
    dimension_length['x'] = x.values.shape[0]
    dimension_length['y'] = y.values.shape[0]
    dimension_length['mid_levels'] = len(mid_levels.values)
    dimension_length['interface_levels'] = len(interface_levels.values)

    mid_level_coords = np.arange(dimension_length['mid_levels'])
    int_level_coords = np.arange(dimension_length['interface_levels'])
    full_coords = {}
    full_coords['x'] = x
    full_coords['y'] = y
    full_coords['mid_levels'] = mid_level_coords
    full_coords['interface_levels'] = int_level_coords

    for dimension in additional_dimensions.keys():
        dimension_length[dimension] = additional_dimensions[dimension].shape[0]
        full_coords[dimension] = additional_dimensions[dimension].values

    quantity_dims = list(description['dims'])

    array_dims = [dimension_length[dimension] for dimension in quantity_dims]
    quantity_coords = [full_coords[dimension] for dimension in quantity_dims]

    quantity_dims = [x.label if elem is 'x' else elem for elem in quantity_dims]
    quantity_dims = [y.label if elem is 'y' else elem for elem in quantity_dims]

    if 'dtype' in description:
        dtype = description['dtype']
    else:
        dtype = 'float64'

    if 'default_value' in description:
        if dtype in ['a100', object]:
            quantity_array = np.ones(array_dims, dtype=dtype)
            quantity_array[:] = description['default_value']
        else:
            quantity_array = np.ones(array_dims, order='F', dtype=dtype)*description['default_value']
    elif 'init_func' in description:
        quantity_array = description['init_func'](array_dims, quantity_description, initial_state)
    else:
        raise ValueError(
            'Malformed description for quantity {}:\
            must contain default_value or init_func'.format(quantity_name))

    return DataArray(quantity_array, dims=quantity_dims,
                     coords=quantity_coords, attrs={'units': description['units']})
