from sympl import (DataArray, set_dimension_names,
                   ensure_no_shared_keys)
import numpy as np


def init_mid_level_pressures(array_dims, quantity_description):
    # We assume the vertical level is the last dimension
    vert_levels = array_dims[-1]

    p_surf = quantity_description['surface_air_pressure']['init_value']
    spacing = np.linspace(0.995, 0.001, vert_levels)
    single_column = p_surf*spacing
    single_column = single_column[np.newaxis, np.newaxis, :]

    return single_column*np.ones(array_dims, order='F')


def init_interface_level_pressures(array_dims, quantity_description):

    vert_levels = array_dims[-1]
    p_surf = quantity_description['surface_air_pressure']['init_value']
    spacing = np.linspace(0.995, 0.001, vert_levels-1)
    midlevel = p_surf*spacing

    interface = np.zeros(vert_levels)
    interface[1:-1] = 0.5*(midlevel[:-1] + midlevel[1:])
    interface[0] = p_surf
    interface[-1] = 0.0005*p_surf

    return interface*np.ones(array_dims, order='F')


def init_interface_level_sigma(array_dims, quantity_description):

    vert_levels = array_dims[-1]
    spacing = np.linspace(0.995, 0.001, vert_levels-1)
    midlevel = spacing

    interface = np.zeros(vert_levels)
    interface[1:-1] = 0.5*(midlevel[:-1] + midlevel[1:])
    interface[0] = 1.
    interface[-1] = 0.0005

    return interface*np.ones(array_dims, order='F')


def init_mid_level_sigma(array_dims, quantity_description):

    vert_levels = array_dims[-1]
    spacing = np.linspace(0.995, 0.001, vert_levels)
    midlevel = spacing

    return midlevel*np.ones(array_dims, order='F')


_quantity_descriptions = {
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
    'surface_air_pressure': {
        'dims': ['x', 'y'],
        'units': 'Pa',
        'init_value': 1.e5
    },
    'air_temperature': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'degK',
        'init_value': 290.
    },
    'surface_temperature': {
        'dims': ['x', 'y'],
        'units': 'degK',
        'init_value': 300.
    },
    'northward_wind': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'm s^-1',
        'init_value': 0.
    },
    'eastward_wind': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'm s^-1',
        'init_value': 0.
    },
    'surface_longwave_emissivity': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'init_value': 1.
    },
    'specific_humidity': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'g/kg',
        'init_value': 0.
    },
    'mole_fraction_of_ozone_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'mole_fraction_of_carbon_dioxide_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 330e-6
    },
    'mole_fraction_of_methane_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'mole_fraction_of_nitrous_oxide_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'mole_fraction_of_oxygen_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.21
    },
    'mole_fraction_of_nitrogen_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.78
    },
    'mole_fraction_of_hydrogen_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 500e-9
    },
    'mole_fraction_of_cfc11_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'mole_fraction_of_cfc12_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'mole_fraction_of_cfc22_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'mole_fraction_of_carbon_tetrachloride_in_air': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'cloud_area_fraction_in_atmosphere_layer': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_value': 0.
    },
    'atmosphere_optical_thickness_due_to_aerosol': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_value': 0.
    },
    'mass_content_of_cloud_ice_in_atmosphere_layer': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'kg m^-2',
        'init_value': 0.
    },
    'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'kg m^-2',
        'init_value': 0.
    },
    'cloud_ice_particle_size': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'micrometer',
        'init_value': 20.
    },
    'cloud_water_droplet_radius': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'micrometer',
        'init_value': 10.
    },
    'atmosphere_optical_thickness_due_to_cloud': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_value': 0.
    },
    'longwave_heating_rate': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'init_value': 0.
    },
    'longwave_heating_rate_assuming_clear_sky': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'init_value': 0.
    },
    'upwelling_longwave_flux_in_air': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'upwelling_longwave_flux_in_air_assuming_clear_sky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'downwelling_longwave_flux_in_air': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'downwelling_longwave_flux_in_air_assuming_clear_sky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'precipitation_amount': {
        'dims': ['x', 'y'],
        'units': 'kg m^-2',
        'init_value': 0.
    },
    'convective_precipitation_amount': {
        'dims': ['x', 'y'],
        'units': 'kg m^-2',
        'init_value': 0.
    }
}


def get_default_state(component_list, x={}, y={}, z={}, input_state={}):
    """
    Return a state dictionary required to run the model.

    The model comprises of components in :code:`component_list`. If coordinate
    values in :code:`x`, :code:`y` and :code:`z` are not provided, a single column
    centered at 0 degrees south, 0 degrees east, with 30 vertical levels is used.

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

        z (dict,optional): A dictionary containing keys :code:`label`, :code:`values`,
            and :code:`units`. :code:`label` refers to the name the coordinate axis will assume.
            :code:`values` refers to the array of coordinate values. :code:`units` refers to the
            units of the coordinate values.
            If :code:`z` is an empty dictionary, 30 levels of arbitrary units are used.

        input_state (dict, optional): A dictionary containing some quantities that will
            also be added to the final output state. Must not contain any quantities that
            the output state will overwrite.


    Returns:
        default_state (dict):
            A state dictionary containing the requested
            quantities using the provided coordinate state.

    Raises:
        ValueError:
            if :code:`component_list` is empty or the shape of :code:`x['values']` and
            :code:`y['values']` is not the same.
    """

    if len(component_list) == 0:
        raise ValueError('Component list must contain at least one component')

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

    if len(z.keys()) == 0:
        z_coordinate_values = np.arange(30)
        z_coordinate_label = 'mid_levels'
        z_coordinate_units = ''
    else:
        z_coordinate_values = z['values']
        z_coordinate_label = z['label']
        z_coordinate_units = z['units']

    if not x_coordinate_values.ndim == y_coordinate_values.ndim:
        raise ValueError('x and y coordinates must have the same shape')

    use_2d_coordinate = False
    x_physical_coordinate_values = None
    x_physical_coordinate_label = None
    x_physical_coordinate_units = None

    y_physical_coordinate_values = None
    y_physical_coordinate_label = None
    y_physical_coordinate_units = None

    if x_coordinate_values.ndim == 2:
        if not x_coordinate_values.shape == y_coordinate_values.shape:
            raise ValueError(
                'If x and y are 2d coordinates, they must have the same shape')

        x_physical_coordinate_values = x_coordinate_values
        x_physical_coordinate_label = x_coordinate_label
        x_physical_coordinate_units = x_coordinate_units

        y_physical_coordinate_values = y_coordinate_values
        y_physical_coordinate_label = y_coordinate_label
        y_physical_coordinate_units = y_coordinate_units

        x_coordinate_values = np.arange(x_physical_coordinate_values.shape[0])
        y_coordinate_values = np.arange(x_physical_coordinate_values.shape[1])

        x_coordinate_label = 'logical_x_coordinate'
        y_coordinate_label = 'logical_y_coordinate'

        x_coordinate_units = ''
        y_coordinate_units = ''

        use_2d_coordinate = True

    x_coordinate = DataArray(
        x_coordinate_values,
        dims=('x_coordinate_label'),
        attrs={'units': x_coordinate_units, 'label': x_coordinate_label})

    output_state[x_coordinate_label] = x_coordinate
    set_dimension_names(x=x_coordinate_label)

    y_coordinate = DataArray(
        y_coordinate_values,
        dims=('y_coordinate_label'),
        attrs={'units': y_coordinate_units, 'label': y_coordinate_label})

    output_state[y_coordinate_label] = y_coordinate
    set_dimension_names(y=y_coordinate_label)

    z_coordinate = DataArray(
        z_coordinate_values,
        dims=('z_coordinate_label',),
        attrs={'units': z_coordinate_units, 'label': z_coordinate_label})

    output_state[z_coordinate_label] = z_coordinate
    set_dimension_names(
        z=list({'mid_levels', 'interface_levels'}.union(z_coordinate_label)))

    output_state['x'] = x_coordinate
    output_state['y'] = y_coordinate
    output_state['z'] = z_coordinate

    quantity_list = set()
    temporary_description = _quantity_descriptions.copy()
    additional_dimensions = {}
    additional_descriptions = {}

    for component in component_list:
        quantity_list = quantity_list.union(set(component.inputs))

        if hasattr(component, 'extra_dimensions'):
            ensure_no_shared_keys(additional_dimensions,
                                  component.extra_dimensions)

            for dimension in component.extra_dimensions.keys():
                output_state[dimension] = DataArray(
                    component.extra_dimensions[dimension], dims=(dimension,))
                additional_dimensions[dimension] = output_state[dimension]

        if hasattr(component, 'quantity_descriptions'):
            ensure_no_shared_keys(additional_descriptions,
                                  component.quantity_descriptions)
            additional_descriptions.update(component.quantity_descriptions)

    temporary_description.update(additional_descriptions)
    for name in quantity_list:
        output_state[name] = get_default_values(
            name, x_coordinate, y_coordinate,
            z_coordinate, temporary_description,
            additional_dimensions)

        if use_2d_coordinate:

            output_state[name].coords[x_physical_coordinate_label] = (
                (x_coordinate.label, y_coordinate.label), x_physical_coordinate_values)

            output_state[name][x_physical_coordinate_label].attrs['units'] = x_physical_coordinate_units

            output_state[name].coords[y_physical_coordinate_label] = (
                (x_coordinate.label, y_coordinate.label), y_physical_coordinate_values)

            output_state[name][y_physical_coordinate_label].attrs['units'] = y_physical_coordinate_units

    ensure_no_shared_keys(input_state, output_state)
    output_state.update(input_state)
    return output_state


def get_default_values(quantity_name, x, y, z,
                       quantity_description, additional_dimensions={}):
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
    dimension_length['mid_levels'] = len(z.values)
    dimension_length['interface_levels'] = dimension_length['mid_levels']+1

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

    if 'init_value' in description:
        quantity_array = np.ones(array_dims, order='F')*description['init_value']
    elif 'init_func' in description:
        quantity_array = description['init_func'](array_dims, quantity_description)
    else:
        raise ValueError('Malformed description for quantity {}'.format(quantity_name))

    return DataArray(quantity_array, dims=quantity_dims,
                     coords=quantity_coords, attrs={'units': description['units']})
