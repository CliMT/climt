from sympl import (DataArray, set_dimension_names,
                   ensure_no_shared_keys)
import numpy as np

def init_mid_level_pressures(array_dims, quantity_description):
    #We assume the vertical level is the last dimension
    vert_levels = array_dims[-1]

    p_surf = quantity_description['surface_pressure']['init_value']
    spacing = np.linspace(0.995, 0.001, vert_levels)
    single_column = p_surf*spacing
    single_column = single_column[np.newaxis, np.newaxis, :]
    
    return single_column*np.ones(array_dims)

def init_interface_level_pressures(array_dims, quantity_description):

    vert_levels = array_dims[-1]
    p_surf = quantity_description['surface_pressure']['init_value']
    spacing = np.linspace(0.995, 0.001, vert_levels-1)
    midlevel = p_surf*spacing
    
    interface = np.zeros(vert_levels)
    interface[1:-1] = 0.5*(midlevel[:-1] + midlevel[1:])
    interface[0] = p_surf
    interface[-1] = 0.0005*p_surf

    return interface*np.ones(array_dims)

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
    'surface_pressure': {
        'dims': ['x', 'y'],
        'units': 'Pa',
        'init_value': 1.e5
    },
    'air_temperature': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'degK',
        'init_value': 290.
    },
    'air_temperature_on_interface_levels': {
        'dims': ['x', 'y', 'interface_levels'],
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
    'surface_emissivity': {
        'dims': ['x', 'y'],
        'units': 'dimensionless',
        'init_value': 1.
    },
    'specific_humidity': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'g/kg',
        'init_value': 0.
    },
    'ozone_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'carbon_dioxide_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 330e-6
    },
    'methane_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'nitrous_oxide_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'oxygen_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.21
    },
    'nitrogen_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.78
    },
    'hydrogen_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 500e-9
    },
    'cfc11_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'cfc12_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'cfc22_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'ccl4_mixing_ratio': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'mole/mole',
        'init_value': 0.
    },
    'cloud_fraction': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_value': 0.
    },
    'cloud_optical_depth': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_value': 0.
    },
    'cloud_ice_water_path': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'g m^-2',
        'init_value': 0.
    },
    'cloud_liquid_water_path': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'g m^-2',
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
    'aerosol_optical_depth': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'dimensionless',
        'init_value': 0.
    },
    'longwave_heating_rate': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'init_value': 0.
    },
    'longwave_heating_rate_clearsky': {
        'dims': ['x', 'y', 'mid_levels'],
        'units': 'K day^-1',
        'init_value': 0.
    },
    'upward_longwave_flux': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'upward_longwave_flux_clearsky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'downward_longwave_flux': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
    'downward_longwave_flux_clearsky': {
        'dims': ['x', 'y', 'interface_levels'],
        'units': 'W m^-2',
        'init_value': 0.
    },
}

def get_default_state(component_list, x={}, y={}, z={}, input_state={}):
    """
    Returns a state dictionary with the required arrays to run a model
    which comprises of components in ``component_list``. If coordinate
    values in ``x``,``y`` and ``z`` are not provided, a single column
    centered at 0 degrees south, 0 degrees east of 30 levels is used.

    Args:
        component_list (iterable): The _components for which a default
            state is required, in the order that they are called.
            It is assumed diagnostic outputs are updated into
            the state and passed into the next component, and that the state
            from the previous component is passed into the next component. The
            function will attempt to determine required inputs from the series
            of _components.

        x (dict,optional): A dictionary containing keys ``label``, ``values``,
            and ``units``. ``label`` refers to the name the coordinate axis will assume.
            ``values`` refers to the array of coordinate values. ``units`` refers to the
            units of the coordinate values.
            If ``x`` is an empty dictionary, a single default value of 0 degrees longitude
            is used.

        y (dict,optional): A dictionary containing keys ``label``, ``values``,
            and ``units``. ``label`` refers to the name the coordinate axis will assume.
            ``values`` refers to the array of coordinate values. ``units`` refers to the
            units of the coordinate values.
            If ``y`` is an empty dictionary, a single default value of 0 degrees latitude
            is used.

        z (dict,optional): A dictionary containing keys ``label``, ``values``,
            and ``units``. ``label`` refers to the name the coordinate axis will assume.
            ``values`` refers to the array of coordinate values. ``units`` refers to the
            units of the coordinate values.
            If ``z`` is an empty dictionary, 30 levels of arbitrary units are used.

        input_state (dict, optional): A dictionary containing some quantities that will
            also be added to the final output state. Must not contain any quantities that
            the output state will overwrite.


    Returns:
        default_state (dict): A state dictionary containing the requested
            quantities using the provided coordinate state.
    """
    
    if len(component_list) == 0:
        raise ValueError('Component list must contain at least one component')

    output_state = {}

    if len(x.keys()) == 0:
        x_coordinate = DataArray(
            np.zeros((1,)), dims=('x',), attrs={'units':'degrees_north', 'label': 'longitude'})

        output_state['longitude'] = x_coordinate
        set_dimension_names(x='longitude')
    else:
        x_coordinate = DataArray(
            x['values'], dims=('x',), attrs={'units': x['units'], 'label': x['label']})

        output_state[x['label']] = x_coordinate
        set_dimension_names(x=x['label'])

    if len(y.keys()) == 0:
        y_coordinate = DataArray(
            np.zeros((1,)), dims=('y',), attrs={'units':'degrees_east', 'label': 'latitude'})

        output_state['latitude'] = y_coordinate
        set_dimension_names(y='latitude')
    else:
        y_coordinate = DataArray(
            x['values'], dims=('y',), attrs={'units': y['units'], 'label': y['label']})

        output_state[y['label']] = y_coordinate
        set_dimension_names(y=y['label'])

    if len(z.keys()) == 0:
        z_coordinate = DataArray(
            np.arange(30), dims=('mid_levels',), attrs={'units': '', 'label': 'mid_levels'})

        output_state['mid_levels'] = z_coordinate
    else:
        z_coordinate = DataArray(
            z['values'], dims=('mid_levels',), attrs={'units': z['units'], 'label': z['label']})

        output_state[z['label']] = z_coordinate
        #set_dimension_names(z=z['label'])

    quantity_list = set()
    temporary_description = _quantity_descriptions.copy()
    additional_dimensions = {}
    additional_descriptions = {}

    for component in component_list:
        quantity_list = quantity_list.union(set(component.inputs))

        if hasattr(component, 'extra_dimensions'):
            ensure_no_shared_keys(additional_dimensions,
                                  component.extra_dimensions)
            additional_dimensions.update(component.extra_dimensions)
            for dimension in component.extra_dimensions.keys():
                output_state[dimension] = DataArray(
                    component.extra_dimensions[dimension], dims=(dimension,))

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

    ensure_no_shared_keys(input_state, output_state)
    output_state.update(input_state)
    return output_state


def get_default_values(quantity_name, x, y, z, 
                       quantity_description, additional_dimensions={}):
    """
    Returns default values for individual quantities.

    Args:
    quantity_name (string): string corresponding to a quantity in the pre-defined
        dictionary of names in ``quantity_description``.

    x (DataArray): DataArray containing the definition of the x coordinates.

    y (DataArray): DataArray containing the definition of the y coordinates.

    z (DataArray): DataArray containing the definition of the mid-level z coordinates.
        Quantities on interface levels will have a vertical dimension of size ``len(z)+1``.

    quantity_description (dict): Description of dimensions, units and default values of
        quantities to be used to create DataArrays.

    additional_dimensions (dict): Dictionary containing names and coordinates of dimensions
        that are not x,y,z (used for spectral bands in radiative codes, for example).

    """


    description = quantity_description[quantity_name]

    dimension_length = {}
    dimension_length['x'] = len(x.values)
    dimension_length['y'] = len(y.values)
    dimension_length['mid_levels'] = len(z.values)
    dimension_length['interface_levels'] = dimension_length['mid_levels']+1

    for dimension in additional_dimensions.keys():
        dimension_length[dimension] = len(additional_dimensions[dimension])

    quantity_dims = description['dims']

    array_dims = [dimension_length[dimension] for dimension in quantity_dims]

    quantity_dims = [x.label if elem == 'x' else elem for elem in quantity_dims]
    quantity_dims = [y.label if elem == 'y' else elem for elem in quantity_dims]

    if 'init_value' in description:
        quantity_array = np.ones(array_dims, order='F')*description['init_value']
    elif 'init_func' in description:
        quantity_array = description['init_func'](array_dims, quantity_description)
    else:
        raise ValueError('Malformed description for quantity {}'.format(quantity_name))

    return DataArray(quantity_array, dims= quantity_dims, attrs={'units': description['units']})
