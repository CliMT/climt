import pint

unit_registry = pint.UnitRegistry()
unit_registry.define('degrees_north = degree_north = degree_N = degrees_N = degreeN = degreesN')
unit_registry.define('degrees_east = degree_east = degree_E = degrees_E = degreeE = degreesE')
unit_registry.define('percent = 0.01*count = %')


def to_units(value, units):
    if not hasattr(value, 'attrs') or 'units' not in value.attrs:
        raise TypeError(
            'Cannot retrieve units from type {}'.format(type(value)))
    elif unit_registry(value.attrs['units']) != unit_registry(units):
        attrs = value.attrs
        value = (unit_registry(value.attrs['units'])*value).to(units).magnitude
        attrs['units'] = units
        value.attrs = attrs
    return value
