from sympl import set_constant, get_constants_string


def list_available_constants():
    """
    Prints all the constants currently registered with sympl.
    """
    print(get_constants_string())


def set_constants_from_dict(constant_descriptions):
    """
    Modify/Add constants in the library.

    Args:

       constant_descriptions (dict):
           Dictionary containing the description of the constants.
           The key should be the name of the constant, and the value
           should be a dictionary containing the following keys:

           * value (float):
               The value assigned.

           * units (string):
                   The units of the value, e.g, m/s, J/kg.

    """

    for name in constant_descriptions.keys():

        value = constant_descriptions[name]['value']
        units = constant_descriptions[name]['units']

        set_constant(name, value, units)
