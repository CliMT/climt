from sympl import Diagnostic, DiagnosticComposite

_quantity_func_dictionary = {}


def get_default_state(coordinate_state, quantity_list=(), component_list=()):
    """

    Args:
        coordinate_state (dict): A state dictionary containing coordinate-like
            quantities (as DataArrays) necessary to calculate the needed
            state quantities.
        quantity_list (iterable, optional): Names of quantities to include in
            the returned default state.
        component_list (iterable, optional): The _components for which a default
            state is required, in the order that they are called.
            It is assumed diagnostic outputs are updated into
            the state and passed into the next component, and that the state
            from the previous component is passed into the next component. The
            function will attempt to determine required inputs from the series
            of _components.

    Returns:
        default_state (dict): A state dictionary containing the requested
            quantities using the provided coordinate state.
    """
    quantity_list = set(quantity_list)
    last_outputs = ()
    for component in component_list:
        quantity_list.update(set(component.inputs).difference(last_outputs))
        if isinstance(component, (Diagnostic, DiagnosticComposite)):
            last_outputs = component.diagnostics
        else:
            last_outputs = component.outputs
    quantity_list.difference_update(coordinate_state.keys())
    output_state = coordinate_state.copy()
    for name in quantity_list:
        output_state[name] = _quantity_func_dictionary[name](coordinate_state)
