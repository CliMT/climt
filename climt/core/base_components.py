import abc
from .exceptions import SharedKeyException


class Implicit(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def ensure_state_is_valid_input(self, state):
        """
        Args:
            state (dict): A model state dictionary.

        Raises:
            InvalidStateException: if state is not a valid input for the
                Implicit instance.

        Returns:
            None
        """
        return

    @abc.abstractmethod
    def __call__(self, state, timestep):
        """
        Gets diagnostics from the current model state and steps the state
        forward in time according to the timestep.

        Args:
            state (dict): A model state dictionary.

        Returns:
            current_diagnostics (dict): A dicitonary whose keys are strings
                indicating state quantities and values are the value of those
                quantities at the time of the input state.
            next_state (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the timestep after input state.
        """
        return


class Prognostic(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def ensure_state_is_valid_input(self, state):
        """
        Args:
            state (dict): A model state dictionary.

        Raises:
            InvalidStateException: if state is not a valid input for the
                Prognostic instance.

        Returns:
            None
        """
        return

    @abc.abstractmethod
    def __call__(self, state):
        """
        Gets tendencies and diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            tendencies (dict): A dictionary whose keys are strings indicating
                state quantities and values are the time derivative of those
                quantities in units/second at the time of the input state.
            diagnostics (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.
        """
        return


class Diagnostic(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def ensure_state_is_valid_input(self, state):
        """
        Args:
            state (dict): A model state dictionary.

        Raises:
            InvalidStateException: if state is not a valid input for the
                Diagnostic instance.

        Returns:
            None
        """
        return

    @abc.abstractmethod
    def __call__(self, state):
        """
        Gets diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            diagnostics (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.
        """
        return


class Monitor(object):
    __metaclass__ = abc.ABCMeta

    def ensure_state_is_valid_input(self, state):
        """
        Args:
            state (dict): A model state dictionary.

        Raises:
            InvalidStateException: if state is not a valid input for the
                Diagnostic instance.

        Returns:
            None
        """
        return  # by default all states are valid

    @abc.abstractmethod
    def store(self, state):
        """
        Stores the given state in the Monitor and performs class-specific
        actions.

        Args:
            state (dict): A model state dictionary.

        Returns:
            None
        """
        return


class ComponentCollection(object):

    component_class = None

    def __init__(self, component_list):
        """
        Args:
            component_list (iterable): the components that should be
                wrapped by this object.
        """
        if self.component_class is not None:
            ensure_components_have_class(component_list, self.component_class)
        self._components = component_list

    def ensure_state_is_valid_input(self, state):
        for component in self._components:
            component.ensure_state_is_valid_input(state)


def ensure_components_have_class(components, component_class):
    for component in components:
        if not isinstance(component, component_class):
            raise ValueError(
                "Components must be of class {}".format(component_class))


class PrognosticCollection(ComponentCollection):

    component_class = Prognostic

    def __call__(self, state):
        """
        Gets tendencies and diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            tendencies (dict): A dictionary whose keys are strings indicating
                state quantities and values are the time derivative of those
                quantities in units/second at the time of the input state.
            diagnostics (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.

        Raises:
            SharedKeyException: if multiple Prognostic objects contained in the
            collection return the same diagnostic quantity.
        """
        return_tendencies = {}
        return_diagnostics = {}
        for prognostic in self._components:
            tendencies, diagnostics = prognostic(state)
            add_dicts_inplace(return_tendencies, tendencies)
            # ensure we won't overwrite already-existing diagnostics
            ensure_no_shared_keys(return_diagnostics, diagnostics)
            return_diagnostics.update(diagnostics)
        return return_tendencies, return_diagnostics


def add_dicts_inplace(dict1, dict2):
    """
    Takes two dictionaries. For any keys in both dictionary, it adds the value
    in dict2 to the value in dict1. This is done in-place if the values are
    array-like, to avoid data copying. None is returned.
    """
    for key in dict2.keys():
        if key not in dict1:
            dict1[key] = dict2[key]
        else:
            try:
                # works for array-like objects, in-place
                dict1[key][:] += dict2[key][:]
            except TypeError:
                dict1[key] += dict2[key]
    return  # not returning anything emphasizes that this is in-place


def ensure_no_shared_keys(dict1, dict2):
    """
    Raises SharedKeyException if there exists a key present in both
    dictionaries.
    """
    if len(set(dict1.keys()).intersection(dict2.keys())) > 0:
        raise SharedKeyException()


class DiagnosticCollection(ComponentCollection):

    component_class = Diagnostic

    def __call__(self, state):
        """
        Gets diagnostics from the passed model state.

        Args:
            state (dict): A model state dictionary.

        Returns:
            diagnostics (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.

        Raises:
            SharedKeyException: if multiple Diagnostic objects contained in the
            collection return the same diagnostic quantity.
        """
        return_diagnostics = {}
        for diagnostic_component in self._components:
            diagnostics = diagnostic_component(state)
            # ensure we won't overwrite already-existing diagnostics
            ensure_no_shared_keys(return_diagnostics, diagnostics)
            return_diagnostics.update(diagnostics)
        return return_diagnostics


class MonitorCollection(ComponentCollection):

    component_class = Monitor

    def store(self, state):
        """
        Stores the given state in the Monitor and performs class-specific
        actions.

        Args:
            state (dict): A model state dictionary.

        Returns:
            None
        """
        for monitor in self._components:
            monitor.store(state)
