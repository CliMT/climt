import abc
from .util import ensure_no_shared_keys, add_dicts_inplace


class Implicit(object):
    """
    Attributes:
        inputs (tuple of str): The quantities required in the state when the
            object is called.
        outputs (tuple of str): The quantities for which values for the new
            state are returned when the object is called.
        diagnostics (tuple of str): The diagnostic quantities that are inserted
            into the old state when the object is called.
    """
    __metaclass__ = abc.ABCMeta

    inputs = ()
    outputs = ()
    diagnostics = ()

    @abc.abstractmethod
    def __call__(self, state, timestep):
        """
        Gets diagnostics from the current model state and steps the state
        forward in time according to the timestep.

        Args:
            state (dict): A model state dictionary. Will be updated with any
                diagnostic quantities produced by this object for the time of
                the input state.

        Returns:
            next_state (dict): A dicitonary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the timestep after input state.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Implicit instance for other reasons.
        """


class Prognostic(object):
    """
    Attributes:
        inputs (tuple of str): The quantities required in the state when the
            object is called.
        tendencies (tuple of str): The quantities for which tendencies are
            returned when the object is called.
        diagnostics (tuple of str): The diagnostic quantities returned when
            the object is called.
    """
    __metaclass__ = abc.ABCMeta

    inputs = ()
    tendencies = ()
    diagnostics = ()

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
            diagnostics (dict): A dictionary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the time of the input state.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Prognostic instance.
        """


class Diagnostic(object):
    """
    Attributes:
        inputs (tuple of str): The quantities required in the state when the
            object is called.
        diagnostics (tuple of str): The diagnostic quantities returned when
            the object is called.
    """
    __metaclass__ = abc.ABCMeta

    inputs = ()
    diagnostics = ()

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

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Prognostic instance.
        """

    def update_state(self, state):
        """
        Gets diagnostics from the passed model state and updates the state with
        those diagnostics (in place).

        Args:
            state (dict): A model state dictionary.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Diagnostic instance.
        """
        diagnostics = self(state)
        ensure_no_shared_keys(state, diagnostics)
        state.update(diagnostics)


class Monitor(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def store(self, state):
        """
        Stores the given state in the Monitor and performs class-specific
        actions.

        Args:
            state (dict): A model state dictionary.

        Raises:
            InvalidStateException: if state is not a valid input for the
                Diagnostic instance.
        """


class ComponentComposite(object):

    component_class = None

    def __init__(self, component_list):
        """
        Args:
            component_list (iterable): The components that should be
                wrapped by this object.
        """
        if self.component_class is not None:
            ensure_components_have_class(component_list, self.component_class)
        self._components = component_list

    def _combine_attribute(self, attr):
        return_attr = []
        for component in self._components:
            return_attr.extend(getattr(component, attr))
        return tuple(set(return_attr))  # set to deduplicate

    @property
    def inputs(self):
        return self._combine_attribute('inputs')

    @property
    def diagnostics(self):
        return self._combine_attribute('diagnostics')


def ensure_components_have_class(components, component_class):
    for component in components:
        if not isinstance(component, component_class):
            raise TypeError(
                "require components of type {}".format(component_class))


class PrognosticComposite(ComponentComposite):

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
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for a
                Prognostic instance.
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

    @property
    def tendencies(self):
        return self._combine_attribute('tendencies')


class DiagnosticComposite(ComponentComposite):

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
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for a
                Diagnostic instance.
        """
        return_diagnostics = {}
        for diagnostic_component in self._components:
            diagnostics = diagnostic_component(state)
            # ensure we won't overwrite already-existing diagnostics
            ensure_no_shared_keys(return_diagnostics, diagnostics)
            return_diagnostics.update(diagnostics)
        return return_diagnostics

    def update_state(self, state):
        """
        Gets diagnostics from the passed model state and updates the state with those
        diagnostics (in place).

        Args:
            state (dict): A model state dictionary.

        Raises:
            InvalidStateException: if state already includes any diagnostics being output.
        """
        for diagnostic_component in self._components:
            diagnostic_component.update_state(state)


class MonitorComposite(ComponentComposite):

    component_class = Monitor

    def store(self, state):
        """
        Stores the given state in the Monitor and performs class-specific
        actions.

        Args:
            state (dict): A model state dictionary.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for a
                Monitor instance.
        """
        for monitor in self._components:
            monitor.store(state)
