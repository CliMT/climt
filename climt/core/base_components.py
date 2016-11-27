import abc

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
