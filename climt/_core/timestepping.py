from .base_components import PrognosticComposite
import abc


class TimeStepper(object):
    """An object which integrates model state forward in time.

    It uses Prognostic and Diagnostic objects to update the current model state
    with diagnostics, and to return the model state at the next timestep.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, prognostic_list, **kwargs):
        self._prognostic = PrognosticComposite(*prognostic_list)

    @abc.abstractmethod
    def step(self, state, timestep):
        """
        Retrieves any diagnostics and returns a new state corresponding
        to the next timestep.

        Args:
            state (dict): The current model state.
            timestep (timedelta): The amount of time to step forward.

        Returns:
            diagnostics (dict): Diagnostics from the timestep of the input
                state.
            new_state (dict): The model state at the next timestep.
        """

    def _copy_untouched_quantities(self, old_state, new_state):
        for key in old_state.keys():
            if key not in new_state:
                new_state[key] = old_state[key]

    @property
    def inputs(self):
        return self._prognostic.inputs

    @property
    def outputs(self):
        return self._prognostic.tendencies


class AdamsBashforth(TimeStepper):
    """A TimeStepper using the Adams-Bashforth scheme."""

    def __init__(self, prognostic_list, order=3):
        """
        Initialize an Adams-Bashforth time stepper.

        Args:
            prognostic_list (iterable of Prognostic): Objects used to get
                tendencies for time stepping.
            diagnostic_list (iterable of Diagnostic): Objects used to get
                diagnostics before time stepping.
            order (int, optional): The order of accuracy to use. Must be between
                1 and 4. 1 is the same as the Euler method. Default is 3.
        """
        if isinstance(order, float) and order.is_integer():
            order = int(order)
        if not isinstance(order, int):
            raise TypeError('order must be an integer')
        if not 1 <= order <= 4:
            raise ValueError('order must be between 1 and 4')
        self._order = order
        self._timestep = None
        self._tendencies_list = []
        super(AdamsBashforth, self).__init__(prognostic_list)

    def step(self, state, timestep):
        """
        Updates the input state dictionary and returns a new state corresponding
        to the next timestep.

        Args:
            state (dict): The current model state. Will be updated in-place by
                the call with any diagnostics from the current timestep.
            timestep (timedelta): The amount of time to step forward.

        Returns:
            new_state (dict): The model state at the next timestep.

        Raises:
            ValueError: If the timestep is not the same as the last time
                step() was called on this instance of this object.
        """
        self._ensure_constant_timestep(timestep)
        state = state.copy()
        tendencies, diagnostics = self._prognostic(state)
        self._tendencies_list.append(tendencies)
        new_state = self._perform_step(state, timestep)
        self._copy_untouched_quantities(state, new_state)
        if len(self._tendencies_list) == self._order:
            self._tendencies_list.pop(0)  # remove the oldest entry
        return diagnostics, new_state

    def _perform_step(self, state, timestep):
        # if we don't have enough previous tendencies built up, use lower order
        order = min(self._order, len(self._tendencies_list))
        if order == 1:
            new_state = step_forward_euler(
                state, self._tendencies_list[-1], timestep)
        elif order == 2:
            new_state = second_bashforth(state, self._tendencies_list, timestep)
        elif order == 3:
            new_state = third_bashforth(state, self._tendencies_list, timestep)
        elif order == 4:
            new_state = fourth_bashforth(state, self._tendencies_list, timestep)
        else:
            # the following should never happen, if it is there's a bug
            raise RuntimeError('order should be integer between 1 and 4')
        return new_state

    def _ensure_constant_timestep(self, timestep):
        if self._timestep is None:
            self._timestep = timestep
        elif self._timestep != timestep:
            raise ValueError(
                'timestep must be constant for Adams-Bashforth time stepping')


class Leapfrog(TimeStepper):
    """A TimeStepper using the Leapfrog scheme.

    This scheme calculates the
    values at time $t_{n+1}$ using the derivatives at $t_{n}$ and values at
    $t_{n-1}$. Following the step, an Asselin filter is applied to damp the
    computational mode that results from the scheme and maintain stability. The
    Asselin filter brings the values at $t_{n}$ (and optionally the values at
    $t_{n+1}$, according to Williams (2009)) closer to the mean of the values
    at $t_{n-1}$ and $t_{n+1}$."""

    def __init__(
            self, prognostic_list, asselin_strength=0.05,
            alpha=0.5):
        """
        Initialize a Leapfrog time stepper.

        Args:
            prognostic_list (iterable of Prognostic): Objects used to get
                tendencies for time stepping.
            diagnostic_list (iterable of Diagnostic): Objects used to get
                diagnostics before time stepping.
            asselin_strength (float, optional): The filter parameter used to
                determine the strength of the Asselin filter. Default is 0.05.
            alpha (float, optional): Constant from Williams (2009), where the
                midpoint is shifted by alpha*influence, and the right point is
                shifted by (1-alpha)*influence. If alpha is 1 then the behavior
                is that of the classic Robert-Asselin time filter, while if it
                is 0.5 the filter will conserve the three-point mean.
                Default is 0.5.

        References:
            Williams, P., 2009: A Proposed Modification to the Robert-Asselin
                Time Filter. Mon. Wea. Rev., 137, 2538--2546,
                doi: 10.1175/2009MWR2724.1.
        """
        self._old_state = None
        self._asselin_strength = asselin_strength
        self._timestep = None
        self._alpha = alpha
        super(Leapfrog, self).__init__(prognostic_list)

    def step(self, state, timestep):
        """
        Updates the input state dictionary and returns a new state corresponding
        to the next timestep.

        Args:
            state (dict): The current model state. Will be updated in-place by
                the call due to the Robert-Asselin-Williams filter.
            timestep (timedelta): The amount of time to step forward.

        Returns:
            new_state (dict): The model state at the next timestep.

        Raises:
            SharedKeyException: If a Diagnostic object has an output that is
                already in the state at the start of the timestep.
            ValueError: If the timestep is not the same as the last time
                step() was called on this instance of this object.
        """
        original_state = state
        state = state.copy()
        self._ensure_constant_timestep(timestep)
        tendencies, diagnostics = self._prognostic(state)
        if self._old_state is None:
            new_state = step_forward_euler(state, tendencies, timestep)
        else:
            state, new_state = step_leapfrog(
                self._old_state, state, tendencies, timestep,
                asselin_strength=self._asselin_strength, alpha=self._alpha)
        self._copy_untouched_quantities(state, new_state)
        self._old_state = state
        for key in original_state.keys():
            original_state[key] = state[key]  # allow filtering to be applied
        return diagnostics, new_state

    def _ensure_constant_timestep(self, timestep):
        if self._timestep is None:
            self._timestep = timestep
        elif self._timestep != timestep:
            raise ValueError(
                'timestep must be constant for Leapfrog time stepping')


def step_leapfrog(
        old_state, state, tendencies, timestep, asselin_strength=0.05,
        alpha=1.):
    """
    Steps the model state forward in time using the given tendencies and the
    leapfrog time scheme, with a Robert-Asselin time filter.

    Args:
        old_state (dict): Model state at the last timestep.
        state (dict): Model state at the current timestep. May be modified by
            this function call, specifically by the Asselin filter.
        tendencies (dict): Time derivatives at the current timestep in
            units/second.
        timestep (timedelta): The amount of time to step forward.
        asselin_strength (float, optional): Asselin filter strength.
            Default is 0.05.
        alpha (float, optional): Constant from Williams (2009), where the
            midpoint is shifted by alpha*influence, and the right point is
            shifted by (alpha-1)*influence. If alpha is 1 then the behavior
            is that of the classic Robert-Asselin time filter, while if it
            is 0.5 the filter will conserve the three-point mean.
            Default is 1.

    Returns:
        state (dict): The input state, modified in place.
        new_state (dict): Model state at the next timestep.
    """
    new_state = {}
    for key in tendencies.keys():
        new_state[key] = (
            old_state[key] + 2*tendencies[key]*timestep.total_seconds())
        filter_influence = 0.5*asselin_strength*(
            old_state[key] - 2*state[key] + new_state[key])
        state[key] += alpha * filter_influence
        if alpha != 1.:
            new_state[key] += (alpha - 1.) * filter_influence
    return state, new_state


def step_forward_euler(state, tendencies, timestep):
    return_state = {}
    for key in tendencies.keys():
        return_state[key] = state[key] + tendencies[key]*timestep.total_seconds()
    return return_state


def second_bashforth(state, tendencies_list, timestep):
    """Return the new state using second-order Adams-Bashforth. tendencies_list
    should be a list of dictionaries whose values are tendencies in
    units/second (from oldest to newest), and timestep should be a timedelta
    object. The dictionaries in tendencies_list should all have the same keys.
    """
    return_state = {}
    for key in tendencies_list[0].keys():
        return_state[key] = state[key] + timestep.total_seconds() * (
            1.5*tendencies_list[-1][key] - 0.5*tendencies_list[-2][key]
        )
    return return_state


def third_bashforth(state, tendencies_list, timestep):
    """Return the new state using third-order Adams-Bashforth. tendencies_list
    should be a list of dictionaries whose values are tendencies in
    units/second (from oldest to newest), and timestep should be a timedelta
    object."""
    return_state = {}
    for key in tendencies_list[0].keys():
        return_state[key] = state[key] + timestep.total_seconds() * (
            23./12*tendencies_list[-1][key] - 4./3*tendencies_list[-2][key] +
            5./12*tendencies_list[-3][key]
        )
    return return_state


def fourth_bashforth(state, tendencies_list, timestep):
    """Return the new state using fourth-order Adams-Bashforth. tendencies_list
    should be a list of dictionaries whose values are tendencies in
    units/second (from oldest to newest), and timestep should be a timedelta
    object."""
    return_state = {}
    for key in tendencies_list[0].keys():
        return_state[key] = state[key] + timestep.total_seconds() * (
            55./24*tendencies_list[-1][key] - 59./24*tendencies_list[-2][key] +
            37./24*tendencies_list[-3][key] - 3./8*tendencies_list[-4][key]
        )
    return return_state
