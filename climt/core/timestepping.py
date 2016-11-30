from .base_components import PrognosticCollection, DiagnosticCollection
from .util import ensure_no_shared_keys
import abc


class TimeStepper(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, prognostic_list, diagnostic_list=(), **kwargs):
        self._prognostic = PrognosticCollection(prognostic_list)
        self._diagnostic = DiagnosticCollection(diagnostic_list)

    @abc.abstractmethod
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
        """


class AdamsBashforth(TimeStepper):

    def __init__(self, prognostic_list, diagnostic_list=(), order=2):
        """
        Initialize an Adams-Bashforth time stepper.

        Args:
            prognostic_list (iterable of Prognostic): Objects used to get
                tendencies for time stepping.
            diagnostic_list (iterable of Diagnostic): Objects used to get
                diagnostics before time stepping.
            order (int, optional): The order of accuracy to use. Must be between
                1 and 4. 1 is the same as the Euler method. Default is 2.
        """
        self._order = order
        self._timestep = None
        self._tendencies_list = []
        super(AdamsBashforth, self).__init__(prognostic_list, diagnostic_list)

    def step(self, state, timestep):
        self._ensure_constant_timestep(timestep)
        self._diagnostic.update_state(state)
        tendencies, diagnostics = self._prognostic(state)
        ensure_no_shared_keys(state, diagnostics)
        state.update(diagnostics)
        self._tendencies_list.append(tendencies)
        # if we don't have enough previous tendencies built up, use lower order
        order = min(self._order, len(self._tendencies_list))
        if order == 1:
            new_state = step_forward_euler(state, tendencies, timestep)
        elif order == 2:
            new_state = second_bashforth(state, self._tendencies_list, timestep)
        elif order == 3:
            new_state = third_bashforth(state, self._tendencies_list, timestep)
        elif order == 4:
            new_state = fourth_bashforth(state, self._tendencies_list, timestep)
        if len(self._tendencies_list) == self._order:
            self._tendencies_list.pop(0)  # remove the oldest entry
        return new_state

    def _ensure_constant_timestep(self, timestep):
        if self._timestep is None:
            self._timestep = timestep
        elif self._timestep != timestep:
            raise ValueError('timestep must be constant for Leapfrog time stepping')


class Leapfrog(TimeStepper):

    def __init__(
            self, prognostic_list, diagnostic_list=(), asselin_strength=0.05,
            alpha=1.):
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
                Default is 1.

        References:
            Williams, P., 2009: A Proposed Modification to the Robert-Asselin
                Time Filter. Mon. Wea. Rev., 137, 2538--2546,
                doi: 10.1175/2009MWR2724.1.
        """
        self._old_state = None
        self._asselin_strength = asselin_strength
        self._timestep = None
        self._alpha = alpha
        super(Leapfrog, self).__init__(prognostic_list, diagnostic_list)

    def step(self, state, timestep):
        self._ensure_constant_timestep(timestep)
        self._diagnostic.update_state(state)
        tendencies, diagnostics = self._prognostic(state)
        ensure_no_shared_keys(state, diagnostics)
        state.update(diagnostics)
        if self._old_state is None:
            new_state = step_forward_euler(state, tendencies, timestep)
        else:
            state, new_state = step_leapfrog(
                self._old_state, state, tendencies, timestep,
                asselin_strength=self._asselin_strength, alpha=self._alpha)
        self._old_state = state
        return new_state

    def _ensure_constant_timestep(self, timestep):
        if self._timestep is None:
            self._timestep = timestep
        elif self._timestep != timestep:
            raise ValueError('timestep must be constant for Leapfrog time stepping')


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
