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
        return state, {}


class AdamsBashforth(TimeStepper):
    pass


class Leapfrog(TimeStepper):

    def __init__(
            self, prognostic_list, diagnostic_list=(), asselin_strength=0.05,
            alpha=1.):
        """
        Initialize a Leapfrog time stepper.

        Args:
            prognostic_list (iterable of Prognostic): objects used to get
                tendencies for time stepping
            diagnostic_list (iterable of Diagnostic): objects used to get
                diagnostics before time stepping
            asselin_strength (float, optional): the filter parameter used to
                determine the strength of the Asselin filter. Default is 0.05.
            alpha (float, optional): constant from Williams (2009), where the
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
        state = state.copy()  # we don't want to modify the input state
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
        return state, new_state

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
        old_state (dict): model state at the last timestep
        state (dict): model state at the current timestep. May be modified by
            this function call, specifically by the Asselin filter.
        tendencies (dict): time derivatives at the current timestep in
            units/second
        timestep (timedelta): length of the timestep
        asselin_strength (float, optional): Asselin filter strength.
            Default is 0.05.
        alpha (float, optional): constant from Williams (2009), where the
            midpoint is shifted by alpha*influence, and the right point is
            shifted by (1-alpha)*influence. If alpha is 1 then the behavior
            is that of the classic Robert-Asselin time filter, while if it
            is 0.5 the filter will conserve the three-point mean.
            Default is 1.

    Returns:
        state (dict): the input state, modified in place
        new_state (dict): model state at the next timestep
    """
    new_state = {}
    for key in tendencies.keys():
        new_state[key] = (
            old_state[key] + 2*tendencies[key]*timestep.total_seconds())
        filter_influence = 0.5*asselin_strength*(
            old_state[key] - 2*state[key] + new_state[key])
        state[key] += alpha * filter_influence
        if alpha != 1.:
            new_state[key] += (1 - alpha) * filter_influence
    return state, new_state


def step_forward_euler(state, tendencies, timestep):
    return_state = {}
    for key in tendencies.keys():
        return_state[key] = state[key] + tendencies[key]*timestep.total_seconds()
    return return_state
