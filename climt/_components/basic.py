from .._core.base_components import Prognostic
from .._core.array import DataArray
from .._core.units import unit_registry as ureg


class ConstantPrognostic(Prognostic):

    def __init__(self, tendencies, diagnostics=None):
        self._tendencies = tendencies
        if diagnostics is not None:
            self._diagnostics = diagnostics
        else:
            self._diagnostics = {}

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
        """
        return self._tendencies.copy(), self._diagnostics.copy()


class RelaxationPrognostic(Prognostic):

    def __init__(self, quantity_name, equilibrium_value=None,
                 relaxation_timescale=None):
        self._quantity_name = quantity_name
        self._equilibrium_value = equilibrium_value
        self._tau = relaxation_timescale

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
        """
        value = state[self._quantity_name].values
        if self._equilibrium_value is None:
            equilibrium = state['equilibrium_' + self._quantity_name].to_units(
                value.attrs['units']).values
        else:
            equilibrium = self._equilibrium_value.to_units(
                value.attrs['units']).values
        if self._tau is None:
            tau = state[
                self._quantity_name + '_relaxation_timescale'].to_units(
                's').values
        else:
            tau = self._tau.to_units('s').values
        tendency_unit_string = str(
            ureg(state[self._quantity_name].attrs['units']) / ureg('s'))
        tendencies = {
            self._quantity_name: DataArray(
                (equilibrium - value)/tau,
                dims=value.dims,
                attrs={'units': tendency_unit_string}
            )
        }
        return tendencies, {}
