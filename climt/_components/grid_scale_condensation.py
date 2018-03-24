from sympl import (
    Implicit, DataArray,
    get_numpy_array, combine_dimensions)
from .._core import bolton_q_sat, bolton_dqsat_dT, get_constant
import numpy as np


class GridScaleCondensation(Implicit):
    """
    Calculate condensation due to supersaturation of water.

    Condenses supersaturated water at the grid scale, assuming all
    condensed water falls as precipitation.
    """

    inputs = (
        'air_temperature', 'specific_humidity', 'air_pressure',
        'air_pressure_on_interface_levels',
    )
    diagnostics = (
        'precipitation_amount',
    )
    outputs = (
        'air_temperature', 'specific_humidity',
    )

    def __init__(self):
        """
        Initialise component.
        """

        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
        self._Lv = get_constant('latent_heat_of_condensation', 'J/kg')
        self._Rd = get_constant('gas_constant_of_dry_air', 'J/kg/degK')
        self._Rh2O = get_constant('gas_constant_of_vapor_phase', 'J/kg/degK')
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._rhow = get_constant('density_of_liquid_phase', 'kg/m^3')

    def __call__(self, state, timestep):
        """
        Gets diagnostics from the current model state and steps the state
        forward in time according to the timestep.

        Args:
            state (dict): A model state dictionary. Will be updated with any
                diagnostic quantities produced by this object for the time of
                the input state.

        Returns:
            next_state (dict): A dictionary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the timestep after input state.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a valid input for the
                Implicit instance for other reasons.
        """
        T = get_numpy_array(
            state['air_temperature'].to_units('degK'),
            out_dims=('x', 'y', 'z'))
        q = get_numpy_array(
            state['specific_humidity'].to_units('kg/kg'),
            out_dims=('x', 'y', 'z'))
        p = get_numpy_array(
            state['air_pressure'].to_units('Pa'),
            out_dims=('x', 'y', 'z'))
        p_interface = get_numpy_array(
            state['air_pressure_on_interface_levels'].to_units('Pa'),
            out_dims=('x', 'y', 'z'))

        q_sat = bolton_q_sat(T, p, self._Rd, self._Rh2O)
        saturated = q > q_sat
        dqsat_dT = bolton_dqsat_dT(
            T[saturated], self._Lv, self._Rh2O, q_sat[saturated])

        condensed_q = np.zeros_like(q)
        condensed_q[saturated] = (
            q[saturated] - q_sat[saturated])/(
            1 + self._Lv/self._Cpd * dqsat_dT)
        new_q = q.copy()
        new_T = T.copy()
        new_q[saturated] -= condensed_q[saturated]
        new_T[saturated] += self._Lv/self._Cpd * condensed_q[saturated]
        mass = (p_interface[:, :, 1:] - p_interface[:, :, :-1])/(
            self._g*self._rhow)
        precipitation = np.sum(condensed_q * mass, axis=2)

        dims_3d = combine_dimensions(
            [state['air_temperature'], state['specific_humidity'],
             state['air_pressure']],
            out_dims=('x', 'y', 'z'))
        dims_2d = dims_3d[:-1]
        diagnostics = {
            'column_integrated_precipitation_rate': DataArray(
                precipitation / timestep.total_seconds(),
                dims=dims_2d, attrs={'units': 'kg/s'}).squeeze()
        }
        new_state = {
            'air_temperature': DataArray(
                new_T, dims=dims_3d,
                attrs=state['air_temperature'].attrs).squeeze(),
            'specific_humidity': DataArray(
                new_q, dims=dims_3d,
                attrs=state['specific_humidity'].attrs).squeeze(),
        }
        return new_state, diagnostics
