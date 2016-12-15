from .._core.base_components import Implicit
from .._core.array import DataArray
from .._core.util import jit, replace_none_with_default
from .._core.util import get_3d_numpy_array as to_3d_array
import numpy as np


@jit(nopython=True)
def bolton_q_sat(T, p, Rd, Rh20):
    es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    epsilon = Rd/Rh20
    return epsilon*es/(p - (1 - epsilon)*es)


@jit(nopython=True)
def bolton_dqsat_dT(T, Lv, Rh20, q_sat):
    """Uses the assumptions of equation 12 in Reed and Jablonowski, 2012. In
    particular, assumes d(qsat)/dT is approximately epsilon/p*d(es)/dT"""
    return Lv*q_sat/(Rh20*T**2)


class GridScaleCondensation(Implicit):
    """Condenses supersaturated water at the grid scale, assuming all
    condensed water falls as precipitation."""

    inputs = (
        'air_temperature', 'specific_humidity', 'air_pressure', 'air_mass',
    )
    diagnostic_outputs = (
        'column_integrated_precipitation_rate',
    )
    tendency_outputs = (
        'air_temperature', 'specific_humidity',
    )

    def __init__(self,
                 gas_constant_of_dry_air=None,
                 gas_constant_of_water_vapor=None,
                 heat_capacity_of_dry_air_at_constant_pressure=None,
                 latent_heat_of_vaporization_of_water=None,):
        """

        Args:
            gas_constant_of_dry_air: Value in $J kg^{-1} K^{-1}$.
                Default taken from climt.default_constants.
            gas_constant_of_water_vapor: Value in $J kg^{-1} K^{-1}$
                Default taken from climt.default_constants.
            heat_capacity_of_dry_air_at_constant_pressure: Value in
                $J kg^{-1} K^{-1}$
                Default taken from climt.default_constants.
            latent_heat_of_vaporization_of_water: Value in $J kg^{-1}$.
                Default taken from climt.default_constants.
        """

        self._Cpd = replace_none_with_default(
            'heat_capacity_of_dry_air_at_constant_pressure',
            heat_capacity_of_dry_air_at_constant_pressure)
        self._Lv = replace_none_with_default(
            'latent_heat_of_vaporization_of_water',
            latent_heat_of_vaporization_of_water)
        self._Rd = replace_none_with_default('gas_constant_of_dry_air',
                                             gas_constant_of_dry_air)
        self._Rh20 = replace_none_with_default('gas_constant_of_water_vapor',
                                               gas_constant_of_water_vapor)
        self._q_sat = bolton_q_sat
        self._dqsat_dT = bolton_dqsat_dT

    def __call__(self, state, time_step):
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
        T = to_3d_array(state['air_temperature'].to_units('degK'))
        q = to_3d_array(state['specific_humidity'].to_units('kg/kg'))
        p = to_3d_array(state['air_pressure'].to_units('Pa'))
        mass = to_3d_array(state['air_mass'].to_units('kg'))
        q_sat = self._q_sat(T, p, self._Rd, self._Rh20)
        saturated = q > q_sat
        dqsat_dT = self._dqsat_dT(
            T[saturated], self._Lv, self._Rh20, q_sat[saturated])
        condensed_q = np.zeros_like(q)
        condensed_q[saturated] = (
            q[saturated] - q_sat[saturated])/(
            1 + self._Lv/self._Cpd * dqsat_dT)
        new_q = q.copy()
        new_T = T.copy()
        new_q[saturated] -= condensed_q[saturated]
        new_T[saturated] += self._Lv/self._Cpd * condensed_q[saturated]
        precipitation = np.sum(condensed_q * mass, axis=2)
        state['column_integrated_precipitation_rate'] = DataArray(
            precipitation/time_step.total_seconds(),
            dims=('x', 'y'), attrs={'units': 'kg/s'})
        new_state = {
            'air_temperature': DataArray(
                new_T, dims=state['air_temperature'].dims,
                attrs=state['air_temperature'].attrs),
            'specific_humidity': DataArray(
                new_q, dims=state['specific_humidity'].dims,
                attrs=state['specific_humidity'].attrs),
        }
        return new_state
