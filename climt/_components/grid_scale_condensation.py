from sympl import Stepper, get_constant
from .._core import bolton_q_sat, bolton_dqsat_dT
import numpy as np


class GridScaleCondensation(Stepper):
    """
    Calculate condensation due to supersaturation of water.

    Condenses supersaturated water at the grid scale, assuming all
    condensed water falls as precipitation.
    """

    input_properties = {
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'kg/kg',
        },
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa'
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'Pa',
        }
    }

    diagnostic_properties = {
        'precipitation_amount': {
            'dims': ['*'],
            'units': 'kg m^-2',
        }
    }

    output_properties = {
        'air_temperature': {'units': 'degK'},
        'specific_humidity': {'units': 'kg/kg'},
    }

    def __init__(self, **kwargs):
        self._update_constants()
        super(GridScaleCondensation, self).__init__(**kwargs)

    def _update_constants(self):
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
        self._Lv = get_constant('latent_heat_of_condensation', 'J/kg')
        self._Rd = get_constant('gas_constant_of_dry_air', 'J/kg/degK')
        self._Rh2O = get_constant('gas_constant_of_vapor_phase', 'J/kg/degK')
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._rhow = get_constant('density_of_liquid_phase', 'kg/m^3')

    def array_call(self, raw_state, timestep):
        """
        Gets diagnostics from the current model state and steps the state
        forward in time according to the timestep.

        Args:
            raw_state (dict): A model state dictionary of numpy arrays
                satisfying input_properties. Will be updated with any
                diagnostic quantities produced by this object for the time of
                the input state.

        Returns:
            next_state (dict): A dictionary whose keys are strings indicating
                state quantities and values are the value of those quantities
                at the timestep after input state.

        Raises:
            KeyError: If a required quantity is missing from the state.
            InvalidStateException: If state is not a value input for the
                Stepper instance for other reasons.
        """
        self._update_constants()
        T = raw_state['air_temperature']
        q = raw_state['specific_humidity']
        p = raw_state['air_pressure']
        p_interface = raw_state['air_pressure_on_interface_levels']

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
        mass = (p_interface[1:, :] - p_interface[:-1, :])/(
            self._g*self._rhow)
        precipitation = np.sum(condensed_q * mass, axis=0)

        diagnostics = {
            'precipitation_amount': precipitation,
        }
        outputs = {
            'air_temperature': new_T,
            'specific_humidity': new_q,
        }

        return diagnostics, outputs
