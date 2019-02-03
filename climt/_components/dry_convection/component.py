from sympl import (
    get_constant, Stepper, initialize_numpy_arrays_with_properties
)
import numpy as np


class DryConvectiveAdjustment(Stepper):
    """
    A conservative scheme to keep the temperature profile close to the
    dry adiabat if it is super-adiabatic.
    """
    input_properties = {
        'air_temperature': {
            'units': 'degK',
            'dims': ['mid_levels', '*'],
        },
        'air_pressure': {
            'units': 'Pa',
            'dims': ['mid_levels', '*'],
        },
        'air_pressure_on_interface_levels': {
            'units': 'Pa',
            'dims': ['interface_levels', '*'],
            'alias': 'P_int',
        },
        'specific_humidity': {
            'units': 'kg/kg',
            'dims': ['mid_levels', '*'],
        },
    }

    output_properties = {
        'air_temperature': {
            'units': 'degK',
        },
        'specific_humidity': {
            'units': 'kg/kg',
        },
    }

    diagnostic_properties = {}

    def array_call(self, state, time_step):

        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
        self._Cvap = get_constant('heat_capacity_of_vapor_phase', 'J/kg/K')
        self._Rdair = get_constant('gas_constant_of_dry_air', 'J/kg/degK')
        self._Pref = get_constant('reference_air_pressure', 'Pa')
        self._Rv = get_constant('gas_constant_of_vapor_phase', 'J/kg/K')

        q = state['specific_humidity']

        output_arrays = initialize_numpy_arrays_with_properties(
            self.output_properties, state, self.input_properties
        )

        output_temperature = output_arrays['air_temperature']
        output_temperature[:] = state['air_temperature']

        output_q = output_arrays['specific_humidity']
        output_q[:] = q

        rd_cp = self.gas_constant(q)/self.heat_capacity(q)
        theta = state['air_temperature']*(self._Pref/state['air_pressure'])**rd_cp

        num_levels = q.shape[0]

        pdiff = state['P_int'][:-1, :] - state['P_int'][1:, :]
        theta_q = theta*(1 + output_q*self._Rv/self._Rdair - output_q)

        for column in range(q.shape[-1]):
            for level in range(num_levels-1, -1, -1):

                dp = pdiff[:, column]
                theta_sum = np.cumsum(theta_q[level::, column])
                divisor = np.arange(1, num_levels - level+1)

                theta_avg = (theta_sum/divisor)[1::]

                theta_lesser = (theta_avg > theta_q[level+1::, column])
                if np.sum(theta_lesser) == 0:
                    continue

                convect_to_level = len(theta_lesser) - np.argmax(theta_lesser[::-1])

                if level == 0:
                    convect_to_level = max(convect_to_level, 1)

                if convect_to_level == 0:
                    continue
                stable_level = level + convect_to_level

                q_conv = output_q[level:stable_level, column]
                t_conv = output_temperature[level:stable_level, column]
                dp_conv = dp[level:stable_level]
                p_conv_high = state['P_int'][level, column]
                p_conv_low = state['P_int'][stable_level, column]

                enthalpy = self.heat_capacity(q_conv)*t_conv
                integral_enthalpy = np.sum(enthalpy*dp_conv)
                mean_conv_q = np.sum(q_conv*dp_conv)/(p_conv_high - p_conv_low)

                output_q[level:stable_level, column] = mean_conv_q

                rdcp_conv = self.gas_constant(mean_conv_q)/self.heat_capacity(mean_conv_q)

                theta_coeff = (
                    state['air_pressure'][level:stable_level, column]/self._Pref)**rdcp_conv

                integral_theta_den = np.sum(self.heat_capacity(q_conv)*theta_coeff*dp_conv)

                mean_theta = integral_enthalpy/integral_theta_den

                output_temperature[level:stable_level, column] = mean_theta*theta_coeff

        return {}, output_arrays

    def heat_capacity(self, q):
        """
        Calculate heat capacity based on amount of q
        """

        return self._Cpd*(1-q) + self._Cvap*q

    def gas_constant(self, q):
        """
        Calculate gas constant based on amount of q
        """

        return self._Rdair*(1-q) + self._Rv*q
