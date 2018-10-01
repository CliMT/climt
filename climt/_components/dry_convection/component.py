from sympl import (
    get_constant, Stepper, initialize_numpy_arrays_with_properties
)
import logging
try:
    from ._dry_convection import calculate_dry_adjustment
except ImportError:
    logging.warning(
        'Import failed. Dry Convection is likely not compiled and will not be '
        'available.'
    )


class DryConvectiveAdjustment(Stepper):
    """
    A conservative scheme to keep the temperature profile close to the
    dry adiabat if it is super-adiabatic.
    """
    input_properties = {
        'air_temperature': {
            'units': 'degK',
            'dims': ['*', 'mid_levels'],
        },
        'air_pressure': {
            'units': 'Pa',
            'dims': ['*', 'mid_levels'],
        },
        'air_pressure_on_interface_levels': {
            'units': 'Pa',
            'dims': ['*', 'interface_levels'],
            'alias': 'P_int',
        },
        'specific_humidity': {
            'units': 'kg/kg',
            'dims': ['*', 'mid_levels'],
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

        output_temp = output_arrays['air_temperature']
        output_temp[:] = state['air_temperature']

        output_q = output_arrays['specific_humidity']
        output_q[:] = q

        rd_cp = gas_constant(q, self._Rdair, self._Rv)/heat_capacity(q, self._Cpd, self._Cvap)
        theta = state['air_temperature']*(self._Pref/state['air_pressure'])**rd_cp
        theta_q = theta*(1 + output_q*self._Rv/self._Rdair - output_q)

        # print('initial theta:', theta)

        calculate_dry_adjustment(
            state['P_int'], state['air_pressure'],
            theta_q, output_q, output_temp,
            self._Cpd, self._Cvap, self._Rdair, self._Rv,
            self._Pref, q.shape[0], q.shape[1])

        return {}, output_arrays


def heat_capacity(q, dry_constant, moist_constant):

    return dry_constant*(1 - q) + moist_constant*q


def gas_constant(q, dry_constant, moist_constant):

    return dry_constant*(1 - q) + moist_constant*q
