from sympl import DiagnosticComponent, initialize_numpy_arrays_with_properties
import logging
import numpy as np
try:
    from . import _dcmip
except ImportError:
    logging.warning(
        'Import Failed. DCMIP initial conditions will not be available.')


class DcmipInitialConditions(DiagnosticComponent):
    """
        Climt interface to the DCMIP initial conditions.
        Currently only provides interfaces to tests 4 and 5.
    """

    input_properties = {
        'latitude': {
            'dims': ['*'],
            'units': 'degrees_north',
        },
        'longitude': {
            'dims': ['*'],
            'units': 'degrees_east',
        },
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa',
        },
    }

    diagnostic_properties = {
        'eastward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'northward_wind': {
            'dims': ['mid_levels', '*'],
            'units': 'm s^-1',
        },
        'air_temperature': {
            'dims': ['mid_levels', '*'],
            'units': 'degK',
        },
        'surface_geopotential': {
            'dims': ['*'],
            'units': 'm^2 s^-2',
        },
        'surface_air_pressure': {
            'dims': ['*'],
            'units': 'Pa',
        },
        'specific_humidity': {
            'dims': ['mid_levels', '*'],
            'units': 'g/g',
        },
    }


    def __init__(self,
                 condition_type='baroclinic_wave',
                 add_perturbation=True,
                 moist=False,
                 **kwargs):
        """
        Initialize the DCMIP module.

        Args:
            condition_type (str, optional):
                The type of initial conditions desired. Can be
                one of :code:`'baroclinic_wave'` or
                :code:`'tropical_cyclone'`.
            add_perturbation (bool, optional):
                Whether a perturbation must be added. Only applies
                to the baroclinic wave test.
        """
        if condition_type not in ['baroclinic_wave', 'tropical_cyclone']:
            raise ValueError("type_of_output has to be one \
                             of 'baroclinic_wave' or 'tropical_cyclone'")
        if condition_type is 'tropical_cyclone' and moist is False:
            raise ValueError("moist_simulation must be True for tropical cyclone test")
        self._condition_type = condition_type
        self._add_perturbation = add_perturbation
        self._moist = moist
        super(DcmipInitialConditions, self).__init__(**kwargs)

    def array_call(self, state):

        longitude = np.radians(state['longitude'])
        latitude = np.radians(state['latitude'])

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        if self._condition_type is 'baroclinic_wave':
            u, v, t, q, p_surface, phi_surface = _dcmip.get_baroclinic_wave_ics(
                state['air_pressure'],
                longitude, latitude,
                perturb=self._add_perturbation,
                moist_sim=self._moist)
        elif self._condition_type is 'tropical_cyclone':
            u, v, t, q, p_surface, phi_surface = _dcmip.get_tropical_cyclone_ics(
                state['air_pressure'],
                longitude, latitude,
                perturb=self._add_perturbation,
                moist_sim=self._moist)

        diagnostics['eastward_wind'][:] = u
        diagnostics['northward_wind'][:] = v
        diagnostics['air_temperature'][:] = t
        diagnostics['surface_geopotential'][:] = phi_surface
        diagnostics['specific_humidity'][:] = q
        diagnostics['surface_air_pressure'][:] = p_surface

        return diagnostics
