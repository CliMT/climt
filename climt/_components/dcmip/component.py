from sympl import (
    DiagnosticComponent, initialize_numpy_arrays_with_properties,
    get_constant)
import logging
import numpy as np
try:
    from . import _dcmip
except ImportError as error:
    logging.warning(
        'Import Failed. DCMIP initial conditions will not be available.')
    print(error)


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
        'atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'dimensionless',
            'alias': 'ak'
        },
        'atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'dimensionless',
            'alias': 'bk'
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
        'air_pressure': {
            'dims': ['mid_levels', '*'],
            'units': 'Pa',
        },
        'air_pressure_on_interface_levels': {
            'dims': ['interface_levels', '*'],
            'units': 'Pa',
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
        if condition_type == 'tropical_cyclone' and moist is False:
            raise ValueError("moist_simulation must be True for tropical cyclone test")
        self._condition_type = condition_type
        self._add_perturbation = add_perturbation
        self._moist = moist
        super(DcmipInitialConditions, self).__init__(**kwargs)

    def array_call(self, state):
        toa_pressure = get_constant('top_of_model_pressure', 'Pa')
        rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1')

        longitude = np.radians(state['longitude'])
        latitude = np.radians(state['latitude'])

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        if self._condition_type == 'baroclinic_wave':
            u, v, t, q, p_surface, phi_surface = _dcmip.get_baroclinic_wave_ics(
                state['air_pressure'],
                longitude, latitude,
                perturb=self._add_perturbation,
                moist_sim=self._moist)
        elif self._condition_type == 'tropical_cyclone':
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
        p_interface = (
            state['ak'] + state['bk']*(p_surface - toa_pressure))
        delta_p = p_interface[1:, :] - p_interface[:-1, :]
        rk = rd/cpd

        diagnostics['air_pressure_on_interface_levels'][:] = p_interface
        diagnostics['air_pressure'][:] = (
            (p_interface[1:, :]**(rk+1) - p_interface[:-1, :]**(rk+1)) / (
                (rk+1) * delta_p
            )
        ) ** (1./rk)

        return diagnostics
