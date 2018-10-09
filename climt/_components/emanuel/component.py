# -*- coding: utf-8 -*-
from ..._core import bolton_q_sat, ensure_contiguous_state
from sympl import (
    ImplicitTendencyComponent, get_constant, initialize_numpy_arrays_with_properties)
import numpy as np
import logging
try:
    from . import _emanuel_convection
except ImportError as error:
    logging.warning(
        'Import failed. Emanuel Convection is likely not compiled and will not '
        'be available.'
    )
    print(error)


class EmanuelConvection(ImplicitTendencyComponent):
    """
    The Emanuel convection scheme from `[Emanuel and Zivkovic-Rothman]`_

    .. _[Emanuel and Zivkovic-Rothman]:
            http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1999)056%3C1766%3ADAEOAC%3E2.0.CO%3B2

    """

    input_properties = {
        'air_temperature': {
            'dims': ['*', 'mid_levels'],
            'units': 'degK',
        },
        'specific_humidity': {
            'dims': ['*', 'mid_levels'],
            'units': 'kg/kg',
        },
        'eastward_wind': {
            'dims': ['*', 'mid_levels'],
            'units': 'm s^-1',
        },
        'northward_wind': {
            'dims': ['*', 'mid_levels'],
            'units': 'm s^-1',
        },
        'air_pressure': {
            'dims': ['*', 'mid_levels'],
            'units': 'mbar',
        },
        'air_pressure_on_interface_levels': {
            'dims': ['*', 'interface_levels'],
            'units': 'mbar',
        },
        'cloud_base_mass_flux': {
            'dims': ['*'],
            'units': 'kg m^-2 s^-1',
        },
    }

    diagnostic_properties = {
        'convective_state': {
            'dims': ['*'],
            'units': 'dimensionless',
            'dtype': np.int32,
        },
        'convective_precipitation_rate': {
            'dims': ['*'],
            'units': 'mm day^-1',
        },
        'convective_downdraft_velocity_scale': {
            'dims': ['*'],
            'units': 'm s^-1',
        },
        'convective_downdraft_temperature_scale': {
            'dims': ['*'],
            'units': 'degK',
        },
        'convective_downdraft_specific_humidity_scale': {
            'dims': ['*'],
            'units': 'kg/kg',
        },
        'cloud_base_mass_flux': {
            'dims': ['*'],
            'units': 'kg m^-2 s^-1',
        },
        'atmosphere_convective_available_potential_energy': {
            'dims': ['*'],
            'units': 'J kg^-1',
        },
        'air_temperature_tendency_from_convection': {
            'dims': ['*', 'mid_levels'],
            'units': 'degK day^-1',
        }
    }

    tendency_properties = {
        'air_temperature': {'units': 'degK s^-1'},
        'specific_humidity': {'units': 'kg/kg s^-1'},
        'eastward_wind': {'units': 'm s^-2'},
        'northward_wind': {'units': 'm s^-2'},
    }

    def __init__(self,
                 minimum_convecting_layer=1,
                 autoconversion_water_content_threshold=0.0011,
                 autoconversion_temperature_threshold=-55,
                 entrainment_mixing_coefficient=1.5,
                 downdraft_area_fraction=0.05,
                 precipitation_fraction_outside_cloud=0.12,
                 speed_water_droplets=50.0,
                 speed_snow=5.5,
                 rain_evaporation_coefficient=1.0,
                 snow_evaporation_coefficient=0.8,
                 convective_momentum_transfer_coefficient=0.7,
                 downdraft_surface_velocity_coefficient=10.0,
                 convection_bouyancy_threshold=0.9,
                 mass_flux_relaxation_rate=0.1,
                 mass_flux_damping_rate=0.1,
                 reference_mass_flux_timescale=300.,
                 **kwargs):
        """

        Args:

            minimum_convecting_layer (int, optional):
                The least model level from which convection can be initiated. Normally set
                to :code:`1` if using bulk PBL schemes. Else, it should be set to the first
                model level at which the temperature is defined.

            autoconversion_water_content_threshold (float, optional):
                The amount of water vapour in :math:`kg/kg`
                above which condensation occurs in warm (above freezing point of water)
                clouds. This value linearly reduces to zero between the freezing point and
                the :code:`autoconversion_temperature_threshold`.

            autoconversion_temperature_threshold (float, optional):
                The temperature in :math:`^\circ C` below which
                all water vapour is converted to rain/snow.

            entrainment_mixing_coefficient (float, optional):
                The coefficient of mixing for entrainment of environmental air into
                the cloud.

            downdraft_area_fraction (float, optional):
                The fractional area covered by unsaturated downdrafts.

            precipitation_fraction_outside_cloud (float, optional):
                The fraction of precipitation falling outside the cloud.

            speed_water_droplets (float, optional):
                The speed of descent of water droplets in :math:`Pa/s`.

            speed_snow (float, optional):
                The speed of descent of snow in :math:`Pa/s`.

            rain_evaporation_coefficient (float, optional):
                Coefficient governing the rate of evaporation of rain.

            snow_evaporation_coefficient (float, optional):
                Coefficient governing the rate of evaporation of snow.

            convective_momentum_transfer_coefficient (float, optional):
                Coefficient **between 0 and 1** governing momentum transport
                by clouds. A value of 1 **shuts off** momentum transport.

            downdraft_surface_velocity_coefficient (float, optional):
                Coefficient mulitplying the downdraft mass flux to calculate
                the downdraft velocity scale.

            convection_bouyancy_threshold (float, optional):
                The maximum negative temperature perturbation in :math:`degK` a parcel can
                have below the temperature at its level of free convection.
                If difference is greater, and previous cloud base mass flux
                is zero, there is no convection.

            mass_flux_relaxation_rate (float, optional):
                Coefficient governing the rate of relaxation to subcloud-layer
                quasi-equilibrium.

            mass_flux_damping_rate (float, optional):
                Coefficient which damps the currently calculated mass flux towards
                the value from the previous time step.

            reference_mass_flux_timescale (float, optional):
                Timescale used to calculate the actual damping coefficient along with
                :code:`mass_flux_damping_rate` and the current time step.

        """

        if (convective_momentum_transfer_coefficient < 0 or
                convective_momentum_transfer_coefficient > 1):
            raise ValueError("Momentum transfer coefficient must be between 0 and 1.")

        if (downdraft_area_fraction < 0 or
                downdraft_area_fraction > 1):
            raise ValueError("Downdraft fraction must be between 0 and 1.")

        if (precipitation_fraction_outside_cloud < 0 or
                precipitation_fraction_outside_cloud > 1):
            raise ValueError("Outside cloud precipitation fraction must be between 0 and 1.")

        self._con_mom_txfr = convective_momentum_transfer_coefficient
        self._downdraft_area_frac = downdraft_area_fraction
        self._precip_frac_outside_cloud = precipitation_fraction_outside_cloud
        self._min_conv_layer = minimum_convecting_layer
        self._crit_humidity = autoconversion_water_content_threshold
        self._crit_temp = autoconversion_temperature_threshold
        self._entrain_coeff = entrainment_mixing_coefficient
        self._droplet_speed = speed_water_droplets
        self._snow_speed = speed_snow
        self._rain_evap = rain_evaporation_coefficient
        self._snow_evap = snow_evaporation_coefficient
        self._beta = downdraft_surface_velocity_coefficient
        self._dtmax = convection_bouyancy_threshold
        self._mf_damp = mass_flux_damping_rate
        self._alpha = mass_flux_relaxation_rate
        self._mf_timescale = reference_mass_flux_timescale
        self._ntracers = 0
        self._set_fortran_constants()
        super(EmanuelConvection, self).__init__(**kwargs)

    def _set_fortran_constants(self):
        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/degK')
        self._Cpv = get_constant('heat_capacity_of_vapor_phase', 'J/kg/degK')
        self._Rdair = get_constant('gas_constant_of_dry_air', 'J/kg/degK')
        self._Rcond = get_constant('gas_constant_of_vapor_phase', 'J/kg/degK')
        self._Lv = get_constant('latent_heat_of_condensation', 'J/kg')
        self._rho_condensible = get_constant('density_of_liquid_phase', 'kg/m^3')
        self._Cl = get_constant('specific_enthalpy_of_vapor_phase', 'J/kg')
        _emanuel_convection.init_emanuel_convection(
            self._min_conv_layer,
            self._crit_humidity,
            self._crit_temp, self._entrain_coeff,
            self._downdraft_area_frac,
            self._precip_frac_outside_cloud,
            self._droplet_speed, self._snow_speed,
            self._rain_evap, self._snow_evap,
            self._con_mom_txfr, self._dtmax,
            self._beta, self._alpha, self._mf_damp,
            self._Cpd, self._Cpv, self._Cl,
            self._Rcond, self._Rdair,
            self._Lv, self._g,
            self._rho_condensible, self._mf_timescale)

    @ensure_contiguous_state
    def array_call(self, raw_state, timestep):
        """
        Get convective heating and moistening.

        Args:

            raw_state (dict):
                The state dictionary of numpy arrays satisfying this
                component's input properties.

        Returns:

            tendencies (dict), diagnostics (dict):
                * The heating and moistening tendencies
                * Any diagnostics associated.

        """
        self._set_fortran_constants()

        num_cols, num_levs = raw_state['air_temperature'].shape

        max_conv_level = num_levs - 3

        tendencies = initialize_numpy_arrays_with_properties(
            self.tendency_properties, raw_state, self.input_properties
        )
        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )

        q_sat = bolton_q_sat(
            raw_state['air_temperature'],
            raw_state['air_pressure'] * 100,
            self._Cpd, self._Cpv
        )

        _emanuel_convection.convect(
            num_levs,
            num_cols,
            max_conv_level,
            self._ntracers,
            timestep.total_seconds(),
            raw_state['air_temperature'],
            raw_state['specific_humidity'],
            q_sat,
            raw_state['eastward_wind'],
            raw_state['northward_wind'],
            raw_state['air_pressure'],
            raw_state['air_pressure_on_interface_levels'],
            diagnostics['convective_state'],
            diagnostics['convective_precipitation_rate'],
            diagnostics['convective_downdraft_velocity_scale'],
            diagnostics['convective_downdraft_temperature_scale'],
            diagnostics['convective_downdraft_specific_humidity_scale'],
            raw_state['cloud_base_mass_flux'],
            diagnostics['atmosphere_convective_available_potential_energy'],
            tendencies['air_temperature'],
            tendencies['specific_humidity'],
            tendencies['eastward_wind'],
            tendencies['northward_wind'])

        diagnostics['air_temperature_tendency_from_convection'][:] = (
            tendencies['air_temperature'] * 86400.)
        diagnostics['cloud_base_mass_flux'][:] = raw_state['cloud_base_mass_flux']
        return tendencies, diagnostics
