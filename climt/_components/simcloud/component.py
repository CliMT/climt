import numpy as np
from collections import namedtuple
from sympl import DiagnosticComponent, Stepper, get_constant, initialize_numpy_arrays_with_properties, get_tracer_names
from scipy.special import lambertw

# Bundle of tunable parameters + physical constants passed to the framework-agnostic
# core. Keeping it immutable means the core has no hidden dependency on a component
# instance and can be unit-tested with plain numpy arrays.
SimCloudParams = namedtuple(
    'SimCloudParams',
    [
    'linear_a_surf',
    'linear_a_top',
    'linear_power',
    'p_ref',
    'qv_surface',
    'freezedry_power',
    'park_a',
    'park_b',
    'T_min',
    'T_max',
    'reff_liq',
    'reff_ice',
    'qcl_val',
    'do_add_stratocumulus',
    'g', 'Rd', 'Rv', 'Cpd', 'L', 'epsilon', 'kappa',
    'c_vv',   
    'c_vl',   
    'c_pv',   
    'E0v',    
    'T_trip',
    ])

def simcloud_core(air_pressure, interface_pressure, temperature, q_vapor,
                  surface_pressure, P):
    """
    Framework-agnostic SimCloud physics.

    All inputs are plain numpy arrays shaped (n_layer, n_profile), except
    ``surface_pressure`` which is shaped (n_profile,). ``P`` is a SimCloudParams.

    Returns grid-mean fields as a dict:
        'cloud_fraction'  : total cloud area fraction (0-1)
        'qcl'             : grid-mean cloud liquid mixing ratio (kg/kg)
        'qci'             : grid-mean cloud ice mixing ratio (kg/kg)
        'liquid_fraction' : liquid fraction of condensate (0-1)
    """
    air_pressure = np.asarray(air_pressure)
    temperature = np.asarray(temperature)
    q_vapor = np.asarray(q_vapor)
    surface_pressure = np.asarray(surface_pressure)

    n_layer, n_profile = air_pressure.shape
    surface_pressure_2d = surface_pressure[np.newaxis, :]

    # kelvin to celsius
    tc = temperature - 273.15
    tc_clip_w = np.maximum(tc, -100.0)
    tc_clip_i = np.maximum(tc, -100.0)

    es_w = np.exp(34.4942 - 4924.99 / (tc_clip_w + 237.1)) / (tc_clip_w + 105.0)**1.57
    es_i = np.exp(43.4942 - 6545.8 / (tc_clip_i + 278.0)) / (tc_clip_i + 868.0)**2.0

    liquid_fraction = np.clip((temperature - P.T_min) / (P.T_max - P.T_min), 0.0, 1.0)
    sat_vapor_pressure = liquid_fraction * es_w + (1.0 - liquid_fraction) * es_i
    sat_specific_humidity = P.epsilon * sat_vapor_pressure / (air_pressure - (1.0 - P.epsilon) * sat_vapor_pressure)
    relative_humidity = np.clip(q_vapor / sat_specific_humidity, 0.0, 1.0)

    # large scale cloud fraction
    exp_arg = np.clip(1.0 - (surface_pressure_2d / air_pressure)**P.linear_power, -50, 1.0)
    coeff_a = P.linear_a_top + (P.linear_a_surf - P.linear_a_top) * np.exp(exp_arg)
    large_scale_cloud_fraction = np.clip(coeff_a * (relative_humidity - 1.0) + 1.0, 0.0, 1.0)

    # freeze dry adjustment
    qv_crit = (air_pressure / surface_pressure_2d)**P.freezedry_power * P.qv_surface
    f_fd = np.clip(q_vapor / qv_crit, 0.15, 1.0)
    large_scale_cloud_fraction = large_scale_cloud_fraction * f_fd

    # marine stratus clouds
    stratocumulus_cloud_fraction = np.zeros_like(large_scale_cloud_fraction)
    if P.do_add_stratocumulus:
        # lcl
        # Safeguard surface temperature and relative humidity from unphysical inputs (like 0 K or negative RH)
        t_s = np.maximum(temperature[0, :], 200.0)
        rh_s = np.maximum(relative_humidity[0, :], 1e-4)
        rm_s = (1.0 - q_vapor[0, :]) * P.Rd + q_vapor[0, :] * P.Rv
        cpm_s = (1.0 - q_vapor[0, :]) * P.Cpd + q_vapor[0, :] * (P.Rv / P.epsilon)

        a_r = (cpm_s / rm_s) + (P.c_vl - P.c_pv) / P.Rv
        b_r = -(P.E0v - (P.c_vv - P.c_vl) * P.T_trip) / (P.Rv * t_s)
        c_r = b_r / a_r
        arg_w = (rh_s + 1e-10)**(1.0 / a_r) * c_r * np.exp(c_r)

        tlcl = (c_r / lambertw(arg_w, k=-1).real) * t_s
        zlcl = cpm_s * (t_s - tlcl) / P.g
        plcl = surface_pressure * (tlcl / t_s)**(cpm_s / rm_s)

        # elf proxies
        theta = temperature * (P.p_ref / air_pressure)**P.kappa
        k700 = np.argmin(np.abs(air_pressure - 0.7 * P.p_ref), axis=0)
        t700 = np.array([temperature[k700[i], i] for i in range(n_profile)])
        p700 = np.array([air_pressure[k700[i], i] for i in range(n_profile)])
        theta_700 = t700 * (P.p_ref / p700)**P.kappa
        lts = theta_700 - theta[0, :]

        z700 = - (P.Rd * 0.5 * (temperature[0, :] + t700) / P.g) * np.log(p700 / surface_pressure)

        def gamma_moist(T, pres):
            tc_m = T - 273.15
            tc_c = np.maximum(tc_m, -100)
            es_m = np.exp(34.4942 - 4924.99 / (tc_c + 237.1)) / (tc_c + 105)**1.57
            qs_m = P.epsilon * es_m / (pres - (1.0 - P.epsilon) * es_m)

            # as per Woods and Bretherton (2006); what is given here is the moist-adiabatic potential temperature gradien
            R = (1 + P.L*qs_m/(P.Rd*T)) / (1 + P.L**2*qs_m/(P.Cpd*P.Rv*T**2))
            return (P.g / P.Cpd) * (1.0 - R)

        g_lcl = gamma_moist(tlcl, plcl)
        g_700 = gamma_moist(t700, p700)
        zinv = np.clip(-lts/g_700 + z700 + 2750.0*(g_lcl/g_700), zlcl, zlcl + 2750.0)
        elf = f_fd[0, :] * (1.0 - np.sqrt(zinv * np.maximum(zlcl, 0)) / 2750.0)  # 2750 m is a constant scale height
        low_cld = np.clip(P.park_a * elf + P.park_b, 0.0, 1.0)

        dth = np.diff(theta, axis=0)
        dp = np.diff(air_pressure, axis=0)
        dth_dp = dth / dp  # potential temperature lapse rate

        for i in range(n_profile):
            p_col = air_pressure[:, i]
            mask = (p_col > 0.7 * P.p_ref)[:-1] & (dth_dp[:, i] * 100.0 < -0.05)
            if np.any(mask):
                idx = np.argmin(np.where(mask, dth_dp[:, i], 1e10))
                stratocumulus_cloud_fraction[idx, i] = low_cld[i]

    total_cloud_fraction = np.maximum(large_scale_cloud_fraction, stratocumulus_cloud_fraction)
    total_cloud_fraction = np.where(air_pressure < 10000.0, 0.0, total_cloud_fraction)

    qcl_val_kg = P.qcl_val / 1000.0
    in_cloud_water_mixing_ratio = np.minimum(
        qcl_val_kg,
        np.maximum(3e-7, qcl_val_kg * (temperature - 220.0) / 60.0)
    )
    grid_mean_condensate = in_cloud_water_mixing_ratio * total_cloud_fraction

    return {
        'cloud_fraction': total_cloud_fraction,
        'qcl': grid_mean_condensate * liquid_fraction,
        'qci': grid_mean_condensate * (1.0 - liquid_fraction),
        'liquid_fraction': liquid_fraction,
    }


class SimCloud(DiagnosticComponent):
    """
    SimCloud version 1.0: A simple diagnostic cloud scheme for idealized climate models.
    Based on Liu et al. (2021), GMD, and the ISCA reference implementation.
    LCL calculation follows the exact expression from Romps (2017), JAS.
    Saturation vapor pressure follows Huang (2018), JAMC.
    """

    input_properties = {
        'air_pressure': {'dims': ['mid_levels', '*'], 'units': 'Pa'},
        'air_pressure_on_interface_levels': {'dims': ['interface_levels', '*'], 'units': 'Pa'},
        'air_temperature': {'dims': ['mid_levels', '*'], 'units': 'degK'},
        'specific_humidity': {'dims': ['mid_levels', '*'], 'units': 'kg/kg'},
        'surface_air_pressure': {'dims': ['*'], 'units': 'Pa'},
    }

    diagnostic_properties = {
        'cloud_area_fraction_in_atmosphere_layer': {'dims': ['mid_levels', '*'], 'units': 'dimensionless'},
        'mass_content_of_cloud_liquid_water_in_atmosphere_layer': {'dims': ['mid_levels', '*'], 'units': 'kg m^-2'},
        'mass_content_of_cloud_ice_in_atmosphere_layer': {'dims': ['mid_levels', '*'], 'units': 'kg m^-2'},
        'cloud_water_droplet_radius': {'dims': ['mid_levels', '*'], 'units': 'micrometer'},
        'cloud_ice_particle_size': {'dims': ['mid_levels', '*'], 'units': 'micrometer'},
    }

    def __init__(
        self,
        linear_a_surf=42.0,
        linear_a_top=12.0,
        linear_power=8.5,
        p_ref=1e5,
        qv_surface=0.003,
        freezedry_power=2.5,
        park_a=1.272,
        park_b=-0.366,
        T_min=233.15,
        T_max=268.15,
        reff_liq=14.0,
        reff_ice=25.0,
        qcl_val=0.2,
        do_add_stratocumulus=True,
        **kwargs
    ):
        self.linear_a_surf = linear_a_surf
        self.linear_a_top = linear_a_top
        self.linear_power = linear_power
        self.p_ref = p_ref
        self.qv_surface = qv_surface
        self.freezedry_power = freezedry_power
        self.park_a = park_a
        self.park_b = park_b
        self.T_min = T_min
        self.T_max = T_max
        self.reff_liq = reff_liq
        self.reff_ice = reff_ice
        self.qcl_val = qcl_val
        self.do_add_stratocumulus = do_add_stratocumulus

        self._g = get_constant('gravitational_acceleration', 'm/s^2')
        self._Rd = get_constant('gas_constant_of_dry_air', 'J kg^-1 K^-1')
        self._Rv = get_constant('gas_constant_of_water_vapor', 'J kg^-1 K^-1')
        self._Cpd = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J/kg/K')
        self._L = get_constant('latent_heat_of_condensation', 'J/kg')
        self._epsilon = self._Rd / self._Rv
        self._kappa = self._Rd / self._Cpd

        self.params = SimCloudParams(
            linear_a_surf=linear_a_surf, 
            linear_a_top=linear_a_top,
            linear_power=linear_power, 
            p_ref=p_ref, 
            qv_surface=qv_surface,
            freezedry_power=freezedry_power, 
            park_a=park_a, 
            park_b=park_b,
            T_min=T_min, 
            T_max=T_max, 
            reff_liq=reff_liq, 
            reff_ice=reff_ice,
            qcl_val=qcl_val, 
            do_add_stratocumulus=do_add_stratocumulus,
            g=self._g, Rd=self._Rd, Rv=self._Rv, 
            Cpd=self._Cpd, L=self._L,
            epsilon=self._epsilon, 
            kappa=self._kappa,
            c_vv=1418.0,      
            c_vl=4119.0,      
            c_pv=1879.0,      
            E0v=2.3740e6,     
            T_trip=273.16,   
        )

        super(SimCloud, self).__init__(**kwargs)

    def array_call(self, state):
        air_pressure = np.asarray(state['air_pressure'])
        interface_pressure = np.asarray(state['air_pressure_on_interface_levels'])
        temperature = np.asarray(state['air_temperature'])
        q_vapor = np.asarray(state['specific_humidity'])
        surface_pressure = np.asarray(state['surface_air_pressure'])

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, state, self.input_properties
        )

        core = simcloud_core(air_pressure, interface_pressure, temperature,
                             q_vapor, surface_pressure, self.params)

        # grid-mean condensate (kg/kg) -> column-integrated path (kg m^-2)
        dp = np.abs(np.diff(interface_pressure, axis=0))
        diagnostics['cloud_area_fraction_in_atmosphere_layer'] = core['cloud_fraction']
        diagnostics['mass_content_of_cloud_liquid_water_in_atmosphere_layer'] = core['qcl'] * dp / self._g
        diagnostics['mass_content_of_cloud_ice_in_atmosphere_layer'] = core['qci'] * dp / self._g
        diagnostics['cloud_water_droplet_radius'] = self.reff_liq * np.ones_like(air_pressure)
        diagnostics['cloud_ice_particle_size'] = self.reff_ice * np.ones_like(air_pressure)

        return diagnostics

class SimCloudCondensation(Stepper):
    """
    Prognostic-diagnostic version of SimCloud that conserves total water (vapor + liquid + ice)
    and includes latent heating.
    """

    input_properties = dict(SimCloud.input_properties)

    output_properties = {
        'air_temperature': {'units': 'degK'},
        'specific_humidity': {'units': 'kg/kg'}
    }

    diagnostic_properties = {
        **SimCloud.diagnostic_properties,
        'precipitation_amount': {'dims': ['*'], 'units': 'kg m^-2'},
    }

    def __init__(self, simcloud_instance=None, **kwargs):
        self._simcloud = simcloud_instance or SimCloud(**kwargs)

        self.params = self._simcloud.params
        self._g = self.params.g
        self._Cpd = self.params.Cpd
        self._L = self.params.L

        try:
            self._Lf = get_constant('latent_heat_of_fusion', 'J/kg')
        except Exception:
            self._Lf = 333550.0

        try:
            self._rhow = get_constant('density_of_liquid_phase', 'kg/m^3')
        except Exception:
            self._rhow = 1000.0

        self.input_properties = self.input_properties.copy()
        self.output_properties = self.output_properties.copy()

        tracer_names = get_tracer_names()
        self.use_prognostic_tracers = False
        if 'cloud_liquid_water_mixing_ratio' in tracer_names and 'cloud_ice_mixing_ratio' in tracer_names:
            self.use_prognostic_tracers = True
            self.input_properties['cloud_liquid_water_mixing_ratio'] = {'dims': ['mid_levels', '*'], 'units': 'kg/kg'}
            self.input_properties['cloud_ice_mixing_ratio'] = {'dims': ['mid_levels', '*'], 'units': 'kg/kg'}
            self.output_properties['cloud_liquid_water_mixing_ratio'] = {'units': 'kg/kg'}
            self.output_properties['cloud_ice_mixing_ratio'] = {'units': 'kg/kg'}

        super(SimCloudCondensation, self).__init__(**kwargs)

    def array_call(self, state, timestep):
        T = state['air_temperature']
        q_vapor = state['specific_humidity']
        p = np.asarray(state['air_pressure'])
        p_interface = np.asarray(state['air_pressure_on_interface_levels'])
        surface_pressure = state['surface_air_pressure']

        dp = np.abs(np.diff(p_interface, axis=0))
        g = self._g

        if self.use_prognostic_tracers:
            q_water = state['cloud_liquid_water_mixing_ratio']
            q_ice = state['cloud_ice_mixing_ratio']
        else:
            q_water = np.zeros_like(q_vapor)
            q_ice = np.zeros_like(q_vapor)

        q_total = q_vapor + q_water + q_ice

        # diagnose clouds against total water (vapor + condensate)
        core = simcloud_core(p, p_interface, T, q_total, surface_pressure, self.params)

        q_water_diag = core['qcl']
        q_ice_diag = core['qci']

        # column-integrated paths reported as diagnostics (before any clipping)
        liq_path = q_water_diag * dp / g
        ice_path = q_ice_diag * dp / g

        # Ensure specific humidity does not go negative by clipping and scaling cloud water/ice
        q_condensed = q_water_diag + q_ice_diag
        insufficient_water = q_total < q_condensed

        scale = np.where(insufficient_water, q_total / (q_condensed + 1e-20), 1.0)
        q_water_diag = q_water_diag * scale
        q_ice_diag = q_ice_diag * scale
        q_vapor_new = np.maximum(0.0, q_total - q_water_diag - q_ice_diag)

        dq_water = q_water_diag - q_water
        dq_ice = q_ice_diag - q_ice

        # heating rate (J/kg / Cpd -> degK)
        dT = (dq_water * self._L + dq_ice * (self._L + self._Lf)) / self._Cpd
        T_new = T + dT

        # precipitation
        if not self.use_prognostic_tracers:
            precipitation = np.sum((q_water_diag + q_ice_diag) * dp / g, axis=0)
        else:
            precipitation = np.zeros(q_vapor.shape[1:])

        diagnostics = {
            'cloud_area_fraction_in_atmosphere_layer': core['cloud_fraction'],
            'mass_content_of_cloud_liquid_water_in_atmosphere_layer': liq_path,
            'mass_content_of_cloud_ice_in_atmosphere_layer': ice_path,
            'cloud_water_droplet_radius': self.params.reff_liq * np.ones_like(p),
            'cloud_ice_particle_size': self.params.reff_ice * np.ones_like(p),
            'precipitation_amount': precipitation
        }

        outputs = {
            'air_temperature': T_new,
            'specific_humidity': q_vapor_new
        }

        if self.use_prognostic_tracers:
            outputs['cloud_liquid_water_mixing_ratio'] = q_water_diag
            outputs['cloud_ice_mixing_ratio'] = q_ice_diag

        return diagnostics, outputs
