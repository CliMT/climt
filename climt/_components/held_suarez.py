# -*- coding: utf-8 -*-
from sympl import TendencyComponent, get_constant
import numpy as np
from typing import NamedTuple
from .._core.backend import jit_compile, prange

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x: x

class HeldSuarezParams(NamedTuple):
    sigma_b: float; k_f: float; k_a: float; k_s: float; delta_T_y: float; delta_theta_z: float; p0: float; kappa: float

class HeldSuarez(TendencyComponent):
    input_properties = {
        "eastward_wind": {"dims": ["*", "mid_levels"], "units": "m s^-1"},
        "northward_wind": {"dims": ["*", "mid_levels"], "units": "m s^-1"},
        "air_temperature": {"dims": ["*", "mid_levels"], "units": "degK"},
        "air_pressure": {"dims": ["*", "mid_levels"], "units": "Pa"},
        "surface_air_pressure": {"dims": ["*"], "units": "Pa"},
        "latitude": {"dims": ["*"], "units": "degrees_north"},
    }
    tendency_properties = {
        "eastward_wind": {"units": "m s^-2"},
        "northward_wind": {"units": "m s^-2"},
        "air_temperature": {"units": "degK s^-1"},
    }
    diagnostic_properties = {}

    def __init__(self, sigma_boundary_layer_top=0.7, k_f=1/86400.0, k_a=1/40.0/86400.0, k_s=1/4.0/86400.0, equator_pole_temperature_difference=60, delta_theta_z=10, **kwargs):
        self._sigma_b = sigma_boundary_layer_top; self._k_f = k_f; self._k_a = k_a; self._k_s = k_s; self._delta_T_y = equator_pole_temperature_difference; self._delta_theta_z = delta_theta_z
        self._update_constants(); super(HeldSuarez, self).__init__(**kwargs)

    def _update_constants(self):
        self._p0 = get_constant("reference_air_pressure", "Pa")
        self._Cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK")
        self._R_d = get_constant("gas_constant_of_dry_air", "J/kg/degK")
        self._kappa = self._R_d / self._Cpd
        self._params = HeldSuarezParams(sigma_b=float(self._sigma_b), k_f=float(self._k_f), k_a=float(self._k_a), k_s=float(self._k_s), delta_T_y=float(self._delta_T_y), delta_theta_z=float(self._delta_theta_z), p0=float(self._p0), kappa=float(self._kappa))

    def array_call(self, raw_state):
        t = raw_state["air_temperature"]
        u = raw_state["eastward_wind"]; v = raw_state["northward_wind"]; p = raw_state["air_pressure"]; ps = raw_state["surface_air_pressure"]; lat = raw_state["latitude"]
        tend_u, tend_v, tend_t = _held_suarez_kernel_np(u, v, t, p, ps, lat, self._params)
        return {"eastward_wind": tend_u, "northward_wind": tend_v, "air_temperature": tend_t}, {}

@njit
def _held_suarez_kernel_np(u, v, t, p, ps, lat, params):
    # ps and lat are (ncol,)
    # p, u, v, t are (ncol, nlev)
    ncol, nlev = u.shape
    tend_u = np.zeros(u.shape); tend_v = np.zeros(v.shape); tend_t = np.zeros(t.shape)
    
    deg_to_rad = np.pi / 180.0
    
    for i in prange(ncol):
        ps_val = ps[i]; lat_val = lat[i]
        lat_rad = lat_val * deg_to_rad
        sin_lat_sq = np.sin(lat_rad)**2
        cos_lat_sq = np.cos(lat_rad)**2
        cos_lat_q = cos_lat_sq**2
        
        for j in range(nlev):
            sigma = p[i, j] / ps_val; p_norm = p[i, j] / params.p0
            Teq = max(200.0, (315.0 - params.delta_T_y * sin_lat_sq - params.delta_theta_z * np.log(p_norm) * cos_lat_sq) * p_norm**params.kappa)
            sigma_fac = max(0.0, (sigma - params.sigma_b) / (1.0 - params.sigma_b))
            k_t = params.k_a + (params.k_s - params.k_a) * sigma_fac * cos_lat_q
            k_v = params.k_f * sigma_fac
            tend_u[i, j] = -k_v * u[i, j]; tend_v[i, j] = -k_v * v[i, j]; tend_t[i, j] = -k_t * (t[i, j] - Teq)
    return tend_u, tend_v, tend_t

