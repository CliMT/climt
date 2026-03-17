# -*- coding: utf-8 -*-
from sympl import TendencyComponent, DiagnosticComponent, get_constant
import numpy as np
from .._core.backend import jit_compile, prange

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x: x

class GrayLongwaveRadiation(TendencyComponent):

    input_properties = {
        "longwave_optical_depth_on_interface_levels": {
            "dims": ["interface_levels", "*"],
            "units": "dimensionless",
            "alias": "tau",
        },
        "air_temperature": {
            "dims": ["mid_levels", "*"],
            "units": "degK",
            "alias": "sl",
        },
        "surface_temperature": {
            "dims": ["*"],
            "units": "degK",
            "alias": "T_surface",
        },
        "air_pressure": {
            "dims": ["mid_levels", "*"],
            "units": "Pa",
            "alias": "p",
        },
        "air_pressure_on_interface_levels": {
            "dims": ["interface_levels", "*"],
            "units": "Pa",
            "alias": "p_interface",
        },
    }

    diagnostic_properties = {
        "downwelling_longwave_flux_in_air": {
            "dims": ["interface_levels", "*"],
            "units": "W m^-2",
            "alias": "lw_down",
        },
        "upwelling_longwave_flux_in_air": {
            "dims": ["interface_levels", "*"],
            "units": "W m^-2",
            "alias": "lw_up",
        },
        "longwave_heating_rate": {"dims": ["mid_levels", "*"], "units": "degK day^-1"},
    }

    tendency_properties = {
        "air_temperature": {"units": "degK s^-1"},
    }

    def array_call(self, state):
        t = getattr(state["sl"], "data", state["sl"])
        tau = getattr(state["tau"], "data", state["tau"])
        t_surf = getattr(state["T_surface"], "data", state["T_surface"])
        p_int = getattr(state["p_interface"], "data", state["p_interface"])

        orig_shape_t = t.shape
        orig_shape_p_int = p_int.shape

        t_flat = np.reshape(t, (t.shape[0], -1))
        tau_flat = np.reshape(tau, (tau.shape[0], -1))
        t_surf_flat = np.reshape(t_surf, (-1,))
        p_int_flat = np.reshape(p_int, (p_int.shape[0], -1))

        sigma_sb = get_constant("stefan_boltzmann_constant", "W/m^2/K^4")
        g = get_constant("gravitational_acceleration", "m/s^2")
        cpd = get_constant("heat_capacity_of_dry_air_at_constant_pressure", "J/kg/K")

        t_flat = np.asarray(t_flat); tau_flat = np.asarray(tau_flat); t_surf_flat = np.asarray(t_surf_flat); p_int_flat = np.asarray(p_int_flat)
        downward_flux, upward_flux = _gray_lw_kernel_np(
            t_flat, p_int_flat, t_surf_flat, tau_flat, sigma_sb
        )

        # Compute tendency in numpy so that float32 fluxes are upcast to float64
        # (matching the original pure-Python implementation's output dtype).
        net_lw_flux = upward_flux - downward_flux
        lw_temperature_tendency = (
            g / cpd
            * (net_lw_flux[1:, :] - net_lw_flux[:-1, :])
            / (p_int_flat[1:, :] - p_int_flat[:-1, :])
        )

        downward_flux = np.reshape(downward_flux, orig_shape_p_int)
        upward_flux = np.reshape(upward_flux, orig_shape_p_int)
        lw_temperature_tendency = np.reshape(lw_temperature_tendency, orig_shape_t)

        return {"sl": lw_temperature_tendency}, {"lw_down": downward_flux, "lw_up": upward_flux, "longwave_heating_rate": lw_temperature_tendency * 86400.0}


class Frierson06LongwaveOpticalDepth(DiagnosticComponent):

    input_properties = {
        "air_pressure_on_interface_levels": {"dims": ["interface_levels", "*"], "units": "Pa"},
        "surface_air_pressure": {"dims": ["*"], "units": "Pa"},
        "latitude": {"dims": ["*"], "units": "degrees_N"},
    }
    diagnostic_properties = {"longwave_optical_depth_on_interface_levels": {"dims": ["interface_levels", "*"], "units": "dimensionless"}}

    def __init__(self, linear_optical_depth_parameter=0.1, longwave_optical_depth_at_equator=6, longwave_optical_depth_at_poles=1.5, **kwargs):
        self._fl = linear_optical_depth_parameter; self._tau0e = longwave_optical_depth_at_equator; self._tau0p = longwave_optical_depth_at_poles; super(Frierson06LongwaveOpticalDepth, self).__init__(**kwargs)

    def array_call(self, state):
        p_int = state["air_pressure_on_interface_levels"]
        ps = getattr(state["surface_air_pressure"], "data", state["surface_air_pressure"])
        lat = getattr(state["latitude"], "data", state["latitude"])
        p_int_raw = getattr(p_int, "data", p_int)
        orig_shape = p_int_raw.shape
        p_int_flat = np.reshape(p_int_raw, (p_int_raw.shape[0], -1))
        ps_flat = np.reshape(ps, (-1,)); lat_flat = np.reshape(lat, (-1,))
        p_int_flat = np.asarray(p_int_flat); ps_flat = np.asarray(ps_flat); lat_flat = np.asarray(lat_flat)
        tau = _frierson_tau_kernel_np(lat_flat, p_int_flat, ps_flat, self._tau0e, self._tau0p, self._fl)
        return {"longwave_optical_depth_on_interface_levels": np.reshape(tau, orig_shape)}

@njit
def _gray_lw_kernel_np(T, p_interface, T_surface, tau, sigma):
    nlev, ncol = T.shape
    upward_flux = np.zeros((nlev + 1, ncol))
    downward_flux = np.zeros((nlev + 1, ncol))

    # Pre-calculate sigma * T**4
    T4 = sigma * T**4

    for i in prange(ncol):
        # Upward flux (surface to TOA)
        upward_flux[0, i] = sigma * T_surface[i]**4
        for k in range(1, nlev + 1):
            dtau = tau[k, i] - tau[k - 1, i]
            trans = np.exp(-dtau)
            upward_flux[k, i] = upward_flux[k - 1, i] * trans + T4[k - 1, i] * (1.0 - trans)

        # Downward flux (TOA to surface)
        downward_flux[nlev, i] = 0.0
        for k in range(nlev - 1, -1, -1):
            dtau = tau[k + 1, i] - tau[k, i]
            trans = np.exp(-dtau)
            downward_flux[k, i] = downward_flux[k + 1, i] * trans + T4[k, i] * (1.0 - trans)

    return downward_flux, upward_flux


@njit
def _frierson_tau_kernel_np(latitude, p_interface, ps, tau0e, tau0p, fl):
    nlev_plus_1, ncol = p_interface.shape
    tau = np.zeros((nlev_plus_1, ncol))
    deg_to_rad = np.pi / 180.0
    
    for i in prange(ncol):
        lat_rad = latitude[i] * deg_to_rad
        tau_0 = tau0e + (tau0p - tau0e) * np.sin(lat_rad) ** 2
        ps_val = ps[i]
        for k in range(nlev_plus_1):
            sigma = p_interface[k, i] / ps_val
            tau[k, i] = tau_0 * (1.0 - (fl * sigma + (1.0 - fl) * sigma**4))
    return tau

