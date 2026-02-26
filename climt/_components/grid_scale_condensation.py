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
        "air_temperature": {
            "dims": ["mid_levels", "*"],
            "units": "degK",
        },
        "specific_humidity": {
            "dims": ["mid_levels", "*"],
            "units": "kg/kg",
        },
        "air_pressure": {"dims": ["mid_levels", "*"], "units": "Pa"},
        "air_pressure_on_interface_levels": {
            "dims": ["interface_levels", "*"],
            "units": "Pa",
        },
    }

    diagnostic_properties = {
        "precipitation_amount": {
            "dims": ["*"],
            "units": "kg m^-2",
        }
    }

    output_properties = {
        "air_temperature": {"units": "degK"},
        "specific_humidity": {"units": "kg/kg"},
    }

    def __init__(self, **kwargs):
        self._update_constants()
        super(GridScaleCondensation, self).__init__(**kwargs)

from sympl import Stepper, get_constant
from .._core import bolton_q_sat, bolton_dqsat_dT
import numpy as np
from typing import NamedTuple
from .._core.backend import get_array_namespace, jit_compile, prange

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    njit = lambda x, **kwargs: x

class GSCParams(NamedTuple):
    Cpd: float; Lv: float; Rd: float; Rh2O: float; g: float; rhow: float

class GridScaleCondensation(Stepper):
    """
    Calculate condensation due to supersaturation of water.

    Condenses supersaturated water at the grid scale, assuming all
    condensed water falls as precipitation.
    """

    input_properties = {
        "air_temperature": {
            "dims": ["mid_levels", "*"],
            "units": "degK",
        },
        "specific_humidity": {
            "dims": ["mid_levels", "*"],
            "units": "kg/kg",
        },
        "air_pressure": {"dims": ["mid_levels", "*"], "units": "Pa"},
        "air_pressure_on_interface_levels": {
            "dims": ["interface_levels", "*"],
            "units": "Pa",
        },
    }

    diagnostic_properties = {
        "precipitation_amount": {
            "dims": ["*"],
            "units": "kg m^-2",
        }
    }

    output_properties = {
        "air_temperature": {"units": "degK"},
        "specific_humidity": {"units": "kg/kg"},
    }

    def __init__(self, **kwargs):
        self._update_constants()
        super(GridScaleCondensation, self).__init__(**kwargs)

    def _update_constants(self):
        self._Cpd = get_constant(
            "heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"
        )
        self._Lv = get_constant("latent_heat_of_condensation", "J/kg")
        self._Rd = get_constant("gas_constant_of_dry_air", "J/kg/degK")
        self._Rh2O = get_constant("gas_constant_of_vapor_phase", "J/kg/degK")
        self._g = get_constant("gravitational_acceleration", "m/s^2")
        self._rhow = get_constant("density_of_liquid_phase", "kg/m^3")
        self._params = GSCParams(
            Cpd=float(self._Cpd), Lv=float(self._Lv), Rd=float(self._Rd),
            Rh2O=float(self._Rh2O), g=float(self._g), rhow=float(self._rhow)
        )

    def array_call(self, state, timestep):
        t = state["air_temperature"]; q = state["specific_humidity"]
        p = state["air_pressure"]; p_int = state["air_pressure_on_interface_levels"]
        
        xp = get_array_namespace(t)
        orig_shape = t.shape
        
        t_flat = xp.reshape(t, (t.shape[0], -1))
        q_flat = xp.reshape(q, (q.shape[0], -1))
        p_flat = xp.reshape(p, (p.shape[0], -1))
        p_int_flat = xp.reshape(p_int, (p_int.shape[0], -1))
        
        if xp is np:
            t_new, q_new, precip = _gsc_kernel_np(t_flat, q_flat, p_flat, p_int_flat, self._params)
        else:
            t_new, q_new, precip = _gsc_kernel_jax(t_flat, q_flat, p_flat, p_int_flat, self._params)
            
        return {"precipitation_amount": xp.reshape(precip, orig_shape[1:])}, {
            "air_temperature": xp.reshape(t_new, orig_shape),
            "specific_humidity": xp.reshape(q_new, orig_shape)
        }

@njit
def _gsc_kernel_np(T, q, p, p_int, params):
    nlev, ncol = T.shape
    new_T = T.copy()
    new_q = q.copy()
    precip = np.zeros(ncol)
    
    # Pre-calculate q_sat
    # es = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    # epsilon = params.Rd / params.Rh2O
    # q_sat = epsilon * es / (p - (1 - epsilon) * es)
    
    for i in prange(ncol):
        p_col = p[:, i]
        p_int_col = p_int[:, i]
        T_col = T[:, i]
        q_col = q[:, i]
        
        col_precip = 0.0
        for k in range(nlev):
            # Bolton q_sat logic
            es = 611.2 * np.exp(17.67 * (T_col[k] - 273.15) / (T_col[k] - 29.65))
            eps_val = params.Rd / params.Rh2O
            q_sat_k = eps_val * es / (p_col[k] - (1.0 - eps_val) * es)
            
            if q_col[k] > q_sat_k:
                dqsat_dT = params.Lv * q_sat_k / (params.Rh2O * T_col[k]**2)
                condensed_q = (q_col[k] - q_sat_k) / (1.0 + params.Lv / params.Cpd * dqsat_dT)
                
                new_q[k, i] = q_col[k] - condensed_q
                new_T[k, i] = T_col[k] + params.Lv / params.Cpd * condensed_q
                
                dp = p_int_col[k] - p_int_col[k+1]
                mass = dp / (params.g * params.rhow)
                col_precip += condensed_q * mass
        
        precip[i] = col_precip
        
    return new_T, new_q, precip

def _gsc_kernel_jax(T, q, p, p_int, params):
    import jax.numpy as jnp
    
    q_sat = bolton_q_sat(T, p, params.Rd, params.Rh2O)
    dqsat_dT = bolton_dqsat_dT(T, params.Lv, params.Rh2O, q_sat)
    
    condensed_q = jnp.where(q > q_sat, (q - q_sat) / (1.0 + params.Lv / params.Cpd * dqsat_dT), 0.0)
    
    new_q = q - condensed_q
    new_T = T + params.Lv / params.Cpd * condensed_q
    
    dp = p_int[:-1, :] - p_int[1:, :]
    mass = dp / (params.g * params.rhow)
    precip = jnp.sum(condensed_q * mass, axis=0)
    
    return new_T, new_q, precip
