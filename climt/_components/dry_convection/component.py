from typing import NamedTuple

import numpy as np
from sympl import Stepper, get_constant

from ..._core.backend import jit_compile, prange


class DryAdjParams(NamedTuple):
    Cpd: float
    Cvap: float
    Rdair: float
    Pref: float
    Rv: float


class DryConvectiveAdjustment(Stepper):
    """
    A conservative scheme to keep the temperature profile close to the
    dry adiabat if it is super-adiabatic.
    """

    input_properties = {
        "air_temperature": {
            "units": "degK",
            "dims": ["mid_levels", "*"],
        },
        "air_pressure": {
            "units": "Pa",
            "dims": ["mid_levels", "*"],
        },
        "air_pressure_on_interface_levels": {
            "units": "Pa",
            "dims": ["interface_levels", "*"],
            "alias": "P_int",
        },
        "specific_humidity": {
            "units": "kg/kg",
            "dims": ["mid_levels", "*"],
        },
    }

    output_properties = {
        "air_temperature": {
            "units": "degK",
        },
        "specific_humidity": {
            "units": "kg/kg",
        },
    }

    diagnostic_properties = {}

    def array_call(self, state, time_step):
        t = state["air_temperature"]
        q = state["specific_humidity"]
        p = state["air_pressure"]
        p_int = state["P_int"]

        orig_shape = t.shape
        t_flat = np.reshape(t, (t.shape[0], -1))
        q_flat = np.reshape(q, (q.shape[0], -1))
        p_flat = np.reshape(p, (p.shape[0], -1))
        p_int_flat = np.reshape(p_int, (p_int.shape[0], -1))

        params = DryAdjParams(
            Cpd=get_constant(
                "heat_capacity_of_dry_air_at_constant_pressure", "J/kg/degK"
            ),
            Cvap=get_constant("heat_capacity_of_vapor_phase", "J/kg/K"),
            Rdair=get_constant("gas_constant_of_dry_air", "J/kg/degK"),
            Pref=get_constant("reference_air_pressure", "Pa"),
            Rv=get_constant("gas_constant_of_vapor_phase", "J/kg/K"),
        )

        t_new, q_new = _dry_adj_kernel_np(t_flat, q_flat, p_flat, p_int_flat, params)

        return {}, {
            "air_temperature": np.reshape(t_new, orig_shape),
            "specific_humidity": np.reshape(q_new, orig_shape),
        }


@jit_compile
def _dry_adj_kernel_np(T, q, p, p_int, params):
    nlev, ncol = T.shape
    T_new = T.copy()
    q_new = q.copy()

    eps = params.Rv / params.Rdair

    for i in prange(ncol):
        p_col = p[:, i]
        p_int_col = p_int[:, i]
        pdiff = np.zeros(nlev)
        for m in range(nlev):
            pdiff[m] = p_int_col[m] - p_int_col[m + 1]

        # TOA to Surface
        for k in range(nlev - 1, -1, -1):
            rd_cp = np.zeros(nlev)
            theta_q = np.zeros(nlev)
            for m in range(nlev):
                rd_cp[m] = (
                    params.Rdair * (1.0 - q_new[m, i]) + params.Rv * q_new[m, i]
                ) / (params.Cpd * (1.0 - q_new[m, i]) + params.Cvap * q_new[m, i])
                theta_q[m] = (
                    T_new[m, i]
                    * (params.Pref / p_col[m]) ** rd_cp[m]
                    * (1.0 + q_new[m, i] * eps - q_new[m, i])
                )

            current_theta_sum = 0.0
            max_unstable_idx = -1

            for m in range(k, nlev):
                current_theta_sum += theta_q[m]
                theta_avg = current_theta_sum / (m - k + 1)

                if m > k and theta_avg > theta_q[m]:
                    max_unstable_idx = m

            if max_unstable_idx != -1:
                # Mix from k to max_unstable_idx
                p_high = p_int_col[k]
                p_low = p_int_col[max_unstable_idx + 1]

                sum_enthalpy = 0.0
                sum_q_dp = 0.0
                for m in range(k, max_unstable_idx + 1):
                    cp_m = params.Cpd * (1.0 - q_new[m, i]) + params.Cvap * q_new[m, i]
                    sum_enthalpy += cp_m * T_new[m, i] * pdiff[m]
                    sum_q_dp += q_new[m, i] * pdiff[m]

                mean_q = sum_q_dp / (p_high - p_low)
                cp_mixed = params.Cpd * (1.0 - mean_q) + params.Cvap * mean_q
                rd_mixed = params.Rdair * (1.0 - mean_q) + params.Rv * mean_q
                rdcp_mixed = rd_mixed / cp_mixed

                sum_theta_den = 0.0
                for m in range(k, max_unstable_idx + 1):
                    theta_coeff = (p_col[m] / params.Pref) ** rdcp_mixed
                    sum_theta_den += cp_mixed * theta_coeff * pdiff[m]

                mean_theta = sum_enthalpy / sum_theta_den

                for m in range(k, max_unstable_idx + 1):
                    q_new[m, i] = mean_q
                    theta_coeff = (p_col[m] / params.Pref) ** rdcp_mixed
                    T_new[m, i] = mean_theta * theta_coeff

    return T_new, q_new
