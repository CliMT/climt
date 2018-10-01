cimport cython
cimport numpy as cnp
import numpy as np

@cython.boundscheck(False)
cpdef void calculate_dry_adjustment(
    cnp.ndarray [double, ndim=2] P_int,
    cnp.ndarray [double, ndim=2] P,
    cnp.ndarray [double, ndim=2] theta_q,
    cnp.ndarray [double, ndim=2] output_q,
    cnp.ndarray [double, ndim=2] output_temp,
    double cpd,
    double cvap,
    double rdair,
    double rv,
    double pref,
    int num_columns, int num_levels) except *:

    cdef int column, level

    for column in range(num_columns):
        for level in range(num_levels-1, -1, -1):

            dp = P_int[column, :-1] - P_int[column, 1:]
            theta_sum = np.cumsum(theta_q[column, level::])
            divisor = np.arange(1, num_levels - level+1)

            theta_avg = (theta_sum/divisor)[1::]

            theta_lesser = (theta_avg > theta_q[column, level+1::])
            if np.sum(theta_lesser) == 0:
                continue

            convect_to_level = len(theta_lesser) - np.argmax(theta_lesser[::-1])

            if level == 0:
                convect_to_level = max(convect_to_level, 1)

            if convect_to_level == 0:
                continue
            stable_level = level + convect_to_level

            q_conv = output_q[column, level:stable_level]
            t_conv = output_temp[column, level:stable_level]
            dp_conv = dp[level:stable_level]
            p_conv_high = P_int[column, level]
            p_conv_low = P_int[column, stable_level]

            enthalpy = heat_capacity(q_conv, cpd, cvap)*t_conv
            integral_enthalpy = np.sum(enthalpy*dp_conv)
            mean_conv_q = np.sum(q_conv*dp_conv)/(p_conv_high - p_conv_low)

            output_q[column, level:stable_level] = mean_conv_q

            rdcp_conv = gas_constant(mean_conv_q, rdair, rv)/heat_capacity(
                mean_conv_q, cpd, cvap)

            theta_coeff = (
                P[column, level:stable_level]/pref)**rdcp_conv

            integral_theta_den = np.sum(heat_capacity(q_conv, cpd,
                                                       cvap)*theta_coeff*dp_conv)

            mean_theta = integral_enthalpy/integral_theta_den

            output_temp[column, level:stable_level] = mean_theta*theta_coeff


ctypedef fused float_or_arr:
    double
    cnp.ndarray

cpdef double sum_arr(double [:] arr) nogil:

    cdef int len_arr = arr.shape[0]
    cdef int index
    cdef double result

    for index in range(len_arr):
        result += arr[index]

    return result


@cython.boundscheck(False)
def heat_capacity(
    float_or_arr q,
    double cpd,
    double cvap):
    """
    Calculate heat capacity based on amount of q
    """

    return cpd*(1-q) + cvap*q

@cython.boundscheck(False)
def gas_constant(
    float_or_arr q,
    double rdair,
    double rv):
    """
    Calculate gas constant based on amount of q
    """

    return rdair*(1-q) + rv*q
