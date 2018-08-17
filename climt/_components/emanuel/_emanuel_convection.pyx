cimport numpy as cnp
import numpy as np
cimport cython
from cython.parallel import prange

'''
Code signatures for the emanuel convection fortran code
'''

cdef extern:
    void init_emanuel_convection_fortran(
        cnp.int32_t *pbl, cnp.int32_t *least_conv_level,
        double *thresh_water_level,
        double *crit_temp, double *entrain_coeff,
        double *downdraft_frac_area,
        double *precip_frac_outside_cloud,
        double *rain_speed, double *snow_speed,
        double *rain_evap_coeff, double *snow_evap_coeff,
        double *mom_tran_coeff, double *max_neg_temp_pert,
        double *beta, double *alpha, double *damp_amp,
        double *Cpd, double *Cpv, double *Cl,
        double *gas_const_vapour, double *gas_const_air,
        double *lat_heat, double *grav, double *density_water,
        double *reference_mass_flux_timescale)

cdef extern nogil:
    void emanuel_convection(double *temp, double *q, double *qs,
                            double *u, double *v,
                            double *pmid, double *pint,
                            cnp.int32_t *nlevs,
                            cnp.int32_t *max_conv_lev,
                            cnp.int32_t *num_tracers,
                            double *dt, cnp.int32_t *conv_state,
                            double *dtemp, double *dq, double *du,
                            double *dv, double *precip,
                            double *downdraft_vel_scale,
                            double *downdraft_temp_scale,
                            double *downdraft_q_scale,
                            double *cloud_base_mass_flux,
                            double *cape, double *tracers,
                            double *dtracers)


cdef int initialised
'''
Arrays and variables to pass to the fortran code
'''

def init_emanuel_convection(
    cnp.int32_t py_least_conv_level,
    double py_thresh_water_level,
    double py_crit_temp, double py_entrain_coeff,
    double py_downdraft_frac_area,
    double py_precip_frac_outside_cloud,
    double py_rain_speed, double py_snow_speed,
    double py_rain_evap_coeff, double py_snow_evap_coeff,
    double py_mom_tran_coeff, double py_max_neg_temp_pert,
    double py_beta, double py_alpha, double py_damp_amp,
    double py_Cpd, double py_Cpv, double py_Cl,
    double py_gas_const_vapour, double py_gas_const_air,
    double py_lat_heat, double py_grav,
    double py_density_water,
    double py_reference_timescale):

    global initialised

    # The dry convective PBL adjustment modifies inputs, so it has been moved
    # to a separate component. Here we disable it.
    cdef cnp.int32_t py_pbl = 0

    initialised = 0

    init_emanuel_convection_fortran(
        &py_pbl,
        &py_least_conv_level,
        &py_thresh_water_level,
        &py_crit_temp, &py_entrain_coeff,
        &py_downdraft_frac_area,
        &py_precip_frac_outside_cloud,
        &py_rain_speed, &py_snow_speed,
        &py_rain_evap_coeff, &py_snow_evap_coeff,
        &py_mom_tran_coeff, &py_max_neg_temp_pert,
        &py_beta, &py_alpha, &py_damp_amp,
        &py_Cpd, &py_Cpv, &py_Cl,
        &py_gas_const_vapour, &py_gas_const_air,
        &py_lat_heat, &py_grav,
        &py_density_water,
        &py_reference_timescale)


    # max_conv_lev = nlevs - 3 #Max value = nlevs - 1

    initialised = 1

@cython.boundscheck(False)
def convect(
    cnp.int32_t nlevs,
    cnp.int32_t ncols,
    cnp.int32_t max_conv_lev,
    cnp.int32_t num_tracers,
    double dt,
    cnp.double_t[:, ::1] temp,
    cnp.double_t[:, ::1] q,
    cnp.double_t[:, ::1] qs,
    cnp.double_t[:, ::1] u,
    cnp.double_t[:, ::1] v,
    cnp.double_t[:, ::1] pmid,
    cnp.double_t[:, ::1] pint,
    cnp.int32_t[::1] conv_state,
    cnp.double_t[::1] precip,
    cnp.double_t[::1] downdraft_vel_scale,
    cnp.double_t[::1] downdraft_temp_scale,
    cnp.double_t[::1] downdraft_q_scale,
    cnp.double_t[::1] cloud_base_mass_flux,
    cnp.double_t[::1] cape,
    cnp.double_t[:, ::1] dtemp,
    cnp.double_t[:, ::1] dq,
    cnp.double_t[:, ::1] du,
    cnp.double_t[:, ::1] dv,
    cnp.double_t[:, :, ::1] dtracers=None,
    cnp.double_t[:, :, ::1] tracers=None):


    global initialised
    cdef Py_ssize_t col, lon
    cdef cnp.double_t[::1] cy_temp, cy_q, cy_qs, cy_u,\
            cy_v, cy_pmid, cy_pint, cy_dtemp, cy_dq,\
            cy_du, cy_dv
    cdef cnp.double_t[:, ::1] cy_tracer, cy_dtracer

    cy_temp = np.zeros(nlevs)
    cy_q = np.zeros(nlevs)
    cy_qs = np.zeros(nlevs)
    cy_u = np.zeros(nlevs)
    cy_v = np.zeros(nlevs)
    cy_pmid = np.zeros(nlevs)
    cy_pint = np.zeros(nlevs+1)
    cy_dtemp = np.zeros(nlevs)
    cy_dq = np.zeros(nlevs)
    cy_du = np.zeros(nlevs)
    cy_dv = np.zeros(nlevs)
    cy_tracer = np.zeros((num_tracers, nlevs))
    cy_dtracer = np.zeros((num_tracers, nlevs))

    if dtracers is None:
        dtracers = np.zeros((ncols, 1, nlevs))

    if tracers is None:
        tracers = np.zeros((ncols, 1, nlevs))

    if initialised == 0:
        raise ValueError('Emanuel scheme not initialised.')

    for col in range(ncols):
    # for lat in prange(nlats, nogil=True):

        cy_temp[:] = temp[col, :]
        cy_q[:] = q[col, :]
        cy_qs[:] = qs[col, :]
        cy_u[:] = u[col, :]
        cy_v[:] = v[col, :]
        cy_pmid[:] = pmid[col, :]
        cy_pint[:] = pint[col, :]
        cy_tracer[:] = tracers[col, :, :]

        # Call fortran code
        emanuel_convection(
            <double *>&cy_temp[0],
            <double *>&cy_q[0],
            <double *>&cy_qs[0],
            <double *>&cy_u[0],
            <double *>&cy_v[0],
            <double *>&cy_pmid[0],
            <double *>&cy_pint[0],

            &nlevs, &max_conv_lev, &num_tracers,
            &dt,

            <cnp.int32_t *>&conv_state[col],
            <double *>&cy_dtemp[0],
            <double *>&cy_dq[0],
            <double *>&cy_du[0],
            <double *>&cy_dv[0],

            <double *>&precip[col],
            <double *>&downdraft_vel_scale[col],
            <double *>&downdraft_temp_scale[col],
            <double *>&downdraft_q_scale[col],
            <double *>&cloud_base_mass_flux[col],
            <double *>&cape[col],
            <double *>&cy_tracer[0, 0],
            <double *>&cy_dtracer[0, 0])

        dtemp[col, :] = cy_dtemp[:]
        dq[col, :] = cy_dq[:]
        du[col, :] = cy_du[:]
        dv[col, :] = cy_dv[:]
        #dtracers[col, :, :] = cy_dtracer[:, :]
