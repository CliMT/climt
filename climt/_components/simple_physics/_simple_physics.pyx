cimport numpy as cnp
import numpy as np

# Simple physics package fortran code signature
# For a description, see Reed and Jablonowski (2012)
cdef extern:
   void simple_physics(
       cnp.int32_t *numColumns, cnp.int32_t *numLevels,
       double *dt, double *latitude,
       double *t, double *q,
       double *u, double* v,
       double *pressure, double *ifacePressure,
       double *layerThickness, double *invThickness,
       double *surfPressure, double *precipRate,
       cnp.int32_t *test, cnp.int32_t *do_lsc, cnp.int32_t *do_pbl,
       cnp.int32_t *do_surf_flux, cnp.int32_t *use_ts_ext, double *t_surf,
       cnp.int32_t *use_qsurf_ext, double *q_surf,
       double *sens_ht_flx, double *lat_ht_flx)

cdef extern:
    void set_fortran_constants(
        double *grav, double *cpd, double *r_air,
        double *latent_heat, double *r_condensible,
        double *planetary_radius, double *rotation_rate,
        double *density_condensible, double *top_pbl,
        double *pbl_decay, double *drag_coeff_sens_lat,
        double *Cd0_ext, double *Cd1_ext, double *Cm_ext)

#Arrays to pass to fortran. The dimensions are (levels, num_columns) for 3d and (num_columns) for 2d
cdef cnp.double_t[:] temp, q, u, v, p, p_iface, thickness, inv_thickness
cdef cnp.double_t[:] out_temp, out_q, out_u, out_v
cdef double surf_press, precip, latitude, surf_temp
cdef cnp.int32_t num_cols, num_levs, initialised
cdef double time_step
cdef cnp.int32_t cyclone, lsc, pbl, lhf, ext_ts, use_ext_qsurf


def set_physical_constants(double grav, double cpd, double r_air, double latent_heat,
                           double r_cond, double radius, double rotation,
                           double rho_cond, double pbltop, double pbldecay,
                           double drag_latent_sens, double drag_base_mom,
                           double drag_varying_mom, double drag_max_mom):
    '''
    Pass the model physical constants to the simple physics fortran code

    Args:
    grav (float): value of acceleration due to gravity.

    cpd (float): The gas constant of dry air.

    r_air (float): The specific heat of dry air.

    latent_heat (float): The latent heat of condensation of the condensible.

    r_cond (float): The gas constant of the condensible substance.

    radius (float): The radius of the planet.

    rotation (float): The rotation rate of the planet.
    '''
    global initialised

    set_fortran_constants(&grav, &cpd, &r_air, &latent_heat, &r_cond, &radius, &rotation,
                          &rho_cond, &pbltop, &pbldecay, &drag_latent_sens,
                          &drag_base_mom, &drag_varying_mom, &drag_max_mom)

    initialised = 1


def init_simple_physics(int sim_cyclone, int sim_lsc, int sim_pbl,
                        int sim_lhf, int use_ext_ts, int use_ext_qsrf):

    global cyclone, lsc, pbl, lhf, ext_ts, use_ext_qsurf

    cyclone = sim_cyclone
    lsc = sim_lsc
    pbl = sim_pbl
    lhf = sim_lhf
    ext_ts = use_ext_ts
    use_ext_qsurf = use_ext_qsrf


#Returns a list of tendencies
def get_new_state(cnp.ndarray[cnp.double_t, ndim=2] u_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] v_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] T_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] p_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] p_iface_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] q_ext,
                    cnp.ndarray[cnp.double_t, ndim=1] surf_press_ext,
                    cnp.ndarray[cnp.double_t, ndim=1] ts_ext,
                    cnp.ndarray[cnp.double_t, ndim=1] qsurf_ext,
                    cnp.ndarray[cnp.double_t, ndim=1] latitude_ext,
                    double time_step):
    '''
    This function takes fetches new state from the Simple Physics code.

    '''

    if not initialised:
        print 'Module not initialised'
        return

    num_cols, num_levs = u_ext.shape

    u_ext = np.ascontiguousarray(u_ext[:, ::-1])
    v_ext = np.ascontguousarray(v_ext[:, ::-1])
    q_ext = np.ascontiguousarray(q_ext[:, ::-1])
    T_ext = np.ascontiguousarray(T_ext[:, ::-1])
    p_ext = np.ascontiguousarray(p_ext[:, ::-1])
    p_iface_ext = np.ascontiguousarray(p_iface_ext[:, ::-1])

    thickness = np.ascontiguousarray(p_iface_ext[1::, :] - p_iface_ext[0:-1, :])
    inv_thickness = np.ascontiguousarray(1./thickness)

    precip_out = np.zeros([num_cols])
    sh_flx_out = np.zeros([num_cols])
    lh_flx_out = np.zeros([num_cols])

    do_simple_physics(
        num_cols, num_levs, time_step,
        u_ext, v_ext, T_ext, p_ext, p_iface_ext,
        thickness, inv_thickness,
        q_ext, surf_press_ext, ts_ext, qsurf_ext, latitude_ext,
        precip_out, sh_flx_out, lh_flx_out)

    u_ext = u_ext[:, ::-1]
    v_ext = v_ext[:, ::-1]
    q_ext = q_ext[:, ::-1]
    T_ext = T_ext[:, ::-1]
    return T_ext, u_ext, v_ext, q_ext, precip_out, sh_flx_out, lh_flx_out


def do_simple_physics(
    cnp.int32_t num_cols,
    cnp.int32_t num_levs,
    cnp.double_t time_step,
    cnp.double_t[:, ::1] u_ext,
    cnp.double_t[:, ::1] v_ext,
    cnp.double_t[:, ::1] temp_ext,
    cnp.double_t[:, ::1] p_ext,
    cnp.double_t[:, ::1] p_iface_ext,
    cnp.double_t[:, ::1] thickness_ext,
    cnp.double_t[:, ::1] inv_thickness_ext,
    cnp.double_t[:, ::1] q_ext,
    cnp.double_t[::1] surf_press_ext,
    cnp.double_t[::1] ts_ext,
    cnp.double_t[::1] qsurf_ext,
    cnp.double_t[::1] latitude_ext,
    cnp.double_t[::1] precip_out,
    cnp.double_t[::1] sh_flx_out,
    cnp.double_t[::1] lh_flx_out):

    global lsc, pbl, lhf, ext_ts, cyclone, use_ext_qsurf
    cdef cnp.int32_t one = 1

    for i_col in range(num_cols):
        #Call fortran code with these arguments
        simple_physics(
            &one, &num_levs,
            &time_step,
            &latitude_ext[i_col],
            <double *>&temp_ext[i_col, 0],
            <double *>&q_ext[i_col, 0],
            <double *>&u_ext[i_col, 0],
            <double *>&v_ext[i_col, 0],
            <double *>&p_ext[i_col, 0],
            <double *>&p_iface_ext[i_col, 0],
            <double *>&thickness_ext[i_col, 0],
            <double *>&inv_thickness_ext[i_col, 0],
            &surf_press_ext[i_col],
            &precip_out[i_col],
            &cyclone,
            &lsc, &pbl, &lhf,
            &ext_ts,
            &ts_ext[i_col],
            &use_ext_qsurf,
            &qsurf_ext[i_col],
            &sh_flx_out[i_col],
            &lh_flx_out[i_col])
