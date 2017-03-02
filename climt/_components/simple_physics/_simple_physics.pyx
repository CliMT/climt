cimport numpy as cnp
import numpy as np

# Simple physics package fortran code signature
# For a description, see Reed and Jablonowski (2012)
cdef extern:
   void simple_physics(int *numColumns, int *numLevels,
                           double *dt, double *latitude,
                           double *t, double *q,
                           double *u, double* v,
                           double *pressure, double *ifacePressure,
                           double *layerThickness, double *invThickness,
                           double *surfPressure, double *precipRate,
                           int *test, int *do_lsc, int *do_pbl,
                           int *do_surf_flux, int *use_ts_ext, double *t_surf)

cdef extern:
    void set_fortran_constants(double *grav, double *cpd, double *r_air,
                                double *latent_heat, double *r_condensible,
                                double *planetary_radius, double *rotation_rate)

#Arrays to pass to fortran. The dimensions are (levels, num_columns) for 3d and (num_columns) for 2d
cdef cnp.double_t[::1] temp, q, u, v, p, p_iface, thickness, inv_thickness
cdef double surf_press, precip, latitude, surf_temp
cdef int num_cols, num_levs, initialised
cdef double time_step
cdef int cyclone, lsc, pbl, lhf, ext_ts

def set_physical_constants(double grav, double cpd, double r_air, double latent_heat, 
                           double r_cond, double radius, double rotation):
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

    set_fortran_constants(&grav, &cpd, &r_air, &latent_heat, &r_cond, &radius, &rotation)

    initialised = 1

def init_simple_physics(int sim_cyclone, int sim_lsc, int sim_pbl, int sim_lhf, int use_ext_ts):

    global cyclone, lsc, pbl, lhf, ext_ts

    cyclone = sim_cyclone
    lsc = sim_lsc
    pbl = sim_pbl
    lhf = sim_lhf
    ext_ts = use_ext_ts


def init_arrays(int columns, int lons, int lats, int levels, double dt, int simulate_cyclone):
    '''
    Initialise global arrays and some constants
    '''

    global temp, q, u, v, p, p_iface, thickness, inv_thickness,\
        surf_press, precip, latitude, num_cols, num_levs, initialised, time_step

    temp = np.empty((levels), dtype=np.double, order='F')
    q = np.empty((levels), dtype=np.double, order='F')
    u = np.empty((levels), dtype=np.double, order='F')
    v = np.empty((levels), dtype=np.double, order='F')
    p = np.empty((levels), dtype=np.double, order='F')
    thickness = np.empty((levels), dtype=np.double, order='F')
    inv_thickness = np.empty((levels), dtype=np.double, order='F')

    p_iface = np.empty((levels+1), dtype=np.double, order='F')

    test = np.zeros((lons,lats,levels), dtype=np.double, order='F')

    #surf_press = np.empty(columns, dtype=np.double, order='F')
    #precip = np.empty(columns, dtype=np.double, order='F')
    #latitude = np.empty(columns, dtype=np.double, order='F')

    #TODO Until we are happy and stable, only one column
    num_cols = 1
    num_levs = levels
    time_step = dt
    cyclone_simulation = simulate_cyclone

    initialised = 1

#Returns a list of tendencies
def get_new_state(cnp.ndarray[cnp.double_t, ndim=3] u_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] v_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] temp_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] p_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] p_iface_ext,
                    cnp.ndarray[cnp.double_t, ndim=3] q_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] surf_press_ext,
                    cnp.ndarray[cnp.double_t, ndim=2] ts_ext,
                    cnp.ndarray[cnp.double_t, ndim=1] latitude_ext,
                    double time_step):
    '''
    This function takes fetches tendencies from the Simple Physics code.    

    '''

    global temp, q, u, v, p, p_iface, latitude, thickness, inv_thickness,\
        surf_press, precip, num_cols, num_levs, initialised, \
        lsc, pbl, lhf, ext_ts, cyclone, surf_temp

    cdef  cnp.double_t[::1] ptr
    if not initialised: 
        print 'Module not initialised'
        return

    nlons, nlats, nlevs = np.asfortranarray(u_ext).shape

    num_cols = 1
    num_levs = nlevs

    temp = np.empty((nlevs), dtype=np.double, order='F')
    q = np.empty((nlevs), dtype=np.double, order='F')
    u = np.empty((nlevs), dtype=np.double, order='F')
    v = np.empty((nlevs), dtype=np.double, order='F')
    p = np.empty((nlevs), dtype=np.double, order='F')
    thickness = np.empty((nlevs), dtype=np.double, order='F')
    inv_thickness = np.empty((nlevs), dtype=np.double, order='F')

    p_iface = np.empty((nlevs+1), dtype=np.double, order='F')

    u_out = np.zeros((nlons,nlats,nlevs),order='F')
    v_out = np.zeros((nlons,nlats,nlevs),order='F')
    t_out = np.zeros((nlons,nlats,nlevs),order='F')
    q_out = np.zeros((nlons,nlats,nlevs),order='F')
    precip_out = np.zeros((nlons,nlats),order='F')

    for lon in range(nlons):
        for lat in range(nlats):

            ptr = np.asfortranarray(u_ext[lon,lat,::-1])
            u[:] = ptr[:]

            ptr = np.asfortranarray((p_iface_ext[lon,lat,0:-1] - p_iface_ext[lon,lat,1::])[::-1])
            thickness[:] = ptr[:]
            
            ptr = np.asfortranarray(temp_ext[lon,lat,::-1])
            temp[:] = ptr[:]

            ptr = np.asfortranarray(v_ext[lon,lat,::-1])
            v[:] = ptr[:]

            ptr = np.asfortranarray(q_ext[lon,lat,::-1])
            q[:] = ptr[:]

            ptr = np.asfortranarray(p_ext[lon,lat,::-1])
            p[:] = ptr[:]

            ptr = np.asfortranarray(p_iface_ext[lon,lat,::-1])
            p_iface[:] = ptr[:]

            ptr = 1./np.asfortranarray(thickness)
            inv_thickness[:] = ptr[:] 

            surf_press = surf_press_ext[lon,lat]

            latitude = latitude_ext[lat]

            surf_temp = ts_ext[lon,lat]


            #Call fortran code with these arguments
            simple_physics(&num_cols, &num_levs,
                       &time_step, &latitude,
                       <double *>&temp[0], <double *>&q[0],
                       <double *>&u[0], <double *>&v[0],
                       <double *>&p[0], <double *>&p_iface[0],
                       <double *>&thickness[0], <double *>&inv_thickness[0],
                       &surf_press, &precip, &cyclone,
                        &lsc, &pbl, &lhf,
                        &ext_ts, &surf_temp)


            t_out[lon,lat,:] = temp[::-1]
            u_out[lon,lat,:] = u[::-1]
            v_out[lon,lat,:] = v[::-1]
            q_out[lon,lat,:] = q[::-1]
            precip_out[lon,lat] = precip
    

    return t_out, u_out, v_out, q_out, precip_out
