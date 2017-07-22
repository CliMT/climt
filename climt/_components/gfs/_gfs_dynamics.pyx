cimport numpy as cnp
import numpy as np
from numpy import empty


# Typedef for function pointer returning void and taking no arguments (for now)

ctypedef void (*pyPhysicsCallback)(double *, double *, double *, double *,
                                   double *, double *, double *,
                              double *, double *, double *, double *,
                                   double *, double *, double *)


# Variables that control the behaviour of the model
# bint is Cython's boolean type
cdef extern:
    bint dry, adiabatic, heldsuarez

#Function to set global mass of dry air
cdef extern:
    void gfs_set_dry_pressure(double *pdry_init)

#Function to set time step
cdef extern:
    void gfs_set_time_step(double *dt)

#Function to get longitude/latitude map from shtns
cdef extern:
    void gfs_get_lon_lat(double *longitudes, double *latitudes)

#Function to set physical constants in the dycore
cdef extern:
    void gfs_set_constants(double *planetary_radius,
                           double *omega,
                           double *gas_constant,
                           double *gas_constant_dry_air,
                           double *gas_constant_condensible,
                           double *gravity,
                           double *specific_heat_dry_air,
                           double *specific_heat_condensible)

#Function to set grid dimensions in dycore
cdef extern:
    void gfs_set_model_dimensions(int *num_lats,
                                  int *num_lons,
                                  int *num_levs,
                                  int *truncation,
                                  int *spectral_dim,
                                  int *num_tracers)

#Function to set model time (might be eliminated if I have time!
#since it is not really used)
cdef extern:
    void gfs_set_model_time(double *time)

# Function to init dynamics
cdef extern:
    void gfs_init_dynamics(int *num_damp_levs, double *tau_damping)

# Function to init physics
cdef extern:
    void gfs_init_physics()

#Function to calculate pressure fields after other fields have been updated
cdef extern:
    void gfs_calculate_pressure()

# Function to step fields by one dt
cdef extern:
    void gfs_take_one_step()

#Function to register an external callback
cdef extern:
    void gfsRegisterPhysicsCallback(pyPhysicsCallback callback)

#Function to add u,v tendency to vrt,div tendency
cdef extern:
    void gfs_uv_to_vrtdiv(
        double *pyUg,
        double *pyVg,
        double *pyVrtg,
        double *pyDivg)

#Function to calculate tendencies within the fortran code (For testing)
cdef extern:
    void gfs_calculate_tendencies(
        double *pyVrtg,
        double *pyDivg,
        double *pyVirtTempg,
        double *pyPressGrid,
        double *pySurfPressure,
        double *pyTracerg,
        double t,
        double dt)

#Function to set tendencies within the fortran code, which will be called by the
# dynamical core after the semi-implicit step
cdef extern:
    void gfs_set_tendencies(
        double *pyVrtTend,
        double *pyDivTend,
        double *pyVirtTempTend,
        double *pyQTend,
        double *pyLnpsTend,
        double *pyTracerTend)

#Function to set model topography
cdef extern:
    void gfs_set_topography(double *surface_geopotential)

#Function to set initial atmospheric mass
cdef extern:
    void gfs_set_dry_pressure(double *dry_pressure)

#Function to deallocate arrays in physics, etc.,
cdef extern:
    void gfs_finalise()

#Function to convert the input grid arrays to spectral arrays
# The data must be available in py{Ug,Vg,VirtTempg,Tracerg,Lnpsg}
cdef extern:
    void gfs_convert_to_spectral()

#Function to convert the spectral arrays to grid arrays
# The result will be available in py{Ug,Vg,VirtTempg,Tracerg,Lnpsg}
cdef extern:
    void gfs_convert_to_grid()

#Function to initialise the arrays for computation
# These are allocated in Python and passed in
cdef extern:
    void gfs_initialise_spectral_arrays(
        complex *pyVrtSpec,
        complex *pyDivSpec,
        complex *pyVirtTempSpec,
        complex *pyTracerSpec,
        complex *pyTopoSpec,
        complex *pyLnPsSpec,
        complex *pyDissSpec,
        complex *pyDmpProf,
        complex *pyDiffProf)

cdef extern:
    void gfs_assign_grid_arrays(
        double *pyUg,
        double *pyVg,
        double *pyVirtTempg,
        double *pyLnPsg,
        double *pyTracerg,
        double *pyVrtg,
        double *pyDivg)

cdef extern:
    void gfs_initialise_grid_arrays(
        double *pyDlnpdtg,
        double *pyEtaDotg,
        double *pyPhisg,
        double *pyDPhisdxg,
        double *pyDPhisdyg,
        double *pyDlnpsdtg)

cdef extern:
    void gfs_assign_pressure_arrays(
        double *pySurfPressure,
        double *pyPressGrid,
        double *pyInterfacePressure)

cdef extern:
    void gfs_init_coordinate_arrays(
        double *ak, double *bk, double *ck,
        double *dbk, double *si, double *sl,
        double *bkl, double *alfa, double *rlnp, double *dpk)

cdef extern:
    void gfs_initialise_spectral_physics_arrays(
        complex *pyVrtSpecTend,
        complex *pyDivSpecTend,
        complex *pyVirtTempSpecTend,
        complex *pyTracerSpecTend,
        complex *pyLnPsSpec)

cdef extern:
    void gfs_initialise_grid_physics_arrays(
        double *pyPwat)

cdef extern:
    void gfs_initialise_timestepper_arrays(
        double *pyDmhyb, double *pyAmhyb, double *pyBmhyb,
        double *pyTref, double *pyPkref, double *pyDpkref,
        double *pyAlfaref, double *pySvhyb, double *pyTorhyb)

cdef extern:
    void gfs_assign_sht_arrays(
        double *pyLap, double *pyInvLap, double *pyGaussWeights,
        double *pyLons, double *pyLats, double *pyAreaWeights,
        int *pyDegree, int *pyOrder, int *pyNlons, int *pyNlats,
        int *pyNlm)

cdef void testFunc():
    print 'a'



# Grid space arrays, cython for initialisation
cdef cnp.double_t[::1,:,:,:] pyTracerg
cdef cnp.double_t[::1,:,:] pyUg, pyVg, pyVrtg, pyDivg,\
        pyVirtTempg, pyDlnpdtg, pyEtaDotg
cdef cnp.double_t[::1,:] pyLnPsg, pyPhis, pyDPhisdx,\
        pyDPhisdy, pyDlnpsdt, pyPwat

# Pressure arrays, in grid space
cdef cnp.double_t[::1,:] pySurfPressure
cdef cnp.double_t[::1,:,:] pyPressGrid
cdef cnp.double_t[::1,:,:] pyInterfacePressure

# Vertical Coordinate arrays
cdef cnp.double_t[::1] pyAk, pyBk, pyCk, pyDbk, pySi, pySl, pyBkl
cdef cnp.double_t[::1,:,:] pyAlfa, pyRlnp, pyDpk

# Timestepper arrays
cdef cnp.double_t[::1, :] pyAmhyb, pyBmhyb
cdef cnp.double_t[::1, :, :, :] pyDmhyb
cdef cnp.double_t[::1] pyTref, pyPkref, pyDpkref, pyAlfaref
cdef cnp.double_t[::1] pySvhyb, pyTorhyb

# Spectral arrays, using cython for declaration
cdef cnp.complex128_t[::1,:,:] pyTracerSpec

cdef cnp.complex128_t[::1,:] pyVrtSpec, pyDivSpec, pyVirtTempSpec

cdef cnp.complex128_t[:] pyTopoSpec, pyLnPsSpec,\
            pyDissSpec, pyDmpProf, pyDiffProf

cdef cnp.complex128_t[::1,:,:] pyTracerSpecTend

cdef cnp.complex128_t[::1,:] pyVrtSpecTend, pyDivSpecTend, pyVirtTempSpecTend

cdef cnp.complex128_t[:] pyLnPsSpecTend

# Temporary arrays for setting tendency terms
cdef cnp.double_t[::1,:,:] tempVrtTend, tempDivTend, tempVirtTempTend, tempUTend,tempVTend
cdef cnp.double_t[::1,:,:,:] tempTracerTend
cdef cnp.double_t[::1,:] tempLnpsTend

# Shtns arrays
cdef cnp.double_t[::1, :] pyLats, pyLons, pyAreaWeights
cdef cnp.double_t[::1] pyLap, pyInvLap, pyGaussWeights
cdef cnp.int_t[::1] pyDegree, pyOrder


# Grid size
#cdef public int numLats, numLons, numTrunc, numLevs, spectralDim, numTracers
cdef int nlm

# Model state
cdef int modelIsInitialised, physicsEnabled
cdef double zero_model_time=0

# Physics subroutine pointer
cdef object physicsCallback

#TODO get constant from physcons.f90
cdef double fv

def set_constants(double radius, double omega,
                  double R, double Rd, double Rv,
                  double g, double Cp, double C_cond):

    gfs_set_constants(&radius, &omega, &R, &Rd, &Rv, &g, &Cp, &C_cond)

def set_model_grid(int nlats, int nlons,
                   int nlevs, int ntrunc, int ndimspec,
                   int ntracers):

    global pyDmhyb, pyAmhyb, pyBmhyb,\
            pyTref, pyPkref, pyDpkref,\
            pyAlfaref, pySvhyb, pyTorhyb,\
            pyAk, pyBk, pyCk, pySi, pySl, pyDbk,\
            pyBkl, pyAlfa, pyRlnp, pyDpk,\
            pyLap, pyInvLap, pyGaussWeights, pyLats, pyLons, pyAreaWeights,\
            pyDegree, pyOrder, nlm


    nlm = int((ntrunc+1)*(ntrunc+2)/2)

    pyLap = np.zeros(nlm, dtype=np.double, order='F')
    pyInvLap = np.zeros(nlm, dtype=np.double, order='F')
    pyGaussWeights = np.zeros(nlats, dtype=np.double, order='F')

    pyDegree = np.zeros(nlm, dtype=np.int, order='F')
    pyOrder = np.zeros(nlm, dtype=np.int, order='F')

    pyLats = np.zeros((nlons, nlats), dtype=np.double, order='F')
    pyLons = np.zeros((nlons, nlats), dtype=np.double, order='F')
    pyAreaWeights = np.zeros((nlons, nlats), dtype=np.double, order='F')

    pyDmhyb = np.zeros((nlevs, nlevs, ntrunc+1, 3),
                       dtype=np.double, order='F')

    pyAmhyb = np.zeros((nlevs, nlevs), dtype=np.double, order='F')
    pyBmhyb = np.zeros((nlevs, nlevs), dtype=np.double, order='F')

    pyTref = np.zeros((nlevs), dtype=np.double, order='F')
    pyPkref = np.zeros((nlevs+1), dtype=np.double, order='F')
    pyDpkref = np.zeros((nlevs), dtype=np.double, order='F')
    pyAlfaref = np.zeros((nlevs), dtype=np.double, order='F')
    pySvhyb = np.zeros((nlevs), dtype=np.double, order='F')
    pyTorhyb = np.zeros((nlevs), dtype=np.double, order='F')

    pyAk = np.zeros((nlevs+1), dtype=np.double, order='F')
    pyBk = np.zeros((nlevs+1), dtype=np.double, order='F')
    pySi = np.zeros((nlevs+1), dtype=np.double, order='F')
    pyCk = np.zeros((nlevs), dtype=np.double, order='F')
    pyDbk = np.zeros((nlevs), dtype=np.double, order='F')
    pySl = np.zeros((nlevs), dtype=np.double, order='F')
    pyBkl = np.zeros((nlevs), dtype=np.double, order='F')

    pyAlfa = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
    pyRlnp = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
    pyDpk = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')

    gfs_set_model_dimensions(&nlats, &nlons, &nlevs, &ntrunc,
                             &ndimspec, &ntracers)

    gfs_initialise_timestepper_arrays(
        <double *>&pyDmhyb[0,0,0,0],
        <double *>&pyAmhyb[0,0],
        <double *>&pyBmhyb[0,0],
        <double *>&pyTref[0],
        <double *>&pyPkref[0],
        <double *>&pyDpkref[0],
        <double *>&pyAlfaref[0],
        <double *>&pySvhyb[0],
        <double *>&pyTorhyb[0])

    gfs_init_coordinate_arrays(
        <double *>&pyAk[0],
        <double *>&pyBk[0],
        <double *>&pyCk[0],
        <double *>&pyDbk[0],
        <double *>&pySi[0],
        <double *>&pySl[0],
        <double *>&pyBkl[0],
        <double *>&pyAlfa[0,0,0],
        <double *>&pyRlnp[0,0,0],
        <double *>&pyDpk[0,0,0])

    gfs_assign_sht_arrays(
        <double *>&pyLap[0],
        <double *>&pyInvLap[0],
        <double *>&pyGaussWeights[0],
        <double *>&pyLons[0, 0],
        <double *>&pyLats[0, 0],
        <double *>&pyAreaWeights[0, 0],
        <int *>&pyDegree[0],
        <int *>&pyOrder[0],
        &nlons, &nlats, &nlm)

    init_spectral_arrays(ndimspec, nlevs, ntracers)
    init_grid_arrays(nlons, nlats, nlevs, ntracers)

# Initialise dynamics and physics
def init_model(double dry_pressure, int num_damp_levels, double tau_damping):

    # Now that the arrays are initialised, call dynamics and physics init


    gfs_set_dry_pressure(&dry_pressure)
    gfs_set_model_time(&zero_model_time)
    gfs_init_dynamics(&num_damp_levels, &tau_damping)
    #gfs_get_lon_lat(<double *>&__longitudes[0,0], <double *>&__latitudes[0,0])
    longitudes = np.ascontiguousarray(pyLons).copy()
    latitudes = np.ascontiguousarray(pyLats).copy()
    area_weights = np.ascontiguousarray(pyAreaWeights).copy()
    gaussian_weights = np.ascontiguousarray(pyGaussWeights).copy()
    sigma_levels = np.ascontiguousarray(pySl).copy()
    sigma_inteface_levels = np.ascontiguousarray(pySi).copy()

    return gaussian_weights, area_weights, latitudes, longitudes, sigma_levels, sigma_inteface_levels



# Create the spectral arrays (defined in spectral_data.f90)
def init_spectral_arrays(spectral_dim, num_levs, num_tracers):

    global pyTracerSpec, pyTracerSpecTend, pyTopoSpec, \
            pyLnPsSpec, pyLnPsSpecTend,\
            pyDissSpec, pyDmpProf, pyDiffProf,\
            pyVrtSpec, pyVrtSpecTend, pyDivSpec, pyDivSpecTend,\
            pyVirtTempSpec, pyVirtTempSpecTend

    pyTracerSpec = np.zeros(
        (spectral_dim, num_levs, num_tracers),dtype=complex, order='F')

    pyVrtSpec = np.zeros((spectral_dim, num_levs),dtype=complex, order='F')
    pyVirtTempSpec = np.zeros(
        (spectral_dim, num_levs),dtype=complex, order='F')

    pyDivSpec = np.zeros((spectral_dim, num_levs),dtype=complex,order='F')


    pyTopoSpec = np.zeros(spectral_dim,dtype=complex,order='F')
    pyLnPsSpec = np.zeros(spectral_dim,dtype=complex,order='F')
    pyDissSpec = np.zeros(spectral_dim,dtype=complex,order='F')

    pyDmpProf = np.zeros(num_levs,dtype=complex,order='F')
    pyDiffProf = np.zeros(num_levs,dtype=complex,order='F')

    pyTracerSpecTend = np.zeros(
        (spectral_dim, num_levs, num_tracers),dtype=complex, order='F')

    pyVrtSpecTend = np.zeros(
        (spectral_dim, num_levs),dtype=complex, order='F')

    pyVirtTempSpecTend = np.zeros(
        (spectral_dim, num_levs),dtype=complex, order='F')

    pyDivSpecTend = np.zeros(
        (spectral_dim, num_levs),dtype=complex,order='F')

    pyLnPsSpecTend = np.zeros(spectral_dim,dtype=complex,order='F')

    gfs_initialise_spectral_arrays(
        <double complex *>&pyVrtSpec[0,0],
        <double complex *>&pyDivSpec[0,0],
        <double complex *>&pyVirtTempSpec[0,0],
        <double complex *>&pyTracerSpec[0,0,0],
        <double complex *>&pyTopoSpec[0],
        <double complex *>&pyLnPsSpec[0],
        <double complex *>&pyDissSpec[0],
        <double complex *>&pyDmpProf[0],
        <double complex *>&pyDiffProf[0])

    gfs_initialise_spectral_physics_arrays(
        <double complex *>&pyVrtSpecTend[0,0],
        <double complex *>&pyDivSpecTend[0,0],
        <double complex *>&pyVirtTempSpecTend[0,0],
        <double complex *>&pyTracerSpecTend[0,0,0],
        <double complex *>&pyLnPsSpecTend[0])



# Create the grid arrays (defined in grid_data.f90)
def init_grid_arrays(num_lons, num_lats, num_levs, num_tracers):

    global pyDlnpdtg, pyEtaDotg, pyPhis, pyDPhisdx, pyDPhisdy, \
            pyDlnpsdt, pyPwat, tempVrtTend, tempDivTend

    pyDlnpdtg = np.zeros((num_lons, num_lats, num_levs),
                         dtype=np.double, order='F')
    pyEtaDotg = np.zeros((num_lons, num_lats, num_levs+1),
                         dtype=np.double, order='F')


    pyPhis = np.zeros((num_lons, num_lats), dtype=np.double, order='F')
    pyDPhisdx = np.zeros((num_lons, num_lats), dtype=np.double, order='F')
    pyDPhisdy = np.zeros((num_lons, num_lats), dtype=np.double, order='F')
    pyDlnpsdt = np.zeros((num_lons, num_lats), dtype=np.double, order='F')
    pyPwat = np.zeros((num_lons, num_lats), dtype=np.double, order='F')



    tempVrtTend = np.zeros((num_lons, num_lats, num_levs),
                           dtype=np.double, order='F')
    tempDivTend = np.zeros((num_lons, num_lats, num_levs),
                           dtype=np.double, order='F')


    gfs_initialise_grid_arrays(
        <double *>&pyDlnpdtg[0,0,0],
        <double *>&pyEtaDotg[0,0,0],
        <double *>&pyPhis[0,0],
        <double *>&pyDPhisdx[0,0],
        <double *>&pyDPhisdy[0,0],
        <double *>&pyDlnpsdt[0,0])

    gfs_initialise_grid_physics_arrays(
        <double *>&pyPwat[0,0])

#Intialise pressure arrays
def assign_pressure_arrays(
    cnp.double_t[::1, :] surface_pressure,
    cnp.double_t[::1, :, :] midlevel_pressure,
    cnp.double_t[::1, :, :] interface_pressure):


    gfs_assign_pressure_arrays(
        <double *>&surface_pressure[0,0],
        <double *>&midlevel_pressure[0,0,0],
        <double *>&interface_pressure[0,0,0])

# Set tendencies for dynamical core to use in physics
def assign_tendencies(
    cnp.double_t[::1, :, :] u_tend,
    cnp.double_t[::1, :, :] v_tend,
    cnp.double_t[::1, :, :] virtemp_tend,
    cnp.double_t[::1, :, :] q_tend,
    cnp.double_t[::1, :] lnps_tend,
    cnp.double_t[::1, :, :, :] tracer_tend):

    global tempVrtTend, tempDivTend

    gfs_uv_to_vrtdiv(
        <double *>&u_tend[0,0,0],
        <double *>&v_tend[0,0,0],
        <double *>&tempVrtTend[0,0,0],
        <double *>&tempDivTend[0,0,0])

    gfs_set_tendencies(
        <double *>&tempVrtTend[0,0,0],
        <double *>&tempDivTend[0,0,0],
        <double *>&virtemp_tend[0,0,0],
        <double *>&q_tend[0,0,0],
        <double *>&lnps_tend[0,0],
        <double *>&tracer_tend[0,0,0,0])

'''
# Take one step
def oneStepForward():

    gfs_convert_to_grid()
    gfs_calc_pressure()


    if physicsEnabled:
        #TODO don't call physics callback directly. use a helper function which will remove
        #TODO individual fields from tracerg and assign them to q, ozone, etc., and then
        #TODO call physics routines
        #physicsCallback(pyUg,
        #                     pyVg,
        tendList = physicsCallback(
            pyUg,
            pyVg,
            pyVirtTempg,
            pyPressGrid,
            pySurfPressure,
            pyTracerg,
            __latitudes)

        setTendencies(tendList)

    else:
        calculate_tendencies(
            <double *>&pyVrtg[0,0,0],
            <double *>&pyDivg[0,0,0],
            <double *>&pyVirtTempg[0,0,0],
            <double *>&pyPressGrid[0,0,0],
            <double *>&pySurfPressure[0,0],
            <double *>0,
            t,
            dt)

        gfs_take_one_step()

'''
# Register a callback which calculates the physics (to be used in stand-alone
# mode only)

cdef setPhysicsCallback(physicsFnPtr):

    physicsCallback = physicsFnPtr
    physicsEnabled = True
    #        gfsRegisterPhysicsCallback(testFunc)


def assign_grid_arrays(
    cnp.double_t[::1, :, :] u,
    cnp.double_t[::1, :, :] v,
    cnp.double_t[::1, :, :] tvirt,
    cnp.double_t[::1, :] lnpsg,
    cnp.double_t[::1, :, :, :] tracers,
    cnp.double_t[::1, :, :] vorticity,
    cnp.double_t[::1, :, :] divergence):

    gfs_assign_grid_arrays(
        <double *>&u[0,0,0],
        <double *>&v[0,0,0],
        <double *>&tvirt[0,0,0],
        <double *>&lnpsg[0,0],
        <double *>&tracers[0,0,0,0],
        <double *>&vorticity[0,0,0],
        <double *>&divergence[0,0,0])

'''
def getResult():

    gfs_convert_to_grid()
    gfs_calc_pressure()

    outputList = []

    theta = np.asfortranarray(pyVirtTempg)
    q = np.asfortranarray(pyTracerg[:,:,:,0])
    q[q<0] = 0

    temp = theta.copy()
    if(np.any(q > 0)):
        temp = theta/(1+fv*q)

    outputList.append(np.asarray(pyUg).copy(order='F'))
    outputList.append(np.asarray(pyVg).copy(order='F'))
    outputList.append(temp)
    outputList.append(q)
    outputList.append(np.asarray(pySurfPressure).copy(order='F'))
    outputList.append(np.asarray(pyPressGrid).copy(order='F'))
    iface_press = np.zeros(np.asarray(pyInterfacePressure).shape, dtype=np.double, order='F')
    iface_press[:] = np.asfortranarray(pyInterfacePressure)[:,:,::-1]
    outputList.append(iface_press)

    return(outputList)
'''

def update_spectral_arrays():

    gfs_convert_to_spectral()


'''
def integrateFields(field_list, increment_list):
    #Only to be used in CLIMT mode

    if climt_mode:

        q = np.asfortranarray(pyTracerg[:,:,:,0])
        q[q<0] = 0

        temptrac = np.zeros((nlons,nlats,nlevs,ntrac),dtype=np.double,order='F')
        uTend,vTend, tempTend,qTend,psTend = increment_list

        virtTempTend = tempTend
        temperature = pyVirtTempg/(1+fv*q)

        #Additional tendency term while converting to virtual temp -- see Pg. 12 in doc/gfsModelDoc.pdf
        if(np.any(q > 0)):
            virtTempTend = tempTend*(1+fv*q) + fv*temperature*qTend

        psTend[psTend==0] = 1.
        lnpsTend = np.log(psTend)

        u,v,virtemp,q,ps = field_list


        #CliMT gives increments; convert to tendencies
        uTend /= dt

        vTend /= dt

        virtTempTend /= dt

        qTend /= dt

        lnpsTend /= dt

        temptrac[:,:,:,0] = qTend

        increment_list = uTend,vTend,virtTempTend,lnpsTend,temptrac
        setTendencies(increment_list)


        #Step forward in time
        take_one_step()

        #Convert back to grid space
        convert_to_grid()

        # only ln(Ps) is calculated in the dynamics. This calculates
        # the values on the full grid
        calculate_pressure()

        ug = np.asfortranarray(pyUg.copy())
        vg = np.asfortranarray(pyVg.copy())
        virtempg = np.asfortranarray(pyVirtTempg.copy())
        qg = np.asfortranarray(pyTracerg[:,:,:,0].copy())
        psg = np.asfortranarray(pySurfPressure.copy())
        press = np.asfortranarray(pyPressGrid.copy())
        iface_press = np.zeros(np.array(pyInterfacePressure).shape, dtype=np.double, order='F')
        iface_press[:] = np.asfortranarray(pyInterfacePressure.copy())[:,:,::-1]

        qg[qg<0] = 0

        if(np.any(qg > 0)):
            #output temperature, not virtual temp
            virtempg = virtempg/(1+fv*qg)

        return(ug,vg,virtempg,qg,psg,press,iface_press)
'''

def take_one_step():
    gfs_take_one_step()

def convert_to_grid():
    gfs_convert_to_grid()

def calculate_pressure():
    gfs_calculate_pressure()

def set_topography(cnp.double_t[::1, :] topography):

    gfs_set_topography(<double *>&topography[0,0])

def set_time_step(double time_step):

    gfs_set_time_step(&time_step)

def shut_down_model():
    #Remember to set time to zero!
    gfs_finalise()
    gfs_set_model_time(&zero_model_time)
