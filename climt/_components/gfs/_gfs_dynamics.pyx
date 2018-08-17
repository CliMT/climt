cimport numpy as cnp
import numpy as np
import cython

# Typedef for function pointer returning void and taking no arguments (for now)

ctypedef void (*pyPhysicsCallback)(double *, double *, double *, double *,
                                   double *, double *, double *,
                              double *, double *, double *, double *,
                                   double *, double *, double *)

# Variables that control the behaviour of the model
# bint is Cython's boolean type
cdef extern:
    bint dry, adiabatic, heldsuarez

cdef extern:
    void gfs_lnps_to_spectral()
    void gfs_vrt_div_to_spectral()
    void gfs_virtemp_to_spectral()
    void gfs_tracer_to_spectral()

#Function to set global mass of dry air
cdef extern:
    void gfs_set_dry_pressure(double *pdry_init)

#Function to set top of model pressure
cdef extern:
    void gfs_set_model_top_pressure(double *p_toa)

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
#cdef extern:
#    void gfsRegisterPhysicsCallback(pyPhysicsCallback callback)

#Function to add u,v tendency to vrt,div tendency
cdef extern:
    void gfs_uv_to_vrtdiv(
        double *pyUg,
        double *pyVg,
        double *pyVrtg,
        double *pyDivg)

#Function to set tendencies within the fortran code, which will be called by the
# dynamical core after the semi-implicit step
cdef extern:
    void gfs_set_tendencies(
        double *pyVrtTend,
        double *pyDivTend,
        double *pyVirtTempTend,
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


# Grid space arrays, cython for initialisation
cdef cnp.double_t[:,:,:,:] pyTracerg
cdef cnp.double_t[:,:,:] pyUg, pyVg, pyVrtg, pyDivg,\
        pyVirtTempg, pyDlnpdtg, pyEtaDotg
cdef cnp.double_t[:,:] pyLnPsg, pyPhis, pyDPhisdx,\
        pyDPhisdy, pyDlnpsdt, pyPwat

# Pressure arrays, in grid space
cdef cnp.double_t[:,:] pySurfPressure
cdef cnp.double_t[:,:,:] pyPressGrid
cdef cnp.double_t[:,:,:] pyInterfacePressure

# Vertical Coordinate arrays
cdef cnp.double_t[:] pyAk, pyBk, pyCk, pyDbk, pySi, pySl, pyBkl
cdef cnp.double_t[:,:,:] pyAlfa, pyRlnp, pyDpk

# Timestepper arrays
cdef cnp.double_t[:, :] pyAmhyb, pyBmhyb
cdef cnp.double_t[:, :, :, :] pyDmhyb
cdef cnp.double_t[:] pyTref, pyPkref, pyDpkref, pyAlfaref
cdef cnp.double_t[:] pySvhyb, pyTorhyb

# Spectral arrays, using cython for declaration
cdef cnp.complex128_t[:,:,:] pyTracerSpec

cdef cnp.complex128_t[:,:] pyVrtSpec, pyDivSpec, pyVirtTempSpec

cdef cnp.complex128_t[:] pyTopoSpec, pyLnPsSpec,\
            pyDissSpec, pyDmpProf, pyDiffProf

cdef cnp.complex128_t[:,:,:] pyTracerSpecTend

cdef cnp.complex128_t[:,:] pyVrtSpecTend, pyDivSpecTend, pyVirtTempSpecTend

cdef cnp.complex128_t[:] pyLnPsSpecTend

# Temporary arrays for setting tendency terms
cdef cnp.double_t[:,:,:] tempVrtTend, tempDivTend, tempVirtTempTend, tempUTend,tempVTend
cdef cnp.double_t[:,:,:,:] tempTracerTend
cdef cnp.double_t[:,:] tempLnpsTend

# Shtns arrays
cdef cnp.double_t[:, :] pyLats, pyLons, pyAreaWeights
cdef cnp.double_t[:] pyLap, pyInvLap, pyGaussWeights
cdef cnp.int_t[:] pyDegree, pyOrder


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

def lnps_to_spectral():
    gfs_lnps_to_spectral()

def vrt_div_to_spectral():
    gfs_vrt_div_to_spectral()

def virtemp_to_spectral():
    gfs_virtemp_to_spectral()

def tracer_to_spectral():
    gfs_tracer_to_spectral()

def set_constants(double radius, double omega,
                  double R, double Rd, double Rv,
                  double g, double Cp, double C_cond):

    gfs_set_constants(&radius, &omega, &R, &Rd, &Rv, &g, &Cp, &C_cond)

#@cython.boundscheck(False)
def set_model_grid(
        int nlats, int nlons,
        int nlevs, int ntrunc, int ndimspec,
        int ntracers,
        cnp.double_t[:] ak_in,
        cnp.double_t[:] bk_in,
    ):

    global pyDmhyb, pyAmhyb, pyBmhyb,\
            pyTref, pyPkref, pyDpkref,\
            pyAlfaref, pySvhyb, pyTorhyb,\
            pyAk, pyBk, pyCk, pySi, pySl, pyDbk,\
            pyBkl, pyAlfa, pyRlnp, pyDpk,\
            pyLap, pyInvLap, pyGaussWeights, pyLats, pyLons, pyAreaWeights,\
            pyDegree, pyOrder, nlm


    nlm = int((ntrunc+1)*(ntrunc+2)/2)

    pyLap = np.zeros(nlm, dtype=np.double)
    pyInvLap = np.zeros(nlm, dtype=np.double)
    pyGaussWeights = np.zeros(nlats, dtype=np.double)

    pyDegree = np.zeros(nlm, dtype=np.int)
    pyOrder = np.zeros(nlm, dtype=np.int)

    pyLats = np.zeros((nlats, nlons), dtype=np.double)
    pyLons = np.zeros((nlats, nlons), dtype=np.double)
    pyAreaWeights = np.zeros((nlats, nlons), dtype=np.double)

    pyDmhyb = np.zeros((3, ntrunc+1, nlevs, nlevs),
                       dtype=np.double)

    pyAmhyb = np.zeros((nlevs, nlevs), dtype=np.double)
    pyBmhyb = np.zeros((nlevs, nlevs), dtype=np.double)

    pyTref = np.zeros((nlevs), dtype=np.double)
    pyPkref = np.zeros((nlevs+1), dtype=np.double)
    pyDpkref = np.zeros((nlevs), dtype=np.double)
    pyAlfaref = np.zeros((nlevs), dtype=np.double)
    pySvhyb = np.zeros((nlevs), dtype=np.double)
    pyTorhyb = np.zeros((nlevs), dtype=np.double)

    pyAk = np.zeros((nlevs+1), dtype=np.double)
    pyBk = np.zeros((nlevs+1), dtype=np.double)
    cdef int i
    for i in range(nlevs+1):
        pyAk[i] = ak_in[i]
        pyBk[i] = bk_in[i]
    pySi = np.zeros((nlevs+1), dtype=np.double)
    pyCk = np.zeros((nlevs), dtype=np.double)
    pyDbk = np.zeros((nlevs), dtype=np.double)
    pySl = np.zeros((nlevs), dtype=np.double)
    pyBkl = np.zeros((nlevs), dtype=np.double)

    pyAlfa = np.zeros((nlevs, nlats, nlons), dtype=np.double)
    pyRlnp = np.zeros((nlevs, nlats, nlons), dtype=np.double)
    pyDpk = np.zeros((nlevs, nlats, nlons), dtype=np.double)

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

    init_grid_arrays(nlons, nlats, nlevs, ntracers)
    init_spectral_arrays(ndimspec, nlevs, ntracers)


# Initialise dynamics and physics
def init_model(
    double dry_pressure, int num_damp_levels, double tau_damping,
    double model_top_pressure):

    # Now that the arrays are initialised, call dynamics and physics init


    gfs_set_dry_pressure(&dry_pressure)
    gfs_set_model_top_pressure(&model_top_pressure)
    gfs_set_model_time(&zero_model_time)
    gfs_init_dynamics(&num_damp_levels, &tau_damping)
    #gfs_get_lon_lat(<double *>&__longitudes[0,0], <double *>&__latitudes[0,0])

    longitudes = np.ascontiguousarray(pyLons).copy()
    latitudes = np.ascontiguousarray(pyLats).copy()
    area_weights = np.ascontiguousarray(pyAreaWeights).copy()
    gaussian_weights = np.ascontiguousarray(pyGaussWeights).copy()
    # sigma_levels = np.ascontiguousarray(pySl).copy()
    # sigma_inteface_levels = np.ascontiguousarray(pySi).copy()

    return gaussian_weights, area_weights, latitudes, longitudes#, sigma_levels, sigma_interface_levels


# Create the spectral arrays (defined in spectral_data.f90)
@cython.boundscheck(False)
def init_spectral_arrays(spectral_dim, num_levs, num_tracers):

    global pyTracerSpec, pyTracerSpecTend, pyTopoSpec, \
            pyLnPsSpec, pyLnPsSpecTend,\
            pyDissSpec, pyDmpProf, pyDiffProf,\
            pyVrtSpec, pyVrtSpecTend, pyDivSpec, pyDivSpecTend,\
            pyVirtTempSpec, pyVirtTempSpecTend

    pyTracerSpec = np.zeros(
        (num_tracers, num_levs, spectral_dim), dtype=complex)

    pyVrtSpec = np.zeros((num_levs, spectral_dim), dtype=complex)
    pyVirtTempSpec = np.zeros(
        (spectral_dim, num_levs), dtype=complex)

    pyDivSpec = np.zeros((num_levs, spectral_dim), dtype=complex)


    pyTopoSpec = np.zeros(spectral_dim, dtype=complex)
    pyLnPsSpec = np.zeros(spectral_dim, dtype=complex)
    pyDissSpec = np.zeros(spectral_dim, dtype=complex)

    pyDmpProf = np.zeros(num_levs, dtype=complex)
    pyDiffProf = np.zeros(num_levs, dtype=complex)

    pyTracerSpecTend = np.zeros(
        (num_tracers, num_levs, spectral_dim), dtype=complex)

    pyVrtSpecTend = np.zeros(
        (num_levs, spectral_dim), dtype=complex)

    pyVirtTempSpecTend = np.zeros(
        (num_levs, spectral_dim), dtype=complex)

    pyDivSpecTend = np.zeros(
        (num_levs, spectral_dim), dtype=complex)

    pyLnPsSpecTend = np.zeros(spectral_dim, dtype=complex)

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

    pyDlnpdtg = np.zeros((num_levs, num_lats, num_lons),
                         dtype=np.double)
    pyEtaDotg = np.zeros((num_levs+1, num_lats, num_lons),
                         dtype=np.double)


    pyPhis = np.zeros((num_lats, num_lons), dtype=np.double)
    pyDPhisdx = np.zeros((num_lats, num_lons), dtype=np.double)
    pyDPhisdy = np.zeros((num_lats, num_lons), dtype=np.double)
    pyDlnpsdt = np.zeros((num_lats, num_lons), dtype=np.double)
    pyPwat = np.zeros((num_lats, num_lons), dtype=np.double)



    tempVrtTend = np.zeros((num_levs, num_lats, num_lons),
                           dtype=np.double)
    tempDivTend = np.zeros((num_levs, num_lats, num_lons),
                           dtype=np.double)


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
    cnp.double_t[:, :] surface_pressure,
    cnp.double_t[:, :, :] midlevel_pressure,
    cnp.double_t[:, :, :] interface_pressure):


    gfs_assign_pressure_arrays(
        <double *>&surface_pressure[0,0],
        <double *>&midlevel_pressure[0,0,0],
        <double *>&interface_pressure[0,0,0])

# Set tendencies for dynamical core to use in physics
@cython.boundscheck(False)
def assign_tendencies(
    cnp.double_t[:, :, :] u_tend,
    cnp.double_t[:, :, :] v_tend,
    cnp.double_t[:, :, :] virtemp_tend,
    cnp.double_t[:, :] lnps_tend,
    cnp.double_t[:, :, :, :] tracer_tend):

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
        <double *>&lnps_tend[0,0],
        <double *>&tracer_tend[0,0,0,0])

# Register a callback which calculates the physics (to be used in stand-alone
# mode only)

#cdef setPhysicsCallback(physicsFnPtr):

#    physicsCallback = physicsFnPtr
#    physicsEnabled = True
    #        gfsRegisterPhysicsCallback(testFunc)


@cython.boundscheck(False)
def assign_grid_arrays(
    cnp.double_t[:, :, :] u,
    cnp.double_t[:, :, :] v,
    cnp.double_t[:, :, :] tvirt,
    cnp.double_t[:, :] lnpsg,
    cnp.double_t[:, :, :, :] tracers,
    cnp.double_t[:, :, :] vorticity,
    cnp.double_t[:, :, :] divergence):

    gfs_assign_grid_arrays(
        <double *>&u[0,0,0],
        <double *>&v[0,0,0],
        <double *>&tvirt[0,0,0],
        <double *>&lnpsg[0,0],
        <double *>&tracers[0,0,0,0],
        <double *>&vorticity[0,0,0],
        <double *>&divergence[0,0,0])


def update_spectral_arrays():
    gfs_convert_to_spectral()


def take_one_step():
    gfs_take_one_step()

def convert_to_grid():
    gfs_convert_to_grid()

def calculate_pressure():
    gfs_calculate_pressure()

def set_topography(cnp.double_t[:, :] topography):

    gfs_set_topography(<double *>&topography[0,0])

def set_time_step(double time_step):

    gfs_set_time_step(&time_step)

def shut_down_model():
    #Remember to set time to zero!
    gfs_finalise()
    gfs_set_model_time(&zero_model_time)
