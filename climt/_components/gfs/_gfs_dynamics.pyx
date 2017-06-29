cimport numpy as cnp
import numpy as np
from numpy import empty


# Typedef for function pointer returning void and taking no arguments (for now)

ctypedef void (*pyPhysicsCallback)(double *, double *, double *, double *,\
                                   double *, double *, double *,\
                              double *, double *, double *, double *,\
                                   double *, double *, double *)

# Variables for grid sizes from the library
cdef extern:
    int nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

# Variables for simulation time control
cdef extern:
    int fhmax, fhout

# time step and full time from library
cdef extern:
    double dt, deltim, t

# Variables that control the behaviour of the model
# bint is Cython's boolean type
cdef extern:
    bint dry, adiabatic, heldsuarez

#Variable for initial pressure
cdef extern:
    double pdryini

# Function definitions to do our work
cdef extern:
    void gfsReadNamelist()

#Function to get longitude/latitude map from shtns
cdef extern:
    void get_lon_lat(double *longitudes, double *latitudes)

# Function to init dynamics
cdef extern:
    void gfsInitDynamics()

# Function to init physics
cdef extern:
    void gfsInitPhysics()

#Function to calculate pressure fields after other fields have been updated
cdef extern:
    void gfsCalcPressure()

# Function to step fields by one dt
cdef extern:
    void gfsTakeOneStep()

#Function to register an external callback
cdef extern:
    void gfsRegisterPhysicsCallback(pyPhysicsCallback callback)

#Function to add u,v tendency to vrt,div tendency
cdef extern:
    void gfs_uv_to_vrtdiv(\
                          double *pyUg,\
                          double *pyVg,\
                          double *pyVrtg,\
                          double *pyDivg)

#Function to calculate tendencies within the fortran code (For testing)
cdef extern:
    void calculate_tendencies(\
                              double *pyVrtg,\
                              double *pyDivg,\
                              double *pyVirtTempg,\
                              double *pyPressGrid,\
                              double *pySurfPressure,\
                              double *pyTracerg,\
                              double t,\
                              double dt)

#Function to set tendencies within the fortran code, which will be called by the 
# dynamical core after the semi-implicit step
cdef extern:
    void set_tendencies(\
                              double *pyVrtTend,\
                              double *pyDivTend,\
                              double *pyVirtTempTend,\
                              double *pyTracerTend,\
                              double *pyLnpsTend,\
                              double t,\
                              double dt)


cdef extern:
    void gfsSetTopography(double *surface_geopotential)
#Function to deallocate arrays in physics, etc.,

cdef extern:
    void gfsFinalise()

#Function to convert the input grid arrays to spectral arrays
# The data must be available in py{Ug,Vg,VirtTempg,Tracerg,Lnpsg}
cdef extern:
    void gfsConvertToSpec()

#Function to convert the spectral arrays to grid arrays
# The result will be available in py{Ug,Vg,VirtTempg,Tracerg,Lnpsg}
cdef extern:
    void gfsConvertToGrid()

#Function to initialise the arrays for computation
# These are allocated in Python and passed in
cdef extern: 
    void initialiseSpectralArrays(\
        complex *pyVrtSpec, \
        complex *pyDivSpec,\
        complex *pyVirtTempSpec,\
        complex *pyTracerSpec,\
        complex *pyTopoSpec,\
        complex *pyLnPsSpec,\
        complex *pyDissSpec,\
        complex *pyDmpProf,\
        complex *pyDiffProf)

cdef extern: 
    void initialiseGridArrays(\
        double *pyUg,\
        double *pyVg,\
        double *pyVrtg,\
        double *pyDivg,\
        double *pyVirtTempg,\
        double *pyTracerg,\
        double *pyDlnpdtg,\
        double *pyEtaDotg,\
        double *pyLnPsg,\
        double *pyPhis,\
        double *pyDPhisdx,\
        double *pyDPhisdy,\
        double *pyDlnpsdt)

cdef extern:
    void initialisePressureArrays(\
        double *pySurfPressure,\
        double *pyPressGrid,\
        double *pyInterfacePressure)

cdef void testFunc():
    print 'a'


cdef class _gfs_dynamics:

# Grid space arrays, cython for initialisation
    cdef public cnp.double_t[::1,:,:,:] pyTracerg
    cdef public cnp.double_t[::1,:,:] pyUg, pyVg, pyVrtg, pyDivg,\
               pyVirtTempg, pyDlnpdtg, pyEtaDotg
    cdef public cnp.double_t[::1,:] pyLnPsg, pyPhis, pyDPhisdx,\
               pyDPhisdy, pyDlnpsdt
    
# Pressure arrays, in grid space
    cdef public cnp.double_t[::1,:] pySurfPressure
    cdef public cnp.double_t[::1,:,:] pyPressGrid
    cdef public cnp.double_t[::1,:,:] pyInterfacePressure


# Spectral arrays, using cython for declaration
    cdef public cnp.complex128_t[::1,:,:] pyTracerSpec

    cdef public cnp.complex128_t[::1,:] pyVrtSpec, pyDivSpec, pyVirtTempSpec
    
    cdef public cnp.complex128_t[:] pyTopoSpec, pyLnPsSpec, \
                pyDissSpec, pyDmpProf, pyDiffProf

# Arrays for storing lon/lat
    cdef public cnp.ndarray longitudes, latitudes

# Temporary arrays for setting tendency terms
    cdef cnp.double_t[::1,:,:] \
        tempVrtTend, tempDivTend, tempVirtTempTend, tempUTend,tempVTend
    cdef cnp.double_t[::1,:,:,:] tempTracerTend
    cdef cnp.double_t[::1,:] tempLnpsTend, __latitudes, __longitudes

# Grid size
    cdef public int numLats, numLons, numTrunc, numLevs, spectralDim, numTracers

# Model state
    cdef int modelIsInitialised, physicsEnabled

# Physics subroutine pointer
    cdef object physicsCallback

# Are we running inside CliMT?
    cdef int climt_mode

#TODO get constant from physcons.f90
    cdef double fv

    def __init__(self, numLons=192, numLats=94, \
                    simTimeHours=24, timestep=1200.0,\
                    useHeldSuarez=True, dryPressure=1e5,\
                    numTracers=1,physics=None,climt_mode=False):

        global adiabatic, dry, nlats, nlons, nlevs,\
            ntrunc, ndimspec, ntrac, fhmax, deltim,\
            heldsuarez, dt, pdryini, fhout, lats
# Read in the namelist
#        gfsReadNamelist()


# Read in the grid sizes (mainly to improve readability)
        if(numLats):
            nlats = <int>numLats
        
        if(numLons):
            nlons = <int>numLons
            ntrunc = <int>numLons/3 - 2
            ndimspec = (ntrunc+1)*(ntrunc+2)/2

        if(timestep):
            deltim = <double>timestep
            dt = <double>timestep

        if(physics):
            self.physicsCallback = physics
            self.physicsEnabled = True

        self.climt_mode = climt_mode

        self.fv = (4.6150e+2/2.8705e+2) - 1

        nlevs = 28
        pdryini = <double> dryPressure
        ntrac = <int>numTracers
        heldsuarez = useHeldSuarez
        fhmax = 9600
        fhout = 24
        adiabatic = False
        dry = True

        self.numLats = nlats
        self.numLons = nlons
        self.numLevs = nlevs
        self.numTrunc = ntrunc
        self.spectralDim = ndimspec
        self.numTracers = ntrac

        print 'Lats, lons, levs, trunc, dims, tracers', nlats, nlons, nlevs,\
            ntrunc, ndimspec, ntrac


# method to reconfigure model after instantiation. HAVE to call init model
# afterwards            
    def configureModel(self, numLons=None, numLats=None):
        
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

        if(numLats):
            self.numLats = <int>numLats
            nlats = self.numLats

        if(numLons):
            print self.numLons
            self.numLons = <int>numLons
            nlons = self.numLons
            self.numTrunc = self.numLons/3
            self.spectralDim = (self.numTrunc+1)*(self.numTrunc+2)/2
            
            ntrunc = self.numTrunc
            ndimspec = self.spectralDim

        print 'Current Lats, lons, trunc, dims', nlats, nlons, ntrunc, ndimspec
        print 'Current Lats, lons, trunc, dims', self.numLats, self.numLons,\
        self.numTrunc, self.spectralDim

# Initialise arrays and dynamics and physics
    def initModel(self):


        if(self.modelIsInitialised):
            self.shutDownModel()

        self.initSpectralArrays()
        self.initGridArrays()
        self.initPressureArrays()

# Now that the arrays are initialised, call dynamics and physics init

        gfsInitDynamics()
        gfsInitPhysics()
        get_lon_lat(<double *>&self.__longitudes[0,0], <double *>&self.__latitudes[0,0])
        self.longitudes = np.asfortranarray(self.__longitudes)
        self.latitudes = np.asfortranarray(self.__latitudes)

        self.modelIsInitialised = 1

# Create the spectral arrays (defined in spectral_data.f90)

    def initSpectralArrays(self):
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

        self.pyTracerSpec = np.zeros((ndimspec, nlevs, ntrac),dtype=complex, order='F')

        self.pyVrtSpec = np.zeros((ndimspec, nlevs),dtype=complex, order='F')
        self.pyVirtTempSpec = np.zeros((ndimspec, nlevs),dtype=complex, order='F')
        self.pyDivSpec = np.zeros((ndimspec, nlevs),dtype=complex,order='F')

        
        self.pyTopoSpec = np.zeros(ndimspec,dtype=complex,order='F')
        self.pyLnPsSpec = np.zeros(ndimspec,dtype=complex,order='F')
        self.pyDissSpec = np.zeros(ndimspec,dtype=complex,order='F')
        
        self.pyDmpProf = np.zeros(nlevs,dtype=complex,order='F')
        self.pyDiffProf = np.zeros(nlevs,dtype=complex,order='F')

        if(ntrac > 0):
            initialiseSpectralArrays(\
                 <double complex *>&self.pyVrtSpec[0,0], \
                 <double complex *>&self.pyDivSpec[0,0],\
                 <double complex *>&self.pyVirtTempSpec[0,0],\
                 <double complex *>&self.pyTracerSpec[0,0,0],\
                 <double complex *>&self.pyTopoSpec[0],\
                 <double complex *>&self.pyLnPsSpec[0],\
                 <double complex *>&self.pyDissSpec[0],\
                 <double complex *>&self.pyDmpProf[0],\
                 <double complex *>&self.pyDiffProf[0])
        else:
            initialiseSpectralArrays(\
                 <double complex *>&self.pyVrtSpec[0,0], \
                 <double complex *>&self.pyDivSpec[0,0],\
                 <double complex *>&self.pyVirtTempSpec[0,0],\
                 <double complex *>0,\
                 <double complex *>&self.pyTopoSpec[0],\
                 <double complex *>&self.pyLnPsSpec[0],\
                 <double complex *>&self.pyDissSpec[0],\
                 <double complex *>&self.pyDmpProf[0],\
                 <double complex *>&self.pyDiffProf[0])


# Create the grid arrays (defined in grid_data.f90)

    def initGridArrays(self):
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac


        self.pyTracerg = np.zeros((nlons, nlats, nlevs, ntrac),\
                dtype=np.double, order='F')

        self.pyUg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyVg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyVrtg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyDivg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyVirtTempg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyDlnpdtg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyEtaDotg = np.zeros((nlons, nlats, nlevs+1), dtype=np.double, order='F')


        self.pyLnPsg = np.zeros((nlons, nlats), dtype=np.double, order='F')
        self.pyPhis = np.zeros((nlons, nlats), dtype=np.double, order='F')
        self.pyDPhisdx = np.zeros((nlons, nlats), dtype=np.double, order='F')
        self.pyDPhisdy = np.zeros((nlons, nlats), dtype=np.double, order='F')
        self.pyDlnpsdt = np.zeros((nlons, nlats), dtype=np.double, order='F')


        self.tempTracerTend = np.zeros((nlons, nlats, nlevs, ntrac),\
                dtype=np.double, order='F')

        self.tempUTend = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.tempVTend = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.tempVrtTend = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')

        self.tempDivTend = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.tempVirtTempTend = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')

        self.tempLnpsTend = np.zeros((nlons, nlats), dtype=np.double, order='F')

        self.__latitudes = np.zeros((nlons, nlats), dtype=np.double, order='F')
        self.__longitudes = np.zeros((nlons, nlats), dtype=np.double, order='F')


        if(ntrac > 0):
            initialiseGridArrays(\
                 <double *>&self.pyUg[0,0,0],\
                 <double *>&self.pyVg[0,0,0],\
                 <double *>&self.pyVrtg[0,0,0],\
                 <double *>&self.pyDivg[0,0,0],\
                 <double *>&self.pyVirtTempg[0,0,0],\
                 <double *>&self.pyTracerg[0,0,0,0],\
                 <double *>&self.pyDlnpdtg[0,0,0],\
                 <double *>&self.pyEtaDotg[0,0,0],\
                 <double *>&self.pyLnPsg[0,0],\
                 <double *>&self.pyPhis[0,0],\
                 <double *>&self.pyDPhisdx[0,0],\
                 <double *>&self.pyDPhisdy[0,0],\
                 <double *>&self.pyDlnpsdt[0,0])
        else:
            initialiseGridArrays(\
                 <double *>&self.pyUg[0,0,0],\
                 <double *>&self.pyVg[0,0,0],\
                 <double *>&self.pyVrtg[0,0,0],\
                 <double *>&self.pyDivg[0,0,0],\
                 <double *>&self.pyVirtTempg[0,0,0],\
                 <double *>0,\
                 <double *>&self.pyDlnpdtg[0,0,0],\
                 <double *>&self.pyEtaDotg[0,0,0],\
                 <double *>&self.pyLnPsg[0,0],\
                 <double *>&self.pyPhis[0,0],\
                 <double *>&self.pyDPhisdx[0,0],\
                 <double *>&self.pyDPhisdy[0,0],\
                 <double *>&self.pyDlnpsdt[0,0])

#Intialise pressure arrays
    def initPressureArrays(self):
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

        self.pySurfPressure = np.zeros((nlons, nlats), dtype=np.double, order='F')
        self.pyPressGrid = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        self.pyInterfacePressure = np.zeros((nlons, nlats, nlevs+1), dtype=np.double, order='F')

        initialisePressureArrays(\
                <double *>&self.pySurfPressure[0,0],\
                <double *>&self.pyPressGrid[0,0,0],\
                <double *>&self.pyInterfacePressure[0,0,0])

# Set tendencies for dynamical core to use in physics
    def setTendencies(self,tendency_list):

        global nlons,nlats,nlevs,ntrac

        cdef cnp.double_t[::1,:,:] tempvrt,tempdiv,tempvt,tempu,tempv
        cdef cnp.double_t[::1,:,:,:] temptracer
        cdef cnp.double_t[::1,:] templnps

        uTend,vTend,virtTempTend,lnpsTend,tracerTend = tendency_list
        #vrtTend,divTend,virtTempTend,lnpsTend,tracerTend = tendency_list
        '''
        if virtTempTend is None:
            tempvt = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        else:
            tempvt = np.asfortranarray(virtTempTend)

        if lnpsTend is None:
            templnps = np.zeros((nlons, nlats), dtype=np.double, order='F')
        else:
            templnps = np.asfortranarray(lnpsTend)

        if tracerTend is None:
            temptracer = np.zeros((nlons, nlats, nlevs, ntrac), dtype=np.double, order='F')
        else:
            temptracer = np.asfortranarray(tracerTend)

        if uTend is None:
            tempu = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        else:
            tempu = np.asfortranarray(uTend)

        if vTend is None:
            tempv = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F')
        else:
            tempv = np.asfortranarray(vTend)
        '''

        tempvt = np.asfortranarray(virtTempTend)

        templnps = np.asfortranarray(lnpsTend)

        temptracer = np.asfortranarray(tracerTend)

        tempu = np.asfortranarray(uTend)

        tempv = np.asfortranarray(vTend)


        #if (uTend.any() and vTend.any()):

        #print
        #print 'adding wind tendencies'
        #print

        tempu = np.asfortranarray(uTend)
        tempv = np.asfortranarray(vTend)
        self.tempUTend[:] = tempu
        self.tempVTend[:] = tempv

        gfs_uv_to_vrtdiv(\
                         <double *>&self.tempUTend[0,0,0],\
                         <double *>&self.tempVTend[0,0,0],\
                         <double *>&self.tempVrtTend[0,0,0],\
                         <double *>&self.tempDivTend[0,0,0],\
                         )

            #tempvrt = np.asfortranarray(self.tempVrtTend) + np.asfortranarray(tempvrt) 
            #tempdiv = np.asfortranarray(self.tempDivTend) + np.asfortranarray(tempdiv)


        #self.tempVrtTend[:] = tempvrt
        #self.tempDivTend[:] = tempdiv
        self.tempVirtTempTend[:] = tempvt
        self.tempLnpsTend[:] = templnps
        self.tempTracerTend[:] = temptracer

        #print 'In set_tend: ', abs(np.asfortranarray(self.tempUTend)[:,:,8]).max(), abs(np.asfortranarray(self.tempVTend)).max(),\
        #    abs(np.asfortranarray(self.tempVirtTempTend)[:,:,8]).max()
        #print 'Vrt, div: ', abs(np.asfortranarray(self.tempVrtTend)[:,:,8]).max(), abs(np.asfortranarray(self.tempDivTend)).max()

        set_tendencies(\
                       <double *>&self.tempVrtTend[0,0,0],\
                       <double *>&self.tempDivTend[0,0,0],\
                       <double *>&self.tempVirtTempTend[0,0,0],\
                       <double *>&self.tempLnpsTend[0,0],\
                       <double *>&self.tempTracerTend[0,0,0,0],\
                       t,dt)

# Take one step
    def oneStepForward(self):

        gfsConvertToGrid()
        gfsCalcPressure()


        if self.physicsEnabled:
            #TODO don't call physics callback directly. use a helper function which will remove
            #TODO individual fields from tracerg and assign them to q, ozone, etc., and then
            #TODO call physics routines
                #self.physicsCallback(self.pyUg,\
                #                     self.pyVg,\
            tendList = \
                self.physicsCallback(self.pyUg,\
                                     self.pyVg,\
                                     self.pyVirtTempg,\
                                     self.pyPressGrid,\
                                     self.pySurfPressure,\
                                     self.pyTracerg,\
                                     self.__latitudes)

            self.setTendencies(tendList)

        else:
            calculate_tendencies(\
                                 <double *>&self.pyVrtg[0,0,0],\
                                 <double *>&self.pyDivg[0,0,0],\
                                 <double *>&self.pyVirtTempg[0,0,0],\
                                 <double *>&self.pyPressGrid[0,0,0],\
                                 <double *>&self.pySurfPressure[0,0],\
                                 <double *>0,\
                                 t,\
                                 dt)
        gfsTakeOneStep()

# Register a callback which calculates the physics (to be used in stand-alone
# mode only)

    cdef setPhysicsCallback(self, physicsFnPtr):

        self.physicsCallback = physicsFnPtr
        self.physicsEnabled = True
#        gfsRegisterPhysicsCallback(testFunc)
    '''
    Does not work!
    def setInitialConditions(self, inputList):
    
    myug,myvg,myvirtempg,mytracerg,mylnpsg = inputList
    
    self.pyUg[:] = myug[:]
    self.pyVg[:] = myvg[:]
    self.pyVirtTempg[:] = myvirtempg[:]
    self.pyTracerg[:] = mytracerg[:]
    self.pyLnPsg[:] = mylnpsg[:]
    
    gfsConvertToSpec()
    '''
    def getResult(self):
    
        gfsConvertToGrid()
        gfsCalcPressure()

        outputList = []

        theta = np.asfortranarray(self.pyVirtTempg)
        q = np.asfortranarray(self.pyTracerg[:,:,:,0])
        q[q<0] = 0

        temp = theta.copy()
        if(np.any(q > 0)):
            temp = theta/(1+self.fv*q)

        outputList.append(np.asarray(self.pyUg).copy(order='F'))
        outputList.append(np.asarray(self.pyVg).copy(order='F'))
        outputList.append(temp)
        outputList.append(q)
        outputList.append(np.asarray(self.pySurfPressure).copy(order='F'))
        outputList.append(np.asarray(self.pyPressGrid).copy(order='F'))
        iface_press = np.zeros(np.asarray(self.pyInterfacePressure).shape, dtype=np.double, order='F')
        iface_press[:] = np.asfortranarray(self.pyInterfacePressure)[:,:,::-1]
        outputList.append(iface_press)

        return(outputList)

    # method to override the parent class (Component) method (to be used in CliMT
    # mode only)

    #!!DO NOT USE!! This function is preserved for historical purposes only
    #!!Use only for testing, to get increments given a set of fields.
    def driver(self, myug, myvg, myvirtempg, myqg, mylnpsg, mypress, double simTime=-1.):


        cdef cnp.double_t[::1,:,:] tempug,tempvg,tempvtg
        cdef cnp.double_t[::1,:,:] tempqg
        cdef cnp.double_t[::1,:] templnpsg
    
    
        if(simTime >= 0):
            global t
            t = simTime
    
        #myug,myvg,myvirtempg,myqg,mylnpsg,mypress = inputArgs
    
        myug = np.asfortranarray(myug)
        myvg = np.asfortranarray(myvg)
        myvirtempg = np.asfortranarray(myvirtempg)
        myqg = np.asfortranarray(myqg)
        mylnpsg = np.asfortranarray(mylnpsg)
    
        #Convert to memory view so that assignment can be made to arrays
        tempug = myug
        tempvg = myvg
        tempvtg = myvirtempg
    
        tempqg = myqg
        templnpsg = mylnpsg
    
        print np.amax(self.pyVirtTempg - myvirtempg)
        #Assign to model arrays
        self.pyUg[:] = tempug
        self.pyVg[:] = tempvg
        self.pyVirtTempg[:] = tempvtg
        self.pyTracerg[:,:,:,0] = tempqg
        self.pyLnPsg[:] = templnpsg
    
    #Convert to spectral space
        gfsConvertToSpec()
    
    #Step forward in time
        self.oneStepForward()
    
    #Convert back to grid space
        gfsConvertToGrid()
    
    # only ln(Ps) is calculated in the dynamics. This calculates
    # the values on the full grid
        gfsCalcPressure()
    
        '''
        ug = np.asfortranarray(self.pyUg.copy())
        vg = np.asfortranarray(self.pyVg.copy())
        virtempg = np.asfortranarray(self.pyVirtTempg.copy())
        tracerg = np.asfortranarray(self.pyTracerg[:,:,:,0].copy())
        lnpsg = np.asfortranarray(mylnpsg.copy())
        press = np.asfortranarray(mypress.copy())
    
        return(ug,vg,virtempg,tracerg,lnpsg,press)
        '''
    
        ugInc = np.ascontiguousarray(self.pyUg - myug)
        vgInc = np.ascontiguousarray(self.pyVg - myvg)
        virtempgInc = np.ascontiguousarray(self.pyVirtTempg - myvirtempg)
        tracergInc = np.ascontiguousarray(self.pyTracerg[:,:,:,0] - myqg)
        lnpsgInc = np.ascontiguousarray(self.pyLnPsg - mylnpsg)
        pressInc = np.ascontiguousarray(self.pyPressGrid - mypress)
    
    
        return(ugInc,vgInc,virtempgInc,tracergInc,lnpsgInc,pressInc)
        #return(ugInc,vgInc,virtempgInc,tracergInc,lnpsgInc,pressInc\
               #,ug,vg,virtempg,tracerg,lnpsg,press)
    
    def printTimes(self):
        global dt,t
        print 'Timestep: ',dt, 'Total Time:', t
    
    def get_nlat(self):
        return self.numLats
    
    def get_nlon(self):
        return self.numLons
    
    def get_nlev(self):
        return self.numLevs

    def initial_conditions(self, ug, vg, virtempg, psg, qg=None):
            #Only to be used in CLIMT mode
        global ntrac,nlons,nlats,nlevs,dt
     
        cdef cnp.double_t[::1,:,:] tempug,tempvg,tempvtg
        cdef cnp.double_t[::1,:,:] tempqg
        cdef cnp.double_t[::1,:] templnpsg

        myug = np.asfortranarray(ug)
        myvg = np.asfortranarray(vg)
        myvirtempg = np.asfortranarray(virtempg)
        #Convert to ln(ps)
        mylnpsg = np.asfortranarray(np.log(psg))
        myqg = np.asfortranarray(np.zeros(ug.shape))
        if qg is not None:
            myqg = np.asfortranarray(qg)


        #Obtain memory view so that assignment can be made to model arrays
        tempug = myug
        tempvg = myvg
        tempvtg = myvirtempg

        templnpsg = mylnpsg
        tempqg = myqg

        #Assign to model arrays
        #print np.amax(np.abs(myug - self.pyUg))
        self.pyUg[:] = tempug[:]
        self.pyVg[:] = tempvg[:]
        self.pyVirtTempg[:] = tempvtg[:]
        self.pyTracerg[:,:,:,0] = tempqg
        self.pyLnPsg[:] = templnpsg[:]
        #self.tempUTend[:] = tempug
        #self.tempVTend[:] = tempvg

        gfsConvertToSpec()

 
    def integrateFields(self,field_list,increment_list):
        #Only to be used in CLIMT mode
        global ntrac,nlons,nlats,nlevs,dt
 
       
        if self.climt_mode:

            q = np.asfortranarray(self.pyTracerg[:,:,:,0])
            q[q<0] = 0
            
            temptrac = np.zeros((nlons,nlats,nlevs,ntrac),dtype=np.double,order='F')
            uTend,vTend, tempTend,qTend,psTend = increment_list
 
            virtTempTend = tempTend
            temperature = self.pyVirtTempg/(1+self.fv*q)

            #Additional tendency term while converting to virtual temp -- see Pg. 12 in doc/gfsModelDoc.pdf
            if(np.any(q > 0)):
                virtTempTend = tempTend*(1+self.fv*q) + self.fv*temperature*qTend

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
            self.setTendencies(increment_list)
 

            #Step forward in time
            self.take_one_step()
 
            #Convert back to grid space
            self.convert_to_grid()
 
            # only ln(Ps) is calculated in the dynamics. This calculates
            # the values on the full grid
            self.calculate_pressure()
 
            ug = np.asfortranarray(self.pyUg.copy())
            vg = np.asfortranarray(self.pyVg.copy())
            virtempg = np.asfortranarray(self.pyVirtTempg.copy())
            qg = np.asfortranarray(self.pyTracerg[:,:,:,0].copy())
            psg = np.asfortranarray(self.pySurfPressure.copy())
            press = np.asfortranarray(self.pyPressGrid.copy())
            iface_press = np.zeros(np.array(self.pyInterfacePressure).shape, dtype=np.double, order='F')
            iface_press[:] = np.asfortranarray(self.pyInterfacePressure.copy())[:,:,::-1]

            qg[qg<0] = 0

            if(np.any(qg > 0)):
                #output temperature, not virtual temp
                virtempg = virtempg/(1+self.fv*qg)

            return(ug,vg,virtempg,qg,psg,press,iface_press)
            
 
    def take_one_step(self):
        gfsTakeOneStep()

    def convert_to_grid(self):
        gfsConvertToGrid()

    def calculate_pressure(self):
        gfsCalcPressure()

    def set_topography(self, topography):
        cdef cnp.double_t[::1,:] temp_topo

        topography = np.asfortranarray(topography)

        temp_topo = topography

        self.pyPhis[:] = temp_topo[:]

        gfsSetTopography(<double *>&self.pyPhis[0,0])

    def shutDownModel(self):

        global t
    
        if self.modelIsInitialised:
            t = 0
            gfsFinalise()
            self.modelIsInitialised = False

    def __dealloc__(self):

        self.shutDownModel()

