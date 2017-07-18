cimport numpy as cnp
import numpy as np

#DCMIP function in the fortran code which returns
#the initial conditions for the baroclinic wave test
cdef extern:
    void dcmipBaroclinicWave(int *moist,double *X,\
                             double *lon, double *lat,\
                             double *p, double *z,\
                             int *zcoords,\
                             double * u, double * v,\
                             double * w, double * t,\
                             double * phis, double *ps,\
                             double * rho, double *q,\
                             double * q1, double *q2, int* perturb)

#the initial conditions for the tropical cyclone test
cdef extern:
    void dcmipTropicalCyclone(double *lon, double *lat,\
                             double *p, double *z,\
                             int *zcoords,\
                             double * u, double * v,\
                             double * w, double * t,\
                             double * phis, double *ps,\
                             double * rho, double *q)


#What we are sending to the code
cdef double latitude, longitude, pressure
cdef int moist
cdef int add_perturbation

#What we want from the code
cdef double u, v, t, phis, ps, q

#What we dont care about
cdef double dummy, X
cdef int dummy_int, zcoords

def get_baroclinic_wave_ics(press, lon, lat, perturb=False, moist_sim=False):
    '''
    Provides the initial conditions to run the DCMIP baroclinic
    wave test

    '''

    global longitude, latitude, pressure, moist, zcoords, add_perturbation, moist

    num_lons, num_lats, num_levs = press.shape

    zonal_vel = np.zeros(press.shape, dtype=np.double, order='F')
    merid_vel = np.zeros(press.shape, dtype=np.double, order='F')
    temperature = np.zeros(press.shape, dtype=np.double, order='F')
    specific_humidity = np.zeros(press.shape, dtype=np.double, order='F')
    surf_geop = np.zeros((num_lons,num_lats), dtype=np.double, order='F')
    surf_pressure = np.zeros((num_lons,num_lats), dtype=np.double, order='F')

    # This makes the code evaluate eta using the pressure.
    zcoords = 0

    if moist_sim:
        moist = 1
    else:
        moist = 0

    add_perturbation = 0

    if perturb is True:
        add_perturbation = 1

    #Scale factor, not used in code, but present in argument list!
    X = 1.

    # lat-lon in RADIANS!!

    for i in range(num_lons):
        for j in range(num_lats):
            for k in range(num_levs):
                longitude = lon[i,j]
                latitude = lat[i,j]
                pressure = press[i,j,k]
                dcmipBaroclinicWave(&moist, &X,\
                                    &longitude, &latitude,\
                                    &pressure, &dummy,\
                                    &zcoords, \
                                    &u, &v,\
                                    &dummy, &t,\
                                    &phis, &ps,\
                                    &dummy, &q,\
                                    &dummy, &dummy, &add_perturbation)
                zonal_vel[i,j,k] = u
                merid_vel[i,j,k] = v
                temperature[i,j,k] = t
                specific_humidity[i,j,k] = q
                surf_geop[i,j] = phis
                surf_pressure[i,j] = ps

    return zonal_vel, merid_vel, temperature, specific_humidity, surf_pressure, surf_geop

def get_tropical_cyclone_ics(press, lon, lat, **kwargs):
    '''
    Provides the initial conditions to run the DCMIP baroclinic
    wave test

    '''

    global longitude, latitude, pressure, zcoords

    num_lons, num_lats, num_levs = press.shape

    zonal_vel = np.zeros(press.shape, dtype=np.double, order='F')
    merid_vel = np.zeros(press.shape, dtype=np.double, order='F')
    temperature = np.zeros(press.shape, dtype=np.double, order='F')
    vapour = np.zeros(press.shape, dtype=np.double, order='F')
    surf_geop = np.zeros((num_lons,num_lats), dtype=np.double, order='F')
    surf_pressure = np.zeros((num_lons,num_lats), dtype=np.double, order='F')

    # This makes the code evaluate eta using the pressure.
    zcoords = 0

    # lat-lon in RADIANS!!

    for i in range(num_lons):
        for j in range(num_lats):
            for k in range(num_levs):
                longitude = lon[i,j]
                latitude = lat[i,j]
                pressure = press[i,j,k]
                dcmipTropicalCyclone(&longitude, &latitude,\
                                    &pressure, &dummy,\
                                    &zcoords, \
                                    &u, &v,\
                                    &dummy, &t,\
                                    &phis, &ps,\
                                    &dummy, &q)
                zonal_vel[i,j,k] = u
                merid_vel[i,j,k] = v
                temperature[i,j,k] = t
                vapour[i,j,k] = q
                surf_geop[i,j] = phis
                surf_pressure[i,j] = ps

    return zonal_vel, merid_vel, temperature, vapour, surf_pressure, surf_geop
