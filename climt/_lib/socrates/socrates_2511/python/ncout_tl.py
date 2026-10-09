# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Program to create netCDF file of temperature on levels (.tl file).
This is extrapolated from the temperature in layers (.t file),
and the surface temperature (.tstar) if available.
'''

import numpy as np
import nctools as nc
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os.path, sys

def ncout_tl(basename):
    CONST = 287. * 250. / (9.80665 * 1000.0)
    
    t_data = Dataset(basename + '.t')
    
    lon = t_data.variables['lon'][:]
    lat = t_data.variables['lat'][:]
    p = t_data.variables['plev'][:]
    t = t_data.variables['t'][:]
    
    n_lon = len(lon)
    n_lat = len(lat)
    layers = len(p)
    
    tlev = np.zeros((layers + 1)*n_lat*n_lon).reshape(layers + 1,n_lat,n_lon)
    
    order = np.argsort(p)[::-1]
    p = p[order]
    t = t[order, :, :]
    
    if (os.path.exists(basename + '.tstar')):
        tstar_data = Dataset(basename + '.tstar')
        tstar = tstar_data.variables['tstar'][:]
        pstar = tstar_data.variables['plev'][:]
    else :
        pstar = np.array([p[0]])
      
    pl = np.zeros(layers + 1)
    pl[0] = pstar
    pl[1:layers] = 0.5 * (p[0:layers-1] + p[1:layers])
    pl[layers] = p[layers-1] * pl[layers-1] / p[layers-2]
    
    log_p = np.log(p)
    log_pl = np.log(pl)
    
    
    for i in range(n_lon):
        for j in range(n_lat):
# for points outside the domain, np.interp uses the y value for the minimum 
# x-value in the domain, so these values need to be calculated manually
            left = None
            right = None
            if (np.max(log_pl) > log_p[0]):
                gradient = (t[1] - t[0]) / (log_p[1] - log_p[0])
                intercept = t[0] - gradient*log_p[0]
                right = gradient * np.max(log_pl) + intercept
            if (np.min(log_pl) < log_p[-1]):
                gradient = (t[-1] - t[-2]) / (log_p[-1] - log_p[-2])
                intercept = t[-1] - gradient*log_p[-1]
                left = gradient * np.min(log_pl) + intercept
            tlev[:,j,i] = np.interp(log_pl, log_p[::-1], t[::-1,j,i], 
                                    left = left, right = right)
    
    if (os.path.exists(basename + '.tstar')):
# By default use the bottom layer temperature for the temperature
# just above the surface:
        tlev[0,:,:] = t[0,:,:]
# Uncomment the following to use the surface temperature instead
#        tlev[0,:,:] = tstar
    else:
        nc.ncout3d(basename + '.tstar', lon, lat, pstar, tlev[0,:,:],
                   longname = 'Surface Temperature', units = 'K')
    
    plt.figure()
    plt.xlabel('Temperature (K)')
    plt.ylabel('Approx height (km)')
    plt.title('Interpolation of temperature to levels')
    plt.plot(t[:,0,0], -np.log(p/pl[0])*CONST)
    plt.plot(tlev[:,0,0], -np.log(pl/pl[0])*CONST, 'r+')
    plt.show()
    
    nc.ncout3d(basename + '.tl', lon, lat, pl, tlev, 
               longname = 'Temperature on levels', units = 'K')
    
    print('ncout_tl - file: ',basename + '.tl')
    
if __name__ == '__main__':
    if (len(sys.argv) > 1):
        basename = sys.argv[1]
    else:
        raise RuntimeError('please enter a file base name')
    ncout_tl(basename)
