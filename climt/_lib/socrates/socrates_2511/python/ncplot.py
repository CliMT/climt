# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots a mean profile of the variable in the supplied file
against height (calculated from pressure assuming isothermal atmos.)
If channel dimension exists, spectrum at boundaries is also plotted.
'''

CONST = 287. * 250. / (9.80665 * 1000.0) #

from netCDF4 import Dataset
import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    if (len(sys.argv) > 4):
        filename4 = sys.argv[4]
        name4 = filename4[filename4.find('.') + 1:]
        dgs4 = Dataset(filename4)
        var4 = dgs4.variables[name4][:]
    if (len(sys.argv) > 3):
        filename3 = sys.argv[3]
        name3 = filename3[filename3.find('.') + 1:]
        dgs3 = Dataset(filename3)
        var3 = dgs3.variables[name3][:]
    if (len(sys.argv) > 2):
        filename2 = sys.argv[2]
        name2 = filename2[filename2.find('.') + 1:]
        dgs2 = Dataset(filename2)
        var2 = dgs2.variables[name2][:]
    if (len(sys.argv) > 1):
        filename = sys.argv[1]
        name = filename[filename.find('.') + 1:]
        dgs = Dataset(filename)
        var = dgs.variables[name][:]
    else:
        raise RuntimeError('please enter a file name')

lon = dgs.variables['lon'][:]
lat = dgs.variables['lat'][:]
p = dgs.variables['plev'][:]

n_lon = len(lon)
n_lat = len(lat)
layers= len(p)
e_factor=2.718

if ('ph_rate' in name):
    e_factor=1.0e2
    xlim = 0.25

if (len(sys.argv) > 2):
    var = var - var2

try:
    width     = dgs.variables['bandwidth'][:]
    wl_short  = dgs.variables['wl_short'][:]
    wl_long   = dgs.variables['wl_long'][:]
    n_channel = len(width)
except:
    vmean = np.zeros(layers)
    for i in np.arange(layers):
        vmean[i] = np.sum(var[i, :, :]) / (n_lon * n_lat)
    ax1 = plt.figure().add_subplot(111)
    ax1.plot(vmean, -np.log(p/max(p))*CONST, label=name)
    ax1.set_title(dgs.variables[name].title)
    if ('ph_rate' in name):
        ax1.set_xlabel('J rate (s$^{-1}$)')
    elif (name == 'aflx'):
        ax1.set_xlabel('Actinic Flux (W m$^{-2}$)')
    elif (name == 'hrts'):
        ax1.set_xlabel('Heating rate (K day$^{-1}$)')
    else:
        ax1.set_xlabel('Flux (W m$^{-2}$)')
    ax1.set_ylabel('Approx height (km)')
else:
    fig=plt.figure()
    vmean = np.zeros(layers)
    vtop = np.sum(var[:, 0, :, :]) / (n_lon * n_lat)
    for i in np.arange(layers):
        vmean[i] = np.sum(var[:, i, :, :]) / (n_lon * n_lat)
        if (vmean[i] > vtop/e_factor):
          e_layer = i
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(vmean, -np.log(p/max(p))*CONST)
    ax1.set_title(dgs.variables[name].title)
    ax1.set_xlabel(name)
    ax1.set_ylabel('Approx height (km)')
    wn = 0.5e-2/wl_short + 0.5e-2/wl_long
    wl = 0.5e6*(wl_short + wl_long)
    toa_spec = np.zeros(n_channel)
    surf_spec = np.zeros(n_channel)
    for ch in range(0, n_channel):
        toa_spec[ch]  = np.sum(var[ch,0, :,:])/(width[ch]*n_lon*n_lat)
        surf_spec[ch] = np.sum(var[ch,layers-1,:,:])/(width[ch]*n_lon*n_lat)
    ax2.plot(wl, toa_spec, color='blue', label='TOA')
    ax2.set_xlabel('Wavelength (micron)')
    if ('ph_rate' in name):
        ax2.set_ylabel('J rate (s$^{-1}$ m$^{-1}$)')
        ax2.set_yscale('symlog')
        e_height = -np.log(p[e_layer]/max(p))*CONST
        ax1.plot([min(vmean),max(vmean)], [e_height,e_height], color='green')
        ax1.set_xscale('log')
        mid_spec = np.zeros(n_channel)
        for ch in range(0, n_channel):
            mid_spec[ch]  = np.sum(var[ch,e_layer, :,:])/(width[ch]*n_lon*n_lat)
        ax2.plot(wl, mid_spec, color='green', label='Mid atmos')
        ax2.set_title('Top & Mid atmosphere spectra')
        ax2.set_xlim(right=xlim)
    elif (name == 'aflx'):
        ax2.set_ylabel('Actinic Flux (W m$^{-2}$ m$^{-1}$)')
        ax2.set_yscale('symlog')
        e_layer = int(layers/2)
        e_height = -np.log(p[e_layer]/max(p))*CONST
        ax1.plot([min(vmean),max(vmean)], [e_height,e_height], color='green')
        mid_spec = np.zeros(n_channel)
        for ch in range(0, n_channel):
            mid_spec[ch]  = np.sum(var[ch,e_layer, :,:])/(width[ch]*n_lon*n_lat)
        ax2.plot(wl, mid_spec, color='green', label='Mid atmos')
        ax2.set_title('Top & Mid atmosphere spectra')
    elif (name == 'hrts'):
        ax2.set_ylabel('Heating rate (K day$^{-1}$ m$^{-1}$)')
        e_layer = int(layers/2)
        e_height = -np.log(p[e_layer]/max(p))*CONST
        ax1.plot([min(vmean),max(vmean)], [e_height,e_height], color='green')
        mid_spec = np.zeros(n_channel)
        for ch in range(0, n_channel):
            mid_spec[ch]  = np.sum(var[ch,e_layer, :,:])/(width[ch]*n_lon*n_lat)
        ax2.plot(wl, mid_spec, color='green', label='Mid atmos')
        ax2.set_title('Top & Mid atmosphere spectra')
    else:
        ax2.plot(wl, surf_spec, color='green', label='Surface')
        ax2.set_title('TOA & surface spectrum')
        ax2.set_ylabel('Flux (Wm-2m-1)')
        ax2.set_xscale('symlog')
    plt.legend()

plt.tight_layout()
plt.show()
