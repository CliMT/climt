# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Converts the text files from http://circ.gsfc.nasa.gov
into netCDF input files for the ES code. Takes the case number
as an argument and should be run from within the dirctory containing
the circ files.
Make sure you have sourced set_rad_env first.
'''

import nctools as nc
import numpy as np
import sys, math, subprocess

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        arg1 = int(sys.argv[1])
    else:
        raise RuntimeError('please enter a CIRC case number '
                           '(e.g. python circ2nc.py 2)')

rad_data = subprocess.run(['echo $RAD_DATA'], stdout=subprocess.PIPE, universal_newlines=True, shell=True).stdout[:-1]
SPECTRAL_FILE = [rad_data + '/spectra/ga9/sp_sw_ga9'
                 , rad_data + '/spectra/ga9_ref/sp_sw_260_jm3']
SP_FILE_NAME = ['ga9', '260']

nspfs = len(SPECTRAL_FILE)
ncase = str(arg1)
Tsfc_sza_nlev = 'Tsfc_sza_nlev_case' + ncase + '.txt'
aerosol_input = 'aerosol_input_case' + ncase + '.txt'
cloud_input = 'cloud_input_case'  +  ncase  +  '.txt'
layer_input = 'layer_input_case' + ncase + '.txt'
level_input = 'level_input_case' + ncase + '.txt'
sfcalbedo_input = 'sfcalbedo_input_case' + ncase + '.txt'
SW_charts = 'SW_charts_1cm-1_case' + ncase + '.txt'
LW_lblrtm_bb = 'LW_lblrtm_bb_case' + ncase + '.txt'

with open(Tsfc_sza_nlev, 'r') as file:
    contents = file.readlines()
nlevs = int(contents[5])
tstar = float(contents[6])
szen = float(contents[7])
stoa = float(contents[8])

print('Number of levels: ', nlevs)
print('Surface temperature: ', tstar)
print('Solar zenith angle: ', szen)
print('Solar irradiance: ', stoa)
stoa=stoa/math.cos(szen*math.pi/180.)
print('Solar irradiance (corrected): ', stoa)

with open(aerosol_input, 'r') as file:
    contents = file.readlines()
alpha = float(contents[3])
str = (''.join(contents[6:])).split()
aerosol = [float(i) for i in str]
aerosol = np.array(aerosol).reshape(nlevs - 1, 4)

with open(cloud_input, 'r') as file:
    contents = file.readlines()
str = (''.join(contents[3:])).split()
cloud = [float(i) for i in str]
cloud = np.array(cloud).reshape(nlevs - 1, 6)

with open(layer_input, 'r') as file:
    contents = file.readlines()
legend = contents[2].split()
str = (''.join(contents[3:])).split()
layer = [float(i) for i in str]
layer = np.array(layer).reshape(nlevs - 1, 13)

with open(level_input, 'r') as file:
    contents = file.readlines()
str = (''.join(contents[1:])).split()
level = [float(i) for i in str]
level = np.array(level).reshape(nlevs, 4)

with open(sfcalbedo_input, 'r') as file:
    contents = file.readlines()
str = (''.join(contents[6:])).split()
albedo = [float(i) for i in str]
albedo = np.array(albedo).reshape(49180, 4)
wvl = albedo[:, 0]

with open(SW_charts, 'r') as file:
    contents = file.readlines()
str = (''.join(contents[5:])).split()
charts = [float(i) for i in str]
charts = np.array(charts).reshape(49180, 5)

with open(LW_lblrtm_bb, 'r') as file:
    contents = file.readlines()
str = (''.join(contents[4:])).split()
lwout = [float(i) for i in str]
lwout = np.array(lwout).reshape(nlevs, 5)

# -------------------------------------

basename = 'case' + ncase
lon = np.zeros(1)
lat = np.zeros(1)
pl = np.zeros(nlevs)
tl = np.zeros(nlevs)
z = np.zeros(nlevs) # height on levels in m

lwflxup = np.zeros(nlevs)
lwflxdn = np.zeros(nlevs)

lwalb   = 0.0
basis   = 1

pstar = np.zeros(1) + level[0, 2]
for i in np.arange(nlevs):
    z[i] = level[i, 1] * 1000.0
    pl[i] = level[i, 2] * 100.0
    tl[i] = level[i, 3]
    lwflxup[i] = lwout[nlevs - 1 - i, 2]
    lwflxdn[i] = lwout[nlevs - 1 - i, 3]

lwHR = np.zeros(nlevs - 1)

p = np.zeros(nlevs - 1)
t = np.zeros(nlevs - 1)
h2o = np.zeros(nlevs - 1)
co2 = np.zeros(nlevs - 1)
o3 = np.zeros(nlevs - 1)
n2o = np.zeros(nlevs - 1)
co = np.zeros(nlevs - 1)
ch4 = np.zeros(nlevs - 1)
o2 = np.zeros(nlevs - 1)
ccl4 = np.zeros(nlevs - 1)
cfc11 = np.zeros(nlevs - 1)
cfc12 = np.zeros(nlevs - 1)

clfr = np.zeros(nlevs - 1)
lwc = np.zeros(nlevs - 1)
iwc = np.zeros(nlevs - 1)
re = np.zeros(nlevs - 1)
ire = np.zeros(nlevs - 1)

for i in np.arange(nlevs - 1):
    lwHR[i] = lwout[nlevs - 1 - i, 4]
    p[i] = layer[i, 1] * 100.0
    t[i] = layer[i, 2]
    h2o[i] = layer[i, 3] * 18.0154 / 28.97
    co2[i] = layer[i, 4] * 44.0098 / 28.97
    o3[i] = layer[i, 5] * 47.9982 / 28.97
    n2o[i] = layer[i, 6] * 44.0128 / 28.97
    co[i] = layer[i, 7] * 28.0101 / 28.97
    ch4[i] = layer[i, 8] * 16.0428 / 28.97
    o2[i] = layer[i, 9] * 31.9988 / 28.97
    ccl4[i] = layer[i, 10] * 153.8215 / 28.97
    cfc11[i] = layer[i, 11] * 137.370 / 28.97
    cfc12[i] = layer[i, 12] * 120.910 / 28.97

    clfr[i] = cloud[i, 1]
    lwc[i] = cloud[i, 2] / ((z[i + 1] - z[i] ) * 1000.0)
    iwc[i] = cloud[i, 3] / ((z[i + 1] - z[i] ) * 1000.0)
  
    re[i] = cloud[i, 4] * 1.0e-6
    ire[i] = cloud[i, 5] * 1.0e-6

re[np.where(re == 0.0)] = 6.0e-6
ire[np.where(ire == 0.0)] = 30.0e-6

# Temperature:
nc.ncout3d(basename + '.tstar', lon, lat, pstar, tstar, longname= 'Surface temperature', units= 'K')
nc.ncout3d(basename + '.t', lon, lat, p, t, longname= 'Temperature', units= 'K')
nc.ncout3d(basename + '.tl', lon, lat, pl, tl, longname= 'Temperature on levels', units= 'K')

# Solar beam:
nc.ncout2d(basename + '.stoa', lon, lat, stoa, longname= 'Solar Irradiance', units= 'WM-2')
nc.ncout2d(basename + '.szen', lon, lat, szen, longname= 'Solar Zenith Angle', units= 'degrees')

# Surface albedos:
# Note: <basename>.surf should be linked to <basename>.surfsw or
# <basename>.surflw  for the SW and LW calls to the radiation code
# respectively.
nc.ncout_surf(basename + '.surflw', lon, lat, basis, lwalb)

# Aerosols:
ammr = aerosol[0:nlevs - 1, 1] / (z[1:nlevs] - z[0:nlevs - 1])

sel = np.where(aerosol[:, 1] > 0.0)
if (sel[0].size > 0):
    omega= aerosol[sel[0][0], 2]
    asym= aerosol[sel[0][0], 3]
else:
    omega= 0.0
    asym= 0.0

for i in range(nspfs):
    bands= 0
    with open(SPECTRAL_FILE[i], 'r') as file:
        contents = file.readlines()
        for (index, line) in enumerate(contents[0:99]):
            if (line[1:24] == 'umber of spectral bands' or 
                line[1:24] == 'UMBER OF SPECTRAL BANDS'):
                bands = int(line[26:])
            if (line[1:35] == 'pecification of spectral intervals' or 
                line[1:35] == 'PECIFICATION OF SPECTRAL INTERVALS'):
                    str = (''.join(contents[index + 3:index + bands + 3])).split()
                    band_limits = [float(j) for j in str]
                    band_limits = np.array(band_limits).reshape(bands, 3)

    kabs =  np.zeros(bands)
    kscat =  np.zeros(bands)
    absp =  np.zeros((nlevs - 1)*bands).reshape(bands, nlevs - 1)
    scat =  np.zeros((nlevs - 1)* bands).reshape(bands, nlevs - 1)

    swalb =  np.zeros(bands)
    wswalb =  np.zeros(bands)
    w2swalb =  np.zeros(bands)
    sp_file = SP_FILE_NAME[i] 
    for band in range(bands):
        ll= max(1.0 / (100.0*band_limits[band, 1]), 820.5)
        ul= min(1.0 / (100.0*band_limits[band, 2]), 49998.5)
        sel= np.where((wvl >= ul) & (wvl <= ll))
        swalb[band]= np.sum(albedo[sel, 1]*albedo[sel, 3]) / np.sum(albedo[sel, 3])
        wswalb[band]= np.sum(charts[sel, 4]) / max(np.sum(charts[sel, 1]),1.0e-10)
        w2swalb[band]= np.sum(albedo[sel, 1]*charts[sel, 1]) / max(np.sum(charts[sel,1]),1.0e-10)

        kext= np.sum((1.0e4 / wvl[sel])**(-alpha)) / sel[0].size
        if (i == 0):
            print(band + 1, ' Tau: ', kext*aerosol[0, 1], kext*aerosol[1, 1], end=' ')  
            print(kext*aerosol[2, 1], kext*aerosol[3, 1], kext*aerosol[4, 1], end=' ') 
            print(kext*aerosol[5, 1])
        kscat[band]= omega*kext
        kabs[band]= kext - kscat[band]
        absp[band, :]= ammr*kabs[band]
        scat[band, :]= ammr*kscat[band]

    nc.ncout_spectral_surf(basename + '.surfsw_' + sp_file, lon, lat, bands, swalb)  
    nc.ncout_spectral_surf(basename + '.surfwsw_' + sp_file, lon, lat, bands, wswalb  )
    nc.ncout_spectral_surf(basename + '.surfw2sw_' + sp_file, lon, lat, bands, w2swalb)
    print('bands=', bands)
    print('absp=', absp)
    nc.ncout_opt_prop(basename + '.op_soot_' + sp_file, lon, lat, p, bands, absp, scat, asym)

# Spectrally constant albedos:
if (ncase == '1'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.196)
if (ncase == '2'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.188)
if (ncase == '3'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.171)
if (ncase == '4'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.670)
if (ncase == '5'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.670)
if (ncase == '6'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.136)
if (ncase == '7'): nc.ncout_surf(basename + '.surfsw_con', lon, lat, 1, 0.164)

# Main Gases:
nc.ncout3d(basename + '.q',   lon, lat, p, h2o, longname = 'Specific Humidity')
nc.ncout3d(basename + '.co2', lon, lat, p, co2, longname = 'Carbon Dioxide MMR')
nc.ncout3d(basename + '.ch4', lon, lat, p, ch4, longname = 'Methane MMR')
nc.ncout3d(basename + '.n2o', lon, lat, p, n2o, longname = 'Dinitrogen oxide MMR')
nc.ncout3d(basename + '.o2',  lon, lat, p, o2,  longname = 'Oxygen MMR')
nc.ncout3d(basename + '.o3',  lon, lat, p, o3,  longname = 'Ozone MMR')

# CFCs:
nc.ncout3d(basename + '.cfc11',   lon, lat, p, cfc11,   longname = 'CFC-11 MMR')
nc.ncout3d(basename + '.cfc12',   lon, lat, p, cfc12,   longname = 'CFC-12 MMR')

# Gases currently unused:
#nc.ncout3d(basename + '.co',  lon, lat, p, h2o, longname = 'Carbon Monoxide MMR')
#nc.ncout3d(basename + '.ccl4',lon, lat, p, ccl4,longname = 'CCL4 MMR')

# Cloud fields:
nc.ncout3d(basename + '.clfr',lon, lat, p, clfr,longname = 'Cloud fraction')
nc.ncout3d(basename + '.lwc', lon, lat, p, lwc, longname = 'Liquid water content', units = 'KGM-3')
nc.ncout3d(basename + '.iwc', lon, lat, p, iwc, longname = 'Ice water content', units = 'KGM-3')
nc.ncout3d(basename + '.re',  lon, lat, p, re,  longname = 'Droplet effective radius', units = 'M')
nc.ncout3d(basename + '.ire', lon, lat, p, ire, longname = 'Ice effective radius', units = 'M')

subprocess.call(['Cdentomix','-q',''+ basename + '.q','-t','' + basename + '.t','-o','' + basename + '.lwm','' + basename + '.lwc'])
subprocess.call(['Cdentomix','-q','' + basename + '.q','-t','' + basename + '.t','-o','' + basename + '.iwm','' + basename + '.iwc'])

# Convective cloud fields:
#nc.ncout3d(basename + '.ccfr',lon, lat, p, 0.0, longname = 'Convective Cloud fraction')
#nc.ncout3d(basename +'.lwmcv',lon, lat, p, 0.0, longname = 'Convective liquid water content', units = 'KGM-3')
#nc.ncout3d(basename +'.iwmcv',lon, lat, p, 0.0, longname = 'Convective ice water content', units = 'KGM-3')
#nc.ncout3d(basename + '.recv',lon, lat, p, re,  longname = 'Droplet effective radius', units = 'M')
#nc.ncout3d(basename +'.irecv',lon, lat, p, ire, longname = 'Ice effective radius', units = 'M')

# LW output files:
nc.ncout3d(basename + '_lwref.dflx', lon, lat, pl, lwflxdn, longname = 'downward flux', units = 'WM-2')
nc.ncout3d(basename + '_lwref.uflx', lon, lat, pl, lwflxup, longname = 'upward flux', units = 'WM-2')
nc.ncout3d(basename + '_lwref.hrts', lon, lat, p, lwHR, longname = 'heating rates', units = 'K.D-1')
