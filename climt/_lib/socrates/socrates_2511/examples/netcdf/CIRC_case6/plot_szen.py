# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots values against solar zenith angle
'''

from netCDF4 import Dataset
import sys
import numpy as np
import matplotlib.pyplot as plt

hrl=8
szen_list_sg=np.arange(75,99)
nszen_sg=len(szen_list_sg)
surfsw_sg=np.zeros(nszen_sg)
szen_sg=np.zeros(nszen_sg)
strat_hr_sg=np.zeros(nszen_sg)
for i in np.arange(nszen_sg):
    szen_sg[i]=szen_list_sg[i]
    filename='case6_sw_sg_ga7_'+str(szen_list_sg[i])+'.vflx'
    dgs = Dataset(filename)
    p = dgs.variables['plev'][:]
    layers= len(p)
    surfsw_sg[i]= dgs.variables['vflx'][layers-1]
    filename='case6_sw_sg_ga7_'+str(szen_list_sg[i])+'.hrts'
    dgs = Dataset(filename)
    p = dgs.variables['plev'][:]
    layers= len(p)
    strat_hr_sg[i]= dgs.variables['hrts'][hrl]
    pressure=str(p[hrl])

szen_list=np.arange(75,90)
nszen=len(szen_list)
surfsw=np.zeros(nszen+1)
szen=np.zeros(nszen+1)
strat_hr=np.zeros(nszen+1)
surfsw[nszen]=0.0
szen[nszen]=90.0
strat_hr[nszen]=0.0
for i in np.arange(nszen):
    szen[i]=szen_list_sg[i]
    filename='case6_sw_ga7_'+str(szen_list_sg[i])+'.vflx'
    dgs = Dataset(filename)
    p = dgs.variables['plev'][:]
    layers= len(p)
    surfsw[i]= dgs.variables['vflx'][layers-1]
    filename='case6_sw_ga7_'+str(szen_list_sg[i])+'.hrts'
    dgs = Dataset(filename)
    strat_hr[i]= dgs.variables['hrts'][hrl] 
                                            
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.plot(szen_sg, surfsw_sg)
ax1.plot(szen, surfsw)
ax1.set_xlabel('Solar zenith angle (degrees)')
ax1.set_ylabel('Surface SW Flux (W m-2)')
ax2 = fig.add_subplot(122)
ax2.plot(szen_sg, strat_hr_sg)
ax2.plot(szen, strat_hr)
ax2.set_xlabel('Solar zenith angle (degrees)')
ax2.set_ylabel('Heating rates (K day-1) at P='+pressure)

plt.tight_layout()
plt.show()
