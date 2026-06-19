# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Plots output of RFMIP calculations
'''

from netCDF4 import Dataset
import nctools as nc
import numpy as np
import matplotlib.pyplot as plt
import sys, math, subprocess, time

def plot_rfmip_out():
    physics_index_expt=1
    physics_index_cont=2
    rfdat = Dataset('multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc')
    lw_u_data = Dataset('rlu_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_expt))
    lw_d_data = Dataset('rld_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_expt))
    sw_u_data = Dataset('rsu_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_expt))
    sw_d_data = Dataset('rsd_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_expt))
    lw_u_data_2 = Dataset('rlu_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_cont))
    lw_d_data_2 = Dataset('rld_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_cont))
    sw_u_data_2 = Dataset('rsu_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_cont))
    sw_d_data_2 = Dataset('rsd_Efx_HadGEM3-GC31-LL_rad-irf_r1i1p{0}f3_gn.nc'.format(physics_index_cont))
    nlevs = len(rfdat.dimensions['level'])
    nprof = len(rfdat.dimensions['site'])
    nexpt = len(rfdat.dimensions['expt'])

    expt_label    = rfdat.variables['expt_label'][:] # expt
    weight        = rfdat.variables['profile_weight'][:] # profile

    pl            = lw_u_data.variables['plev'][:,:] # profile, levels
    lw_down       = lw_d_data.variables['rld'][:,:,:] # expt, profile, levels
    lw_up         = lw_u_data.variables['rlu'][:,:,:] # expt, profile, levels
    sw_down       = sw_d_data.variables['rsd'][:,:,:] # expt, profile, levels
    sw_up         = sw_u_data.variables['rsu'][:,:,:] # expt, profile, levels
    
    lw_down2      = lw_d_data_2.variables['rld'][:,:,:] # expt, profile, levels
    lw_up2        = lw_u_data_2.variables['rlu'][:,:,:] # expt, profile, levels
    sw_down2      = sw_d_data_2.variables['rsd'][:,:,:] # expt, profile, levels
    sw_up2        = sw_u_data_2.variables['rsu'][:,:,:] # expt, profile, levels


    label=expt_label.copy()
    label[0]='PD'
    label[1]='PI'
    label[3]='future'
    label[8]='PI_CO2'
    label[9]='PI_CH4'
    label[10]='PI_N2O'
    label[11]='PI_O3'
    label[12]='PI_HCs'
    label[13]='p4K'
    label[14]='p4KconRH'
    label[15]='PI_all'
    label[16]='future_all'
    label[17]='LGM'

    height = np.zeros((nprof, nlevs))
    for i in np.arange(nprof):
        height[i,:] = -np.log(pl[i,:]/max(pl[i,:]))*287.0*250.0/(9.80665*1000.0)

    lw_down_mean = np.zeros((nexpt, nlevs))
    lw_up_mean   = np.zeros((nexpt, nlevs))
    sw_down_mean = np.zeros((nexpt, nlevs))
    sw_up_mean   = np.zeros((nexpt, nlevs))
    for j in np.arange(nlevs):
        lw_down_mean[0,j] = np.sum(lw_down[0,:,j]*weight[:])
        lw_up_mean[0,j]   = np.sum(lw_up[0,:,j]*weight[:])
        sw_down_mean[0,j] = np.sum(sw_down[0,:,j]*weight[:])
        sw_up_mean[0,j]   = np.sum(sw_up[0,:,j]*weight[:])
    for i in np.arange(nexpt-1)+1:
        for j in np.arange(nlevs):
            lw_down_mean[i,j] = np.sum(lw_down[i,:,j]*weight[:]) - np.sum(lw_down[0,:,j]*weight[:])
            lw_up_mean[i,j]   = np.sum(lw_up[i,:,j]*weight[:])   - np.sum(lw_up[0,:,j]*weight[:])
            sw_down_mean[i,j] = np.sum(sw_down[i,:,j]*weight[:]) - np.sum(sw_down[0,:,j]*weight[:])
            sw_up_mean[i,j]   = np.sum(sw_up[i,:,j]*weight[:])   - np.sum(sw_up[0,:,j]*weight[:])

    lw_down2_mean = np.zeros((nexpt, nlevs))
    lw_up2_mean   = np.zeros((nexpt, nlevs))
    sw_down2_mean = np.zeros((nexpt, nlevs))
    sw_up2_mean   = np.zeros((nexpt, nlevs))
    for j in np.arange(nlevs):
        lw_down2_mean[0,j] = np.sum(lw_down2[0,:,j]*weight[:])
        lw_up2_mean[0,j]   = np.sum(lw_up2[0,:,j]*weight[:])
        sw_down2_mean[0,j] = np.sum(sw_down2[0,:,j]*weight[:])
        sw_up2_mean[0,j]   = np.sum(sw_up2[0,:,j]*weight[:])
    for i in np.arange(nexpt-1)+1:
        for j in np.arange(nlevs):
            lw_down2_mean[i,j] = np.sum(lw_down2[i,:,j]*weight[:]) - np.sum(lw_down2[0,:,j]*weight[:])
            lw_up2_mean[i,j]   = np.sum(lw_up2[i,:,j]*weight[:])   - np.sum(lw_up2[0,:,j]*weight[:])
            sw_down2_mean[i,j] = np.sum(sw_down2[i,:,j]*weight[:]) - np.sum(sw_down2[0,:,j]*weight[:])
            sw_up2_mean[i,j]   = np.sum(sw_up2[i,:,j]*weight[:])   - np.sum(sw_up2[0,:,j]*weight[:])

    fig=plt.figure()
    params = plt.subplots_adjust(left=.03, right=.99, top=.98, bottom=.02)
    for expt in np.arange(18):
        ax = fig.add_subplot(5,4,expt+1)
#    for expt in [9]:
#        ax = fig.add_subplot(1,1,1)
        ax.plot(sw_down_mean[expt,:], height[expt,:], color='blue', label='SW down')
        ax.plot(sw_up_mean[expt,:], height[expt,:], color='green', label='SW up')
        ax.plot(lw_down_mean[expt,:], height[expt,:], color='red', label='LW down')
        ax.plot(lw_up_mean[expt,:], height[expt,:], color='purple', label='LW up')
        ax.plot(sw_down2_mean[expt,:], height[expt,:], color='blue', linestyle='--', label='SW down (control)')
        ax.plot(sw_up2_mean[expt,:], height[expt,:], color='green', linestyle='--', label='SW up (control)')
        ax.plot(lw_down2_mean[expt,:], height[expt,:], color='red', linestyle='--', label='LW down (control)')
        ax.plot(lw_up2_mean[expt,:], height[expt,:], color='purple', linestyle='--', label='LW up (control)')
        ax.set_title(expt_label[expt])
        ax.set_ylabel('Approx height (km)')

    leg=plt.legend()
    leg.set_draggable(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_rfmip_out()
