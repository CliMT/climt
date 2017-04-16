cimport numpy as cnp
import numpy as np


#Flags
cdef cnp.int32_t rrtm_cloud_overlap_method
cdef cnp.int32_t rrtm_calculate_change_up_flux
cdef cnp.int32_t rrtm_cloud_props_flag
cdef cnp.int32_t rrtm_ice_props_flag
cdef cnp.int32_t rrtm_liq_droplet_flag
#Definitions of the functions we want to use within
#CliMT: initialisation and radiative transfer
cdef extern:
    void rrtmg_longwave_init(double *cp_d_air)

cdef extern:
    void rrtmg_set_constants(double *pi, double *grav, double *planck,
                            double *boltz, double *clight, double *avogadro,
                            double *loschmidt, double *gas_const,
                            double *stef_boltz, double *secs_per_day)

cdef extern:
    void rrtmg_longwave(cnp.int32_t *num_cols, cnp.int32_t *num_layers,
                       cnp.int32_t *cld_overlap_method, cnp.int32_t *calculate_change_up_flux,
                       double *layer_pressure, double *interface_pressure,
                       double *layer_temp, double *interface_temp, double *surface_temp,
                       double *h2o_vmr, double *o3_vmr, double *co2_vmr,
                       double *ch4_vmr, double *n2o_vmr, double *o2_vmr,
                       double *cfc11_vmr, double *cfc12_vmr, double *cfc22_vmr, double *ccl4_vmr,
                       double *surf_emissivity, cnp.int32_t *cloud_props_flag, cnp.int32_t *ice_props_flag,
                       cnp.int32_t *liq_droplet_flag, double *cloud_fraction,
                       double *cloud_tau, double *cloud_ice_path, double *cloud_liq_path,
                       double *cloud_ice_eff_size, double *cloud_droplet_eff_radius, double *aerosol_tau,
                       double *upward_longwave_flux, double *downward_longwave_flux,
                       double *longwave_heating_rate, double *up_lw_flux_clearsky,
                       double *down_lw_flux_clearsky, double *lw_heating_rate_clearsky,
                       double *duflx_dt, double *duflxc_dt)


#Use random cloud overlap model as default.
def initialise_rrtm_radiation(double cp_d_air, cnp.int32_t cloud_overlap_method=1,
                              cnp.int32_t calculate_change_up_flux=0,
                              cnp.int32_t cloud_props_flag=0, cnp.int32_t ice_props_flag=0,
                              cnp.int32_t liq_droplet_flag=0):

    global rrtm_cloud_overlap_method, rrtm_calculate_change_up_flux,\
        rrtm_cloud_props_flag, rrtm_ice_props_flag, rrtm_liq_droplet_flag

    rrtmg_longwave_init(&cp_d_air)

    rrtm_cloud_overlap_method = cloud_overlap_method
    rrtm_calculate_change_up_flux = calculate_change_up_flux
    rrtm_cloud_props_flag = cloud_props_flag
    rrtm_ice_props_flag = ice_props_flag
    rrtm_liq_droplet_flag = liq_droplet_flag

def set_constants(double pi, double grav, double planck, double boltz,
                  double clight, double avogadro, double loschmidt, double gas_const,
                  double stef_boltz, double secs_per_day):

    rrtmg_set_constants(&pi, &grav, &planck,
                        &boltz, &clight, &avogadro,
                        &loschmidt, &gas_const,
                        &stef_boltz, &secs_per_day)

def rrtm_calculate_longwave_fluxes(
    cnp.int32_t num_cols, cnp.int32_t num_layers,
    cnp.double_t[::1,:] layer_pressure,
    cnp.double_t[::1,:] interface_pressure,
    cnp.double_t[::1,:] layer_temp,
    cnp.double_t[::1,:] interface_temp,
    cnp.double_t[::1] surface_temp,
    cnp.double_t[::1,:] h2o_vmr,
    cnp.double_t[::1,:] o3_vmr,
    cnp.double_t[::1,:] co2_vmr,
    cnp.double_t[::1,:] ch4_vmr,
    cnp.double_t[::1,:] n2o_vmr,
    cnp.double_t[::1,:] o2_vmr,
    cnp.double_t[::1,:] cfc11_vmr,
    cnp.double_t[::1,:] cfc12_vmr,
    cnp.double_t[::1,:] cfc22_vmr,
    cnp.double_t[::1,:] ccl4_vmr,
    cnp.double_t[::1,:] surf_emissivity,
    cnp.double_t[::1,:] cloud_fraction,
    cnp.double_t[::1,:,:] aerosol_tau,
    cnp.double_t[::1,:] upward_longwave_flux,
    cnp.double_t[::1,:] downward_longwave_flux,
    cnp.double_t[::1,:] longwave_heating_rate,
    cnp.double_t[::1,:] up_lw_flux_clearsky,
    cnp.double_t[::1,:] down_lw_flux_clearsky,
    cnp.double_t[::1,:] lw_heating_rate_clearsky,
    cnp.double_t[::1,:,:] cloud_tau=None,
    cnp.double_t[::1,:] cloud_ice_path=None,
    cnp.double_t[::1,:] cloud_liq_path=None,
    cnp.double_t[::1,:] cloud_ice_eff_size=None,
    cnp.double_t[::1,:] cloud_droplet_eff_radius=None,
    cnp.double_t[::1,:] duflx_dt=None,
    cnp.double_t[::1,:] duflxc_dt=None):

    global rrtm_cloud_overlap_method, rrtm_calculate_change_up_flux,\
        rrtm_cloud_props_flag, rrtm_ice_props_flag, rrtm_liq_droplet_flag

    if rrtm_calculate_change_up_flux == 0:
        duflx_dt = np.zeros((1,1))
        duflxc_dt = np.zeros((1,1))

    rrtmg_longwave(
        &num_cols, &num_layers,
        &rrtm_cloud_overlap_method,
        &rrtm_calculate_change_up_flux,
        <double *>&layer_pressure[0,0],
        <double *>&interface_pressure[0,0],
        <double *>&layer_temp[0,0],
        <double *>&interface_temp[0,0],
        <double *>&surface_temp[0],
        <double *>&h2o_vmr[0,0],
        <double *>&o3_vmr[0,0],
        <double *>&co2_vmr[0,0],
        <double *>&ch4_vmr[0,0],
        <double *>&n2o_vmr[0,0],
        <double *>&o2_vmr[0,0],
        <double *>&cfc11_vmr[0,0],
        <double *>&cfc12_vmr[0,0],
        <double *>&cfc22_vmr[0,0],
        <double *>&ccl4_vmr[0,0],
        <double *>&surf_emissivity[0,0],
        &rrtm_cloud_props_flag,
        &rrtm_ice_props_flag,
        &rrtm_liq_droplet_flag,
        <double *>&cloud_fraction[0,0],
        <double *>&cloud_tau[0,0,0],
        <double *>&cloud_ice_path[0,0],
        <double *>&cloud_liq_path[0,0],
        <double *>&cloud_ice_eff_size[0,0],
        <double *>&cloud_droplet_eff_radius[0,0],
        <double *>&aerosol_tau[0,0,0],
        <double *>&upward_longwave_flux[0,0],
        <double *>&downward_longwave_flux[0,0],
        <double *>&longwave_heating_rate[0,0],
        <double *>&up_lw_flux_clearsky[0,0],
        <double *>&down_lw_flux_clearsky[0,0],
        <double *>&lw_heating_rate_clearsky[0,0],
        <double *>&duflx_dt[0,0],
        <double *>&duflxc_dt[0,0])

    return (upward_longwave_flux, downward_longwave_flux,
            longwave_heating_rate, up_lw_flux_clearsky,
            down_lw_flux_clearsky, lw_heating_rate_clearsky)
            #TODO ignore for now duflx_dt, duflxc_dt)
