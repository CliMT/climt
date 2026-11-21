"""cffi build script: one compiled module per Fortran library.

Run by setuptools through ``cffi_modules`` in ``setup.py`` (API, out-of-line
mode). A build-time file only: it lives in ``fortran/`` next to
``_build_config.py``, outside the installed package. Each ``*_ffi`` object
below compiles ``climt_fortran._<name>_cffi``, a small extension linked
against the static library ``_lib/<system>/lib*.a``. The prototypes are those
of the routines' ``bind(c)`` names; every argument is passed by pointer, as
Fortran expects.
"""

import os

import cffi
from _build_config import build_config  # fortran/, put on sys.path by setup.py

_CONFIG = build_config()


def _ffi(module_name, prototypes, static_library, extra_flags=()):
    """An FFI declaring ``prototypes`` that builds ``climt_fortran.<module_name>``
    linked against ``static_library`` from ``_lib/<system>/``. ``extra_flags``
    go to both the compiler and the linker (``-fopenmp`` for RRTMG, as the
    1.0 Cython extensions had)."""
    ffi = cffi.FFI()
    ffi.cdef(prototypes)
    ffi.set_source(
        "climt_fortran." + module_name,
        "#include <stdint.h>\n" + prototypes,
        # Relink when the Fortran library changes, not only when the C text does.
        depends=[os.path.join(_CONFIG["lib_dir"], static_library)],
        libraries=_CONFIG["libraries"],
        library_dirs=_CONFIG["library_dirs"],
        include_dirs=_CONFIG["include_dirs"],
        extra_compile_args=_CONFIG["extra_compile_args"] + list(extra_flags),
        extra_link_args=([os.path.join(_CONFIG["lib_dir"], static_library)] + list(extra_flags)
                         + _CONFIG["extra_link_args"]),
    )
    return ffi


EMANUEL_PROTOTYPES = """
void init_emanuel_convection_fortran(
    int32_t *pbl, int32_t *least_conv_level,
    double *thresh_water_level,
    double *crit_temp, double *entrain_coeff,
    double *downdraft_frac_area,
    double *precip_frac_outside_cloud,
    double *rain_speed, double *snow_speed,
    double *rain_evap_coeff, double *snow_evap_coeff,
    double *mom_tran_coeff, double *max_neg_temp_pert,
    double *beta, double *alpha, double *damp_amp,
    double *Cpd, double *Cpv, double *Cl,
    double *gas_const_vapour, double *gas_const_air,
    double *lat_heat, double *grav, double *density_water,
    double *reference_mass_flux_timescale);

void emanuel_convection(
    double *temp, double *q, double *qs,
    double *u, double *v,
    double *pmid, double *pint,
    int32_t *nlevs,
    int32_t *max_conv_lev,
    int32_t *num_tracers,
    double *dt, int32_t *conv_state,
    double *dtemp, double *dq, double *du,
    double *dv, double *precip,
    double *downdraft_vel_scale,
    double *downdraft_temp_scale,
    double *downdraft_q_scale,
    double *cloud_base_mass_flux,
    double *cape, double *tracers,
    double *dtracers);
"""

emanuel_ffi = _ffi("_emanuel_convection_cffi", EMANUEL_PROTOTYPES, "libemanuel.a")

RRTMG_LW_PROTOTYPES = """
void rrtmg_set_constants(
    double *pi, double *grav, double *planck, double *boltz,
    double *clight, double *avogadro, double *loschmidt,
    double *gas_const, double *stef_boltz, double *secs_per_day);

void rrtmg_lw_ini_wrapper(double *cp_d_air);

void mcica_subcol_lw_wrapper(
    int32_t *col_lon, int32_t *num_cols, int32_t *num_layers,
    int32_t *cld_overlap_method, int32_t *permute_seed,
    int32_t *random_number_generator, double *layer_pressure,
    double *cloud_fraction, double *cloud_ice_path, double *cloud_liq_path,
    double *cloud_ice_eff_size, double *cloud_droplet_eff_radius,
    double *cloud_tau,
    double *cloud_fraction_mcica, double *cloud_ice_path_mcica,
    double *cloud_liq_path_mcica, double *cloud_ice_eff_size_mcica,
    double *cloud_droplet_eff_radius_mcica, double *cloud_tau_mcica);

void rrtmg_lw_mcica_wrapper(
    int32_t *num_cols, int32_t *num_layers,
    int32_t *cld_overlap_method, int32_t *calculate_change_up_flux,
    double *layer_pressure, double *interface_pressure,
    double *layer_temp, double *interface_temp, double *surface_temp,
    double *h2o_vmr, double *o3_vmr, double *co2_vmr,
    double *ch4_vmr, double *n2o_vmr, double *o2_vmr,
    double *cfc11_vmr, double *cfc12_vmr, double *cfc22_vmr,
    double *ccl4_vmr, double *surf_emissivity,
    int32_t *cloud_props_flag, int32_t *ice_props_flag,
    int32_t *liq_droplet_flag, double *cloud_fraction_mcica,
    double *cloud_tau_mcica, double *cloud_ice_path_mcica,
    double *cloud_liq_path_mcica, double *cloud_ice_eff_size_mcica,
    double *cloud_droplet_eff_radius_mcica, double *aerosol_tau,
    double *upward_longwave_flux, double *downward_longwave_flux,
    double *longwave_heating_rate, double *up_lw_flux_clearsky,
    double *down_lw_flux_clearsky, double *lw_heating_rate_clearsky,
    double *duflx_dt, double *duflxc_dt);

void rrtmg_lw_nomcica_wrapper(
    int32_t *num_cols, int32_t *num_layers,
    int32_t *cld_overlap_method, int32_t *calculate_change_up_flux,
    double *layer_pressure, double *interface_pressure,
    double *layer_temp, double *interface_temp, double *surface_temp,
    double *h2o_vmr, double *o3_vmr, double *co2_vmr,
    double *ch4_vmr, double *n2o_vmr, double *o2_vmr,
    double *cfc11_vmr, double *cfc12_vmr, double *cfc22_vmr,
    double *ccl4_vmr, double *surf_emissivity,
    int32_t *cloud_props_flag, int32_t *ice_props_flag,
    int32_t *liq_droplet_flag, double *cloud_fraction,
    double *cloud_tau, double *cloud_ice_path,
    double *cloud_liq_path, double *cloud_ice_eff_size,
    double *cloud_droplet_eff_radius, double *aerosol_tau,
    double *upward_longwave_flux, double *downward_longwave_flux,
    double *longwave_heating_rate, double *up_lw_flux_clearsky,
    double *down_lw_flux_clearsky, double *lw_heating_rate_clearsky,
    double *duflx_dt, double *duflxc_dt);
"""

rrtmg_lw_ffi = _ffi("_rrtmg_lw_cffi", RRTMG_LW_PROTOTYPES, "librrtmg_lw.a", extra_flags=["-fopenmp"])

RRTMG_SW_PROTOTYPES = """
void rrtmg_sw_set_constants(
    double *pi, double *grav, double *planck, double *boltz,
    double *clight, double *avogadro, double *loschmidt,
    double *gas_const, double *stef_boltz, double *secs_per_day);

void rrtmg_sw_ini_wrapper(double *cp_d_air);

void mcica_subcol_sw_wrapper(
    int32_t *col_lon, int32_t *num_cols, int32_t *num_layers,
    int32_t *cld_overlap_method, int32_t *permute_seed,
    int32_t *rnd_gen_flag, double *layer_pressure,
    double *cloud_fraction, double *cloud_ice_path, double *cloud_liq_path,
    double *cloud_ice_eff_size, double *cloud_droplet_eff_radius,
    double *cloud_tau, double *cloud_single_scat_albedo,
    double *cloud_asym, double *cloud_fwd_scat_frac,
    double *cloud_fraction_mcica, double *cloud_ice_path_mcica,
    double *cloud_liq_path_mcica, double *cloud_ice_eff_size_mcica,
    double *cloud_droplet_eff_radius_mcica, double *cloud_tau_mcica,
    double *cloud_single_scat_albedo_mcica, double *cloud_asym_mcica,
    double *cloud_fwd_scat_frac_mcica);

void rrtmg_sw_mcica_wrapper(
    int32_t *num_cols, int32_t *num_layers,
    int32_t *cld_overlap_method, int32_t *aerosol_input_format,
    double *layer_pressure, double *interface_pressure, double *layer_temp,
    double *interface_temp, double *surface_temp, double *h2o_vmr,
    double *o3_vmr, double *co2_vmr, double *ch4_vmr, double *n2o_vmr,
    double *o2_vmr, double *surface_direct_sw_albedo,
    double *surface_diffuse_sw_albedo, double *surface_direct_nir_albedo,
    double *surface_diffuse_nir_albedo, double *cosine_zenith_angle,
    double *flux_adj_earth_sun_dist, int32_t *day_of_year,
    double *solar_constant, int32_t *solar_variability_type,
    int32_t *cloud_props_flag, int32_t *ice_props_flag,
    int32_t *liq_droplet_flag, double *cloud_fraction_mcica,
    double *cloud_tau_mcica, double *cloud_single_scat_albedo_mcica,
    double *cloud_asym_mcica, double *cloud_fwd_scat_frac_mcica,
    double *cloud_ice_path_mcica, double *cloud_liq_path_mcica,
    double *cloud_ice_eff_size_mcica,
    double *cloud_droplet_eff_radius_mcica, double *aerosol_tau,
    double *aer_single_scat_albedo, double *aer_asym,
    double *aer_tau_at_55_micron,
    double *upward_shortwave_flux, double *downward_shortwave_flux,
    double *shortwave_heating_rate, double *up_sw_flux_clearsky,
    double *down_sw_flux_clearsky, double *sw_heating_rate_clearsky,
    double *sol_var_scale_factors, double *fac_sunspot_ampl,
    double *solar_cycle_frac);

void rrtmg_sw_nomcica_wrapper(
    int32_t *num_cols, int32_t *num_layers,
    int32_t *cld_overlap_method, int32_t *aerosol_input_format,
    double *layer_pressure, double *interface_pressure, double *layer_temp,
    double *interface_temp, double *surface_temp, double *h2o_vmr,
    double *o3_vmr, double *co2_vmr, double *ch4_vmr, double *n2o_vmr,
    double *o2_vmr, double *surface_direct_sw_albedo,
    double *surface_diffuse_sw_albedo, double *surface_direct_nir_albedo,
    double *surface_diffuse_nir_albedo, double *cosine_zenith_angle,
    double *flux_adj_earth_sun_dist, int32_t *day_of_year,
    double *solar_constant, int32_t *solar_variability_type,
    int32_t *cloud_props_flag, int32_t *ice_props_flag,
    int32_t *liq_droplet_flag, double *cloud_fraction,
    double *cloud_tau, double *cloud_single_scat_albedo,
    double *cloud_asym, double *cloud_fwd_scat_frac,
    double *cloud_ice_path, double *cloud_liq_path,
    double *cloud_ice_eff_size, double *cloud_droplet_eff_radius,
    double *aerosol_tau, double *aer_single_scat_albedo, double *aer_asym,
    double *aer_tau_at_55_micron,
    double *upward_shortwave_flux, double *downward_shortwave_flux,
    double *shortwave_heating_rate, double *up_sw_flux_clearsky,
    double *down_sw_flux_clearsky, double *sw_heating_rate_clearsky,
    double *sol_var_scale_factors, double *fac_sunspot_ampl,
    double *solar_cycle_frac);
"""

rrtmg_sw_ffi = _ffi("_rrtmg_sw_cffi", RRTMG_SW_PROTOTYPES, "librrtmg_sw.a", extra_flags=["-fopenmp"])
