! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Run the Socrates radiative transfer code

module socrates_runes

use rad_pcf, only: &
  ip_source_illuminate                      => ip_solar, &
  ip_source_thermal                         => ip_infra_red, &
  ip_cloud_representation_off               => ip_cloud_off, &
  ip_cloud_representation_ice_water         => ip_cloud_ice_water, &
  ip_cloud_representation_combine_ice_water => ip_cloud_combine_ice_water, &
  ip_cloud_representation_csiw              => ip_cloud_csiw, &
  ip_cloud_representation_split_ice_water   => ip_cloud_split_ice_water, &
  ip_overlap_max_random                     => ip_max_rand, &
  ip_overlap_random                         => ip_rand, &
  ip_overlap_exponential_random             => ip_exponential_rand, &
  ip_inhom_homogeneous                      => ip_homogeneous, &
  ip_inhom_scaling                          => ip_scaling, &
  ip_inhom_mcica                            => ip_mcica, &
  ip_inhom_cairns                           => ip_cairns, &
  ip_inhom_tripleclouds_2019                => ip_tripleclouds_2019, &
  ip_cloud_entrapment_zero                  => ip_zero_entrapment, &
  ip_cloud_entrapment_max                   => ip_max_entrapment, &
  ip_droplet_re_external                    => ip_re_external, &
  ip_droplet_re_constant                    => ip_re_constant, &
  ip_droplet_re_liu                         => ip_re_liu, &
  ip_droplet_re_default                     => ip_re_default, &
  ip_scatter_full, &
  ip_scatter_none                           => ip_no_scatter_abs, &
  ip_scatter_approx, &
  ip_scatter_hybrid, &
  ip_gas_overlap_hybrid                     => ip_overlap_hybrid, &
  ip_gas_overlap_random                     => ip_overlap_random, &
  ip_gas_overlap_exact_major                => ip_overlap_exact_major, &
  ip_gas_overlap_k_eqv                      => ip_overlap_k_eqv_scl

use socrates_def_diag, only: StrDiag

implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_RUNES'
contains

subroutine runes(n_profile, n_layer, diag, &
  spectrum, spectrum_name, mcica_data, &
  profile_list, n_layer_stride, n_level_stride, &
  n_cloud_layer, n_aer_mode, n_aer_layer, n_tile, n_subcol_gen, &
  p_layer, t_layer, p_level, t_level, r_layer, r_level, mass, density, &
  h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
  cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
  he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
  n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
  hdo, hcl, hf, cosso, tosso, yosos, &
  ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
  hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2, &
  p_layer_1d, t_layer_1d, p_level_1d, t_level_1d, r_layer_1d, r_level_1d, &
  mass_1d, density_1d, &
  h2o_1d, co2_1d, o3_1d, n2o_1d, co_1d, ch4_1d, o2_1d, no_1d, so2_1d, no2_1d, &
  nh3_1d, hno3_1d, n2_1d, cfc11_1d, cfc12_1d, cfc113_1d, hcfc22_1d, hfc125_1d, &
  hfc134a_1d, cfc114_1d, tio_1d, vo_1d, h2_1d, he_1d, ocs_1d, na_1d, k_1d, &
  feh_1d, crh_1d, li_1d, rb_1d, cs_1d, ph3_1d, c2h2_1d, hcn_1d, h2s_1d, ar_1d, &
  o_1d, n_1d, no3_1d, n2o5_1d, hono_1d, ho2no2_1d, h2o2_1d, c2h6_1d, ch3_1d, &
  h2co_1d, ho2_1d, hdo_1d, hcl_1d, hf_1d, cosso_1d, tosso_1d, yosos_1d, &
  ch3cho_1d, ch3ooh_1d, ch3coch3_1d, ch3cocho_1d, chocho_1d, c2h5cho_1d, &
  hoch2cho_1d, c2h5coch3_1d, mvk_1d, macr_1d, pan_1d, ch3ono2_1d, &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, co_mix_ratio, &
  ch4_mix_ratio, o2_mix_ratio, no_mix_ratio, so2_mix_ratio, no2_mix_ratio, &
  nh3_mix_ratio, hno3_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, &
  cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, hfc125_mix_ratio, &
  hfc134a_mix_ratio, cfc114_mix_ratio, tio_mix_ratio, vo_mix_ratio, &
  h2_mix_ratio, he_mix_ratio, ocs_mix_ratio, na_mix_ratio, k_mix_ratio, &
  feh_mix_ratio, crh_mix_ratio, li_mix_ratio, rb_mix_ratio, cs_mix_ratio, &
  ph3_mix_ratio, c2h2_mix_ratio, hcn_mix_ratio, h2s_mix_ratio, ar_mix_ratio, &
  o_mix_ratio, n_mix_ratio, no3_mix_ratio, n2o5_mix_ratio, hono_mix_ratio, &
  ho2no2_mix_ratio, h2o2_mix_ratio, c2h6_mix_ratio, ch3_mix_ratio, &
  h2co_mix_ratio, ho2_mix_ratio, hdo_mix_ratio, hcl_mix_ratio, hf_mix_ratio, &
  cosso_mix_ratio, tosso_mix_ratio, yosos_mix_ratio, ch3cho_mix_ratio, &
  ch3ooh_mix_ratio, ch3coch3_mix_ratio, ch3cocho_mix_ratio, chocho_mix_ratio, &
  c2h5cho_mix_ratio, hoch2cho_mix_ratio, c2h5coch3_mix_ratio, mvk_mix_ratio, &
  macr_mix_ratio, pan_mix_ratio, ch3ono2_mix_ratio, &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_co_well_mixed, l_ch4_well_mixed, l_o2_well_mixed, l_no_well_mixed, &
  l_so2_well_mixed, l_no2_well_mixed, l_nh3_well_mixed, l_hno3_well_mixed, &
  l_n2_well_mixed, l_cfc11_well_mixed, l_cfc12_well_mixed, &
  l_cfc113_well_mixed, l_hcfc22_well_mixed, l_hfc125_well_mixed, &
  l_hfc134a_well_mixed, l_cfc114_well_mixed, l_tio_well_mixed, &
  l_vo_well_mixed, l_h2_well_mixed, l_he_well_mixed, l_ocs_well_mixed, &
  l_na_well_mixed, l_k_well_mixed, l_feh_well_mixed, l_crh_well_mixed, &
  l_li_well_mixed, l_rb_well_mixed, l_cs_well_mixed, l_ph3_well_mixed, &
  l_c2h2_well_mixed, l_hcn_well_mixed, l_h2s_well_mixed, l_ar_well_mixed, &
  l_o_well_mixed, l_n_well_mixed, l_no3_well_mixed, l_n2o5_well_mixed, &
  l_hono_well_mixed, l_ho2no2_well_mixed, l_h2o2_well_mixed, &
  l_c2h6_well_mixed, l_ch3_well_mixed, l_h2co_well_mixed, l_ho2_well_mixed, &
  l_hdo_well_mixed, l_hcl_well_mixed, l_hf_well_mixed, l_cosso_well_mixed, &
  l_tosso_well_mixed, l_yosos_well_mixed, l_ch3cho_well_mixed, &
  l_ch3ooh_well_mixed, l_ch3coch3_well_mixed, l_ch3cocho_well_mixed, &
  l_chocho_well_mixed, l_c2h5cho_well_mixed, l_hoch2cho_well_mixed, &
  l_c2h5coch3_well_mixed, l_mvk_well_mixed, l_macr_well_mixed, &
  l_pan_well_mixed, l_ch3ono2_well_mixed, &
  l_flux_ground, t_ground, flux_ground, flux_ground_1d, &
  cos_zenith_angle, solar_irrad, l_orog, orog_corr, &
  l_grey_albedo, grey_albedo, albedo_diff, albedo_dir, &
  albedo_diff_1d, albedo_dir_1d, &
  l_tile, l_flux_tile, &
  frac_tile, t_tile, flux_tile, albedo_diff_tile, albedo_dir_tile, &
  frac_tile_1d, t_tile_1d, flux_tile_1d, &
  albedo_diff_tile_1d, albedo_dir_tile_1d, &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  liq_rsd, ice_rsd, liq_conv_rsd, ice_conv_rsd, &
  liq_nc, ice_nc, liq_conv_nc, ice_conv_nc, &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d, &
  liq_rsd_1d, ice_rsd_1d, liq_conv_rsd_1d, ice_conv_rsd_1d, &
  liq_nc_1d, ice_nc_1d, liq_conv_nc_1d, ice_conv_nc_1d, &
  liq_dim_constant, liq_conv_dim_constant, &
  cloud_vertical_decorr, conv_vertical_decorr, &
  cloud_horizontal_rsd, &
  liq_dim_aparam, liq_dim_bparam, &
  layer_heat_capacity, layer_heat_capacity_1d, &
  i_source, i_scatter_method, i_gas_overlap, &
  i_cloud_representation, i_overlap, i_inhom, &
  i_cloud_entrapment, i_mcica_sampling, &
  i_st_water, i_cnv_water, i_st_ice, i_cnv_ice, i_drop_re, &
  rand_seed, &
  l_rayleigh, l_mixing_ratio, l_aerosol_mode, &
  aer_mix_ratio, aer_absorption, aer_scattering, aer_asymmetry, &
  aer_mix_ratio_1d, aer_absorption_1d, aer_scattering_1d, aer_asymmetry_1d, &
  mean_rel_humidity, mean_rel_humidity_1d, &
  l_water_soluble, water_soluble, water_soluble_1d, &
  l_dust_like, dust_like, dust_like_1d, &
  l_oceanic, oceanic, oceanic_1d, &
  l_soot, soot, soot_1d, &
  l_ash, ash, ash_1d, &
  l_sulphuric, sulphuric, sulphuric_1d, &
  l_ammonium_sulphate, ammonium_sulphate, ammonium_sulphate_1d, &
  l_saharan_dust, saharan_dust, saharan_dust_1d, &
  l_accum_sulphate, accum_sulphate, accum_sulphate_1d, &
  l_aitken_sulphate, aitken_sulphate, aitken_sulphate_1d, &
  l_fresh_soot, fresh_soot, fresh_soot_1d, &
  l_aged_soot, aged_soot, aged_soot_1d, &
  l_sodium_chloride, sodium_chloride, sodium_chloride_1d, &
  l_seasalt_film, seasalt_film, seasalt_film_1d, &
  l_seasalt_jet, seasalt_jet, seasalt_jet_1d, &
  l_dust_div1, dust_div1, dust_div1_1d, &
  l_dust_div2, dust_div2, dust_div2_1d, &
  l_dust_div3, dust_div3, dust_div3_1d, &
  l_dust_div4, dust_div4, dust_div4_1d, &
  l_dust_div5, dust_div5, dust_div5_1d, &
  l_dust_div6, dust_div6, dust_div6_1d, &
  l_biomass_1, biomass_1, biomass_1_1d, &
  l_biomass_2, biomass_2, biomass_2_1d, &
  l_biogenic, biogenic, biogenic_1d, &
  l_ocff_fresh, ocff_fresh, ocff_fresh_1d, &
  l_ocff_aged, ocff_aged, ocff_aged_1d, &
  l_delta, delta, delta_1d, &
  l_murk, murk, murk_1d, &
  l_nitrate, nitrate, nitrate_1d, &
  l_twobindust_1, twobindust_1, twobindust_1_1d, &
  l_twobindust_2, twobindust_2, twobindust_2_1d, &
  l_invert, l_profile_last, l_debug, i_profile_debug)


use def_spectrum, only: StrSpecData
use def_mcica,    only: StrMcica
use def_control,  only: StrCtrl,  deallocate_control
use def_dimen,    only: StrDim
use def_atm,      only: StrAtm,   deallocate_atm
use def_bound,    only: StrBound, deallocate_bound
use def_cld,      only: StrCld,   deallocate_cld, deallocate_cld_prsc, &
                                  deallocate_cld_mcica
use def_aer,      only: StrAer,   deallocate_aer, deallocate_aer_prsc
use def_out,      only: StrOut,   deallocate_out

use socrates_set_spectrum, only: spectrum_array_name, spectrum_array, &
                                 mcica_spectrum_name, mcica_data_array

use socrates_set_control,   only: set_control
use socrates_set_dimen,     only: set_dimen
use socrates_set_atm,       only: set_atm
use socrates_set_bound,     only: set_bound
use socrates_set_cld,       only: set_cld
use socrates_set_cld_dim,   only: set_cld_dim
use socrates_set_cld_mcica, only: set_cld_mcica
use socrates_set_aer,       only: set_aer
use socrates_set_diag,      only: set_diag

use realtype_rd, only: RealExt
use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength
use rad_pcf, only: i_normal, i_err_fatal

implicit none

! Output diagnostic fields
type(StrDiag), intent(inout) :: diag

! Spectral data:
type (StrSpecData), intent(in), target, optional :: spectrum
character(len=*), intent(in), optional :: spectrum_name

! Mcica data
type (StrMcica), intent(in), target, optional :: mcica_data

integer, intent(in) :: n_profile
!   Number of columns to operate on
integer, intent(in) :: n_layer
!   Number of layers for radiation
integer, intent(in), optional :: profile_list(:)
!   List of profiles to use from input fields
integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays
integer, intent(in), optional :: n_level_stride
!   Number of levels in input 1d arrays
integer, intent(in), optional :: n_tile
!   Number of surface tiles
integer, intent(in), optional :: n_cloud_layer
!   Number of potentially cloudy layers
integer, intent(in), optional :: n_aer_mode
!   Number of aerosol modes
integer, intent(in), optional :: n_aer_layer
!   Number of aerosol layers in 1d arrays
integer, intent(in), optional :: n_subcol_gen
!   Number of sub-columns to generate for cloud sampling

real(RealExt), intent(in), optional :: p_layer(:, :)
real(RealExt), intent(in), optional :: p_layer_1d(:)
!   Pressure at layer centres
real(RealExt), intent(in), optional :: t_layer(:, :)
real(RealExt), intent(in), optional :: t_layer_1d(:)
!   Temperature at layer centres
real(RealExt), intent(in), optional :: p_level(:, :)
real(RealExt), intent(in), optional :: p_level_1d(:)
!   Pressure at layer boundaries
real(RealExt), intent(in), optional :: t_level(:, :)
real(RealExt), intent(in), optional :: t_level_1d(:)
!   Temperature at layer boundaries
real(RealExt), intent(in), optional :: r_layer(:, :)
real(RealExt), intent(in), optional :: r_layer_1d(:)
!   Radius (height from centre of planet) at layer centres
real(RealExt), intent(in), optional :: r_level(:, :)
real(RealExt), intent(in), optional :: r_level_1d(:)
!   Radius (height from centre of planet) at layer boundaries
real(RealExt), intent(in), optional :: mass(:, :)
real(RealExt), intent(in), optional :: mass_1d(:)
!   Mass of layer (kg m-2)
real(RealExt), intent(in), optional :: density(:, :)
real(RealExt), intent(in), optional :: density_1d(:)
!   Density of layer (kg m-3)
real(RealExt), intent(in), optional :: h2o(:, :)
real(RealExt), intent(in), optional :: h2o_1d(:)
!   Mass mixing ratio of water vapour
real(RealExt), intent(in), optional :: co2(:, :)
real(RealExt), intent(in), optional :: co2_1d(:)
!   Mass mixing ratio of carbon dioxide
real(RealExt), intent(in), optional :: o3(:, :)
real(RealExt), intent(in), optional :: o3_1d(:)
!   Mass mixing ratio of ozone
real(RealExt), intent(in), optional :: n2o(:, :)
real(RealExt), intent(in), optional :: n2o_1d(:)
!   Mass mixing ratio of dinitrogen oxide
real(RealExt), intent(in), optional :: co(:, :)
real(RealExt), intent(in), optional :: co_1d(:)
!   Mass mixing ratio of carbon monoxide
real(RealExt), intent(in), optional :: ch4(:, :)
real(RealExt), intent(in), optional :: ch4_1d(:)
!   Mass mixing ratio of methane
real(RealExt), intent(in), optional :: o2(:, :)
real(RealExt), intent(in), optional :: o2_1d(:)
!   Mass mixing ratio of oxygen
real(RealExt), intent(in), optional :: no(:, :)
real(RealExt), intent(in), optional :: no_1d(:)
!   Mass mixing ratio of nitrogen monoxide
real(RealExt), intent(in), optional :: so2(:, :)
real(RealExt), intent(in), optional :: so2_1d(:)
!   Mass mixing ratio of sulphur dioxide
real(RealExt), intent(in), optional :: no2(:, :)
real(RealExt), intent(in), optional :: no2_1d(:)
!   Mass mixing ratio of nitrogen dioxide
real(RealExt), intent(in), optional :: nh3(:, :)
real(RealExt), intent(in), optional :: nh3_1d(:)
!   Mass mixing ratio of ammonia
real(RealExt), intent(in), optional :: hno3(:, :)
real(RealExt), intent(in), optional :: hno3_1d(:)
!   Mass mixing ratio of nitric acid
real(RealExt), intent(in), optional :: n2(:, :)
real(RealExt), intent(in), optional :: n2_1d(:)
!   Mass mixing ratio of nitrogen
real(RealExt), intent(in), optional :: cfc11(:, :)
real(RealExt), intent(in), optional :: cfc11_1d(:)
!   Mass mixing ratio of CFC11 (CFCl3)
real(RealExt), intent(in), optional :: cfc12(:, :)
real(RealExt), intent(in), optional :: cfc12_1d(:)
!   Mass mixing ratio of CFC12 (CF2Cl2)
real(RealExt), intent(in), optional :: cfc113(:, :)
real(RealExt), intent(in), optional :: cfc113_1d(:)
!   Mass mixing ratio of CFC113 (CF2ClCFCl2)
real(RealExt), intent(in), optional :: hcfc22(:, :)
real(RealExt), intent(in), optional :: hcfc22_1d(:)
!   Mass mixing ratio of HCFC22 (CHF2Cl)
real(RealExt), intent(in), optional :: hfc125(:, :)
real(RealExt), intent(in), optional :: hfc125_1d(:)
!   Mass mixing ratio of HFC125 (C2HF5)
real(RealExt), intent(in), optional :: hfc134a(:, :)
real(RealExt), intent(in), optional :: hfc134a_1d(:)
!   Mass mixing ratio of HFC134A (CF4CFH2)
real(RealExt), intent(in), optional :: cfc114(:, :)
real(RealExt), intent(in), optional :: cfc114_1d(:)
!   Mass mixing ratio of CFC114 (C2Cl2F4)
real(RealExt), intent(in), optional :: tio(:, :)
real(RealExt), intent(in), optional :: tio_1d(:)
!   Mass mixing ratio of titanium oxide
real(RealExt), intent(in), optional :: vo(:, :)
real(RealExt), intent(in), optional :: vo_1d(:)
!   Mass mixing ratio of vanadium oxide
real(RealExt), intent(in), optional :: h2(:, :)
real(RealExt), intent(in), optional :: h2_1d(:)
!   Mass mixing ratio of hydrogen
real(RealExt), intent(in), optional :: he(:, :)
real(RealExt), intent(in), optional :: he_1d(:)
!   Mass mixing ratio of helium
real(RealExt), intent(in), optional :: ocs(:, :)
real(RealExt), intent(in), optional :: ocs_1d(:)
!   Mass mixing ratio of carbonyl sulphide
real(RealExt), intent(in), optional :: na(:, :)
real(RealExt), intent(in), optional :: na_1d(:)
!   Mass mixing ratio of sodium
real(RealExt), intent(in), optional :: k(:, :)
real(RealExt), intent(in), optional :: k_1d(:)
!   Mass mixing ratio of potassium
real(RealExt), intent(in), optional :: feh(:, :)
real(RealExt), intent(in), optional :: feh_1d(:)
!   Mass mixing ratio of iron hydride
real(RealExt), intent(in), optional :: crh(:, :)
real(RealExt), intent(in), optional :: crh_1d(:)
!   Mass mixing ratio of chromium hydride
real(RealExt), intent(in), optional :: li(:, :)
real(RealExt), intent(in), optional :: li_1d(:)
!   Mass mixing ratio of lithium
real(RealExt), intent(in), optional :: rb(:, :)
real(RealExt), intent(in), optional :: rb_1d(:)
!   Mass mixing ratio of rubidium
real(RealExt), intent(in), optional :: cs(:, :)
real(RealExt), intent(in), optional :: cs_1d(:)
!   Mass mixing ratio of cesium
real(RealExt), intent(in), optional :: ph3(:, :)
real(RealExt), intent(in), optional :: ph3_1d(:)
!   Mass mixing ratio of phosphine
real(RealExt), intent(in), optional :: c2h2(:, :)
real(RealExt), intent(in), optional :: c2h2_1d(:)
!   Mass mixing ratio of acetylene
real(RealExt), intent(in), optional :: hcn(:, :)
real(RealExt), intent(in), optional :: hcn_1d(:)
!   Mass mixing ratio of hydrogen cyanide
real(RealExt), intent(in), optional :: h2s(:, :)
real(RealExt), intent(in), optional :: h2s_1d(:)
!   Mass mixing ratio of hydrogen sulphide
real(RealExt), intent(in), optional :: ar(:, :)
real(RealExt), intent(in), optional :: ar_1d(:)
!   Mass mixing ratio of argon
real(RealExt), intent(in), optional :: o(:, :)
real(RealExt), intent(in), optional :: o_1d(:)
!   Mass mixing ratio of atomic oxygen
real(RealExt), intent(in), optional :: n(:, :)
real(RealExt), intent(in), optional :: n_1d(:)
!   Mass mixing ratio of atomic nitrogen
real(RealExt), intent(in), optional :: no3(:, :)
real(RealExt), intent(in), optional :: no3_1d(:)
!   Mass mixing ratio of nitrate radical
real(RealExt), intent(in), optional :: n2o5(:, :)
real(RealExt), intent(in), optional :: n2o5_1d(:)
!   Mass mixing ratio of dinitrogen pentoxide
real(RealExt), intent(in), optional :: hono(:, :)
real(RealExt), intent(in), optional :: hono_1d(:)
!   Mass mixing ratio of nitrous acid
real(RealExt), intent(in), optional :: ho2no2(:, :)
real(RealExt), intent(in), optional :: ho2no2_1d(:)
!   Mass mixing ratio of peroxynitric acid
real(RealExt), intent(in), optional :: h2o2(:, :)
real(RealExt), intent(in), optional :: h2o2_1d(:)
!   Mass mixing ratio of hydrogen peroxide
real(RealExt), intent(in), optional :: c2h6(:, :)
real(RealExt), intent(in), optional :: c2h6_1d(:)
!   Mass mixing ratio of ethane
real(RealExt), intent(in), optional :: ch3(:, :)
real(RealExt), intent(in), optional :: ch3_1d(:)
!   Mass mixing ratio of methyl radical
real(RealExt), intent(in), optional :: h2co(:, :)
real(RealExt), intent(in), optional :: h2co_1d(:)
!   Mass mixing ratio of formaldehyde
real(RealExt), intent(in), optional :: ho2(:, :)
real(RealExt), intent(in), optional :: ho2_1d(:)
!   Mass mixing ratio of hydroperoxy radical
real(RealExt), intent(in), optional :: hdo(:, :)
real(RealExt), intent(in), optional :: hdo_1d(:)
!   Mass mixing ratio of semiheavy water
real(RealExt), intent(in), optional :: hcl(:, :)
real(RealExt), intent(in), optional :: hcl_1d(:)
!   Mass mixing ratio of hydrogen chloride
real(RealExt), intent(in), optional :: hf(:, :)
real(RealExt), intent(in), optional :: hf_1d(:)
!   Mass mixing ratio of hydrogen fluoride
real(RealExt), intent(in), optional :: cosso(:, :)
real(RealExt), intent(in), optional :: cosso_1d(:)
!   Mass mixing ratio of cis-OSSO
real(RealExt), intent(in), optional :: tosso(:, :)
real(RealExt), intent(in), optional :: tosso_1d(:)
!   Mass mixing ratio of trans-OSSO
real(RealExt), intent(in), optional :: yosos(:, :)
real(RealExt), intent(in), optional :: yosos_1d(:)
!   Mass mixing ratio of OSO-S
real(RealExt), intent(in), optional :: ch3cho(:, :)
real(RealExt), intent(in), optional :: ch3cho_1d(:)
!   Mass mixing ratio of acetaldehyde
real(RealExt), intent(in), optional :: ch3ooh(:, :)
real(RealExt), intent(in), optional :: ch3ooh_1d(:)
!   Mass mixing ratio of methylhydroperoxide
real(RealExt), intent(in), optional :: ch3coch3(:, :)
real(RealExt), intent(in), optional :: ch3coch3_1d(:)
!   Mass mixing ratio of acetone
real(RealExt), intent(in), optional :: ch3cocho(:, :)
real(RealExt), intent(in), optional :: ch3cocho_1d(:)
!   Mass mixing ratio of methylglyoxal (MGLY)
real(RealExt), intent(in), optional :: chocho(:, :)
real(RealExt), intent(in), optional :: chocho_1d(:)
!   Mass mixing ratio of glyoxal
real(RealExt), intent(in), optional :: c2h5cho(:, :)
real(RealExt), intent(in), optional :: c2h5cho_1d(:)
!   Mass mixing ratio of propanal
real(RealExt), intent(in), optional :: hoch2cho(:, :)
real(RealExt), intent(in), optional :: hoch2cho_1d(:)
!   Mass mixing ratio of glycolaldehyde
real(RealExt), intent(in), optional :: c2h5coch3(:, :)
real(RealExt), intent(in), optional :: c2h5coch3_1d(:)
!   Mass mixing ratio of methyl ethyl ketone (MEK)
real(RealExt), intent(in), optional :: mvk(:, :)
real(RealExt), intent(in), optional :: mvk_1d(:)
!   Mass mixing ratio of methyl vinyl ketone (MVK)
real(RealExt), intent(in), optional :: macr(:, :)
real(RealExt), intent(in), optional :: macr_1d(:)
!   Mass mixing ratio of methacrolein (MACR)
real(RealExt), intent(in), optional :: pan(:, :)
real(RealExt), intent(in), optional :: pan_1d(:)
!   Mass mixing ratio of peroxyacetyl nitrate (PAN)
real(RealExt), intent(in), optional :: ch3ono2(:, :)
real(RealExt), intent(in), optional :: ch3ono2_1d(:)
!   Mass mixing ratio of methylnitrate

real(RealExt), intent(in), optional :: &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, co_mix_ratio, &
  ch4_mix_ratio, o2_mix_ratio, no_mix_ratio, so2_mix_ratio, no2_mix_ratio, &
  nh3_mix_ratio, hno3_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, &
  cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, hfc125_mix_ratio, &
  hfc134a_mix_ratio, cfc114_mix_ratio, tio_mix_ratio, vo_mix_ratio, &
  h2_mix_ratio, he_mix_ratio, ocs_mix_ratio, na_mix_ratio, k_mix_ratio, &
  feh_mix_ratio, crh_mix_ratio, li_mix_ratio, rb_mix_ratio, cs_mix_ratio, &
  ph3_mix_ratio, c2h2_mix_ratio, hcn_mix_ratio, h2s_mix_ratio, ar_mix_ratio, &
  o_mix_ratio, n_mix_ratio, no3_mix_ratio, n2o5_mix_ratio, hono_mix_ratio, &
  ho2no2_mix_ratio, h2o2_mix_ratio, c2h6_mix_ratio, ch3_mix_ratio, &
  h2co_mix_ratio, ho2_mix_ratio, hdo_mix_ratio, hcl_mix_ratio, hf_mix_ratio, &
  cosso_mix_ratio, tosso_mix_ratio, yosos_mix_ratio, ch3cho_mix_ratio, &
  ch3ooh_mix_ratio, ch3coch3_mix_ratio, ch3cocho_mix_ratio, chocho_mix_ratio, &
  c2h5cho_mix_ratio, hoch2cho_mix_ratio, c2h5coch3_mix_ratio, mvk_mix_ratio, &
  macr_mix_ratio, pan_mix_ratio, ch3ono2_mix_ratio
!   Well mixed mass mixing ratios

logical, intent(in), optional :: &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_co_well_mixed, l_ch4_well_mixed, l_o2_well_mixed, l_no_well_mixed, &
  l_so2_well_mixed, l_no2_well_mixed, l_nh3_well_mixed, l_hno3_well_mixed, &
  l_n2_well_mixed, l_cfc11_well_mixed, l_cfc12_well_mixed, &
  l_cfc113_well_mixed, l_hcfc22_well_mixed, l_hfc125_well_mixed, &
  l_hfc134a_well_mixed, l_cfc114_well_mixed, l_tio_well_mixed, &
  l_vo_well_mixed, l_h2_well_mixed, l_he_well_mixed, l_ocs_well_mixed, &
  l_na_well_mixed, l_k_well_mixed, l_feh_well_mixed, l_crh_well_mixed, &
  l_li_well_mixed, l_rb_well_mixed, l_cs_well_mixed, l_ph3_well_mixed, &
  l_c2h2_well_mixed, l_hcn_well_mixed, l_h2s_well_mixed, l_ar_well_mixed, &
  l_o_well_mixed, l_n_well_mixed, l_no3_well_mixed, l_n2o5_well_mixed, &
  l_hono_well_mixed, l_ho2no2_well_mixed, l_h2o2_well_mixed, &
  l_c2h6_well_mixed, l_ch3_well_mixed, l_h2co_well_mixed, l_ho2_well_mixed, &
  l_hdo_well_mixed, l_hcl_well_mixed, l_hf_well_mixed, l_cosso_well_mixed, &
  l_tosso_well_mixed, l_yosos_well_mixed, l_ch3cho_well_mixed, &
  l_ch3ooh_well_mixed, l_ch3coch3_well_mixed, l_ch3cocho_well_mixed, &
  l_chocho_well_mixed, l_c2h5cho_well_mixed, l_hoch2cho_well_mixed, &
  l_c2h5coch3_well_mixed, l_mvk_well_mixed, l_macr_well_mixed, &
  l_pan_well_mixed, l_ch3ono2_well_mixed
!   Flag to use the well mixed ratios

logical, intent(in), optional :: l_flux_ground
!   Set effective surface emission over whole grid-box
real(RealExt), intent(in), optional :: t_ground(:)
!   Effective radiative temperature over whole grid-box
real(RealExt), intent(in), optional :: flux_ground(:, :)
!   Effective surface emission over whole grid-box (n_profile, n_band)
real(RealExt), intent(in), optional :: flux_ground_1d(:)
!   1d Effective surface emission over whole grid-box (n_band)
real(RealExt), intent(in), optional :: cos_zenith_angle(:)
!   Cosine of solar zenith angle
real(RealExt), intent(in), optional :: solar_irrad(:)
!   Solar irradiance at top-of-atmosphere (mean over timestep)

logical, intent(in), optional :: l_orog
!   Apply orographic correction
real(RealExt), intent(in), optional :: orog_corr(:)
!   Orographic correction factor

logical, intent(in), optional :: l_grey_albedo
!   Set a single grey albedo / emissivity for the surface
real(RealExt), intent(in), optional :: grey_albedo
!   Grey surface albedo

real(RealExt), intent(in), optional :: albedo_diff(:, :)
!   Spectral diffuse albedo (n_profile, n_band)
real(RealExt), intent(in), optional :: albedo_dir(:, :)
!   Spectral direct albedo (n_profile, n_band)
real(RealExt), intent(in), optional :: albedo_diff_1d(:)
!   1d spectral diffuse albedo (n_band)
real(RealExt), intent(in), optional :: albedo_dir_1d(:)
!   1d spectral direct albedo (n_band)

logical, intent(in), optional :: l_tile
!   Use tiled surface properties
logical, intent(in), optional :: l_flux_tile(:)
!   Set effective surface emission for selected tiles
real(RealExt), intent(in), optional :: frac_tile(:, :)
!   Tile fractions (n_profile, n_tile)
real(RealExt), intent(in), optional :: t_tile(:, :)
!   Tile temperatures (n_profile, n_tile)
real(RealExt), intent(in), optional :: flux_tile(:, :, :)
!   Tile emissions (n_profile, n_tile, n_band)
real(RealExt), intent(in), optional :: albedo_diff_tile(:, :, :)
!   Diffuse tile albedo (n_profile, n_tile, n_band)
real(RealExt), intent(in), optional :: albedo_dir_tile(:, :, :)
!   Direct tile albedo (n_profile, n_tile, n_band)
real(RealExt), intent(in), optional :: frac_tile_1d(:)
!   1d tile fractions (n_tile)
real(RealExt), intent(in), optional :: t_tile_1d(:)
!   1d tile temperatures (n_tile)
real(RealExt), intent(in), optional :: flux_tile_1d(:)
!   Tile emissions (n_tile*n_band)
real(RealExt), intent(in), optional :: albedo_diff_tile_1d(:)
!   1d diffuse tile albedo (n_tile*n_band)
real(RealExt), intent(in), optional :: albedo_dir_tile_1d(:)
!   1d direct tile albedo (n_tile*n_band)

real(RealExt), intent(in), dimension (:, :), optional :: &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  liq_rsd, ice_rsd, liq_conv_rsd, ice_conv_rsd, &
  liq_nc, ice_nc, liq_conv_nc, ice_conv_nc
real(RealExt), intent(in), dimension (:), optional :: &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d, &
  liq_rsd_1d, ice_rsd_1d, liq_conv_rsd_1d, ice_conv_rsd_1d, &
  liq_nc_1d, ice_nc_1d, liq_conv_nc_1d, ice_conv_nc_1d
real(RealExt), intent(in), optional :: &
  liq_dim_constant, liq_conv_dim_constant
!   Liquid and ice cloud fractions, gridbox mean mixing ratios,
!   effective dimensions, relative standard deviation of condensate,
!   and number concentration

real(RealExt), intent(in), optional :: cloud_vertical_decorr
!   Decorrelation pressure scale for cloud vertical overlap
real(RealExt), intent(in), optional :: conv_vertical_decorr
!   Decorrelation pressure scale for convective cloud vertical overlap
real(RealExt), intent(in), optional :: cloud_horizontal_rsd
!   Relative standard deviation of sub-grid cloud condensate
real(RealExt), intent(in), optional :: liq_dim_aparam
real(RealExt), intent(in), optional :: liq_dim_bparam
!   Parameters for calculating cloud droplet effective radius

real(RealExt), intent(in), optional :: layer_heat_capacity(:, :)
real(RealExt), intent(in), optional :: layer_heat_capacity_1d(:)
!   Heat capacity of layer

integer, intent(in), optional :: i_source
!   Select source of radiation
integer, intent(in), optional :: i_scatter_method
!   Select scattering parametrisation
integer, intent(in), optional :: i_gas_overlap
!   Select treatment of overlapping gas absorption
integer, intent(in), optional :: &
  i_cloud_representation, i_overlap, i_inhom, i_cloud_entrapment, &
  i_mcica_sampling, i_st_water, i_st_ice, i_cnv_water, i_cnv_ice, i_drop_re
!   Select treatment of cloud
integer, intent(in), optional :: rand_seed(:)
!   Random seed for cloud generator

logical, intent(in), optional :: l_rayleigh
!   Include Rayleigh scattering

logical, intent(in), optional :: l_mixing_ratio
!   Assume mass mixing ratios are with respect to dry mass

logical, intent(in), optional :: l_aerosol_mode
!   Include aerosol optical properties specified by mode

real(RealExt), intent(in), optional :: aer_mix_ratio(:, :, :)
!   MODE aerosol mass-mixing ratio (n_profile, n_layer, n_mode)

real(RealExt), intent(in), optional :: aer_mix_ratio_1d(:)
!   1d MODE aerosol mass-mixing ratio (n_aer_layer*n_mode)

real(RealExt), intent(in), dimension(:, :, :, :), optional :: &
  aer_absorption, aer_scattering, aer_asymmetry
!   MODE aerosol optical properties (n_profile, n_layer, n_mode, n_band)

real(RealExt), intent(in), dimension(:), optional :: &
  aer_absorption_1d, aer_scattering_1d, aer_asymmetry_1d
!   1d MODE aerosol optical properties (n_aer_layer*n_mode*n_band)

real(RealExt), intent(in), optional :: mean_rel_humidity(:, :)
real(RealExt), intent(in), optional :: mean_rel_humidity_1d(:)
!   Mean relative humidity applicable for CLASSIC aerosols (clear-sky)

logical, intent(in), optional :: &
  l_water_soluble, l_dust_like, l_oceanic, l_soot, l_ash, l_sulphuric, &
  l_ammonium_sulphate, l_saharan_dust, &
  l_accum_sulphate, l_aitken_sulphate, &
  l_fresh_soot, l_aged_soot, &
  l_sodium_chloride, l_seasalt_film, l_seasalt_jet, &
  l_dust_div1, l_dust_div2, l_dust_div3, &
  l_dust_div4, l_dust_div5, l_dust_div6, &
  l_biomass_1, l_biomass_2, &
  l_biogenic, &
  l_ocff_fresh, l_ocff_aged, &
  l_delta, l_murk, &
  l_nitrate, &
  l_twobindust_1, l_twobindust_2
! Flags to include CLASSIC aerosols

real(RealExt), intent(in), dimension(:, :), optional :: &
  water_soluble, dust_like, oceanic, soot, ash, sulphuric, &
  ammonium_sulphate, saharan_dust, &
  accum_sulphate, aitken_sulphate, &
  fresh_soot, aged_soot, &
  sodium_chloride, seasalt_film, seasalt_jet, &
  dust_div1, dust_div2, dust_div3, &
  dust_div4, dust_div5, dust_div6, &
  biomass_1, biomass_2, &
  biogenic, &
  ocff_fresh, ocff_aged, &
  delta, murk, &
  nitrate, &
  twobindust_1, twobindust_2
! CLASSIC aerosol mass mixing ratios

real(RealExt), intent(in), dimension(:), optional :: &
  water_soluble_1d, dust_like_1d, oceanic_1d, soot_1d, ash_1d, sulphuric_1d, &
  ammonium_sulphate_1d, saharan_dust_1d, &
  accum_sulphate_1d, aitken_sulphate_1d, &
  fresh_soot_1d, aged_soot_1d, &
  sodium_chloride_1d, seasalt_film_1d, seasalt_jet_1d, &
  dust_div1_1d, dust_div2_1d, dust_div3_1d, &
  dust_div4_1d, dust_div5_1d, dust_div6_1d, &
  biomass_1_1d, biomass_2_1d, &
  biogenic_1d, &
  ocff_fresh_1d, ocff_aged_1d, &
  delta_1d, murk_1d, &
  nitrate_1d, &
  twobindust_1_1d, twobindust_2_1d
! 1d CLASSIC aerosol mass mixing ratios

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields and diagnostics

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information


! Spectral data:
type(StrSpecData), pointer :: spec

! Mcica data:
type(StrMcica), pointer :: mcica
type(StrMcica), target :: mcica_dummy

! Controlling options:
type(StrCtrl) :: control

! Dimensions:
type(StrDim) :: dimen

! Atmospheric properties:
type(StrAtm) :: atm

! Boundary conditions:
type(StrBound) :: bound

! Cloud properties:
type(StrCld) :: cld

! Aerosol properties:
type(StrAer) :: aer

! Output fields from core radiation code:
type(StrOut) :: radout

! Controlling options for diagnostic calls:
type (StrCtrl) :: control_clean

! Output fields from diagnostic calls:
type(StrOut) :: radout_clean

integer :: id_spec, id_mcica
!   Loop variables

integer :: ierr = i_normal
character (len=errormessagelength) :: cmessage
character (len=*), parameter :: RoutineName = 'RUNES'

nullify(spec, mcica)

if (present(spectrum_name)) then
  do id_spec=1, size(spectrum_array)
    if (spectrum_array_name(id_spec) == spectrum_name) exit
    if (id_spec == size(spectrum_array)) then
      cmessage = 'Spectrum name not found.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  end do
  spec => spectrum_array(id_spec)
  mcica => mcica_dummy
  if (present(i_cloud_representation).and.present(i_inhom)) then
    if ( (i_cloud_representation /= ip_cloud_representation_off) .and. &
         (i_inhom == ip_inhom_mcica) ) then
      if (allocated(mcica_data_array).and.allocated(mcica_spectrum_name)) then
        do id_mcica=1, size(mcica_data_array)
          if (mcica_spectrum_name(i_source, id_mcica) == spectrum_name) exit
          if (id_mcica == size(mcica_data_array)) then
            cmessage = 'Spectrum name not associated with MCICA data.'
            ierr=i_err_fatal
            call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
          end if
        end do
        mcica => mcica_data_array(id_mcica)
      else
        cmessage = 'MCICA data has not been read in correctly.'
        ierr=i_err_fatal
        call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      end if
    end if
  end if
else if (present(spectrum)) then
  spec => spectrum
  mcica => mcica_dummy
  if (present(i_cloud_representation).and.present(i_inhom)) then
    if ( (i_cloud_representation /= ip_cloud_representation_off) .and. &
         (i_inhom == ip_inhom_mcica) ) then
      if (present(mcica_data)) then
        mcica => mcica_data
      else
        cmessage = 'No mcica_data object has been passed to socrates_runes.'
        ierr=i_err_fatal
        call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      end if
    end if
  end if
else
  cmessage = 'No spectrum name or object supplied.'
  ierr=i_err_fatal
  call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
end if

call set_control(control, diag, spec, &
  isolir                 = i_source, &
  l_rayleigh             = l_rayleigh, &
  l_orog                 = l_orog, &
  l_mixing_ratio         = l_mixing_ratio, &
  l_aerosol_mode         = l_aerosol_mode, &
  l_tile                 = l_tile, &
  l_flux_ground          = l_flux_ground, &
  l_flux_tile            = l_flux_tile, &
  n_tile                 = n_tile, &
  n_cloud_layer          = n_cloud_layer, &
  n_aer_mode             = n_aer_mode, &
  i_scatter_method       = i_scatter_method, &
  i_gas_overlap          = i_gas_overlap, &
  i_cloud_representation = i_cloud_representation, &
  i_overlap              = i_overlap, &
  i_inhom                = i_inhom, &
  i_cloud_entrapment     = i_cloud_entrapment, &
  i_mcica_sampling       = i_mcica_sampling, &
  i_st_water             = i_st_water, &
  i_cnv_water            = i_cnv_water, &
  i_st_ice               = i_st_ice, &
  i_cnv_ice              = i_cnv_ice, &
  i_drop_re              = i_drop_re, &
  l_set_defaults         = .true.)

call set_dimen(dimen, control, n_profile, n_layer, &
  mcica_data    = mcica, &
  n_tile        = n_tile, &
  n_cloud_layer = n_cloud_layer, &
  n_aer_mode    = n_aer_mode, &
  n_subcol_gen  = n_subcol_gen )

call set_atm(atm, control, dimen, spec, n_profile, n_layer, &
  profile_list, n_layer_stride, n_level_stride, &
  p_layer, t_layer, mass, density, p_level, t_level, r_layer, r_level, &
  p_layer_1d, t_layer_1d, mass_1d, density_1d, p_level_1d, t_level_1d, &
  r_layer_1d, r_level_1d, &
  h2o, co2, o3, n2o, co, ch4, o2, no, so2, no2, nh3, hno3, n2, &
  cfc11, cfc12, cfc113, hcfc22, hfc125, hfc134a, cfc114, tio, vo, h2, &
  he, ocs, na, k, feh, crh, li, rb, cs, ph3, c2h2, hcn, h2s, ar, o, n, no3, &
  n2o5, hono, ho2no2, h2o2, c2h6, ch3, h2co, ho2, &
  hdo, hcl, hf, cosso, tosso, yosos, &
  ch3cho, ch3ooh, ch3coch3, ch3cocho, chocho, c2h5cho, &
  hoch2cho, c2h5coch3, mvk, macr, pan, ch3ono2, &
  h2o_1d, co2_1d, o3_1d, n2o_1d, co_1d, ch4_1d, o2_1d, no_1d, so2_1d, no2_1d, &
  nh3_1d, hno3_1d, n2_1d, cfc11_1d, cfc12_1d, cfc113_1d, hcfc22_1d, hfc125_1d, &
  hfc134a_1d, cfc114_1d, tio_1d, vo_1d, h2_1d, he_1d, ocs_1d, na_1d, k_1d, &
  feh_1d, crh_1d, li_1d, rb_1d, cs_1d, ph3_1d, c2h2_1d, hcn_1d, h2s_1d, ar_1d, &
  o_1d, n_1d, no3_1d, n2o5_1d, hono_1d, ho2no2_1d, h2o2_1d, c2h6_1d, ch3_1d, &
  h2co_1d, ho2_1d, hdo_1d, hcl_1d, hf_1d, cosso_1d, tosso_1d, yosos_1d, &
  ch3cho_1d, ch3ooh_1d, ch3coch3_1d, ch3cocho_1d, chocho_1d, c2h5cho_1d, &
  hoch2cho_1d, c2h5coch3_1d, mvk_1d, macr_1d, pan_1d, ch3ono2_1d, &
  h2o_mix_ratio, co2_mix_ratio, o3_mix_ratio, n2o_mix_ratio, co_mix_ratio, &
  ch4_mix_ratio, o2_mix_ratio, no_mix_ratio, so2_mix_ratio, no2_mix_ratio, &
  nh3_mix_ratio, hno3_mix_ratio, n2_mix_ratio, cfc11_mix_ratio, &
  cfc12_mix_ratio, cfc113_mix_ratio, hcfc22_mix_ratio, hfc125_mix_ratio, &
  hfc134a_mix_ratio, cfc114_mix_ratio, tio_mix_ratio, vo_mix_ratio, &
  h2_mix_ratio, he_mix_ratio, ocs_mix_ratio, na_mix_ratio, k_mix_ratio, &
  feh_mix_ratio, crh_mix_ratio, li_mix_ratio, rb_mix_ratio, cs_mix_ratio, &
  ph3_mix_ratio, c2h2_mix_ratio, hcn_mix_ratio, h2s_mix_ratio, ar_mix_ratio, &
  o_mix_ratio, n_mix_ratio, no3_mix_ratio, n2o5_mix_ratio, hono_mix_ratio, &
  ho2no2_mix_ratio, h2o2_mix_ratio, c2h6_mix_ratio, ch3_mix_ratio, &
  h2co_mix_ratio, ho2_mix_ratio, hdo_mix_ratio, hcl_mix_ratio, hf_mix_ratio, &
  cosso_mix_ratio, tosso_mix_ratio, yosos_mix_ratio, ch3cho_mix_ratio, &
  ch3ooh_mix_ratio, ch3coch3_mix_ratio, ch3cocho_mix_ratio, chocho_mix_ratio, &
  c2h5cho_mix_ratio, hoch2cho_mix_ratio, c2h5coch3_mix_ratio, mvk_mix_ratio, &
  macr_mix_ratio, pan_mix_ratio, ch3ono2_mix_ratio, &
  l_h2o_well_mixed, l_co2_well_mixed, l_o3_well_mixed, l_n2o_well_mixed, &
  l_co_well_mixed, l_ch4_well_mixed, l_o2_well_mixed, l_no_well_mixed, &
  l_so2_well_mixed, l_no2_well_mixed, l_nh3_well_mixed, l_hno3_well_mixed, &
  l_n2_well_mixed, l_cfc11_well_mixed, l_cfc12_well_mixed, &
  l_cfc113_well_mixed, l_hcfc22_well_mixed, l_hfc125_well_mixed, &
  l_hfc134a_well_mixed, l_cfc114_well_mixed, l_tio_well_mixed, &
  l_vo_well_mixed, l_h2_well_mixed, l_he_well_mixed, l_ocs_well_mixed, &
  l_na_well_mixed, l_k_well_mixed, l_feh_well_mixed, l_crh_well_mixed, &
  l_li_well_mixed, l_rb_well_mixed, l_cs_well_mixed, l_ph3_well_mixed, &
  l_c2h2_well_mixed, l_hcn_well_mixed, l_h2s_well_mixed, l_ar_well_mixed, &
  l_o_well_mixed, l_n_well_mixed, l_no3_well_mixed, l_n2o5_well_mixed, &
  l_hono_well_mixed, l_ho2no2_well_mixed, l_h2o2_well_mixed, &
  l_c2h6_well_mixed, l_ch3_well_mixed, l_h2co_well_mixed, l_ho2_well_mixed, &
  l_hdo_well_mixed, l_hcl_well_mixed, l_hf_well_mixed, l_cosso_well_mixed, &
  l_tosso_well_mixed, l_yosos_well_mixed, l_ch3cho_well_mixed, &
  l_ch3ooh_well_mixed, l_ch3coch3_well_mixed, l_ch3cocho_well_mixed, &
  l_chocho_well_mixed, l_c2h5cho_well_mixed, l_hoch2cho_well_mixed, &
  l_c2h5coch3_well_mixed, l_mvk_well_mixed, l_macr_well_mixed, &
  l_pan_well_mixed, l_ch3ono2_well_mixed, &
  l_invert, l_profile_last, l_debug, i_profile_debug)

call set_bound(bound, control, dimen, spec, n_profile, &
  profile_list        = profile_list, &
  n_tile              = n_tile, &
  t_ground            = t_ground, &
  flux_ground         = flux_ground, &
  flux_ground_1d      = flux_ground_1d, &
  cos_zenith_angle    = cos_zenith_angle, &
  solar_irrad         = solar_irrad, &
  orog_corr           = orog_corr, &
  l_grey_albedo       = l_grey_albedo, &
  grey_albedo         = grey_albedo, &
  albedo_diff         = albedo_diff, &
  albedo_dir          = albedo_dir, &
  albedo_diff_1d      = albedo_diff_1d, &
  albedo_dir_1d       = albedo_dir_1d, &
  frac_tile           = frac_tile, &
  t_tile              = t_tile, &
  flux_tile           = flux_tile, &
  albedo_diff_tile    = albedo_diff_tile, &
  albedo_dir_tile     = albedo_dir_tile, &
  frac_tile_1d        = frac_tile_1d, &
  t_tile_1d           = t_tile_1d, &
  flux_tile_1d        = flux_tile_1d, &
  albedo_diff_tile_1d = albedo_diff_tile_1d, &
  albedo_dir_tile_1d  = albedo_dir_tile_1d, &
  l_profile_last      = l_profile_last, &
  l_debug             = l_debug, &
  i_profile_debug     = i_profile_debug )

call set_cld(cld, control, dimen, spec, atm, &
  profile_list          = profile_list, &
  n_layer_stride        = n_layer_stride, &
  cloud_frac            = cloud_frac, &
  conv_frac             = conv_frac, &
  liq_frac              = liq_frac, &
  ice_frac              = ice_frac, &
  liq_conv_frac         = liq_conv_frac, &
  ice_conv_frac         = ice_conv_frac, &
  liq_mmr               = liq_mmr, &
  ice_mmr               = ice_mmr, &
  liq_conv_mmr          = liq_conv_mmr, &
  ice_conv_mmr          = ice_conv_mmr, &
  liq_rsd               = liq_rsd, &
  ice_rsd               = ice_rsd, &
  liq_conv_rsd          = liq_conv_rsd, &
  ice_conv_rsd          = ice_conv_rsd, &
  cloud_frac_1d         = cloud_frac_1d, &
  conv_frac_1d          = conv_frac_1d, &
  liq_frac_1d           = liq_frac_1d, &
  ice_frac_1d           = ice_frac_1d, &
  liq_conv_frac_1d      = liq_conv_frac_1d, &
  ice_conv_frac_1d      = ice_conv_frac_1d, &
  liq_mmr_1d            = liq_mmr_1d, &
  ice_mmr_1d            = ice_mmr_1d, &
  liq_conv_mmr_1d       = liq_conv_mmr_1d, &
  ice_conv_mmr_1d       = ice_conv_mmr_1d, &
  liq_rsd_1d            = liq_rsd_1d, &
  ice_rsd_1d            = ice_rsd_1d, &
  liq_conv_rsd_1d       = liq_conv_rsd_1d, &
  ice_conv_rsd_1d       = ice_conv_rsd_1d, &
  cloud_vertical_decorr = cloud_vertical_decorr, &
  conv_vertical_decorr  = conv_vertical_decorr, &
  cloud_horizontal_rsd  = cloud_horizontal_rsd, &
  l_invert              = l_invert, &
  l_profile_last        = l_profile_last, &
  l_debug               = l_debug, &
  i_profile_debug       = i_profile_debug )

call set_cld_dim(cld, control, dimen, spec, atm, &
  profile_list          = profile_list, &
  n_layer_stride        = n_layer_stride, &
  liq_nc                = liq_nc, &
  ice_nc                = ice_nc, &
  liq_conv_nc           = liq_conv_nc, &
  ice_conv_nc           = ice_conv_nc, &
  liq_dim               = liq_dim, &
  ice_dim               = ice_dim, &
  liq_conv_dim          = liq_conv_dim, &
  ice_conv_dim          = ice_conv_dim, &
  liq_nc_1d             = liq_nc_1d, &
  ice_nc_1d             = ice_nc_1d, &
  liq_conv_nc_1d        = liq_conv_nc_1d, &
  ice_conv_nc_1d        = ice_conv_nc_1d, &
  liq_dim_1d            = liq_dim_1d, &
  ice_dim_1d            = ice_dim_1d, &
  liq_conv_dim_1d       = liq_conv_dim_1d, &
  ice_conv_dim_1d       = ice_conv_dim_1d, &
  liq_dim_constant      = liq_dim_constant, &
  liq_conv_dim_constant = liq_conv_dim_constant, &
  liq_dim_aparam        = liq_dim_aparam, &
  liq_dim_bparam        = liq_dim_bparam, &
  l_invert              = l_invert, &
  l_profile_last        = l_profile_last, &
  l_debug               = l_debug, &
  i_profile_debug       = i_profile_debug )

call set_cld_mcica(cld, mcica, control, dimen, spec, atm, &
  profile_list = profile_list, &
  rand_seed    = rand_seed )

call set_aer(aer, control, dimen, spec, &
  n_profile, n_layer, n_aer_mode, profile_list, n_layer_stride, n_aer_layer, &
  aer_mix_ratio, aer_absorption, aer_scattering, aer_asymmetry, &
  aer_mix_ratio_1d, aer_absorption_1d, aer_scattering_1d, aer_asymmetry_1d, &
  mean_rel_humidity, mean_rel_humidity_1d, &
  l_water_soluble, water_soluble, water_soluble_1d, &
  l_dust_like, dust_like, dust_like_1d, &
  l_oceanic, oceanic, oceanic_1d, &
  l_soot, soot, soot_1d, &
  l_ash, ash, ash_1d, &
  l_sulphuric, sulphuric, sulphuric_1d, &
  l_ammonium_sulphate, ammonium_sulphate, ammonium_sulphate_1d, &
  l_saharan_dust, saharan_dust, saharan_dust_1d, &
  l_accum_sulphate, accum_sulphate, accum_sulphate_1d, &
  l_aitken_sulphate, aitken_sulphate, aitken_sulphate_1d, &
  l_fresh_soot, fresh_soot, fresh_soot_1d, &
  l_aged_soot, aged_soot, aged_soot_1d, &
  l_sodium_chloride, sodium_chloride, sodium_chloride_1d, &
  l_seasalt_film, seasalt_film, seasalt_film_1d, &
  l_seasalt_jet, seasalt_jet, seasalt_jet_1d, &
  l_dust_div1, dust_div1, dust_div1_1d, &
  l_dust_div2, dust_div2, dust_div2_1d, &
  l_dust_div3, dust_div3, dust_div3_1d, &
  l_dust_div4, dust_div4, dust_div4_1d, &
  l_dust_div5, dust_div5, dust_div5_1d, &
  l_dust_div6, dust_div6, dust_div6_1d, &
  l_biomass_1, biomass_1, biomass_1_1d, &
  l_biomass_2, biomass_2, biomass_2_1d, &
  l_biogenic, biogenic, biogenic_1d, &
  l_ocff_fresh, ocff_fresh, ocff_fresh_1d, &
  l_ocff_aged, ocff_aged, ocff_aged_1d, &
  l_delta, delta, delta_1d, &
  l_murk, murk, murk_1d, &
  l_nitrate, nitrate, nitrate_1d, &
  l_twobindust_1, twobindust_1, twobindust_1_1d, &
  l_twobindust_2, twobindust_2, twobindust_2_1d, &
  l_invert, l_profile_last)

! DEPENDS ON: radiance_calc
if (control%l_clean .or. control%l_clear_clean) then
  ! Calculate radiative transfer for clean-air diagnostics
  control_clean = control
  control_clean%l_aerosol = .false.
  control_clean%l_aerosol_mode = .false.
  control_clean%l_clear = control%l_clear_clean

  control_clean%l_flux_direct_band = associated(diag%flux_direct_clean_band)
  control_clean%l_flux_down_band = associated(diag%flux_down_clean_band)
  control_clean%l_flux_up_band = associated(diag%flux_up_clean_band)
  control_clean%l_flux_direct_clear_band = associated( &
             diag%flux_direct_clear_clean_band)
  control_clean%l_flux_down_clear_band = associated( &
             diag%flux_down_clear_clean_band)
  control_clean%l_flux_up_clear_band = associated( &
             diag%flux_up_clear_clean_band)
  control_clean%l_photolysis_rate = control%l_photolysis_rate_clean
  control_clean%l_photolysis_rate_clear = control%l_photolysis_rate_clear_clean

  control_clean%l_aerosol_absorption_band = .false.
  control_clean%l_aerosol_scattering_band = .false.
  control_clean%l_aerosol_asymmetry_band  = .false.
  control_clean%l_cloud_extinction        = .false.
  control_clean%l_cloud_absorptivity      = .false.
  call radiance_calc(control_clean, dimen, spec, atm, cld, aer, bound, &
                     radout_clean)
end if

if (associated(diag%heating_rate) .or. &
    associated(diag%flux_direct) .or. &
    associated(diag%flux_down) .or. &
    associated(diag%flux_up) .or. &
    associated(diag%flux_direct_surf) .or. &
    associated(diag%flux_down_surf) .or. &
    associated(diag%flux_up_surf) .or. &
    associated(diag%flux_direct_toa) .or. &
    associated(diag%flux_down_toa) .or. &
    associated(diag%flux_up_toa) .or. &
    associated(diag%flux_up_tile) .or. &
    control%l_clear .or. &
    control%l_blue_flux_surf .or. &
    control%l_aerosol_absorption_band .or. &
    control%l_aerosol_scattering_band .or. &
    control%l_aerosol_asymmetry_band .or. &
    control%l_cloud_absorptivity .or. &
    control%l_cloud_extinction) then
  ! Calculate radiative transfer
  call radiance_calc(control, dimen, spec, atm, cld, aer, bound, radout)
end if

call set_diag(diag, &
  control, dimen, spec, atm, cld, mcica, aer, bound, radout, radout_clean, &
  n_profile, n_layer, &
  profile_list           = profile_list, &
  n_layer_stride         = n_layer_stride, &
  n_tile                 = n_tile, &
  layer_heat_capacity    = layer_heat_capacity, &
  layer_heat_capacity_1d = layer_heat_capacity_1d, &
  l_invert               = l_invert, &
  l_profile_last         = l_profile_last)

call deallocate_out(radout)
call deallocate_aer_prsc(aer)
call deallocate_aer(aer)
call deallocate_cld_mcica(cld)
call deallocate_cld_prsc(cld)
call deallocate_cld(cld)
call deallocate_bound(bound)
call deallocate_atm(atm)
call deallocate_control(control)

end subroutine runes
end module socrates_runes
