! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Define the Socrates diagnostic output type

module socrates_def_diag

use realtype_rd, only: RealExt
use gas_list_pcf, only: npd_pathway

implicit none

type :: StrPhotol
  real(RealExt), pointer :: rate(:,:) => null()
  ! Photolysis rate, s-1 (n_profile, n_layer)
  real(RealExt), pointer :: rate_clear(:,:) => null()
  ! Clear-sky photolysis rate, s-1 (n_profile, n_layer)
  real(RealExt), pointer :: rate_clean(:,:) => null()
  ! Clean-air photolysis rate, s-1 (n_profile, n_layer)
  real(RealExt), pointer :: rate_clear_clean(:,:) => null()
  ! Clear-clean photolysis rate, s-1 (n_profile, n_layer)
end type StrPhotol

type :: StrDiag

type(StrPhotol) :: photolysis(npd_pathway)
! Photolysis diagnostics

real(RealExt), pointer :: photolysis_rate(:,:,:) => null()
! Photolysis rate, s-1 (n_profile, n_layer, n_pathway)

real(RealExt), pointer :: heating_rate(:,:) => null()
! Heating rate, Ks-1 (n_profile, n_layer)

real(RealExt), pointer :: flux_direct(:,:) => null()
! Direct (unscattered) downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_down(:,:) => null()
! Downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up(:,:) => null()
! Upwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_direct_clear(:,:) => null()
! Clear-sky direct (unscattered) downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_down_clear(:,:) => null()
! Clear-sky downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up_clear(:,:) => null()
! Clear-sky upwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_direct_clean(:,:) => null()
! Clean-air (aerosol free) direct downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_down_clean(:,:) => null()
! Clean-air downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up_clean(:,:) => null()
! Clean-air upwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_direct_clear_clean(:,:) => null()
! Clear-clean direct downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_down_clear_clean(:,:) => null()
! Clear-clean downwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up_clear_clean(:,:) => null()
! Clear-clean upwards flux, Wm-2 (n_profile, 0:n_layer)

real(RealExt), pointer :: flux_up_tile(:,:) => null()
! Upwards flux on tiles, Wm-2 (n_profile, n_tile)

real(RealExt), pointer :: flux_up_blue_tile(:,:) => null()
! Upwards blue flux on tiles, Wm-2 (n_profile, n_tile)

real(RealExt), pointer :: flux_direct_surf(:) => null()
! Direct (unscattered) downwards flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_surf(:) => null()
! Downwards flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_surf(:) => null()
! Upwards flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_clear_surf(:) => null()
! Clear-sky direct downwards surface flux, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_clear_surf(:) => null()
! Clear-sky downwards surface flux, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_clear_surf(:) => null()
! Clear-sky upwards surface flux, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_clean_surf(:) => null()
! Clean-air direct flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_clean_surf(:) => null()
! Clean-air downwards flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_clean_surf(:) => null()
! Clean-air upwards flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_clear_clean_surf(:) => null()
! Clear-clean direct surface flux, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_clear_clean_surf(:) => null()
! Clear-clean downwards surface flux, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_clear_clean_surf(:) => null()
! Clear-clean upwards surface flux, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_blue_surf(:) => null()
! Direct blue flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_blue_surf(:) => null()
! Total downward blue flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_uv_surf(:) => null()
! Direct UV flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_uv_surf(:) => null()
! Downwards UV flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_uv_surf(:) => null()
! Upwards UV flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_uv_clear_surf(:) => null()
! Clear-sky direct UV flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_uv_clear_surf(:) => null()
! Clear-sky downwards UV flux at the surface, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_uv_clear_surf(:) => null()
! Clear-sky upwards UV flux at the surface, Wm-2 (n_profile)

real(RealExt) :: wavelength_uv = 3.2e-07_RealExt
! Upper wavelength boundary for the UV flux diagnostics, defaults to 320nm

real(RealExt), pointer :: flux_direct_toa(:) => null()
! Direct downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_toa(:) => null()
! Downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_toa(:) => null()
! Upwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_clear_toa(:) => null()
! Clear-sky direct downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_clear_toa(:) => null()
! Clear-sky downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_clear_toa(:) => null()
! Clear-sky upwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_clean_toa(:) => null()
! Clean-air direct downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_clean_toa(:) => null()
! Clean-air downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_clean_toa(:) => null()
! Clean-air upwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_clear_clean_toa(:) => null()
! Clear-clean direct downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_down_clear_clean_toa(:) => null()
! Clear-clean downwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_up_clear_clean_toa(:) => null()
! Clear-clean upwards flux at top-of-atmosphere, Wm-2 (n_profile)

real(RealExt), pointer :: flux_direct_band(:,:,:) => null()
! Direct downwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_down_band(:,:,:) => null()
! Downwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_up_band(:,:,:) => null()
! Upwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_direct_clear_band(:,:,:) => null()
! Clear-sky direct flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_down_clear_band(:,:,:) => null()
! Clear-sky downwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_up_clear_band(:,:,:) => null()
! Clear-sky upwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_direct_clean_band(:,:,:) => null()
! Clean-air direct flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_down_clean_band(:,:,:) => null()
! Clean-air downwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_up_clean_band(:,:,:) => null()
! Clean-air upwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_direct_clear_clean_band(:,:,:) => null()
! Clear-clean direct flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_down_clear_clean_band(:,:,:) => null()
! Clear-clean downwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: flux_up_clear_clean_band(:,:,:) => null()
! Clear-clean upwards flux per band, Wm-2 (n_profile, 0:n_layer, n_band)

real(RealExt), pointer :: total_cloud_cover(:) => null()
! Total cloud cover (n_profile)

real(RealExt), pointer :: total_cloud_fraction(:,:) => null()
! Total cloud fraction in layers (n_profile, n_layer)

real(RealExt), pointer :: liq_frac(:,:) => null()
! Liquid cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: liq_conv_frac(:,:) => null()
! Liquid convective cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: liq_part_frac(:,:) => null()
! Liquid cloud fraction of the total cloud (n_profile, n_layer)

real(RealExt), pointer :: liq_conv_part_frac(:,:) => null()
! Liquid convective cloud fraction of the total cloud (n_profile, n_layer)

real(RealExt), pointer :: liq_incloud_mmr(:,:) => null()
! Liquid in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: liq_mmr(:,:) => null()
! Liquid gridbox mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: liq_inconv_mmr(:,:) => null()
! Liquid convective in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: liq_conv_mmr(:,:) => null()
! Liquid convective gridbox mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: liq_dim(:,:) => null()
! Cloud droplet effective dimension (n_profile, n_layer)

real(RealExt), pointer :: liq_conv_dim(:,:) => null()
! Convective cloud droplet effective dimension (n_profile, n_layer)

real(RealExt), pointer :: liq_path(:) => null()
! Liquid mass path in grid column (n_profile)

real(RealExt), pointer :: liq_conv_path(:) => null()
! Convective liquid mass path in grid column (n_profile)

real(RealExt), pointer :: ice_frac(:,:) => null()
! Ice cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_frac(:,:) => null()
! Ice convective cloud fraction (n_profile, n_layer)

real(RealExt), pointer :: ice_part_frac(:,:) => null()
! Ice cloud fraction of the total cloud (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_part_frac(:,:) => null()
! Ice convective cloud fraction of the total cloud (n_profile, n_layer)

real(RealExt), pointer :: ice_incloud_mmr(:,:) => null()
! Ice in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: ice_mmr(:,:) => null()
! Ice gridbox mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: ice_inconv_mmr(:,:) => null()
! Ice convective in-cloud mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_mmr(:,:) => null()
! Ice convective gridbox mean mixing ratio (n_profile, n_layer)

real(RealExt), pointer :: ice_dim(:,:) => null()
! Cloud ice-crystal effective dimension (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_dim(:,:) => null()
! Convective cloud ice-crystal effective dimension (n_profile, n_layer)

real(RealExt), pointer :: ice_re(:,:) => null()
! Cloud ice-crystal effective radius (n_profile, n_layer)

real(RealExt), pointer :: ice_conv_re(:,:) => null()
! Convective cloud ice-crystal effective radius (n_profile, n_layer)

real(RealExt), pointer :: ice_path(:) => null()
! Ice mass path in grid column (n_profile)

real(RealExt), pointer :: ice_conv_path(:) => null()
! Convective ice mass path in grid column (n_profile)

real(RealExt), pointer :: cloud_top_liq_dim(:) => null()
! Cloud droplet effective radius at cloud top weighted by cloud fraction

real(RealExt), pointer :: cloud_top_liq_weight(:) => null()
! Weight for liquid cloud fraction at cloud top

real(RealExt), pointer :: cloud_top_warm_liq_dim(:) => null()
! Warm cloud droplet effective radius at cloud top weighted by cloud fraction

real(RealExt), pointer :: cloud_top_warm_liq_weight(:) => null()
! Weight for warm liquid cloud fraction at cloud top

integer, pointer :: n_subcol_cloud(:) => null()
! Number of generated sub-columns containing cloud

real(RealExt), pointer :: liq_subcol_scaling(:,:,:) => null()
! Scaling factor for liquid condensate in each sub-column
! (n_profile, n_layer, n_subcol_gen)

real(RealExt), pointer :: ice_subcol_scaling(:,:,:) => null()
! Scaling factor for ice condensate in each sub-column
! (n_profile, n_layer, n_subcol_gen)

real(RealExt), pointer :: cloud_absorptivity(:, :) => null()
! Absorptivity of cloud weighted by cloud fraction
! and upward clear-sky infra-red flux

real(RealExt), pointer :: cloud_weight_absorptivity(:, :) => null()
! Weight for cloud_absorptivity:
! cloud fraction * upward clear-sky infra-red flux

real(RealExt), pointer :: cloud_extinction(:, :) => null()
! Cloud extinction weighted by cloud fraction
! and downward clear-sky solar flux

real(RealExt), pointer :: cloud_weight_extinction(:, :) => null()
! Weight for cloud_extinction:
! cloud fraction * downward clear-sky solar flux

real(RealExt), pointer :: cloud_thermal_absorptivity(:, :) => null()
! Absorptivity of cloud averaged using the Planckian flux
! in each band for the local temperature, or at a particular
! wavelength if a value > 0 is specified:
real(RealExt) :: cloud_absorptivity_wavelength = -1.0_RealExt

real(RealExt), pointer :: cloud_solar_extinction(:, :) => null()
! Cloud extinction averaged using solar spectrum in each band

real(RealExt), pointer :: aerosol_optical_depth(:,:,:) => null()
! Total aerosol optical depth (n_profile, n_layer, n_band)

real(RealExt), pointer :: aerosol_scat_optical_depth(:,:,:) => null()
! Total aerosol scattering optical depth (n_profile, n_layer, n_band)

real(RealExt), pointer :: aerosol_asymmetry_scat(:,:,:) => null()
! Total aerosol asymmetry weighted by scattering optical depth

end type StrDiag

end module socrates_def_diag
