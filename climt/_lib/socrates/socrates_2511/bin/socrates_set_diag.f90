! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the diagnostics to be output from the Socrates runes interface
!
!------------------------------------------------------------------------------
module socrates_set_diag
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_DIAG'
contains

subroutine set_diag(diag, control, dimen, spectrum, &
  atm, cld, mcica_data, aer, bound, radout, radout_clean, &
  n_profile, n_layer, profile_list, n_layer_stride, n_tile, &
  layer_heat_capacity, layer_heat_capacity_1d, &
  l_invert, l_profile_last)

use socrates_def_diag, only: StrDiag

use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use def_cld,      only: StrCld
use def_mcica,    only: StrMcica
use def_aer,      only: StrAer
use def_bound,    only: StrBound
use def_out,      only: StrOut

use realtype_rd, only: RealK, RealExt
use rad_pcf, only: ip_cloud_off, ip_mcica, &
  ip_clcmp_st_water, ip_clcmp_st_ice, ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
  ip_ice_adt, ip_ice_agg_de, ip_ice_agg_de_sun, ip_ice_baran

use socrates_cloud_abs_diag, only: cloud_abs_diag
use socrates_cloud_ext_diag, only: cloud_ext_diag
use socrates_cloud_level_diag, only: cloud_level_diag

implicit none

! Output diagnostic fields
type(StrDiag),   intent(inout) :: diag

! Control options:
type(StrCtrl),      intent(in) :: control

! Dimensions:
type(StrDim),       intent(in) :: dimen

! Spectral data:
type (StrSpecData), intent(in) :: spectrum

! Atmospheric properties:
type(StrAtm),       intent(in) :: atm

! Cloud properties:
type(StrCld),       intent(in) :: cld

! MCICA data:
type(StrMcica),     intent(in) :: mcica_data

! Aerosol properties:
type(StrAer),       intent(in) :: aer

! Boundary conditions:
type(StrBound),     intent(in) :: bound

! Output fields from core radiation code:
type(StrOut),       intent(in) :: radout, radout_clean

integer, intent(in) :: n_profile
!   Number of columns to operate on
integer, intent(in) :: n_layer
!   Number of layers for radiation
integer, intent(in), optional :: profile_list(:)
!   List of profiles to fill in output fields
integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays
integer, intent(in), optional :: n_tile
!   Number of surface tiles

real(RealExt), intent(in), optional :: layer_heat_capacity(:, :)
real(RealExt), intent(in), optional :: layer_heat_capacity_1d(:)
!   Heat capacity of layer

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields and diagnostics


integer :: i, ii, l, ll, j, k, list(n_profile)
!   Loop variables
integer :: layer_offset, level_offset, stride_layer
!   Offset to loop counters to allow indexing in inverted order
logical :: l_last
!   Local flag to loop over profiles last in 1d fields and diagnostics
real(RealExt) :: flux_divergence(n_profile, n_layer)
!   Flux divergence across layer (Wm-2)


if (present(profile_list)) then
  list = profile_list(1:n_profile)
else
  do l=1, n_profile
    list(l) = l
  end do
end if

layer_offset = 0
level_offset = 0
if (present(l_invert)) then
  if (l_invert) then
    ! The layer is indexed using an inverted loop counter
    layer_offset = n_layer + 1
    level_offset = n_layer
  end if
end if

if (present(l_profile_last)) then
  l_last = l_profile_last
else
  l_last = .false.
end if

! Set the number of layers in the 1d arrays
if (present(n_layer_stride)) then
  stride_layer = n_layer_stride
else
  stride_layer = n_layer
end if


!------------------------------------------------------------------------------
! Heating rates
!------------------------------------------------------------------------------
if (associated(diag%heating_rate)) then
  if (control%l_flux_div) then
    do i=1, n_layer
      ii = abs(layer_offset-i)
      do l=1, n_profile
        flux_divergence(l, ii) = real( &
          sum(radout%flux_div(l, i, 1:control%n_channel)), RealExt)
      end do
    end do
  else
    do i=1, n_layer
      ii = abs(layer_offset-i)
      do l=1, n_profile
        flux_divergence(l, ii) = real( &
          sum(radout%flux_down(l, i-1, 1:control%n_channel)) - &
          sum(radout%flux_down(l, i,   1:control%n_channel)) + &
          sum(radout%flux_up(  l, i,   1:control%n_channel)) - &
          sum(radout%flux_up(  l, i-1, 1:control%n_channel)), RealExt)
      end do
    end do
  end if
  if (present(layer_heat_capacity)) then
    if (l_last) then
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(i, list(l)) = &
            flux_divergence(l, i) / layer_heat_capacity(i, list(l))
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(list(l), i) = &
            flux_divergence(l, i) / layer_heat_capacity(list(l), i)
        end do
      end do
    end if
  else if (present(layer_heat_capacity_1d)) then
    if (l_last) then
      do l=1, n_profile
        do i=1, n_layer
          ll = stride_layer*(list(l)-1) + i
          diag%heating_rate(i, list(l)) = &
            flux_divergence(l, i) / layer_heat_capacity_1d(ll)
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          ll = n_profile*(i-1) + list(l)
          diag%heating_rate(list(l), i) = &
            flux_divergence(l, i) / layer_heat_capacity_1d(ll)
        end do
      end do
    end if
  else
    ! Just return the flux_divergence if no heat capacity supplied
    if (l_last) then
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(i, list(l)) = flux_divergence(l, i)
        end do
      end do
    else
      do i=1, n_layer
        do l=1, n_profile
          diag%heating_rate(list(l), i) = flux_divergence(l, i)
        end do
      end do
    end if
  end if
  ! Fill extra layers with zeros
  if (l_last) then
    do l=1, n_profile
      do i=lbound(diag%heating_rate, 1), 0
        diag%heating_rate(i, list(l)) = 0.0_RealExt
      end do
      do i=n_layer+1, ubound(diag%heating_rate, 1)
        diag%heating_rate(i, list(l)) = 0.0_RealExt
      end do
    end do
  else
    do i=lbound(diag%heating_rate, 2), 0
      do l=1, n_profile
        diag%heating_rate(list(l), i) = 0.0_RealExt
      end do
    end do
    do i=n_layer+1, ubound(diag%heating_rate, 2)
      do l=1, n_profile
        diag%heating_rate(list(l), i) = 0.0_RealExt
      end do
    end do
  end if
end if


!------------------------------------------------------------------------------
! Fluxes
!------------------------------------------------------------------------------
call sum_flux_channels(diag%flux_direct, radout%flux_direct, &
                  surf=diag%flux_direct_surf, &
                   toa=diag%flux_direct_toa)
call sum_flux_channels(diag%flux_down, radout%flux_down, &
                  surf=diag%flux_down_surf, &
                   toa=diag%flux_down_toa)
call sum_flux_channels(diag%flux_up, radout%flux_up, &
                  surf=diag%flux_up_surf, &
                   toa=diag%flux_up_toa)
call sum_flux_channels(diag%flux_direct_clear, radout%flux_direct_clear, &
                  surf=diag%flux_direct_clear_surf, &
                   toa=diag%flux_direct_clear_toa)
call sum_flux_channels(diag%flux_down_clear, radout%flux_down_clear, &
                  surf=diag%flux_down_clear_surf, &
                   toa=diag%flux_down_clear_toa)
call sum_flux_channels(diag%flux_up_clear, radout%flux_up_clear, &
                  surf=diag%flux_up_clear_surf, &
                   toa=diag%flux_up_clear_toa)
call sum_flux_channels(diag%flux_direct_clean, radout_clean%flux_direct, &
                  surf=diag%flux_direct_clean_surf, &
                   toa=diag%flux_direct_clean_toa)
call sum_flux_channels(diag%flux_down_clean, radout_clean%flux_down, &
                  surf=diag%flux_down_clean_surf, &
                   toa=diag%flux_down_clean_toa)
call sum_flux_channels(diag%flux_up_clean, radout_clean%flux_up, &
                  surf=diag%flux_up_clean_surf, &
                   toa=diag%flux_up_clean_toa)
call sum_flux_channels(diag%flux_direct_clear_clean, &
               radout_clean%flux_direct_clear, &
                  surf=diag%flux_direct_clear_clean_surf, &
                   toa=diag%flux_direct_clear_clean_toa)
call sum_flux_channels(diag%flux_down_clear_clean, &
               radout_clean%flux_down_clear, &
                  surf=diag%flux_down_clear_clean_surf, &
                   toa=diag%flux_down_clear_clean_toa)
call sum_flux_channels(diag%flux_up_clear_clean, &
               radout_clean%flux_up_clear, &
                  surf=diag%flux_up_clear_clean_surf, &
                   toa=diag%flux_up_clear_clean_toa)

call sum_tile_channels(diag%flux_up_tile, radout%flux_up_tile)
call sum_tile_channels(diag%flux_up_blue_tile, radout%flux_up_blue_tile)
if (associated(diag%flux_direct_blue_surf)) then
  do l=1, n_profile
    diag%flux_direct_blue_surf(list(l)) &
      = real(radout%flux_direct_blue_surf(l), RealExt)
  end do
end if
if (associated(diag%flux_down_blue_surf)) then
  do l=1, n_profile
    diag%flux_down_blue_surf(list(l)) &
      = real(radout%flux_down_blue_surf(l), RealExt)
  end do
end if
call set_flux_bands(diag%flux_direct_band, radout%flux_direct_band)
call set_flux_bands(diag%flux_down_band, radout%flux_down_band)
call set_flux_bands(diag%flux_up_band, radout%flux_up_band)
call set_flux_bands(diag%flux_direct_clear_band, radout%flux_direct_clear_band)
call set_flux_bands(diag%flux_down_clear_band, radout%flux_down_clear_band)
call set_flux_bands(diag%flux_up_clear_band, radout%flux_up_clear_band)
call set_flux_bands(diag%flux_direct_clean_band, radout_clean%flux_direct_band)
call set_flux_bands(diag%flux_down_clean_band, radout_clean%flux_down_band)
call set_flux_bands(diag%flux_up_clean_band, radout_clean%flux_up_band)
call set_flux_bands(diag%flux_direct_clear_clean_band, &
            radout_clean%flux_direct_clear_band)
call set_flux_bands(diag%flux_down_clear_clean_band, &
            radout_clean%flux_down_clear_band)
call set_flux_bands(diag%flux_up_clear_clean_band, &
            radout_clean%flux_up_clear_band)
call sum_surf_weight_bands(diag%flux_direct_uv_surf, radout%flux_direct_band)
call sum_surf_weight_bands(diag%flux_down_uv_surf, radout%flux_down_band)
call sum_surf_weight_bands(diag%flux_up_uv_surf, radout%flux_up_band)
call sum_surf_weight_bands(diag%flux_direct_uv_clear_surf, &
                           radout%flux_direct_clear_band)
call sum_surf_weight_bands(diag%flux_down_uv_clear_surf, &
                           radout%flux_down_clear_band)
call sum_surf_weight_bands(diag%flux_up_uv_clear_surf, &
                           radout%flux_up_clear_band)

!------------------------------------------------------------------------------
! Photolysis rates
!------------------------------------------------------------------------------
do k=1, spectrum%photol%n_pathway
  call sum_photol_channels(diag%photolysis(k)%rate, &
                           radout%photolysis_rate, k)
  call sum_photol_channels(diag%photolysis(k)%rate_clear, &
                           radout%photolysis_rate_clear, k)
  call sum_photol_channels(diag%photolysis(k)%rate_clean, &
                           radout_clean%photolysis_rate, k)
  call sum_photol_channels(diag%photolysis(k)%rate_clear_clean, &
                           radout_clean%photolysis_rate_clear, k)
end do
call sum_photol_channels_all(diag%photolysis_rate, radout%photolysis_rate)

!------------------------------------------------------------------------------
! Cloud diagnostics
!------------------------------------------------------------------------------
if (associated(diag%total_cloud_cover)) then
  if (control%i_cloud_representation == ip_cloud_off) then
    do l=1, n_profile
      diag%total_cloud_cover(list(l)) = 0.0_RealExt
    end do
  else
    if (control%i_inhom == ip_mcica) then
      do l=1, n_profile
        diag%total_cloud_cover(list(l)) &
          = real(cld%n_subcol_cld(l), RealExt) &
          / real(dimen%nd_subcol_gen, RealExt)
      end do
    else
      do l=1, n_profile
        diag%total_cloud_cover(list(l)) &
          = real(radout%tot_cloud_cover(l), RealExt)
      end do
    end if
  end if
end if

if (associated(diag%n_subcol_cloud)) then
  if (control%i_cloud_representation /= ip_cloud_off .and. &
      control%i_inhom == ip_mcica) then
    do l=1, n_profile
      diag%n_subcol_cloud(list(l)) = cld%n_subcol_cld(l)
    end do
  else
    do l=1, n_profile
      diag%n_subcol_cloud(list(l)) = 0
    end do
  end if
end if

if (associated(diag%total_cloud_fraction)) then
  call set_layer_scalar(diag%total_cloud_fraction, 0.0_RealExt)
  if (control%i_cloud_representation /= ip_cloud_off) then
    call set_cloud_prop(diag%total_cloud_fraction, cld%w_cloud)
  end if
end if

! Cloud component diagnostics
call set_layer_scalar(diag%liq_dim, 0.0_RealExt)
call set_layer_scalar(diag%liq_incloud_mmr, 0.0_RealExt)
call set_layer_scalar(diag%liq_mmr, 0.0_RealExt)
call set_layer_scalar(diag%liq_frac, 0.0_RealExt)
call set_layer_scalar(diag%liq_part_frac, 0.0_RealExt)
call set_layer_scalar(diag%ice_dim, 0.0_RealExt)
call set_layer_scalar(diag%ice_re, 0.0_RealExt)
call set_layer_scalar(diag%ice_incloud_mmr, 0.0_RealExt)
call set_layer_scalar(diag%ice_mmr, 0.0_RealExt)
call set_layer_scalar(diag%ice_frac, 0.0_RealExt)
call set_layer_scalar(diag%ice_part_frac, 0.0_RealExt)
call set_layer_scalar(diag%liq_conv_dim, 0.0_RealExt)
call set_layer_scalar(diag%liq_inconv_mmr, 0.0_RealExt)
call set_layer_scalar(diag%liq_conv_mmr, 0.0_RealExt)
call set_layer_scalar(diag%liq_conv_frac, 0.0_RealExt)
call set_layer_scalar(diag%liq_conv_part_frac, 0.0_RealExt)
call set_layer_scalar(diag%ice_conv_dim, 0.0_RealExt)
call set_layer_scalar(diag%ice_conv_re, 0.0_RealExt)
call set_layer_scalar(diag%ice_inconv_mmr, 0.0_RealExt)
call set_layer_scalar(diag%ice_conv_mmr, 0.0_RealExt)
call set_layer_scalar(diag%ice_conv_frac, 0.0_RealExt)
call set_layer_scalar(diag%ice_conv_part_frac, 0.0_RealExt)
call set_layer_scalar(diag%cloud_solar_extinction, 0.0_RealExt)
call set_layer_scalar(diag%cloud_thermal_absorptivity, 0.0_RealExt)
call set_list_scalar(diag%liq_path, 0.0_RealExt)
call set_list_scalar(diag%ice_path, 0.0_RealExt)
call set_list_scalar(diag%liq_conv_path, 0.0_RealExt)
call set_list_scalar(diag%ice_conv_path, 0.0_RealExt)
call set_3d_scalar(diag%liq_subcol_scaling, 0.0_RealExt)
call set_3d_scalar(diag%ice_subcol_scaling, 0.0_RealExt)
if (control%i_cloud_representation /= ip_cloud_off) then
  do k=1, cld%n_condensed
    select case (cld%type_condensed(k))
    case (ip_clcmp_st_water)
      call set_cloud_prop(diag%liq_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%liq_incloud_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%liq_mmr, cld%condensed_mix_ratio(:,:,k), &
                          cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%liq_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%liq_part_frac, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_path(diag%liq_path, atm%mass, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)), &
                          cld%condensed_mix_ratio(:,:,k))
      if (control%i_inhom == ip_mcica) &
        call set_cloud_subcol(diag%liq_subcol_scaling, &
                              cld%c_sub(:, :, :, cld%i_cloud_type(k)))
    case (ip_clcmp_st_ice)
      call set_cloud_prop(diag%ice_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%ice_incloud_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%ice_mmr, cld%condensed_mix_ratio(:,:,k), &
                          cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%ice_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%ice_part_frac, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_path(diag%ice_path, atm%mass, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)), &
                          cld%condensed_mix_ratio(:,:,k))
      if (control%i_inhom == ip_mcica) &
        call set_cloud_subcol(diag%ice_subcol_scaling, &
                              cld%c_sub(:, :, :, cld%i_cloud_type(k)))
      if (associated(diag%ice_re)) then
        call set_cloud_prop(diag%ice_re, cld%condensed_dim_char(:,:,k))
        select case (cld%i_condensed_param(ip_clcmp_st_ice))
        case (ip_ice_adt, ip_ice_agg_de, ip_ice_agg_de_sun, ip_ice_baran)
          ! Here the effective dimension represents the diameter
          if (l_last) then
            diag%ice_re(:, list) = 0.5_RealExt * diag%ice_re(:, list)
          else
            diag%ice_re(list, :) = 0.5_RealExt * diag%ice_re(list, :)
          end if
        end select
      end if
    case (ip_clcmp_cnv_water)
      call set_cloud_prop(diag%liq_conv_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%liq_inconv_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%liq_conv_mmr, cld%condensed_mix_ratio(:,:,k), &
                          cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%liq_conv_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%liq_conv_part_frac, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_path(diag%liq_conv_path, atm%mass, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)), &
                          cld%condensed_mix_ratio(:,:,k))
    case (ip_clcmp_cnv_ice)
      call set_cloud_prop(diag%ice_conv_dim, cld%condensed_dim_char(:,:,k))
      call set_cloud_prop(diag%ice_inconv_mmr, cld%condensed_mix_ratio(:,:,k))
      call set_cloud_prop(diag%ice_conv_mmr, cld%condensed_mix_ratio(:,:,k), &
                          cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%ice_conv_frac, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_prop(diag%ice_conv_part_frac, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)))
      call set_cloud_path(diag%ice_conv_path, atm%mass, cld%w_cloud, &
                          cld%frac_cloud(:, :, cld%i_cloud_type(k)), &
                          cld%condensed_mix_ratio(:,:,k))
      if (associated(diag%ice_conv_re)) then
        call set_cloud_prop(diag%ice_conv_re, cld%condensed_dim_char(:,:,k))
        select case (cld%i_condensed_param(ip_clcmp_cnv_ice))
        case (ip_ice_adt, ip_ice_agg_de, ip_ice_agg_de_sun, ip_ice_baran)
          ! Here the effective dimension represents the diameter
          if (l_last) then
            diag%ice_conv_re(:, list) = 0.5_RealExt * diag%ice_conv_re(:, list)
          else
            diag%ice_conv_re(list, :) = 0.5_RealExt * diag%ice_conv_re(list, :)
          end if
        end select
      end if
    end select
  end do
end if

! Observed effective radius diagnostics
if (associated(diag%cloud_top_liq_dim) .and. &
    associated(diag%cloud_top_liq_weight)) then
  call cloud_level_diag(control, dimen, atm, cld, .true., &
    list, diag%cloud_top_liq_dim, diag%cloud_top_liq_weight)
end if
if (associated(diag%cloud_top_warm_liq_dim) .and. &
    associated(diag%cloud_top_warm_liq_weight)) then
  call cloud_level_diag(control, dimen, atm, cld, .false., &
    list, diag%cloud_top_warm_liq_dim, diag%cloud_top_warm_liq_weight)
end if

! Mean cloud absorptivity and extinction diagnostics
if (associated(diag%cloud_absorptivity)) &
  call set_layer_prop(diag%cloud_absorptivity, radout%cloud_absorptivity)

if (associated(diag%cloud_weight_absorptivity)) &
  call set_layer_prop(diag%cloud_weight_absorptivity, &
                      radout%cloud_weight_absorptivity)

if (associated(diag%cloud_extinction)) &
  call set_layer_prop(diag%cloud_extinction, radout%cloud_extinction)

if (associated(diag%cloud_weight_extinction)) &
  call set_layer_prop(diag%cloud_weight_extinction, &
                      radout%cloud_weight_extinction)

if (associated(diag%cloud_solar_extinction)) &
  call cloud_ext_diag(diag%cloud_solar_extinction, &
                      control, dimen, spectrum, atm, cld, &
                      list, layer_offset, l_last)

if (associated(diag%cloud_thermal_absorptivity)) &
  call cloud_abs_diag(diag%cloud_thermal_absorptivity, &
                      diag%cloud_absorptivity_wavelength, &
                      control, dimen, spectrum, atm, cld, &
                      list, layer_offset, l_last)

!------------------------------------------------------------------------------
! Aerosol diagnostics
!------------------------------------------------------------------------------
call set_aerosol_prop(diag%aerosol_optical_depth, &
                      radout%aerosol_absorption_band, &
                      radout%aerosol_scattering_band)
call set_aerosol_prop(diag%aerosol_scat_optical_depth, &
                      radout%aerosol_scattering_band)
call set_aerosol_prop(diag%aerosol_asymmetry_scat, &
                      radout%aerosol_asymmetry_band)

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
subroutine sum_flux_channels(field, field_channels, surf, toa)

  implicit none
  
  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in), allocatable :: field_channels(:, :, :)
  real(RealExt), intent(inout), optional, pointer :: surf(:)
  real(RealExt), intent(inout), optional, pointer :: toa(:)
  
  if (associated(field)) then
    if (l_last) then
      do i=0, n_layer
        ii = abs(level_offset-i)
        do l=1, n_profile
          field(ii, list(l)) &
            = real(sum(field_channels(l, i, 1:control%n_channel)), RealExt)
        end do
      end do
      if (present(surf)) then
        if (associated(surf)) then
          do l=1, n_profile
            surf(list(l)) = field(abs(level_offset-n_layer), list(l))
          end do
        end if
      end if
      if (present(toa)) then
        if (associated(toa)) then
          do l=1, n_profile
            toa(list(l)) = field(abs(level_offset), list(l))
          end do
        end if
      end if
    else
      do i=0, n_layer
        ii = abs(level_offset-i)
        do l=1, n_profile
          field(list(l), ii) &
            = real(sum(field_channels(l, i, 1:control%n_channel)), RealExt)
        end do
      end do
      if (present(surf)) then
        if (associated(surf)) then
          do l=1, n_profile
            surf(list(l)) = field(list(l), abs(level_offset-n_layer))
          end do
        end if
      end if
      if (present(toa)) then
        if (associated(toa)) then
          do l=1, n_profile
            toa(list(l)) = field(list(l), abs(level_offset))
          end do
        end if
      end if
    end if
  else
    if (present(surf)) then
      if (associated(surf)) then
        do l=1, n_profile
          surf(list(l)) &
            = real(sum(field_channels(l, n_layer, 1:control%n_channel)),RealExt)
        end do
      end if
    end if
    if (present(toa)) then
      if (associated(toa)) then
        do l=1, n_profile
          toa(list(l)) &
            = real(sum(field_channels(l, 0, 1:control%n_channel)), RealExt)
        end do
      end if
    end if
  end if

end subroutine sum_flux_channels


!------------------------------------------------------------------------------
subroutine sum_tile_channels(field, field_channels)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in), allocatable :: field_channels(:, :, :)

  integer :: j

  if (associated(field).and.present(n_tile)) then
    if (l_last) then
      do l=1, n_profile
        do i=1, n_tile
          field(i, list(l)) = 0.0_RealExt
        end do
      end do
    else
      do i=1, n_tile
        do l=1, n_profile
          field(list(l), i) = 0.0_RealExt
        end do
      end do
    end if
    if (control%l_tile) then
      if (l_last) then
        do i=1, n_tile
          do j=1, bound%n_point_tile
            l = bound%list_tile(j)
            field(i, list(l)) &
              = real(sum(field_channels(j, i, 1:control%n_channel)), RealExt)
          end do
        end do
      else
        do i=1, n_tile
          do j=1, bound%n_point_tile
            l = bound%list_tile(j)
            field(list(l), i) &
              = real(sum(field_channels(j, i, 1:control%n_channel)), RealExt)
          end do
        end do
      end if
    end if
  end if

end subroutine sum_tile_channels


!------------------------------------------------------------------------------
subroutine set_flux_bands(field, field_bands)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealK), intent(in), allocatable :: field_bands(:, :, :)
  integer :: k_lower, k_upper, k_dim
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
      k_dim = 2
    else
      i_dim = 2
      k_dim = 3
    end if
    k_lower = max(lbound(field, k_dim), 1)
    k_upper = min(ubound(field, k_dim), spectrum%basic%n_band)
    i_lower = max(lbound(field, i_dim), 0)
    i_upper = min(ubound(field, i_dim), n_layer)

    if (l_last) then
      do l=1, n_profile
        do k=k_lower, k_upper
          do i=lbound(field, i_dim), i_lower-1
            field(i, k, list(l)) = 0.0_RealExt
          end do
          do i=i_upper+1, ubound(field, i_dim)
            field(i, k, list(l)) = 0.0_RealExt
          end do
        end do
      end do
    else
      do k=k_lower, k_upper
        do i=lbound(field, i_dim), i_lower-1
          do l=1, n_profile
            field(list(l), i, k) = 0.0_RealExt
          end do
        end do
        do i=i_upper+1, ubound(field, i_dim)
          do l=1, n_profile
            field(list(l), i, k) = 0.0_RealExt
          end do
        end do
      end do
    end if
    if (l_last) then
      do k=k_lower, k_upper
        do i=i_lower, i_upper
          ii = abs(level_offset-i)
          do l=1, n_profile
            field(i, k, list(l)) = real(field_bands(l, ii, k), RealExt)
          end do
        end do
      end do
    else
      do k=k_lower, k_upper
        do i=i_lower, i_upper
          ii = abs(level_offset-i)
          do l=1, n_profile
            field(list(l), i, k) = real(field_bands(l, ii, k), RealExt)
          end do
        end do
      end do
    end if
  end if

end subroutine set_flux_bands


!------------------------------------------------------------------------------
subroutine sum_surf_weight_bands(field, field_bands)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:)
  real(RealK), intent(in), allocatable :: field_bands(:, :, :)

  if (associated(field)) then
    do l=1, n_profile
      field(list(l)) = 0.0_RealExt
    end do
    do k=1, spectrum%basic%n_band
      if (control%weight_diag(k) > 0.0) then
        do l=1, n_profile
          field(list(l)) = field(list(l)) &
            + real(control%weight_diag(k) * field_bands(l, n_layer, k), RealExt)
        end do
      end if
    end do
  end if

end subroutine sum_surf_weight_bands


!------------------------------------------------------------------------------
subroutine sum_photol_channels(field, field_channels, i_pathway)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in), allocatable :: field_channels(:, :, :, :)
  integer, intent(in) :: i_pathway

  real(RealK), allocatable :: photol(:, :)

  if (associated(field)) then
    allocate(photol(n_profile, n_layer))
    photol=0.0_RealK
    do j=1, control%n_channel
      do i=1, n_layer
        do l=1, n_profile
          photol(l, i) = photol(l, i) + field_channels(l, i, i_pathway, j)
        end do
      end do
    end do
    call set_layer_prop(field, photol)
    deallocate(photol)
  end if

end subroutine sum_photol_channels


!------------------------------------------------------------------------------
subroutine sum_photol_channels_all(field, field_channels)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealK), intent(in), allocatable :: field_channels(:, :, :, :)

  real(RealK), allocatable :: photol(:, :, :)

  if (associated(field)) then
    allocate(photol(n_profile, n_layer, spectrum%photol%n_pathway))
    photol=0.0_RealK
    do j=1, control%n_channel
      do k=1, spectrum%photol%n_pathway
        do i=1, n_layer
          do l=1, n_profile
            photol(l, i, k) = photol(l, i, k) + field_channels(l, i, k, j)
          end do
        end do
      end do
    end do
    call set_3d_prop(field, photol)
    deallocate(photol)
  end if

end subroutine sum_photol_channels_all


!------------------------------------------------------------------------------
subroutine set_cloud_prop(field, cld_field1, cld_field2, cld_field3)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in) :: cld_field1(:, dimen%id_cloud_top:)
  real(RealK), intent(in), optional :: cld_field2(:, dimen%id_cloud_top:)
  real(RealK), intent(in), optional :: cld_field3(:, dimen%id_cloud_top:)
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
    else
      i_dim = 2
    end if
    if (layer_offset == n_layer + 1) then
      ! Field output is inverted
      i_lower = max(n_layer + 1 - ubound(field, i_dim), dimen%id_cloud_top)
      i_upper = min(n_layer + 1 - lbound(field, i_dim), n_layer)
    else
      i_lower = max(lbound(field, i_dim), dimen%id_cloud_top)
      i_upper = min(ubound(field, i_dim), n_layer)
    end if
    ! Fill diagnostic between requested layers
    if (present(cld_field2).and.present(cld_field3)) then
      if (l_last) then
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(ii, list(l)) = real( &
              cld_field1(l, i) * cld_field2(l, i) * cld_field3(l, i), RealExt)
          end do
        end do
      else
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), ii) = real( &
              cld_field1(l, i) * cld_field2(l, i) * cld_field3(l, i), RealExt)
          end do
        end do
      end if
    else if (present(cld_field2)) then
      if (l_last) then
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(ii, list(l)) = real( &
              cld_field1(l, i) * cld_field2(l, i), RealExt)
          end do
        end do
      else
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), ii) = real( &
              cld_field1(l, i) * cld_field2(l, i), RealExt)
          end do
        end do
      end if
    else
      if (l_last) then
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(ii, list(l)) = real(cld_field1(l, i), RealExt)
          end do
        end do
      else
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), ii) = real(cld_field1(l, i), RealExt)
          end do
        end do
      end if
    end if
  end if

end subroutine set_cloud_prop


!------------------------------------------------------------------------------
subroutine set_cloud_path(field, d_mass, cld_field1, cld_field2, cld_field3)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:)
  real(RealK), intent(in) :: d_mass(:, :)
  real(RealK), intent(in) :: cld_field1(:, dimen%id_cloud_top:)
  real(RealK), intent(in) :: cld_field2(:, dimen%id_cloud_top:)
  real(RealK), intent(in) :: cld_field3(:, dimen%id_cloud_top:)

  if (associated(field)) then
    do i=dimen%id_cloud_top, n_layer
      do l=1, n_profile
        field(list(l)) = field(list(l)) + real( d_mass(l, i) &
          * cld_field1(l, i) * cld_field2(l, i) * cld_field3(l, i), RealExt)
      end do
    end do
  end if

end subroutine set_cloud_path


!------------------------------------------------------------------------------
subroutine set_cloud_subcol(field, cld_field)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealK), intent(in) :: cld_field(:, dimen%id_cloud_top:, :)
  integer :: j_lower, j_upper, j_dim, j
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
      j_dim = 2
    else
      i_dim = 2
      j_dim = 3
    end if
    j_lower = max(lbound(field, j_dim), 1)
    j_upper = min(ubound(field, j_dim), dimen%nd_subcol_gen)
    if (layer_offset == n_layer + 1) then
      ! Field output is inverted
      i_lower = max(n_layer + 1 - ubound(field, i_dim), dimen%id_cloud_top)
      i_upper = min(n_layer + 1 - lbound(field, i_dim), n_layer)
    else
      i_lower = max(lbound(field, i_dim), dimen%id_cloud_top)
      i_upper = min(ubound(field, i_dim), n_layer)
    end if
    ! Fill diagnostic between requested layers
    if (l_last) then
      do j=j_lower, j_upper
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(ii, j, list(l)) = real(cld_field(l, i, j), RealExt)
          end do
        end do
      end do
    else
      do j=j_lower, j_upper
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), ii, j) = real(cld_field(l, i, j), RealExt)
          end do
        end do
      end do
    end if
  end if

end subroutine set_cloud_subcol


!------------------------------------------------------------------------------
subroutine set_layer_prop(field, field_in)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealK), intent(in), allocatable :: field_in(:, :)
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
    else
      i_dim = 2
    end if
    i_lower = max(lbound(field, i_dim), 1)
    i_upper = min(ubound(field, i_dim), n_layer)
    ! Fill diagnostic between requested layers
    if (l_last) then
      do l=1, n_profile
        do i=lbound(field, i_dim), i_lower-1
          field(i, list(l)) = 0.0_RealExt
        end do
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          field(i, list(l)) = real(field_in(l, ii), RealExt)
        end do
        do i=i_upper+1, ubound(field, i_dim)
          field(i, list(l)) = 0.0_RealExt
        end do
      end do
    else
      do i=lbound(field, i_dim), i_lower-1
        do l=1, n_profile
          field(list(l), i) = 0.0_RealExt
        end do
      end do
      do i=i_lower, i_upper
        ii = abs(layer_offset-i)
        do l=1, n_profile
          field(list(l), i) = real(field_in(l, ii), RealExt)
        end do
      end do
      do i=i_upper+1, ubound(field, i_dim)
        do l=1, n_profile
          field(list(l), i) = 0.0_RealExt
        end do
      end do
    end if
  end if

end subroutine set_layer_prop


!------------------------------------------------------------------------------
subroutine set_3d_prop(field, field_in)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealK), intent(in), allocatable :: field_in(:, :, :)
  integer :: k_lower, k_upper, k_dim
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
      k_dim = 2
    else
      i_dim = 2
      k_dim = 3
    end if
    k_lower = max(lbound(field, k_dim), 1)
    k_upper = min(ubound(field, k_dim), size(field_in, 3))
    i_lower = max(lbound(field, i_dim), 1)
    i_upper = min(ubound(field, i_dim), n_layer)

    if (l_last) then
      do l=1, n_profile
        do k=k_lower, k_upper
          do i=lbound(field, i_dim), i_lower-1
            field(i, k, list(l)) = 0.0_RealExt
          end do
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
              field(i, k, list(l)) = real(field_in(l, ii, k), RealExt)
          end do
          do i=i_upper+1, ubound(field, i_dim)
            field(i, k, list(l)) = 0.0_RealExt
          end do
        end do
      end do
    else
      do k=k_lower, k_upper
        do i=lbound(field, i_dim), i_lower-1
          do l=1, n_profile
            field(list(l), i, k) = 0.0_RealExt
          end do
        end do
        do i=i_lower, i_upper
          ii = abs(layer_offset-i)
          do l=1, n_profile
            field(list(l), i, k) = real(field_in(l, ii, k), RealExt)
          end do
        end do
        do i=i_upper+1, ubound(field, i_dim)
          do l=1, n_profile
            field(list(l), i, k) = 0.0_RealExt
          end do
        end do
      end do
    end if
  end if

end subroutine set_3d_prop


!------------------------------------------------------------------------------
subroutine set_list_scalar(field, scalar)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:)
  real(RealExt), intent(in) :: scalar

  if (associated(field)) then
    do l=1, n_profile
      field(list(l)) = scalar
    end do
  end if

end subroutine set_list_scalar


!------------------------------------------------------------------------------
subroutine set_layer_scalar(field, scalar)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :)
  real(RealExt), intent(in) :: scalar

  if (associated(field)) then
    ! Fill diagnostic between requested layers
    if (l_last) then
      do l=1, n_profile
        do i=lbound(field, 1), ubound(field, 1)
          field(i, list(l)) = scalar
        end do
      end do
    else
      do i=lbound(field, 2), ubound(field, 2)
        do l=1, n_profile
          field(list(l), i) = scalar
        end do
      end do
    end if
  end if

end subroutine set_layer_scalar


!------------------------------------------------------------------------------
subroutine set_3d_scalar(field, scalar)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealExt), intent(in) :: scalar
  integer :: j

  if (associated(field)) then
    ! Fill diagnostic between requested layers
    if (l_last) then
      do l=1, n_profile
        do j=lbound(field, 2), ubound(field, 2)
          do i=lbound(field, 1), ubound(field, 1)
            field(i, j, list(l)) = scalar
          end do
        end do
      end do
    else
      do j=lbound(field, 3), ubound(field, 3)
        do i=lbound(field, 2), ubound(field, 2)
          do l=1, n_profile
            field(list(l), i, j) = scalar
          end do
        end do
      end do
    end if
  end if

end subroutine set_3d_scalar


!------------------------------------------------------------------------------
subroutine set_aerosol_prop(field, aer_field1, aer_field2)

  implicit none

  real(RealExt), intent(inout), pointer :: field(:, :, :)
  real(RealK), intent(in), allocatable :: aer_field1(:, :, :)
  real(RealK), intent(in), allocatable, optional :: aer_field2(:, :, :)
  integer :: k_lower, k_upper, k_dim
  integer :: i_lower, i_upper, i_dim

  if (associated(field)) then
    if (l_last) then
      i_dim = 1
      k_dim = 2
    else
      i_dim = 2
      k_dim = 3
    end if
    k_lower = max(lbound(field, k_dim), 1)
    k_upper = min(ubound(field, k_dim), spectrum%basic%n_band)
    i_lower = max(lbound(field, i_dim), 1)
    i_upper = min(ubound(field, i_dim), n_layer)

    if (l_last) then
      do l=1, n_profile
        do k=k_lower, k_upper
          do i=lbound(field, i_dim), i_lower-1
            field(i, k, list(l)) = 0.0_RealExt
          end do
          do i=i_upper+1, ubound(field, i_dim)
            field(i, k, list(l)) = 0.0_RealExt
          end do
        end do
      end do
    else
      do k=k_lower, k_upper
        do i=lbound(field, i_dim), i_lower-1
          do l=1, n_profile
            field(list(l), i, k) = 0.0_RealExt
          end do
        end do
        do i=i_upper+1, ubound(field, i_dim)
          do l=1, n_profile
            field(list(l), i, k) = 0.0_RealExt
          end do
        end do
      end do
    end if
    if (present(aer_field2)) then
      if (l_last) then
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(i, k, list(l)) = real(( aer_field1(l, ii, k) &
                + aer_field2(l, ii, k) ) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      else
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(list(l), i, k) = real(( aer_field1(l, ii, k) &
                + aer_field2(l, ii, k) ) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      end if
    else
      if (l_last) then
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(i, k, list(l)) = real( &
                aer_field1(l, ii, k) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      else
        do k=k_lower, k_upper
          do i=i_lower, i_upper
            ii = abs(layer_offset-i)
            do l=1, n_profile
              field(list(l), i, k) = real( &
                aer_field1(l, ii, k) * atm%mass(l, ii), RealExt)
            end do
          end do
        end do
      end if
    end if
  end if

end subroutine set_aerosol_prop

end subroutine set_diag
end module socrates_set_diag
