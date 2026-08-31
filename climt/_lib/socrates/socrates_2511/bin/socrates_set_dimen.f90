! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates dimension type
!
!------------------------------------------------------------------------------
module socrates_set_dimen
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_DIMEN'
contains

subroutine set_dimen(dimen, control, n_profile, n_layer, mcica_data, &
  n_channel, n_tile, n_cloud_layer, n_aer_mode, n_subcol_gen, &
  n_direction, n_viewing_level, n_brdf_basis_fnc, n_brdf_trunc, &
  n_profile_aerosol_prsc, n_profile_cloud_prsc, &
  n_opt_level_aerosol_prsc, n_opt_level_cloud_prsc)

use def_dimen,   only: StrDim
use def_control, only: StrCtrl
use def_mcica,   only: StrMcica, ip_mcica_full_sampling, &
  ip_mcica_single_sampling, ip_mcica_optimal_sampling
use rad_pcf,     only: &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, ip_cloud_csiw, &
  ip_cloud_combine_homogen, ip_cloud_combine_ice_water, &
  ip_cloud_split_homogen, ip_cloud_split_ice_water, &
  ip_cloud_column_max, ip_solver_mix_app_scat, ip_solver_mix_direct, &
  ip_solver_mix_direct_hogan, ip_solver_triple_app_scat, ip_solver_triple, &
  ip_solver_triple_hogan, ip_two_stream, ip_spherical_harmonic, &
  ip_sph_mode_flux, ip_trunc_triangular, ip_trunc_azim_sym, &
  i_normal, i_err_fatal
use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none

! Dimensions:
type(StrDim), intent(inout) :: dimen

! Control options:
type(StrCtrl), intent(in) :: control

! Mcica data:
type(StrMcica), intent(in), optional :: mcica_data

integer, intent(in) :: n_profile
integer, intent(in) :: n_layer
integer, intent(in), optional :: n_channel, n_tile
integer, intent(in), optional :: n_cloud_layer, n_aer_mode
integer, intent(in), optional :: n_subcol_gen
integer, intent(in), optional :: n_direction, n_viewing_level
integer, intent(in), optional :: n_brdf_basis_fnc, n_brdf_trunc
integer, intent(in), optional :: n_profile_aerosol_prsc
integer, intent(in), optional :: n_profile_cloud_prsc
integer, intent(in), optional :: n_opt_level_aerosol_prsc
integer, intent(in), optional :: n_opt_level_cloud_prsc

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_DIMEN'
character (len=errormessagelength) :: cmessage


dimen%nd_profile = n_profile
dimen%nd_layer = n_layer

if (present(n_channel)) then
  dimen%nd_channel = n_channel
else
  dimen%nd_channel = 1
end if


! Cloud
select case (control%i_cloud_representation)
case (ip_cloud_homogen, ip_cloud_ice_water, &
      ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
  dimen%nd_cloud_type      = 2
  dimen%nd_cloud_component = 2
case (ip_cloud_conv_strat, ip_cloud_csiw, &
      ip_cloud_split_homogen, ip_cloud_split_ice_water)
  dimen%nd_cloud_type      = 4
  dimen%nd_cloud_component = 4
case default
  dimen%nd_cloud_type      = 1
  dimen%nd_cloud_component = 1
end select
select case(control%i_cloud)
case (ip_cloud_column_max)
  if (present(n_cloud_layer)) then
    dimen%nd_column = 3 * n_cloud_layer + 2
  else
    dimen%nd_column = 3 * n_layer + 2
  end if
case default
  dimen%nd_column = 1
end select
if (present(mcica_data)) then
  if (present(n_subcol_gen)) then
    dimen%nd_subcol_gen = n_subcol_gen
  else
    dimen%nd_subcol_gen = mcica_data%n_subcol_gen
  end if
  select case (control%i_mcica_sampling)
  case (ip_mcica_full_sampling)
    dimen%nd_subcol_req = dimen%nd_subcol_gen
  case (ip_mcica_single_sampling)
    dimen%nd_subcol_req = mcica_data%n_subcol_req_single
  case (ip_mcica_optimal_sampling)
    dimen%nd_subcol_req = mcica_data%n_subcol_req_optimal
  case default
    dimen%nd_subcol_req = 1
  end select
else
  dimen%nd_subcol_gen = 1
  dimen%nd_subcol_req = 1
end if

if (control%l_cloud) then
  if (present(n_cloud_layer)) then
    dimen%id_cloud_top = dimen%nd_layer + 1 - max(1,n_cloud_layer)
    dimen%nd_layer_clr = dimen%nd_layer - max(0,n_cloud_layer)
  else
    dimen%id_cloud_top = 1
    dimen%nd_layer_clr = dimen%nd_layer
  end if
else
  dimen%id_cloud_top = dimen%nd_layer
  dimen%nd_layer_clr = dimen%nd_layer
end if


! Aerosol
if (present(n_aer_mode)) then
  dimen%nd_aerosol_mode = max(1,n_aer_mode)
else
  dimen%nd_aerosol_mode = 1
end if

! Arrays for prescribed optical properties
if (present(n_profile_aerosol_prsc)) then
  dimen%nd_profile_aerosol_prsc   = n_profile_aerosol_prsc
else
  dimen%nd_profile_aerosol_prsc   = 1
end if
if (present(n_profile_cloud_prsc)) then
  dimen%nd_profile_cloud_prsc     = n_profile_cloud_prsc
else
  dimen%nd_profile_cloud_prsc     = 1
end if
if (present(n_opt_level_aerosol_prsc)) then
  dimen%nd_opt_level_aerosol_prsc = n_opt_level_aerosol_prsc
else
  dimen%nd_opt_level_aerosol_prsc = 1
end if
if (present(n_opt_level_cloud_prsc)) then
  dimen%nd_opt_level_cloud_prsc   = n_opt_level_cloud_prsc
else
  dimen%nd_opt_level_cloud_prsc   = 1
end if


! Tiled surface.
if (control%l_tile) then
  dimen%nd_point_tile = max(1,n_profile)
  if (present(n_tile)) then
    dimen%nd_tile = n_tile
  else
    dimen%nd_tile = 1
  end if
else
  dimen%nd_point_tile = 1
  dimen%nd_tile       = 1
end if

! Solver dependent dimensions
if (control%i_angular_integration == ip_two_stream) then

  dimen%nd_viewing_level    = 1
  dimen%nd_radiance_profile = 1
  dimen%nd_j_profile        = 1
  dimen%nd_direction        = 1
  dimen%nd_brdf_basis_fnc   = 2
  dimen%nd_brdf_trunc       = 1
  dimen%nd_flux_profile     = dimen%nd_profile
  dimen%nd_2sg_profile      = dimen%nd_profile
  dimen%nd_source_coeff     = 2
  dimen%nd_max_order        = 1
  dimen%nd_sph_coeff        = 1

  select case(control%i_solver)
  case(ip_solver_mix_app_scat, ip_solver_mix_direct, ip_solver_mix_direct_hogan)
    dimen%nd_overlap_coeff=8
    dimen%nd_region=2
  case(ip_solver_triple_app_scat, ip_solver_triple, ip_solver_triple_hogan)
    dimen%nd_overlap_coeff=18
    dimen%nd_region=3
  case default
    dimen%nd_overlap_coeff=1
    dimen%nd_region=2
  end select

else if (control%i_angular_integration == ip_spherical_harmonic) then

  if (present(n_direction)) then
    dimen%nd_direction = n_direction
  else
    dimen%nd_direction = 1
  end if
  if (present(n_brdf_basis_fnc)) then
    dimen%nd_brdf_basis_fnc = n_brdf_basis_fnc
  else
    dimen%nd_brdf_basis_fnc = 2
  end if
  if (present(n_brdf_trunc)) then
    dimen%nd_brdf_trunc = n_brdf_trunc
  else
    dimen%nd_brdf_trunc = 1
  end if
  dimen%nd_radiance_profile = dimen%nd_profile
  dimen%nd_j_profile        = 1

  if (control%i_sph_mode == ip_sph_mode_flux) then
    dimen%nd_viewing_level = n_layer+1
    dimen%nd_flux_profile  = dimen%nd_profile
  else
    if (present(n_viewing_level)) then
      dimen%nd_viewing_level = n_viewing_level
    else
      dimen%nd_viewing_level = 1
    end if
    dimen%nd_flux_profile = 1
  end if

  dimen%nd_2sg_profile   = 1
  dimen%nd_source_coeff  = 1
  dimen%nd_overlap_coeff = 1
  dimen%nd_region        = 1
  dimen%nd_max_order     = control%ls_global_trunc + 2
  if (control%i_truncation == ip_trunc_triangular) then
    dimen%nd_sph_coeff = &
      (control%ls_global_trunc+3)*(control%ls_global_trunc+4)/2
  else if (control%i_truncation == ip_trunc_azim_sym) then
    dimen%nd_sph_coeff = control%ls_global_trunc+2
  else
    cmessage = 'Illegal truncation'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if

end if

end subroutine set_dimen
end module socrates_set_dimen
