! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates control type
!
!------------------------------------------------------------------------------
module socrates_set_control
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CONTROL'
contains

subroutine set_control(control, diag, spectrum, l_set_defaults, &
  l_rayleigh, l_gas, l_continuum, l_cont_gen, l_orog, l_solvar, &
  l_rescale, l_ir_source_quad, l_mixing_ratio, &
  l_aerosol, l_aerosol_mode, l_aerosol_ccn, &
  l_tile, l_flux_ground, l_flux_tile, n_tile, n_cloud_layer, n_aer_mode, &
  isolir, i_scatter_method, i_gas_overlap, &
  i_cloud_representation, i_overlap, i_inhom, &
  i_cloud_entrapment, i_mcica_sampling, &
  i_st_water, i_cnv_water, i_st_ice, i_cnv_ice, i_drop_re )

use def_control, only: StrCtrl, allocate_control
use socrates_def_diag, only: StrDiag
use def_spectrum, only: StrSpecData
use realtype_rd, only: RealK, RealExt
use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength
use missing_data_mod, only: imdi
use rad_pcf, only: &
  ip_solar, ip_pifm80, ip_scatter_full, ip_infra_red, ip_elsasser, &
  ip_two_stream, ip_ir_gauss, ip_spherical_harmonic, ip_overlap_k_eqv_scl, &
  ip_cloud_off, ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, &
  ip_cloud_csiw, ip_cloud_combine_homogen, ip_cloud_combine_ice_water, &
  ip_cloud_split_homogen, ip_cloud_split_ice_water, &
  ip_max_rand, ip_rand, ip_exponential_rand, ip_homogeneous, &
  ip_scaling, ip_mcica, ip_cairns, ip_cloud_ice_water, ip_cloud_mcica, &
  ip_zero_entrapment, ip_max_entrapment, &
  ip_no_scatter_abs, ip_no_scatter_ext, ip_solver_no_scat, &
  ip_solver_homogen_direct, ip_scatter_approx, ip_solver_mix_app_scat, &
  ip_solver_homogen_direct, ip_solver_mix_direct, ip_solver_mix_direct_hogan, &
  ip_cloud_mix_max, ip_cloud_part_corr, ip_cloud_mix_random, &
  ip_solver_triple_app_scat, ip_solver_triple, &
  ip_solver_triple_hogan, ip_cloud_triple, ip_cloud_part_corr_cnv, &
  ip_cloud_clear, ip_scale_ses2, ip_overlap_mix_ses2, ip_re_external, &
  i_normal, i_err_fatal
use gas_list_pcf, only: npd_gases, npd_products
use def_mcica, only: ip_mcica_optimal_sampling
use socrates_set_spectrum, only: set_weight_blue

implicit none

! Control options:
type(StrCtrl), intent(inout) :: control

! Diagnostic fields
type(StrDiag), intent(in) :: diag

! Spectral data:
type(StrSpecData), intent(in), optional :: spectrum

logical, intent(in), optional :: l_set_defaults, &
  l_rayleigh, l_gas, l_continuum, l_cont_gen, l_orog, l_solvar, &
  l_rescale, l_ir_source_quad, l_mixing_ratio, &
  l_aerosol, l_aerosol_mode, l_aerosol_ccn, &
  l_tile, l_flux_ground, l_flux_tile(:)

integer, intent(in), optional :: n_tile, n_cloud_layer, n_aer_mode

integer, intent(in), optional :: isolir, i_scatter_method, i_gas_overlap, &
  i_cloud_representation, i_overlap, i_inhom, &
  i_cloud_entrapment, i_mcica_sampling, &
  i_st_water, i_cnv_water, i_st_ice, i_cnv_ice, i_drop_re

! Local variables
integer :: i, k
integer :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CONTROL_DEFAULTS'
character (len=errormessagelength) :: cmessage


! Logical options
if (present(l_rayleigh)) control%l_rayleigh = l_rayleigh
if (present(l_gas)) control%l_gas = l_gas
if (present(l_continuum)) control%l_continuum = l_continuum
if (present(l_cont_gen)) control%l_cont_gen = l_cont_gen
if (present(l_orog)) then
  control%l_orog = l_orog
  control%l_orog_fix = l_orog
  if (l_orog) control%l_flux_div = .true.
end if
if (present(l_solvar)) control%l_solvar = l_solvar
if (present(l_rescale)) control%l_rescale = l_rescale
if (present(l_ir_source_quad)) control%l_ir_source_quad = l_ir_source_quad
if (present(l_mixing_ratio)) control%l_mixing_ratio = l_mixing_ratio
if (present(l_aerosol)) control%l_aerosol = l_aerosol
if (present(l_aerosol_mode)) control%l_aerosol_mode = l_aerosol_mode
if (present(l_aerosol_ccn)) control%l_aerosol_ccn = l_aerosol_ccn
if (present(l_tile)) then
  control%l_tile = l_tile
  control%l_tile_emissivity = l_tile
end if
if (present(l_flux_ground)) control%l_flux_ground = l_flux_ground


! Integer options
if (present(isolir)) control%isolir = isolir
if (present(i_scatter_method)) control%i_scatter_method = i_scatter_method
if (present(i_gas_overlap)) control%i_gas_overlap = i_gas_overlap
if (present(i_cloud_representation)) &
  control%i_cloud_representation = i_cloud_representation
if (present(i_overlap)) control%i_overlap = i_overlap
if (present(i_inhom)) control%i_inhom = i_inhom
if (present(i_cloud_entrapment)) &
  control%i_cloud_entrapment = i_cloud_entrapment
if (present(i_mcica_sampling)) control%i_mcica_sampling = i_mcica_sampling
if (present(i_st_water)) control%i_st_water = i_st_water
if (present(i_cnv_water)) control%i_cnv_water = i_cnv_water
if (present(i_st_ice)) control%i_st_ice = i_st_ice
if (present(i_cnv_ice)) control%i_cnv_ice = i_cnv_ice
if (present(i_drop_re)) control%i_drop_re = i_drop_re


! Diagnostic options
control%l_flux_direct_band = associated(diag%flux_direct_band)
control%l_flux_down_band = associated(diag%flux_down_band)
control%l_flux_up_band = associated(diag%flux_up_band)
control%l_flux_direct_clear_band = associated(diag%flux_direct_clear_band)
control%l_flux_down_clear_band = associated(diag%flux_down_clear_band)
control%l_flux_up_clear_band = associated(diag%flux_up_clear_band)

control%l_photolysis_rate = associated(diag%photolysis_rate)
do k=1, spectrum%photol%n_pathway
  if (associated(diag%photolysis(k)%rate)) &
    control%l_photolysis_rate = .true.
  if (associated(diag%photolysis(k)%rate_clear)) &
    control%l_photolysis_rate_clear = .true.
  if (associated(diag%photolysis(k)%rate_clean)) &
    control%l_photolysis_rate_clean = .true.
  if (associated(diag%photolysis(k)%rate_clear_clean)) &
    control%l_photolysis_rate_clear_clean = .true.
end do

if (associated(diag%flux_direct_clear)      .or. &
    associated(diag%flux_down_clear)        .or. &
    associated(diag%flux_up_clear)          .or. &
    associated(diag%flux_direct_clear_surf) .or. &
    associated(diag%flux_down_clear_surf)   .or. &
    associated(diag%flux_up_clear_surf)     .or. &
    associated(diag%flux_direct_clear_toa)  .or. &
    associated(diag%flux_down_clear_toa)    .or. &
    associated(diag%flux_up_clear_toa)      .or. &
    control%l_flux_direct_clear_band        .or. &
    control%l_flux_down_clear_band          .or. &
    control%l_flux_up_clear_band            .or. &
    control%l_photolysis_rate_clear) then
  control%l_clear = .true.
end if

if (associated(diag%flux_direct_clean)      .or. &
    associated(diag%flux_down_clean)        .or. &
    associated(diag%flux_up_clean)          .or. &
    associated(diag%flux_direct_clean_surf) .or. &
    associated(diag%flux_down_clean_surf)   .or. &
    associated(diag%flux_up_clean_surf)     .or. &
    associated(diag%flux_direct_clean_toa)  .or. &
    associated(diag%flux_down_clean_toa)    .or. &
    associated(diag%flux_up_clean_toa)      .or. &
    associated(diag%flux_direct_clean_band) .or. &
    associated(diag%flux_down_clean_band)   .or. &
    associated(diag%flux_up_clean_band)     .or. &
    control%l_photolysis_rate_clean) then
  control%l_clean = .true.
end if

if (associated(diag%flux_direct_clear_clean)      .or. &
    associated(diag%flux_down_clear_clean)        .or. &
    associated(diag%flux_up_clear_clean)          .or. &
    associated(diag%flux_direct_clear_clean_surf) .or. &
    associated(diag%flux_down_clear_clean_surf)   .or. &
    associated(diag%flux_up_clear_clean_surf)     .or. &
    associated(diag%flux_direct_clear_clean_toa)  .or. &
    associated(diag%flux_down_clear_clean_toa)    .or. &
    associated(diag%flux_up_clear_clean_toa)      .or. &
    associated(diag%flux_direct_clear_clean_band) .or. &
    associated(diag%flux_down_clear_clean_band)   .or. &
    associated(diag%flux_up_clear_clean_band)     .or. &
    control%l_photolysis_rate_clear_clean) then
  control%l_clear_clean = .true.
end if

if (associated(diag%flux_up_blue_tile)     .or. &
    associated(diag%flux_direct_blue_surf) .or. &
    associated(diag%flux_down_blue_surf) ) then
  control%l_blue_flux_surf = .true.
end if

if (associated(diag%flux_direct_uv_surf)) &
          control%l_flux_direct_band = .true.
if (associated(diag%flux_down_uv_surf)) &
          control%l_flux_down_band = .true.
if (associated(diag%flux_up_uv_surf)) &
          control%l_flux_up_band = .true.
if (associated(diag%flux_direct_uv_clear_surf)) &
          control%l_flux_direct_clear_band = .true.
if (associated(diag%flux_down_uv_clear_surf)) &
          control%l_flux_down_clear_band = .true.
if (associated(diag%flux_up_uv_clear_surf)) &
          control%l_flux_up_clear_band = .true.

if (associated(diag%aerosol_optical_depth)) then
  control%l_aerosol_absorption_band = .true.
  control%l_aerosol_scattering_band = .true.
end if
if (associated(diag%aerosol_scat_optical_depth)) then
  control%l_aerosol_scattering_band = .true.
end if
if (associated(diag%aerosol_asymmetry_scat)) then
  control%l_aerosol_asymmetry_band = .true.
end if

if (associated(diag%cloud_absorptivity) .or. &
    associated(diag%cloud_weight_absorptivity)) then
  control%l_cloud_absorptivity = .true.
end if
if (associated(diag%cloud_extinction) .or. &
    associated(diag%cloud_weight_extinction)) then
  control%l_cloud_extinction = .true.
end if

! Defaults and checking
if (present(l_set_defaults)) then
  if (l_set_defaults) then

    ! Two-stream is the default
    call set_int_default(control%i_angular_integration, ip_two_stream)

    select case (control%i_angular_integration)
    case(ip_two_stream)

      ! Source-specific options
      select case(control%isolir)
      case(ip_solar)
        call set_int_default(control%i_2stream,        ip_pifm80)
        call set_int_default(control%i_scatter_method, ip_scatter_full)
      case(ip_infra_red)
        call set_int_default(control%i_2stream,        ip_elsasser)
        call set_int_default(control%i_scatter_method, ip_scatter_full)
        if (.not.present(l_ir_source_quad)) control%l_ir_source_quad = .true.
      case default
        cmessage = 'Radiative source option unrecognised.'
        ierr=i_err_fatal
        CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)
      end select

      ! Set general default options
      if (.not.present(l_rescale)) control%l_rescale = .true.
      if (.not.present(l_gas)) control%l_gas = .true.
      if (.not.present(l_continuum)) control%l_continuum = .true.
      if (.not.present(l_cont_gen)) control%l_cont_gen = .true.
      if (.not.present(l_aerosol)) control%l_aerosol = .true.
      call set_int_default(control%first_band, 1)
      if (present(spectrum)) then
        call set_int_default(control%last_band, spectrum%basic%n_band)
      end if
      call set_int_default(control%n_channel, 1)
      call set_int_default(control%n_order_forward, 2)
      call set_int_default(control%i_gas_overlap, ip_overlap_k_eqv_scl)

      ! Prevent unnecessary calculation of equivalent extinction when
      ! there is only one term in the LW
      if (control%i_scatter_method /= ip_scatter_full) then
        ! Not done for full scattering to maintain LFRic KGO
        control%l_grey_single = .true.
      end if

      ! Consistent tiling options
      if (.not.present(l_tile)) then
        control%l_tile = .true.
        ! Bug-fix changes KGO:
        ! control%l_tile_emissivity = .true.
      end if
      if (present(n_tile)) then
        if (n_tile < 1) control%l_tile = .false.
      else
        control%l_tile = .false.
      end if

      ! Consistent cloud options and defaults
      call set_int_default(control%i_cloud_representation, ip_cloud_off)
      call set_int_default(control%i_overlap, ip_max_rand)
      call set_int_default(control%i_inhom, ip_homogeneous)
      call set_int_default(control%i_cloud_entrapment, ip_zero_entrapment)
      call set_int_default(control%i_mcica_sampling, ip_mcica_optimal_sampling)
      call set_int_default(control%i_drop_re, ip_re_external)
      if (present(n_cloud_layer)) then
        control%l_global_cloud_top = .true.
        if (n_cloud_layer < 1) control%i_cloud_representation = ip_cloud_off
      end if
      select case(control%i_cloud_representation)
      case(ip_cloud_homogen, ip_cloud_ice_water, &
           ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
        control%l_cloud = .true.
        if (.not.control%l_drop .and. .not.control%l_ice) then
          control%l_drop = .true.
          control%l_ice  = .true.
        end if
        if (control%i_inhom == ip_mcica) then
          control%i_cloud = ip_cloud_mcica
          if ( (control%i_scatter_method == ip_no_scatter_abs) .or. &
               (control%i_scatter_method == ip_no_scatter_ext) ) then
            control%i_solver       = ip_solver_no_scat
            control%i_solver_clear = ip_solver_no_scat
          else
            control%i_solver       = ip_solver_homogen_direct
            control%i_solver_clear = ip_solver_homogen_direct
          end if
        else
          if (control%i_scatter_method == ip_scatter_approx) then
            control%i_solver       = ip_solver_mix_app_scat
          else if (control%i_cloud_entrapment == ip_max_entrapment) then
            control%i_solver       = ip_solver_mix_direct
          else
            control%i_solver       = ip_solver_mix_direct_hogan
          end if
          control%i_solver_clear = ip_solver_homogen_direct
          if (control%i_overlap == ip_max_rand) then
            control%i_cloud = ip_cloud_mix_max
          else if (control%i_overlap == ip_exponential_rand) then
            control%i_cloud = ip_cloud_part_corr
          else
            control%i_cloud = ip_cloud_mix_random
          end if
        end if
      case(ip_cloud_conv_strat, ip_cloud_csiw, &
           ip_cloud_split_homogen, ip_cloud_split_ice_water)
        ! Not compatible with control%i_inhom == ip_mcica
        if (control%i_inhom == ip_mcica) control%i_inhom = ip_homogeneous
        control%l_cloud = .true.
        if (.not.control%l_drop .and. .not.control%l_ice) then
          control%l_drop = .true.
          control%l_ice  = .true.
        end if
        if (control%i_scatter_method == ip_scatter_approx) then
          control%i_solver       = ip_solver_triple_app_scat
        else if (control%i_cloud_entrapment == ip_max_entrapment) then
          control%i_solver       = ip_solver_triple
        else
          control%i_solver       = ip_solver_triple_hogan
        end if
        control%i_solver_clear = ip_solver_homogen_direct
        if (control%i_overlap == ip_max_rand) then
          control%i_cloud = ip_cloud_triple
        else if (control%i_overlap == ip_exponential_rand) then
          control%i_cloud = ip_cloud_part_corr_cnv
        else
          control%i_cloud = ip_cloud_mix_random
        end if
      case(ip_cloud_off)
        ! No treatment of cloud
        control%l_cloud        = .false.
        control%l_drop         = .false.
        control%l_ice          = .false.
        control%l_microphysics = .false.
        control%i_cloud        = ip_cloud_clear
        if ( (control%i_scatter_method == ip_no_scatter_abs) .or. &
             (control%i_scatter_method == ip_no_scatter_ext) ) then
          control%i_solver       = ip_solver_no_scat
          control%i_solver_clear = ip_solver_no_scat
        else
          control%i_solver       = ip_solver_homogen_direct
          control%i_solver_clear = ip_solver_homogen_direct
        end if
      end select

      ! Consistent aerosol options
      if (.not.present(l_aerosol_mode)) control%l_aerosol_mode = .true.
      if (present(n_aer_mode)) then
        if (n_aer_mode < 1) control%l_aerosol_mode = .false.
      else
        control%l_aerosol_mode = .false.
      end if

    case(ip_ir_gauss, ip_spherical_harmonic)

      ! Currently no checking done or defaults set

    case default

      cmessage = 'Method of angular integration unrecognised.'
      ierr=i_err_fatal
      CALL ereport(ModuleName//':'//RoutineName, ierr, cmessage)

    end select

  end if
end if


if (present(spectrum)) then
  ! Allocate band-by-band control options
  call allocate_control(control, spectrum)
  
  ! Set properties for individual bands.
  if (control%n_channel == control%last_band - control%first_band + 1) then
    do i = 1, control%n_channel
      control%map_channel(control%first_band + i-1) = i
    end do
  else
    do i = 1, spectrum%basic%n_band
      control%map_channel(i) = 1
    end do
  end if
  do i = 1, spectrum%basic%n_band
    control%weight_band(i)           = 1.0_RealK
    control%weight_diag(i)           = 0.0_RealK
    control%i_scatter_method_band(i) = control%i_scatter_method
    control%i_gas_overlap_band(i)    = control%i_gas_overlap
    if (any(spectrum%gas%i_scale_fnc(i,:) == ip_scale_ses2)) then
      ! If SES2 scaling is used in the band the overlap must also use SES2:
      control%i_gas_overlap_band(i)  = ip_overlap_mix_ses2
    end if
  end do

  ! UV flux diagnostics
  if (control%isolir == ip_solar) then
    if (associated(diag%flux_direct_uv_surf) .or. &
        associated(diag%flux_down_uv_surf) .or. &
        associated(diag%flux_up_uv_surf) .or. &
        associated(diag%flux_direct_uv_clear_surf) .or. &
        associated(diag%flux_down_uv_clear_surf) .or. &
        associated(diag%flux_up_uv_clear_surf)) then
      call set_weight_blue(spectrum, control%weight_diag, diag%wavelength_uv)
    end if
    if (associated(diag%flux_direct_uv_clear_surf) .or. &
        associated(diag%flux_down_uv_clear_surf) .or. &
        associated(diag%flux_up_uv_clear_surf)) then
      do i = 1, spectrum%basic%n_band
        control%l_clear_band(i) = control%weight_diag(i) > 0.0
      end do
    end if
  end if

end if

if (present(n_tile)) then
  call allocate_control(control=control, n_tile=n_tile)
  if (present(l_flux_tile)) control%l_flux_tile = l_flux_tile
end if

contains

  subroutine set_int_default(i_option, i_default)
    implicit none
    integer, intent(inout) :: i_option
    integer, intent(in) :: i_default
    if (i_option == imdi) then
      i_option = i_default
    end if
  end subroutine set_int_default

end subroutine set_control
end module socrates_set_control
