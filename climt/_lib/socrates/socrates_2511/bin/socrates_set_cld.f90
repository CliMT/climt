! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the variables in the Socrates cloud type
!
!------------------------------------------------------------------------------
module socrates_set_cld
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CLD'
contains

subroutine set_cld(cld, control, dimen, spectrum, atm, &
  profile_list, n_layer_stride, cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_rsd, ice_rsd, liq_conv_rsd, ice_conv_rsd, &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_rsd_1d, ice_rsd_1d, liq_conv_rsd_1d, ice_conv_rsd_1d, &
  cloud_vertical_decorr, conv_vertical_decorr, cloud_horizontal_rsd, &
  l_invert, l_profile_last, l_debug, i_profile_debug)

use def_cld,      only: StrCld, allocate_cld, allocate_cld_prsc
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use realtype_rd,  only: RealK, RealExt
use rad_pcf,      only: &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat, ip_cloud_csiw, &
  ip_cloud_combine_homogen, ip_cloud_combine_ice_water, &
  ip_cloud_split_homogen, ip_cloud_split_ice_water, &
  ip_clcmp_st_water, ip_clcmp_st_ice, ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
  ip_phase_water, ip_phase_ice, ip_cloud_type_homogen, &
  ip_cloud_type_water, ip_cloud_type_ice, &
  ip_cloud_type_strat, ip_cloud_type_conv, &
  ip_cloud_type_sw, ip_cloud_type_si, ip_cloud_type_cw, ip_cloud_type_ci, &
  ip_drop_unparametrized, ip_ice_unparametrized, &
  ip_scaling, ip_cairns, ip_mcica, ip_tripleclouds_2019, &
  i_normal, i_err_fatal
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none


! Cloud properties:
type(StrCld),      intent(out) :: cld

! Control options:
type(StrCtrl),     intent(in)  :: control

! Dimensions:
type(StrDim),      intent(in)  :: dimen

! Spectral data:
type(StrSpecData), intent(in)  :: spectrum

! Atmospheric properties:
type(StrAtm),      intent(in)  :: atm

integer, intent(in), optional :: profile_list(:)
!   List of profiles to use from input fields

integer, intent(in), optional :: n_layer_stride
!   Number of layers in input 1d arrays

real(RealExt), intent(in), dimension(:, :), optional :: &
  cloud_frac, conv_frac, &
  liq_frac, ice_frac, liq_conv_frac, ice_conv_frac, &
  liq_mmr, ice_mmr, liq_conv_mmr, ice_conv_mmr, &
  liq_rsd, ice_rsd, liq_conv_rsd, ice_conv_rsd
!   Liquid and ice cloud fractions, gridbox mean mixing ratios,
!   and relative standard deviation of condensate

real(RealExt), intent(in), dimension(:), optional :: &
  cloud_frac_1d, conv_frac_1d, &
  liq_frac_1d, ice_frac_1d, liq_conv_frac_1d, ice_conv_frac_1d, &
  liq_mmr_1d, ice_mmr_1d, liq_conv_mmr_1d, ice_conv_mmr_1d, &
  liq_rsd_1d, ice_rsd_1d, liq_conv_rsd_1d, ice_conv_rsd_1d
!   Liquid and ice cloud fractions, gridbox mean mixing ratios,
!   and relative standard deviation of condensate input as 1d fields

real(RealExt), intent(in), optional :: cloud_vertical_decorr
!   Decorrelation pressure scale for cloud vertical overlap
real(RealExt), intent(in), optional :: conv_vertical_decorr
!   Decorrelation pressure scale for convective cloud vertical overlap
real(RealExt), intent(in), optional :: cloud_horizontal_rsd
!   Relative standard deviation of sub-grid cloud condensate

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information


! Local variables
integer :: i, j, k, kk, l, ll, list(dimen%nd_profile)
!   Loop variables
integer :: i_phase, i_param_type, n_cloud_parameter, i_thin, i_thick
integer :: layer_offset, stride_layer
!   Working variables
real(RealK), dimension(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer) :: &
  cond_mmr, cond_rsd, frac, frac_liq, frac_ice, frac_thin, ratio_thin
!   Working arrays

logical :: l_last
!   Local flag to loop over profiles last in 1d fields
logical :: l_combine
!   Combine stratiform and convective cloud
logical :: l_split
!   Split cloud into optically thick and thin regions

real(RealK) :: eps = EPSILON(1.0_RealK)
real(RealK) :: min_cloud_fraction = 0.001_RealK
real(RealK) :: frac_c, frac_m

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CLD'
character (len=errormessagelength) :: cmessage

! Functions called
integer, external :: set_n_cloud_parameter


! Allocate structure for the core radiation code interface
call allocate_cld(cld, dimen, spectrum)
call allocate_cld_prsc(cld, dimen, spectrum)

if (.not.control%l_cloud) then
  return
end if

if (present(profile_list)) then
  list(1:atm%n_profile) = profile_list(1:atm%n_profile)
else
  do l=1, atm%n_profile
    list(l) = l
  end do
end if

layer_offset = 0
if (present(l_invert)) then
  if (l_invert) then
    ! The layer is indexed using an inverted loop counter
    layer_offset = atm%n_layer + 1
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
  stride_layer = atm%n_layer
end if

!------------------------------------------------------------------------------
! Set properties of condensed components
!------------------------------------------------------------------------------

if (control%l_ice .and. control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water, &
        ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_st_ice
  case (ip_cloud_conv_strat, ip_cloud_csiw, &
        ip_cloud_split_homogen, ip_cloud_split_ice_water)
    cld%n_condensed = 4
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_st_ice
    cld%type_condensed(3) = ip_clcmp_cnv_water
    cld%type_condensed(4) = ip_clcmp_cnv_ice
  end select
else if (control%l_ice .and. .not.control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water, &
        ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
    cld%n_condensed = 1
    cld%type_condensed(1) = ip_clcmp_st_ice
  case (ip_cloud_conv_strat, ip_cloud_csiw, &
        ip_cloud_split_homogen, ip_cloud_split_ice_water)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_ice
    cld%type_condensed(2) = ip_clcmp_cnv_ice
  end select
else if (.not.control%l_ice .and. control%l_drop) then
  select case (control%i_cloud_representation)
  case (ip_cloud_homogen, ip_cloud_ice_water, &
        ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
    cld%n_condensed = 1
    cld%type_condensed(1) = ip_clcmp_st_water
  case (ip_cloud_conv_strat, ip_cloud_csiw, &
        ip_cloud_split_homogen, ip_cloud_split_ice_water)
    cld%n_condensed = 2
    cld%type_condensed(1) = ip_clcmp_st_water
    cld%type_condensed(2) = ip_clcmp_cnv_water
  end select
else
  cmessage = 'Cloud on, but no condensed components included.'
  ierr=i_err_fatal
  call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
end if

select case (control%i_cloud_representation)
case (ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
  ! Combine stratiform and convective cloud into a single region
  l_combine = .true.
  l_split = .false.
case (ip_cloud_split_homogen, ip_cloud_split_ice_water)
  ! Stratiform and convective cloud is combined
  l_combine = .true.
  ! Total cloud is then split into optically thick and thin regions
  l_split = .true.
case default
  l_combine = .false.
  l_split = .false.
end select

! Initialise convective cloud
cld%c_cloud = 0.0_RealK
cld%c_ratio = 0.0_RealK

do i=1, cld%n_condensed
  select case (cld%type_condensed(i))
  case (ip_clcmp_st_water)
    i_phase = ip_phase_water
    i_param_type = control%i_st_water
  case (ip_clcmp_st_ice)
    i_phase = ip_phase_ice
    i_param_type = control%i_st_ice
  case (ip_clcmp_cnv_water)
    i_phase = ip_phase_water
    i_param_type = control%i_cnv_water
  case (ip_clcmp_cnv_ice)
    i_phase = ip_phase_ice
    i_param_type = control%i_cnv_ice
  end select

  select case (i_phase)
  case (ip_phase_water)
    if (i_param_type <= 0) then
      cld%i_condensed_param(i) = ip_drop_unparametrized
      cmessage = 'Prescribed liquid cloud not yet implemented.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (i_param_type > spectrum%dim%nd_drop_type) then
      cmessage = 'Liquid cloud type outside allowed range.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (spectrum%drop%l_drop_type(i_param_type)) then
      ! Take parametrisation from spectral file
      cld%i_condensed_param(i) = spectrum%drop%i_drop_parm(i_param_type)
      cld%condensed_n_phf(i) = spectrum%drop%n_phf(i_param_type)
      ! DEPENDS ON: set_n_cloud_parameter
      n_cloud_parameter = set_n_cloud_parameter( cld%i_condensed_param(i), &
        cld%type_condensed(i), cld%condensed_n_phf(i) )
      do j=1, spectrum%basic%n_band
        do k=1, n_cloud_parameter
          cld%condensed_param_list(k, i, j) &
            = spectrum%drop%parm_list(k, j, i_param_type)
        end do
      end do

      ! Assign liquid mass mixing ratio, effective dimension,
      ! and relative standard deviation
      select case (cld%type_condensed(i))

      case (ip_clcmp_st_water)
        if (present(liq_mmr).or.present(liq_mmr_1d)) then
          ! Set the liquid condensate mass mixing ratio
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             liq_mmr, liq_mmr_1d)
          if (l_combine .and. &
              (present(liq_conv_mmr).or.present(liq_conv_mmr_1d))) then
            ! Add in the convective condensate if using combined cloud
            call set_cld_field(cond_mmr, liq_conv_mmr, liq_conv_mmr_1d)
            cld%condensed_mix_ratio(:, :, i) = &
              cld%condensed_mix_ratio(:, :, i) + cond_mmr
            ! Increment the convective condensate for the ratio calculation
            cld%c_ratio = cld%c_ratio + cond_mmr
          end if
        else
          cmessage = 'Liquid MMR not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
        if (present(liq_rsd).or.present(liq_rsd_1d)) then
          ! Set the relative standard deviation (RSD)
          call set_cld_field(cond_rsd, liq_rsd, liq_rsd_1d)
        else if (present(cloud_horizontal_rsd)) then
          cond_rsd = real(cloud_horizontal_rsd, RealK)
        else
          cond_rsd = 0.0_RealK
        end if

      case (ip_clcmp_cnv_water)
        if (present(liq_conv_mmr).or.present(liq_conv_mmr_1d)) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             liq_conv_mmr, liq_conv_mmr_1d)
        else if (l_split) then
          ! If splitting cloud then convective mmr may not have been provided
          cld%condensed_mix_ratio(:, :, i) = 0.0_RealK
        else
          cmessage = 'Convective liquid MMR not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
        if (l_combine .and. (present(liq_mmr).or.present(liq_mmr_1d))) then
          ! Add in the stratiform condensate if using combined cloud
          call set_cld_field(cond_mmr, liq_mmr, liq_mmr_1d)
          cld%condensed_mix_ratio(:, :, i) = &
            cld%condensed_mix_ratio(:, :, i) + cond_mmr
        end if
        if (present(liq_conv_rsd).or.present(liq_conv_rsd_1d)) then
          ! Set the relative standard deviation (RSD)
          call set_cld_field(cond_rsd, liq_conv_rsd, liq_conv_rsd_1d)
        else if (present(cloud_horizontal_rsd)) then
          cond_rsd = real(cloud_horizontal_rsd, RealK)
        else
          cond_rsd = 0.0_RealK
        end if

      end select

    else
      cmessage = 'Liquid cloud type not in spectral file.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if

  case (ip_phase_ice)
    if (i_param_type <= 0) then
      cld%i_condensed_param(i) = ip_ice_unparametrized
      cmessage = 'Prescribed ice cloud not yet implemented.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (i_param_type > spectrum%dim%nd_ice_type) then
      cmessage = 'Ice cloud type outside allowed range.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (spectrum%ice%l_ice_type(i_param_type)) then
      ! Take parametrisation from spectral file
      cld%i_condensed_param(i) = spectrum%ice%i_ice_parm(i_param_type)
      cld%condensed_n_phf(i) = spectrum%ice%n_phf(i_param_type)
      n_cloud_parameter = set_n_cloud_parameter( cld%i_condensed_param(i), &
        cld%type_condensed(i), cld%condensed_n_phf(i) )
      do j=1, spectrum%basic%n_band
        do k=1, n_cloud_parameter
          cld%condensed_param_list(k, i, j) &
            = spectrum%ice%parm_list(k, j, i_param_type)
        end do
      end do

      ! Assign ice mass mixing ratio, effective dimension,
      ! and relative standard deviation
      select case (cld%type_condensed(i))

      case (ip_clcmp_st_ice)
        if (present(ice_mmr).or.present(ice_mmr_1d)) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             ice_mmr, ice_mmr_1d)
          if (l_combine .and. &
              (present(ice_conv_mmr).or.present(ice_conv_mmr_1d))) then
            ! Add in the convective condensate if using combined cloud
            call set_cld_field(cond_mmr, ice_conv_mmr, ice_conv_mmr_1d)
            cld%condensed_mix_ratio(:, :, i) = &
              cld%condensed_mix_ratio(:, :, i) + cond_mmr
            ! Increment the convective condensate for the ratio calculation
            cld%c_ratio = cld%c_ratio + cond_mmr
          end if
        else
          cmessage = 'Ice MMR not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
        if (present(ice_rsd).or.present(ice_rsd_1d)) then
          ! Set the relative standard deviation (RSD)
          call set_cld_field(cond_rsd, ice_rsd, ice_rsd_1d)
        else if (present(cloud_horizontal_rsd)) then
          cond_rsd = real(cloud_horizontal_rsd, RealK)
        else
          cond_rsd = 0.0_RealK
        end if

      case (ip_clcmp_cnv_ice)
        if (present(ice_conv_mmr).or.present(ice_conv_mmr_1d)) then
          call set_cld_field(cld%condensed_mix_ratio(:, :, i), &
                             ice_conv_mmr, ice_conv_mmr_1d)
        else if (l_split) then
          ! If splitting cloud then convective mmr may not have been provided
          cld%condensed_mix_ratio(:, :, i) = 0.0_RealK
        else
          cmessage = 'Convective ice MMR not provided.'
          ierr=i_err_fatal
          call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
        end if
        if (l_combine .and. (present(ice_mmr).or.present(ice_mmr_1d))) then
          ! Add in the stratiform condensate if using combined cloud
          call set_cld_field(cond_mmr, ice_mmr, ice_mmr_1d)
          cld%condensed_mix_ratio(:, :, i) = &
            cld%condensed_mix_ratio(:, :, i) + cond_mmr
        end if
        if (present(ice_conv_rsd).or.present(ice_conv_rsd_1d)) then
          ! Set the relative standard deviation (RSD)
          call set_cld_field(cond_rsd, ice_conv_rsd, ice_conv_rsd_1d)
        else if (present(cloud_horizontal_rsd)) then
          cond_rsd = real(cloud_horizontal_rsd, RealK)
        else
          cond_rsd = 0.0_RealK
        end if

      end select

    else
      cmessage = 'Ice cloud type not in spectral file.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    end if
  end select

  ! Treatment of condensate inhomogeneity
  select case (control%i_inhom)
  case (ip_scaling)
    ! Scale the mixing ratio directly
    ! (equivalent to eqn 11a from Cairns et al. 2000).
    cld%condensed_mix_ratio(:, :, i) = &
      cld%condensed_mix_ratio(:, :, i) / (1.0_RealK + cond_rsd**2)
  case (ip_cairns)
    ! Set the relative variance as the square of the
    ! relative standard deviation. This will be used
    ! to renormalise the single scattering properties.
    cld%condensed_rel_var_dens(:, :, i) = cond_rsd**2
  case (ip_mcica, ip_tripleclouds_2019)
    ! Use the relative variation field to store the RSD
    cld%condensed_rel_var_dens(:, :, i) = cond_rsd
  end select

end do ! over condensed components


! Set the decorrelation scalings for cloud vertical overlap
if (present(cloud_vertical_decorr)) then
  cld%dp_corr_strat = real(cloud_vertical_decorr, RealK)
else
  cld%dp_corr_strat = 0.0_RealK
end if
if (present(conv_vertical_decorr)) then
  cld%dp_corr_conv  = real(conv_vertical_decorr, RealK)
else
  cld%dp_corr_conv = 0.0_RealK
end if


!------------------------------------------------------------------------------
! Set cloud amounts and convert mixing ratios to in-cloud values
!------------------------------------------------------------------------------

! Set cloud fractions
select case (control%i_cloud_representation)

case (ip_cloud_homogen, ip_cloud_combine_homogen)
  cld%n_cloud_type = 1
  do i = 1, cld%n_condensed
    cld%i_cloud_type(i) = ip_cloud_type_homogen
  end do
  if (present(cloud_frac).or.present(cloud_frac_1d)) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_homogen), &
                       cloud_frac, cloud_frac_1d)
  else
    cmessage = 'Cloud fraction not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
  if (l_combine .and. ( present(conv_frac).or.present(conv_frac_1d) )) then
    ! Add in the convective fraction if using combined cloud
    call set_cld_field(cld%c_cloud, conv_frac, conv_frac_1d)
    cld%frac_cloud(:, :, ip_cloud_type_homogen) = &
      cld%frac_cloud(:, :, ip_cloud_type_homogen) + cld%c_cloud
  end if

case (ip_cloud_ice_water, ip_cloud_combine_ice_water)
  cld%n_cloud_type = 2
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      cld%i_cloud_type(i) = ip_cloud_type_water
    case (ip_clcmp_st_ice)
      cld%i_cloud_type(i) = ip_cloud_type_ice
    end select
  end do
  if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
      (present(ice_frac).or.present(ice_frac_1d)) .and. &
      (present(cloud_frac).or.present(cloud_frac_1d))) then
    call set_cld_field(frac_liq, liq_frac, liq_frac_1d)
    call set_cld_field(frac_ice, ice_frac, ice_frac_1d)
    call set_cld_field(frac, cloud_frac, cloud_frac_1d)
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        if (frac_liq(l, k) + frac_ice(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          call split_mix_phase(ip_cloud_type_water, ip_cloud_type_ice)
          cld%frac_cloud(l, k, ip_cloud_type_water) = &
            frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_ice) = &
            frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_water) = 0.0_RealK
          cld%frac_cloud(l, k, ip_cloud_type_ice) = 0.0_RealK
        end if
      end do
    end do
  else if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
           (present(ice_frac).or.present(ice_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_water), &
                       liq_frac, liq_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_ice), &
                       ice_frac, ice_frac_1d)
  else
    cmessage = 'Liquid and ice cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
  if (l_combine .and. &
      (present(liq_conv_frac).or.present(liq_conv_frac_1d)) .and. &
      (present(ice_conv_frac).or.present(ice_conv_frac_1d))) then
    ! Add in the convective fraction if using combined cloud
    call set_cld_field(frac_liq, liq_conv_frac, liq_conv_frac_1d)
    call set_cld_field(frac_ice, ice_conv_frac, ice_conv_frac_1d)
    if (present(conv_frac).or.present(conv_frac_1d)) then
      call set_cld_field(frac, conv_frac, conv_frac_1d)
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          if (frac_liq(l, k) + frac_ice(l, k) > eps) then
            ! Split mixed phase fraction between ice and liquid
            call split_mix_phase(ip_cloud_type_water, ip_cloud_type_ice)
            frac_liq(l, k) = &
              frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
            frac_ice(l, k) = &
              frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          end if
        end do
      end do
    end if
    cld%frac_cloud(:, :, ip_cloud_type_water) = &
      cld%frac_cloud(:, :, ip_cloud_type_water) + frac_liq
    cld%frac_cloud(:, :, ip_cloud_type_ice) = &
      cld%frac_cloud(:, :, ip_cloud_type_ice) + frac_ice
    cld%c_cloud = frac_liq + frac_ice
  end if

case (ip_cloud_conv_strat, ip_cloud_split_homogen)
  cld%n_cloud_type = 2
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      cld%i_cloud_type(i) = ip_cloud_type_strat
    case (ip_clcmp_st_ice)
      cld%i_cloud_type(i) = ip_cloud_type_strat
    case (ip_clcmp_cnv_water)
      cld%i_cloud_type(i) = ip_cloud_type_conv
    case (ip_clcmp_cnv_ice)
      cld%i_cloud_type(i) = ip_cloud_type_conv
    end select
  end do
  if ((present(cloud_frac).or.present(cloud_frac_1d)) .and. &
      (present(conv_frac).or.present(conv_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_strat), &
                       cloud_frac, cloud_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_conv), &
                       conv_frac, conv_frac_1d)
  else if (l_split .and. &
           (present(cloud_frac).or.present(cloud_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_strat), &
                       cloud_frac, cloud_frac_1d)
    cld%frac_cloud(:, :, ip_cloud_type_conv) = 0.0_RealK
  else
    cmessage = 'Cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if

  if (l_split) then
    ! Combine the stratiform and convective fractions and then split
    ! into optically thick and thin regions
    frac = cld%frac_cloud(:, :, ip_cloud_type_strat) &
         + cld%frac_cloud(:, :, ip_cloud_type_conv)
    select case (control%i_inhom)
    case (ip_tripleclouds_2019)
      ! Use equations A1 and A2 from Hogan et al. (2019) to set the
      ! fractions and mixing ratios in each region of the split cloud.
      ! A single RSD field is needed so use the weighted mean of liq and ice.
      cond_rsd = 0.0_RealK
      cond_mmr = 0.0_RealK
      do i = 1, cld%n_condensed
        if (cld%i_cloud_type(i) == ip_cloud_type_strat) then
          cond_rsd = cond_rsd + cld%condensed_rel_var_dens(:, :, i) &
                              * cld%condensed_mix_ratio(:, :, i)
          cond_mmr = cond_mmr + cld%condensed_mix_ratio(:, :, i)
        end if
      end do
      cond_rsd = cond_rsd / max(cond_mmr, epsilon(cond_mmr))
      frac_m = 0.4_RealK/2.25_RealK
      frac_c = 0.525_RealK/2.25_RealK
      frac_thin = max(0.5_RealK, min(0.9_RealK, cond_rsd*frac_m + frac_c))
      cld%frac_cloud(:, :, ip_cloud_type_strat) = frac*frac_thin
      cld%frac_cloud(:, :, ip_cloud_type_conv) = frac*(1.0_RealK-frac_thin)
      ratio_thin = 0.025_RealK + 0.975_RealK &
        * exp(-cond_rsd-cond_rsd**2/2.0_RealK-cond_rsd**3/4.0_RealK)
      do i = 1, cld%n_condensed
        select case (cld%i_cloud_type(i))
        case (ip_cloud_type_strat)
          ! Optically thin region
          cld%condensed_mix_ratio(:, :, i) = &
            cld%condensed_mix_ratio(:, :, i) * ratio_thin * frac_thin
        case (ip_cloud_type_conv)
          ! Optically thick region
          cld%condensed_mix_ratio(:, :, i) = &
            cld%condensed_mix_ratio(:, :, i) * (1.0_RealK-ratio_thin*frac_thin)
        end select
      end do
    case default
      cld%frac_cloud(:, :, ip_cloud_type_strat) = 0.5_RealK*frac
      cld%frac_cloud(:, :, ip_cloud_type_conv) = 0.5_RealK*frac
      cld%condensed_mix_ratio(:, :, 1:cld%n_condensed) = &
        0.5_RealK*cld%condensed_mix_ratio(:, :, 1:cld%n_condensed)
    end select
  end if

case (ip_cloud_csiw, ip_cloud_split_ice_water)
  cld%n_cloud_type = 4
  do i = 1, cld%n_condensed
    select case (cld%type_condensed(i))
    case (ip_clcmp_st_water)
      cld%i_cloud_type(i) = ip_cloud_type_sw
    case (ip_clcmp_st_ice)
      cld%i_cloud_type(i) = ip_cloud_type_si
    case (ip_clcmp_cnv_water)
      cld%i_cloud_type(i) = ip_cloud_type_cw
    case (ip_clcmp_cnv_ice)
      cld%i_cloud_type(i) = ip_cloud_type_ci
    end select
  end do
  if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
      (present(ice_frac).or.present(ice_frac_1d)) .and. &
      (present(cloud_frac).or.present(cloud_frac_1d))) then
    call set_cld_field(frac_liq, liq_frac, liq_frac_1d)
    call set_cld_field(frac_ice, ice_frac, ice_frac_1d)
    call set_cld_field(frac, cloud_frac, cloud_frac_1d)
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        if (frac_liq(l, k) + frac_ice(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          call split_mix_phase(ip_cloud_type_sw, ip_cloud_type_si)
          cld%frac_cloud(l, k, ip_cloud_type_sw) = &
            frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_si) = &
            frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_sw) = 0.0_RealK
          cld%frac_cloud(l, k, ip_cloud_type_si) = 0.0_RealK
        end if
      end do
    end do
  else if ((present(liq_frac).or.present(liq_frac_1d)) .and. &
           (present(ice_frac).or.present(ice_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_sw), &
                       liq_frac, liq_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_si), &
                       ice_frac, ice_frac_1d)
  else
    cmessage = 'Liquid and ice cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if
  if ((present(liq_conv_frac).or.present(liq_conv_frac_1d)) .and. &
      (present(ice_conv_frac).or.present(ice_conv_frac_1d)) .and. &
      (present(conv_frac).or.present(conv_frac_1d))) then
    call set_cld_field(frac_liq, liq_conv_frac, liq_conv_frac_1d)
    call set_cld_field(frac_ice, ice_conv_frac, ice_conv_frac_1d)
    call set_cld_field(frac, conv_frac, conv_frac_1d)
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        if (frac_liq(l, k) + frac_ice(l, k) > eps) then
          ! Split mixed phase fraction between ice and liquid
          call split_mix_phase(ip_cloud_type_cw, ip_cloud_type_ci)
          cld%frac_cloud(l, k, ip_cloud_type_cw) = &
            frac(l, k)*frac_liq(l, k) / (frac_liq(l, k)+frac_ice(l, k))
          cld%frac_cloud(l, k, ip_cloud_type_ci) = &
            frac(l, k)*frac_ice(l, k) / (frac_liq(l, k)+frac_ice(l, k))
        else
          cld%frac_cloud(l, k, ip_cloud_type_cw) = 0.0_RealK
          cld%frac_cloud(l, k, ip_cloud_type_ci) = 0.0_RealK
        end if
      end do
    end do
  else if ((present(liq_conv_frac).or.present(liq_conv_frac_1d)) .and. &
           (present(ice_conv_frac).or.present(ice_conv_frac_1d))) then
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_cw), &
                       liq_conv_frac, liq_conv_frac_1d)
    call set_cld_field(cld%frac_cloud(:, :, ip_cloud_type_ci), &
                       ice_conv_frac, ice_conv_frac_1d)
  else if (l_split) then
    cld%frac_cloud(:, :, ip_cloud_type_cw) = 0.0_RealK
    cld%frac_cloud(:, :, ip_cloud_type_ci) = 0.0_RealK
  else
    cmessage = 'Liquid and ice convective cloud fractions not provided.'
    ierr=i_err_fatal
    call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
  end if

  if (l_split) then
    ! Combine the stratiform and convective fractions and then split
    ! into optically thick and thin regions
    select case (control%i_inhom)
    case (ip_tripleclouds_2019)
      ! Use equations A1 and A2 from Hogan et al. (2019) to set the
      ! fractions and mixing ratios in each region of the split cloud.
      do i_phase = 1, 2
        ! Liquid and ice fractions are split according to their RSD.
        select case (i_phase)
        case (ip_phase_water)
          i_thin = ip_cloud_type_sw
          i_thick = ip_cloud_type_cw
        case (ip_phase_ice)
          i_thin = ip_cloud_type_si
          i_thick = ip_cloud_type_ci
        end select
        ! Select the RSD for this phase (or zero if this phase is turned off)
        cond_rsd = 0.0_RealK
        do i = 1, cld%n_condensed
          if (cld%i_cloud_type(i) == i_thin) then
            cond_rsd = cld%condensed_rel_var_dens(:, :, i)
          end if
        end do
        ! First combine the strat and conv cloud fractions for this phase
        frac = cld%frac_cloud(:, :, i_thin) + cld%frac_cloud(:, :, i_thick)
        frac_m = 0.4_RealK/2.25_RealK
        frac_c = 0.525_RealK/2.25_RealK
        ! Optically thin fraction set between 0.5 and 0.9 of total
        frac_thin = max(0.5_RealK, min(0.9_RealK, cond_rsd*frac_m + frac_c))
        cld%frac_cloud(:, :, i_thin) = frac*frac_thin
        ! Optically thick fraction set between 0.1 and 0.5 of total
        cld%frac_cloud(:, :, i_thick) = frac*(1.0_RealK-frac_thin)
        ! Condensate in optically thin region reduced according to equations
        ! A1 and A2 in Hogan et al 2019. The condensate in the optically 
        ! thick region is adjusted to preserve the correct gridbox mean.
        ratio_thin = 0.025_RealK + 0.975_RealK &
          * exp(-cond_rsd-cond_rsd**2/2.0_RealK-cond_rsd**3/4.0_RealK)
        do i = 1, cld%n_condensed
          if (cld%i_cloud_type(i) == i_thin) then
            ! Optically thin region
            cld%condensed_mix_ratio(:, :, i) = &
              cld%condensed_mix_ratio(:, :, i) * ratio_thin * frac_thin
          else if (cld%i_cloud_type(i) == i_thick) then
            ! Optically thick region
            cld%condensed_mix_ratio(:, :, i) = &
              cld%condensed_mix_ratio(:, :, i)*(1.0_RealK-ratio_thin*frac_thin)
          end if
        end do
      end do
    case default
      ! The default case splits the cloud into two equal regions.
      frac_liq = cld%frac_cloud(:, :, ip_cloud_type_sw) &
               + cld%frac_cloud(:, :, ip_cloud_type_cw)
      frac_ice = cld%frac_cloud(:, :, ip_cloud_type_si) &
               + cld%frac_cloud(:, :, ip_cloud_type_ci)
      cld%frac_cloud(:, :, ip_cloud_type_sw) = 0.5_RealK*frac_liq
      cld%frac_cloud(:, :, ip_cloud_type_cw) = 0.5_RealK*frac_liq
      cld%frac_cloud(:, :, ip_cloud_type_si) = 0.5_RealK*frac_ice
      cld%frac_cloud(:, :, ip_cloud_type_ci) = 0.5_RealK*frac_ice
      cld%condensed_mix_ratio(:, :, 1:cld%n_condensed) = &
        0.5_RealK*cld%condensed_mix_ratio(:, :, 1:cld%n_condensed)
    end select
  end if

end select

if (l_combine) then
  ! Set the convective cloud contribution to combined cloud
  do k = dimen%id_cloud_top, atm%n_layer
    do l = 1, atm%n_profile
      ! Temporarily set c_ratio to the in-cloud convective condensate mmr
      cld%c_ratio(l, k) = cld%c_ratio(l, k) / max(cld%c_cloud(l, k), eps)
    end do
  end do
  cond_mmr = 0.0_RealK
  do i = 1, cld%n_condensed
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        cond_mmr(l, k) = cond_mmr(l, k) + cld%condensed_mix_ratio(l, k, i)
      end do
    end do
  end do
end if

! Convert mass mixing ratios to in-cloud values
do i = 1, cld%n_condensed
  do k = dimen%id_cloud_top, atm%n_layer
    do l = 1, atm%n_profile
      if (cld%frac_cloud(l, k, cld%i_cloud_type(i)) > eps) then
        cld%condensed_mix_ratio(l, k, i) = cld%condensed_mix_ratio(l, k, i) &
          / cld%frac_cloud(l, k, cld%i_cloud_type(i))
      else
        cld%condensed_mix_ratio(l, k, i) = 0.0_RealK
      end if
    end do
  end do
end do

! Normalise the cloud fractions
do k = dimen%id_cloud_top, atm%n_layer
  do l = 1, atm%n_profile
    cld%w_cloud(l, k) = sum(cld%frac_cloud(l, k, 1:cld%n_cloud_type))
    if (cld%w_cloud(l, k) > min_cloud_fraction) then
      do j=1, cld%n_cloud_type
        cld%frac_cloud(l, k, j) = cld%frac_cloud(l, k, j) / cld%w_cloud(l, k)
      end do
      if (l_combine) then
        ! Set c_ratio to the in-cloud convective condensate MMR
        ! divided by the total in-cloud condensate MMR
        cld%c_ratio(l, k) = &
          cld%c_ratio(l, k) * cld%w_cloud(l, k) / max(cond_mmr(l, k), eps)
      end if
    else
      cld%w_cloud(l, k) = 0.0_RealK
      cld%frac_cloud(l, k, 1:cld%n_cloud_type) = 0.0_RealK
      if (l_combine) then
        cld%c_cloud(l, k) = 0.0_RealK
        cld%c_ratio(l, k) = 0.0_RealK
      end if
    end if
    if (cld%w_cloud(l, k) > 1.0_RealK + min_cloud_fraction) then
      cmessage = 'Cloud fraction greater than 1.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)
    else if (cld%w_cloud(l, k) > 1.0_RealK) then
      cld%w_cloud(l, k) = 1.0_RealK
    end if
  end do
end do


contains

  subroutine set_cld_field(out_field, full_field, oned_field)
    implicit none

    real(RealK), intent(inout) :: out_field(:, dimen%id_cloud_top:)
!     Output field
    real(RealExt), intent(in), optional :: full_field(:, :)
!     Full field variable
    real(RealExt), intent(in), optional :: oned_field(:)
!     One-dimensional variable

    if (present(full_field)) then
      if (l_last) then
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            out_field(l, k) = max(real(full_field(kk, list(l)),RealK),0.0_RealK)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            out_field(l, k) = max(real(full_field(list(l), kk),RealK),0.0_RealK)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            ll = stride_layer*(list(l)-1) + kk
            out_field(l, k) = max(real(oned_field(ll), RealK), 0.0_RealK)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            ll = atm%n_profile*(kk-1) + list(l)
            out_field(l, k) = max(real(oned_field(ll), RealK), 0.0_RealK)
          end do
        end do
      end if
    else
      cmessage = 'The required cloud fields have not been provided.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)      
    end if
  end subroutine set_cld_field


  subroutine split_mix_phase(i_cloud_type_water, i_cloud_type_ice)
    implicit none

    integer, intent(in) :: i_cloud_type_water, i_cloud_type_ice
    real(RealK) :: frac_mix, liq_mix_mmr, ice_mix_mmr

    ! Split mixed phase fraction between ice and liquid
    ! using the ratio of their in-cloud MMRs

    frac_mix = frac_liq(l, k) + frac_ice(l, k) - frac(l, k)
    if (frac_mix > eps .and. &
        frac_liq(l, k) > eps .and. frac_ice(l, k) > eps) then
      liq_mix_mmr = 0.0_RealK
      ice_mix_mmr = 0.0_RealK
      do i = 1, cld%n_condensed
        if (cld%i_cloud_type(i) == i_cloud_type_water) then
          liq_mix_mmr = liq_mix_mmr &
            + cld%condensed_mix_ratio(l, k, i) / frac_liq(l, k)
        else if (cld%i_cloud_type(i) == i_cloud_type_ice) then
          ice_mix_mmr = ice_mix_mmr &
            + cld%condensed_mix_ratio(l, k, i) / frac_ice(l, k)
        end if
      end do
      if (liq_mix_mmr + ice_mix_mmr > eps) then
        frac_liq(l, k) = max(frac_liq(l, k) - frac_mix &
          * ice_mix_mmr / (liq_mix_mmr + ice_mix_mmr), 0.0_RealK)
        frac_ice(l, k) = max(frac_ice(l, k) - frac_mix &
          * liq_mix_mmr / (liq_mix_mmr + ice_mix_mmr), 0.0_RealK)
      end if
    end if
  end subroutine split_mix_phase

end subroutine set_cld
end module socrates_set_cld
