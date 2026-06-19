! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the droplet and ice crystal dimension in the Socrates cloud type
!
!------------------------------------------------------------------------------
module socrates_set_cld_dim
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_SET_CLD_DIM'
contains

subroutine set_cld_dim(cld, control, dimen, spectrum, atm, &
  profile_list, n_layer_stride, &
  liq_nc, ice_nc, liq_conv_nc, ice_conv_nc, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim, &
  liq_nc_1d, ice_nc_1d, liq_conv_nc_1d, ice_conv_nc_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d, &
  liq_dim_constant, liq_conv_dim_constant, &
  liq_dim_aparam, liq_dim_bparam, &
  l_invert, l_profile_last, l_debug, i_profile_debug)

use def_cld,      only: StrCld
use def_control,  only: StrCtrl
use def_dimen,    only: StrDim
use def_spectrum, only: StrSpecData
use def_atm,      only: StrAtm
use realtype_rd,  only: RealK, RealExt
use rad_pcf,      only: &
  ip_cloud_combine_homogen, ip_cloud_combine_ice_water, &
  ip_cloud_split_homogen, ip_cloud_split_ice_water, &
  ip_clcmp_st_water, ip_clcmp_st_ice, ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
  ip_re_external, ip_re_constant, &
  i_normal, i_err_fatal
use ereport_mod,  only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none


! Cloud properties:
type(StrCld),      intent(inout) :: cld

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
  liq_nc, ice_nc, liq_conv_nc, ice_conv_nc, &
  liq_dim, ice_dim, liq_conv_dim, ice_conv_dim
real(RealExt), intent(in), dimension(:), optional :: &
  liq_nc_1d, ice_nc_1d, liq_conv_nc_1d, ice_conv_nc_1d, &
  liq_dim_1d, ice_dim_1d, liq_conv_dim_1d, ice_conv_dim_1d
real(RealExt), intent(in), optional :: &
  liq_dim_constant, liq_conv_dim_constant
!   Liquid and ice number concentration and effective dimensions

real(RealExt), intent(in), optional :: liq_dim_aparam
real(RealExt), intent(in), optional :: liq_dim_bparam
!   Parameters for calculating cloud droplet effective radius

logical, intent(in), optional :: l_invert
!   Flag to invert fields in the vertical
logical, intent(in), optional :: l_profile_last
!   Loop over profiles is last in input fields

logical, intent(in), optional :: l_debug
integer, intent(in), optional :: i_profile_debug
!   Options for outputting debugging information


! Local variables
integer :: i, k, kk, l, ll, list(dimen%nd_profile)
!   Loop variables
integer :: i_param_type, layer_offset, stride_layer
!   Working variable
logical :: l_last
!   Local flag to loop over profiles last in 1d fields
logical :: l_combine
!   Combine stratiform and convective cloud
logical :: l_split
!   Split cloud into optically thick and thin regions

integer                      :: ierr = i_normal
character (len=*), parameter :: RoutineName = 'SET_CLD_DIM'
character (len=errormessagelength) :: cmessage


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


select case (control%i_cloud_representation)
case (ip_cloud_combine_homogen, ip_cloud_combine_ice_water)
  ! Combine stratiform and convective cloud into a single region
  l_combine = .true.
  l_split = .false.
case (ip_cloud_split_homogen, ip_cloud_split_ice_water)
  ! Stratiform and convective cloud is combined and the
  ! total cloud is split into optically thick and thin regions
  l_combine = .true.
  l_split = .true.
case default
  l_combine = .false.
  l_split = .false.
end select


do i=1, cld%n_condensed
  if (cld%type_condensed(i) == ip_clcmp_st_water .or. &
     (cld%type_condensed(i) == ip_clcmp_cnv_water .and. l_split)) then
    if (control%i_drop_re == ip_re_external .and. &
        (present(liq_dim).or.present(liq_dim_1d))) then
      ! Set the droplet effective radius
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         liq_dim, liq_dim_1d)
    else if (control%i_drop_re == ip_re_constant .and. &
             present(liq_dim_constant)) then
      ! Set the droplet effective radius to a constant
      cld%condensed_dim_char(:, :, i) = real(liq_dim_constant, RealK)
    else
      call set_liq_dim(liq_nc, liq_nc_1d)
    end if
    i_param_type = control%i_st_water

  else if (cld%type_condensed(i) == ip_clcmp_cnv_water) then
    if (control%i_drop_re == ip_re_external .and. &
        (present(liq_conv_dim).or.present(liq_conv_dim_1d))) then
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         liq_conv_dim, liq_conv_dim_1d)
    else if (control%i_drop_re == ip_re_constant .and. &
             present(liq_conv_dim_constant)) then
      cld%condensed_dim_char(:, :, i) = real(liq_conv_dim_constant, RealK)
    else
      call set_liq_dim(liq_conv_nc, liq_conv_nc_1d)
    end if
    i_param_type = control%i_cnv_water

  else if (cld%type_condensed(i) == ip_clcmp_st_ice .or. &
          (cld%type_condensed(i) == ip_clcmp_cnv_ice .and. l_split)) then
    if (present(ice_dim).or.present(ice_dim_1d)) then
      ! Set the ice-crystal effective dimension
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         ice_dim, ice_dim_1d)
    else if (l_combine) then
      call set_ice_dim(ice_nc, ice_nc_1d, ice_conv_nc, ice_conv_nc_1d)
    else
      call set_ice_dim(ice_nc, ice_nc_1d)
    end if
    i_param_type = control%i_st_ice

  else if (cld%type_condensed(i) == ip_clcmp_cnv_ice) then
    if (present(ice_conv_dim).or.present(ice_conv_dim_1d)) then
      call set_cld_field(cld%condensed_dim_char(:, :, i), &
                         ice_conv_dim, ice_conv_dim_1d)
    else
      call set_ice_dim(ice_conv_nc, ice_conv_nc_1d)
    end if
    i_param_type = control%i_cnv_ice

  end if

  select case (cld%type_condensed(i))
  case (ip_clcmp_st_water, ip_clcmp_cnv_water)
    ! Constrain effective radius to be within parametrisation bounds
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        cld%condensed_dim_char(l, k, i) = min( max( &
          cld%condensed_dim_char(l, k, i), &
          spectrum%drop%parm_min_dim(i_param_type) ), &
          spectrum%drop%parm_max_dim(i_param_type) )
      end do
    end do
  case (ip_clcmp_st_ice, ip_clcmp_cnv_ice)
    ! Constrain effective dimension to be within parametrisation bounds
    do k = dimen%id_cloud_top, atm%n_layer
      do l = 1, atm%n_profile
        cld%condensed_dim_char(l, k, i) = min( max( &
          cld%condensed_dim_char(l, k, i), &
          spectrum%ice%parm_min_dim(i_param_type) ), &
          spectrum%ice%parm_max_dim(i_param_type) )
      end do
    end do
  end select

end do ! over condensed components



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
            out_field(l, k) = real(full_field(kk, list(l)), RealK)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            out_field(l, k) = real(full_field(list(l), kk), RealK)
          end do
        end do
      end if
    else if (present(oned_field)) then
      if (l_last) then
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            ll = stride_layer*(list(l)-1) + kk
            out_field(l, k) = real(oned_field(ll), RealK)
          end do
        end do
      else
        do k = dimen%id_cloud_top, atm%n_layer
          kk = abs(layer_offset-k)
          do l=1, atm%n_profile
            ll = atm%n_profile*(kk-1) + list(l)
            out_field(l, k) = real(oned_field(ll), RealK)
          end do
        end do
      end if
    else
      cmessage = 'The required cloud fields have not been provided.'
      ierr=i_err_fatal
      call ereport(ModuleName//':'//RoutineName, ierr, cmessage)      
    end if
  end subroutine set_cld_field


  subroutine set_liq_dim(full_nc, oned_nc)

    use rad_pcf, only: ip_re_liu
    use rad_ccf, only: pi, rho_water
    implicit none

    real(RealExt), intent(in), optional :: full_nc(:, :)
!     Full field cloud droplet number concentration
    real(RealExt), intent(in), optional :: oned_nc(:)
!     One-dimensional cloud droplet number concentration

    real(RealK) :: cdnc(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer)
!     Working cloud droplet number concentration
    real(RealK), parameter :: eps = epsilon(1.0_RealK)

    ! Parameters for Liu spectral dispersion
    real(RealK) :: aparam
    real(RealK) :: bparam
    real(RealK) :: beta

    if (present(liq_dim_aparam)) then
      aparam = real(liq_dim_aparam, RealK)
    else
      aparam = 0.07_RealK
    end if

    if (present(liq_dim_bparam)) then
      bparam = real(liq_dim_bparam, RealK)
    else
      bparam = -0.14_RealK
    end if

    select case (control%i_drop_re)
    case (ip_re_liu)
      call set_cld_field(cdnc, full_nc, oned_nc)
      ! Apply Liu spectral dispersion
      do k = dimen%id_cloud_top, atm%n_layer
        do l=1, atm%n_profile
          beta = aparam * ((MAX(eps, &
            cld%condensed_mix_ratio(l, k, i)) * atm%density(l, k) &
            * 1.0e-3_RealK / (cdnc(l, k) * 1e-06_RealK))**(bparam))
          cld%condensed_dim_char(l, k, i) = MAX(0.0_RealK, &
            3.0_RealK * cld%condensed_mix_ratio(l, k, i) * atm%density(l, k) &
            / (4.0_RealK * pi * rho_water * (beta**(-3)) * cdnc(l, k))) &
            **(1.0_RealK/3.0_RealK)
        end do
      end do
    case default
      ! A default value of 7-microns is assumed.
      cld%condensed_dim_char(:, :, i) = 7.0e-6_RealK
    end select

  end subroutine set_liq_dim


  subroutine set_ice_dim(full_nc, oned_nc, full_conv_nc, oned_conv_nc)

    use rad_pcf, only: &
      ip_ice_adt, ip_ice_agg_de, ip_ice_agg_de_sun, ip_ice_baran, &
      ip_ice_pade_2_phf
    use rad_ccf, only: pi
    implicit none

    real(RealExt), intent(in), optional :: full_nc(:, :), full_conv_nc(:, :)
!     Full field cloud ice number concentration
    real(RealExt), intent(in), optional :: oned_nc(:), oned_conv_nc(:)
!     One-dimensional cloud ice number concentration

    real(RealK) :: cinc(dimen%nd_profile, dimen%id_cloud_top:dimen%nd_layer)
    real(RealK) :: cinc_conv(dimen%nd_profile,dimen%id_cloud_top:dimen%nd_layer)
    real(RealK) :: cinc_incloud
!     Working cloud ice number concentration
    real(RealK), parameter :: eps = epsilon(1.0_RealK)
    real(RealK), parameter :: min_cinc = tiny(1.0_RealK)
!     Tolerance to prevent divide by zero

    ! Parameters for the aggregate parametrization.
    real (RealK), parameter :: a0_agg_cold = 7.5094588e-04_RealK
    real (RealK), parameter :: b0_agg_cold = 5.0830326e-07_RealK
    real (RealK), parameter :: a0_agg_warm = 1.3505403e-04_RealK
    real (RealK), parameter :: b0_agg_warm = 2.6517429e-05_RealK
    real (RealK), parameter :: t_switch    = 216.208_RealK
    real (RealK), parameter :: t0_agg      = 279.5_RealK
    real (RealK), parameter :: s0_agg      = 0.05_RealK
    real (RealK), parameter :: dge2de      = &
      (3.0_RealK/2.0_RealK)*(3.0_RealK/(2.0_RealK*sqrt(3.0_RealK)))

    ! Parameters for Baran's diagnostic parametrisation of effective size.
    ! m and n are the slope and intercept (positive sign) in units of metres.
    real (RealK), parameter :: m_ice_baran = 1.868e-6_RealK
    real (RealK), parameter :: n_ice_baran = 353.613e-6_RealK

    ! Parameters for the calculation of equivalent spherical radius
    real (RealK), parameter :: rho_ice = 9.17e+02_RealK

    select case (cld%i_condensed_param(i))

    case (ip_ice_adt)
      ! This parametrization is based on the mean maximum
      ! dimension of the crystal, determined as a function of
      ! the local temperature. The size is limited to its value
      ! at the freezing level.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          cld%condensed_dim_char(l, k, i) = min( 7.198755e-04_RealK, &
            exp(5.522e-02_RealK * (atm%t(l, k) - 2.7965e+02_RealK)) &
            / 9.702e+02_RealK )
        end do
      end do

    case (ip_ice_agg_de, ip_ice_agg_de_sun)
      ! Aggregate parametrization based on effective dimension.
      ! The fit provided here is based on Stephan Havemann's fit of
      ! Dge with temperature, consistent with David Mitchell's treatment
      ! of the variation of the size distribution with temperature. The
      ! parametrization of the optical properties is based on De
      ! (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
      ! (=(2*SQRT(3)/3)*volume/projected area), which explains the
      ! conversion factor. The fit to Dge is in two sections, because
      ! Mitchell's relationship predicts a cusp at 216.208 K. Limits
      ! of 8 and 124 microns are imposed on Dge: these are based on this
      ! relationship and should be reviewed if it is changed. Note also
      ! that the relationship given here is for polycrystals only.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          ! Preliminary calculation of Dge.
          if (atm%t(l, k) < t_switch) then
            cld%condensed_dim_char(l, k, i) = &
              a0_agg_cold*exp(s0_agg*(atm%t(l, k)-t0_agg))+b0_agg_cold
          else
            cld%condensed_dim_char(l, k, i) = &
              a0_agg_warm*exp(s0_agg*(atm%t(l, k)-t0_agg))+b0_agg_warm
          end if
          ! Limit and convert to De.
          cld%condensed_dim_char(l, k, i) = dge2de * min( 1.24e-04_RealK, &
            max(8.0e-06_RealK, cld%condensed_dim_char(l, k, i)) )
        end do
      end do

    case (ip_ice_baran)
      ! Diagnostic De=De(T) parameterisation developed by A. Baran.
      ! Fit to data based on Field et al., JAS, 2007, 10.1175/2007JAS2344.1.
      ! The weighting of the ensemble model follows that assumed in
      ! experiment 4 of Baran et al. (2014 and 2016). The weights applied
      ! at each size bin are 0.50,0.20,0.30, 0.0, 0.0, 0.0 (that is
      ! hexagons+rosettes account for 70% of the mass at each size bin).
      ! De is defined as 1.50*(mass of PSD/area of PSD).
      ! The size is limited to its value at the freezing level.
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          cld%condensed_dim_char(l, k, i) &
            = m_ice_baran * atm%t(l, k) - n_ice_baran
        end do
      end do

    case (ip_ice_pade_2_phf)
      ! Pade fits based on the equivalent spherical radius.
      call set_cld_field(cinc, full_nc, oned_nc)
      if (present(full_conv_nc).or.present(oned_conv_nc)) then
        call set_cld_field(cinc_conv, full_conv_nc, oned_conv_nc)
        do k = dimen%id_cloud_top, atm%n_layer
          do l = 1, atm%n_profile
            ! Combined cloud requires total ice number  
            cinc(l, k) = cinc(l, k) + cinc_conv(l, k)
          end do
        end do
      end if
      do k = dimen%id_cloud_top, atm%n_layer
        do l = 1, atm%n_profile
          ! Convert grid-box mean values to in-cloud values
          cinc_incloud = max( min_cinc, cinc(l, k) / max( eps, &
            cld%w_cloud(l, k) * cld%frac_cloud(l, k, cld%i_cloud_type(i)) ) )
          ! Calculate equivalent spherical radius of ice crystals
          cld%condensed_dim_char(l, k, i) = max( 0.0_RealK, &
            3.0_RealK * cld%condensed_mix_ratio(l, k, i) * atm%density(l, k) &
            / (4.0_RealK * pi * rho_ice * cinc_incloud) ) &
            **(1.0_RealK/3.0_RealK)
        end do
      end do

    case default
      ! A default value of 30-microns is assumed.
      cld%condensed_dim_char(:, :, i) = 30.0e-6_RealK

    end select

  end subroutine set_ice_dim

end subroutine set_cld_dim
end module socrates_set_cld_dim
