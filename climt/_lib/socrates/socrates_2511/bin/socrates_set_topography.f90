! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Set the slope and horizon angles for topographic corrections

module socrates_set_topography

implicit none
private
public :: set_topography

character(len=*), parameter :: ModuleName='SOCRATES_SET_TOPOGRAPHY'

contains

! For further details see: Manners et al (2012), Radiative transfer over
! resolved topographic features for high-resolution weather prediction,
! https://doi.org/10.1002/qj.956

subroutine set_topography(n_profile, pi, n_horiz_ang, n_horiz_layer, &
  grad_x, grad_y, grad_aspect, horizon_angle, horizon_aspect, &
  slope_angle, slope_aspect, skyview)

use realtype_rd, only: RealExt
use errormessagelength_mod, only: errormessagelength
use ereport_mod, only: ereport
use rad_pcf, only: i_normal, i_err_fatal
use missing_data_mod, only: rmdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Dimensions
integer, intent(in) :: n_profile
real(RealExt), intent(in) :: pi
integer, intent(in), optional :: n_horiz_ang, n_horiz_layer

! Input ancillary data
real(RealExt), intent(in), optional :: &
  grad_x(n_profile), grad_y(n_profile), grad_aspect(n_profile)
real(RealExt), intent(inout), optional :: &
  horizon_angle(:), horizon_aspect(:)

! Output variabes
real(RealExt), intent(out), optional :: &
  slope_angle(n_profile), slope_aspect(n_profile), skyview(n_profile)

! Local variables
integer :: k, l, ll, n_ang, n_layer
real(RealExt) :: work(n_profile)
real(RealExt), parameter :: eps = epsilon(1.0_RealExt)

!integer :: ierr = i_normal
!character(len=errormessagelength) :: cmessage
character(len=*), parameter :: RoutineName='SET_TOPOGRAPHY'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Slope Angle
if (present(slope_angle)) then
  if (present(grad_x) .and. present(grad_y)) then
    ! Find slope angles and aspects from x & y gradients:
    slope_angle = atan( (grad_x**2 + grad_y**2)**0.5_RealExt )
  else
    slope_angle = 0.0_RealExt
  end if
end if

! Slope Aspect
if (present(slope_aspect)) then
  if (present(grad_x) .and. present(grad_y)) then
    work = grad_x
    where (abs(grad_x) < eps) work = sign(eps, grad_x)
    slope_aspect = pi - atan(grad_y/work) + sign(pi/2.0_RealExt, work)
  else
    slope_aspect = 0.0_RealExt
  end if
  if (present(grad_aspect)) then
    ! Add bearing of grid north so aspects are relative to true North:
    slope_aspect = modulo(slope_aspect + grad_aspect, pi*2.0_RealExt)
  end if
end if

! Horizon Angles
if (present(horizon_angle) .and. present(horizon_aspect)) then
  if (present(n_horiz_ang)) then
    n_ang = n_horiz_ang
  else
    ! Default to 16 horizon angles
    n_ang = 16
  end if
  if (present(n_horiz_layer)) then
    n_layer = n_horiz_layer
  else
    ! Surface horizon angles only by default
    n_layer = 1
  end if

  if (present(slope_angle) .and. present(slope_aspect)) then
    ! The surface horizon angle must include self-shadowing from the slope
    ! of the grid-box in each direction.
    ll = 0
    do l = 1, n_profile
      do k = 1, n_ang
        ll = ll + 1
        horizon_angle(ll) = min( horizon_angle(ll), atan2( 1.0_RealExt, &
          -tan(slope_angle(l)) * cos(slope_aspect(l) - horizon_aspect(ll)) ) )
        ll = ll + n_layer-1
      end do
    end do
  end if
end if

if (present(skyview)) then
  ! Note calculation of the horizon angle in the inclined reference frame
  ! is incorrectly formulated in Manners et al (2012), eqn 13. A corrected
  ! formulation is used here.
  ll = 0
  do l = 1, n_profile
    skyview(l) = 0.0_RealExt
    do k = 1, n_ang
      ll = ll + 1
      skyview(l) = skyview(l) + 1.0_RealExt &
        - ( cos(slope_angle(l)) * cos(horizon_angle(ll)) &
          + sin(slope_angle(l)) * sin(horizon_angle(ll)) &
          * cos(slope_aspect(l) - horizon_aspect(ll)) )**2
    end do
    if (n_ang > 0) then
      skyview(l) = skyview(l) / ( real(n_ang, RealExt) * cos(slope_angle(l)) )
    end if
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine set_topography
end module socrates_set_topography
