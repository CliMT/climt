! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Correction for the angle of solar incidence on sloping terrain.

module solinc_mod

implicit none
character(len=*), parameter, private :: ModuleName = 'SOLINC_MOD'
contains

! Calculate the orographic correction to be applied to the direct SW flux
! at the surface using the angle of solar incidence on sloping terrain.
! Shading from surrounding terrain may optionally be included.
! For further details see: Manners et al (2012), Radiative transfer over
! resolved topographic features for high-resolution weather prediction,
! https://doi.org/10.1002/qj.956
subroutine solinc(n_profile, cos_zenith_angle, sol_azimuth, &
  slope_aspect, slope_angle, orog_corr, l_shading, &
  n_horiz_ang, horiz_aspect, horiz_angle, cosz_beg, cosz_end)

use realtype_rd, only: RealK
use rad_ccf, only: pi
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

implicit none


integer, intent(in) :: n_profile
!   Number of columns

real(RealK), intent(in) :: cos_zenith_angle(n_profile)
!   Cosine of the solar zenith angle
real(RealK), intent(in) :: sol_azimuth(n_profile)
!   Solar azimuth angle in radians clockwise from grid north
real(RealK), intent(in) :: slope_aspect(n_profile)
!   Azimuthal angle of slope normal with respect to North
real(RealK), intent(in) :: slope_angle(n_profile)
!   Zenith angle of slope normal

real(RealK), intent(out) :: orog_corr(n_profile)
!   Orographic correction factor

logical, intent(in), optional :: l_shading
!   Include effect of shading from surrounding terrain
integer, intent(in), optional :: n_horiz_ang
!   Number of horizon angles
real(RealK), intent(in), optional :: horiz_aspect(:, :)
!   Azimuthal angle of horizon in radians clockwise from true north
!   (n_horiz_ang, n_profile)
real(RealK), intent(in), optional :: horiz_angle(:, :, :)
!   Zenith angle of horizon (obscuring terrain) in radians
!   (n_horiz_layer, n_horiz_ang, n_profile)
real(RealK), intent(in), optional :: cosz_beg(n_profile), cosz_end(n_profile)
!   Cosine of the solar zenith angle at the beginning and end of timestep

! Local variables
integer :: l
integer :: max_sub(n_profile)
integer :: min_sub(n_profile)

real(RealK), allocatable :: rel_horiz_aspect(:, :)
real(RealK) :: cosz_hor(n_profile)
real(RealK) :: cosz_lit(n_profile)
real(RealK) :: lit_frac(n_profile)

integer(KIND=jpim), parameter :: zhook_in  = 0
integer(KIND=jpim), parameter :: zhook_out = 1
real(KIND=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SOLINC'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! orog_corr is equal to the cosine of the angle between the incoming
! solar insolation and the normal to the mean slope, divided by
! (cosine of the solar zenith angle x cosine of the slope angle).

where (cos_zenith_angle /= 0.0_RealK .and. slope_angle /= 0.0_RealK)

  orog_corr = 1.0_RealK + tan( acos(cos_zenith_angle) ) * &
     tan(slope_angle) * cos( sol_azimuth - slope_aspect )

  orog_corr = max(orog_corr, epsilon(orog_corr))

elsewhere
  orog_corr = 1.0_RealK
end where

if (present(l_shading)) then
if (l_shading) then
  ! Calculate cosine of terrain horizon angle
  allocate(rel_horiz_aspect(n_horiz_ang, n_profile))
  rel_horiz_aspect = modulo( horiz_aspect - &
    spread(sol_azimuth, 1, n_horiz_ang), pi*2.0_RealK)
  !PGI needs the mask for unknown reasons (see ticket um:#2159). The result of
  !the modulo above should be +ve definite anyway since pi*2.0 is positive.
  min_sub=minloc(rel_horiz_aspect, dim=1, mask=rel_horiz_aspect >= 0.0_RealK)
  max_sub=maxloc(rel_horiz_aspect, dim=1, mask=rel_horiz_aspect >= 0.0_RealK)
  do l=1, n_profile
    cosz_hor(l) = cos( (horiz_angle(1, min_sub(l), l) &
      * (pi*2.0_RealK - rel_horiz_aspect(max_sub(l), l)) &
      + horiz_angle(1, max_sub(l), l) * rel_horiz_aspect(min_sub(l), l)) &
      / (pi*2.0_RealK - rel_horiz_aspect(max_sub(l), l) &
      + rel_horiz_aspect(min_sub(l), l)) )
  end do
  deallocate(rel_horiz_aspect)

  ! Adjust orographic correction for times when sun is below horizon
  ! (assumes change of zenith angle is roughly linear over the timestep)
  where (cos_zenith_angle /= 0.0_RealK)
    where (cosz_hor > cosz_beg .and. cosz_hor > cosz_end)
      ! Sun is below the horizon
      orog_corr = epsilon(orog_corr)
    elsewhere (cosz_hor > cosz_beg .and. cosz_hor < cosz_end)
      ! Sun rises during timestep
      lit_frac = (cosz_end - cosz_hor) / (cosz_end - cosz_beg)
      cosz_lit = (cosz_end + cosz_hor) * 0.5_RealK
      orog_corr = lit_frac * ( cosz_lit + sin(acos(cosz_lit)) * &
        tan(slope_angle) * cos(sol_azimuth - slope_aspect) ) / cos_zenith_angle
      orog_corr = max(orog_corr, epsilon(orog_corr))
    elsewhere (cosz_hor < cosz_beg .and. cosz_hor > cosz_end)
      ! Sun sets during timestep
      lit_frac = (cosz_beg - cosz_hor) / (cosz_beg - cosz_end)
      cosz_lit = (cosz_beg + cosz_hor) * 0.5_RealK
      orog_corr = lit_frac * ( cosz_lit + sin(acos(cosz_lit)) * &
        tan(slope_angle) * cos(sol_azimuth - slope_aspect) ) / cos_zenith_angle
      orog_corr = max(orog_corr, epsilon(orog_corr))
    end where
  end where
end if
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine solinc
end module solinc_mod
