! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Calculation of horizon angle along a 1D slice

module socrates_horizon_1d

implicit none
private
public :: horizon_1d

character(len=*), parameter :: ModuleName='SOCRATES_HORIZON_1D'

contains

subroutine horizon_1d(orog, nd_points, n_points, horiz_limit, &
                      sep_ang, horiz_ang1, horiz_ang2)

use realtype_rd, only: RealK
use rad_ccf, only: pi, planet_radius

implicit none

! Description:
!   Determines horizon angles in the forwards and backwards directions
!
! Method:
!   Uses a method similar to that outlined in Dozier and Frew (1990
!   IEEE Trans. Geosci. Remote Sens., 28, 963-969) adapted for use
!   within a limited horizon and allowing for the curvature of the
!   Earth.

integer, intent(in) :: nd_points, n_points, horiz_limit

real(RealK), intent(in) :: orog(1-horiz_limit:nd_points+horiz_limit)
real(RealK), intent(in) :: sep_ang(1-horiz_limit:nd_points+horiz_limit-1)

real(RealK), intent(out) :: horiz_ang1(nd_points), horiz_ang2(nd_points)

! Local constants
real(RealK), parameter :: pi_over_2 = pi/2.0_RealK

! Local variables
integer :: i, k
integer :: horiz_point(1-horiz_limit:n_points+horiz_limit)
real(RealK) :: ang, h_ang


horiz_ang1(:) = 0.0_RealK
horiz_point(:) = -horiz_limit

! First proceed forwards along the strip calculating horizon angles in
! the backwards direction
do i=1,n_points
  k=1
  ! Set maximum horizon angle to the horizontal
  horiz_ang1(i) = pi_over_2
  ! Initially each point is its own horizon
  horiz_point(i) = i
  do
    ! Loop backwards from each point until we reach the horizon limit
    if (k > horiz_limit) exit
    ! Only consider points with higher orography
    if ( orog(i-k) > orog(horiz_point(i)) ) then
      ang = sum(sep_ang(i-k:i-1))
      ! Horizon angle considers curvature of the earth
      h_ang = atan2( (planet_radius+orog(i-k))*sin(ang), &
        orog(i-k)*cos(ang) - orog(i) - planet_radius*(1.0_RealK-cos(ang)) )
      if (h_ang < horiz_ang1(i)) then
        ! Record the new horizon angle and the point which forms it
        horiz_ang1(i) = h_ang
        horiz_point(i) = i-k
      end if
      ! If the horizon angle to this point is smaller (zenith to horizon)
      ! than the further horizon angle viewed from this point then we
      ! need not consider any further points.
      if (i-k >= 1) then
        if (h_ang < horiz_ang1(i-k)) exit
      end if
    end if
    ! If this point is its own horizon we can stop
    if (horiz_point(i-k) == i-k) then
      exit
      ! If this point's own horizon is within our horizon limit then we
      ! can jump straight to the further horizon point
    else if (horiz_point(i-k) >= i-horiz_limit) then
      k=i-horiz_point(i-k)
      ! Otherwise we look to the next point along
    else
      k=k+1
    end if
  end do
end do

horiz_ang2(:) = 0.0_RealK
horiz_point(:) = n_points+horiz_limit+1

! Now proceed backwards along the strip calculating horizon angles in
! the forwards direction
do i=n_points,1,-1
  k=1
  horiz_ang2(i) = pi_over_2
  horiz_point(i) = i
  do
    if (k > horiz_limit) exit
    if ( orog(i+k) > orog(horiz_point(i)) ) then
      ang = sum(sep_ang(i:i+k-1))
      h_ang = atan2( (planet_radius+orog(i+k))*sin(ang), &
        orog(i+k)*cos(ang) - orog(i) - planet_radius*(1.0_RealK-cos(ang)) )
      if (h_ang < horiz_ang2(i)) then
        horiz_ang2(i) = h_ang
        horiz_point(i) = i+k
      end if
      if (i+k <= n_points) then
        if (h_ang < horiz_ang2(i+k)) exit
      end if
    end if
    if (horiz_point(i+k) == i+k) then
      exit
    else if (horiz_point(i+k) <= i+horiz_limit) then
      k=horiz_point(i+k)-i
    else
      k=k+1
    end if
  end do
end do

end subroutine horizon_1d
end module socrates_horizon_1d
