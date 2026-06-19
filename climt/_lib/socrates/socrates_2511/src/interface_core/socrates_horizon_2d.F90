! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Calculation of 16 horizon angles for a regular 2D grid

module socrates_horizon_2d

implicit none
private
public :: horizon_2d

character(len=*), parameter :: ModuleName='SOCRATES_HORIZON_2D'

contains

subroutine horizon_2d(n_x, n_y, n_horiz_layer, n_horiz_ang, horiz_limit, &
                      latitude, longitude, orog, horiz_ang, horiz_aspect)

use realtype_rd, only: RealK
use rad_ccf, only: pi
use socrates_horizon_1d, only: horizon_1d

implicit none

integer, intent(in) :: n_x, n_y, n_horiz_layer, n_horiz_ang, horiz_limit

real(RealK), intent(in) :: latitude(1-horiz_limit:n_x+horiz_limit, &
                                    1-horiz_limit:n_y+horiz_limit)
real(RealK), intent(in) :: longitude(1-horiz_limit:n_x+horiz_limit, &
                                     1-horiz_limit:n_y+horiz_limit)
real(RealK), intent(in) :: orog(1-horiz_limit:n_x+horiz_limit, &
                                1-horiz_limit:n_y+horiz_limit)

real(RealK), intent(out) :: horiz_ang(n_x, n_y, n_horiz_layer, n_horiz_ang)
real(RealK), intent(out) :: horiz_aspect(n_x, n_y, n_horiz_ang)

! Local constants

integer, parameter :: ip_N   = 1
integer, parameter :: ip_NE  = 2
integer, parameter :: ip_E   = 3
integer, parameter :: ip_SE  = 4
integer, parameter :: ip_S   = 5
integer, parameter :: ip_SW  = 6
integer, parameter :: ip_W   = 7
integer, parameter :: ip_NW  = 8

integer, parameter :: ip_NNE = 9
integer, parameter :: ip_ENE = 10
integer, parameter :: ip_ESE = 11
integer, parameter :: ip_SSE = 12
integer, parameter :: ip_SSW = 13
integer, parameter :: ip_WSW = 14
integer, parameter :: ip_WNW = 15
integer, parameter :: ip_NNW = 16

real(RealK), parameter :: pi_over_2 = pi/2.0_RealK
real(RealK), parameter :: twopi = pi*2.0_RealK

! Local variables

integer :: i,ii,j,jj,k
integer :: nd_points

real(RealK) :: sinlat(1-horiz_limit:n_x+horiz_limit, &
                      1-horiz_limit:n_y+horiz_limit)
real(RealK) :: coslat(1-horiz_limit:n_x+horiz_limit, &
                      1-horiz_limit:n_y+horiz_limit)
real(RealK) :: max_ang

real(RealK) :: sep_ang(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit-1)
real(RealK) :: lon_1d(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit-1)
real(RealK) :: dlon(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit-1)
real(RealK) :: cos_dlon(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit-1)
real(RealK) :: sinlat_1d(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit-1)
real(RealK) :: coslat_1d(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit-1)
real(RealK) :: orog_1d(1-horiz_limit:2*MAX(n_x,n_y)+horiz_limit)
real(RealK) :: horiz_ang1(2*MAX(n_x,n_y))
real(RealK) :: horiz_ang2(2*MAX(n_x,n_y))


nd_points=2*max(n_x,n_y)

! Calculate arrays used more than once
sinlat=sin(latitude)
coslat=cos(latitude)

! Max separation angle to catch errors where 1D slice wraps around grid
max_ang = pi_over_2/real(horiz_limit, RealK)

! Assume grid directions approximate to great circles from the local
! bearings. Horizon angles:
!  8   1   2
!   \  |  /
!    \ | /
! 7---   ---3
!    / | \
!   /  |  \
!  6   5   4


! Calculate horizon angles along 1D slices.

! Horizon angles W and E
! ----------------------
do j=1, n_y
  i = 1-horiz_limit
  orog_1d(i) = orog(i,j)
  do i=2-horiz_limit, n_x+horiz_limit
    orog_1d(i) = orog(i,j)
    dlon(i-1) = longitude(i,j)-longitude(i-1,j)
    cos_dlon(i-1) = cos(dlon(i-1))
    sep_ang(i-1) = min( max_ang, acos( sinlat(i-1,j) * sinlat(i,j) &
      + coslat(i-1,j) * coslat(i,j) * cos_dlon(i-1) ) )
  end do

  call horizon_1d(orog_1d, nd_points, n_x, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  do i=1, n_x
    horiz_ang(i,j,1,ip_W) = horiz_ang1(i)
    horiz_ang(i,j,1,ip_E) = horiz_ang2(i)
    horiz_aspect(i,j,ip_W) = modulo( atan2( &
      - coslat(i-1,j) * sin(dlon(i-1)), sinlat(i-1,j) * coslat(i,j) &
      - coslat(i-1,j) * sinlat(i,j) * cos_dlon(i-1) ), twopi )
    horiz_aspect(i,j,ip_E) = modulo( atan2( &
        coslat(i+1,j) * sin(dlon(i)), sinlat(i+1,j) * coslat(i,j) &
      - coslat(i+1,j) * sinlat(i,j) * cos_dlon(i) ), twopi )
  end do
end do

! Horizon angles S and N
! ----------------------
do i=1, n_x
  j = 1-horiz_limit
  orog_1d(j) = orog(i,j)
  do j=2-horiz_limit, n_y+horiz_limit
    orog_1d(j) = orog(i,j)
    dlon(j-1) = longitude(i,j)-longitude(i,j-1)
    cos_dlon(j-1) = cos(dlon(j-1))
    sep_ang(j-1) = min( max_ang, acos( sinlat(i,j-1) * sinlat(i,j) &
      + coslat(i,j-1) * coslat(i,j) * cos_dlon(j-1) ) )
  end do

  call horizon_1d(orog_1d, nd_points, n_y, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  do j=1, n_y
    horiz_ang(i,j,1,ip_S) = horiz_ang1(j)
    horiz_ang(i,j,1,ip_N) = horiz_ang2(j)
    horiz_aspect(i,j,ip_S) = modulo( atan2( &
      - coslat(i,j-1) * sin(dlon(j-1)), sinlat(i,j-1) * coslat(i,j) &
      - coslat(i,j-1) * sinlat(i,j) * cos_dlon(j-1) ), twopi )
    horiz_aspect(i,j,ip_N) = modulo( atan2( &
        coslat(i,j+1) * sin(dlon(j)), sinlat(i,j+1) * coslat(i,j) &
      - coslat(i,j+1) * sinlat(i,j) * cos_dlon(j) ), twopi )
  end do
end do


! Slices along the diagonals use a sheared array.

! Horizon angles NW and SE
! ------------------------
do j=1-horiz_limit, n_y+horiz_limit
  i = 1-horiz_limit
  orog_1d(i) = orog(i,j)
  sinlat_1d(i) = sinlat(i,j)
  coslat_1d(i) = coslat(i,j)
  lon_1d(i) = longitude(i,j)
  ! Map the diagonals onto 1D slices, wrapping around the grid to
  ! minimise the number of loops required
  do i=2-horiz_limit, n_x+horiz_limit
    k = modulo(j-i, n_y+2*horiz_limit) + 1-horiz_limit
    orog_1d(i) = orog(i,k)
    sinlat_1d(i) = sinlat(i,k)
    coslat_1d(i) = coslat(i,k)
    lon_1d(i) = longitude(i,k)
    dlon(i-1) = lon_1d(i) - lon_1d(i-1)
    cos_dlon(i-1) = cos(dlon(i-1))
    sep_ang(i-1)=min(max_ang, acos( sinlat_1d(i) * sinlat_1d(i-1) &
      + coslat_1d(i) * coslat_1d(i-1) * cos_dlon(i-1) ) )
  end do

  call horizon_1d(orog_1d, nd_points, n_x, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  ! Map the slices back on to the 2D grid
  do i=1, n_x
    k = modulo(j-i, n_y+2*horiz_limit) + 1-horiz_limit
    if (k >= 1 .and. k <= n_y) then
      horiz_ang(i,k,1,ip_NW) = horiz_ang1(i)
      horiz_ang(i,k,1,ip_SE) = horiz_ang2(i)
      horiz_aspect(i,k,ip_NW) = modulo( atan2( &
        - coslat_1d(i-1) * sin(dlon(i-1)), sinlat_1d(i-1) * coslat_1d(i) &
        - coslat_1d(i-1) * sinlat_1d(i) * cos_dlon(i-1) ), twopi )
      horiz_aspect(i,k,ip_SE) = modulo( atan2( &
          coslat_1d(i+1) * sin(dlon(i)), sinlat_1d(i+1) * coslat_1d(i) &
        - coslat_1d(i+1) * sinlat_1d(i) * cos_dlon(i) ), twopi )
   end if
  end do
end do

! Horizon angles SW and NE
! ------------------------
do j=1-horiz_limit, n_y+horiz_limit
  i = 1-horiz_limit
  orog_1d(i) = orog(i,j)
  sinlat_1d(i) = sinlat(i,j)
  coslat_1d(i) = coslat(i,j)
  lon_1d(i) = longitude(i,j)
  do i=2-horiz_limit, n_x+horiz_limit
    k = modulo(j+i-2*(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    orog_1d(i) = orog(i,k)
    sinlat_1d(i) = sinlat(i,k)
    coslat_1d(i) = coslat(i,k)
    lon_1d(i) = longitude(i,k)
    dlon(i-1) = lon_1d(i) - lon_1d(i-1)
    cos_dlon(i-1) = cos(dlon(i-1))
    sep_ang(i-1)=min(max_ang, acos( sinlat_1d(i) * sinlat_1d(i-1) &
      + coslat_1d(i) * coslat_1d(i-1) * cos_dlon(i-1) ) )
  end do

  call horizon_1d(orog_1d, nd_points, n_x, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  do i=1, n_x
    k = modulo(j+i-2*(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    if (k >= 1 .and. k <= n_y) then
      horiz_ang(i,k,1,ip_SW) = horiz_ang1(i)
      horiz_ang(i,k,1,ip_NE) = horiz_ang2(i)
      horiz_aspect(i,k,ip_SW) = modulo( atan2( &
        - coslat_1d(i-1) * sin(dlon(i-1)), sinlat_1d(i-1) * coslat_1d(i) &
        - coslat_1d(i-1) * sinlat_1d(i) * cos_dlon(i-1) ), twopi )
      horiz_aspect(i,k,ip_NE) = modulo( atan2( &
          coslat_1d(i+1) * sin(dlon(i)), sinlat_1d(i+1) * coslat_1d(i) &
        - coslat_1d(i+1) * sinlat_1d(i) * cos_dlon(i) ), twopi )
    end if
  end do
end do


! The intermediate diagonals require a shear of 2 grid-boxes, with
! points in-between given by the mean of the 2 points bisected by the
! line of sight.

! Horizon angles NNW and SSE
! --------------------------
do j=1-horiz_limit, n_y+horiz_limit
  i = 1 - horiz_limit/2
  k = modulo(j-2*i+(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
  orog_1d(2*i) = orog(i,k)
  sinlat_1d(i) = sinlat(i,k)
  coslat_1d(i) = coslat(i,k)
  lon_1d(i) = longitude(i,k)
  ! 1D slices are created by looping along the slowest changing
  ! direction, setting up both the main grid-points and the previous
  ! bisected point on each loop
  do i=2-horiz_limit/2, n_x+horiz_limit/2
    k = modulo(j-2*i+(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    sinlat_1d(i) = sinlat(i,k)
    coslat_1d(i) = coslat(i,k)
    lon_1d(i) = longitude(i,k)
    dlon(i-1) = lon_1d(i) - lon_1d(i-1)
    cos_dlon(i-1) = cos(dlon(i-1))
    ! The separation angle between main points is calculated and then
    ! halved to give the angle to the bisected points
    sep_ang(2*i-2)=min(max_ang, acos( sinlat_1d(i) * sinlat_1d(i-1) &
      + coslat_1d(i) * coslat_1d(i-1) * cos_dlon(i-1) ) * 0.5_RealK)
    sep_ang(2*i-1) = sep_ang(2*i-2)
    orog_1d(2*i) = orog(i,k)
    jj = modulo(k+1-(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    ! The orographic height of the bisected point is interpolated
    orog_1d(2*i-1) = (orog(i,jj)+orog(i-1,jj))*0.5_RealK
  end do
  ! First point in array is beyond the horizon so use arbitrary values
  orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
  sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

  call horizon_1d(orog_1d, nd_points, 2*n_x, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  ! The horizon angles calculated for the main grid-points are then
  ! mapped back on to the 2D grid
  do i=1, n_x
    k = modulo(j-2*i+(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    if (k >= 1 .and. k <= n_y) then
      horiz_ang(i,k,1,ip_NNW) = horiz_ang1(2*i)
      horiz_ang(i,k,1,ip_SSE) = horiz_ang2(2*i)
      horiz_aspect(i,k,ip_NNW) = modulo( atan2( &
        - coslat_1d(i-1) * sin(dlon(i-1)), sinlat_1d(i-1) * coslat_1d(i) &
        - coslat_1d(i-1) * sinlat_1d(i) * cos_dlon(i-1) ), twopi )
      horiz_aspect(i,k,ip_SSE) = modulo( atan2( &
          coslat_1d(i+1) * sin(dlon(i)), sinlat_1d(i+1) * coslat_1d(i) &
        - coslat_1d(i+1) * sinlat_1d(i) * cos_dlon(i) ), twopi )
    end if
  end do
end do

! Horizon angles SSW and NNE
! --------------------------
do j=1-horiz_limit, n_y+horiz_limit
  i = 1 - horiz_limit/2
  k = modulo(j+2*i-3*(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
  orog_1d(2*i) = orog(i,k)
  sinlat_1d(i) = sinlat(i,k)
  coslat_1d(i) = coslat(i,k)
  lon_1d(i) = longitude(i,k)
  do i=2-horiz_limit/2, n_x+horiz_limit/2
    k = modulo(j+2*i-3*(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    sinlat_1d(i) = sinlat(i,k)
    coslat_1d(i) = coslat(i,k)
    lon_1d(i) = longitude(i,k)
    dlon(i-1) = lon_1d(i) - lon_1d(i-1)
    cos_dlon(i-1) = cos(dlon(i-1))
    sep_ang(2*i-2)=min(max_ang, acos( sinlat_1d(i) * sinlat_1d(i-1) &
      + coslat_1d(i) * coslat_1d(i-1) * cos_dlon(i-1) ) * 0.5_RealK)
    sep_ang(2*i-1) = sep_ang(2*i-2)
    orog_1d(2*i) = orog(i,k)
    jj = modulo(k-1-(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    orog_1d(2*i-1) = (orog(i,jj)+orog(i-1,jj))*0.5_RealK
  end do
  orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
  sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

  call horizon_1d(orog_1d, nd_points, 2*n_x, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  do i=1, n_x
    k = modulo(j+2*i-3*(1-horiz_limit), n_y+2*horiz_limit) + 1-horiz_limit
    if (k >= 1 .and. k <= n_y) then
      horiz_ang(i,k,1,ip_SSW) = horiz_ang1(2*i)
      horiz_ang(i,k,1,ip_NNE) = horiz_ang2(2*i)
      horiz_aspect(i,k,ip_SSW) = modulo( atan2( &
        - coslat_1d(i-1) * sin(dlon(i-1)), sinlat_1d(i-1) * coslat_1d(i) &
        - coslat_1d(i-1) * sinlat_1d(i) * cos_dlon(i-1) ), twopi )
      horiz_aspect(i,k,ip_NNE) = modulo( atan2( &
          coslat_1d(i+1) * sin(dlon(i)), sinlat_1d(i+1) * coslat_1d(i) &
        - coslat_1d(i+1) * sinlat_1d(i) * cos_dlon(i) ), twopi )
    end if
  end do
end do


! Horizon angles ESE and WNW
! --------------------------
do i=1-horiz_limit, n_x+horiz_limit
  j = 1 - horiz_limit/2
  k = modulo(i-2*j+(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
  orog_1d(2*j) = orog(k,j)
  sinlat_1d(j) = sinlat(k,j)
  coslat_1d(j) = coslat(k,j)
  lon_1d(j) = longitude(k,j)
  do j=2-horiz_limit/2, n_y+horiz_limit/2
    k = modulo(i-2*j+(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
    sinlat_1d(j) = sinlat(k,j)
    coslat_1d(j) = coslat(k,j)
    lon_1d(j) = longitude(k,j)
    dlon(j-1) = lon_1d(j) - lon_1d(j-1)
    cos_dlon(j-1) = cos(dlon(j-1))
    sep_ang(2*j-2)=min(max_ang, acos( sinlat_1d(j) * sinlat_1d(j-1) &
      + coslat_1d(j) * coslat_1d(j-1) * cos_dlon(j-1) ) * 0.5_RealK)
    sep_ang(2*j-1) = sep_ang(2*j-2)
    orog_1d(2*j) = orog(k,j)
    ii = modulo(k+1-(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
    orog_1d(2*j-1) = (orog(ii,j)+orog(ii,j-1))*0.5_RealK
  END DO
  orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
  sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

  call horizon_1d(orog_1d, nd_points, 2*n_y, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  do j=1, n_y
    k = modulo(i-2*j+(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
    if (k >= 1 .and. k <= n_x) then
      horiz_ang(k,j,1,ip_ESE) = horiz_ang1(2*j)
      horiz_ang(k,j,1,ip_WNW) = horiz_ang2(2*j)
      horiz_aspect(k,j,ip_ESE) = modulo( atan2( &
        - coslat_1d(j-1) * sin(dlon(j-1)), sinlat_1d(j-1) * coslat_1d(j) &
        - coslat_1d(j-1) * sinlat_1d(j) * cos_dlon(j-1) ), twopi )
      horiz_aspect(k,j,ip_WNW) = modulo( atan2( &
          coslat_1d(j+1) * sin(dlon(j)), sinlat_1d(j+1) * coslat_1d(j) &
        - coslat_1d(j+1) * sinlat_1d(j) * cos_dlon(j) ), twopi )
    end if
  end do
end do

! Horizon angles WSW and ENE
! --------------------------
do i=1-horiz_limit, n_x+horiz_limit
  j = 1 - horiz_limit/2
  k = modulo(i+2*j-3*(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
  orog_1d(2*j) = orog(k,j)
  sinlat_1d(j) = sinlat(k,j)
  coslat_1d(j) = coslat(k,j)
  lon_1d(j) = longitude(k,j)
  do j=2-horiz_limit/2, n_y+horiz_limit/2
    k = modulo(i+2*j-3*(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
    sinlat_1d(j) = sinlat(k,j)
    coslat_1d(j) = coslat(k,j)
    lon_1d(j) = longitude(k,j)
    dlon(j-1) = lon_1d(j) - lon_1d(j-1)
    cos_dlon(j-1) = cos(dlon(j-1))
    sep_ang(2*j-2)=min(max_ang, acos( sinlat_1d(j) * sinlat_1d(j-1) &
      + coslat_1d(j) * coslat_1d(j-1) * cos_dlon(j-1) ) * 0.5_RealK)
    sep_ang(2*j-1) = sep_ang(2*j-2)
    orog_1d(2*j) = orog(k,j)
    ii = modulo(k-1-(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
    orog_1d(2*j-1) = (orog(ii,j)+orog(ii,j-1))*0.5_RealK
  end do
  orog_1d(1-horiz_limit) = orog_1d(2-horiz_limit)
  sep_ang(1-horiz_limit) = sep_ang(2-horiz_limit)

  call horizon_1d(orog_1d, nd_points, 2*n_y, horiz_limit, &
                  sep_ang, horiz_ang1, horiz_ang2)

  do j=1, n_y
    k = modulo(i+2*j-3*(1-horiz_limit), n_x+2*horiz_limit) + 1-horiz_limit
    if (k >= 1 .and. k <= n_x) then
      horiz_ang(k,j,1,ip_WSW) = horiz_ang1(2*j)
      horiz_ang(k,j,1,ip_ENE) = horiz_ang2(2*j)
      horiz_aspect(k,j,ip_WSW) = modulo( atan2( &
        - coslat_1d(j-1) * sin(dlon(j-1)), sinlat_1d(j-1) * coslat_1d(j) &
        - coslat_1d(j-1) * sinlat_1d(j) * cos_dlon(j-1) ), twopi )
      horiz_aspect(k,j,ip_ENE) = modulo( atan2( &
          coslat_1d(j+1) * sin(dlon(j)), sinlat_1d(j+1) * coslat_1d(j) &
        - coslat_1d(j+1) * sinlat_1d(j) * cos_dlon(j) ), twopi )
    end if
  end do
end do

end subroutine horizon_2d
end module socrates_horizon_2d
