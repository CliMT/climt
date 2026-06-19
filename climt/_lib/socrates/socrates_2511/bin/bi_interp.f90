! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to interpolate irregular 2D data
SUBROUTINE bi_interp(x_in,y_in,z_in,rmdi,x_out,y_out,z_out)

  USE realtype_rd

  IMPLICIT NONE

  REAL (RealK), INTENT(IN)  :: x_in(4), y_in(4), z_in(4), rmdi
  REAL (RealK), INTENT(IN)  :: x_out, y_out
  REAL (RealK), INTENT(OUT) :: z_out

!   in(3)------------ in(1)
!        |          |
! Inc. ^ |    out   |
!  y   | |     .    |
!        |          |
!   in(4)------------ in(2)
!             ->
!        Increasing x


! Local variables
  REAL (RealK) :: fx(2),fy(2),x,y,d

! If only one good value then no need for interpolation
  IF (COUNT(z_in == rmdi) == 3) THEN
    z_out = MAXVAL(z_in,z_in /= rmdi)
    RETURN
  END IF

! Calculate fractional contributions
  IF (z_in(2) == rmdi) THEN
    fx(1) = 0.0_RealK
  ELSE IF (z_in(4) == rmdi) THEN
    fx(1) = 1.0_RealK
  ELSE
    fx(1) = (x_out-x_in(4))/(x_in(2)-x_in(4))
  END IF

  IF (z_in(1) == rmdi) THEN
    fx(2) = 0.0_RealK
  ELSE IF (z_in(3) == rmdi) THEN
    fx(2) = 1.0_RealK
  ELSE
    fx(2) = (x_out-x_in(3))/(x_in(1)-x_in(3))
  END IF

  IF (z_in(3) == rmdi) THEN
    fy(1) = 0.0_RealK
  ELSE IF (z_in(4) == rmdi) THEN
    fy(1) = 1.0_RealK
  ELSE
    fy(1) = (y_out-y_in(4))/(y_in(3)-y_in(4))
  END IF

  IF (z_in(1) == rmdi) THEN
    fy(2) = 0.0_RealK
  ELSE IF (z_in(2) == rmdi) THEN
    fy(2) = 1.0_RealK
  ELSE
    fy(2) = (y_out-y_in(2))/(y_in(1)-y_in(2))
  END IF

  IF (ANY(z_in == rmdi)) THEN
!   If there is missing data use simple weighted average:
    d = fx(2)*fy(2) + fx(1)*(1.0_RealK-fy(2)) &
      + (1.0_RealK-fx(2))*fy(1) + (1.0_RealK-fx(1))*(1.0_RealK-fy(1))
    z_out = ( z_in(1)*fx(2)*fy(2) +z_in(2)*fx(1)*(1.0_RealK-fy(2)) &
      + z_in(3)*(1.0_RealK-fx(2))*fy(1) &
      + z_in(4)*(1.0_RealK-fx(1))*(1.0_RealK-fy(1)) ) / d
  ELSE
!   Affine interpolation:
    d = fy(1)*fx(1) - fy(1)*fx(2) - fy(2)*fx(1) + fy(2)*fx(2) - 1.0_RealK
    x = ( fy(1)*fx(1) - fx(1) - fy(1)*fx(2) ) / d
    y = ( fy(1)*fx(1) - fy(1) - fy(2)*fx(1) ) / d
    z_out = z_in(1)*x*y + z_in(2)*x*(1.0_RealK-y) &
      + z_in(3)*(1.0_RealK-x)*y + z_in(4)*(1.0_RealK-x)*(1.0_RealK-y)
  END IF

END SUBROUTINE bi_interp
