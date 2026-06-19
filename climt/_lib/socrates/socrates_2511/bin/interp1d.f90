! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE interp1d_mod
IMPLICIT NONE
CONTAINS

! One-dimensional linear interpolation function.
!
! Description:
!   Out-of-bounds values are set equal to the nearest edge value.
!
!------------------------------------------------------------------------------
FUNCTION interp1d(x ,y, xi, n) RESULT(yi)

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: &
      n
!         Size of x and y arrays
  REAL(RealK), INTENT(IN) :: &
      x(n) &
!         Array with x-values
    , y(n) &
!         Array with y-values
    , xi
!         Value at which the y-coordinate is wanted

  REAL(RealK) :: &
      yi
!         Value at xi

  INTEGER :: &
      i
!         Loop index

  IF (xi <= x(1)) THEN
    yi = y(1)
    RETURN
  ELSE IF (xi >= x(n)) THEN
    yi = y(n)
    RETURN
  ELSE
    DO i=2,n
      IF (xi <= x(i)) THEN
        yi = (y(i) - y(i-1))/(x(i) - x(i-1))*(xi - x(i-1)) + y(i-1)
        RETURN
      END IF
    END DO
  END IF

END FUNCTION interp1d
END MODULE interp1d_mod
