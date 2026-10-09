! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to perform trapezoidal integration.
!
FUNCTION trapezoid_90(x, y) RESULT(int)
!
! Description:
!   This routine performs (non-uniform) trapezoidal integration
!   using specified absscissae and ordinates.
!
! Method:
!   Straightforward.
!
!
!
!
! Modules used
  USE realtype_rd
!
!
  IMPLICIT NONE
!
!
! Dummy arguments
  REAL  (RealK), Intent(IN) :: x(:)
!   Abscissae
  REAL  (RealK), Intent(IN) :: y(:)
!   Ordinates
  REAL  (RealK) :: int
!   Integral
!
! Local arguments
  INTEGER :: i
!   Loop variable
!
!
!
! Enter the integrating loop.
  int=0.0_RealK
  DO i=1, SIZE(x)-1
    int = int + (x(i+1) - x(i)) * (y(i+1) + y(i))
  ENDDO
  int = 0.5_RealK * int
!
!
!
  RETURN
END FUNCTION trapezoid_90
