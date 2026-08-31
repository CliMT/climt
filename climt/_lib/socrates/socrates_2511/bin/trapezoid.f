! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to perform trapezoidal integration.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      function trapezoid(n_points, x, y)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
      REAL  (RealK) ::
     &    trapezoid
!           Value of integral
      INTEGER, Intent(IN) ::
     &    n_points
!           Number of points
      REAL  (RealK), Intent(IN) ::
     &    x(n_points)
!           Abscissae
     &  , y(n_points)
!           Ordinates
!
!     Local arguments
      INTEGER
     &    i
!           Loop variable
!
!     Enter the integrating loop.
      trapezoid=0.0_RealK
      DO i=1, n_points-1
        trapezoid=trapezoid+(x(i+1)-x(i))*(y(i+1)+y(i))
      ENDDO
      trapezoid=0.5_RealK*trapezoid
!
!
!
      RETURN
      END
