! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set coefficients for a spline fit.
!
! Method:
!   The second derivatives of the cubic spline are found by solving
!   a matrix equation using natural boundary conditions.
!
!- ---------------------------------------------------------------------
MODULE spline_fit_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE spline_fit(n_data, x, y, y2)

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE


! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      n_data
!       Number of data points
  REAL (RealK), INTENT(IN) ::                                           &
      x(n_data)                                                         &
!       Absicissae of data
    , y(n_data)
!       Ordinates of data
  REAL (RealK), INTENT(OUT) ::                                          &
      y2(n_data)
!       Calculated second derivatives

! Local arguments.
  INTEGER ::                                                            &
      i
!       Loop variable
  REAL (RealK) ::                                                       &
      sigma                                                             &
!       Temporary variable
    , rho                                                               &
!       Temporary variable
    , u(n_data)
!       Temporary working variable


! No fit is possible with one term.
  IF (n_data == 1) RETURN

! This routine is essentially tridiagonal solution of a matrix
! equation. the starting point is the natural boundary condition.
  y2(1)=0.0_RealK
  u(1)=0.0_RealK
  DO i=2, n_data-1
    sigma=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    rho=sigma*y2(i-1)+2.0_RealK
    y2(i)=(sigma-1.0_RealK)/rho
    u(i)=(6.0_RealK*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))          &
      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sigma*u(i-1))/rho
  END DO
  y2(n_data)=0.0_RealK

! Backsubstitution.
  DO i=n_data-1, 1, -1
    y2(i)=y2(i)*y2(i+1)+u(i)
  END DO

END SUBROUTINE spline_fit
END MODULE spline_fit_mod
