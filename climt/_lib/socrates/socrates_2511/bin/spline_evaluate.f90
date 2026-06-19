! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to evaluate a cubic spline.
!
!- ---------------------------------------------------------------------
MODULE spline_evaluate_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE spline_evaluate(ierr, n_data, x, y, y2, x_point, y_point)

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_range

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_data
!       Number of data points
  REAL (RealK), INTENT(IN) ::                                           &
      x(n_data)                                                         &
!       Abscissae of data
    , y(n_data)                                                         &
!       Ordinates of data
    , y2(n_data)                                                        &
!       Second derivative
    , x_point
!       Point for evaluation
  REAL (RealK), INTENT(OUT) ::                                          &
      y_point
!       Value of spline

! Local arguments.
  INTEGER ::                                                            &
      i_low                                                             &
!       Lower limit of bracketing pair
    , i_high                                                            &
!       Upper limit of bracketing pair
    , i
!       Index of new point
  REAL (RealK) ::                                                       &
      delta                                                             &
!       Width of interval
    , a                                                                 &
!       Weight of lower point
    , b                                                                 &
!       Weight of upper point
    , c
!       Temporary store


  IF (n_data == 1) THEN
!   Return the only available value.
    y_point=y(1)
  ELSE IF ( (x_point > x(n_data)).OR.(x_point < x(1)) ) THEN
!   Outside the range: return an error flag for corrective action.
    ierr=i_err_range
    RETURN
  ELSE
!   Find the bracketing points by a binary chop.
    i_low=1
    i_high=n_data
    DO
      IF ((i_high-i_low) <=  1) EXIT
      i=(i_low+i_high)/2
      IF (x_point > x(i)) THEN
        i_low=i
      ELSE
        i_high=i
      END IF
    END DO
    delta=x(i_high)-x(i_low)
    a=(x(i_high)-x_point)/delta
    b=(x_point-x(i_low))/delta
    c=(delta**2/6.0_RealK)
    y_point=a*(y(i_low)+(a*a-1)*y2(i_low)*c)                            &
      +b*(y(i_high)+(b*b-1.0_RealK)*y2(i_high)*c)
  END IF

END SUBROUTINE spline_evaluate
END MODULE spline_evaluate_mod
