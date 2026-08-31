! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to fit a single grey term to transmission data.
!
SUBROUTINE exponent_fit_90 &
!
(n, u, trans, k, ierr)
!
! Description:
!
! Calculates a k-value for grey transmissions.
!
! Method:
!
! Newton-Raphson iteration is used. This routone is similar to
! optimal_k.
!
! Modules used:
  USE realtype_rd
  USE ck_parm_acf
  USE def_std_io_icf
  USE error_pcf
!
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(IN) :: n
!   Number of data points
  REAL  (RealK), Intent(IN), Dimension(:) :: u
!   Amount of absoprber
  REAL  (RealK), Intent(IN), Dimension(:) :: trans
!   Transmissions to be fitted
!
  REAL  (RealK), Intent(OUT) :: k
!   Calculated absorption coefficient
!
!
!
! Local variables.
  INTEGER :: iteration
!   Number of iterations
  INTEGER :: i
!   Loop variable
  REAL  (RealK) :: error
!   Mean square error in the fitted transmission
  REAL  (RealK) :: d_sq_error
!   First derivative of the mean square error in the fitted transmission
  REAL  (RealK) :: d2_sq_error
!   Second derivative of the mean square error in the fitted transmission
  REAL  (RealK) :: k_prime
!   Provisional value of the absorption coefficient
  REAL  (RealK) :: x
!   Auxiliary variable
!
!
!
! Start from k=0.
  k = 0.0_RealK
!
  iteration = 0
!
! Iterate to improve value of k.
  DO
    error       = 0.0_RealK
    d_sq_error  = 0.0_RealK
    d2_sq_error = 0.0_RealK
    DO i = 1, n
      x = exp(-k * u(i))
      error = error + (x - trans(i))**2
      d_sq_error = d_sq_error - 2.0_RealK * u(i) * x * (x - trans(i))
      d2_sq_error = d2_sq_error + &
           2.0_RealK * x * u(i) * u(i) * &
           (2.0_RealK * x - trans(i))
    ENDDO
    error       = SQRT((1.0_RealK / REAL(n, RealK)) * error)
    d_sq_error  = (2.0_RealK / REAL(n, RealK)) * d_sq_error
    d2_sq_error = (2.0_RealK / REAL(n, RealK)) * d2_sq_error
!
!   Check for convergence.
    IF ( (ABS(d_sq_error) < 2048.0*EPSILON(d_sq_error)) .OR. &
         (ABS(d_sq_error/d2_sq_error) < EPSILON(d_sq_error)) ) THEN
      EXIT
    ELSE IF (iteration > np_kopt_max_iter) THEN
      WRITE(iu_err, '(/A)') &
        'Failure to converge in Newton-Raphson iteration.'
      ierr=i_abort_calculation
      RETURN
    ELSE IF (ABS(d2_sq_error) > EPSILON(d2_sq_error)) THEN
      iteration = iteration + 1
      IF (d2_sq_error > 0.0_RealK) THEN
        k_prime = k - d_sq_error/d2_sq_error
        IF (k_prime < 0.0_RealK) THEN
!         If the second derivative is such that k becomes negative.
!         Set k to half the previous value of k.
          k = k * 0.5_RealK
        ELSE
          k = k_prime
        ENDIF
      ELSE
!       If started from a high initial value k diverges
!       because the error asymptotes to a constant value for large k
!       and this convex approach to the asymptote appears to
!       Newton-Raphson iteration as an approach to a minimum. The
!       region where this occurs will be distinguished by a negative
!       second derivative.
        k = k * 0.5_RealK
      ENDIF
    ELSE
      WRITE(iu_err, '(/A)') &
        'Ill-conditioned division in Newton-Raphson iteration.'
      ierr=i_abort_calculation
      RETURN
    ENDIF
  ENDDO
!
!
!
  RETURN
!
!
!
END SUBROUTINE exponent_fit_90
