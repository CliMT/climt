! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the optimal k-value for each term.
!
SUBROUTINE optimal_k &
!
(n_nu, nu_inc, k, wgt, integ_wgt, k_mean, tol, &
 l_k_wgt, k_wgt, l_wgt_scale_sqrt, u_wgt_scale, k_opt, error, ierr)
!
! Description:
!
! Calculates optimal k using Newton-Raphson method
!
! Method:
!
! Uses standard formula
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE error_pcf
  USE ck_parm_acf
!
!
  IMPLICIT NONE
!
!
!
! Subroutine arguments
  INTEGER, Intent(IN) :: n_nu
!   Number of frequencies 
  INTEGER, Intent(INOUT) :: ierr
!   Error status
  REAL  (RealK), Intent(IN) :: nu_inc
!   Frequencies increment
  REAL  (RealK), Intent(IN) :: k(:)
!   Monochromatic absorption coefficients
  REAL  (RealK), Intent(IN) :: wgt(:)
!   Weightings for monochromatic absorption
  REAL  (RealK), Intent(IN) :: integ_wgt
!   Integral of the weighting function
  REAL  (RealK), Intent(IN) :: k_mean
!   Simple mean absorption coefficient for the band
  REAL  (RealK), Intent(IN) :: tol
!   Tolerance required of the fit
  REAL  (RealK), Intent(IN) :: k_wgt(:)
!   Monochromatic absorption coefficients to use in weighting
  REAL  (RealK), Intent(IN) :: u_wgt_scale
!   Factor to scale continuum column mass to gas column mass
  REAL  (RealK), Intent(OUT) :: k_opt
!   Optimal absorption coefficient
  REAL  (RealK), Intent(OUT) :: error
!   Root mean square error in the fitted transmission
  LOGICAL, Intent(IN) :: l_k_wgt
!   Use k_wgt as additional weights in transmissions
  LOGICAL, Intent(IN) :: l_wgt_scale_sqrt
!   If true gas column mass is scaled using the square root of the continuum
!   column mass. Otherwise a linear scaling is used. 
!
!
! Local variables
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: iter, kmap(n_nu)
!   Counter for iterations
  REAL  (RealK) :: umin
!   Minimum pathlength for the absorber
  REAL  (RealK) :: umax
!   Maximum pathlength for the absorber
  REAL  (RealK) :: u(n_path_kopt_default)
!   Path lengths of absorber
  REAL  (RealK) :: u_wgt
!   Path length of gas used as weight in continuum transmissions
  REAL  (RealK), Allocatable :: wgt_trans(:)
!   Product of weightings and transmissions along defined paths
  REAL  (RealK) :: trans_exact(n_path_kopt_default)
!   Exact transmissions along the path
  REAL  (RealK) :: exku
!   Variable for efficiency
  REAL  (RealK) :: d_sq_error
!   First derivative of the mean square error in the fitted transmission
  REAL  (RealK) :: d2_sq_error
!   Second derivative of the mean square error in the fitted transmission
  REAL  (RealK) :: k_prime
!   Provisional value of the absorption coefficient
  REAL  (RealK) :: dk
!   Last step size in the Newton-Raphson iteration.
!
  INTERFACE
    SUBROUTINE map_heap_func(a, map)
      USE realtype_rd
      REAL  (RealK), Intent(IN), Dimension(:) :: a
      INTEGER, Intent(OUT), Dimension(:) :: map
    END SUBROUTINE map_heap_func
  END INTERFACE
!   
!- End of header
!
!
!
! Set limits on the range of amounts of absorber to cover the
! key region of the curve of growth.
  umax=umax_kopt
  IF ( MINVAL(k(1:n_nu)) > 0.0_RealK ) &
    umax = MIN(umax_kopt, -LOG(tol)/MINVAL(k(1:n_nu)))
  umin=umin_kopt
  IF ( MAXVAL(k(1:n_nu)) > 0.0_RealK ) &
    umin = MAX(umin_kopt, -LOG(1.0_RealK-tol)/MAXVAL(k(1:n_nu)))
! Weak absorption may occasionally give transposed lines: in this
! case the grey extinction will be optimal.
  IF (umin > umax) then
    k_opt = k_mean
    error = 0.0_RealK
    RETURN
  ENDIF
!
! Set a suitable range of amounts of absorber.
  DO i=1, n_path_kopt_default
    u(i)=umin*EXP( REAL(i-1, RealK) * &
      LOG(umax/umin) / REAL((n_path_kopt_default-1), RealK))
  ENDDO
!
! Calculate the exact transmission over these paths.
  ALLOCATE(wgt_trans(n_nu))
  IF (l_k_wgt) THEN
    DO i = 1, n_path_kopt_default
      IF (l_wgt_scale_sqrt) THEN
        u_wgt=u_wgt_scale*SQRT(u(i))
      ELSE
        u_wgt=u_wgt_scale*u(i)
      END IF
      DO j=1, n_nu
        wgt_trans(j)=wgt(j)*EXP(-k(j)*u(i)-k_wgt(j)*u_wgt)
      ENDDO
      trans_exact(i)=SUM(wgt_trans(1:n_nu)) &
        / (SUM(wgt(1:n_nu)*EXP(-k_wgt(1:n_nu)*u_wgt)))
    ENDDO
  ELSE
    DO i = 1, n_path_kopt_default
      DO j=1, n_nu
        wgt_trans(j)=wgt(j)*EXP(-k(j)*u(i))
      ENDDO
      trans_exact(i)=nu_inc * SUM(wgt_trans(1:n_nu)) / integ_wgt
    ENDDO
  END IF
  DEALLOCATE(wgt_trans)
!
! Carry out Newton-Raphson iteration to mimimize the squared
! error in the fit, initializing with the lowest absorption.
  k_opt=MINVAL(k(1:n_nu))
  iter=0

! Calculate mid-point k value:
!  CALL map_heap_func(k(1:n_nu), kmap)
!  k_opt=k(kmap((n_nu+1)/2))

  DO
!
!   Calculate the mean square error and its first two derivates at 
!   the estimated k.
    error       = 0.0_RealK
    d_sq_error  = 0.0_RealK
    d2_sq_error = 0.0_RealK
    DO i=1, n_path_kopt_default
      exku=EXP(-k_opt*u(i))
      error=error+(exku-trans_exact(i))**2
      d_sq_error=d_sq_error-u(i)*exku*(exku-trans_exact(i))
      d2_sq_error=d2_sq_error+u(i)**2*exku* &
        (2.0_RealK*exku-trans_exact(i))
    ENDDO
    error=SQRT((1.0_RealK/REAL(n_path_kopt_default, RealK))*error)
    d_sq_error=(2.0_RealK/REAL(n_path_kopt_default, RealK))*d_sq_error
    d2_sq_error=(2.0_RealK/REAL(n_path_kopt_default, RealK))*d2_sq_error
!
!   Check for convergence.
    IF ( ((error < tol).AND.(iter > np_kopt_max_iter)) .OR. &
         (ABS(d_sq_error) < 2048.0*EPSILON(d_sq_error)) .OR. &
         (ABS(d_sq_error/d2_sq_error) < EPSILON(d_sq_error)) ) THEN
      ierr=i_normal
      EXIT
    ELSE IF (iter > np_kopt_max_iter) THEN
      IF (ierr == i_normal) THEN
        ! Attempt convergence starting from k_mean
        iter = 0
        k_opt = k_mean
        ierr = i_warning
      ELSE
        WRITE(iu_err, '(/A)') &
          'Error: Failure to converge in Newton-Raphson iteration.'
        ! If failed to converge again then abort
        ierr = i_abort_calculation
        RETURN
      END IF
    ELSE IF (ABS(d2_sq_error) > EPSILON(d2_sq_error)) THEN
      iter=iter+1
      IF (d2_sq_error > 0.0_RealK) THEN
        k_prime = k_opt - d_sq_error/d2_sq_error
        IF (k_prime < 0.0_RealK) THEN
!         If the second derivative is such that k_opt becomes negative.
!         Set k_k_opt to half the previous value of k_opt.
          k_opt = k_opt * 0.5_RealK
          dk = -k_opt
        ELSE IF ( iter > 1 .AND.                                        &
                  dk*d_sq_error/d2_sq_error > 0.0_RealK .AND.           &
                  ABS(2.0_RealK*d_sq_error) > ABS(dk*d2_sq_error) ) THEN
!         To prevent the iteration oscillating around a solution we
!         halve the step size when the increment changes sign and
!         convergence is slow.
          dk = -d_sq_error/(2.0_RealK*d2_sq_error)
          k_opt = k_opt + dk
        ELSE
          dk = k_prime - k_opt
          k_opt = k_prime
        ENDIF
      ELSE
!       If started from a high initial value k_opt diverges
!       because the error asymptotes to a constant value for large k
!       and this convex approach to the asymptote appears to
!       Newton-Raphson iteration as an approach to a minimum. The
!       region where this occurs will be distinguished by a negative
!       second derivative.
        k_opt = k_opt * 0.5_RealK
        dk = -k_opt
      ENDIF
    ELSE
      IF (ierr == i_normal) THEN
        ! Attempt convergence starting from k_mean
        iter = 0
        k_opt = k_mean
        ierr = i_warning
      ELSE
        WRITE(iu_err, '(/A,A)') 'Warning: ', &
          'Ill-conditioned division in Newton-Raphson iteration.'
        ! If failed to converge again then use k_mean
        k_opt = k_mean
        error=0.0_RealK
        ierr=i_normal
        EXIT
      END IF
    ENDIF
  ENDDO

END SUBROUTINE optimal_k
       
