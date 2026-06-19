! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to perform scalar minimization in a fixed direction.
!
SUBROUTINE line_search_cloud_90(ierr, iu_monitor, &
  fit_species, &
  n_data, radius_effective, vol_fraction, data, &
  i_fit_type, property, n_parameter, parm, &
  h, h_l2_norm, search_distance  &
  )
!
! Description:
!   The value of a parameter, S, giving the distance travelled
!   along a direction, H, is adjusted to minimize a function F.
!   The minimum is first bracketed and then Brent's method is
!   applied to find the minimum. This subroutine applies 
!   specifically to the fitting of cloud data.
!
!
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE cloud_fit_parm_acf
  USE cloud_fitting
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
! Include header files.
!
! Dummy arguments
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(INOUT) :: iu_monitor
!   Unit number for monitoring information
  CHARACTER  (LEN=5), Intent(IN) :: fit_species
!   Phase of condensate to be fitted
  INTEGER, Intent(IN) :: n_data
!   Number of data points
  INTEGER, Intent(IN) :: i_fit_type
!   Type of fit required
  CHARACTER  (LEN=10), Intent(IN) :: property
!   Single-scattering property to be fitted
  INTEGER, Intent(IN) :: n_parameter
!   Number of parameters
  REAL  (RealK), Intent(IN), Dimension(n_data) :: vol_fraction
!   Volume fraction of absorber
  REAL  (RealK), Intent(IN), Dimension(n_data) :: radius_effective
!   Effective radii of particles
  REAL  (RealK), Intent(IN), Dimension(n_data) :: data
!   Data to be fitted
  REAL  (RealK), Intent(IN), Dimension(n_parameter) :: h
!   Searching direction
  REAL  (RealK), Intent(IN) :: h_l2_norm
!   L2-norm of H
  REAL  (RealK), Intent(INOUT), Dimension(n_parameter) :: parm
!   Fitting parameters
  REAL  (RealK), Intent(INOUT) :: search_distance
!   Searching distance
!
!
!
! Local arguments.
  INTEGER :: i
!   Loop variable
  INTEGER :: iteration
!   Number of iteration
  LOGICAL :: l_fit
!   Variable for fitting
  REAL  (RealK) :: s_low
!   Lower limit of enclosing interval
  REAL  (RealK) :: s_high
!   Upper limit of enclosing interval
  REAL  (RealK) :: s_middle
!   Mid-point of interval
  REAL  (RealK) :: s_first
!   Abscissa of lowest point so far
  REAL  (RealK) :: s_second
!   Abscissa of sceond lowest point
  REAL  (RealK) :: s_third
!   Absicissa of third lowest
  REAL  (RealK) :: s_new
!   Abscissa of new point
  REAL  (RealK) :: s_initial
!   Initial searching distance
  REAL  (RealK) :: f_first
!   Lowest point so far
  REAL  (RealK) :: f_second
!   Sceond lowest point
  REAL  (RealK) :: f_third
!   Third lowest
  REAL  (RealK) :: f_new
!   New point
  REAL  (RealK) :: step
!   Step in s
  REAL  (RealK) :: penult_step
!   Previous step
  REAL  (RealK) :: tol_scaled
!   Scaled tolerance for minimum
  REAL  (RealK) :: s0
!   Position of minimum of quadratic
  REAL  (RealK), Dimension(n_parameter) :: x1
!   Search position
  REAL  (RealK), Dimension(n_parameter) :: x2
!   Search position
  REAL  (RealK), Dimension(n_parameter) :: x3
!   Search position
  REAL  (RealK), Dimension(n_parameter) :: x_new
!   New search position
  REAL  (RealK), Parameter :: golden_ratio = 6.1803399E-01_RealK
!   Golden ratio for search
!
!
! Subroutines called:
  EXTERNAL &
      fit_parabola_90
!
! Variables related to the treatment of ill-conditioning
  REAL  (RealK) :: tol_conv
!   The tolerance used to assess convergence
  REAL  (RealK) :: sq_tol_conv
!   The square root of the above
!
!
!
! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  tol_conv=EPSILON(s_new)
  sq_tol_conv=SQRT(tol_conv)
!
! Record the initial searching distance.
  s_initial=search_distance/h_l2_norm
  s_new=s_initial
!
! Find three points enclosing the minimum.
  DO
    x1(:)=parm(:)
    x2(:)=parm(:)-s_new*h(:)
    x3(:)=parm(:)+s_new*h(:)
!
!   Find the residuals at these points.
    f_first=cloud_residual(fit_species, radius_effective, data, &
      i_fit_type, property, x1 &
      )
    f_second=cloud_residual(fit_species, radius_effective, data, &
      i_fit_type, property, x2 &
      )
    f_third=cloud_residual(fit_species, radius_effective, data, &
      i_fit_type, property, x3 &
      )
!
!   Increase s_new if minimum is not enclosed.
    IF ( (f_first >= f_second) .OR. (f_first.ge.f_third) ) THEN
      IF ( (f_second < f_third) .AND. (f_first >= f_second) ) THEN
        parm=x1
      ELSE IF (f_first >= f_third) THEN
        parm=x3
      ENDIF
      IF (abs(s_new/s_initial) < s_max_ratio) THEN
        s_new=1.5_RealK*s_new
      ELSE
        WRITE(iu_monitor, '(/a)') &
          '*** Warning: Line search range has been exceeded.'
        ierr=i_abort_calculation
        RETURN
      ENDIF
    ELSE
      EXIT
    ENDIF
!
  ENDDO
!
! The minimum is bracketed between -s and +s. Locate it 
! using Brent's method.
!
! Initialize.
  s_high=s_new
  s_low=-s_new
  s_first=0.0_RealK
  s_second=s_first
  s_third=s_first
  s_new=0.0_RealK
  f_second=f_first
  f_third=f_first
!
! Set the iteration count to 0.
  iteration=0
!
! Test for convergence.
  DO
!
    tol_scaled=(1.0_RealK/h_l2_norm+abs(s_first))*sq_tol_conv
  IF ( (s_first < s_low+tol_scaled) .AND. &
       (s_first > s_high-tol_scaled) ) THEN
!   If convergence has occurred we reset the minimum if different
!   from the original point and finish.
    IF (iteration >= 1) THEN
!     Set the new values of the parameters.
      parm=x_new
!     Return the new step as the searching distance.
      search_distance=abs(s_new*h_l2_norm)
    ENDIF
    EXIT
    ELSE IF (iteration >= np_max_line_search) THEN
      WRITE(iu_monitor, '(/a)') &
        '*** Warning in subroutine line_search_cloud_90.'
      WRITE(iu_monitor, '(a)') &
        'Too many iterations have occurred.'
      ierr=i_abort_calculation
!     Set the new value of the scaling array.
      parm=x_new
!     Return the new step as the searching distance.
      search_distance=ABS(s_new*h_l2_norm)
      EXIT
    ENDIF
!
!   Another iteration if convergence has not occurred.
    iteration=iteration+1
    IF (iteration == 1) THEN
!     On the first step set a large value so as not to reject
!     parabolic interpolation unfairly.
      penult_step=s_high-s_low
    ELSE
      penult_step=step
    ENDIF
!
!   Find mid-point
    s_middle=0.5_RealK*(s_low+s_high)
!
!   Fit a parabola to the lowest three points.
    CALL fit_parabola_90(s_first, s_second, s_third, &
      f_first, f_second, f_third, &
      s0, l_fit)
!
!   If the fit succeeded test that it satisfies the conditions.
    IF ( (l_fit) .AND.           &
         (s_low < s0) .AND.      &
         (s_high > s0) .AND.     &
         (abs(s0-s_first) < 0.5_RealK*abs(penult_step)) ) THEN
!     The step is acceptable.
      s_new=s0
      step=s_new-s_first
    ELSE
!     A golden section search is used.
      IF (s_first >= s_middle) THEN
        step=s_low-s_first
      ELSE
        step=s_high-s_first
      ENDIF
      step=golden_ratio*step
      s_new=s_first+step
    ENDIF
!
!   Functional evaluation at the new point.
    x_new=parm+s_new*h
    f_new=cloud_residual(fit_species, radius_effective, data, &
      i_fit_type, property, x_new &
      )
!
!   Rearrange the points depending on the value of f_new.
    IF (f_new <= f_first) THEN
!     There is a new sequence of the lowest points, and the
!     size of the enclosing interval can be reduced.
      IF (s_new >= s_first) THEN
        s_low=s_first
      ELSE
        s_high=s_first
      ENDIF
      s_third=s_second
      f_third=f_second
      s_second=s_first
      f_second=f_first
      s_first=s_new
      f_first=f_new
!
    ELSE
!     The interval can at least be refined.
      IF (s_new < s_first) THEN
        s_low=s_new
      ELSE
        s_high=s_new
      ENDIF
      IF ( (f_new <= f_second) .OR. &
           (ABS(s_second-s_first) < sq_tol_conv) ) THEN
        s_third=s_second
        f_third=f_second
        s_second=s_new
        f_second=f_new
      ELSE IF ( (f_new <= f_third) .OR. &
                (ABS(s_third-s_first) < sq_tol_conv) .OR. &
                (ABS(s_third-s_second) < sq_tol_conv) ) THEN
        s_third=s_new
        f_third=f_new
      ENDIF
    ENDIF
!
  ENDDO
!
!
!
END SUBROUTINE line_search_cloud_90
