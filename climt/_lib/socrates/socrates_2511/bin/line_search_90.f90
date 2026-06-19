! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to perform scalar minimization in a fixed direction. 
!
SUBROUTINE line_search_90(ierr, iu_monitor, &
  n_set, length_set, &
  u_set, trans_set, weight_point, &
  i_type_residual, scaling, &
  n_k, w_k, k_abs, &
  i_scale_function, ScalePrm, p_array, t_array, &
  h, h_l2_norm, search_distance &
  )
!
! Description:
!   The value of a parameter, S, giving the distance travelled 
!   along a direction, H, is adjusted to minimize a function F.
!   The minimum is first bracketed and then Brent's method is 
!   applied to find the minimum. This routine is specific to the
!   minimization of errors in fitting gaseous transmissions.
!
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE scale_parameters_acf
  USE rad_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(IN) :: iu_monitor
!   Unit number for output of detailed monitoring information
  INTEGER, Intent(IN) :: n_set
!   Number of sets
  INTEGER, Intent(IN), Dimension(:) :: length_set
!   Lengths of sets
  INTEGER, Intent(IN) :: n_k
!   Number of k-terms
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  INTEGER, Intent(IN) :: i_type_residual
!   Fit to the scaling function
  REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
!   Amounts of absorber
  REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
!   Transmissivities
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
!   Scalings required
  REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
!   Weights for data points
  REAL  (RealK), Intent(IN), Dimension(:) :: w_k
!   Weights for k-terms
  REAL  (RealK), Intent(IN), Dimension(:) :: k_abs
!   k-terms
  REAL  (RealK), Intent(IN), Dimension(:) :: p_array
!   Pressure array
  REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!   Temperature array
  REAL  (RealK), Intent(IN), Dimension(:) :: h
!   Direction of searching
  REAL  (RealK), Intent(IN) :: h_l2_norm
!   L2 norm of h.
  REAL  (RealK), Intent(INOUT) :: search_distance
!   Esitmated searching distance
  REAL  (RealK), Intent(INOUT), Dimension(:) :: ScalePrm
!   Parameters for scale function
!
! Local variables.
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
  REAL  (RealK), Allocatable, Dimension(:) :: x1
!   Search position
  REAL  (RealK), Allocatable, Dimension(:) :: x2
!   Search position
  REAL  (RealK), Allocatable, Dimension(:) :: x3
!   Search position
  REAL  (RealK), Allocatable, Dimension(:) :: x_new
!   New search position
!
  REAL  (RealK), Parameter :: golden_ratio = 0.61803399_RealK
!   Golden ratio for search
!
! Variables related to the treatment of ill-conditioning
  REAL  (RealK) :: tol_conv
!   The tolerance used to assess convergence
  REAL  (RealK) :: sq_tol_conv
!   The square root of the above
!
! Subroutines called:
  EXTERNAL fit_parabola_90
!
! Functions called:
  INTERFACE
!
    FUNCTION residual_trans_90 (n_set, length_set, u_set, trans_set, &
      weight_point, i_type_residual, scaling, &
      n_term, w_k, k_abs, i_scale_function, ScalePrm, p_array, t_array) &
    RESULT(residual)
!
    USE realtype_rd
!
    INTEGER, Intent(IN) :: n_set
    INTEGER, Intent(IN) :: n_term
    INTEGER, Intent(IN), Dimension(:) :: length_set
    INTEGER, Intent(IN) :: i_scale_function
    INTEGER, Intent(IN) :: i_type_residual
    REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
    REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
    REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
    REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
    REAL  (RealK), Intent(IN), Dimension(:) :: w_k
    REAL  (RealK), Intent(IN), Dimension(:) :: k_abs
    REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
    REAL  (RealK), Intent(IN), Dimension(:) :: p_array
    REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!
    REAL  (RealK) :: residual
!
    END FUNCTION residual_trans_90
!
  END INTERFACE
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
! Set aside local space.
  ALLOCATE(x1(n_scale_variable(i_scale_function)))
  ALLOCATE(x2(n_scale_variable(i_scale_function)))
  ALLOCATE(x3(n_scale_variable(i_scale_function)))
  ALLOCATE(x_new(n_scale_variable(i_scale_function)))
!
! Find three points enclosing the minimum.
  DO
    x1(:)=ScalePrm(:)
    x2(:)=ScalePrm(:)-s_new*h(:)
    x3(:)=ScalePrm(:)+s_new*h(:)
!
!   Find the residuals at these points.
    f_first=residual_trans_90(n_set, length_set, &
      u_set, trans_set, weight_point, &
      i_type_residual, scaling, &
      n_k, w_k, k_abs, &
      i_scale_function, x1, p_array, t_array)
    f_second=residual_trans_90(n_set, length_set, &
      u_set, trans_set, weight_point, &
      i_type_residual, scaling, &
      n_k, w_k, k_abs, &
      i_scale_function, x2, p_array, t_array)
    f_third=residual_trans_90(n_set, length_set, &
      u_set, trans_set, weight_point, &
      i_type_residual, scaling, &
      n_k, w_k, k_abs, &
      i_scale_function, x3, p_array, t_array)
!
!   Increase s_new if minimum is not enclosed.
    IF ( (f_first >= f_second) .OR. (f_first >= f_third) ) THEN
      IF ( (f_second < f_third) .AND. (f_first >= f_second) ) THEN
        ScalePrm=x1
      ELSE IF (f_first >= f_third) THEN
        ScalePrm=x3
      ENDIF
      IF (ABS(s_new/s_initial) < s_max_ratio) THEN
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
!
! The minimum is bracketed between -s and +s. locate it 
! using Brent's method.
!
! Initialize.
  s_high=s_new
  s_low=-s_new
  s_first=0.0_RealK
  s_second=s_first
  s_third=s_first
  s_new=0.0e+00_RealK
  f_second=f_first
  f_third=f_first
!
! Set the iteration count to 0.
  iteration=0
!
! Test for convergence.
  DO
!
    tol_scaled=(1.0_realk/h_l2_norm+abs(s_first))*sq_tol_conv
    IF ( (s_first < s_low+tol_scaled) .AND. &
         (s_first > s_high-tol_scaled) ) THEN
!     If convergence has occurred we reset the minimum if different
!     from the original point and finish.
      IF (iteration >= 1) THEN
!       Set the new value of ScalePrm.
        ScalePrm=x_new
!       Return the new step as the searching distance.
        search_distance=ABS(s_new*h_l2_norm)
      ENDIF
      EXIT
    ELSE IF (iteration >= np_max_line_search) THEN
      WRITE(iu_monitor, '(/a)') '*** Warning in subroutine line_search.'
      WRITE(iu_monitor, '(a)') 'Too many iterations have occurred.'
      ierr=i_abort_calculation
!     Set the new value of ScalePrm.
      ScalePrm=x_new
!     Return the new step as the searching distance.
      search_distance=ABS(s_new*h_l2_norm)
      EXIT
    ENDIF
!
!   Another iteration if convergence has not occurred.
    iteration=iteration+1
    IF (iteration == 1) THEN
!     Set this to a large value to avoid unfair rejection of
!     parabolic fitting later.
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
         (ABS(s0-s_first) < 0.5_RealK*abs(penult_step)) ) THEN
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
    x_new=ScalePrm+s_new*h
    f_new=residual_trans_90(n_set, length_set, &
      u_set, trans_set, weight_point, &
      i_type_residual, scaling, &
      n_k, w_k, k_abs, &
      i_scale_function, x_new, p_array, t_array)
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
           (ABS(s_second-s_first) < tol_conv) ) THEN
        s_third=s_second
        f_third=f_second
        s_second=s_new
        f_second=f_new
      ELSE IF ( (f_new <= f_third) .OR.                &
                (ABS(s_third-s_first) < tol_conv) .OR. &
                (ABS(s_third-s_second) < tol_conv) ) THEN
        s_third=s_new
        f_third=f_new
      ENDIF
    ENDIF
!
  ENDDO
!
  DEALLOCATE(x1)
  DEALLOCATE(x2)
  DEALLOCATE(x3)
  DEALLOCATE(x_new)
!
!
!
END
