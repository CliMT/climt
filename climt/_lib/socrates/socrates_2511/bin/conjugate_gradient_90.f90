! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to determine scaling by conjugate gradients.
!
SUBROUTINE conjugate_gradient_90 &
!
(ierr, iu_monitor, l_re_entry, &
 n_set, length_set, u_set, trans_set, &
 i_type_residual, scaling, weight_point, &
 n_term, w_esft, k_esft, &
 i_scale_function, scale_vector_long, p_array, t_array, &
 rms_residual)
!
! Description:
!   The Polak-Ribiere form of the conjugate gradient algorithm
!   is implemented to find a scaling function for the amount
!   of absorber.
!
!
! Modules used:
  USE realtype_rd
  USE def_std_io_icf
  USE dimensions_spec_ucf
  USE dimensions_pp_ucf
  USE scale_parameters_acf
  USE rad_pcf
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
  INTEGER, Intent(IN) :: iu_monitor
!   Unit number for output to file of detailed information
  INTEGER, Intent(IN) :: n_set
!   Number of sets of data
  INTEGER, Intent(IN) :: n_term
!   Number of ESFT terms
  INTEGER, Intent(IN), Dimension(:) :: length_set
!   Lengths of sets
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  INTEGER, Intent(IN) :: i_type_residual
!   Type of residual used
  LOGICAL, Intent(IN) :: l_re_entry
!   Flag for re-entry
  REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
!   Amounts of absorber
  REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
!   Transmissions
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
!   Scalings required
  REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
!   Weights for data points
  REAL  (RealK), Intent(IN), Dimension(:) :: w_esft
!   ESFT weights
  REAL  (RealK), Intent(IN), Dimension(:) :: k_esft
!   ESFT exponents
  REAL  (RealK), Intent(IN), Dimension(:) :: p_array
!   Pressure array
  REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!   Temperature array
  REAL  (RealK), Intent(INOUT), Dimension(:, :) :: scale_vector_long
!   Parameters for scaling
  REAL  (RealK), Intent(OUT) :: rms_residual
!   R.m.s. residual error in transmission
!
! Local arguments.
  INTEGER :: iteration
!   Iteration count
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: n_var
!   Number of scaling variables
  REAL  (RealK), Allocatable, Dimension(:) :: gradient
!   Gradient of residual
  REAL  (RealK), Allocatable, Dimension(:) :: scale_vector
!   Parameters for scaling
  REAL  (RealK) :: residual_old
!   Old value of residual
  REAL  (RealK) :: residual_new
!   New value of residual
  REAL  (RealK) :: residual_typical
!   Typical value of residual
  REAL  (RealK), Allocatable, Dimension(:) :: g
!   Reversed gradient
  REAL  (RealK), Allocatable, Dimension(:) :: h
!   Direction of search
  REAL  (RealK) :: h_l2_norm
!   L2 norm of h
  REAL  (RealK) :: dgi
!   Temporary constant
  REAL  (RealK) :: gi2
!   Temporary constant
  REAL  (RealK) :: search_distance
!   Searching distance
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
    REAL  (RealK) :: residual_trans
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
!
    SUBROUTINE residual_gradient_90 &
!
      (n_set, length_set, &
       u_set, trans_set, weight_point, &
       i_type_residual, scaling, &
       n_term, w_esft, k_esft, &
       i_scale_function, ScalePrm, p_array, t_array, &
       gradient)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: n_set
      INTEGER, Intent(IN) :: n_term
      INTEGER, Intent(IN), Dimension(:) :: length_set
      INTEGER, Intent(IN) :: i_scale_function
      INTEGER, Intent(IN) :: i_type_residual
      REAL  (RealK), Intent(OUT), Dimension(:) :: gradient
      REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
      REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
      REAL  (RealK), Intent(IN), Dimension(:) :: w_esft  
      REAL  (RealK), Intent(IN), Dimension(:) :: k_esft
      REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
      REAL  (RealK), Intent(IN), Dimension(:) :: p_array
      REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!
    END SUBROUTINE residual_gradient_90
!
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
      USE realtype_rd
!
      INTEGER, Intent(INOUT) :: ierr
      INTEGER, Intent(IN) :: iu_monitor
      INTEGER, Intent(IN) :: n_set
      INTEGER, Intent(IN), Dimension(:) :: length_set
      INTEGER, Intent(IN) :: n_k
      INTEGER, Intent(IN) :: i_scale_function
      INTEGER, Intent(IN) :: i_type_residual
      REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
      REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
      REAL  (RealK), Intent(IN), Dimension(:) :: w_k  
      REAL  (RealK), Intent(IN), Dimension(:) :: k_abs
      REAL  (RealK), Intent(IN), Dimension(:) :: p_array
      REAL  (RealK), Intent(IN), Dimension(:) :: t_array
      REAL  (RealK), Intent(IN), Dimension(:) :: h
      REAL  (RealK), Intent(IN) :: h_l2_norm
      REAL  (RealK), Intent(INOUT) :: search_distance
      REAL  (RealK), Intent(INOUT), Dimension(:) :: ScalePrm
!
    END SUBROUTINE line_search_90
!
!
    SUBROUTINE terminate_scale_90 &
!
      (iu_monitor, n_set, length_set, &
       u_set, trans_set, weight_point, &
       i_type_residual, scaling, &
       n_k, w_k, k_abs, &
       i_scale_function, scale_vector, p_array, t_array, &
       scale_vector_long &
       )
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: iu_monitor
      INTEGER, Intent(IN) :: n_set
      INTEGER, Intent(IN), Dimension(:) :: length_set
      INTEGER, Intent(IN) :: n_k
      INTEGER, Intent(IN) :: i_scale_function
      INTEGER, Intent(IN) :: i_type_residual
      REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
      REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
      REAL  (RealK), Intent(IN), Dimension(:) :: w_k
      REAL  (RealK), Intent(IN), Dimension(:) :: k_abs
      REAL  (RealK), Intent(IN), Dimension(:) :: scale_vector
      REAL  (RealK), Intent(IN), Dimension(:) :: p_array
      REAL  (RealK), Intent(IN), Dimension(:) :: t_array
      REAL  (RealK), Intent(OUT), Dimension(:, :) :: scale_vector_long
!
    END SUBROUTINE terminate_scale_90
!
!
  END INTERFACE
!
!
!
  n_var=SIZE(scale_vector_long, 1)
  ALLOCATE(scale_vector(n_var))
  ALLOCATE(gradient(n_var))
  ALLOCATE(g(n_var))
  ALLOCATE(h(n_var))
!
! Initialize the scaling vector.
  scale_vector(:)=scale_vector_long(:, 1)
!
! Find a typical value of the residual
  residual_typical=residual_trans_90(n_set, length_set, &
    u_set, trans_set, weight_point, &
    i_type_residual, scaling, &
    n_term, w_esft, k_esft, &
    i_scale_function, &
    SPREAD(0.0_RealK, 1, n_var), &
    p_array, t_array)
!
  WRITE(iu_monitor, '(//a)') 'Attempting to find optimal scaling '
  WRITE(iu_monitor, '(a/)') &
    'for pressure and temperature using the conjugate gradient algorithm,'
!
! Initialize the iterating loop.
  iteration=0
  residual_old=residual_trans_90(n_set, length_set, &
    u_set, trans_set, weight_point, &
    i_type_residual, scaling, &
    n_term, w_esft, k_esft, &
    i_scale_function, scale_vector, p_array, t_array)
  CALL residual_gradient_90(n_set, length_set, &
     u_set, trans_set, weight_point, &
     i_type_residual, scaling, &
     n_term, w_esft, k_esft, &
     i_scale_function, scale_vector, p_array, t_array, &
     gradient)
  g=-gradient
  h=g
  rms_residual=SQRT(residual_old)
!
  DO
    iteration=iteration+1
!
    h_l2_norm=0.0_RealK
    DO i=1, SIZE(h)
      h_l2_norm=h_l2_norm+h(i)*h(i)
    ENDDO
    h_l2_norm=SQRT(h_l2_norm)
!
!   If this is very small further minimization is pointless.
    IF (h_l2_norm < SQRT(epsilon(h_l2_norm))*residual_typical) THEN
      CALL terminate_scale_90(iu_monitor, n_set, length_set, &
        u_set, trans_set, weight_point, &
        i_type_residual, scaling, &
        n_term, w_esft, k_esft, &
        i_scale_function, scale_vector, p_array, t_array, &
        scale_vector_long)
      IF (.NOT.l_re_entry .AND. (iteration == 1)) THEN
!       If the program exits without performing a minimization
!       and this is not a run to refine an existing set of
!       scaling parameters they are set to 0.
        DO i=1, n_scale_variable(i_scale_function)
          scale_vector_long(i, 1:n_term)=SPREAD(0.0_RealK, 1, n_term)
        ENDDO
        rms_residual=SQRT(residual_typical)
      ENDIF
      EXIT
    ENDIF
!
    IF (iteration == 1) THEN
      search_distance=1.0_RealK
    ENDIF
!
!   Search along H for a minimum.
    search_distance=1.0_RealK
    CALL line_search_90(ierr, iu_monitor, &
      n_set, length_set, u_set, trans_set, &
      weight_point, &
      i_type_residual, scaling, &
      n_term, w_esft, k_esft, &
      i_scale_function, scale_vector, p_array, t_array, &
      h, h_l2_norm, search_distance)
    IF (ierr /= i_normal) THEN
      IF (ierr == i_abort_calculation) THEN
        WRITE(*, '(/a)') 'Attempt to find scaling abandoned.'
        WRITE(iu_monitor, '(/a)') 'Attempt to find scaling abandoned.'
        rms_residual=SQRT(residual_typical)
      ENDIF
      EXIT
    ENDIF
!
!   Assess the condition for normal terminiation
    residual_new=residual_trans_90(n_set, length_set, &
      u_set, trans_set, weight_point, &
      i_type_residual, scaling, &
      n_term, w_esft, k_esft, &
      i_scale_function, scale_vector, p_array, t_array)
    rms_residual=SQRT(residual_new)
    WRITE(iu_monitor, '(/a, i4)') 'Iteration = ', iteration
    WRITE(iu_monitor, '(a, 1pe16.9)') &
      'Mean square error = ', residual_new
    WRITE(iu_monitor, '(a, (t22, 3(1pe16.9, 2x)))') &
      'scaling parameters = ', &
      scale_vector(1:n_scale_variable(i_scale_function))
    IF (ABS(residual_new-residual_old) < &
        EPSILON(residual_old)*residual_typical) THEN
!     We are done.
      CALL terminate_scale_90(iu_monitor, n_set, length_set, &
        u_set, trans_set, weight_point, &
        i_type_residual, scaling, &
        n_term, w_esft, k_esft, &
        i_scale_function, scale_vector, p_array, t_array, &
        scale_vector_long)
      RETURN
    ELSE
      residual_old=residual_new
    ENDIF
!
!   Proceed to determine the conjugate direction
!   find the gradient at the new minimum.
    CALL residual_gradient_90(n_set, length_set, &
      u_set, trans_set, weight_point, &
      i_type_residual, scaling, &
      n_term, w_esft, k_esft, &
      i_scale_function, scale_vector, p_array, t_array, &
      gradient)
    gi2=0.0_RealK
    dgi=0.0_RealK
    DO i=1, n_scale_variable(i_scale_function)
      gi2=gi2+g(i)*g(i)
      dgi=dgi+gradient(i)*(gradient(i)+g(i))
    ENDDO
    IF (gi2 < EPSILON(residual_typical)*residual_typical) THEN
      CALL terminate_scale_90(iu_monitor, n_set, length_set, &
        u_set, trans_set, weight_point, &
        i_type_residual, scaling, &
        n_term, w_esft, k_esft, &
        i_scale_function, scale_vector, p_array, t_array, &
        scale_vector_long)
      EXIT
    ENDIF
!
    g(:)=-gradient(:)
!   In the case of non-quadratic minima the direction of
!   searching should be re-initialized after a number of steps
!   equal to the number of variables, since otherwise the
!   searching directions may drift away from the optimal ones.
    IF (MOD(iteration, n_scale_variable(i_scale_function)) == 0) THEN
      search_distance=1.0_RealK
      h(:)=g(:)
    ELSE
      h(:)=g(:)+(dgi/gi2)*h(:)
    ENDIF
!
    IF (iteration >= np_max_iteration_cg) THEN
      WRITE(*, '(/a)') &
        'The conjugate gradient scaling has not converged:'
      WRITE(*, '(/a)') &
        'Current values of the scaling parameters will be given.'
      WRITE(iu_monitor, '(/a)') 'Too many iterations without convergence.'
      WRITE(iu_monitor, '(a)') 'Current data will be used.'
      CALL terminate_scale_90(iu_monitor, n_set, length_set, &
        u_set, trans_set, weight_point, &
        i_type_residual, scaling, &
        n_term, w_esft, k_esft, &
        i_scale_function, scale_vector, p_array, t_array, &
        scale_vector_long)
      EXIT
    ENDIF
!
  ENDDO
!
  DEALLOCATE(scale_vector)
  DEALLOCATE(gradient)
  DEALLOCATE(g)
  DEALLOCATE(h)
!
!
!
END SUBROUTINE conjugate_gradient_90
