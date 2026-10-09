! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to fit to cloud properties using conjugate gradients.
!
SUBROUTINE conjugate_gradient_cloud_90(ierr, iu_monitor, &
  fit_species, &
  n_data, d, vol_fraction, actual, &
  i_fit, property, &
  n_parameter, parm &
  )
!
! Method:
!      The Polak-Ribiere from of the conjugate gradient algorithm
!      is implemented to find a fit to cloudy properties.
!
!
!
! Modules used
  USE realtype_rd
  USE cloud_fit_parm_acf
  USE cloud_fitting
  USE def_std_io_icf
  USE error_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(INOUT) :: iu_monitor
!   Unit number for monitoring information
!
  CHARACTER  (LEN=5), Intent(IN)  :: fit_species
!   Phase of condensate to be fitted
!
  INTEGER, Intent(IN) :: n_data
!   Number of data points
  INTEGER, Intent(IN) :: i_fit
!   Type of fit required
  CHARACTER  (LEN=10), Intent(IN) :: property
!   Single-scattering property to be fitted
  INTEGER, Intent(IN) :: n_parameter
!   Number of fitting parameters
  REAL  (RealK), Intent(IN), Dimension(n_data) :: actual
!   Data to be fitted
  REAL  (RealK), Intent(IN), Dimension(n_data) :: vol_fraction
!   Volume fraction of absorber
  REAL  (RealK), Intent(IN), Dimension(n_data) :: d
!   Effective radius of particles
  REAL  (RealK), Intent(INOUT), Dimension(n_parameter) :: parm
!   Fitting parameters
!
!
!
! Local arguments
  INTEGER :: iteration
!   Number of CG iteration
  INTEGER :: i
!   Loop variable
  REAL  (RealK) :: residual_old
!   Former value of residual
  REAL  (RealK) :: residual_new
!   Newer value of residual
  REAL  (RealK) :: residual_best
!   Best value of residual
  REAL  (RealK) :: residual_typical
!   Typical value of residual
  REAL  (RealK), Dimension(n_parameter) :: parm_best
!   Best fitting parameters
  REAL  (RealK), Allocatable, Dimension(:) :: gradient
!   Gradient
  REAL  (RealK), Allocatable, Dimension(:) :: g
!   Reversed gradient
  REAL  (RealK), Allocatable, Dimension(:) :: h
!   Search direction
  REAL  (RealK) :: h_l2_norm
!   L2-norm of h
  REAL  (RealK) :: search_distance
!   Estimated search distance
  REAL  (RealK) :: gi2
!   Temporary constant
  REAL  (RealK) :: dgi
!   Temporary constant
!
! Subroutines called:
  EXTERNAL &
      line_search_cloud_90
!
!
!
  ALLOCATE(gradient(n_parameter))
  ALLOCATE(g(n_parameter))
  ALLOCATE(h(n_parameter))
!
! Calculate the initial value of the residual.
  residual_old=cloud_residual(fit_species, d, actual, &
    i_fit, property, parm &
    )

  residual_best=residual_old
  parm_best=parm

! Take the initial residual as an indication of the typical size.
  residual_typical=residual_old
!
! Initialize the gradient and the direction of searching.
  CALL cloud_residual_gradient(fit_species, d, actual, &
    i_fit, property, parm, gradient &
    )
  g(1:n_parameter)=-gradient(1:n_parameter)
  h(1:n_parameter)=g(1:n_parameter)
!
!
! Initialize the count of iterations.
  iteration=0
  DO
    iteration=iteration+1
!
    h_l2_norm=0.0_RealK
    DO i=1, n_parameter
      h_l2_norm=h_l2_norm+h(i)*h(i)
    ENDDO
    h_l2_norm=SQRT(h_l2_norm)
!
!   Return if h is virtually 0.
    IF (h_l2_norm < SQRT(epsilon(h_l2_norm))*residual_typical) EXIT
!
    IF (iteration == 1) THEN
!     Estimate a suitable searching distance.
      search_distance=1.0_RealK
    ENDIF
!
!   Search along the direction H.
    search_distance=1.0_RealK
    CALL line_search_cloud_90(ierr, iu_monitor, &
      fit_species, &
      n_data, d, vol_fraction, actual, &
      i_fit, property, n_parameter, parm, &
      h, h_l2_norm, search_distance  &
      )
    IF (ierr /= i_normal) THEN
      IF (ierr == i_abort_calculation) THEN
!       Line search was aborted: we recover from this.
        ierr=i_normal
      ELSE
        EXIT
      ENDIF
    ENDIF
!
!   Assess the condition for normal termination.
    residual_new=cloud_residual(fit_species, d, actual, &
      i_fit, property, parm &
      )
    IF (residual_new < residual_best) THEN
      residual_best=residual_new
      parm_best=parm
    END IF
    WRITE(iu_monitor, '(/a, i4)') 'Iteration = ', iteration
    WRITE(iu_monitor, '(a, 1pe16.9)') &
      'Mean square error = ', residual_new
    WRITE(iu_monitor, '(a, (2(1x, 1pe16.9)))') &
      'Scaling parameters = ', (parm(i), i=1, n_parameter)
    IF (ABS(residual_new-residual_old) < &
        EPSILON(residual_typical)*residual_typical) THEN
      parm=parm_best
      EXIT
    ELSE
      residual_old=residual_new
    ENDIF
!
!   Now search in a direction conjugate to the current one.
    CALL cloud_residual_gradient(fit_species, d, actual, &
      i_fit, property, parm, gradient &
      )
    gi2=0.0_RealK
    dgi=0.0_RealK
    DO  i=1, n_parameter
      gi2=gi2+g(i)*g(i)
      dgi=dgi+gradient(i)*(gradient(i)+g(i))
    ENDDO
    IF (gi2 < EPSILON(residual_typical)*residual_typical) EXIT
!
    DO  i=1, n_parameter
      g(i)=-gradient(i)
!     In the case of non-quadratic minima the direction of
!     searching should be re-initialized after a number of steps
!     equal to the number of variables, since otherwise the
!     searching directions may drift away from the optimal ones.
      IF (MOD(iteration, n_parameter) == 0) THEN
        h(i)=g(i)
        search_distance=1.0_RealK
      ELSE
        h(i)=g(i)+(dgi/gi2)*h(i)
      ENDIF
    ENDDO
!
    IF (iteration > np_max_iteration_cg) THEN
      WRITE(iu_stdout, '(/a)') &
        'Too many iterations in the conjugate gradient routine ' // &
        'without convergence.'
      WRITE(iu_stdout, '(a)') 'Best data will be used.'
      parm=parm_best
      EXIT
    ENDIF
!
  ENDDO
!
  DEALLOCATE(gradient)
  DEALLOCATE(g)
  DEALLOCATE(h)
!
!
!
END SUBROUTINE conjugate_gradient_cloud_90
