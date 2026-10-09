! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate scaling for a k-distribution.
SUBROUTINE scale_ck_fit_90 &
!
(ierr, iu_monitor, &
 n_other, &
 n_path, u_set, trans_set, p_set, t_set, &
 p_ref, t_ref, k_ref, &
 i_type_residual, i_scale_function, scale_vector, &
 rms_residual)
!
! Description:
!   An initial guess is made for the scaling parameters and
!   scaled independent variables are defined for numerical 
!   convenience. A conjugate gradient algorithm is then called.
!
!
!
! Modules used.
  USE realtype_rd
  USE def_std_io_icf
  USE rad_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
  INTEGER, Intent(IN) ::  iu_monitor
!   Unit number for output of detailed monitoring information
  INTEGER, Intent(IN) ::  i_scale_function
!   Type of scaling function
  INTEGER, Intent(IN) :: i_type_residual
!   Type of residual used
  INTEGER, Intent(IN) :: n_other
!   Number pairs of p and T, excluding the reference one
  INTEGER, Intent(IN) :: n_path
!   Number of path-lengths
!
  REAL  (RealK), Intent(IN) :: p_ref
!   Reference pressure
  REAL  (RealK), Intent(IN) :: t_ref
!   Reference temperature
  REAL  (RealK), Intent(IN) :: k_ref
!   Reference k-value
  REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
!   Amounts of absorbers
  REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
!   Transmissions
  REAL  (RealK), Intent(IN), Dimension(:) :: p_set
!   Pressures at which data are given
  REAL  (RealK), Intent(IN), Dimension(:) :: t_set
!   Temperatures at which data are given
  REAL  (RealK), Intent(OUT), Dimension(:) :: scale_vector
!   Scaling parameters
  REAL  (RealK), Intent(OUT) :: rms_residual
!   R.m.s. residual error in transmission
!
! Local variables.
  INTEGER :: n_point
!   Total number values of transmission for all pressures, 
!   temperatures and amounts
  INTEGER :: length_set(n_other)
!   Number of data points in each set of tranmissions
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  REAL  (RealK) :: p_norm(n_other)
!   Normalized and converted pressures
  REAL  (RealK) :: t_norm(n_other)
!   Normalized and converted temperatures
  REAL  (RealK) :: weight_point(n_other, n_path)
!   Weights for data points
  REAL  (RealK) :: scaling(n_other, n_path)
!   Scalings required
  REAL  (RealK) :: range_t
!   Range of temperatures encountered
  REAL  (RealK), Allocatable, Dimension(:, :) :: scale_vector_long
!   Scaling parameters: a local vector of dimensions to match the
!   lower routine
!
! Subroutines called:
  INTERFACE
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
      USE realtype_rd
!
      INTEGER, Intent(INOUT) :: ierr
      INTEGER, Intent(IN) :: iu_monitor
      INTEGER, Intent(IN) :: n_set
      INTEGER, Intent(IN) :: n_term
      INTEGER, Intent(IN), Dimension(:) :: length_set
      INTEGER, Intent(IN) :: i_scale_function
      INTEGER, Intent(IN) :: i_type_residual
      LOGICAL, Intent(IN) :: l_re_entry
      REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
      REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
      REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
      REAL  (RealK), Intent(IN), Dimension(:) :: w_esft
      REAL  (RealK), Intent(IN), Dimension(:) :: k_esft
      REAL  (RealK), Intent(IN), Dimension(:) :: p_array
      REAL  (RealK), Intent(IN), Dimension(:) :: t_array
      REAL  (RealK), Intent(INOUT), Dimension(:, :) :: scale_vector_long
      REAL  (RealK), Intent(OUT) :: rms_residual
!
    END SUBROUTINE conjugate_gradient_90
!
  END INTERFACE
!
!
!
! Normalize the pressure and the temperature to a convenient form 
! for fitting.
  IF ( (i_scale_function == IP_scale_power_quad) .OR. &
       (i_scale_function == IP_scale_doppler_quad) ) THEN
    range_t=0.0_RealK
    DO i=1, n_other
      range_t=MAX(range_t, ABS(t_set(i)-t_ref))
    ENDDO
  ENDIF
  DO i=1, n_other
    IF (i_scale_function == IP_scale_power_law) THEN
      p_norm(i)=LOG(p_set(i)/p_ref)
      t_norm(i)=LOG(t_set(i)/t_ref)
    ELSE IF (i_scale_function == IP_scale_power_quad) THEN
      p_norm(i)=LOG(p_set(i)/p_ref)
      t_norm(i)=(t_set(i)-t_ref)/range_t
    ELSE IF (i_scale_function == IP_scale_doppler_quad) THEN
      p_norm(i)=p_set(i)/p_ref
      t_norm(i)=(t_set(i)-t_ref)/range_t
    ENDIF
  ENDDO
!
! Initialize the scaling parameters
  IF (i_scale_function == IP_scale_power_law) THEN
    scale_vector(1:2) = (/ 1.0_RealK, 0.5_RealK /)
  ELSE IF (i_scale_function == IP_scale_power_quad) THEN
    scale_vector(1:3) = (/ 1.0_RealK, 0.0_RealK, 0.0_RealK /)
  ELSE IF (i_scale_function == IP_scale_doppler_quad) THEN
    scale_vector(1:4) = (/ 1.0_RealK, -2.0_RealK, 0.0_RealK, 0.0_RealK /)
  ELSE
    WRITE(iu_err, '(/a)') &
      '*** Error: An unsupported scaling function has been requested.'
    ierr=i_err_fatal
    RETURN
  ENDIF
!
! Arrays for the conjugate gradient algorithm: this is more general
! than is required at this point because of its earlier heritage,
! but this functionality may again become useful.
! Here all points are weighted equally.
  n_point=0
  DO i=1, n_other
    length_set(i)=n_path
    n_point=n_point+length_set(i)
  ENDDO
  DO i=1, n_other
    DO j=1, length_set(i)
      weight_point(i, j)=1.0_RealK/real(n_point, RealK)
    ENDDO
  ENDDO
!
! Pass to the conjugate gradient algorithm. This retains the functionality
! to operate over several k-terms, so scale_vector is embedded in a 
! larger array to match.
  ALLOCATE(scale_vector_long(SIZE(scale_vector),1))
  scale_vector_long(:, 1) = scale_vector(:)
  CALL conjugate_gradient_90(ierr, iu_monitor, .FALSE., &
    n_other, length_set, u_set, trans_set, &
    i_type_residual, scaling, weight_point, &
    1, (/ 1.0_RealK /), (/ k_ref /), &
    i_scale_function, scale_vector_long, p_norm, t_norm, &
    rms_residual)
  IF ( (ierr /= i_normal) .AND. &
       (ierr /= i_abort_calculation) ) RETURN
  scale_vector(:)=scale_vector_long(:, 1)
  DEALLOCATE(scale_vector_long)
!
!
  IF (ierr == i_abort_calculation) THEN
!   If the calculation has failed we set all scaling 
!   parameters to 0.
    IF (i_scale_function == IP_scale_power_law) THEN
      scale_vector(1:2) = 0.0_RealK
    ELSE IF (i_scale_function == IP_scale_power_quad) THEN
      scale_vector(1:3) = 0.0_RealK
    ELSE IF (i_scale_function == IP_scale_doppler_quad) THEN
!     The second element is later to be exponentiated 
!     so 0 is inappropriate.
      scale_vector(1:4) = &
        (/ 0.0_RealK, -1.0e+02_RealK, 0.0_RealK, 0.0_RealK /)
    ENDIF
    RETURN
  ENDIF
!
!
! Normalize the scaling variables to the correct form
! for the radiation code.
  IF (i_scale_function == IP_scale_power_quad) THEN
    scale_vector(2)=t_ref*scale_vector(2)/range_t
    scale_vector(3)=scale_vector(3)*(t_ref/range_t)**2
  ENDIF
  IF (i_scale_function == IP_scale_doppler_quad) THEN
    scale_vector(2)=EXP(scale_vector(2))*p_ref
    scale_vector(3)=t_ref*scale_vector(3)/range_t
    scale_vector(4)=scale_vector(4)*(t_ref/range_t)**2
  ENDIF
!
!
!
  RETURN
END SUBROUTINE scale_ck_fit_90
