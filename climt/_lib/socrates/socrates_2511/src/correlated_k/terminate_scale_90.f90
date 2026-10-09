! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to print information at end of iteration.
SUBROUTINE terminate_scale_90 &
!
(iu_monitor, &
 n_set, length_set, u_set, trans_set, weight_point, &
 i_type_residual, scaling, &
 n_k, w_k, k_abs, &
 i_scale_function, scale_vector, p_array, t_array, &
 scale_vector_long &
 )
!
! Description:
!   Straightforward.
!
!
! Modules used:
  USE realtype_rd
  USE rad_pcf
  USE def_std_io_icf
!
!
  IMPLICIT NONE
!
!
! Include header files.
!
! Dummy arguments.
  INTEGER, Intent(IN) :: iu_monitor
!   Unit number for output of detailed monitoring information
  INTEGER, Intent(IN) :: n_set
!   Number of sets of data
  INTEGER, Intent(IN), Dimension(:) :: length_set
!   Lengths of sets
  INTEGER, Intent(IN) :: n_k
!   Number of terms of fit
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  INTEGER, Intent(IN) :: i_type_residual
!   Fit to the scaling function
  REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
!   Amounts of absorbers
  REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
!   Transmissivities
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
!   Scalings required
  REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
!   Weights for transmissivities
  REAL  (RealK), Intent(IN), Dimension(:) :: w_k
!   Weights for k-terms
  REAL  (RealK), Intent(IN), Dimension(:) :: k_abs
!   Exponents for k-terms
  REAL  (RealK), Intent(IN), Dimension(:) :: scale_vector
!   Scaling parameters
  REAL  (RealK), Intent(IN), Dimension(:) :: p_array
!   Pressure array
  REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!   Temperature array
  REAL  (RealK), Intent(OUT), Dimension(:, :) :: scale_vector_long
!   Output scaling parameters
!
! Local variables.
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  REAL  (RealK), Allocatable, Dimension(:) :: scale_vector_test
!   Temporary scaling variables
  REAL  (RealK) :: residual_test
!   Value of residual at point
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
  ALLOCATE(scale_vector_test(SIZE(scale_vector)))
!
! Display rms residual errors for different choices of scaling
! parameters.
  WRITE(iu_monitor, '(/a)') 'Final scaling calculated:'
  SELECT CASE(i_scale_function)
!
    CASE (IP_scale_power_law) 
      WRITE(iu_monitor, '(6x, a, 1pe10.3, /, 6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector(1), &
        'Exponent of temperature = ', scale_vector(2)
!
    CASE (ip_scale_power_quad)
      WRITE(iu_monitor, '(6x, a, 1pe10.3, / 6x, a, 2(1pe10.3, 2x))') &
        'Exponent of pressure = ', scale_vector(1), &
        'Coefficients of quadratic = ', scale_vector(2:3)
!
    CASE (ip_scale_doppler_quad)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Offset to normalized pressure = ', scale_vector(2)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector(3:4)
!
  END SELECT
!
  residual_test=residual_trans_90(n_set, length_set, &
    u_set, trans_set, weight_point, &
    i_type_residual, scaling, &
    n_k, w_k, k_abs, &
    i_scale_function, scale_vector, p_array, t_array)
  residual_test=SQRT(residual_test)
  WRITE(iu_monitor, '(6x, a, 1pe12.5)') &
    'R.m.s. residual error in transmission = ', residual_test
!
!
  WRITE(iu_monitor, '(/a)') 'Pressure scaling alone:'
!
  SELECT CASE(i_scale_function)
    CASE (ip_scale_power_law)
      scale_vector_test(1:2) = (/ scale_vector(1), 0.0_RealK /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of temperature = ', scale_vector_test(2)
!
    CASE (ip_scale_power_quad)
      scale_vector_test(1:3) = (/ scale_vector(1), 0.0_RealK, 0.0_RealK /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector_test(2:3)
!
    CASE (ip_scale_doppler_quad)
      scale_vector_test(1:4)= (/ scale_vector(1:2), 0.0_RealK, 0.0_RealK /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Offset to normalized pressure = ', scale_vector_test(2)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector_test(3:4)
!
  END SELECT
!
  residual_test=residual_trans_90(n_set, length_set, &
    u_set, trans_set, weight_point, &
    i_type_residual, scaling, &
    n_k, w_k, k_abs, &
    i_scale_function, scale_vector_test, p_array, t_array)
  residual_test=sqrt(residual_test)
  WRITE(iu_monitor, '(6x, a, 1pe12.5)') &
    'R.m.s. residual error in transmission = ', residual_test
!
!
  WRITE(iu_monitor, '(/a)') 'Temperature scaling alone:'
!
  SELECT CASE(i_scale_function)
!
    CASE (ip_scale_power_law)
      scale_vector_test(1:2) = (/ 0.0_RealK, scale_vector(2) /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of temperature = ', scale_vector_test(2)
!
    CASE (ip_scale_power_quad)
      scale_vector_test(1:3) = (/ 0.0e+00_RealK, scale_vector(2:3) /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector_test(2:3)
!
    CASE (ip_scale_doppler_quad)
      scale_vector_test(1:4) = (/ 0.0_RealK, 0.0_RealK, scale_vector(3:4) /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
          'Offset to normalized pressure = ', scale_vector_test(2)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector_test(3:4)
!
  END SELECT
!
  residual_test=residual_trans_90(n_set, length_set, &
    u_set, trans_set, weight_point, &
    i_type_residual, scaling, &
    n_k, w_k, k_abs, &
    i_scale_function, scale_vector_test, p_array, t_array)
  residual_test=sqrt(residual_test)
  WRITE(iu_monitor, '(6x, a, 1pe12.5)') &
    'R.m.s. residual error in transmission = ', residual_test
!
  WRITE(iu_monitor, '(/a)') 'No scaling at all:'
  SELECT CASE(i_scale_function)
!
    CASE (ip_scale_power_law)
      scale_vector_test(1:2) = (/ 0.0_RealK, 0.0_RealK /)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of temperature = ', scale_vector_test(2)
!
    CASE (ip_scale_power_quad)
      scale_vector_test(1:3)=SPREAD(0.0_RealK, 1, 3)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector_test(2:3)
!
    CASE (ip_scale_doppler_quad)
      scale_vector_test(1:4)=SPREAD(0.0_RealK, 1, 4)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Exponent of pressure = ', scale_vector_test(1)
      WRITE(iu_monitor, '(6x, a, 1pe10.3)') &
        'Offset to normalized pressure = ', scale_vector_test(2)
      WRITE(iu_monitor, '(6x, a, 3(1pe10.3, 2x))') &
        'Coefficients of quadratic = ', scale_vector_test(3:4)
!
  END SELECT
!
  residual_test=residual_trans_90(n_set, length_set, &
    u_set, trans_set, weight_point, &
    i_type_residual, scaling, &
    n_k, w_k, k_abs, &
    i_scale_function, scale_vector_test, p_array, t_array)
  residual_test=sqrt(residual_test)
  WRITE(iu_monitor, '(6x, a, 1pe12.5)') &
    'R.m.s. residual error in transmission = ', residual_test
!
! Copy the short scaling vector into the long array.
  DO i=1, n_k
    scale_vector_long(:, i)=scale_vector(:)
  ENDDO
!
  DEALLOCATE(scale_vector_test)
!
!
!
  RETURN
END SUBROUTINE terminate_scale_90
