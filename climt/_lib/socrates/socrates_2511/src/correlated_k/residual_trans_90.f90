! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the residuals in the transmission.
!
FUNCTION residual_trans_90 &
!
(n_set, length_set, &
 u_set, trans_set, weight_point, &
 i_type_residual, scaling, &
 n_term, w_k, k_abs, &
 i_scale_function, ScalePrm, p_array, t_array &
 ) &
!
RESULT(residual)
!
! Description:
!   Straightforward.
!
!
! Modules used:
  USE realtype_rd
  USE type_residual_pcf
  USE offset_residual_trans_acf
!
!
  IMPLICIT NONE
!
!
! Include header files.
!
! Dummy arguments.
  INTEGER, Intent(IN) :: n_set
!   Number of sets of data
  INTEGER, Intent(IN) :: n_term
!   Number of terms of fit
  INTEGER, Intent(IN), Dimension(:) :: length_set
!   Lengths of sets
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  INTEGER, Intent(IN) :: i_type_residual
!   Type of residual used
  REAL  (RealK) :: residual_trans
!   Residual in transmissivity
  REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
!   Amounts of absorbers
  REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
!   Transmissivities
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
!   Scalings required
  REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
!   Weights for transmissivities
  REAL  (RealK), Intent(IN), Dimension(:) :: w_k
!   ESFT wieghts
  REAL  (RealK), Intent(IN), Dimension(:) :: k_abs
!   ESFT exponents
  REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
!   Scaling parameters
  REAL  (RealK), Intent(IN), Dimension(:) :: p_array
!   Pressure array
  REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!   Temperature array
!
  REAL  (RealK) :: residual
!   Calculated residual in transmission or scaling function
!
!
! Local variables.
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: k
!   Loop variable
  REAL  (RealK) :: trans_fit
!   Fitted transmittance
  REAL  (RealK) :: residual_set
!   Residual term of one set
  REAL  (RealK), Allocatable, Dimension(:) :: scale_factor
!   Scaling factor for each set
!
! Functions called:
  INTERFACE
!
    FUNCTION func_scale_90(p_fnc, t_fnc, i_scale_function, ScalePrm) &
      RESULT(f)
!
    USE realtype_rd
!
    INTEGER, Intent(IN) :: i_scale_function
    REAL  (RealK), Intent(IN) :: p_fnc
    REAL  (RealK), Intent(IN) :: t_fnc
    REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
!
    REAL  (RealK) :: f
!
    END FUNCTION func_scale_90
!
  END INTERFACE
!
!
!
  ALLOCATE(scale_factor(n_set))
  DO i=1, n_set
    scale_factor(i) = func_scale_90(p_array(i), t_array(i), &
                                    i_scale_function, ScalePrm)
    scale_factor(i) = MAX(scale_factor(i), 0.0_RealK)
  ENDDO
  residual=0.0_RealK
!
  SELECT CASE(i_type_residual)
!
    CASE (IP_scale_trans_residual)
      DO i=1, n_set
        residual_set=0.0_RealK
        DO j=1, length_set(i)
          trans_fit=0.0_RealK
          DO k=1, n_term
            trans_fit=trans_fit + &
              w_k(k)*exp(-k_abs(k)*u_set(j, i)*scale_factor(i))
          ENDDO
          residual_set=residual_set + &
            weight_point(i, j)*(trans_set(j, i)-trans_fit)**2
        ENDDO
        residual=residual+residual_set
      ENDDO
!
    CASE (IP_scale_full_scl_res)
      DO i=1, n_set
        residual_set=0.0_RealK
        DO j=1, length_set(i)
          residual_set=residual_set + &
            weight_point(i, j) * &
            (LOG((offset_residual+scale_factor(i)) / &
            (offset_residual+scaling(j, i))))**2
        ENDDO
        residual=residual+residual_set
      ENDDO
!
    CASE(IP_scale_mean_scl_res)
      DO i=1, n_set
        residual=residual+ &
          (((offset_residual+scale_factor(i)) / &
            (offset_residual+scaling(1, i)))-1.0_RealK)**2 + &
          (scale_factor(i)-scaling(1, i))**2
      ENDDO
!
  END SELECT
  DEALLOCATE(scale_factor)
!
!
!
  RETURN
  end
