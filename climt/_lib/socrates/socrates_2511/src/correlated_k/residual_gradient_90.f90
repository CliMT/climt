! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate gradient of residuals in transmission.
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
! Description:
!   Straightforward.
!
!
! Modules used:
  USE realtype_rd
  USE offset_residual_trans_acf
  USE type_residual_pcf
  USE rad_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments.
  INTEGER, Intent(IN) :: n_set
!   Number of sets
  INTEGER, Intent(IN) :: n_term
!   Number of terms of fit
  INTEGER, Intent(IN), Dimension(:) :: length_set
!   Lengths of sets
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  INTEGER, Intent(IN) :: i_type_residual
!   Type of residual used
  REAL  (RealK), Intent(OUT), Dimension(:) :: gradient
!   Gradient of residual
  REAL  (RealK), Intent(IN), Dimension(:, :) :: u_set
!   Amounts of absorbers
  REAL  (RealK), Intent(IN), Dimension(:, :) :: trans_set
!   Transmissivities
  REAL  (RealK), Intent(IN), Dimension(:, :) :: scaling
!   Scalings required
  REAL  (RealK), Intent(IN), Dimension(:, :) :: weight_point
!   Weights for transmissivities
  REAL  (RealK), Intent(IN), Dimension(:) :: w_esft
!   ESFT wieghts
  REAL  (RealK), Intent(IN), Dimension(:) :: k_esft
!   ESFT exponents
  REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
!   Parameters for scaling
  REAL  (RealK), Intent(IN), Dimension(:) :: p_array
!   Pressure array
  REAL  (RealK), Intent(IN), Dimension(:) :: t_array
!   Temperature array
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
  REAL  (RealK), Allocatable, Dimension(:) :: d_residual_by_d_scale
!   Derivative of residual wrt scale
  REAL  (RealK) :: sub_term
!   Sub-term of transmittance
  REAL  (RealK) :: derivative_term
!   Term entering derivative
  REAL  (RealK) :: ss
!   Temporary variable
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
    FUNCTION func_scale_derivative_90(p_fnc, t_fnc, &
                                      i_scale_function, ScalePrm) &
      RESULT(df)
!
    USE realtype_rd
!
    INTEGER, Intent(IN) :: i_scale_function
    REAL  (RealK), Intent(IN) :: p_fnc
    REAL  (RealK), Intent(IN) :: t_fnc
    REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
!
    REAL  (RealK), Pointer, Dimension(:) :: df
!
    END FUNCTION func_scale_derivative_90
!
  END INTERFACE
!
!
!
  ALLOCATE(d_residual_by_d_scale(n_set))
!
  SELECT CASE(i_type_residual)
!
    CASE (ip_scale_trans_residual)
      d_residual_by_d_scale(:)=0.0_RealK
      DO i=1, n_set
        DO j=1, length_set(i)
          derivative_term=0.0_RealK
          trans_fit=0.0_RealK
          DO k=1, n_term
            sub_term=w_esft(k)*EXP(-k_esft(k)*u_set(j, i) * &
              func_scale_90(p_array(i), t_array(i), &
                i_scale_function, ScalePrm))
            trans_fit=trans_fit+sub_term
            derivative_term=derivative_term - &
              k_esft(k)*u_set(j, i)*sub_term
          ENDDO
          d_residual_by_d_scale(i)=d_residual_by_d_scale(i) - &
            2.0_RealK*weight_point(i, j) * &
            (trans_set(j, i)-trans_fit)*derivative_term
        ENDDO
      ENDDO
!
    CASE (ip_scale_full_scl_res)
      d_residual_by_d_scale(:)=0.0_RealK
      DO i=1, n_set
        DO j=1, length_set(i)
          ss=func_scale_90(p_array(i), t_array(i), &
            i_scale_function, ScalePrm)
          d_residual_by_d_scale(i)=d_residual_by_d_scale(i) + &
            2.0_RealK*weight_point(i, j)* &
            (LOG((offset_residual+ss)/(offset_residual+scaling(j, i)))) * &
            (1.0_RealK/(1.0_RealK+ss))
        ENDDO
      ENDDO
    CASE (ip_scale_mean_scl_res)
      DO i=1, n_set
        ss=func_scale_90(p_array(i), t_array(i), &
          i_scale_function, ScalePrm)
        d_residual_by_d_scale(i) = 2.0_RealK * &
          (((offset_residual+ss)/(offset_residual+scaling(1, i))) - &
          1.0_RealK) * &
          (1.0_RealK/(1.0_realk+scaling(1, i))) + &
          2.0_RealK*(ss-scaling(1, i))
      ENDDO
!
  END SELECT
!
!
! Find the gradient for the pressure and temperature.
  gradient(:)=0.0_RealK
  DO j=1, n_set
    gradient(:) = gradient(:) + d_residual_by_d_scale(j) * &
      func_scale_derivative_90(p_array(j), t_array(j), &
                               i_scale_function, ScalePrm)
  ENDDO
!
  DEALLOCATE(d_residual_by_d_scale)
!
!
!
  RETURN
END SUBROUTINE residual_gradient_90
