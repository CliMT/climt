! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the derivative of the scaling function.
!
FUNCTION func_scale_derivative_90 &
!
(p_fnc, t_fnc, i_scale_function, ScalePrm) &
!
RESULT(df)
!
! Description:
!   Straightforward.
!
!
! Modules used:
  USE realtype_rd
  USE rad_pcf
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments
  INTEGER, Intent(IN) :: i_scale_function
!   Type of scaling function
  REAL  (RealK), Intent(IN) :: p_fnc
!   Function of pressure
  REAL  (RealK), Intent(IN) :: t_fnc
!   Function of temperature
  REAL  (RealK), Intent(IN), Dimension(:) :: ScalePrm
!   Scaling parameters
!
  REAL  (RealK), Pointer, Dimension(:) :: df
!   Derivative of the scaling function
!
! Local variables:
  REAL  (RealK) :: f
!   Temporary variable (scaling function)
  REAL  (RealK) :: v1
!   Temporary variable
  REAL  (RealK) :: v2
!   Temporary variable
  REAL  (RealK) :: v3
!   Temporary variable
!
!
!
  ALLOCATE(df(SIZE(ScalePrm)))
  SELECT CASE(i_scale_function)
    CASE (IP_scale_power_law)
      f = exp(ScalePrm(1)*p_fnc+ScalePrm(2)*t_fnc)
      df(1:2) = (/  p_fnc * f, t_fnc * f /)
    CASE (IP_scale_power_quad)
      df(1:3) = (/ p_fnc * EXP(ScalePrm(1)*p_fnc) * &
                           (1.0_RealK+ScalePrm(2)*t_fnc + &
                           ScalePrm(3)*t_fnc*t_fnc), &
                   t_fnc * EXP(ScalePrm(1)*p_fnc), &
                   t_fnc * t_fnc * EXP(ScalePrm(1)*p_fnc) &
                /)
    CASE (IP_scale_doppler_quad) 
      v1=exp(ScalePrm(2))
      v2=log((p_fnc+v1)/(1.0_RealK+v1))
      v3=(1.0_RealK+ScalePrm(3)*t_fnc+ScalePrm(4)*t_fnc*t_fnc)
      df(1:4) = (/ v2 * EXP(ScalePrm(1)*v2) * v3, &
                   ScalePrm(1) * ( (1.0_RealK-p_fnc) / &
                     (1.0_RealK+v1)**2 ) * &
                     v1 * EXP((ScalePrm(1)-1.0_RealK)*v2) * v3, &
                   t_fnc * EXP(ScalePrm(1)*v2), &
                   t_fnc * t_fnc * EXP(ScalePrm(1)*v2) &
                /)
  END SELECT
!
!
!
  RETURN
END
