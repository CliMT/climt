! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the scaling function for gaseous absorption.
!
FUNCTION func_scale_90 &
!
(p_fnc, t_fnc, i_scale_function, ScalePrm) &
!
RESULT(f)
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
  REAL  (RealK) :: f
!   Value of scaling function
!
!
!
  SELECT CASE(i_scale_function)
    CASE (ip_scale_power_law)
      f=exp(ScalePrm(1)*p_fnc+ScalePrm(2)*t_fnc)
    CASE (ip_scale_power_quad)
      f=exp(ScalePrm(1)*p_fnc) &
        *(1.0e+00_RealK+ScalePrm(2)*t_fnc &
        +ScalePrm(3)*t_fnc*t_fnc)
    CASE (ip_scale_doppler_quad)
      f=exp(ScalePrm(1) &
        *log((p_fnc+exp(ScalePrm(2))) &
        /(1.0e+00_RealK+exp(ScalePrm(2))))) &
        *(1.0e+00_RealK+ScalePrm(3)*t_fnc &
        +ScalePrm(4)*t_fnc*t_fnc)
  END SELECT
!
!
!
  RETURN
END
