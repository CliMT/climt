! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the derivative Planckian function.
!
SUBROUTINE d_planck_90(nu, t, db)
!
!
! Description:
!   This routine returns an array of values of the derivative of
!   the Planckian radiance with respect to temperature calculated 
!   at a range of supplied frequencies.
!
! Method:
!  Straightforward.
!
!
!
! Modules used
  USE realtype_rd
  USE rad_ccf, ONLY: pi, h_planck, c_light, k_boltzmann
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  REAL  (RealK), Intent(IN) :: nu(:)
!   Frequencies
  REAL  (RealK), Intent(IN) :: T
!   Temperature
  REAL  (RealK), Intent(OUT) :: db(:)
!   Derivative of the Planckian function
!
!
!
! Evaluate a negative exponential to ensure conditioning.
  db = EXP(-( h_planck * c_light * nu / ( k_boltzmann * T) ) )
  db = 2.0_RealK * h_planck * c_light**2 * db * (nu**3) / &
      (1.0_RealK - db)**2
!
!
!
END SUBROUTINE d_planck_90
