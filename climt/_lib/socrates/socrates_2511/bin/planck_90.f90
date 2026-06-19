! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the Planckian function.
!
SUBROUTINE planck_90(nu, t, b)
!
!
! Description:
!   This routine returns an array of values of the Planckian
!   radiance calculated at a range of supplied frquencies.
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
  REAL  (RealK), Intent(OUT) :: b(:)
!   Calculated Planckian radiances
!
!
!
! Evaluate a negative exponential to ensure conditioning.
  b = EXP(-( h_planck * c_light * nu / ( k_boltzmann * T) ) )
  b = 2.0_RealK * h_planck * c_light**2 * b * (nu**3) / &
      (1.0_RealK - b)
!
!
!
END SUBROUTINE planck_90
