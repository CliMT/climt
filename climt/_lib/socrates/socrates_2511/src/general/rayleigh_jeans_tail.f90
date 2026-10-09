! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the Rayleigh-Jeans approximation.
!
! Method:
!   The standard Rayleigh-Jeans formula is used.
!
! Note:
!   This routine is used to provide a tail to the Planck function.
!
!- ---------------------------------------------------------------------
FUNCTION rayleigh_jeans_tail(Sol, lambda)

  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: astronomical_unit, c_light, k_boltzmann
  USE def_solarspec, ONLY: StrSolarSpec

  IMPLICIT NONE


  TYPE (StrSolarSpec), INTENT(IN) :: Sol
!   Solar spectrum
  REAL (RealK), INTENT(IN) :: lambda
!   Initial wavelength
  REAL (RealK) :: rayleigh_jeans_tail
!   Tail irradiance.


! Evaluate black body flux.
  rayleigh_jeans_tail=2.0_RealK*c_light*k_boltzmann &
     *Sol%t_effective/(3.0_RealK*lambda**3)
! Scale to get the flux at 1 AU.
  rayleigh_jeans_tail=rayleigh_jeans_tail &
     *(Sol%radius/astronomical_unit)**2

END FUNCTION rayleigh_jeans_tail
