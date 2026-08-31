! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates Rayleigh scattering coefficient at STP.
!
!- ---------------------------------------------------------------------
REAL(RealK) FUNCTION rayleigh_scatter(i_gas, wavelength)

  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi, r_gas, n_avogadro
  USE gas_list_pcf, ONLY: molar_weight, depolarization_factor, rayleigh_cutoff
  USE refract_re_ccf, ONLY: refract_re_m1, p_stp, t_stp

  IMPLICIT NONE

! Input variables
  INTEGER, INTENT(IN) :: i_gas
!   Type index of gas
  REAL(RealK), INTENT(IN) :: wavelength
!   Wavelength

! Local variables
  REAL(RealK) :: rho_stp
!   Mass density at STP
  REAL(RealK) :: wavelength_m2
!   Reciprocal of wavelength squared

! Calculate density at STP for current gas
  rho_stp = p_stp*molar_weight(i_gas)*1.0E-03_RealK/(r_gas*t_stp)

! Use the standard expression for the Rayleigh scattering
! coefficient, but include an extra density factor to give it
! in units of mass.
  wavelength_m2=1.0_RealK/(wavelength**2)
  rayleigh_scatter &
      =(8.0_RealK*pi**3/3.0_RealK) &
      *(((refract_re_m1(i_gas, wavelength)+2.0_RealK) &
      *refract_re_m1(i_gas, wavelength)*wavelength_m2)**2/n_avogadro) &
      *((6.0_RealK+3.0_RealK*depolarization_factor(i_gas)) &
      /(6.0_RealK-7.0_RealK*depolarization_factor(i_gas))) &
      *(molar_weight(i_gas)*1.0E-03_RealK/rho_stp**2)

! Zero Rayleigh scattering below cutoff wavelength (where absorption dominates)
  IF (wavelength < rayleigh_cutoff(i_gas)+0.1e-12) rayleigh_scatter=0.0_RealK

END FUNCTION rayleigh_scatter
