! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Function to calculate the Rayleigh scattering coefficient at S.T.P.
!
!-----------------------------------------------------------------------
FUNCTION rayleigh_scatter_h2he(lambda, wavelength_refract_index_H2,            &
  refract_index_H2, n_points)

  USE realtype_rd,  ONLY: RealK
  USE rad_ccf,      ONLY: pi, n_avogadro, A_H, A_He, mol_weight_h2he,          &
                          rho_h2he_stp, rho_n_h2he
  USE interp1d_mod, ONLY: interp1d

  IMPLICIT NONE


! Dummy variables.
  INTEGER, Intent(IN) :: &
      n_points

  REAL  (RealK) :: &
      rayleigh_scatter_h2he
!           Name of function
  REAL  (RealK), Intent(IN) :: &
      lambda
!           Wavelength
  REAL (RealK), Intent(IN) :: &
      wavelength_refract_index_H2(1:n_points)
!           Wavelength at which H2 refractive index data is evaluated
  REAL (RealK), Intent(IN) :: &
      refract_index_H2(1:n_points)
!           Refractive index of H2

! Local variables.
  REAL  (RealK) :: &
      refract_index_m1
!           Refractive index at 0C less 1.
  REAL  (RealK) :: &
      lambda_m2
!           Reciprocal of wavelength squared
  REAL  (RealK) :: &
      refract_index_H2_lambda
!           H2 refractive index evaluated at lambda
  REAL  (RealK) :: &
      refract_index_He_lambda
!           He refractive index evaluated at lambda
  REAL  (RealK) :: &
      temp
!           Temporary storage of RHS in Lorentz-Lorentz equation


! Refractive index at 0C for H2 and He
  lambda_m2=1.0_RealK/(lambda*lambda)
  refract_index_H2_lambda = interp1d(wavelength_refract_index_H2, &
    refract_index_H2, lambda, n_points)
  refract_index_He_lambda = 1.0_RealK + &
    0.01470091_RealK/(423.98_RealK-(lambda*1e+6_RealK)**(-2.0_RealK))

! Calculate refractive index of mixture using the Lorentz-Lorentz equation
  temp = (A_H/2.0_RealK)/(A_H/2.0_RealK + A_He)* &
    (refract_index_H2_lambda**2.0_RealK - 1.0_RealK)/ &
    (refract_index_H2_lambda**2.0_RealK + 2.0_RealK) + &
    (A_He)/(A_H/2.0_RealK + A_He)* &
    (refract_index_He_lambda**2.0_RealK - 1.0_RealK)/ &
    (refract_index_He_lambda**2.0_RealK + 2.0_RealK)
  refract_index_m1 = &
    sqrt((2.0_RealK*temp + 1.0_RealK)/(1_RealK - temp)) - 1.0_RealK

! Alternative, simplified way of calculating the refractive index
!  refract_index_m1 = &
!      refract_index_H2_lambda*(A_H/2.0_RealK)/(A_H/2.0_RealK + A_He) + &
!      refract_index_He_lambda*(A_He)/(A_H/2.0_RealK + A_He) - 1.0_RealK
!
!     Use the standard expression for the Rayleigh scattering
!     coefficient, but include an extra density factor to give it
!     in units of mass.

  rayleigh_scatter_h2he &
     =(8.0_RealK*pi**3/3.0_RealK) &
     *(((refract_index_m1+2.0_RealK) &
     *refract_index_m1*lambda_m2)**2/n_avogadro) &
     *((6.0_RealK+3.0_RealK*rho_n_h2he) &
     /(6.0_RealK-7.0_RealK*rho_n_h2he)) &
     *(mol_weight_h2he/rho_h2he_stp**2)

END FUNCTION rayleigh_scatter_h2he
