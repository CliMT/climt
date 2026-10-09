! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Function to calculate the Rayleigh scattering coefficient at S.T.P.
!
!-----------------------------------------------------------------------
      FUNCTION rayleigh_scatter_air(lambda)

      USE realtype_rd
      USE rad_ccf

      IMPLICIT NONE


!     Dummy variables.
      REAL  (RealK) ::
     &    rayleigh_scatter_air
!           Name of function
      REAL  (RealK), Intent(IN) ::
     &    lambda
!           Wavelength

!     Local variables.
      REAL  (RealK) ::
     &    refract_index_m1
!           Refractive index at 0C less 1.
     &  , lambda_m2
!           Reciprocal of wavelength squared
     &  , lambda_m2_limit
!           Reciprocal of wavelength squared limited to lambda > 175nm


!     Calculate the refractive index using Edlen's formula converted to
!     apply at 273.16K. (Limit at short wavelengths where error in fit
!     becomes large.)
      lambda_m2_limit=1.0_RealK/(MAX(lambda,0.175e-6_RealK)**2)
      refract_index_m1=6.78606e-05_RealK+3.11180e+10_RealK
     &   /(1.46e+14_RealK-lambda_m2_limit)
     &   +2.69425e+08_RealK/(4.1e+13_RealK-lambda_m2_limit)

!     Use the standard expression for the Rayleigh scattering
!     coefficient, but include an extra density factor to give it
!     in units of mass.
      lambda_m2=1.0_RealK/(lambda*lambda)
      rayleigh_scatter_air
     &   =(8.0_RealK*pi**3/3.0_RealK)
     &   *(((refract_index_m1+2.0_RealK)
     &   *refract_index_m1*lambda_m2)**2/n_avogadro)
     &   *((6.0_RealK+3.0_RealK*rho_n)
     &   /(6.0_RealK-7.0_RealK*rho_n))
     &   *(mol_weight_air/rho_air_stp**2)

!     Zero Rayleigh scattering below 175nm
!     (Rayleigh scattering is interpolated from the band edges so a small
!     tolerance is used for correct treatment of bands with a limit at 175nm)
      IF (lambda < 0.1750001e-6) rayleigh_scatter_air=0.0_RealK

      END
