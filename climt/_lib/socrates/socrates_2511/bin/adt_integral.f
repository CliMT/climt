! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate scattering properties using ADT.
!
! Method:
!   The standard formulae for the extinction and scattering
!   cross-sections, integrated across a modified-gamma
!   distribution, given by ADT are implemented.
!
!- ---------------------------------------------------------------------
      SUBROUTINE adt_integral(alpha_mg, size, n_relative
     &   , q_ext, q_scatter)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE bna_factor_ccf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      REAL  (RealK), Intent(IN) ::
     &    alpha_mg
!           Parameter of modified gamma distribution
     &  , size
!           Size parameter
      complex  (RealK), Intent(IN) ::
     &    n_relative
!           Relative refractive index
      REAL  (RealK), Intent(OUT) ::
     &    q_ext
!           Extinction efficiency
     &  , q_scatter
!           Scattering efficiency
!
!     Local variables.
      REAL  (RealK) ::
     &    u
!           Dummy variables for adt
     &  , v
!              "
     &  , r
!              "
      complex  (RealK) ::
     &    z
!              "
!
!
      v=size*aimag(n_relative)
      r=4.0_RealK*v
      u=size*(REAL(n_relative, RealK)-1.0_RealK)
      z=2.0_RealK*cmplx(v, u)
!
!     Evaluation of q_ext is split into two parts 
!     to avoid ill-conditioning.
      IF (abs(z) > 1.0e-04_RealK) THEN
         q_ext=4.0_RealK*REAL(0.5_RealK+(1.0_RealK
     &      /(alpha_mg*(alpha_mg+1.0_RealK)))
     &      *(z*alpha_mg+(1.0_RealK+z)
     &      *(1.0_RealK-EXP(alpha_mg*LOG(z+1.0_RealK))))
     &      /(z*z*EXP((1.0_RealK+alpha_mg)*LOG(z+1.0_RealK))))
      ELSE
         q_ext=4.0_RealK*REAL((alpha_mg+2.0_RealK)*z
     &      *(1.0_RealK/3.0_RealK-(alpha_mg+3.0_RealK)
     &      *z*0.125_RealK*(1.0_RealK
     &      -4.0_RealK*z*(alpha_mg+4.0_RealK)
     &      /15.0_RealK)))
      ENDIF
!
!     Rescale R with the BNA factor.
      r=factor_bna*r
      IF (abs(r) > 1.0e-04_RealK) THEN
         q_scatter=q_ext-2.0_RealK*(0.5_RealK+(1.0_RealK
     &      /(alpha_mg*(alpha_mg+1.0_RealK)))
     &      *(r*alpha_mg+(1.0_RealK+r)
     &      *(1.0_RealK-EXP(alpha_mg*LOG(r+1.0_RealK))))
     &      /(r*r*EXP((2.0_RealK+alpha_mg)*LOG(r+1.0_RealK))))
      ELSE
         q_scatter=q_ext-2.0_RealK*(alpha_mg+2.0_RealK)*r
     &      *(1.0_RealK/3.0_RealK-(alpha_mg+3.0_RealK)*r
     &      *0.125_RealK*(1.0_RealK
     &      -4.0_RealK*r*(alpha_mg+4.0_RealK)/15.0_RealK))
      ENDIF
!
!
!
      RETURN
      END
