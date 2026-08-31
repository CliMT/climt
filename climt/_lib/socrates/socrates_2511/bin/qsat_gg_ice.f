! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the saturated mixing ratio for ice.
!
! Method:
!   The standard Goff-Gratsch formula is implemented.
!
!- ---------------------------------------------------------------------
      SUBROUTINE qsat_gg_ice(qs, t, p)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments
      REAL  (RealK), Intent(IN) ::
     &    t
!           Temperature
     &  , p
!           Pressure
      REAL  (RealK), Intent(OUT) ::
     &    qs
!           Saturation mixing ratio
!
!     Local variables.
      REAL  (RealK) ::
     &    x
     &  , u1
     &  , u2
     &  , u3
!           Temporary variables
     &  , ew
!           Saturation vapour pressure of water vapour
     &  , r
!           Humidity mixing ratio
!
!
!
      x=2.7316e+02_RealK/t
      u1=-9.09718_RealK*(x-1.0_RealK)
      u2=-3.56654_RealK*log10(x)
      u3=8.76793_RealK*(1.0_RealK-1.0_RealK/x)
      ew=6.1071e+02_RealK*1.0e+01_RealK**(u1+u2+u3)
      r=6.2197e-01_RealK*ew/(p-ew)
      qs=r/(1.0_RealK+r)
!
!
!
      RETURN
      END
