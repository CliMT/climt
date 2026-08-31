! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find saturated mixing ratio of water.
!
! Method:
!   The standard Goff-Gratsch formula is implemented.
!
!- ---------------------------------------------------------------------
      SUBROUTINE qsat_gg(qs, t, p)
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
     &  , u4
!           Temporary variables
     &  , ew
!           Saturation vapour pressure of water vapour
     &  , r
!           Humidity mixing ratio
!
!
!
      x=3.7316e+02_RealK/t
      u1=-7.90298_RealK*(x-1.0_RealK)
      u2=5.02808_RealK*log10(x)
      u3=-1.3816e-07_RealK*(1.0e+01_RealK**(1.1344e+01_RealK
     &  *(1.0_RealK-1.0_RealK/x))-1.0_RealK)
      u4=8.1328e-03_RealK*(1.0e+01_RealK**(-3.49149_RealK
     &  *(x-1.0_RealK))-1.0_RealK)
      ew=1.013246e+05_RealK*1.0e+01_RealK**(u1+u2+u3+u4)
      r=6.2197e-01_RealK*ew/(p-ew)
      qs=r/(1.0_RealK+r)
!
!
!
      RETURN
      END
