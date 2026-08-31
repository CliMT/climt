! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the planck function.
!
! Method:
!  Straightforward.
!
!- ---------------------------------------------------------------------
      FUNCTION planck(t, lambda)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE rad_ccf

      IMPLICIT NONE


!     Dummy argyuments.
      REAL  (RealK) ::
     &    planck
!           Name of function
      REAL  (RealK), Intent(IN) ::
     &    t
!           Temperature
     &  , lambda
!           Wavelength
!
!     Local variables.
      REAL  (RealK) ::
     &    exponential
!           Inverse exponential factor
!
!
!
!     Evaluate a negative exponential to ensure conditioning.
      exponential=exp(-(h_planck*c_light/k_boltzmann)/(lambda*t))
      planck=2.0_RealK*h_planck*c_light**2*exponential
     &   /((lambda**5)*(1.0_RealK-exponential))
!
!
!
      RETURN
      END
