! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate the extinction of an SRA profile.
!
! Method:
!   A representation is given and the extinction is calculated
!   appropriately.
!
!- ---------------------------------------------------------------------
      FUNCTION extinction_profile(z
     &  , i_representation, n_parm, parm
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE aerosol_representation_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      INTEGER, Intent(IN) ::
     &    i_representation
!           Type of height variation
     &  , n_parm
!           Size of parameter array
      REAL  (RealK), Intent(IN) ::
     &    z
!           Height of level
     &  , parm(n_parm)
!           Functional parameters
      REAL  (RealK) ::
     &    extinction_profile
!           Evaluated extinction
!
!     Local variables.
      REAL  (RealK) ::
     &    lambda
!           Fraction of linear interval
!
!
      IF (i_representation == i_constant) THEN
        extinction_profile=parm(1)
      ELSE IF (i_representation == i_linear) THEN
        lambda=(z-parm(3))/(parm(4)-parm(3))
        extinction_profile=(parm(1)-parm(2))*lambda+parm(2)
      ELSE IF (i_representation == i_exponential) THEN
        extinction_profile=parm(1)*exp(-z/parm(2))
      ENDIF
!
!
!
      RETURN
      END
