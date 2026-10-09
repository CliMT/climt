! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to find whether an aerosol component is in a spectral file.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      FUNCTION l_find_component(n_aerosol, type_aerosol, i_component)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_spec_ucf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      INTEGER, Intent(IN) ::
     &    n_aerosol
!           Number of aerosols
     &  , type_aerosol(npd_aerosol_species)
!           Types of aerosols
     &  , i_component
!           Aerosol component
      LOGICAL ::
     &    l_find_component
!           Function to find component
!
!     Local variables.
      INTEGER
     &    i
!           Loop variable
!
!
      l_find_component=.false.
      i=0
1     i=i+1
      IF (type_aerosol(i) == i_component) THEN
        l_find_component=.true.
        RETURN
      ELSE IF (i < n_aerosol) THEN
        goto 1
      ELSE
        RETURN
      ENDIF
!
!
!
      END
