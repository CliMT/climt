! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a SAMSON file.
!
! Method:
!   The input file is opened and the standard reading routine
!   is called. This is a restricted version for conversion of
!   formats. Only fields on pressure levels can be read.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_genln2_flux(ierr
     &  , file_name
     &  , n_profile, n_level, p, field
     &  , nd_profile, nd_layer
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
!
!     Sizes of arrays:
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Size allocated for profiles
     &  , nd_layer
!           Size allocated for layers
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
!
      INTEGER, Intent(OUT) ::
     &    n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
      REAL  (RealK), Intent(OUT) ::
     &    p(nd_layer+1)
!           Pressures of the input field
     &  , field(nd_profile, nd_layer+1)
!           Field read in
!
!     Local variables
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
     &  , iunit
!           Unit number for reading 
      CHARACTER
     &    file_title*256
!           Title of field
      REAL  (RealK) ::
     &    zero
!           Dummy variable
!
!     Subroutines called:
      EXTERNAL
     &    get_free_unit
!
!
!
      CALL get_free_unit(ierr, iunit)
      IF (ierr /= i_normal) RETURN
      OPEN(unit=iunit, file=file_name)
!
      READ(iunit, '(a)') file_title
      READ(iunit, '(a)') file_title
      READ(iunit, '(a)') file_title
      n_profile=6
      READ(iunit, *)  file_title, file_title, n_level
      DO i=1, n_level
        READ(iunit, *)  p(i), (field(l, i), l=1, n_profile)
      ENDDO
!
      CLOSE(iunit)
!
!
!
      RETURN
      END
