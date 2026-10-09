! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to generate a CDL-file of surface fields.
!
! Method:
!   A set of coordinates and albedos are read and written in
!   the appropriate format.
!
!- ---------------------------------------------------------------------
      PROGRAM gen_surf_cdl
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE dimensions_spec_ucf
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
     &  , l
!           Loop variable
     &  , k
!           Loop variable
!
!     Field to be created:
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_basis
!           Number of basis functions
      CHARACTER
     &    file_surf*80
!           Name of null file
     &  , name_surf*24
!           Name of the horizontal field
     &  , unit_surf*20
!           Units of the horizontal field
     &  , long_surf*40
!           Long name of the horizontal field
!
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , weight(npd_profile, npd_brdf_basis_fnc)
!           Null field generated
!
      data ierr/i_normal/
!     Subroutines called:
      EXTERNAL
     &    output_surf_cdl
!
!
!
!
!     Horizontal coordinates:
!
      WRITE(iu_stdout, '(a)') 'enter the number of longitudes.'
1     read(iu_stdin, *, iostat=ios) n_longitude
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 1
        ENDIF
      ENDIF
      WRITE(iu_stdout, '(a)') 'enter the number of latitudes.'
2     read(iu_stdin, *, iostat=ios) n_latitude
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 2
        ENDIF
      ENDIF
      n_profile=n_longitude*n_latitude
      IF ( (n_profile > npd_profile).OR.
     &     (n_latitude > npd_latitude).OR.
     &     (n_longitude > npd_longitude) ) THEN
        WRITE(iu_err, '(/a)') '+++ too many profiles:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(/a)') 'please reduce.'
          goto 1
        ENDIF
      ENDIF
!
      WRITE(iu_stdout, '(a)') 'specify the longitudes.'
      READ(iu_stdin, *, iostat=ios) (longitude(k), k=1, n_longitude)
      WRITE(iu_stdout, '(a)') 'specify the latitudes.'
      READ(iu_stdin, *, iostat=ios) (latitude(k), k=1, n_latitude)
!
!     Set up the basis functions
      WRITE(iu_stdout, '(a)') 'enter the number of basis functions.'
3     read(iu_stdin, *, iostat=ios) n_basis
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 3
        ENDIF
      ENDIF
      IF (n_basis > npd_brdf_basis_fnc) THEN
        WRITE(iu_err, '(a)') '+++ too many basis functions:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 3
        ENDIF
      ENDIF
!
!
!     Prepare the correct field.
      DO k=1, n_basis
        WRITE(iu_stdout, '(/a, i5, a)') 
     &    'specify the weights for basis function ' 
     &    , k, ' with longitude changing most rapidly.'
        READ(iu_stdin, *) (weight(l, k), l=1, n_profile)
      ENDDO
!
!     Prepare the field for output.
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_surf
      WRITE(iu_stdout, '(/a)') 'enter the name of the output field.'
      READ(iu_stdin, '(a)') name_surf
      WRITE(iu_stdout, '(/a)') 'enter the units of the output field.'
      READ(iu_stdin, '(a)') unit_surf
      WRITE(iu_stdout, '(/a)') 'enter the long-name 
     &                                       of the output field.'
      READ(iu_stdin, '(a)') long_surf
!
      CALL output_surf_cdl(ierr
     &  , file_surf
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_basis
     &  , name_surf, 'float', unit_surf, long_surf, weight
     &  , npd_profile, npd_latitude, npd_longitude, npd_brdf_basis_fnc
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!
!
      STOP
      END
