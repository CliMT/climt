! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to generate a CDL-file of a 2-D field.
!
! Method:
!   The horizontal coordinates and values of the field
!   are read and passed to teh output routines.
!
!- ---------------------------------------------------------------------
      PROGRAM gen_horiz_cdl
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
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
     &  , k
!           Loop variable
     &  , l
!           Loop variable
!
!     Field to be created:
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of pressure levels
      CHARACTER
     &    file_horiz*80
!           Name of null file
     &  , name_horiz*24
!           Name of the horizontal field
     &  , unit_horiz*20
!           Units of the horizontal field
     &  , long_horiz*40
!           Long name of the horizontal field
!
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , field_horiz(npd_profile)
!           Null field generated
      data ierr/i_normal/
!
!     Subroutines called:
      EXTERNAL
     &    output_horiz_cdl
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
!
!     Prepare the correct field.
      WRITE(iu_stdout, '(/a)') 'enter the values of the output '
     &  //'field with longitude changing most rapidly.'
      READ(iu_stdin, *) (field_horiz(l), l=1, n_profile)
!
!     Obtain the properties of the output field.
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_horiz
      WRITE(iu_stdout, '(/a)') 'enter the name of the output field.'
      READ(iu_stdin, '(a)') name_horiz
      WRITE(iu_stdout, '(/a)') 'enter the units of the output field.'
      READ(iu_stdin, '(a)') unit_horiz
      WRITE(iu_stdout, '(/a)') 'enter the long-name
     &                                     of the output field.'
      READ(iu_stdin, '(a)') long_horiz
!
      CALL output_horiz_cdl(ierr
     &  , file_horiz
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , name_horiz, 'float', unit_horiz, long_horiz, field_horiz
     &  , npd_profile, npd_latitude, npd_longitude
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!
!
      STOP
      END
