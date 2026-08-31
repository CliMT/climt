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
!   are read and passed to the output routines.
!
!- ---------------------------------------------------------------------
      PROGRAM gen_view_cdl
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
     &  , i
!           Loop variable
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
     &  , n_viewing_level
!           Number of viewing levels
     &  , n_direction
!           Number of viewing directions
      CHARACTER
     &    file_view*80
!           Name of file of vieing information
!
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , viewing_level(npd_layer+1)
!           Viewing levels
     &  , direction(npd_profile, npd_direction, 2)
!           Viewing directions
!
      data ierr/i_normal/
!     Subroutines called:
      EXTERNAL
     &    output_view_cdl
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
!     Vertical coordinates.
      WRITE(iu_stdout, '(a)') 'enter the number of viewing levels.'
3     read(iu_stdin, *, iostat=ios) n_viewing_level
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 3
        ENDIF
      ENDIF
      IF (n_viewing_level > npd_layer) THEN
        WRITE(iu_err, '(/a)') '+++ too many viewing levels:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(/a)') 'please reduce.'
          goto 3
        ENDIF
      ENDIF
      WRITE(iu_stdout, '(a)') 'specify the viewing levels.'
      READ(iu_stdin, *, iostat=ios) 
     &  (viewing_level(k), k=1, n_viewing_level)
!
      WRITE(iu_stdout, '(a)') 'specify the number 
     &                                   of viewing directions.'
      READ(iu_stdin, *, iostat=ios) n_direction
      WRITE(iu_stdout, '(a, /a)') 
     &  'enter the polar viewing angles with longitude changing most '
     &  , 'rapidly, THEN latitude, then direction.'
      READ(iu_stdin, *, iostat=ios) 
     &  ( ( (direction(l+n_longitude*(i-1), k, 1), l=1, n_longitude)
     &  , i=1, n_latitude), k=1, n_direction)
      WRITE(iu_stdout, '(a, /a)') 
     &  'enter the azimuthal viewing angles with longitude changing '
     &  , 'most rapidly, THEN latitude, then direction.'
      READ(iu_stdin, *, iostat=ios) 
     &  ( ( (direction(l+n_longitude*(i-1), k, 2), l=1, n_longitude)
     &  , i=1, n_latitude), k=1, n_direction)
!
!
!     Obtain the properties of the output field.
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_view
!
      CALL output_view_cdl(ierr
     &  , file_view
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_viewing_level, viewing_level
     &  , n_direction, direction
     &  , npd_profile, npd_latitude, npd_longitude
     &  , (npd_layer+1), npd_direction
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!
!
      STOP
      END
