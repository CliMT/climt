! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to generate a null CDL-file.
!
! Method:
!   Pressure levels are read in and a CDL-file with
!   a null field at these levels is written.
!
!- ---------------------------------------------------------------------
      PROGRAM gen_null
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE def_std_io_icf
      USE gas_list_pcf
      USE rad_pcf
      USE def_data_in_icf
      USE input_head_pcf
      USE shell_sort_mod, ONLY: shell_sort
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
     &  , i
!           Loop variable
     &  , j
!           Loop variable
     &  , k
!           Loop variable
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of pressure levels
     &  , n_level
!           Number of pressure levels
     &  , pointer(npd_layer+1)
!           Pointer for sorting
      CHARACTER(LEN=len_col_header) ::
     &    unit
!           Unit for pressure
      CHARACTER(LEN=80) ::
     &    file_null
!           Name of null file
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      REAL  (RealK) ::
     &    p_level(npd_layer+1)
!           Specified levels
     &  , latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , key(npd_layer+1)
!           Sorting key
     &  , factor
!           Input conversion factor
     &  , field_null(npd_profile, npd_layer+1)
!           Null field generated
!
      data ierr/i_normal/
!     Subroutines called:
      EXTERNAL  output_vert_cdl
!
!
!
!
!
!
!     Vertical coordinates:
!
!     Read from standard input the required pressure levels
      WRITE(iu_stdout, '(/a)')
     &  'enter the unit for the following pressure levels:'
1     read(iu_stdin, '(a)') unit
      k=0
2     k=k+1
      IF (unit == name_unit(k)) THEN
        factor=factor_unit(k)
      ELSE IF (k < npd_unit) THEN
        goto 2
      ELSE
        WRITE(iu_err, '(a)') '+++ unrecognized unit:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter'
          goto 1
        ENDIF
      ENDIF
      WRITE(iu_stdout, '(a)') 'enter the number of pressure levels.'
3     read(iu_stdin, *, iostat=ios) n_level
      IF  (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 3
        ENDIF
      ENDIF
      IF (n_level > npd_layer+1) THEN
        WRITE(iu_err, '(/a)') 'too many levels:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please reduce.'
          goto 3
        ENDIF
      ENDIF
      WRITE(iu_stdout, '(a)') 'specify pressures at these levels.'
      READ(iu_stdin, *, iostat=ios) (p_level(k), k=1, n_level)
!
!
!     Horizontal coordinates:
!
      WRITE(iu_stdout, '(a)') 'enter the number of longitudes.'
4     read(iu_stdin, *, iostat=ios) n_longitude
      IF  (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 4
        ENDIF
      ENDIF
      WRITE(iu_stdout, '(a)') 'enter the number of latitudes.'
5     read(iu_stdin, *, iostat=ios) n_latitude
      IF  (ios /= 0) THEN
        WRITE(iu_err, '(a)') 'illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 5
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
          goto 4
        ENDIF
      ENDIF
      WRITE(iu_stdout, '(a)') 'specify the longitudes.'
      READ(iu_stdin, *, iostat=ios) (longitude(k), k=1, n_longitude)
      WRITE(iu_stdout, '(a)') 'specify the latitudes.'
      READ(iu_stdin, *, iostat=ios) (latitude(k), k=1, n_latitude)
!
!
!
!     Obtain the name of the output file.
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_null
!
!     Sort the pressures
      DO i=1, n_level
        key(i)=p_level(i)
        pointer(i)=i
      ENDDO
      CALL shell_sort(n_level, pointer, key)
      DO i=1, n_level
        key(i)=p_level(pointer(i))
      ENDDO
      DO i=1, n_level
        p_level(i)=factor*key(i)
      ENDDO
!
!     Prepare the correct field.
      DO i=1, n_level
        DO j=1, n_profile
          field_null(j, i)=0.0_RealK
        ENDDO
      ENDDO
!
      CALL output_vert_cdl(ierr
     &  , file_null
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level, 'plev', 4, p_level
     &  , 'null', 'float', 'none', 'null field', field_null
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!
!
      STOP
      END
