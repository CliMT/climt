! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to evaluate the value of a CDL-field at a given pressure.
!
! Method:
!   The field is read in and the value at the selected pressure
!   is calculated by interpolation according to the selected
!   scheme.
!
!- ---------------------------------------------------------------------
      PROGRAM fval
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE error_pcf
      USE interpolate_p_mod
!
!
      IMPLICIT NONE
!
!
!     Declaration of variables.
      CHARACTER
     &    file*80
!           Name of input file
     &  , name_vert_coord*24
!           Name of the vetical coordinate
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
     &  , n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of input profiles
     &  , n_level
!           Number of input levels
     &  , i_mode
!           Mode of interpolation
     &  , i
!           Loop variable
     &  , l
!           Loop variable
      LOGICAL
     &    l_splined
!           True when spline is done
     &  , l_vert_coord
!           Flag asserting that vertical coordinate is set
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , field(npd_profile, npd_layer+1)
!           Input field
     &  , p(npd_profile, npd_layer+1)
!           Pressures in input field
     &  , p_evaluate
!           Pressure for evaluation
     &  , field_evaluated(npd_profile)
!           Evaluated field
     &  , a(npd_layer+1)
!           Gathered interpolation field
     &  , x(npd_layer+1)
!           Abscissa for interpolation
     &  , y(npd_layer+1)
!           Ordinate for interpolation
     &  , y2(npd_layer+1)
!             2nd derivative of ordinate
     &  , p_int(npd_layer+1)
!           Pressure for interpolation
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl
!
!
      data 
     &    l_splined/.false./
     &  , l_vert_coord/.false./
      data ierr/i_normal/
      data n_latitude/0/
      data n_longitude/0/
!
!
!
!     Read in the basic field.
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the field to be evaluated.'
      READ(iu_stdin, '(a)') file
      CALL assign_input_vert_cdl(ierr
     &  , file, 'input field', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, field
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!     Read in the pressure at which the evaluation is to be done.
      WRITE(iu_stdout, '(/a)')
     &  'enter the pressure at which the field is to be evaluated.'
1     read(iu_stdin, *, iostat=ios) p_evaluate
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ unrecognized response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 1
        ENDIF
      ENDIF
!
      WRITE(iu_stdout, '(/a)') 'enter the mode of interpolation.'
2     read(iu_stdin, *, iostat=ios) i_mode
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ unrecognized response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 2
        ENDIF
      ENDIF
      IF ( (i_mode /= IP_1_lin_lin).AND.
     &     (i_mode /= IP_1_log_lin).AND.
     &     (i_mode /= IP_1_lin_log).AND.
     &     (i_mode /= IP_1_log_log).AND.
     6     (i_mode /= IP_3_lin_lin).AND.
     &     (i_mode /= IP_3_log_lin).AND.
     &     (i_mode /= IP_3_lin_log).AND.
     &     (i_mode /= IP_3_log_log) ) THEN
        WRITE(iu_err, '(a)') '+++ unsupported mode:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 2
        ENDIF
      ENDIF
!
!     Perform the interpolation.
      DO l=1, n_profile
        DO i=1, n_level
          a(i)=field(l, i)
          p_int(i)=p(l, i)
        ENDDO
        l_splined=.false.
        CALL interpolate_p(n_level, p_int
     &    , a, x, y, y2
     &    , p_evaluate, field_evaluated(l)
     &    , i_mode, l_splined
     &    )
      ENDDO
!
!     Print the final field.
      WRITE(iu_stdout, '(4(4x, 1pe14.7))') (field_evaluated(l)
     &  , l=1, n_profile)
!
!
!
      END
