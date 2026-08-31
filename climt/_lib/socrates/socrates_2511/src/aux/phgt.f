! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to find the height of a pressure level.
!
! Method:
!   A CDL-file of heights is read and the user supplies
!   a pressure level. The height of this level is interpolated.
!
!- ---------------------------------------------------------------------
      PROGRAM phgt
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
!
!     Declaration of variables.
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
!
      CHARACTER
     &    file*80
!           Name of input file
     &  , name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
!
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of levels
     &  , n_level
!           Number of levels
      INTEGER
     &    i_mode
!           Mode of interpolation
     &  , l
!           Loop variable
     &  , i
!           Loop variable
      LOGICAL
     &    l_splined
!           True when spline is done
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , p(npd_profile, npd_layer+1)
!           Pressures
     &  , field_hgt(npd_profile, npd_layer+1)
!           Height field
     &  , height
!           Height for evaluation
     &  , p_evaluated(npd_profile)
!           Evaluated pressure
     &  , height_field(npd_layer+1)
!           Reversed height field
     &  , p_reversed(npd_layer+1)
!           Reversed pressure field
     &  , x(npd_layer+1)
!           Abscissa for interpolation
     &  , y(npd_layer+1)
!           Ordinate for interpolation
     &  , y2(npd_layer+1)
!             2nd derivative of ordinate
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl
!
!
      data 
     &    l_vert_coord/.false./
     &  , l_splined/.false./
     &  , i_mode/IP_1_lin_log/
      data ierr/i_normal/
!
!
!
!     Read the file of heights.
      WRITE(iu_stdout, '(a)')
     &  'enter the name of the file containing the heights.'
      READ(iu_stdin, '(a)') file
      CALL assign_input_vert_cdl(ierr
     &  , file, 'heights', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, field_hgt
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
!
!     Read in the height at which the pressure is to be found.
      WRITE(iu_stdout, '(/a)')
     &  'enter the height at which pressures are to be found.'
1     read(iu_stdin, *, iostat=ios) height
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ unrecognized response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_err, '(a)') 'please re-enter.'
          goto 1
        ENDIF
      ENDIF
!
!
!     Reverse the height and pressure arrays so that we fit on an array
!     of ascending values.
      DO l=1, n_profile
!
        DO i=1, n_level
          height_field(n_level+1-i)=field_hgt(l, i)
          p_reversed(n_level+1-i)=p(l, i)
        ENDDO
        l_splined=.false.
        CALL interpolate_p(n_level, height_field
     &    , p_reversed, x, y, y2
     &    , height, p_evaluated(l)
     &    , i_mode, l_splined
     &    )
!
      ENDDO
!
!     Print the final field.
      WRITE(iu_stdout, '(1x, 1pe14.7)') 
     &  (p_evaluated(l), l=1, n_profile)
!
!
!
      END
