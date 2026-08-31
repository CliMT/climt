! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to interpolate a netCDF-field on to a new grid.
!
! Method:
!   Two netCDF files are read, one specifying the data and the
!   other the new grid. A method of interpolation is also 
!   supplied and the program used this method to interpolate
!   the data to the final grid.
!   
!
!- ---------------------------------------------------------------------
      PROGRAM interp_cdf
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE error_pcf
      USE filenamelength_mod, ONLY: filenamelength
      USE interpolate_p_mod
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
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
!
!     Input field:
      CHARACTER (LEN=filenamelength) :: file_in
!           Name of input file
      INTEGER
     &    n_latitude_in
!           Number of latitudes in the input file
     &  , n_longitude_in
!           Number of longitudes in the input file
     &  , n_profile_in
!           Number of profiles in the input file
     &  , n_level_in
!           Number of levels in the input file
      REAL  (RealK) ::
     &    latitude_in(npd_latitude)
!           Latitude in the input file
     &  , longitude_in(npd_longitude)
!           Longitude in the input file
     &  , p_in(npd_profile, npd_layer+1)
!           Pressure levels in the input file
     &  , field_in(npd_profile, npd_layer+1)
!           Input field
!
!
!     Grid field:
      CHARACTER (LEN=filenamelength) :: file_grid
!           Name of grid file
      INTEGER
     &    n_latitude_grid
!           Number of latitudes in the grid file
     &  , n_longitude_grid
!           Number of longitudes in the grid file
     &  , n_profile_grid
!           Number of profiles in the grid file
     &  , n_level_grid
!           Number of levels in the grid file
      REAL  (RealK) ::
     &    latitude_grid(npd_latitude)
!           Latitude in the grid file
     &  , longitude_grid(npd_longitude)
!           Longitude in the grid file
     &  , p_grid(npd_profile, npd_layer+1)
!           Pressure levels in the grid file
     &  , field_grid(npd_profile, npd_layer+1)
!           Grid field
      REAL  (RealK) ::
     &    p_int(npd_layer+1)
!           Pressure levels for interpolation
!
!
!     Output field:
      CHARACTER (LEN=filenamelength) :: file_out
!           Name of output file
      CHARACTER
     &    name_out*24
!           Name of the output field
     &  , unit_out*20
!           Units of the output field
     &  , long_out*40
!           Long name of the output field
      REAL  (RealK) ::
     &    field_out(npd_profile, npd_layer+1)
!           Input field
!
!
      CHARACTER
     &    name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
!
      INTEGER
     &    i_mode
!           Mode of interpolation
     &  , i
!           Loop variable
     &  , l
!           Loop variable
      LOGICAL
     &    l_splined
!           True when splining is done
      REAL  (RealK) ::
     &    a_in(npd_layer+1)
!           Gathered interpolation field
     &  , x_in(npd_layer+1)
!           Abscissa for interpolation
     &  , y_in(npd_layer+1)
!           Ordinate for interpolation
     &  , y2_in(npd_layer+1)
!             2nd derivative of ordinate
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdf
     &  , output_vert_cdf
!
      data 
     &    l_splined/.false./
     &  , l_vert_coord/.false./
      data ierr/i_normal/
      data n_latitude_in/0/
      data n_longitude_in/0/
      data n_latitude_grid/0/
      data n_longitude_grid/0/
!
!
!
!     Read in the basic field.
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the field to be interpolated.'
      READ(iu_stdin, '(a)') file_in
      CALL assign_input_vert_cdf(ierr
     &  , file_in, 'input field', l_vert_coord, name_vert_coord
     &  , .true., n_level_in, .NOT.l_vert_coord
     &  , n_latitude_in, latitude_in, n_longitude_in, longitude_in
     &  , 1
     &  , n_profile_in, n_level_in
     &  , p_in, field_in
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
!     Read in the grid onto which interpolation is done.
      WRITE(iu_stdout, '(/a)')
     &  'enter the name of the file containing the grid.'
      READ(iu_stdin, '(a)') file_grid
      CALL assign_input_vert_cdf(ierr
     &  , file_grid, 'grid field', l_vert_coord, name_vert_coord
     &  , .false., n_level_in, .true.
     &  , n_latitude_grid, latitude_grid
     &  , n_longitude_grid, longitude_grid
     &  , 1
     &  , n_profile_grid, n_level_grid
     &  , p_grid, field_grid
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
!
      WRITE(iu_stdout, '(/a)') 'enter the mode of interpolation.'
1     read(iu_stdin, *, iostat=ios) i_mode
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ unrecognized response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'please re-enter.'
          goto 1
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
          WRITE(iu_err, '(a)') 'please re-enter.'
          goto 1
        ENDIF
      ENDIF
!
!     Perform the interpolation.
      DO l=1, n_profile_in
        DO i=1, n_level_in
          a_in(i)=field_in(l, i)
          p_int(i)=p_in(l, i)
        ENDDO
        l_splined=.false.
        DO i=1, n_level_grid
          CALL interpolate_p(n_level_in, p_int
     &      , a_in, x_in, y_in, y2_in
     &      , p_grid(l, i), field_out(l, i)
     &      , i_mode, l_splined
     &      )
        ENDDO
      ENDDO
!
      WRITE(iu_stdout, '(/a)') 'enter the name of the output file.'
      READ(iu_stdin, '(a)') file_out
      WRITE(iu_stdout, '(/a)') 'enter the name of the output field.'
      READ(iu_stdin, '(a)') name_out
      WRITE(iu_stdout, '(/a)') 'enter the units of the output field.'
      READ(iu_stdin, '(a)') unit_out
      WRITE(iu_stdout, '(/a)') 'enter the long-name of the 
     &                                                 output field.'
      READ(iu_stdin, '(a)') long_out
!
      CALL output_vert_cdf(ierr
     &  , file_out
     &  , n_latitude_in, latitude_in, n_longitude_in, longitude_in
     &  , n_profile_in
     &  , n_level_grid
     &  , trim(name_vert_coord), len(trim(name_vert_coord))
     &  , p_grid
     &  , name_out, 'float', unit_out, long_out, field_out
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!
!
      END
