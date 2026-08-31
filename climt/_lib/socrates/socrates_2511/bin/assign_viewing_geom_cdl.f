! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign the viewing geometry from a CDL-file.
!
! Purpose:
!   This subroutine reads a single CDL-file containing a set
!   of polar and azimuthal viewing angles (in degrees) and a
!   set of levels where radiances are required. The levels are
!   specified as real numbers which will be associated with the
!   actual integral numbers of levels; hence 0.0 is the top, 0.5
!   is halfway through the first layer, 1.0 is at the base of the
!   first layer and 3.4 is 40% of the way through the fourth layer.
!
! Method:
!   Straightforward, but note the use of dimensions: ND_PROFILE is
!   used throughout since this routine will be called only when the
!   number of profiles where radiances are required will be equal to
!   the full number.
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_viewing_geom_cdl(ierr
     &  , file_name
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile
     &  , n_direction, direction
     &  , n_viewing_level, viewing_level
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_direction, nd_viewing_level
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE rad_ccf, ONLY: pi
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!     INCLUDE HEADER FILES.
!
!
!
!     Declaration of variables:
!
      INTEGER   !,Intent(OUT)
     &    ierr
!            Error flag
!
      INCLUDE 'cdl_struc.finc'
!
!     Sizes of arrays
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Allowed size for atmospheric profiles
     &  , nd_latitude
!           Allowed size for latitudes
     &  , nd_longitude
!           Allowed size for longitudes
     &  , nd_direction
!           Allowed size for viewing directions
     &  , nd_viewing_level
!           Allowed size for levels where radiances are required
!
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
!
!     Atmospheric structure
      INTEGER, Intent(INOUT) ::
     &    n_profile
!           Number of atmospheric profiles
     &  , n_latitude
!           Number of atmospheric profiles
     &  , n_longitude
!           Number of atmospheric profiles
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
!     Viewing directions
      INTEGER, Intent(OUT) ::
     &    n_direction
!           Number of directions
      REAL  (RealK), Intent(OUT) ::
     &    direction(nd_profile, nd_direction, 2)
!           Viewing directions (converted to the cosine of the polar
!           angle and the azimuthal angles in radians: the direction
!           is towards the viewing point)
!
!     Viewing levels
      INTEGER, Intent(OUT) ::
     &    n_viewing_level
!           Number of levels where radiances are required
      REAL  (RealK), Intent(OUT) ::
     &    viewing_level(nd_viewing_level)
!           Levels where radiances are required
!
!
!     Local Variables
!
      INTEGER
     &    i_lat
!           Loop variable
     &  , i_long
!           Loop variable
     &  , i_dir
!           Loop variable
     &  , k
!           Derived index
     &  , l
!           Derived index
      INTEGER
     &    id_lat
!           CDL `index' for latitudes
     &  , id_long
!           CDL `index' for longitudes
     &  , id_dir
!           CDL `index' for viewing directions
     &  , id_vlev
!           CDL `index' for viewing levels
     &  , iv_pol
!           Number of the variable holding polar angles
     &  , iv_azim
!           Number of the variable holding azimuthal angles
     &  , iv_rlev
!           Number of the variable holding viewing levels
     &  , stride_lat
!           Stride through latitudes
     &  , stride_long
!           Stride through longitudes
     &  , stride_dir
!           Stride through directions
     &  , i_dummy(nd_cdl_dimen_size)
!           Dummy integer array
      REAL  (RealK) ::
     &    dummy(nd_cdl_dimen_size)
!           Dummy real array
!
!
!     Functions invoked:
      INTEGER
     &    calc_cdl_stride
!           Function to calculate the stride in a CDL-array
      EXTERNAL
     &    calc_cdl_stride
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     read_cdl, assign_horiz_cdl, find_dimen_cdl, find_var_cdl
!
!
!
!     Read in the file of viewing angles.
      CALL read_cdl(ierr
     &  , file_name
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit, dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Check the contents of the file.
!
!     Three variables are required here.
!
      IF (n_var /= 3) THEN
        WRITE(iu_err, '(3(/a))')
     &    '*** Error: The file '
     &    , file_name
     &    , 'is incomplete.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Check the horizontal directions:
      CALL assign_horiz_cdl(ierr, file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_fl
     &  , id_lat, n_latitude, latitude, id_long, n_longitude, longitude
     &  , n_profile
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!     Viewing Directions:
!
!     We presume that no directions have been assigned. No actual
!     values are required here, merely the number of directions.
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'direction', 9, 'int', 3
     &  , .false., n_direction, .true., .true.
     &  , id_dir, i_dummy, dummy
     &  , nd_direction
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!
!     Viewing Levels:
!
!     We presume that no viewing levels have been assigned. (This
!     is a dimension, the list of the number of levels, not the
!     actual levels themselves).
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'level', 5, 'int', 3
     &  , .false., n_viewing_level, .true., .true.
     &  , id_vlev, i_dummy, dummy
     &  , nd_viewing_level
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Polar viewing angles:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'pol', 3, 'float', 5
     &  , .true.
     &  , iv_pol
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Azimuthal viewing angles:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'azim', 4, 'float', 5
     &  , .true.
     &  , iv_azim
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!     Levels for calculating radiances:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'rlev', 4, 'float', 5
     &  , .true.
     &  , iv_rlev
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!     The data in the CDL file will be stored with the last
!     `index' changing most rapidly (C-convention), but there
!     will be a common stride between values at the same
!     point. Strides are calculated for each dimension.
!
      stride_lat=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_long
     &   , dimension_size, nd_cdl_dimen)
      stride_dir=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_dir
     &   , dimension_size, nd_cdl_dimen)
!
!
!     Assign the data, converting to the cosine of the polar angle
!     and the azimuthal angle in radians.
      DO i_lat=1, n_latitude
        DO i_long=1, n_longitude
          l=i_long+(i_lat-1)*n_longitude
          DO i_dir=1, n_direction
            k=1+stride_lat*(i_lat-1)+stride_long*(i_long-1)
     &        +stride_dir*(i_dir-1)
            direction(l, i_dir, 1)
     &        =cos(pi*data_fl(k, iv_pol)/1.8e+02_RealK)
            direction(l, i_dir, 2)
     &        =pi*data_fl(k, iv_azim)/1.8e+02_RealK
          ENDDO
        ENDDO
      ENDDO
!
      DO l=1, n_viewing_level
        viewing_level(l)=data_fl(l, iv_rlev)
      ENDDO
!
!
!
      RETURN
      END
