! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign a single-level field in netCDF format.
!
! Purpose:
!   This subroutine reads a single netCDF-file containing a field on
!   one level and assigns it to an array of the code.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_input_novert_cdf(ierr
     &  , file_name, field_name
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, field
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
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
!           Allowed size for profiles
     &  , nd_latitude
!           Allowed size for latitudes
     &  , nd_longitude
!           Allowed size for longitudes
!
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
     &  , field_name*(*)
!           Name of input field
!
      INTEGER, Intent(INOUT) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
!
!
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
      REAL  (RealK), Intent(OUT) ::
     &    field(nd_profile)
!           Field of data
!
!
!     Local Variables
!
      INTEGER
     &    l
!           Loop variable
     &  , id_lat
!           CDL `index' for latitude
     &  , id_long
!           CDL `index' for longitude
     &  , i_lat
!           Loop variable
     &  , i_long
!           Loop variable
     &  , stride_lat
!           Stride through latitudes
     &  , stride_long
!           Stride through longitudes
!
!
!     Functions invoked:
      INTEGER
     &    calc_cdl_stride
!           Function to calculate the stride in a CDL-array
      EXTERNAL
     &    calc_cdl_stride
!
!     Subroutines called:
      EXTERNAL
     &     read_cdf, assign_horiz_cdl
!
!
!
!     Read the CDL-file.
      CALL read_cdf(ierr
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
!     Only one variable is permitted here.
!
      IF (n_var /= 1) THEN
        WRITE(iu_err, '(2(/a))')
     &    '*** Error: Incorrect number of variables in the file '
     &    , file_name
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Check the horizontal dimensions
      CALL assign_horiz_cdl(ierr, file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_fl
     &  , id_lat, n_latitude, latitude, id_long, n_longitude, longitude
     &  , n_profile
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!
!     The data in the CDL file will be stored with the last
!     `index' changing most rapidly (C-convention), but there
!     will be a common stride between values at the same
!     point.
!
      stride_lat=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_long
     &   , dimension_size, nd_cdl_dimen)
!
!     Assign the data.
      DO i_lat=1, n_latitude
        DO i_long=1, n_longitude
          l=i_long+(i_lat-1)*n_longitude
          field(l)
     &      =data_fl(1+stride_lat*(i_lat-1)+stride_long*(i_long-1), 1)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END
