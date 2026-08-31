! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to prepare a CDL-file of viewing fields.
!
! Method:
!   Straightforward; but note the convention that the fields
!   within the program are stored with longitude changing most
!   rapidly.
!
!- ---------------------------------------------------------------------
      SUBROUTINE output_view_cdl(ierr
     &  , file_name
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_viewing_level, viewing_level
     &  , n_direction, direction
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_viewing_level, nd_direction
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
     &  , nd_viewing_level
!           Allowed size for viewing levels
     &  , nd_direction
!           Allowed size for viewing directions
!
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of the file to be written
!
!     Dimensions:
      INTEGER, Intent(IN) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_viewing_level
!           Number of viewing levels
     &  , n_direction
!           Number of viewing directions
!
      REAL  (RealK), Intent(IN) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
     &  , viewing_level(nd_viewing_level)
!           Viewing levels
     &  , direction(nd_profile, nd_direction, 2)
!           Viewing directions
!
!
!     Local Variables
!
      INTEGER
     &    l
!           Loop variable
     &  , k
!           Loop variable
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     write_cdl
!
!
!
!     Set up the dimensions.
      n_dimension=4
      dimension_name(1)='lon'
      dimension_name(2)='lat'
      dimension_name(3)='level'
      dimension_name(4)='direction'
      dimension_type(1)='float'
      dimension_type(2)='float'
      dimension_type(3)='int'
      dimension_type(4)='int'
      dimension_unit(1)='degree'
      dimension_unit(2)='degree'
      dimension_unit(3)='none'
      dimension_unit(4)='none'
      dimension_long(1)='longitude'
      dimension_long(2)='latitude'
      dimension_long(3)='level'
      dimension_long(4)='direction'
      dimension_size(1)=n_longitude
      dimension_size(2)=n_latitude
      dimension_size(3)=n_viewing_level
      dimension_size(4)=n_direction
      DO l=1, n_longitude
        dimension_array_fl(l, 1)=longitude(l)
      ENDDO
      DO l=1, n_latitude
        dimension_array_fl(l, 2)=latitude(l)
      ENDDO
      DO k=1, n_viewing_level
        dimension_array_int(k, 3)=k
      ENDDO
      DO k=1, n_direction
        dimension_array_int(k, 4)=k
      ENDDO
!
!     Set common properties of the variables.
      n_var=3
!
      n_dimension_var(1)=3
      list_dimension_var(1, 1)=1
      list_dimension_var(2, 1)=2
      list_dimension_var(3, 1)=4
      var_name(1)='pol'
      var_type(1)='float'
      var_unit(1)='degree'
      var_long(1)='polar viewing angle'
      n_data(1)=n_profile*n_direction
      DO k=1, n_direction
        DO l=1, n_profile
          data_fl(l+(k-1)*n_profile, 1)=direction(l, k, 1)
        ENDDO
      ENDDO
!
      n_dimension_var(2)=3
      list_dimension_var(1, 2)=1
      list_dimension_var(2, 2)=2
      list_dimension_var(3, 2)=4
      var_name(2)='azim'
      var_type(2)='float'
      var_unit(2)='degree'
      var_long(2)='azimuthal viewing angle'
      n_data(2)=n_profile*n_direction
      DO k=1, n_direction
        DO l=1, n_profile
          data_fl(l+(k-1)*n_profile, 2)=direction(l, k, 2)
        ENDDO
      ENDDO
!
      n_dimension_var(3)=1
      list_dimension_var(1, 3)=3
      var_name(3)='rlev'
      var_type(3)='float'
      var_unit(3)='none'
      var_long(3)='viewing level'
      n_data(3)=n_viewing_level
      DO k=1, n_viewing_level
        data_fl(k, 3)=viewing_level(k)
      ENDDO
!
      CALL write_cdl(ierr
     &  , file_name
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  , n_dimension, dimension_name, dimension_type
     &  , dimension_unit
     &  , dimension_long, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , n_var, var_name, var_type, var_unit, var_long
     &  , n_dimension_var, list_dimension_var
     &  , n_data, data_int, data_fl
     &  )
      IF (ierr /= i_normal) RETURN
!
!
!
      RETURN
      END
