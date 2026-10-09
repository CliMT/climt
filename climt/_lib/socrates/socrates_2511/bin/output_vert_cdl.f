! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to control the writing of a CDL-file of a field.
!
! Purpose:
!   This subroutine receives a 3-D array and sets up a call to
!   the routines to write it to a CDL-file.
!
! Method:
!   Straightforward; but note the convention that the fields
!   within the program are stored with longitude changing most
!   rapidly.
!
!- ---------------------------------------------------------------------
      SUBROUTINE output_vert_cdl(ierr
     &  , file_name
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_level, name_vert_coord, len_vert_coord, p_level
     &  , field_name, field_type, field_unit, field_long, field
     &  , nd_profile, nd_latitude, nd_longitude, nd_top, nd_bottom
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
     &  , nd_top
!           Declared top level of the field
     &  , nd_bottom
!           Declared bottom level of the field
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
     &  , n_level
!           Number of levels
!
      REAL  (RealK), Intent(IN) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
      CHARACTER !, Intent(IN)
     &    name_vert_coord*(*)
!           Name of vertical coordinate
      INTEGER, Intent(IN) ::
     &    len_vert_coord
!           Length of name of vertical coordinate
!
      REAL  (RealK), Intent(IN) ::
     &    p_level(nd_profile, nd_top: nd_bottom)
!           Atmsopheric pressures
!
!     The field:
      CHARACTER !, Intent(IN)
     &    field_name*(*)
!           Name of the field to be written
     &  , field_type*(*)
!           Type of the field to be written
     &  , field_unit*(*)
!           Unit of the field to be written
     &  , field_long*(*)
!           Long title of the field to be written
      REAL  (RealK), Intent(IN) ::
     &    field(nd_profile, nd_top: nd_bottom)
!           Field to be written out
!
!
!     Local Variables
!
      INTEGER
     &    i
!           Loop variable
     &  , l
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
      n_dimension=3
      dimension_name(1)='lon'
      dimension_name(2)='lat'
      dimension_type(1)='float'
      dimension_type(2)='float'
      dimension_unit(1)='degree'
      dimension_unit(2)='degree'
      dimension_long(1)='longitude'
      dimension_long(2)='latitude'
      dimension_size(1)=n_longitude
      dimension_size(2)=n_latitude
      DO l=1, n_longitude
        dimension_array_fl(l, 1)=longitude(l)
      ENDDO
      DO l=1, n_latitude
        dimension_array_fl(l, 2)=latitude(l)
      ENDDO
!
!     Vertical coordinate:
      IF (name_vert_coord(1:len_vert_coord) == 'plev') THEN
        dimension_name(3)='plev'
        dimension_type(3)='float'
        dimension_unit(3)='pa'
        dimension_long(3)='pressure'
        dimension_size(3)=n_level
        DO i=1, n_level
          dimension_array_fl(i, 3)=p_level(1, i-1+nd_top)
        ENDDO
      ELSE IF (name_vert_coord(1:len_vert_coord) == 'level') THEN
        dimension_name(3)='level'
        dimension_type(3)='int'
        dimension_unit(3)='none'
        dimension_long(3)='level'
        dimension_size(3)=n_level
        DO i=1, n_level
          dimension_array_int(i, 3)=i-1+nd_top
        ENDDO
      ELSE
        WRITE(iu_err, '(/a)')
     &    '*** Error: An illegal vertical coordinate was specified'
     &    , 'for the output field.'
      ENDIF
!
!     Set common properties of the variables.
      IF (name_vert_coord(1:len_vert_coord) == 'plev') THEN
        n_var=1
        n_dimension_var=3
        list_dimension_var(1, 1)=1
        list_dimension_var(2, 1)=2
        list_dimension_var(3, 1)=3
      ELSE IF (name_vert_coord(1:len_vert_coord) == 'level') THEN
        n_var=2
        n_dimension_var=3
        list_dimension_var(1, 1:2)=1
        list_dimension_var(2, 1:2)=2
        list_dimension_var(3, 1:2)=3
        var_name(1)='plev'
        var_type(1)='float'
        var_unit(1)='pa'
        var_long(1)='pressure'
        n_data(1)=n_profile*n_level
        DO i=1, n_level
          DO l=1, n_profile
            data_fl(l+(i-1)*n_profile, 1)=p_level(l, i-1+nd_top)
          ENDDO
        ENDDO
      ENDIF
!
!
      var_name(n_var)=field_name
      var_type(n_var)=field_type
      var_unit(n_var)=field_unit
      var_long(n_var)=field_long
      n_data(n_var)=n_profile*n_level
      DO i=1, n_level
        DO l=1, n_profile
          data_fl(l+(i-1)*n_profile, n_var)=field(l, i-1+nd_top)
        ENDDO
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
