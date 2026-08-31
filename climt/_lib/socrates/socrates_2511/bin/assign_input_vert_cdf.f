! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign a netCDF field with vertical structure.
!
! Purpose:
!   This subroutine checks a single netCDF field with vertical structure
!   which is read in and puts it into an array.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_input_vert_cdf(ierr
     &  , file_name, field_name, l_vert_coord, name_vert_coord
     &  , l_check_vert_dimen, n_vert_expected, l_vert_assignable
     &  , n_latitude, latitude, n_longitude, longitude
     &  , i_start
     &  , n_profile, n_level
     &  , p, field
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
!           Top vertical dimension of field to be filled
     &  , nd_bottom
!           Bottom vertical dimension of field to be filled
!
      INTEGER, Intent(IN) ::
     &    n_vert_expected
!           Expected number of vertical levels
     &  , i_start
!           Starting point in array
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
     &  , field_name*(*)
!           Name of input field.
      CHARACTER !, Intent(INOUT)
     &    name_vert_coord*(*)
!           Name of vertical coordinate
      LOGICAL, Intent(IN) ::
     &    l_check_vert_dimen
!           Flag for to check vertical dimension
     &  , l_vert_assignable
!           Flag to allow vertical coordinate to be assigned 
!           by this field
!
      LOGICAL, Intent(INOUT) ::
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
!
      INTEGER, Intent(INOUT) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
!
!
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
     &  , p(nd_profile, nd_top: nd_bottom)
!           Atmsopheric pressures
!
      REAL  (RealK), Intent(INOUT) ::
     &    field(nd_profile, nd_top: nd_bottom)
!           Field of input data
!
!
!     Local Variables
!
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
     &  , i_lat
!           Loop variable
     &  , i_long
!           Loop variable
     &  , id_lat
!           CDL index for latitude
     &  , id_long
!           CDL index for longitude
     &  , iv_plev
!           CDL index for pressure levels
     &  , iv_fld
!           CDL index for field
     &  , idv
!           CDL index
     &  , stride_i
!           Stride over index I
     &  , stride_lat
!           Stride through latitudes
     &  , stride_long
!           Stride through longitudes
!
      LOGICAL
     &    l_loop
!           Logical controlling WHILE loops: this is required
!           to expand the evaluation of double logical tests 
!           which can result in out-of-bounds problems on some
!           machines.
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
!     Read the netCDF-field.
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
!     Check the horizontal dimensions.
      CALL assign_horiz_cdl(ierr, file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_fl
     &  , id_lat, n_latitude, latitude, id_long, n_longitude
     &  , longitude
     &  , n_profile
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
      IF (ierr /= i_normal) RETURN
!
!
!     Vertical dimensions.
!
      IF (l_vert_coord) THEN
!
!       Check that the vertical coordinate in the file matches
!       that assumed.
        idv=1
        l_loop=(idv <= n_dimension)
        IF (l_loop) l_loop=(l_loop.AND.
     &             (dimension_name(idv)(1:len(name_vert_coord)) /= 
     &              name_vert_coord(1:len(name_vert_coord))) )
        DO WHILE (l_loop)
          idv=idv+1
          l_loop=(idv <= n_dimension)
          IF (l_loop) l_loop=(l_loop.AND.
     &               (dimension_name(idv)(1:len(name_vert_coord)) /= 
     &                name_vert_coord(1:len(name_vert_coord))) )
        ENDDO
!
        IF (idv > n_dimension) THEN
          WRITE(iu_err, '(/a, /a, a, a1)')
     &      '*** Error: Mismatch in type of vertical coordinate '
     &      , '           while reading '//field_name//'.'
          ierr=i_err_fatal
          RETURN
        ELSE IF (l_check_vert_dimen) THEN
!         Check that the vertical dimension is as expected
          IF (dimension_size(idv) /= n_vert_expected) THEN
            WRITE(iu_err, '(/a, /a)')
     &        '*** Error: Mismatch in number of vertical data '
     &        , '           while reading '//field_name//'.'
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDIF
!
!       Prepare to assign the vetical coordinates if required. 
        IF (l_vert_assignable) THEN
          n_level=dimension_size(idv)
          iv_plev=idv
        ENDIF
!
!
      ELSE
!
!       Find what coordinate exists in the file if this field may 
!       be used to set vertical coordinates.
        IF (l_vert_assignable) THEN
          idv=1
          DO WHILE ( (idv <= n_dimension).AND.
     &               (.NOT.l_vert_coord) )
!
            IF (dimension_name(idv)(1:4) == 'hlev') THEN
              WRITE(iu_err, '(/a)') 
     &          '*** Sorry, hybrid coordinates are not yet available.'
              RETURN
!             NAME_VERT_COORD(1:LEN(DIMENSION_NAME(IDV)))
!    &          =DIMENSION_NAME(IDV)(1:LEN(DIMENSION_NAME(IDV)))
!             L_VERT_COORD=.TRUE.
!             IF (DIMENSION_TYPE(IDV)(1:5).NE.'float') THEN
!               WRITE(IU_ERR, '(/A)')
!    &            '*** Error: Illegal type for hybrid coordinates '
!    &            //'           while reading '//FIELD_NAME//'.'
!               IERR=I_ERR_FATAL
!               RETURN
!             ENDIF
!
            ELSE IF (dimension_name(idv)(1:4) == 'plev') THEN
              name_vert_coord(1:len(dimension_name(idv)))
     &          =dimension_name(idv)(1:len(dimension_name(idv)))
              l_vert_coord=.true.
              IF (dimension_type(idv)(1:5) /= 'float') THEN
                WRITE(iu_err, '(/a)')
     &            '*** Error: Illegal type for pressure coordinates '
     &            //'           while reading '//field_name
                ierr=i_err_fatal
                RETURN
              ENDIF
!
!             Set the levels, but record the pointer to the position
!             in the input array for later assignment of actual pressures
!             when we know what the striide will be.
              n_level=dimension_size(idv)
              iv_plev=idv
!
            ELSE IF (dimension_name(idv)(1:5) == 'layer') THEN
              name_vert_coord(1:len(dimension_name(idv)))
     &          =dimension_name(idv)(1:len(dimension_name(idv)))
              l_vert_coord=.true.
              IF (dimension_type(idv)(1:3) /= 'int') THEN
                WRITE(iu_err, '(/a)')
     &            '*** Error: Illegal type for layer number '
     &            //'           while reading '//field_name
                ierr=i_err_fatal
                RETURN
              ENDIF
!
!             No values of coordinates are explicitly passed 
!             back in this case.
              n_level=dimension_size(idv)
!
            ELSE IF (dimension_name(idv)(1:5) == 'level') THEN
              name_vert_coord(1:len(dimension_name(idv)))
     &          =dimension_name(idv)(1:len(dimension_name(idv)))
              l_vert_coord=.true.
              IF (dimension_type(idv)(1:3) /= 'int') THEN
                WRITE(iu_err, '(/a)')
     &            '*** Error: Illegal type for level number '
     &            //'           while reading '//field_name
                ierr=i_err_fatal
                RETURN
              ENDIF
!
!             No values of coordinates are explicitly passed 
!             back in this case.
              n_level=dimension_size(idv)
!
            ELSE
!
              idv=idv+1
!
            ENDIF
!
          ENDDO
!
        ELSE
!
          WRITE(iu_err, '(/a, /a)')
     &      '*** Error: The field '//field_name
     &      , '           cannot be used to '
     &      //'assign vertical coordinates.'
          ierr=i_err_fatal
          RETURN
!
        ENDIF
!
      ENDIF
!
!
!     Read the data. Two basic forms of file are possible: either
!     plev is a coordinate, or it may be amongst the data (this
!     must be the case if we have multiple profiles on different
!     pressure levels).
      IF (n_var == 2) THEN
        CALL find_var_cdl(ierr
     &    , file_name
     &    , n_var, var_name, var_type
     &    , 'plev', 4, 'float', 5
     &    , .true.
     &    , iv_plev
     &    , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &    )
        IF (ierr /= i_normal) RETURN
!       This is a trick which works with just two variables. If the
!       code is extended, we will need a more generic relationship.
        iv_fld=3-iv_plev
      ELSE
        iv_fld=1
      ENDIF
!
!     The data in the CDL file will be stored with the last
!     index changing most rapidly (C-convention), but there
!     will be a common stride between values at the same
!     point.
!
      stride_lat=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), id_long
     &   , dimension_size, nd_cdl_dimen)
      stride_i=calc_cdl_stride(n_dimension_var(1)
     &   , list_dimension_var(1, 1), idv
     &   , dimension_size, nd_cdl_dimen)
!
!     Assign the pressure levels if required.
      IF (l_vert_assignable) THEN
        IF (n_var == 1) THEN
          DO i=1, n_level
            DO l=1, n_profile
              p(l, i_start+i-1)
     &          =dimension_array_fl(i, iv_plev)
            ENDDO
          ENDDO
        ELSE IF (n_var == 2) THEN
          DO i=1, n_level
            DO i_lat=1, n_latitude
              DO i_long=1, n_longitude
                l=i_long+(i_lat-1)*n_longitude
                p(l, i-1+i_start)
     &            =data_fl(1+stride_i*(i-1)
     &            +stride_lat*(i_lat-1)
     &            +stride_long*(i_long-1), iv_plev)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
!
!           
!     Assign the data. 
      DO i=1, n_level
        DO i_lat=1, n_latitude
          DO i_long=1, n_longitude
            l=i_long+(i_lat-1)*n_longitude
            field(l, i-1+i_start)
     &        =data_fl(1+stride_i*(i-1)
     &        +stride_lat*(i_lat-1)+stride_long*(i_long-1), iv_fld)
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END
