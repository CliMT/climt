! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign the horizontal dimensions from a CDL-file.
!
! Purpose:
!   This subroutine checks the latitide and longitude coordinates,
!   reads from a CDL file and assigns the latitudes and longitudes
!   and the total number of points, if required.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_horiz_cdl(ierr, file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_fl
     &  , id_lat, n_latitude, latitude, id_long, n_longitude, longitude
     &  , n_profile
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
!           Name of current file
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
      INTEGER, Intent(OUT) ::
     &    id_lat
!           CDL `index' for latitude
     &  , id_long
!           CDL `index' for longitude
!
!
!     Local Variables
!
      INTEGER
     &    l
!           Loop variable
      LOGICAL
     &    l_found_lat
!           Flag for latitude 
     &  , l_found_long
!           Flag for longitude
!
!
!
!     Check the latitude coordinate:
!
      l_found_lat=.false.
      id_lat=1
      DO WHILE ( (id_lat <= n_dimension).AND.(.NOT.l_found_lat) )
        IF (dimension_name(id_lat)(1:3) == 'lat') THEN
!
          l_found_lat=.true.
          IF (dimension_type(id_lat)(1:5) /= 'float') THEN
            WRITE(iu_err, '(/a, /a)')
     &        '*** Error: Illegal type for latitude while reading '
     &        , file_name//'.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!         Check the size of the fields.
          IF (n_latitude > 0) THEN
            IF (n_latitude /= dimension_size(id_lat)) THEN
              WRITE(iu_err, '(/a, /a)')
     &          '*** Error: Number of latitudes is not '
     &          //'consistent with existing value after '
     &          , 'reading '//file_name//'.'
              ierr=i_err_fatal
              RETURN
            ENDIF
          ELSE
            n_latitude=dimension_size(id_lat)
            DO l=1, n_latitude
              latitude(l)=dimension_array_fl(l, id_lat)
            ENDDO
          ENDIF
!
        ELSE
!
          id_lat=id_lat+1
!
        ENDIF
      ENDDO
!
      IF (.NOT.l_found_lat) THEN
        WRITE(iu_err, '(/a, /a)')
     &    '*** Error: No latitude coordinate found while reading '
     &    , file_name//'.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!
!     Longitude
!
      l_found_long=.false.
      id_long=1
      DO WHILE ( (id_long <= n_dimension).AND.(.NOT.l_found_long) )
        IF (dimension_name(id_long)(1:3) == 'lon') THEN
!
          l_found_long=.true.
          IF (dimension_type(id_long)(1:5) /= 'float') THEN
            WRITE(iu_err, '(/a)')
     &        '*** Error: Illegal type for longitude '
     &        //'while reading '//file_name//'.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!         Check the size of the fields.
          IF (n_longitude > 0) THEN
            IF (n_longitude /= dimension_size(id_long)) THEN
              WRITE(iu_err, '(/a, /a)')
     &          '*** Error: Number of longitudes is not '
     &          //'consistent with existing value after '
     &          , 'reading '//file_name//'.'
              ierr=i_err_fatal
              RETURN
            ENDIF
          ELSE
!           Assign the longitudes.
            n_longitude=dimension_size(id_long)
            DO l=1, n_longitude
              longitude(l)=dimension_array_fl(l, id_long)
            ENDDO
          ENDIF
!
        ELSE
!
          id_long=id_long+1
!
        ENDIF
      ENDDO
!
      IF (.NOT.l_found_long) THEN
        WRITE(iu_err, '(/a, /a)')
     &    '*** Error: No longitude coordinate found while reading '
     &    , file_name//'.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Now set the overall horizontal size.
      n_profile =  n_latitude * n_longitude
      IF (n_profile > nd_profile) THEN
        WRITE(iu_err, '(a)') 
     &    '*** Error: Too many atmospheric profiles.'
        ierr = i_err_fatal
        RETURN
      ENDIF
!
!
!
      RETURN
      END
