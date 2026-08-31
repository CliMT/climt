! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to control the writing of netCDF-files of radiances.
!
! Purpose:
!   This subroutine receives radiances as input and calls
!   a routine to write them to netCDF-files.
!
! Method:
!   Straightforward, but note the use of dimensions: ND_PROFILE is
!   used throughout since this routine will be called only when the
!   number of profiles where radiances are required will be equal to
!   the full number.
!
!- ---------------------------------------------------------------------
      SUBROUTINE output_radiance_cdf(ierr
     &  , base_name, length_name
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, n_direction, direction, azim_0
     &  , n_viewing_level, viewing_level
     &  , n_channel, radiance
     &  , nd_profile, nd_latitude, nd_longitude, nd_layer
     &  , nd_viewing_level, nd_direction, nd_channel
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE rad_ccf, ONLY: pi
      USE gas_list_pcf
      USE rad_pcf
      USE input_head_pcf
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
     &  , nd_layer
!           Allowed size for layers
     &  , nd_viewing_level
!           Allowed size for levels where radiances are known
     &  , nd_direction
!           Allowed size for directions
     &  , nd_channel
!           Allowed size for spectral channels
!
      INTEGER, Intent(IN) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_direction
!           Number of directions
     &  , n_viewing_level
!           Number of levels where the radiance is calculated
     &  , n_channel
!           Number of channels used
!
      REAL  (RealK), Intent(IN) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
     &  , direction(nd_profile, nd_direction, 2)
!           Viewing directions
     &  , viewing_level(nd_viewing_level)
!           Viewing levels
     &  , azim_0(nd_profile)
!           Solar azimuthal angles
!
!
      CHARACTER !, Intent(IN)
     &    base_name*(*)
!           Base name of input file
      INTEGER, Intent(IN) ::
     &    length_name
!           Length of basename
!
      REAL  (RealK), Intent(IN) ::
     &    radiance(nd_profile, nd_viewing_level
     &      , nd_direction, nd_channel)
!           Radiances
!
!
!     Local Variables
!
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
     &  , k
!           Loop variable
     &  , ic
!           Loop variable (spectral channels)
     &  , point
!           Point in CDL data-array
!
      CHARACTER !, Intent(IN)
     &    file_name*80
!           Name of input file
!
!
!     SUBROUTINES CALLED:
      EXTERNAL
     &     write_cdf
!
!
!
!
!     There will always be four dimensions. An extra fifth will
!     be presented if multispectral radiances are produced.
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
      dimension_long(3)='radiance level'
      dimension_long(4)='viewing direction'
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
      DO i=1, n_viewing_level
        dimension_array_int(i, 3)=i
      ENDDO
      DO k=1, n_direction
        dimension_array_int(k, 4)=k
      ENDDO
!
!     If multispectral output is to be generated add a dimension
!     for channels.
      IF (n_channel > 1) THEN
        n_dimension=5
        dimension_name(5)='channel'
        dimension_type(5)='int'
        dimension_unit(5)='none'
        dimension_long(5)='spectral channel'
        dimension_size(5)=n_channel
        DO k=1, n_channel
          dimension_array_int(k, 5)=k
        ENDDO
      ENDIF
!
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_radiance)
!
      n_var=4
!
      var_name(1)='pol'
      var_type(1)='float'
      var_unit(1)='degrees'
      var_long(1)='polar viewing angles'
      n_dimension_var(1)=3
      list_dimension_var(1, 1)=1
      list_dimension_var(2, 1)=2
      list_dimension_var(3, 1)=4
      n_data(1)=n_direction*n_profile
      DO k=1, n_direction
        DO l=1, n_profile
          data_fl(l+(k-1)*n_profile, 1)
     &      =(1.8e+02_RealK/pi)*acos(direction(l, k, 1))
        ENDDO
      ENDDO
      var_name(2)='azim'
      var_type(2)='float'
      var_unit(2)='degrees'
      var_long(2)='azimuthal viewing angles'
      n_dimension_var(2)=3
      list_dimension_var(1, 2)=1
      list_dimension_var(2, 2)=2
      list_dimension_var(3, 2)=4
      n_data(2)=n_direction
!     Rotate the grid back from the solar frame to the user's.
      DO k=1, n_direction
        DO l=1, n_profile
          data_fl(l+(k-1)*n_profile, 2)
     &      =(1.8e+02_RealK/pi)*(direction(l, k, 2)+azim_0(l))
        ENDDO
      ENDDO
      var_name(3)='rlev'
      var_type(3)='float'
      var_unit(3)='none'
      var_long(3)='radiance levels'
      n_dimension_var(3)=1
      list_dimension_var(1, 3)=3
      n_data(3)=n_viewing_level
      DO i=1, n_viewing_level
        data_fl(i, 3)=viewing_level(i)
      ENDDO
!
      var_name(4)='radiance'
      var_type(4)='float'
      var_unit(4)='wm-2.str-1'
      var_long(4)='radiance'
      list_dimension_var(1, 4)=1
      list_dimension_var(2, 4)=2
      list_dimension_var(3, 4)=3
      list_dimension_var(4, 4)=4
!
      IF (n_channel == 1) THEN
        n_dimension_var(4)=4
        n_data(4)=n_profile*n_viewing_level*n_direction
        DO k=1, n_direction
          DO i=1, n_viewing_level
            DO l=1, n_profile
              data_fl(l+(i-1+(k-1)*n_viewing_level)*n_profile, 4)
     &          =radiance(l, i, k, 1)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        n_dimension_var(4)=5
        n_data(4)=n_profile*n_viewing_level*n_direction*n_channel
        list_dimension_var(5, 4)=5
        DO ic=1, n_channel
          DO k=1, n_direction
            DO i=1, n_viewing_level
              DO l=1, n_profile
                point=l+(i-1+(k-1+(ic-1)*n_direction)
     &            *n_viewing_level)*n_profile
                data_fl(point, 4)=radiance(l, i, k, ic)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      CALL write_cdf(ierr
     &  , file_name(1: length_name+1+len_file_suffix)
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
