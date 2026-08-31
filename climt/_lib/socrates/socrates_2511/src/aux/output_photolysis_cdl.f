! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to control the writing of CDL-files of mean radiances.
!
! Purpose:
!   This subroutine receives mean radiances as input and calls
!   a routines to write them to CDL-files.
!
! Method:
!   Straightforward, but note the use of dimensions: ND_PROFILE is
!   used throughout since this routine will be called only when the
!   number of profiles where mean radiances are required will be equal
!   to the full number.
!
!- ---------------------------------------------------------------------
      SUBROUTINE output_photolysis_cdl(ierr
     &  , base_name, length_name
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, n_layer
     &  , n_viewing_level, viewing_level
     &  , n_channel, photolysis
     &  , nd_profile, nd_latitude, nd_longitude, nd_layer
     &  , nd_viewing_level, nd_channel
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
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
     &  , n_layer
!           Number of layers
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
     &  , viewing_level(nd_viewing_level)
!           Viewing levels
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
     &    photolysis(nd_profile, nd_viewing_level, nd_channel)
!           Rates of photolysis
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
     &     write_cdl
!
!
!
!
!     There will always be at least three dimensions. An extra fourth
!     will be present if multispectral radiances are produced.
      n_dimension=3
      dimension_name(1)='lon'
      dimension_name(2)='lat'
      dimension_name(3)='level'
      dimension_type(1)='float'
      dimension_type(2)='float'
      dimension_type(3)='int'
      dimension_unit(1)='degree'
      dimension_unit(2)='degree'
      dimension_unit(3)='none'
      dimension_long(1)='longitude'
      dimension_long(2)='latitude'
      dimension_long(3)='viewing level'
      dimension_size(1)=n_longitude
      dimension_size(2)=n_latitude
      dimension_size(3)=n_viewing_level
      DO l=1, n_longitude
        dimension_array_fl(l, 1)=longitude(l)
      ENDDO
      DO l=1, n_latitude
        dimension_array_fl(l, 2)=latitude(l)
      ENDDO
      DO i=1, n_viewing_level
        dimension_array_int(i, 3)=i
      ENDDO
!
!     If multispectral output is to be generated add a dimension
!     for channels.
      IF (n_channel > 1) THEN
        n_dimension=4
        dimension_name(4)='channel'
        dimension_type(4)='int'
        dimension_unit(4)='none'
        dimension_long(4)='spectral channel'
        dimension_size(4)=n_channel
        DO k=1, n_channel
          dimension_array_int(k, 4)=k
        ENDDO
      ENDIF
!
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_photolysis)
!
      n_var=2
!
      var_name(1)='rlev'
      var_type(1)='float'
      var_unit(1)='none'
      var_long(1)='viewing levels'
      n_dimension_var(1)=1
      list_dimension_var(1, 1)=3
      n_data(1)=n_viewing_level
      DO i=1, n_viewing_level
        data_fl(i, 1)=viewing_level(i)
      ENDDO
!
      var_name(2)='jrad'
      var_type(2)='float'
      var_unit(2)='wm-2.str-1'
      var_long(2)='mean radiance'
      list_dimension_var(1, 2)=1
      list_dimension_var(2, 2)=2
      list_dimension_var(3, 2)=3
!
      IF (n_channel == 1) THEN
        n_dimension_var(2)=3
        n_data(2)=n_profile*n_viewing_level
        DO i=1, n_viewing_level
          DO l=1, n_profile
            data_fl(l+(i-1)*n_profile, 2)
     &        =photolysis(l, i, 1)
          ENDDO
        ENDDO
      ELSE
        n_dimension_var(2)=4
        n_data(2)=n_profile*n_viewing_level*n_channel
        list_dimension_var(4, 2)=4
        DO ic=1, n_channel
          DO i=1, n_viewing_level
            DO l=1, n_profile
              point=l+(i-1+(k-1+(ic-1))
     &          *n_viewing_level)*n_profile
              data_fl(point, 2)=photolysis(l, i, ic)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      CALL write_cdl(ierr
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
