! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign the optical properties from a CDL-file.
!
! Purpose:
!   This subroutine reads a single CDL-file containing a set of
!   prescribed optical properties on a range of pressure levels.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_input_opt_cdl(ierr
     &  , file_name
     &  , n_band
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_p_level, p_level, n_phase_term
     &  , absorption, scattering, phase_function
     &  , nd_profile, nd_latitude, nd_longitude
     &  , nd_band, nd_phase_term, nd_opt_level
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
     &  , nd_band
!           Allowed size for spectral bands
     &  , nd_phase_term
!           Allowed size for terms in the phase function
     &  , nd_opt_level
!           Allowed size for vertical levels of data
!
!
      CHARACTER !, Intent(IN)
     &    file_name*(*)
!           Name of input file
!
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of spectral bands
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
     &    n_phase_term
!           Number of terms in the phase function
     &  , n_p_level
!           Number of vertical levels
      REAL  (RealK), Intent(OUT) ::
     &    p_level(nd_profile, nd_opt_level)
!           Pressure levels where data are given
     &  , absorption(nd_profile, nd_opt_level, nd_band)
!           Absorption coefficients
     &  , scattering(nd_profile, nd_opt_level, nd_band)
!           Scattering coefficients
     &  , phase_function(nd_profile, nd_opt_level
     &      , nd_phase_term, nd_band)
!           Scattering coefficients
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
     &  , id_band
!           CDL `index' for spectral band
     &  , id_mom
!           CDL `index' for moments of the phase function
     &  , id_vlev
!           CDL `index' for vertical levels
     &  , iv_abs
!           CDL `index' for absorption
     &  , iv_scat
!           CDL `index' for absorption
     &  , iv_phf
!           CDL `index' for absorption
     &  , i_lat
!           Loop variable
     &  , i_long
!           Loop variable
     &  , i_band
!           Loop variable
     &  , i_mom
!           Loop variable
     &  , i_vlev
!           Loop variable
     &  , stride_lat_as
!           Stride through latitudes (absorption and scattering)
     &  , stride_long_as
!           Stride through longitudes (absorption and scattering)
     &  , stride_band_as
!           Stride through spectral bands (absorption and scattering)
     &  , stride_vlev_as
!           Stride through vertical levels (absorption and scattering)
     &  , stride_lat_phf
!           Stride through latitudes (phase function)
     &  , stride_long_phf
!           Stride through longitudes (phase function)
     &  , stride_band_phf
!           Stride through spectral bands (phase function)
     &  , stride_mom_phf
!           Stride through moments of the phase function
     &  , stride_vlev_phf
!           Stride through vertical levels (phase function)
     &  , offset_as
!           Offset in array of data for absorption and scattering
     &  , offset_phf
!           Offset in array of data for the phase function
!
      LOGICAL
     &    l_assigned
!           FLag for assigned sizes of dimensions
      REAL  (RealK)
     &    dummy(nd_cdl_dimen_size)
!           Dummy integer array
!
!     Functions invoked:
      INTEGER
     &    calc_cdl_stride
!           Function to calculte the stride in a CDL-array
      EXTERNAL
     &    calc_cdl_stride
!
!     Subroutines called:
      EXTERNAL
     &     read_cdl, find_dimen_cdl, find_var_cdl
!
!
!
!     Read in the file of optical properties.
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
!     Three variables are expected.
!
      IF (n_var /= 3) THEN
        WRITE(iu_err, '(3(/a))')
     &    '*** Error: The file'
     &    , file_name
     &    , 'does not contain the correct number of variables.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!
!     Check the dimensions
!
!
!     Latitude
!
      l_assigned=(n_latitude > 0)
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'lat', 3, 'float', 5
     &  , l_assigned, n_latitude, .true., .true.
     &  , id_lat, dummy, latitude
     &  , nd_latitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!     Longitude
!
      l_assigned=(n_longitude > 0)
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'lon', 3, 'float', 5
     &  , l_assigned, n_longitude, .true., .true.
     &  , id_long, dummy, longitude
     &  , nd_latitude
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
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
!     Pressure levels
!
!     The number of levels allowed here is variable.
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'plev', 4, 'float', 5
     &  , .false., n_p_level, .true., .true.
     &  , id_vlev, dummy, p_level(1,:)
     &  , nd_opt_level
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )

      p_level=SPREAD(p_level(1,:),1,nd_profile)

      IF (n_p_level < 2) THEN
        WRITE(iu_err, '(3(/a))')
     &    '*** Error: The file', file_name
     &    //'contains too few vertical levels.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Spectral band:
!
!     No actual values are required in this case.
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'band', 4, 'int', 3
     &  , .true., n_band, .false., .true.
     &  , id_band, dummy, dummy
     &  , nd_band
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!     Moments of the phase function:
!
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'mom', 3, 'int', 3
     &  , .false., n_phase_term, .true., .true.
     &  , id_mom, dummy, dummy
     &  , nd_phase_term
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!     Optical variables:
!
!     Absorption:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'abs', 3, 'float', 5
     &  , .true.
     &  , iv_abs
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!     Scattering:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'scat', 4, 'float', 5
     &  , .true.
     &  , iv_scat
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!     Phase function:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'phf', 3, 'float', 5
     &  , .true.
     &  , iv_phf
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     The data in the CDL file will be stored with the last
!     `index' changing most rapidly (C-convention), but there
!     will be a common stride between values at the same
!     point. The phase function will have different strides
!     because it loops over moments.
!
      stride_lat_as=calc_cdl_stride(n_dimension_var(iv_abs)
     &   , list_dimension_var(1, iv_abs), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long_as=calc_cdl_stride(n_dimension_var(iv_abs)
     &   , list_dimension_var(1, iv_abs), id_long
     &   , dimension_size, nd_cdl_dimen)
      stride_vlev_as=calc_cdl_stride(n_dimension_var(iv_abs)
     &   , list_dimension_var(1, iv_abs), id_vlev
     &   , dimension_size, nd_cdl_dimen)
      stride_band_as=calc_cdl_stride(n_dimension_var(iv_abs)
     &   , list_dimension_var(1, iv_abs), id_band
     &   , dimension_size, nd_cdl_dimen)
!
      stride_lat_phf=calc_cdl_stride(n_dimension_var(iv_phf)
     &   , list_dimension_var(1, iv_phf), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long_phf=calc_cdl_stride(n_dimension_var(iv_phf)
     &   , list_dimension_var(1, iv_phf), id_long
     &   , dimension_size, nd_cdl_dimen)
      stride_vlev_phf=calc_cdl_stride(n_dimension_var(iv_phf)
     &   , list_dimension_var(1, iv_phf), id_vlev
     &   , dimension_size, nd_cdl_dimen)
      stride_band_phf=calc_cdl_stride(n_dimension_var(iv_phf)
     &   , list_dimension_var(1, iv_phf), id_band
     &   , dimension_size, nd_cdl_dimen)
      stride_mom_phf=calc_cdl_stride(n_dimension_var(iv_phf)
     &   , list_dimension_var(1, iv_phf), id_mom
     &   , dimension_size, nd_cdl_dimen)
!
!     Assign the data.
      DO i_vlev=1, n_p_level
        DO i_band=1, n_band
          DO i_lat=1, n_latitude
            DO i_long=1, n_longitude
!
              l=i_long+(i_lat-1)*n_longitude
              offset_as=stride_vlev_as*(i_vlev-1)
     &                 +stride_band_as*(i_band-1)
     &                 +stride_lat_as*(i_lat-1)
     &                 +stride_long_as*(i_long-1)
!
              absorption(l, i_vlev, i_band)
     &          =data_fl(1+offset_as, iv_abs)
              scattering(l, i_vlev, i_band)
     &          =data_fl(1+offset_as, iv_scat)
!
              DO i_mom=1, n_phase_term
                offset_phf=stride_vlev_phf*(i_vlev-1)
     &                    +stride_band_phf*(i_band-1)
     &                    +stride_lat_phf*(i_lat-1)
     &                    +stride_long_phf*(i_long-1)
     &                    +stride_mom_phf*(i_mom-1)
                phase_function(l, i_vlev, i_mom, i_band)
     &            =data_fl(1+offset_phf, iv_phf)
              ENDDO
!
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END
