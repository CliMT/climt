! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign single scattering properties from a CDL-file.
!
! Purpose:
!   This subroutine reads a single CDL-file containing a set of
!   single scattering properties.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE assign_input_ss_cdl(ierr
     &  , file_name
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_layer, n_phase_term
     &  , tau, omega, phase_function
     &  , l_forward, forward_scatter
     &  , l_forward_solar, forward_solar
     &  , nd_profile, nd_latitude, nd_longitude, nd_layer
     &  , nd_phase_term
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
      INTEGER, Intent(INOUT) ::
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
!           Allowed size for vertical layers
     &  , nd_phase_term
!           Allowed size for terms in the phase function
!
!
      CHARACTER, Intent(IN) ::
     &    file_name*(*)
!           Name of input file
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
      REAL  (RealK), Intent(OUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
      LOGICAL, Intent(OUT) ::
     &    l_forward
!           Flag for forward scattering
     &  , l_forward_solar
!           Flag for forward scattering for the solar beam
      INTEGER, Intent(OUT) ::
     &    n_phase_term
!           Number of terms in the phase function
     &  , n_layer
!           Number of vertical levels
      REAL  (RealK), Intent(OUT) ::
     &    tau(nd_profile, nd_layer)
!           Optical depths
     &  , omega(nd_profile, nd_layer)
!           Albedos of single scattering
     &  , forward_scatter(nd_profile, nd_layer)
!           Forward scattering fraction
     &  , forward_solar(nd_profile, nd_layer)
!           Forward scattering fraction for the solar beam
     &  , phase_function(nd_profile, nd_layer, nd_phase_term)
!           Phase function
!
!
!     Local Variables
!
      INTEGER
     &    v_level(nd_profile, nd_layer)
!           Numbers of layers (required only for consistency)
      INTEGER
     &    l
!           Loop variable
     &  , id_lat
!           CDL `index' for latitude
     &  , id_long
!           CDL `index' for longitude
     &  , id_mom
!           CDL `index' for moments of the phase function
     &  , id_vlev
!           CDL `index' for vertical levels
     &  , iv_tau
!           CDL `index' for absorption
     &  , iv_omega
!           CDL `index' for the albedo of single scattering
     &  , iv_frwd
!           CDL `index' for forward scattering
     &  , iv_sfrwd
!           CDL `index' for solar forward scattering
     &  , iv_phf
!           CDL `index' for absorption
     &  , i_lat
!           Loop variable
     &  , i_long
!           Loop variable
     &  , i_mom
!           Loop variable
     &  , i_vlev
!           Loop variable
     &  , stride_lat_to
!           Stride through latitudes (optical depth and albedo)
     &  , stride_long_to
!           Stride through longitudes (optical depth and albedo)
     &  , stride_vlev_to
!           Stride through vertical levels (optical depth and albedo)
     &  , stride_lat_phf
!           Stride through latitudes (phase function)
     &  , stride_long_phf
!           Stride through longitudes (phase function)
     &  , stride_mom_phf
!           Stride through moments of the phase function
     &  , stride_vlev_phf
!           Stride through vertical levels (phase function)
     &  , offset_to
!           Offset in array of data for optical depth and albedo
     &  , offset_phf
!           Offset in array of data for the phase function
!
      LOGICAL
     &    l_assigned
!           FLag for assigned sizes of dimensions
      INTEGER
     &    i_dummy(nd_cdl_dimen_size)
!           Dummy integer array
      REAL  (RealK)
     &    dummy(nd_cdl_dimen_size)
!           Dummy real array
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
!     Three, four or five variables are expected 
!     (the forward scattering fractions may or may not be
!     present).
      l_forward=(n_var == 4).OR.(n_var == 5)
      l_forward_solar=(n_var == 5)
!
      IF ( (n_var < 3).OR.(n_var > 5) ) THEN
        WRITE(iu_err, '(3(/a))')
     &    '*** Error: The file'
     &    , file_name
     &    , 'does not contain the a valid number of variables.'
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
     &  , id_lat, i_dummy, latitude
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
     &  , id_long, i_dummy, longitude
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
!     Vertical levels
!
!     The number of levels allowed here is variable.
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'layer', 5, 'int', 3
     &  , .false., n_layer, .true., .true.
     &  , id_vlev, v_level, dummy
     &  , nd_layer
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
      IF (n_layer < 1) THEN
        WRITE(iu_err, '(3(/a))')
     &    '*** Error: The file', file_name
     &    //'contains too few vertical layers.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Moments of the phase function:
!
      CALL find_dimen_cdl(ierr
     &  , file_name
     &  , n_dimension, dimension_name, dimension_type, dimension_size
     &  , dimension_array_int, dimension_array_fl
     &  , 'mom', 3, 'int', 3
     &  , .false., n_phase_term, .true., .true.
     &  , id_mom, i_dummy, dummy
     &  , nd_phase_term
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!     Single scattering Properties:
!
!     Optical Depth:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'tau', 3, 'float', 5
     &  , .true.
     &  , iv_tau
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!     Albedo of Single Scattering:
!
      CALL find_var_cdl(ierr
     &  , file_name
     &  , n_var, var_name, var_type
     &  , 'omega', 5, 'float', 5
     &  , .true.
     &  , iv_omega
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
!     The forward scattering may be present.
      IF (l_forward) THEN
         CALL find_var_cdl(ierr
     &     , file_name
     &     , n_var, var_name, var_type
     &     , 'frwd', 4, 'float', 5
     &     , .true.
     &     , iv_frwd
     &     , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &     )
      ENDIF
!
!     The forward scattering for the solar beam may be present.
      IF (l_forward_solar) THEN
         CALL find_var_cdl(ierr
     &     , file_name
     &     , n_var, var_name, var_type
     &     , 'sfrwd', 5, 'float', 5
     &     , .true.
     &     , iv_sfrwd
     &     , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &     )
      ENDIF
!
!
!
!     The data in the CDL file will be stored with the last
!     `index' changing most rapidly (C-convention), but there
!     will be a common stride between values at the same
!     point. The phase function will have different strides
!     because it loops over moments.
!
      stride_lat_to=calc_cdl_stride(n_dimension_var(iv_tau)
     &   , list_dimension_var(1, iv_tau), id_lat
     &   , dimension_size, nd_cdl_dimen)
      stride_long_to=calc_cdl_stride(n_dimension_var(iv_tau)
     &   , list_dimension_var(1, iv_tau), id_long
     &   , dimension_size, nd_cdl_dimen)
      stride_vlev_to=calc_cdl_stride(n_dimension_var(iv_tau)
     &   , list_dimension_var(1, iv_tau), id_vlev
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
      stride_mom_phf=calc_cdl_stride(n_dimension_var(iv_phf)
     &   , list_dimension_var(1, iv_phf), id_mom
     &   , dimension_size, nd_cdl_dimen)
!
!     Assign the data.
      DO i_vlev=1, n_layer
        DO i_lat=1, n_latitude
          DO i_long=1, n_longitude
!
            l=i_long+(i_lat-1)*n_longitude
            offset_to=stride_vlev_to*(i_vlev-1)
     &               +stride_lat_to*(i_lat-1)
     &               +stride_long_to*(i_long-1)
!
            tau(l, i_vlev)
     &        =data_fl(1+offset_to, iv_tau)
            omega(l, i_vlev)
     &        =data_fl(1+offset_to, iv_omega)
!
            IF (l_forward) THEN
              forward_scatter(l, i_vlev)
     &          =data_fl(1+offset_to, iv_frwd)
            ENDIF
!
            IF (l_forward_solar) THEN
              forward_solar(l, i_vlev)
     &          =data_fl(1+offset_to, iv_sfrwd)
            ENDIF
!
            DO i_mom=1, n_phase_term
              offset_phf=stride_vlev_phf*(i_vlev-1)
     &                  +stride_lat_phf*(i_lat-1)
     &                  +stride_long_phf*(i_long-1)
     &                  +stride_mom_phf*(i_mom-1)
              phase_function(l, i_vlev, i_mom)
     &          =data_fl(1+offset_phf, iv_phf)
            ENDDO
!
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END
