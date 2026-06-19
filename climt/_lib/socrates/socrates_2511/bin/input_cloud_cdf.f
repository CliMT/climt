! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to control the reading of input for clouds.
!
! Purpose:
!   This subroutine provides overall control of the input of cloudy
!   fields.
!
! Method:
!   To be done.
!
!- ---------------------------------------------------------------------
      SUBROUTINE input_cloud_cdf(ierr
     &  , n_band
     &  , l_drop_type, i_drop_parametrization
     &  , n_drop_phf_term, drop_parameter_list
     &  , drop_min_dim, drop_max_dim
     &  , l_ice_type, i_ice_parametrization
     &  , n_ice_phf_term, ice_parameter_list
     &  , ice_min_dim, ice_max_dim
     &  , base_name, length_name
     &  , l_drop, l_ice
     &  , i_cloud, l_cloud, i_cloud_representation
     &  , n_condensed, type_condensed
     &  , i_condensed_param, condensed_n_phf, condensed_param_list
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , l_vert_coord, name_vert_coord, n_layer, p
     &  , condensed_mix_ratio, condensed_dim_char
     &  , w_cloud, n_cloud_type, i_cloud_type, frac_cloud
     &  , n_opt_level_drop_prsc, n_phase_term_drop_prsc
     &  , drop_pressure_prsc, drop_absorption_prsc
     &  , drop_scattering_prsc, drop_phase_fnc_prsc
     &  , n_opt_level_ice_prsc, n_phase_term_ice_prsc
     &  , ice_pressure_prsc, ice_absorption_prsc
     &  , ice_scattering_prsc, ice_phase_fnc_prsc
     &  , nd_profile, nd_latitude, nd_longitude, nd_layer
     &  , nd_band, nd_phase_term
     &  , nd_drop_type, nd_ice_type, nd_cloud_parameter
     &  , nd_cloud_component, nd_cloud_type
     &  , nd_profile_cloud_prsc, nd_opt_level_cloud_prsc
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE rad_pcf
      USE realtype_rd
      USE def_std_io_icf
      USE gas_list_pcf
      USE input_head_pcf
!
!
      IMPLICIT NONE
!
!
!     INCLUDE HEADER FILES.
!
!
!
!     Dummy Arguments:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
!
!     Dimensions of arrays:
!
!     Physical dimensions:
!
      INTEGER, Intent(IN) ::
     &    nd_latitude
!           Allowed size for latitudes
     &  , nd_longitude
!           Allowed size for longitudes
     &  , nd_profile
!           Maximum number of atmospheric profiles
     &  , nd_layer
!           Maximum number of atmospheric layers
!
!     Spectral dimensions:
      INTEGER, Intent(IN) ::
     &    nd_band
!           Maximum number of spectral bands
     &  , nd_phase_term
!           Maximum number of terms in the phase function
     &  , nd_drop_type
!           Maximum number of types of droplets
     &  , nd_ice_type
!           Maximum number of types of ice crystals
     &  , nd_cloud_parameter
!           Maximum number of parameters for clouds
     &  , nd_cloud_component
!           Maximum number of cloud components
     &  , nd_cloud_type
!           Maximum number of cloud types
!
!     Dimensions for fields of prescribed optical properties:
      INTEGER, Intent(IN) ::
     &    nd_profile_cloud_prsc
!           Size allocated for profiles of prescribed
!           cloudy optical properties
     &  , nd_opt_level_cloud_prsc
!           Size allocated for levels of prescribed
!           cloudy optical properties
!     Dimensions for CDL fields:
      INTEGER, Intent(IN) ::
     &    nd_cdl_dimen
!           Maximum number of CDL dimensions
     &  , nd_cdl_dimen_size
!           Maximum number of components in each dimension
     &  , nd_cdl_data
!           Maximum number of elements of CDL data
     &  , nd_cdl_var
!           Maximum number of elements of CDL arrays
!
!     Characteristics of file-names
      CHARACTER !, Intent(IN)
     &    base_name*(*)
!           Basename of input files
      INTEGER, Intent(IN) ::
     &    length_name
!           Length of input filename
!
!     Vertical coordinates:
      LOGICAL, Intent(INOUT) ::
     &    l_vert_coord
!           Flag for vertical coordinate
      CHARACTER
     &    name_vert_coord*(*)
!           Name of vertical coordinate
!
!     General atmospheric properties
      INTEGER, Intent(INOUT) ::
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_layer
!           Number of layers
!
!     Spectral variables
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of spectral bands
      LOGICAL, Intent(IN) ::
     &    l_drop_type(nd_drop_type)
!           Flags for the presence of data for types of droplets
     &  , l_ice_type(nd_ice_type)
!           Flags for the presence of data for types of ice crystals
      INTEGER, Intent(IN) ::
     &    i_drop_parametrization(nd_drop_type)
!           Type of parametrization for each type of droplet
     &  , n_drop_phf_term(nd_drop_type)
!           Number of moments in the phase function for droplets
     &  , i_ice_parametrization(nd_ice_type)
!           Type of parametrization for each type of ice crystal
     &  , n_ice_phf_term(nd_ice_type)
!           Number of moments in the phase function for ice crystals
      REAL  (RealK), Intent(IN) ::
     &    drop_parameter_list(nd_cloud_parameter
     &      , nd_band, nd_drop_type)
!           Parameters for droplets
     &  , ice_parameter_list(nd_cloud_parameter
     &      , nd_band, nd_ice_type)
!           Parameters for ice crystals
     &  , drop_min_dim(nd_drop_type), drop_max_dim(nd_drop_type)
!           Range of allowed dimensions for spectral file droplet fit
     &  , ice_min_dim(nd_ice_type), ice_max_dim(nd_ice_type)
!           Range of allowed dimensions for spectral file ice crystal fit
!
      LOGICAL, Intent(IN) ::
     &    l_drop
!           Flag for inclusion of droplets
     &  , l_ice
!           Flag for inclusion of ice crystals
      LOGICAL, Intent(OUT) ::
     &    l_cloud
!           Flag for inclusion of clouds
!
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
     &  , p(nd_profile, nd_layer)
!           Pressures in layers
!
!     Cloudy fields
      INTEGER, Intent(OUT) ::
     &    i_cloud
!           Cloud scheme
     &  , n_condensed
!           Number of condensed phases in clouds
     &  , n_cloud_type
!           Number of types of clouds
     &  , i_cloud_type(nd_cloud_component)
!           Cloud type for each component
     &  , type_condensed(nd_cloud_component)
!           Types of condensed phases in clouds
     &  , i_condensed_param(nd_cloud_component)
!           Parametrization schemes for components
     &  , i_cloud_representation
!           Representation of the mixing of components in clouds
     &  , condensed_n_phf(nd_cloud_component)
!           Number of terms in the phase function for condensed
!           components
      REAL  (RealK), Intent(OUT) ::
     &    w_cloud(nd_profile, nd_layer)
!          Amounts of cloud
     &  , frac_cloud(nd_profile, nd_layer, nd_cloud_type)
!          Fractions of types of clouds
     &  , condensed_mix_ratio(nd_profile, nd_layer
     &     , nd_cloud_component)
!          Mixing ratios of condensed components
     &  , condensed_dim_char(nd_profile, nd_layer
     &     , nd_cloud_component)
!          Characteristic dimensions of condensed components
     &  , condensed_param_list(nd_cloud_parameter
     &     , nd_cloud_component, nd_band)
!          Coefficients in parametrizations of condensed components
!
!     Fields for prescribed optical properties of droplets
!
      INTEGER, Intent(OUT) ::
     &    n_opt_level_drop_prsc
!           Number of levels of prescribed
!           optical properties of droplets
     &  , n_phase_term_drop_prsc
!           Number of terms in the prescribed phase function for
!           droplets
      REAL  (RealK), Intent(OUT) ::
     &    drop_pressure_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc)
!           Pressures at which optical properties of
!           droplets are prescribed
     &  , drop_absorption_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc, nd_band)
!           Prescribed absorption by droplets
     &  , drop_scattering_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc, nd_band)
!           Prescribed scattering by droplets
     &  , drop_phase_fnc_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc, nd_phase_term, nd_band)
!           Prescribed phase function of droplets
!
!     Fields for prescribed optical properties of ice crystals
      INTEGER, Intent(OUT) ::
     &    n_opt_level_ice_prsc
!           Number of levels of prescribed
!           optical properties of ice crystals
     &  , n_phase_term_ice_prsc
!           Number of terms in the prescribed phase function for 
!           ice crystals
      REAL  (RealK), Intent(OUT) ::
     &    ice_pressure_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc)
!           Pressures at which optical properties of
!           ice crystals are prescribed
     &  , ice_absorption_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc, nd_band)
!           Prescribed absorption by ice crystals
     &  , ice_scattering_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc, nd_band)
!           Prescribed scattering by ice crystals
     &  , ice_phase_fnc_prsc(nd_profile_cloud_prsc
     &      , nd_opt_level_cloud_prsc, nd_phase_term, nd_band)
!           Prescribed phase functions of ice crystals
!
!
!     Local Variables
!
      CHARACTER
     &    file_name*80
!           Name of input file
     &  , text*80
!           Text for I/O
     &  , suffix_mr*12
!           Suffix for mixing ratios
     &  , suffix_dim_char*12
!           Suffix for characteristic dimensions
     &  , suffix_cl(nd_cloud_component)*12
!           Suffix for clouds
!
!     Cloudy fields
      INTEGER
     &    i_phase
!           Phase of component in cloud
     &  , i_parametrization_type
!           Type of condensed component
     &  , n_cloud_parameter
!           Number of parameters for the parametrization scheme selected
     &  , i
!           Loop variable
     &  , j
!           Loop variable
     &  , k
!           Loop variable
     &  , l
!           Loop variable
!
      LOGICAL
     &    l_exist_opt
!           Flag for existence of file of optical properties
     &  , l_drop_prsc_set
!           Flag for prescribed optical properties of droplets
     &  , l_ice_prsc_set
!           Flag for prescribed optical properties of ice crystals
     &  , l_vert_assignable
!           Flag to permit the assigning of the values of vertical
!           coordinates
!
!
!     Functions called:
      INTEGER
     &    set_n_cloud_parameter
!           Function to find the number of cloud parameters
      EXTERNAL
     &    set_n_cloud_parameter
!
!     Subroutines called:
      EXTERNAL
     &     assign_input_vert_cdf, assign_input_opt_cdf
!
!
!
!
!
      WRITE(iu_stdout, '(/a)') 'Enter cloud scheme.'
      READ(iu_stdin, *) i_cloud
!
!     Set the cloud flag from consideration of the value of I_CLOUD
!     and determine all the components present.
!
      IF ( (i_cloud == IP_cloud_mix_max).OR.
     &     (i_cloud == IP_cloud_mix_random).OR.
     &     (i_cloud == IP_cloud_triple).OR.
     &     (i_cloud == IP_cloud_part_corr).OR.
     &     (i_cloud == IP_cloud_part_corr_cnv).OR.
     &     (i_cloud == IP_cloud_column_max).OR.
     &     (i_cloud == IP_cloud_mcica) ) THEN
!
!       Clouds are included.
        l_cloud=.true.
!
        WRITE(iu_stdout, '(a)') 'Enter representation of clouds.'
        READ(iu_stdin, *) i_cloud_representation
!
        IF ( l_ice.AND.l_drop) THEN
!
!         Both phases are present.
          IF ( (i_cloud_representation == IP_cloud_homogen).OR.
     &         (i_cloud_representation == IP_cloud_ice_water) ) THEN
            n_condensed=2
            type_condensed(1)=IP_clcmp_st_water
            type_condensed(2)=IP_clcmp_st_ice
          ELSE IF ( (i_cloud_representation
     &             == IP_cloud_conv_strat).OR.
     &            (i_cloud_representation == IP_cloud_csiw) ) THEN
            n_condensed=4
            type_condensed(1)=IP_clcmp_st_water
            type_condensed(2)=IP_clcmp_st_ice
            type_condensed(3)=IP_clcmp_cnv_water
            type_condensed(4)=IP_clcmp_cnv_ice
          ENDIF
!
        ELSE IF ( l_ice.AND.(.NOT.l_drop) ) THEN
!
!         Only ice is present.
          IF ( (i_cloud_representation == IP_cloud_homogen).OR.
     &         (i_cloud_representation == IP_cloud_ice_water) ) THEN
            n_condensed=1
            type_condensed(1)=IP_clcmp_st_ice
          ELSE IF ( (i_cloud_representation
     &             == IP_cloud_conv_strat).OR.
     &            (i_cloud_representation == IP_cloud_csiw) ) THEN
            n_condensed=2
            type_condensed(1)=IP_clcmp_st_ice
            type_condensed(2)=IP_clcmp_cnv_ice
          ENDIF
!
        ELSE IF ( (.NOT.l_ice).AND.l_drop ) THEN
!
!         Only droplets are present.
          IF ( (i_cloud_representation == IP_cloud_homogen).OR.
     &         (i_cloud_representation == IP_cloud_ice_water) ) THEN
            n_condensed=1
            type_condensed(1)=IP_clcmp_st_water
          ELSE IF ( (i_cloud_representation
     &             == IP_cloud_conv_strat).OR.
     &            (i_cloud_representation == IP_cloud_csiw) ) THEN
            n_condensed=2
            type_condensed(1)=IP_clcmp_st_water
            type_condensed(2)=IP_clcmp_cnv_water
          ENDIF
!
        ELSE
!
          WRITE(iu_err, '(/a)') 
     &      '*** Error: Clouds are specified, but no cloudy '
     &      //'proceses are included.'
          ierr=i_err_fatal
          RETURN
!
        ENDIF
!
!       Read in any prescribed optical properties:
!
        l_drop_prsc_set=.false.
        IF (l_drop) THEN
!         Determine whether a file of prescribed optical properties
!         exists.
          file_name(1: length_name+1+len_file_suffix)
     &      =base_name(1: length_name)
     &      //'.'//phys_suffix(IP_optical_water)
          INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix))
     &      , exist=l_exist_opt)
!
          IF (l_exist_opt) THEN
            CALL assign_input_opt_cdf(ierr
     &        , file_name(1: length_name+len_file_suffix)
     &        , n_band
     &        , n_latitude, latitude, n_longitude, longitude
     &        , n_profile
     &        , n_opt_level_drop_prsc, drop_pressure_prsc
     &        , n_phase_term_drop_prsc, drop_absorption_prsc
     &        , drop_scattering_prsc, drop_phase_fnc_prsc
     &        , nd_profile_cloud_prsc, nd_latitude, nd_longitude
     &        , nd_band, nd_phase_term, nd_opt_level_cloud_prsc
     &        , nd_cdl_dimen, nd_cdl_dimen_size
     &        , nd_cdl_data*n_band*nd_phase_term, nd_cdl_var
     &        )
          ENDIF
          l_drop_prsc_set=.true.
        ENDIF
!
        l_ice_prsc_set=.false.
        IF (l_ice) THEN
!         Determine whether a file of prescribed optical properties
!         exists.
          file_name(1: length_name+1+len_file_suffix)
     &      =base_name(1: length_name)
     &      //'.'//phys_suffix(IP_optical_ice)
          INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix))
     &      , exist=l_exist_opt)
!
          IF (l_exist_opt) THEN
            CALL assign_input_opt_cdf(ierr
     &        , file_name(1: length_name+len_file_suffix)
     &        , n_band
     &        , n_latitude, latitude, n_longitude, longitude
     &        , n_profile
     &        , n_opt_level_ice_prsc, ice_pressure_prsc
     &        , n_phase_term_ice_prsc, ice_absorption_prsc
     &        , ice_scattering_prsc, ice_phase_fnc_prsc
     &        , nd_profile_cloud_prsc, nd_latitude, nd_longitude
     &        , nd_band, nd_phase_term, nd_opt_level_cloud_prsc
     &        , nd_cdl_dimen, nd_cdl_dimen_size
     &        , nd_cdl_data*n_band*nd_phase_term, nd_cdl_var
     &        )
          ENDIF
          l_ice_prsc_set=.true.
        ENDIF
!
!
!       Read in the data for all components.
        DO i=1, n_condensed
!
          IF (type_condensed(i) == IP_clcmp_st_water) THEN
            text='enter type of (stratiform) droplets'
            i_phase=IP_phase_water
            suffix_mr(1:len_file_suffix)
     &        =phys_suffix(IP_lwm)(1:len_file_suffix)
            suffix_dim_char(1:len_file_suffix)
     &        =phys_suffix(IP_re)(1:len_file_suffix)
          ELSE IF (type_condensed(i) == IP_clcmp_st_ice) THEN
            text='enter type of (stratiform) ice crystals'
            i_phase=IP_phase_ice
            suffix_mr(1:len_file_suffix)
     &        =phys_suffix(IP_iwm)(1:len_file_suffix)
            suffix_dim_char(1:len_file_suffix)
     &        =phys_suffix(IP_ire)(1:len_file_suffix)
          ELSE IF (type_condensed(i) == IP_clcmp_cnv_water) THEN
            text='enter type of convective droplets'
            i_phase=IP_phase_water
            suffix_mr(1:len_file_suffix)
     &        =phys_suffix(IP_lwmcv)(1:len_file_suffix)
            suffix_dim_char(1:len_file_suffix)
     &        =phys_suffix(IP_recv)(1:len_file_suffix)
          ELSE IF (type_condensed(i) == IP_clcmp_cnv_ice) THEN
            text='enter type of convective ice crystals'
            i_phase=IP_phase_ice
            suffix_mr(1:len_file_suffix)
     &        =phys_suffix(IP_iwmcv)(1:len_file_suffix)
            suffix_dim_char(1:len_file_suffix)
     &        =phys_suffix(IP_irecv)(1:len_file_suffix)
          ENDIF
!
          WRITE(iu_stdout, '(a)') text
          READ(iu_stdin, *) i_parametrization_type
!
!
          IF (i_parametrization_type <= 0) THEN
!
            IF (i_phase == IP_phase_water) THEN
              IF (.NOT.l_drop_prsc_set) THEN
                WRITE(iu_err, '(/a, i3)')
     &            '*** Error: No prescribed optical data for '
     &            //'droplets have been supplied.'
                ierr=i_err_fatal
                RETURN
              ENDIF
              i_condensed_param(i)=IP_drop_unparametrized
            ELSE IF (i_phase == IP_phase_ice) THEN
              IF (.NOT.l_ice_prsc_set) THEN
                WRITE(iu_err, '(/a, i3)')
     &            '*** Error: No prescribed optical data for '
     &            //'ice crystals have been supplied.'
                ierr=i_err_fatal
                RETURN
              ENDIF
              i_condensed_param(i)=IP_ice_unparametrized
            ENDIF
!
          ELSE
!
!           This component is parametrized, so mixing ratios and sizes
!           are required.
!
            IF (i_phase == IP_phase_water) THEN
!
              IF (.NOT.l_drop_type(i_parametrization_type)) THEN
                WRITE(iu_err, '(/a, i3)')
     &            '*** Error: The spectral file contains no data '
     &            //'for the type of droplet', i_parametrization_type
                ierr=i_err_fatal
                RETURN
              ENDIF
!
              i_condensed_param(i)
     &          =i_drop_parametrization(i_parametrization_type)
              condensed_n_phf(i)
     &          =n_drop_phf_term(i_parametrization_type)
!
              n_cloud_parameter=set_n_cloud_parameter(
     &           i_condensed_param(i)
     &           , type_condensed(i), condensed_n_phf(i))
              IF (n_cloud_parameter > nd_cloud_parameter) THEN
                WRITE(iu_err, '(/a, /a)')
     &            '*** Error: The number of cloud parameters '
     &            //'exceeds the dimension. ' 
     &            , 'recompile with a larger value '
     &            //'of npd_cloud_parameter.'
                ierr=i_err_fatal
                RETURN
              ENDIF
!
              DO j=1, n_band
                DO k=1, n_cloud_parameter
                  condensed_param_list(k, i, j)
     &              =drop_parameter_list(k, j
     &              , i_parametrization_type)
                ENDDO
              ENDDO
!
            ELSE IF (i_phase == IP_phase_ice) THEN
!
              IF (.NOT.l_ice_type(i_parametrization_type)) THEN
                WRITE(iu_err, '(/a, i3)')
     &            '*** Error: The spectral file contains no data '
     &            //'for the type of ice crystal'
     &            , i_parametrization_type
                STOP
              ENDIF
!
              i_condensed_param(i)
     &          =i_ice_parametrization(i_parametrization_type)
              condensed_n_phf(i)
     &          =n_ice_phf_term(i_parametrization_type)
!
              n_cloud_parameter=set_n_cloud_parameter(
     &           i_condensed_param(i)
     &           , type_condensed(i), condensed_n_phf(i))
              IF (n_cloud_parameter > nd_cloud_parameter) THEN
                WRITE(iu_err, '(/a, /a)')
     &            '*** Error: The number of cloud parameters '
     &            //'exceeds the dimension. ' 
     &            , 'recompile with a larger value '
     &            //'of npd_cloud_parameter.'
                ierr=i_err_fatal
                RETURN
              ENDIF
!
              DO j=1, n_band
                DO k=1, n_cloud_parameter
                  condensed_param_list(k, i, j)
     &              =ice_parameter_list
     &              (k, j, i_parametrization_type)
                ENDDO
              ENDDO
!
            ENDIF
!
!
            file_name(1: length_name+1+len_file_suffix)
     &        =base_name(1: length_name)//'.'//suffix_mr
            l_vert_assignable=.NOT.l_vert_coord
            CALL assign_input_vert_cdf(ierr
     &        , file_name(1: length_name+1+len_file_suffix)
     &        , 'condensed mixing ratios'
     &        , l_vert_coord, name_vert_coord
     &        , .true., n_layer, l_vert_assignable
     &        , n_latitude, latitude, n_longitude, longitude
     &        , 1
     &        , n_profile, n_layer
     &        , p, condensed_mix_ratio(1, 1, i)
     &        , nd_profile, nd_latitude, nd_longitude
     &        , 1, nd_layer
     &        , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &        )
            IF (ierr /= i_normal) RETURN
!
            file_name(1: length_name+1+len_file_suffix)
     &        =base_name(1: length_name)
     &        //'.'//suffix_dim_char
            l_vert_assignable=.NOT.l_vert_coord
            CALL assign_input_vert_cdf(ierr
     &        , file_name(1: length_name+1+len_file_suffix)
     &        , 'condensed CHARACTERistic dimensions'
     &        , l_vert_coord, name_vert_coord
     &        , .true., n_layer, l_vert_assignable
     &        , n_latitude, latitude, n_longitude, longitude
     &        , 1
     &        , n_profile, n_layer
     &        , p, condensed_dim_char(1, 1, i)
     &        , nd_profile, nd_latitude, nd_longitude, 1, nd_layer
     &        , nd_cdl_dimen, nd_cdl_dimen_size
     &        , nd_cdl_data, nd_cdl_var
     &        )
            IF (ierr /= i_normal) RETURN

            IF (i_phase == IP_phase_water) THEN
              DO k = 1, n_layer
                DO l = 1, n_profile
                  condensed_dim_char(l, k, i) = MIN( MAX( 
     &              condensed_dim_char(l, k, i),
     &              drop_min_dim(i_parametrization_type) ),
     &              drop_max_dim(i_parametrization_type) )
                END DO
              END DO
            ELSE IF (i_phase == IP_phase_ice) THEN
              DO k = 1, n_layer
                DO l = 1, n_profile
                  condensed_dim_char(l, k, i) = MIN( MAX( 
     &              condensed_dim_char(l, k, i),
     &              ice_min_dim(i_parametrization_type) ),
     &              ice_max_dim(i_parametrization_type) )
                END DO
              END DO
            END IF

          ENDIF
!
        ENDDO
!
!
!
!       Setting of amounts of cloud.
!
!       Set the number of types of cloud and their properties 
!       from the representation.
!
        IF (i_cloud_representation == IP_cloud_homogen) THEN
          n_cloud_type=1
          DO i=1, n_condensed
            i_cloud_type(i) = ip_cloud_type_homogen
          END DO
          suffix_cl(1)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction)(1:len_file_suffix)
        ELSE IF (i_cloud_representation == IP_cloud_ice_water) THEN
          n_cloud_type=2
          DO i = 1, n_condensed
            SELECT CASE (type_condensed(i))
            CASE (ip_clcmp_st_water)
              i_cloud_type(i) = ip_cloud_type_water
            CASE (ip_clcmp_st_ice)
              i_cloud_type(i) = ip_cloud_type_ice
            CASE DEFAULT
              WRITE(iu_err, '(a)')
     &          '*** Error in input_cloud_cdf: type_condensed '
     &          //'incorrectly set for ip_cloud_ice_water'
              ierr=i_err_fatal
              RETURN
            END SELECT
          END DO
          suffix_cl(1)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_w)(1:len_file_suffix)
          suffix_cl(2)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_i)(1:len_file_suffix)
        ELSE IF (i_cloud_representation == IP_cloud_conv_strat) THEN
          n_cloud_type=2
          DO i = 1, n_condensed
            SELECT CASE (type_condensed(i))
            CASE (ip_clcmp_st_water)
              i_cloud_type(i) = ip_cloud_type_strat
            CASE (ip_clcmp_st_ice)
              i_cloud_type(i) = ip_cloud_type_strat
            CASE (ip_clcmp_cnv_water)
              i_cloud_type(i) = ip_cloud_type_conv
            CASE (ip_clcmp_cnv_ice)
              i_cloud_type(i) = ip_cloud_type_conv
            CASE DEFAULT
              WRITE(iu_err, '(a)')
     &          '*** Error in input_cloud_cdf: type_condensed '
     &          //'incorrectly set for ip_cloud_conv_strat'
              ierr=i_err_fatal
              RETURN
            END SELECT
          END DO
          suffix_cl(1)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction)(1:len_file_suffix)
          suffix_cl(2)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_cv)(1:len_file_suffix)
        ELSE IF (i_cloud_representation == IP_cloud_csiw) THEN
          n_cloud_type=4
          DO i = 1, n_condensed
            SELECT CASE (type_condensed(i))
            CASE (ip_clcmp_st_water)
              i_cloud_type(i) = ip_cloud_type_sw
            CASE (ip_clcmp_st_ice)
              i_cloud_type(i) = ip_cloud_type_si
            CASE (ip_clcmp_cnv_water)
              i_cloud_type(i) = ip_cloud_type_cw
            CASE (ip_clcmp_cnv_ice)
              i_cloud_type(i) = ip_cloud_type_ci
            CASE DEFAULT
              WRITE(iu_err, '(a)')
     &          '*** Error in input_cloud_cdf: type_condensed '
     &          //'incorrectly set for ip_cloud_csiw'
              ierr=i_err_fatal
              RETURN
            END SELECT
          END DO
          suffix_cl(1)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_w)(1:len_file_suffix)
          suffix_cl(2)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_i)(1:len_file_suffix)
          suffix_cl(3)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_w_cv)(1:len_file_suffix)
          suffix_cl(4)(1:len_file_suffix)
     &      =phys_suffix(IP_cloud_fraction_i_cv)(1:len_file_suffix)
        ENDIF
!
        DO i=1, n_cloud_type
          file_name(1: length_name+1+len_file_suffix)
     &        =base_name(1: length_name)
     &      //'.'//suffix_cl(i)
          l_vert_assignable=.NOT.l_vert_coord
          CALL assign_input_vert_cdf(ierr
     &      , file_name(1: length_name+1+len_file_suffix)
     &      , 'cloud fractions'
     &      , l_vert_coord, name_vert_coord
     &      , .true., n_layer, l_vert_assignable
     &      , n_latitude, latitude, n_longitude, longitude
     &      , 1
     &      , n_profile, n_layer
     &      , p, frac_cloud(1, 1, i)
     &      , nd_profile, nd_latitude, nd_longitude, 1, nd_layer
     &      , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &      )
          IF (ierr /= i_normal) RETURN
        ENDDO
!
!       Normalize the cloud fractions.
        DO l=1, n_profile
          DO i=1, n_layer
            w_cloud(l, i)=0.0_RealK
            DO j=1, n_cloud_type
              w_cloud(l, i)=w_cloud(l, i)+frac_cloud(l, i, j)
            ENDDO
            IF (w_cloud(l, i) > 0.0_RealK) THEN
              DO j=1, n_cloud_type
                frac_cloud(l, i, j)=frac_cloud(l, i, j)/w_cloud(l, i)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!
      ELSE
!
!       Clouds are not included.
        l_cloud=.false.
!
      ENDIF
!
!
!
      RETURN
      END
