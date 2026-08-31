! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to control the reading of input for aerosols.
!
! Purpose:
!   This subroutine provides overall control of the input of aerosol
!   fields.
!
! Method:
!   To be done.
!
!- ---------------------------------------------------------------------
      SUBROUTINE input_aerosol_cdl(ierr
     &  , base_name, length_name
     &  , n_band, index_water
     &  , l_vert_coord, name_vert_coord, l_q_unread, l_opt_overlay
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_layer, p
     &  , n_aerosol, type_aerosol, aerosol_mix_ratio
     &  , i_aerosol_parametrization, water_mix_ratio
     &  , n_opt_level_aerosol_prsc, n_phase_term_aerosol_prsc
     &  , aerosol_pressure_prsc, aerosol_absorption_prsc
     &  , aerosol_scattering_prsc, aerosol_phase_fnc_prsc
     &  , nd_profile, nd_latitude, nd_longitude, nd_layer
     &  , nd_band, nd_aerosol_species, nd_phase_term
     &  , nd_profile_prsc, nd_level_prsc
     &  , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE rad_pcf
      USE gas_list_pcf
      USE input_head_pcf
!
!
      IMPLICIT NONE
!
!
!
!
!     Dummy arguments:
      INTEGER
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
     &  , nd_level_prsc
!           Allowed size for levels of optical properties
     &  , nd_profile_prsc
!           Allowed size for prescribed atmospheric profiles
!
!     Spectral dimensions:
      INTEGER, Intent(IN) ::
     &    nd_band
!           Size allocated for spectral bands
     &  , nd_aerosol_species
!           Size allocated for species of aerosols
     &  , nd_phase_term
!           Size allocated for terms in the phase function
!
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
!
!     Characteristics of file-names
      CHARACTER !, Intent(IN)
     &    base_name*(*)
!           Basename of input files
      INTEGER, Intent(IN) ::
     &    length_name
!           Length of input filename
      INTEGER, Intent(IN) ::
     &    n_band
!           Number of levels
     &  , index_water
!           Index of water vapour
      LOGICAL, Intent(INOUT) ::
     &    l_q_unread
!           Flag for profile of specific humidity
     &  , l_vert_coord
!           Flag for vertical coordinate
      LOGICAL, Intent(IN) ::
     &    l_opt_overlay
!           Flag to overlay optical properties
      CHARACTER !, Intent(INOUT)
     &    name_vert_coord*(*)
!           Name of the vertical coordinate
!
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
!     Horizontal structure
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitude
     &  , longitude(nd_longitude)
!           Longitude
!     Vertical structure
      REAL  (RealK), Intent(INOUT) ::
     &    p(nd_profile, nd_layer)
!           Pressures
!
!     Aerosol Fields
      INTEGER, Intent(IN) ::
     &    n_aerosol
!           Number of aerosols
     &  , type_aerosol(nd_aerosol_species)
!           Types of aerosols required
      INTEGER, Intent(INOUT) ::
     &    i_aerosol_parametrization(nd_aerosol_species)
!           Parametrizations adopted for aerosols
!
      INTEGER, Intent(OUT) ::
     &    n_opt_level_aerosol_prsc(nd_aerosol_species)
!           Number of levels of optical properties for each aerosol
     &  , n_phase_term_aerosol_prsc(nd_aerosol_species)
!           Number of terms in the phase function for each species
!
      REAL  (RealK), Intent(OUT) ::
     &    aerosol_mix_ratio(nd_profile, nd_layer, nd_aerosol_species)
!           Mixing ratios of aserosols
     &  , water_mix_ratio(nd_profile, nd_layer)
!           Mixing ratios of water vapour
!
!     Prescribed optical properties
      REAL  (RealK), Intent(OUT) ::
     &    aerosol_pressure_prsc(nd_profile_prsc, nd_level_prsc
     &      , nd_aerosol_species)
!           Pressures of prescribed data
     &  , aerosol_absorption_prsc(nd_profile_prsc, nd_level_prsc
     &      , nd_aerosol_species, nd_band)
!           Prescribed absorption
     &  , aerosol_scattering_prsc(nd_profile_prsc, nd_level_prsc
     &      , nd_aerosol_species, nd_band)
!           Prescribed scattering
     &  , aerosol_phase_fnc_prsc(nd_profile_prsc, nd_level_prsc
     &      , nd_phase_term, nd_aerosol_species, nd_band)
!           Prescribed Phase Functions
!
!
!
!     Local variables:
      CHARACTER
     &    file_name*80
!           Name of current file
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
     &  , il
!           Loop variable
     &  , ib
!           Loop variable
     &  , im
!           Loop variable
      LOGICAL
     &    l_exist_opt
!           Flag for the existence of a file of optical properties
     &  , l_exist
!           Flag for the existence of a file of mass mixing ratios
     &  , l_vert_assignable
!           Flag to permit the assignment of values of the vertical
!           coordinate
      REAL  (RealK) ::
     &    absorption_local(nd_profile_prsc, nd_level_prsc, nd_band)
!           Local array of absorption
     &  , scattering_local(nd_profile_prsc, nd_level_prsc, nd_band)
!           Local array of absorption
     &  , phase_function_local(nd_profile_prsc, nd_level_prsc
     &      , nd_phase_term, nd_band)
!           Local array of absorption
!
!
!
!
!
!
      DO i=1, n_aerosol
!
!       If prescribed optical properties are not supplied mixing
!       ratios must be specified.
!
        IF (l_opt_overlay) THEN
          n_phase_term_aerosol_prsc(i)=0
!         Determine whether a file of prescribed optical properties
!         exists.
          file_name(1: length_name+1+len_file_suffix)
     &      =base_name(1: length_name)
     &      //'.'//aerosol_opt_suffix(type_aerosol(i))
          INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix))
     &      ,exist=l_exist_opt)
!
          IF (l_exist_opt) THEN
            i_aerosol_parametrization(i) = IP_aerosol_unparametrized
            CALL assign_input_opt_cdl(ierr
     &        , file_name(1: length_name+len_file_suffix) 
     &        , n_band
     &        , n_latitude, latitude, n_longitude, longitude, n_profile
     &        , n_opt_level_aerosol_prsc(i)
     &        , aerosol_pressure_prsc(1, 1, i)
     &        , n_phase_term_aerosol_prsc(i)
     &        , absorption_local, scattering_local, phase_function_local
     &        , nd_profile_prsc, nd_latitude, nd_longitude
     &        , nd_band, nd_phase_term, nd_level_prsc
     &        , nd_cdl_dimen, nd_cdl_dimen_size
     &        , nd_profile_prsc*nd_level_prsc*nd_phase_term*n_band
     &        , nd_cdl_var
     &        )
!
!           Transpose the local optical properties into the full arrays.
            DO ib=1, n_band
              DO il=1, n_opt_level_aerosol_prsc(i)
                DO l=1, n_profile
                  aerosol_absorption_prsc(l, il, i, ib)
     &              =absorption_local(l, il, ib)
                  aerosol_scattering_prsc(l, il, i, ib)
     &              =scattering_local(l, il, ib)
                ENDDO
                DO im=1, n_phase_term_aerosol_prsc(i)
                  DO l=1, n_profile
                    aerosol_phase_fnc_prsc(l, il, im, i, ib)
     &                =phase_function_local(l, il, im, ib)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
!
          ENDIF
!          
        ENDIF
!
        IF (.NOT.(l_opt_overlay.AND.l_exist_opt)) THEN
          file_name(1: length_name+1+len_file_suffix)
     &      =base_name(1: length_name)
     &      //'.'//aerosol_suffix(type_aerosol(i))
          INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix)) 
     &        ,exist=l_exist)
          IF (l_exist) THEN
            l_vert_assignable=.NOT.l_vert_coord
            CALL assign_input_vert_cdl(ierr
     &        , file_name(1: length_name+len_file_suffix)
     &        , 'aerosol mixing ratios'
     &        , l_vert_coord, name_vert_coord
     &        , .true., n_layer, l_vert_assignable
     &        , n_latitude, latitude, n_longitude, longitude
     &        , 1
     &        , n_profile, n_layer
     &        , p, aerosol_mix_ratio(1, 1, i)
     &        , nd_profile, nd_latitude, nd_longitude, 1, nd_layer
     &        , nd_cdl_dimen, nd_cdl_dimen_size
     &        , nd_cdl_data, nd_cdl_var
     &        )
            IF (ierr /= i_normal) STOP
          ELSE
            aerosol_mix_ratio(:, :, i) = 0.0_RealK
          ENDIF
!
!         To compute the radiation field for a moist aerosol
!         the profile of water vapour is needed.
          IF ( (i_aerosol_parametrization(i) == IP_aerosol_param_moist)
     &         .AND.l_q_unread) THEN
            IF (index_water <= 0) THEN
              WRITE(iu_err, '(/a)')
     &           '*** Error: Water vapour must be included in '
     &           //'the spectral file to treat moist aerosols.'
              ierr=i_err_fatal
              STOP
            ENDIF
            file_name(1: length_name+1+len_file_suffix)
     &         =base_name(1: length_name)
     &         //'.'//gas_suffix(IP_h2o)(1:len_file_suffix)
            l_vert_assignable=.NOT.l_vert_coord
            CALL assign_input_vert_cdl(ierr
     &        , file_name(1: length_name+1+len_file_suffix)
     &        , 'mixing ratios of water'
     &        , l_vert_coord, name_vert_coord
     &        , .true., n_layer, l_vert_assignable
     &        , n_latitude, latitude, n_longitude, longitude
     &        , 1
     &        , n_profile, n_layer
     &        , p, water_mix_ratio(1, 1)
     &        , nd_profile, nd_latitude, nd_longitude, 1, nd_layer
     &        , nd_cdl_dimen, nd_cdl_dimen_size
     &        , nd_cdl_data, nd_cdl_var
     &        )
            IF (ierr /= i_normal) STOP
            l_q_unread=.false.
          ENDIF
        ENDIF
      ENDDO
!
!
!
      RETURN
      END
