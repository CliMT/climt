! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to set aerosol mixing ratios for SRA profiles.
!
! Method:
!   The pressure, temperatures and specific humidity are read
!   in as CDL-files. These quantities enable the height at
!   each pressure level to be calculated: the SRA profiles are
!   specified in terms of heights. For the profile selected
!   the extinction at 550 nm for each component is calculated
!   and the volume fraction is set to be consistent with the
!   extinction.
!
!- ---------------------------------------------------------------------
      PROGRAM sra_aerosol_profile
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_spec_ucf
      USE dimensions_cdl_ucf
      USE rad_ccf, ONLY: grav_acc
      USE gas_list_pcf
      USE rad_pcf
      USE aerosol_profile_pcf
      USE aerosol_model_pcf
      USE aerosol_representation_pcf
      USE def_std_io_icf
      USE input_head_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables.
      CHARACTER
     &    base_name*80
!           Base names of input files
     &  , file_name*80
!           Name of current file
     &  , name_vert_coord*24
!           Name of vertical coordinate
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
!
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_level
!           Number of levels
     &  , i_aerosol_profile
!           Index of aerosol profile
      LOGICAL
     &    l_vert_coord
!           Flag asserting that vertical coordinate is set
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
     &  , z(npd_layer+1)
!           Heights of levels
     &  , p(npd_profile, npd_layer+1)
!           Pressures
     &  , t(npd_profile, npd_layer+1)
!           Temperatures
     &  , q(npd_profile, npd_layer+1)
!           Specific humidity
     &  , density(0: npd_layer)
!           Densities at levels
     &  , aerosol_mix_ratio(npd_profile, npd_layer+1
     &    , npd_aerosol_species)
!           Number density of aerosol
     &  , pstar(npd_profile)
!           Surface pressure
     &  , tstar(npd_profile)
!           Surface air temperature
     &  , qstar(npd_profile)
!           Surface humidity
     &  , density_star
!           Surface density
      INTEGER
     &    type_aerosol(npd_aerosol_component)
!           Types of aerosols
     &  , i
!           Loop `index'
     &  , j
!           Loop `index'
     &  , l
!           Loop `index'
     &  , i_model
!           Index of aerosol model
     &  , i_aerosol_layer
!           Index of aerosol layer
     &  , length_name
!           Length of base name
     &  , i_component
!           Temporary `index'
      REAL  (RealK) ::
     &    volume_fraction(npd_layer+1, npd_aerosol_component)
!           Calculated aerosol fractions
     &  , volume_fraction_level(npd_aerosol_component)
!           Fractions of components at levels
     &  , ext_550nm_overall
!           Overall extinction at 550nm
!
!     Numerical precision:
      REAL  (RealK) ::
     &    tol_z
!           Tolerance for equality of heights
!
!     Functions called:
      REAL  (RealK) ::
     &    extinction_profile
!           Function to calculate extinction
     &  , fnc_density
!           Function to calculate density
      EXTERNAL
     &    extinction_profile, fnc_density
!
!     Subroutines called:
      EXTERNAL
     &    assign_input_vert_cdl, calc_volume_fraction
     &  , output_vert_cdl

!     Include data for aerosols
      INCLUDE 'aerosol_component.finc'


      data
     &    l_vert_coord/.false./
      data ierr/i_normal/
!
!
!
!
!
!     Set the components of the standard profiles.
      data
     &    type_aerosol(1)/IP_water_soluble/
     &  , type_aerosol(2)/IP_dust_like/
     &  , type_aerosol(3)/IP_oceanic/
     &  , type_aerosol(4)/IP_soot/
     &  , type_aerosol(5)/IP_ash/
     &  , type_aerosol(6)/IP_sulphuric/
!
!
!
!     Stop the program because conversion aerosol_profile_pcf is
!     not yet complete.
      WRITE(iu_err, '(/A)')
     &  '*** Error: This program has not been made operational at F90.'
      IF (EPSILON(tol_z) < 1.0) STOP
!
!     Set the numerical precision
      tol_z=1.0e+03_RealK*epsilon(tol_z)
!
!     Obtain the files containing the grid on which to set the profile.
      WRITE(iu_stdout, '(/a)')
     &  'Enter the base-names of the files to be used.'
      READ(iu_stdin, '(a)') base_name
      j=len(base_name)
1     IF (base_name(j:j) == ' ') THEN
        j=j-1
        IF (j == 0) THEN
          WRITE(iu_err, '(a)') 'No name was supplied.'
          ierr=i_err_fatal
          STOP
        ELSE
          goto 1
        ENDIF
      ENDIF
      length_name=j
!
!     Assign the concentration values for each type of aerosol.
      WRITE(iu_stdout, '(/a)')
     &  'Enter the type of the sra-profile to be set.'
2     READ(iu_stdin, *, iostat=ios) i_aerosol_profile
      IF (ios /= 0) THEN
        WRITE(iu_err, '(a)') '+++ Illegal response:'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'Please re-type.'
          goto 2
        ENDIF
      ENDIF
!
      IF ( (i_aerosol_profile <= 0).OR.
     &     (i_aerosol_profile > npd_aerosol_profile) ) THEN
        WRITE(iu_err, '(/a)') 
     &    '*** Error: Unrecognized aerosol profile.'
        IF (lock_code(.true.)) THEN
          INCLUDE 'batch_error_main.finc'
        ELSE
          WRITE(iu_stdout, '(a)') 'Please re-enter.'
          goto 2
        ENDIF
      ENDIF
!
!
!     Read the fields required to calculate the heights.
!
!     Surface pressure:
!
      file_name(1: length_name+7)=base_name(1: length_name)
     &  //'.'//phys_suffix(IP_pressure_ground)
      CALL assign_input_novert_cdl(ierr
     &  , file_name(1: length_name+7)
     &  , 'surface pressures'
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, pstar
     &  , npd_profile, npd_latitude, npd_longitude
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!     Surface temperature:
!
      file_name(1: length_name+7)=base_name(1: length_name)
     &  //'.'//phys_suffix(IP_temperature_ground)
      CALL assign_input_novert_cdl(ierr
     &  , file_name(1: length_name+7)
     &  , 'surface temperatures'
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, tstar
     &  , npd_profile, npd_latitude, npd_longitude
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!     Temperatures:
!
      file_name(1: length_name+7)=base_name(1: length_name)
     &  //'.'//phys_suffix(IP_temperature)
      CALL assign_input_vert_cdl(ierr
     &  , file_name(1: length_name+7) 
     &  , 'temperatures', l_vert_coord, name_vert_coord
     &  , .true., n_level, .NOT.l_vert_coord
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, t
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
!
!     Specific humidities:
!
      file_name(1: length_name+7)=base_name(1: length_name)
     &  //'.'//phys_suffix(IP_spec_humidity)
      CALL assign_input_vert_cdl(ierr
     &  , file_name(1: length_name+7)
     &  , 'specific humidities', l_vert_coord, name_vert_coord
     &  , .true., n_level, .true.
     &  , n_latitude, latitude, n_longitude, longitude
     &  , 1
     &  , n_profile, n_level
     &  , p, q
     &  , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &  , npd_cdl_dimen, npd_cdl_dimen_size
     &  , npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!     Integrate through the atmosphere to find the heights at each
!     specified level. The specific humidity at the surface is
!     extrapolated using the temperature.
!
      DO l=1, n_profile
        qstar(l)=q(l, n_level)
     &    +((q(l, n_level)-q(l, n_level-1))
     &    /(p(l, n_level)-p(l, n_level-1)))
     &    *((pstar(l)-p(l, n_level-1))
     &    /(p(l, n_level)-p(l, n_level-1)))
        density_star=fnc_density(pstar(l), tstar(l), qstar(l))
        density(n_level)=fnc_density(p(l, n_level), t(l, n_level)
     &    , q(l, n_level))
        z(n_level)=(pstar(l)-p(l, n_level))/(grav_acc*0.5_RealK
     &    *(density_star+density(n_level)))
        DO i=n_level-1, 1, -1
          density(i)=fnc_density(p(l, i), t(l, i), q(l, i))
          z(i)=z(i+1)+(p(l, i+1)-p(l, i))
     &      /(grav_acc*0.5_RealK*(density(i+1)+density(i)))
        ENDDO
!
!       Set all volume fractions to 0. A DO-loop is preferred to
!       a DATA statement in case of multiple calls.
        DO i=1, n_level
          DO j=1, npd_aerosol_component
            volume_fraction(i, j)=0.0_RealK
          ENDDO
        ENDDO
!
        DO i=1, n_level
!         Find the layer of aerosol at this height
          i_aerosol_layer=1
3         if (z(i).gt.z_aerosol_layer(i_aerosol_layer
     &      , i_aerosol_profile)-tol_z) THEN
            i_aerosol_layer=i_aerosol_layer+1
          IF (i_aerosol_layer < 
     &      n_aerosol_layer(i_aerosol_profile)) goto 3
          ENDIF
!
          IF (i_aerosol_layer < 
     &      n_aerosol_layer(i_aerosol_profile)) THEN
!           Calculate the overall extinction at 550nm for this profile.
!           above the top of the profiles the volume fraction is
!           left at 0.0.
            ext_550nm_overall=extinction_profile(z(i)
     &        , i_extinction_rep(i_aerosol_layer, i_aerosol_profile)
     &        , npd_extinction_parameter
     &        , extinction_parameter(1, i_aerosol_layer
     &        , i_aerosol_profile)
     &        )
!
!           Work out the appropriate volume fractions for the model
!           applicable in this layer.
            i_model=i_aerosol_model(i_aerosol_layer, i_aerosol_profile)
            CALL calc_volume_fraction(npd_aerosol_component
     &        , fraction_component(1, i_model)
     &        , ext_550nm_component, ext_550nm_overall
     &        , volume_fraction_level
     &        )
            DO j=1, npd_aerosol_component
              volume_fraction(i, j)=volume_fraction_level(j)
            ENDDO
          ENDIF
        ENDDO
!
!       Convert the volume fraction to a mixing ratio.
        DO i=1, n_level
          DO j=1, 6
            i_component=type_aerosol(j)
            aerosol_mix_ratio(l, i, j)
     &        =volume_fraction(i, i_component)
     &        *density_component(i_component)
     &        /density(i)
          ENDDO
        ENDDO
!
      ENDDO
!
!
!
!     Write out the files of aerosol components.
      DO j=1, 6
        i_component=type_aerosol(j)
        file_name(1: length_name+7)=base_name(1: length_name)
     &    //'.'//aerosol_suffix(i_component)
        CALL output_vert_cdl(ierr
     &    , file_name
     &    , n_latitude, latitude, n_longitude, longitude, n_profile
     &    , n_level, trim(name_vert_coord), len(trim(name_vert_coord))
     &    , p, 'mixrat', 'float', 'none', 'mixing ratio'
     &    , aerosol_mix_ratio
     &    , npd_profile, npd_latitude, npd_longitude, 1, npd_layer+1
     &    , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &    )
      ENDDO
!
!
!
      STOP
      END
