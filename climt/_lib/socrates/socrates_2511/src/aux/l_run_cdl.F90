! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to run the radiation code using input from CDL files.
!
! Method:
!       Options are requested, files are read in and RADIANCE_CALC
!       is called. This driver works with fields given in layers.
!
!- ---------------------------------------------------------------------
PROGRAM l_run_cdl

! Modules to set types of variables:
  USE realtype_rd
  USE def_spectrum
  USE def_mcica,   ONLY: StrMcica, read_mcica_data, ip_mcica_full_sampling, &
                         ip_mcica_single_sampling, ip_mcica_optimal_sampling
  USE def_dimen,   ONLY: StrDim
  USE def_control, ONLY: StrCtrl,  allocate_control,   deallocate_control
  USE def_atm,     ONLY: StrAtm,   allocate_atm,       deallocate_atm
  USE def_cld,     ONLY: StrCld,   allocate_cld,       deallocate_cld, &
                                   allocate_cld_prsc,  deallocate_cld_prsc, &
                                   allocate_cld_mcica, deallocate_cld_mcica
  USE def_aer,     ONLY: StrAer,   allocate_aer,       deallocate_aer, &
                                   allocate_aer_prsc,  deallocate_aer_prsc
  USE def_bound,   ONLY: StrBound, allocate_bound,     deallocate_bound
  USE def_out,     ONLY: StrOut,                       deallocate_out
  USE dimensions_field_ucf
  USE dimensions_cdl_ucf
  USE dimensions_fixed_pcf
  USE dimensions_spec_ucf, ONLY: npd_phase_term
  USE def_std_io_icf
  USE rad_pcf
  USE gas_list_pcf
  USE input_head_pcf
  USE rad_ccf, ONLY: seconds_per_day, grav_acc, cp_air_dry, r_gas_dry, &
                     mol_weight_air, pi
  USE socrates_set_spectrum, only: set_spectrum
  USE socrates_set_cld_mcica, only: set_cld_mcica
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

  CHARACTER (LEN=*), PARAMETER :: RoutineName = 'L_RUN_CDL'

! Declaration of variables.

! Treatment of errors:
  INTEGER :: ierr = i_normal
!       Error flag
  INTEGER :: ios
!       I/O error flag
  CHARACTER (LEN=errormessagelength) :: cmessage
!       Error message
  LOGICAL :: l_interactive
!       Switch for interactive use

! Dimensions:
  TYPE(StrDim) :: dimen

! Control options:
  TYPE(StrCtrl) :: control

! Spectral information:
  TYPE(StrSpecData) :: spectrum

! MCICA data
  TYPE(StrMcica) :: mcica

! Atmospheric properties:
  TYPE(StrAtm) :: atm

! Cloud properties:
  TYPE(StrCld) :: cld

! Aerosol properties:
  TYPE(StrAer) :: aer

! Boundary conditions:
  TYPE(StrBound) :: bound

! Output fields:
  TYPE(StrOut) :: radout

! Input files:
  CHARACTER  (LEN=80) :: base_name
!       Base name of input files
  CHARACTER  (LEN=80) :: file_name
!       Full name of current input file
  CHARACTER  (LEN=1)  :: char_yn
!       Character response variable
  INTEGER :: i_gas
!       Species of gas to be read
  INTEGER :: length_name
!       Length of name of input file
  LOGICAL :: l_q_unread = .FALSE.
!       Flag for profile of specific humidity
  LOGICAL :: l_vert_coord = .FALSE.
!       Flag for vertical coordinate
  LOGICAL :: l_vert_assignable = .FALSE.
!       Flag Permitting the values of the vertical coordinates to be set
  LOGICAL :: l_vert_coord_level = .FALSE.
!       Flag for vertical coordinate on levels

! Physical processes included
  CHARACTER  (LEN=10) :: process_flag
!       String of process letters

! Ouput properties:
  INTEGER :: n_channel
!       Number of channels used for output diagnostics

! General atmospheric properties
  INTEGER :: n_latitude = 0
!       Number of latitudes
  INTEGER :: n_longitude = 0
!       Number of longitudes
  INTEGER :: n_level
!       Number of levels

! Solar Fields
  LOGICAL :: l_azim_0 = .FALSE.
!       Flag recording the presence of solar azimuthal angles

! Cloudy fields
  REAL  (RealK) :: rad_mcica_sigma
!       Normalized cloud condensate std. dev

! Aerosol fields
  LOGICAL :: l_opt_overlay
!       Flag to overlay optical properties

! Fluxes and radiances calculated
  REAL  (RealK), ALLOCATABLE :: flux_diffuse_down(:,:,:)
!       Diffuse downward flux
  REAL  (RealK), ALLOCATABLE :: flux_net(:,:,:)
!       Net flux
  REAL  (RealK), ALLOCATABLE :: heating_rate(:,:,:)
!       Heating rates
  REAL  (RealK), ALLOCATABLE :: contrib_func_i(:,:,:)
!       Contribution function (to the upwards intensity at TOA)
  REAL  (RealK), ALLOCATABLE :: contrib_func_f(:,:,:)
!       Contribution function (to the outgoing flux at TOA)
  
! Controlling variables:
  INTEGER :: i, j, l, ic, ll
!       Loop variables
  INTEGER :: b1
!       First band read in
  INTEGER :: b2
!       Second band read in
  INTEGER :: path_end
!       location of last / in location of spectral_file
  CHARACTER  (LEN=24) :: name_vert_coord
!       Name of vertical coordinate
  CHARACTER (LEN=200) :: MCICA_DATA
!             Path to McICA data file

  LOGICAL :: l_exist
!       Existence flag for file

  LOGICAL :: l_moist_aerosol
!       Flag for moist aerosol
  INTEGER :: i_pointer_water
!       Pointer to water vapour
  INTEGER, ALLOCATABLE :: i_humidity_pointer(:,:)
!       Pointer to look-up table of humidities for aerosols
  REAL  (RealK) :: delta_humidity
!       Increment in look-up table for hum.

! Dummy variables for use in calls to subroutines
! (these are generally declared with the right shape 
! for bounds checking):
  REAL  (RealK), ALLOCATABLE :: dummy_vert_coord(:,:)


! External functions:
  LOGICAL, EXTERNAL :: set_interactive
!   Function to set the flag for interactive operation


! ------------------------------------------------------------------
! Initialisation:
! ------------------------------------------------------------------
! Set the flag for interactive operation
  l_interactive=set_interactive()

  control%i_cloud_representation = IP_cloud_type_homogen
  atm%n_profile = 0
  cld%n_condensed = 0

! ------------------------------------------------------------------
! Spectral Data:
! ------------------------------------------------------------------
! Read in the spectral file.
  WRITE(iu_stdout, "(a)") "Enter the name of the spectral file."
  READ(iu_stdin, "(a)") control%spectral_file
! Derive path for McICA data file
  path_end=SCAN(control%spectral_file,'/',.TRUE.)
  mcica_data(:PATH_END)=control%spectral_file(:PATH_END)
  mcica_data(PATH_END+1:)='mcica_data'

! Read spectral file
  CALL set_spectrum( &
    spectrum      = spectrum, &
    spectral_file = control%spectral_file, &
    l_all_gases   = .true. )

  WRITE(iu_stdout, "(a,i0,a)") &
    "Enter the number of channels for output (1 - ", &
    Spectrum%Basic%n_band, "):"
  READ(iu_stdin, *) n_channel
  IF ((n_channel < 1) .OR. (n_channel > Spectrum%Basic%n_band)) THEN
    WRITE(*, "(a)") "Number of channels out of range, setting to 1."
    n_channel = 1
  END IF

! ------------------------------------------------------------------
! Input Files:
! ------------------------------------------------------------------
! All input fields are supplied in CDL format. The suffix
! determines the contents.
  WRITE(iu_stdout, '(a)') 'Enter the base-name of the cdl files.'
  READ(iu_stdin, '(a)') base_name
  j=len(base_name)
  DO
    IF (j == 0) THEN
      WRITE(iu_err, '(a)') 'No name was supplied.'
      ierr=i_err_fatal
      STOP
    ENDIF
    IF (base_name(j:j) == ' ') THEN
      j=j-1
    ELSE
      EXIT
    ENDIF
  ENDDO
  length_name=j


! ------------------------------------------------------------------
! Diagnostics
! ------------------------------------------------------------------
  control%l_blue_flux_surf = .FALSE.
  control%l_cloud_absorptivity = .FALSE.
  control%l_cloud_extinction = .FALSE.
  control%l_ls_cloud_absorptivity = .FALSE.
  control%l_ls_cloud_extinction = .FALSE.
  control%l_cnv_cloud_absorptivity = .FALSE.
  control%l_cnv_cloud_extinction = .FALSE.
  control%l_flux_direct_band = .FALSE.
  control%l_flux_direct_div_band = .FALSE.
  control%l_flux_direct_sph_band = .FALSE.
  control%l_flux_down_band = .FALSE.
  control%l_flux_up_band = .FALSE.
  control%l_flux_direct_clear_band = .FALSE.
  control%l_flux_direct_clear_div_band = .FALSE.
  control%l_flux_direct_clear_sph_band = .FALSE.
  control%l_flux_down_clear_band = .FALSE.
  control%l_flux_up_clear_band = .FALSE.
  control%l_actinic_flux = .FALSE.
  control%l_actinic_flux_clear = .FALSE.
  control%l_actinic_flux_band = .FALSE.
  control%l_actinic_flux_clear_band = .FALSE.
  control%l_aerosol_absorption_band = .FALSE.
  control%l_aerosol_scattering_band = .FALSE.
  control%l_aerosol_asymmetry_band = .FALSE.
  control%l_spherical_path_diag = .FALSE.
  control%l_contrib_func = .FALSE.
  control%l_contrib_func_band = .FALSE.

! ------------------------------------------------------------------
! Set array sizes and allocate arrays:
! ------------------------------------------------------------------

  dimen%nd_profile                = npd_profile
  dimen%nd_flux_profile           = npd_profile
  dimen%nd_2sg_profile            = npd_profile
  dimen%nd_radiance_profile       = npd_profile
  dimen%nd_j_profile              = 1
  dimen%nd_layer                  = npd_layer
  dimen%nd_layer_clr              = npd_layer
  dimen%id_cloud_top              = 1
  dimen%nd_channel                = n_channel
  dimen%nd_column                 = npd_column
  dimen%nd_max_order              = npd_max_order
  dimen%nd_direction              = npd_direction
  dimen%nd_viewing_level          = npd_layer
  dimen%nd_brdf_basis_fnc         = npd_brdf_basis_fnc
  dimen%nd_brdf_trunc             = npd_brdf_trunc
  dimen%nd_profile_aerosol_prsc   = npd_profile_aerosol_prsc
  dimen%nd_profile_cloud_prsc     = npd_profile_cloud_prsc
  dimen%nd_opt_level_aerosol_prsc = npd_opt_level_aerosol_prsc
  dimen%nd_opt_level_cloud_prsc   = npd_opt_level_cloud_prsc
  dimen%nd_cloud_component        = npd_cloud_component
  dimen%nd_cloud_type             = npd_cloud_type
  dimen%nd_overlap_coeff          = npd_overlap_coeff
  dimen%nd_source_coeff           = npd_source_coeff
  dimen%nd_region                 = npd_region
  dimen%nd_point_tile             = 1
  dimen%nd_tile                   = 1
  dimen%nd_subcol_gen             = 1
  dimen%nd_subcol_req             = 1
  dimen%nd_aerosol_mode           = 1

  CALL allocate_control(control, spectrum)
  CALL allocate_atm(atm, dimen, spectrum)
  CALL allocate_cld(cld, dimen, spectrum)
  CALL allocate_aer(aer, dimen, spectrum)
  CALL allocate_bound(bound, dimen, spectrum)

  ALLOCATE( flux_diffuse_down(dimen%nd_profile, 0:dimen%nd_layer,       &
              dimen%nd_channel)                                       )
  ALLOCATE( flux_net(dimen%nd_profile, 0: dimen%nd_layer,               &
              dimen%nd_channel)                                       )
  ALLOCATE( heating_rate(dimen%nd_profile, dimen%nd_layer,              &
              dimen%nd_channel)                                       )
  IF (control%l_contrib_func) THEN
    ALLOCATE( contrib_func_i(dimen%nd_profile, dimen%nd_layer,          &
                dimen%nd_channel)                                     )
    ALLOCATE( contrib_func_f(dimen%nd_profile, dimen%nd_layer,          &
                dimen%nd_channel)                                     )
  ELSE
    ALLOCATE( contrib_func_i(1, 1, 1) )
    ALLOCATE( contrib_func_f(1, 1, 1) )
  END IF


! ------------------------------------------------------------------
! Spectral Region:
! ------------------------------------------------------------------
  WRITE(iu_stdout, '(/a)') 'Enter the spectral region.'
2 read(iu_stdin, *) control%isolir

  IF (control%isolir == IP_solar) THEN

    IF (.NOT.Spectrum%Basic%l_present(2)) THEN
      WRITE(iu_err, '(/a)')                                             &
        '*** Error: The spectral file contains ' //                     &
        'no solar spectral data.'
      STOP
    ENDIF

!   Assign the solar zenith angles from the input file.
!   They will be converted to trigonometric functions later.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' //                               &
      phys_suffix(IP_solar_zenith_angle)
    CALL assign_input_novert_cdl(ierr,                                  &
      file_name(1: length_name+1+len_file_suffix),                      &
      'solar zenith angles',                                            &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, bound%zen_0,                                       &
      dimen%nd_profile, npd_latitude, npd_longitude,                    &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )

    IF (ierr /= i_normal) STOP
!   The solar azimuthal angles.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' //                               &
      phys_suffix(IP_solar_azimuth)
    CALL assign_input_novert_cdl(ierr,                                  &
      file_name(1: length_name+1+len_file_suffix),                      &
      'solar azimuthal angles',                                         &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, bound%azim_0,                                      &
      dimen%nd_profile, npd_latitude, npd_longitude,                    &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )

    IF (ierr == i_normal) THEN
!     The azimuthal angles have been read and a radiance 
!     calculation will be feasible.
      l_azim_0=.true.
    ELSE IF (ierr == i_err_exist) THEN
!     It will not be possible to carry out a radiance calcaultaion
!     without azimuthal angles, but a two-stream calculation will
!     be possible: the error flag is reset to allow us to continue
!     but L_AZIM_0 will be checked later.
      ierr=i_normal
    ELSE 
!     Other errors are fatal.
      STOP
    ENDIF
!   The file of solar irradiances.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' //                               &
      phys_suffix(IP_solar_toa)
    CALL assign_input_novert_cdl(ierr,                                  &
      file_name(1: length_name+1+len_file_suffix),                      &
      'solar irradiances',                                              &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, bound%solar_irrad,                                 &
      dimen%nd_profile, npd_latitude, npd_longitude,                    &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )

    IF (ierr /= i_normal) STOP


  ELSE IF (control%isolir == IP_infra_red) THEN
    IF (.NOT.Spectrum%Basic%l_present(6)) THEN
      WRITE(iu_err, '(/a)')                                             &
        '*** Error: The spectral file contains no data for the ' //     &
        'Planckian function.' 
      STOP
    ENDIF

    IF (Spectrum%Basic%l_present(2)) THEN
      control%l_solar_tail_flux = .TRUE.
    ENDIF
  ELSE 

    WRITE(iu_err, '(a)') '+++ Unknown spectral region:'
    IF (l_interactive) THEN
      WRITE(iu_err, '(a)') 'Please reenter.'
      goto 2
    ELSE
      STOP
    ENDIF

  ENDIF


! ------------------------------------------------------------------
! Range of Bands:
! ------------------------------------------------------------------
! Restrict the calculation to the range of active bands.
  WRITE(iu_stdout, '(/a)') 'Enter range of active bands.'
3 read(iu_stdin, *, iostat=ios) b1, b2
  IF (ios /= 0) THEN
    WRITE(iu_err, '(a)') '+++ unrecognized response: '
    IF (l_interactive) THEN
      WRITE(iu_err, '(a)') 'Please reenter.'
      goto 3
    ELSE
      STOP
    ENDIF
  ENDIF
  control%last_band=max(b1, b2)
  control%first_band=min(b1, b2)
  IF ( (control%last_band  > Spectrum%Basic%n_band) .OR.                &
       (control%first_band < 1) ) THEN
    WRITE(iu_err, '(a)') 'Bands out of range: '
    IF (l_interactive) THEN
      WRITE(iu_err, '(a)') 'Please reenter.'
      goto 3
    ELSE
      STOP
    ENDIF
  ENDIF

! Map spectral bands into output channels
  IF (n_channel == 1) THEN
    control%map_channel(1:Spectrum%Basic%n_band)=1
  ELSE IF (n_channel == control%last_band - control%first_band + 1) THEN
    DO i = 1, n_channel
      control%map_channel(control%first_band + i-1)=i
    END DO
  ELSE
    WRITE(iu_stdout, '(a)') 'Enter channel map.'
    READ(iu_stdin, *) control%map_channel(control%first_band:control%last_band)
  END IF

! Calculate the weighting for the bands.
  control%weight_band = 1.0_RealK
  WRITE(iu_stdout, '(a)')                                               &
    'Weight the fluxes with a filter function? (y/n)'
  READ(iu_stdin, '(a)') char_yn
  IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
    CALL filter_function(ierr,                                          &
      Spectrum%Basic%n_band,                                            &
      Spectrum%Basic%wavelength_short,                                  &
      Spectrum%Basic%wavelength_long,                                   &
      control%weight_band,                                              &
      Spectrum%Dim%nd_band)
    IF (ierr /= i_normal) STOP
  ENDIF

! ------------------------------------------------------------------
! Scaling of optical depth for direct flux
! ------------------------------------------------------------------
  WRITE(iu_stdout, '(/a)')                                              &
    'Entre treatment of optical depth for direct solar flux (0/1/2)'                            
  WRITE(iu_stdout, '(/a)')                                              &
    '0: no scaling; 1: delta-scaling; 2: circumsolar scaling'
  READ(iu_stdin, *) control%i_direct_tau 
  IF (control%i_direct_tau == 2 ) then
    WRITE(iu_stdout, '(/a)')                                            &
    'Entre half angle of FOV, between 0.25 and 5 degree'
    READ(iu_stdin, *) control%half_angle
  ENDIF

! ------------------------------------------------------------------
! Options:
! ------------------------------------------------------------------
! Obtain the various option flags.
  WRITE(iu_stdout, '(/a)') 'Enter process flag.'
  READ(iu_stdin, '(a)') process_flag
  DO j=1, len(process_flag)
    IF (process_flag(j: j) == 'r') THEN
      control%l_rayleigh=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(3)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no rayleigh scattering data.'
          STOP
      ENDIF

    ELSE IF (process_flag(j: j) == 'a') THEN
      control%l_aerosol=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(11)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no aerosol data.'
        STOP
      ENDIF

    ELSE IF (process_flag(j: j) == 'g') THEN
      control%l_gas=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(5)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no gaseous absorption data.'
        STOP
      ENDIF

    ELSE IF (process_flag(j: j) == 'c') THEN
      control%l_continuum=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(9)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no continuum absorption data.'
        STOP
      ENDIF

    ELSE IF (process_flag(j: j) == 'u') THEN
      control%l_cont_gen=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(19)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no generalised continuum absorption data.'
        STOP
      ENDIF

    ELSE IF (process_flag(j: j) == 'd') THEN
      control%l_drop=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(10)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no data for water droplets.'
        STOP
      ENDIF

    ELSE IF (process_flag(j: j) == 'i') THEN
      control%l_ice=.TRUE.

      IF (.NOT.Spectrum%Basic%l_present(12)) THEN
        write (iu_err, '(/a)')                                          &
          '*** Error: The spectral file contains ' //                   &
          'no data for ice crystals.'
        STOP
      ENDIF

    ENDIF
  ENDDO


! ------------------------------------------------------------------
! Gaseous Absorption:
! ------------------------------------------------------------------
  IF (control%l_gas) THEN

    WRITE(iu_stdout, '(a)')                                             &
      'Enter the treatment of gaseous overlaps.'
    READ(iu_stdin, *) control%i_gas_overlap
    DO j=control%first_band, control%last_band
      control%i_gas_overlap_band(j) = control%i_gas_overlap
    ENDDO
    IF (control%i_gas_overlap == IP_overlap_single) THEN
      WRITE(iu_stdout, '(a)') 'Enter type number of gas.'
      READ(iu_stdin, *) control%i_gas
      IF (.NOT.ANY(Spectrum%Gas%type_absorb(1:Spectrum%Gas%n_absorb) &
                   == control%i_gas)) THEN
        cmessage = 'The selected gas is not in the spectral file: '// &
          TRIM(gas_suffix(Spectrum%Gas%type_absorb(control%i_gas)))
        ierr=i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    ELSE IF (control%i_gas_overlap ==                                   &
      IP_overlap_random_resort_rebin) THEN
      control%n_esft_red = -1
      DO WHILE (control%n_esft_red < 0)
        WRITE(iu_stdout, '(a)') 'Enter number of rebinned k-terms.'
        READ(iu_stdin, *) control%n_esft_red
        IF (control%n_esft_red < 0) THEN
          WRITE(iu_err, '(/a)')                                         &
            '+++ Number of rebinned k-coefficient cannot be ' //        &
            'smaller than 1.'
          WRITE(iu_err, '(a)') 'Please reenter.'
        END IF
      END DO
      IF (control%n_esft_red /= 0) THEN
        control%gpnt_split = -1.0_RealK
        DO WHILE (control%gpnt_split < 0.0_RealK .OR.                   &
                  control%gpnt_split > 1.0_RealK)
          WRITE(iu_stdout, '(a)') 'Enter g-coordinate for ' //          &
            'splitting (0 or 1 = no split).'
          READ(iu_stdin, *) control%gpnt_split
          IF (control%gpnt_split < 0 .OR.                               &
                  control%gpnt_split > 1.0_RealK) THEN
            WRITE(iu_err, '(/a)')                                       &
              '+++ g-coordinate of split must be a number ' //          &
              'between 0 and 1.'
            WRITE(iu_err, '(a)') 'Please reenter.'
          END IF
        END DO
      END IF
    ENDIF

  ENDIF


  DO i_gas=1, Spectrum%Gas%n_absorb

!   Read in each file of mixing ratios. The values of the vertical
!   coordinates are set only on the first occasion when a
!   vertical coordinate is found.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' //                               &
      gas_suffix(Spectrum%Gas%type_absorb(i_gas))(1:len_file_suffix)
    INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix)),      &
            exist=l_exist)
    IF (l_exist) THEN
      l_vert_assignable=.NOT.l_vert_coord
      CALL assign_input_vert_cdl(ierr,                                  &
        file_name(1: length_name+1+len_file_suffix),                    &
        'gaseous mixing ratios',                                        &
        l_vert_coord, name_vert_coord,                                  &
        .true., atm%n_layer, l_vert_assignable,                         &
        n_latitude, atm%lat, n_longitude, atm%lon,                      &
        1,                                                              &
        atm%n_profile, atm%n_layer,                                     &
        atm%p, atm%gas_mix_ratio(1, 1, i_gas),                          &
        dimen%nd_profile, npd_latitude, npd_longitude, 1,               &
        dimen%nd_layer,                                                 &
        npd_cdl_dimen, npd_cdl_dimen_size,                              &
        npd_cdl_data, npd_cdl_var )
      IF (ierr /= i_normal) STOP
    ELSE
      atm%gas_mix_ratio(:, :, i_gas) = 0.0_RealK
    ENDIF

!   Water vapour may be required for other elsewhere (as in
!   the case of moist aerosols), so a check is made to see
!   if it has been read.
    l_q_unread=l_q_unread.OR.                                           &
      (.NOT.(Spectrum%Gas%type_absorb(i_gas) == IP_h2o))

  ENDDO


! ------------------------------------------------------------------
! Aerosol Processes:
! ------------------------------------------------------------------
  IF (control%l_aerosol) THEN
!   On occasion it is desirable to use specified profiles of optical
!   properties for aerosols. To do this an option to overlay 
!   observational data is provided. If this is not set the
!   code will not check for the presence of files of optical
!   properties; but if it is set observed propeties will override any
!   set in the spectral file.
    WRITE(iu_stdout, '(a)')                                             &
      'Overlay profiles of optical properties ? (y/n)'
    READ(iu_stdin, '(a)') char_yn
    IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
      l_opt_overlay=.true.
      dimen%nd_profile_aerosol_prsc   = dimen%nd_profile
      dimen%nd_phf_term_aerosol_prsc  = npd_phase_term
    ELSE
      l_opt_overlay=.false.
      dimen%nd_profile_aerosol_prsc   = 1
      dimen%nd_opt_level_aerosol_prsc = 1
      dimen%nd_phf_term_aerosol_prsc  = 1
    ENDIF
    aer%mr_source = ip_aersrc_classic_ron
  ELSE
    dimen%nd_profile_aerosol_prsc   = 1
    dimen%nd_opt_level_aerosol_prsc = 1
    dimen%nd_phf_term_aerosol_prsc  = 1
    aer%mr_source = ip_aersrc_classic_roff
  ENDIF

  CALL allocate_aer_prsc(aer, dimen, spectrum)
  DO i = 1, Spectrum%Dim%nd_aerosol_species
    aer%mr_type_index(i)=i
  END DO

  IF (control%l_aerosol) THEN
    CALL input_aerosol_cdl(ierr,                                        &
      base_name, length_name,                                           &
      Spectrum%Basic%n_band, Spectrum%Cont%index_water,                 &
      l_vert_coord, name_vert_coord, l_q_unread, l_opt_overlay,         &
      n_latitude, atm%lat, n_longitude, atm%lon, atm%n_profile,         &
      atm%n_layer, atm%p,                                               &
      Spectrum%Aerosol%n_aerosol, Spectrum%Aerosol%type_aerosol,        &
      aer%mix_ratio, Spectrum%Aerosol%i_aerosol_parm,                   &
      atm%gas_mix_ratio(1, 1, Spectrum%Cont%index_water),               &
      aer%n_opt_level_prsc, aer%n_phase_term_prsc,                      &
      aer%pressure_prsc, aer%absorption_prsc,                           &
      aer%scattering_prsc, aer%phase_fnc_prsc,                          &
      dimen%nd_profile, npd_latitude, npd_longitude, dimen%nd_layer,    &
      Spectrum%Dim%nd_band, Spectrum%Dim%nd_aerosol_species,            &
      dimen%nd_phf_term_aerosol_prsc,                                   &
      dimen%nd_profile_aerosol_prsc, dimen%nd_opt_level_aerosol_prsc,   &
      npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP
  ENDIF


! ------------------------------------------------------------------
! Clouds:
! ------------------------------------------------------------------
  l_exist = .FALSE.
  IF (control%l_drop) THEN
!   Determine whether a file of prescribed optical properties exists.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' // phys_suffix(IP_optical_water)
    INQUIRE(file=TRIM(file_name(1:length_name+1+len_file_suffix)),      &
      exist=l_exist)
  ENDIF
  IF (.NOT.l_exist .AND. control%l_ice) THEN
!   Determine whether a file of prescribed optical properties exists.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' // phys_suffix(IP_optical_ice)
    INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix)),      &
      exist=l_exist)
  ENDIF
  IF (l_exist) THEN
    dimen%nd_phf_term_cloud_prsc  = npd_phase_term
  ELSE
    dimen%nd_profile_cloud_prsc   = 1
    dimen%nd_opt_level_cloud_prsc = 1
    dimen%nd_phf_term_cloud_prsc  = 1
  ENDIF

  CALL allocate_cld_prsc(cld, dimen, spectrum)

  CALL input_cloud_cdl(ierr,                                            &
    Spectrum%Basic%n_band,                                              &
    Spectrum%Drop%l_drop_type, Spectrum%Drop%i_drop_parm,               &
    Spectrum%Drop%n_phf, Spectrum%Drop%parm_list,                       &
    Spectrum%Drop%parm_min_dim, Spectrum%Drop%parm_max_dim,             &
    Spectrum%Ice%l_ice_type, Spectrum%Ice%i_ice_parm,                   &
    Spectrum%Ice%n_phf, Spectrum%Ice%parm_list,                         &
    Spectrum%Ice%parm_min_dim, Spectrum%Ice%parm_max_dim,               &
    base_name, length_name,                                             &
    control%l_drop, control%l_ice,                                      &
    control%i_cloud, control%l_cloud, control%i_cloud_representation,   &
    cld%n_condensed, cld%type_condensed, cld%i_condensed_param,         &
    cld%condensed_n_phf, cld%condensed_param_list,                      &
    n_latitude, atm%lat, n_longitude, atm%lon, atm%n_profile,           &
    l_vert_coord, name_vert_coord, atm%n_layer, atm%p,                  &
    cld%condensed_mix_ratio, cld%condensed_dim_char,                    &
    cld%w_cloud, cld%n_cloud_type, cld%i_cloud_type, cld%frac_cloud,    &
    cld%n_opt_level_drop_prsc, cld%n_phase_term_drop_prsc,              &
    cld%drop_pressure_prsc, cld%drop_absorption_prsc,                   &
    cld%drop_scattering_prsc, cld%drop_phase_fnc_prsc,                  &
    cld%n_opt_level_ice_prsc, cld%n_phase_term_ice_prsc,                &
    cld%ice_pressure_prsc, cld%ice_absorption_prsc,                     &
    cld%ice_scattering_prsc, cld%ice_phase_fnc_prsc,                    &
    dimen%nd_profile, npd_latitude, npd_longitude, dimen%nd_layer,      &
    Spectrum%Dim%nd_band, dimen%nd_phf_term_cloud_prsc,                 &
    Spectrum%Dim%nd_drop_type, Spectrum%Dim%nd_ice_type,                &
    Spectrum%Dim%nd_cloud_parameter,                                    &
    dimen%nd_cloud_component, dimen%nd_cloud_type,                      &
    dimen%nd_profile_cloud_prsc, dimen%nd_opt_level_cloud_prsc,         &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
  IF (ierr /= i_normal) STOP

  IF ( (control%i_cloud == IP_cloud_part_corr).OR.                      &
       (control%i_cloud == IP_cloud_part_corr_cnv).OR.                  &
       (control%i_cloud == IP_cloud_mcica) ) THEN
    WRITE(iu_stdout, '(/a)') 'Enter decorrelation pressure scale ' //   &
      'for the overlap of large scale cloud.'
    READ(iu_stdin, *) cld%dp_corr_strat
  ENDIF
    
  IF (control%i_cloud == IP_cloud_part_corr_cnv) THEN
    WRITE(iu_stdout, '(/a)') 'Enter decorrelation pressure scale ' //   &
      'for the overlap of convective cloud.'
    READ(iu_stdin, *) cld%dp_corr_conv
  ENDIF

  IF (control%i_cloud == IP_cloud_mcica) THEN
    WRITE(iu_stdout, '(/a)') 'Enter relative standard deviation ' //  &
      '(std. dev / mean) of cloud water content.'
    READ(iu_stdin, *) rad_mcica_sigma
    WRITE(iu_stdout, '(/a)') 'Enter number for the sampling method ' // &
      'used in McICA,'
    READ(iu_stdin, *) control%i_mcica_sampling
    WRITE (iu_stdout,'(/a)') 'Reading McICA data file: ',mcica_data
    IF (control%isolir == ip_solar) THEN
      CALL read_mcica_data(mcica, mcica_data, sp_sw=spectrum)
    ELSE
      CALL read_mcica_data(mcica, mcica_data, sp_lw=spectrum)
    END IF
  ENDIF

! ------------------------------------------------------------------
! Angular Integration:
! ------------------------------------------------------------------
  CALL angular_control(ierr,                                            &
    base_name, length_name,                                             &
    n_latitude, atm%lat, n_longitude, atm%lon,                          &
    atm%n_profile, atm%n_layer, control%isolir, bound%zen_0,            &
    l_azim_0, bound%azim_0, control%i_angular_integration,              &
    control%l_rescale, control%n_order_forward,                         &
    control%i_2stream, control%i_solver, control%n_order_gauss,         &
    control%i_truncation, control%ls_global_trunc,                      &
    control%ms_min, control%ms_max, control%accuracy_adaptive,          &
    control%euler_factor, control%l_henyey_greenstein_pf,               &
    control%l_lanczos, control%ls_brdf_trunc,                           &
    control%i_sph_algorithm, control%n_order_phase_solar,               &
    atm%n_direction, atm%direction,                                     &
    atm%n_viewing_level, atm%viewing_level, control%i_sph_mode,         &
!                 Derived dimensions
    dimen%nd_sph_coeff, dimen%nd_2sg_profile, dimen%nd_flux_profile,    &
    dimen%nd_radiance_profile, dimen%nd_j_profile,                      &
    dimen%nd_viewing_level,                                             &
!                 Sizes of arrays
    npd_latitude, npd_longitude, dimen%nd_profile,                      &
    dimen%nd_layer, dimen%nd_direction, dimen%nd_max_order,             &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
  IF (ierr /= i_normal) STOP

! Reset dimen%nd_max_order to reduce memory requirements
  IF ( (control%i_angular_integration == IP_two_stream).OR.             &
       (control%i_angular_integration == IP_ir_gauss) ) THEN
    dimen%nd_max_order = 1
  END IF


! ------------------------------------------------------------------
! Surface properties:
! ------------------------------------------------------------------
  file_name(1: length_name+1+len_file_suffix) =                         &
    base_name(1: length_name)//'.'//phys_suffix(IP_surface_char)
  CALL assign_surface_char_cdl(ierr,                                    &
    control%i_angular_integration, control%isolir,                      &
    Spectrum%Basic%n_band, file_name(1: length_name+1+len_file_suffix), &
    'surface characteristics', bound%zen_0,                             &
    Spectrum%Basic%wavelength_short, Spectrum%Basic%wavelength_long,    &
    n_latitude, atm%lat, n_longitude, atm%lon, atm%n_profile,           &
    bound%n_brdf_basis_fnc, control%ls_brdf_trunc, bound%rho_alb,       &
    bound%f_brdf, dimen%nd_profile, npd_latitude, npd_longitude,        &
    Spectrum%Dim%nd_band, dimen%nd_brdf_basis_fnc, dimen%nd_brdf_trunc, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
  IF (ierr /= i_normal) STOP


! ------------------------------------------------------------------
! Treatment of scattering
! ------------------------------------------------------------------
  WRITE(iu_stdout, '(/a)') 'Enter treatment of scattering.'
  READ(iu_stdin, *) control%i_scatter_method
  DO i=control%first_band, control%last_band
    control%i_scatter_method_band(i)=control%i_scatter_method
  ENDDO


! ------------------------------------------------------------------
! Temperatures
! ------------------------------------------------------------------
  file_name(1: length_name+1+len_file_suffix) =                         &
    base_name(1: length_name)//'.'//phys_suffix(IP_temperature)
  l_vert_assignable=.NOT.l_vert_coord
  CALL assign_input_vert_cdl(ierr,                                      &
    file_name(1: length_name+1+len_file_suffix), 'temperatures',        &
    l_vert_coord, name_vert_coord,                                      &
    .true., atm%n_layer, l_vert_assignable,                             &
    n_latitude, atm%lat, n_longitude, atm%lon, 1,                       &
    atm%n_profile, atm%n_layer,                                         &
    atm%p, atm%t,                                                       &
    dimen%nd_profile, npd_latitude, npd_longitude, 1, dimen%nd_layer,   &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
  IF (ierr /= i_normal) STOP

! In the infra-red region the temperatures at the edges of levels
! are required: the pressures at which these temperatures are set
! `define' the edges of the layers and in the solar, only the
! pressures are of concern.
  file_name(1: length_name+1+len_file_suffix) =                         &
    base_name(1: length_name)//'.'//phys_suffix(IP_t_level)
  CALL assign_input_vert_cdl(ierr,                                      &
    file_name(1: length_name+1+len_file_suffix),                        &
    'temperatures on levels',                                           &
    l_vert_coord_level, name_vert_coord,                                &
    .true., atm%n_layer+1, .true.,                                      &
    n_latitude, atm%lat, n_longitude, atm%lon, 0,                       &
    atm%n_profile, n_level,                                             &
    atm%p_level, atm%t_level,                                           &
    dimen%nd_profile, npd_latitude, npd_longitude, 0, dimen%nd_layer,   &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
  IF (ierr /= i_normal) STOP
! Check that the number of levels in this field is consistent.
  IF (n_level /= (atm%n_layer+1)) THEN
    WRITE(iu_err, '(/a)') '*** Error: The number of temperatures ' //   &
      'at layer boundaries is not consistent.'
  ENDIF


! ------------------------------------------------------------------
! Pressures
! ------------------------------------------------------------------

  file_name(1: length_name+1+len_file_suffix) =                         &
    base_name(1: length_name)//'.'//phys_suffix(IP_pressure)
  INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix)),        &
          exist=l_exist)
  IF (l_exist) THEN
    ALLOCATE( dummy_vert_coord(dimen%nd_profile, dimen%nd_layer) )
    CALL assign_input_vert_cdl(ierr,                                    &
      file_name(1: length_name+1+len_file_suffix), 'pressures',         &
      l_vert_coord, name_vert_coord,                                    &
      .true., atm%n_layer, .false.,                                     &
      n_latitude, atm%lat, n_longitude, atm%lon, 1,                     &
      atm%n_profile, atm%n_layer,                                       &
      dummy_vert_coord, atm%p,                                          &
      dimen%nd_profile, npd_latitude, npd_longitude, 1, dimen%nd_layer, &
      npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP
    DEALLOCATE( dummy_vert_coord )
  ENDIF

  file_name(1: length_name+1+len_file_suffix) =                         &
    base_name(1: length_name)//'.'//phys_suffix(IP_p_level)
  INQUIRE(file=Trim(file_name(1:length_name+1+len_file_suffix)),        &
          exist=l_exist)
  IF (l_exist) THEN
    ALLOCATE( dummy_vert_coord(dimen%nd_profile, 0:dimen%nd_layer) )
    CALL assign_input_vert_cdl(ierr,                                    &
      file_name(1: length_name+1+len_file_suffix),                      &
      'pressures on levels',                                            &
      l_vert_coord_level, name_vert_coord,                              &
      .true., atm%n_layer+1, .false.,                                   &
      n_latitude, atm%lat, n_longitude, atm%lon, 0,                     &
      atm%n_profile, n_level,                                           &
      dummy_vert_coord, atm%p_level,                                    &
      dimen%nd_profile, npd_latitude, npd_longitude, 0, dimen%nd_layer, &
      npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP
!   Check that the number of levels in this field is consistent.
    IF (n_level /= (atm%n_layer+1)) THEN
      WRITE(iu_err, '(/a)') '*** Error: The number of pressures ' //    &
        'at layer boundaries is not consistent.'
    ENDIF
    DEALLOCATE( dummy_vert_coord )
  ENDIF


! ------------------------------------------------------------------
! Surface Temperatures
! ------------------------------------------------------------------
  IF (control%isolir == IP_infra_red) THEN
!   Surface Temperatures are required.
    file_name(1: length_name+1+len_file_suffix) =                       &
      base_name(1: length_name) // '.' //                               &
      phys_suffix(IP_temperature_ground)
    CALL assign_input_novert_cdl(ierr,                                  &
      file_name(1: length_name+1+len_file_suffix),                      &
      'surface temperatures',                                           &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, bound%t_ground,                                    &
      dimen%nd_profile, npd_latitude, npd_longitude,                    &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP
  ENDIF


! ------------------------------------------------------------------
! Variation of the temperature within layers.
! ------------------------------------------------------------------
  IF (control%isolir == IP_infra_red) THEN
    WRITE(iu_stdout, '(/a)') 'Is the ir-source function to be ' //      &
      'taken as linear or quadratic in tau? (l/q)'
6   read(iu_stdin, '(a)') char_yn
    IF ( (char_yn == 'L').OR.(char_yn == 'l') ) THEN
      control%l_ir_source_quad=.FALSE.
    ELSE IF ( (char_yn == 'Q').OR.(char_yn == 'q') ) THEN
      control%l_ir_source_quad=.TRUE.
    ELSE
      WRITE(iu_err, '(a)') '+++ Illegal response'
      IF (l_interactive) THEN
        WRITE(iu_err, '(a)') 'Please reenter.'
        goto 6
      ELSE
        STOP
      ENDIF
    ENDIF
  ENDIF


! ------------------------------------------------------------------
! Determination of the mass and density in each layer.
! ------------------------------------------------------------------
  IF (l_vert_coord_level) THEN
    DO i=atm%n_layer, 1, -1
      DO l=1, atm%n_profile
        atm%mass(l, i)=(atm%p_level(l, i)-atm%p_level(l, i-1))/grav_acc
      ENDDO
    ENDDO
  ELSE
    WRITE(iu_err, '(/a)')                                               &
      '*** Error: No coordinates on the boundaries have been set.'
    STOP
  ENDIF
  IF (Spectrum%Cont%index_water > 0) THEN
    DO i=1, atm%n_layer
      DO l=1, atm%n_profile
        atm%density(l, i)=atm%p(l, i)/(r_gas_dry*atm%t(l, i)*(1.0e+00_RealK &
          + (mol_weight_air/(molar_weight(ip_h2o)*1.0E-03_RealK) - 1.0_RealK) &
          * atm%gas_mix_ratio(l, i, Spectrum%Cont%index_water)))
      END DO
    END DO
  ELSE
    DO i=1, atm%n_layer
      DO l=1, atm%n_profile
        atm%density(l, i)=atm%p(l, i)/(r_gas_dry*atm%t(l, i))
      END DO
    END DO   
  END IF


! ------------------------------------------------------------------
! Subgrid cloud generator
! ------------------------------------------------------------------
  IF (control%i_cloud == IP_cloud_mcica) THEN
    dimen%nd_subcol_gen = mcica%n_subcol_gen
    select case (control%i_mcica_sampling)
    case (ip_mcica_full_sampling)
      dimen%nd_subcol_req = mcica%n_subcol_gen
    case (ip_mcica_single_sampling)
      dimen%nd_subcol_req = mcica%n_subcol_req_single
    case (ip_mcica_optimal_sampling)
      dimen%nd_subcol_req = mcica%n_subcol_req_optimal
    end select
    control%i_overlap=2
    cld%c_cloud = 0.0_RealK
    cld%c_ratio = 0.0_RealK
    cld%condensed_rel_var_dens = rad_mcica_sigma
    CALL set_cld_mcica(cld, mcica, control, dimen, spectrum, atm)
  END IF


! ------------------------------------------------------------------
! Calculate the clear-sky mean relative humidity for moist aerosols.
! ------------------------------------------------------------------
l_moist_aerosol=.FALSE.
DO j=1, spectrum%aerosol%n_aerosol
  l_moist_aerosol=l_moist_aerosol.OR.                                          &
    (spectrum%aerosol%i_aerosol_parm(j) == ip_aerosol_param_moist).OR.         &
    (spectrum%aerosol%i_aerosol_parm(j) == ip_aerosol_param_phf_moist)
END DO

IF (l_moist_aerosol) THEN
  i_pointer_water=MAX(spectrum%cont%index_water, 1)
  ALLOCATE( i_humidity_pointer(dimen%nd_profile, dimen%nd_layer) )
  ! DEPENDS ON: set_moist_aerosol_properties
  CALL set_moist_aerosol_properties(ierr,                                      &
    atm%n_profile, atm%n_layer,                                                &
    spectrum%aerosol%n_aerosol, spectrum%aerosol%i_aerosol_parm,               &
    spectrum%aerosol%nhumidity,                                                &
    atm%gas_mix_ratio(1, 1, i_pointer_water), atm%t, atm%p, cld%w_cloud,       &
    delta_humidity, aer%mean_rel_humidity, i_humidity_pointer,                 &
    dimen%nd_profile, dimen%nd_layer, dimen%id_cloud_top,                      &
    spectrum%dim%nd_aerosol_species)
  DEALLOCATE( i_humidity_pointer )
END IF


! ------------------------------------------------------------------
! Calculation of radiances or irradiances.
! ------------------------------------------------------------------
  CALL radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)


! ------------------------------------------------------------------
! Processing of the output.
! ------------------------------------------------------------------
  IF ( (control%i_angular_integration == IP_two_stream).OR.             &
       (control%i_angular_integration == IP_ir_gauss).OR.               &
       ( (control%i_angular_integration == IP_spherical_harmonic).AND.  &
         (control%i_sph_mode == IP_sph_mode_flux) ) ) THEN

    DO ic=1, n_channel
      DO i=0, atm%n_layer
        DO l=1, atm%n_profile
          flux_net(l, i, ic) =                                          &
            radout%flux_down(l, i, ic)-radout%flux_up(l, i, ic)
        ENDDO
      ENDDO
      DO i=1, atm%n_layer
        DO l=1, atm%n_profile
          heating_rate(l, i, ic) =                                      &
            (flux_net(l, i-1, ic)-flux_net(l, i, ic)) /                 &
            (atm%mass(l, i)*cp_air_dry)
!         Convert heating rates to the conventional units of K/day.
          heating_rate(l, i, ic) =                                      &
            seconds_per_day*heating_rate(l, i, ic)
        ENDDO
      ENDDO
      IF (control%isolir == IP_solar) THEN
        DO i=0, atm%n_layer
          DO l=1, atm%n_profile
            flux_diffuse_down(l, i, ic) =                               &
              radout%flux_down(l, i, ic)-radout%flux_direct(l, i, ic)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    IF (control%l_contrib_func) THEN
      contrib_func_i = radout%contrib_funci
      contrib_func_f = radout%contrib_funcf
    END IF

!   Write the output files.
    CALL output_flux_cdl(ierr,                                          &
      control, spectrum,                                                &
      base_name, length_name,                                           &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, atm%n_layer,                                       &
      trim(name_vert_coord), len(trim(name_vert_coord)),                &
      atm%p, atm%p_level,                                               &
      n_channel,                                                        &
      radout%flux_down, flux_diffuse_down,                              &
      radout%flux_up, flux_net,                                         &
      radout%flux_direct, heating_rate,                                 &
      contrib_func_i, contrib_func_f,                                   &
      dimen%nd_profile, npd_latitude, npd_longitude, dimen%nd_layer,    &
      dimen%nd_channel,                                                 &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP

  ENDIF

  IF ( (control%i_angular_integration == IP_spherical_harmonic).AND.    &
       (control%i_sph_mode == IP_sph_mode_j) ) THEN

!   Rates of photolysis:
    CALL output_photolysis_cdl(ierr,                                    &
      base_name, length_name,                                           &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, atm%n_layer,                                       &
      atm%n_viewing_level, atm%viewing_level,                           &
      n_channel, radout%photolysis,                                     &
      dimen%nd_profile, npd_latitude, npd_longitude, dimen%nd_layer,    &
      dimen%nd_viewing_level, dimen%nd_channel,                         &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP

  ENDIF

  IF ( (control%i_angular_integration == IP_spherical_harmonic).AND.    &
       (control%i_sph_mode == IP_sph_mode_rad) ) THEN

!   Radiances:
    CALL output_radiance_cdl(ierr,                                      &
      base_name, length_name,                                           &
      n_latitude, atm%lat, n_longitude, atm%lon,                        &
      atm%n_profile, atm%n_direction, atm%direction, bound%azim_0,      &
      atm%n_viewing_level, atm%viewing_level,                           &
      n_channel, radout%radiance,                                       &
      dimen%nd_profile, npd_latitude, npd_longitude, dimen%nd_layer,    &
      dimen%nd_viewing_level, dimen%nd_direction, dimen%nd_channel,     &
      npd_cdl_dimen, npd_cdl_dimen_size,                                &
      npd_cdl_data, npd_cdl_var )
    IF (ierr /= i_normal) STOP

  ENDIF


! ------------------------------------------------------------------
! Deallocate arrays:
! ------------------------------------------------------------------
  CALL deallocate_cld_mcica(cld)
  CALL deallocate_atm(atm)
  CALL deallocate_cld(cld)
  CALL deallocate_cld_prsc(cld)
  CALL deallocate_aer(aer)
  CALL deallocate_aer_prsc(aer)
  CALL deallocate_bound(bound)
  CALL deallocate_out(radout)

  DEALLOCATE( flux_diffuse_down         )
  DEALLOCATE( flux_net                  )
  DEALLOCATE( heating_rate              )
  DEALLOCATE( contrib_func_i            )
  DEALLOCATE( contrib_func_f            )


END PROGRAM l_run_cdl
