! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to run DISORT using input from CDL files.
!
! Method:
!    This program follows the logical the monochromatic driver and is
!    intended to run DISORT in parallel with the radiance version of
!    the code. Changes to run_mono should be mirrored here. 
!
!    IMPORTANT:
!    The direct beam in DISORT is defined without rescaling, whereas
!    the normal convention in this package is that the direct flux is
!    rescaled. This means that the direct and diffuse downward fluxes
!    are not separately comparable, though the total downward flux 
!    should be. To compare the individual components of the downward
!    flux, the spherical harmonic code must be run twice, once without
!    rescaling to get the direct beam only, and once with it to get 
!    the total downward flux. The difference between these two fluxes
!    will then be comparable with the diffuse downward flux as 
!    calculated by DISORT.
!
!    DISORT must be compiled with static storage as some variables
!    (such as TOL, the convergence criterion for Legendre polynomials
!    are not reinitialized on every call to a subroutine).
!
!- ---------------------------------------------------------------------
      PROGRAM run_disort
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE dimensions_field_ucf
      USE dimensions_cdl_ucf
      USE dimensions_spec_ucf
      USE dimensions_fixed_pcf
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
!     Declaration of variables.
!
!     Treatment of errors:
      INTEGER
     &    ierr
!           Error flag
      LOGICAL
     &    lock_code
!           Logical to forbid interactive looping
      EXTERNAL
     &    lock_code
!
!     Derived dimensions:
      INTEGER
     &    nd_sph_coeff
!           Size of Array for spherical coefficients
     &  , nd_2sg_profile
!           Number of profiles in arrays of fluxes
     &  , nd_flux_profile
!           Number of profiles in arrays of output fluxes
     &  , nd_radiance_profile
!           Number of profiles in arrays of radiances
     &  , nd_viewing_level
!           Number of layers in arrays of radiances
!
!
!
!     Input files:
      CHARACTER
     &    base_name*80
!           Base name of input files
     &  , file_name*80
!           Full name of current input file
     &  , char_yn*1
!           Character response variable
      INTEGER
     &    length_name
!           Length of name of input file
      LOGICAL
     &    l_non_net
!           Flag for non-net fluxes
     &  , l_vert_coord_level
!           Flag for vertical coordinate on levels
!
!
!     Optical Properties:
      INTEGER
     &    n_phase_term
!           Number of terms in the phase function
      LOGICAL
     &    l_forward
!           Flag indicating the presence of forward scattering data
      REAL  (RealK) ::
     &    tau(npd_profile, npd_layer)
!           Optical Depths
     &  , omega(npd_profile, npd_layer)
!           Albedos of single scattering
     &  , phase_function(npd_profile, npd_layer, npd_phase_term)
!           Phase functions
     &  , forward_scatter(npd_profile, npd_layer)
!           Forward scattering fractions
!
!     Angular integration
      INTEGER
     &    i_angular_integration
!           Type of angular integration
     &  , i_2stream
!           Type of two-stream scheme
     &  , n_order_gauss
!           Order of Gaussian integration (IR only)
     &  , i_sph_mode
!           Mode in which the spherical harmonic solver is being used
     &  , i_truncation
!           Type of truncation
     &  , ls_global_trunc
!           Overall order of truncation with spherical harmonics
     &  , ls_brdf_trunc
!           Order of truncation applied to BRDFs
     &  , n_order_forward
!           Order of term used to `define' the 
!           forward scattering fraction
     &  , i_sph_algorithm
!           Algorithm used for spherical harmonic calculation
     &  , ls_trunc_full
!           Full order of truncation used with an iterative solution
      LOGICAL
     &    l_rescale
!           Flag for rescaling of the phase function
     &  , l_ir_source_quad
!           Flag to use a quadratic source function in the IR
     &  , l_henyey_greenstein_pf
!           Flag for Henyey-Greenstein phase functions
     &  , l_lanczos
!           Flag to use Lanczos smoothing of solar phf
!
!     Viewing geometry
      INTEGER
     &    n_direction
!           Number of viewing directions
      REAL  (RealK) ::
     &    direction(npd_direction, 2)
!           Viewing angles
!     Levels of viewing
      INTEGER
     &    n_viewing_level
!           Number of levels where radiances are required
      REAL  (RealK) ::
     &    viewing_level(npd_layer+1)
!           List of required levels: in this driver radiances are
!           calculated at all possible levels in the atmosphere
!
!     Options for solvers:
      INTEGER
     &    i_solver
!           Type of main solver
!
!     General atmospheric properties
      INTEGER
     &    n_latitude
!           Number of latitudes
     &  , n_longitude
!           Number of longitudes
     &  , n_profile
!           Number of profiles
     &  , n_layer
!           Number of layers
     &  , n_level
!           Number of levels
!     Horizontal structure
      REAL  (RealK) ::
     &    latitude(npd_latitude)
!           Latitude
     &  , longitude(npd_longitude)
!           Longitude
      REAL  (RealK) ::
     &    t(npd_profile, npd_layer)
!           Temperatures in layers
     &  , t_level(npd_profile, 0: npd_layer)
!           Temperatures at the boundaries of layers
!
!     Spectral region
      INTEGER
     &    isolir
!
!     Solar Fields
      LOGICAL
     &    l_azim_0
!           Flag recording the presence of solar azimuthal angles
      REAL  (RealK) ::
     &    zen_0(npd_profile)
!           Secants or cosines of solar zenith angles
     &  , azim_0(npd_profile)
!           Azimuthal solar angles
     &  , solar_toa(npd_profile)
!           Solar irradiance at the top of the atmosphere
!
!     Infra-red fields:
      REAL  (RealK) ::
     &    wavelength_short
!           Shorter wavelength of IR region
     &  , wavelength_long
!           Longer wavelength of IR region
     &  , w_temp
!           Sorting variable
!
!     Isotropic incident radiation:
      LOGICAL
     &    l_iso_inc
!           Flag for isotropic incident radiation
     &  , iso_inc(npd_profile)
!           Isotropic incident radiance (this is isotropic over
!           the downward hemisphere, so the corresponding incident
!           flux is pi times this)
!
!
!
!     Surface properties
      INTEGER
     &    n_brdf_basis_fnc
!           Number of basis functions for BRDFs
      REAL  (RealK) ::
     &    rho_alb(npd_profile, npd_brdf_basis_fnc, npd_band)
!           Weights for basis functions of the BRDFs
     &  , f_brdf(npd_brdf_basis_fnc, 0: npd_brdf_trunc/2
     &      , 0: npd_brdf_trunc/2, 0: npd_brdf_trunc)
!           Array of BRDF basis terms
     &  , t_ground(npd_profile)
!           Temperature of the ground
!
!
!     Fluxes and radiances calculated
      REAL  (RealK) ::
     &    flux_direct(npd_profile, 0: npd_layer)
!           Direct flux
     &  , flux_down(npd_profile, 0: npd_layer)
!           Totol downward flux
     &  , flux_diffuse_down(npd_profile, 0: npd_layer)
!           Diffuse downward flux
     &  , flux_up(npd_profile, 0: npd_layer)
!           Upward flux
     &  , flux_net(npd_profile, 0: npd_layer)
!           Net flux
     &  , radiance(npd_profile, npd_layer+1, npd_direction)
!           Radiances calculated
     &  , heating_rate(npd_profile, npd_layer)
!           Heating rates
!
!     Controlling variables:
      INTEGER
     &    i
!           Loop variable
     &  , j
!           Loop variable
     &  , l
!           Loop variable
      CHARACTER
     &    name_vert_coord*24
!           Name of vertical coordinate
      LOGICAL
     &    l_exist
!           Existence flag for files
!
!     Variables required for consistency with called routines:
      INTEGER
     &    p
     &  , p_level
!
!
!
      data l_azim_0/.false./
      data l_henyey_greenstein_pf/.false./
      data l_vert_coord_level/.false./
      data l_forward/.false./
!
!
!
!
!     ------------------------------------------------------------------
!     Input Files:
!     ------------------------------------------------------------------
!
!     All input fields are supplied in CDL format. The suffix
!     determines the contents.
      WRITE(iu_stdout, '(a)') 'enter the base-name of the cdl files.'
      READ(iu_stdin, '(a)') base_name
      j=len(base_name)
1     if (base_name(j:j) == ' ') then
        j=j-1
        IF (j == 0) THEN
          WRITE(iu_err, '(a)') 'no name was supplied.'
          ierr=i_err_fatal
          STOP
        ELSE
          goto 1
        ENDIF
      ENDIF
      length_name=j
!
!
!     ------------------------------------------------------------------
!     Optical Properties:
!     ------------------------------------------------------------------
!
!     Optical Properties:
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)
     &  //'.'//phys_suffix(IP_optical_ss)(1:len_file_suffix)
      CALL assign_input_ss_cdl(ierr
     &  , file_name(1: length_name+1+len_file_suffix)
     &  , n_latitude, latitude, n_longitude, longitude, n_profile
     &  , n_layer, n_phase_term
     &  , tau, omega, phase_function, l_forward, forward_scatter
     &  , npd_profile, npd_latitude, npd_longitude, npd_layer
     &  , npd_phase_term
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!     ------------------------------------------------------------------
!     Spectral Region:
!     ------------------------------------------------------------------
!
      WRITE(iu_stdout, '(/a)') 'enter the spectral region.'
2     read(iu_stdin, *) isolir
!
      IF (isolir == IP_solar) THEN
!
!       Assign the solar zenith angles from the input file.
!       They will be converted to trigonometric functions later.
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)//'.'
     &    //phys_suffix(IP_solar_zenith_angle)
        CALL assign_input_novert_cdl(ierr
     &    , file_name(1: length_name+1+len_file_suffix)
     &    , 'solar zenith angles'
     &    , n_latitude, latitude, n_longitude, longitude
     &    , n_profile, zen_0
     &    , npd_profile, npd_latitude, npd_longitude
     &    , npd_cdl_dimen, npd_cdl_dimen_size
     &    , npd_cdl_data, npd_cdl_var
     &    )
        IF (ierr /= i_normal) STOP
!       The solar azimuthal angles.
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)
     &    //'.'//phys_suffix(IP_solar_azimuth)
        INQUIRE(file=file_name(1: length_name+1+len_file_suffix)
     &    , exist=l_azim_0)
        IF (l_azim_0) THEN
          CALL assign_input_novert_cdl(ierr
     &      , file_name(1: length_name+1+len_file_suffix)
     &      , 'solar azimuthal angles'
     &      , n_latitude, latitude, n_longitude, longitude
     &      , n_profile, azim_0
     &      , npd_profile, npd_latitude, npd_longitude
     &      , npd_cdl_dimen, npd_cdl_dimen_size
     &      , npd_cdl_data, npd_cdl_var
     &      )
        ENDIF
!       The the file of solar irradiances.
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)
     &    //'.'//phys_suffix(IP_solar_toa)
        CALL assign_input_novert_cdl(ierr
     &    , file_name(1: length_name+1+len_file_suffix)
     &    , 'solar irradiances'
     &    , n_latitude, latitude, n_longitude, longitude
     &    , n_profile, solar_toa
     &    , npd_profile, npd_latitude, npd_longitude
     &    , npd_cdl_dimen, npd_cdl_dimen_size
     &    , npd_cdl_data, npd_cdl_var
     &    )
        IF (ierr /= i_normal) STOP
!
!
      ELSE IF (isolir == IP_infra_red) THEN
!
        WRITE(iu_stdout, '(/a)') 'enter the range of wavelengths (m)'
        READ(iu_stdin, *) wavelength_short, wavelength_long
        IF (wavelength_long < wavelength_short) THEN
          w_temp=wavelength_long
          wavelength_short=wavelength_long
          wavelength_short=w_temp
        ENDIF
!
!       ----------------------------------------------------------------
!       Temperatures
!       ----------------------------------------------------------------
!
        IF (isolir == IP_infra_red) THEN
!         In the infra-red region the temperatures at the edges of
!         levels are required.
          file_name(1: length_name+1+len_file_suffix)
     &       =base_name(1: length_name)//'.'//phys_suffix(IP_t_level)
          CALL assign_input_vert_cdl(ierr
     &      , file_name(1: length_name+1+len_file_suffix)
     &      , 'temperatures on levels'
     &      , l_vert_coord_level, name_vert_coord
     &      , .true., n_layer+1, .NOT.l_vert_coord_level
     &      , n_latitude, latitude, n_longitude, longitude, 0
     &      , n_profile, n_level
     &      , p_level, t_level
     &      , npd_profile, npd_latitude, npd_longitude, 0, npd_layer
     &      , npd_cdl_dimen, npd_cdl_dimen_size
     &      , npd_cdl_data, npd_cdl_var
     &      )
          IF (ierr /= i_normal) STOP
!         Check that the number of levels in this field is consistent.
          IF (n_level /= (n_layer+1)) THEN
            WRITE(iu_err, '(/a)') 
     &        '*** Error: The number of temperatures '
     &        //'at layer boundaries is not consistent.'
          ENDIF
!
!         Surface Temperatures are also required.
          file_name(1: length_name+1+len_file_suffix)
     &      =base_name(1: length_name)
     &      //'.'//phys_suffix(IP_temperature_ground)
          CALL assign_input_novert_cdl(ierr
     &      , file_name(1: length_name+1+len_file_suffix)
     &      , 'surface temperatures'
     &      , n_latitude, latitude, n_longitude, longitude
     &      , n_profile, t_ground
     &      , npd_profile, npd_latitude, npd_longitude
     &      , npd_cdl_dimen, npd_cdl_dimen_size
     &      , npd_cdl_data, npd_cdl_var
     &      )
          IF (ierr /= i_normal) STOP
!
!         --------------------------------------------------------------
!         Variation of the temperature within layers.
!         --------------------------------------------------------------
          WRITE(iu_stdout, '(/a)') 'is the ir-source function to be '
     &       //'taken as linear or quadratic in tau? (l/q)'
6         read(iu_stdin, '(a)') char_yn
          IF ( (char_yn == 'L').OR.(char_yn == 'l') ) THEN
            l_ir_source_quad=.false.
          ELSE IF ( (char_yn == 'Q').OR.(char_yn == 'q') ) THEN
            l_ir_source_quad=.true.
          ELSE
            WRITE(iu_err, '(a)') '+++ illegal response'
            IF (lock_code(.true.)) THEN
              STOP
            ELSE
              WRITE(iu_err, '(a)') 'please reenter.'
              goto 6
            ENDIF
          ENDIF
        ENDIF
!
      ELSE 
!
        WRITE(iu_err, '(a)') '+++ unknown spectral region:'
        IF (lock_code(.true.)) THEN
          STOP
        ELSE
          WRITE(iu_err, '(a)') 'please reenter.'
          goto 2
        ENDIF
!
      ENDIF
!
!     ------------------------------------------------------------------
!     Isotropic sources
!     ------------------------------------------------------------------
      WRITE(iu_stdout, '(/a)') 'do you want to INCLUDE isotropic '
     &   //'incident radiation? (y/n)'
16    read(iu_stdin, '(a)') char_yn
      IF ( (char_yn == 'N').OR.(char_yn == 'n') ) THEN
        l_iso_inc=.false.
      ELSE IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
        l_iso_inc=.true.
!       The isotropic source at the top.
        file_name(1: length_name+1+len_file_suffix)
     &    =base_name(1: length_name)
     &    //'.'//phys_suffix(IP_iso_inc)
        INQUIRE(file=file_name(1: length_name+1+len_file_suffix)
     &    , exist=l_exist)
        IF (l_exist) THEN
          CALL assign_input_novert_cdl(ierr
     &      , file_name(1: length_name+1+len_file_suffix)
     &      , 'isotropic incident radiation'
     &      , n_latitude, latitude, n_longitude, longitude
     &      , n_profile, iso_inc
     &      , npd_profile, npd_latitude, npd_longitude
     &      , npd_cdl_dimen, npd_cdl_dimen_size
     &      , npd_cdl_data, npd_cdl_var
     &      )
        ELSE
          WRITE(iu_err, '(/a)')
     &      '*** Error: A file of isotropic incident radiation ;'
     &      //'does not exist.'
          STOP
        ENDIF
      ELSE
        WRITE(iu_err, '(a)') '+++ illegal response'
        IF (lock_code(.true.)) THEN
          STOP
        ELSE
          WRITE(iu_err, '(a)') 'please reenter.'
          goto 16
        ENDIF
      ENDIF
!
!     ------------------------------------------------------------------
!     Angular Integration:
!     ------------------------------------------------------------------
!
      CALL angular_control(ierr
     &  , base_name, length_name
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, n_layer, isolir, zen_0, l_azim_0, azim_0
     &  , i_angular_integration, l_rescale, n_order_forward
     &  , i_2stream, i_solver, l_non_net
     &  , n_order_gauss
     &  , i_truncation, ls_global_trunc, l_henyey_greenstein_pf
     &  , l_lanczos, ls_brdf_trunc
     &  , i_sph_algorithm, ls_trunc_full
     &  , n_direction, direction
     &  , n_viewing_level, viewing_level, i_sph_mode
!                       Derived dimensions
     &  , nd_sph_coeff
     &  , nd_2sg_profile, nd_flux_profile
     &  , nd_radiance_profile, nd_viewing_level
!                       Sizes of arrays
     &  , npd_latitude, npd_longitude
     &  , npd_profile, npd_layer, npd_direction
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
!
!     ------------------------------------------------------------------
!     Surface properties:
!     ------------------------------------------------------------------
!
!     Set surface properties.
      file_name(1: length_name+1+len_file_suffix)
     &  =base_name(1: length_name)//'.'//phys_suffix(IP_surface_char)
      CALL assign_surface_char_cdl(ierr
     &  , i_angular_integration, isolir, 1
     &  , file_name(1: length_name+1+len_file_suffix)
     &  , 'surface characteristics'
     &  , zen_0, wavelength_short, wavelength_long
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, n_brdf_basis_fnc, ls_brdf_trunc, rho_alb, f_brdf
     &  , npd_profile, npd_latitude, npd_longitude
     &  , npd_band, npd_brdf_basis_fnc, npd_brdf_trunc
     &  , npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var
     &  )
      IF (ierr /= i_normal) STOP
!
!
!
!     ------------------------------------------------------------------
!     Calculation of radiances.
!     ------------------------------------------------------------------
!
      CALL disort_interface(ierr
     &  , i_angular_integration, l_rescale, n_order_forward
     &  , i_2stream, i_solver, n_order_gauss
     &  , i_truncation, ls_global_trunc, ls_brdf_trunc
     &  , l_henyey_greenstein_pf
     &  , isolir, zen_0, azim_0, solar_toa
     &  , t_level, l_ir_source_quad, t, t_ground
     &  , wavelength_short, wavelength_long
     &  , l_iso_inc, iso_inc
     &  , n_profile, n_layer
     &  , tau, omega, phase_function, n_phase_term, forward_scatter
     &  , n_direction, direction
     &  , n_viewing_level, viewing_level
     &  , n_brdf_basis_fnc, f_brdf, rho_alb
     &  , flux_direct, flux_down, flux_up, radiance
     &  , npd_profile, npd_layer, npd_phase_term, nd_sph_coeff
     &  , nd_flux_profile, nd_radiance_profile
     &  , nd_viewing_level, npd_direction
     &  , npd_brdf_basis_fnc, npd_brdf_trunc
     &  , npd_source_coeff
     &  )
      IF (ierr /= i_normal) STOP
!
!
!     ------------------------------------------------------------------
!     Processing of the output.
!     ------------------------------------------------------------------
      IF ( (i_angular_integration == IP_two_stream).OR.
     &     (i_angular_integration == IP_ir_gauss).OR.
     &     ( (i_angular_integration == IP_spherical_harmonic).AND.
     &       (i_sph_flux_mode == IP_sph_mode_flux) ) ) THEN
!
        IF (l_non_net) THEN
          DO i=0, n_layer
            DO l=1, n_profile
              flux_net(l, i)=flux_down(1, i)-flux_up(1, i)
            ENDDO
          ENDDO
        ENDIF
        IF (isolir == IP_solar) THEN
          DO i=0, n_layer
            DO l=1, n_profile
              flux_diffuse_down(1, i)=flux_down(1, i)-flux_direct(1, i)
            ENDDO
          ENDDO
        ENDIF
!
!       Write the output files.
!
        CALL output_flux_cdl(ierr
     &    , base_name, length_name
     &    , isolir, l_non_net, i_angular_integration
     &    , n_latitude, latitude, n_longitude, longitude
     &    , n_profile, n_layer, 'level', 5, p, p_level
     &    , flux_down, flux_diffuse_down, flux_up
     &    , flux_net, flux_direct, heating_rate
     &    , npd_profile, npd_latitude, npd_longitude, npd_layer
     &    , npd_cdl_dimen, npd_cdl_dimen_size
     &    , npd_cdl_data, npd_cdl_var
     &    )
!
      ENDIF
!
      IF ( (i_angular_integration == IP_spherical_harmonic).AND.
     &     (i_sph_mode == IP_sph_mode_rad) ) THEN
!
!       Radiances:
!
        CALL output_radiance_cdl(ierr
     &    , base_name, length_name
     &    , n_latitude, latitude, n_longitude, longitude
     &    , n_profile, n_direction, direction, azim_0
     &    , n_viewing_level, viewing_level
     &    , radiance
     &    , npd_profile, npd_latitude, npd_longitude, npd_layer
     &    , nd_viewing_level, npd_direction
     &    , npd_cdl_dimen, npd_cdl_dimen_size
     &    , npd_cdl_data, npd_cdl_var
     &    )
!
      ENDIF
!
      IF ( (i_angular_integration == IP_spherical_harmonic).AND.
     &     (i_sph_mode == IP_sph_mode_j) ) THEN
!
!       Rates of photolysis:
!
        CALL output_photolysis_cdl(ierr
     &    , base_name, length_name
     &    , n_latitude, latitude, n_longitude, longitude
     &    , n_profile, n_layer
     &    , n_viewing_level, viewing_level
     &    , radiance
     &    , npd_profile, npd_latitude, npd_longitude, npd_layer
     &    , nd_viewing_level
     &    , npd_cdl_dimen, npd_cdl_dimen_size
     &    , npd_cdl_data, npd_cdl_var
     &    )
!
      ENDIF
!
!
!
      STOP
      END
