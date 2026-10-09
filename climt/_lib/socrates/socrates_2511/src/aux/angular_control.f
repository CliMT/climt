! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the control for angular integration.
!
! Purpose:
!   This subroutine sets the controlling options for angular 
!   integration and also assigns the viewing geometry.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE angular_control(ierr
     &  , base_name, length_name
     &  , n_latitude, latitude, n_longitude, longitude
     &  , n_profile, n_layer, isolir, zen_0, l_azim_0, azim_0
     &  , i_angular_integration, l_rescale, n_order_forward
     &  , i_2stream, i_solver
     &  , n_order_gauss
     &  , i_truncation, ls_global_trunc, ms_min, ms_max
     &  , accuracy_adaptive, euler_factor
     &  , l_henyey_greenstein_pf, l_lanczos, ls_brdf_trunc
     &  , i_sph_algorithm, n_order_phase_solar
     &  , n_direction, direction
     &  , n_viewing_level, viewing_level, i_sph_mode
!           Derived dimensions
     &  , nd_sph_coeff
     &  , nd_2sg_profile, nd_flux_profile
     &  , nd_radiance_profile, nd_j_profile, nd_viewing_level
!           Sizes of arrays
     &  , nd_latitude, nd_longitude
     &  , nd_profile, nd_layer, nd_direction, nd_max_order
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
      USE rad_ccf, ONLY: pi
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
!     Treatment of errors:
      INTEGER   !,Intent(OUT)
     &    ierr
!            Error flag
!
      INCLUDE 'cdl_struc.finc'
!
!     Sizes of dummy arrays
      INTEGER, Intent(IN) ::
     &    nd_latitude
!           Allowed size for latitudes
     &  , nd_longitude
!           Allowed size for longitudes
     &  , nd_profile
!           Maximum number of profiles
     &  , nd_layer
!           Maximum number of layers
     &  , nd_direction
!           Maximum number of dierctions of radiance
     &  , nd_max_order
!           Allowed size for spherical harmonics
!
!
!     Input files:
      CHARACTER !, Intent(IN)
     &    base_name*80
!           Base name of input files
      INTEGER, Intent(IN) ::
     &    length_name
!           Length of name of input file
!
!     Horizontal grid
      INTEGER, Intent(INOUT) ::
     &    n_latitude
!           Number of atmospheric profiles
     &  , n_longitude
!           Number of atmospheric profiles
      REAL  (RealK), Intent(INOUT) ::
     &    latitude(nd_latitude)
!           Latitudes
     &  , longitude(nd_longitude)
!           Longitudes
!
!     General atmospheric properties
      INTEGER, Intent(INOUT) ::
     &    n_profile
!           Number of profiles
      INTEGER, Intent(IN) ::
     &    n_layer
!           Number of layers
!
!     Spectral region
      INTEGER, Intent(IN) ::
     &    isolir
!
!     Solar Fields
      LOGICAL, Intent(IN) ::
     &    l_azim_0
!           Flag recording the presence of solar azimuthal angles
      REAL  (RealK), Intent(INOUT) ::
     &    zen_0(nd_profile)
!           Secants or cosines of solar zenith angles
     &  , azim_0(nd_profile)
!           Azimuthal solar angles
!
!
!     Angular integration
      INTEGER, Intent(OUT) ::
     &    i_angular_integration
!           Type of angular integration
     &  , i_2stream
!           Type of two-stream scheme
     &  , n_order_gauss
!           Order of Gaussian integration (IR only)
     &  , i_truncation
!           Type of truncation
     &  , ls_global_trunc
!           Overall order of truncation with spherical harmonics
     &  , ms_min
!           Lowest azimuthal order calculated
     &  , ms_max
!           Highest azimuthal order calculated
     &  , ls_brdf_trunc
!           Order of truncation applied to BRDFs
     &  , n_order_forward
!           Order of term used to `define' the 
!           forward scattering fraction
     &  , i_sph_mode
!           Mode in which the spherical harmonic solver is to be used
     &  , i_sph_algorithm
!           Algorithm used for spherical harmonic calculation
     &  , n_order_phase_solar
!           Number of terms in the solar phase function retained
!           in an iterative solution
      LOGICAL, Intent(INOUT) ::
     &    l_rescale
!           Flag for rescaling of the phase function
     &  , l_henyey_greenstein_pf
!           Flag for Henyey-Greenstein phase functions
     &  , l_lanczos
!           Flag to use Lanczos smoothing of solar phf
      REAL  (RealK), Intent(OUT) ::
     &    accuracy_adaptive
!           Accuracy used in adaptive truncation
     &  , euler_factor
!           Factor applied to the last term of an alternating series
!
!     Derived dimensions:
      INTEGER, Intent(OUT) ::
     &    nd_sph_coeff
!           Size of Array for spherical coefficients
     &  , nd_2sg_profile
!           Number of profiles in arrays of fluxes
     &  , nd_flux_profile
!           Number of profiles in arrays of output fluxes
     &  , nd_radiance_profile
!           Number of profiles in arrays of radiances
     &  , nd_j_profile
!           Number of profiles in arrays of photolysis
     &  , nd_viewing_level
!           Number of layers in arrays of radiances
!
!
!
!     Options for solvers:
      INTEGER, Intent(INOUT) ::
     &    i_solver
!           Type of main solver
!
!
!     Viewing geometry
      INTEGER, Intent(OUT) ::
     &    n_direction
!           Number of viewing directions
      REAL  (RealK), Intent(OUT) ::
     &    direction(nd_profile, nd_direction, 2)
!           Viewing angles
!
!     Levels of viewing
      INTEGER, Intent(OUT) ::
     &    n_viewing_level
!           Number of levels where radiances are required
      REAL  (RealK), Intent(OUT) ::
     &    viewing_level(nd_layer+1)
!           List of required levels: in this driver radiances are
!           calculated at all possible levels in the atmosphere
!
!
!     Local variables:
      CHARACTER
     &    file_name*80
!           Full name of current input file
     &  , char_yn*1
!           Character response variable
!
!     Controlling variables:
      INTEGER
     &    id
!           Loop variable
     &  , i
!           Loop variable
     &  , l
!           Loop variable
!
!
!
!     Subroutines called:
      EXTERNAL
     &    assign_viewing_geom_cdl
!
!
!
!
!
      WRITE(iu_stdout, '(a)') 'Enter type of angular integration'
      READ(iu_stdin, *) i_angular_integration
!
      IF (i_angular_integration == IP_two_stream) THEN
!
        WRITE(iu_stdout, '(a)') 'Enter two-stream scheme.'
        READ(iu_stdin, *) i_2stream
        WRITE(iu_stdout, '(a)') 'Enter rescaling flag.'
        READ(iu_stdin, *) l_rescale
        n_order_forward=2
!       The standard two-stream rescaling is based the square of the
!       asymmetry parameter which is really an assumption that the
!       phase function is of the Henyey-Greenstein form
        IF (l_rescale) l_henyey_greenstein_pf=.true.
!
        WRITE(iu_stdout, '(a)') 'Enter solver.'
        READ(iu_stdin, *) i_solver
!
!       Arrays of fluxes must be of the full size.
        nd_2sg_profile=nd_profile
        nd_flux_profile=nd_profile
        nd_radiance_profile=1
        nd_j_profile=1
        nd_viewing_level=1
        nd_sph_coeff=1
!
!       Convert the zenith angles to secants.
        IF (isolir == IP_solar) THEN
          DO l=1, n_profile
            zen_0(l)=1.0_RealK/cos(pi*zen_0(l)/180.0_RealK)
          ENDDO
        ENDIF
!
      ELSE IF (i_angular_integration == IP_ir_gauss) THEN
!
        WRITE(iu_stdout, '(a)') 'Enter order of gaussian quadrature.'
        READ(iu_stdin, *) n_order_gauss
!
!       Arrays of fluxes must be of the full size.
        nd_2sg_profile=nd_profile
        nd_flux_profile=nd_profile
        nd_j_profile=1
        nd_radiance_profile=1
        nd_viewing_level=1
        nd_sph_coeff=1
!
      ELSE IF (i_angular_integration == IP_spherical_harmonic) THEN
!
!       When the reduced method of calculation is used this refers to
!       the direct spherical harmonic core of the algorithm.
        WRITE(iu_stdout, '(a)') 'Enter type of truncation.'
        READ(iu_stdin, *) i_truncation
        DO WHILE ( (i_truncation /= IP_trunc_triangular).AND.
     &             (i_truncation /= IP_trunc_rhombohedral).AND.
     &             (i_truncation /= IP_trunc_azim_sym).AND.
     &             (i_truncation /= IP_trunc_adaptive) )
          WRITE(iu_err, '(/a)')
     &      '+++ This number does not represent a supported truncation.'
          WRITE(iu_err, '(a)') 'Please reenter.'
          READ(iu_stdin, *) i_truncation
        ENDDO
!
!       Non-axially symmetric truncations should not be used in the IR.
        IF ( (isolir == IP_infra_red) .AND. 
     &       (i_truncation /= IP_trunc_azim_sym) ) THEN
          WRITE(iu_err, '(/a, /a)')
     &      '+++ Note: An axially symmetric truncation should be  
     &       used in the infra-red. The truncation has been reset.' 
          i_truncation = IP_trunc_azim_sym
        ENDIF
!
        WRITE(iu_stdout, '(a)') 'Enter rescaling flag.'
        READ(iu_stdin, *) l_rescale
!
        WRITE(iu_stdout, '(a)') 'Enter order of truncation.'
        READ(iu_stdin, *) ls_global_trunc
!
!       Check the dimension of angular arrays to see that it will
!       be sufficient, noting that we need an extra order above the
!       nominal truncation.
        DO WHILE (nd_max_order < ls_global_trunc+1)
          WRITE(iu_err, '(/a, /a, i5, /a)') 
     &      '+++ Space is not available for this order of truncation.'
     &      , 'the maximum is one less than ', nd_max_order
     &      , 'Please enter a lower number.'
          READ(iu_stdin, *) ls_global_trunc
        ENDDO
!
        IF ( (i_truncation == IP_trunc_triangular).OR.
     &       (i_truncation == IP_trunc_rhombohedral).OR.
     &       (i_truncation == IP_trunc_adaptive) ) THEN
          WRITE(iu_stdout, '(a)') 'Enter lowest azimuthal order.'
          READ(iu_stdin, *) ms_min
!
          WRITE(iu_stdout, '(a)') 'Enter highest azimuthal order.'
          READ(iu_stdin, *) ms_max
!
          IF ( (ms_max < ms_min).OR.
     &         (ms_min < 0) ) THEN
            WRITE(iu_err, '(/a)')
     &        '*** Error: Illegal range of azimuthal orders.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!
          IF (ms_max > ls_global_trunc) THEN
            WRITE(iu_err, '(/a, /a)')
     &        '*** Warning: Maximum azimuthal orders exceeds order'
     &        , 'of truncation: it will be reset.'
            ms_max=ls_global_trunc
          ENDIF
!
        ELSE IF (i_truncation == IP_trunc_azim_sym) THEN
!
!         Only the zeroth azimuthal order can be relevant.
          ms_min=0
          ms_max=0
!
        ENDIF
!
        IF (i_truncation == IP_trunc_adaptive) THEN
          READ(iu_stdin, *) accuracy_adaptive
        ENDIF
!
!       The simplest form of Euler's transformation can be applied
!       to sum the alternating series of spherical harmonics, that is,
!       only half the final order is added.
        WRITE(iu_stdout, '(a)')
     &    'Should Euler''s transformation be applied? (Y/N)'
        READ(iu_stdin, '(a)') char_yn
        DO WHILE ( (char_yn /= 'y').AND.
     &             (char_yn /= 'Y').AND.
     &             (char_yn /= 'n').AND.
     &             (char_yn /= 'N') ) 
          READ(iu_stdin, '(a)') char_yn
          WRITE(iu_err, '(/a)')
     &      '+++ invalid response: please reenter.'
        ENDDO
        IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
          euler_factor=0.5_RealK
        ELSE
          euler_factor=1.0_RealK
        ENDIF
!
        IF (l_rescale) THEN
!         The forward scattering fraction is defined by rescaling the
!         phase function to zero the term just beyond the truncation.
          n_order_forward=ls_global_trunc+1
        ENDIF
!
        WRITE(iu_stdout, '(a)') 
     &    'Enter order of truncation applied to brdfs.'
        READ(iu_stdin, *) ls_brdf_trunc
!
        l_henyey_greenstein_pf=.false.
        WRITE(iu_stdout, '(a)') 
     &    'Are henyey-greenstein phase functions to be used? (y/n)'
64      read(iu_stdin, *) char_yn
!
        IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
!
!         Set the higher order terms of the phase function.
          l_henyey_greenstein_pf=.true.
!
        ELSE IF ( (char_yn /= 'N').AND.(char_yn /= 'n') ) THEN
          WRITE(iu_err, '(/a)') 
     &      '+++ Illegal response: please reenter.'
          goto 64
        ENDIF
!
        WRITE(iu_stdout, '(a, /a)') 
     &    'Enter "R", "F" or "J" to calculate radiances, hemispheric'
     &    , 'fluxes or mean radiances.'
66      read(iu_stdin, *) char_yn
!
        IF ( (char_yn == 'F').OR.(char_yn == 'f') ) THEN
!
          i_sph_mode=IP_sph_mode_flux
!         The zeroth azimuthal order is required.
          IF (ms_min > 0) THEN
            WRITE(iu_err, '(/a)')
     &        '*** Error: The zeroth azimuthal order is required'
     &        //'to obtain hemispheric fluxes.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!
!         When fluxes are calculated a direct spherical harmonic
!         solution should be used.
          i_sph_algorithm=IP_sph_direct
!
        ELSE IF ( (char_yn == 'R').OR.(char_yn == 'r') ) THEN
!
          i_sph_mode=IP_sph_mode_rad
!
          IF (isolir == IP_solar) THEN
            IF (.NOT.l_azim_0) THEN
              WRITE(iu_err, '(/a)')
     &          '*** Error: A calculation of radiances is not possible '
     &          , 'without solar azimuthal angles.'
              RETURN
            ENDIF
          ENDIF
!
          WRITE(iu_stdout, '(a)') 
     &      'Enter the spherical harmonic algorithm.'
          READ(iu_stdin, *) i_sph_algorithm
          DO WHILE ( (i_sph_algorithm /= IP_sph_direct).AND.
     &               (i_sph_algorithm /= IP_sph_reduced_iter) )
            WRITE(iu_err, '(/a)')
     &        '+++ This number does not represent a valid algorithm.'
            WRITE(iu_err, '(a)') 'Please reenter.'
            READ(iu_stdin, *) i_sph_algorithm
          ENDDO
!
!         The test used to change the order of truncation with
!         adaptive truncation is currently formulated in terms 
!         of the amplitude of the spherical harmonics: these are
!         calculated explicitly only if the direct algorithm is
!         used.
          IF ( (i_sph_algorithm /= IP_sph_reduced_iter).AND.
     &         (i_truncation == IP_trunc_adaptive) ) THEN
            WRITE(iu_err, '(/a, /a, //a)')
     &        '*** Error: The adaptive truncation is not currently'
     &        , 'available with the fast spherical harmonic algorithm.'
     &        , 'the direct algorithm will be used instead.'
            i_sph_algorithm=IP_sph_direct
          ENDIF
!
          IF ( (i_sph_algorithm == IP_sph_reduced_iter).AND.
     &         (isolir == IP_solar) ) THEN
            WRITE(iu_stdout, '(/a)')
     &        'Enter the order of truncation of the solar beam.'
            READ(iu_stdin, *) n_order_phase_solar
            IF (n_order_phase_solar < ls_global_trunc) THEN
              WRITE(iu_stdout, '(a, /a)')
     &          'The solar order of truncation is less than the '
     &          , 'order for direct calculation: it will be reset.'
              n_order_phase_solar=ls_global_trunc
            ENDIF

            l_lanczos=.false.
            WRITE(iu_stdout, '(a)') 
     &        'Use Lanczos smoothing of the solar phase function? (y/n)'
            READ(iu_stdin, *) char_yn
            IF ( (char_yn == 'Y').OR.(char_yn == 'y') ) THEN
              l_lanczos=.true.
            ENDIF
          ENDIF
!
        ELSE IF ( (char_yn == 'J').OR.(char_yn == 'j') ) THEN
!
          nd_j_profile=nd_profile
!
          i_sph_mode=IP_sph_mode_j
!         The zeroth azimuthal order is required.
          IF (ms_min > 0) THEN
            WRITE(iu_err, '(/a)')
     &        '*** Error: The zeroth azimuthal order is required'
     &        //'to obtain mean radiances.'
            ierr=i_err_fatal
            RETURN
          ENDIF
!
!         When mean radiances are calculated a direct spherical
!         harmonic solution should be used.
          i_sph_algorithm=IP_sph_direct
!
        ELSE 
!
          WRITE(iu_err, '(/a)') 
     &      '+++ Illegal response: please reenter.'
          goto 66
!
        ENDIF
!
        IF (isolir == IP_solar) THEN
!         Convert the zenith angles to cosines.
          DO l=1, n_profile
            zen_0(l)=cos(pi*zen_0(l)/180.0_RealK)
          ENDDO
!         Convert azimuthal angles to radians.
          IF (l_azim_0) THEN
            DO l=1, n_profile
              azim_0(l)=(pi/1.8e+02_RealK)*azim_0(l)
            ENDDO
          ENDIF
        ENDIF
!
!
!       Set the derived dimensions:
!
!       Arrays of radiances must be of the full size.
        IF (i_sph_mode == IP_sph_mode_flux) THEN
          nd_flux_profile=nd_profile
        ELSE IF (i_sph_mode == IP_sph_mode_rad) THEN
          nd_flux_profile=1
        ENDIF
        nd_2sg_profile=1
        nd_radiance_profile=nd_profile
        nd_j_profile=1
        nd_viewing_level=nd_layer+1
!
        IF ( (i_truncation == IP_trunc_triangular).OR.
     &       (i_truncation == IP_trunc_adaptive) ) THEN
!
!         We need to consider coefficients of harmonics covering
!         the range of orders 0<=l<=LS_GLOBAL_TRUNC+2 and 
!         MS_MIN<=m<=MS_MAX since coefficients for negative m can 
!         be obtained from those for positive m by symmetry. The 
!         addition of 2 to the polar order is necessary to obtain 
!         an even number of polar orders for each azimuthal order, 
!         as required by the method of solution below and to carry 
!         yet one more order for the solar source. This gives l+3-m 
!         terms for each m.
!
!         In the adaptive case we must reserve equally many terms
!         for safety.
!
          nd_sph_coeff=(ms_max-ms_min+1)*(ls_global_trunc+3)
     &      +(ms_min*(ms_min-1)-ms_max*(ms_max+1))/2
!
        ELSE IF (i_truncation == IP_trunc_rhombohedral) THEN
!
!         Reset the maximum azimuthal order to a value consistent
!         with the truncation if too high a value has been supplied.
          ms_max
     &      =min(ms_max, (ls_global_trunc+mod(ms_min, 2)+ms_min-1)/2)
!         Harmonics are counted roughly as above: for each valid
!         m, with a little redundancy, there are L'+4+M1-2m harmonics.
!         (If refining this be careful to note that the number of
!         orders does not fall off by 2 with every order because of
!         the extra term we need to retain when m is odd to get an
!         even number of harmonics).
          nd_sph_coeff=(ms_max-ms_min+1)*(ls_global_trunc+4+ms_min)
     &      +ms_min*(ms_min-1)-ms_max*(ms_max+1)
!
        ELSE IF (i_truncation == IP_trunc_azim_sym) THEN
!
!         Here, only terms with m=0 are considered, so there is
!         only one term for each allowed l, including l=0.
          nd_sph_coeff=ls_global_trunc+2
!
        ELSE
!
          WRITE(iu_err, '(/a)')
     &      '*** Error: Illegal truncation.'
          RETURN
!
        ENDIF
!
        IF (i_sph_mode == IP_sph_mode_rad) THEN
!
!         The viewing geometry:
          file_name(1: length_name+1+len_file_suffix)
     &      =base_name(1: length_name)
     &      //'.'//phys_suffix(IP_view_geom)(1:len_file_suffix)
          CALL assign_viewing_geom_cdl(ierr
     &      , file_name(1: length_name+1+len_file_suffix)
     &      , n_latitude, latitude, n_longitude, longitude
     &      , n_profile
     &      , n_direction, direction
     &      , n_viewing_level, viewing_level
     &      , nd_profile, nd_latitude, nd_longitude
     &      , nd_direction, nd_viewing_level
     &      , nd_cdl_dimen, nd_cdl_dimen_size, nd_cdl_data, nd_cdl_var
     &      )
          IF (ierr /= i_normal) RETURN
!
!         Check that the viewing levels are in range.
          DO i=1, n_viewing_level
            IF ( (viewing_level(i) < 0.0_RealK) .OR.
     &           (viewing_level(i) > n_layer) ) THEN
              WRITE(iu_err, '(/a)')
     &          '*** Error: Viewing levels are out of range.'
              ierr=i_err_fatal
              RETURN
            ENDIF
          ENDDO
!
          IF (isolir == IP_infra_red) THEN
!           Dummy solar azimuthal angles are initialized here to allow
!           the use of unified code later. This choice is equivalent
!           to sunlight travelling in the direction phi=0.
            DO l=1, n_profile
              azim_0(l)=0
            ENDDO 
          ENDIF
!
!         Rotate the grid about the vertical to bring to define the
!         zero of azimuthal angles as the solar direction.
          DO id=1, n_direction
            DO l=1, n_profile
              direction(l, id, 2)=direction(l, id, 2)-azim_0(l)
            ENDDO
          ENDDO
!
        ELSE
!
!         The radiation field is required on all levels.
          n_viewing_level=n_layer+1
          DO i=1, n_layer+1
            viewing_level(i)=REAL(i-1, RealK)
          ENDDO
        ENDIF
!
      ELSE
!
        WRITE(iu_err, '(/a)') 
     &    '*** Error: illegal type of angular integration.'
        ierr=i_err_fatal
        RETURN
!
      ENDIF
!
!
!
      RETURN
      END
