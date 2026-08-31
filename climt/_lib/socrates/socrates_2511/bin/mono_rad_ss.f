! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calcaulate fluxes or radiances.
!
! Purpose:
!   To calculate monochromatic radiances or fluxes from given
!   single-scattering properties.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE mono_rad_ss(ierr
     &  , control, bound
     &  , i_angular_integration, l_rescale, n_order_forward
     &  , i_2stream, i_solver, n_order_gauss
     &  , i_truncation, ls_global_trunc, ms_min, ms_max
     &  , accuracy_adaptive, euler_factor
     &  , ls_brdf_trunc, i_sph_mode
     &  , i_sph_algorithm, n_order_phase_solar
     &  , l_henyey_greenstein_pf
     &  , isolir, zen_0, solar_irrad
     &  , t_level, l_ir_source_quad, t, t_ground
     &  , wavelength_short, wavelength_long
     &  , l_iso_inc, iso_inc
     &  , n_profile, n_layer
     &  , tau, omega, phase_fnc, n_phase_term
     &  , l_forward, forward_scatter
     &  , l_forward_solar, forward_solar
     &  , n_direction, direction
     &  , n_viewing_level, viewing_level
     &  , n_brdf_basis_fnc, f_brdf, rho_alb
     &  , flux_direct, flux_down, flux_up, radiance, photolysis
     &  , nd_profile, nd_layer
     &  , nd_phase_term, nd_max_order, nd_sph_coeff
     &  , nd_flux_profile, nd_radiance_profile, nd_j_profile
     &  , nd_viewing_level, nd_direction
     &  , nd_brdf_basis_fnc, nd_brdf_trunc
     &  , nd_source_coeff
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_control, ONLY: StrCtrl
      USE def_bound,   ONLY: StrBound
      USE def_spherical_geometry, ONLY: StrSphGeo
      USE def_std_io_icf
      USE rad_pcf
      USE rad_ccf, ONLY: pi
      USE calc_brdf_mod, ONLY: calc_brdf
      USE calc_cg_coeff_mod, ONLY: calc_cg_coeff
      USE calc_uplm_sol_mod, ONLY: calc_uplm_sol
      USE calc_uplm_zero_mod, ONLY: calc_uplm_zero
      USE gauss_angle_mod, ONLY: gauss_angle
      USE rescale_phase_fnc_mod, ONLY: rescale_phase_fnc
      USE rescale_tau_omega_mod, ONLY: rescale_tau_omega
      USE set_rad_layer_mod, ONLY: set_rad_layer
      USE set_truncation_mod, ONLY: set_truncation
      USE sol_scat_cos_mod, ONLY: sol_scat_cos
      USE sph_solver_mod, ONLY: sph_solver
      USE two_stream_mod, ONLY: two_stream
!
!
      IMPLICIT NONE
!
!
!
!
!
!     Declaration of variables.
!
!     Treatment of errors:
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag

!     Control options:
      TYPE(StrCtrl), INTENT(IN) :: control

!     Boundary conditions:
      TYPE(StrBound), INTENT(IN) :: bound

!     Sizes of arrays:
      INTEGER, Intent(IN) ::
     &    nd_profile
!           Size allocated for atmospheric profiles
     &  , nd_layer
!           Size allocated for atmospheric layers
     &  , nd_phase_term
!           Size allocated for terms in the phase function
     &  , nd_max_order
!           Size allocated for orders of spherical harmonics
     &  , nd_sph_coeff
!           Size of array for spherical coefficients
     &  , nd_flux_profile
!           Number of profiles in arrays of output fluxes
     &  , nd_direction
!           Size allowed for directional radiances
     &  , nd_brdf_basis_fnc
!           Size allowed for BRDF basis functions
     &  , nd_brdf_trunc
!           Size allowed for orders of BRDFs
     &  , nd_radiance_profile
!           Number of profiles in arrays of radiances
     &  , nd_j_profile
!           Number of profiles in arrays of rates of photolysis
     &  , nd_viewing_level
!           Number of layers in arrays of radiances
     &  , nd_source_coeff
!           Size allocated for two-stream source coefficients
!
!     Angular integration
      INTEGER, Intent(IN) ::
     &    i_angular_integration
!           Type of angular integration
     &  , i_2stream
!           Type of two-stream scheme
     &  , i_solver
!           Type of solver
     &  , n_order_gauss
!           Order of Gaussian integration (IR only)
     &  , i_truncation
!           Type of spherical truncation
     &  , ls_global_trunc
!           Overall order of spherical truncation
     &  , ls_brdf_trunc
!           Order of truncation applied to BRDFs
     &  , n_order_forward
!           Order of term used to `define' the 
!           forward scattering fraction
     &  , ms_min
!           Lowest azimuthal order calculated
     &  , ms_max
!           Highest azimuthal order calculated
     &  , i_sph_mode
!           Mode in which the spherical harmonic solver is being used
     &  , i_sph_algorithm
!           Algorithm used for spherical harmonic calculation
     &  , n_order_phase_solar
!           Number of terms retained in the phase function for the
!           solar beam
      LOGICAL, Intent(IN) ::
     &    l_rescale
!           Flag for rescaling of the phase function
     &  , l_forward
!           Flag for explicit specification of forward scattering
     &  , l_forward_solar
!           Flag for explicit specification of rescaling for the solar 
!           beam
     &  , l_ir_source_quad
!           Flag to use a quadratic source function in the IR
     &  , l_henyey_greenstein_pf
!           Flag for Henyey-Greenstein phase functions
      REAL  (RealK), Intent(IN) ::
     &    accuracy_adaptive
!           Accuracy for adaptive truncation
     &  , euler_factor
!           Factor applied to the last term of the alternating series
!
!
!     Optical Properties:
      INTEGER, Intent(IN) ::
     &    n_phase_term
!           Order of terms in the phase function
!     Optical properties will be altered by rescaling and resetting for
!     Henyey-Greenstein phase functions.
      REAL  (RealK), Intent(INOUT) ::
     &    tau(nd_profile, nd_layer)
!           Optical depths
     &  , omega(nd_profile, nd_layer)
!           Albedos of single scattering
     &  , phase_fnc(nd_profile, nd_layer, nd_phase_term)
!           Optical depths
     &  , forward_scatter(nd_profile, nd_layer)
!           Forward-scattering fractions
     &  , forward_solar(nd_profile, nd_layer)
!           Forward-scattering for the solar beam
!
!     Specification of radiances required:
      INTEGER, Intent(IN) ::
     &    n_direction
!           Number of directions at which to calculate radiances
     &  , n_viewing_level
!           Number of levels where radiances are required
      REAL  (RealK), Intent(IN) ::
     &    direction(nd_radiance_profile, nd_direction, 2)
!           Directions in which to calculate radiances
     &  , viewing_level(nd_layer+1)
!           List of required levels: in this driver radiances are
!           calculated at all possible levels in the atmosphere
!
!     General atmospheric properties
      INTEGER, Intent(IN) ::
     &    n_profile
!           Number of profiles
     &  , n_layer
!           Number of layers
      REAL  (RealK), Intent(IN) ::
     &    t(nd_profile, nd_layer)
!           Temperatures at the mid-points of layers
     &  , t_level(nd_profile, 0: nd_layer)
!           Temperatures at the edges of layers
     &  , t_ground(nd_profile)
!           Temperature of the surface
!
!     Spectral region
      INTEGER, Intent(IN) ::
     &    isolir
!           Switch for spectral region
!
!     Solar Fields
      REAL  (RealK), Intent(IN) ::
     &    zen_0(nd_profile)
!           Secants or cosines of solar zenith angles
     &  , solar_irrad(nd_profile)
!           Solar irradiance at the top of the atmosphere
!
!     Infra-red properties:
      REAL  (RealK), Intent(IN) ::
     &    wavelength_short
!           Shorter wavelength of region (IR only)
     &  , wavelength_long
!           Longer wavelength of region (IR only)
!
!     Isotropic incident radiation:
      LOGICAL, Intent(IN) ::
     &    l_iso_inc
!           Flag for incident isotropic radiation
      REAL  (RealK), Intent(IN) ::
     &    iso_inc(nd_profile)
!           Incident isotropic radiance
!
!     Surface properties:
      INTEGER, Intent(IN) ::
     &    n_brdf_basis_fnc
!           Number of basis functions for BRDFs
      REAL  (RealK), Intent(IN) ::
     &    rho_alb(nd_profile, nd_brdf_basis_fnc)
!           Weights of the basis functions
     &  , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!           Array of BRDF basis terms
!
!
!     Fluxes and radiances calculated:
      REAL  (RealK), Intent(OUT) ::
     &    flux_direct(nd_flux_profile, 0: nd_layer)
!           Direct flux
     &  , flux_down(nd_flux_profile, 0: nd_layer)
!           Totol downward flux
     &  , flux_up(nd_flux_profile, 0: nd_layer)
!           Totol downward flux
     &  , radiance(nd_radiance_profile, nd_viewing_level, nd_direction)
!           Calculated radiances
     &  , photolysis(nd_j_profile, nd_viewing_level)
!           Rate of photolysis


!     Local variables:

!     Spherical geometry fields
      TYPE(StrSphGeo) :: sph

!     Spherical variables:
      INTEGER
     &    ls_max_order
!           Maximum order of terms required
     &  , ls_local_trunc(0: nd_max_order)
!           Actual truncation for each particular value of m
     &  , ms_trunc(0: nd_max_order)
!           Maximum azimuthal quantum number for each order
     &  , ia_sph_mm(0: nd_max_order)
!           Address of spherical coefficient of (m, m) for each m
     &  , n_order_phase
!           Number of terms in the phase function to be treated
!           explicitly: conventionally rescaling is carried out
!           to zero the first omitted term, so we may require a
!           term one order higher to set the rescaling parameter
     &  , n_order_required
!           Order of terms which the phase function must contain
!           for the options selected
!
!     Additional variables for angular integration:
      LOGICAL
     &    l_solar_phf
!           Logical to specify a separate treatment of the singly
!           scattered solar beam
     &  , l_rescale_solar_phf
!           Logical to apply rescaling to the singly scattered
!           solar phase function
      REAL  (RealK) ::
     &    cos_sol_view(nd_radiance_profile, nd_direction)
!           Cosines of the angles between the solar direction and
!           the viewing direction
     &  , tau_noscal(nd_profile, nd_layer)
!           Unscaled optical depths
     &  , phase_fnc_solar(nd_radiance_profile, nd_layer, nd_direction)
!           The phase function evaluated at in the viewing directions
!           for scattering from the solar beam
!     Legendre polynomials:
      REAL  (RealK) ::
     &    cnst1
!           Constant in recurrence for Legendre polynomials
     &  , p_legendre_ls(nd_radiance_profile)
!           Legendre polynomial at the current order
     &  , p_legendre_ls_m1(nd_radiance_profile)
!           Legendre polynomial at the previous order
     &  , p_legendre_tmp(nd_radiance_profile)
!           Temporary Legendre polynomial
!
      REAL  (RealK) ::
     &    cg_coeff(nd_sph_coeff)
!           Clebsch-Gordan coefficients
     &  , uplm_zero(nd_sph_coeff)
!           Values of spherical harmonics at polar angles pi/2
     &  , uplm_sol(nd_radiance_profile, nd_sph_coeff)
!           Values of spherical harmonics in the solar direction
     &  , brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)
!           The BRDF evaluated for scattering from the solar
!           beam into the viewing direction
     &  , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!           The BRDF evaluated for scattering from isotropic
!           radiation into the viewing direction
!
!     Controlling variables
      INTEGER
     &    i
!           Loop variable
     &  , l
!           Loop variable
     &  , ls
!           Loop variable
     &  , id
!           Loop variable
!
!     Derived dimensions:
      INTEGER
     &    nd_red_eigensystem
!           Size allowed for the reduced eigensystem
     &  , nd_sph_equation
!           Size allowed for spherical harmonic equations
     &  , nd_sph_diagonal
!           Size allowed for diagonals of the spherical harmonic
!           matrix
     &  , nd_sph_cf_weight
!           Size allowed for application of weights of the C. F.
     &  , nd_sph_u_range
!           Size allowed for range of values of u^+|- contributing
!           on any viewing level
!
!     Planckian variables:
      REAL  (RealK) ::
     &    diff_planck(nd_profile, nd_layer)
!           Difference in pi*Planckian function
     &  , diff_planck_2(nd_profile, nd_layer)
!           2x2nd differences of Planckian
     &  , d_planck_flux_surface(nd_profile)
!           Differential Planckian flux from the surface
     &  , planck_flux(nd_profile, 0: nd_layer)
!           Total Planckian flux (this has the full dimensions here
!           as it is used as a dummy array below)
     &  , planck_radiance(nd_radiance_profile, nd_viewing_level)
!           Total Planckian radiance
!
      INTEGER
     &    i_rad_layer(nd_viewing_level)
!           Layers in which radiances are calculated
      REAL  (RealK) ::
     &    frac_rad_layer(nd_viewing_level)
!           Fractions below the tops of the layers
!
      REAL  (RealK) ::
     &    dummy(nd_profile, nd_layer)
!           Dummy real variable
!
!     Machine tolerances
      REAL  (RealK) ::
     &    eps_r
!           Smallest number such that 1.0-EPS_R is not 1.0 to the
!           machine's precision
!       
!
!     Upper boundary conditions:
      REAL  (RealK) ::
     &    flux_inc_down(nd_profile)
!           Incident total flux
     &  , flux_inc_direct(nd_profile)
!           Incident direct flux
!
!     Radiances and fluxes:
      REAL  (RealK) ::
     &    i_direct(nd_radiance_profile, 0: nd_layer)
!           Direct radiance (the first row contains the incident
!           solar radiance: the other rows are calculated)
      REAL  (RealK) ::
     &    flux_total(nd_flux_profile, 2*nd_layer+2)
!           Total fluxes
!
!
!     Subroutines called:
      EXTERNAL planck_ss_source
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      eps_r=epsilon(zen_0(1))
!
!     Preliminary calculations for spherical harmonics (we need to do
!     these tasks before rescaling which is under a common IF-test:
!     other operations can wait).
      IF (i_angular_integration == IP_spherical_harmonic) THEN
!
!       Set limits on ranges of harmonics and set pointers to arrays.
        CALL set_truncation(ierr
     &    , i_truncation, ls_global_trunc
     &    , ls_max_order, ls_local_trunc
     &    , ms_min, ms_max, ms_trunc
     &    , ia_sph_mm, n_order_phase
     &    , nd_max_order
     &    )
!
!       Determine whether special treatment of the solar
!       beam is required.
        l_solar_phf=(isolir == IP_solar).AND.
     &              (i_sph_algorithm == IP_sph_reduced_iter)
        l_rescale_solar_phf=l_rescale.AND.l_solar_phf
!       Calculate the solar scattering angles if treating the
!       solar beam separately.
        IF (l_solar_phf) THEN
          CALL sol_scat_cos(n_profile, n_direction
     &      , zen_0, direction, cos_sol_view
     &      , nd_profile, nd_direction)
!
          DO id=1, n_direction
!
            DO l=1, n_profile
!             Initialize the Legendre polynomials at the zeroth and
!             first orders.
              p_legendre_ls_m1(l)=1.0_RealK
              p_legendre_ls(l)=cos_sol_view(l, id)
            ENDDO
            DO i=1, n_layer
              DO l=1, n_profile
                phase_fnc_solar(l, i, id)=1.0_RealK
     &            +3.0_RealK*p_legendre_ls(l)*phase_fnc(l, i, 1)
              ENDDO
            ENDDO
!
            DO ls=2, n_order_phase_solar
!             Calculate higher orders by recurrences.
              cnst1=1.0_RealK-1.0_RealK/real(ls, RealK)
              DO l=1, n_profile
                p_legendre_tmp(l)=p_legendre_ls(l)
                p_legendre_ls(l)
     &            =(1.0_RealK+cnst1)*p_legendre_ls(l)
     &            *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
                p_legendre_ls_m1(l)=p_legendre_tmp(l)
              ENDDO
              IF (l_henyey_greenstein_pf) THEN
                DO i=1, n_layer
                  DO l=1, n_profile
                    phase_fnc_solar(l, i, id)=phase_fnc_solar(l, i, id)
     &                +real(2*ls+1, RealK)*phase_fnc(l, i, 1)**ls
     &                *p_legendre_ls(l)
                  ENDDO
                ENDDO
              ELSE
                DO i=1, n_layer
                  DO l=1, n_profile
                    phase_fnc_solar(l, i, id)=phase_fnc_solar(l, i, id)
     &                +real(2*ls+1, RealK)*phase_fnc(l, i, ls)
     &                *p_legendre_ls(l)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
!
!
!       Check the information on the phase function.
        IF (l_henyey_greenstein_pf) THEN
          n_order_required=1
        ELSE
          IF (l_rescale) THEN
            n_order_required=max(n_order_phase, n_order_forward)
          ELSE
            n_order_required=n_order_phase
          ENDIF
          IF (l_solar_phf) THEN
            n_order_required=max(n_order_phase, n_order_phase_solar)
          ENDIF
        ENDIF
        IF (n_order_required > n_phase_term) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: Too many terms are required '
     &      //'from the phase function.'
          ierr=i_err_fatal
          RETURN
        ENDIF
!
!       Determine which layers will be required to give radiances.
        CALL set_rad_layer(ierr
     &    , n_layer, n_viewing_level, viewing_level
     &    , i_rad_layer, frac_rad_layer
     &    , nd_viewing_level
     &    )
        if(ierr /= i_normal) RETURN
!
      ELSE IF (i_angular_integration == IP_two_stream) THEN
!
        n_order_phase=1
!
      ENDIF
!
      IF ( (i_angular_integration == IP_two_stream).OR.
     &     (i_angular_integration == IP_spherical_harmonic) ) THEN
!
!       Override the phase function supplied if Henyey-Greenstein
!       phase functions are used.
        IF (l_henyey_greenstein_pf) THEN
          DO i=1, n_layer
            DO ls=1, n_order_phase
              DO l=1, n_profile
                phase_fnc(l, i, ls)=phase_fnc(l, i, 1)**ls
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!

        IF (control%i_direct_tau == ip_direct_noscaling) THEN
          DO i=1, n_layer
            DO l=1, n_profile
              tau_noscal(l, i) = tau(l, i)
            END DO
          END DO
        END IF

        IF (l_rescale) THEN
!
!         Set the forward scattering fraction if it is not explicitly
!         provided.
          IF (.NOT.l_forward) THEN
!
            IF (l_henyey_greenstein_pf) THEN
              DO i=1, n_layer
                DO l=1, n_profile
                  forward_scatter(l, i)
     &              =phase_fnc(l, i, 1)**n_order_forward
                ENDDO
              ENDDO
            ELSE
              DO i=1, n_layer
                DO l=1, n_profile
                  forward_scatter(l, i)
     &              =phase_fnc(l, i, n_order_forward)
                ENDDO
              ENDDO
            ENDIF
!
          ENDIF
!
!         Repeat for the solar beam.
!
          IF (l_rescale_solar_phf) THEN
!
!           If explicit data are not provided they must be calculated.
            IF (.NOT.l_forward_solar) THEN
!
              IF (l_henyey_greenstein_pf) THEN
                DO i=1, n_layer
                  DO l=1, n_profile
                    forward_solar(l, i)
     &                =phase_fnc(l, i, 1)**(n_order_phase_solar+1)
                  ENDDO
                ENDDO
              ELSE
                DO i=1, n_layer
                  DO l=1, n_profile
                    forward_solar(l, i)
     &                =phase_fnc(l, i, n_order_phase_solar+1)
                  ENDDO
                ENDDO
              ENDIF
!
            ENDIF
!
          ENDIF
!
          CALL rescale_phase_fnc(n_profile, 1, n_layer
     &      , n_direction, cos_sol_view
     &      , n_order_phase
     &      , phase_fnc, forward_scatter, forward_solar
     &      , l_rescale_solar_phf, n_order_phase_solar, phase_fnc_solar
     &      , nd_profile, nd_radiance_profile, nd_layer, 1
     &      , nd_direction, nd_max_order
     &      )
!
          CALL rescale_tau_omega(n_profile, 1, n_layer
     &      , tau, omega, forward_scatter
     &      , nd_profile, nd_layer, 1
     &      )
!
        ENDIF
      ENDIF
!
!
!     Limit TAU and OMEGA to avoid arithmetic exceptions.
      DO i=1, n_layer
        DO l=1, n_profile
          tau(l, i)=max(tau(l, i), sqrt(tiny(tau(1,1))))
          omega(l, i)=min(omega(l, i), 1.0_RealK-eps_r)
        ENDDO
      ENDDO
      IF (control%i_direct_tau == ip_direct_noscaling) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            tau_noscal(l, i)=                                            
     &        max(tau_noscal(l, i), sqrt(tiny(tau_noscal(1,1))))
          END DO
        END DO
      END IF
!
!     Calculate Planckian terms:
      IF (isolir == IP_infra_red) THEN
        CALL planck_ss_source(ierr, n_profile, n_layer
     &    , i_angular_integration
     &    , t_level, l_ir_source_quad, t, t_ground
     &    , wavelength_short, wavelength_long
     &    , n_viewing_level, i_rad_layer, frac_rad_layer
     &    , diff_planck, diff_planck_2, d_planck_flux_surface
     &    , planck_flux, planck_radiance
     &    , nd_profile, nd_layer
     &    , nd_radiance_profile, nd_viewing_level
     &    )
        IF (ierr /= i_normal) STOP
      ENDIF
!
!
!
      IF (i_angular_integration == IP_two_stream) THEN
!
!       Check that the solver is supported in this application.
        IF ( (i_solver /= IP_solver_pentadiagonal).AND.
     &       (i_solver /= IP_solver_homogen_direct) ) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: The two-stream solver selected is not '
     &      //'valid here.'
          ierr=i_err_fatal
          RETURN
        ENDIF
!
        IF (isolir == IP_solar) THEN
          DO l=1, n_profile
            flux_inc_direct(l)=solar_irrad(l)/zen_0(l)
            flux_inc_down(l)=solar_irrad(l)/zen_0(l)
          ENDDO
        ELSE
          DO l=1, n_profile
            flux_inc_direct(l)=0.0_RealK
            flux_inc_down(l)=-planck_flux(l, 0)
          ENDDO
        ENDIF
!
        IF (l_iso_inc) THEN
          DO l=1, n_profile
            flux_inc_down(l)=flux_inc_down(l)+pi*iso_inc(l)
          ENDDO
        ENDIF
!
        CALL two_stream(ierr
     &    , control, bound
!                       Atmospheric properties
     &    , n_profile, n_layer
!                       Two-stream scheme
     &    , i_2stream
!                       Options for solver
     &    , i_solver, control%i_scatter_method
!                       Options for equivalent extinction
     &    , .false., dummy
!                       Spectral region
     &    , isolir
!                       Infra-red properties
     &    , diff_planck
     &    , l_ir_source_quad, diff_planck_2
!                       Conditions at TOA
     &    , flux_inc_down, flux_inc_direct, zen_0
!                       Surface conditions
     &    , rho_alb(1, IP_surf_alb_diff)
     &    , rho_alb(1, IP_surf_alb_dir), d_planck_flux_surface
!                       Spherical geometry
     &    , sph
!                       Single scattering properties
     &    , tau_noscal, tau, omega, phase_fnc(1, 1, 1)
!                       Fluxes calculated
     &    , flux_direct, flux_total
!                       Sizes of arrays
     &    , nd_profile, nd_layer, nd_source_coeff
     &    )
!
      ELSE IF (i_angular_integration == IP_spherical_harmonic) THEN
!
!       Calculate Clebsch-Gordan coefficients once and for all.
        CALL calc_cg_coeff(ls_max_order
     &    , ia_sph_mm, ms_min, ms_trunc
     &    , cg_coeff
     &    , nd_max_order, nd_sph_coeff)
!
!       Calculate spherical harmonics at polar angles of pi/2 for
!       use in Marshak's boundary conditions.
        CALL calc_uplm_zero(ms_min, ms_max, ia_sph_mm
     &    , ls_local_trunc, uplm_zero
     &    , nd_max_order, nd_sph_coeff)
!
!
        IF (isolir == IP_solar) THEN
!         Calculate the spherical harmonics
          CALL calc_uplm_sol(n_profile, ms_min, ms_max, ia_sph_mm
     &      , ls_local_trunc, zen_0, uplm_sol
     &      , nd_profile, nd_max_order, nd_sph_coeff)
!         Set the radiances at the top of the atmosphere.
          DO l=1, n_profile
            i_direct(l, 0)=solar_irrad(l)
          ENDDO
        ENDIF
!
!       Precalcuate arrays of terms for the BRDF.
        IF (i_sph_algorithm == IP_sph_reduced_iter) THEN
!         Precalcuate arrays of terms for the BRDF.
          CALL calc_brdf(isolir, ms_min, ms_max, ia_sph_mm
     &      , uplm_sol, uplm_zero
     &      , n_brdf_basis_fnc, ls_brdf_trunc, f_brdf
     &      , n_profile, n_direction, direction
     &      , brdf_sol, brdf_hemi
     &      , nd_profile, nd_radiance_profile, nd_direction
     &      , nd_max_order, nd_sph_coeff
     &      , nd_brdf_basis_fnc, nd_brdf_trunc)
        ENDIF
!
!
!
!
!       Set the incident diffuse flux.
        IF (isolir == IP_solar) THEN
          DO l=1, n_profile
            flux_inc_down(l)=0.0_RealK
          ENDDO
        ELSE
          DO l=1, n_profile
            flux_inc_down(l)=-planck_flux(l, 0)
          ENDDO
        ENDIF
!
        IF (l_iso_inc) THEN
          DO l=1, n_profile
            flux_inc_down(l)=flux_inc_down(l)+pi*iso_inc(l)
          ENDDO
        ENDIF
!
!       Precalculate dimensions for the dynamically allocated
!       arrays.
        nd_red_eigensystem=(nd_max_order+1)/2
        nd_sph_equation=2*nd_layer*nd_red_eigensystem
        nd_sph_diagonal=6*nd_red_eigensystem
        IF (i_sph_algorithm == IP_sph_direct) THEN
          nd_sph_cf_weight=nd_max_order+1
          nd_sph_u_range=2*nd_red_eigensystem
        ELSE IF (i_sph_algorithm == IP_sph_reduced_iter) THEN
          nd_sph_cf_weight=nd_direction
          nd_sph_u_range=nd_sph_equation
        ENDIF
!
        CALL sph_solver(     
!                       Atmospheric sizes
     &      n_profile, n_layer
!                       Angular integration
     &    , ms_min, ms_max, i_truncation, ls_local_trunc
     &    , cg_coeff, uplm_zero, ia_sph_mm
     &    , accuracy_adaptive, euler_factor
     &    , i_sph_algorithm, i_sph_mode
!                       Spectral Region
     &    , isolir
!                       Options for Equivalent Extinction
     &    , .false., dummy
!                       Solar Fields
     &    , i_direct, zen_0, uplm_sol
!                       Infra-red Properties
     &    , diff_planck, flux_inc_down
     &    , l_ir_source_quad, diff_planck_2
!                       Optical properies
     &    , tau, omega, phase_fnc, phase_fnc_solar
!                       Surface Conditions
     &    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb
     &    , f_brdf, brdf_sol, brdf_hemi
     &    , d_planck_flux_surface
!                       Levels for calculating radiances
     &    , n_viewing_level, i_rad_layer, frac_rad_layer
!                       Viewing Geometry
     &    , n_direction, direction
!                       Calculated radiances or fluxes
     &    , flux_direct, flux_total, radiance, photolysis
!                       Dimensions of arrays
     &    , nd_profile, nd_layer
     &    , nd_flux_profile, nd_radiance_profile, nd_j_profile
     &    , nd_max_order, nd_sph_coeff
     &    , nd_brdf_basis_fnc, nd_brdf_trunc
     &    , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal
     &    , nd_sph_cf_weight, nd_sph_u_range
     &    , nd_viewing_level, nd_direction
     &    )
!
      ELSE IF (i_angular_integration == IP_ir_gauss) THEN
!
!       Full angular resolution using Gaussian integration.
        CALL gauss_angle(n_profile, n_layer
     &    , n_order_gauss
     &    , tau
     &    , flux_inc_down
     &    , diff_planck, d_planck_flux_surface
     &    , rho_alb(1, IP_surf_alb_diff)
     &    , flux_total
     &    , nd_profile, nd_layer
     &    )
        IF (ierr /= i_normal) RETURN
!
      ENDIF
!
      IF ( (i_angular_integration == IP_two_stream).OR.
     &     (i_angular_integration == IP_ir_gauss).OR.
     &   ( (i_angular_integration == IP_spherical_harmonic).AND.
     &     (i_sph_mode == IP_sph_mode_flux) ) ) THEN
!       Copy the total fluxes across to the output arrays.
        DO i=0, n_layer
          DO l=1, n_profile
            flux_up(l, i)=flux_total(l, 2*i+1)
            flux_down(l, i)=flux_total(l, 2*i+2)
          ENDDO
        ENDDO
      ENDIF
!
!     Make any adjustments to fluxes and radiances to convert
!     to actual values.
      IF (isolir == IP_infra_red) THEN
        IF ( (i_angular_integration == ip_two_stream).OR.
     &       (i_angular_integration == ip_ir_gauss) ) THEN
          DO i=0, n_layer
            DO l=1, n_profile
              flux_up(l, i)=flux_up(l, i)+planck_flux(l, i)
              flux_down(l, i)=flux_down(l, i)+planck_flux(l, i)
            END DO
          END DO
        ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN
!         Planckian radiances are always used with spherical harmonics,
!         even when calculating fluxes. The number of levels should
!         be set appropriately above.
          IF (i_sph_mode == ip_sph_mode_flux) THEN
            DO i=0, n_layer
              DO l=1, n_profile
                flux_up(l, i)=flux_up(l, i)+pi*planck_radiance(l, i+1)
                flux_down(l, i)=flux_down(l, i)
     &            +pi*planck_radiance(l, i+1)
              END DO
            END DO
          ELSE IF (i_sph_mode == ip_sph_mode_rad) THEN
            DO id=1, n_direction
              DO i=1, n_viewing_level
                DO l=1, n_profile
                  radiance(l, i, id)=radiance(l, i, id)
     &              +planck_radiance(l, i)
                END DO
              END DO
            END DO
          END IF
        END IF
      END IF

      END
