! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve for radiances in harmonics.
!
! Method:
!   After setting the basic properties for the radiance solver
!   a matrix is built and solved.
!
!- ---------------------------------------------------------------------
MODULE sph_solver_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SPH_SOLVER_MOD'
CONTAINS
SUBROUTINE sph_solver(                                                  &
!                 Atmospheric sizes
      n_profile, n_layer                                                &
!                 Angular integration
    , ms_min, ms_max, i_truncation, ls_local_trunc                      &
    , cg_coeff, uplm_zero, ia_sph_mm                                    &
    , accuracy_adaptive, euler_factor                                   &
    , i_sph_algorithm, i_sph_mode                                       &
!                 Spectral Region
    , isolir                                                            &
!                   Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Solar Fields
    , i_direct, mu_0, uplm_sol                                          &
!                 Infra-red Properties
    , diff_planck, flux_inc_down                                        &
    , l_ir_source_quad, diff_planck_2                                   &
!                   Optical properies
    , tau, omega, phase_fnc, phase_fnc_solar                            &
!                 Surface Conditions
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
    , d_planck_flux_surface                                             &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                   Viewing Geometry
    , n_direction, direction                                            &
!                 Calculated Radiances or Fluxes
    , flux_direct, flux_total, radiance_mono, photolysis                &
!                 Dimensions of arrays
    , nd_profile, nd_layer                                              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc                                  &
    , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal              &
    , nd_sph_cf_weight, nd_sph_u_range                                  &
    , nd_viewing_level, nd_direction                                    &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_solar, ip_sph_direct, ip_sph_mode_flux,         &
                     ip_sph_mode_rad, ip_sph_reduced_iter,              &
                     ip_trunc_adaptive
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE build_sph_matrix_mod, ONLY: build_sph_matrix
  USE cg_kappa_ms_mod, ONLY: cg_kappa_ms
  USE eval_uplm_mod, ONLY: eval_uplm
  USE hemi_sph_integ_mod, ONLY: hemi_sph_integ
  USE increment_rad_cf_mod, ONLY: increment_rad_cf
  USE single_scat_sol_mod, ONLY: single_scat_sol
  USE sph_matrix_solver_mod, ONLY: sph_matrix_solver

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_flux_profile                                                   &
!       Size allocated for profiles where fluxes are calculated
    , nd_radiance_profile                                               &
!       Size allocated for profiles where radiances are calculated
    , nd_j_profile                                                      &
!       Size allocated for profiles where mean radiances
!       are calculated
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_direction                                                      &
!       Size allowed for viewing directions
    , nd_viewing_level                                                  &
!       Size allowed for levels where the radiance is calculated
    , nd_max_order                                                      &
!       Size allowed for orders of spherical harmonics
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for orders of BRDFs
    , nd_sph_coeff                                                      &
!       Size allowed for spherical harmonic coefficients
    , nd_red_eigensystem                                                &
!       Size allowed for the reduced eigensystem
    , nd_sph_equation                                                   &
!       Size allowed for spherical harmonic equations
    , nd_sph_diagonal                                                   &
!       Size allowed for diagonals of the spherical harmonic
!       matrix
    , nd_sph_cf_weight                                                  &
!       Size allowed for application of weights of the C. F.
    , nd_sph_u_range
!       Size allowed for the range of u^+|- contributing on any
!       viewing level


  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , isolir
!       Spectral region
!                 Angular integration
  INTEGER, INTENT(IN) ::                                                &
      ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                   &
!       Orders of truncation at each azimuthal order
    , i_truncation                                                      &
!       Type of truncation used
    , i_sph_mode                                                        &
!       Mode in which the spherical solver is used
    , i_sph_algorithm
!       Algorithm used to solve the spherical system
  REAL (RealK), INTENT(IN) ::                                           &
      cg_coeff(nd_sph_coeff)                                            &
!       Clebsch-Gordan coefficients
    , uplm_zero(nd_sph_coeff)                                           &
!       Values of spherical harmonics at polar angles pi/2
    , uplm_sol(nd_profile, nd_sph_coeff)                                &
!       Values of spherical harmonics in the solar direction
    , accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last polar order
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Use quadratic source term

!                   Variables for equivalent extinction
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Apply scaling to solar flux
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment of solar beam with equivalent extinction

  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)                                         &
!       Optical depth
    , omega(nd_profile, nd_layer)                                       &
!       Albedos of single scattering
    , phase_fnc(nd_profile, nd_layer, nd_max_order)                     &
!       Moments of the phase function
    , phase_fnc_solar(nd_radiance_profile, nd_layer, nd_direction)
!       The phase function evaluated for scattering from
!       the solar beam into the viewing directions

  REAL(RealK), INTENT(IN) ::                                            &
      mu_0(nd_profile)
!       Cosines of solar zenith angles

  REAL(RealK), INTENT(IN) ::                                            &
      diff_planck(nd_profile, nd_layer)                                 &
!       Difference in the FLUX Planckian function
    , diff_planck_2(nd_profile, nd_layer)                               &
!       2x2nd differences of Planckian
    , flux_inc_down(nd_profile)
!       Incident downward flux (in real calculations this is used
!       only in the IR where it is assumed to be Planckian, but
!       in may be used in the solar in idealized test cases)
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_profile, 0: nd_layer)
!       Direct solar radiance (only the first row is set on input)
  INTEGER, INTENT(IN) ::                                                &
      ls_brdf_trunc                                                     &
!       Order of trunation of BRDFs
    , n_brdf_basis_fnc
!       Number of BRDF basis functions
  REAL (RealK), INTENT(IN) ::                                           &
      d_planck_flux_surface(nd_profile)                                 &
!       Differential Planckian flux from the surface
    , rho_alb(nd_profile, nd_brdf_basis_fnc)                            &
!       Weights of the basis functions
    , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                         &
!       Array of BRDF basis terms
    , brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)             &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction

!                   Viewing Geometry
  INTEGER, INTENT(IN) ::                                                &
      n_direction                                                       &
!       Number of viewing directions
    , n_viewing_level                                                   &
!       Number of levels where radiances are calculated
    , i_rad_layer(nd_viewing_level)
!       Layers in which radiances are calculated
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)                   &
!       Viewing directions
    , frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers

!                   Radiances calculated
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct Flux
    , flux_total(nd_flux_profile, 2*nd_layer+2)                         &
!       Total Fluxes
    , radiance_mono(nd_radiance_profile                                 &
        , nd_viewing_level, nd_direction)                               &
!       Radiances
    , photolysis(nd_j_profile, nd_viewing_level)
!       Rates of photolysis


! Local variabales.
  INTEGER                                                               &
      ms                                                                &
!       Azimuthal order
    , lsr                                                               &
!       Reduced polar order
    , n_red_eigensystem                                                 &
!       Size of the reduced eigensystem
    , n_equation                                                        &
!       Number of equations
    , i                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable (directions)
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      kappa(nd_max_order/2, nd_max_order/2)                             &
!       Integrals of pairs of spherical harmonics over the downward
!       hemisphere
    , cgk(nd_brdf_trunc/2+1, nd_max_order)                              &
!       Products of the Clebsch-Gordan coefficients and the
!       hemispheric integrals
    , up_lm(nd_profile, nd_max_order+1, nd_direction)                   &
!       Polar parts of spherical harmonics
    , weight_u(nd_profile, nd_viewing_level                             &
        , nd_sph_cf_weight, nd_sph_u_range)                             &
!       Weights to be applied to the vector U containing the
!       complementary functions
    , c_ylm(nd_profile, nd_viewing_level, nd_max_order+1)               &
!       Coefficients in the expansion of the radiance
!       in spherical harmonics
    , a(nd_profile, nd_sph_equation, nd_sph_diagonal)                   &
!       Matrix on the LHS of the equation for spherical
!       harmonics
    , b(nd_profile, nd_sph_equation)                                    &
!       RHS of matrix equation
    , upm(nd_profile, nd_sph_equation)                                  &
!       Variables u+|-
    , azim_factor(nd_profile, nd_direction)
!       Azimuthal factors
  INTEGER                                                               &
      ls_trunc_calc                                                     &
!       Order of truncation required in calculations
    , ls_brdf_trunc_calc                                                &
!       Order of truncation of BRDFs required in calculations
    , ls_significant
!       Maximum significant polar order

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SPH_SOLVER'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the direct radiances which are independent of the
! azimuthal order.
  IF (isolir == ip_solar) THEN

    IF (l_scale_solar) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          i_direct(l, i)=i_direct(l, i-1)                               &
            *adjust_solar_ke(l, i)*EXP(-tau(l, i)/mu_0(l))
        END DO
      END DO
    ELSE
      DO i=1, n_layer
        DO l=1, n_profile
          i_direct(l, i)=i_direct(l, i-1)*EXP(-tau(l, i)/mu_0(l))
        END DO
      END DO
    END IF

  END IF

! Initialize the monochromatic radiances if they are required.
  IF (i_sph_mode == ip_sph_mode_rad) THEN
    DO id=1, n_direction
      DO i=1, n_viewing_level
        DO l=1, n_profile
          radiance_mono(l, i, id)=0.0e+00_RealK
        END DO
      END DO
    END DO
  END IF



! For each non-negative azimuthal order a matrix of coefficients
! is built and solved.
  DO ms=ms_min, ms_max

!   Set the order of truncation required for the calculation.
!   If an adaptive truncation is used we may not calculate all
!   the declared polar orders, though space is reserved for them
!   and is initialized (this is necessary as a subsequent call to
!   this routine with different optical properties may require a
!   different number of terms). On entry we have no estimate of the
!   number of terms required, other than the assigned truncation.
    IF (i_truncation == ip_trunc_adaptive) THEN
      IF (ms == ms_min) THEN
        ls_trunc_calc=ls_local_trunc(ms)
      ELSE
!       The polar order must exceed the azimuthal order by at least
!       one, and must give an even number of harmonics.
        ls_trunc_calc=MAX((ms+1), ls_significant)
        ls_trunc_calc=ls_trunc_calc+MOD((ls_trunc_calc+ms+1), 2)
        ls_trunc_calc=MIN(ls_trunc_calc, ls_local_trunc(ms))
      END IF
    ELSE
      ls_trunc_calc=ls_local_trunc(ms)
    END IF

!   Azimuthal factors are required if calculating radiances.
    IF (i_sph_mode == ip_sph_mode_rad) THEN
!     The recalculation of azimuthal factors for each monochromatic
!     radiance is inefficient, but avoids the need for the storage
!     of arrays with different azimuthal orders. The equation of
!     transfer is solved only for positive azimuthal orders. Hence,
!     if the azimuthal order is 0 we add in the term, but if the
!     azimuthal order is non-zero it also represents the negative
!     order and acquires a weighting factor of 2.
      IF (ms == 0) THEN
        DO id=1, n_direction
          DO l=1, n_profile
            azim_factor(l, id)=1.0e+00_RealK
          END DO
        END DO
      ELSE
        DO id=1, n_direction
          DO l=1, n_profile
            azim_factor(l, id)                                          &
            =2.0e+00_RealK*COS(REAL(ms,RealK)*direction(l,id,2))
          END DO
        END DO
      END IF
!     Calculate spherical harmonics in the viewing directions.
!     A judgement about storage has been made here. It is
!     deemed acceptable to store values for a fixed azimuthal
!     order, but deemed that too much storage would be required
!     for all azmuthal orders to be held at once.
      DO id=1, n_direction
        CALL eval_uplm(ms, ls_trunc_calc                                &
          , n_profile, direction(1, id, 1), up_lm(1, 1, id)             &
          , nd_profile)
      END DO
    END IF

!   Calculate integrals of products of spherical harmonics for
!   use in Marshak's boundary conditions. These arrays are
!   recalculated each time to save storage.
    CALL hemi_sph_integ(ls_trunc_calc, ms, uplm_zero(ia_sph_mm(ms))     &
      , kappa                                                           &
      , nd_max_order                                                    &
      )
!   Perform preliminary calculations for the BRDF: we need to check
!   that the order of truncation of BRDFs does not exceed the
!   order of calculation: it must also be even.
    ls_brdf_trunc_calc                                                  &
      =MIN(ls_brdf_trunc, ls_trunc_calc-MOD(ls_trunc_calc, 2))
    CALL cg_kappa_ms(ms, ls_trunc_calc, ls_brdf_trunc_calc              &
      , cg_coeff, kappa                                                 &
      , cgk                                                             &
      , nd_max_order, nd_brdf_trunc                                     &
      )

!   Initialize the spherical harmonics at this azimuthal order.
    IF ( (i_sph_mode == ip_sph_mode_flux).OR.                           &
         (i_sph_algorithm == ip_sph_direct) ) THEN
      DO i=1, n_viewing_level
        DO lsr=1, ls_local_trunc(ms)+1-ms
          DO l=1, n_profile
            c_ylm(l, i, lsr)=0.0e+00_RealK
          END DO
        END DO
      END DO
    END IF

!   In the infra-red region differential quantities are used. The
!   incident flux will then be determined from the Plankian flux
!   (which is used for consistency with two-stream calculations)
!   by dividing by pi: we must also multiply by sqrt(4.pi) to get
!   the weighting of the zeroth spherical harmonic. Other elements
!   of the array were zeroed before.

    n_red_eigensystem=(ls_trunc_calc+1-ms)/2

    CALL build_sph_matrix(i_sph_algorithm, euler_factor                 &
!                   Basic sizes
      , n_profile, n_layer, ls_trunc_calc                               &
      , ms, n_red_eigensystem                                           &
!                   Numerical arrays of spherical terms
      , cg_coeff(ia_sph_mm(ms)), kappa, up_lm                           &
!                   Solar variables
      , isolir, i_direct, mu_0, uplm_sol(1, ia_sph_mm(ms))              &
      , azim_factor                                                     &
!                   Infra-red variables
      , diff_planck, l_ir_source_quad, diff_planck_2                    &
!                   Isotropic incident flux
      , flux_inc_down                                                   &
!                   Optical properties
      , tau, omega, phase_fnc                                           &
!                   Surface Fields
      , ls_brdf_trunc_calc, n_brdf_basis_fnc, rho_alb                   &
      , f_brdf, brdf_sol, brdf_hemi, cgk                                &
      , d_planck_flux_surface                                           &
!                   Levels where radiances are calculated
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!                   Viewing Geometry
      , n_direction, direction(1, 1, 1)                                 &
!                   Output variables
      , a, b, c_ylm, weight_u, radiance_mono                            &
!                   Dimensions
      , nd_profile, nd_radiance_profile                                 &
      , nd_layer, nd_viewing_level, nd_direction                        &
      , nd_max_order, nd_brdf_basis_fnc, nd_brdf_trunc                  &
      , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal            &
      , nd_sph_cf_weight, nd_sph_u_range                                &
      )

!   Apply standard Gaussian elimination to obtain coefficients
!   u^+|-.

    n_equation=2*n_layer*n_red_eigensystem

    CALL sph_matrix_solver(n_profile, n_layer, n_red_eigensystem        &
      , a, b                                                            &
      , upm                                                             &
      , nd_profile, nd_sph_equation, nd_sph_diagonal                    &
      )

!   Increment the radiances with the contributions from
!   the complementary function.
    CALL increment_rad_cf(n_profile                                     &
      , n_direction, azim_factor                                        &
      , n_viewing_level, i_rad_layer                                    &
      , i_sph_mode, i_sph_algorithm                                     &
      , ms, ls_trunc_calc, euler_factor                                 &
      , isolir, mu_0, kappa, up_lm                                      &
      , n_red_eigensystem, n_equation, weight_u, upm                    &
      , i_direct, c_ylm, flux_direct, flux_total                        &
      , radiance_mono, photolysis                                       &
      , nd_profile, nd_flux_profile                                     &
      , nd_radiance_profile, nd_j_profile                               &
      , nd_layer, nd_direction, nd_viewing_level                        &
      , nd_max_order, nd_sph_equation, nd_sph_cf_weight                 &
      , nd_sph_u_range                                                  &
      )

    IF (i_truncation == ip_trunc_adaptive) THEN

!     Reduce the polar order of truncation if higher polar orders
!     make an insignificant contribution to the radiance field.
!     At least two polar orders are required.
      ls_significant=ms+1
      DO lsr=ms+2, ls_trunc_calc+1-ms
        DO i=1, n_viewing_level
          DO l=1, n_profile
            IF (ABS(c_ylm(l, i, lsr)) >  accuracy_adaptive              &
              *ABS(c_ylm(l, i, 1))) ls_significant=lsr+ms-1
          END DO
        END DO
      END DO
    END IF

  END DO

  IF ( (i_sph_algorithm == ip_sph_reduced_iter).AND.                    &
       (isolir == ip_solar) ) THEN
!   Add in the singly scattered solar beam using the
!   potentially higher order of truncation.
    CALL single_scat_sol(n_profile, n_layer                             &
      , n_direction, direction                                          &
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
      , i_direct, mu_0                                                  &
      , tau, omega, phase_fnc_solar                                     &
      , radiance_mono                                                   &
      , nd_profile, nd_radiance_profile                                 &
      , nd_layer, nd_direction, nd_viewing_level                        &
      )
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sph_solver
END MODULE sph_solver_mod
