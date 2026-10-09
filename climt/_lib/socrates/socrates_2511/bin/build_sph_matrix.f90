! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to build up the matrix for radiances.
!
! Purpose:
!   This routine assembles the stepped matrix to solve the equation
!   of transfer for a specific value of the azimuthal quantum
!   number.
!
! Method:
!
!   Labelling of variables and equations:
!     The variables are u_{ik}^{-|+} where i runs over all the layers
!   1,...,N_LAYER and k runs over positive eigenvalues of the reduced
!   eigensystem 1,...,N_RED_EIGENSYSTEM: there are thus 2n_e variables
!   describing the radiance in each layer, so the number of a variable
!   is
!      IV=2n_e(i-1)+k+n_e(1+|-1)/2
!   (Note that u_{ik}^- preceeds u_{ik}^+ by N_RED_EIGENSYSTEM).
!   At the top of the atmosphere (L'+1-m)/2 conditions are applied by
!   Marshak's conditions, where l'=m+1,...,L' in steps of 2, so for
!   this boundary
!      IE=(l'+1-m)/2
!   At the i'th interior boundary a condition of continuity is applied
!   to I_{lm}, where l=m,...,L'. To match the numbering of the equations
!   at the boundary values of l=m, m+2,...,L'-1 in steps of 2 are taken
!   first followed by those with l=m+1,...,L', so the number of the
!   equation is
!      IE=n_e(2i-1)+(l-m)/2+1,           l=m, m+2,...,L'-1
!      IE=n_e(2i-1)+(l+1-m)/2+n_e,       l=m+1,...,L',
!   allowing for n_e conditions at the top of the model and 2n_e
!   conditions at higher interfaces. At the bottom of the atmosphere
!   Marshak's condition is imposed using the harmonics l'=m+1,...,L'
!   in steps of 2, so the numbering of equations is
!      IE=n_e(2N_LAYER-1)+(l'+1-m)/2
!     Each of these equations couples together u_{ik}^{+|-} in the
!   layers above and below the interface. Hence, each equation
!   IE=(2i-1)n_e+1,...,(2i+1)n_e involves the variables IV=(2i-1)n_e+1,
!   ...,(2i+1)n_e, producing a stepped diagonal matrix which can be
!   encoded in an array of 4n_e columns with IE indexing the rows.
!   3n_e-1 sub-diagonals. The mapping is:
!     (IE, IV) --> (IE, IV-2*N_RED_EIGENSYSTEM*(I-1))
!
!- ---------------------------------------------------------------------
MODULE build_sph_matrix_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BUILD_SPH_MATRIX_MOD'
CONTAINS
SUBROUTINE build_sph_matrix(i_sph_algorithm, euler_factor               &
!                 Basic sizes
    , n_profile, n_layer, ls_trunc, ms, n_red_eigensystem               &
!                 Numerical arrays of spherical terms
    , cg_coeff, kappa, up_lm                                            &
!                 Solar variables
    , isolir, i_direct, mu_0, uplm_sol, azim_factor                     &
!                 Infra-red variables
    , diff_planck, l_ir_source_quad, diff_planck_2                      &
!                 Diffuse incident field
    , flux_down_inc                                                     &
!                 Optical properies
    , tau, omega, phase_fnc                                             &
!                 Surface Fields
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi, cgk                                  &
    , d_planck_flux_surface                                             &
!                 Levels where radiances are calculated
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                   Viewing Geometry
    , n_direction, mu_v                                                 &
!                 Output variables
    , a, b, c_ylm, weight_u, radiance                                   &
!                 Dimensions
    , nd_profile, nd_radiance_profile, nd_layer                         &
    , nd_viewing_level, nd_direction                                    &
    , nd_max_order, nd_brdf_basis_fnc, nd_brdf_trunc                    &
    , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal              &
    , nd_sph_cf_weight, nd_sph_u_range                                  &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_solar, ip_infra_red, ip_sph_direct,             &
                     ip_sph_reduced_iter 
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE calc_surf_rad_mod, ONLY: calc_surf_rad
  USE calc_top_rad_mod, ONLY: calc_top_rad
  USE eig_sys_mod, ONLY: eig_sys
  USE layer_part_integ_mod, ONLY: layer_part_integ
  USE set_dirn_weights_mod, ONLY: set_dirn_weights
  USE set_level_weights_mod, ONLY: set_level_weights

  IMPLICIT NONE

! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_radiance_profile                                               &
!       Size allocated for atmospheric profiles where radiances
!       are calculated
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_viewing_level                                                  &
!       Allocated size for levels where radiances are calculated
    , nd_direction                                                      &
!       Allocated size for viewing directions
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_brdf_basis_fnc                                                 &
!       Size allocated for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allocated for orders in BRDFs
    , nd_red_eigensystem                                                &
!       Size allocated for the reduced eigensystem
    , nd_sph_equation                                                   &
!       Size allocated for spherical harmonic equations
    , nd_sph_diagonal                                                   &
!       Size allocated for diagonals in matrix for harmonics
    , nd_sph_cf_weight                                                  &
!       Size allocated for enetities to be incremented by the
!       complementary function
    , nd_sph_u_range
!       Range of values of u^+|- contributing on any viewing
!       level


! Dummy arguments
! Atmospheric structrure:
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric layers
    , n_layer
!       Number of atmospheric layers

! Spherical harmonic structure:
  INTEGER, INTENT(IN) ::                                                &
      i_sph_algorithm                                                   &
!       Algorithm for the spherical harmonic solution
    , ls_trunc                                                          &
!       The truncating order of the system of equations
    , ms                                                                &
!       Azimuthal order
    , n_red_eigensystem
!       Size of the reduced eigensystem
  REAL (RealK), INTENT(IN) ::                                           &
      euler_factor
!       Factor applied to the last term of an alternating series

  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Flag for spectral region

! Optical properties:
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)                                         &
!       Optical depths of the layers
    , omega(nd_profile, nd_layer)                                       &
!       Albedos of single scattering of the layers
    , phase_fnc(nd_profile, nd_layer, nd_max_order)
!       Phase functions of the layers

! Solar Fields:
  REAL (RealK), INTENT(IN) ::                                           &
      mu_0(nd_profile)                                                  &
!       Cosine of solar zenith angle
    , i_direct(nd_profile, 0: nd_layer)                                 &
!       The direct solar radiance
    , uplm_sol(nd_profile, ls_trunc+2-ms)
!       Spherical harmonics of the solar angle

! Infra-red quantities:
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Flag for quadratic source function in the IR
  REAL (RealK), INTENT(IN) ::                                           &
      diff_planck(nd_profile, nd_layer)                                 &
!       Differences in the hemispheric Planckian FLUX (bottom-top)
!       across the layer
    , diff_planck_2(nd_profile, nd_layer)
!       Twice the second differences in the hemispheric Planckian
!       FLUX
  REAL (RealK), INTENT(IN) ::                                           &
      cg_coeff(ls_trunc+1-ms)                                           &
!       Clebsch-Gordan coefficients
    , kappa(nd_max_order/2, nd_max_order/2)                             &
!       Integrals of pairs of spherical harmonics over the downward
!       hemisphere
    , cgk(nd_brdf_trunc/2+1, nd_max_order)                              &
!       Products of the Clebsch-Gordan coefficients and the
!       hemispheric integrals
    , up_lm(nd_profile, nd_max_order+1, nd_direction)                   &
!       Polar parts of spherical harmonics in viewing directions
    , flux_down_inc(nd_profile)
!       Diffuse hemispherically isotropic incident flux

! Surface Fields:
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
      mu_v(nd_profile, nd_direction)                                    &
!       Cosines of polar viewing directions
    , azim_factor(nd_profile, nd_direction)                             &
!       Azimuthal factors
    , frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers


  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_radiance_profile                                      &
        , nd_viewing_level, nd_direction)
!       Radiances to be incremented (note that at lower
!       levels this is declared withh ND_PROFILE, but this
!       is fine since those routines will be called only
!       when the two sizes are equal)

  REAL (RealK), INTENT(OUT) ::                                          &
      a(nd_profile, nd_sph_equation, nd_sph_diagonal)                   &
!       Matrix in the LHS of the equations
    , b(nd_profile, nd_sph_equation)                                    &
!       Vector of forcings for the matrix equation
    , weight_u(nd_profile, nd_viewing_level, nd_sph_cf_weight           &
        , nd_sph_u_range)
!       Weights to be applied to the vector U containing the
!       complementary functions

!                 Radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      c_ylm(nd_profile, nd_viewing_level, ls_trunc+1-ms)
!       Coefficients for radiances


! Local variables
  INTEGER                                                               &
      ie                                                                &
!       Number of the equation
    , ivma                                                              &
!       Index for u^- in the layer above the interface
!       This and the next three variables are used in two forms:
!       with an offset in indexing the matrix A and without
!       an offset in indexing the array WEIGHT_U.
    , ivmb                                                              &
!       Index for u^- in the layer below the interface
    , ivpa                                                              &
!       Index for u^+ in the layer above the interface
    , ivpb                                                              &
!       Index for u^+ in the layer below the interface
    , i_above                                                           &
!       Index for layer above in two-dimensional arrays
    , i_below                                                           &
!       Index for layer below in two-dimensional arrays
    , i_assign_level
!       Level where a radiance is to be assigned
  INTEGER                                                               &
      ls_p                                                              &
!       Primed polar order
    , ls                                                                &
!       Polar order
    , lsr_p                                                             &
!       Reduced primed polar order (LSR_P is MS-1 less than
!       LS_P to facilitate addressing of arrays which do not
!       hold redundant space for m>l')
    , lsr                                                               &
!       Reduced polar order
    , ls_d                                                              &
!       Dummy polar order
    , lsr_d                                                             &
!       Reduced dummy polar order
    , ls_dd                                                             &
!       Dummy polar order
    , lsr_dd                                                            &
!       Reduced dummy polar order
    , i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l
!       Loop variable
  LOGICAL                                                               &
      l_assign
!       Controlling logical for assigning levels
  REAL (RealK) ::                                                       &
      ss(nd_profile, 0: nd_max_order)                                   &
!       S-coefficients for the current layer
    , ssrt(nd_profile, 0: nd_max_order)
!       Square roots of S-coefficients
  REAL (RealK) ::                                                       &
      mu(nd_profile, nd_red_eigensystem, 2)                             &
!       Eigenvaluse of the reduced system
    , eig_vec(nd_profile, 2*nd_red_eigensystem                          &
        , nd_red_eigensystem, 2)                                        &
!       Eigenvectors of the full systems for positive eigenvalues
!       (these are scaled by the s-coefficients in the routine
!       EIG_SYS)
    , theta(nd_profile, nd_red_eigensystem, 2)                          &
!       Array of exponentials of optical depths along slant paths
    , source_top(nd_profile, ls_trunc+1-ms, 2)                          &
!       Source function at the top of the layer
    , source_bottom(nd_profile, ls_trunc+1-ms, 2)
!       Source function at the bottom of the layer
  REAL (RealK) ::                                                       &
      surface_term(nd_profile, ls_trunc+1-ms)                           &
!       Surface terms involving BRDFs
    , b_factor(nd_profile)                                              &
!       Contribution to the RHS of the equations
    , ksi                                                               &
!       Expression involving the BRDF
    , phi                                                               &
!       Expression involving the BRDF
    , phi_d                                                             &
!       Expression involving the BRDF
    , lambda                                                            &
!       Expression involving the BRDF
    , lambda_d
!       Expression involving the BRDF
  REAL (RealK) ::                                                       &
      z_sol(nd_profile, ls_trunc+1-ms)
!       Coefficient of the solar source function at the top of
!       the layer
  REAL (RealK) ::                                                       &
      q_0(nd_profile)                                                   &
!       Term for thermal particular integral
    , q_1(nd_profile)
!       Term for thermal particular integral
  INTEGER                                                               &
      k_sol(nd_profile)
!       Index of eigenvalue closest to the cosine of the solar
!       zenith angle
  REAL (RealK) ::                                                       &
      upm_c(nd_profile, 2*nd_red_eigensystem)
!       Weights for exponentials in conditioning term

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='BUILD_SPH_MATRIX'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialize the matrix.
  DO ie=1, 2*n_layer*n_red_eigensystem
    DO k=1, 6*n_red_eigensystem
      DO l=1, n_profile
        a(l, ie, k)=0.0e+00_RealK
      END DO
    END DO
  END DO

! To keep track of the layers in which radiances are required
! I_ASSIGN_LEVEL is used: we search for the layer containing this
! level, as indicated by I_RAD_LAYER and set the elements of
! WEIGHT_U for later use with the vector giving the complementary
! function. The terms of the particular integral are assigned to
! C_YLM. Initialize to look for the first level.
  i_assign_level=1

! I_BELOW and I_ABOVE hold variables for the layers below and
! above the current interface. They are flipped to enable us
! to use arrays with a dimension of 2, without the need to copy
! lots of data.
  i_below=1



! Begin by determining the properties of the top layer.
  IF (ms == 0) THEN
    DO l=1, n_profile
      ss(l, ms)=1.0e+00_RealK-omega(l, 1)
      ssrt(l, ms)=SQRT(ss(l, ms))
    END DO
  END IF
  DO ls=MAX(1, ms), ls_trunc
    DO l=1, n_profile
      ss(l, ls)=1.0e+00_RealK-omega(l, 1)*phase_fnc(l, 1, ls)
      ssrt(l, ls)=SQRT(ss(l, ls))
    END DO
  END DO

! Calculate the eigenvalues and eigenvectors for this layer.
  CALL eig_sys(n_profile, ls_trunc, ms, n_red_eigensystem               &
    , cg_coeff, ssrt(1, 0)                                              &
    , mu(1, 1, i_below), eig_vec(1, 1, 1, i_below)                      &
    , nd_profile, nd_red_eigensystem, nd_max_order                      &
    )

! Calculate the exponential terms for this layer
  DO k=1, n_red_eigensystem
    DO l=1, n_profile
      theta(l, k, i_below)=EXP(-tau(l, 1)/mu(l, k, i_below))
    END DO
  END DO

! Find the particular integral in this layer.
  CALL layer_part_integ(                                                &
      n_profile, ls_trunc, ms, n_red_eigensystem                        &
    , cg_coeff, mu(1, 1, i_below)                                       &
    , eig_vec(1, 1, 1, i_below), theta(1, 1, i_below)                   &
    , isolir, i_direct(1, 0), mu_0, uplm_sol                            &
    , diff_planck(1, 1), l_ir_source_quad, diff_planck_2(1, 1)          &
    , tau(1, 1), ss(1, 0)                                               &
    , source_top(1, 1, i_below), source_bottom(1, 1, i_below)           &
    , upm_c, k_sol, z_sol, q_0, q_1                                     &
    , nd_profile, nd_max_order, nd_red_eigensystem                      &
    )


! Impose Marshak's boundary conditions at the top of the atmosphere.
! For each allowed order of l' (LS_P, or LSR_P in the reduced
! notation), those with odd parity, the integral of Y_l'^m and the
! boundary condition on the radiance is formed and integrated over
! the downward hemisphere.

  DO lsr_p=2, ls_trunc+1-ms, 2

    ie=lsr_p/2

!   Begin with the exceptional case in which l=l' and KAPPA is 1/2.
    DO l=1, n_profile
      b(l, ie)=-0.5e+00_RealK*source_top(l, lsr_p, i_below)
    END DO
!   For other values of l, which must be odd when l' is even and
!   vice versa, the precalculated values are used. A hemispherically
!   isotropic incident radiance may exist if l=m=0, so we this
!   case exceptionally, adjusting the beginning of the loop.
    IF (ms == 0) THEN
      DO l=1, n_profile
        b(l, ie)=b(l, ie)+kappa(lsr_p/2, 1)                             &
          *(2.0e+00_RealK*flux_down_inc(l)/SQRT(pi)                     &
          -source_top(l, 1, i_below))
      END DO
    END IF
    DO lsr=MAX(3-2*ms, 1), ls_trunc-ms, 2
      DO l=1, n_profile
        b(l, ie)=b(l, ie)                                               &
          -kappa(lsr_p/2, (lsr+1)/2)*source_top(l, lsr, i_below)
      END DO
    END DO

!   Now calculate the coefficients of the matrix of unknowns,
!   u_{mik}^{+|-}.
    DO k=1, n_red_eigensystem
!     Variable numbers:
!     To accord with the general structure of the compressed matrix
!     the equations for the top boundary conditions are
!     right-justified by advancing the column by
!     2*N_RED_EIGENSYSTEM.
      ivmb=k+2*n_red_eigensystem
      ivpb=ivmb+n_red_eigensystem
!     In Marshak's procedure, l'+m will be odd, so apart from
!     the term where l'=l, l+m will be even in all the non-zero
!     terms of the sum over l, so it is easy to obtain
!     A(L, IE, IVMB) by a simple subtraction from A(L, IE, IVPB).
!     Begin with the term l=l'.
      DO l=1, n_profile
        a(l, ie, ivpb)=0.5e+00_RealK*eig_vec(l, lsr_p, k, i_below)
      END DO
      DO ls=ms, ls_trunc-1, 2
        lsr=ls+1-ms
        DO l=1, n_profile
          a(l, ie, ivpb)                                                &
            =a(l, ie, ivpb)+kappa(lsr_p/2, (lsr+1)/2)                   &
            *eig_vec(l, lsr, k, i_below)
        END DO
      END DO
      DO l=1, n_profile
        a(l, ie, ivmb)                                                  &
          =a(l, ie, ivpb)-eig_vec(l, lsr_p, k, i_below)
        a(l, ie, ivpb)=a(l, ie, ivpb)*theta(l, k, i_below)
      END DO
    END DO
  END DO

! Set the weightings to be applied to the solution of the
! linear system of equations.
  IF (i_sph_algorithm == ip_sph_direct) THEN
!   If we solve the problem directly the weightings will
!   apply to coefficients of the spherical harmonics.

!   The next test is done in two parts to ensure that it reamains
!   within bounds on I_RAD_LAYER.
    l_assign=(i_assign_level <= n_viewing_level)
    IF (l_assign) l_assign=(i_rad_layer(i_assign_level) == 1)

    CALL set_level_weights(1, n_profile, ls_trunc                       &
      , ms, n_red_eigensystem                                           &
      , cg_coeff, mu(1, 1, i_below), eig_vec(1, 1, 1, i_below)          &
      , isolir, z_sol(1, 1), mu_0                                       &
      , q_0, l_ir_source_quad, q_1                                      &
      , upm_c, k_sol                                                    &
      , tau, ss                                                         &
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
      , l_assign, i_assign_level                                        &
      , c_ylm, weight_u(1, 1, 1, 1)                                     &
      , nd_profile, nd_viewing_level                                    &
      , nd_max_order                                                    &
      , nd_red_eigensystem, nd_sph_cf_weight                            &
      )
  ELSE IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
!   Here the weights couple directly to radiances in
!   particular directions.
    CALL set_dirn_weights(n_profile                                     &
      , ms, ls_trunc, up_lm                                             &
      , n_direction, mu_v, azim_factor                                  &
      , n_viewing_level, i_rad_layer, frac_rad_layer, 1                 &
      , n_red_eigensystem                                               &
      , mu(1, 1, i_below), eig_vec(1, 1, 1, i_below)                    &
      , isolir, z_sol(1, 1), mu_0                                       &
      , l_ir_source_quad, diff_planck                                   &
      , upm_c, k_sol                                                    &
      , tau, omega, phase_fnc                                           &
      , weight_u(1, 1, 1, 1), radiance                                  &
      , nd_profile, nd_layer, nd_direction, nd_viewing_level            &
      , nd_red_eigensystem, nd_max_order                                &
      )
  END IF



! For each interior level, 1,.., N_LAYER-1, continuity is imposed
! on the (LS,MS)th component of the radiance field, for LS in the
! range MS,..., LS_TRUNC (2*N_RED_EIGENSYSTEM orders). At each
! stage we need information about the layers above and below the
! layer. Arrays such as THETA therefore have an extra dimension of
! size 2 to hold both values without the need to declare storage
! for the whole column. To avoid copying values this last `index' is
! accessed using the variables I_ABOVE and I_BELOW which are
! flipped as we pass through each layer. To follow the indexing
! note that when the loop variable is I we are looking at the
! I-1st interface.

  DO i=2, n_layer

!   Flip the indices for the layer above and below.
    i_above=i_below
    i_below=3-i_below

!   Calculate the condensed optical properties of the
!   current layer.
    IF (ms == 0) THEN
      DO l=1, n_profile
        ss(l, ms)=1.0e+00_RealK-omega(l, i)
        ssrt(l, ms)=SQRT(ss(l, ms))
      END DO
    END IF
    DO ls=MAX(1, ms), ls_trunc
      DO l=1, n_profile
        ss(l, ls)=1.0e+00_RealK-omega(l, i)*phase_fnc(l, i, ls)
        ssrt(l, ls)=SQRT(ss(l, ls))
      END DO
    END DO

!   Calculate the eigenvalues and eigenvectors for the current
!   layer which is that below the interface.
    CALL eig_sys(n_profile, ls_trunc, ms, n_red_eigensystem             &
      , cg_coeff, ssrt(1, 0)                                            &
      , mu(1, 1, i_below), eig_vec(1, 1, 1, i_below)                    &
      , nd_profile, nd_red_eigensystem, nd_max_order                    &
      )


!   Calculate the exponential terms for this layer
    DO k=1, n_red_eigensystem
      DO l=1, n_profile
         theta(l, k, i_below)=EXP(-tau(l, i)/mu(l, k, i_below))
      END DO
    END DO

!   Find the particular integral in this layer.
    CALL layer_part_integ(                                              &
        n_profile, ls_trunc, ms, n_red_eigensystem                      &
      , cg_coeff, mu(1, 1, i_below)                                     &
      , eig_vec(1, 1, 1, i_below), theta(1, 1, i_below)                 &
      , isolir, i_direct(1, i-1), mu_0, uplm_sol                        &
      , diff_planck(1, i), l_ir_source_quad, diff_planck_2(1, i)        &
      , tau(1, i), ss(1, 0)                                             &
      , source_top(1, 1, i_below), source_bottom(1, 1, i_below)         &
      , upm_c, k_sol, z_sol, q_0, q_1                                   &
      , nd_profile, nd_max_order, nd_red_eigensystem                    &
      )

!   Loop over the permitted orders of LS, compressing entries
!   into the matrix.
    DO lsr=1, 2*n_red_eigensystem

!     Number the equation:
      IF (MOD(lsr, 2) == 1) THEN
        ie=n_red_eigensystem*(2*i-3)+(lsr+1)/2
      ELSE IF (MOD(lsr, 2) == 0) THEN
        ie=n_red_eigensystem*(2*i-2)+lsr/2
      END IF

!     Loop over eigenvalues.
      DO k=1, n_red_eigensystem
!       Assign number to the variables in the equation
        ivma=k
        ivpa=ivma+n_red_eigensystem
        ivmb=ivpa+n_red_eigensystem
        ivpb=ivmb+n_red_eigensystem
        DO l=1, n_profile
          a(l, ie, ivma)=eig_vec(l, lsr, k, i_above)                    &
            *theta(l, k, i_above)*REAL(1-2*MOD(lsr-1, 2), RealK)
          a(l, ie, ivpa)=eig_vec(l, lsr, k, i_above)
          a(l, ie, ivmb)=-eig_vec(l, lsr, k, i_below)                   &
            *REAL(1-2*MOD(lsr-1, 2), RealK)
          a(l, ie, ivpb)=-eig_vec(l, lsr, k, i_below)                   &
            *theta(l, k, i_below)
        END DO
      END DO

      DO l=1, n_profile
        b(l, ie)=source_top(l, lsr, i_below)                            &
          -source_bottom(l, lsr, i_above)
      END DO

    END DO

    IF (i_sph_algorithm == ip_sph_direct) THEN
!     If we solve the problem directly the weightings will
!     apply to coefficients of the spherical harmonics.
!     An assignment is required only if there are remaining
!     viewing levels and we are in the right layer.

!     The next test is done in two parts to ensure that it reamains
!     within bounds on I_RAD_LAYER.
      l_assign=(i_assign_level <= n_viewing_level)
      IF (l_assign) l_assign=(i_rad_layer(i_assign_level) == i)

!     The different indexing of WEIGHT_U in the following two
!     calls is intentional. In the first case we interpolate
!     the radiance in one layer, so the final index runs only over
!     the eigensystem for that layer. In the second case, there
!     are contributions to the radiance at a particular level
!     from all layers, so the final index must be allow for
!     contributions from all layers.

      CALL set_level_weights(i, n_profile, ls_trunc                     &
        , ms, n_red_eigensystem                                         &
        , cg_coeff, mu(1, 1, i_below), eig_vec(1, 1, 1, i_below)        &
        , isolir, z_sol(1, 1), mu_0                                     &
        , q_0, l_ir_source_quad, q_1                                    &
        , upm_c, k_sol                                                  &
        , tau(1, i), ss                                                 &
        , n_viewing_level, i_rad_layer, frac_rad_layer                  &
        , l_assign, i_assign_level                                      &
        , c_ylm, weight_u(1, 1, 1, 1)                                   &
        , nd_profile, nd_viewing_level                                  &
        , nd_max_order                                                  &
        , nd_red_eigensystem, nd_sph_cf_weight                          &
        )
    ELSE IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
!     Here the weights couple directly to radiances in
!     particular directions.
      CALL set_dirn_weights(n_profile                                   &
        , ms, ls_trunc, up_lm                                           &
        , n_direction, mu_v, azim_factor                                &
        , n_viewing_level, i_rad_layer, frac_rad_layer, i               &
        , n_red_eigensystem                                             &
        , mu(1, 1, i_below), eig_vec(1, 1, 1, i_below)                  &
        , isolir, z_sol(1, 1), mu_0                                     &
        , l_ir_source_quad, diff_planck                                 &
        , upm_c, k_sol                                                  &
        , tau, omega, phase_fnc                                         &
        , weight_u(1, 1, 1, 1+2*n_red_eigensystem*(i-1)), radiance      &
        , nd_profile, nd_layer, nd_direction, nd_viewing_level          &
        , nd_red_eigensystem, nd_max_order                              &
        )
    END IF


  END DO




! Impose the surface boundary condition using the appropriate
! bidirectional reflection functions with Marshak's conditions.

! Flip the `index' to the layer above the interface.
  i_above=i_below

  DO lsr_p=2, ls_trunc+1-ms, 2

    ls_p=lsr_p+ms-1

    ie=n_red_eigensystem*(2*n_layer-1)+lsr_p/2

!   Initialize the RHS of the equations.
    DO l=1, n_profile
      b(l, ie)=0.0e+00_RealK
    END DO

!   Terms in this equation fall into two groups: those which
!   involve the BRDF and those which do not. These latter
!   terms, which arise directly from Marshak's conditions
!   are treated first.

!   Begin with the exceptional case where l=l' so KAPPA is 1/2.
    DO l=1, n_profile
      b(l,ie)=b(l,ie)-REAL(1-2*MOD(ls_p,2),RealK)*0.5e+00_RealK         &
        *source_bottom(l, lsr_p, i_above)
    END DO
    IF ( (isolir == ip_infra_red).AND.(ms == 0).AND.                    &
         (lsr_p == 1) ) THEN
!     The Planckian flux is used instead of the radiance for
!     consistency with the two-stream equations.
      DO l=1, n_profile
        b(l, ie)=b(l, ie)+d_planck_flux_surface(l)/SQRT(pi)
      END DO
    END IF
!   For other values of l, which must be odd when l' is even and
!   vice versa, the precalculated values are used.
    DO lsr=1, ls_trunc-ms, 2
      DO l=1, n_profile
        b(l, ie)=b(l, ie)                                               &
          -REAL(1-2*MOD(lsr+ms-1, 2), RealK)                            &
          *kappa(lsr_p/2, (lsr+1)/2)                                    &
          *source_bottom(l, lsr, i_above)
      END DO
      IF ( (isolir == ip_infra_red).AND.(ms == 0).AND.                  &
           (lsr == 1) ) THEN
        DO l=1, n_profile
          b(l, ie)=b(l, ie)+kappa(lsr_p/2, (lsr+1)/2)                   &
            *d_planck_flux_surface(l)*2.0e+00_RealK/SQRT(pi)
        END DO
      END IF
    END DO

!   Now calculate the coefficients of the matrix of unknowns,
!   u_{mik}^{+|-}.
    DO k=1, n_red_eigensystem
!     Variable numbers:
      ivma=k
      ivpa=ivma+n_red_eigensystem
!     KAPPA has not been calculated for those values of l for
!     which it is 0: we therefore add the terms in two groups,
!     firstly those for l=l' and then those for other values of
!     l where KAPPA is non-zero. As at the top, of the atmosphere
!     it is possible to evaluate A(L, IE, IVMA) from
!     A(L, IE, IVPA)
      ls_p=lsr_p+ms-1
      DO l=1, n_profile
        a(l, ie, ivpa)=a(l, ie, ivpa)+REAL(1-2*MOD(ls_p, 2), RealK)     &
          *0.5e+00_RealK*eig_vec(l, lsr_p, k, i_above)
      END DO
      DO ls=ms, ls_trunc-1, 2
        lsr=ls+1-ms
        DO l=1, n_profile
          a(l, ie, ivpa)=a(l, ie, ivpa)+REAL(1-2*MOD(ls, 2), RealK)     &
            *kappa(lsr_p/2, (lsr+1)/2)*eig_vec(l, lsr, k, i_above)
        END DO
      END DO

      DO l=1, n_profile
        a(l, ie, ivma)=a(l, ie, ivpa)+REAL(1-2*MOD(ms, 2), RealK)       &
            *eig_vec(l, lsr_p, k, i_above)
      END DO

    END DO


!   The second group of terms involves the BRDF.
!   There will be no contribution from orders
!   above the order of trunction of the BRDF.
    IF (ms <= ls_brdf_trunc) THEN
!     Add in the solar or infra-red contributions involving the
!     BRDF basis functions which do not involve terms in KAPPA.
      IF (isolir == ip_solar) THEN

        DO j=1, n_brdf_basis_fnc
          DO l=1, n_profile
            b_factor(l)=0.0e+00_RealK
          END DO
          DO ls=ms, ls_brdf_trunc-MOD(ms, 2), 2
            lsr=ls+1-ms
            ksi=kappa(lsr_p/2, 1)*f_brdf(j, 0, ls/2, ms)
            DO ls_d=ms+2, ls_brdf_trunc-MOD(ms, 2), 2
              lsr_d=ls_d-ms+1
              ksi=ksi+kappa(lsr_p/2, (lsr_d+1)/2)                       &
                *f_brdf(j, ls_d/2, ls/2, ms)
            END DO
            DO l=1, n_profile
              b_factor(l)=b_factor(l)+ksi*uplm_sol(l, lsr)
            END DO
          END DO
          DO l=1, n_profile
            b(l, ie)=b(l, ie)+i_direct(l, n_layer)*mu_0(l)              &
              *REAL(1-2*MOD(ms, 2), RealK)                              &
              *rho_alb(l, j)*b_factor(l)
          END DO
        END DO

      ELSE IF (isolir == ip_infra_red) THEN
        IF (ms == 0) THEN
          DO j=1, n_brdf_basis_fnc
            lambda=0.0e+00_RealK
            DO ls_d=0, ls_brdf_trunc, 2
              lsr_d=ls_d+1
              lambda_d=0.0e+00_RealK
              DO ls_dd=0, ls_brdf_trunc, 2
                lsr_dd=ls_dd+1
                lambda_d=lambda_d+kappa(lsr_p/2, (lsr_dd+1)/2)          &
                  *f_brdf(j, ls_dd/2, ls_d/2, ms)
              END DO
              lambda=lambda+kappa(1, (lsr_d+1)/2)*lambda_d
            END DO
            DO l=1, n_profile
              b(l, ie)=b(l, ie)                                         &
                +rho_alb(l, j)*lambda                                   &
                *SQRT(4.0e+00_RealK*pi/3.0e+00_RealK)                   &
                *d_planck_flux_surface(l)/pi
            END DO
          END DO
        END IF
      END IF

      DO ls=ms, ls_trunc

        lsr=ls+1-ms

        DO l=1, n_profile
          surface_term(l, lsr)=0.0e+00_RealK
        END DO
        DO j=1, n_brdf_basis_fnc
          phi=0.0e+00_RealK
          DO ls_d=ms, ls_brdf_trunc-MOD(ms, 2), 2
            lsr_d=ls_d-ms+1
            phi_d=0.0e+00_RealK
            DO ls_dd=ms, ls_brdf_trunc-MOD(ms, 2), 2
              lsr_dd=ls_dd-ms+1
              phi_d=phi_d+cgk((lsr_dd+1)/2, lsr)                        &
                *f_brdf(j, ls_d/2, ls_dd/2, ms)
            END DO
            phi=phi+kappa(lsr_p/2, (lsr_d+1)/2)*phi_d
          END DO
          DO l=1, n_profile
            surface_term(l, lsr)=surface_term(l, lsr)                   &
              +rho_alb(l, j)*phi*REAL(1-2*MOD(ms, 2), RealK)
          END DO
        END DO

!       Add on the contribution to the RHS.
        DO l=1, n_profile
          b(l, ie)=b(l, ie)                                             &
            -source_bottom(l, lsr, i_above)*surface_term(l, lsr)
        END DO
!       Add in the contributions to the matrix on the LHS.
        DO k=1, n_red_eigensystem
!         Variable numbers:
          ivma=k
          ivpa=ivma+n_red_eigensystem
          DO l=1, n_profile
            a(l, ie, ivma)=a(l, ie, ivma)                               &
              +surface_term(l, lsr)*REAL(1-2*MOD(lsr-1, 2), RealK)      &
              *eig_vec(l, lsr, k, i_above)
            a(l, ie, ivpa)=a(l, ie, ivpa)                               &
              +surface_term(l, lsr)*eig_vec(l, lsr, k, i_above)
          END DO

        END DO

      END DO

    END IF

    DO k=1, n_red_eigensystem
      ivma=k
      DO l=1, n_profile
        a(l, ie, ivma)=a(l, ie, ivma)*theta(l, k, i_above)
      END DO
    END DO

  END DO


  IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
!   Calculate the contribution of radiation reflected from the
!   surface.

    CALL calc_surf_rad(n_profile, n_layer, tau                          &
      , ms, ls_trunc, euler_factor                                      &
      , isolir, i_direct(1, n_layer), mu_0, d_planck_flux_surface       &
      , n_brdf_basis_fnc, ls_brdf_trunc, f_brdf                         &
      , rho_alb, brdf_sol, brdf_hemi, cgk                               &
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
      , n_direction, mu_v, up_lm, azim_factor                           &
      , n_red_eigensystem, eig_vec(1, 1, 1, i_above)                    &
      , theta(1, 1, i_above), source_bottom(1, 1, i_above)              &
      , radiance                                                        &
      , weight_u(1, 1, 1, 1+2*n_red_eigensystem*(n_layer-1))            &
      , nd_profile, nd_layer, nd_direction, nd_viewing_level            &
      , nd_red_eigensystem, nd_max_order, nd_brdf_basis_fnc             &
      , nd_brdf_trunc                                                   &
      )

!   Isotropic incident fluxes are permitted (and required in the
!   differential formulation of the IR).
    IF (ms == 0) THEN
      CALL calc_top_rad(n_profile, tau                                  &
        , n_viewing_level, i_rad_layer, frac_rad_layer                  &
        , n_direction, mu_v                                             &
        , flux_down_inc                                                 &
        , radiance                                                      &
        , nd_profile, nd_layer, nd_direction, nd_viewing_level          &
        )
    END IF
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE build_sph_matrix
END MODULE build_sph_matrix_mod
