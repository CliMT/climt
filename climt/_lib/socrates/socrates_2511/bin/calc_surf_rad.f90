! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set weights for the C.F. at the surface
!
! Purpose:
!   The contribution to the radiance of radiation reflected from the
!   surface is evaluated.
!
! Method:
!   The iterated expression for the radiance involves a contribution
!   from the radiance reflected from the surface. In principle, this
!   could be provided by the upward radiance, but in practice this
!   would be of low accuracy since the spherical harmonic series for
!   the radiance will be noisy at low orders of truncation. It is
!   better to evaluate the reflected radiance using the BRDFs, even
!   though it is more expensive to do so; this ensures that no
!   radiation will appear to be reflected from a non-reflecting
!   surface. Given these constraints the algorithm is essentially
!   straightforward.
!
!- ---------------------------------------------------------------------
MODULE calc_surf_rad_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_SURF_RAD_MOD'
CONTAINS
SUBROUTINE calc_surf_rad(n_profile, n_layer, tau                        &
    , ms, ls_trunc, euler_factor                                        &
    , isolir, i_direct_surf, mu_0, d_planck_flux_surface                &
    , n_brdf_basis_fnc, ls_brdf_trunc, f_brdf                           &
    , rho_alb, brdf_sol, brdf_hemi, cgk                                 &
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
    , n_direction, mu_v, up_lm, azim_factor                             &
    , n_red_eigensystem, eig_vec, theta, source_base                    &
    , radiance, weight_u                                                &
    , nd_profile, nd_layer, nd_direction, nd_viewing_level              &
    , nd_red_eigensystem, nd_max_order, nd_brdf_basis_fnc               &
    , nd_brdf_trunc                                                     &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_solar, ip_infra_red
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_red_eigensystem                                                &
!       Size allocated for the reduced eigensystem
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_brdf_basis_fnc                                                 &
!       Size allocated for basis functions of BRDFs
    , nd_brdf_trunc
!       Size allocated for orders in basis functions of BRDFs


! The atmosphere:
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric profiles
    , n_layer
!       Number of atmospheric layers
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, nd_layer)
!       Optical depths

! Controlling spherical orders:
  INTEGER, INTENT(IN) ::                                                &
      ms                                                                &
!       Current azimuthal order
    , ls_trunc
!       Order of polar truncation
  REAL (RealK), INTENT(IN) ::                                           &
      euler_factor
!       Factor applied to the last term of the series

! Variables for solar or thermal sources
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Spectral region
  REAL (RealK), INTENT(IN) ::                                           &
      i_direct_surf(nd_profile)                                         &
!       The direct solar radiance at the surface
    , mu_0(nd_profile)
!       Cosines of the solar zenith angle
  REAL (RealK), INTENT(IN) ::                                           &
      d_planck_flux_surface(nd_profile)
!       Differential Planckian flux at the surface

! Variables related to the BRDFs
  INTEGER, INTENT(IN) ::                                                &
      n_brdf_basis_fnc                                                  &
!       Number of basis functions used in BRDFs
    , ls_brdf_trunc
!       Order of polar truncation applied to BRDFs
  REAL (RealK), INTENT(IN) ::                                           &
      f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                         &
!       BRDF basis functions
    , rho_alb(nd_profile, nd_brdf_basis_fnc)                            &
!       Weights of applied to the basis functions of the BRDF
    , brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)             &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction
  REAL (RealK), INTENT(IN) ::                                           &
      cgk(nd_brdf_trunc/2+1, nd_max_order)
!       Integrals of pairs of spherical harmonics over the downward
!       hemisphere


! Viewing geometry:
  INTEGER, INTENT(IN) ::                                                &
      n_direction                                                       &
!       Number of directions
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , i_rad_layer(nd_viewing_level)
!       Indices of layers containing viewing levels
  REAL (RealK), INTENT(IN) ::                                           &
      frac_rad_layer(nd_viewing_level)                                  &
!       Fraction optical depth into its layer of the
!       radiance level
    , mu_v(nd_profile, nd_direction)                                    &
!       Cosines of polar viewing angles
    , azim_factor(nd_profile, nd_direction)                             &
!       Azimuthal factors
    , up_lm(nd_profile, nd_max_order+1, nd_direction)
!       Spherical harmonics at a fixed azimuthal order

  INTEGER, INTENT(IN) ::                                                &
      n_red_eigensystem
!       Size of the reduced eigensystem
  REAL (RealK), INTENT(IN) ::                                           &
      eig_vec(nd_profile, 2*nd_red_eigensystem                          &
        , nd_red_eigensystem)                                           &
!       Eigenvalues of the full eigensystem scaled by
!       the s-parameters
    , theta(nd_profile, nd_red_eigensystem)                             &
!       Array of exponentials of optical depths along slant paths
    , source_base(nd_profile, ls_trunc+1-ms)
!       Source function at the bottom of the layer


  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_profile, nd_viewing_level, nd_direction)
!       Radiances (to be incremented by the contribution of
!       the particular integral)
  REAL (RealK), INTENT(INOUT) ::                                        &
      weight_u(nd_profile, nd_viewing_level                             &
        , nd_direction, 2*nd_red_eigensystem)
!       Weights for the coefficients in the complementary
!       function


! Local variables
  INTEGER                                                               &
      l                                                                 &
!       Loop variable (points)
    , ir                                                                &
!       Loop variable (viewing levels)
    , id                                                                &
!       Loop variable (directions)
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , ll
!       Loop variable
  INTEGER                                                               &
      ls                                                                &
!       Loop variable (polar orders)
    , lsr                                                               &
!       Loop variable (reduced polar orders)
    , ls_d                                                              &
!       Loop variable (polar orders)
    , lsr_d                                                             &
!       Loop variable (reduced polar orders)
    , ls_dd                                                             &
!       Loop variable (polar orders)
    , lsr_dd
!       Loop variable (reduced polar orders)
  INTEGER                                                               &
      n_list_up                                                         &
!       Number of points where the viewing direction is upward
    , list_up(nd_profile)
!       List of points where the viewing direction is upward
  REAL (RealK) ::                                                       &
      trans(nd_profile)                                                 &
!       Tranmission along the line of sight from the surface
!       to the viewing level
    , x
!       Temporary variable
! Working arrays realated to the BRDF:
  REAL (RealK) ::                                                       &
      xy(nd_profile)                                                    &
!       Product of (Clebsch-Gordan coefficient * kappa) and
!       the spherical harmonic at the current order in the
!       viewing direction
    , ryx(nd_profile, ls_trunc-ms+1)                                    &
!       Sum over basis functions of products of the above
!       and albedo weights
    , rvyx_m(nd_profile, nd_red_eigensystem)                            &
!       Sum over polar orders of product of RYX and elements
!       of the eigenvalue for each eigenvalue for application
!       to terms in negative exponentials
    , rvyx_p(nd_profile, nd_red_eigensystem)                            &
!       Sum over polar orders of product of RYX and elements
!       of the eigenvalue for each eigenvalue for application
!       to terms in positive exponentials
    , rsyx(nd_profile)                                                  &
!       Sum over polar orders of product of RYX and elements
!       of the source function
    , brdf_full(nd_profile)
!       Full BRDF weighted and summed over all basis functions

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_SURF_RAD'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! For each direction and observing level we calculate the
! contribution of the particular integral to the radiance
! from the surface and appropriate weightings for the
! complementary function.

  DO id=1, n_direction

!   Collect upward directions.
    n_list_up=0
    DO l=1, n_profile
      IF (mu_v(l, id) >  0.0e+00_RealK) THEN
        n_list_up=n_list_up+1
        list_up(n_list_up)=l
      END IF
    END DO

!   Calculate the angular arrays related to the BRDF. At higher
!   azimuthal orders there will be contributions because all
!   terms of the BRDF that would contribute are beyond the
!   level of truncation and so are zero.

    IF (ms <= ls_brdf_trunc-MOD(ms, 2)) THEN

      DO j=1, n_brdf_basis_fnc
        DO ls=ms, ls_trunc
          lsr=ls-ms+1
          DO ls_d=ms, ls_brdf_trunc-MOD(ms, 2), 2
            lsr_d=ls_d-ms+1
            x=0.0e+00_RealK
            DO ls_dd=ms, ls_brdf_trunc-MOD(ms, 2), 2
              lsr_dd=ls_dd-ms+1
              x=x-cgk((lsr_dd+1)/2, lsr)*f_brdf(j, ls_d, ls_dd, ms)
            END DO
            IF (ls_d == ms) THEN
!             Initialize this time.
              DO l=1, n_profile
                xy(l)=x*up_lm(l, lsr_d, id)
              END DO
            ELSE
!             Now add the increments.
              DO l=1, n_profile
                xy(l)=xy(l)+x*up_lm(l, lsr_d, id)
              END DO
            END IF
          END DO
          IF (j  == 1) THEN
!           Initialize this time.
            DO l=1, n_profile
              ryx(l, lsr)=rho_alb(l, 1)*xy(l)
            END DO
          ELSE
!           Increment for subsequent basis functions.
            DO l=1, n_profile
              ryx(l, lsr)=ryx(l, lsr)+rho_alb(l, j)*xy(l)
            END DO
          END IF
        END DO
      END DO

      DO k=1, n_red_eigensystem
        DO l=1, n_profile
          x=euler_factor*ryx(l, ls_trunc-ms+1)                          &
            *eig_vec(l, ls_trunc-ms+1, k)
          rvyx_m(l, k)=x
          rvyx_p(l, k)=x
        END DO
        DO lsr= ls_trunc-ms, 1, -1
          DO l=1, n_profile
            x=ryx(l, lsr)*eig_vec(l, lsr, k)
            rvyx_m(l, k)=rvyx_m(l, k)                                   &
              +x*REAL(1-2*MOD(lsr-1, 2), RealK)
            rvyx_p(l, k)=rvyx_p(l, k)+x
          END DO
        END DO
        DO l=1, n_profile
          rvyx_m(l, k)=rvyx_m(l, k)*theta(l, k)
        END DO
      END DO

      DO l=1, n_profile
        rsyx(l)=euler_factor*ryx(l, ls_trunc-ms+1)                      &
          *source_base(l, ls_trunc-ms+1)
      END DO
      DO lsr= ls_trunc-ms, 1, -1
        DO l=1, n_profile
          rsyx(l)=rsyx(l)+ryx(l, lsr)*source_base(l, lsr)
        END DO
      END DO

    END IF


    DO ir=1, n_viewing_level

!     Calculate minus the slantwise transmission from the
!     surface to the level in question. TRANS is used is
!     hold intermediate results.
      DO ll=1, n_list_up
        l=list_up(ll)
        trans(l)                                                        &
          =(1.0e+00_RealK                                               &
          -frac_rad_layer(ir))*tau(l, i_rad_layer(ir))
      END DO
      DO i=i_rad_layer(ir)+1, n_layer
        DO ll=1, n_list_up
          l=list_up(ll)
          trans(l)=trans(l)+tau(l, i)
        END DO
      END DO
      DO ll=1, n_list_up
        l=list_up(ll)
        trans(l)                                                        &
          =EXP(-trans(l)/mu_v(l, id))
      END DO

!     Add in the terms from the BRDF if in range.
      IF (ms <= ls_brdf_trunc-MOD(ms,2)) THEN

        DO ll=1, n_list_up
          l=list_up(ll)
!         Add the contribution from the source function at the
!         base of the layer.
          radiance(l, ir, id)=radiance(l, ir, id)                       &
            +trans(l)*rsyx(l)*azim_factor(l, id)
        END DO

!       Increment the weights applied to the complementary function.
        DO k=1, n_red_eigensystem
          DO ll=1, n_list_up
            l=list_up(ll)
            weight_u(l, ir, id, k)=weight_u(l, ir, id, k)               &
              +trans(l)*rvyx_m(l, k)
            weight_u(l, ir, id, k+n_red_eigensystem)                    &
              =weight_u(l, ir, id, k+n_red_eigensystem)                 &
              +trans(l)*rvyx_p(l, k)
          END DO
        END DO

      END IF


!     Add the direct solar or thermal contributions to the radiance.
!     The azimuthal dependencies are included in the solar and
!     hemispheric parts of the BRDF, so the should be added in just
!     once, most naturally at the zeroth order.

      IF (ms == 0) THEN
        IF (isolir == ip_solar) THEN
          DO ll=1, n_list_up
            l=list_up(ll)
            brdf_full(l)=rho_alb(l, 1)*brdf_sol(l, 1, id)
          END DO

          DO j=2, n_brdf_basis_fnc
            DO ll=1, n_list_up
              l=list_up(ll)
              brdf_full(l)=brdf_full(l)                                 &
                +rho_alb(l, j)*brdf_sol(l, j, id)
            END DO


          END DO
          DO ll=1, n_list_up
            l=list_up(ll)
            radiance(l, ir, id)=radiance(l, ir, id)                     &
              +trans(l)*i_direct_surf(l)*mu_0(l)                        &
              *brdf_full(l)
          END DO

        ELSE IF (isolir == ip_infra_red) THEN
          DO ll=1, n_list_up
            l=list_up(ll)
            brdf_full(l)=rho_alb(l, 1)*brdf_hemi(l, 1, id)
          END DO
          DO j=2, n_brdf_basis_fnc
            DO ll=1, n_list_up
              l=list_up(ll)
              brdf_full(l)=brdf_full(l)                                 &
                +rho_alb(l, j)*brdf_hemi(l, j, id)
            END DO
          END DO
          DO ll=1, n_list_up
            l=list_up(ll)
            radiance(l, ir, id)=radiance(l, ir, id)                     &
              +trans(l)                                                 &
              *(1.0e+00_RealK-brdf_full(l))                             &
              *d_planck_flux_surface(l)/pi
          END DO
        END IF
      END IF

    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_surf_rad
END MODULE calc_surf_rad_mod
