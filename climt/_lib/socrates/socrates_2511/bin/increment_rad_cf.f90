! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to increment radiances at a given azimuthal order.
!
! Method:
!   The weights of the terms in the complementary function
!   of the direct solution by spherical harmonics, u_{imk}^+-
!   are now available. For each viewing level and direction
!   we multiply by the precalculated coefficients and the
!   factor representing the azimuthal dependence to complete the
!   calculation of the radiance.
!
!- ---------------------------------------------------------------------
MODULE increment_rad_cf_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INCREMENT_RAD_CF_MOD'
CONTAINS
SUBROUTINE increment_rad_cf(n_profile                                   &
    , n_direction, azim_factor                                          &
    , n_viewing_level, i_rad_layer                                      &
    , i_sph_mode, i_sph_algorithm, ms, ls_trunc, euler_factor           &
    , isolir, mu_0, kappa, up_lm                                        &
    , n_red_eigensystem, n_equation, weight_u, upm                      &
    , i_direct, c_ylm, flux_direct, flux_total                          &
    , radiance, j_radiance                                              &
    , nd_profile, nd_flux_profile                                       &
    , nd_radiance_profile, nd_j_profile                                 &
    , nd_layer, nd_direction, nd_viewing_level                          &
    , nd_max_order, nd_sph_equation, nd_sph_cf_weight                   &
    , nd_sph_u_range                                                    &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_solar, ip_sph_direct, ip_sph_mode_flux,         &
                     ip_sph_mode_j, ip_sph_mode_rad, ip_sph_reduced_iter
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for points where radiances are calculated
    , nd_flux_profile                                                   &
!       Size allocated for profiles where fluxes are calculated
    , nd_radiance_profile                                               &
!       Size allocated for profiles where radiances are calculated
    , nd_j_profile                                                      &
!       Size allocated for profiles where mean radiances
!       are calculated
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_direction                                                      &
!       Size allocated for order of spherical calculation
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_max_order                                                      &
!       Size allocated for orders of direct calculation of
!       spherical harmonics
    , nd_sph_equation                                                   &
!       Size allocated for spherical equations
    , nd_sph_cf_weight                                                  &
!       Size allocated for entities to be weighted by the C. F.
    , nd_sph_u_range
!       Size allowed for range of values of u^+|- contributing
!       on any viewing level


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile
!       Number of profiles

! Spectral decomposition
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Spectral region

! Viewing geometry:
  INTEGER, INTENT(IN) ::                                                &
      n_direction                                                       &
!       Number of directions in which radiances are calculated
    , n_viewing_level                                                   &
!       Number of levels where the radiance is calculated
    , i_rad_layer(nd_viewing_level)
!       Layers of the atmosphere in which viewing levels fall
  REAL (RealK), INTENT(IN) ::                                           &
      azim_factor(nd_profile, nd_direction)                             &
!       Factors for representing the azimuthal dependence
    , mu_0(nd_profile)                                                  &
!       Cosines of the solar zenith angles
    , kappa(nd_max_order/2, nd_max_order/2)                             &
!       Integrals of Y_l^m*.Y_l^m over the downward hemisphere
    , up_lm(nd_profile, nd_max_order+1, nd_direction)
!       Polar parts of spherical harmonics

! Angular integration:
  INTEGER, INTENT(IN) ::                                                &
      i_sph_mode                                                        &
!       Mode in which the spherical harmonic code is called
    , i_sph_algorithm                                                   &
!       Algorithm used to solve spherical harmonic problem
    , ms                                                                &
!       Azimuthal order of truncation
    , ls_trunc
!       Polar order of truncation
  REAL (RealK), INTENT(IN) ::                                           &
      euler_factor
!       Factor weighting the last term of the series

! Components of the solution of the linear system
  INTEGER, INTENT(IN) ::                                                &
      n_red_eigensystem                                                 &
!       Size of the reduced eigensystem
    , n_equation
!       Number of spherical equations
  REAL (RealK), INTENT(IN) ::                                           &
      weight_u(nd_profile, nd_viewing_level                             &
        , nd_sph_cf_weight, nd_sph_u_range)                             &
!       Weights for coefficients in equations
    , upm(nd_profile, nd_sph_equation)
!       Variables u+|-

  REAL (RealK), INTENT(IN) ::                                           &
      i_direct(nd_profile, 0: nd_layer)
!       Direct radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      c_ylm(nd_profile, nd_viewing_level, ls_trunc+1-ms)
!       Spherical harmonic coefficients for radiances

! Calculated radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_radiance_profile, nd_viewing_level, nd_direction)
!       Radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct fluxes
    , flux_total(nd_flux_profile, 2*nd_layer+2)
!       Total fluxes

!       Mean radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      j_radiance(nd_j_profile, nd_viewing_level)
!       Mean radiances


! Local arguments.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , iv                                                                &
!       Loop variable (viewing level)
    , ie                                                                &
!       Loop variable
    , ls                                                                &
!       Polar order
    , lsr                                                               &
!       Reduced polar order
    , offset_u
!       Offset applied to the elements of u^+|- to move to
!       elements relevant to the current layer
  REAL (RealK) ::                                                       &
      contribution                                                      &
!       Contribution of the current order to the flux
    , cnst_ls
!       Constant term involving the polar order

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='INCREMENT_RAD_CF'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (i_sph_algorithm == ip_sph_direct) THEN

!   Radiances or fluxes are calculated directly from the
!   spherical harmonics.

!   Determine the coefficients of the spherical harmonics
!   from the solution of the eigenproblem.
    DO iv=1, n_viewing_level
      offset_u=2*n_red_eigensystem*(i_rad_layer(iv)-1)
      DO k=1, 2*n_red_eigensystem
        DO lsr=1, ls_trunc+1-ms
          DO l=1, n_profile
            c_ylm(l, iv, lsr)=c_ylm(l, iv, lsr)                         &
              +weight_u(l, iv, lsr, k)*upm(l, k+offset_u)
          END DO
        END DO
      END DO
    END DO

    IF (i_sph_mode == ip_sph_mode_flux) THEN

!     Although this routine is called to increment radiances over
!     angular orders, when run to calculate fluxes it should
!     only be called once during each monochromatic calculation.

      DO iv=1, n_viewing_level
        DO l=1, n_profile
          contribution=c_ylm(l, iv, 2)*SQRT(pi/3.0e+00_RealK)
!         Upward flux:
          flux_total(l, 2*iv-1)=contribution
!         Downward flux:
          flux_total(l, 2*iv)=-contribution
        END DO
      END DO
      DO ls=0, ls_trunc, 2
        cnst_ls=2.0e+00_RealK*kappa(1, (ls+2)/2)                        &
          *SQRT(pi/3.0e+00_RealK)
        DO iv=1, n_viewing_level
          DO l=1, n_profile
            contribution=cnst_ls*c_ylm(l, iv, ls+1)
            flux_total(l, 2*iv-1)                                       &
              =flux_total(l, 2*iv-1)-contribution
            flux_total(l, 2*iv)                                         &
              =flux_total(l, 2*iv)-contribution
          END DO
        END DO
      END DO

      IF (isolir == ip_solar) THEN
        DO iv=1, n_viewing_level
          DO l=1, n_profile
            flux_direct(l, iv-1)=i_direct(l, iv-1)*mu_0(l)
            flux_total(l, 2*iv)=flux_total(l, 2*iv)                     &
              +flux_direct(l, iv-1)
          END DO
        END DO
      END IF

    ELSE IF (i_sph_mode == ip_sph_mode_j) THEN

!     Although this routine is called to increment radiances over
!     angular orders, when run to calculate mean radiances it should
!     be called only once during each monochromatic calculation.

      DO iv=1, n_viewing_level
        DO l=1, n_profile
          j_radiance(l, iv)=c_ylm(l, iv, 2)*SQRT(4.0e+00_RealK*pi)
        END DO
      END DO

      IF (isolir == ip_solar) THEN
        DO iv=1, n_viewing_level
          DO l=1, n_profile
            j_radiance(l, iv)=j_radiance(l, iv)                         &
              +i_direct(l, iv)
          END DO
        END DO
      END IF

    ELSE IF (i_sph_mode == ip_sph_mode_rad) THEN

!     Determine the radiances directly from the amplitudes of
!     the harmonics.
      DO id=1, n_direction

!       Add in the contributions on each viewing level. To improve
!       convergence of the alternating series the contribution
!       from the last term may be reduced in size.
        DO iv=1, n_viewing_level
          DO lsr=1, ls_trunc-ms
            DO l=1, n_profile
              radiance(l, iv, id)=radiance(l, iv, id)                   &
                +azim_factor(l, id)*c_ylm(l, iv, lsr)                   &
                *up_lm(l, lsr, id)
            END DO
          END DO
          DO l=1, n_profile
            radiance(l, iv, id)=radiance(l, iv, id)+euler_factor        &
              *azim_factor(l, id)*c_ylm(l, iv, ls_trunc+1-ms)           &
              *up_lm(l, ls_trunc+1-ms, id)
          END DO
        END DO

      END DO
    END IF

  ELSE IF (i_sph_algorithm == ip_sph_reduced_iter) THEN

    DO id=1, n_direction
      DO iv=1, n_viewing_level
        DO ie=1, n_equation
          DO l=1, n_profile
            radiance(l, iv, id)=radiance(l, iv, id)                     &
              +azim_factor(l, id)                                       &
              *weight_u(l, iv, id, ie)*upm(l, ie)
          END DO
        END DO
      END DO
    END DO


  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE increment_rad_cf
END MODULE increment_rad_cf_mod
