! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate monochromatic radiances using IPA.
!
! Method:
!   In this subroutine a long vector for radiance calculations
!   is set up using the information on the types of cloud present.
!
!- ---------------------------------------------------------------------
MODULE calc_radiance_ipa_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_RADIANCE_IPA_MOD'
CONTAINS
SUBROUTINE calc_radiance_ipa(                                           &
!                 Atmospheric Properties
      n_profile, n_layer, n_cloud_top                                   &
!                 Angular Integration
    , n_order_phase, ms_min, ms_max, ls_local_trunc                     &
    , i_truncation, accuracy_adaptive, euler_factor                     &
    , i_sph_algorithm, i_sph_mode, l_rescale                            &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , diff_planck                                                       &
    , l_ir_source_quad, diff_planck_2                                   &
!                 Conditions at TOA
    , flux_inc_down, zen_0                                              &
!                 Conditions at Surface
    , d_planck_flux_surface                                             &
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
!                 Optical Properties
    , ss_prop                                                           &
!                 Cloud Geometry
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                 Viewing Geometry
    , n_direction, direction                                            &
!                 Calculated fluxes or radiances
    , flux_direct, flux_total, i_direct, radiance, j_radiance           &
!                 Dimensions of Arrays
    , nd_profile, nd_layer                                              &
    , nd_column                                                         &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc                                  &
    , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal              &
    , nd_sph_cf_weight, nd_sph_u_range                                  &
    , nd_viewing_level, nd_direction                                    &
    , nd_profile_column                                                 &
    )


  USE realtype_rd, ONLY: RealK
  USE def_ss_prop, ONLY: str_ss_prop
  USE rad_pcf, ONLY: ip_infra_red, ip_solar, ip_sph_mode_flux,          &
                     ip_sph_reduced_iter, ip_sph_mode_rad
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE sph_solver_mod, ONLY: sph_solver

  IMPLICIT NONE

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_flux_profile                                                   &
!       Size allocated for profiles of output fluxes
    , nd_radiance_profile                                               &
!       Size allocated for profiles of radiances
    , nd_j_profile                                                      &
!       Size allocated for profiles of photolysis rates
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_column                                                         &
!       Size allocated for columns at a grid-point
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
!       Size allowed for the spherical harmonic eigensystem
    , nd_sph_equation                                                   &
!       Size allowed for spherical harmonic equations
    , nd_sph_diagonal                                                   &
!       Size allowed for diagonals of the spherical harmonic
!       matrix
    , nd_sph_cf_weight                                                  &
!       Size allowed for application of weights of the C. F.
    , nd_sph_u_range                                                    &
!       Size allowed for range of values of u^+|- contributing
!       on any viewing level
    , nd_direction                                                      &
!       Size allocated for viewing dierctions
    , nd_profile_column
!       Number of profiles of subcolumns considered at once


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_order_phase
!       Number of orders retained in the phase function

!                 Spherical arrays
  INTEGER, INTENT(IN) ::                                                &
      ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , i_truncation                                                      &
!       Type of speherical truncation
    , i_sph_mode                                                        &
!       Mode in which the spherical harmonic solver is used
    , i_sph_algorithm                                                   &
!       Spherical harmonic algorithm
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)
!       Orders of truncation at each azimuthal order
  REAL (RealK) ::                                                       &
      cg_coeff(nd_sph_coeff)                                            &
!       Clebsch-Gordan coefficients
    , uplm_zero(nd_sph_coeff)                                           &
!       Values of spherical harmonics at polar angles pi/2
    , uplm_sol(nd_profile, nd_sph_coeff)
!       Values of spherical harmonics in the solar direction
  REAL (RealK), INTENT(IN) ::                                           &
      accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Spectral region
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar                                                     &
!       Scale solar beam
    , l_ir_source_quad                                                  &
!       Use a quadratic source term
    , l_rescale
!       Flag for rescaling of the optical properties

! Fields for equivalent extinction
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment of solar beam with equivalent extinction

! Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

! Planckian terms:
  REAL (RealK), INTENT(IN) ::                                           &
      diff_planck(nd_profile, nd_layer)                                 &
!       Change in Planckian function
    , diff_planck_2(nd_profile, nd_layer)                               &
!       Twice 2nd differences in Planckian
    , d_planck_flux_surface(nd_profile)
!       Differential Planckian flux from the surface

! Conditions at TOA
  REAL (RealK), INTENT(IN) ::                                           &
      zen_0(nd_profile)                                                 &
!       Secant of zenith angle
    , flux_inc_down(nd_profile)
!       Incident total flux

! Conditions at surface
  INTEGER, INTENT(IN) ::                                                &
      ls_brdf_trunc                                                     &
!       Order of trunation of BRDFs
    , n_brdf_basis_fnc
!       Number of BRDF basis functions
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb(nd_profile, nd_brdf_basis_fnc)                            &
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

! Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_column_slv(nd_profile)                                          &
!       Number of columns to be solved in each profile
    , list_column_slv(nd_profile, nd_column)                            &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(nd_profile, nd_column)                              &
!       Layer in the current column to change
    , i_clm_cld_typ(nd_profile, nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK), INTENT(IN) ::                                           &
      area_column(nd_profile, nd_column)
!       Area of each column

!                 Levels at which radiances will be calculated
  INTEGER, INTENT(IN) ::                                                &
      n_viewing_level                                                   &
!       Number of levels where radiances are calculated
    , i_rad_layer(nd_viewing_level)
!       Layers in which radiances are calculated
  REAL (RealK), INTENT(IN) ::                                           &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers

!                   Viewing Geometry
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)
!       Viewing directions

!                   Calculated Fluxes or Radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_profile, 0: nd_layer)
!       Direct radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct Flux
    , flux_total(nd_flux_profile, 2*nd_layer+2)                         &
!       Total Fluxes
    , radiance(nd_radiance_profile, nd_viewing_level, nd_direction)     &
!       Radiances
    , j_radiance(nd_j_profile, nd_viewing_level)
!       Photolysis rates



! Local variables.
  INTEGER                                                               &
      ls                                                                &
!       Polar order of harmonic
    , ms                                                                &
!       Azimuthal order of harmonic
    , i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , js                                                                &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , lp                                                                &
!       Index of current real grid-point during assignments
    , ll                                                                &
!       Index in the long array of columns to be taken in one go
    , ll_copy                                                           &
!       Index of column to be copied
    , icl                                                               &
!       Index of notional sub-column
    , ics                                                               &
!       Index of current sub-column where a solution is required
    , icc                                                               &
!       Temporary variable listing the layer in the current column
!       where a change is required
    , ict
!       Temporary variable listing the type of optical region moved
!       into be the current change
  INTEGER                                                               &
      n_long                                                            &
!       Length of long vector
    , target(nd_profile_column)
!       Actual target grid-point for point in the long array
  REAL (RealK) ::                                                       &
      weight_column(nd_profile_column)
!       Weight applied to each column in the sum
  LOGICAL                                                               &
      l_new
!       Flag to consider a new grid-point

! Properties of vectors of subcolumns
  REAL (RealK) ::                                                       &
      tau_long(nd_profile_column, nd_layer)                             &
!       Long vector of optical depth
    , omega_long(nd_profile_column, nd_layer)                           &
!       Long vector of albedo of single scattering
    , phase_fnc_long(nd_profile_column, nd_layer, nd_max_order)         &
!       Long vector of phase functions
    , phase_fnc_solar_long(nd_profile_column                            &
        , nd_layer, nd_direction)                                       &
!       Long vector of solar phase functions
    , forward_scatter_long(nd_profile_column, nd_layer)                 &
!       Long vector of forward scattering fractions
    , adjust_solar_ke_long(nd_profile_column, nd_layer)                 &
!       Long vector of solar scalings
    , zen_0_long(nd_profile_column)                                     &
!       Long vector of cosines of the solar zenith angle
    , uplm_sol_long(nd_profile_column, nd_sph_coeff)                    &
!       Long vector of spherical harmonics at the solar angle
    , diff_planck_long(nd_profile_column, nd_layer)                     &
!       Long vector of differences in the Planckian
    , diff_planck_2_long(nd_profile_column, nd_layer)                   &
!       Long vector of second differences in the Planckian
    , flux_inc_down_long(nd_profile_column)                             &
!       Long vector of incident downward fluxes
    , d_planck_flux_surface_long(nd_profile_column)                     &
!       Long vector of differential Planckian fluxes
!       at the surface
    , rho_alb_long(nd_profile_column, nd_brdf_basis_fnc)                &
!       Long vector of weightings of BRDF basis functions
    , brdf_sol_long(nd_profile_column, nd_brdf_basis_fnc                &
        , nd_direction)                                                 &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi_long(nd_profile_column, nd_brdf_basis_fnc               &
        , nd_direction)                                                 &
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction
    , direction_long(nd_profile_column, nd_direction, 2)
!       Viewing directions

! Calculated Fluxes or Radiances in subcolumns
  REAL (RealK) ::                                                       &
      flux_direct_column(nd_profile_column, 0: nd_layer)                &
!       Direct Flux
    , flux_total_column(nd_profile_column, 2*nd_layer+2)                &
!       Total Fluxes
    , i_direct_column(nd_profile_column, 0: nd_layer)                   &
!       Direct radiances
    , radiance_column(nd_profile_column, nd_viewing_level               &
        , nd_direction)                                                 &
!       Radiances
    , photolysis_column(nd_profile_column, nd_viewing_level)
!       Photolysis rates

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_RADIANCE_IPA'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Zero the output arrays ready for incrementing.
  IF (i_sph_mode == ip_sph_mode_flux) THEN

    DO i=1, 2*n_layer+2
      DO l=1, n_profile
        flux_total(l, i)=0.0e+00_RealK
      END DO
    END DO

    IF (isolir == ip_solar) THEN
      DO i=0, n_layer
        DO l=1, n_profile
          flux_direct(l, i)=0.0e+00_RealK
        END DO
      END DO
    END IF

  ELSE

    DO id=1, n_direction
      DO i=1, n_viewing_level
        DO l=1, n_profile
          radiance(l, i, id)=0.0e+00_RealK
        END DO
      END DO
    END DO

    DO i=1, n_viewing_level
      DO l=1, n_profile
        j_radiance(l, i)=0.0e+00_RealK
      END DO
    END DO

    IF (isolir == ip_solar) THEN
!     The top level contains the input: other values are zeroed
!     to allow incrementing.
      DO i=1, n_layer
        DO l=1, n_profile
          i_direct(l, i)=0.0e+00_RealK
        END DO
      END DO
    END IF

  END IF


! Start feeding points into the long array. This is
! not written to vectorize as that is quite complicated.

  lp=1
  l_new=.TRUE.

  DO WHILE (lp <= n_profile)

    ll=0

    DO WHILE ( (ll <  nd_profile_column).AND.(lp <= n_profile) )

      ll=ll+1
      target(ll)=lp

      IF (l_new) THEN

!       We consider a new grid-point and so must set the first
!       notional column which is contains no cloud.
        icl=1
        ics=1
        DO i=1, n_cloud_top-1
          tau_long(ll, i)=ss_prop%tau_clr(lp, i)
          omega_long(ll, i)=ss_prop%omega_clr(lp, i)
          DO ls=1, n_order_phase
            phase_fnc_long(ll, i, ls)=ss_prop%phase_fnc_clr(lp, i, ls)
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          tau_long(ll, i)=ss_prop%tau(lp, i, 0)
          omega_long(ll, i)=ss_prop%omega(lp, i, 0)
          DO ls=1, n_order_phase
            phase_fnc_long(ll, i, ls)=ss_prop%phase_fnc(lp, i, ls, 0)
          END DO
        END DO
        IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
          DO id=1, n_direction
            DO i=1, n_cloud_top-1
              phase_fnc_solar_long(ll, i, id)                           &
                =ss_prop%phase_fnc_solar_clr(lp, i, id)
            END DO
            DO i=n_cloud_top, n_layer
              phase_fnc_solar_long(ll, i, id)                           &
                =ss_prop%phase_fnc_solar(lp, i, id, 0)
            END DO
          END DO
          IF (l_rescale) THEN
            DO i=1, n_cloud_top-1
              forward_scatter_long(ll, i)                               &
                =ss_prop%forward_scatter_clr(lp, i)
            END DO
            DO i=n_cloud_top, n_layer
              forward_scatter_long(ll, i)                               &
                =ss_prop%forward_scatter(lp, i, 0)
            END DO
          END IF
        END IF


        l_new=.FALSE.


      ELSE

!       Copy the previous column over. Normally this will be the
!       previous one, but if we are starting a new batch it will
!       be the one at the end of the previous batch.
        IF (ll >  1) THEN
          ll_copy=ll-1
        ELSE
          ll_copy=n_long
        END IF

        DO i=1, n_layer
          tau_long(ll, i)=tau_long(ll_copy, i)
          omega_long(ll, i)=omega_long(ll_copy, i)
          DO ls=1, n_order_phase
            phase_fnc_long(ll, i, ls)                                   &
              =phase_fnc_long(ll_copy, i, ls)
          END DO
        END DO
        IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
          DO id=1, n_direction
            DO i=1, n_layer
              phase_fnc_solar_long(ll, i, id)                           &
                =phase_fnc_solar_long(ll_copy, i, id)
            END DO
          END DO
          IF (l_rescale) THEN
            DO i=1, n_layer
              forward_scatter_long(ll, i)                               &
                =forward_scatter_long(ll_copy, i)
            END DO
          END IF
        END IF

      END IF

!     Move through the notional columns at this grid-point
!     adjusting individiual layers until we find one where the
!     equations are to be solved.
      DO WHILE (icl <  list_column_slv(lp, ics))
        icc=i_clm_lyr_chn(lp, icl)
        ict=i_clm_cld_typ(lp, icl)

        tau_long(ll, icc)=ss_prop%tau(lp, icc, ict)
        omega_long(ll, icc)=ss_prop%omega(lp, icc, ict)
        DO ls=1, n_order_phase
          phase_fnc_long(ll, icc, ls)                                   &
            =ss_prop%phase_fnc(lp, icc, ls, ict)
        END DO
        IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
          DO id=1, n_direction
            phase_fnc_solar_long(ll, icc, id)                           &
              =ss_prop%phase_fnc_solar(lp, icc, id, ict)
          END DO
          IF (l_rescale) THEN
            forward_scatter_long(ll, icc)                               &
              =ss_prop%forward_scatter(lp, icc, ict)
          END IF
        END IF

        icl=icl+1
      END DO


!     Set arrays which are independent of cloud changes.
      IF (isolir == ip_solar) THEN

        IF (l_scale_solar) THEN
          DO i=1, n_layer
            adjust_solar_ke_long(ll, i)=adjust_solar_ke(lp, i)
          END DO
        END IF

        zen_0_long(ll)=zen_0(lp)
        i_direct_column(ll, 0)=i_direct(lp, 0)
        DO ms=ms_min, ms_max
          DO ls=ms, ls_local_trunc(ms)+1
            js=ia_sph_mm(ms)+ls-ms
            uplm_sol_long(ll, js)=uplm_sol(lp, js)
          END DO
        END DO

        IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
          DO id=1, n_direction
            DO j=1, n_brdf_basis_fnc
              brdf_sol_long(ll, j, id)=brdf_sol(lp, j, id)
            END DO
          END DO
        END IF

      ELSE IF (isolir == ip_infra_red) THEN

        d_planck_flux_surface_long(ll)                                  &
          =d_planck_flux_surface(lp)
        DO i=1, n_layer
          diff_planck_long(ll, i)=diff_planck(lp, i)
        END DO
        IF (l_ir_source_quad) THEN
          DO i=1, n_layer
            diff_planck_2_long(ll, i)=diff_planck_2(lp, i)
          END DO
        END IF
        IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
          DO id=1, n_direction
            DO j=1, n_brdf_basis_fnc
              brdf_hemi_long(ll, j, id)=brdf_hemi(lp, j, id)
            END DO
          END DO
        END IF

      END IF

!     Set the viewing directions.
      IF (i_sph_mode == ip_sph_mode_rad) THEN
        DO id=1, n_direction
          direction_long(ll, id, 1)=direction(lp, id, 1)
          direction_long(ll, id, 2)=direction(lp, id, 2)
        END DO
      END IF

      flux_inc_down_long(ll)=flux_inc_down(lp)
      DO js=1, n_brdf_basis_fnc
        rho_alb_long(ll, js)=rho_alb(lp, js)
      END DO

!     The curent notional column will contain the fraction of
!     the grid-box required for incrementing.
      weight_column(ll)=area_column(lp, icl)


!     Prepare for the next column, moving on to the next grid-point
!     as required.
      ics=ics+1
      IF (ics >  n_column_slv(lp)) THEN
        lp=lp+1
        l_new=.TRUE.
      END IF

    END DO

!   Set N_LONG which will be required for the next batch after LL
!   has been reset.
    n_long=ll


    CALL sph_solver(                                                    &
!                   Atmospheric sizes
        n_long, n_layer                                                 &
!                   Angular integration
      , ms_min, ms_max, i_truncation, ls_local_trunc                    &
      , cg_coeff, uplm_zero, ia_sph_mm                                  &
      , accuracy_adaptive, euler_factor                                 &
      , i_sph_algorithm, i_sph_mode                                     &
!                   Spectral Region
      , isolir                                                          &
!                   Options for Equivalent Extinction
      , l_scale_solar, adjust_solar_ke_long                             &
!                   Solar Fields
      , i_direct_column, zen_0_long, uplm_sol_long                      &
!                   Infra-red Properties
      , diff_planck_long, flux_inc_down_long                            &
      , l_ir_source_quad, diff_planck_2_long                            &
!                   Optical properies
      , tau_long, omega_long, phase_fnc_long                            &
      , phase_fnc_solar_long                                            &
!                   Surface Conditions
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb_long                   &
      , f_brdf, brdf_sol_long, brdf_hemi_long                           &
      , d_planck_flux_surface_long                                      &
!                   Levels for calculating radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!                   Viewing Geometry
      , n_direction, direction_long                                     &
!                   Calculated Radiances or Fluxes
      , flux_direct_column, flux_total_column, radiance_column          &
      , photolysis_column                                               &
!                   Dimensions of arrays
      , nd_profile_column, nd_layer                                     &
      , nd_profile_column, nd_profile_column, nd_profile_column         &
      , nd_max_order, nd_sph_coeff                                      &
      , nd_brdf_basis_fnc, nd_brdf_trunc                                &
      , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal            &
      , nd_sph_cf_weight, nd_sph_u_range                                &
      , nd_viewing_level, nd_direction                                  &
      )



!   Scatter the calculated fluxes or radiances back to their
!   appropriate grid-points.
    IF (i_sph_mode == ip_sph_mode_flux) THEN

      DO i=1, 2*n_layer+2
        DO ll=1, n_long
          l=target(ll)
          flux_total(l, i)=flux_total(l, i)                             &
            +weight_column(ll)*flux_total_column(ll, i)
        END DO
      END DO

      IF (isolir == ip_solar) THEN
        DO i=0, n_layer
          DO ll=1, n_long
            l=target(ll)
            flux_direct(l, i)=flux_direct(l, i)                         &
              +weight_column(ll)*flux_direct_column(ll, i)
          END DO
        END DO
      END IF

    ELSE


      DO id=1, n_direction
        DO i=1, n_viewing_level
          DO ll=1, n_long
            l=target(ll)
            radiance(l, i, id)=radiance(l, i, id)                       &
              +weight_column(ll)*radiance_column(ll, i, id)
          END DO
        END DO
      END DO

      DO i=1, n_viewing_level
        DO ll=1, n_long
          l=target(ll)
          j_radiance(l, i)=j_radiance(l, i)                             &
            +weight_column(ll)*photolysis_column(ll, i)
        END DO
      END DO

      IF (isolir == ip_solar) THEN
        DO i=1, n_layer
          DO ll=1, n_long
            l=target(ll)
            i_direct(l, i)=i_direct(l, i)                               &
              +weight_column(ll)*i_direct_column(ll, i)
          END DO
        END DO
      END IF

    END IF


  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_radiance_ipa
END MODULE calc_radiance_ipa_mod
