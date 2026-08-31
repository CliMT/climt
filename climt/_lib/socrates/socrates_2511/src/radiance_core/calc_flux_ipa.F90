! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate monochromatic fluxes using IPA.
!
! Method:
!   In this subroutine a long vector for two-stream flux calculations
!   is set up using the information on the types of cloud present.
!
!- ---------------------------------------------------------------------
MODULE calc_flux_ipa_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_FLUX_IPA_MOD'
CONTAINS
SUBROUTINE calc_flux_ipa(ierr                                           &
    , control, bound                                                    &
!                 Atmospheric Properties
    , n_profile, n_layer, n_cloud_top                                   &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Algorithmic options
    , i_2stream, i_solver, i_scatter_method                             &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , diff_planck                                                       &
    , l_ir_source_quad, diff_planck_2                                   &
!                 Conditions at TOA
    , flux_inc_down, flux_inc_direct, sec_0                             &
!                 Conditions at Surface
    , d_planck_flux_surface, rho_alb                                    &
!                 Spherical geometry
    , sph                                                               &
!                 Optical Properties
    , ss_prop                                                           &
!                 Cloud Geometry
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                   Calculated fluxes
    , flux_direct, flux_total                                           &
!                 Options for clear-sky fluxes
    , l_clear, i_solver_clear                                           &
!                   Calculated fluxes
    , flux_direct_clear, flux_total_clear                               &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_profile_column, nd_source_coeff                                &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_direct_csr_scaling, ip_direct_noscaling,        &
                     ip_infra_red, ip_solar, ip_surf_alb_diff,          &
                     ip_surf_alb_dir
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE copy_clr_full_mod, ONLY: copy_clr_full
  USE two_stream_mod, ONLY: two_stream

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_layer_clr                                                      &
!       Size allocated for totally clear atmospheric layers
    , nd_column                                                         &
!       Size allocated for columns at a grid-point
    , nd_profile_column                                                 &
!       Number of profiles of subcolumns considered at once
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_source_coeff
!       Number of coefficients in the source function


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top
!       Topmost cloudy layer

  INTEGER, INTENT(IN) ::                                                &
      isolir                                                            &
!       Spectral region
    , i_2stream                                                         &
!       Two-stream scheme selected
    , i_solver                                                          &
!       Solver selected
    , i_scatter_method
!       Method of treating scattering
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar                                                     &
!       Scale solar beam
    , l_ir_source_quad
!       Use a quadratic source term
!       the singly scattered solar beam

! Fields for equivalent extinction
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment of solar beam with equivalent extinction

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!       Spherical geometry fields

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
      sec_0(nd_profile)                                                 &
!       Secant of zenith angle
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)
!       Incident total flux

! Conditions at surface
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb(nd_profile, 2)
!       Weights of the basis functions

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

!                   Calculated Fluxes
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct Flux
    , flux_total(nd_profile, 2*nd_layer+2)
!       Total Fluxes

!                 Options for clear-sky fluxes
  LOGICAL                                                               &
      l_clear
!       Flag for clear-sky fluxes
  INTEGER                                                               &
      i_solver_clear
!       Solver selected for clear-sky fluxes
!                   Calculated clear-sky fluxes
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct_clear(nd_profile, 0: nd_layer)                        &
!       Direct Clear-sky Flux
    , flux_total_clear(nd_profile, 2*nd_layer+2)
!       Total Clear-sky Flux



! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l                                                                 &
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
      weight_long(nd_profile_column)
!       Weight applied to each column in the sum
  LOGICAL                                                               &
      l_new
!       Flag to consider a new grid-point

! Properties of vectors of subcolumns
  REAL (RealK) ::                                                       &
      tau_long(nd_profile_column, nd_layer)                             &
!       Long vector of optical depth
    , tau_long_dir(nd_profile_column, nd_layer)                         &
!       Unscaled Long vector of optical depth
    , omega_long(nd_profile_column, nd_layer)                           &
!       Long vector of albedo of single scattering
    , asymmetry_long(nd_profile_column, nd_layer)                       &
!       Long vector of asymmetries
    , adjust_solar_ke_long(nd_profile_column, nd_layer)                 &
!       Long vector of solar scalings
    , sec_0_long(nd_profile_column)                                     &
!       Long vector of cosines of the solar zenith angle
    , diff_planck_long(nd_profile_column, nd_layer)                     &
!       Long vector of differences in the Planckian
    , diff_planck_2_long(nd_profile_column, nd_layer)                   &
!       Long vector of second differences in the Planckian
    , flux_inc_direct_long(nd_profile_column)                           &
!       Long vector of incident direct downward fluxes
    , flux_inc_down_long(nd_profile_column)                             &
!       Long vector of incident downward fluxes
    , d_planck_flux_surface_long(nd_profile_column)                     &
!       Long vector of differential Planckian fluxes
!       at the surface
    , rho_alb_long(nd_profile_column, 2)
!       Long vector of weightings of BRDF basis functions

! Calculated Fluxes in subcolumns
  REAL (RealK) ::                                                       &
      flux_direct_long(nd_profile_column, 0: nd_layer)                  &
!       Direct Flux
    , flux_total_long(nd_profile_column, 2*nd_layer+2)
!       Total Fluxes

! Clear-sky optical properties of the whole column
  REAL  (RealK), ALLOCATABLE ::                                         &
      tau_clr_f(:, :)                                                   &
!       Clear-sky optical depth for the whole column
    , tau_clr_dir_f(:, :)                                               &
!       Unscaled clear-sky optical depth
    , omega_clr_f(:, :)                                                 &
!       Clear-sky albedos of single scattering for the whole column
    , phase_fnc_clr_f(:, :, :)
!       Moments of the clear-sky phase function for the whole column

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_FLUX_IPA'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Zero the output arrays ready for incrementing.
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
          asymmetry_long(ll, i)=ss_prop%phase_fnc_clr(lp, i, 1)
        END DO
        DO i=n_cloud_top, n_layer
          tau_long(ll, i)=ss_prop%tau(lp, i, 0)
          omega_long(ll, i)=ss_prop%omega(lp, i, 0)
          asymmetry_long(ll, i)=ss_prop%phase_fnc(lp, i, 1, 0)
        END DO

        IF (control%i_direct_tau == ip_direct_noscaling .OR.            &
            control%i_direct_tau == ip_direct_csr_scaling) THEN
          DO i=1, n_cloud_top-1
            tau_long_dir(ll, i)=ss_prop%tau_clr_dir(lp, i)
          END DO
          DO i=n_cloud_top, n_layer
            tau_long_dir(ll, i)=ss_prop%tau_dir(lp, i, 0)
          END DO
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
          asymmetry_long(ll, i)                                         &
            =asymmetry_long(ll_copy, i)
        END DO
        IF (control%i_direct_tau == ip_direct_noscaling .OR.            &
            control%i_direct_tau == ip_direct_csr_scaling) THEN
          DO i=1, n_layer
            tau_long_dir(ll, i)=tau_long_dir(ll_copy, i)
          END DO
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
        asymmetry_long(ll, icc)                                         &
          =ss_prop%phase_fnc(lp, icc, 1, ict)

        icl=icl+1
      END DO
      IF (control%i_direct_tau == ip_direct_noscaling .OR.              &
          control%i_direct_tau == ip_direct_csr_scaling) THEN
        DO WHILE (icl <  list_column_slv(lp, ics))
          icc=i_clm_lyr_chn(lp, icl)
          ict=i_clm_cld_typ(lp, icl)
          tau_long_dir(ll, icc)=ss_prop%tau_dir(lp, icc, ict)
          icl=icl+1
        END DO
      END IF


!     Set arrays which are independent of cloud changes.
      IF (isolir == ip_solar) THEN

        IF (l_scale_solar) THEN
          DO i=1, n_layer
            adjust_solar_ke_long(ll, i)=adjust_solar_ke(lp, i)
          END DO
        END IF

        sec_0_long(ll)=sec_0(lp)
        flux_inc_direct_long(ll)=flux_inc_direct(lp)
        d_planck_flux_surface_long(ll)=0.0e+00_RealK

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

      END IF

      flux_inc_down_long(ll)=flux_inc_down(lp)
      rho_alb_long(ll, ip_surf_alb_dir)                                 &
        =rho_alb(lp, ip_surf_alb_dir)
      rho_alb_long(ll, ip_surf_alb_diff)                                &
        =rho_alb(lp, ip_surf_alb_diff)

!     The curent notional column will contain the fraction of
!     the grid-box required for incrementing.
      weight_long(ll)=area_column(lp, icl)


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


!   N.B. The clear-sky option cannot be used here.
    CALL two_stream(ierr                                                &
      , control, bound                                                  &
!                   Atmospheric properties
      , n_long, n_layer                                                 &
!                   Two-stream scheme
      , i_2stream                                                       &
!                   Options for solver
      , i_solver, i_scatter_method                                      &
!                   Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke_long                             &
!                   Spectral region
      , isolir                                                          &
!                   Infra-red properties
      , diff_planck_long                                                &
      , l_ir_source_quad, diff_planck_2_long                            &
!                   Conditions at TOA
      , flux_inc_down_long, flux_inc_direct_long, sec_0_long            &
!                   Surface conditions
      , rho_alb_long(1, ip_surf_alb_diff)                               &
      , rho_alb_long(1, ip_surf_alb_dir), d_planck_flux_surface_long    &
!                   Spherical geometry
      , sph                                                             &
!                   Single scattering properties
      , tau_long_dir, tau_long, omega_long, asymmetry_long(1, 1)     &
!                   Fluxes calculated
      , flux_direct_long, flux_total_long                               &
!                   Sizes of arrays
      , nd_profile_column, nd_layer, nd_source_coeff                    &
      )



!   Scatter the calculated fluxes back to their
!   appropriate grid-points.

    DO i=1, 2*n_layer+2
      DO ll=1, n_long
        l=target(ll)
        flux_total(l, i)=flux_total(l, i)                               &
          +weight_long(ll)*flux_total_long(ll, i)
      END DO
    END DO

    IF (isolir == ip_solar) THEN
      DO i=0, n_layer
        DO ll=1, n_long
          l=target(ll)
          flux_direct(l, i)=flux_direct(l, i)                           &
            +weight_long(ll)*flux_direct_long(ll, i)
        END DO
      END DO
    END IF


  END DO

! Calculate the clear-sky fluxes if required.
  IF (l_clear) THEN

!   Set aside space for the clear optical properties and copy
!   them across.
    ALLOCATE(tau_clr_f(nd_profile, nd_layer))
    ALLOCATE(tau_clr_dir_f(nd_profile, nd_layer))
    ALLOCATE(omega_clr_f(nd_profile, nd_layer))
    ALLOCATE(phase_fnc_clr_f(nd_profile, nd_layer, 1))

    CALL copy_clr_full(n_profile, n_layer, n_cloud_top, control, 1      &
      , ss_prop%tau_clr, ss_prop%tau_clr_dir                            &
      , ss_prop%omega_clr, ss_prop%phase_fnc_clr                        &
      , ss_prop%tau, ss_prop%tau_dir                                    &
      , ss_prop%omega, ss_prop%phase_fnc                                &
      , tau_clr_f, tau_clr_dir_f, omega_clr_f, phase_fnc_clr_f       &
!                   Sizes of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, 1                    &
      )

    CALL two_stream(ierr                                                &
      , control, bound                                                  &
!                   Atmospheric properties
      , n_profile, n_layer                                              &
!                   Two-stream scheme
      , i_2stream                                                       &
!                   Options for solver
      , i_solver_clear, i_scatter_method                                &
!                   Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                   Spectral region
      , isolir                                                          &
!                   Infra-red properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                   Conditions at TOA
      , flux_inc_down, flux_inc_direct, sec_0                           &
!                   Surface conditions
      , rho_alb(1, ip_surf_alb_diff)                                    &
      , rho_alb(1, ip_surf_alb_dir), d_planck_flux_surface              &
!                   Spherical geometry
      , sph                                                             &
!                   Single scattering properties
      , tau_clr_dir_f, tau_clr_f, omega_clr_f                           &
      , phase_fnc_clr_f(1, 1, 1)                                        &
!                   Fluxes calculated
      , flux_direct_clear, flux_total_clear                             &
!                   Sizes of arrays
      , nd_profile, nd_layer, nd_source_coeff                           &
      )

!   Remove the arrays that are no longer required.
    DEALLOCATE(tau_clr_f)
    DEALLOCATE(tau_clr_dir_f)
    DEALLOCATE(omega_clr_f)
    DEALLOCATE(phase_fnc_clr_f)

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_flux_ipa
END MODULE calc_flux_ipa_mod
