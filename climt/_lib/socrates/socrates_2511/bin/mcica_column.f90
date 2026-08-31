! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine to solve the two-stream equations for completely
!  overcast columns and completely clear-sky columns
!
! Method:
!   The two-stream coefficients are calculated for the clear-sky
!   case and the overcast case. From these clear and cloudy
!   transmission and reflection coefficients are determined. For
!   each case (cloudy and clear) a suitable solver is called.
!
!- ---------------------------------------------------------------------
MODULE mcica_column_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MCICA_COLUMN_MOD'
CONTAINS
SUBROUTINE mcica_column(ierr                                            &
  , control, cld, bound                                                 &
!                   Atmospheric properties
  , n_profile, n_layer                                                  &
!                   Two-stream scheme
  , i_2stream                                                           &
!                   Options for solver
  , i_solver, i_scatter_method                                          &
!                   Options for equivalent extinction
  , l_scale_solar, adjust_solar_ke                                      &
!                   Spectral region
  , isolir                                                              &
!                   Infra-red properties
  , diff_planck                                                         &
  , l_ir_source_quad, diff_planck_2                                     &
!                   Conditions at TOA
  , flux_inc_down, flux_inc_direct, sec_0                               &
!                   Conditions at surface
  , diffuse_albedo, direct_albedo, d_planck_flux_surface                &
!                   Spherical geometry
  , sph                                                                 &
!                   Optical Properties
  , ss_prop                                                             &
!                   Cloud geometry
  , n_cloud_top, index_subcol                                           &
!                   Calculated fluxes
  , flux_direct, flux_total, l_actinic, actinic_flux                    &
!                   Flags for clear-sky calculations
  , l_clear                                                             &
!                   Calculated clear-sky fluxes
  , flux_direct_clear, flux_total_clear, actinic_flux_clear             &
!                   Dimensions of arrays
  , nd_profile, nd_layer, nd_layer_clr, id_ct                           &
  , nd_source_coeff                                                     &
  , nd_cloud_type                                                       &
  )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_cld,     ONLY: StrCld
  USE def_bound,   ONLY: StrBound
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_direct_csr_scaling, ip_direct_noscaling,        &
                     ip_infra_red, ip_no_scatter_abs, ip_no_scatter_ext,&
                     ip_scatter_approx, ip_solar,  ip_scatter_full
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  USE set_n_source_coeff_mod, ONLY: set_n_source_coeff
  USE calc_actinic_flux_mod, ONLY: calc_actinic_flux
  USE column_solver_mod, ONLY: column_solver
  USE ir_source_mod, ONLY: ir_source
  USE two_coeff_fast_lw_mod, ONLY: two_coeff_fast_lw
  USE two_coeff_mod, ONLY: two_coeff

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Cloud properties:
  TYPE(StrCld),       INTENT(IN)    :: cld

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
    nd_profile                                                          &
!       Size allocated for atmospheric profiles
  , nd_layer                                                            &
!       Size allocated for atmospheric layers
  , nd_layer_clr                                                        &
!       Size allocated for completely clear layers
  , id_ct                                                               &
!       Topmost declared cloudy layer
  , nd_source_coeff                                                     &
!       Size allocated for coefficients in the source function
  , nd_cloud_type
!       Size allocated for types of clouds


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
    n_profile                                                           &
!       Number of profiles
  , n_layer                                                             &
!       Number of layers
  , n_cloud_top                                                         &
!       Top cloudy layer
  , index_subcol                                                        &
!       Index of current sub-grid cloud column
  , isolir                                                              &
!       Spectral region
  , i_2stream                                                           &
!       Two-stream scheme
  , i_solver                                                            &
!       Solver used
  , i_scatter_method
!       Method of treating scattering
  INTEGER, INTENT(INOUT) ::                                             &
    ierr
!       Error flag
  LOGICAL, INTENT(IN) ::                                                &
    l_clear                                                             &
!       Calculate clear-sky fluxes
  , l_scale_solar                                                       &
!       Flag to scale solar
  , l_ir_source_quad
!       Use quadratic source term

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!       Spherical geometry fields

! Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

  REAL (RealK), INTENT(IN) ::                                           &
    sec_0(nd_profile)                                                   &
!       Secant of solar zenith angle
  , diffuse_albedo(nd_profile)                                          &
!       Diffuse albedo
  , direct_albedo(nd_profile)                                           &
!       Direct albedo
  , flux_inc_down(nd_profile)                                           &
!       Incident total flux
  , flux_inc_direct(nd_profile)                                         &
!       Incident direct flux
  , diff_planck(nd_profile, nd_layer)                                   &
!       Change in Planckian function
  , d_planck_flux_surface(nd_profile)                                   &
!       Flux from surface
  , adjust_solar_ke(nd_profile, nd_layer)                               &
!       Adjustment of solar beam with equivalent extinction
  , diff_planck_2(nd_profile, nd_layer)
!         2x2nd difference of Planckian

! Fluxes calculated
  REAL (RealK), INTENT(OUT) ::                                          &
    flux_direct(nd_profile, 0: nd_layer)                                &
!       Direct flux
  , flux_total(nd_profile, 2*nd_layer+2)                                &
!       Long flux vector
  , actinic_flux(nd_profile, nd_layer)                                  &
!       Actinic flux
  , flux_direct_clear(nd_profile, 0: nd_layer)                          &
!       Clear direct flux
  , flux_total_clear(nd_profile, 2*nd_layer+2)                          &
!       Clear total flux
  , actinic_flux_clear(nd_profile, nd_layer)
!       Clear actinic flux
  LOGICAL, INTENT(IN) :: l_actinic
!       Actinic fluxes calculated


! Local variabales.
  INTEGER ::                                                            &
    n_source_coeff                                                      &
!       Number of source coefficients
  , i, ii, j, k, l, ll
!       Loop variables


! Clear-sky coefficients:
  REAL (RealK) ::                                                       &
    trans(nd_profile, nd_layer)                                         &
!       Free transmission of layer
  , reflect(nd_profile, nd_layer)                                       &
!       Free reflectance of layer
  , trans_0(nd_profile, nd_layer)                                       &
!       Free direct transmission of layer
  , trans_0_dir(nd_profile, nd_layer)                                   &
!       Free direct transmission of layer using direct tau
  , source_coeff(nd_profile, nd_layer, nd_source_coeff)                 &
!       Free source coefficients
  , s_down(nd_profile, nd_layer)                                        &
!       Free downward source
  , s_up(nd_profile, nd_layer)
!       Free upward source


! Coefficients in the two-stream equations:
  REAL (RealK) ::                                                       &
    trans_temp(nd_profile, 1)                                           &
!       Temporary diffuse transmission coefficient
  , reflect_temp(nd_profile, 1)                                         &
!       Temporary diffuse reflection coefficient
  , trans_0_temp(nd_profile, 1)                                         &
!       Temporary direct transmission coefficient
  , trans_0_temp_dir(nd_profile, 1)                                     &
!       Temporary unscaled direct transmission coefficient
  , source_coeff_temp(nd_profile, 1, nd_source_coeff)
!       Temporary source coefficients in two-stream equations

! Variables for gathering:
  INTEGER                                                               &
    n_list(nd_layer, nd_cloud_type)                                     &
!       Number of points in list
  , l_list(nd_profile, nd_layer, nd_cloud_type)                         &
!       List of collected points
  , n_list_sph(nd_layer, nd_cloud_type)                                 &
!       Number of points in list for spherical geometry
  , l_list_sph(nd_profile, nd_layer, nd_cloud_type)
!       List of collected points for spherical geometry
  REAL (RealK) ::                                                       &
    tau_gathered(nd_profile, 1)                                         &
!       Gathered optical depth
  , tau_gathered_dir(nd_profile, 1)                                     &
!       Optical depth for direct flux
  , omega_gathered(nd_profile, 1)                                       &
!       Gathered alebdo of single scattering
  , asymmetry_gathered(nd_profile, 1)                                   &
!       Gathered asymmetry
  , sec_0_gathered(nd_profile)                                          &
!       Gathered secant of the solar zenith angle
  , path_div_gathered(nd_profile, 1)
!       Gathered path scaling for calculating direct flux divergence

! Variables for actinic flux calculation
  LOGICAL :: l_mask(nd_profile, nd_layer)
!       Mask of cloudy points over all cloud types
  REAL(RealK), ALLOCATABLE :: tau_abs(:, :)
!       Temporary array for the optical depth due to absorption

  INTEGER :: path_base

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'MCICA_COLUMN'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the number of source coefficients for the approximation
  n_source_coeff=set_n_source_coeff(isolir, l_ir_source_quad)

! Calculate the transmission and reflection coefficients and
! source terms for the clear sky
  IF ( (i_scatter_method == ip_scatter_full) .OR.                       &
       (i_scatter_method == ip_scatter_approx) ) THEN
    CALL two_coeff(ierr, control                                        &
    , n_profile, 1, n_cloud_top-1                                       &
    , i_2stream                                                         &
    , ss_prop%phase_fnc_clr                                             &
    , ss_prop%omega_clr, ss_prop%tau_clr_dir, ss_prop%tau_clr           &
    , isolir, sec_0, sph%common%path_div                                &
    , trans, reflect, trans_0_dir, trans_0                              &
    , source_coeff                                                      &
    , nd_profile, 1, nd_layer_clr, 1, nd_layer, nd_source_coeff         &
    )
    CALL two_coeff(ierr, control                                        &
    , n_profile, n_cloud_top, n_layer                                   &
    , i_2stream                                                         &
    , ss_prop%phase_fnc(:, :, :, 0)                                     &
    , ss_prop%omega(:, :, 0)                                            &
    , ss_prop%tau_dir(:, :, 0), ss_prop%tau(:, :, 0)                    &
    , isolir, sec_0, sph%common%path_div                                &
    , trans, reflect, trans_0_dir, trans_0                              &
    , source_coeff                                                      &
    , nd_profile, id_ct, nd_layer, 1, nd_layer, nd_source_coeff         &
    )
  ELSE IF ( (i_scatter_method == ip_no_scatter_abs) .OR.                &
            (i_scatter_method == ip_no_scatter_ext) ) THEN
    CALL two_coeff_fast_lw(n_profile, 1, n_cloud_top-1                  &
      , l_ir_source_quad, ss_prop%tau_clr                               &
      , trans, source_coeff                                             &
      , nd_profile, nd_layer, 1, nd_layer_clr, nd_source_coeff)
    CALL two_coeff_fast_lw(n_profile, n_cloud_top, n_layer              &
      , l_ir_source_quad, ss_prop%tau(:, :, 0)                          &
      , trans, source_coeff                                             &
      , nd_profile, nd_layer, id_ct, nd_layer, nd_source_coeff)
    DO i=1, n_layer
      DO l=1, n_profile
        reflect(l, i)=0.0e+00_RealK
      END DO
    END DO
  END IF

! Calculate the optical depth to absorption for the actinic flux
  IF (l_actinic) THEN
    ALLOCATE(tau_abs(nd_profile, nd_layer))
    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        tau_abs(l, i) = &
          ss_prop%tau_clr(l, i) * ( 1.0_RealK - ss_prop%omega_clr(l, i) )
      END DO
    END DO
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        tau_abs(l, i) = &
          ss_prop%tau(l, i, 0) * ( 1.0_RealK - ss_prop%omega(l, i, 0) )
      END DO
    END DO
  END IF

! Calculate clear-sky fluxes before optical properties are overwritten
! for cloud points.
  IF (l_clear) THEN
    IF (isolir == ip_infra_red) THEN
      CALL ir_source(n_profile, 1, n_layer                              &
      , source_coeff, diff_planck                                       &
      , l_ir_source_quad, diff_planck_2                                 &
      , s_down, s_up                                                    &
      , nd_profile, nd_layer, nd_source_coeff)
    END IF

    CALL column_solver(ierr, control, bound, sph%common, sph%clear      &
    , n_profile, n_layer                                                &
    , i_scatter_method, i_solver                                        &
    , trans, reflect, trans_0_dir, trans_0, source_coeff                &
    , isolir, flux_inc_direct, flux_inc_down                            &
    , s_down, s_up                                                      &
    , diffuse_albedo, direct_albedo                                     &
    , d_planck_flux_surface                                             &
    , l_scale_solar, adjust_solar_ke                                    &
    , flux_direct_clear, flux_total_clear                               &
    , nd_profile, nd_layer, nd_source_coeff)

!   Calculate the clear-sky actinic flux
    IF (l_actinic) THEN
      CALL calc_actinic_flux(control, sph%clear, sph%common, &
        n_profile, n_layer, tau_abs, &
        flux_total_clear, flux_direct_clear, sec_0, &
        l_scale_solar, adjust_solar_ke, &
        actinic_flux_clear, &
        nd_profile, nd_layer)
    END IF
  END IF


! Repeat the calculation for cloudy regions.
! (Clouds are indexed beginning with index 1 in the last
! dimension of arrays of optical properties.)

! Zero points where cloud exists so they can be overwritten by the sum
! of the cloudy contributions (which already include clear-sky)
  l_mask = .FALSE.
  IF (isolir == ip_solar) THEN
    DO k=1, cld%n_cloud_type
      DO i=n_cloud_top, n_layer
        n_list(i,k)=0
        DO l=1, n_profile
          IF (cld%c_sub(l,i,index_subcol,k)*cld%frac_cloud(l,i,k)       &
              > 0.0_RealK) THEN
            n_list(i,k)=n_list(i,k)+1
            l_list(n_list(i,k),i,k)=l
            trans(l, i)=0.0_RealK
            reflect(l, i)=0.0_RealK
            trans_0(l, i)=0.0_RealK
            trans_0_dir(l, i)=0.0_RealK
            IF (l_actinic) THEN
              l_mask(l, i) = .TRUE.
              tau_abs(l, i)=0.0_RealK
            END IF
            DO j=1, n_source_coeff
              source_coeff(l, i, j)=0.0_RealK
            END DO
          END IF
        END DO
      END DO
    END DO
  ELSE
    DO k=1, cld%n_cloud_type
      DO i=n_cloud_top, n_layer
        n_list(i,k)=0
        DO l=1, n_profile
          IF (cld%c_sub(l,i,index_subcol,k)*cld%frac_cloud(l,i,k)       &
              > 0.0_RealK) THEN
            n_list(i,k)=n_list(i,k)+1
            l_list(n_list(i,k),i,k)=l
            trans(l, i)=0.0_RealK
            reflect(l, i)=0.0_RealK
            IF (l_actinic) THEN
              l_mask(l, i) = .TRUE.
              tau_abs(l, i)=0.0_RealK
            END IF
            DO j=1, n_source_coeff
              source_coeff(l, i, j)=0.0_RealK
            END DO
          END IF
        END DO
      END DO
    END DO
  END IF

! Calculate the transmission and reflection coefficients for
! each type of cloud and increment the totals, weighting with
! the cloud fraction.
  DO k=1, cld%n_cloud_type
    DO i=n_cloud_top, n_layer
      IF (n_list(i,k) >  0) THEN

!       Gather the optical properties.
!       Here we must consider one layer at a time. To reduce
!       storage the temporary arrays are only one layer thick,
!       but they will be passed to the subroutine where they
!       will be declared as running from the Ith to the Ith layer
!       to make the code more readable at the lower level.
        IF ( (i_scatter_method == ip_scatter_full) .OR.                 &
             (i_scatter_method == ip_scatter_approx) ) THEN
          DO l=1, n_list(i,k)
            ll=l_list(l,i,k)
            tau_gathered(l, 1)=ss_prop%tau(ll, i, k)
            omega_gathered(l, 1)=ss_prop%omega(ll, i, k)
            asymmetry_gathered(l, 1)=ss_prop%phase_fnc(ll, i, 1, k)
          END DO
          IF (control%i_direct_tau == ip_direct_noscaling .OR.          &
              control%i_direct_tau == ip_direct_csr_scaling) THEN
            DO l=1, n_list(i,k)
              ll=l_list(l,i,k)
              tau_gathered_dir(l, 1)=ss_prop%tau_dir(ll, i, k)
            END DO
          END IF

          IF (isolir == ip_solar) THEN
            IF (control%l_spherical_solar) THEN
              DO l=1, n_list(i,k)
                ll=l_list(l,i,k)
                path_div_gathered(l,1)=sph%common%path_div(ll,i)
              END DO
            ELSE
              DO l=1, n_list(i,k)
                ll=l_list(l,i,k)
                sec_0_gathered(l)=sec_0(ll)
              END DO
            END IF
          END IF

          CALL two_coeff(ierr, control                                  &
            , n_list(i,k), i, i                                         &
            , i_2stream                                                 &
            , asymmetry_gathered, omega_gathered                        &
            , tau_gathered_dir, tau_gathered                            &
            , isolir, sec_0_gathered, path_div_gathered                 &
            , trans_temp, reflect_temp                                  &
            , trans_0_temp_dir, trans_0_temp                            &
            , source_coeff_temp                                         &
            , nd_profile, i, i, i, i, nd_source_coeff                   &
            )

          DO l=1, n_list(i,k)
            ll=l_list(l,i,k)
            trans(ll, i)=trans(ll, i)                                   &
              +cld%frac_cloud(ll, i, k)*trans_temp(l, 1)
            reflect(ll, i)=reflect(ll, i)                               &
              +cld%frac_cloud(ll, i, k)*reflect_temp(l, 1)
          END DO
          DO j=1, n_source_coeff
            DO l=1, n_list(i,k)
              ll=l_list(l,i,k)
              source_coeff(ll, i, j)=source_coeff(ll, i, j)             &
                +cld%frac_cloud(ll, i, k)                               &
                *source_coeff_temp(l, 1, j)
            END DO
          END DO
          IF (isolir == ip_solar) THEN
            DO l=1, n_list(i,k)
              ll=l_list(l,i,k)
              trans_0(ll, i)=trans_0(ll, i)                             &
                +cld%frac_cloud(ll, i, k)                               &
                *trans_0_temp(l, 1)
            END DO
            IF (control%i_direct_tau == ip_direct_noscaling .OR.        &
                control%i_direct_tau == ip_direct_csr_scaling) THEN
              DO l=1, n_list(i,k)
                ll=l_list(l,i,k)
                trans_0_dir(ll, i)=trans_0_dir(ll, i)                   &
                  +cld%frac_cloud(ll, i, k)                             &
                  *trans_0_temp_dir(l, 1)
              END DO
            END IF
          END IF
          IF (l_actinic) THEN
            DO l=1, n_list(i,k)
              ll=l_list(l,i,k)
              ! Calculate effective optical depth to absorption by
              ! assuming the flux divergence is distributed
              ! proportional to 1-trans
              tau_abs(ll, i) = tau_abs(ll,i) + cld%frac_cloud(ll, i, k) &
                * (1.0_RealK - trans_temp(l, 1)) / ( tau_gathered(l, 1) &
                * (1.0_RealK - omega_gathered(l, 1)) )
            END DO
          END IF
        ELSE IF ( (i_scatter_method == ip_no_scatter_abs) .OR.          &
                  (i_scatter_method == ip_no_scatter_ext) ) THEN
          DO l=1, n_list(i,k)
            tau_gathered(l, 1)=ss_prop%tau(l_list(l,i,k), i, k)
          END DO
          CALL two_coeff_fast_lw(n_list(i,k), 1, 1                      &
            , l_ir_source_quad, tau_gathered                            &
            , trans_temp, source_coeff_temp                             &
            , nd_profile, 1, 1, 1, nd_source_coeff)
          DO l=1, n_list(i,k)
            ll=l_list(l,i,k)
            trans(ll, i)=trans(ll, i)                                   &
              +cld%frac_cloud(ll, i, k)*trans_temp(l, 1)
          END DO
          DO j=1, n_source_coeff
            DO l=1, n_list(i,k)
              ll=l_list(l,i,k)
              source_coeff(ll, i, j)=source_coeff(ll, i, j)             &
                +cld%frac_cloud(ll, i, k)*source_coeff_temp(l, 1, j)
            END DO
          END DO
          IF (l_actinic) THEN
            DO l=1, n_list(i,k)
              ll=l_list(l,i,k)
              tau_abs(ll, i) = tau_abs(ll,i) + cld%frac_cloud(ll, i, k) &
                * (1.0_RealK - trans_temp(l, 1)) / tau_gathered(l, 1)
            END DO
          END IF
        END IF

      END IF

    END DO
  END DO
  IF (l_actinic) THEN
    WHERE (l_mask)
      ! Convert weighted sum to absorption optical depth for actinic flux
      tau_abs = (1.0_RealK - trans) / tau_abs
    END WHERE
  END IF

! Calculate the transmission through cloudy layers for spherical geometry
  IF (control%l_spherical_solar) THEN
    DO ii=0, n_layer+1
      ! First zero the clear-sky transmission where there is cloud
      ! and where the spherical path touches the layer
      DO k=1, cld%n_cloud_type
        DO i=n_cloud_top, n_layer
          n_list_sph(i,k)=0
          DO l=1, n_list(i,k)
            ll=l_list(l,i,k)
            IF (sph%common%path(ll,i,ii) > 0.0_RealK) THEN
              n_list_sph(i,k)=n_list_sph(i,k)+1
              l_list_sph(n_list_sph(i,k),i,k)=ll
              sph%common%trans_0_cloud(ll,i,ii)=0.0_RealK
            END IF
          END DO
        END DO
      END DO
      ! Now add transmission through cloud fractions (note that
      ! frac_cloud adds up to 1 for MCICA sub-columns)
      DO k=1, cld%n_cloud_type
        DO i=n_cloud_top, n_layer
          DO l=1, n_list_sph(i,k)
            ll=l_list_sph(l,i,k)
            !!!!! Optimise later through use of gathering and exp_v
            sph%common%trans_0_cloud(ll,i,ii) &
              = sph%common%trans_0_cloud(ll,i,ii) &
              + cld%frac_cloud(ll, i, k) &
              *EXP(-sph%common%path(ll,i,ii)*ss_prop%tau(ll,i,k))
          END DO
        END DO
      END DO
      ! Calculate total transmission through cloudy layers and the
      ! layers above cloud top (calculated in spherical_trans_coeff).
      DO l=1, n_profile
        path_base = sph%common%path_base(l, ii)
        sph%allsky%trans_0(l, ii) = sph%allsky%trans_0(l, ii) * PRODUCT( &
          sph%common%trans_0_cloud(l,n_cloud_top:path_base,ii) )
      END DO
    END DO
  END IF


! Now calculate the fluxes including clouds
  IF (isolir == ip_infra_red) THEN
    CALL ir_source(n_profile, 1, n_layer                                &
    , source_coeff, diff_planck                                         &
    , l_ir_source_quad, diff_planck_2                                   &
    , s_down, s_up                                                      &
    , nd_profile, nd_layer, nd_source_coeff)
  END IF

  CALL column_solver(ierr, control, bound, sph%common, sph%allsky       &
    , n_profile, n_layer                                                &
    , i_scatter_method, i_solver                                        &
    , trans, reflect, trans_0_dir, trans_0, source_coeff                &
    , isolir, flux_inc_direct, flux_inc_down                            &
    , s_down, s_up                                                      &
    , diffuse_albedo, direct_albedo                                     &
    , d_planck_flux_surface                                             &
    , l_scale_solar, adjust_solar_ke                                    &
    , flux_direct, flux_total                                           &
    , nd_profile, nd_layer, nd_source_coeff)

  IF (l_actinic) THEN
    CALL calc_actinic_flux(control, sph%allsky, sph%common, &
      n_profile, n_layer, tau_abs, flux_total, flux_direct, sec_0, &
      l_scale_solar, adjust_solar_ke, &
      actinic_flux, &
      nd_profile, nd_layer)
    DEALLOCATE(tau_abs)
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE mcica_column
END MODULE mcica_column_mod
