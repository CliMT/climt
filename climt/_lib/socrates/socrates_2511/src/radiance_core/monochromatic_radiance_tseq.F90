! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to solve for the monochromatic radiances.
!
! Method:
!   The final single scattering properties are calculated
!   and rescaled. An appropriate subroutine is called to
!   calculate the radiances depending on the treatment of
!   cloudiness.
!
!- ---------------------------------------------------------------------
MODULE monochromatic_radiance_tseq_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MONOCHROMATIC_RADIANCE_TSEQ_MOD'
CONTAINS
SUBROUTINE monochromatic_radiance_tseq(ierr                             &
    , control, cld, bound                                               &
!                 Atmospheric Propertries
    , n_profile, n_layer                                                &
!                 Options for Solver
    , i_2stream, i_solver, i_scatter_method                             &
!                 Optical Properties
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , diff_planck, l_ir_source_quad, diff_planck_2                      &
!                 Conditions at TOA
    , sec_0, flux_inc_direct, flux_inc_down                             &
!                 Surface Properties
    , d_planck_flux_surface                                             &
    , rho_alb                                                           &
!                 Spherical geometry
    , sph                                                               &
!                 Optical Properties
    , ss_prop                                                           &
!                 Cloudy Properties
    , i_cloud                                                           &
!                 Cloud Geometry
    , n_cloud_top, index_subcol                                         &
    , n_region, k_clr, i_region_cloud, frac_region                      &
    , w_free, cloud_overlap                                             &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                 Fluxes Calculated
    , flux_direct, flux_total, l_actinic, actinic_flux                  &
!                 Flags for Clear-sky Calculation
    , l_clear, i_solver_clear                                           &
!                 Clear-sky Fluxes Calculated
    , flux_direct_clear, flux_total_clear, actinic_flux_clear           &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_cloud_type, nd_region, nd_overlap_coeff                        &
    , nd_source_coeff, nd_max_order                                     &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_cld,     ONLY: StrCld
  USE def_bound,   ONLY: StrBound
  USE def_ss_prop, ONLY: str_ss_prop
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: ip_cloud_clear, ip_surf_alb_diff, ip_surf_alb_dir, &
                     ip_cloud_mcica, ip_cloud_triple, ip_solar,         &
                     ip_cloud_mix_max, ip_cloud_mix_random,             &
                     ip_cloud_part_corr, ip_cloud_part_corr_cnv,        &
                     ip_cloud_column_max
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE calc_actinic_flux_mod, ONLY: calc_actinic_flux
  USE calc_flux_ipa_mod, ONLY: calc_flux_ipa
  USE copy_clr_full_mod, ONLY: copy_clr_full
  USE mcica_column_mod, ONLY: mcica_column
  USE mix_column_mod, ONLY: mix_column
  USE triple_column_mod, ONLY: triple_column
  USE two_stream_mod, ONLY: two_stream

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Cloud properties:
  TYPE(StrCld),       INTENT(IN)    :: cld

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_layer_clr                                                      &
!       Size allocated for completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_column                                                         &
!       Number of columns per point
    , nd_cloud_type                                                     &
!       Maximum number of types of cloud
    , nd_region                                                         &
!       Maximum number of cloudy regions
    , nd_overlap_coeff                                                  &
!       Maximum number of overlap coefficients
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_max_order
!       Size allocated for spherical harmonics


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

!                 Atmospheric properties
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers

!                 Angular integration
  INTEGER, INTENT(IN) ::                                                &
      i_2stream                                                         &
!       Two-stream scheme
    , i_scatter_method
!       Method of treating scattering

!                 Options for solver
  INTEGER, INTENT(IN) ::                                                &
      i_solver
!       Solver used

!                 Variables for equivalent extinction
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Apply scaling to solar flux
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment of solar beam with equivalent extinction

!                 Spectral region
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Visible or IR

!                 Infra-red properties
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Flag for quadratic IR-source
  REAL (RealK), INTENT(IN) ::                                           &
      diff_planck(nd_profile, nd_layer)                                 &
!       DIfferences in the Planckian function across layers
    , diff_planck_2(nd_profile, nd_layer)
!       Twice the second differences of Planckian source function

!                 Conditions at TOA
  REAL (RealK), INTENT(IN) ::                                           &
      sec_0(nd_profile)                                                 &
!       Secants (two-stream) or cosines (spherical harmonics)
!       of the solar zenith angles
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)
!       Incident downward flux

!                 Surface properties
  REAL (RealK), INTENT(IN) ::                                           &
      d_planck_flux_surface(nd_profile)
!       Differential Planckian flux from the surface
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb(nd_profile, 2)
!       Weights of the basis functions

  TYPE(StrSphGeo), INTENT(INOUT) :: sph
!       Spherical geometry fields

!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

!                 Cloudy properties
  INTEGER, INTENT(IN) ::                                                &
      i_cloud
!       Cloud scheme used

!                 Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top                                                       &
!       Topmost cloudy layer
    , index_subcol                                                      &
!       Index of current sub-grid cloud column
    , n_region                                                          &
!       Number of cloudy regions
    , i_region_cloud(nd_cloud_type)                                     &
!       Regions in which types of clouds fall
    , k_clr
!       Index of clear-sky region

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
      w_free(nd_profile, id_ct: nd_layer)                               &
!       Clear-sky fraction
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)                                             &
!       Coefficients for energy transfer at interfaces
    , area_column(nd_profile, nd_column)                                &
!       Areas of columns
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region


!                 Calculated Fluxes
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_total(nd_profile, 2*nd_layer+2)                              &
!       Total flux
    , actinic_flux(nd_profile, nd_layer)
!       Actinic flux
  LOGICAL, INTENT(IN) :: l_actinic
!       Flag for calculation of actinic flux

!                 Flags for clear-sky calculations
  LOGICAL, INTENT(IN) ::                                                &
      l_clear
!       Calculate clear-sky properties
  INTEGER, INTENT(IN) ::                                                &
      i_solver_clear
!       Clear solver used

!                 Clear-sky fluxes calculated
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct_clear(nd_profile, 0: nd_layer)                        &
!       Clear-sky direct flux
    , flux_total_clear(nd_profile, 2*nd_layer+2)                        &
!       Clear-sky total flux
    , actinic_flux_clear(nd_profile, nd_layer)
!       Clear-sky actinic flux



! Local variables.
  INTEGER                                                               &
      nd_profile_column
!       Size allowed for profiles taken simultaneously in a
!       decomposition into columns
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , i
!       Loop variable
! Full clear-sky single-scattering properties
  REAL (RealK), ALLOCATABLE ::                                          &
      tau_clr_f(:, :)                                                   &
!       Clear-sky optical depth
    , tau_clr_dir_f(:, :)                                               &
!       Unscaled clear-sky optical depth
    , omega_clr_f(:, :)                                                 &
!       Clear-sky albedo of single scattering
    , phase_fnc_clr_f(:, :, :)
!       Clear-sky phase function
  REAL(RealK), ALLOCATABLE :: tau_abs(:, :)
!       Temporary array for the optical depth due to absorption

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MONOCHROMATIC_RADIANCE_TSEQ'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Choose an appropriate routine to calculate the fluxes as
! determined by the cloud scheme selected.

  IF (i_cloud == ip_cloud_clear) THEN
!   Allocate and set dynamic arrays.
    ALLOCATE(tau_clr_f(nd_profile, nd_layer))
    ALLOCATE(tau_clr_dir_f(nd_profile, nd_layer))
    ALLOCATE(omega_clr_f(nd_profile, nd_layer))
    ALLOCATE(phase_fnc_clr_f(nd_profile, nd_layer, 1))

    CALL copy_clr_full(n_profile, n_layer, n_cloud_top                  &
      , control, 1                                                      &
      , ss_prop%tau_clr, ss_prop%tau_clr_dir, ss_prop%omega_clr         &
      , ss_prop%phase_fnc_clr                                           &
      , ss_prop%tau, ss_prop%tau_dir                                    &
      , ss_prop%omega, ss_prop%phase_fnc                                &
      , tau_clr_f, tau_clr_dir_f, omega_clr_f, phase_fnc_clr_f          &
!                   Sizes of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, 1                    &
      )

!   A two-stream scheme with no clouds.
    CALL two_stream(ierr                                                &
      , control, bound                                                  &
!                 Atmospheric properties
      , n_profile, n_layer                                              &
!                 Two-stream scheme
      , i_2stream                                                       &
!                 Options for solver
      , i_solver, i_scatter_method                                      &
!                 Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                 Spectral region
      , isolir                                                          &
!                 Infra-red properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                 Conditions at TOA
      , flux_inc_down, flux_inc_direct, sec_0                           &
!                 Surface conditions
      , rho_alb(1, ip_surf_alb_diff)                                    &
      , rho_alb(1, ip_surf_alb_dir), d_planck_flux_surface              &
!                 Spherical geometry
      , sph                                                             &
!                 Single scattering properties
      , tau_clr_dir_f, tau_clr_f                                        &
      , omega_clr_f, phase_fnc_clr_f(1, 1, 1)                           &
!                 Fluxes calculated
      , flux_direct, flux_total                                         &
!                 Sizes of arrays
      , nd_profile, nd_layer, nd_source_coeff                           &
      )

!   Calculate the actinic flux
    IF (l_actinic) THEN
      ALLOCATE(tau_abs(nd_profile, nd_layer))
      DO i=1, n_layer
        DO l=1, n_profile
          tau_abs(l, i) = tau_clr_f(l, i)*(1.0_RealK-omega_clr_f(l, i))
        END DO
      END DO
      CALL calc_actinic_flux(control, sph%allsky, sph%common, &
        n_profile, n_layer, tau_abs, flux_total, flux_direct, sec_0, &
        l_scale_solar, adjust_solar_ke, &
        actinic_flux, &
        nd_profile, nd_layer)
      DEALLOCATE(tau_abs)
    END IF

!   Release temporary storage.
    DEALLOCATE(tau_clr_f)
    DEALLOCATE(tau_clr_dir_f)
    DEALLOCATE(omega_clr_f)
    DEALLOCATE(phase_fnc_clr_f)

    IF (l_clear) THEN
!     The clear fluxes here can be copied directly without
!     any further calculation.
      IF (isolir == ip_solar) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            flux_direct_clear(l, i)=flux_direct(l, i)
          END DO
        END DO
        IF (control%l_spherical_solar) THEN
          DO i=0, n_layer+1
            DO l=1, n_profile
              sph%clear%flux_direct(l, i) &
                = sph%allsky%flux_direct(l, i)
            END DO
          END DO
          DO i=1, n_layer
            DO l=1, n_profile
              sph%clear%flux_direct_div(l, i) &
                = sph%allsky%flux_direct_div(l, i)
            END DO
          END DO
        END IF
      END IF
      DO i=1, 2*n_layer+2
        DO l=1, n_profile
          flux_total_clear(l, i)=flux_total(l, i)
        END DO
      END DO
      IF (l_actinic) THEN
        DO i=1, n_layer
          DO l=1, n_profile
            actinic_flux_clear(l, i)=actinic_flux(l, i)
          END DO
        END DO
      END IF
    END IF

  ELSE IF (i_cloud == ip_cloud_mcica) THEN


    CALL mcica_column(ierr                                              &
      , control, cld, bound                                             &
!                 Atmospheric properties
      , n_profile, n_layer                                              &
!                 Two-stream scheme
      , i_2stream                                                       &
!                 Options for solver
      , i_solver, i_scatter_method                                      &
!                 Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                 Spectral region
      , isolir                                                          &
!                 Infra-red properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                 Conditions at TOA
      , flux_inc_down, flux_inc_direct, sec_0                           &
!                 Conditions at surface
      , rho_alb(1, ip_surf_alb_diff), rho_alb(1, ip_surf_alb_dir)       &
      , d_planck_flux_surface                                           &
!                 Spherical geometry
      , sph                                                             &
!                 Single scattering properties
      , ss_prop                                                         &
!                 Cloud geometry
      , n_cloud_top, index_subcol                                       &
!                 Fluxes calculated
      , flux_direct, flux_total, l_actinic, actinic_flux                &
!                 Flags for clear-sky calculations
      , l_clear                                                         &
!                 Clear-sky fluxes calculated
      , flux_direct_clear, flux_total_clear, actinic_flux_clear         &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct                       &
      , nd_source_coeff                                                 &
      , nd_cloud_type                                                   &
      )

  ELSE IF ((i_cloud == ip_cloud_mix_max).OR.                            &
           (i_cloud == ip_cloud_mix_random).OR.                         &
           ( (i_cloud == ip_cloud_part_corr).AND.                       &
             (n_region == 2) ) ) THEN

!   Clouds are treated using the coupled overlaps originally
!   introduced by Geleyn and Hollingsworth.

    CALL mix_column(ierr                                                &
      , control, bound                                                  &
!                 Atmospheric properties
      , n_profile, n_layer, k_clr                                       &
!                 Two-stream scheme
      , i_2stream                                                       &
!                 Options for solver
      , i_solver, i_scatter_method                                      &
!                 Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                 Spectral region
      , isolir                                                          &
!                 Infra-red properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                 Conditions at TOA
      , flux_inc_down, flux_inc_direct, sec_0                           &
!                 Conditions at surface
      , rho_alb(1, ip_surf_alb_diff), rho_alb(1, ip_surf_alb_dir)       &
      , d_planck_flux_surface                                           &
!                 Spherical geometry
      , sph                                                             &
!                 Single scattering properties
      , ss_prop                                                         &
!                 Cloud geometry
      , n_cloud_top                                                     &
      , cld%n_cloud_type, cld%frac_cloud                                &
      , w_free, cld%w_cloud                                             &
      , cloud_overlap                                                   &
!                 Fluxes calculated
      , flux_direct, flux_total                                         &
!                 Flags for clear-sky calculations
      , l_clear, i_solver_clear                                         &
!                 Clear-sky fluxes calculated
      , flux_direct_clear, flux_total_clear                             &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct                       &
      , nd_max_order, nd_source_coeff                                   &
      , nd_cloud_type, nd_overlap_coeff                                 &
      )

  ELSE IF ((i_cloud == ip_cloud_triple).OR.                             &
           ( (i_cloud == ip_cloud_part_corr_cnv).AND.                   &
             (n_region == 3) ) ) THEN

!   Clouds are treated using a decomposition of the column
!   into clear-sky, stratiform and convective regions.

    CALL triple_column(ierr                                             &
      , control, bound                                                  &
!                 Atmospheric properties
      , n_profile, n_layer                                              &
!                 Two-stream scheme
      , i_2stream                                                       &
!                 Options for solver
      , i_solver, i_scatter_method                                      &
!                 Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                 Spectral region
      , isolir                                                          &
!                 Infra-red properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                 Conditions at TOA
      , flux_inc_down, flux_inc_direct, sec_0                           &
!                 Conditions at surface
      , rho_alb(1, ip_surf_alb_diff), rho_alb(1, ip_surf_alb_dir)       &
      , d_planck_flux_surface                                           &
!                 Spherical geometry
      , sph                                                             &
!                 Single scattering properties
      , ss_prop                                                         &
!                 Cloud geometry
      , n_cloud_top                                                     &
      , cld%n_cloud_type, cld%frac_cloud                                &
      , n_region, i_region_cloud, frac_region                           &
      , w_free, cld%w_cloud                                             &
      , cloud_overlap                                                   &
!                 Fluxes calculated
      , flux_direct, flux_total                                         &
!                 Flags for clear-sky calculations
      , l_clear, i_solver_clear                                         &
!                 Clear-sky fluxes calculated
      , flux_direct_clear, flux_total_clear                             &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct                       &
      , nd_max_order, nd_source_coeff                                   &
      , nd_cloud_type, nd_region, nd_overlap_coeff                      &
      )

  ELSE IF (i_cloud == ip_cloud_column_max) THEN
!   Clouds are treated on the assumption of maximum overlap
!   in a column model.

!   Set a dimension to allow the subcolumns of several profiles
!   to be considered at once.
    nd_profile_column=MAX(1, n_profile)
    DO l=1, n_profile
      nd_profile_column=MAX(nd_profile_column, n_column_slv(l))
    END DO

    CALL calc_flux_ipa(ierr                                             &
      , control, bound                                                  &
!                 Atmospheric properties
      , n_profile, n_layer, n_cloud_top                                 &
!                 Options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                 Algorithmic options
      , i_2stream, i_solver, i_scatter_method                           &
!                 Spectral region
      , isolir                                                          &
!                 Infra-red properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                 Conditions at TOA
      , flux_inc_down, flux_inc_direct, sec_0                           &
!                 Conditions at surface
      , d_planck_flux_surface, rho_alb                                  &
!                 Spherical geometry
      , sph                                                             &
!                 Single scattering properties
      , ss_prop                                                         &
!                 Cloud geometry
      , n_column_slv, list_column_slv                                   &
      , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!                 Fluxes calculated
      , flux_direct, flux_total                                         &
!                 Flags for clear-sky calculations
      , l_clear, i_solver_clear                                         &
!                 Clear-sky fluxes calculated
      , flux_direct_clear, flux_total_clear                             &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column            &
      , nd_profile_column, nd_source_coeff                              &
      )

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE monochromatic_radiance_tseq
END MODULE monochromatic_radiance_tseq_mod
