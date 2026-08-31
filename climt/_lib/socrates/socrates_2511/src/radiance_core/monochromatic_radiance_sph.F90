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
MODULE monochromatic_radiance_sph_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MONOCHROMATIC_RADIANCE_SPH_MOD'
CONTAINS
SUBROUTINE monochromatic_radiance_sph(ierr                              &
!                 Atmospheric Propertries
    , control, n_profile, n_layer                                       &
!                 Angular Integration
    , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc       &
    , accuracy_adaptive, euler_factor, i_sph_algorithm                  &
    , i_sph_mode, l_rescale                                             &
!                   Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , diff_planck, l_ir_source_quad, diff_planck_2                      &
!                 Conditions at TOA
    , zen_0, flux_inc_down                                              &
    , i_direct                                                          &
!                 Surface Properties
    , d_planck_flux_surface                                             &
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
!                 Optical Properties
    , ss_prop                                                           &
!                 Cloudy Properties
    , i_cloud                                                           &
!                 Cloud Geometry
    , n_cloud_top                                                       &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                   Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                   Viewing geometry
    , n_direction, direction                                            &
!                 Calculated Fluxes
    , flux_direct, flux_total                                           &
!                   Calculated radiances
    , radiance                                                          &
!                   Calculated mean radiances
    , j_radiance                                                        &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                &
    , nd_direction                                                      &
    )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_ss_prop, ONLY: str_ss_prop
  USE rad_pcf, ONLY: i_err_fatal, ip_cloud_clear, ip_solar,             &
                     ip_sph_direct, ip_sph_reduced_iter,                &
                     ip_cloud_mix_max, ip_cloud_mix_random,             &
                     ip_cloud_triple, ip_cloud_part_corr,               &
                     ip_cloud_part_corr_cnv, ip_cloud_column_max
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE calc_radiance_ipa_mod, ONLY: calc_radiance_ipa
  USE copy_clr_full_mod, ONLY: copy_clr_full
  USE copy_clr_sol_mod, ONLY: copy_clr_sol
  USE sph_solver_mod, ONLY: sph_solver

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_layer_clr                                                      &
!       Maximum number of completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_flux_profile                                                   &
!       Maximum number of profiles in arrays of fluxes
    , nd_radiance_profile                                               &
!       Maximum number of profiles in arrays of radiances
    , nd_j_profile                                                      &
!       Maximum number of profiles in arrays of mean radiances
    , nd_column                                                         &
!       Number of columns per point
    , nd_max_order                                                      &
!       Maximum order of spherical harmonics used
    , nd_sph_coeff                                                      &
!       Allocated size for spherical coefficients
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for orders of BRDFs
    , nd_viewing_level                                                  &
!       Allocated size for levels where radiances are calculated
    , nd_direction
!       Allocated size for viewing directions

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
      n_order_phase                                                     &
!       Highest order retained in the phase function
    , i_truncation                                                      &
!       Type of spherical truncation adopted
    , ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                   &
!       Orders of truncation at each azimuthal order
    , i_sph_mode                                                        &
!       Mode in which the spherical harmonic solver is being used
    , i_sph_algorithm
!       Algorithm used for spherical harmonic calculation
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale
!       Flag for rescaling of the optical properties
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
      zen_0(nd_profile)                                                 &
!       Secants (two-stream) or cosines (spherical harmonics)
!       of the solar zenith angles
    , flux_inc_down(nd_profile)
!       Incident downward flux
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_profile, 0: nd_layer)
!       Direct radiance (the first row contains the incident
!       solar radiance: the other rows are calculated)

!                 Surface properties
  REAL (RealK), INTENT(IN) ::                                           &
      d_planck_flux_surface(nd_profile)
!       Differential Planckian flux from the surface
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

!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

!                 Cloudy properties
  INTEGER, INTENT(IN) ::                                                &
      i_cloud
!       Cloud scheme used

!                 Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top
!       Topmost cloudy layer

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
!       Areas of columns


!                 Levels where radiance are calculated
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

!                 Calculated Fluxes
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct flux
    , flux_total(nd_flux_profile, 2*nd_layer+2)
!       Total flux

!                 Calculated radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      radiance(nd_radiance_profile, nd_viewing_level, nd_direction)
!       Radiances

!                 Calculated mean radiances
  REAL (RealK), INTENT(OUT) ::                                          &
      j_radiance(nd_j_profile, nd_viewing_level)
!       Mean radiances


! Local variables.
  INTEGER                                                               &
      nd_red_eigensystem                                                &
!       Size allowed for the reduced eigensystem
    , nd_sph_equation                                                   &
!       Size allowed for spherical harmonic equations
    , nd_sph_diagonal                                                   &
!       Size allowed for diagonals of the spherical harmonic
!       matrix
    , nd_sph_cf_weight                                                  &
!       Size allowed for entities to be incremented by the
!       complementary function of the linear system
    , nd_sph_u_range                                                    &
!       Size allowed for range of values of u^+|- contributing
!       on any viewing level
    , nd_profile_column
!       Size allowed for profiles taken simultaneously in a
!       decomposition into columns
  INTEGER                                                               &
      l
!       Loop variable

  REAL  (RealK), ALLOCATABLE ::                                         &
      tau_clr_f(:, :)                                                   &
!       Clear-sky optical depth for the whole column
    , tau_clr_dir_f(:, :)                                               &
!       Unscaled clear-sky optical depth
    , omega_clr_f(:, :)                                                 &
!       Clear-sky albedo of single scattering for the whole column
    , phase_fnc_clr_f(:, :, :)                                          &
!       Clear-sky phase function for the whole column
    , forward_scatter_clr_f(:, :)                                       &
!       Clear-sky forward scattering for the whole column
    , phase_fnc_solar_clr_f(:, :, :)
!       Clear-sky solar phase function in viewing directions
!       for the whole column


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'MONOCHROMATIC_RADIANCE_SPH'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Split the method of solution according to the cloud scheme.
  IF (i_cloud == ip_cloud_clear) THEN

!   Precalculate dimensions for the dynamically allocated arrays.
    nd_red_eigensystem=(nd_max_order+1)/2
    nd_sph_equation=2*nd_layer*nd_red_eigensystem
    nd_sph_diagonal=6*nd_red_eigensystem
    IF (i_sph_algorithm == ip_sph_direct) THEN
      nd_sph_cf_weight=nd_max_order+1
      nd_sph_u_range=2*nd_red_eigensystem
    ELSE IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
      nd_sph_cf_weight=nd_direction
      nd_sph_u_range=nd_sph_equation
    END IF

!   Allocate and set dynamic arrays.
    ALLOCATE(tau_clr_f(nd_profile, nd_layer))
    ALLOCATE(tau_clr_dir_f(nd_profile, nd_layer))
    ALLOCATE(omega_clr_f(nd_profile, nd_layer))
    ALLOCATE(phase_fnc_clr_f(nd_profile, nd_layer, nd_max_order))
    ALLOCATE(forward_scatter_clr_f(nd_profile, nd_layer))
    ALLOCATE(phase_fnc_solar_clr_f(nd_radiance_profile                  &
      , nd_layer, nd_direction))

    CALL copy_clr_full(n_profile, n_layer, n_cloud_top                  &
      , control, n_order_phase                                          &
      , ss_prop%tau_clr, ss_prop%tau_clr_dir                            &
      , ss_prop%omega_clr, ss_prop%phase_fnc_clr                        &
      , ss_prop%tau, ss_prop%tau_dir                                    &
      , ss_prop%omega, ss_prop%phase_fnc                                &
      , tau_clr_f, tau_clr_dir_f, omega_clr_f, phase_fnc_clr_f          &
!                   Sizes of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_max_order         &
      )
    IF ( (i_sph_algorithm == ip_sph_reduced_iter).AND.                  &
         (isolir == ip_solar) ) THEN
      CALL copy_clr_sol(n_profile, n_layer, n_cloud_top                 &
        , n_direction, l_rescale                                        &
        , ss_prop%forward_scatter_clr                                   &
        , ss_prop%phase_fnc_solar_clr                                   &
        , ss_prop%forward_scatter, ss_prop%phase_fnc_solar              &
        , forward_scatter_clr_f                                         &
        , phase_fnc_solar_clr_f                                         &
!                   Sizes of arrays
        , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_direction       &
        )
    END IF

    CALL sph_solver(                                                    &
!                   Atmospheric sizes
        n_profile, n_layer                                              &
!                   Angular integration
      , ms_min, ms_max, i_truncation, ls_local_trunc                    &
      , cg_coeff, uplm_zero, ia_sph_mm                                  &
      , accuracy_adaptive, euler_factor                                 &
      , i_sph_algorithm, i_sph_mode                                     &
!                   Spectral Region
      , isolir                                                          &
!                 Options for Equivalent Extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                   Solar Fields
      , i_direct, zen_0, uplm_sol                                       &
!                   Infra-red Properties
      , diff_planck, flux_inc_down                                      &
      , l_ir_source_quad, diff_planck_2                                 &
!                   Optical properies
      , tau_clr_f, omega_clr_f, phase_fnc_clr_f                         &
      , phase_fnc_solar_clr_f                                           &
!                   Surface Conditions
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                        &
      , f_brdf, brdf_sol, brdf_hemi                                     &
      , d_planck_flux_surface                                           &
!                   Levels for calculating radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!                   Viewing Geometry
      , n_direction, direction                                          &
!                   Radiances Calculated
      , flux_direct, flux_total, radiance, j_radiance                   &
!                   Dimensions of arrays
      , nd_profile, nd_layer                                            &
      , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
      , nd_max_order, nd_sph_coeff                                      &
      , nd_brdf_basis_fnc, nd_brdf_trunc                                &
      , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal            &
      , nd_sph_cf_weight, nd_sph_u_range                                &
      , nd_viewing_level, nd_direction                                  &
      )

!   Release temporary storage.
    DEALLOCATE(tau_clr_f)
    DEALLOCATE(tau_clr_dir_f)
    DEALLOCATE(omega_clr_f)
    DEALLOCATE(phase_fnc_clr_f)
    DEALLOCATE(forward_scatter_clr_f)
    DEALLOCATE(phase_fnc_solar_clr_f)



  ELSE IF ((i_cloud == ip_cloud_mix_max).OR.                            &
           (i_cloud == ip_cloud_mix_random).OR.                         &
           (i_cloud == ip_cloud_triple).OR.                             &
           (i_cloud == ip_cloud_part_corr).OR.                          &
           (i_cloud == ip_cloud_part_corr_cnv) ) THEN

    cmessage = '*** ERROR: Radiances cannot yet be computed using '     &
      //'coupled overlaps.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  ELSE IF (i_cloud == ip_cloud_column_max) THEN

!   Clouds are treated using the independent pixel approximation,
!   as directed by the decompositional arrays.

!   Set a dimension to allow the subcolumns of several profiles
!   to be considered at once.
    nd_profile_column=MAX(1, n_profile)
    DO l=1, n_profile
      nd_profile_column=MAX(nd_profile_column, n_column_slv(l))
    END DO

!   Precalculate dimensions for the dynamically allocated arrays.
    nd_red_eigensystem=(nd_max_order+1)/2
    nd_sph_equation=2*nd_layer*nd_red_eigensystem
    nd_sph_diagonal=6*nd_red_eigensystem
    IF (i_sph_algorithm == ip_sph_direct) THEN
      nd_sph_cf_weight=nd_max_order+1
      nd_sph_u_range=2*nd_red_eigensystem
    ELSE IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
      nd_sph_cf_weight=nd_direction
      nd_sph_u_range=nd_sph_equation
    END IF


    CALL calc_radiance_ipa(                                             &
!                 Atmospheric Properties
        n_profile, n_layer, n_cloud_top                                 &
!                   Angular Integration
      , n_order_phase, ms_min, ms_max, ls_local_trunc                   &
      , i_truncation, accuracy_adaptive, euler_factor                   &
      , i_sph_algorithm, i_sph_mode, l_rescale                          &
!                   Precalculated angular arrays
      , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                        &
!                 Options for Equivalent Extinction
      , l_scale_solar, adjust_solar_ke                                  &
!                 Spectral Region
      , isolir                                                          &
!                 Infra-red Properties
      , diff_planck                                                     &
      , l_ir_source_quad, diff_planck_2                                 &
!                 Conditions at TOA
      , flux_inc_down, zen_0                                            &
!                 Conditions at Surface
      , d_planck_flux_surface                                           &
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                        &
      , f_brdf, brdf_sol, brdf_hemi                                     &
!                 Optical Properties
      , ss_prop                                                         &
!                 Cloud Geometry
      , n_column_slv, list_column_slv                                   &
      , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!                   Levels for calculating radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!                   Viewing Geometry
      , n_direction, direction                                          &
!                   Calculated fluxes or radiances
      , flux_direct, flux_total, i_direct, radiance, j_radiance         &
!                 Dimensions of Arrays
      , nd_profile, nd_layer                                            &
      , nd_column                                                       &
      , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
      , nd_max_order, nd_sph_coeff                                      &
      , nd_brdf_basis_fnc, nd_brdf_trunc                                &
      , nd_red_eigensystem, nd_sph_equation, nd_sph_diagonal            &
      , nd_sph_cf_weight, nd_sph_u_range                                &
      , nd_viewing_level, nd_direction                                  &
      , nd_profile_column                                               &
      )

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE monochromatic_radiance_sph
END MODULE monochromatic_radiance_sph_mod
