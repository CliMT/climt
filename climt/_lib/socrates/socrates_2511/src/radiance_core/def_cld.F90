! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for cloud data.
!
! Description:
!   This module contains the declaration of the structure
!   used to store cloud data in the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_cld

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrCld

  INTEGER :: n_cloud_type
!   Number of types of cloud
  INTEGER :: n_condensed
!   Number of condensed components in clouds

  INTEGER, ALLOCATABLE :: type_condensed(:)
!   Types of condensed components
  INTEGER, ALLOCATABLE :: i_cloud_type(:)
!   Types of cloud to which each component contributes
  INTEGER, ALLOCATABLE :: i_condensed_param(:)
!   Parametrization schemes for components
  INTEGER, ALLOCATABLE :: condensed_n_phf(:)
!   Number of terms in the phase function

  REAL (RealK) :: dp_corr_strat
!   Decorrelation pressure scale for large scale cloud
  REAL (RealK) :: dp_corr_conv
!   Decorrelation pressure scale for convective cloud

  REAL (RealK), ALLOCATABLE :: w_cloud(:, :)
!   Total cloud area fraction in layers
  REAL (RealK), ALLOCATABLE :: frac_cloud(:, :, :)
!   In-cloud fractions of different types of cloud
  REAL (RealK), ALLOCATABLE :: condensed_mix_ratio(:, :, :)
!   Mass mixing ratios of condensate
  REAL (RealK), ALLOCATABLE :: condensed_dim_char(:, :, :)
!   Characteristic dimensions of condensed species
  REAL (RealK), ALLOCATABLE :: condensed_param_list(:, :, :)
!   Coefficients in parametrizations of condensed phases

  REAL (RealK), ALLOCATABLE :: c_cloud(:, :)
!   Convective cloud area fraction in layers
  REAL (RealK), ALLOCATABLE :: c_ratio(:, :)
!   Ratio of convective cloud condensate to mean condensate

  REAL (RealK), ALLOCATABLE :: condensed_rel_var_dens(:, :, :)
!   Relative variance of cloud density used to correct single scattering
!   parameters for cloud inhomogeneity using ip_cairns

! Prescribed optical properties
  INTEGER :: n_opt_level_drop_prsc
!   Number of levels of prescribed optical properties of droplets
  INTEGER :: n_phase_term_drop_prsc
!   Number of terms in the phase function for prescribed water droplets
  INTEGER :: n_opt_level_ice_prsc
!   Number of levels of prescribed optical properties of ice crystals
  INTEGER :: n_phase_term_ice_prsc
!   Number of terms in the prescribed phase function for ice crystals

  REAL (RealK), ALLOCATABLE :: drop_pressure_prsc(:, :)
!   Pressures at which optical properties of droplets are prescribed
  REAL (RealK), ALLOCATABLE :: drop_absorption_prsc(:, :, :)
!   Prescribed absorption by droplets
  REAL (RealK), ALLOCATABLE :: drop_scattering_prsc(:, :, :)
!   Prescribed scattering by droplets
  REAL (RealK), ALLOCATABLE :: drop_phase_fnc_prsc(:, :, :, :)
!   Prescribed phase function of droplets
  REAL (RealK), ALLOCATABLE :: ice_pressure_prsc(:, :)
!   Pressures at which optical properties of ice crystals are prescribed
  REAL (RealK), ALLOCATABLE :: ice_absorption_prsc(:, :, :)
!   Prescribed absorption by ice crystals
  REAL (RealK), ALLOCATABLE :: ice_scattering_prsc(:, :, :)
!   Prescribed scattering by ice crystals
  REAL (RealK), ALLOCATABLE :: ice_phase_fnc_prsc(:, :, :, :)
!   Prescribed phase functions of ice crystals

! Variables for the Monte Carlo Independent Column Approximation (MCICA)
  REAL (RealK), ALLOCATABLE :: c_sub(:, :, :, :)
!   Scaling factor for condensate in each sub-column
  REAL (RealK), ALLOCATABLE :: frac_cloudy(:)
!   Fraction of the profile which is cloudy
  INTEGER, ALLOCATABLE :: n_subcol_cld(:)
!   Number of sub-columns containing cloud
  INTEGER, ALLOCATABLE :: subcol_k(:, :)
!   Number of sub-columns sampled for each k_term in each band.
  INTEGER, ALLOCATABLE :: first_subcol_k(:, :)
!   The first sub-column sampled for each k-term.
  INTEGER, ALLOCATABLE :: subcol_reorder(:)
!   Order of sub-columns. (Sub-columns are rearranged so that each sub-column
!   is sampled by k-terms of equivalent "importance" in the LW as in the SW.)


END TYPE StrCld


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_cld(cld, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: strspecdata

IMPLICIT NONE

TYPE (StrCld),      INTENT(INOUT) :: cld
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(cld%type_condensed))                                       &
  ALLOCATE(cld%type_condensed       ( dimen%nd_cloud_component               ))

IF (.NOT. ALLOCATED(cld%i_cloud_type))                                         &
  ALLOCATE(cld%i_cloud_type         ( dimen%nd_cloud_component               ))

IF (.NOT. ALLOCATED(cld%i_condensed_param))                                    &
  ALLOCATE(cld%i_condensed_param    ( dimen%nd_cloud_component               ))

IF (.NOT. ALLOCATED(cld%condensed_n_phf))                                      &
  ALLOCATE(cld%condensed_n_phf      ( dimen%nd_cloud_component               ))

IF (.NOT. ALLOCATED(cld%w_cloud))                                              &
  ALLOCATE(cld%w_cloud              ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer    ))

IF (.NOT. ALLOCATED(cld%frac_cloud))                                           &
  ALLOCATE(cld%frac_cloud           ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer,     &
                                      dimen%nd_cloud_type                    ))

IF (.NOT. ALLOCATED(cld%condensed_mix_ratio))                                  &
  ALLOCATE(cld%condensed_mix_ratio  ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer,     &
                                      dimen%nd_cloud_component               ))

IF (.NOT. ALLOCATED(cld%condensed_dim_char))                                   &
  ALLOCATE(cld%condensed_dim_char   ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer,     &
                                      dimen%nd_cloud_component               ))

IF (.NOT. ALLOCATED(cld%condensed_param_list))                                 &
  ALLOCATE(cld%condensed_param_list ( sp%dim%nd_cloud_parameter,               &
                                      dimen%nd_cloud_component,                &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(cld%c_cloud))                                              &
  ALLOCATE(cld%c_cloud              ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer    ))

IF (.NOT. ALLOCATED(cld%c_ratio))                                              &
  ALLOCATE(cld%c_ratio              ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer    ))

IF (.NOT. ALLOCATED(cld%condensed_rel_var_dens))                               &
  ALLOCATE(cld%condensed_rel_var_dens ( dimen%nd_profile,                      &
                                        dimen%id_cloud_top : dimen%nd_layer,   &
                                        dimen%nd_cloud_component             ))

END SUBROUTINE allocate_cld
!------------------------------------------------------------------------------
SUBROUTINE allocate_cld_prsc(cld, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: strspecdata

IMPLICIT NONE

TYPE (StrCld),      INTENT(INOUT) :: cld
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(cld%drop_pressure_prsc))                                   &
  ALLOCATE(cld%drop_pressure_prsc   ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc          ))

IF (.NOT. ALLOCATED(cld%drop_absorption_prsc))                                 &
  ALLOCATE(cld%drop_absorption_prsc ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(cld%drop_scattering_prsc))                                 &
  ALLOCATE(cld%drop_scattering_prsc ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(cld%drop_phase_fnc_prsc))                                  &
  ALLOCATE(cld%drop_phase_fnc_prsc  ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      dimen%nd_phf_term_cloud_prsc,            &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(cld%ice_pressure_prsc))                                    &
  ALLOCATE(cld%ice_pressure_prsc    ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc          ))

IF (.NOT. ALLOCATED(cld%ice_absorption_prsc))                                  &
  ALLOCATE(cld%ice_absorption_prsc  ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(cld%ice_scattering_prsc))                                  &
  ALLOCATE(cld%ice_scattering_prsc  ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(cld%ice_phase_fnc_prsc))                                   &
  ALLOCATE(cld%ice_phase_fnc_prsc   ( dimen%nd_profile_cloud_prsc,             &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      dimen%nd_phf_term_cloud_prsc,            &
                                      sp%dim%nd_band                         ))

END SUBROUTINE allocate_cld_prsc
!------------------------------------------------------------------------------
SUBROUTINE allocate_cld_mcica(cld, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: strspecdata

IMPLICIT NONE

TYPE (StrCld),      INTENT(INOUT) :: cld
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(cld%c_sub))                                                &
  ALLOCATE(cld%c_sub                ( dimen%nd_profile,                        &
                                      dimen%id_cloud_top : dimen%nd_layer,     &
                                      dimen%nd_subcol_gen,                     &
                                      cld%n_cloud_type                       ))

IF (.NOT. ALLOCATED(cld%frac_cloudy))                                          &
  ALLOCATE(cld%frac_cloudy          ( dimen%nd_profile                       ))

IF (.NOT. ALLOCATED(cld%n_subcol_cld))                                         &
  ALLOCATE(cld%n_subcol_cld         ( dimen%nd_profile                       ))

IF (.NOT. ALLOCATED(cld%subcol_k))                                             &
  ALLOCATE(cld%subcol_k             ( sp%dim%nd_band,                          &
                                      sp%dim%nd_k_term                       ))

IF (.NOT. ALLOCATED(cld%first_subcol_k))                                       &
  ALLOCATE(cld%first_subcol_k       ( sp%dim%nd_band,                          &
                                      sp%dim%nd_k_term + 1                   ))

IF (.NOT. ALLOCATED(cld%subcol_reorder))                                       &
  ALLOCATE(cld%subcol_reorder       ( dimen%nd_subcol_req                    ))

END SUBROUTINE allocate_cld_mcica
!------------------------------------------------------------------------------
SUBROUTINE deallocate_cld(cld)

IMPLICIT NONE

TYPE (StrCld), INTENT(INOUT) :: cld

IF (ALLOCATED(cld%condensed_rel_var_dens)) &
  DEALLOCATE(cld%condensed_rel_var_dens)
IF (ALLOCATED(cld%c_ratio))              DEALLOCATE(cld%c_ratio)
IF (ALLOCATED(cld%c_cloud))              DEALLOCATE(cld%c_cloud)
IF (ALLOCATED(cld%condensed_param_list)) DEALLOCATE(cld%condensed_param_list)
IF (ALLOCATED(cld%condensed_dim_char))   DEALLOCATE(cld%condensed_dim_char)
IF (ALLOCATED(cld%condensed_mix_ratio))  DEALLOCATE(cld%condensed_mix_ratio)
IF (ALLOCATED(cld%frac_cloud))           DEALLOCATE(cld%frac_cloud)
IF (ALLOCATED(cld%w_cloud))              DEALLOCATE(cld%w_cloud)
IF (ALLOCATED(cld%condensed_n_phf))      DEALLOCATE(cld%condensed_n_phf)
IF (ALLOCATED(cld%i_condensed_param))    DEALLOCATE(cld%i_condensed_param)
IF (ALLOCATED(cld%i_cloud_type))         DEALLOCATE(cld%i_cloud_type)
IF (ALLOCATED(cld%type_condensed))       DEALLOCATE(cld%type_condensed)

END SUBROUTINE deallocate_cld
!------------------------------------------------------------------------------
SUBROUTINE deallocate_cld_prsc(cld)

IMPLICIT NONE

TYPE (StrCld), INTENT(INOUT) :: cld

IF (ALLOCATED(cld%ice_phase_fnc_prsc))   DEALLOCATE(cld%ice_phase_fnc_prsc)
IF (ALLOCATED(cld%ice_scattering_prsc))  DEALLOCATE(cld%ice_scattering_prsc)
IF (ALLOCATED(cld%ice_absorption_prsc))  DEALLOCATE(cld%ice_absorption_prsc)
IF (ALLOCATED(cld%ice_pressure_prsc))    DEALLOCATE(cld%ice_pressure_prsc)
IF (ALLOCATED(cld%drop_phase_fnc_prsc))  DEALLOCATE(cld%drop_phase_fnc_prsc)
IF (ALLOCATED(cld%drop_scattering_prsc)) DEALLOCATE(cld%drop_scattering_prsc)
IF (ALLOCATED(cld%drop_absorption_prsc)) DEALLOCATE(cld%drop_absorption_prsc)
IF (ALLOCATED(cld%drop_pressure_prsc))   DEALLOCATE(cld%drop_pressure_prsc)

END SUBROUTINE deallocate_cld_prsc
!------------------------------------------------------------------------------
SUBROUTINE deallocate_cld_mcica(cld)

IMPLICIT NONE

TYPE (StrCld), INTENT(INOUT) :: cld

IF (ALLOCATED(cld%subcol_reorder))       DEALLOCATE(cld%subcol_reorder)
IF (ALLOCATED(cld%first_subcol_k))       DEALLOCATE(cld%first_subcol_k)
IF (ALLOCATED(cld%subcol_k))             DEALLOCATE(cld%subcol_k)
IF (ALLOCATED(cld%n_subcol_cld))         DEALLOCATE(cld%n_subcol_cld)
IF (ALLOCATED(cld%frac_cloudy))          DEALLOCATE(cld%frac_cloudy)
IF (ALLOCATED(cld%c_sub))                DEALLOCATE(cld%c_sub)

END SUBROUTINE deallocate_cld_mcica
!------------------------------------------------------------------------------

END MODULE def_cld
