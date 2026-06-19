! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for aerosol data.
!
! Description:
!   This module contains the declaration of the structure
!   used to store aerosol data in the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_aer

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrAer

! CLASSIC aerosols
  INTEGER, ALLOCATABLE :: mr_type_index(:)
!   Index relating aerosol_mix_ratio aerosols to aerosols in
!   the spectral information
  INTEGER, ALLOCATABLE :: mr_source(:)
!   Scheme/source of the aerosol data, to determine use in
!   changing radiative fluxes and use in diagnostics
  REAL (RealK), ALLOCATABLE :: mix_ratio(:, :, :)
!   Mixing ratios of aerosols
  REAL (RealK), ALLOCATABLE :: mean_rel_humidity(:, :)
!   Mean relative humidity applicable for aerosols (clear-sky)

! MODE aerosols
  INTEGER :: n_mode
!   Number of aerosol modes
  REAL (RealK), ALLOCATABLE :: mode_mix_ratio(:, :, :)
!   Total mass-mixing ratio for each mode
  REAL (RealK), ALLOCATABLE :: mode_absorption(:, :, :, :)
!   Waveband-averaged modal aerosol absorption
  REAL (RealK), ALLOCATABLE :: mode_scattering(:, :, :, :)
!   Waveband-averaged modal aerosol scattering
  REAL (RealK), ALLOCATABLE :: mode_asymmetry(:, :, :, :)
!   Waveband-averaged modal aerosol asymmetry

! Prescribed optical properties
  INTEGER, ALLOCATABLE :: n_opt_level_prsc(:)
!   Number of levels of prescribed optical properties of aerosols
  INTEGER, ALLOCATABLE :: n_phase_term_prsc(:)
!   Number of terms in the prescribed phase functions of aerosols
  REAL (RealK), ALLOCATABLE :: pressure_prsc(:, :, :)
!   Pressures at which optical properties of aerosols are prescribed
  REAL (RealK), ALLOCATABLE :: absorption_prsc(:, :, :, :)
!   Prescribed absorption by aerosols
  REAL (RealK), ALLOCATABLE :: scattering_prsc(:, :, :, :)
!   Prescribed scattering by aerosols
  REAL (RealK), ALLOCATABLE :: phase_fnc_prsc(:, :, :, :, :)
!   Prescribed phase functions of aerosols

END TYPE StrAer


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_aer(aer, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrAer),      INTENT(INOUT) :: aer
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(aer%mr_type_index))                                        &
  ALLOCATE(aer%mr_type_index        ( sp%dim%nd_aerosol_mr                   ))

IF (.NOT. ALLOCATED(aer%mr_source))                                            &
  ALLOCATE(aer%mr_source            ( sp%dim%nd_aerosol_mr                   ))

IF (.NOT. ALLOCATED(aer%mix_ratio))                                            &
  ALLOCATE(aer%mix_ratio            ( dimen%nd_profile,                        &
                                      dimen%nd_layer,                          &
                                      sp%dim%nd_aerosol_mr                   ))

IF (.NOT. ALLOCATED(aer%mean_rel_humidity))                                    &
  ALLOCATE(aer%mean_rel_humidity    ( dimen%nd_profile,                        &
                                      dimen%nd_layer                         ))

IF (.NOT. ALLOCATED(aer%mode_mix_ratio))                                       &
  ALLOCATE(aer%mode_mix_ratio       ( dimen%nd_profile,                        &
                                      dimen%nd_layer,                          &
                                      dimen%nd_aerosol_mode                  ))

IF (.NOT. ALLOCATED(aer%mode_absorption))                                      &
  ALLOCATE(aer%mode_absorption      ( dimen%nd_profile,                        &
                                      dimen%nd_layer,                          &
                                      dimen%nd_aerosol_mode,                   &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(aer%mode_scattering))                                      &
  ALLOCATE(aer%mode_scattering      ( dimen%nd_profile,                        &
                                      dimen%nd_layer,                          &
                                      dimen%nd_aerosol_mode,                   &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(aer%mode_asymmetry))                                       &
  ALLOCATE(aer%mode_asymmetry       ( dimen%nd_profile,                        &
                                      dimen%nd_layer,                          &
                                      dimen%nd_aerosol_mode,                   &
                                      sp%dim%nd_band                         ))

END SUBROUTINE allocate_aer
!------------------------------------------------------------------------------
SUBROUTINE allocate_aer_prsc(aer, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrAer),      INTENT(INOUT) :: aer
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(aer%n_opt_level_prsc))                                     &
  ALLOCATE(aer%n_opt_level_prsc     ( sp%dim%nd_aerosol_species              ))

IF (.NOT. ALLOCATED(aer%n_phase_term_prsc))                                    &
  ALLOCATE(aer%n_phase_term_prsc    ( sp%dim%nd_aerosol_species              ))

IF (.NOT. ALLOCATED(aer%pressure_prsc))                                        &
  ALLOCATE(aer%pressure_prsc        ( dimen%nd_profile_aerosol_prsc,           &
                                      dimen%nd_opt_level_cloud_prsc,           &
                                      sp%dim%nd_aerosol_species              ))

IF (.NOT. ALLOCATED(aer%absorption_prsc))                                      &
  ALLOCATE(aer%absorption_prsc      ( dimen%nd_profile_aerosol_prsc,           &
                                      dimen%nd_opt_level_aerosol_prsc,         &
                                      sp%dim%nd_aerosol_species,               &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(aer%scattering_prsc))                                      &
  ALLOCATE(aer%scattering_prsc      ( dimen%nd_profile_aerosol_prsc,           &
                                      dimen%nd_opt_level_aerosol_prsc,         &
                                      sp%dim%nd_aerosol_species,               &
                                      sp%dim%nd_band                         ))

IF (.NOT. ALLOCATED(aer%phase_fnc_prsc))                                       &
  ALLOCATE(aer%phase_fnc_prsc       ( dimen%nd_profile_aerosol_prsc,           &
                                      dimen%nd_opt_level_aerosol_prsc,         &
                                      dimen%nd_phf_term_aerosol_prsc,          &
                                      sp%dim%nd_aerosol_species,               &
                                      sp%dim%nd_band                         ))

END SUBROUTINE allocate_aer_prsc
!------------------------------------------------------------------------------
SUBROUTINE deallocate_aer(aer)

IMPLICIT NONE

TYPE (StrAer), INTENT(INOUT) :: aer

IF (ALLOCATED(aer%mode_asymmetry))    DEALLOCATE(aer%mode_asymmetry)
IF (ALLOCATED(aer%mode_scattering))   DEALLOCATE(aer%mode_scattering)
IF (ALLOCATED(aer%mode_absorption))   DEALLOCATE(aer%mode_absorption)
IF (ALLOCATED(aer%mode_mix_ratio))    DEALLOCATE(aer%mode_mix_ratio)
IF (ALLOCATED(aer%mean_rel_humidity)) DEALLOCATE(aer%mean_rel_humidity)
IF (ALLOCATED(aer%mix_ratio))         DEALLOCATE(aer%mix_ratio)
IF (ALLOCATED(aer%mr_source))         DEALLOCATE(aer%mr_source)
IF (ALLOCATED(aer%mr_type_index))     DEALLOCATE(aer%mr_type_index)

END SUBROUTINE deallocate_aer
!------------------------------------------------------------------------------
SUBROUTINE deallocate_aer_prsc(aer)

IMPLICIT NONE

TYPE (StrAer), INTENT(INOUT) :: aer

IF (ALLOCATED(aer%phase_fnc_prsc))   DEALLOCATE(aer%phase_fnc_prsc)
IF (ALLOCATED(aer%scattering_prsc))  DEALLOCATE(aer%scattering_prsc)
IF (ALLOCATED(aer%absorption_prsc))  DEALLOCATE(aer%absorption_prsc)
IF (ALLOCATED(aer%pressure_prsc))    DEALLOCATE(aer%pressure_prsc)
IF (ALLOCATED(aer%n_opt_level_prsc)) DEALLOCATE(aer%n_opt_level_prsc)
IF (ALLOCATED(aer%n_opt_level_prsc)) DEALLOCATE(aer%n_opt_level_prsc)

END SUBROUTINE deallocate_aer_prsc
!------------------------------------------------------------------------------

END MODULE def_aer
