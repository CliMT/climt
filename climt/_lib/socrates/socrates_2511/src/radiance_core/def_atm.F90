! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for atmospheric data.
!
! Description:
!   This module contains the declaration of the structure
!   used to store atmospheric data in the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_atm

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrAtm

! Grid structure:
  INTEGER :: n_profile
!   Total number of horizontal profiles
  INTEGER :: n_layer
!   Number of atmospheric layers
  REAL (RealK), ALLOCATABLE :: lon(:)
!   Longitude coordinates of points
  REAL (RealK), ALLOCATABLE :: lat(:)
!   Latitude coordinates of points

! Specification of the viewing geometry
  INTEGER :: n_direction
!   Number of directions at which to calculate radiances
  INTEGER :: n_viewing_level
!   Number of levels where the radiance is required
  REAL (RealK), ALLOCATABLE :: direction(:, :, :)
!   Directions in which to calculate radiances
  REAL (RealK), ALLOCATABLE :: viewing_level(:)
!   List of levels where the radiance is required

! Thermodynamic Fields:
  REAL (RealK), ALLOCATABLE :: mass(:, :)
!   Column masses in each layer
  REAL (RealK), ALLOCATABLE :: density(:, :)
!   Representative density in each layer
  REAL (RealK), ALLOCATABLE :: p(:, :)
!   Pressures at the centres of layers
  REAL (RealK), ALLOCATABLE :: p_level(:, :)
!   Pressures at the edges of layers
  REAL (RealK), ALLOCATABLE :: t(:, :)
!   Temperatures at the centres of layers
  REAL (RealK), ALLOCATABLE :: t_level(:, :)
!   Temperatures at the edges of layers

! Heights (from the centre of the planet)
  REAL (RealK), ALLOCATABLE :: r_layer(:, :)
!   Height / radius at the centres of layers
  REAL (RealK), ALLOCATABLE :: r_level(:, :)
!   Height / radius at the edges of layers

! Gaseous Fields:
  REAL (RealK), ALLOCATABLE :: gas_mix_ratio(:, :, :)
!   Gaseous mixing ratios

END TYPE StrAtm


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_atm(atm, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrAtm),      INTENT(INOUT) :: atm
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(atm%lon))                                                  &
  ALLOCATE(atm%lon           ( dimen%nd_profile                              ))

IF (.NOT. ALLOCATED(atm%lat))                                                  &
  ALLOCATE(atm%lat           ( dimen%nd_profile                              ))

IF (.NOT. ALLOCATED(atm%direction))                                            &
  ALLOCATE(atm%direction     ( dimen%nd_radiance_profile,                      &
                               dimen%nd_direction,                             &
                               2                                             ))

IF (.NOT. ALLOCATED(atm%viewing_level))                                        &
  ALLOCATE(atm%viewing_level ( dimen%nd_viewing_level                        ))

IF (.NOT. ALLOCATED(atm%mass))                                                 &
  ALLOCATE(atm%mass          ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(atm%density))                                              &
  ALLOCATE(atm%density       ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(atm%p))                                                    &
  ALLOCATE(atm%p             ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(atm%t))                                                    &
  ALLOCATE(atm%t             ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(atm%p_level))                                              &
  ALLOCATE(atm%p_level       ( dimen%nd_profile,                               &
                               0 : dimen%nd_layer                            ))

IF (.NOT. ALLOCATED(atm%t_level))                                              &
  ALLOCATE(atm%t_level       ( dimen%nd_profile,                               &
                               0 : dimen%nd_layer                            ))

IF (.NOT. ALLOCATED(atm%r_layer))                                              &
  ALLOCATE(atm%r_layer       ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(atm%r_level))                                              &
  ALLOCATE(atm%r_level       ( dimen%nd_profile,                               &
                               0 : dimen%nd_layer                            ))

IF (.NOT. ALLOCATED(atm%gas_mix_ratio))                                        &
  ALLOCATE(atm%gas_mix_ratio ( dimen%nd_profile,                               &
                               dimen%nd_layer,                                 &
                               sp%dim%nd_species                             ))

END SUBROUTINE allocate_atm
!------------------------------------------------------------------------------
SUBROUTINE deallocate_atm(atm)

IMPLICIT NONE

TYPE (StrAtm), INTENT(INOUT) :: atm

IF (ALLOCATED(atm%gas_mix_ratio)) DEALLOCATE(atm%gas_mix_ratio)
IF (ALLOCATED(atm%r_level))       DEALLOCATE(atm%r_level)
IF (ALLOCATED(atm%r_layer))       DEALLOCATE(atm%r_layer)
IF (ALLOCATED(atm%t_level))       DEALLOCATE(atm%t_level)
IF (ALLOCATED(atm%p_level))       DEALLOCATE(atm%p_level)
IF (ALLOCATED(atm%t))             DEALLOCATE(atm%t)
IF (ALLOCATED(atm%p))             DEALLOCATE(atm%p)
IF (ALLOCATED(atm%density))       DEALLOCATE(atm%density)
IF (ALLOCATED(atm%mass))          DEALLOCATE(atm%mass)
IF (ALLOCATED(atm%viewing_level)) DEALLOCATE(atm%viewing_level)
IF (ALLOCATED(atm%direction))     DEALLOCATE(atm%direction)
IF (ALLOCATED(atm%lat))           DEALLOCATE(atm%lat)
IF (ALLOCATED(atm%lon))           DEALLOCATE(atm%lon)

END SUBROUTINE deallocate_atm
!------------------------------------------------------------------------------

END MODULE def_atm
