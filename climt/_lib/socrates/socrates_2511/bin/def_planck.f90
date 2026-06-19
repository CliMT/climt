! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for Planckian emission data.
!
! Description:
!   This module contains the declaration of the structure
!   used to store Planckian emission data in the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_planck

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrPlanck

  REAL (RealK), ALLOCATABLE :: flux(:, :)
!   Planckian flux in band at edges of layers
  REAL (RealK), ALLOCATABLE :: diff(:, :)
!   Difference in the Planckian flux across layers
  REAL (RealK), ALLOCATABLE :: diff_2(:, :)
!   Twice the 2nd difference in the Planckian flux across layers
  REAL (RealK), ALLOCATABLE :: flux_ground(:)
!   Planckian flux at the surface temperature
  REAL (RealK), ALLOCATABLE :: flux_tile(:, :)
!   Local Planckian fluxes on surface tiles
  REAL (RealK), ALLOCATABLE :: radiance(:, :)
!   Planckian radiances in the band at viewing levels

END TYPE StrPlanck


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_planck(planck, dimen)

USE def_dimen, ONLY: strdim

IMPLICIT NONE

TYPE (StrPlanck),   INTENT(INOUT) :: planck
TYPE (StrDim),      INTENT(IN)    :: dimen

IF (.NOT. ALLOCATED(planck%flux))                                              &
  ALLOCATE(planck%flux       ( dimen%nd_profile,                               &
                               0 : dimen%nd_layer                            ))

IF (.NOT. ALLOCATED(planck%diff))                                              &
  ALLOCATE(planck%diff       ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(planck%diff_2))                                            &
  ALLOCATE(planck%diff_2     ( dimen%nd_profile,                               &
                               dimen%nd_layer                                ))

IF (.NOT. ALLOCATED(planck%flux_ground))                                       &
  ALLOCATE(planck%flux_ground( dimen%nd_profile                              ))

IF (.NOT. ALLOCATED(planck%flux_tile))                                         &
  ALLOCATE(planck%flux_tile  ( dimen%nd_point_tile,                            &
                               dimen%nd_tile                                 ))

IF (.NOT. ALLOCATED(planck%radiance))                                          &
  ALLOCATE(planck%radiance   ( dimen%nd_radiance_profile,                      &
                               dimen%nd_viewing_level                        ))

END SUBROUTINE allocate_planck
!------------------------------------------------------------------------------
SUBROUTINE deallocate_planck(planck)

IMPLICIT NONE

TYPE (StrPlanck), INTENT(INOUT) :: planck

IF (ALLOCATED(planck%radiance))    DEALLOCATE(planck%radiance)
IF (ALLOCATED(planck%flux_tile))   DEALLOCATE(planck%flux_tile)
IF (ALLOCATED(planck%flux_ground)) DEALLOCATE(planck%flux_ground)
IF (ALLOCATED(planck%diff_2))      DEALLOCATE(planck%diff_2)
IF (ALLOCATED(planck%diff))        DEALLOCATE(planck%diff)
IF (ALLOCATED(planck%flux))        DEALLOCATE(planck%flux)

END SUBROUTINE deallocate_planck
!------------------------------------------------------------------------------

END MODULE def_planck
