! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for spherical geometry fields
!
! Description:
!   This module contains the declaration of the structure
!   used to store fields related to treatment of spherical geometry
!   for the solar beam.
!
!------------------------------------------------------------------------------
MODULE def_spherical_geometry

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrSphCommon
  INTEGER, ALLOCATABLE :: path_base(:,:)
!   Lowest layer beam passes through  
  REAL (RealK), ALLOCATABLE :: path(:,:,:)
!   Path length through spherical shells
  REAL (RealK), ALLOCATABLE :: path_div(:,:)
!   Path scaling for calculating the flux divergence
  REAL (RealK), ALLOCATABLE :: trans_0_cloud(:,:,:)
!   Transmission through each cloudy layer 
  REAL (RealK), ALLOCATABLE :: flux_inc_direct(:,:)
!   Incident solar flux for each layer
  REAL (RealK), ALLOCATABLE :: adjust_solar_ke(:,:)
!   Adjustment of solar transmission to `include' effects
!   of minor absorbers and take out equivalent extinction
END TYPE StrSphCommon


TYPE StrSphComp
  REAL (RealK), ALLOCATABLE :: trans_0(:,:)
!   Transmission through all spherical layers to current layer
  REAL (RealK), ALLOCATABLE :: flux_direct(:,:)
!   Direct flux arriving at each layer
  REAL (RealK), ALLOCATABLE :: flux_direct_div(:,:)
!   Direct flux divergence across layer
END TYPE StrSphComp


TYPE StrSphGeo
  TYPE (StrSphCommon) :: Common
  TYPE (StrSphComp)   :: Clear
  TYPE (StrSphComp)   :: AllSky
END TYPE StrSphGeo


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_sph(sph, dimen)

USE def_dimen, ONLY: strdim

IMPLICIT NONE

TYPE (StrSphGeo),   INTENT(INOUT) :: sph
TYPE (StrDim),      INTENT(IN)    :: dimen

IF (.NOT. ALLOCATED(sph%common%path_base))                                     &
  ALLOCATE(sph%common%path_base       ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%common%path))                                          &
  ALLOCATE(sph%common%path            ( dimen%nd_profile,                      &
                                        dimen%nd_layer,                        &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%common%path_div))                                      &
  ALLOCATE(sph%common%path_div        ( dimen%nd_profile,                      &
                                        dimen%nd_layer                       ))

IF (.NOT. ALLOCATED(sph%common%trans_0_cloud))                                 &
  ALLOCATE(sph%common%trans_0_cloud   ( dimen%nd_profile,                      &
                                        dimen%id_cloud_top:dimen%nd_layer,     &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%common%flux_inc_direct))                               &
  ALLOCATE(sph%common%flux_inc_direct ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%common%adjust_solar_ke))                               &
  ALLOCATE(sph%common%adjust_solar_ke ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%clear%trans_0))                                        &
  ALLOCATE(sph%clear%trans_0          ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%clear%flux_direct))                                    &
  ALLOCATE(sph%clear%flux_direct      ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%clear%flux_direct_div))                                &
  ALLOCATE(sph%clear%flux_direct_div  ( dimen%nd_profile,                      &
                                        dimen%nd_layer                       ))

IF (.NOT. ALLOCATED(sph%allsky%trans_0))                                       &
  ALLOCATE(sph%allsky%trans_0         ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%allsky%flux_direct))                                   &
  ALLOCATE(sph%allsky%flux_direct     ( dimen%nd_profile,                      &
                                        0:dimen%nd_layer+1                   ))

IF (.NOT. ALLOCATED(sph%allsky%flux_direct_div))                               &
  ALLOCATE(sph%allsky%flux_direct_div ( dimen%nd_profile,                      &
                                        dimen%nd_layer                       ))

END SUBROUTINE allocate_sph
!------------------------------------------------------------------------------
SUBROUTINE deallocate_sph(sph)

IMPLICIT NONE

TYPE (StrSphGeo), INTENT(INOUT) :: sph

IF (ALLOCATED(sph%allsky%flux_direct_div)) &
   DEALLOCATE(sph%allsky%flux_direct_div)
IF (ALLOCATED(sph%allsky%flux_direct)) &
   DEALLOCATE(sph%allsky%flux_direct)
IF (ALLOCATED(sph%allsky%trans_0)) &
   DEALLOCATE(sph%allsky%trans_0)
IF (ALLOCATED(sph%clear%flux_direct_div)) &
   DEALLOCATE(sph%clear%flux_direct_div)
IF (ALLOCATED(sph%clear%flux_direct)) &
   DEALLOCATE(sph%clear%flux_direct)
IF (ALLOCATED(sph%clear%trans_0)) &
   DEALLOCATE(sph%clear%trans_0)
IF (ALLOCATED(sph%common%adjust_solar_ke)) &
   DEALLOCATE(sph%common%adjust_solar_ke)
IF (ALLOCATED(sph%common%flux_inc_direct)) &
   DEALLOCATE(sph%common%flux_inc_direct)
IF (ALLOCATED(sph%common%trans_0_cloud)) &
   DEALLOCATE(sph%common%trans_0_cloud)
IF (ALLOCATED(sph%common%path_div)) &
   DEALLOCATE(sph%common%path_div)
IF (ALLOCATED(sph%common%path)) &
   DEALLOCATE(sph%common%path)
IF (ALLOCATED(sph%common%path_base)) &
   DEALLOCATE(sph%common%path_base)

END SUBROUTINE deallocate_sph
!------------------------------------------------------------------------------

END MODULE def_spherical_geometry
