! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for surface and TOA boundary conditions.
!
! Description:
!   This module contains the declaration of the structure
!   used to store boundary conditions for the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_bound

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrBound

! Solar Fields
  REAL (RealK), ALLOCATABLE :: zen_0(:)
!   Secants or cosines of solar zenith angles
  REAL (RealK), ALLOCATABLE :: azim_0(:)
!   Azimuthal solar angles
  REAL (RealK), ALLOCATABLE :: solar_irrad(:)
!   Solar irradiance at the top of the atmosphere
  REAL (RealK), ALLOCATABLE :: lit(:, :)
!   Lit fraction of timestep for each layer (using spherical geometry)
  REAL (RealK), ALLOCATABLE :: cos_zen(:, :)
!   Cosines of solar zenith angles for each layer (using spherical geometry)

! Surface properties
  INTEGER :: n_brdf_basis_fnc
!   Number of basis functions for BRDFs
  REAL (RealK), ALLOCATABLE :: rho_alb(:, :, :)
!   Weights for basis functions of the BRDFs
  REAL (RealK), ALLOCATABLE :: f_brdf(:, :, :, :)
!   Array of BRDF basis terms
  REAL (RealK), ALLOCATABLE :: t_ground(:)
!   Temperature of ground
  REAL (RealK), ALLOCATABLE :: flux_ground(:, :)
!   Emission of ground
  REAL (RealK), ALLOCATABLE :: orog_corr(:)
!   Correction factor for the direct solar flux
!   reaching the surface for sloping terrain.

! Arrays related to tiling of the surface
  INTEGER :: n_point_tile
!   Number of points to tile
  INTEGER :: n_tile
!   Number of tiles used
  INTEGER, ALLOCATABLE :: list_tile(:)
!   List of points with surface tiling
  REAL (RealK), ALLOCATABLE :: rho_alb_tile(:, :, :, :)
!   Weights for the basis functions of the BRDFs at the tiled points
  REAL (RealK), ALLOCATABLE :: frac_tile(:, :)
!   Fraction of tiled grid-points occupied by each tile
  REAL (RealK), ALLOCATABLE :: t_tile(:, :)
!   Local surface temperatures on individual tiles
  REAL (RealK), ALLOCATABLE :: flux_tile(:, :, :)
!   Local surface flux on individual tiles

END TYPE StrBound


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_bound(bound, dimen, sp)

USE def_dimen, ONLY: strdim
USE def_spectrum, ONLY: strspecdata

IMPLICIT NONE

TYPE (StrBound),    INTENT(INOUT) :: bound
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(bound%zen_0))                                              &
  ALLOCATE(bound%zen_0           ( dimen%nd_profile                          ))

IF (.NOT. ALLOCATED(bound%azim_0))                                             &
  ALLOCATE(bound%azim_0          ( dimen%nd_profile                          ))

IF (.NOT. ALLOCATED(bound%solar_irrad))                                        &
  ALLOCATE(bound%solar_irrad     ( dimen%nd_profile                          ))

IF (.NOT. ALLOCATED(bound%lit))                                                &
  ALLOCATE(bound%lit             ( dimen%nd_profile,                           &
                                   0:dimen%nd_layer+1                        ))

IF (.NOT. ALLOCATED(bound%cos_zen))                                            &
  ALLOCATE(bound%cos_zen         ( dimen%nd_profile,                           &
                                   0:dimen%nd_layer+1                        ))

IF (.NOT. ALLOCATED(bound%rho_alb))                                            &
  ALLOCATE(bound%rho_alb         ( dimen%nd_profile,                           &
                                   dimen%nd_brdf_basis_fnc,                    &
                                   sp%dim%nd_band                            ))

IF (.NOT. ALLOCATED(bound%f_brdf))                                             &
  ALLOCATE(bound%f_brdf          ( dimen%nd_brdf_basis_fnc,                    &
                                   0: dimen%nd_brdf_trunc/2,                   &
                                   0: dimen%nd_brdf_trunc/2,                   &
                                   0: dimen%nd_brdf_trunc                    ))

IF (.NOT. ALLOCATED(bound%t_ground))                                           &
  ALLOCATE(bound%t_ground        ( dimen%nd_profile                          ))

IF (.NOT. ALLOCATED(bound%flux_ground))                                        &
  ALLOCATE(bound%flux_ground     ( dimen%nd_profile,                           &
                                   sp%dim%nd_band                            ))

IF (.NOT. ALLOCATED(bound%orog_corr))                                          &
  ALLOCATE(bound%orog_corr       ( dimen%nd_profile                          ))

IF (.NOT. ALLOCATED(bound%list_tile))                                          &
  ALLOCATE(bound%list_tile       ( dimen%nd_point_tile                       ))

IF (.NOT. ALLOCATED(bound%rho_alb_tile))                                       &
  ALLOCATE(bound%rho_alb_tile    ( dimen%nd_point_tile,                        &
                                   dimen%nd_brdf_basis_fnc,                    &
                                   dimen%nd_tile,                              &
                                   sp%dim%nd_band                            ))

IF (.NOT. ALLOCATED(bound%frac_tile))                                          &
  ALLOCATE(bound%frac_tile       ( dimen%nd_point_tile,                        &
                                   dimen%nd_tile                             ))

IF (.NOT. ALLOCATED(bound%t_tile))                                             &
  ALLOCATE(bound%t_tile          ( dimen%nd_point_tile,                        &
                                   dimen%nd_tile                             ))

IF (.NOT. ALLOCATED(bound%flux_tile))                                          &
  ALLOCATE(bound%flux_tile       ( dimen%nd_point_tile,                        &
                                   dimen%nd_tile,                              &
                                   sp%dim%nd_band                            ))

END SUBROUTINE allocate_bound
!------------------------------------------------------------------------------
SUBROUTINE deallocate_bound(bound)

IMPLICIT NONE

TYPE (StrBound), INTENT(INOUT) :: bound

IF (ALLOCATED(bound%flux_tile))    DEALLOCATE(bound%flux_tile)
IF (ALLOCATED(bound%t_tile))       DEALLOCATE(bound%t_tile)
IF (ALLOCATED(bound%frac_tile))    DEALLOCATE(bound%frac_tile)
IF (ALLOCATED(bound%rho_alb_tile)) DEALLOCATE(bound%rho_alb_tile)
IF (ALLOCATED(bound%list_tile))    DEALLOCATE(bound%list_tile)
IF (ALLOCATED(bound%orog_corr))    DEALLOCATE(bound%orog_corr)
IF (ALLOCATED(bound%flux_ground))  DEALLOCATE(bound%flux_ground)
IF (ALLOCATED(bound%t_ground))     DEALLOCATE(bound%t_ground)
IF (ALLOCATED(bound%f_brdf))       DEALLOCATE(bound%f_brdf)
IF (ALLOCATED(bound%rho_alb))      DEALLOCATE(bound%rho_alb)
IF (ALLOCATED(bound%cos_zen))      DEALLOCATE(bound%cos_zen)
IF (ALLOCATED(bound%lit))          DEALLOCATE(bound%lit)
IF (ALLOCATED(bound%solar_irrad))  DEALLOCATE(bound%solar_irrad)
IF (ALLOCATED(bound%azim_0))       DEALLOCATE(bound%azim_0)
IF (ALLOCATED(bound%zen_0))        DEALLOCATE(bound%zen_0)

END SUBROUTINE deallocate_bound
!------------------------------------------------------------------------------

END MODULE def_bound
