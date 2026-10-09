! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE spherical_path_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SPHERICAL_PATH_MOD'
CONTAINS

! Calculation of path of solar beam using spherical geometry
!
!------------------------------------------------------------------------------
SUBROUTINE spherical_path(dimen, atm, bound, sph)

USE realtype_rd,            ONLY: RealK
USE def_dimen,              ONLY: StrDim
USE def_atm,                ONLY: StrAtm
USE def_bound,              ONLY: StrBound
USE def_spherical_geometry, ONLY: StrSphGeo
USE parkind1,               ONLY: jprb, jpim
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE


! Dimensions:
TYPE(StrDim),    INTENT(IN)    :: dimen

! Atmospheric properties:
TYPE(StrAtm),    INTENT(IN)    :: atm

! Boundary conditions:
TYPE(StrBound),  INTENT(IN)    :: bound

! Spherical geometry fields
TYPE(StrSphGeo), INTENT(INOUT) :: sph


! Local variables
INTEGER :: l, i, ii
!   Loop integers
REAL :: impact_sq(dimen%nd_profile, 0:dimen%nd_layer+1)
!   Impact parameter squared: height at which beam is tangential to the surface
REAL :: r_layer_sq(dimen%nd_profile, dimen%nd_layer)
!   Radius on layers squared
REAL :: r_level_sq(dimen%nd_profile, 0:dimen%nd_layer)
!   Radius on levels squared
REAL :: recip_layer_depth(dimen%nd_profile, dimen%nd_layer)
!   Reciprocal of the layer depth
REAL :: path_upper, path_lower, path_weight
!   Limiting paths and weight when the impact parameter is in the layer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPHERICAL_PATH'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! First deal with the path for beams leaving the top-of-atmosphere (TOA).
! Only consider radiation from below.
ii=0
DO l = 1, atm%n_profile
  r_level_sq(l,0) = atm%r_level(l,0)*atm%r_level(l,0)
  IF (bound%cos_zen(l,ii) < 0.0_RealK) THEN
    impact_sq(l,ii)=r_level_sq(l,0) &
      *(1.0_RealK-bound%cos_zen(l,ii)*bound%cos_zen(l,ii))
  END IF
END DO
DO i=1, atm%n_layer
  DO l = 1, atm%n_profile
    recip_layer_depth(l,i) = 1.0_RealK/(atm%r_level(l,i-1)-atm%r_level(l,i))
    r_level_sq(l,i) = atm%r_level(l,i)*atm%r_level(l,i)
    IF (bound%cos_zen(l,ii) < 0.0_RealK) THEN
      IF (impact_sq(l,ii) >= r_level_sq(l,i-1)) THEN
        ! Beam does not pass through this layer
        sph%common%path(l,i,ii) = 0.0_RealK
      ELSE IF (impact_sq(l,ii) >= r_level_sq(l,i)) THEN
        ! Beam enters top of layer, passes down to impact parameter
        ! then back up to the top of the layer
        sph%common%path(l,i,ii) = &
          2.0_RealK*SQRT(r_level_sq(l,i-1) - impact_sq(l,ii)) &
          * recip_layer_depth(l,i)
      ELSE
        ! Beam passes through layer twice, once down, once up
        sph%common%path(l,i,ii) = &
          2.0_RealK*( SQRT(r_level_sq(l,i-1) - impact_sq(l,ii))   &
                    - SQRT(r_level_sq(l,i)   - impact_sq(l,ii)) ) &
          * recip_layer_depth(l,i)
      END IF
    ELSE
      ! If beam arriving from above it does not pass through any layers
      sph%common%path(l,i,ii) = 0.0_RealK
    END IF
  END DO
END DO

! Now consider the path to each layer.
DO ii=1, atm%n_layer
  DO l = 1, atm%n_profile
    IF (bound%lit(l,ii) > 0.0_RealK) THEN
      r_layer_sq(l,ii)=atm%r_layer(l,ii)*atm%r_layer(l,ii)
      impact_sq(l,ii)=r_layer_sq(l,ii) &
        *(1.0_RealK-bound%cos_zen(l,ii)*bound%cos_zen(l,ii))
    END IF
  END DO
  ! The spherical path taken before arriving at the layer
  DO i=1, atm%n_layer
    DO l = 1, atm%n_profile
      IF (bound%lit(l,ii) > 0.0_RealK) THEN
        IF ( bound%cos_zen(l,ii) >= 0.0_RealK .OR. &
             impact_sq(l,ii) >= r_level_sq(l,ii) ) THEN
          ! Where the solar beam is above the horizontal or where the impact
          ! parameter is within the layer we just consider the path the beam
          ! takes before arriving at the top of the current layer.
          IF (i < ii) THEN
            ! The path length is scaled by the height of the layer so that it
            ! can be used later to scale the optical depth.
            sph%common%path(l,i,ii) = &
              ( SQRT(r_level_sq(l,i-1) - impact_sq(l,ii))   &
              - SQRT(r_level_sq(l,i)   - impact_sq(l,ii)) ) &
              * recip_layer_depth(l,i)
          ELSE
            sph%common%path(l,i,ii) = 0.0_RealK
          END IF
        ELSE
          ! Otherwise the beam will be arriving at the bottom of the current
          ! layer and we need to consider the path from the top of the
          ! atmosphere down to the impact parameter and then upwards from that
          ! point to the bottom of the current layer.
          IF (i <= ii) THEN
            ! Beam passes through upper layers once on the way down
            sph%common%path(l,i,ii) = &
              ( SQRT(r_level_sq(l,i-1) - impact_sq(l,ii))   &
              - SQRT(r_level_sq(l,i)   - impact_sq(l,ii)) ) &
              * recip_layer_depth(l,i)            
          ELSE
            IF (impact_sq(l,ii) >= r_level_sq(l,i-1)) THEN
              ! Beam does not pass through this layer
              sph%common%path(l,i,ii) = 0.0_RealK
            ELSE IF (impact_sq(l,ii) >= r_level_sq(l,i)) THEN
              ! Beam enters top of layer, passes down to impact parameter
              ! then back up to the top of the layer
              sph%common%path(l,i,ii) = &
                2.0_RealK*SQRT(r_level_sq(l,i-1) - impact_sq(l,ii)) &
                * recip_layer_depth(l,i)
            ELSE
              ! Beam passes through layer twice, once down, once up
              sph%common%path(l,i,ii) = &
                2.0_RealK*( SQRT(r_level_sq(l,i-1) - impact_sq(l,ii))   &
                          - SQRT(r_level_sq(l,i)   - impact_sq(l,ii)) ) &
                * recip_layer_depth(l,i)
            END IF
          END IF
        END IF
      ELSE
        sph%common%path(l,i,ii) = 0.0_RealK
      END IF ! If lit
    END DO
  END DO
  ! Calculate the path within the layer, including the scaling for the column
  ! element over which the flux divergence is calculated (path_div), and the
  ! path before arriving at the head of this column element (path).
  DO l = 1, atm%n_profile
    IF (bound%lit(l,ii) > 0.0_RealK) THEN
      IF ( bound%cos_zen(l,ii) >= 0.0_RealK ) THEN
        ! Path to the middle of the layer from the top, scaled by the
        ! depth down to the middle of the layer.
        sph%common%path_div(l,ii) = &
          ( SQRT(r_level_sq(l,ii-1) - impact_sq(l,ii)) &
          - SQRT(r_layer_sq(l,ii)   - impact_sq(l,ii)) ) &
          / ( atm%r_level(l,ii-1) - atm%r_layer(l,ii) )
      ELSE IF ( impact_sq(l,ii) >= r_level_sq(l,ii) ) THEN
        ! Path down to impact parameter and back up to the middle of the layer.
        sph%common%path(l,ii,ii) = &
          SQRT(r_level_sq(l,ii-1) - impact_sq(l,ii)) + &
          SQRT(r_layer_sq(l,ii)   - impact_sq(l,ii))
        ! Here we need to partition this path between the part contained in
        ! the column element and the path before arriving at the element.
        ! This should transition smoothly between the known limits for the
        ! path when the impact parameter coincides with the middle of the
        ! layer and when it coincides with the lower boundary of the layer.
        path_upper = SQRT(r_level_sq(l,ii-1) - r_layer_sq(l,ii))
        path_lower = SQRT(r_layer_sq(l,ii)   - r_level_sq(l,ii))
        path_weight = ( atm%r_layer(l,ii) - SQRT(impact_sq(l,ii)) ) &
                    / ( atm%r_layer(l,ii) - atm%r_level(l,ii) )
        ! Subtract the path within the column element and scale the
        ! remaining path by the depth of the layer.
        sph%common%path(l,ii,ii) = ( sph%common%path(l,ii,ii) &
          - path_upper*(1.0_RealK - path_weight) - path_lower*path_weight ) &
          * recip_layer_depth(l,ii)
        ! The scaling for the column element weights the contributions from
        ! the upper and lower path lengths by the corresponding depths of the
        ! upper and lower portions of the layer.
        sph%common%path_div(l,ii) &
          = path_upper*(1.0_RealK - path_weight) &
          / ( atm%r_level(l,ii-1) - atm%r_layer(l,ii) ) &
          + path_lower*path_weight &
          / ( atm%r_layer(l,ii) - atm%r_level(l,ii) )
      ELSE
        ! Path to the middle of the layer from the bottom, scaled by the
        ! depth up to the middle of the layer.
        sph%common%path_div(l,ii) = &
          ( SQRT(r_layer_sq(l,ii) - impact_sq(l,ii)) &
          - SQRT(r_level_sq(l,ii) - impact_sq(l,ii)) ) &
          / ( atm%r_layer(l,ii) - atm%r_level(l,ii) )
      END IF
    ELSE
      sph%common%path_div(l,ii) = 1.0_RealK
    END IF
  END DO
END DO

! Finally deal with the path to the surface which is only lit when the
! beam arrives from above
ii=atm%n_layer+1
DO l = 1, atm%n_profile
  IF (bound%lit(l,ii) > 0.0_RealK) THEN
    impact_sq(l,ii) = r_level_sq(l,atm%n_layer) &
      * (1.0_RealK-bound%cos_zen(l,ii)*bound%cos_zen(l,ii))
  END IF
END DO
DO i=1, atm%n_layer
  DO l = 1, atm%n_profile
    IF (bound%lit(l,ii) > 0.0_RealK) THEN
      sph%common%path(l,i,ii) = &
        ( SQRT(r_level_sq(l,i-1) - impact_sq(l,ii))   &
        - SQRT(r_level_sq(l,i)   - impact_sq(l,ii)) ) &
        * recip_layer_depth(l,i)
    ELSE
      sph%common%path(l,i,ii) = 0.0_RealK
    END IF
  END DO
END DO

! Set the lowest layer beam passes through
sph%common%path_base(:, :) = 1
DO ii = 0, atm%n_layer+1
  DO i = 1, atm%n_layer
    DO l = 1, atm%n_profile
      IF (sph%common%path(l, i, ii) > 0.0_RealK) THEN
        sph%common%path_base(l, ii) = i
      END IF
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE spherical_path
END MODULE spherical_path_mod
