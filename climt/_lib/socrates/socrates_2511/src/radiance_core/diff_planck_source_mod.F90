! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate differences in source functions.
!
! Method:
!   The Planck fluxes for the band at the boundaries of layers
!   are found and differences across layers are determined. If
!   the Planckian is being taken to vary quadratically across
!   the layer second differences are found.
!
!------------------------------------------------------------------------------
MODULE diff_planck_source_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DIFF_PLANCK_SOURCE_MOD'
CONTAINS

SUBROUTINE diff_planck_source(control, dimen, spectrum, atm, bound, &
  i_band, i_rad_layer, frac_rad_layer, planck)

  USE def_control,  ONLY: StrCtrl
  USE def_dimen,    ONLY: StrDim
  USE def_spectrum, ONLY: StrSpecData
  USE def_atm,      ONLY: StrAtm
  USE def_bound,    ONLY: StrBound
  USE def_planck,   ONLY: StrPlanck
  USE realtype_rd,  ONLY: RealK
  USE rad_pcf,      ONLY: ip_spherical_harmonic, ip_surf_alb_diff
  USE rad_ccf,      ONLY: pi
  USE yomhook,      ONLY: lhook, dr_hook
  USE parkind1,     ONLY: jprb, jpim

  USE planck_flux_band_mod, ONLY: planck_flux_band

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),     INTENT(IN)  :: control

! Dimensions:
  TYPE(StrDim),      INTENT(IN)  :: dimen

! Spectral data:
  TYPE(StrSpecData), INTENT(IN)  :: spectrum

! Atmospheric properties:
  TYPE(StrAtm),      INTENT(IN)  :: atm

! Boundary conditions:
  TYPE(StrBound),    INTENT(IN)  :: bound

  INTEGER, INTENT(IN) :: &
    i_band, &
!     Spectral band
    i_rad_layer(dimen%nd_viewing_level)
!     Layers in which to intercept radiances
  REAL (RealK), INTENT(IN) :: &
    frac_rad_layer(dimen%nd_viewing_level)
!     Fractions below the tops of the layers

! Planckian emission fields
  TYPE(StrPlanck), INTENT(INOUT) :: planck


! Local variables.
  INTEGER :: i, k, l
!     Loop variables
  REAL (RealK), PARAMETER :: frac_tile_tol = TINY(1.0_RealK)
!     Tolerance for tile fraction when setting Planckian to zero

  LOGICAL :: l_flux_tile = .FALSE.
!     Local flag to use tile emission rather than temperature

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='DIFF_PLANCK_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (control%i_angular_integration == ip_spherical_harmonic) THEN
!   Calculate the Planckian radiance on viewing levels.
    DO i=1, atm%n_viewing_level
      DO l=1, atm%n_profile
        planck%radiance(l, i) = planck_flux_band( spectrum, i_band, &
          atm%t_level(l, i_rad_layer(i)-1) &
          +(atm%t_level(l, i_rad_layer(i))-atm%t_level(l, i_rad_layer(i)-1)) &
          *frac_rad_layer(i) )
        planck%radiance(l, i)=planck%radiance(l, i)/pi
      END DO
    END DO
  END IF

! Calculate the change in the Planckian flux across each layer.
  DO l=1, atm%n_profile
    planck%flux(l, 0) = planck_flux_band(spectrum, i_band, atm%t_level(l, 0))
  END DO
  DO i=1, atm%n_layer
    DO l=1, atm%n_profile
      planck%flux(l,i) = planck_flux_band(spectrum, i_band, atm%t_level(l,i))
      planck%diff(l,i) = planck%flux(l, i) - planck%flux(l, i-1)
    END DO
  END DO

  IF (control%l_ir_source_quad) THEN
    ! Calculate the second difference if required.
    DO i=1, atm%n_layer
      DO l=1, atm%n_profile
        planck%diff_2(l, i) = planck_flux_band(spectrum, i_band, atm%t(l, i))
        planck%diff_2(l, i) = 2.0e+00_RealK*(planck%flux(l, i) &
          + planck%flux(l, i-1) - 2.0e+00_RealK * planck%diff_2(l, i))
      END DO
    END DO
  END IF

! Planckian flux at the surface.
! Calculated from surface emission or skin temperature.
  IF (control%l_flux_ground) THEN
    DO l=1, atm%n_profile
      planck%flux_ground(l) &
        = bound%flux_ground(l, i_band) &
        /(1.0_RealK - bound%rho_alb(l, ip_surf_alb_diff, i_band))
    END DO
  ELSE
    DO l=1, atm%n_profile
      planck%flux_ground(l) &
        = planck_flux_band(spectrum, i_band, bound%t_ground(l))
    END DO
  END IF

! Local Planckian fluxes will be required on tiled surfaces.
! Furthermore, the overall Planckian will be calculated as a
! weighted sum of the individual components: this allows for
! variations in the Planckian between spectral bands more
! satisfactorily than the use of an equivalent temperature can.
  IF (control%l_tile) THEN

    ! Fluxes for tiles calculated from surface emission or skin temperature.
    DO k=1, bound%n_tile
      IF (allocated(control%l_flux_tile)) l_flux_tile = control%l_flux_tile(k)
      DO l=1, bound%n_point_tile
        IF (l_flux_tile .AND. bound%frac_tile(l, k) > frac_tile_tol) THEN
          planck%flux_tile(l, k) &
            = bound%flux_tile(l, k, i_band) &
            /(1.0_RealK - bound%rho_alb_tile(l, ip_surf_alb_diff, k, i_band))
        ELSE IF (bound%frac_tile(l, k) > frac_tile_tol) THEN
          planck%flux_tile(l, k) &
            = planck_flux_band(spectrum, i_band, bound%t_tile(l, k))
        ELSE
          planck%flux_tile(l, k) = 0.0_RealK
        END IF
      END DO
    END DO

    IF (control%l_tile_emissivity) THEN
      DO l=1, bound%n_point_tile
        planck%flux_ground(bound%list_tile(l)) &
          =bound%frac_tile(l, 1)*planck%flux_tile(l, 1) &
          *(1.0_RealK - bound%rho_alb_tile(l, ip_surf_alb_diff, 1, i_band))
      END DO
      DO k=2, bound%n_tile
        DO l=1, bound%n_point_tile
          planck%flux_ground(bound%list_tile(l)) &
            =planck%flux_ground(bound%list_tile(l)) &
            +bound%frac_tile(l, k)*planck%flux_tile(l, k) &
            *(1.0_RealK - bound%rho_alb_tile(l, ip_surf_alb_diff, k, i_band))
        END DO
      END DO
      DO l=1, bound%n_point_tile
        planck%flux_ground(bound%list_tile(l)) &
          =planck%flux_ground(bound%list_tile(l)) / ( 1.0_RealK &
          -bound%rho_alb(bound%list_tile(l), ip_surf_alb_diff, i_band) )
      END DO
    ELSE
      DO l=1, bound%n_point_tile
        planck%flux_ground(bound%list_tile(l)) &
          =bound%frac_tile(l, 1)*planck%flux_tile(l, 1)
      END DO
      DO k=2, bound%n_tile
        DO l=1, bound%n_point_tile
          planck%flux_ground(bound%list_tile(l)) &
            =planck%flux_ground(bound%list_tile(l)) &
            +bound%frac_tile(l, k)*planck%flux_tile(l, k)
        END DO
      END DO
    END IF

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diff_planck_source
END MODULE diff_planck_source_mod
