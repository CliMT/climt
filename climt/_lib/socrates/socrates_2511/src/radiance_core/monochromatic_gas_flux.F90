! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate fluxes including only gaseous absorption.
!
! Method:
!   Transmission coefficients for each layer are calculated
!   from the gaseous absorption alone. fluxes are propagated
!   upward or downward through the column using these
!   coefficients and source terms.
!
!- ---------------------------------------------------------------------
MODULE monochromatic_gas_flux_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MONOCHROMATIC_GAS_FLUX_MOD'
CONTAINS
SUBROUTINE monochromatic_gas_flux(n_profile, n_layer                    &
    , tau_gas                                                           &
    , isolir, sec_0, flux_inc_direct, flux_inc_down                     &
    , diff_planck, d_planck_flux_surface                                &
    , diffuse_albedo, direct_albedo                                     &
    , diffusivity_factor                                                &
    , flux_direct, flux_diffuse                                         &
    , nd_profile, nd_layer                                              &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_infra_red, ip_solar
  USE vectlib_mod, ONLY : exp_v
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer
!       Maximum number of layers


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , isolir
!       Spectral region
  REAL (RealK), INTENT(IN) ::                                           &
      tau_gas(nd_profile, nd_layer)                                     &
!       Gaseous optical depths
    , sec_0(nd_profile)                                                 &
!       Secant of zenith angle
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)                                         &
!       Incident diffuse flux
    , d_planck_flux_surface(nd_profile)                                 &
!       Difference in Planckian fluxes between the surface
!       and the overlying air
    , diff_planck(nd_profile, nd_layer)                                 &
!       Difference in Planckian function
    , diffuse_albedo(nd_profile)                                        &
!       Diffuse surface albedo
    , direct_albedo(nd_profile)                                         &
!       Direct surface albedo
    , diffusivity_factor
!       Diffusivity factor
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_diffuse(nd_profile, 2*nd_layer+2)
!       Diffuse flux

! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      trans(n_profile, n_layer)                                         &
    , trans_in(n_profile, n_layer)                                      &
!       Transmissivities
    , source_up(nd_profile, nd_layer)                                   &
!       Upward source function
    , source_down(nd_profile, nd_layer)
!       Downward source function

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r                                                             &
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision
    , sq_eps_r
!       The square root of the above

  REAL (RealK) :: temp(n_profile)
  REAL (RealK) :: temp_in(n_profile)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MONOCHROMATIC_GAS_FLUX'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(tau_gas(1, 1))
  sq_eps_r=SQRT(eps_r)

  DO i=1, n_layer
     DO l=1, n_profile
        trans_in(l, i)= -diffusivity_factor*tau_gas(l, i)
     END DO
  END DO
  CALL exp_v(n_layer*n_profile,trans_in,trans)

  IF (isolir == ip_solar) THEN
    DO i=1, n_layer
      DO l=1, n_profile
        source_up(l, i)=0.0e+00_RealK
        source_down(l, i)=0.0e+00_RealK
      END DO
    END DO
  ELSE IF (isolir == ip_infra_red) THEN
    DO i=1, n_layer
      DO l=1, n_profile
        source_up(l, i)=(1.0e+00_RealK-trans(l, i)+sq_eps_r)            &
          *diff_planck(l, i)                                            &
          /(diffusivity_factor*tau_gas(l, i)+sq_eps_r)
        source_down(l, i)=-source_up(l, i)
      END DO
    END DO
  END IF

! The direct flux.
  IF (isolir == ip_solar) THEN
    DO l=1, n_profile
      flux_direct(l, 0)=flux_inc_direct(l)
    END DO
    DO i=1, n_layer
      DO l=1, n_profile
        temp_in(l)=-tau_gas(l, i)*sec_0(l)
      END DO
      CALL exp_v(n_profile, temp_in, temp)
      DO l=1, n_profile
        flux_direct(l, i)=flux_direct(l, i-1)*temp(l)
      END DO
    END DO
  END IF

! Downward fluxes.
  DO l=1, n_profile
    flux_diffuse(l, 2)=flux_inc_down(l)
  END DO
  DO i=1, n_layer
    DO l=1, n_profile
      flux_diffuse(l, 2*i+2)=trans(l, i)*flux_diffuse(l, 2*i)           &
        +source_down(l, i)
    END DO
  END DO

! Upward fluxes.
  IF (isolir == ip_solar) THEN
    DO l=1, n_profile
      flux_diffuse(l, 2*n_layer+1)=                                     &
        +diffuse_albedo(l)*flux_diffuse(l, 2*n_layer+2)                 &
        +direct_albedo(l)*flux_direct(l, n_layer)
    END DO
  ELSE
    DO l=1, n_profile
      flux_diffuse(l, 2*n_layer+1)=d_planck_flux_surface(l)             &
        +diffuse_albedo(l)*flux_diffuse(l, 2*n_layer+2)
    END DO
  END IF
  DO i=n_layer, 1, -1
    DO l=1, n_profile
      flux_diffuse(l, 2*i-1)=trans(l, i)*flux_diffuse(l, 2*i+1)         &
        +source_up(l, i)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE monochromatic_gas_flux
END MODULE monochromatic_gas_flux_mod
