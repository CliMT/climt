! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the solar flux and source terms.
!
!- ---------------------------------------------------------------------
MODULE solar_source_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLAR_SOURCE_MOD'
CONTAINS
SUBROUTINE solar_source(control, bound, n_profile, n_layer              &
     , flux_inc_direct                                                  &
     , trans_0_dir, trans_0, source_coeff                               &
     , l_scale_solar, adjust_solar_ke                                   &
     , flux_direct                                                      &
     , s_down, s_up                                                     &
     , nd_profile, nd_layer, nd_source_coeff                            &
     )


  USE realtype_rd, ONLY: RealK
  USE def_control, ONLY: StrCtrl
  USE def_bound,   ONLY: StrBound
  USE rad_pcf,     ONLY: ip_direct_csr_scaling,ip_direct_noscaling,     &
                         ip_scf_solar_down,ip_scf_solar_up
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Boundary conditions:
  TYPE(StrBound),     INTENT(IN)    :: bound

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_source_coeff
!       Size allocated for coefficients in the source terms


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers

  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Scaling applied to solar beam

  REAL (RealK), INTENT(IN) ::                                           &
      flux_inc_direct(nd_profile)                                       &
!       Incident solar flux
    , trans_0(nd_profile, nd_layer)                                     &
!       Direct transmission coefficient
    , trans_0_dir(nd_profile, nd_layer)                                 &
!       Direct transmission using direct tau
    , source_coeff(nd_profile, nd_layer, nd_source_coeff)               &
!       Reflection coefficient
    , adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment to solar flux


  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , s_down(nd_profile, nd_layer)                                      &
!       Downward source function
    , s_up(nd_profile, nd_layer)
!       Upward source function


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) ::                                                       &
      flux_direct_dir(nd_profile, 0: nd_layer)
!       Direct flux without scaling


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLAR_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO l=1, n_profile
    flux_direct(l, 0)=flux_inc_direct(l)
  END DO
  IF (control%i_direct_tau == ip_direct_noscaling .OR.                  &
      control%i_direct_tau == ip_direct_csr_scaling) THEN
    DO l=1, n_profile
      flux_direct_dir(l, 0)=flux_inc_direct(l)
    END DO
  END IF

  IF (l_scale_solar) THEN
!   The solar flux may be multiplied by a scaling factor if an
!   equivalent extinction is used.
    DO i=1, n_layer
      DO l=1, n_profile
!       The delta-scaled direct flux is used for the source terms
!       and the total downwards flux.
        flux_direct(l, i)                                               &
          =flux_direct(l, i-1)*trans_0(l, i)                            &
          *adjust_solar_ke(l, i)
        s_up(l, i)=source_coeff(l, i, ip_scf_solar_up)                  &
          *flux_direct(l, i-1)
        s_down(l, i)=(source_coeff(l, i, ip_scf_solar_down)             &
          -trans_0(l, i))*flux_direct(l, i-1)                           &
          +flux_direct(l, i)
      END DO
    END DO
    IF (control%i_direct_tau == ip_direct_noscaling .OR.                &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      DO i=1, n_layer
        DO l=1, n_profile
!         The unscaled direct flux is calculated separately.
          flux_direct_dir(l, i)                                         &
            =flux_direct_dir(l, i-1)*trans_0_dir(l, i)                  &
            *adjust_solar_ke(l, i)
        END DO
      END DO
    END IF

  ELSE

    DO i=1, n_layer
      DO l=1, n_profile
        flux_direct(l, i)                                               &
          =flux_direct(l, i-1)*trans_0(l, i)
        s_up(l, i)=source_coeff(l, i, ip_scf_solar_up)                  &
          *flux_direct(l, i-1)
        s_down(l, i)=source_coeff(l, i, ip_scf_solar_down)              &
          *flux_direct(l, i-1)
      END DO
    END DO
    IF (control%i_direct_tau == ip_direct_noscaling .OR.                &
        control%i_direct_tau == ip_direct_csr_scaling) THEN
      DO i=1, n_layer
        DO l=1, n_profile
          flux_direct_dir(l, i)                                         &
           =flux_direct_dir(l, i-1)*trans_0_dir(l, i)
        END DO
      END DO
    END IF
  END IF


  IF (control%i_direct_tau == ip_direct_noscaling .OR.                  &
      control%i_direct_tau == ip_direct_csr_scaling) THEN
!   From this point, use the unscaled direct flux as the direct component.
    DO i= 0, n_layer
      DO l=1, n_profile
        flux_direct(l, i)=flux_direct_dir(l, i)
      END DO
    END DO
  END IF


! Correct the direct flux at the ground for sloping terrain
  IF (control%l_orog) THEN
    IF (control%l_orog_fix) THEN
!     The total flux at the surface must also add the correction made
!     to the direct flux.
      s_down(1:n_profile, n_layer) =                                    &
        s_down(1:n_profile, n_layer) +                                  &
        flux_direct(1:n_profile, n_layer) *                             &
        (bound%orog_corr(1:n_profile) - 1.0_RealK)

      flux_direct(1:n_profile, n_layer) =                               &
        flux_direct(1:n_profile, n_layer) *                             &
        bound%orog_corr(1:n_profile)
    ELSE
      flux_direct(1:n_profile, n_layer) =                               &
        flux_direct(1:n_profile, n_layer) *                             &
        bound%orog_corr(1:n_profile)

!     The total flux at the surface must also add the correction made
!     to the direct flux.
      s_down(1:n_profile, n_layer) =                                    &
        s_down(1:n_profile, n_layer) +                                  &
        flux_direct(1:n_profile, n_layer) *                             &
        (bound%orog_corr(1:n_profile) - 1.0_RealK)
    END IF
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solar_source
END MODULE solar_source_mod
