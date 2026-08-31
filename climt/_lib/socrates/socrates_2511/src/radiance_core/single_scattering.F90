! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to find the optical depth and single scattering albedo.
!
! Method:
!   Depending on the treatment of scattering, the optical and
!   and single scattering albedo are determined from the
!   extinctions supplied.
!
!- ---------------------------------------------------------------------
MODULE single_scattering_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SINGLE_SCATTERING_MOD'
CONTAINS
SUBROUTINE single_scattering(i_scatter_method_band                      &
    , n_profile, i_first_layer, i_last_layer                            &
    , d_mass                                                            &
    , k_grey_tot, k_ext_scat, k_gas_abs                                 &
    , tau, omega                                                        &
    , nd_profile, nd_layer, id_lt, id_lb                                &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_no_scatter_abs, ip_no_scatter_ext,              &
                     ip_scatter_approx, ip_scatter_full
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , id_lt                                                             &
!       Topmost declared layer for optical properties
    , id_lb
!       Bottom declared layer for optical properties


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method_band
!       Treatment of scattering in this band

!                 Atmospheric properties
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_first_layer                                                     &
!       First layer to consider
    , i_last_layer
!       Last layer to consider
  REAL (RealK), INTENT(IN) ::                                           &
      d_mass(nd_profile, nd_layer)
!       Mass thickness of each layer

!                 Optical properties
  REAL (RealK), INTENT(IN) ::                                           &
      k_grey_tot(nd_profile, id_lt: id_lb)                              &
!       Absorptive extinction
    , k_ext_scat(nd_profile, id_lt: id_lb)                              &
!       Scattering extinction
    , k_gas_abs(nd_profile, nd_layer)
!       Gaseous extinction

!                 Single scattering properties
  REAL (RealK), INTENT(OUT) ::                                          &
      tau(nd_profile, id_lt: id_lb)                                     &
!       Optical depth
    , omega(nd_profile, id_lt: id_lb)
!       Single scattering albedo



! Local variables.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , i
!       Loop variable

  REAL (RealK) ::                                                       &
      k_total
!       Total extinction including gaseous contributions

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SINGLE_SCATTERING'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(tau(1, id_lt))

! The machine tolerance is added to the denominator in the
! expression for omega to prevent division by zero: this is
! significant only if the total extinction is small, and thus
! will not sensibly affect any physical results.

  IF ((i_scatter_method_band == ip_scatter_full) .OR.                   &
      (i_scatter_method_band == ip_scatter_approx)) THEN

    DO i=i_first_layer, i_last_layer
      DO l=1, n_profile
        k_total=k_grey_tot(l, i)+k_gas_abs(l, i)
        tau(l, i)=k_total*d_mass(l, i)
        omega(l, i)=k_ext_scat(l, i)/(k_total+TINY(omega))
        omega(l, i)                                                     &
          =MIN(omega(l, i), 1.0e+00_RealK-3.2e+01_RealK*eps_r)
      END DO
    END DO

  ELSE IF (i_scatter_method_band == ip_no_scatter_abs) THEN

!   The scattering extinction is ignored completely, so
!   only the absorptive contributions to the single
!   scattering properties are included. If full scattering
!   is not to be used in the IR this is normally the appropriate
!   approximation as scattering is still dominated by the
!   forward peak.

    DO i=i_first_layer, i_last_layer
      DO l=1, n_profile
        tau(l, i)=(k_grey_tot(l, i)+k_gas_abs(l, i)                     &
          -k_ext_scat(l, i))*d_mass(l, i)
        omega(l, i)=0.0e+00_RealK
      END DO
    END DO

  ELSE IF (i_scatter_method_band == ip_no_scatter_ext) THEN

!   The scattering extinction is added on to the absorption.
!   This option is usually a bad approximation to the effects
!   of scattering in the IR, but may occasionally be appropriate
!   if the asymmetry is low.

    DO i=i_first_layer, i_last_layer
      DO l=1, n_profile
        tau(l, i)=(k_grey_tot(l, i)+k_gas_abs(l, i))                    &
          *d_mass(l, i)
        omega(l, i)=0.0e+00_RealK
      END DO
    END DO

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE single_scattering
END MODULE single_scattering_mod
