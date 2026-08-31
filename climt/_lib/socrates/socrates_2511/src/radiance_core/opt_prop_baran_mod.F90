! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate optical properties for the baran parametrisation
!
!------------------------------------------------------------------------------
MODULE opt_prop_baran_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE opt_prop_baran(id_ct, &
  first_layer, last_layer, n_cloud_profile, i_cloud_profile, &
  cloud_parameter, cond_mass_frac, dim_char, t, &
  k_ext_tot, k_ext_scat, asymmetry, phase_fnc_weighted, forward_scat_weighted)

USE realtype_rd, ONLY: RealK

! Sizes of arrays
  INTEGER, INTENT(IN) :: id_ct
!       Topmost declared cloudy layer

! Dummy variables.
  INTEGER, INTENT(IN) :: first_layer, last_layer
!       First and last layers to loop over
  INTEGER, INTENT(IN) :: n_cloud_profile(id_ct:)
!       Number of cloudy profiles
  INTEGER, INTENT(IN) :: i_cloud_profile(:, id_ct:)
!       Profiles containing clouds

  REAL(RealK), INTENT(IN) :: cloud_parameter(:)
!       Cloud parameters
  REAL(RealK), INTENT(IN) :: cond_mass_frac(:, id_ct:)
!       Condensate mass mixing ratio
  REAL(RealK), INTENT(IN) :: dim_char(:, id_ct:)
!       Characteristic dimension
  REAL(RealK), INTENT(IN) :: t(:, :)
!       Temperature

  REAL(RealK), INTENT(OUT) :: k_ext_tot(:, id_ct:)
!       Total extinction
  REAL(RealK), INTENT(OUT), OPTIONAL :: k_ext_scat(:, id_ct:)
!       Scattering extinction
  REAL(RealK), INTENT(OUT), OPTIONAL :: asymmetry(:, id_ct:)
!       Asymmetry
  REAL(RealK), INTENT(OUT), OPTIONAL :: phase_fnc_weighted(:, id_ct:, :)
!       Phase function weighted by scattering extinction
  REAL(RealK), INTENT(OUT), OPTIONAL :: forward_scat_weighted(:, id_ct:)
!       Forward scattering weighted by scattering extinction

! Local variables.
  INTEGER :: l, ll, i
  REAL(RealK) :: x, omega

  IF (PRESENT(forward_scat_weighted) .AND. &
      PRESENT(phase_fnc_weighted) .AND. &
      PRESENT(asymmetry) .AND. &
      PRESENT(k_ext_scat)) THEN
    DO i=first_layer, last_layer
      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        call calc_k_ext_tot()
        call calc_k_ext_scat()
        call calc_asymmetry()
        call calc_forward_scat()
      END DO
    END DO
  ELSE IF (PRESENT(phase_fnc_weighted) .AND. &
           PRESENT(asymmetry) .AND. &
           PRESENT(k_ext_scat)) THEN
    DO i=first_layer, last_layer
      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        call calc_k_ext_tot()
        call calc_k_ext_scat()
        call calc_asymmetry()
      END DO
    END DO
  ELSE IF (PRESENT(k_ext_scat)) THEN
    DO i=first_layer, last_layer
      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        call calc_k_ext_tot()
        call calc_k_ext_scat()
      END DO
    END DO
  ELSE
    DO i=first_layer, last_layer
      DO ll=1, n_cloud_profile(i)
        l=i_cloud_profile(ll, i)
        call calc_k_ext_tot()
      END DO
    END DO
  END IF

CONTAINS

SUBROUTINE calc_k_ext_tot()
IMPLICIT NONE
  x = MIN(cond_mass_frac(l, i) + 1.0e-12_RealK, 1.0e-01_RealK)
  k_ext_tot(l, i) = cloud_parameter(1) * x / (t(l, i)**4)
END SUBROUTINE calc_k_ext_tot

SUBROUTINE calc_k_ext_scat()
IMPLICIT NONE
  x = MIN(cond_mass_frac(l, i) + 1.0e-12_RealK, &
          1.0e-03_RealK - EPSILON(1.0_RealK))
  omega = MIN(cloud_parameter(2) &
           + (cloud_parameter(3) * x * t(l, i)), &
              1.0_RealK - EPSILON(1.0_RealK))
  k_ext_scat(l, i) = k_ext_tot(l, i) * omega
END SUBROUTINE calc_k_ext_scat

SUBROUTINE calc_asymmetry()
IMPLICIT NONE
  asymmetry(l, i) = MIN(cloud_parameter(4) &
                     + (cloud_parameter(5) * x * t(l, i)), &
                        1.0_RealK - EPSILON(1.0_RealK))
  phase_fnc_weighted(l, i, 1) = k_ext_scat(l, i) * asymmetry(l, i)
END SUBROUTINE calc_asymmetry

SUBROUTINE calc_forward_scat()
IMPLICIT NONE
  forward_scat_weighted(l, i) = phase_fnc_weighted(l, i, 1) * asymmetry(l, i)
END SUBROUTINE calc_forward_scat

END SUBROUTINE opt_prop_baran
END MODULE opt_prop_baran_mod
