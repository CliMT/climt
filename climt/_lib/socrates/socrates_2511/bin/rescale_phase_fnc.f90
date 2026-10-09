! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to rescale the phase function.
!
! Method:
!   The standard rescaling of the phase function is used.
!
!- ---------------------------------------------------------------------
MODULE rescale_phase_fnc_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'RESCALE_PHASE_FNC_MOD'
CONTAINS
SUBROUTINE rescale_phase_fnc(n_profile                                  &
    , i_layer_first, i_layer_last, n_direction, cos_sol_view            &
    , n_order_phase, phase_fnc, forward_scatter, forward_solar          &
    , l_rescale_solar_phf, n_order_phase_solar, phase_fnc_solar         &
    , nd_profile, nd_radiance_profile, nd_layer, id_1                   &
    , nd_direction, nd_max_order                                        &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_radiance_profile                                               &
!       Size allocated for atmospheric profiles used specifically
!       for calculating radiances
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , id_1                                                              &
!       Topmost declared layer for optical properties
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_direction
!       Size allocated for viewing directions

! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to rescale
    , i_layer_last                                                      &
!       Last layer to rescale
    , n_direction
!       Number of directions
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale_solar_phf
!       Flag to rescale the singly scattered solar phase function
  REAL (RealK), INTENT(IN) ::                                           &
      forward_scatter(nd_profile, id_1: nd_layer)                       &
!       Forward scattering
    , forward_solar(nd_profile, id_1: nd_layer)                         &
!       Forward scattering for the solar beam
    , cos_sol_view(nd_radiance_profile, nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing directions
  INTEGER, INTENT(IN) ::                                                &
      n_order_phase                                                     &
!       Order of terms in the phase function to be retained
    , n_order_phase_solar
!       Order of terms retained in treating singly scattered
!       solar radiation
  REAL (RealK), INTENT(INOUT) ::                                        &
      phase_fnc(nd_profile, id_1: nd_layer, nd_max_order)               &
!       Phase function
    , phase_fnc_solar(nd_radiance_profile, id_1: nd_layer               &
        , nd_direction)
!       The phase function for single scattered solar radiation

! Local variables
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , ls
!       Loop variable

! Legendre polynomials:
  REAL (RealK) ::                                                       &
      cnst1                                                             &
!       Constant in recurrence for Legendre polynomials
    , p_legendre_ls(nd_radiance_profile)                                &
!       Legendre polynomial at the current order
    , p_legendre_ls_m1(nd_radiance_profile)                             &
!       Legendre polynomial at the previous order
    , p_legendre_tmp(nd_radiance_profile)
!       Temporary Legendre polynomial

  REAL (RealK) ::                                                       &
      peak(nd_profile)
!       Forward scattering peak

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='RESCALE_PHASE_FNC'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO k=1, n_order_phase
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        phase_fnc(l, i, k)                                              &
          =(phase_fnc(l, i, k)-forward_scatter(l, i))                   &
          /(1.0e+00_RealK-forward_scatter(l, i))
      END DO
    END DO
  END DO


  IF (l_rescale_solar_phf) THEN

    DO id=1, n_direction

!     As usual we do not store Legendre polynomials:
      DO l=1, n_profile
        p_legendre_ls_m1(l)=1.0e+00_RealK
        p_legendre_ls(l)=cos_sol_view(l, id)
        peak(l)=1.0e+00_RealK+p_legendre_ls(l)*REAL(2*1+1, RealK)
      END DO

      DO ls=2, n_order_phase_solar
!       Calculate higher orders by recurrences.
        cnst1=1.0e+00_RealK-1.0e+00_RealK/REAL(ls, RealK)
        DO l=1, n_profile
          p_legendre_tmp(l)=p_legendre_ls(l)
          p_legendre_ls(l)                                              &
            =(1.0e+00_RealK+cnst1)*p_legendre_ls(l)                     &
            *cos_sol_view(l, id)-cnst1*p_legendre_ls_m1(l)
          p_legendre_ls_m1(l)=p_legendre_tmp(l)
          peak(l)=peak(l)+REAL(2*ls+1, RealK)*p_legendre_ls(l)
        END DO
      END DO

!     This is not precisely a rescaling because we do not conserve
!     the forward peak, but what is calculated is what contributes
!     to scattered radiation outside the aureole.
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          phase_fnc_solar(l, i, id)=(phase_fnc_solar(l, i, id)          &
            -forward_solar(l, i)*peak(l))                               &
            /(1.0e+00_RealK-forward_scatter(l, i))
        END DO
      END DO

    END DO

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE rescale_phase_fnc
END MODULE rescale_phase_fnc_mod
