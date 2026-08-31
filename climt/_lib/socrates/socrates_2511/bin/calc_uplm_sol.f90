! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate Upsilon_l^m(0) for the solar direction.
!
! Purpose:
!   This routine is called to determine the values of spherical
!   harmonics in the solar direction.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used.
!
!- ---------------------------------------------------------------------
MODULE calc_uplm_sol_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_UPLM_SOL_MOD'
CONTAINS
SUBROUTINE calc_uplm_sol(n_profile, ms_min, ms_max, ia_sph_mm           &
    , ls_local_trunc, zen_0, uplm_sol                                   &
    , nd_profile, nd_max_order, nd_sph_coeff)


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_sph_coeff                                                      &
!       Number of spherical coefficients
    , nd_max_order
!       Maximum order of calculation
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of atmospheric profiles
    , ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical oefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)
!       Local truncation at this order
  REAL (RealK), INTENT(IN) ::                                           &
      zen_0(nd_profile)
!       Cosines of solar zenith angles
  REAL (RealK), INTENT(OUT) ::                                          &
      uplm_sol(nd_profile, nd_sph_coeff)
!       Array of Upsilon_l^m evaluated in the solar direction.


! Local variables
  INTEGER                                                               &
      ls                                                                &
!       Order of harmonic
    , ms                                                                &
!       Azimuthal quantum number of harmonic
    , j                                                                 &
!       Temporary address
    , k                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      product                                                           &
!       Factorial terms of Y_lm
    , lsr                                                               &
!       Real polar order of harmonic
    , msr
!       Real azimuthal order of harmonic

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_UPLM_SOL'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Note here that ZEN_0 holds the cosine of the zenith angle, so
! the cosine of the solar direction is actually -ZEN_0.
  DO ms=ms_min, ms_max

!   Calculate Upsilon_m^m(n_sol) to start the recurrence.
    j=ia_sph_mm(ms)

    product=1.0e+00_RealK
    msr=REAL(ms, RealK)
    IF (ms >  0) THEN
      DO k=1, ms
        product=(1.0e+00_RealK-5.0e-01_RealK/REAL(k, RealK))*product
      END DO
      DO l=1, n_profile
        uplm_sol(l, j)=(-1.0e+00_RealK)**ms                             &
          *SQRT((1.0e+00_RealK-zen_0(l)*zen_0(l))**ms*product           &
          *(2.0e+00_RealK*msr+1.0e+00_RealK)/(4.0e+00_RealK*pi))
      END DO
    ELSE
      DO l=1, n_profile
        uplm_sol(l, j)=1.0e+00_RealK/SQRT(4.0e+00_RealK*pi)
      END DO
    END IF

!   Calculate the next polar order to enable the recurrence to
!   start.
    IF (ms <= ls_local_trunc(ms)+1) THEN
      DO l=1, n_profile
        uplm_sol(l, j+1)                                                &
          =-zen_0(l)*SQRT(2.0e+00_RealK*msr                             &
          +3.0e+00_RealK)*uplm_sol(l, j)
      END DO
    END IF

!   Complete the recurrence on l.
    DO ls=ms+2, ls_local_trunc(ms)+1
      j=ia_sph_mm(ms)+ls-ms
      lsr=REAL(ls, RealK)
      DO l=1, n_profile
        uplm_sol(l, j)                                                  &
          =SQRT(((2.0e+00_RealK*lsr-1.0e+00_RealK)                      &
          *(2.0e+00_RealK*lsr+1.0e+00_RealK))                           &
          /((lsr+msr)*(lsr-msr)))                                       &
          *(-zen_0(l))*uplm_sol(l, j-1)                                 &
          -SQRT(((2.0e+00_RealK*lsr+1.0e+00_RealK)                      &
          *(lsr-1.0e+00_RealK-msr)*(lsr-1.0e+00_RealK+msr))             &
          /((2.0e+00_RealK*lsr-3.0e+00_RealK)                           &
          *(lsr-msr)*(lsr+msr)))*uplm_sol(l, j-2)
      END DO
    END DO

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_uplm_sol
END MODULE calc_uplm_sol_mod
