! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate spherical harmonics excluding expoential.
!
! Purpose:
!   Spherical harmonics, Upsilon_lm, are calculated for given directions
!   for all values of l at a fixed value of m.
!
! Method:
!   Y_mm is known so upward recurrence on l is used.
!
!- ---------------------------------------------------------------------
MODULE eval_uplm_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'EVAL_UPLM_MOD'
CONTAINS
SUBROUTINE eval_uplm(ms, n_max_order, n_direction, x                    &
     , up_lm                                                            &
     , nd_direction                                                     &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
       nd_direction
!         Maximum number of directions

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      ms                                                                &
!       Azimuthal quantum number of spherical harmonic
    , n_max_order
!       Maximum order of harmonics to calculate
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of directions
  REAL (RealK), INTENT(IN) ::                                           &
      x(nd_direction)
!       Cosines of polar angels of viewing directions
  REAL (RealK), INTENT(OUT) ::                                          &
      up_lm(nd_direction, n_max_order+1-ms)
!       Non-azimuthal parts of spherical harmonics


! Local variables
  INTEGER                                                               &
      ls                                                                &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k
!       Loop variable
  REAL (RealK) ::                                                       &
      l                                                                 &
!       Principal quantum number of harmonic
    , m                                                                 &
!       Azimuthal quantum number of harmonic
    , product
!       Factorial terms in Y_lm

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='EVAL_UPLM'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Start the recurrence for Y_mm.
  product=1.0e+00_RealK
  m=REAL(ms, RealK)
  IF (ms >  0) THEN
    DO j=1, ms
      product=(1.0e+00_RealK-5.0e-01_RealK/REAL(j, RealK))*product
    END DO
    DO k=1, n_direction
      up_lm(k, 1)=(-1.0e+00_RealK)**ms                                  &
        *SQRT((1.0e+00_RealK-x(k)*x(k))**ms*product                     &
        *(2.0e+00_RealK*m+1.0e+00_RealK)/(4.0e+00_RealK*pi))
    END DO
  ELSE
    DO k=1, n_direction
      up_lm(k, 1)=1.0e+00_RealK/SQRT(4.0e+00_RealK*pi)
    END DO
  END IF


! Calculate Y_(m+1),m if it is within bounds.
  IF (ms <  n_max_order) THEN
    DO k=1, n_direction
      up_lm(k, 2)=x(k)*SQRT(2.0e+00_RealK*m+3.0e+00_RealK)              &
        *up_lm(k, 1)
    END DO
  END IF


! Complete the recurrence on l.
  DO ls=ms+2, n_max_order
    l=REAL(ls, RealK)
    DO k=1, n_direction
      up_lm(k, ls+1-ms)=x(k)                                            &
        *SQRT(((2.0e+00_RealK*l-1.0e+00_RealK)                          &
        *(2.0e+00_RealK*l+1.0e+00_RealK))                               &
        /((l+m)*(l-m)))*up_lm(k, ls-ms)                                 &
        -SQRT(((2.0e+00_RealK*l+1.0e+00_RealK)                          &
        *(l-1.0e+00_RealK-m)*(l-1.0e+00_RealK+m))                       &
        /((2.0e+00_RealK*l-3.0e+00_RealK)*(l+m)*(l-m)))                 &
        *up_lm(k, ls-1-ms)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eval_uplm
END MODULE eval_uplm_mod
