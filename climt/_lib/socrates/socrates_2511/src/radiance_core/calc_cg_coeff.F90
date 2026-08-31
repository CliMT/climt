! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate spherical Clebsch-Gordan coefficients.
!
! Purpose:
!   The routine yields the Clebsch-Gordan coefficients between the
!   spherical harmonics Y_l^m and Y_1^0, c_{lm}^+ in the notation of
!   the description of the algorithm, or <l+1,m|1,0,l,m> in standard
!   notation. These are stored in one array with addressing determined
!   by the truncation.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used. Only values for m>0 are required.
!
!- ---------------------------------------------------------------------
MODULE calc_cg_coeff_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_CG_COEFF_MOD'
CONTAINS
SUBROUTINE calc_cg_coeff(ls_max_order                                   &
    , ia_sph_mm, ms_min, ms_trunc                                       &
    , cg_coeff                                                          &
    , nd_max_order, nd_sph_coeff)


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_sph_coeff
!       Size of array of spherical coefficients
  INTEGER, INTENT(IN) ::                                                &
      ls_max_order                                                      &
!       Maximum order of harmonics required
    , ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_trunc(0: nd_max_order)                                         &
!       Truncation in MS for this order
    , ia_sph_mm(0: nd_max_order)
!       Position of Clebsh-Gordan coefficient with m=0 for the
!       given value of l.
  REAL (RealK), INTENT(OUT) ::                                          &
      cg_coeff(nd_sph_coeff)


! Local variables
  INTEGER                                                               &
      ls                                                                &
!      Order of harmonic
    , ms
!       Azimuthal order of harmonic
  REAL (RealK) ::                                                       &
      inv
!       l-dependent denominator

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_CG_COEFF'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO ls=0, ls_max_order
    inv=1.0e+00_RealK/REAL((2*ls+1)*(2*ls+3), RealK)
    DO ms=ms_min, ms_trunc(ls)
      cg_coeff(ia_sph_mm(ms)+ls-ms)                                     &
        =SQRT(REAL((ls+1-ms)*(ls+1+ms), RealK)*inv)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_cg_coeff
END MODULE calc_cg_coeff_mod
