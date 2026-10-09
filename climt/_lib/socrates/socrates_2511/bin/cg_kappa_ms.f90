! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the product of CG and KAPPA.
!
! Purpose:
!   This routine calculates a sum of products of hemispheric
!   integrals and Clebsch-Gordan coefficients used in the
!   specification of BRDFs. These terms might be stored, but
!   this could involve the use of too much memory.
!
! Method:
!   Indexing is a bit complicated. The BRDF is truncated at an
!   even order and l'+m and l+m must both be even, so the effective
!   truncation when m is odd is one order lower. Nominally,
!   CGK is CGK(l',l) where l' takes alternate values though l takes
!   consecutive values, so we map into the actual array as
!       (l', l) --> ( (l'-m+2)/2, (l-m+1) )
!
!- ---------------------------------------------------------------------
MODULE cg_kappa_ms_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CG_KAPPA_MS_MOD'
CONTAINS
SUBROUTINE cg_kappa_ms(ms, ls_trunc, ls_brdf_trunc                      &
    , cg_coeff, kappa                                                   &
    , cgk                                                               &
    , nd_max_order, nd_brdf_trunc                                       &
    )


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_brdf_trunc
!       Size allocated for orders of spherical harmonics
!       in BRDFs

! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      ms                                                                &
!       Azimuthal order
    , ls_trunc                                                          &
!       The order of truncation applied to spherical harmonics
    , ls_brdf_trunc
!       The order of truncation applied to the BRDF
  REAL (RealK), INTENT(IN) ::                                           &
      cg_coeff(ls_trunc-ms+1)                                           &
!       Clebsch-Gordan coefficients
    , kappa(nd_max_order/2, nd_max_order/2)
!       Integrals of pairs of spherical harmonics over the downward
!       hemisphere

  REAL (RealK), INTENT(OUT) ::                                          &
      cgk(nd_brdf_trunc/2+1, nd_max_order)
!       Integrals of pairs of spherical harmonics over the downward
!       hemisphere

! Local variables:
  INTEGER                                                               &
      lsr_p                                                             &
!       Reduced primed polar order
    , lsr
!       Reduced polar order

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CG_KAPPA_MS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Consider first the case where l+m is even. In this case the
! documented formula is applied directly, with the omission
! of a term in the first element of the array where l' would
! be out of range.
  DO lsr=1, ls_trunc+1-ms, 2
    cgk(1, lsr)=cg_coeff(1)*kappa(1, (lsr+1)/2)
    DO lsr_p=3, ls_brdf_trunc-ms+1-MOD(ms, 2), 2
      cgk((lsr_p+1)/2, lsr)                                             &
        =cg_coeff(lsr_p)*kappa((lsr_p+1)/2, (lsr+1)/2)                  &
        +cg_coeff(lsr_p-1)*kappa((lsr_p-1)/2, (lsr+1)/2)
    END DO
  END DO
! In the case where l+m is odd the array is generally zero, so
! we initialize all such values and calculate exceptional cases
! later. Note that KAPPA does not appear in these loops because
! in the compressed format these trivially evaluated values are
! not stored.
  DO lsr=2, ls_trunc+1-ms, 2
    DO lsr_p=1, ls_brdf_trunc-ms+1-MOD(ms, 2), 2
      cgk((lsr_p+1)/2, lsr)=0.0e+00_RealK
    END DO
  END DO
! The case l=l'+1:
  DO lsr=2, ls_brdf_trunc-ms-MOD(ms, 2)+2, 2
    cgk(lsr/2, lsr)=cg_coeff(lsr-1)*0.5e+00_RealK
  END DO
! The case l=l'-1:
  DO lsr=2, ls_brdf_trunc-ms-MOD(ms, 2), 2
    cgk(lsr/2+1, lsr)=cg_coeff(lsr)*0.5e+00_RealK
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cg_kappa_ms
END MODULE cg_kappa_ms_mod
