! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate hemispheric spherical integrals.
!
! Purpose:
!   This routine calculates hemispheric integrals of products
!   of spherical harmonics for a fixed azimuthal order for use
!   in Marshak's boundary conditions.
!
! Method:
!   We require the integral of Y_l'^m* Y_l^m over the downward
!   hemisphere for those l' such that l'+m is odd. If l=l' the
!   integral is trivially 1/2, but otherwise it will be zero
!   unless l+l' is odd. To reduce storage we omit the case l=l'
!   here and then map
!   (l', l) --> ( (l'-m+1)/2, (l-m+2)/2)
!   in the stored array.
!
!   The integrals are evaluated from the values of spherical
!   harmonics and their derivatives at a polar angle of pi/2.
!
!- ---------------------------------------------------------------------
MODULE hemi_sph_integ_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'HEMI_SPH_INTEG_MOD'
CONTAINS
SUBROUTINE hemi_sph_integ(ls_trunc, ms, uplm_zero                       &
    , kappa                                                             &
    , nd_max_order                                                      &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_max_order
!       Size allocated for orders of spherical harmonics

! Dummy arguments
  INTEGER, INTENT(IN) ::                                                &
      ls_trunc                                                          &
!       The truncating order of the system of equations
    , ms
!       Azimuthal order
  REAL (RealK), INTENT(IN) ::                                           &
      uplm_zero(ls_trunc+1-ms)
!       Spherical harmonics and derivatives at a polar angle of
!       pi/2

  REAL (RealK), INTENT(OUT) ::                                          &
      kappa(nd_max_order/2, nd_max_order/2)
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

  CHARACTER(LEN=*), PARAMETER :: RoutineName='HEMI_SPH_INTEG'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The outer loop is over l' where l'+m is odd. Indexing is done
! using the reduced indices l'+1-m and l+1-m.
  DO lsr_p=2, ls_trunc+1-ms, 2
    DO lsr=1, ls_trunc-ms, 2
      kappa(lsr_p/2, (lsr+1)/2)=2.0e+00_RealK*pi                        &
        *REAL(1-2*MOD(lsr_p, 2), RealK)                                 &
        *uplm_zero(lsr)*uplm_zero(lsr_p)                                &
        /REAL((lsr-lsr_p)*(lsr+lsr_p-1+2*ms), RealK)
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE hemi_sph_integ
END MODULE hemi_sph_integ_mod
