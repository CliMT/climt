! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate basis functions for the diffuse albedo.
!
! Purpose:
!   This routine takes the BRDF supplied and calculates a diffuse
!   albedo for isotropic radiation, which is used to define an
!   equivalent extinction.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used. See calc_brdf.f for a note on
!   the symmetries of the BRDF and storage.
!
!- ---------------------------------------------------------------------
MODULE diff_albedo_basis_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DIFF_ALBEDO_BASIS_MOD'
CONTAINS
SUBROUTINE diff_albedo_basis(n_brdf_basis_fnc                           &
    , ls_brdf_trunc, f_brdf                                             &
    , uplm_zero                                                         &
    , diffuse_alb_basis                                                 &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_sph_coeff                    &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_brdf_basis_fnc                                                 &
!       Size allocated for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allocated for truncation of BRDFs
    , nd_sph_coeff
!       Size allocated for spherical coefficients (dimensioning
!       as ND_BRDF_TRUNC+1 might seem natural, but this can
!       lead to problems at low orders of truncation if
!       ND_BRDF_TRUNC is set too large higher in the program.


! Dummy arguments.
  REAL (RealK), INTENT(IN) ::                                           &
      uplm_zero(nd_sph_coeff)
!       Array of Upsilon_l^m and derivatives at polar angles of pi/2
  INTEGER, INTENT(IN) ::                                                &
      n_brdf_basis_fnc                                                  &
!       Number of basis functions for BRDFs
    , ls_brdf_trunc
!       Order of truncation applied to BRDFs
  REAL (RealK), INTENT(IN) ::                                           &
      f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!       Array of moments of BRDF basis functions

  REAL (RealK), INTENT(OUT) ::                                          &
      diffuse_alb_basis(nd_brdf_basis_fnc)
!       The diffuse albedo for isotropic radiation calculated
!       from the appropriate BRDF basis function


! Local variables
  INTEGER                                                               &
      lsr                                                               &
!       Reduced polar order of harmonic
    , lsr_p                                                             &
!       Reduced polar order of harmonic
    , j
!       Loop variable

  REAL (RealK) ::                                                       &
      factor(nd_brdf_trunc+1)                                           &
!       Term involving a sum over l' calculated for speed.
    , sum_p(nd_brdf_trunc+1)
!       Sum of products of the BRDF and factors over l'

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='DIFF_ALBEDO_BASIS'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO j=1, n_brdf_basis_fnc
    diffuse_alb_basis(j)=0.0e+00_RealK

    DO lsr_p=1, ls_brdf_trunc+1, 2
      factor(lsr_p)=uplm_zero(lsr_p)                                    &
        *REAL(1-2*MOD(lsr_p-1, 2), RealK)                               &
        /REAL((lsr_p-2)*(lsr_p+1), RealK)
    END DO

    DO lsr=1, ls_brdf_trunc+1, 2
      sum_p(lsr)=0.0e+00_RealK
      DO lsr_p=1, ls_brdf_trunc+1, 2
        sum_p(lsr)=sum_p(lsr)+factor(lsr_p)                             &
          *f_brdf(j, (lsr-1)/2, (lsr_p-1)/2, 0)
      END DO
      diffuse_alb_basis(j)=diffuse_alb_basis(j)                         &
        +sum_p(lsr)*uplm_zero(lsr)                                      &
        /REAL((lsr-2)*(lsr+1), RealK)
    END DO

    diffuse_alb_basis(j)=diffuse_alb_basis(j)*4.0e+00_RealK*pi

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diff_albedo_basis
END MODULE diff_albedo_basis_mod
