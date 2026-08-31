! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate Upsilon_l^m(0) or its derivative.
!
! Purpose:
!   This routine is called to determine the value of a spherical
!   harmonic with theta=pi/2 and phi=0 or the derivative for
!   alternate orders. This minimizes storage.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used.
!
!- ---------------------------------------------------------------------
MODULE calc_uplm_zero_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_UPLM_ZERO_MOD'
CONTAINS
SUBROUTINE calc_uplm_zero(ms_min, ms_max, ia_sph_mm                     &
    , ls_local_trunc, uplm_zero                                         &
    , nd_max_order, nd_sph_coeff)


  USE realtype_rd, ONLY: RealK
  USE rad_ccf, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      nd_sph_coeff                                                      &
!       Number of spherical coefficients
    , nd_max_order
!       Maximum order of calculation
  INTEGER, INTENT(IN) ::                                                &
      ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)
!       Local truncation at this order
  REAL (RealK), INTENT(OUT) ::                                          &
      uplm_zero(nd_sph_coeff)
!       Array of Upsilon_l^m and derivatives at polar angles of pi/2


! Local variables
  INTEGER                                                               &
      ls                                                                &
!       Order of harmonic
    , ms                                                                &
!       Azimuthal quantum number of harmonic
    , j                                                                 &
!       Temporary address
    , k
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_UPLM_ZERO'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO ms=ms_min, ms_max

!   Calculate Upsilon_m^m(0) to start the recurrence.
    j=ia_sph_mm(ms)
    uplm_zero(j)=1.0e+00_RealK/(4.0e+00_RealK*pi)
    DO k=3, 2*ms+1, 2
       uplm_zero(j)=uplm_zero(j)*REAL(k, RealK)/REAL(k-1, RealK)
    END DO
    uplm_zero(j)=REAL(1-2*MOD(ms, 2), RealK)*SQRT(uplm_zero(j))


!   Calculate DUpsilon_{m+1}^m(0) to start the recurrence for the
!   derivatives.
    j=j+1
    uplm_zero(j)=3.0e+00_RealK/(4.0e+00_RealK*pi)
    DO k=5, 2*ms+3, 2
       uplm_zero(j)=uplm_zero(j)*REAL(k, RealK)/REAL(k-3, RealK)
    END DO
    uplm_zero(j)=REAL(1-2*MOD(ms, 2), RealK)*SQRT(uplm_zero(j))

!   Now apply the recurrence formulae:
    DO ls=ms+2, ls_local_trunc(ms)-1, 2
      j=ia_sph_mm(ms)+ls-ms

!     Recurrence for Upsilon_l^m.
      uplm_zero(j)=-uplm_zero(j-2)                                      &
        *SQRT(REAL((2*ls+1)*(ls+ms-1)*(ls-ms-1), RealK)                 &
        /REAL((2*ls-3)*(ls+ms)*(ls-ms), RealK))

!     Recurrence for the derivative Upsilon_(l+1)^m.
      uplm_zero(j+1)=-uplm_zero(j-1)                                    &
        *SQRT(REAL((2*ls+3)*(ls+1+ms)*(ls+1-ms), RealK)                 &
        /REAL((2*ls-1)*(ls+ms)*(ls-ms), RealK))

    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_uplm_zero
END MODULE calc_uplm_zero_mod
