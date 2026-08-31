! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set arrays describing the spherical truncation.
!
! Purpose:
!   This routine sets an arrays of pointers to control the allocation
!   of memory.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
MODULE set_truncation_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_TRUNCATION_MOD'
CONTAINS
SUBROUTINE set_truncation(ierr                                          &
    , i_truncation, ls_global_trunc                                     &
    , ls_max_order, ls_local_trunc                                      &
    , ms_min, ms_max, ms_trunc                                          &
    , ia_sph_mm, n_order_phase                                          &
    , nd_max_order                                                      &
    )

  USE rad_pcf, ONLY: i_err_fatal, ip_trunc_adaptive, ip_trunc_azim_sym, &
                     ip_trunc_rhombohedral, ip_trunc_triangular
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


  INTEGER, INTENT(IN) ::                                                &
      nd_max_order
!       Size allocated for orders of spherical harmonics

! Dummy arguments
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      i_truncation                                                      &
!       Type of spherical truncation
    , ls_global_trunc                                                   &
!       Global order of truncation
    , ms_min                                                            &
!       Lowest value of of the azimuthal order calculated
    , ms_max
!       Highest value of of the azimuthal order calculated
  INTEGER, INTENT(OUT) ::                                               &
      ls_max_order                                                      &
!       Maximum order of spherical harmonic terms required
    , ls_local_trunc(0: nd_max_order)                                   &
!       Truncating order for individual azimuthal quantum numbers
    , ms_trunc(0: nd_max_order)                                         &
!       Maximum azimuthal quantum number for each order
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficients of (m, m) for each m
    , n_order_phase
!       Order of terms in the phase function to be used in
!       direct calculation of spherical harmonics


! Local Variables
  INTEGER                                                               &
      ls                                                                &
!       Order of spherical harmonic
    , ms
!       Azimuthal order of spherical harmonic

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SET_TRUNCATION'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Carry out a preliminary check that the truncation is appropriate.
  IF ( (i_truncation == ip_trunc_azim_sym).AND.                         &
     (ms_max >  0) ) THEN

    cmessage = '*** Error: An azimuthally symmetric truncation is not ' &
      //'appropriate if MS_MAX > 0.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END IF

  IF ( (i_truncation == ip_trunc_triangular).OR.                        &
     (i_truncation == ip_trunc_adaptive) ) THEN

!   Set the characteristics for a triangular truncation:
!   azimuthal orders from MS_MIN to MS_MAX are needed.
!   In order to keep an even number of harmonics for all
!   azimuthal orders, the maximum order must be set
!   1 greater than the (odd) order of the truncation. In
!   addition, an extra order beyond the truncation for the
!   particular m is required for the case of solar radiation.
!   ("+4" appears in the loop below because MS is one greater
!   then the order for which space is being reserved.) Note
!   finally that space is allocated for the unused harmonic
!   LS_GLOBAL_TRUNC+2 for even values of MS just to keep the
!   programming simple.

!   The adaptive truncation comes here as well since the maximum
!   conceivable number of harmonics might be required for each
!   azimuthal order.
    ls_max_order=ls_global_trunc+1
    ms_trunc(ms_min)=ms_min
    DO ls=ms_min+1, ls_max_order
      ms_trunc(ls)=MIN(MIN(ms_max, ls), ls_global_trunc)
    END DO
    ia_sph_mm(ms_min)=1
    DO ms=ms_min+1, ms_max
      ia_sph_mm(ms)=ia_sph_mm(ms-1)+ls_global_trunc+4-ms
    END DO
!   For each MS an even number of terms must be calculated. The
!   global truncation will be odd.
    DO ms=ms_min, ms_max
      ls_local_trunc(ms)=ls_global_trunc+MOD(ms, 2)
    END DO

  ELSE IF (i_truncation == ip_trunc_rhombohedral) THEN

!   Set the characteristics for a rhombohedral truncation.
!   If calculation begins with an odd azimuthal order, one
!   extra order will be required to ensure that even of polar
!   orders are calculated.
    ls_max_order=ls_global_trunc+MOD(ms_min, 2)
!   N.B. If LS_MAX_ORDER is ever changed, make sure that the
!   following code is still valid. Here LS_MAX_ORDER logically
!   means LS_GLOBAL_TRUNC+MOD(MS_MIN, 2).

!   Reset the maximum azimuthal order if it has been set too high
    ms_trunc(ms_min)=ms_min
    DO ls=ms_min+1, ls_max_order
      ms_trunc(ls)=MIN(ms_max, MIN(ls_max_order-ls+ms_min, ls))
    END DO
    ia_sph_mm(ms_min)=1
!   The "+4" rather than "+3" below allows for the fact that the
!   requisite number of harmonics does not fall off exactly
!   linearly with the azimuthal order, but does so in steps.
    DO ms=ms_min+1, ms_max
      ia_sph_mm(ms)=ia_sph_mm(ms-1)                                     &
        +ls_global_trunc+4+ms_min-2*ms+1
    END DO
!   For each MS an even number of terms must be calculated. The
!   global truncation will be odd.
    DO ms=ms_min, ms_max
      ls_local_trunc(ms)=ls_global_trunc                                &
        +MOD(ms_min, 2)-(ms-ms_min)
    END DO

  ELSE IF (i_truncation == ip_trunc_azim_sym) THEN

!   Set the characteristics for an azimuthally symmetric truncation.
!   This will be the normal case in the infra-red region.
    ls_max_order=ls_global_trunc
    ms_trunc(0)=0
    ia_sph_mm(0)=1
    DO ls=1, ls_max_order
      ms_trunc(ls)=0
    END DO
!   Set the address of one extra order for use with solar sources.
    ls_local_trunc(0)=ls_global_trunc

  ELSE

    cmessage = '***Error: An illegal truncation has been requested.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END IF

! Calculate enough terms of the phase function to satsify the
! truncation at every azimuthal order.
  n_order_phase=1
  DO ms=ms_min, ms_max
    n_order_phase=MAX(n_order_phase, ls_local_trunc(ms))
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_truncation
END MODULE set_truncation_mod
