! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate temperature interpolation factors for generalised
! continuum k-terms.
!
!-----------------------------------------------------------------------

MODULE inter_t_lookup_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTER_T_LOOKUP_MOD'
CONTAINS
SUBROUTINE inter_t_lookup(nd_profile, nd_layer, nd_t_lookup_cont &
     , n_profile, n_layer, t, t_lookup_cont &
     , wt_ct, jt_ct)

  USE realtype_rd, ONLY: RealK
  USE yomhook,     ONLY: lhook, dr_hook
  USE parkind1,    ONLY: jprb, jpim

  IMPLICIT NONE


  INTEGER, INTENT(IN) :: &
       nd_profile &
!       Max number of profile
     , nd_layer &
!       Max number of layer
     , nd_t_lookup_cont &
!       Number of lookup temperatures
     , n_layer &
!       Number of layers
     , n_profile
!       Number of profiles

  REAL (RealK), INTENT(IN) :: &
       t(nd_profile, nd_layer) &
!       Layer temperature
     , t_lookup_cont(nd_t_lookup_cont)
!       Lookup table temperatures


  REAL (RealK), INTENT(OUT) :: &
       wt_ct(nd_profile, nd_layer)
!       Weight of jt_low-term in temperature interpolation.

  INTEGER, INTENT(OUT) :: &
       jt_ct(nd_profile, nd_layer)
!       Index of reference temperature at level I such that the actual
!       temperature is between JT and JT+1

! Local variables
  REAL (RealK) :: &
       t_layer
!        Inter-medium variable

  INTEGER :: &
       i, l
!       Vertical and horizontal loop index

  REAL, PARAMETER :: eps = EPSILON(1.0_RealK)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='INTER_T_LOOKUP'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


  DO i=1, n_layer
    DO l=1, n_profile
!     Find the reference temperature on the lower side of the
!     layer temperature for each layer. Store these indices in jt_ct.
!     Store in wt_ct the fraction of the way between the layer
!     temperature and the next highest reference temperature.
      t_layer=MIN(MAX(t(l,i), t_lookup_cont(1)*(1.0_RealK + eps)), &
                  t_lookup_cont(nd_t_lookup_cont)*(1.0_RealK - eps))

      jt_ct(l,i) = MINLOC(t_layer - t_lookup_cont(:), 1, &
                          t_layer - t_lookup_cont(:) >= 0.0_RealK)

      wt_ct(l,i) = 1.0_RealK - (t_layer - t_lookup_cont(jt_ct(l,i)))/ &
                  (t_lookup_cont(jt_ct(l,i)+1) - t_lookup_cont(jt_ct(l,i)))
    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE inter_t_lookup
END MODULE inter_t_lookup_mod
