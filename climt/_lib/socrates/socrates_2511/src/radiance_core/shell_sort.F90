! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to perform a shell sort.
!
! Method:
!   The standard shell sorting algorithm is implemented.
!
!- ---------------------------------------------------------------------
MODULE shell_sort_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHELL_SORT_MOD'
CONTAINS
SUBROUTINE shell_sort(n, pointer, key)


  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n
!       Number of elements
  INTEGER, INTENT(INOUT) ::                                             &
      pointer(n)
!       Pointer to succeeding elements
  REAL (RealK), INTENT(IN) ::                                           &
      key(n)
!       Key for sorting

! Local variables.
  INTEGER                                                               &
      gap                                                               &
!       Searching interval
    , pointer_temp                                                      &
!       Temporary value of pointer
    , j                                                                 &
!       Loop variable
    , k
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='SHELL_SORT'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (n == 1) THEN
    pointer(1)=1
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  gap=n
  DO WHILE(gap >= 2)
    gap=gap/2
    DO j=gap, n-1
      DO k=j-gap+1, 1, -gap
        IF (key(pointer(k)) >  key(pointer(k+gap))) THEN
                pointer_temp=pointer(k)
          pointer(k)=pointer(k+gap)
          pointer(k+gap)=pointer_temp
        END IF
      END DO
    END DO
  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE shell_sort
END MODULE shell_sort_mod
