! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to get a free unit number
!
!- ---------------------------------------------------------------------
SUBROUTINE get_free_unit(ierr, iunit)

  USE rad_pcf, ONLY: i_normal, i_err_io
  USE ereport_mod, ONLY: ereport

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: ierr  ! Error flag
  INTEGER, INTENT(OUT)   :: iunit ! Unit number
  LOGICAL :: l_open               ! Flag for open unit

  CHARACTER (LEN=80)            :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'get_free_unit'

  ierr=i_normal
  iunit=20
  INQUIRE(unit=iunit, opened=l_open)
  DO WHILE ( (l_open).AND.(iunit < 100) )
    IF (l_open) iunit=iunit+1
    INQUIRE(unit=iunit, opened=l_open)
  END DO

  IF (iunit > 100) THEN
    cmessage = 'No free units are available for i/o.'
    ierr=i_err_io
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

END SUBROUTINE get_free_unit
