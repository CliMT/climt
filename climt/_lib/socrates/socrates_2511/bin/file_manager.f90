! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Module containing the file management types and routines
!              which handle the assignment and release of file units
!              and storage of metadata relating to files
!              (These are simplified versions of the routines required
!               from the UM.) 
!------------------------------------------------------------------------------
MODULE file_manager

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE assign_file_unit(filename, iunit, handler)
IMPLICIT NONE

CHARACTER(*),           INTENT(IN)  :: filename
INTEGER,                INTENT(OUT) :: iunit
CHARACTER(*),           INTENT(IN)  :: handler

INTEGER :: errorstatus
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'file_manager:assign_file_unit'
LOGICAL :: l_open

iunit=20
INQUIRE(UNIT=iunit, OPENED=l_open)
DO WHILE ( (l_open).AND.(iunit < 100) )
  IF (l_open) iunit=iunit+1
  INQUIRE(UNIT=iunit, OPENED=l_open)
END DO

IF (iunit > 100) THEN
  cmessage = 'No free units are available for i/o.'
  errorstatus = 1
  CALL ereport(RoutineName, errorstatus, cmessage)
END IF

END SUBROUTINE assign_file_unit

!------------------------------------------------------------------------------

SUBROUTINE release_file_unit(iunit, handler)
IMPLICIT NONE

INTEGER,      INTENT(IN) :: iunit
CHARACTER(*), INTENT(IN) :: handler

! Dummy routine

END SUBROUTINE release_file_unit

!------------------------------------------------------------------------------

END MODULE file_manager
