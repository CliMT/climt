! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set unit numbers for standard I/O.

MODULE def_std_io_icf

  IMPLICIT NONE

  INTEGER, PARAMETER  ::  iu_stdin = 5
!                           Unit number for standard input
  INTEGER, PARAMETER  ::  iu_stdout = 6
!                           Unit number for standard output
  INTEGER, PARAMETER  ::  iu_err = 0
!                           Unit number for error messages

END MODULE def_std_io_icf
