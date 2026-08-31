! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set unit numbers for the input of data

MODULE def_data_in_icf

  IMPLICIT NONE

  INTEGER, Parameter :: iu_raw_in   = 70
!           Unit number for raw data 
  INTEGER, Parameter :: iu_data_out = 71
!           Unit number for output
  INTEGER, Parameter :: iu_filter   = 72
!           Unit number for input of the filter fuction 

END MODULE def_data_in_icf  
