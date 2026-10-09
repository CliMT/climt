! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set parameters for conjugate gradient algorithm.
!
MODULE scale_parameters_acf
!
! Description:
!
! This module contains the settings of parameters used in the
! conjugate gradient algorithm which defines scaling functions.
!
!- End of header
!
! Modules used:
  USE realtype_rd
!
!
!
  INTEGER, Parameter :: NP_max_iteration_cg = 100
!   Maximum number of conjugate gradient cycles
  INTEGER, Parameter :: NP_max_line_search  = 100
!   Maximum number of line searches
  REAL  (RealK) :: s_max_ratio = 2.0E+04_RealK
!   Maximum line search parameter
!
END MODULE scale_parameters_acf
