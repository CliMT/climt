! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set parameters for conjugate gradient parametrization.
!
MODULE cloud_fit_parm_acf
!
! Description:
!
! This module defines parameters used in the conjugate gradient
! fitting of cloud parametrizations.
!
!
! Modules used:
  USE realtype_rd
!
!
  IMPLICIT NONE
!
!
!- End of header
!
!
  INTEGER, Parameter :: np_max_iteration_cg = 10000
!   Maximum number of conjugate-gradient iterations
  INTEGER, Parameter :: np_max_line_search  = 100
!   Maximum number of line searches
!
  REAL  (RealK), Parameter :: s_max_ratio   = 1.0E+03_RealK
!   Maximum line search ratio
!
!
END MODULE cloud_fit_parm_acf
