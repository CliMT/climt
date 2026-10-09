! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set algorithmic parameters for correlated-k fits.
!
MODULE ck_parm_acf
!
! Description:
!   This module defines algorithmis parameters used in deriving
!   the correlated-k fits.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd

  IMPLICIT NONE


  INTEGER, PARAMETER :: n_path_kopt_default = 20
!   Default number of absorptive paths
  INTEGER, PARAMETER :: np_kopt_max_iter = 100
!   Maximum number of iterations of the Newton-Raphson algorithm
!   to find the optimal value of k 
  REAL (RealK), PARAMETER :: umin_kopt_default = 1.0E-09_RealK
!   Minimum absorptive pathlength (kg/m2)
  REAL (RealK), PARAMETER :: umax_kopt_default = 1.0E+05_RealK
!   Maximum absorptive pathlength (kg/m2)
  REAL (RealK) :: umin_kopt = umin_kopt_default
!   Minimum absorptive pathlength (kg/m2)
  REAL (RealK) :: umax_kopt = umax_kopt_default
!   Maximum absorptive pathlength (kg/m2)

END MODULE ck_parm_acf
