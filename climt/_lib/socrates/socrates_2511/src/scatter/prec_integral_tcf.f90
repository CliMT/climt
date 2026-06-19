! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set tolerances for integration routines.
!
MODULE prec_integral_tcf
!
! Description:
!
! This module defines the tolerances used to assess whether
! inegrals over size distributions have converged.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd
!
  REAL  (RealK), Parameter :: tol_panel      = 1.0E-04_RealK
!         Threshold for neglect of panel
  REAL  (RealK), Parameter :: tol_refinement = 1.0E-04_RealK
!         Tolerance on Romberg step
!
END MODULE prec_integral_tcf
