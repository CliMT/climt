! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define types of residual for fitting transmissions.
!
MODULE type_residual_pcf
!
! Description:
!   This module defines the permitted types of residual for minimzing
! errors in the scaling function.
!
!- End of header
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, Parameter :: ip_scale_trans_residual = 1
!   Minimize the residual in the transmission
  INTEGER, Parameter :: ip_scale_full_scl_res = 2
!   Minimize the residuals in the scaling
  INTEGER, Parameter :: ip_scale_mean_scl_res = 3
!   Minimize the residuals in the mean scaling
!
!
!
END MODULE type_residual_pcf
