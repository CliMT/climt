! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set offset for the residual in the scaling.
!
MODULE offset_residual_trans_acf
!
! Description:
!   This module defines the offset to the scaling of the residual.
! This parameter may require adjustment.
!
!- End of header
!
!
  USE realtype_rd
!
  IMPLICIT NONE
!
!
!
  REAL  (RealK), Parameter :: offset_residual = 1.0E-02_RealK
!   Offset to scaling in residual to avoid excessive weight being
!   applied to errors when the scaling function is close to 0.
!
!
!
END MODULE offset_residual_trans_acf
