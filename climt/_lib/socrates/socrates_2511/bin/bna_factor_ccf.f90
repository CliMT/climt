! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set the BNA factor for ADT calculations.
!
MODULE bna_factor_ccf
!
! Description:
!
! This module defines the BNA factor used to correct ADT.
!
!- End of header
!
!
!  Modules used
   USE realtype_rd
!
!
   IMPLICIT NONE
!
!
!
   REAL  (RealK), Parameter :: factor_bna = 1.7_RealK
!
!
!
END MODULE bna_factor_ccf
