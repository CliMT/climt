! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module setting maximum permitted size of scatterers
!
MODULE max_size_acf
!
! Description:
!
! This module defines the maximum dimension of scattering
! object used in the code.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd
!
  REAL  (RealK), Parameter :: max_dimen_scatterer = 4.0E-03_RealK
!           Maximum size of scatterer
!
END MODULE max_size_acf
