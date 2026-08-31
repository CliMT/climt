! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set methods of weighting optical properties.
!
MODULE method_weight_pcf
!
!
! Description:
!
! This module contains the default sizes for spectral arrays. In the
! course of time it should become redundant.
!
!- End of header
!
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, Parameter :: IP_weight_thin       = 1
!   Weighting in thin limit
  INTEGER, Parameter :: IP_weight_thick      = 2
!   Weighting in thick limit
  INTEGER, Parameter :: IP_weight_thin_thick = 3
!   Combined thin and thick weighting
!
END MODULE method_weight_pcf
