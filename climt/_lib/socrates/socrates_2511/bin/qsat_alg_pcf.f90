! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define algorithms for saturated mixing ratios.
!
MODULE qsat_alg_pcf
!
! Description:
!
!   This module defines the allowed algorithms for calculating
!   saturated mixing ratios.
!
!- End of header
!
!
  IMPLICIT NONE
!
!
  INTEGER, Parameter :: IP_UM_look_up  = 1
!   Look-up tables from the Unified Model
  INTEGER, Parameter :: IP_goff_gratsch = 2
!   Full Goff-Gratsch formulae
!
!
!
END MODULE qsat_alg_pcf
