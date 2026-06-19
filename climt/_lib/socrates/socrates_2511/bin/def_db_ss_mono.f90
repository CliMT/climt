! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module defining the geometry of ice crystals.
!
MODULE def_db_ss_mono
!
! Description:
!   This module defines the structure for monochromatic 
!   single scattering information in a database.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd
!
!
!
  TYPE STR_db_ss_mono
!
    REAL (RealK), Pointer, Dimension(:) :: dm
!     Maximum dimension of crystal
    REAL (RealK), Pointer, Dimension(:, :) :: ss
!     Square of maximum dimension of crystal
    REAL (RealK), Pointer, Dimension(:, :) :: d2_ss
!     Square of maximum dimension of crystal
!
  END TYPE STR_db_ss_mono
!
!
!
END MODULE def_db_ss_mono
