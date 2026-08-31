! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module defining the geometry of ice crystals.
!
MODULE def_db_crystal_geometry
!
! Description:
!   This module defines the geometrical information on crystal
!   geometry used with external databases.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd
!
!
!
  TYPE STR_db_cryst_geom
!
    INTEGER :: n_geom
!     Number of entries in the geometrical file
!
    REAL (RealK), Pointer :: dm(:)
!     Maximum dimension of crystal
    REAL (RealK), Pointer :: dm2(:)
!     Square of maximum dimension of crystal
    REAL (RealK), Pointer :: dm3(:)
!     Cube of maximum dimension of crystal
    REAL (RealK), Pointer :: proj_area(:)
!     Projected area of crystal
    REAL (RealK), Pointer :: d2_proj_area(:)
!     2nd derivative of projected area of crystal
    REAL (RealK), Pointer :: volume(:)
!     Projected area of crystal
    REAL (RealK), Pointer :: d2_volume(:)
!     2nd derivative of projected area of crystal
!
!
  END TYPE STR_db_cryst_geom
!
!
!
END MODULE def_db_crystal_geometry
