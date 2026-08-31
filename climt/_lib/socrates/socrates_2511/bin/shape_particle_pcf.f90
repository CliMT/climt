! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+    Module to set the shape of particles.
!
MODULE shape_particle_pcf
!
! Description:
!
! This module defines the recognised shapes of ice crystals.
!
!- End of header
!
!
      INTEGER, Parameter :: npd_shape                = 11
!       Number of permitted shapes
!
      INTEGER, Parameter :: ip_shape_null            =  0
!       Unassigned shape
!
      INTEGER, Parameter :: ip_shape_sphere          =  1
!       Spherical particles
!
!     Obsolete hexagonal cylinders
      INTEGER, Parameter :: ip_shape_hexcyl          =  2
!       Hexagonal cylinders
!
!     Shapes used in David Mitchell's schemes
      INTEGER, Parameter :: ip_shape_polycrystal     =  3
!       Polycrystals
      INTEGER, Parameter :: ip_shape_plate           =  4
!       Plates
      INTEGER, Parameter :: ip_shape_rosette         =  5
!       Rosettes
      INTEGER, Parameter :: ip_shape_column          =  6
!       Hexagonal columns
!
!     Shapes available in Baran's databases
      INTEGER, Parameter :: ip_shape_db_aggregate    =  7
!       Baran's aggregate particle
      INTEGER, Parameter :: ip_shape_db_undefined    =  8
!       Undefined shape for use with databases:
!       this allows the use of observational size 
!       distributions with various shapes
      INTEGER, Parameter :: ip_shape_db_column       =  9
!       Baran's aggregate particle
      INTEGER, Parameter :: ip_shape_db_aggregate_py = 10
!       Baran's aggregate particle
      INTEGER, Parameter :: ip_shape_db_dendrite     = 11
!       Baran's aggregate particle
!
END MODULE shape_particle_pcf
