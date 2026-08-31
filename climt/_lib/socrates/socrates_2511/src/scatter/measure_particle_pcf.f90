! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+    Module to set measures of particle sizes.
!
MODULE measure_particle_pcf
!
! Description:
!
! This defines the ways in which the sizes of ice crystals may
! be defined.
!
!- End of header
!
!
      INTEGER, Parameter :: ip_measure_radius = 1
!           Radius of the particle (principally for spheres)
      INTEGER, Parameter :: ip_measure_max_dimen = 2
!           Maximum dimension of the particle
      INTEGER, Parameter :: ip_measure_proj_radius = 3
!           Radius of sphere of equivalent projected area
!
END MODULE measure_particle_pcf
