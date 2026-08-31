! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set permitted types of scattering algorithm.
!
MODULE scatter_algorithm_pcf
!
! Description:
!
! This module defines the recognised algorithms used to 
! calculate electromagnetic scattering.
!
!- End of header
!
!
  INTEGER, Parameter :: npd_algorithm_scatter = 5
!           Number of algorithms
!
  INTEGER, Parameter :: ip_algorithm_mie = 1
!           Identification number for full Mie calculation
  INTEGER, Parameter :: ip_algorithm_adt = 2
!           Identification number for ADT
  INTEGER, Parameter :: ip_adtmitchell96 = 4
!           Identification number for Mitchell's enhanced ADT
  INTEGER, Parameter :: ip_scat_database = 5
!           Identification number for calculation using an
!           external database
!
  CHARACTER (LEN=40), Parameter, Dimension(npd_algorithm_scatter) :: &
    name_scatter_algorithm = (/ &
    'Full Mie scattering calculation         ', &
    'Anomalous diffraction approximation     ', &
    'Obsolete scheme: removed                ', &
    'Modified ADT of Mitchell (1996)         ', &
    'External database                       ' /)
!
END MODULE scatter_algorithm_pcf
