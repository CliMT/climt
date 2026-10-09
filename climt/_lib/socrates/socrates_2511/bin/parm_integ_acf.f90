! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set parameters for infinite integrals.
!
MODULE parm_integ_acf
!
! Description:
!
! This defines the parameters used to carry out the integration
! over size distributions.
!
!- End of header
!
!
! Modules used:
   USE realtype_rd
!
!
!!    npd_refinement, npd_panel_point and npd_point are not
!!    really independent. ideally, npd_point=1+(npd_panel_point-1)
!!    *(2**(npd_refinement-1)). Other setting will not make full use of
!!    the declared space.
!
  INTEGER, Parameter :: npd_panel_point = 5
!           Number of points in panel
  INTEGER, Parameter :: npd_integral    = 20
!           Number of integrals available
  INTEGER, Parameter :: npd_point       = 16387
!           Max. no. of points in interval
  INTEGER, Parameter :: npd_panel       = 20
!           Maximum number of panels
  INTEGER, Parameter :: npd_refinement  = 12
!           Maximum number of refinements
!
  REAL  (RealK), Parameter :: p_panel_ratio = 3.2_RealK
!           Size of each panel
!
END MODULE parm_integ_acf
