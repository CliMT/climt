! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+    Module to set identification of distribution types.
!
MODULE distribution_pcf
!
! Description:
!   This module defines the available size distributions.
!
!- End of header
!
!
  IMPLICIT NONE
!
!
!
  INTEGER, Parameter  ::  npd_distribution = 7
!                           Number of permitted distributions
!
!
  INTEGER, Parameter  ::  ip_external = 1
!                           External distribution with dN/dr
!                           specified explicitly
  INTEGER, Parameter  ::  ip_log_normal = 2
!                           Log normal distribution
  INTEGER, Parameter  ::  ip_modified_gamma = 3
!                           Modified gamma distribution
  INTEGER, Parameter  ::  ip_heymsfield_platt = 4
!                           Heymsfield and Platt's distribution
  INTEGER, Parameter  ::  ip_mitchell_97 = 5
!                           Mitchell's ice distribution of 1997
!                           (See D. L. Mitchell, J. M. Edwards and
!                           P. N. Francis (1997) "GCM Sensitivity
!                           of globally averaged albedo and OLR to
!                           ice crystal shape", 7th ARM Science Team
!                           Meeting, San Antonio, TX.)
  INTEGER, Parameter  ::  ip_mitchell_trop_00 = 6
!                           Mitchell's ice distribution of 2000
!                           (See D. L. Mitchell, D. Ivanova, A. Macke
!                           and G. M. McFarquhar (2000) "A GCM 
!                           parameterization of bimodal spectra for
!                           ice clouds", 9th ARM Science Team
!                           Meeting.)
  INTEGER, Parameter  ::  ip_ivanova_mlat_00 = 7
!                           Mitchell's ice distribution of 2000
!                           (See D. Ivanova, D. L. Mitchell, W. P. Annott,
!                           and W. Poellot (2000) "A GCM 
!                           parameterization of bimodal spectra for
!                           mid-latitude ice clouds", 13th Int. Conf. on
!                           Clouds and Precipitation.)
!
  CHARACTER (LEN=30), Parameter :: &
                              name_distribution(npd_distribution) = (/ &
                                     "External                      ", &
                                     "Log-normal                    ", &
                                     "Modified Gamma                ", &
                                     "Heymsfield and Platt          ", &
                                     "Mitchell (1996)               ", &
                                     "Mitchell (1998) Tropical      ", &
                                     "Mitchell (2000) Mid-latitude  "  &
                                    /)
!
END MODULE distribution_pcf
