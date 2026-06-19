! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set the possible types of correlated-k fits.
!
MODULE ck_fit_pcf
!
! Description:
!   This module defines the ways in which the correlated-k fits
!   may be derived.
!
!- End of header

  IMPLICIT NONE


  INTEGER, Parameter  ::  IP_ck_none = 0
!   No calculation of k-terms
  INTEGER, Parameter  ::  IP_ck_tol = 1
!   The correlated-k fits are done to a specified tolerance
!   for each term                        
  INTEGER, Parameter  ::  IP_ck_fixed_n = 2
!   A fixed specified number of terms is used for each fit
  INTEGER, Parameter  ::  IP_ck_bin = 3
!   Absorption coefficients are split into a number of bins
!   depending on their scaling behaviour. The number of terms
!   in each bin depends on a specified tolerance.

END MODULE ck_fit_pcf
