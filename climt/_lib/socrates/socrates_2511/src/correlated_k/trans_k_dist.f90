! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate transmission over a k-distribution.
FUNCTION trans_k_dist &
!
(n_nu, k, nu_inc, wgt, integ_wgt, n_path, u) &
!
RESULT (trans)
!
! Description:
!   This routine returns an array of tramissions calculated
!   over a k-distribution for a range of amounts of absorber.
!
! Method:
!   Straightforward.
!
!
!
! Modules used.
  USE realtype_rd
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  INTEGER, Intent(IN) :: n_nu
!   Number of frequencies
  INTEGER, Intent(IN) :: n_path
!   Number of paths
  REAL  (RealK), Intent(IN), Dimension(n_nu) :: k
!   Monochromatic absorption coefficients
  REAL  (RealK), Intent(IN) :: nu_inc
!   Frequency increment
  REAL  (RealK), Intent(IN), Dimension(n_nu) :: wgt
!   Weighting function
  REAL  (RealK), Intent(IN) :: integ_wgt
!   Integral of the weighting function
  REAL  (RealK), Intent(IN), Dimension(n_path) :: u
!   Pathlengths of absorber
!
  REAL  (RealK), Dimension(n_path) :: trans
!   Calculated transmissions
!
! Local variables.
  INTEGER :: i
!   Loop variable
  REAL  (RealK), Dimension(n_nu) :: wgt_trans
!
!
!- End of Header
!
!
!
  DO i=1, n_path
    wgt_trans(1:n_nu)=wgt(1:n_nu)*EXP(-k(1:n_nu)*u(i))
    trans(i)=nu_inc * SUM(wgt_trans(1:n_nu))/integ_wgt
  ENDDO
!
!
!
  RETURN
END FUNCTION trans_k_dist
