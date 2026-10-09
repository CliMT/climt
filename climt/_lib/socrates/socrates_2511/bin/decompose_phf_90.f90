! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to decompose a phase function into moments.
!
Subroutine decompose_phf_90 &
!
(n_weight, mu_angle, weight, &
 n_phf_term, i_stokes, phase_fnc_term, &
 ierr &
)
!
! Description:
!   This subroutine receives a distribution of intensities (the first
!   Stokes parmeter) specified at discrete angles on a Gaussian grid
!   and calculates the moments of the phase function from them.
!
! Method:
!   Each moment of the phase function is the an inner product with
!   the appropriate Legendre polynomial. The Legendre polynomials
!   are evaluated by upward recurrence, while the integrals are 
!   evaluated using Gaussian quadrature with precalculated weights.
!
!
!
! Modules used
  USE realtype_rd
  USE error_pcf
  USE def_std_io_icf
!
!
  IMPLICIT NONE
!
!
!
! Dummy varaiables:
  INTEGER, Intent(IN) :: n_weight
!           Number of weighting points
  INTEGER, Intent(IN) :: n_phf_term
!           Order of truncation 
!
  REAL  (RealK), Intent(IN) :: mu_angle(n_weight)
!           Polar angles used in quadrature
  REAL  (RealK), Intent(IN) :: weight(n_weight)
!           Weights at quadrature points
  REAL  (RealK), Intent(IN) :: i_stokes(n_weight)
!           Stokes's first parameter at the weighting points
!
!
  REAL  (RealK), Intent(OUT) :: phase_fnc_term(n_phf_term)
!           Moments of the phase function
!
  INTEGER, Intent(INOUT) :: ierr
!           Error flag
!
!
!
! Local variables
  INTEGER :: order
!           Order of term
  INTEGER :: i
!           Loop variable
  REAL  (RealK) :: lp_n(n_weight)
!           Legendre polynomials of order N at the weighting points
  REAL  (RealK) :: lp_nm1(n_weight)
!           Legendre polynomials of order N-1 at the weighting points
  REAL  (RealK) :: lp_nm2
!           Legendre polynomial of order N-2 at the weighting point
!           (there is no need to store values at all points for this)
  REAL  (RealK) :: ls
!           Real variable holding the current order
  REAL  (RealK) :: normalization
!           Normalization for the phase function
!
!    
!
! Check that the intended order of quadrature will be sufficient
! to yield the moments reliably.
  IF (n_phf_term > n_weight) THEN
    WRITE(iu_err, '(/A)') &
      '*** Error: The order of Gaussian quadrature is too low.'
    ierr=i_err_fatal
    RETURN
  ENDIF
!
!
! Initialize the unsed zeroth order to start the recurrence.
  normalization = 0.0_RealK
  DO i = 1, n_weight
    lp_n(i) = 1.0_RealK
    normalization = normalization + weight(I) * i_stokes(i)
  ENDDO
  order=1
!
  DO ; IF (order > n_phf_term) EXIT
!
    phase_fnc_term(order) = 0.0_RealK
    DO i=1, n_weight
!
!     Use a recurrence relation to determine the Legendre
!     polynomials.
      IF (order >= 2) THEN
        lp_nm2    = lp_nm1(i)
        lp_nm1(i) = lp_n(i)
        ls        = REAL(order, RealK)
        lp_n(i)   = ((2.0_RealK * ls - 1.0_RealK) * &
          mu_angle(i) * lp_nm1(i) - &
          (ls-1.0_RealK) * lp_nm2) / ls
      ELSE IF (order == 1) THEN
        lp_nm1(i) = lp_n(i)
        lp_n(i)   = mu_angle(i)
      ENDIF
!
      phase_fnc_term(order)=phase_fnc_term(order) + &
        weight(i) * lp_n(i) * i_stokes(i)
!
    ENDDO
!
    phase_fnc_term(order)=phase_fnc_term(order)/normalization
    order=order+1
!
  ENDDO
!
!
!
  RETURN
END
