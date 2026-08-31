! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the g quadrature points for correlated-k.
!
SUBROUTINE set_g_point_90 &
!
(n_nu, nu_inc, kabs, wgt, integ_wgt, &
 i_ck_fit, tol, max_path, &
 l_kabs_wgt, kabs_wgt, l_wgt_scale_sqrt, u_wgt_scale, &
 nd_k_term, iu_monitor, &
 n_k, w_k, k_opt, k_ave, &
 ig, ierr)
!
! Description:
!   This routine is supplied with the table of monochromatic
!   absorption coefficients, determines an appropriate number 
!   of terms for the fit, if required, and returns the appropriate
!   k-values and weights.
!
!
! Modules used:
  USE realtype_rd
  USE error_pcf
  USE def_std_io_icf
  USE ck_fit_pcf
  USE ck_parm_acf
!
!
  IMPLICIT NONE
!
!
!
! Dummy variables.
  INTEGER, Intent(IN) :: n_nu
!   Number of frequencies covering the band
  INTEGER, Intent(IN) :: iu_monitor
!   Unit number for monitoring output
  INTEGER, Intent(IN) :: nd_k_term
!   Size allocated for k-terms
  INTEGER, Intent(IN) :: i_ck_fit
!   Type of correlated-k fit required
  INTEGER, Intent(INOUT) :: ig(0:nd_k_term)
!   Indices of the last quadrature point lying in the interval
!   associated with the current term
  INTEGER, Intent(INOUT) :: n_k
!   Number of terms in the fit
  INTEGER, Intent(INOUT) :: ierr
!   Error flag
!
  REAL  (RealK), Intent(IN) :: nu_inc
!   Increment of wavenumber across through the monochromatic 
!   absorption ocefficients
  REAL  (RealK), Intent(IN) :: kabs(:)
!   Monochromatic absorption coefficients
  REAL  (RealK), Intent(IN) :: wgt(:)
!   Monochromatic weightings
  REAL  (RealK), INTENT(IN) :: integ_wgt
!   Integral of the weighting across the band
  REAL  (RealK), Intent(IN) :: tol
!   Tolerance for the c-k fit
  REAL  (RealK), Intent(IN) :: max_path
!   Maximum pathlength to be considered for the absorber
  REAL  (RealK), Intent(IN) :: kabs_wgt(:)
!   Monochromatic absorption coefficients to use in weighting
  REAL  (RealK), Intent(IN) :: u_wgt_scale
!   Factor to scale continuum column mass to gas column mass
  REAL  (RealK), Intent(OUT) :: w_k(:)
!   Weights for the k-terms
  REAL  (RealK), Intent(OUT) :: k_opt(:)
!   Optimal absorption coefficients for each k-term
  REAL  (RealK), Intent(OUT) :: k_ave(:)
!   Mean absorption coefficients for each k-term
!
  LOGICAL, Intent(IN) :: l_kabs_wgt
!   Use k_wgt as additional weights in transmissions
  LOGICAL, Intent(IN) :: l_wgt_scale_sqrt
!   If true gas column mass is scaled using the square root of the continuum
!   column mass. Otherwise a linear scaling is used. 
!   
!
! Local Variables:
  INTEGER :: i, j
!   Loop variable
  INTEGER :: il
!   Index of monochromatic absorption point at the left of the
!   chopping interval
  INTEGER :: ir
!   Index of monochromatic absorption point at the right of the
!   chopping interval
  INTEGER :: im
!   Index of monochromatic absorption point in the middle of the
!   chopping interval
  INTEGER :: n_nu_k
!   Number of frequencies in the interval of this k-term
  INTEGER :: n_div
!   Number of divisions in the log of absorption
  INTEGER:: i0
!   Pointer for clarity of coding
  INTEGER :: i1
!   Pointer for clarity of coding

  REAL  (RealK) :: k_max_grey
!   Maximum value of absorption which can be considered grey
  REAL  (RealK) :: dlogk
!   Interval of log(k) across the whole band
  REAL  (RealK) :: dlogk_inc
!   Interval of log(k) across the interval of the current term
  REAL  (RealK) :: integ_k
!   Integral of the weighting function across the region represented
!   by any k-term
  REAL  (RealK) :: ktop(nd_k_term)
!   Maximum absorption in any interval
  REAL  (RealK), Dimension(n_nu) :: wgt_k
!   Product of the absorption coefficient and the weighting function
  REAL  (RealK) :: km
!   Absorption coefficient in the middle of the chopping interval
  REAL  (RealK) :: err_norm_k
!   Error norm for each individual k-term
  REAL  (RealK) :: err_norm
!   Error norm for fitting
!   
! Functions called:
!
  INTERFACE
!
    SUBROUTINE optimal_k &
      (n_nu, nu_inc, k, wgt, integ_k, k_mean, tol, &
       l_k_wgt, k_wgt, l_wgt_scale_sqrt, u_wgt_scale, &
       k_opt, error, ierr)
!
      USE realtype_rd
!
      INTEGER, Intent(IN) :: n_nu
      INTEGER, Intent(INOUT) :: ierr
      REAL  (RealK), Intent(IN) :: nu_inc
      REAL  (RealK), Intent(IN) :: k(:)
      REAL  (RealK), Intent(IN) :: wgt(:)
      REAL  (RealK), Intent(IN) :: integ_k
      REAL  (RealK), Intent(IN) :: k_mean
      REAL  (RealK), Intent(IN) :: tol
      REAL  (RealK), Intent(IN) :: k_wgt(:)
      REAL  (RealK), Intent(IN) :: u_wgt_scale
      REAL  (RealK), Intent(OUT) :: k_opt
      REAL  (RealK), Intent(OUT) :: error
      LOGICAL, Intent(IN) :: l_k_wgt
      LOGICAL, Intent(IN) :: l_wgt_scale_sqrt
!
    END SUBROUTINE
!
  END INTERFACE
!
!   
!
! Care is needed lest kabs should have entries equal to 0. Ideally, we
! would split into equal increments in log(k), but if kabs is very small
! or 0 this will generate many small terms, so we actually need to find
! when k is small enough to be taken as grey.
! 
! Initially, determine the limit at which k can be taken as grey.
! A power series expansion of the exponential suggests that the 
! error should be proportional to the inverse second power of the
! number of terms.
  SELECT CASE(i_ck_fit)
    CASE(IP_ck_tol,IP_ck_bin)
      k_max_grey = 0.5_RealK * SQRT(tol) / max_path
    CASE(IP_ck_fixed_n)
!     Experience shows that about 25 terms give an error of about
!     10^-6, so we are a little more stringent.
      k_max_grey = 0.01_RealK / (max_path * REAL(n_k, RealK) )
  END SELECT
!
  IF (k_max_grey < kabs(1)) THEN
!   A straight division in dlogk is desirable.
    dlogk = LOG(kabs(n_nu) / kabs(1))
  ELSE
    dlogk = LOG(kabs(n_nu) / k_max_grey)
  ENDIF
  IF (dlogk < EPSILON(dlogk) ) THEN
!   Only one term is needed.
    ig(0)=0
    ig(1)=n_nu
    w_k(1)=1.0_RealK
    k_ave(1)=kabs(n_nu)
    k_opt(1)=kabs(n_nu)
  ENDIF
!   
!
  IF (i_ck_fit == IP_ck_tol .OR. i_ck_fit == IP_ck_bin) THEN
!
    WRITE(iu_monitor, '(/a/)') &
      "Trial fits with increasing numbers of terms."
!   Start with a fit of one term.
    n_k=1
  ENDIF

  WRITE(iu_monitor, '(a)') &
    '-------------------------------------------------------'
  n_div = n_k
  DO
    dlogk_inc=dlogk/REAL(n_div, RealK)
!   Split into equal increments in dlog(k).
    ig(0)=0
    err_norm=0.0_RealK
    WRITE(iu_monitor, '(a)') &
      '-------------------------------------------------------'
!   The number of k-terms may be less than the number of intervals in the
!   case where some intervals contain no absorption coefficients.
!   We reset the number of terms to 1 and increment where absorption
!   coefficients are present in the interval.
    n_k = 1
    DO j=1, n_div
!     ktop marks the upper boundary of the interval: all points
!     with absorption greater than ktop fall in the next interval.
      ktop(n_k)=kabs(n_nu)/EXP(REAL(n_div-j, RealK)*dlogk_inc)

!     Determine the index of the last quadrature point which falls
!     in this interval by using a binary chop.
!     Set the left and right-hand ends of the interval
      il=ig(n_k-1)
      ir=n_nu
      DO
        IF (j == n_div) THEN
          il=n_nu
          EXIT
        ENDIF
        IF (il >= ir-1) EXIT
        im=(il+ir)/2
        km=kabs(im)
        IF (km > ktop(n_k)) THEN
          ir=im
        ELSE IF (km <= ktop(n_k)) THEN
          il=im
        ENDIF
      ENDDO
!     Check whether there are any coefficients in the interval.
      IF (il > ig(n_k-1)) THEN
        ig(n_k)=il
        IF (j < n_div) n_k=n_k+1
      END IF
    END DO
    DO i=1, n_k
      i0 = ig(i-1)+1
      i1 = ig(i)
      n_nu_k=i1-i0+1
!     Calculate the simple mean k-value across this interval.
      integ_k=nu_inc * SUM(wgt(i0:i1))
      wgt_k(i0:i1)=kabs(i0:i1)*wgt(i0:i1)
      k_ave(i)=nu_inc * SUM(wgt_k(i0:i1))/integ_k
      w_k(i)=integ_k/integ_wgt
      CALL optimal_k(n_nu_k, nu_inc, kabs(i0:i1), wgt(i0:i1), &
         integ_k, k_ave(i), tol, &
         l_kabs_wgt, kabs_wgt(i0:i1), l_wgt_scale_sqrt, u_wgt_scale, &
         k_opt(i), err_norm_k, ierr)
!     The error norm is constructed by weighting the squared error 
!     for each k-term by the weight applied to the term. It's not
!     clear that this is necessarily the best measure for applications.
      err_norm=err_norm+w_k(i)*err_norm_k*err_norm_k
      WRITE(iu_monitor,'(a,i4,f10.6,e14.6)') ' k-term, g-interval, error:', &
        i,REAL(i1)/REAL(n_nu), w_k(i)*err_norm_k*err_norm_k
    END DO

    err_norm=SQRT(err_norm)
    WRITE(iu_monitor, '(a, i5, a6, a, 1pe12.5)') &
      'Trial fit contains ', n_k, ' terms:', &
      ' Error = ', err_norm
!
    IF (i_ck_fit == IP_ck_fixed_n) THEN
!     We have a fit with the prescribed number of terms and so
!     pass on to the next stage.
        WRITE(iu_monitor, "(a)") 'Fit accepted.'
      EXIT
    ELSE
      IF (err_norm < tol) THEN
        WRITE(iu_monitor, "(a)") 'Fit to tolerance has been obtained.'
        EXIT
      ELSE IF (n_k > nd_k_term-1) THEN
        WRITE(iu_monitor, "(a)") 'Fit to tolerance has not been obtained.'
        WRITE(iu_err, '(/A)') &
          '*** Warning: Too many terms are required for a ' &
          //'fit to the requested accuracy'
        EXIT
      ELSE
        n_div=n_div+1
      ENDIF
    ENDIF
!
  ENDDO

END SUBROUTINE set_g_point_90
