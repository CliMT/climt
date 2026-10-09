! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate mie scattering cross-sections and asymmetry.
!
! Method:
!   This routine is closely based on that given in "Absorption 
!   and Scattering of Light by Small Particles" by C. F. Bohren 
!   and Huffman. It uses downward recurrence.
!
!- ---------------------------------------------------------------------
      SUBROUTINE mie_scatter(ierr, x, n_relative
     &  , q_scatter, q_ext, asymmetry
     &  , l_phase, n_angle, mu_angle, s1, s2
     &  , nd_scatt_angle, nd_log_deriv
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, Intent(IN) ::
     &    nd_scatt_angle
!           Size allocated for angles where the phase function is
     &  , nd_log_deriv
!           Size allocated for arrays of logarithmic derivatives
!           calculated
!
!
!     Dummy arguments.
      complex  (RealK), Intent(IN) ::
     &    n_relative
!           Relative refractive `index'
      complex  (RealK), Intent(OUT) ::
     &    s1(nd_scatt_angle)
!           Element of scattering matrix
     &  , s2(nd_scatt_angle)
!           Element of scattering matrix
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      INTEGER, Intent(IN) ::
     &    n_angle
!           Number of scattering angles
      LOGICAL, Intent(IN) ::
     &    l_phase
!           Calculate phase function
      REAL  (RealK), Intent(IN) ::
     &    x
!           Size parameter
     &  , mu_angle(nd_scatt_angle)
!           Cosines of scattering angles
      REAL  (RealK), Intent(OUT) ::
     &    q_scatter
!           Scattering efficiency
     &  , q_ext
!           Extinction efficiency
     &  , asymmetry
!           Asymmetry factor
!
!     Local arguments.
      complex  (RealK) ::
     &    d(nd_log_deriv)
!           Logarithmic derivatives
     &  , c
!           Complex number of term
     &  , y
!           Size parameter*refractive `index'
     &  , ksi_n
!           Riccati-bessel function at n
     &  , ksi_nm1
!           Riccati-bessel function at n-1
     &  , ksi_nm2
!           Riccati-bessel function at n-2
     &  , a
!           Coefficient in series
     &  , b
!           Coefficient in series
     &  , a_old
!           Previous coefficient
     &  , b_old
!           Previous coefficient
      INTEGER
     &    n_term
!           Number of terms computed
     &  , n_down
!           Start of downward recurrence
     &  , n_y
!           Integral part of y
     &  , i_term
!           Index of term
     &  , i_angle
!           Anglular `index'
      REAL  (RealK) ::
     &    r
!           Numerical factor
     &  , rx
!           Numerical factor
     &  , factor
!           Numerical factor
     &  , x2
!           Numerical factor
     &  , psi_n
!           Riccati-bessel function at n
     &  , psi_nm1
!           Riccati-bessel function at n-1
     &  , psi_nm2
!           Riccati-bessel function at n-2
     &  , chi_n
!           Riccati-bessel function at n
     &  , chi_nm1
!           Riccati-bessel function at n-1
     &  , chi_nm2
!           Riccati-bessel function at n-2
     &  , pi_n(nd_scatt_angle)
!           Riccati-bessel function at n
     &  , pi_nm1(nd_scatt_angle)
!           Riccati-bessel function at n-1
     &  , pi_nm2(nd_scatt_angle)
!           Riccati-bessel function at n-2
     &  , tau_n(nd_scatt_angle)
!           Riccati-bessel function at n
!
!
      y=x*n_relative
!
!     Set the number of terms required in the scattering series.
      n_term=int(2.0_RealK+4.0_RealK*x**3.333e-01_RealK+x)
!     Set the term at which downward recurrence begins to be
!     large enough.
      n_y=1.1_RealK*abs(y)
      n_down=max(n_term, n_y)+15
!
!     Check array bounds.
      IF (n_down > nd_log_deriv) THEN
        WRITE(iu_err, '(/a)') '*** error: too many logarithmic '
     &    //'derivatives are required.'
        WRITE(iu_err, '(a, 1pe12.5)') 'size parameter = ', x
        WRITE(iu_err, '(a, i8)') 'recuirrence starts at ', n_down
        WRITE(iu_err, '(a)')
     &    'increase nd_log_deriv and re-compile.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Calculate the logarithmic derivatives, D, by downward recurrence.
!     D(N_DOWN) is initialized to 0.
      d(n_down)=CMPLX(0.0_RealK, 0.0_RealK)
      DO i_term=n_down-1, 1, -1
         c=(REAL(i_term, RealK)+1.0_RealK)/y
         d(i_term)=c-1.0_RealK/(d(i_term+1)+c)
      ENDDO
!
!     Main loop of series follows: Riccati-Bessel functions 
!     are calculated by upward recurrence.
!
!     Initializtion: this starts with at -1, i.e. n-2 for n=1 below.
      psi_nm2=cos(x)
      psi_nm1=sin(x)
      chi_nm2=-sin(x)
      chi_nm1=cos(x)
      ksi_nm2=CMPLX(psi_nm2, -chi_nm2)
      ksi_nm1=CMPLX(psi_nm1, -chi_nm1)
      q_scatter=0.0_RealK
      q_ext=0.0_RealK
      asymmetry=0.0_RealK
!     Initialization of angularly-dependent terms.
      IF (l_phase) THEN
        DO i_angle=1, n_angle
          pi_nm2(i_angle)=0.0_RealK
          pi_nm1(i_angle)=0.0_RealK
          s1(i_angle)=(0.0_RealK, 0.0_RealK)
          s2(i_angle)=(0.0_RealK, 0.0_RealK)
        ENDDO
      ENDIF
!
      i_term=0
!
100     i_term=i_term+1
        r=REAL(i_term, RealK)
        rx=r/x
        factor=(2.0_RealK*r+1.0_RealK)/(r*(r+1.0_RealK))
!       Compute new psi and ksi by recurrence.
        psi_n=(2.0_RealK*r-1.0_RealK)*psi_nm1/x-psi_nm2
        chi_n=(2.0_RealK*r-1.0_RealK)*chi_nm1/x-chi_nm2
        ksi_n=CMPLX(psi_n, -chi_n)
        a=((d(i_term)/n_relative+rx)*psi_n-psi_nm1)
     &    /((d(i_term)/n_relative+rx)*ksi_n-ksi_nm1)
        b=((n_relative*d(i_term)+rx)*psi_n-psi_nm1)
     &    /((n_relative*d(i_term)+rx)*ksi_n-ksi_nm1)
!
        q_scatter=q_scatter+(2.0_RealK*r+1.0_RealK)
     &    *((abs(a))**2+(abs(b))**2)
        q_ext=q_ext+(2.0_RealK*r+1.0_RealK)*REAL(a+b, RealK)
!
        IF (i_term > 1) THEN
          asymmetry=asymmetry+(r*r-1.0_RealK)
     &      *REAL(a_old*conjg(a)+b_old*conjg(b), RealK)
     &      /r+(2.0_RealK*r-1.0_RealK)
     &      *REAL(a_old*conjg(b_old), RealK)
     &      /(r*(r-1.0_RealK))
        ENDIF
!
        IF (l_phase) THEN
!         Functions pi and tau are required for elements of scattering
!         matrix.
          DO i_angle=1, n_angle
            IF (i_term == 1) THEN
              pi_n(i_angle)=1.0_RealK
            ELSE
              pi_n(i_angle)=((2.0_RealK*r-1.0_RealK)
     &          /(r-1.0_RealK))*mu_angle(i_angle)
     &          *pi_nm1(i_angle)-(r/(r-1.0_RealK))*pi_nm2(i_angle)
            ENDIF
            tau_n(i_angle)=r*mu_angle(i_angle)*pi_n(i_angle)
     &        -(r+1.0_RealK)*pi_nm1(i_angle)
            s1(i_angle)=s1(i_angle)+factor*(a*pi_n(i_angle)
     &        +b*tau_n(i_angle))
            s2(i_angle)=s2(i_angle)+factor*(a*tau_n(i_angle)
     &        +b*pi_n(i_angle))
          ENDDO
        ENDIF
!
!       Test for next term of series and loop around if required.
        IF (i_term < n_term) THEN
!         Advance stored terms for recurrence.
          a_old=a
          b_old=b
          psi_nm2=psi_nm1
          psi_nm1=psi_n
          chi_nm2=chi_nm1
          chi_nm1=chi_n
          ksi_nm1=ksi_n
          IF (l_phase) THEN
            DO i_angle=1, n_angle
              pi_nm2(i_angle)=pi_nm1(i_angle)
              pi_nm1(i_angle)=pi_n(i_angle)
            ENDDO
          ENDIF
          goto 100
        ENDIF
!
      x2=1.0_RealK/x**2
      q_scatter=2.0_RealK*q_scatter*x2
      q_ext=2.0_RealK*q_ext*x2
      asymmetry=4.0_RealK*asymmetry*x2/q_scatter
!
!
!
      RETURN
      END
