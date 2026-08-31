! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate specific humidities from Gill's formula.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE qsat_gill(sat_spec_hum, t, p, n_profile)
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     This routine computes the saturated specific humidity
!     at temperature T and pressure p, using the formulae given in 
!     Appendix four of Adrian Gill's book.
!
!     Note that the formulae work with pressures in hectopascals and
!     temperatures in degrees celsius. the conversions are made
!     inside this routine and should have no impact on the rest
!     of the code.
!
!     This routine was perpetrated by D L Roberts (12/8/93).
!
!     Modified to cope with very low pressures by DLR (27/10/93).
!
      INTEGER, Intent(IN) :: n_profile
!           Number of profiles.
!
      REAL  (RealK), Intent(IN) :: p(n_profile), t(n_profile)
!           Pressure in pascals, Temperature in kelvin
!
      REAL  (RealK), Intent(OUT) :: sat_spec_hum(n_profile)
!           Saturated specific humidity at T and p.
!
!
!     Local variables.
!
      INTEGER :: l
!
      REAL  (RealK) :: press ! pressure in hectopascals.
      REAL  (RealK) :: temp ! temperature in celsius.
      REAL  (RealK) :: a ! a temporary holding variable.
      REAL  (RealK) :: ew ! saturation vapour pressure of 
!                         ! PURE WATER VAPOUR.
      REAL  (RealK) :: ewdash ! sat vap pressure of water vapour in air.
      REAL  (RealK) :: fw ! the ratio between ewdash and ew.
      REAL  (RealK) :: zero_degc ! kelvin equivalent of zero celsius.
!
!     The value assigned is that used in v3.1 of the unified model.
!
      parameter( zero_degc=273.15e+00_RealK )
!
      REAL  (RealK) :: epsilon ! the ratio of the 
!                              ! MOLECULAR MASS OF WATER
!                            to that of dry air.
!     The value assigned is that used in v3.1 of the unified model.
!
      parameter( epsilon=0.62198e+00_RealK )
!
      REAL  (RealK) :: eta ! one minus epsilon.
      parameter( eta = 1.0e+00_RealK - epsilon )
!
!
!     Loop over all points.
!
      DO l = 1,n_profile
!
!       Convert to local units for temperature and pressure.
!       
        temp = t(l) - zero_degc
        press = p(l)*0.01e+00_RealK
!       
!       Equation (A4.7) of Gill's book.
!       
        fw = 1.0e+00_RealK + 1.0e-06_realk*press*( 4.5e+00_realk &
          + 6.0e-04_RealK*temp*temp )
!       
!       Equation (A4.5) of Gill.
!       
        a = ( 7.859e-01_RealK + 3.477e-02_realk*temp ) &
           /( 1.0e+00_RealK + 4.12e-03_realk*temp )
        ew = 1.0e+01_RealK**a
!       
!       Equation (A4.6) of Gill.
!       
        ewdash = fw*ew
!       
!       The next equation is a rearrangement of Gill's (A4.3),
!       with w subscripts added because we require saturation.
!       
!       Note that at very low pressures a fix has to be applied,
!       to avoid a singularity.
!       
        IF (press  >  ewdash) THEN
          sat_spec_hum(l) = (epsilon*ewdash)/(press-ewdash*eta)
        ELSE
          sat_spec_hum(l) = 1.0e+00_RealK
        ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE qsat_gill
