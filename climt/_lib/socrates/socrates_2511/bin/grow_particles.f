! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to treat the growth of particles with humidity.
!
! Method:
!
!   THIS ROUTINE CHANGES THE LOGARITHMS OF THE MEAN AND STANDARD
!   DEVIATION OF A LOG-NORMAL DISTRIBUTION OF PARTICLE SIZES, TO
!   ACCOUNT FOR GROWTH DUE TO HUMIDITY.
!
!   AUTHOR: D.L.ROBERTS.
!
!   THIS IS VERSION 6 OF THE ROUTINE. (OCT 15TH 1998)
!   VERSION 1 DID NOT ALTER THE STANDARD DEVIATION, I.E. ALL
!   PARTICLES WERE ASSUMED TO GROW AT THE SAME RATE. THIS
!   ASSUMPTION WAS RELAXED AT VERSION 2. HOWEVER, THIS
!   INTRODUCED A SUBTLE BUG: WHEN FITZGERALD'S BETA PARAMETER
!   IS NOT UNITY, ONE HAS TO TAKE INTO ACCOUNT THE FACT THAT
!   HIS FORMULA FOR ALPHA WAS STATED IN A FORM DERIVED FOR
!   PARTICLE RADII IN MICRONS. SINCE THIS PROGRAM WORKS IN
!   METRES, ALPHA HAS TO BE INCREASED BY A FACTOR WHICH TURNS
!   OUT TO BE 10.0**( 6.0*(BETA-1) ). VERSION 3 DIFFERED FROM
!   VERSION 2 BY INCLUDING THIS FACTOR.
!   VERSION 3, HOWEVER, PRODUCED A DISCONTINUITY IN THE SIZE
!   OF THE PARTICLES AT RELATIVE HUMIDITY 0.81. VERSION 4
!   CORRECTS THIS BY INTERPOLATING BETA BETWEEN HUMIDITIES
!   0.3 AND 0.81, IN THE SAME WAY AS ALPHA HAS ALWAYS BEEN
!   INTERPOLATED. VERSION 5 IS ALMOST THE SAME AS VERSION 4,
!   BUT WITH THE UNIT CONVERSION CORRECTION ONLY APPLIED AT
!   HUMIDITIES GREATER THAN OR EQUAL TO 0.81, AS IN VERSION 3.
!   VERSION 6 DIFFERS FROM VERSION 5 IN THAT SODIUM CHLORIDE
!   IS NOW HANDLED. THE TREATMENT OF AMMONIUM SULPHATE IS
!   UNALTERED FROM VERSION 5.
!
!   ABOVE THE DELIQUESCENCE POINT, TAKEN AS 0.81,
!   THE SCHEME IS THE ONE DUE TO J.W.FITZGERALD (1975).
!   IT IS NOT VERY GOOD BUT WILL DO FOR THE TIME BEING.
!   THE HYSTERESIS EFFECT IS MODELLED BY A REGION WHERE GROWTH
!   IS LINEAR. THE JUSTIFICATION FOR THIS IS THAT NEAR THE
!   DELIQUESCENCE POINT, THE PROBABILITY THAT A PARTICLE HAS
!   PREVIOUSLY BEEN EXPOSED TO A HUMIDITY IN EXCESS OF THAT
!   NEEDED FOR DELIQUESCENCE IS HIGH, SO THAT ON AVERAGE THE
!   UPPER BRANCH OF THE GROWTH CURVE SHOULD RECEIVE A HIGHER
!   WEIGHTING THAN THE LOWER BRANCH. AT LOW HUMIDITIES THE
!   CONVERSE APPLIES.
!   NOTE THAT A CEILING IS PLACED AT A HUMIDITY OF 0.995,
!   BECAUSE OF THE HIGHLY NONLINEAR GROWTH AS SATURATION
!   IS APPROACHED. IN ANY CASE, AT SUCH HIGH HUMIDITIES THE
!   IMPACT OF AEROSOL WILL BE SWAMPED BY THE EFFECT OF CLOUD.
!
!   NOTE THAT FOR THE PURPOSE OF DILUTING THE REFRACTIVE INDEX
!   OF A PARTICLE, IT IS STILL ASSUMED THAT ALL PARTICLES GROW
!   BY THE SAME FACTOR. (THIS IS RETURNED AS GROWTH_FACTOR.)
!   THIS IS NOT A SERIOUS INCONSISTENCY, BECAUSE BY THE TIME
!   BETA IS SIGNIFICANTLY DIFFERENT FROM UNITY, THE PARTICLES
!   ARE VERY DILUTE, SO THAT SMALL ERRORS IN THE AMOUNT OF WATER
!   DO NOT MAKE A LARGE DIFFERENCE TO THE REFRACTIVE INDEX.
!
!   EXTRA COMMENTS REGARDING SODIUM CHLORIDE.
!   -----------------------------------------
!
!   THE CODE FOR SODIUM CHLORIDE DIFFERS FROM THE AMMONIUM
!   SULPHATE CASE IN THE FOLLOWING WAYS.
!   THE DELIQUESCENCE POINT IS TAKEN AS 0.75.
!   THE EFFLORESCENCE POINT IS TAKEN AS 0.42.
!   (SEE TANG ET AL 1977, J.AEROSOL SCI. 8, 149-159.)
!   ALPHA FOR SODIUM CHLORIDE IS 1.35 TIMES ITS VALUE FOR
!   AMMONIUM SULPHATE. (BETA IS THE SAME VALUE FOR BOTH
!   SUBSTANCES.)
!   NOTE THAT WE ASSUME IT IS OK TO EXTEND THE FORMULAE FOR
!   ALPHA AND BETA DOWN FROM 0.81 TO 0.75. IN VIEW OF THE
!   GRAPHS OF THESE FUNCTIONS IN FIGS 3 AND 4 OF FITZGERALD'S
!   PAPER, THIS IS PROBABLY A HARMLESS ASSUMPTION. 
!
!   EXTRA COMMENTS REGARDING AMMONIUM NITRATE.
!   -----------------------------------------
!   THE CODE FOR AMMONIUM NITRATE DIFFERS FROM THE AMMONIUM
!   SULPHATE CASE IN THE FOLLOWING WAYS.
!   THE DELIQUESCENCE POINT IS TAKEN AS 0.61.
!
!   EXTRA COMMENTS REGARDING BIOMASS, BIOGENIC, AND FOSSIL-FUEL
!   ORGANIC CARBON AEROSOLS:
!   HYGROSCOPIC GROWTH FOR THOSE AEROSOL TYPE IS PRESCRIBED
!   USING ARRAYS OF GROWTH FACTORS AS A FUNCTION OF RELATIVE
!   HUMIDITY. FOSSIL-FUEL ORGANIC CARBON AEROSOLS USE THE SAME
!   FACTORS AS BIOMASS-BURNING.
!
!   WARNING.
!   --------
!
!   A SPECIAL PATCH OF CODE HAS BEEN ADDED TO GENERATE VALUES AT
!   A HUMIDITY OF 1.0. IT SHOULD BE NOTED THAT THIS IS OUTSIDE 
!   THE ALLOWABLE RANGE OF FITZGERALD'S SCHEME. THE VALUES 
!   RETURNED HAVE NO RELATION TO WHAT ACTUALLY HAPPENS AT 
!   SATURATION IN REALITY. THE POINT OF THIS PATCH IS THAT VALUES 
!   AT A HUMIDITY OF UNITY ARE NECESSARY FOR THE LOOK-UP TABLE.
!
!   THE VALUES OF ALPHA AND BETA AT HUMIDITY 1.0 ARE CHOSEN BY
!   EXTRAPOLATING LINEARLY FROM THE VALUES AT 0.995 AND A LOWER
!   HUMIDITY, LASTHUM, THAT SHOULD BE SET TO THE LAST HUMIDITY
!   VALUE BEFORE 1.0 IN THE LOOKUP TABLE. THE EFFECT OF THIS
!   IS THAT INTERPOLATION BETWEEN HUMIDITY VALUES LASTHUM AND 1.0
!   WILL GIVE THE VALUES OF ALPHA AND BETA THAT WOULD BE OBTAINED
!   BY LINEAR INTERPOLATION BETWEEN LASTHUM AND 0.995.
!
!   THIS IS STILL NOT EXACTLY IDEAL, BUT IT IS THE BEST I CAN DO
!   FOR THE TIME BEING.
!
!
!- ---------------------------------------------------------------------
      SUBROUTINE grow_particles(humidity,ln_r0_ln,ln_sigma_ln
     &          ,ichem_type,growth_factor,lasthum,ierr)
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE rad_pcf
!
!
      IMPLICIT NONE
!
!
!
      REAL  (RealK) ::
     &  humidity,
!              relative humidity
     &  ln_r0_ln,
!              log of mean of log-normal distn.
     &  ln_sigma_ln,
!              log of st.dev. of log-normal distn.
     &  growth_factor,
!              growth factor for dilution.
     &  alpha,
!              fitzgerald's alpha parameter
     &  beta,
!              fitzgerald's beta parameter
     &  phi,
!              fitzgerald's phi parameter
     &  alpha81,
!              value of alpha at humidity 0.81.
     &  beta81,                
!              value of beta at humidity 0.81.
     &  power
!              exponent for units conversion factor.
!
!     declarations for special patch.
!
      REAL  (RealK) ::
     & lasthum,
!              last lookup humidity value below 1.0
     & alpha_last,
!              value of alpha at lasthum
     & beta_last
!              value of beta at lasthum
!
      INTEGER ierr             ! error flag
      INTEGER ichem_type       ! indicator for chemical type.
!     allowed values are:  ip_ammonium_sulphate.
!                          ip_accum_sulphate.
!                          ip_aitken_sulphate.
!                          ip_sodium_chloride.
!                          ip_seasalt_jet.
!                          ip_seasalt_film.
!                          ip_biomass_1
!                          ip_biomass_2
!                          ip_biogenic
!                          ip_ocff_fresh
!                          ip_ocff_aged
!                          ip_nitrate
!                          ip_delta
!
!     Growth factors for biomass-burning aerosols as a function of relative
!     humidity. They were computed to match water uptake properties observed 
!     during the SAFARI campaign [Magi and Hobbs, 2003]. They are given for
!     21 values of relative humidity, ranging from 0 to 100% with a step of
!     5%.
!
      INTEGER, PARAMETER :: nhum = 21
      INTEGER i_hum
      REAL   (RealK), PARAMETER :: humidity_step = 0.05_RealK
!
      REAL   (RealK), PARAMETER ::
     &    biom_gf1(nhum) = 
!           fresh biomass
     &     (/ 1.0000_RealK, 1.0000_RealK, 1.0000_RealK, 
     &        1.0000_RealK, 1.0000_RealK, 1.0025_RealK, 
     &        1.0045_RealK, 1.0070_RealK, 1.0085_RealK, 
     &        1.0100_RealK, 1.0180_RealK, 1.0280_RealK,
     &        1.0430_RealK, 1.0600_RealK, 1.0850_RealK, 
     &        1.1100_RealK, 1.1450_RealK, 1.1870_RealK, 
     &        1.2350_RealK, 1.2850_RealK, 1.343_RealK /)
     
      REAL   (RealK), PARAMETER ::
     &    biom_gf2(nhum) = 
!           aged biomass
     &     (/ 1.000_RealK, 1.000_RealK, 1.000_RealK, 
     &        1.000_RealK, 1.000_RealK, 1.005_RealK, 
     &        1.015_RealK, 1.025_RealK, 1.040_RealK, 
     &        1.060_RealK, 1.080_RealK, 1.115_RealK,
     &        1.150_RealK, 1.185_RealK, 1.230_RealK, 
     &        1.280_RealK, 1.333_RealK, 1.390_RealK, 
     &        1.447_RealK, 1.507_RealK, 1.573_RealK /)

!
!     Growth factors for biogenic aerosols as a function of relative humidity.
!     They are taken from Varutbangkul et al.: Hygroscopicity of secondary 
!     organic aerosols formed by oxidation of cycloalkenes, monoterpenes, 
!     sesquiterpenes, and related compounds, Atmos. Chem. Phys., 6, 2367-2388, 
!     doi:10.5194/acp-6-2367-2006, 2006. They are given for 21 values of 
!     relative humidity, ranging from 0 to 100% with a step of 5%.
!
      REAL    (RealK) ::
     &    biogenic_gf(nhum) =
     &     (/ 1.00000_RealK, 1.00055_RealK, 1.00113_RealK,
     &        1.00186_RealK, 1.00281_RealK, 1.00401_RealK,
     &        1.00553_RealK, 1.00744_RealK, 1.00983_RealK,
     &        1.01279_RealK, 1.01646_RealK, 1.02100_RealK,
     &        1.02663_RealK, 1.03366_RealK, 1.04257_RealK,
     &        1.05409_RealK, 1.06959_RealK, 1.09195_RealK,
     &        1.12896_RealK, 1.21524_RealK, 1.69764_RealK /)
     
!     statement function specifications.
!
!
!     select type of nucleus.
!
      IF ( ( ichem_type  ==  IP_ammonium_sulphate ) .OR.
     &     ( ichem_type  ==  IP_accum_sulphate ) .OR.
     &     ( ichem_type  ==  IP_aitken_sulphate ) ) THEN
!
!     ammonium sulphate section.
!
        alpha = 1.0_RealK
        beta = 1.0_RealK
        growth_factor = 1.0_RealK
!
!     specify alpha.
!     no growth takes place below a humidity of 0.30.
!
        IF ( humidity  >=  0.3_RealK ) THEN
!
          IF ( humidity  <  0.81_RealK ) THEN
!
!     we have to be careful here to use the actual value of alpha
!     at 0.81, taking into account the unit conversion factor.
!     otherwise the function would become discontinuous at 0.81.
!
            power = 6.0_RealK*( bfunc(0.81_RealK) 
     &        - 1.0_RealK )
            alpha81 = afunc1(0.81_RealK)*(10.0_RealK**power)
            alpha = 1.0_RealK + (humidity-0.3_RealK)
     &         *(alpha81-1.0_RealK)/0.51_RealK
          ELSE IF ( humidity  <=  0.97_RealK ) THEN
            alpha = afunc1(humidity)
          ELSE IF ( humidity  <=  0.995_RealK ) THEN
            phi = pfunc(humidity)
            alpha = afunc2(humidity,phi)
          ELSE IF ( humidity  <  0.99999_RealK ) THEN
!
!     the maximum value of 0.995 for the humidity leads to
!     this value of alpha.
!
            phi = pfunc(0.995_RealK)
            alpha = afunc2(0.995_RealK,phi)
          ELSE
!
!     special patch for unit humidity. (see health warning.)
!     note that we again have to put in the unit conversion
!     factor, in order to get the extrapolation right.
!
          IF ( lasthum  <=  0.3_RealK ) THEN
              alpha_last = 1.0_RealK
          ELSE IF ( lasthum  <  0.81_RealK ) THEN
              power = 6.0_RealK*( bfunc(0.81_RealK) 
     &          - 1.0_RealK )
              alpha81 = afunc1(0.81_RealK)*(10.0_RealK**power)
              alpha_last = 1.0_RealK + (lasthum-0.3_RealK)
     &           *(alpha81-1.0_RealK)/0.51_RealK
          ELSE IF ( lasthum  <=  0.97_RealK ) THEN
              power = 6.0_RealK*( bfunc(lasthum) - 1.0_RealK )
              alpha_last = afunc1(lasthum)*(10.0_RealK**power)
          ELSE IF ( lasthum  <  0.995_RealK ) THEN
              phi = pfunc(lasthum)
              power = 6.0_RealK*( bfunc(lasthum) - 1.0_RealK )
              alpha_last = afunc2(lasthum,phi)*(10.0_RealK**power)
          ELSE
              WRITE(iu_err,998)
              ierr = i_err_fatal
              RETURN
          ENDIF
!
!     set alpha to its value at 0.995 before extrapolating to 1.
!     note that we again have to put in the unit conversion
!     factor, in order to get the extrapolation right.
!
            phi = pfunc(0.995_RealK)
            power = 6.0_RealK*( bfunc(0.995_RealK) 
     &              - 1.0_RealK )
            alpha = afunc2(0.995_RealK,phi)*(10.0_RealK**power)
            alpha = ( (1.0_RealK-lasthum)*alpha 
     &              - 0.005_RealK*alpha_last )
     &              / (0.995_RealK - lasthum)
          ENDIF
        ENDIF
!
!     specify beta.
!     below the deliquescence point, we interpolate beta
!     in the same way as alpha. (this is a change from
!     version 3, in which it was assumed that all
!     particles grew by the same factor, so beta=1.0.)
!
        IF ( humidity  >=  0.3_RealK .AND. 
     &       humidity  <  0.81_RealK ) THEN
          beta81 = bfunc(0.81_RealK)
          beta = 1.0_RealK + (humidity-0.3_RealK)
     &      *(beta81-1.0_RealK)/0.51_RealK
        ENDIF
        IF ( humidity  >=  0.81_RealK ) THEN
          beta = bfunc(humidity)
        ENDIF
        IF ( humidity  >  0.995_RealK ) THEN
          beta = bfunc(0.995_RealK)
        ENDIF
!
!     special patch for unit humidity. (see health warning.)
!
        IF ( humidity  >=  0.99999_RealK ) THEN
          IF ( lasthum  <=  0.3_RealK ) THEN
            beta_last = 1.0_RealK
          ELSE IF ( lasthum  <  0.81_RealK ) THEN
            beta81 = bfunc(0.81_RealK)
            beta_last = 1.0_RealK 
     &        + (lasthum-0.3_RealK)
     &        *(beta81-1.0_RealK)/0.51_RealK
          ELSE IF ( lasthum  <  0.995_RealK ) THEN
            beta_last = bfunc(lasthum)
          ELSE
            WRITE(iu_err,998)
            ierr = i_err_fatal
            RETURN
          ENDIF
!
!     now extrapolate to unit humidity.
!
          beta = ( (1.0_RealK-lasthum)*beta 
     &              - 0.005_RealK*beta_last )
     &              / (0.995_RealK - lasthum)
        ENDIF
!
!     increase logs of distribution mean and standard deviation.
!
        ln_r0_ln = beta*ln_r0_ln + log(alpha)
        ln_sigma_ln = beta*ln_sigma_ln
        growth_factor = alpha
!
!     now correct for the unit conversion problem when humidity
!     is at or above 0.81. this arises because fitzgerald's
!     formula is derived for particle radii in microns.
!     we do not make this correction below humidity=0.81,
!     because the true values of alpha (i.e. already corrected)
!     were used in the interpolation. this is the only
!     significant change from version 4 to version 5.
!
        IF ( humidity >= 0.81_RealK .AND. 
     &       humidity < 9.9999e-01_RealK ) THEN
          power = 6.0_RealK*(beta - 1.0_RealK)
          ln_r0_ln = ln_r0_ln + power*log(10.0_RealK)
          growth_factor = growth_factor*(10.0_RealK**power)
        ENDIF
!
! -----------------------------------------------------------------
!     sodium chloride section.
!
!     caution: in this section the variables alpha81 and beta81
!     are used for values at relative humidity 0.75.
!
      ELSE IF ( ( ichem_type  ==  IP_sodium_chloride ) .OR.
     &          ( ichem_type  ==  IP_seasalt_jet ) .OR.
     &          ( ichem_type  ==  IP_seasalt_film ) ) THEN
! 
        alpha = 1.0_RealK
        beta = 1.0_RealK
        growth_factor = 1.0_RealK
!
!     specify alpha.
!     no growth takes place below a humidity of 0.42.
!
        IF ( humidity  >=  0.42_RealK ) THEN
!
          IF ( humidity  <  0.75_RealK ) THEN
!
!     we have to be careful here to use the actual value of alpha
!     at 0.75, taking into account the unit conversion factor.
!     otherwise the function would become discontinuous at 0.75.
!
            power = 6.0_RealK*( bfunc(0.75_RealK) 
     &        - 1.0_RealK )
            alpha81 = cfunc1(0.75_RealK)*(10.0_RealK**power)
            alpha = 1.0_RealK + (humidity-0.42_RealK)
     &         *(alpha81-1.0_RealK)/0.33_RealK
          ELSE IF ( humidity  <=  0.97_RealK ) THEN
            alpha = cfunc1(humidity)
          ELSE IF ( humidity  <=  0.995_RealK ) THEN
            phi = pfunc(humidity)
            alpha = cfunc2(humidity,phi)
          ELSE IF ( humidity  <  0.99999_RealK ) THEN
!
!     the maximum value of 0.995 for the humidity leads to
!     this value of alpha.
!
            phi = pfunc(0.995_RealK)
            alpha = cfunc2(0.995_RealK,phi)
          ELSE
!
!     special patch for unit humidity. (see health warning.)
!     note that we again have to put in the unit conversion
!     factor, in order to get the extrapolation right.
!
          IF ( lasthum  <=  0.42_RealK ) THEN
              alpha_last = 1.0_RealK
          ELSE IF ( lasthum  <  0.75_RealK ) THEN
              power = 6.0_RealK*( bfunc(0.75_RealK) 
     &          - 1.0_RealK )
              alpha81 = cfunc1(0.75_RealK)*(10.0_RealK**power)
              alpha_last = 1.0_RealK + (lasthum-0.42_RealK)
     &           *(alpha81-1.0_RealK)/0.33_RealK
          ELSE IF ( lasthum  <=  0.97_RealK ) THEN
              power = 6.0_RealK*( bfunc(lasthum) - 1.0_RealK )
              alpha_last = cfunc1(lasthum)*(10.0_RealK**power)
          ELSE IF ( lasthum  <  0.995_RealK ) THEN
              phi = pfunc(lasthum)
              power = 6.0_RealK*( bfunc(lasthum) - 1.0_RealK )
              alpha_last = cfunc2(lasthum,phi)*(10.0_RealK**power)
          ELSE
              WRITE(iu_err,998)
              ierr = i_err_fatal
              RETURN
          ENDIF
!
!     set alpha to its value at 0.995 before extrapolating to 1.
!     note that we again have to put in the unit conversion
!     factor, in order to get the extrapolation right.
!
            phi = pfunc(0.995_RealK)
            power = 6.0_RealK*( bfunc(0.995_RealK) 
     &        - 1.0_RealK )
            alpha = cfunc2(0.995_RealK,phi)*(10.0_RealK**power)
            alpha = ( (1.0_RealK-lasthum)*alpha 
     &              - 0.005_RealK*alpha_last )
     &              / (0.995_RealK - lasthum)
          ENDIF
        ENDIF
!
!     specify beta.
!     below the deliquescence point, we interpolate beta
!     in the same way as alpha. (this is a change from
!     version 3, in which it was assumed that all
!     particles grew by the same factor, so beta=1.0.)
!
        IF ( humidity  >=  0.42_RealK .AND. 
     &       humidity  <  0.75_RealK ) THEN
          beta81 = bfunc(0.75_RealK)
          beta = 1.0_RealK + (humidity-0.42_RealK)
     &      *(beta81-1.0_RealK)/0.33_RealK
        ENDIF
        IF ( humidity  >=  0.75_RealK ) THEN
          beta = bfunc(humidity)
        ENDIF
        IF ( humidity  >  0.995_RealK ) THEN
          beta = bfunc(0.995_RealK)
        ENDIF
!
!     special patch for unit humidity. (see health warning.)
!
        IF ( humidity  >=  0.99999_RealK ) THEN
          IF ( lasthum  <=  0.42_RealK ) THEN
            beta_last = 1.0_RealK
          ELSE IF ( lasthum  <  0.75_RealK ) THEN
            beta81 = bfunc(0.75_RealK)
            beta_last = 1.0_RealK 
     &        + (lasthum-0.42_RealK)*(beta81-1.0_RealK)
     &        /0.33_RealK
          ELSE IF ( lasthum  <  0.995_RealK ) THEN
            beta_last = bfunc(lasthum)
          ELSE
            WRITE(iu_err,998)
            ierr = i_err_fatal
            RETURN
          ENDIF
!
!     now extrapolate to unit humidity.
!
          beta = ( (1.0_RealK-lasthum)*beta 
     &              - 0.005_RealK*beta_last )
     &              / (0.995_RealK - lasthum)
        ENDIF
!
!     increase logs of distribution mean and standard deviation.
!
        ln_r0_ln = beta*ln_r0_ln + log(alpha)
        ln_sigma_ln = beta*ln_sigma_ln
        growth_factor = alpha
!
!     now correct for the unit conversion problem when humidity
!     is at or above 0.75. this arises because fitzgerald's
!     formula is derived for particle radii in microns.
!     we do not make this correction below humidity=0.75,
!     because the true values of alpha (i.e. already corrected)
!     were used in the interpolation. this is the only
!     significant change from version 4 to version 5.
!
        IF ( humidity >= 0.75_RealK .AND. 
     &       humidity < 9.9999e-01_RealK ) THEN
          power = 6.0_RealK*(beta - 1.0_RealK)
          ln_r0_ln = ln_r0_ln + power*log(10.0_RealK)
          growth_factor = growth_factor*(10.0_RealK**power)
        ENDIF

!
! Ammonium nitrate (taken from Fitzgerald, 1975).
!
      ELSE IF ( ( ichem_type == IP_NITRATE ) .OR.
     &          ( ichem_type == IP_DELTA ) ) THEN
      
        alpha = 1.0_RealK
        beta = 1.0_RealK
        growth_factor = 1.0_RealK
        
!
! Specify alpha.
! No growth takes place below a humidity of 0.30
! (the efflorescence point varies depending on
!  conditions, and might not exist at all)
! The deliquescence point is taken at 61% RH.
!
        IF ( humidity  >=  0.30_RealK ) THEN
        
          IF ( humidity  <  0.61_RealK ) THEN
!
!     we have to be careful here to use the actual value of alpha
!     at 0.61, taking into account the unit conversion factor.
!     otherwise the function would become discontinuous at 0.61.
!     (Note that the 0.31 is simply RH(deliq) - RH(effl)).
!
            power = 6.0_RealK*( bfunc(0.61_RealK) 
     &        - 1.0_RealK )
            alpha81 = nfunc1(0.61_RealK)*(10.0_RealK**power)
            alpha = 1.0_RealK + (humidity-0.3_RealK)
     &         *(alpha81-1.0_RealK)/0.31_RealK
          ELSE IF ( humidity  <=  0.97_RealK ) THEN
            alpha = nfunc1(humidity)
          ELSE IF ( humidity  <=  0.995_RealK ) THEN
            phi = pfunc(humidity)
            alpha = nfunc2(humidity,phi)
          ELSE IF ( humidity  <  0.99999_RealK ) THEN
!
!     the maximum value of 0.995 for the humidity leads to
!     this value of alpha.
!
            phi = pfunc(0.995_RealK)
            alpha = nfunc2(0.995_RealK,phi)
          ELSE
!
!     special patch for unit humidity. (see health warning.)
!     note that we again have to put in the unit conversion
!     factor, in order to get the extrapolation right.
!
          IF ( lasthum  <=  0.3_RealK ) THEN
              alpha_last = 1.0_RealK
          ELSE IF ( lasthum  <  0.61_RealK ) THEN
              power = 6.0_RealK*( bfunc(0.61_RealK) 
     &          - 1.0_RealK )
              alpha81 = nfunc1(0.61_RealK)*(10.0_RealK**power)
              alpha_last = 1.0_RealK + (lasthum-0.3_RealK)
     &           *(alpha81-1.0_RealK)/0.31_RealK
          ELSE IF ( lasthum  <=  0.97_RealK ) THEN
              power = 6.0_RealK*( bfunc(lasthum) - 1.0_RealK )
              alpha_last = nfunc1(lasthum)*(10.0_RealK**power)
          ELSE IF ( lasthum  <  0.995_RealK ) THEN
              phi = pfunc(lasthum)
              power = 6.0_RealK*( bfunc(lasthum) - 1.0_RealK )
              alpha_last = nfunc2(lasthum,phi)*(10.0_RealK**power)
          ELSE
              WRITE(iu_err,998)
              ierr = i_err_fatal
              RETURN
          ENDIF
!
!     set alpha to its value at 0.995 before extrapolating to 1.
!     note that we again have to put in the unit conversion
!     factor, in order to get the extrapolation right.
!
            phi = pfunc(0.995_RealK)
            power = 6.0_RealK*( bfunc(0.995_RealK) 
     &        - 1.0_RealK )
            alpha = nfunc2(0.995_RealK,phi)*(10.0_RealK**power)
            alpha = ( (1.0_RealK-lasthum)*alpha 
     &              - 0.005_RealK*alpha_last )
     &              / (0.995_RealK - lasthum)
          ENDIF
        ENDIF
!
!     specify beta.
!     below the deliquescence point, we interpolate beta
!     in the same way as alpha. (this is a change from
!     version 3, in which it was assumed that all
!     particles grew by the same factor, so beta=1.0.)
!
        IF ( humidity  >=  0.3_RealK .AND. 
     &       humidity  <  0.61_RealK ) THEN
          beta81 = bfunc(0.61_RealK)
          beta = 1.0_RealK + (humidity-0.3_RealK)
     &      *(beta81-1.0_RealK)/0.31_RealK
        ENDIF
        IF ( humidity  >=  0.61_RealK ) THEN
          beta = bfunc(humidity)
        ENDIF
        IF ( humidity  >  0.995_RealK ) THEN
          beta = bfunc(0.995_RealK)
        ENDIF
!
!     special patch for unit humidity. (see health warning.)
!
        IF ( humidity  >=  0.99999_RealK ) THEN
          IF ( lasthum  <=  0.3_RealK ) THEN
            beta_last = 1.0_RealK
          ELSE IF ( lasthum  <  0.61_RealK ) THEN
            beta81 = bfunc(0.61_RealK)
            beta_last = 1.0_RealK 
     &        + (lasthum-0.3_RealK)
     &        *(beta81-1.0_RealK)/0.31_RealK
          ELSE IF ( lasthum  <  0.995_RealK ) THEN
            beta_last = bfunc(lasthum)
          ELSE
            WRITE(iu_err,998)
            ierr = i_err_fatal
            RETURN
          ENDIF
!
!     now extrapolate to unit humidity.
!
          beta = ( (1.0_RealK-lasthum)*beta 
     &              - 0.005_RealK*beta_last )
     &              / (0.995_RealK - lasthum)
        ENDIF
!
!     increase logs of distribution mean and standard deviation.
!
        ln_r0_ln = beta*ln_r0_ln + log(alpha)
        ln_sigma_ln = beta*ln_sigma_ln
        growth_factor = alpha
!
!     now correct for the unit conversion problem when humidity
!     is at or above 0.61. this arises because fitzgerald''s
!     formula is derived for particle radii in microns.
!     we do not make this correction below humidity=0.61,
!     because the true values of alpha (i.e. already corrected)
!     were used in the interpolation. this is the only
!     significant change from version 4 to version 5.
!
        IF ( humidity >= 0.61_RealK .AND. 
     &       humidity < 9.9999e-01_RealK ) THEN
          power = 6.0_RealK*(beta - 1.0_RealK)
          ln_r0_ln = ln_r0_ln + power*log(10.0_RealK)
          growth_factor = growth_factor*(10.0_RealK**power)
        ENDIF
!
! Biomass-burning aerosol by simply using an pre-defined array
! of growth factors (see declaration of arrays biom_gf).
! We suppose that the standard deviation of the lognormal does not vary 
! with humidity (i.e. all particle sizes grow at the same rate).
! Note that parameters ln_r0_ln and ln_sigma_ln are always for dry particles.
!
! Fossil-fuel organic carbon aerosols use the same growth factors as biomass
! burning aerosols.
!
      ELSE IF ( ( ichem_type == IP_BIOMASS_1 ) .OR.
     &          ( ichem_type == IP_BIOMASS_2 ) .OR.
     &          ( ichem_type == IP_OCFF_FRESH ) .OR.
     &          ( ichem_type == IP_OCFF_AGED  ) ) THEN

!
! get the index in the growth factor array
! (we verify that 0.0 <= humidity <= 1.0)
!

          IF( humidity .LT. 0.0_RealK .OR. 
     &        humidity .GT. 1.0_RealK) THEN
            WRITE(iu_err,*) 
     &         'Humidity ', humidity, ' has an unexpected value.'
            ierr = i_err_fatal
          ELSE
            i_hum = NINT(humidity / humidity_step) + 1
            IF (ichem_type == IP_BIOMASS_1 .OR. 
     &          ichem_type == IP_OCFF_FRESH) THEN
              growth_factor = biom_gf1(i_hum)
            ELSE
              growth_factor = biom_gf2(i_hum)
            ENDIF
            ln_r0_ln = LOG( EXP(ln_r0_ln) * growth_factor)
            ln_sigma_ln = ln_sigma_ln
          ENDIF 
!
! Biogenic aerosol. Like for biomass-burning, pre-defined growth factors are 
! used.
!
      ELSE IF ( ichem_type == IP_BIOGENIC ) THEN

!
! get the index in the growth factor array
! (we verify that 0.0 <= humidity <= 1.0)
!

          IF( humidity .LT. 0.0_RealK .OR.
     &        humidity .GT. 1.0_RealK) THEN
            WRITE(iu_err,*)
     &         'Humidity ', humidity, ' has an unexpected value.'
            ierr = i_err_fatal
          ELSE
            i_hum = NINT(humidity / humidity_step) + 1
            growth_factor = biogenic_gf(i_hum)
            ln_r0_ln = LOG( EXP(ln_r0_ln) * growth_factor)
            ln_sigma_ln = ln_sigma_ln  
          ENDIF
          
      ELSE
!
!     trap incorrect types here.
!
        WRITE(iu_err,999) ichem_type
        ierr = i_err_fatal
!
      ENDIF
!
      RETURN
!
998   format(1x,'lasthum should be less than 0.995')
999   format(1x,'unavailable aerosol composition type',i5)
!
!
!
      contains
!
      FUNCTION afunc1(h) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: r
        r = 1.2_RealK
     &    *EXP( (0.066_RealK*h)/(1.058_RealK - h) )
      END FUNCTION
!
      FUNCTION bfunc(h) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: r
        r = EXP( (7.7e-04_RealK*h)/(1.009_RealK - h) )
      END FUNCTION
!
      FUNCTION pfunc(h) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: r
        r = 1.058_RealK - ( (0.0155_RealK*(h 
     &                      - 0.97_RealK))
     &                      / (1.02_RealK - h**1.4_RealK) )
      END FUNCTION
!
      FUNCTION afunc2(h, p) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: p
        REAL  (RealK) :: r
        r = 1.2_RealK*EXP( (0.066_RealK*h)/(p - h) )
      END FUNCTION
!
      FUNCTION cfunc1(h) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: r
        r = afunc1(h)*1.35_RealK
      END FUNCTION
!
      FUNCTION cfunc2(h, p) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: p
        REAL  (RealK) :: r
        r = afunc2(h,p)*1.35_RealK
      END FUNCTION
!
      FUNCTION nfunc1(h) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: r
        r = afunc1(h)*1.06_RealK
      END FUNCTION
!
      FUNCTION nfunc2(h, p) RESULT(r)
        REAL  (RealK) :: h
        REAL  (RealK) :: p
        REAL  (RealK) :: r
        r = afunc2(h,p)*1.06_RealK
      END FUNCTION
!
!
      END
