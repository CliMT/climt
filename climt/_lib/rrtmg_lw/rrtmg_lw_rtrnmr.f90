!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/11/12 20:52:26 $
!
      module rrtmg_lw_rtrnmr

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, ngs
      use rrlw_tbl, only: tblint, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrrtx, hnamrtx

      implicit none

      contains

!-----------------------------------------------------------------------------
      subroutine rtrnmr(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                        cldfrac, taucloud, planklay, planklev, plankbnd, &
                        pwvcm, fracs, taut, & 
                        totuflux, totdflux, fnet, htr, &
                        totuclfl, totdclfl, fnetc, htrc, &
                        totuflux_bnd, totdflux_bnd, fnet_bnd, htr_bnd, &
                        totuclfl_bnd, totdclfl_bnd, fnetc_bnd, htrc_bnd, &
                        idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with a maximum-random cloud overlap method.
!  This subroutine also provides the optional capability to calculate
!  the derivative of upward flux respect to surface temperature using
!  the pre-tabulated derivative of the Planck function with respect to 
!  temperature integrated over each spectral band.
!***************************************************************************

! ------- Declarations -------

! ----- Input -----
      integer(kind=im), intent(in) :: nlayers         ! total number of layers
      integer(kind=im), intent(in) :: istart          ! beginning band of calculation
      integer(kind=im), intent(in) :: iend            ! ending band of calculation
      integer(kind=im), intent(in) :: iout            ! output option flag

! Atmosphere
      real(kind=rb), intent(in) :: pz(0:)             ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=rb), intent(in) :: semiss(:)          ! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: planklay(:,:)      ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=rb), intent(in) :: planklev(0:,:)     ! 
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(in) :: plankbnd(:)        ! 
                                                      !    Dimensions: (nbndlw)
      real(kind=rb), intent(in) :: fracs(:,:)         ! 
                                                      !    Dimensions: (nlayers,ngptw)
      real(kind=rb), intent(in) :: taut(:,:)          ! gaseous + aerosol optical depths
                                                      !    Dimensions: (nlayers,ngptlw)

! Clouds
      integer(kind=im), intent(in) :: ncbands         ! number of cloud spectral bands
      real(kind=rb), intent(in) :: cldfrac(:)         ! layer cloud fraction
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: taucloud(:,:)      ! layer cloud optical depth
                                                      !    Dimensions: (nlayers,nbndlw)
      integer(kind=im), intent(in) :: idrv            ! flag for calculation of dF/dt from 
                                                      ! Planck derivative [0=off, 1=on]
      real(kind=rb), intent(in) :: dplankbnd_dt(:)    ! derivative of Planck function wrt temp
                                                      !    Dimensions: (nbndlw)

! ----- Output -----
      real(kind=rb), intent(out) :: totuflux(0:)      ! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdflux(0:)      ! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: fnet(0:)          ! net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: htr(0:)           ! longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuclfl(0:)      ! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totdclfl(0:)      ! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: fnetc(0:)         ! clear sky net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: htrc(0:)          ! clear sky longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: totuflux_bnd(0:,:)! upward longwave flux by band (w/m2)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: totdflux_bnd(0:,:)! downward longwave flux by band (w/m2)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: fnet_bnd(0:,:)    ! net longwave flux by band (w/m2)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: htr_bnd(0:,:)     ! longwave heating rate by band (k/day)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: totuclfl_bnd(0:,:)! clear sky upward longwave flux by band (w/m2)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: totdclfl_bnd(0:,:)! clear sky downward longwave flux by band (w/m2)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: fnet_bnd(0:,:)    ! clear sky net longwave flux by band (w/m2)
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=rb), intent(out) :: htrc_bnd(0:,:)    ! clear sky longwave heating rate by band (k/day)
                                                      !    Dimensions: (0:nlayers, nbndlw)
      real(kind=rb), intent(out) :: dtotuflux_dt(0:)  ! change in upward longwave flux (w/m2/k)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (0:nlayers)
      real(kind=rb), intent(out) :: dtotuclfl_dt(0:)  ! change in upward longwave flux (w/m2/k)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (0:nlayers)

! ----- Local -----
! Declarations for radiative transfer
      real(kind=rb) :: abscld(nlayers,nbndlw)
      real(kind=rb) :: atot(nlayers)
      real(kind=rb) :: atrans(nlayers)
      real(kind=rb) :: bbugas(nlayers)
      real(kind=rb) :: bbutot(nlayers)
      real(kind=rb) :: clrurad(0:nlayers)
      real(kind=rb) :: clrdrad(0:nlayers)
      real(kind=rb) :: efclfrac(nlayers,nbndlw)
      real(kind=rb) :: uflux(0:nlayers)
      real(kind=rb) :: dflux(0:nlayers)
      real(kind=rb) :: urad(0:nlayers)
      real(kind=rb) :: drad(0:nlayers)
      real(kind=rb) :: uclfl(0:nlayers)
      real(kind=rb) :: dclfl(0:nlayers)
      real(kind=rb) :: odcld(nlayers,nbndlw)

      real(kind=rb) :: secdiff(nbndlw)                 ! secant of diffusivity angle
      real(kind=rb) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)! diffusivity angle adjustment coefficients
      real(kind=rb) :: wtdiff, rec_6
      real(kind=rb) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=rb) :: odepth, odtot, odepth_rec, odtot_rec, gassrc, ttot
      real(kind=rb) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=rb) :: rad0, reflect, radlu, radclru

      real(kind=rb) :: duflux_dt(0:nlayers)
      real(kind=rb) :: duclfl_dt(0:nlayers)
      real(kind=rb) :: d_urad_dt(0:nlayers)
      real(kind=rb) :: d_clrurad_dt(0:nlayers)
      real(kind=rb) :: d_rad0_dt, d_radlu_dt, d_radclru_dt

      integer(kind=im) :: icldlyr(nlayers)             ! flag for cloud in layer
      integer(kind=im) :: ibnd, ib, iband, lay, lev, l ! loop indices
      integer(kind=im) :: igc                          ! g-point interval counter
      integer(kind=im) :: iclddn                       ! flag for cloud in down path
      integer(kind=im) :: ittot, itgas, itr            ! lookup table indices
      integer(kind=im) :: ipat(16,0:2)

! Declarations for cloud overlap adjustment
      real(kind=rb) :: faccld1(nlayers+1),faccld2(nlayers+1)
      real(kind=rb) :: facclr1(nlayers+1),facclr2(nlayers+1)
      real(kind=rb) :: faccmb1(nlayers+1),faccmb2(nlayers+1)
      real(kind=rb) :: faccld1d(0:nlayers),faccld2d(0:nlayers)
      real(kind=rb) :: facclr1d(0:nlayers),facclr2d(0:nlayers)
      real(kind=rb) :: faccmb1d(0:nlayers),faccmb2d(0:nlayers)

      real(kind=rb) :: fmax, fmin, rat1, rat2
      real(kind=rb) :: clrradd, cldradd, clrradu, cldradu, oldclr, oldcld
      real(kind=rb) :: rad, cldsrc, radmod

      integer(kind=im) :: istcld(nlayers+1),istcldd(0:nlayers)

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer
!    d_radlu_dt                   ! spectrally summed upward radiance 
!    d_radclru_dt                 ! spectrally summed clear sky upward radiance 
!    d_urad_dt                    ! upward radiance by layer
!    d_clrurad_dt                 ! clear sky upward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)
!    totuflux_bnd                 ! upward longwave flux by band (w/m2)
!    totdflux_bnd                 ! downward longwave flux by band (w/m2)
!    fnet_bnd                     ! net longwave flux by band (w/m2)
!    htr_bnd                      ! longwave heating rate by band (k/day)
!    totuclfl_bnd                 ! clear sky upward longwave flux by band (w/m2)
!    totdclfl_bnd                 ! clear sky downward longwave flux by band (w/m2)
!    fnetc_bnd                    ! clear sky net longwave flux by band (w/m2)
!    htrc_bnd                     ! clear sky longwave heating rate by band (k/day)
!    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature
!    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature


! These arrays indicate the spectral 'region' (used in the 
! calculation of ice cloud optical depths) corresponding
! to each spectral band.  See cldprop.f for more details.
      data ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
                 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

! This secant and weight corresponds to the standard diffusivity 
! angle.  This initial value is redefined below for some bands.
      data wtdiff /0.5_rb/
      data rec_6 /0.166667_rb/

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.
      data a0 / 1.66_rb,  1.55_rb,  1.58_rb,  1.66_rb, &
                1.54_rb, 1.454_rb,  1.89_rb,  1.33_rb, &
               1.668_rb,  1.66_rb,  1.66_rb,  1.66_rb, &
                1.66_rb,  1.66_rb,  1.66_rb,  1.66_rb /
      data a1 / 0.00_rb,  0.25_rb,  0.22_rb,  0.00_rb, &
                0.13_rb, 0.446_rb, -0.10_rb,  0.40_rb, &
              -0.006_rb,  0.00_rb,  0.00_rb,  0.00_rb, &
                0.00_rb,  0.00_rb,  0.00_rb,  0.00_rb /
      data a2 / 0.00_rb, -12.0_rb, -11.7_rb,  0.00_rb, &
               -0.72_rb,-0.243_rb,  0.19_rb,-0.062_rb, &
               0.414_rb,  0.00_rb,  0.00_rb,  0.00_rb, &
                0.00_rb,  0.00_rb,  0.00_rb,  0.00_rb /

      do ibnd = 1,nbndlw
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff(ibnd) = 1.66_rb
         else
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm)
           if (secdiff(ibnd) .gt. 1.80_rb) secdiff(ibnd) = 1.80_rb
           if (secdiff(ibnd) .lt. 1.50_rb) secdiff(ibnd) = 1.50_rb
         endif
      enddo

      hvrrtx = '$Revision: 1.7 $'

      urad(0) = 0.0_rb
      drad(0) = 0.0_rb
      totuflux(0) = 0.0_rb
      totdflux(0) = 0.0_rb
      clrurad(0) = 0.0_rb
      clrdrad(0) = 0.0_rb
      totuclfl(0) = 0.0_rb
      totdclfl(0) = 0.0_rb
! Initialize arrays at surface by band
      do ibnd = 1,nbndlw
         totuflux_bnd(0,ibnd) = 0.0_rb
         totdflux_bnd(0,ibnd) = 0.0_rb
         totuclfl_bnd(0,ibnd) = 0.0_rb
         totdclfl_bnd(0,ibnd) = 0.0_rb
      enddo
      if (idrv .eq. 1) then
         d_urad_dt(0) = 0.0_rb
         d_clrurad_dt(0) = 0.0_rb
         dtotuflux_dt(0) = 0.0_rb
         dtotuclfl_dt(0) = 0.0_rb
      endif

      do lay = 1, nlayers
         urad(lay) = 0.0_rb
         drad(lay) = 0.0_rb
         totuflux(lay) = 0.0_rb
         totdflux(lay) = 0.0_rb
         clrurad(lay) = 0.0_rb
         clrdrad(lay) = 0.0_rb
         totuclfl(lay) = 0.0_rb
         totdclfl(lay) = 0.0_rb
! Initialize arrays at levels by band
         do ibnd = 1,nbndlw
            totuflux_bnd(lay,ibnd) = 0.0_rb
            totdflux_bnd(lay,ibnd) = 0.0_rb
            totuclfl_bnd(lay,ibnd) = 0.0_rb
            totdclfl_bnd(lay,ibnd) = 0.0_rb
         enddo
         if (idrv .eq. 1) then
            d_urad_dt(lay) = 0.0_rb
            d_clrurad_dt(lay) = 0.0_rb
            dtotuflux_dt(lay) = 0.0_rb
            dtotuclfl_dt(lay) = 0.0_rb
         endif

         do ib = 1, ncbands
            if (cldfrac(lay) .ge. 1.e-6_rb) then
               odcld(lay,ib) = secdiff(ib) * taucloud(lay,ib)
               icldlyr(lay) = 1
            else
               odcld(lay,ib) = 0.0_rb
               icldlyr(lay) = 0
            endif
         enddo
      enddo

! Maximum/Random cloud overlap parameter

      istcld(1) = 1
      istcldd(nlayers) = 1
      do lev = 1, nlayers

         if (icldlyr(lev).eq.1) then
! Maximum/random cloud overlap
            istcld(lev+1) = 0
            if (lev .eq. nlayers) then
               faccld1(lev+1) = 0._rb
               faccld2(lev+1) = 0._rb
               facclr1(lev+1) = 0._rb
               facclr2(lev+1) = 0._rb
               faccmb1(lev+1) = 0._rb
               faccmb2(lev+1) = 0._rb
            elseif (cldfrac(lev+1) .ge. cldfrac(lev)) then
               faccld1(lev+1) = 0._rb
               faccld2(lev+1) = 0._rb
               if (istcld(lev) .eq. 1) then
                  facclr1(lev+1) = 0._rb
                  facclr2(lev+1) = 0._rb
                  if (cldfrac(lev) .lt. 1._rb) facclr2(lev+1) = &
                     (cldfrac(lev+1)-cldfrac(lev))/(1._rb-cldfrac(lev))
                  facclr2(lev) = 0._rb
                  faccld2(lev) = 0._rb
               else
                  fmax = max(cldfrac(lev),cldfrac(lev-1))
                  if (cldfrac(lev+1) .gt. fmax) then
                     facclr1(lev+1) = rat2
                     facclr2(lev+1) = (cldfrac(lev+1)-fmax)/(1._rb-fmax)
                  elseif (cldfrac(lev+1) .lt. fmax) then
                     facclr1(lev+1) = (cldfrac(lev+1)-cldfrac(lev))/ &
                        (cldfrac(lev-1)-cldfrac(lev))
                     facclr2(lev+1) = 0._rb
                  else
                     facclr1(lev+1) = rat2
                     facclr2(lev+1) = 0._rb
                  endif
               endif
               if (facclr1(lev+1).gt.0._rb .or. facclr2(lev+1).gt.0._rb) then
                  rat1 = 1._rb
                  rat2 = 0._rb
               else
                  rat1 = 0._rb
                  rat2 = 0._rb
               endif
            else
               facclr1(lev+1) = 0._rb
               facclr2(lev+1) = 0._rb
               if (istcld(lev) .eq. 1) then
                  faccld1(lev+1) = 0._rb
                  faccld2(lev+1) = (cldfrac(lev)-cldfrac(lev+1))/cldfrac(lev)

                  facclr2(lev) = 0._rb
                  faccld2(lev) = 0._rb
               else
                  fmin = min(cldfrac(lev),cldfrac(lev-1))
                  if (cldfrac(lev+1) .le. fmin) then
                     faccld1(lev+1) = rat1
                     faccld2(lev+1) = (fmin-cldfrac(lev+1))/fmin
                  else
                     faccld1(lev+1) = (cldfrac(lev)-cldfrac(lev+1))/(cldfrac(lev)-fmin)
                     faccld2(lev+1) = 0._rb
                  endif
               endif
               if (faccld1(lev+1).gt.0._rb .or. faccld2(lev+1).gt.0._rb) then
                  rat1 = 0._rb
                  rat2 = 1._rb
               else
                  rat1 = 0._rb
                  rat2 = 0._rb
               endif
            endif
            faccmb1(lev+1) = facclr1(lev+1) * faccld2(lev) * cldfrac(lev-1) 
            faccmb2(lev+1) = faccld1(lev+1) * facclr2(lev) * (1._rb - cldfrac(lev-1)) 
         else
            istcld(lev+1) = 1
         endif
      enddo

      do lev = nlayers, 1, -1
         if (icldlyr(lev).eq.1) then
            istcldd(lev-1) = 0
            if (lev .eq. 1) then
               faccld1d(lev-1) = 0._rb
               faccld2d(lev-1) = 0._rb
               facclr1d(lev-1) = 0._rb
               facclr2d(lev-1) = 0._rb
               faccmb1d(lev-1) = 0._rb
               faccmb2d(lev-1) = 0._rb
            elseif (cldfrac(lev-1) .ge. cldfrac(lev)) then
               faccld1d(lev-1) = 0._rb
               faccld2d(lev-1) = 0._rb
               if (istcldd(lev) .eq. 1) then
                  facclr1d(lev-1) = 0._rb
                  facclr2d(lev-1) = 0._rb
                  if (cldfrac(lev) .lt. 1._rb) facclr2d(lev-1) = &
                     (cldfrac(lev-1)-cldfrac(lev))/(1._rb-cldfrac(lev))
                  facclr2d(lev) = 0._rb
                  faccld2d(lev) = 0._rb
               else
                  fmax = max(cldfrac(lev),cldfrac(lev+1))
                  if (cldfrac(lev-1) .gt. fmax) then
                     facclr1d(lev-1) = rat2
                     facclr2d(lev-1) = (cldfrac(lev-1)-fmax)/(1._rb-fmax)
                  elseif (cldfrac(lev-1) .lt. fmax) then
                     facclr1d(lev-1) = (cldfrac(lev-1)-cldfrac(lev))/ &
                        (cldfrac(lev+1)-cldfrac(lev))
                     facclr2d(lev-1) = 0.
                  else
                     facclr1d(lev-1) = rat2
                     facclr2d(lev-1) = 0._rb
                  endif
               endif
               if (facclr1d(lev-1).gt.0._rb .or. facclr2d(lev-1).gt.0._rb)then
                  rat1 = 1._rb
                  rat2 = 0._rb
               else
                  rat1 = 0._rb
                  rat2 = 0._rb
               endif
            else
               facclr1d(lev-1) = 0._rb
               facclr2d(lev-1) = 0._rb
               if (istcldd(lev) .eq. 1) then
                  faccld1d(lev-1) = 0._rb
                  faccld2d(lev-1) = (cldfrac(lev)-cldfrac(lev-1))/cldfrac(lev)
                  facclr2d(lev) = 0._rb
                  faccld2d(lev) = 0._rb
               else
                  fmin = min(cldfrac(lev),cldfrac(lev+1))
                  if (cldfrac(lev-1) .le. fmin) then
                     faccld1d(lev-1) = rat1
                     faccld2d(lev-1) = (fmin-cldfrac(lev-1))/fmin
                  else
                     faccld1d(lev-1) = (cldfrac(lev)-cldfrac(lev-1))/(cldfrac(lev)-fmin)
                     faccld2d(lev-1) = 0._rb
                  endif
               endif
               if (faccld1d(lev-1).gt.0._rb .or. faccld2d(lev-1).gt.0._rb)then
                  rat1 = 0._rb
                  rat2 = 1._rb
               else
                  rat1 = 0._rb
                  rat2 = 0._rb
               endif
            endif
            faccmb1d(lev-1) = facclr1d(lev-1) * faccld2d(lev) * cldfrac(lev+1) 
            faccmb2d(lev-1) = faccld1d(lev-1) * facclr2d(lev) * (1._rb - cldfrac(lev+1))
         else
            istcldd(lev-1) = 1
         endif
      enddo

      igc = 1
! Loop over frequency bands.
      do iband = istart, iend

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1
         if (ncbands .eq. 1) then
            ib = ipat(iband,0)
         elseif (ncbands .eq.  5) then
            ib = ipat(iband,1)
         elseif (ncbands .eq. 16) then
            ib = ipat(iband,2)
         endif

! Loop over g-channels.
 1000    continue

! Radiative transfer starts here.
         radld = 0._rb
         radclrd = 0._rb
         iclddn = 0

! Downward radiative transfer loop.  

         do lev = nlayers, 1, -1
               plfrac = fracs(lev,igc)
               blay = planklay(lev,iband)
               dplankup = planklev(lev,iband) - blay
               dplankdn = planklev(lev-1,iband) - blay
               odepth = secdiff(iband) * taut(lev,igc)

               if (odepth .lt. 0.0_rb) odepth = 0.0_rb
! Cloudy layer
               if (icldlyr(lev).eq.1) then
                  iclddn = 1
                  odtot = odepth + odcld(lev,ib)
                  if (odtot .lt. 0.06_rb) then
                     atrans(lev) = odepth - 0.5_rb*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     atot(lev) =  odtot - 0.5_rb*odtot*odtot
                     odtot_rec = rec_6*odtot
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                  
                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)
                  elseif (odepth .le. 0.06_rb) then
                     atrans(lev) = odepth - 0.5_rb*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     odtot = odepth + odcld(lev,ib)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_rb
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     atot(lev) = 1._rb - exp_tbl(ittot)

                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  else
                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5_rb
                     odepth = tau_tbl(itgas)
                     atrans(lev) = 1._rb - exp_tbl(itgas)
                     tfacgas = tfn_tbl(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odcld(lev,ib)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_rb
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1._rb - exp_tbl(ittot)

                     bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  endif

                  if (istcldd(lev) .eq. 1) then
                     cldradd = cldfrac(lev) * radld
                     clrradd = radld - cldradd
                     oldcld = cldradd
                     oldclr = clrradd
                     rad = 0._rb
                  endif
                  ttot = 1._rb - atot(lev)
                  cldsrc = bbdtot * atot(lev)
                  cldradd = cldradd * ttot + cldfrac(lev) * cldsrc
                  clrradd = clrradd * (1._rb-atrans(lev)) + (1._rb - cldfrac(lev)) * gassrc
                  radld = cldradd + clrradd
                  drad(lev-1) = drad(lev-1) + radld

                  radmod = rad * &
                       (facclr1d(lev-1) * (1.-atrans(lev)) + &
                       faccld1d(lev-1) *  ttot) - &
                       faccmb1d(lev-1) * gassrc + &
                       faccmb2d(lev-1) * cldsrc

                  oldcld = cldradd - radmod
                  oldclr = clrradd + radmod
                  rad = -radmod + facclr2d(lev-1)*oldclr - faccld2d(lev-1)*oldcld
                  cldradd = cldradd + rad
                  clrradd = clrradd - rad
! Clear layer
               else
                  if (odepth .le. 0.06_rb) then
                     atrans(lev) = odepth-0.5_rb*odepth*odepth
                     odepth = rec_6*odepth
                     bbd = plfrac*(blay+dplankdn*odepth)
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)
                  else
                     tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5_rb
                     transc = exp_tbl(itr)
                     atrans(lev) = 1._rb-transc
                     tausfac = tfn_tbl(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                  endif   
                  radld = radld + (bbd-radld)*atrans(lev)
                  drad(lev-1) = drad(lev-1) + radld
                endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                 if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev) 
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd
                  else
                     radclrd = radld
                     clrdrad(lev-1) = drad(lev-1)
                  endif
            enddo

! Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.
!  Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
!  they are defined in subroutine setcoef. 

         rad0 = fracs(1,igc) * plankbnd(iband)
         if (idrv .eq. 1) then
            d_rad0_dt = fracs(1,igc) * dplankbnd_dt(iband)
         endif

!  Add in reflection of surface downward radiance.
         reflect = 1._rb - semiss(iband)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd

! Upward radiative transfer loop.

         urad(0) = urad(0) + radlu
         clrurad(0) = clrurad(0) + radclru
         if (idrv .eq. 1) then
            d_radlu_dt = d_rad0_dt
            d_urad_dt(0) = d_urad_dt(0) + d_radlu_dt
            d_radclru_dt = d_rad0_dt
            d_clrurad_dt(0) = d_clrurad_dt(0) + d_radclru_dt
         endif

         do lev = 1, nlayers
! Cloudy layer
            if (icldlyr(lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               if (istcld(lev) .eq. 1) then
                  cldradu = cldfrac(lev) * radlu
                  clrradu = radlu - cldradu
                  oldcld = cldradu
                  oldclr = clrradu
                  rad = 0._rb
               endif
               ttot = 1._rb - atot(lev)
               cldsrc = bbutot(lev) * atot(lev)
               cldradu = cldradu * ttot + cldfrac(lev) * cldsrc
               clrradu = clrradu * (1.0_rb-atrans(lev)) + (1._rb - cldfrac(lev)) * gassrc
! Total sky radiance
               radlu = cldradu + clrradu
               urad(lev) = urad(lev) + radlu
               radmod = rad * &
                   (facclr1(lev+1)*(1.0_rb-atrans(lev))+ &
                   faccld1(lev+1) *  ttot) - &
                   faccmb1(lev+1) * gassrc + &
                   faccmb2(lev+1) * cldsrc
               oldcld = cldradu - radmod
               oldclr = clrradu + radmod
               rad = -radmod + facclr2(lev+1)*oldclr - faccld2(lev+1)*oldcld
               cldradu = cldradu + rad
               clrradu = clrradu - rad
               if (idrv .eq. 1) then
                  d_radlu_dt = d_radlu_dt * cldfrac(lev) * (1.0_rb - atot(lev)) + &
                         d_radlu_dt * (1.0_rb - cldfrac(lev)) * (1.0_rb - atrans(lev))
                  d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
               endif
! Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               urad(lev) = urad(lev) + radlu
               if (idrv .eq. 1) then
                  d_radlu_dt = d_radlu_dt * (1.0_rb - atrans(lev))
                  d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
               endif
            endif
!  Set clear sky stream to total sky stream as long as all layers
!  are clear (iclddn=0).  Streams must be calculated separately at 
!  all layers when a cloud is present (iclddn=1), because surface 
!  reflectance is different for each stream.
               if (iclddn.eq.1) then
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev) 
                  clrurad(lev) = clrurad(lev) + radclru
               else
                  radclru = radlu
                  clrurad(lev) = urad(lev)
               endif
               if (idrv .eq. 1) then
                  if (iclddn.eq.1) then
                     d_radclru_dt = d_radclru_dt * (1.0_rb - atrans(lev))
                     d_clrurad_dt(lev) = d_clrurad_dt(lev) + d_radclru_dt
                  else
                     d_radclru_dt = d_radlu_dt
                     d_clrurad_dt(lev) = d_urad_dt(lev)
                  endif
               endif
         enddo

! Increment g-point counter
         igc = igc + 1
! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs(iband)) go to 1000

! Process longwave output from band.
! Calculate upward, downward, and net flux.
         do lev = nlayers, 0, -1
            uflux(lev) = urad(lev)*wtdiff
            dflux(lev) = drad(lev)*wtdiff
            urad(lev) = 0.0_rb
            drad(lev) = 0.0_rb
! First calculate the fluxes by band.
            totuflux_bnd(lev,iband) = uflux(lev) * delwave(iband)
            totdflux_bnd(lev,iband) = dflux(lev) * delwave(iband)
! Later accumulate them.
            totuflux(lev) = totuflux(lev) + totuflux_bnd(lev,iband)
            totdflux(lev) = totdflux(lev) + totdflux_bnd(lev,iband)
            uclfl(lev) = clrurad(lev)*wtdiff
            dclfl(lev) = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_rb
            clrdrad(lev) = 0.0_rb
            totuclfl_bnd(lev,iband) = uclfl(lev) * delwave(iband)
            totdclfl_bnd(lev,iband) = dclfl(lev) * delwave(iband)
            totuclfl(lev) = totuclfl(lev) + totuclfl_bnd(lev,iband)
            totdclfl(lev) = totdclfl(lev) + totdclfl_bnd(lev,iband)
         enddo

! Calculate total change in upward flux wrt surface temperature
         if (idrv .eq. 1) then
            do lev = nlayers, 0, -1
               duflux_dt(lev) = d_urad_dt(lev) * wtdiff
               d_urad_dt(lev) = 0.0_rb
               dtotuflux_dt(lev) = dtotuflux_dt(lev) + duflux_dt(lev) * delwave(iband) * fluxfac

               duclfl_dt(lev) = d_clrurad_dt(lev) * wtdiff
               d_clrurad_dt(lev) = 0.0_rb
               dtotuclfl_dt(lev) = dtotuclfl_dt(lev) + duclfl_dt(lev) * delwave(iband) * fluxfac
            enddo
         endif

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(0) = totuflux(0) * fluxfac
      totdflux(0) = totdflux(0) * fluxfac
      fnet(0) = totuflux(0) - totdflux(0)
      totuclfl(0) = totuclfl(0) * fluxfac
      totdclfl(0) = totdclfl(0) * fluxfac
      fnetc(0) = totuclfl(0) - totdclfl(0)
! Calculate fluxes at surface by band
      do iband = 1, nbndlw
         totuflux_bnd(0,iband) = totuflux_bnd(0,iband) * fluxfac
         totdflux_bnd(0,iband) = totdflux_bnd(0,iband) * fluxfac
         fnet_bnd(0,iband) = totuflux_bnd(0,iband) - totdflux_bnd(0,iband)
         totuclfl_bnd(0,iband) = totuclfl_bnd(0,iband) * fluxfac
         totdclfl_bnd(0,iband) = totdclfl_bnd(0,iband) * fluxfac
         fnetc_bnd(0,iband) = totuclfl_bnd(0,iband) - totdclfl_bnd(0,iband)
      enddo

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(lev) = totuflux(lev) * fluxfac
         totdflux(lev) = totdflux(lev) * fluxfac
         fnet(lev) = totuflux(lev) - totdflux(lev)
         totuclfl(lev) = totuclfl(lev) * fluxfac
         totdclfl(lev) = totdclfl(lev) * fluxfac
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
! Calculate fluxes at model levels by band
         do iband = 1, nbndlw
            totuflux_bnd(lev,iband) = totuflux_bnd(lev,iband) * fluxfac
            totdflux_bnd(lev,iband) = totdflux_bnd(lev,iband) * fluxfac
            fnet_bnd(lev,iband) = totuflux_bnd(lev,iband) - totdflux_bnd(lev,iband)
            totuclfl_bnd(lev,iband) = totuclfl_bnd(lev,iband) * fluxfac
            totdclfl_bnd(lev,iband) = totdclfl_bnd(lev,iband) * fluxfac
            fnetc_bnd(lev,iband) = totuclfl_bnd(lev,iband) - totdclfl_bnd(lev,iband)
         enddo
         l = lev - 1

! Calculate heating rates at model layers
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev))
! Calculate heating rates at model layers by band
         do iband = 1, nbndlw
            htr_bnd(l,iband)=heatfac*(fnet_bnd(l,iband)-fnet_bnd(lev,iband))/(pz(l)-pz(lev)) 
            htrc_bnd(l,iband)=heatfac*(fnetc_bnd(l,iband)-fnetc_bnd(lev,iband))/(pz(l)-pz(lev))
         enddo
      enddo

! Set heating rate to zero in top layer
      htr(nlayers) = 0.0_rb
      htrc(nlayers) = 0.0_rb
! Set heating rate to zero in top layer by band
      do iband = 1, nbndlw
         htr_bnd(nlayers,iband) = 0.0_rb
         htrc_bnd(nlayers,iband) = 0.0_rb
      enddo

      end subroutine rtrnmr

      end module rrtmg_lw_rtrnmr

