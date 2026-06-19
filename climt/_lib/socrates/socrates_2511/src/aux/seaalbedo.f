! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE seaalbedo(zenith, wavelength, albedo)
      USE realtype_rd, ONLY: RealK
      USE spline_evaluate_mod, ONLY: spline_evaluate
      IMPLICIT NONE
      
      ! Arguments with intent(in)
      !
      !  zenith     -- solar zenith angle (radians)
      !  wavelength -- wavelength (m)
      !
      ! Argument with intent(out)
      !
      !  albedo     -- ocean albedo
      !
      REAL(RealK) zenith
      REAL(RealK) wavelength
      REAL(RealK) albedo

      REAL(RealK) n_pure_water
      REAL n_sea_water
      REAL salinity, temperature, psi
      REAL vent10, vent41feet, ww, pi
      PARAMETER (pi=3.1415926)
      REAL(RealK) rmorel
      REAL rfresnel, rfres_plat, rmousse
      PARAMETER (salinity=35.0)
      PARAMETER (temperature=20.0) 
      PARAMETER (vent10=7.0)

      psi = SNGL(zenith)
      if (psi .lt. 0.01) psi = 0.01

C     prevent going out of range
      if (wavelength .gt. 5.e-6)  then 
         wavelength = 5.e-6
      ENDIF
      if (wavelength .lt. 0.2e-6) then 
         wavelength = 0.2e-6
      ENDIF

      ! convert to microns
      wavelength = wavelength * 1.e6
 
C-donne l'indice de refraction de l'eau pure a la l.o. wavelength a 25oC
      CALL m_pure_water(REAL(wavelength,RealK), n_pure_water)
C-donne l'indice de refraction de l'eau de mer. 
      CALL m_sea_water(SNGL(n_pure_water), salinity, 
     &                 temperature, n_sea_water)
C-donne le vent aux bonnes altitudes
      CALL vent_z(vent10, 10.0, 12.464, vent41feet)
C-cacule la reflection de Fresnel 
      CALL fresnel(psi, n_sea_water, vent41feet, rfresnel)
      CALL fresnel_plat(psi, n_sea_water, rfres_plat)
      IF (rfresnel.GT.1.0) rfresnel = rfres_plat
      rfresnel=MAX(0.02, rfresnel)
C-calcule la surface des whitecaps et leur albedo
      CALL mousse(wavelength, vent10, ww, rmousse)
C-calcule l'albedo de ce qui remonte 
      CALL morel(wavelength, rmorel)
C-somme les contributions
      albedo = (1.0d0-DBLE(ww))*
     &      (DBLE(rfresnel+rmorel))+DBLE(ww)*DBLE(rmousse)
      albedo = MIN(albedo, 1.0d0)
      
      RETURN
      END
C=================================================================
      SUBROUTINE vent_z(vent1, z1, z2, vent2)
      IMPLICIT NONE
      REAL vent2, vent1, z1, z2
      REAL u_star, ct_karman
      PARAMETER (ct_karman=0.35)
      u_star=0.03*vent1
      vent2=vent1+u_star/ct_karman*log(z2/z1)
      vent2=MAX(0.0, vent2)
      RETURN
      END
C---------------------------------------------------------
      SUBROUTINE morel(wv,rmorel)
      USE realtype_rd, ONLY: RealK
      USE spline_evaluate_mod, ONLY: spline_evaluate
      USE spline_fit_mod, ONLY: spline_fit
      IMPLICIT NONE
      INTEGER ierr
      INTEGER wvmax
      PARAMETER (wvmax=15)
      REAL(RealK) wvi(wvmax), rmoreli(wvmax), rmoreli2(wvmax)
      REAL(RealK) wv, rmorel 
      DATA wvi/0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
     .         0.7, 0.8, 0.9, 1.0,  2.0,  3.0, 4.0/
      DATA rmoreli/0.10, 0.10, 0.10, 0.0898, 0.0424, 0.0137, 
     .             0.0034, 0.0018, 0.0010, 0.0001, 0.0001, 
     .             0.0001, 0.0001, 0.0001, 0.0001/

      CALL spline_fit(wvmax, wvi, rmoreli, rmoreli2)
      CALL spline_evaluate(ierr,wvmax,wvi,rmoreli,rmoreli2,wv,rmorel)
      rmorel=0.52*rmorel
      RETURN 
      END
c----------------------------------------------------------
      SUBROUTINE m_pure_water(wv,m)
      USE realtype_rd, ONLY: RealK
      USE spline_evaluate_mod, ONLY: spline_evaluate
      USE spline_fit_mod, ONLY: spline_fit
      IMPLICIT NONE
C----------------------------------------------------------
C--REFERENCE = Jerlov, Marine Optics, 1976
C            = Hale and Querry, 1973
C--spline interpolation to pure water refractive index as a 
C--function of wavelength
C-----------------------------------------------------------
      INTEGER ierr
      INTEGER wvmax
c-------Jerlov
c      PARAMETER (wvmax=10) 
c-------Hale and Querry 1973
      PARAMETER (wvmax=65)
      REAL(RealK) mi(wvmax), wvi(wvmax)
      REAL(RealK) y2(wvmax)
c-------Jerlov
c      DATA mi/1.3773, 1.3569, 1.3480, 1.3433, 1.3403, 
c     .        1.3371, 1.3330, 1.3289, 1.3247, 1.3210/
c      DATA wvi/0.250, 0.308, 0.359, 0.400, 0.434, 
c     .         0.486, 0.589, 0.768, 1.000, 1.250/
c-------Hale and Querry 1973
      DATA mi/1.396, 1.373, 1.362, 1.354, 1.349, 1.346, 1.343, 
     .        1.341, 1.339, 1.338, 1.337, 1.336, 1.335, 1.334, 
     .        1.333, 1.333, 1.332, 1.332, 1.331, 1.331, 1.331, 
     .        1.330, 1.330, 1.330, 1.329, 1.329, 1.329, 1.328, 
     .        1.328, 1.328, 1.327, 1.327, 1.327, 1.324, 1.321, 
     .        1.317, 1.312, 1.306, 1.296, 1.279, 1.242, 1.219, 
     .        1.188, 1.157, 1.142, 1.149, 1.201, 1.292, 1.371,
     .        1.426, 1.467, 1.483, 1.478, 1.467, 1.450, 1.432, 
     .        1.420, 1.410, 1.400, 1.385, 1.374, 1.364, 1.357, 
     .        1.351, 1.346/ 
      DATA wvi/0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 
     .         0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 
     .         0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 
     .         0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 
     .         0.900, 0.925, 0.950, 0.975, 1.000, 1.200, 1.400,
     .         1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.650, 
     .         2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000,
     .         3.050, 3.100, 3.150, 3.200, 3.250, 3.300, 3.350, 
     .         3.400, 3.450, 3.500, 3.600, 3.700, 3.800, 3.900,
     .         4.000, 4.100/
      REAL(RealK) m, wv

      CALL spline_fit(wvmax, wvi, mi, y2)
      CALL spline_evaluate(ierr, wvmax, wvi, mi, y2, wv, m)
      RETURN
      END
C-----------------------------------------------------------------

      SUBROUTINE mousse(wavelength,vent10,ww,reflectance)
      USE realtype_rd, ONLY: RealK
      USE spline_evaluate_mod, ONLY: spline_evaluate
      USE spline_fit_mod, ONLY: spline_fit
      IMPLICIT NONE
C----------------------------------------------------
C--REFERENCE : Koepke 1984
C--            Whitlock et al.    1982
C--INPUTS 
C--wavelength in um
C--vent10 : 10 meter wind
C--OUTPUTS 
C--ww : surface covered by whitecaps
C--reflectance : reflectance of whitecaps
C----------------------------------------------------
      INTEGER ierr
      REAL vent10, ww 
      INTEGER nwvmax
      PARAMETER (nwvmax=37)
      REAL(RealK) wv(1:nwvmax), whitlock(1:nwvmax)
      REAL(RealK) whitlock2(1:nwvmax)
      DATA wv/0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
     .        1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
     .        1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 
     .        2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1,
     .        3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 
     .        3.9, 4.0/
      DATA whitlock/0.54, 0.54, 0.54, 0.54, 0.53, 0.51, 0.485,
     .              0.46, 0.41, 0.34, 0.24, 0.18, 0.23, 0.26, 
     .              0.24, 0.14, 0.11, 0.13, 0.15, 0.14, 0.11,
     .              0.08, 0.06, 0.01, 0.00, 0.00, 0.00, 0.00, 
     .              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     .              0.00, 0.00/
      REAL(RealK) wavelength, reflec
      REAL xx, sp, f_ef_p, f_ef_s
      PARAMETER (xx=0.5, f_ef_p=0.4, f_ef_s=0.18)
      REAL wp, ws, f_ef, reflectance

      CALL spline_fit(nwvmax, wv, whitlock, whitlock2)
      sp=0.
      IF (vent10.GT.8.75) sp=-2.1 + 0.24 * vent10
      ww=2.95E-6*(vent10**3.52)
      wp=1./(1.+xx*sp)
      ws=sp/(1.+xx*sp)
      f_ef=wp*f_ef_p + ws*f_ef_s
      CALL spline_evaluate(ierr, nwvmax, wv, whitlock, whitlock2,
     .                     wavelength, reflec)
      reflec=MAX(0.0_RealK, reflec)
      reflectance=f_ef*reflec
      ww=ww*(1.+sp)/(1.+xx*sp)
      RETURN
      END
C-------------------------------------------------------------
      SUBROUTINE m_sea_water(n_pure_water, salinity,
     .                       temperature,n_sea_water)
      IMPLICIT NONE 
C------------------------------------------------------
C--REFERENCE Stelson 1990
C------------Jerlov 1976
c------------Ivanoff
C--INPUTS
C--salinity parts per thousands
C--temperature in Celsius
C--n_pure_water : reference pure water refractive index at 
C--               the specific wavelength considered
C--OUTPUTS 
C--n sea water refractive index
C------------------------------------------------------
      REAL r_nacl, r_w, m_nacl, m_w
      PARAMETER (m_nacl=58.5, m_w=18.01534) !--g mol-1
c--r_nacl est fitte sur Jerlov a 20oC
      PARAMETER (r_nacl=8.87)               !--cm3/mol-1
      REAL d, d_0, salinity, temperature, v, xr
      REAL x_w, x_nacl, mol_nacl, mol_w
      REAL n_pure_water, v_w_pure, n_sea_water
c      REAL sigma_0
    
c---determination de r_w pour fitter n_pure_water
      CALL eos(temperature,0.,d_0)
      d_0=d_0/1000.
      v_w_pure=m_w/d_0
      r_w=v_w_pure*(n_pure_water*n_pure_water-1.)/
     .             (n_pure_water*n_pure_water+2)
c--formule Ivanoff---------
c      sigma_0=-0.093+salinity*
c     .        (0.8149+salinity*(-0.000482+0.0000068*salinity))
c      d=(sigma_0/1000.+1.)*d_0           !--density g cm-3
c--formule unesco (lodyc)--
      CALL eos(temperature,salinity,d)
      d=d/1000.
c--------------------------
      mol_w=1000./m_w                    !--mol
      mol_nacl=salinity/m_nacl           !--mol
      x_w=mol_w/(mol_w+mol_nacl)         !--mol/mol
      x_nacl=mol_nacl/(mol_w+mol_nacl)   !--mol/mol
      v=(x_nacl*m_nacl+x_w*m_w)/d        !--cm3 mol-1
      xr=x_nacl*r_nacl+x_w*r_w
      n_sea_water=SQRT( (1.+2.*xr/v)/(1.-xr/v) )
      RETURN
      END
C--------------------------------------------------------------
      SUBROUTINE eos(temperature, salinity, WKZ)
      IMPLICIT NONE
      REAL :: temperature, salinity, WKZ
CC  EOS80 IS VALID FOR S= 0 TO 42 PSU ; T= -2 TO 40 CELCIUS ; P= 0
CC  TO 10000 DECIBARS. CHECK VALUE: RHO = 1059.8204 KG/M**3 FOR
CC  P=10000 DBAR, T = 40 DEG CELCIUS, S=40 PSU
      REAL :: WKX, WKY
      REAL :: ZA, ZAW, ZA1, ZB,ZBW, ZB1, ZC, ZD, ZE, ZH, ZKW, ZK0, ZP
      REAL :: ZQ, ZR1, ZR2, ZR3, ZR4, ZS, ZSR, ZT, ZXK

      ZP = 0.0
      ZH = 0.0
      ZT = temperature
      ZS = salinity
      ZXK= ZH*FSATG(ZS,ZT,ZP)
      ZT = ZT + 0.5*ZXK
      ZQ = ZXK
      ZP = ZP + 0.5*ZH
      ZXK= ZH*FSATG(ZS,ZT,ZP)
      ZT = ZT + 0.29289322*(ZXK-ZQ)
      ZQ = 0.58578644*ZXK + 0.121320344*ZQ
      ZXK= ZH*FSATG(ZS,ZT,ZP)
      ZT = ZT + 1.707106781*(ZXK-ZQ)
      ZQ = 3.414213562*ZXK - 4.121320344*ZQ
      ZP = ZP + 0.5*ZH
      ZXK= ZH*FSATG(ZS,ZT,ZP)
      WKX= ZT + (ZXK-2.0*ZQ)/6.0
C
      ZT = WKX
C ...   square root salinity
      ZSR= SQRT(ABS(ZS))
      WKY=ZSR
C ...   compute density pure water at atm pressure
      ZR1= ((((6.536332E-9*ZT-1.120083E-6)*ZT+1.001685E-4)*ZT
     $              -9.095290E-3)*ZT+6.793952E-2)*ZT+999.842594
C ...   seawater density atm pressure
      ZR2= (((5.3875E-9*ZT-8.2467E-7)*ZT+7.6438E-5)*ZT
     $             -4.0899E-3)*ZT+8.24493E-1
      ZR3= (-1.6546E-6*ZT+1.0227E-4)*ZT-5.72466E-3
      ZR4= 4.8314E-4
      WKZ= (ZR4*ZS + ZR3*ZSR + ZR2)*ZS + ZR1
      ZT = WKX
      ZSR= WKY
      ZH = 0.0
C
      ZE = (9.1697E-11*ZT+2.0816E-9)*ZT-9.9348E-8
      ZBW= (5.2787E-9*ZT-6.12293E-7)*ZT+8.50935E-6
      ZB = ZBW + ZE * ZS
C
      ZD = 1.91075E-4
      ZC = (-1.6078E-6*ZT-1.0981E-5)*ZT+2.2838E-3
      ZAW= ((-5.77905E-7*ZT+1.16092E-4)*ZT+1.43713E-3)*ZT
     $              +3.239908
      ZA = ( ZD*ZSR + ZC)*ZS + ZAW
C
      ZB1= (-5.3009E-3*ZT+1.6483E-1)*ZT+7.944E-1
      ZA1= ((-6.1670E-4*ZT+1.09987E-1)*ZT-6.03459)*ZT+546.746
      ZKW= (((-5.155288E-4*ZT+1.360477E-1)*ZT-23.27105)*ZT
     $                +1484.206)*ZT+196522.1
      ZK0= (ZB1*ZSR + ZA1)*ZS + ZKW
C ...   evaluate pressure polynomial
      WKZ=WKZ/(1.0-ZH/(ZK0+ZH*(ZA+ZB*ZH)))

      CONTAINS
        REAL FUNCTION FSATG(PS,PT,PHP)
          IMPLICIT NONE
          REAL :: PS,PT,PHP
          FSATG = (((-2.1687E-16*PT+1.8676E-14)*PT-4.6206E-13)*PHP
     &      +((2.7759E-12*PT-1.1351E-10)*(PS-35.)+((-5.4481E-14*PT
     &      +8.733E-12)*PT-6.7795E-10)*PT+1.8741E-8))*PHP
     &      +(-4.2393E-8*PT+1.8932E-6)*(PS-35.)
     &      +((6.6228E-10*PT-6.836E-8)*PT+8.5258E-6)*PT+3.5803E-5
        END FUNCTION
      END
C--------------------------------------------------------
      SUBROUTINE fresnel(psi, m, W, rfresnel)
      IMPLICIT NONE
C--------------------------------------------------------
C---REFERENCE = Cox and Munk (1956)
C               Cox and Munk (1954)
C---------------------------------------------------------
C---INPUTS
C---psi = zenith angle (in radians)
C---W = wind at 41 feet
C---m = water refractive index
C---OUTPUTS
C---rfresnel = Fresnel albedo for wind W
C---LOCAL
C---psir = angle of refracted beam (radians)
C---sigma2 = standard deviation of the slope distribution
C---rho_psi = fresnel's albedo for flat sea
C---a, b, c = coefficients (Cox and Munk 1956)
C---gf_psi = F(psi) = rho(psi) * cos(psi)
C---gfprim = derivative of gf_psi
C---gfsec = second derivative of gf_psi
C---erf = error function
C---ik = erf(k)
C---------------------------------------------------------
      REAL psi, psir, rho_psi, gf_psi, f_psi
      REAL gfprim_psi, gfprim_eps, gfsec_psi
      REAL rho_eps, rho_eps2, gf_eps, gf_eps2
      REAL sigma, sigma2, m, W, pi
      PARAMETER (pi=3.1415927)
      REAL a, b, c, k, rmin, rmax, rfresnel
      REAL eps, epsr, eps2, epsr2, epsilon, ik 
      PARAMETER (epsilon=0.001)
      REAL psi_c, rfmax

      sigma2 = 0.003 + 0.00512*W
      sigma=SQRT(sigma2)

      psi_c=psi
      rfmax=0.0
 100  psir=ASIN(SIN(psi_c)/m)
      rho_psi=0.5*( (SIN(psi_c-psir)/SIN(psi_c+psir))**2.
     .             +(TAN(psi_c-psir)/TAN(psi_c+psir))**2.) 
      eps=psi_c+epsilon
      epsr=ASIN(SIN(eps)/m)
      rho_eps=0.5*( (SIN(eps-epsr)/SIN(eps+epsr))**2.
     .             +(TAN(eps-epsr)/TAN(eps+epsr))**2.)
      eps2=eps+epsilon
      epsr2=ASIN(SIN(eps2)/m)
      rho_eps2=0.5*( (SIN(eps2-epsr2)/SIN(eps2+epsr2))**2.
     .              +(TAN(eps2-epsr2)/TAN(eps2+epsr2))**2.)
      gf_psi=rho_psi*COS(psi_c) 
      gf_eps=rho_eps*COS(eps)
      gf_eps2=rho_eps2*COS(eps2)
      gfprim_psi=(gf_eps-gf_psi)/epsilon 
      gfprim_eps=(gf_eps2-gf_eps)/epsilon
      gfsec_psi=(gfprim_eps-gfprim_psi)/epsilon 
      a=-gfprim_psi/gf_psi
      b=0.5+0.5*gfsec_psi/gf_psi
      c=0.5+0.5*gfprim_psi/gf_psi/TAN(psi_c)
      f_psi=0.5*(b+c)
      k=1./sigma/TAN(psi_c)
      ik=erf(k)
      rmin=rho_psi*(0.5*(1.+ik)+0.5/SQRT(pi)*a*sigma*EXP(-k*k)+
     .              0.25*b*sigma2*(1.+ik-2./SQRT(pi)*k*EXP(-k*k))+
     .              0.25*c*sigma2*(1.+ik) ) 
      k=0.5/sigma/TAN(psi_c)
      ik=erf(k)
      rmax=rho_psi*(0.5*(1.+ik)+0.5/SQRT(pi)*a*sigma*EXP(-k*k)+
     .              0.25*b*sigma2*(1.+ik-2./SQRT(pi)*k*EXP(-k*k))+
     .              0.25*c*sigma2*(1.+ik) )
      rfresnel=(rmin+rmax)/2. 
      rfresnel=MAX(0.0,rfresnel)
      rfresnel=MIN(1.0,rfresnel)
      IF (rfmax.LE.rfresnel) THEN 
        rfmax=rfresnel 
        psi_c=psi_c-0.001
        IF (psi_c.LT.0.) GOTO 200
        GOTO 100
      ENDIF
 200  rfresnel=rfmax
      RETURN
      END
c-----------------------------------------------------------
      SUBROUTINE fresnel_plat(psi, m, rfresnel_plat)
      IMPLICIT NONE
C--------------------------------------------------------
C---INPUTS
C---psi = zenith angle (in radians)
C---m = water refractive index
C---OUTPUTS
C---rfresnel = Fresnel albedo for flat sea
C---LOCAL
C---psir = angle of refracted beam (radians)
C---------------------------------------------------------
      REAL psi, psir
      REAL m, rfresnel_plat

      psir=ASIN(SIN(psi)/m)
      rfresnel_plat=0.5*( (SIN(psi-psir)/SIN(psi+psir))**2.
     .                   +(TAN(psi-psir)/TAN(psi+psir))**2.)
      RETURN
      END
