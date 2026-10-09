! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module with data and functions to calculate real refractive indices
!  of gases.
!
! Description:
!   Formulas for the real refractive index is taken from
!   refractiveindex.info, sources for individual gases are provided
!   below.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
MODULE refract_re_ccf

USE realtype_rd, ONLY: RealK
USE gas_list_pcf

IMPLICIT NONE

! Number of gases refractive index data
INTEGER, PARAMETER :: n_ri_gases = 8

! Derived type to keep validity parameters for dispersion formulas
TYPE StrRefractParam
  INTEGER     :: type_gas
!   Gas types (as defined in gas_list_pcf)
  REAL(RealK) :: wavelength_min
!   Minimum wavelength at which dispersion formula is valid
  REAL(RealK) :: wavelength_max
!   Maximum wavelength at which dispersion formula is valid
  REAL(RealK) :: pressure
!   Pressure at which dispersion formula is valid
  REAL(RealK) :: temperature
!   Temperature at which dispersion formula is valid
END TYPE StrRefractParam

! Define dispersion formula validity parameters for each gas
TYPE(StrRefractParam), PARAMETER :: RefractParam(n_ri_gases) = (/ &
  StrRefractParam( &
    type_gas        = ip_n2, &
    wavelength_min  = 0.4679E-06_RealK, &
    wavelength_max  = 2.0587E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 273.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_o2, &
    wavelength_min  = 0.4000E-06_RealK, &
    wavelength_max  = 1.8000E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 293.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_h2o, &
    wavelength_min  = 0.3500E-06_RealK, &
    wavelength_max  = 1.2000E-06_RealK, &
    pressure        = 1333.0_RealK, &
    temperature     = 293.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_ar, &
    wavelength_min  = 0.1404E-06_RealK, &
    wavelength_max  = 0.5677E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 273.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_co2, &
    wavelength_min  = 0.1807E-06_RealK, &
    wavelength_max  = 1.6945E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 273.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_ch4, &
    wavelength_min  = 0.5290E-06_RealK, &
    wavelength_max  = 0.6585E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 293.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_h2, &
    wavelength_min  = 0.1680E-06_RealK, &
    wavelength_max  = 1.6945E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 273.15_RealK), &

  StrRefractParam( &
    type_gas        = ip_he, &
    wavelength_min  = 0.2753E-06_RealK, &
    wavelength_max  = 0.5462E-06_RealK, &
    pressure        = 101325.0_RealK, &
    temperature     = 273.15_RealK) &
  /)

! Standard atmospheric conditions
REAL(RealK), PARAMETER :: &
    p_stp = 101325.0_RealK, &
!     Standard pressure
    t_stp = 273.15_RealK
!     Standard temperature

CONTAINS

! Refractive index of gas less 1 at STP
REAL(RealK) FUNCTION refract_re_m1(i_gas, wavelength)

! Input variables
  INTEGER, INTENT(IN) :: i_gas
!   Type index of gas
  REAL(RealK), INTENT(IN) :: wavelength
!   Wavelengths at which to calculate refractive index

! Local variables
  INTEGER :: i_ri_gas
!   Local index of gas
  INTEGER :: i
!   Loop index

! Find local index of gas
  i_ri_gas = -1
  DO i = 1, n_ri_gases
    IF (i_gas == RefractParam(i)%type_gas) THEN
      i_ri_gas = i
      EXIT
    END IF
  END DO

! If exact refractive index is not available set it to unity
  IF (i_ri_gas == -1) THEN
    refract_re_m1 = 0.0_RealK
    RETURN
  END IF

  IF (wavelength > RefractParam(i_ri_gas)%wavelength_max) THEN
    refract_re_m1 = refract_re_m1_gas(i_gas, &
        RefractParam(i_ri_gas)%wavelength_max)
  ELSE IF (wavelength < RefractParam(i_ri_gas)%wavelength_min) THEN
    refract_re_m1 = refract_re_m1_gas(i_gas, &
        RefractParam(i_ri_gas)%wavelength_min)
  ELSE
    refract_re_m1 = refract_re_m1_gas(i_gas, wavelength)
  END IF

! Convert refractive index to STP
  refract_re_m1 = refract_re_m1* &
      p_stp*RefractParam(i_ri_gas)%temperature/ &
      (RefractParam(i_ri_gas)%pressure*t_stp)

END FUNCTION refract_re_m1

! Refractive index less 1 at measured P,T
REAL(RealK) FUNCTION refract_re_m1_gas(i_gas, wavelength)

! Input variables
  INTEGER, INTENT(IN) :: i_gas
!   Type index of gas
  REAL(RealK), INTENT(IN) :: wavelength
!   Wavelength at which to calculate refractive index

! Local variables
  REAL(RealK) :: wl_um
!   Wavelength in micro metre

  wl_um = wavelength*1.0E+06_RealK

  SELECT CASE (i_gas)

  CASE (ip_n2)
!   E. R. Peck and B. N. Khanna. J. Opt. Soc. Am. 56, 1059-1063 (1963)
    refract_re_m1_gas = 6.8552E-05_RealK + &
        3.243157E-02_RealK/(144.0_RealK - wl_um**(-2))

  CASE (ip_o2)
!   J. Zhang, Z. H. Lu, and L. J. Wang. Appl. Opt. 47, 3143-3151 (2008)
!   P. Křen. Appl. Opt. 50, 6484-6485 (2011)
    refract_re_m1_gas = 1.181494E-04_RealK + &
        9.708931E-03_RealK/(75.4_RealK - wl_um**(-2))

  CASE (ip_h2o)
!   P. E. Ciddor. Appl. Opt. 35, 1566-1573 (1996)
    refract_re_m1_gas = 1.022E-08_RealK*(295.235_RealK + &
        2.6422_RealK*wl_um**(-2) - &
        0.032380_RealK*wl_um**(-4) + &
        0.004028_RealK*wl_um**(-6))

  CASE (ip_ar)
!   A. Bideau-Mehu, Y. Guern, R. Abjean, A. Johannin-Gilles. J. Quant.
!   Spectrosc. Rad. Transfer 25, 395-402 (1981)
    refract_re_m1_gas = &
        2.50141E-03_RealK/(91.012_RealK - wl_um**(-2)) + &
        5.00283E-04_RealK/(87.892_RealK - wl_um**(-2)) + &
        5.22343E-02_RealK/(214.02_RealK - wl_um**(-2))

  CASE (ip_co2)
!   A. Bideau-Mehu, Y. Guern, R. Abjean and A. Johannin-Gilles. Opt. Commun. 9,
!   432-434 (1973)
    refract_re_m1_gas = &
        6.99100E-02_RealK/(166.175_RealK - wl_um**(-2)) + &
        1.44720E-03_RealK/(79.609_RealK - wl_um**(-2)) + &
        6.42941E-05_RealK/(56.3064_RealK - wl_um**(-2)) + &
        5.21306E-05_RealK/(46.0196_RealK - wl_um**(-2)) + &
        1.46847E-06_RealK/(0.0584738_RealK - wl_um**(-2))

  CASE (ip_ch4)
!   S. Loria. Annalen der Physik 334, 605-622 (1909)
    refract_re_m1_gas = 0.00042607_RealK + 6.1396687E-06_RealK/(wl_um**2)

  CASE (ip_h2)
!   E. R. Peck and S. Hung. J. Opt. Soc. Am. 67, 1550-1554 (1977)
    refract_re_m1_gas = &
        0.0148956_RealK/(180.7_RealK - wl_um**(-2)) + &
        0.0049037_RealK/(92.0_RealK - wl_um**(-2))

  CASE (ip_he)
!   C. Cuthbertson and M. Cuthbertson. Proc. R. Soc. London A 135, 40-47 (1936)
!   C. Cuthbertson and M. Cuthbertson. Proc. R. Soc. London A 83, 149-151 (1910)
    refract_re_m1_gas = &
        0.014755297_RealK/(426.2974_RealK - wl_um**(-2))

  CASE DEFAULT
!   No data for gas, refractive index is set to unity
    refract_re_m1_gas = 0.0_RealK

  END SELECT

END FUNCTION refract_re_m1_gas

END MODULE refract_re_ccf