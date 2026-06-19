! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module to set indexing numbers of gaseous absorbing species.
!
! Description:
!   This module defines the identifiers defining the physical types
!   of each molecular absorbing species.
!   The numbering 1-12 agrees with HITRAN.
!
!- ---------------------------------------------------------------------
MODULE gas_list_pcf

USE realtype_rd, ONLY: RealK

IMPLICIT NONE

INTEGER, PRIVATE :: i

INTEGER, PARAMETER :: npd_gases = 75
!   Number of indexed gases

INTEGER, PARAMETER :: IP_h2o = 1
!   Identifier for water vapour
INTEGER, PARAMETER :: IP_co2 = 2
!   Identifier for carbon dioxide
INTEGER, PARAMETER :: IP_o3 = 3
!   Identifier for ozone
INTEGER, PARAMETER :: IP_n2o = 4
!   Identifier for dinitrogen oxide
INTEGER, PARAMETER :: IP_co = 5
!   Identifier for carbon monoxide
INTEGER, PARAMETER :: IP_ch4 = 6
!   Identifier for methane
INTEGER, PARAMETER :: IP_o2 = 7
!   Identifier for oxygen
INTEGER, PARAMETER :: IP_no = 8
!   Identifier for nitrogen monoxide
INTEGER, PARAMETER :: IP_so2 = 9
!   Identifier for sulphur dioxide
INTEGER, PARAMETER :: IP_no2 = 10
!   Identifier for nitrogen dioxide
INTEGER, PARAMETER :: IP_nh3 = 11
!   Identifier for ammonia
INTEGER, PARAMETER :: IP_hno3 = 12
!   Identifier for nitric acid
INTEGER, PARAMETER :: IP_n2 = 13
!   Identifier for nitrogen
INTEGER, PARAMETER :: IP_cfc11 = 14
!   Identifier for CFC11 (CFCl3)
INTEGER, PARAMETER :: IP_cfc12 = 15
!   Identifier for CFC12 (CF2Cl2)
INTEGER, PARAMETER :: IP_cfc113 = 16
!   Identifier for CFC113 (CF2ClCFCl2)
INTEGER, PARAMETER :: IP_hcfc22 = 17
!   Identifier for HCFC22 (CHF2Cl)
INTEGER, PARAMETER :: IP_hfc125 = 18
!   Identifier for HFC125 (C2HF5)
INTEGER, PARAMETER :: IP_hfc134a = 19
!   Identifier for HFC134A (CF3CFH2)
INTEGER, PARAMETER :: IP_cfc114 = 20
!   Identifier for CFC114 (C2Cl2F4)
INTEGER, PARAMETER :: IP_tio = 21
!   Identifier for TiO
INTEGER, PARAMETER :: IP_vo = 22
!   Identifier for VO
INTEGER, PARAMETER :: IP_h2 = 23
!   Identifier for hydrogen
INTEGER, PARAMETER :: IP_he = 24
!   Identifier for helium
INTEGER, PARAMETER :: IP_ocs = 25
!   Identifier for carbonyl sulphide
INTEGER, PARAMETER :: IP_na = 26
!   Identifier for sodium
INTEGER, PARAMETER :: IP_k = 27
!   Identifier for potassium
INTEGER, PARAMETER :: IP_feh = 28
!   Identifier for iron hydride
INTEGER, PARAMETER :: IP_crh = 29
!   Identifier for chromium hydride
INTEGER, PARAMETER :: IP_li = 30
!   Identifier for lithium
INTEGER, PARAMETER :: IP_rb = 31
!   Identifier for rubidium
INTEGER, PARAMETER :: IP_cs = 32
!   Identifier for cesium
INTEGER, PARAMETER :: IP_ph3 = 33
!   Identifier for phosphine
INTEGER, PARAMETER :: IP_c2h2 = 34
!   Identifier for acetylene
INTEGER, PARAMETER :: IP_hcn = 35
!   Identifier for hydrogen cyanide
INTEGER, PARAMETER :: IP_h2s = 36
!   Identifier for hydrogen sulphide
INTEGER, PARAMETER :: IP_ar = 37
!   Identifier for argon
INTEGER, PARAMETER :: IP_air = 38
!   Identifier for all other gases, used by generalised continuum
INTEGER, PARAMETER :: IP_o = 39
!   Identifier for atomic oxygen
INTEGER, PARAMETER :: IP_n = 40
!   Identifier for atomic nitrogen
INTEGER, PARAMETER :: IP_no3  = 41
!   Identifier for nitrate radical
INTEGER, PARAMETER :: IP_n2o5 = 42
!   Identifier for dinitrogen pentoxide
INTEGER, PARAMETER :: IP_hono = 43
!   Identifier for nitrous acid
INTEGER, PARAMETER :: IP_ho2no2 = 44
!   Identifier for peroxynitric acid
INTEGER, PARAMETER :: IP_h2o2 = 45
!   Identifier for hydrogen peroxide
INTEGER, PARAMETER :: IP_c2h6 = 46
!   Identifier for ethane
INTEGER, PARAMETER :: IP_ch3 = 47
!   Identifier for methyl radical
INTEGER, PARAMETER :: IP_h2co = 48
!   Identifier for formaldehyde
INTEGER, PARAMETER :: IP_ho2 = 49
!   Identifier for hydroperoxy radical
INTEGER, PARAMETER :: IP_hdo = 50
!   Identifier for semiheavy water
INTEGER, PARAMETER :: IP_hcl = 51
!   Identifier for hydrogen chloride
INTEGER, PARAMETER :: IP_hf = 52
!   Identifier for hydrogen fluoride
INTEGER, PARAMETER :: IP_cosso = 53
!   Identifier for cis-OSSO
INTEGER, PARAMETER :: IP_tosso = 54
!   Identifier for trans-OSSO
INTEGER, PARAMETER :: IP_yosos = 55
!   Identifier for OSO-S
INTEGER, PARAMETER :: IP_ch3cho = 56
!   Identifier for acetaldehyde
INTEGER, PARAMETER :: IP_ch3ooh = 57
!   Identifier for methylhydroperoxide
INTEGER, PARAMETER :: IP_ch3coch3 = 58
!   Identifier for acetone
INTEGER, PARAMETER :: IP_ch3cocho = 59
!   Identifier for methylglyoxal (MGLY)
INTEGER, PARAMETER :: IP_chocho = 60
!   Identifier for glyoxal
INTEGER, PARAMETER :: IP_c2h5cho = 61
!   Identifier for propanal
INTEGER, PARAMETER :: IP_hoch2cho = 62
!   Identifier for glycolaldehyde
INTEGER, PARAMETER :: IP_c2h5coch3 = 63
!   Identifier for methyl ethyl ketone (MEK)
INTEGER, PARAMETER :: IP_mvk = 64
!   Identifier for methyl vinyl ketone (MVK)
INTEGER, PARAMETER :: IP_macr = 65
!   Identifier for methacrolein (MACR)
INTEGER, PARAMETER :: IP_pan = 66
!   Identifier for peroxyacetyl nitrate (PAN)
INTEGER, PARAMETER :: IP_ch3ono2 = 67
!   Identifier for methylnitrate
INTEGER, PARAMETER :: IP_c2h3 = 68
!   Identifier for vinyl radical
INTEGER, PARAMETER :: IP_c2h4 = 69
!   Identifier for ethylene
INTEGER, PARAMETER :: IP_oh = 70
!   Identifier for hydroxyl radical
INTEGER, PARAMETER :: IP_hco = 71
!   Identifier for formyl radical
INTEGER, PARAMETER :: IP_n2o4 = 72
!   Identifier for dinitrogen tetroxide
INTEGER, PARAMETER :: IP_c2n2 = 73
!   Identifier for cyanogen
INTEGER, PARAMETER :: IP_n2h4 = 74
!   Identifier for hydrazine
INTEGER, PARAMETER :: IP_n2o3 = 75
!   Identifier for dinitrogen trioxide


! Column headers for reading data in raw_input
CHARACTER (LEN=10), PARAMETER :: header_gas(npd_gases) = (/ &
    'H2O       ', 'CO2       ', 'O3        ', 'N2O       ', &
    'CO        ', 'CH4       ', 'O2        ', 'NO        ', &
    'SO2       ', 'NO2       ', 'NH3       ', 'HNO3      ', &
    'N2        ', 'CFC11     ', 'CFC12     ', 'CFC113    ', &
    'HCFC22    ', 'HFC125    ', 'HFC134A   ', 'CFC114    ', &
    'TiO       ', 'VO        ', 'H2        ', 'He        ', &
    'OCS       ', 'Na        ', 'K         ', 'FeH       ', &
    'CrH       ', 'Li        ', 'Rb        ', 'Cs        ', &
    'PH3       ', 'C2H2      ', 'HCN       ', 'H2S       ', &
    'Ar        ', 'AIR       ', 'O         ', 'N         ', &
    'NO3       ', 'N2O5      ', 'HONO      ', 'HO2NO2    ', &
    'H2O2      ', 'C2H6      ', 'CH3       ', 'H2CO      ', &
    'HO2       ', 'HDO       ', 'HCl       ', 'HF        ', &
    'cOSSO     ', 'tOSSO     ', 'yOSOS     ', 'CH3CHO    ', &
    'CH3OOH    ', 'CH3COCH3  ', 'CH3COCHO  ', 'CHOCHO    ', &
    'C2H5CHO   ', 'HOCH2CHO  ', 'C2H5COCH3 ', 'MVK       ', &
    'MACR      ', 'PAN       ', 'CH3ONO2   ', 'C2H3      ', &
    'C2H4      ', 'OH        ', 'HCO       ', 'N2O4      ', &
    'C2N2      ', 'N2H4      ', 'N2O3      '/)

! File suffixes
CHARACTER (LEN=12), PARAMETER :: gas_suffix(npd_gases) = (/ &
    'q           ', 'co2         ', 'o3          ', 'n2o         ', &
    'co          ', 'ch4         ', 'o2          ', 'no          ', &
    'so2         ', 'no2         ', 'nh3         ', 'hno3        ', &
    'n2          ', 'cfc11       ', 'cfc12       ', 'cfc113      ', &
    'hcfc22      ', 'hfc125      ', 'hfc134a     ', 'cfc114      ', &
    'tio         ', 'vo          ', 'h2          ', 'he          ', &
    'ocs         ', 'na          ', 'k           ', 'feh         ', &
    'crh         ', 'li          ', 'rb          ', 'cs          ', &
    'ph3         ', 'c2h2        ', 'hcn         ', 'h2s         ', &
    'ar          ', 'air         ', 'o           ', 'n           ', &
    'no3         ', 'n2o5        ', 'hono        ', 'ho2no2      ', &
    'h2o2        ', 'c2h6        ', 'ch3         ', 'h2co        ', &
    'ho2         ', 'hdo         ', 'hcl         ', 'hf          ', &
    'cosso       ', 'tosso       ', 'yosos       ', 'ch3cho      ', &
    'ch3ooh      ', 'ch3coch3    ', 'ch3cocho    ', 'chocho      ', &
    'c2h5cho     ', 'hoch2cho    ', 'c2h5coch3   ', 'mvk         ', &
    'macr        ', 'pan         ', 'ch3ono2     ', 'c2h3        ', &
    'c2h4        ', 'oh          ', 'hco         ', 'n2o4        ', &
    'c2n2        ', 'n2h4        ', 'n2o3        '/)

! Long names
CHARACTER (LEN=20), PARAMETER :: name_absorb(npd_gases) = (/ &
                                   "Water Vapour        ", &
                                   "Carbon Dioxide      ", &
                                   "Ozone               ", &
                                   "Dinitrogen Oxide    ", &
                                   "Carbon monoxide     ", &
                                   "Methane             ", &
                                   "Oxygen              ", &
                                   "Nitrogen monoxide   ", &
                                   "Sulphur dioxide     ", &
                                   "Nitrogen dioxide    ", &
                                   "Ammonia             ", &
                                   "Nitric acid         ", &
                                   "Nitrogen            ", &
                                   "CFC11               ", &
                                   "CFC12               ", &
                                   "CFC113              ", &
                                   "HCFC22              ", &
                                   "HFC125              ", &
                                   "HFC134A             ", &
                                   "CFC114              ", &
                                   "Titanium oxide      ", &
                                   "Vanadium oxide      ", &
                                   "Hydrogen            ", &
                                   "Helium              ", &
                                   "Carbonyl sulphide   ", &
                                   "Sodium              ", &
                                   "Potassium           ", &
                                   "Iron hydride        ", &
                                   "Chromium hydride    ", &
                                   "Lithium             ", &
                                   "Rubidium            ", &
                                   "Cesium              ", &
                                   "Phosphine           ", &
                                   "Acetylene           ", &
                                   "Hydrogen cyanide    ", &
                                   "Hydrogen sulphide   ", &
                                   "Argon               ", &
                                   "Dry air             ", &
                                   "Atomic oxygen       ", &
                                   "Atomic nitrogen     ", &
                                   "Nitrate radical     ", &
                                   "Dinitrogen pentoxide", &
                                   "Nitrous acid        ", &
                                   "Peroxynitric acid   ", &
                                   "Hydrogen peroxide   ", &
                                   "Ethane              ", &
                                   "Methyl radical      ", &
                                   "Formaldehyde        ", &
                                   "Hydroperoxy radical ", &
                                   "Semiheavy water     ", &
                                   "Hydrogen chloride   ", &
                                   "Hydrogen fluoride   ", &
                                   "cis-OSSO            ", &
                                   "trans-OSSO          ", &
                                   "OSO-S               ", &
                                   "Acetaldehyde        ", &
                                   "Methylhydroperoxide ", &
                                   "Acetone             ", &
                                   "Methylglyoxal       ", &
                                   "Glyoxal             ", &
                                   "Propanal            ", &
                                   "Glycolaldehyde      ", &
                                   "Methyl ethyl ketone ", &
                                   "Methyl vinyl ketone ", &
                                   "Methacrolein        ", &
                                   "Peroxyacetyl nitrate", &
                                   "Methylnitrate       ", &
                                   "Vinyl radical       ", &
                                   "Ethylene            ", &
                                   "Hydroxyl radical    ", &
                                   "Formyl radical      ", &
                                   "Dinitrogen tetroxide", &
                                   "Cyanogen            ", &
                                   "Hydrazine           ", &
                                   "Dinitrogen trioxide "/)


! Molecular weights taken from "General Inorganic Chemistry"
! by J. A. Duffy (1970), Longmans (except where stated).
REAL (RealK), PARAMETER :: molar_weight(npd_gases) = (/ &
  18.0153_RealK,     & ! H2O
  44.0100_RealK,     & ! CO2
  47.9982_RealK,     & ! O3
  44.0128_RealK,     & ! N2O
  28.0106_RealK,     & ! CO
  16.0430_RealK,     & ! CH4
  31.9988_RealK,     & ! O2
  30.0061_RealK,     & ! NO
  64.0628_RealK,     & ! SO2
  46.0055_RealK,     & ! NO2
  17.0306_RealK,     & ! NH3
  63.0129_RealK,     & ! HNO3
  28.0134_RealK,     & ! N2
  137.3686_RealK,    & ! CFC11
  120.9140_RealK,    & ! CFC12
  187.3765_RealK,    & ! CFC113
  86.46892_RealK,    & ! HCFC22
  120.02227_RealK,   & ! HFC125
  102.03184_RealK,   & ! HFC134a
  170.921_RealK,     & ! CFC114 (from NIST)
  63.866_RealK,      & ! TiO (from NIST)
  66.9409_RealK,     & ! VO (from NIST)
  2.01588_RealK,     & ! H2 (from NIST)
  4.002602_RealK,    & ! He (from NIST)
  60.075_RealK,      & ! OCS
  22.98976928_RealK, & ! Na (from NIST)
  39.0983_RealK,     & ! K (from NIST)
  56.853_RealK,      & ! FeH (from NIST)
  53.004_RealK,      & ! CrH (from NIST)
  6.941_RealK,       & ! Li (from NIST)
  85.4678_RealK,     & ! Rb (from NIST)
  132.9054519_RealK, & ! Cs (from NIST)
  33.99758_RealK,    & ! PH3 (from NIST)
  26.0373_RealK,     & ! C2H2 (from NIST)
  27.0253_RealK,     & ! HCN (from NIST)
  34.081_RealK,      & ! H2S (from NIST)
  39.948_RealK,      & ! Ar (from NIST)
  28.966_RealK,      & ! Dry air
  15.9994_RealK,     & ! O (from NIST)
  14.00674_RealK,    & ! N (from NIST)
  63.0128_RealK,     & ! NO3 (from NIST)
  108.0104_RealK,    & ! N2O5 (from NIST)
  47.0134_RealK,     & ! HONO (from NIST)
  79.0122_RealK,     & ! HO2NO2 (from NIST)
  34.0147_RealK,     & ! H2O2 (from NIST)
  30.0690_RealK,     & ! C2H6 (from NIST)
  15.0345_RealK,     & ! CH3  (from NIST)
  30.0260_RealK,     & ! H2CO (from NIST
  33.0067_RealK,     & ! HO2 (from NIST)
  19.0214_RealK,     & ! HDO (from NIST)
  36.461_RealK,      & ! HCl (from NIST)
  20.00689_RealK,    & ! HF (from NIST)
  96.129_RealK,      & ! cis-OSSO (from NIST)
  96.129_RealK,      & ! trans-OSSO (from NIST)
  96.129_RealK,      & ! OSO-S (from NIST)
  44.0526_RealK,     & ! CH3CHO (from NIST)
  48.0413_RealK,     & ! CH3OOH (from NIST)
  58.0791_RealK,     & ! CH3COCH3 (from NIST)
  72.0627_RealK,     & ! CH3COCHO (from NIST)
  58.0361_RealK,     & ! CHOCHO (from NIST)
  58.0791_RealK,     & ! C2H5CHO (from NIST)
  60.0520_RealK,     & ! HOCH2CHO (from NIST)
  72.1057_RealK,     & ! C2H5COCH3 (from NIST)
  70.0898_RealK,     & ! MVK (from NIST)
  70.0898_RealK,     & ! MACR (from NIST)
  121.0489_RealK,    & ! PAN (from NIST)
  77.0394_RealK,     & ! CH3ONO2 (from NIST)
  27.0452_RealK,     & ! C2H3 (from NIST)
  28.0532_RealK,     & ! C2H4 (from NIST)
  17.0073_RealK,     & ! OH (from NIST)
  29.0180_RealK,     & ! HCO (from NIST)
  91.0110_RealK,     & ! N2O4 (from NIST)
  52.0348_RealK,     & ! C2N2 (from NIST)
  32.0452_RealK,     & ! N2H4 (from NIST)
  76.0116_RealK     /) ! N2O3 (from NIST)


! Array of identifiers in HITRAN for each gas in the radiation code.
INTEGER, PARAMETER :: hitran_number(npd_gases) = (/ &
  1,   & ! H2O
  2,   & ! CO2
  3,   & ! O3
  4,   & ! N2O
  5,   & ! CO
  6,   & ! CH4
  7,   & ! O2
  8,   & ! NO
  9,   & ! SO2
  10,  & ! NO2
  11,  & ! NH3
  12,  & ! HNO3
  22,  & ! N2
  0,   & ! CFC11
  0,   & ! CFC12
  0,   & ! CFC113
  0,   & ! HCFC22
  0,   & ! HFC125
  0,   & ! HFC134a
  0,   & ! CFC114
  0,   & ! TiO
  0,   & ! VO
  45,  & ! H2
  0,   & ! He
  19,  & ! OCS
  0,   & ! Na
  0,   & ! K
  0,   & ! FeH
  0,   & ! CrH
  0,   & ! Li
  0,   & ! Rb
  0,   & ! Cs
  28,  & ! PH3
  26,  & ! C2H2
  23,  & ! HCN
  31,  & ! H2S
  0,   & ! Ar
  0,   & ! Dry air
  34,  & ! O
  0,   & ! N
  0,   & ! NO3
  0,   & ! N2O5
  0,   & ! HONO
  0,   & ! HO2NO2
  25,  & ! H2O2
  27,  & ! C2H6
  0,   & ! CH3
  20,  & ! H2CO
  33,  & ! HO2
  1,   & ! HDO
  15,  & ! HCl
  14,  & ! HF
  (0, i=ip_cosso, npd_gases) /)

! Maximum number of specified HITRAN isotopes for a given absorber
INTEGER, PARAMETER :: npd_isotopes = 3

! List of HITRAN isotopes for each absorber (0 for all isotopes)
INTEGER, PARAMETER :: hitran_isotopes(npd_isotopes, npd_gases) &
  = RESHAPE ( [INTEGER :: &
  (0, i=1, npd_isotopes*ip_ho2),            & ! H2O -> HO2
  4, 5, 6, (0, i=1, npd_isotopes-3),        & ! HDO: HD16O, HD18O, HD17O
  (0, i=1, npd_isotopes*(npd_gases-ip_hdo)) & ! HCl -> last gas
  ], shape=[npd_isotopes, npd_gases] )

! Depolarization factors used to compute the Rayleigh scattering coefficients
REAL (RealK), PARAMETER :: depolarization_factor(npd_gases) = (/ &
  0.0_RealK,     & ! H2O
  0.0922_RealK,  & ! CO2 (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0_RealK,     & ! O3
  0.1197_RealK,  & ! N2O (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.08_RealK,    & ! CO (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0_RealK,     & ! CH4
  0.06_RealK,    & ! O2 (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0218_RealK,  & ! NO (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0_RealK,     & ! SO2
  0.0_RealK,     & ! NO2
  0.0_RealK,     & ! NH3
  0.0_RealK,     & ! HNO3
  0.0305_RealK,  & ! N2 (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0_RealK,     & ! CFC11
  0.0_RealK,     & ! CFC12
  0.0_RealK,     & ! CFC113
  0.0_RealK,     & ! HCFC22
  0.0_RealK,     & ! HFC125
  0.0_RealK,     & ! HFC134a
  0.0_RealK,     & ! CFC114
  0.0_RealK,     & ! TiO
  0.0_RealK,     & ! VO
  0.0221_RealK,  & ! H2 (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.025_RealK,   & ! He (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0_RealK,     & ! OCS
  0.0_RealK,     & ! Na
  0.0_RealK,     & ! K
  0.0_RealK,     & ! FeH
  0.0_RealK,     & ! CrH
  0.0_RealK,     & ! Li
  0.0_RealK,     & ! Rb
  0.0_RealK,     & ! Cs
  0.0_RealK,     & ! PH3
  0.0_RealK,     & ! C2H2
  0.0_RealK,     & ! HCN
  0.0_RealK,     & ! H2S
  0.0006_RealK,  & ! Ar (Parthasarathy, Indian J. Phys. 25, 21 (1951))
  0.0279_RealK,  & ! Dry air
  0.0_RealK,     & ! O
  0.0_RealK,     & ! N
 (0.0_RealK, i=ip_no3, npd_gases) /)

! Minimum wavelength to consider Rayleigh scattering
REAL (RealK), PARAMETER :: rayleigh_cutoff(npd_gases) = (/ &
  175.0E-09_RealK, & ! H2O ( threshold for H2O -> O(1D) + H2 )
  167.1E-09_RealK, & ! CO2 ( threshold for CO2 -> CO + O(1D) )
  0.0_RealK,       & ! O3
  0.0_RealK,       & ! N2O
  0.0_RealK,       & ! CO
  0.0_RealK,       & ! CH4
  175.0E-09_RealK, & ! O2 ( threshold for O2 -> O(3P) + O(1D) )
  0.0_RealK,       & ! NO
  0.0_RealK,       & ! SO2
  0.0_RealK,       & ! NO2
  0.0_RealK,       & ! NH3
  0.0_RealK,       & ! HNO3
  79.8E-09_RealK,  & ! N2 ( threshold for N2 -> N2+ )
 (0.0_RealK, i=ip_cfc11, npd_gases) /)

! Maximum number of photolysis products for a given absorber
INTEGER, PARAMETER :: npd_products = 9
INTEGER, PARAMETER :: npd_pathway = npd_gases * (npd_products + 1)

CHARACTER(LEN=56), PARAMETER :: blank = ""
! Description of photolysis products
CHARACTER(LEN=56), PARAMETER :: photol_products(npd_products, npd_gases) &
  = RESHAPE([CHARACTER(LEN=56) ::    &
  "H2O -> O(3P) + H2             ",  &
  "H2O -> OH(X2Pi) + H           ",  &
  "H2O -> O(1D) + H2             ",  &
  "H2O -> OH(A2Sigma+) + H       ",  &
  "H2O -> O(3P) + H + H          ",  &
  (blank, i=1, npd_products-5),      & ! H2O
  "CO2 -> CO + O(3P)             ",  &
  "CO2 -> CO + O(1D)             ",  &
  "CO2 -> CO + O(1S)             ",  &
  "CO2 -> CO(a3Pi) + O(3P)       ",  &
  "CO2 -> CO2+                   ",  &
  "CO2 -> CO + O+                ",  &
  "CO2 -> CO+ + O(3P)            ",  &
  "CO2 -> O2 + C+                ",  &
  (blank, i=1, npd_products-8),      & ! CO2
  "O3 -> O(3P) + O2(X3Sigmag-)   ",  &
  "O3 -> O(3P) + O2(a1Deltag)    ",  &
  "O3 -> O(3P) + O2(b1Sigmag+)   ",  &
  "O3 -> O(1D) + O2(X3Sigmag-)   ",  &
  "O3 -> O(1D) + O2(a1Deltag)    ",  &
  "O3 -> O(1D) + O2(b1Sigmag+)   ",  &
  "O3 -> 3 O(3P)                 ",  &
  "O3 -> O(1S) + O2(a1Deltag)    ",  &
  (blank, i=1, npd_products-8),      & ! O3
  "N2O -> N2 + O(1D)             ",  &
  "N2O -> N2 + O(3P)             ",  &
  "N2O -> N(4S) + NO(2Pi)        ",  &
  "N2O -> N2 + O(1S)             ",  &
  (blank, i=1, npd_products-4),      & ! N2O
  "CO -> C + O(3P)               ",  &
  (blank, i=1, npd_products-1),      & ! CO
  "CH4 -> CH3 + H                ",  &
  "CH4 -> CH2(1) + H2            ",  &
  "CH4 -> CH2(3) + H + H         ",  &
  "CH4 -> CH + H2 + H            ",  &
  (blank, i=1, npd_products-4),      & ! CH4
  "O2 -> O(3P) + O(3P)           ",  &
  "O2 -> O(3P) + O(1D)           ",  &
  "O2 -> O(1D) + O(1D)           ",  &
  "O2 -> O(3P) + O(1S)           ",  &
  "O2 -> O(1D) + O(1S)           ",  &
  "O2 -> O2+                     ",  &
  "O2 -> O+ + O                  ",  &
  (blank, i=1, npd_products-7),      & ! O2
  "NO -> O(3P) + N(4S)           ",  &
  (blank, i=1, npd_products-1),      & ! NO
  "SO2 -> SO + O(3P)             ",  &
  (blank, i=1, npd_products-1),      & ! SO2
  "NO2 -> NO + O(3P)             ",  &
  "NO2 -> NO + O(1D)             ",  &
  (blank, i=1, npd_products-2),      & ! NO2
  "NH3 -> NH2 + H                ",  &
  (blank, i=1, npd_products-1),      & ! NH3
  "HNO3 -> OH + NO2              ",  &
  "HNO3 -> HONO + O(3P)          ",  &
  "HNO3 -> H + NO3               ",  &
  "HNO3 -> OH + NO2*(12B2)       ",  &
  "HNO3 -> HONO + O(1D)          ",  &
  "HNO3 -> HONO (a3A)+ O(3P)     ",  &
  (blank, i=1, npd_products-6),      & ! HNO3
  "N2 -> N + N                   ",  &
  "N2 -> N2+                     ",  &
  "N2 -> N+ + N                  ",  &
  (blank, i=1, npd_products-3),      & ! N2
  (blank, i=1, npd_products),        & ! CFC11
  (blank, i=1, npd_products),        & ! CFC12
  (blank, i=1, npd_products),        & ! CFC113
  (blank, i=1, npd_products),        & ! HCFC22
  (blank, i=1, npd_products),        & ! HFC125
  (blank, i=1, npd_products),        & ! HFC134a
  (blank, i=1, npd_products),        & ! CFC114
  (blank, i=1, npd_products),        & ! TiO
  (blank, i=1, npd_products),        & ! VO
  "H2 -> H + H                   ",  &
  (blank, i=1, npd_products-1),      & ! H2
  (blank, i=1, npd_products),        & ! He
  "OCS -> CO + S(3P)             ",  &
  "OCS -> CO + S(1D)             ",  &
  "OCS -> CO + S(1S)             ",  &
  (blank, i=1, npd_products-3),      & ! OCS
  (blank, i=1, npd_products),        & ! Na
  (blank, i=1, npd_products),        & ! K
  (blank, i=1, npd_products),        & ! FeH
  (blank, i=1, npd_products),        & ! CrH
  (blank, i=1, npd_products),        & ! Li
  (blank, i=1, npd_products),        & ! Rb
  (blank, i=1, npd_products),        & ! Cs
  (blank, i=1, npd_products),        & ! PH3
  "C2H2 -> C2H + H               ",  &
  (blank, i=1, npd_products-1),      & ! C2H2
  "HCN -> CN + H                 ",  &
  (blank, i=1, npd_products-1),      & ! HCN
  (blank, i=1, npd_products),        & ! H2S
  (blank, i=1, npd_products),        & ! Ar
  (blank, i=1, npd_products),        & ! Dry air
  "O -> O+(4S)                   ",  &
  "O -> O+(2D)                   ",  &
  "O -> O+(2P)                   ",  &
  "O -> O+(4Pe)                  ",  &
  "O -> O+(2Pe)                  ",  &
  "O -> O++                      ",  &
  "O -> O+++                     ",  &
  (blank, i=1, npd_products-7),      & ! O
  "N -> N+                       ",  &
  "N -> N++                      ",  &
  (blank, i=1, npd_products-2),      & ! N
  "NO3 -> NO + O2                ",  &
  "NO3 -> NO2 + O(3P)            ",  &
  (blank, i=1, npd_products-2)    ,  & ! NO3
  "N2O5 -> NO3 + NO2             ",  &
  "N2O5 -> NO3 + NO + O(3P)      ",  &
  (blank, i=1, npd_products-2)    ,  & ! N2O5
  "HONO -> OH + NO               ",  &
  "HONO -> H + NO2               ",  &
  (blank, i=1, npd_products-2)    ,  & ! HONO
  "HO2NO2 -> HO2 + NO2           ",  &
  "HO2NO2 -> OH + NO3            ",  &
  "HO2NO2 -> O(3P) + HNO3        ",  &
  "HO2NO2 -> H + NO2 + O2        ",  &
  "HO2NO2 -> HO2 + NO + O(3P)    ",  &
  "HO2NO2 -> OH + NO2 + O(3P)    ",  &
  "HO2NO2 -> H + O(3P) + NO3     ",  &
  "HO2NO2 -> HONO + O2(1Sigma)   ",  &
  "HO2NO2 -> HONO + O2(1Lambda)  ",  &
  (blank, i=1, npd_products-9)    ,  & ! HO2NO2
  "H2O2 -> OH + OH               ",  &
  "H2O2 -> H2O + O(1D)           ",  &
  "H2O2 -> H + HO2               ",  &
  (blank, i=1, npd_products-3)    ,  & ! H2O2
  "C2H6 -> C2H4 + H2             ",  &
  "C2H6 -> C2H4 + H + H          ",  &
  "C2H6 -> C2H2 + H2 + H2        ",  &
  "C2H6 -> CH4 + CH2(1)          ",  &
  "C2H6 -> CH3 + CH3             ",  &
  (blank, i=1, npd_products-5)    ,  & ! C2H6
  "CH3 -> CH2(1) + H             ",  & !
  (blank, i=1, npd_products-1)    ,  & ! CH3
  "H2CO -> H + HCO               ",  &
  "H2CO -> H2 + CO               ",  &
  "H2CO -> H + H + CO            ",  &
  (blank, i=1, npd_products-3)    ,  & ! H2CO
  "HO2 -> OH + O(3P)             ",  &
  "HO2 -> OH + O(1D)             ",  &
  (blank, i=1, npd_products-2)    ,  & ! HO2
  (blank, i=1, npd_products)      ,  & ! HDO
  (blank, i=1, npd_products)      ,  & ! HCl
  (blank, i=1, npd_products)      ,  & ! HF
  (blank, i=1, npd_products)      ,  & ! cis-OSSO
  (blank, i=1, npd_products)      ,  & ! trans-OSSO
  (blank, i=1, npd_products)      ,  & ! OSO-S
  "CH3CHO -> CH3 + HCO           ",  &
  "CH3CHO -> CH4 + CO            ",  &
  "CH3CHO -> CH3CO + H           ",  &
  (blank, i=1, npd_products-3)    ,  & ! CH3CHO
  "CH3OOH -> CH3O + OH           ",  &
  "CH3OOH -> CH3 + HO2           ",  &
  "CH3OOH -> CH3O2 + H           ",  &
  "CH3OOH -> CH3OH + O(3P)       ",  &
  (blank, i=1, npd_products-4)    ,  & ! CH3OOH
  "CH3COCH3 -> CH3CO + CH3       ",  &
  "CH3COCH3 -> 2 CH3 + CO        ",  &
  (blank, i=1, npd_products-2)    ,  & ! CH3COCH3
  "CH3COCHO -> CH3CO + HCO       ",  &
  "CH3COCHO -> CH4 + 2 CO        ",  &
  "CH3COCHO -> CH3CHO + CO       ",  &
  (blank, i=1, npd_products-3)    ,  & ! CH3COCHO
  "CHOCHO -> HCO + HCO           ",  &
  "CHOCHO -> H2 + 2 CO           ",  &
  "CHOCHO -> HCHO + CO           ",  &
  "CHOCHO -> H + CO + HCO        ",  &
  (blank, i=1, npd_products-4)    ,  & ! CHOCHO
  "C2H5CHO -> C2H5 + HCO         ",  &
  "C2H5CHO -> C2H6 + CO          ",  &
  "C2H5CHO -> C2H4 + HCHO        ",  &
  "C2H5CHO -> CH3 + CH2CHO       ",  &
  (blank, i=1, npd_products-4)    ,  & ! C2H5CHO
  "HOCH2CHO -> CH2OH + HCO       ",  &
  "HOCH2CHO -> CH3OH + CO        ",  &
  "HOCH2CHO -> OH + CH2CHO       ",  &
  "HOCH2CHO -> HOCH2CO + H       ",  &
  (blank, i=1, npd_products-4)    ,  & ! HOCH2CHO
  "C2H5COCH3 -> CH3CO + C2H5     ",  &
  "C2H5COCH3 -> C2H5CO + CH3     ",  &
  "C2H5COCH3 -> C2H5 + CO + CH3  ",  &
  (blank, i=1, npd_products-3)    ,  & ! C2H5COCH3
  "CH3C(O)CH=CH2 -> CH3-CH=CH2 + CO  ", &
  "CH3C(O)CH=CH2 -> CH=CH2 + CH3C(O) ", &
  "CH3C(O)CH=CH2 -> CH=CH2C(O) + CH3 ", &
  (blank, i=1, npd_products-3)        , & ! CH3C(O)CH=CH2 (MVK)
  "CH2=C(CH3)CHO -> CH2=CCH3 + CHO   ", &
  "CH2=C(CH3)CHO -> C3H6 + CO        ", &
  "CH2=C(CH3)CHO -> H + CH2=C(CH3)CO ", &
  (blank, i=1, npd_products-3)        , & ! CH2=C(CH3)CHO (MACR)
  "CH3C(O)OONO2 -> CH3C(O)OO + NO2   ", &
  "CH3C(O)OONO2 -> CH3C(O)O + NO3    ", &
  (blank, i=1, npd_products-2)        , & ! CH3C(O)OONO2 (PAN)
  "CH3ONO2 -> CH3O + NO2             ", &
  "CH3ONO2 -> HCHO + HONO            ", &
  "CH3ONO2 -> HCHO + NO + OH         ", &
  "CH3ONO2 -> CH3ONO + O(3P)         ", &
  "CH3ONO2 -> CH3 + NO3              ", &
  "CH3ONO2 -> CH2ONO2 + H            ", &
  "CH3ONO2 -> CH3O + NO + O(3P)      ", &
  "CH3ONO2 -> CH3ONO + O(1D)         ", &
  (blank, i=1, npd_products-8)        , & ! CH3ONO2
  "C2H3 -> C2H2 + H                  ", &
  (blank, i=1, npd_products-1)        , & ! C2H3
  "C2H4 -> C2H2 + H2                 ", & 
  "C2H4 -> C2H2 + H + H              ", &
  (blank, i=1, npd_products-2)        , & ! C2H4
  "OH -> O(1D) + H                   ", &
  (blank, i=1, npd_products-1)        , & ! OH
  "HCO -> H + CO                     ", &
  (blank, i=1, npd_products-1)        , & ! HCO
  "N2O4 -> NO2 + NO2                 ", &
  (blank, i=1, npd_products-1)        , & ! N2O4
  "C2N2 -> C2 + N2                   ", &
  (blank, i=1, npd_products-1)        , & ! C2N2
  "N2H4 -> N2H3 + H                  ", &
  (blank, i=1, npd_products-1)        , & ! N2H4
  "N2O3 -> NO2 + NO                  ", &
  (blank, i=1, npd_products-1)          & ! N2O3
  ], shape=[npd_products, npd_gases] )

! Name used by UKCA for photolysis pathway
CHARACTER(LEN=56), PARAMETER :: photol_fldname(0:npd_products, npd_gases) &
  = RESHAPE([CHARACTER(LEN=56) ::    &
  (blank, i=0, 1),                   & ! H2O
  "jh2o                          ",  & ! H2O -> OH(X2Pi) + H
  (blank, i=3, npd_products),        & ! H2O
  (blank, i=0, 0),                   & ! CO2
  "jco2                          ",  & ! CO2 -> CO + O(3P)
  (blank, i=2, npd_products),        & ! CO2
  (blank, i=0, 0),                   & ! O3
  "jo3b                          ",  & ! O3 -> O(3P) + O2(X3Sigmag-)
  (blank, i=2, 3),                   & ! O3
  "jo3a                          ",  & ! O3 -> O(1D) + O2(X3Sigmag-)
  (blank, i=5, npd_products),        & ! O3
  (blank, i=0, npd_products),        & ! N2O
  (blank, i=0, npd_products),        & ! CO
  (blank, i=0, npd_products),        & ! CH4
  "jo2                           ",  & ! O2 -> Unspecified
  "jo2                           ",  & ! O2 -> O(3P) + O(3P)
  "jo2b                          ",  & ! O2 -> O(3P) + O(1D)
  (blank, i=3, npd_products),        & ! O2
  (blank, i=0, npd_products),        & ! NO
  (blank, i=0, npd_products),        & ! SO2
  (blank, i=0, 0),                   & ! NO2
  "jno2                          ",  & ! NO2 -> NO + O(3P)
  (blank, i=2, npd_products),        & ! NO2
  (blank, i=0, npd_products),        & ! NH3
  (blank, i=0, 0),                   & ! HNO3
  "jhono2                        ",  & ! HNO3 -> OH + NO2
  (blank, i=2, npd_products),        & ! HNO3
  (blank, i=0, npd_products),        & ! N2
  (blank, i=0, npd_products),        & ! CFC11
  (blank, i=0, npd_products),        & ! CFC12
  (blank, i=0, npd_products),        & ! CFC113
  (blank, i=0, npd_products),        & ! HCFC22
  (blank, i=0, npd_products),        & ! HFC125
  (blank, i=0, npd_products),        & ! HFC134a
  (blank, i=0, npd_products),        & ! CFC114
  (blank, i=0, npd_products),        & ! TiO
  (blank, i=0, npd_products),        & ! VO
  (blank, i=0, npd_products),        & ! H2
  (blank, i=0, npd_products),        & ! He
  (blank, i=0, npd_products),        & ! OCS
  (blank, i=0, npd_products),        & ! Na
  (blank, i=0, npd_products),        & ! K
  (blank, i=0, npd_products),        & ! FeH
  (blank, i=0, npd_products),        & ! CrH
  (blank, i=0, npd_products),        & ! Li
  (blank, i=0, npd_products),        & ! Rb
  (blank, i=0, npd_products),        & ! Cs
  (blank, i=0, npd_products),        & ! PH3
  (blank, i=0, npd_products),        & ! C2H2
  (blank, i=0, npd_products),        & ! HCN
  (blank, i=0, npd_products),        & ! H2S
  (blank, i=0, npd_products),        & ! Ar
  (blank, i=0, npd_products),        & ! Dry air
  (blank, i=0, npd_products),        & ! O
  (blank, i=0, npd_products),        & ! N
  (blank, i=0, 0),                   & ! NO3
  "jno3a                         ",  & ! NO3 -> NO + O2
  "jno3b                         ",  & ! NO3 -> NO2 + O(3P)
  (blank, i=3, npd_products),        & ! NO3
  (blank, i=0, 0),                   & ! N2O5
  "jn2o5                         ",  & ! N2O5 -> NO3 + NO2
  (blank, i=2, npd_products),        & ! N2O5
  (blank, i=0, 0),                   & ! HONO
  "jhono                         ",  & ! HONO -> OH + NO
  (blank, i=2, npd_products),        & ! HONO
  "jpna                          ",  & ! HO2NO2 -> Unspecified
  "jpna67                        ",  & ! HO2NO2 -> HO2 + NO2
  "jpna33                        ",  & ! HO2NO2 -> OH + NO3
  (blank, i=3, npd_products),        & ! HO2NO2
  (blank, i=0, 0),                   & ! H2O2
  "jh2o2                         ",  & ! H2O2 -> OH + OH
  (blank, i=2, npd_products),        & ! H2O2
  (blank, i=0, npd_products),        & ! C2H6
  (blank, i=0, npd_products),        & ! CH3
  "jhchoa                        ",  & ! H2CO -> Unspecified
  (blank, i=1, 1),                   & ! H2CO
  "jhchob                        ",  & ! H2CO -> H2 + CO
  (blank, i=3, npd_products),        & ! H2CO
  (blank, i=0, npd_products),        & ! HO2
  (blank, i=0, npd_products),        & ! HDO
  (blank, i=0, npd_products),        & ! HCl
  (blank, i=0, npd_products),        & ! HF
  (blank, i=0, npd_products),        & ! cis-OSSO
  (blank, i=0, npd_products),        & ! trans-OSSO
  (blank, i=0, npd_products),        & ! OSO-S
  "jaceta                        ",  & ! CH3CHO -> Unspecified
  (blank, i=1, 1),                   & ! CH3CHO
  "jacetb                        ",  & ! CH3CHO -> CH4 + CO
  (blank, i=3, npd_products),        & ! CH3CHO
  "jmhp                          ",  & ! CH3OOH -> Unspecified
  (blank, i=1, npd_products),        & ! CH3OOH
  "jaceto                        ",  & ! CH3COCH3 -> Unspecified
  (blank, i=1, npd_products),        & ! CH3COCH3
  "jmkal                         ",  & ! CH3COCHO -> Unspecified
  (blank, i=1, npd_products),        & ! CH3COCHO
  (blank, i=0, npd_products),        & ! CHOCHO
  "jetcho                        ",  & ! C2H5CHO -> Unspecified
  (blank, i=1, npd_products),        & ! C2H5CHO
  (blank, i=0, npd_products),        & ! HOCH2CHO
  (blank, i=0, npd_products),        & ! C2H5COCH3
  "jmvk                          ",  & ! MVK -> Unspecified
  (blank, i=1, npd_products),        & ! MVK
  "jmacr                         ",  & ! MACR -> Unspecified
  (blank, i=1, npd_products),        & ! MACR
  (blank, i=0, 0),                   & ! PAN
  "jpan                          ",  & ! CH3C(O)OONO2 -> CH3C(O)OO + NO2
  (blank, i=2, npd_products),        & ! PAN
  "jmena                         ",  & ! CH3ONO2 -> Unspecified
  (blank, i=1, npd_products),        & ! CH3ONO2
  (blank, i=0, npd_products),        & ! C2H3
  (blank, i=0, npd_products),        & ! C2H4
  (blank, i=0, npd_products),        & ! OH
  (blank, i=0, npd_products),        & ! HCO
  (blank, i=0, npd_products),        & ! N2O4
  (blank, i=0, npd_products),        & ! C2N2
  (blank, i=0, npd_products),        & ! N2H4
  (blank, i=0, npd_products)         & ! N2O3
  ], shape=[npd_products+1, npd_gases] )

! Threshold wavelength defining energy required for photolysis
REAL (RealK), PARAMETER :: threshold_wavelength(npd_products, npd_gases) &
  = RESHAPE ( [REAL(RealK) ::       &
  246.0E-09_RealK,                  & ! H2O -> O(3P) + H2
  242.0E-09_RealK,                  & ! H2O -> OH(X2Pi) + H
  175.0E-09_RealK,                  & ! H2O -> O(1D) + H2
  134.0E-09_RealK,                  & ! H2O -> OH(A2Sigma+) + H
  129.0E-09_RealK,                  & ! H2O -> O(3P) + H + H
  (0.0_RealK, i=1, npd_products-5), & ! H2O
  227.5E-09_RealK,                  & ! CO2 -> CO + O(3P) : Heubner 92
  167.1E-09_RealK,                  & ! CO2 -> CO + O(1D) : Heubner 92
  128.6E-09_RealK,                  & ! CO2 -> CO + O(1S) : Heubner 92
  108.2E-09_RealK,                  & ! CO2 -> CO(a3Pi) + O(3P) : Heubner 92
  89.922E-09_RealK,                 & ! CO2 -> CO2+ : Heubner 92
  65.026E-09_RealK,                 & ! CO2 -> CO + O+ : Heubner 92
  63.693E-09_RealK,                 & ! CO2 -> CO+ + O : Heubner 92
  54.655E-09_RealK,                 & ! CO2 -> O2 + C+ : Heubner 92
  (0.0_RealK, i=1, npd_products-8), & ! CO2
  1180.0E-09_RealK,                 & ! O3 -> O(3P) + O2(X3Sigmag-)
   612.0E-09_RealK,                 & ! O3 -> O(3P) + O2(a1Deltag)
   463.0E-09_RealK,                 & ! O3 -> O(3P) + O2(b1Sigmag+)
   411.0E-09_RealK,                 & ! O3 -> O(1D) + O2(X3Sigmag-)
   310.0E-09_RealK,                 & ! O3 -> O(1D) + O2(a1Deltag)
   267.0E-09_RealK,                 & ! O3 -> O(1D) + O2(b1Sigmag+)
   201.0E-09_RealK,                 & ! O3 -> 3 O(3P)
   196.0E-09_RealK,                 & ! O3 -> O(1S) + O2(a1Deltag)
  (0.0_RealK, i=1, npd_products-8), & ! O3
  336.0E-09_RealK,                  & ! N2O -> N2 + O(1D)
  713.0E-09_RealK,                  & ! N2O -> N2 + O(3P)
  248.0E-09_RealK,                  & ! N2O -> N(4S) + NO(2Pi)
  210.0E-09_RealK,                  & ! N2O -> N2 + O(1S)
  (0.0_RealK, i=1, npd_products-4), & ! N2O
  111.78E-09_RealK,                 & ! CO -> CO + O(3P) :Huebner 92
  (0.0_RealK, i=1, npd_products-1), & ! CO
  277.0E-09_RealK,                  & ! CH4 -> CH3 + H
  237.3E-09_RealK,                  & ! CH4 -> CH2(1) + H2
  132.2E-09_RealK,                  & ! CH4 -> CH2(3) + H + H
  137.0E-09_RealK,                  & ! CH4 -> CH + H2 + H
  (0.0_RealK, i=1, npd_products-4), & ! CH4 :Huebner 92
  242.3E-09_RealK,                  & ! O2 -> O(3P) + O(3P)
  175.0E-09_RealK,                  & ! O2 -> O(3P) + O(1D)
  137.0E-09_RealK,                  & ! O2 -> O(1D) + O(1D)
  132.0E-09_RealK,                  & ! O2 -> O(3P) + O(1S)
  110.0E-09_RealK,                  & ! O2 -> O(1D) + O(1S)
  102.78E-09_RealK,                 & ! O2 -> O2+
   66.2E-09_RealK,                  & ! O2 -> O+ + O
  (0.0_RealK, i=1, npd_products-7), & ! O2
  191.0E-09_RealK,                  & ! NO -> O(3P) + N(4S) :Huebner 92
  (0.0_RealK, i=1, npd_products-1), & ! NO
  218.7E-09_RealK,                  & ! SO2 -> SO + O(3P) : Becker 95
  (0.0_RealK, i=1, npd_products-1), & ! SO2
  398.0E-09_RealK,                  & ! NO2 -> NO + O(3P)
  244.0E-09_RealK,                  & ! NO2 -> NO + O(1D)
  (0.0_RealK, i=1, npd_products-2), & ! NO2
  279.8E-09_RealK,                  & ! NH3 -> NH2 + H :Huebner
  (0.0_RealK, i=1, npd_products-1), & ! NH3
  604.0E-09_RealK,                  & ! HNO3 -> OH + NO2
  393.0E-09_RealK,                  & ! HNO3 -> HONO + O(3P)
  278.0E-09_RealK,                  & ! HNO3 -> H + NO3
  381.0E-09_RealK,                  & ! HNO3 -> OH + NO2*(12B2)
  242.0E-09_RealK,                  & ! HNO3 -> HONO + O(1D)
  216.0E-09_RealK,                  & ! HNO3 -> HONO (a3A)+ O(3P)
  (0.0_RealK, i=1, npd_products-6), & ! HNO3 : IUPAC PNOx2
   98.6E-09_RealK,                  & ! N2 -> N + N
   79.8E-09_RealK,                  & ! N2 -> N2+
   51.0E-09_RealK,                  & ! N2 -> N+ + N
  (0.0_RealK, i=1, npd_products-3), & ! N2
  (0.0_RealK, i=1, npd_products),   & ! CFC11
  (0.0_RealK, i=1, npd_products),   & ! CFC12
  (0.0_RealK, i=1, npd_products),   & ! CFC113
  (0.0_RealK, i=1, npd_products),   & ! HCFC22
  (0.0_RealK, i=1, npd_products),   & ! HFC125
  (0.0_RealK, i=1, npd_products),   & ! HFC134a
  (0.0_RealK, i=1, npd_products),   & ! CFC114
  (0.0_RealK, i=1, npd_products),   & ! TiO
  (0.0_RealK, i=1, npd_products),   & ! VO
  276.885E-09_RealK,                & ! H2 -> H + H :Huebner 92
  (0.0_RealK, i=1, npd_products-1), & ! H2
  (0.0_RealK, i=1, npd_products),   & ! He
  388.0E-09_RealK,                  & ! OCS -> CO + S(3P)
  285.0E-09_RealK,                  & ! OCS -> CO + S(1D)
  209.0E-09_RealK,                  & ! OCS -> CO + S(1S)
  (0.0_RealK, i=1, npd_products-3), & ! OCS
  (0.0_RealK, i=1, npd_products),   & ! Na
  (0.0_RealK, i=1, npd_products),   & ! K
  (0.0_RealK, i=1, npd_products),   & ! FeH
  (0.0_RealK, i=1, npd_products),   & ! CrH
  (0.0_RealK, i=1, npd_products),   & ! Li
  (0.0_RealK, i=1, npd_products),   & ! Rb
  (0.0_RealK, i=1, npd_products),   & ! Cs
  (0.0_RealK, i=1, npd_products),   & ! PH3
  230.6E-09_RealK,                  & ! C2H2 -> C2H + H :Huebner 92
  (0.0_RealK, i=1, npd_products-1), & ! C2H2
  195.0E-09_RealK,                  & ! HCN -> CN + H :Huebner 92
  (0.0_RealK, i=1, npd_products-1), & ! HCN
  (0.0_RealK, i=1, npd_products),   & ! H2S
  (0.0_RealK, i=1, npd_products),   & ! Ar
  (0.0_RealK, i=1, npd_products),   & ! Dry air
   91.25E-09_RealK,                 & ! O -> O+(4S) 
   73.18E-09_RealK,                 & ! O -> O+(2D) 
   66.58E-09_RealK,                 & ! O -> O+(2P) 
   43.50E-09_RealK,                 & ! O -> O+(4Pe)
   31.00E-09_RealK,                 & ! O -> O+(2Pe)
   24.80E-09_RealK,                 & ! O -> O++    
   12.179E-09_RealK,                & ! O -> O+++   
  (0.0_RealK, i=1, npd_products-7), & ! O
   85.92E-09_RealK,                 & ! N -> N+
   28.00E-09_RealK,                 & ! N -> N++
  (0.0_RealK, i=1, npd_products-2), & ! N
   7320.0E-09_RealK,                & ! NO3 -> NO + O2 : JPL 19-5
   574.0E-09_RealK,                 & ! NO3 -> NO2 + O(3P) : JPL 19-5
  (0.0_RealK, i=1, npd_products-2), & ! NO3
   1255.0E-09_RealK,                & ! N2O5 -> NO3 + NO2 : JPL 19-5
   298.0E-09_RealK,                 & ! N2O5 -> NO3 + NO + O(3P) : JPL 19-5
  (0.0_RealK, i=1, npd_products-2), & ! N2O5
   579.0E-09_RealK,                 & ! HONO -> OH + NO : JPL 19-5
   362.0E-09_RealK,                 & ! HONO -> H + NO2 : JPL 19-5
  (0.0_RealK, i=1, npd_products-2), & ! HONO
  1207.0E-09_RealK,                 & ! HO2NO2 -> HO2 + NO2 : JPL 19-5
   726.0E-09_RealK,                 & ! HO2NO2 -> OH + NO3 : JPL 19-5
   713.0E-09_RealK,                 & ! HO2NO2 -> O(3P) + HNO3 : JPL 19-5
   393.0E-09_RealK,                 & ! HO2NO2 -> H + NO2 + O2 : JPL 19-5
   339.0E-09_RealK,                 & ! HO2NO2 -> HO2 + NO + O(3P) : JPL 19-5
   321.0E-09_RealK,                 & ! HO2NO2 -> OH + NO2 + O(3P) : JPL 19-5
   201.0E-09_RealK,                 & ! HO2NO2 -> H + O(3P) + NO3 : JPL 19-5
   911.0E-09_RealK,                 & ! HO2NO2 -> HONO + O2(1Sigma) : JPL 19-5
   1744.0E-09_RealK,                & ! HO2NO2 -> HONO + O2(1Lambda) : JPL 19-5
  (0.0_RealK, i=1, npd_products-9), & ! HO2NO2
   557.0E-09_RealK,                 & ! H2O2 -> OH + OH : JPL 19-5
   359.0E-09_RealK,                 & ! H2O2 -> H2O + O(1D) : JPL 19-5
   324.0E-09_RealK,                 & ! H2O2 -> H + HO2 : JPL 19-5
  (0.0_RealK, i=1, npd_products-3), & ! H2O2
   874.3E-09_RealK,                 & ! C2H6 -> C2H4 + H2 : Huebner 92
   290.0E-09_RealK,                 & ! C2H6 -> C2H4 + H + H :Huebner 92/Lias 70
   428.7E-09_RealK,                 & ! C2H6 -> C2H2 + H2 + H2 :Chang 2020
   272.6E-09_RealK,                 & ! C2H6 -> CH4 + CH2(1) :Huebner 92
   322.0E-09_RealK,                 & ! C2H6 -> CH3 + CH3 :Huebner 92
  (0.0_RealK, i=1, npd_products-5), & ! C2H6
   216.0E-09_RealK,                 & ! CH3 -> CH2(1) + H
  (0.0_RealK, i=1, npd_products-1), & ! CH3- Venot 2012-highest non zero qy
   330.0E-09_RealK,                 & ! H2CO -> H + HCO : JPL 19-5
   361.0E-09_RealK,                 & ! H2CO -> H2 + CO : JPL 19-5
   283.0E-09_RealK,                 & ! H2CO -> H + H + CO : JPL 19-5
  (0.0_RealK, i=1, npd_products-3), & ! H2CO
   438.0E-09_RealK,                 & ! HO2 -> OH + O(3P) : JPL 19-5
   259.0E-09_RealK,                 & ! HO2 -> OH + O(1D) : JPL 19-5
  (0.0_RealK, i=1, npd_products-2), & ! HO2
  (0.0_RealK, i=1, npd_products),   & ! HDO
  (0.0_RealK, i=1, npd_products),   & ! HCl
  (0.0_RealK, i=1, npd_products),   & ! HF
  (0.0_RealK, i=1, npd_products),   & ! cis-OSSO
  (0.0_RealK, i=1, npd_products),   & ! trans-OSSO
  (0.0_RealK, i=1, npd_products),   & ! OSO-S
   335.0E-09_RealK,                 & ! CH3CHO -> CH3 + HCO : JPL 19-5
   1.0_RealK,                       & ! CH3CHO -> CH4 + CO : JPL 19-5
   320.0E-09_RealK,                 & ! CH3CHO -> CH3CO + H : JPL 19-5
  (0.0_RealK, i=1, npd_products-3), & ! CH3CHO
   645.0E-09_RealK,                 & ! CH3OOH -> CH3O + OH : JPL 19-5
   410.0E-09_RealK,                 & ! CH3OOH -> CH3 + HO2 : JPL 19-5
   334.0E-09_RealK,                 & ! CH3OOH -> CH3O2 + H : JPL 19-5
   670.0E-09_RealK,                 & ! CH3OOH -> CH3OH + O(3P) : JPL 19-5
  (0.0_RealK, i=1, npd_products-4), & ! CH3OOH
   338.0E-09_RealK,                 & ! CH3COCH3 -> CH3CO + CH3 : JPL 19-5
   299.0E-09_RealK,                 & ! CH3COCH3 -> 2 CH3 + CO : JPL 19-5
  (0.0_RealK, i=1, npd_products-2), & ! CH3COCH3
   422.0E-09_RealK,                 & ! CH3COCHO -> CH3CO + HCO : JPL 19-5
   1.0_RealK,                       & ! CH3COCHO -> CH4 + 2 CO : JPL 19-5
   1.0_RealK,                       & ! CH3COCHO -> CH3CHO + CO : JPL 19-5
  (0.0_RealK, i=1, npd_products-3), & ! CH3COCHO
   399.0E-09_RealK,                 & ! CHOCHO -> HCO + HCO : JPL 19-5
   1.0_RealK,                       & ! CHOCHO -> H2 + 2 CO : JPL 19-5
   1.0_RealK,                       & ! CHOCHO -> HCHO + CO : JPL 19-5
   329.0E-09_RealK,                 & ! CHOCHO -> H + CO + HCO : JPL 19-5
  (0.0_RealK, i=1, npd_products-4), & ! CHOCHO
   341.0E-09_RealK,                 & ! C2H5CHO -> C2H5 + HCO : JPL 19-5
   1.0_RealK,                       & ! C2H5CHO -> C2H6 + CO : JPL 19-5
   926.0E-09_RealK,                 & ! C2H5CHO -> C2H4 + HCHO : JPL 19-5
   349.0E-09_RealK,                 & ! C2H5CHO -> CH3 + CH2CHO : JPL 19-5
  (0.0_RealK, i=1, npd_products-4), & ! C2H5CHO
   347.0E-09_RealK,                 & ! HOCH2CHO -> CH2OH + HCO
  1910.0E-09_RealK,                 & ! HOCH2CHO -> CH3OH + CO 
   328.0E-09_RealK,                 & ! HOCH2CHO -> OH + CH2CHO
   314.0E-09_RealK,                 & ! HOCH2CHO -> HOCH2CO + H
  (0.0_RealK, i=1, npd_products-4), & ! HOCH2CHO
   342.0E-09_RealK,                 & ! C2H5COCH3 -> CH3CO + C2H5
   339.0E-09_RealK,                 & ! C2H5COCH3 -> C2H5CO + CH3
   303.0E-09_RealK,                 & ! C2H5COCH3 -> C2H5 + CO + CH3
  (0.0_RealK, i=1, npd_products-3), & ! C2H5COCH3 : Zborowska 21
   1.0_RealK,                       & ! CH3C(O)CH=CH2 -> CH3-CH=CH2 + CO
   1.0_RealK,                       & ! CH3C(O)CH=CH2 -> CH=CH2 + CH3C(O)
   1.0_RealK,                       & ! CH3C(O)CH=CH2 -> CH=CH2C(O) + CH3
  (0.0_RealK, i=1, npd_products-3), & ! CH3C(O)CH=CH2 (MVK) : JPL 19-5
   1.0_RealK,                       & ! CH2=C(CH3)CHO -> CH2=CCH3 + CHO
   1.0_RealK,                       & ! CH2=C(CH3)CHO -> C3H6 + CO
   1.0_RealK,                       & ! CH2=C(CH3)CHO -> H + CH2=C(CH3)CO
  (0.0_RealK, i=1, npd_products-3), & ! CH2=C(CH3)CHO (MACR) : JPL 19-5
  1004.0E-09_RealK,                 & ! CH3C(O)OONO2 -> CH3C(O)OO + NO2
   963.0E-09_RealK,                 & ! CH3C(O)OONO2 -> CH3C(O)O + NO3
  (0.0_RealK, i=1, npd_products-2), & ! CH3C(O)OONO2 (PAN) : IUPAC P21
   697.0E-09_RealK,                 & ! CH3ONO2 -> CH3O + NO2       
   1.0_RealK,                       & ! CH3ONO2 -> HCHO + HONO      
   497.0E-09_RealK,                 & ! CH3ONO2 -> HCHO + NO + OH   
   391.0E-09_RealK,                 & ! CH3ONO2 -> CH3ONO + O(3P)   
   344.0E-09_RealK,                 & ! CH3ONO2 -> CH3 + NO3        
   294.0E-09_RealK,                 & ! CH3ONO2 -> CH2ONO2 + H      
   250.0E-09_RealK,                 & ! CH3ONO2 -> CH3O + NO + O(3P)
   241.0E-09_RealK,                 & ! CH3ONO2 -> CH3ONO + O(1D)   
  (0.0_RealK, i=1, npd_products-8), & ! CH3ONO2 : JPL 19-5
   237.0E-09_RealK,                 & ! C2H3 -> C2H2 + H   
  (0.0_RealK, i=1, npd_products-1), & ! C2H3 : Venot 2012-highest non zero qy
   720.0E-09_RealK,                 & ! C2H4 -> C2H2 + H2
   196.0E-09_RealK,                 & ! C2H4 -> C2H2 + H + H   
  (0.0_RealK, i=1, npd_products-2), & ! C2H4 : Huebner 92
   511.4E-09_RealK,                 & ! OH -> O(1D) + H  
  (0.0_RealK, i=1, npd_products-1), & ! OH : Huebner 1992
   229.0E-09_RealK,                 & ! HCO -> H + CO 
  (0.0_RealK, i=1, npd_products-1), & ! HCO : Venot 2012-highest non zero xsc
  453.0E-09_RealK,                  & ! N2O4 -> NO2 + NO2 
  (0.0_RealK, i=1, npd_products-1), & ! N2O4 : Venot 2012-highest non zero xsc
  224.0E-09_RealK,                  & ! C2N2 -> C2 + N2 
  (0.0_RealK, i=1, npd_products-1), & ! C2N2 : Venot 2012-highest non zero xsc
  290.0E-09_RealK,                  & ! N2H4 -> N2H3 + H
  (0.0_RealK, i=1, npd_products-1), & ! N2H4 : Venot 2012-highest non zero xsc
  398.0E-09_RealK,                  & ! N2O3 -> NO2 + NO
  (0.0_RealK, i=1, npd_products-1)  & ! N2O3 : Venot 2012-highest non zero xsc
  ], shape=[npd_products, npd_gases] )

! Unless otherwise stated, data comes from JPL publication No. 15-10:
! Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies
! Other references:
!  * JPL 19-5   : JPL publication No. 19-5
!  * Heubner 92 : Heubner et al (1992, p120) DOI: 10.1007/978-94-017-3023-5_1
!  * Becker 95 : Becker et al (1995) DOI: 10.1016/0301-0104(95)00114-4
!  * Zborowska 21 : Zborowska et al (2021) DOI: 10.1016/j.atmosenv.2021.118383
!  * IUPAC Task Group on Atmospheric Chemical Kinetic Data Evaluation:
!    * http://iupac.pole-ether.fr
!    * HNO3: PNOx2: https://iupac-aeris.ipsl.fr/datasheets/pdf/PNOx2.pdf
!    * PAN: P21: https://iupac-aeris.ipsl.fr/datasheets/pdf/P21.pdf
!  * Venot 2012 : Venot et al (2012) DOI: 10.1051/0004-6361/201219310
!  * Lias 70   : Lias et al (1970) DOI: 10.1063/1.1673226
!  * Chang 2020   : Chang et al (2020) DOI: 10.1039/D0SC01746A 

END MODULE gas_list_pcf
