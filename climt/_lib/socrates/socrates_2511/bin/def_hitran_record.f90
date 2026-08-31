! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to declare a structure for a HITRAN record
!
MODULE def_hitran_record
!
! Description:
!
! This module contains the declaration of the structure
! currently used for elements of the HITRAN spectroscopic
! database.
!
!- End of header
!
! Modules used:
  USE realtype_rd
!
! Format for HITRAN Parameters, 1986 - 2001
  CHARACTER (LEN = 68), PARAMETER :: hitran_record_frmt86 = &
       "(i2,i1,f12.6,2(1pe10.3),2(0pf5.4),f10.4,f4.2,f8.6,2i3,2a9,3i1,3i2)"
!
  TYPE StrHitranRec86
!
    INTEGER           :: mol_num         ! 1=h20, 2=co2, 3=o3, etc
    INTEGER           :: iso_num         ! 1=most abundant, 2=second most, etc
    REAL  (RealK)     :: frequency       ! in cm-1
    REAL  (RealK)     :: intensity       ! in cm-1/(molecule cm-2) at 296K
    REAL  (RealK)     :: trans_momsq     ! in Debye squared
    REAL  (RealK)     :: air_broadhw     ! in cm-1/atm at 296K
    REAL  (RealK)     :: self_broadhw    ! in cm-1/atm at 296K
    REAL  (RealK)     :: lstate_energy   ! in cm-1
    REAL  (RealK)     :: air_broad_coeff ! unitless
    REAL  (RealK)     :: pshift          ! in cm-1/atm at 296K
    INTEGER           :: ustate_globalqi
    INTEGER           :: lstate_globalqi
    CHARACTER (LEN=9) :: ustate_localqi
    CHARACTER (LEN=9) :: lstate_localqi
    INTEGER           :: accindex_frequency
    INTEGER           :: accindex_intensity
    INTEGER           :: accindex_abhw
    INTEGER           :: tabrefs_frequency
    INTEGER           :: tabrefs_intensity
    INTEGER           :: tabrefs_abhw
!
  END TYPE StrHitranRec86

! Format for HITRAN Parameters, Editions after 2001
  CHARACTER (LEN = 74), PARAMETER :: hitran_record_frmt = &
    "(i2,i1,f12.6,2(1pe10.3),2(0pf5.4),f10.4,f4.2,f8.6,4a15,6i1,6i2,a1,2f7.1)"

! Bespoke format for HITRAN Parameters downloaded from hitran.org
  CHARACTER (LEN = 74), PARAMETER :: hitran_bespoke_frmt = &
    "(i3,i3,f12.6,2(1pe10.3),f6.4,f5.3,f10.4,f7.4,f9.6,4a15,6i1,6i2,a1,2f7.1)"

! Pointer to either standard or bespoke format depending on input options
! (default set here, options set in corr_k.f90)
  CHARACTER (LEN = 74) :: hitran_record_format = hitran_record_frmt

  TYPE StrHitranRec

    INTEGER           :: mol_num         ! 1=h20, 2=co2, 3=o3, etc
    INTEGER           :: iso_num         ! 1=most abundant, 2=second most, etc
    REAL  (RealK)     :: frequency       ! in cm-1
    REAL  (RealK)     :: intensity       ! in cm-1/(molecule cm-2) at 296K
    REAL  (RealK)     :: trans_momsq     ! in Debye squared
    REAL  (RealK)     :: air_broadhw     ! in cm-1/atm at 296K
    REAL  (RealK)     :: self_broadhw    ! in cm-1/atm at 296K
    REAL  (RealK)     :: lstate_energy   ! in cm-1
    REAL  (RealK)     :: air_broad_coeff ! unitless
    REAL  (RealK)     :: pshift          ! in cm-1/atm at 296K
    CHARACTER(LEN=15) :: ustate_globalqi
    CHARACTER(LEN=15) :: lstate_globalqi
    CHARACTER(LEN=15) :: ustate_localqi
    CHARACTER(LEN=15) :: lstate_localqi
    INTEGER           :: accindex_frequency
    INTEGER           :: accindex_intensity
    INTEGER           :: accindex_abhw
    INTEGER           :: accindex_sbhw
    INTEGER           :: accindex_abcoeff
    INTEGER           :: accindex_pshift   
    INTEGER           :: tabrefs_frequency
    INTEGER           :: tabrefs_intensity
    INTEGER           :: tabrefs_abhw     
    INTEGER           :: tabrefs_sbhw     
    INTEGER           :: tabrefs_abcoeff  
    INTEGER           :: tabrefs_pshift   
    CHARACTER(LEN=1)  :: line_mix_flag
    REAL  (RealK)     :: ustat_wgt
    REAL  (RealK)     :: lstat_wgt

  END TYPE StrHitranRec

! Format for HITRAN cross-section headers (2000)
  CHARACTER (LEN = 53), PARAMETER :: xsc_header_frmt = &
    "(a20,2f10.4,i7,f7.2,f6.1,1pe10.3,a5,a15,4x,a3,i3)"

  CHARACTER (LEN = 15), PARAMETER :: xsc_data_frmt = &
    "(10(1pe10.3))"

! Bespoke format for UV cross-section headers 
  CHARACTER (LEN = 53), PARAMETER :: uvxsc_header_frmt = &
    "(a20,2f16.4,i7,f7.2,f6.1,1pe16.9,1x,a5,a15,4x,a3,i3)"

  CHARACTER (LEN = 15), PARAMETER :: uvxsc_data_frmt = &
    "(10(1pe16.9))"

! Pointers to either xsc or uvxsc formats depending on input options
! (defaults set here, options set in corr_k.f90)
  CHARACTER (LEN = 53) :: xsc_header_format = xsc_header_frmt
  CHARACTER (LEN = 15) :: xsc_data_format = xsc_data_frmt
  
  TYPE StrXscHead
    CHARACTER(LEN=20)     :: chemical_symbol
    REAL (RealK)          :: wavenumber_min ! cm-1
    REAL (RealK)          :: wavenumber_max
    INTEGER               :: no_pts
    REAL (RealK)          :: temperature    ! K
    REAL (RealK)          :: pressure       ! Torr
    REAL (RealK)          :: max_xsc        ! cm2 molecule-1
    CHARACTER(LEN=5)      :: resolution     ! cm-1 (mA in UV)
    CHARACTER(LEN=15)     :: common_name
    CHARACTER(LEN=3)      :: broadening_gas
    INTEGER               :: re_no
  END TYPE StrXscHead

  TYPE StrXscRec
    TYPE (StrXscHead) :: head
    REAL (RealK), POINTER :: data(:) => NULL()
  END TYPE StrXscRec

! Format for HITRAN CIA headers
  CHARACTER (LEN = 43), PARAMETER :: cia_header_frmt = &
    "(a20,2f10.4,i7,f7.2,e7.3,1pe10.3,a5,a21,i3)"

  CHARACTER (LEN = 16), PARAMETER :: cia_data_frmt = &
    "(f10.4,1x,e10.3)"

  TYPE StrCIAHead
    CHARACTER(LEN=20)     :: chemical_symbol
    REAL (RealK)          :: wavenumber_min ! cm-1
    REAL (RealK)          :: wavenumber_max
    INTEGER               :: no_pts
    REAL (RealK)          :: temperature    ! K
    REAL (RealK)          :: pressure       ! Torr
    REAL (RealK)          :: max_xsc        ! cm2 molecule-1
    CHARACTER(LEN=5)      :: resolution     ! cm-1
    CHARACTER(LEN=21)     :: comments
    INTEGER               :: re_no
  END TYPE StrCIAHead

  TYPE StrCIARec
    TYPE (StrCIAHead) :: head
    REAL (RealK), ALLOCATABLE :: wavenumber(:)
    REAL (RealK), ALLOCATABLE :: data(:)
  END TYPE StrCIARec


END MODULE def_hitran_record
