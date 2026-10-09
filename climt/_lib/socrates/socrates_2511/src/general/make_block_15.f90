! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to make spectral blocks of type 15.
!
! Method:
!   Aerosol optical properties are calculated for the AOD wavelengths
!   using a dummy spectrum structure which is passed to make_block_11.
!   The optical properties are therefore calculated in a consistent
!   manner with those in block 11.
!
!------------------------------------------------------------------------------
SUBROUTINE make_block_15(Sp, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT), TARGET :: Sp
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

  INTEGER :: j, k
!   Loop variables
  TYPE(StrSpecData) :: aod_sp
!   Dummy spectral configuration to hold aerosol optical depths

! Aerosol types which gather together aerosol components
  INTEGER, PARAMETER :: ip_type_sulphate = 1
  INTEGER, PARAMETER :: ip_type_dust     = 2
  INTEGER, PARAMETER :: ip_type_seasalt  = 3
  INTEGER, PARAMETER :: ip_type_soot     = 4
  INTEGER, PARAMETER :: ip_type_biomass  = 5
  INTEGER, PARAMETER :: ip_type_biogenic = 6
  INTEGER, PARAMETER :: ip_type_ocff     = 7
  INTEGER, PARAMETER :: ip_type_delta    = 8
  INTEGER, PARAMETER :: ip_type_nitrate  = 9
  INTEGER, PARAMETER :: ip_type_twobdust = 10

  INTEGER, PARAMETER :: aod_type_component(npd_aerosol_component) = &
    (/ ip_type_sulphate, &       ! "Water soluble       "
       ip_type_dust,     &       ! "Dust-like           "
       ip_type_seasalt,  &       ! "Oceanic             "
       ip_type_soot,     &       ! "Soot                "
       0,                &       ! "Volcanic Ash        "
       ip_type_sulphate, &       ! "Sulphuric Acid      "
       ip_type_sulphate, &       ! "Ammonium Sulphate   "
       0,                &       ! "Uncharacterized     "
       ip_type_dust,     &       ! "Saharan Dust        "
       ip_type_sulphate, &       ! "Accum. Sulphate     "
       ip_type_sulphate, &       ! "Aitken Sulphate     "
       ip_type_soot,     &       ! "Fresh Soot          "
       ip_type_soot,     &       ! "Aged Soot           "
       ip_type_seasalt,  &       ! "Generic NaCl        "
       ip_type_seasalt,  &       ! "NaCl film mode      "
       ip_type_seasalt,  &       ! "NaCl jet mode       "
       ip_type_dust,     &       ! "Dust Division 1     "
       ip_type_dust,     &       ! "Dust Division 2     "
       ip_type_dust,     &       ! "Dust Division 3     "
       ip_type_dust,     &       ! "Dust Division 4     "
       ip_type_dust,     &       ! "Dust Division 5     "
       ip_type_dust,     &       ! "Dust Division 6     "
       ip_type_biomass,  &       ! "Biomass Division 1  "
       ip_type_biomass,  &       ! "Biomass Division 2  "
       ip_type_biogenic, &       ! "Biogenic            "
       ip_type_ocff,     &       ! "Fresh fossil-fuel OC"
       ip_type_ocff,     &       ! "Aged fossil-fuel OC "
       ip_type_delta,    &       ! "Delta aerosol       "
       0,                &       ! "Murk                "
       ip_type_nitrate,  &       ! "Ammonium nitrate    "
       ip_type_twobdust, &       ! "Two-bin Dust Div 1  "
       ip_type_twobdust  /)      ! "Two-bin Dust Div 2  "
  

  IF (.NOT.Sp%Basic%l_present(15)) THEN
!   AOD wavelengths are currently hardwired
    Sp%Aerosol%n_aod_wavel = 6
    Sp%Dim%nd_aod_wavel = Sp%Aerosol%n_aod_wavel
    IF (ALLOCATED(Sp%Aerosol%i_aod_type)) &
       DEALLOCATE(Sp%Aerosol%i_aod_type)
    ALLOCATE(Sp%Aerosol%i_aod_type(Sp%Dim%nd_aerosol_species))
    IF (ALLOCATED(Sp%Aerosol%aod_wavel)) &
       DEALLOCATE(Sp%Aerosol%aod_wavel)
    ALLOCATE(Sp%Aerosol%aod_wavel(Sp%Dim%nd_aod_wavel))
    IF (ALLOCATED(Sp%Aerosol%aod_abs)) &
       DEALLOCATE(Sp%Aerosol%aod_abs)
    ALLOCATE(Sp%Aerosol%aod_abs(Sp%Dim%nd_humidity, &
      Sp%Dim%nd_aerosol_species, Sp%Dim%nd_aod_wavel))
    IF (ALLOCATED(Sp%Aerosol%aod_scat)) &
       DEALLOCATE(Sp%Aerosol%aod_scat)
    ALLOCATE(Sp%Aerosol%aod_scat(Sp%Dim%nd_humidity, &
      Sp%Dim%nd_aerosol_species, Sp%Dim%nd_aod_wavel))
    Sp%Aerosol%aod_wavel(1:6) = 1.0e-9 * &
      (/ 380.0, 440.0, 550.0, 670.0, 865.0, 1020.0 /)
  END IF

! Copy the spectrum structure to include all the relevant aerosol species
  aod_sp = Sp

! Change the number of bands to the number of AOD wavelengths
  aod_sp%basic%n_band = Sp%Aerosol%n_aod_wavel
  aod_sp%dim%nd_band = aod_sp%basic%n_band

! Remove current aerosol data
  aod_sp%Aerosol%l_aero_spec = .FALSE.

  CALL make_block_11(aod_sp, ierr)

! Move the block 11 data to block 15
  DO k = 1, Sp%Aerosol%n_aod_wavel
    DO j = 1, Sp%Aerosol%n_aerosol
      IF (aod_sp%Aerosol%l_aero_spec(j)) THEN
        Sp%Aerosol%aod_abs(:,j,k)  = aod_sp%Aerosol%abs(:,j,k)
        Sp%Aerosol%aod_scat(:,j,k) = aod_sp%Aerosol%scat(:,j,k)
      END IF
    END DO
  END DO

! Set the AOD type
  Sp%Aerosol%i_aod_type = aod_type_component(Sp%Aerosol%type_aerosol)

  Sp%Basic%l_present(15)=.TRUE.

END SUBROUTINE make_block_15
