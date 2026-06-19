! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read wavelengths for ice scattering from a database.
!
SUBROUTINE get_db_wavelengths &
 (iu_db_input, nd_wavelength, nd_size_scat, n_wavelength, wavelength, &
  n_rec_block, db_record, n_angle, ierr)
!
! Method:
!   The database is read to form a list of wavelengths present
!   and to establish the number of records and their locations.
!
!
!
! Modules used:
  USE def_std_io_icf
  USE error_pcf
  USE realtype_rd
  USE def_sct_db
!
!
  IMPLICIT NONE
!
!
!
! Dummy arguments
  INTEGER, Intent(In) :: nd_wavelength
!   Size allocated for array of wavelengths
  INTEGER, Intent(In) :: nd_size_scat
!   Size allocated for the number of scattering 
!   entries at each wavelength
!
  INTEGER, Intent(IN) :: iu_db_input
!   Unit number for input from the database
!
  INTEGER, Intent(Out) :: n_wavelength
!   Number of distinct wavelengths in the database
  REAL (RealK), Intent(Out), Dimension(nd_wavelength) :: &
    wavelength
!   Array of actual wavelengths (set to a fixed passed dimension)
  INTEGER, Intent(Out), Dimension(nd_wavelength) :: n_rec_block
!   Number of blocks of data for each wavelength
  INTEGER, Intent(Out), Dimension(nd_wavelength, nd_size_scat) :: &
    db_record
!   Pointers to the records for each wavelength and each block at
!   that wavelength
!
!
  INTEGER, Intent(Out) :: n_angle
!   Number of angles in each list of the phase function
!
  INTEGER, Intent(Out) :: ierr
!
!
!
! Local variables:
  INTEGER :: n_block
!   Number of blocks of scattering data
  INTEGER :: i_block
!   Loop variable   
  LOGICAL :: new_wavelength
!   Flag that states if a new wavelength
!   has been read from the database
  INTEGER :: j
!   Loop variable
  TYPE (str_sctdb_entry) :: ice_sct
!   Structure holding the ice scattering data
  INTEGER :: j_pt_wl
!   Pointer to the wavelength to which the data apply
!
!
!
! Read the total number of blocks from the unformatted database file.
  READ(iu_db_input, REC=1) n_block
!
! Read each block of data from the unformatted database file.
  i_block      = 1
  n_wavelength = 0
  n_rec_block  = 0
!
  process_block: DO

!   Read the data from the direct access file.
    READ(iu_db_input, REC=i_block+1) &
      ice_sct%dm, ice_sct%wavelength, ice_sct%csca, ice_sct%cext, &
      ice_sct%omega, ice_sct%asymm, ice_sct%n_angle
!
!
!   Test whether this is a new wavelength or close to an existing one.
    j=1
    DO
!
      IF (j > n_wavelength) THEN
!       We have exhausted the existing wavelengths, so this must be
!       a new wavelength.
        n_wavelength             = n_wavelength + 1
        wavelength(n_wavelength) = ice_sct%wavelength
        j_pt_wl                  = n_wavelength
        EXIT
      ELSE IF (ABS(wavelength(j)-ice_sct%wavelength) < &
        1.0E-5*MIN(wavelength(j), &
        REAL(ice_sct%wavelength, RealK)) ) THEN
!       This is an existing wavelength.
        j_pt_wl = j
        EXIT
      ELSE
!       Advance to the next wavelength
        j = j+1
      ENDIF
!
    ENDDO
!
!   Set up the locations for this block of data.
    n_rec_block(j_pt_wl)                     = n_rec_block(j_pt_wl)+1
    db_record(j_pt_wl, n_rec_block(j_pt_wl)) = i_block+1
!
!
!   Check if number of scattering angles is consistent 
!   throughout the records
    IF (i_block == 1) n_angle = ice_sct%n_angle
    IF ( ice_sct%n_angle /= n_angle .AND. i_block > 1 ) THEN
      WRITE(iu_err, '(/A, /A)') '*** Error: The number of ' // &
        'scattering angles is not consistent throughout ', &
        'the database records.'
      ierr = i_err_fatal
      RETURN
    ENDIF
!
!   Advance to the next block
    i_block = i_block+1
!   Stop at the last block.
    IF (i_block > n_block) EXIT
!
  ENDDO process_block
!
!
!
  RETURN
END SUBROUTINE get_db_wavelengths
