! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read in HITRAN line data
!
SUBROUTINE read_hitran ( &
     iu_lbl,          & ! in
     iu_monitor,      & ! in
     molecule,        & ! in
     isotopes,        & ! in
     nd_lines,        & ! in
     wavenum_min,     & ! in
     wavenum_max,     & ! in
     number_of_lines, & ! out
     hitran_data )      ! out
  
  
! Description:
!
! Reads in HITRAN line data for required molecular species across
! required spectral interval into structure defined above
! 
! Method:
!
! Straightforward
!
! Modules used:
  USE realtype_rd, ONLY: RealK
  USE def_hitran_record
  USE hitran_cnst, ONLY: molecule_names, iso_abundance, number_isotopes
  USE gas_list_pcf, ONLY: npd_isotopes
  
  IMPLICIT NONE
  
  ! Subroutine arguments
  
  ! Scalar arguments with intent(in):
  INTEGER, INTENT(IN) :: nd_lines     ! Size allocated for number of lines
  INTEGER, INTENT(IN) :: molecule     ! molecule number
  INTEGER, INTENT(IN) :: isotopes(npd_isotopes) ! isotope numbers
  INTEGER, INTENT(IN) :: iu_lbl       ! unit number for I/O from LbL file
  INTEGER, INTENT(IN) :: iu_monitor   ! unit number for monitoring output
  REAL (RealK), INTENT(IN) :: wavenum_min  ! in cm-1
  REAL (RealK), INTENT(IN) :: wavenum_max  ! in cm-1
  
  ! Scalar arguments with intent(out):
  INTEGER, INTENT(OUT) :: number_of_lines
  
  ! Array arguments with intent(out):
  TYPE(StrHitranRec), DIMENSION(nd_lines), INTENT(OUT) :: hitran_data
  
  ! Local parameters:
  CHARACTER (LEN = *), PARAMETER :: routine_name = "read_hitran"
  
  INTEGER, PARAMETER :: input_unit  = 10
  
  ! Local scalars:
  TYPE(StrHitranRec) :: single_record
  
  
  INTEGER :: i, j, count_lines
  INTEGER :: io_status    ! Error code for file I/O
  
  REAL :: upper_cutoff
  REAL :: lower_cutoff

  REAL (RealK) :: abundance_scaling

!- End of header
  
  WRITE(*,"(a,a6)")          "Gas required: ", molecule_names(molecule)
  WRITE(iu_monitor,"(a,a6)") "Gas required: ", molecule_names(molecule)
  
  WRITE(*,"(a,f12.6,2x,f12.6)") &
    "Band limits: ", wavenum_min, wavenum_max
  WRITE(iu_monitor,"(a,f12.6,2x,f12.6)") &
    "Band limits: ", wavenum_min, wavenum_max
  
  
  
    WRITE(*,"(a)")          "Opened HITRAN data file "
    WRITE(iu_monitor,"(a)") "Opened HITRAN data file "

    upper_cutoff = wavenum_max + (wavenum_max-wavenum_min)/10.0
    lower_cutoff = wavenum_min - (wavenum_max-wavenum_min)/10.0

    number_of_lines = 0
    count_lines     = 0

    abundance_scaling = 0.0_RealK
    j = SUM(number_isotopes(1:molecule-1))
    DO i=1, npd_isotopes
      IF (isotopes(i) > 0) THEN
        abundance_scaling = abundance_scaling + iso_abundance(j+i)
      END IF
    END DO

    DO
      count_lines = count_lines + 1
      io_status   = 0
  
      READ(iu_lbl, hitran_record_format, IOSTAT = io_status) single_record
      ! Correct the isotope number from 0 to 10 (for CO2 in .par format data)
      IF (single_record % iso_num == 0) single_record % iso_num = 10
      IF (io_status == 0) THEN
        IF (single_record % frequency <= upper_cutoff) THEN
          IF (single_record % frequency >= lower_cutoff) THEN
            IF (single_record % mol_num == molecule) THEN
!             The strengths of lines in HITRAN are scaled by the
!             terrestrial abundance. For most applications we read
!             all lines for all isotopes. For some gases a subset of
!             isotopes are specified, in which case the line strengths
!             are rescaled so the isotope abundances sum to one.
              IF (isotopes(1) == 0 .OR. &
                  ANY(isotopes == single_record % iso_num)) THEN              
                number_of_lines              = number_of_lines + 1
                hitran_data(number_of_lines) = single_record
                IF (abundance_scaling > 0.0_RealK) THEN
                  hitran_data(number_of_lines)%intensity &
                    = single_record%intensity / abundance_scaling
                END IF
              ENDIF
            ENDIF
          ENDIF
        ELSE
          EXIT
        ENDIF
      ELSE IF (io_status < 0) THEN
        EXIT
      ELSE
        WRITE(*,"(a,i6)")          &
             "Error reading HITRAN file at line: ", count_lines 
        WRITE(iu_monitor,"(a,i6)") &
             "Error reading HITRAN file at line: ", count_lines 
        EXIT
      ENDIF
  
    ENDDO
  
END SUBROUTINE read_hitran
