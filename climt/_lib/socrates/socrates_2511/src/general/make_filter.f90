PROGRAM make_filter

  USE def_inst_flt, ONLY: StrFiltResp
  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE

  INTEGER :: i
  REAL (RealK) :: wmin, wmax
  INTEGER :: iu_flt = 35
!   Unit number for the filter function file
  CHARACTER  (LEN=20) :: ch_unit
!   Character string holding the definition of the unit
  CHARACTER  (LEN=80) :: filename="UV_index"
  TYPE (StrFiltResp) :: filter
!   Instrumental response function

  filter%satellite = "UV index"
  filter%channel   = "250 - 400nm"
  ch_unit          = "*MICRON"
  filter%n_pts     = 3001
  wmin = 0.25_RealK
  wmax = 0.4_RealK
  ALLOCATE(filter%wavenumber(filter%n_pts))
  ALLOCATE(filter%response(filter%n_pts))
  DO i = 1, filter%n_pts
    filter%wavenumber(i) = wmin + (wmax-wmin)*REAL(i-1,RealK)/(filter%n_pts-1)
    IF (filter%wavenumber(i) < 0.25) THEN
      filter%response(i) = 0.0_RealK
    ELSE IF (filter%wavenumber(i) < 0.298) THEN
      filter%response(i) = 1.0_RealK
    ELSE IF (filter%wavenumber(i) < 0.328) THEN
      filter%response(i) = 10.0_RealK**(0.094_RealK* &
        (298.0_RealK-filter%wavenumber(i)*1000.0_RealK))
    ELSE IF (filter%wavenumber(i) <= 0.4) THEN
      filter%response(i) = 10.0_RealK**(0.015_RealK* &
        (140.0_RealK-filter%wavenumber(i)*1000.0_RealK))
    ELSE
      filter%response(i) = 0.0_RealK
    END IF
  END DO
  filter%response = filter%response*40.0_RealK

  OPEN(unit=iu_flt, file=filename, status='unknown')
  WRITE(iu_flt, '(3(/a/),/,i0)') TRIM(filter%satellite), &
    TRIM(filter%channel), TRIM(ch_unit), filter%n_pts
  DO i = 1, filter%n_pts
    WRITE(iu_flt, '(2(1pe16.9))') filter%wavenumber(i),  filter%response(i)
  END DO
  CLOSE(iu_flt)

END PROGRAM make_filter
