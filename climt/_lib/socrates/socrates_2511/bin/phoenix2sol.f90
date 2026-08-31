PROGRAM phoenix2sol

  USE realtype_rd, ONLY: RealK
  USE def_solarspec, ONLY: StrSolarSpec
  USE rad_ccf, ONLY: solar_radius, astronomical_unit

  IMPLICIT NONE

  TYPE (StrSolarSpec) :: sol
  ! Stellar spectrum

  CHARACTER (LEN=256) :: line, input_file, output_file
  
  INTEGER :: ierr, ios
  ! Error flags
  INTEGER :: i
  ! Loop variables
  INTEGER :: iu_input, iu_output
  ! Unit number for the input and output solar spectrum data files
  REAL (RealK) :: scale_wv, scale_irr
  ! Scaling for wavelength and irradiance to correct units

  CALL get_command_argument(1, input_file)
  CALL get_command_argument(2, output_file)
  
  ! Read data from the input stellar spectrum.
  CALL get_free_unit(ierr, iu_input)
  OPEN(unit=iu_input, file=input_file, iostat=ios, status='UNKNOWN')
      
  ! Read first to find the number of points in the spectrum.
  Sol%n_points = 0
  DO
    READ(iu_input, '(A)', IOSTAT=ios) line
    IF (ios /= 0) THEN
      EXIT
    ELSE
      Sol%n_points=Sol%n_points+1
    ENDIF
  ENDDO

  ALLOCATE(Sol%wavelength(Sol%n_points))
  ALLOCATE(Sol%irrad(     Sol%n_points))

  ! Trappist1
  ! sol%t_effective = 2600.0_RealK
  ! sol%radius = 0.121 * solar_radius

  ! Proxima Centauri
  sol%t_effective = 3000.0_RealK
  sol%radius = 0.1542 * solar_radius

  scale_wv=1.0E-10_RealK ! Wavelength in Angstroms

  ! Convert from Erg/s/cm^2/A to W/m^2/m (x1.0E7)
  ! and from surface of star to 1AU.
  scale_irr=1.0E+07 * sol%radius**2 / astronomical_unit**2

  REWIND(iu_input)
  DO i=1, Sol%n_points
    READ(iu_input, *, IOSTAT=ios) Sol%wavelength(i), Sol%irrad(i)
    IF (ios /= 0) THEN
      print*, 'end of file'
      STOP
    END IF
  END DO

  CLOSE(iu_input)
  
  Sol%wavelength = Sol%wavelength * scale_wv
!  Sol%irrad = 10.0_RealK**(sol%irrad-8.0) * scale_irr
  Sol%irrad = sol%irrad * scale_irr
  sol%l_binned = .FALSE.

  CALL write_solar_spectrum(output_file, sol, ierr)
  
end program phoenix2sol
