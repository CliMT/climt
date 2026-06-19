! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to write averaged scattering properties.
!
SUBROUTINE write_average_90 &
!
(iu_average, n_band, i_block, ss_data, &
 i_output_type, &
 absorption, scattering, n_phf_term, phase_fnc, &
 nd_band, nd_phase_term &
)
!
! Description:
!   Single scattering data averaged in frequency space are written to
! a file.
!
! Modules used
  USE realtype_rd
  USE scatter_pp_pcf
  USE def_s_scat_prop
  USE file_type_pcf
!
!
  IMPLICIT NONE
!
!
! Sizes of dummy arrays
  INTEGER, Intent(IN) :: nd_band
!   Size allocated for spectral bands
  INTEGER, Intent(IN) :: nd_phase_term
!   Size allocated for terms in the phase function
!
! Dummy arguments.
  INTEGER, Intent(IN) :: iu_average
!   Unit number for writing averaged data
  INTEGER, Intent(IN) :: n_band
!   Number of bands
  INTEGER, Intent(IN) :: i_block
!   Index of block
  TYPE  (StrSingleScatSzD), Intent(IN) :: ss_data
!   Single scattering data, including geometrical properties
  INTEGER, Intent(IN) :: i_output_type
!   Output type
  INTEGER, Intent(IN) :: n_phf_term
!   Number of terms in the phase function
  REAL  (RealK), Intent(IN), Dimension(nd_band) :: absorption
!   Absorption extinction
  REAL  (RealK), Intent(IN), Dimension(nd_band) :: scattering
!   Scattering extinction
  REAL  (RealK), Intent(IN), Dimension(nd_band, nd_phase_term) :: phase_fnc
!   Phase function
!
! Local variables.
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
!
!
!
  WRITE(iu_average, '(a13, i5)') '*FILE TYPE = ', i_output_type
  WRITE(iu_average, '(/,a)') &
    'Scattering properties averaged across each band'
  WRITE(iu_average, '(/, 3x, a8, i5)') 'Block = ', i_block
  WRITE(iu_average, '(3x, a20, i5, 3x, a20)') &
    'Type of scatterer = ', ss_data%i_scatter_type, ss_data%name_type
  IF (ss_data%i_scatter_type == ip_type_aerosol) THEN
    WRITE(iu_average, '(3x, a21, i5, 3x, a20)') &
      'Index of component = ', &
      ss_data%i_component, ss_data%name_component
  ENDIF
  IF ( i_output_type  ==  IT_file_ave_mie_dry ) THEN
    WRITE(iu_average, &
      '(/, 3x, a19, 1pe12.5, 1x, a4, /, ' // &
      ' 3x, a19, 1pe12.5, 1x, a2, /, 3x, a19, 1pe12.5)') &
        'Number density   = ', ss_data%number_density, 'm-3:', &
        'Character. Dim.  = ', ss_data%dim_char, 'm:', &
        'Volume fraction  = ', ss_data%vol_frac
  ELSE IF ( i_output_type  ==  IT_file_ave_phf_mie_dry ) THEN
    WRITE(iu_average, &
      '(/, 3x, a19, 1pe12.5, 1x, a4, /, ' // &
      ' 3x, a19, 1pe12.5, 1x, a2, /, 3x, a19, 1pe12.5)') &
        'Number density   = ', ss_data%number_density, 'm-3:', &
        'Character. Dim.  = ', ss_data%dim_char, 'm:', &
        'Volume Fraction  = ', ss_data%vol_frac
    WRITE(iu_average, '(/, 3x, a39, 1x, i3, /)') &
      'Number of terms in the phase function =', n_phf_term
  ELSE IF ( i_output_type  ==  IT_file_ave_mie_humid ) THEN
    WRITE(iu_average, '(/, 3x,a11,f7.5)') 'Humidity = ', ss_data%humidity
    WRITE(iu_average, &
      '(/, 3x, a19, 1pe12.5, 1x, a4, /, ' // &
      ' 3x, a19, 1pe12.5, 1x, a2, /, 3x, a36, 1pe12.5, 1x, a2, /, ' // &
      ' 3x, a19, 1pe12.5, /,  3x, a36, 1pe12.5)') &
        'Number density   = ', ss_data%number_density, 'm-3:', &
        'Effective radius = ', ss_data%dim_char, 'm:', &
        'Effective radius of dry particles = ' , ss_data%re_dry, 'm:', &
        'Volume fraction  = ', ss_data%vol_frac, &
        'Volume fraction of dry particles  = ', ss_data%vol_frac_dry
  ELSE IF ( i_output_type  ==  IT_file_ave_phf_mie_humid ) THEN
    WRITE(iu_average, '(/, 3x,a11,f7.5)') 'Humidity = ', ss_data%humidity
    WRITE(iu_average, '(/, 3x, a19, 1pe12.5, 1x, a4, / ' // &
      ' 3x, a19, 1pe12.5, 1x, a2, /, 3x, a36, 1pe12.5, 1x, a2, / ' // &
      ' 3x, a19, 1pe12.5, /,  3x, a36, 1pe12.5)') &
        'Number density   = ', ss_data%number_density, 'm-3:', &
        'Effective radius = ', ss_data%dim_char, 'm:', &
        'Effective radius of dry particles = ', &
        ss_data%re_dry, 'm:', &
        'Volume fraction  = ', ss_data%vol_frac, &
        'Volume fraction of dry particles  = ', &
        ss_data%vol_frac_dry
    WRITE(iu_average, '(/, 3x, a39, 1x, i3, /)') &
      'Number of terms in the phase function =', n_phf_term
  ELSE IF ( i_output_type  ==  IT_file_ave_scat_mass ) THEN
    WRITE(iu_average, &
      '(/, 3x, a21, 1pe12.5, 1x, a6, /, ' // &
      ' 3x, a21, 1pe12.5, 1x, a3, /, 3x, a21, 1pe12.5, 1x, a6)') &
        'Mass mixing ratio  = ', ss_data%mass_mixing_ratio, 'kg/kg:', &
        'Mean particle mass = ', ss_data%dim_char, 'kg:', &
        'Air density =        ', ss_data%air_density, 'kg/m3:'
    WRITE(iu_average, '(/, 3x, a39, 1x, i3, /)') &
      'Number of terms in the phase function =', n_phf_term
  ENDIF
  WRITE(iu_average, '(/, a4, 8x, a10, 10x, a10, 10x, a10)') &
    'Band', 'Absorption', 'Scattering', 'Phase fnc.'
  WRITE(iu_average, '(12x, a10, 10x, a10)') &
    '  (m-1)   ', '  (m-1)   '
  DO i=1, n_band
    WRITE(iu_average, &
      '(i5,3(4x, 1pe16.9), /, (5x, 3(4x, 1pe16.9)))') &
      i, absorption(i), scattering(i), &
      (phase_fnc(i, j), j=1, n_phf_term)
  ENDDO
  WRITE(iu_average, '(//)')

END SUBROUTINE write_average_90
