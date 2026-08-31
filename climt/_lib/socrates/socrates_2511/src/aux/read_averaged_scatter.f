! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a file of averaged scattering data.
!
! Method:
!   Straightforward.
!
!- ---------------------------------------------------------------------
      SUBROUTINE read_averaged_scatter(ierr, file_in
     &  , n_band, absorption, scattering, n_phf_term, phase_fnc
     &  , nd_band, nd_phase_term
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE file_type_pcf
      USE scatter_pp_pcf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, Intent(IN) ::
     &    nd_band
!           Size allocated for spectral bands
     &  , nd_phase_term
!           Size allocated for terms in the phase function
!            
!
!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
      CHARACTER !, Intent(IN)
     &    file_in*80
!           Name of file to read
      INTEGER, Intent(OUT) ::
     &    n_phf_term
!           Number of terms in the phase function
     &  , n_band
!           Number of spectral bands
      REAL  (RealK), Intent(OUT) ::
     &    absorption(nd_band)
!           Absorption
     &  , scattering(nd_band)
!           Scattering extinction
     &  , phase_fnc(nd_phase_term, nd_band)
!           Phase function 
!
!     Local variables.
      INTEGER
     &    iunit_scat_ave
!           Unit used to read the data
     &  , ios
!           I/O error flag
     &  , i_input_type
!           Type of input file
     &  , i_block
!           Number of block in file
     &  , i_scatter_type
!           Type of scatterer
     &  , i_component
!           Aerosol component
     &  , n_band_rec
!           Dummy number of bands (N_BAND itself cannot be used as this
!           variable may be reset to 0.)
     &  , k
!           Loop variable
      LOGICAL
     &    l_exist
!           Flag for existence of file
      REAL  (RealK) ::
     &    number_total
!           Total number of particles
     &  , radius_eff
!           Effective radius
     &  , volume_fraction
!           Volume fraction of spheres
!
      EXTERNAL
     &    get_free_unit
!
!
!
!     Determine a free unit to read the file.
      CALL get_free_unit(ierr, iunit_scat_ave)
      IF (ierr /= i_normal) RETURN
!
!     Check that the file exists.
      INQUIRE(file=file_in, exist=l_exist)
      IF (.NOT.l_exist) THEN
        WRITE(iu_err, '(3(/a))') 
     &    '*** Error: The file ', file_in, 'does not exist.'
        ierr=i_err_exist
        RETURN
      ENDIF
!
!     Open the file.
      OPEN(unit=iunit_scat_ave, file=file_in, iostat=ios, status='old')
      IF (ios /= 0) THEN
         WRITE(iu_err, '(3(/a))')
     &      '*** Error: The file ', file_in, 'could not be opened.'
         ierr=i_err_io
         RETURN
      ENDIF
!
!     Check the type of the input file.  
      READ(iunit_scat_ave, '(13x, i5)') i_input_type
      IF ( (i_input_type /= it_file_ave_mie_dry).AND.
     &     (i_input_type /= it_file_ave_phf_mie_dry) ) THEN
        WRITE(iu_err, '(3(/a))') 
     &    '*** Error: This file', file_in, 'is of an invalid type.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      READ(iunit_scat_ave, '(///, 11x, i5, /, 23x, i5)')
     &  i_block, i_scatter_type
      IF (i_scatter_type == IP_type_aerosol)
     &  READ(iunit_scat_ave, '(24x, i5)') i_component
!
      READ(iunit_scat_ave, '(/, 3(22x, 1pe12.5, /))')
     &  number_total, radius_eff, volume_fraction
!
      IF (i_input_type == it_file_ave_phf_mie_dry) THEN
        READ(iunit_scat_ave, '(43x, i3 , /)') n_phf_term
      ELSE
        n_phf_term=1
      ENDIF
      READ(iunit_scat_ave, '(//)')
!
      READ(iunit_scat_ave, '(i5)', iostat=ios) n_band_rec
      DO WHILE ( (n_band_rec > 0).AND.(ios == 0) )
        IF (n_band_rec > nd_band) THEN
           WRITE(iu_err, '(2(/a))')
     &        '*** Error: The size allocated for spectral bands '
     &        //'is too small for the file ', file_in
           ierr=i_err_fatal
           RETURN
        ENDIF
        backspace(iunit_scat_ave)
        READ(iunit_scat_ave, *)
     &    n_band, absorption(n_band), scattering(n_band)
     &    , (phase_fnc(k, n_band), k=1, n_phf_term)
        READ(iunit_scat_ave, '(i5)', iostat=ios) n_band_rec
      ENDDO
      READ(iunit_scat_ave, '(/)')
      CLOSE(iunit_scat_ave)
!
!
!
      RETURN
      END
