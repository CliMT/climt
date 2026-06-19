! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to format a distribution of particle sizes.
!
! Method:
!   Directives are inserted into a file of raw data conrtaining
!   size spectra. This program reads the file and formats
!   the output for the scattering program.
!
!- ---------------------------------------------------------------------
      PROGRAM format_size
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE def_std_io_icf
      USE error_pcf
      USE def_data_in_icf
      USE dimensions_pp_ucf
!
!
      IMPLICIT NONE
!
!
!
!     Declaration of variables
      INTEGER
     &    ierr
!           Error flag
     &  , ios
!           I/O error flag
     &  , length_base
!           Length of base file name
     &  , n_distribution
!           Number of distributions
     &  , n_radius
!           Number of radii
     &  , n_size
!           Number of bin sizes
     &  , i
!           Loop variable
     &  , j
!           Loop variable
     &  , i_zero
!           Position of 0 in collating seq.
     &  , n1
!           Work variable
     &  , n2
!           Work variable
      LOGICAL
     &    l_set_radius
!           True if radii are set
     &  , l_set_size
!           True if sizes of bins are set
      REAL  (RealK) ::
     &    radius(npd_size)
!           Radii of bins
     &  , number(npd_size)
!           Numbers within bins
     &  , bin_size(npd_size)
!           Sizes of bins
      CHARACTER
     &    base_name*80
!           Base name for output files
     &  , file_name*80
!           Name of output file
     &  , title*80
!           Title of distribution
     &  , line*80
!           Line of input
     &  , shape*10
!           Shape of particle
!
      data ierr/i_normal/
!     Subroutines called:
      EXTERNAL
     &    open_file_in
!
!
!
!     Find the location of 0 in the compiler's collating sequence.
      i_zero=ichar('0')
!
!     Open input file.
      CALL open_file_in(ierr, iu_raw_in
     &  , 'Enter the name of the input file of size distributions.')
      IF (ierr == i_err_fatal) STOP
!
      WRITE(iu_stdout, '(a)')
     &  'Enter the base name of the output files.'
      READ(iu_stdin, '(a)') base_name
!     Find length of base name.
      j=0
      length_base=0
1     j=j+1
      IF (base_name(j:j) /= ' ') THEN
        length_base=length_base+1
        goto 1
      ENDIF
!
      WRITE(iu_stdout, '(a)')
     &  'Enter the shape of the particles.'
      READ(iu_stdin, '(a)') shape
!
      n_distribution=0
      l_set_radius=.false.
!     Read each line until a directive is found.
2     read(iu_raw_in, '(a)', end=3) line
      IF (line(1:11) == '*BIN_RADIUS') THEN
        READ(iu_raw_in, '(a)') title
        READ(iu_raw_in, *) n_radius
        IF (n_radius > npd_size) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: Number of radii exceeds permitted maximum.'
          STOP
        ENDIF
        READ(iu_raw_in, *) (radius(i), i=1, n_radius)
!       Convert from microns.
        DO i=1, n_radius
          radius(i)=1.0e-06_RealK*radius(i)
        ENDDO
        l_set_radius=.true.
!
      ELSE IF (line(1:11) == '*BIN_SIZE') THEN
        READ(iu_raw_in, '(a)') title
        READ(iu_raw_in, *) n_size
        IF (n_size > npd_size) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: Number of sizes exceeds permitted maximum.'
          STOP
        ENDIF
        READ(iu_raw_in, *) (bin_size(i), i=1, n_size)
!       Convert from microns.
        DO i=1, n_radius
          bin_size(i)=1.0e-06_RealK*bin_size(i)
        ENDDO
        l_set_size=.true.
!
!
      ELSE IF (line(1:11) == '*BIN_NUMBER') THEN
        IF (.NOT.l_set_radius) THEN
          WRITE(iu_err, '(/a)')
     &        '*** Error: Radii for the bins are not set: '
     &      , 'specify radii at the top of the file.'
          STOP
        ENDIF
        IF (.NOT.l_set_size) THEN
          WRITE(iu_err, '(/a)')
     &        '*** Error: Sizes for the bins are not set: '
     &      , 'specifiy sizes at the top of the file.'
          STOP
        ENDIF
        IF (n_radius /= n_size) THEN
          WRITE(iu_err, '(/a)')
     &      '*** Error: Number of radii and number of sizes '
     &      //'do not match.'
          STOP
        ENDIF
        READ(iu_raw_in, '(a)') title
        READ(iu_raw_in, *, iostat=ios) (number(i), i=1, n_radius)
        IF (ios /= 0) THEN
          WRITE(iu_err, '(/a, 1x, i5, a)')
     &        '*** Error: Distribution', n_distribution
     &      , ' could not be read.'
          STOP
        ENDIF
        n_distribution=n_distribution+1
        IF (n_distribution > 99) THEN
          WRITE(iu_err, '(/a)')
     &       '*** Error: Too many distributions for the file-'
     &       //'naming convention.'
          STOP
        ENDIF
!
!       Write out the size distribution in the form acceptable to the
!       Mie code.
        n1=n_distribution/10
        n2=n_distribution-10*n1
        file_name(1: length_base+3)=base_name(1: length_base)
     &    //'_'//char(n1+i_zero)//char(n2+i_zero)
        OPEN(file=file_name(1: length_base+3)
     &    , status='new', unit=iu_data_out)
        WRITE(iu_data_out, '(a)') 'Observational size distribution.'
        WRITE(iu_data_out, '(a)') title
        WRITE(iu_data_out, '(a)') shape
        WRITE(iu_data_out, '(a)') '*BEGIN_DATA'
        DO i=1, n_radius
!         Convert to actual number to a number per unit range of sizes.
          WRITE(iu_data_out, '(2(5x, 1pe13.6))')
     &      radius(i), number(i)/bin_size(i)
        ENDDO
        WRITE(iu_data_out, '(a)') '*END'
        CLOSE(iu_data_out)
!
      ELSE IF (line(1:1) == '*') THEN
        WRITE(iu_err, '(/a)')
     &    '*** Error: Unrecognized *-directive in input.'
        STOP
      ENDIF
      goto 2
!
3     CLOSE(iu_raw_in)
!
!
!
      STOP
      END
