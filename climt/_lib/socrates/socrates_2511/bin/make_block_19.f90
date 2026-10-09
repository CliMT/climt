! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to make spectral blocks of type 19.
!
! Method:
!   Initialy, a transparent grey fit is set for each gas.
!   A file is opened and an ESFT fit is read from the file.
!
!- ---------------------------------------------------------------------
SUBROUTINE make_block_19(Spectrum, ierr)

  USE realtype_rd
  USE def_spectrum
  USE rad_pcf
  USE file_type_pcf

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT), TARGET :: Spectrum
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local arguments.
  CHARACTER (LEN=80) :: line
!   Input line
  INTEGER :: iu_esft
!   Unit number for file of ESFT/k-term data
  INTEGER :: ios
!   IO status
  INTEGER :: i_input_type
!   Type of input file
  INTEGER :: i_index_gas_1, i_index_gas_2
!   Indices of continuum gases
  INTEGER :: i_band
!   Number of band
  INTEGER :: i_index
!   Number of continuum
  INTEGER :: i, j, k, it!, l, ip
!   Loop variables
  LOGICAL :: l_index_band(Spectrum%Dim%nd_band, Spectrum%ContGen%n_cont)
!   Absorbers present
  REAL (RealK) :: t_ref
!   Reference temperature (not used)
  INTEGER i_scale_k
!   Type of scaling applied to each k-term
  INTEGER :: i_scale_fnc
!   Type of scaling function
  INTEGER :: nd_k_term_cont_alloc
!   Previously allocated size for k-terms
  INTEGER, ALLOCATABLE :: arr_tmp_int_3d(:, :, :)
  REAL (RealK), ALLOCATABLE :: arr_tmp_real_3d(:, :, :)
  REAL (RealK), ALLOCATABLE :: arr_tmp_real_4d(:, :, :, :)
!   Temporary arrays used when resizing existing arrays

  LOGICAL, EXTERNAL :: non_blank
!   Function to detect blank lines

! Pointers to dimensions: used to shorten declarations later
  INTEGER, POINTER :: nd_band
!   Size allocated for spectral bands
  INTEGER, POINTER :: nd_k_term_cont
!   Size allocated for continuum k-terms
  INTEGER, POINTER :: nd_cont
!   Size allocated for continua
!   INTEGER, POINTER :: nd_scale_variable
! !   Size allocated for scaling variables

! Alias pointers to dimensions to the actual structure.
  nd_band            => Spectrum%Dim%nd_band
  nd_k_term_cont     => Spectrum%Dim%nd_k_term_cont
  nd_cont            => Spectrum%Dim%nd_cont

! If the block does not exist it is filled with grey null fits.
  IF (.NOT.Spectrum%Basic%l_present(19)) THEN
!   Allocate space for the arrays of k-terms.
    IF (ALLOCATED(Spectrum%ContGen%i_band_k_cont)) &
        DEALLOCATE(Spectrum%ContGen%i_band_k_cont)
    ALLOCATE(Spectrum%ContGen%i_band_k_cont(nd_band, nd_cont))
    IF (ALLOCATED(Spectrum%ContGen%i_scat_cont)) &
        DEALLOCATE(Spectrum%ContGen%i_scat_cont)
    ALLOCATE(Spectrum%ContGen%i_scat_cont(nd_k_term_cont, nd_band, nd_cont))
    IF (ALLOCATED(Spectrum%ContGen%i_cont_overlap_band)) &
        DEALLOCATE(Spectrum%ContGen%i_cont_overlap_band)
    ALLOCATE(Spectrum%ContGen%i_cont_overlap_band(nd_band, nd_cont))
    IF (ALLOCATED(Spectrum%ContGen%k_cont)) &
        DEALLOCATE(Spectrum%ContGen%k_cont)
    ALLOCATE(Spectrum%ContGen%k_cont(nd_k_term_cont, nd_band, nd_cont))
    IF (ALLOCATED(Spectrum%ContGen%w_cont)) &
        DEALLOCATE(Spectrum%ContGen%w_cont)
    ALLOCATE(Spectrum%ContGen%w_cont(nd_k_term_cont, nd_band, nd_cont))
    Spectrum%ContGen%i_scat_cont=0
    Spectrum%ContGen%i_cont_overlap_band=0
    DO i = 1, Spectrum%Basic%n_band
      DO j = 1, Spectrum%ContGen%n_band_cont(i)
        Spectrum%ContGen%i_band_k_cont(i, j) = 1
        Spectrum%ContGen%k_cont(1, i, j) = 0.0_RealK
        Spectrum%ContGen%w_cont(1, i, j) = 1.0_RealK
      ENDDO
    ENDDO
  ENDIF

! Obtain the band data from the prepared file of ESFT terms.
  CALL get_free_unit(ierr, iu_esft)
  CALL open_file_in(ierr, iu_esft, &
    'enter the name of the file of continuum esft data.')
  DO
    READ(iu_esft, '(A)', IOSTAT=ios) line
    IF (ios /= 0) THEN
      WRITE(*, '(/a)') '***error: file type not found.'
      ierr = i_err_fatal
      RETURN
    END IF
    IF (line(1:10) == '*FILE TYPE') THEN
      BACKSPACE iu_esft
      EXIT
    END IF
  END DO
  IF (ierr /= i_normal) RETURN

! Assemble the list of indexing numbers.
  DO i = 1, Spectrum%Basic%n_band
    l_index_band(i,:) = .FALSE.
    DO j = 1, Spectrum%ContGen%n_band_cont(i)
      l_index_band(i, Spectrum%ContGen%index_cont(j, i)) = .TRUE.
    END DO
  END DO

  outer: DO
    inner: DO
      READ(iu_esft, '(A)', IOSTAT=ios) line
      IF (ios < 0) EXIT outer
      IF (line(1:10) == '*FILE TYPE') THEN
        BACKSPACE iu_esft
        EXIT inner
      END IF
    END DO inner
    READ(iu_esft, '(15x, i5, //)', IOSTAT=ios) i_input_type
    IF (ios < 0) EXIT
    IF (i_input_type /= it_file_cont_gen_fit) THEN
      WRITE(*, '(/a)') &
        '***error: the esft data have an invalid file type.'
      ierr=i_err_fatal
      RETURN
    END IF
    READ(iu_esft, '(14x, i5, 18x, i5, 5x, i5)') &
      i_band, i_index_gas_1, i_index_gas_2

!   Find the index of this continuum in the spectral file
    i_index = -1
    DO j = 1, Spectrum%ContGen%n_cont
      IF ((Spectrum%ContGen%index_cont_gas_1(j) == i_index_gas_1 .AND. &
           Spectrum%ContGen%index_cont_gas_2(j) == i_index_gas_2) .OR. &
          (Spectrum%ContGen%index_cont_gas_1(j) == i_index_gas_2 .AND. &
           Spectrum%ContGen%index_cont_gas_2(j) == i_index_gas_1)) THEN
        i_index = j
        EXIT
      END IF
    END DO

!   Check that the continuum was found in the spectral file
    IF (i_index == -1) THEN
      WRITE(*, '(/A/)') '*** Error in subroutine make_block_19'
      WRITE(*,'(a, i4)') 'Continuum not found in spectral file.'
      ierr=i_err_fatal
      RETURN
    END IF

!   Find the position of this datum in the array of gases.
    IF (.NOT.l_index_band(i_band, i_index)) THEN
      WRITE(*, '(/a, i5)') 'Adding continuum to band', i_band
      l_index_band(i_band, i_index) = .TRUE.
      Spectrum%ContGen%n_band_cont(i_band) = &
        Spectrum%ContGen%n_band_cont(i_band) + 1
      Spectrum%ContGen%index_cont(Spectrum%ContGen%n_band_cont(i_band), &
        i_band) = i_index
    END IF

    READ(iu_esft, '(21x, 1pe10.3)') t_ref
    READ(iu_esft, '(//)')

!   Read over the transmission data.
    DO
      READ(iu_esft, '(a)') line
      IF (.NOT.non_blank(line)) EXIT
    END DO

!   Read number of ESFT/k-terms and scaling data
    READ(iu_esft, '(/, 23x, i5, 20x, i5, 21x, i5, //)') &
      Spectrum%ContGen%i_band_k_cont(i_band, i_index), &
      i_scale_k, i_scale_fnc

!   Only a look-up table in temperature is supported at the moment
    IF (i_scale_k /= ip_scale_term .OR. i_scale_fnc /= ip_scale_t_lookup) THEN
      WRITE(*, '(/A/)') '*** Error in subroutine make_block_19'
      WRITE(*,'(a, i4)') 'Invalid scaling for continuum.'
      ierr=i_err_fatal
      RETURN
    END IF

!   Resize arrays if the number of k-terms is greater than that allocated
    nd_k_term_cont_alloc = nd_k_term_cont
    IF (Spectrum%ContGen%i_band_k_cont(i_band, i_index) > &
        nd_k_term_cont_alloc) THEN
      nd_k_term_cont = Spectrum%ContGen%i_band_k_cont(i_band, i_index)
      ALLOCATE(arr_tmp_int_3d(nd_k_term_cont_alloc, nd_band, nd_cont))
      ALLOCATE(arr_tmp_real_3d(nd_k_term_cont_alloc, nd_band, nd_cont))

      arr_tmp_int_3d = Spectrum%ContGen%i_scat_cont
      DEALLOCATE(Spectrum%ContGen%i_scat_cont)
      ALLOCATE(Spectrum%ContGen%i_scat_cont(nd_k_term_cont, nd_band, nd_cont))
      Spectrum%ContGen%i_scat_cont(1:nd_k_term_cont_alloc, :, :) = &
          arr_tmp_int_3d
      Spectrum%ContGen%i_scat_cont(nd_k_term_cont_alloc+1:nd_k_term_cont, &
          :, :)=0

      arr_tmp_real_3d = Spectrum%ContGen%k_cont
      DEALLOCATE(Spectrum%ContGen%k_cont)
      ALLOCATE(Spectrum%ContGen%k_cont(nd_k_term_cont, nd_band, nd_cont))
      Spectrum%ContGen%k_cont(1:nd_k_term_cont_alloc, :, :) = arr_tmp_real_3d

      arr_tmp_real_3d = Spectrum%ContGen%w_cont
      DEALLOCATE(Spectrum%ContGen%w_cont)
      ALLOCATE(Spectrum%ContGen%w_cont(nd_k_term_cont, nd_band, nd_cont))
      Spectrum%ContGen%w_cont(1:nd_k_term_cont_alloc, :, :) = arr_tmp_real_3d

      DEALLOCATE(arr_tmp_int_3d)
      DEALLOCATE(arr_tmp_real_3d)
    END IF

!   Read optical depth = 1 ESFT/k-terms and weights
    DO k=1, Spectrum%ContGen%i_band_k_cont(i_band, i_index)
      READ(iu_esft, '(2(3x, 1pe16.9), (t39, 2(3x, 1pe16.9)))') &
        Spectrum%ContGen%k_cont(k, i_band, i_index), &
        Spectrum%ContGen%w_cont(k, i_band, i_index)
    END DO
    READ(iu_esft, '(/)')

!   Read in look-up table
    READ(iu_esft, '(14x, i4, 12x, i4)') Spectrum%ContGen%n_t_lookup_cont

    IF (Spectrum%ContGen%n_t_lookup_cont > Spectrum%Dim%nd_t_lookup_cont) THEN
      Spectrum%Dim%nd_t_lookup_cont = Spectrum%ContGen%n_t_lookup_cont
      IF (ALLOCATED(Spectrum%ContGen%t_lookup_cont)) &
          DEALLOCATE(Spectrum%ContGen%t_lookup_cont)
      ALLOCATE(Spectrum%ContGen%t_lookup_cont(Spectrum%Dim%nd_t_lookup_cont))
      IF (ALLOCATED(Spectrum%ContGen%k_lookup_cont)) &
          DEALLOCATE(Spectrum%ContGen%k_lookup_cont)
      ALLOCATE(Spectrum%ContGen%k_lookup_cont(Spectrum%Dim%nd_t_lookup_cont, &
          nd_k_term_cont, nd_cont, nd_band))

    ELSE IF (Spectrum%ContGen%i_band_k_cont(i_band, i_index) > &
             nd_k_term_cont_alloc) THEN
!     Resize look-up table as the number of k-terms is greater than that
!     allocated.
      ALLOCATE(arr_tmp_real_4d(Spectrum%Dim%nd_t_lookup_cont, &
          nd_k_term_cont_alloc, nd_cont, nd_band))
      arr_tmp_real_4d = Spectrum%ContGen%k_lookup_cont
      DEALLOCATE(Spectrum%ContGen%k_lookup_cont)
      ALLOCATE(Spectrum%ContGen%k_lookup_cont(Spectrum%Dim%nd_t_lookup_cont, &
        nd_k_term_cont, nd_cont, nd_band))
      Spectrum%ContGen%k_lookup_cont(:, 1:nd_k_term_cont_alloc, :, :) = &
        arr_tmp_real_4d
      DEALLOCATE(arr_tmp_real_4d)
    END IF
    
    READ(iu_esft, '(6(1PE13.6))', IOSTAT=ios) &
      (Spectrum%ContGen%t_lookup_cont(it), &
       it=1, Spectrum%ContGen%n_t_lookup_cont)
    IF (ios /= 0) THEN
      WRITE(*, '(/A/)') '*** Error in subroutine make_block_19'
      WRITE(*,'(a)') 'Error occurred reading T look-up table.'
      ierr=i_err_fatal
      RETURN
    END IF
!     Skip over the headers.
    READ(iu_esft, '(/)')
    DO k=1, Spectrum%ContGen%i_band_k_cont(i_band, i_index)
      READ(iu_esft, '(6(1PE13.6))', IOSTAT=ios) &
        (Spectrum%ContGen%k_lookup_cont(it,k,i_index,i_band), &
         it=1, Spectrum%ContGen%n_t_lookup_cont)
      IF (ios /= 0) THEN
        WRITE(*, '(/A/)') '*** Error in subroutine make_block_19'
        WRITE(*,'(a, 4i4)') 'Error occurred reading ESFT/k-terms.'
        ierr=i_err_fatal
        RETURN
      END IF
    END DO
  END DO outer

  Spectrum%Basic%l_present(19)=.TRUE.
  CLOSE(iu_esft)

END SUBROUTINE make_block_19
