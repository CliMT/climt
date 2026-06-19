! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to convert Fennelly and Torr (1992) data to QY files.
!
! This program works on a digitised version of table I from 
! https://doi.org/10.1016/0092-640X(92)90004-2 which can be found in
! examples/sp_uv/fennelly_and_torr_data.dat
! This contains photoionisation and photoabsorption cross-sections for
! O, N2, O2, and N broken down by the products of each branch. The ratio
! of these branches to the total cross-section for each species provides
! the branching ratio or quantum yield (QY) for each branch.
PROGRAM ft2qy

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE

  INTEGER :: ierr, ios
!   Error flag
  INTEGER :: iu_dat, iu_out
  !   File unit numbers
  INTEGER :: i, j
  INTEGER :: header_length = 7
  INTEGER :: data_length

  CHARACTER (LEN=256) :: infile
  CHARACTER (LEN=16) :: outfiles(15) = (/ &
    'O_1.qy          ', &
    'O_2.qy          ', &
    'O_3.qy          ', &
    'O_4.qy          ', &
    'O_5.qy          ', &
    'O.dat           ', &
    'N2.dat          ', &
    'N2_2.qy         ', &
    'N2_3.qy         ', &
    'N2_1.qy         ', &
    'O2_ft.dat       ', &
    'O2_6.qy         ', &
    'O2_7.qy         ', &
    'O2_2.qy         ', &
    'N.dat           '/)
  CHARACTER (LEN=16) :: outfile
  
  REAL (RealK), ALLOCATABLE :: in_dat(:, :)
  REAL (RealK), ALLOCATABLE :: in_data(:)
  REAL (RealK), ALLOCATABLE :: in_wn(:)


  CALL get_command_ARGUMENT(1, infile)

! Open the cross-section data file
  CALL get_free_unit(ierr, iu_dat)
  OPEN(UNIT=iu_dat, FILE=infile, IOSTAT=ierr, STATUS='OLD')

! Count the number of lines in the file
  data_length = 0
  DO
    READ(iu_dat, *, iostat=ios)
    IF (ios /= 0) EXIT
    data_length = data_length + 1
  END DO

  REWIND(iu_dat)

  ! Skip over header
  DO i=1, header_length
    READ(iu_dat, *)
    data_length = data_length - 1
  END DO

  ALLOCATE ( in_dat(16, data_length)  )
  ALLOCATE ( in_data(data_length)    )
  ALLOCATE ( in_wn(data_length)      )

! Read the data
  READ(iu_dat, *) in_dat
  CLOSE(iu_dat)

! Convert wavelength from Angstrom to nm
  in_wn = in_dat(1,:)*0.1_RealK

! Loop through the columns in the input data to determine
! cross-sections and quantum yields for output
  DO i=2, 16
    SELECT CASE (i)
    CASE(2,3,4,5,6)
      ! Determine quantum yield for each ionisation branch of atomic oxygen
      ! by dividing through by the total photoionisation cross-section
      DO j=1,data_length
        IF (in_dat(7,j) > 0.0_RealK) THEN
          in_data(j) = in_dat(i,j)/in_dat(7,j)
        ELSE
          in_data(j) = 0.0_RealK
        END IF
      END DO
    CASE(7,8)
      ! Total absorption cross-section for O and N2 converted to cm2
      in_data = in_dat(i,:)*1.0E-18_RealK
    CASE(9,10)
      ! Determine quantum yields for N2+ and N+ by dividing through by
      ! total absorption cross-section (or total ionisation cross-section
      ! if the data gives an inconsistently higher value)
      DO j=1,data_length
        IF (MAX(in_dat(8,j), in_dat(11,j)) > 0.0_RealK) THEN
          in_data(j) = in_dat(i,j)/MAX(in_dat(8,j), in_dat(11,j))
        ELSE
          in_data(j) = 0.0_RealK
        END IF
      END DO
    CASE(11)
      ! Determine quantum yield for N2 -> N + N as the only non-ionising
      ! branch for N2 absorption
      DO j=1,data_length
        IF (in_dat(8,j) > in_dat(11,j)) THEN
          in_data(j) = (in_dat(8,j)-in_dat(11,j))/in_dat(8,j)
        ELSE
          in_data(j) = 0.0_RealK
        END IF
      END DO
    CASE(12)
      ! Total absorption cross-section for O2 converted to cm2
      in_data = in_dat(i,:)*1.0E-18_RealK
    CASE(13,14)
      ! Determine quantum yields for O2+ and O+ by dividing through by
      ! total absorption cross-section (or total ionisation cross-section
      ! if the data gives an inconsistently higher value)
      DO j=1,data_length
        IF (MAX(in_dat(12,j), in_dat(15,j)) > 0.0_RealK) THEN
          in_data(j) = in_dat(i,j)/MAX(in_dat(12,j), in_dat(15,j))
        ELSE
          in_data(j) = 0.0_RealK
        END IF
      END DO
    CASE(15)
      ! At these wavelengths, the non-ionising cross-section leads to
      ! O2 -> O(3P) + O(1D) photoionisation
      DO j=1,data_length
        IF (in_dat(12,j) > in_dat(15,j)) THEN
          in_data(j) = (in_dat(12,j)-in_dat(15,j))/in_dat(12,j)
        ELSE
          in_data(j) = 0.0_RealK
        END IF
      END DO
    CASE(16)
      ! Total absorption (all ionising) cross-section for N converted to cm2
      in_data = in_dat(i,:)*1.0E-18_RealK
    END SELECT
    outfile = outfiles(i-1)
    ! Open a file for output
    CALL get_free_unit(ierr, iu_out)
    OPEN(UNIT=iu_out, FILE=outfile, IOSTAT=ierr, STATUS='UNKNOWN')
    SELECT CASE (i)
    CASE(7,8,12,16)
      ! Output total cross-section data files for each species
      DO j=1,data_length
        WRITE(iu_out,'(F10.5, 1PE14.5)') in_wn(j), in_data(j)
      END DO
    CASE DEFAULT
      ! Output quantum yield data files for each branch
      DO j=1,data_length
        WRITE(iu_out,'(2F10.5)') in_wn(j), in_data(j)
      END DO
    END SELECT
    CLOSE(iu_out)
  END DO

  DEALLOCATE ( in_wn      )
  DEALLOCATE ( in_data    )
  DEALLOCATE ( in_dat     )

END PROGRAM ft2qy
