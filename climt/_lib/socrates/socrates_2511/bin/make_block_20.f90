! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to make spectral blocks of type 20.
!
! Description:
!   This routine creates look-up tables of quantum yields for photolysis
!
!------------------------------------------------------------------------------
SUBROUTINE make_block_20(Sp, ierr)

  USE realtype_rd, ONLY: RealK
  USE def_spectrum, ONLY: StrSpecData
  USE def_std_io_icf, ONLY: iu_stdin, iu_err
  USE gas_list_pcf, ONLY: npd_gases, npd_products, &
    name_absorb, photol_products, threshold_wavelength
  USE rad_pcf, ONLY: i_err_fatal

  IMPLICIT NONE

  TYPE (StrSpecData), Intent(INOUT) :: Sp
!   Spectral file to be assigned
  INTEGER, Intent(INOUT) :: ierr
!   Error flag

! Local variables:
  INTEGER :: ios
!   I/O error flag
  INTEGER :: i, j, i_pe, i_qy, i_sub
!   Loop variable
  INTEGER :: i_type, type_index(npd_gases)
!   Types of gases
  INTEGER :: iu_qy
!   Unit number for the quantum yield data file
  INTEGER :: iu_pe
!   Unit number for the photoelectron factors data file
  INTEGER :: data_length_qy, header_length_qy
  INTEGER :: data_length_pe
  INTEGER :: n_qy
  INTEGER :: pathway_absorber(Sp%Dim%nd_pathway)
  INTEGER :: pathway_products(Sp%Dim%nd_pathway)
  INTEGER :: n_t_lookup_photol(Sp%Dim%nd_pathway), n_t_lookup
  INTEGER :: n_wl_lookup_photol(Sp%Dim%nd_pathway)
  INTEGER :: nd_t_lookup_photol, nd_wl_lookup_photol
  REAL(RealK) :: t_lookup_photol(Sp%Dim%nd_t_lookup_photol, Sp%Dim%nd_pathway)
  REAL(RealK) :: wl_lookup_photol(Sp%Dim%nd_wl_lookup_photol,Sp%Dim%nd_pathway)
  REAL(RealK) :: quantum_yield(Sp%Dim%nd_t_lookup_photol, &
                               Sp%Dim%nd_wl_lookup_photol, &
                               Sp%Dim%nd_pathway)
  REAL(RealK) :: threshold_wl(Sp%Dim%nd_pathway)
  LOGICAL :: l_thermalise(Sp%Dim%nd_pathway)

  REAL(RealK), ALLOCATABLE :: qy(:)
  REAL(RealK), ALLOCATABLE :: qy_data(:, :)
  REAL(RealK), ALLOCATABLE :: pe_data(:, :)
  REAL(RealK), ALLOCATABLE :: qy_unique(:, :)
  REAL(RealK), ALLOCATABLE :: qy_tmp(:, :)
  REAL(RealK), ALLOCATABLE :: qy_sub(:, :)
  REAL(RealK), ALLOCATABLE :: t_lookup(:)
  CHARACTER (LEN=256) :: qy_file
  CHARACTER (LEN=256) :: pe_file
  CHARACTER (LEN=1) :: char_in
  CHARACTER (LEN=256) :: line_in
  LOGICAL :: l_exist_qy
  LOGICAL :: l_exist_pe

  ! Add a photolysis pathway.
  IF (Sp%Photol%n_pathway > 0) THEN
    pathway_absorber   = Sp%Photol%pathway_absorber
    pathway_products   = Sp%Photol%pathway_products
    n_t_lookup_photol  = Sp%Photol%n_t_lookup_photol
    n_wl_lookup_photol = Sp%Photol%n_wl_lookup_photol
    t_lookup_photol    = Sp%Photol%t_lookup_photol
    wl_lookup_photol   = Sp%Photol%wl_lookup_photol
    quantum_yield      = Sp%Photol%quantum_yield
    threshold_wl       = Sp%Photol%threshold_wavelength
    l_thermalise       = Sp%Photol%l_thermalise
  END IF
  IF (ALLOCATED(Sp%Photol%pathway_absorber    )) &
     DEALLOCATE(Sp%Photol%pathway_absorber    )
  IF (ALLOCATED(Sp%Photol%pathway_products    )) &
     DEALLOCATE(Sp%Photol%pathway_products    )
  IF (ALLOCATED(Sp%Photol%n_t_lookup_photol   )) &
     DEALLOCATE(Sp%Photol%n_t_lookup_photol   )
  IF (ALLOCATED(Sp%Photol%n_wl_lookup_photol  )) &
     DEALLOCATE(Sp%Photol%n_wl_lookup_photol  )
  IF (ALLOCATED(Sp%Photol%t_lookup_photol     )) &
     DEALLOCATE(Sp%Photol%t_lookup_photol     )
  IF (ALLOCATED(Sp%Photol%wl_lookup_photol    )) &
     DEALLOCATE(Sp%Photol%wl_lookup_photol    )
  IF (ALLOCATED(Sp%Photol%quantum_yield       )) &
     DEALLOCATE(Sp%Photol%quantum_yield       )
  IF (ALLOCATED(Sp%Photol%threshold_wavelength)) &
     DEALLOCATE(Sp%Photol%threshold_wavelength)
  IF (ALLOCATED(Sp%Photol%l_thermalise        )) &
     DEALLOCATE(Sp%Photol%l_thermalise        )

  Sp%Photol%n_pathway = Sp%Photol%n_pathway + 1
  Sp%Dim%nd_pathway = Sp%Photol%n_pathway
  ALLOCATE(Sp%Photol%pathway_absorber(Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%pathway_products(Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%n_t_lookup_photol(Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%n_wl_lookup_photol(Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%threshold_wavelength(Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%l_thermalise(Sp%Dim%nd_pathway))

  IF (Sp%Photol%n_pathway > 1) THEN
    Sp%Photol%pathway_absorber(1:Sp%Photol%n_pathway-1) &
            = pathway_absorber(1:Sp%Photol%n_pathway-1)
    Sp%Photol%pathway_products(1:Sp%Photol%n_pathway-1) &
            = pathway_products(1:Sp%Photol%n_pathway-1)
    Sp%Photol%n_t_lookup_photol(1:Sp%Photol%n_pathway-1) &
            = n_t_lookup_photol(1:Sp%Photol%n_pathway-1)
    Sp%Photol%n_wl_lookup_photol(1:Sp%Photol%n_pathway-1) &
            = n_wl_lookup_photol(1:Sp%Photol%n_pathway-1)
    Sp%Photol%threshold_wavelength(1:Sp%Photol%n_pathway-1) &
            = threshold_wl(1:Sp%Photol%n_pathway-1)
    Sp%Photol%l_thermalise(1:Sp%Photol%n_pathway-1) &
            = l_thermalise(1:Sp%Photol%n_pathway-1)
  END IF

  ! Photolysis absorber
  DO i=1, Sp%Gas%n_absorb
    WRITE(*, '(i3, 7x, a)') &
      Sp%Gas%type_absorb(i), name_absorb(Sp%Gas%type_absorb(i))
    type_index(Sp%Gas%type_absorb(i)) = i
  END DO
  WRITE(*,'(A)') &
    'Enter the gas identifier for the absorbing species:'
  DO j=1, 3
    READ(iu_stdin, *, IOSTAT=ios) i_type
    IF (ios == 0 .AND. ANY(i_type == &
        Sp%Gas%type_absorb(1:Sp%Gas%n_absorb))) THEN
      Sp%Photol%pathway_absorber(Sp%Photol%n_pathway) = type_index(i_type)
      EXIT
    ELSE IF (j < 3) THEN
        WRITE(iu_err, '(A)') &
          '+++ Invalid input: Please re-enter.'
    ELSE
      WRITE(iu_err, '(A)') &
        '*** Error: Invalid gas identifiers.'
      ierr = i_err_fatal
      RETURN
    END IF
  END DO

  ! Products of photolysis
  WRITE(*,'(I3,7x,A)') 0, 'Unspecified'
  DO j=1, npd_products
    IF (LEN(TRIM(photol_products(j, i_type))) > 1) THEN
      WRITE(*,'(I3,7x,A)') j, photol_products(j, i_type)
    END IF
  END DO
  WRITE(*,'(A)') &
    'Enter the reaction products:'
  DO j=1, 3
    READ(iu_stdin, *, IOSTAT=ios) &
      Sp%Photol%pathway_products(Sp%Photol%n_pathway)
    IF (ios == 0 .AND. &
        Sp%Photol%pathway_products(Sp%Photol%n_pathway) <= npd_products) THEN
      EXIT
    ELSE IF (j < 3) THEN
        WRITE(iu_err, '(A)') &
          '+++ Invalid input: Please re-enter.'
    ELSE
      WRITE(iu_err, '(A)') &
        '*** Error: Invalid photolysis products.'
      ierr = i_err_fatal
      RETURN
    END IF
  END DO

  ! Thermalisation indicator
  WRITE(*,'(A)') &
    'Should the energy needed for photolysis be immediately thermalised? (y/n)'
  READ(iu_stdin, '(a)', IOSTAT=ios) char_in
  IF (ios == 0 .AND. (char_in == 'Y' .OR. char_in == 'y')) THEN
    Sp%Photol%l_thermalise(Sp%Photol%n_pathway) = .TRUE.
  ELSE
    Sp%Photol%l_thermalise(Sp%Photol%n_pathway) = .FALSE.
  END IF

  ! Threshold wavelength
  i = Sp%Photol%pathway_products(Sp%Photol%n_pathway)
  IF (i == 0) THEN
    WRITE(*,'(A)') &
      'Enter the threshold wavelength (m) for photolysis:'
    READ(iu_stdin, *, IOSTAT=ios) &
      Sp%Photol%threshold_wavelength(Sp%Photol%n_pathway)
  ELSE
    Sp%Photol%threshold_wavelength(Sp%Photol%n_pathway) &
      = threshold_wavelength(i, i_type)
  END IF

  ! Read in quantum yields from a file
  CALL get_free_unit(ios, iu_qy)
  WRITE(*,'(A)') &
    'Enter filename containing quantum yields:'
  READ(iu_stdin, '(A)', IOSTAT=ios) qy_file
  INQUIRE(FILE=qy_file, EXIST=l_exist_qy)

  IF (l_exist_qy .AND. qy_file /= '') THEN
    OPEN(unit=iu_qy, file=qy_file, iostat=ios, status='old')
    IF (ios /= 0) THEN
      WRITE(iu_err, '(A)') &
        '*** Error: cannot open file'
      ierr = i_err_fatal
      RETURN
    END IF

    ! Count the number of lines in the file
    data_length_qy = 0
    header_length_qy = 0
    DO
      READ(iu_qy, '(A)', iostat=ios) line_in
      IF (ios /= 0) EXIT
      IF (line_in(1:2) == '*T') THEN
        ! Where a *TEMPERATURES directive (or just *T) is present
        ! then the number of temperatures is read followed by a
        ! list of their values. Everything prior to the directive
        ! is considered part of the header.
        READ(line_in(INDEX(line_in,' '):), *) n_t_lookup
        ALLOCATE(t_lookup(n_t_lookup))
        READ(iu_qy, *, iostat=ios) t_lookup
        header_length_qy = data_length_qy + 2
        data_length_qy = 0
      ELSE
        data_length_qy = data_length_qy + 1
      END IF
    END DO
    REWIND(iu_qy)
    DO i=1, header_length_qy
      READ(iu_qy, *)
    END DO
    ! If no *TEMPERATURES directive is found, a single (arbitrary)
    ! temperature is assumed.
    IF (.NOT.ALLOCATED(t_lookup)) THEN
      n_t_lookup = 1
      ALLOCATE(t_lookup(n_t_lookup))
      t_lookup = 0.0_RealK
    END IF
    ! Read in the data (wavelength in nm, quantum yield)
    ALLOCATE(qy_data(1 + n_t_lookup, data_length_qy))
    READ(iu_qy, *, iostat=ios) qy_data
    CLOSE(iu_qy)

    ! Determine the upper wavelength bounds for unique QY values
    ALLOCATE(qy_unique(1 + n_t_lookup, 0:data_length_qy))
    ALLOCATE(qy(n_t_lookup))
    qy = -99.0_RealK
    n_qy = 0
    DO i=1, data_length_qy
      IF (ANY(qy_data(2:,i) /= qy)) THEN
        qy = qy_data(2:,i)
        n_qy = n_qy+1
        qy_unique(2:,n_qy)=qy
      END IF
      ! Upper wavelength bound at mid-point to next data point
      ! (converted from nm to m)
      IF (i == data_length_qy) THEN
        qy_unique(1,n_qy) &
          =(1.5_RealK*qy_data(1,i)-0.5_RealK*qy_data(1,i-1))*1.0E-09_RealK
      ELSE
        qy_unique(1,n_qy)=(qy_data(1,i)+qy_data(1,i+1))*0.5E-09_RealK
      END IF
      IF (ANY(qy > 0.0_RealK)) THEN
        ! Adjust threshold wavelength for last QY > 0
        Sp%Photol%threshold_wavelength(Sp%Photol%n_pathway) &
          = MAX( Sp%Photol%threshold_wavelength(Sp%Photol%n_pathway), &
                 qy_data(1,i)*1.0E-09_RealK )
      END IF
    END DO
    IF (ALL(qy == 0.0_RealK)) THEN
      n_qy = n_qy-1
    END IF
    DEALLOCATE(qy)
    DEALLOCATE(qy_data)
    ! Quantum yield is constrained to zero beyond threshold
    qy_unique(1,n_qy) = MIN(qy_unique(1,n_qy), &
      Sp%Photol%threshold_wavelength(Sp%Photol%n_pathway))
  ELSE
    n_qy=1
    ALLOCATE(qy_unique(2,0:n_qy))
    qy_unique(1,1)=Sp%Photol%threshold_wavelength(Sp%Photol%n_pathway)
    qy_unique(2,1)=1.0_RealK
    n_t_lookup = 1
    ALLOCATE(t_lookup(n_t_lookup))
    t_lookup = 0.0_RealK
  END IF

  ! Read in photoelectron factors from a file
  CALL get_free_unit(ios, iu_pe)
  WRITE(*,'(A)') &
    'Enter filename containing photoelectron factors:'
  READ(iu_stdin, '(A)', IOSTAT=ios) pe_file
  INQUIRE(FILE=pe_file, EXIST=l_exist_pe)
  IF (l_exist_pe .AND. pe_file /= '') THEN
    OPEN(unit=iu_pe, file=pe_file, iostat=ios, status='old')
    IF (ios /= 0) THEN
      WRITE(iu_err, '(A)') &
        '*** Error: cannot open file'
      ierr = i_err_fatal
      RETURN
    END IF

    ! Count the number of lines in the file
    data_length_pe = 0
    DO
      READ(iu_pe, *, iostat=ios)
      IF (ios /= 0) EXIT
      data_length_pe = data_length_pe + 1
    END DO
    REWIND(iu_pe)
    ! Read in the data (wavelength (nm) min, max, photoelectron factor)
    ALLOCATE(pe_data(3,data_length_pe))
    READ(iu_pe, *, iostat=ios) pe_data
    pe_data(1:2,:)=pe_data(1:2,:)*1.0E-09_RealK
    CLOSE(iu_pe)

    ! Reallocate qy_unique to hold extra values
    ALLOCATE(qy_tmp(1 + n_t_lookup, n_qy))
    qy_tmp=qy_unique(:,1:n_qy)
    DEALLOCATE(qy_unique)
    ALLOCATE(qy_unique(1 + n_t_lookup, 0:n_qy+data_length_pe+1))
    i_qy=1
    i_pe=1
    i=1
    DO
      ! Merge in photoelectron factors to quantum yield data
      ! (uncomment the print statements if debugging is required)
      IF (i_qy > n_qy .AND. i_pe > data_length_pe) THEN
        n_qy = i-1
!        PRINT*, '0:',n_qy
        EXIT
      ELSE IF (i_pe > data_length_pe) THEN
        qy_unique(1,i)=qy_tmp(1,i_qy)
        qy_unique(2:,i)=qy_tmp(2:,i_qy)
        i=i+1
        i_qy=i_qy+1
!        PRINT*, '1:',i, i_qy, i_pe, qy_unique(1,i-1), qy_unique(2,i-1)
      ELSE IF (i_qy > n_qy) THEN
        qy_unique(1,i)=pe_data(2,i_pe)
        qy_unique(2:,i)=pe_data(3,i_pe)
        i=i+1
        i_pe=i_pe+1        
!        PRINT*, '2:',i, i_qy, i_pe, qy_unique(1,i-1), qy_unique(2,i-1)
      ELSE IF (qy_tmp(1,i_qy) <= pe_data(2,i_pe)) THEN
        qy_unique(1,i)=qy_tmp(1,i_qy)
        qy_unique(2:,i)=qy_tmp(2:,i_qy)+pe_data(3,i_pe)
        i=i+1
        i_qy=i_qy+1
        IF (qy_tmp(1,i_qy) == pe_data(2,i_pe)) i_pe=i_pe+1
!        PRINT*, '3:',i, i_qy, i_pe, qy_unique(1,i-1), qy_unique(2,i-1)
      ELSE
        qy_unique(1,i)=pe_data(2,i_pe)
        qy_unique(2:,i)=qy_tmp(2:,i_qy)+pe_data(3,i_pe)
        i=i+1
        i_pe=i_pe+1
!        PRINT*, '4:',i, i_qy, i_pe, qy_unique(1,i-1), qy_unique(2,i-1)
      END IF
    END DO
    DEALLOCATE(qy_tmp)
    DEALLOCATE(pe_data)
  END IF

  IF (Sp%Basic%l_present(17)) THEN
    ! Mean Quantum Yield values across the sub-bands.
    ALLOCATE(qy_sub(n_t_lookup, 0:Sp%Var%n_sub_band))
    qy_sub(:, :) = 0.0_RealK
    qy_unique(1,0)=0.0_RealK
    i=1
    outer: DO i_sub=1, Sp%Var%n_sub_band
      inner: DO
        IF (qy_unique(1,i) > sp%var%wavelength_sub_band(1, i_sub)) THEN
          ! Add QY value to weighted mean
          qy_sub(:, i_sub) = qy_sub(:, i_sub) + qy_unique(2:, i) &
             * ( 1.0_RealK / MAX(qy_unique(1,i-1), &
                                 sp%var%wavelength_sub_band(1, i_sub)) &
               - 1.0_RealK / MIN(qy_unique(1,i), &
                                 sp%var%wavelength_sub_band(2, i_sub)) ) &
             / ( 1.0_RealK / sp%var%wavelength_sub_band(1, i_sub) &
               - 1.0_RealK / sp%var%wavelength_sub_band(2, i_sub) )
        END IF
        IF (qy_unique(1,i) > sp%var%wavelength_sub_band(2, i_sub)) EXIT inner
        IF (i == n_qy) EXIT outer
        i = i+1
      END DO inner
    END DO outer
    ! Reduce to unique QY values
    i_qy=0
    qy_sub(:, 0)=-99.9_RealK
    DO i_sub=1, Sp%Var%n_sub_band
      IF (ANY(ABS(qy_sub(:, i_sub)-qy_sub(:, i_sub-1)) > 1.0E-10_RealK)) THEN
        ! If QY is different assign a new unique value
        i_qy=i_qy+1
        n_qy=UBOUND(qy_unique,2)
        IF (i_qy > n_qy) THEN
          ALLOCATE(qy_tmp(1 + n_t_lookup, n_qy))
          qy_tmp=qy_unique(:,1:n_qy)
          DEALLOCATE(qy_unique)
          ALLOCATE(qy_unique(1 + n_t_lookup, i_qy))
          qy_unique(:,1:n_qy)=qy_tmp
          DEALLOCATE(qy_tmp)
        END IF
        qy_unique(1, i_qy) = sp%var%wavelength_sub_band(2, i_sub)
        qy_unique(2:, i_qy) = qy_sub(:, i_sub)
      ELSE
        ! If QY same as previous band, just reset the limiting wavelength
        qy_unique(1, i_qy) = sp%var%wavelength_sub_band(2, i_sub)
      END IF
    END DO
    n_qy=i_qy
    IF (ALL(qy_unique(2:, i_qy) <= 1.0E-10_RealK)) THEN
      n_qy = n_qy-1
    END IF
    DEALLOCATE(qy_sub)
  END IF

  Sp%Photol%n_t_lookup_photol(Sp%Photol%n_pathway) = n_t_lookup
  Sp%Photol%n_wl_lookup_photol(Sp%Photol%n_pathway) = n_qy

  nd_t_lookup_photol  = MAXVAL(Sp%Photol%n_t_lookup_photol)
  nd_wl_lookup_photol = MAXVAL(Sp%Photol%n_wl_lookup_photol)

  ALLOCATE(Sp%Photol%t_lookup_photol(nd_t_lookup_photol, &
                                     Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%wl_lookup_photol(nd_wl_lookup_photol, &
                                      Sp%Dim%nd_pathway))
  ALLOCATE(Sp%Photol%quantum_yield(nd_t_lookup_photol, &
                                   nd_wl_lookup_photol, &
                                   Sp%Dim%nd_pathway))
  IF (Sp%Photol%n_pathway > 1) THEN
    Sp%Photol%t_lookup_photol(1:Sp%Dim%nd_t_lookup_photol, &
                              1:Sp%Photol%n_pathway-1) &
      = t_lookup_photol(1:Sp%Dim%nd_t_lookup_photol, &
                        1:Sp%Photol%n_pathway-1)
    Sp%Photol%wl_lookup_photol(1:Sp%Dim%nd_wl_lookup_photol, &
                               1:Sp%Photol%n_pathway-1) &
      = wl_lookup_photol(1:Sp%Dim%nd_wl_lookup_photol, &
                         1:Sp%Photol%n_pathway-1)
    Sp%Photol%quantum_yield(1:Sp%Dim%nd_t_lookup_photol, &
                            1:Sp%Dim%nd_wl_lookup_photol, &
                            1:Sp%Photol%n_pathway-1) &
      = quantum_yield(1:Sp%Dim%nd_t_lookup_photol, &
                      1:Sp%Dim%nd_wl_lookup_photol, &
                      1:Sp%Photol%n_pathway-1)
  END IF
  Sp%Photol%t_lookup_photol(1:n_t_lookup, Sp%Photol%n_pathway) &
    = t_lookup
  DEALLOCATE(t_lookup)
  Sp%Photol%wl_lookup_photol(1:n_qy, Sp%Photol%n_pathway) &
    = qy_unique(1, 1:n_qy)
  Sp%Photol%quantum_yield(1:n_t_lookup, 1:n_qy, Sp%Photol%n_pathway) &
    = qy_unique(2:, 1:n_qy)
  DEALLOCATE(qy_unique)

  Sp%Dim%nd_t_lookup_photol = nd_t_lookup_photol
  Sp%Dim%nd_wl_lookup_photol = nd_wl_lookup_photol
  Sp%Basic%l_present(20)=.TRUE.

END SUBROUTINE make_block_20
