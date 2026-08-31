! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set mapping between height and pressure.
!
! Method:
!   A reference profile is passed to the subroutine to define
!   the atmospheric state. This is checked to see that it
!   contains both pressure and height data. If not an error
!   is deemed to have occurrred. If it does, the reference 
!   state is set.
!
!- ---------------------------------------------------------------------
      SUBROUTINE set_state(ierr
     &  , n_level_profile, n_column_profile, i_data_type, profile
     &  , l_reference, n_level_ref, p_ref, z_ref
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
      USE input_head_pcf
      USE dimensions_field_ucf
      USE def_std_io_icf
      USE error_pcf
!
!
      IMPLICIT NONE
!
!
!
!     Dummy arguments.
      INTEGER, Intent(INOUT) ::
     &    ierr
!           Error flag
     &  , n_level_ref
!           Number of levels in reference
      INTEGER, Intent(IN) ::
     &    n_column_profile
!           Number of columns in profile
     &  , n_level_profile
!           Number of levels in profile
     &  , i_data_type(npd_data_column)
!           Types of data present
      LOGICAL, Intent(OUT) ::
     &    l_reference
!           True if reference state present
      REAL  (RealK), Intent(IN) ::
     &    profile(npd_layer+1, npd_data_column)
!           Proposed reference profile
      REAL  (RealK), Intent(OUT) ::
     &     p_ref(npd_layer+1)
!          Reference pressure
     &  , z_ref(npd_layer+1)
!           Reference height
!
!     Local variables.
      INTEGER
     &    i_column
!           Column index
     &  , i_column_pressure
!           Column containing pressure
     &  , i_column_height
!           Column containing height
     &  , i
!           Loop variable
!
!
!     Search for the columns containing the height and the pressure.
      i_column_pressure=0
      i_column_height=0
      i_column=0
1     i_column=i_column+1
      IF (i_data_type(i_column) == IP_pressure) THEN
        i_column_pressure=i_column
      ELSE IF (i_column < n_column_profile) THEN
        goto 1
      ELSE
        WRITE(iu_err, '(/a)') 
     &    '*** Error: Proposed reference profile contains '
     &    //'no pressure data.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
      i_column=0
2     i_column=i_column+1
      IF (i_data_type(i_column) == IP_height) THEN
        i_column_height=i_column
      ELSE IF (i_column < n_column_profile) THEN
        goto 2
      ELSE
        WRITE(iu_err, '(/a)') 
     &    '*** Error: Proposed reference profile contains '
     &    //'no height data.'
        ierr=i_err_fatal
        RETURN
      ENDIF
!
!     Data are present: pass the relevant values 
!     into the reference arrays.
      l_reference=.true.
      n_level_ref=n_level_profile
      DO i=1, n_level_profile
        p_ref(i)=profile(i, i_column_pressure)
        z_ref(i)=profile(i, i_column_height)
      ENDDO
!
!
!
      RETURN
      END
