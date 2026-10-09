! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to form composite array of pressures.
!
! Method:
!   Pressures in the background array are used as levels,
!   except in the region where the merged data are to be taken.
!
!- ---------------------------------------------------------------------
      SUBROUTINE merge_pressure(ierr
     &  , n_level_front, p_front, n_level_back, p_back
     &  , n_level_composite, p_composite
     &  , p_front_low, p_front_high
     &  , l_top_only
     &  )
!
!
!
!     Modules to set types of variables:
      USE realtype_rd
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
      INTEGER, Intent(IN) ::
     &    n_level_front
!           Number of levels merging array
     &  , n_level_back
!           Number of levels in background
      LOGICAL, Intent(IN) ::
     &    l_top_only
!           Logical to merge only upper levels
      REAL  (RealK), Intent(IN) ::
     &    p_front(npd_layer+1)
!           Merging pressures
     &  , p_back(npd_layer+1)
!           Background pressures
     &  , p_front_low
!           Low pressure of merging field
     &  , p_front_high
!           High pressure of merging field
!
      INTEGER, Intent(OUT) ::
     &    n_level_composite
!           Number of composite levels
      REAL  (RealK), Intent(OUT) ::
     &    p_composite(npd_layer+1)
!           Composite array
!
!     Local variables.
      INTEGER
     &    n_back
!           Temporary number of background levs
     &  , n_new
!           Putative size of merged array
     &  , i
!           Loop variable
!
!     Numerical precision:
      REAL  (RealK) ::
     &    tol_p
!           Tolerance for equality of pressures
!
!
!
!     Set the numerical precision
      tol_p=1.0e+02_RealK*epsilon(tol_p)
!
!     Pressures form the array to be merged are inserted into the
!     background array replacing old levels in this range.
      n_back=1
1     if (p_back(n_back).lt.p_front_low*(1.0_RealK-tol_p)) then
        n_back=n_back+1
        IF (n_back < n_level_back) goto 1
      ENDIF
      n_back=n_back-1
      DO i=1, n_back
        p_composite(i)=p_back(i)
      ENDDO
      n_level_composite=n_back
      n_new=n_level_composite+n_level_front
      IF (n_new > npd_layer+1) THEN
        WRITE(iu_err, '(/a)')
     &    '** Error: Merged array exceeds available size: '
     &    //'increase npd_layer and recompile.'
        ierr=i_err_fatal
        RETURN
      ENDIF
      DO i=1, n_level_front
        p_composite(n_back+i)=p_front(i)
      ENDDO
      n_level_composite=n_new
!
!     Return now unless the bottom of the background file is needed.
      IF (l_top_only) RETURN
!
!     P(N_BAKC) must satisfy the test below, but we need to increment
!     explicitly since n_back may still be 0.
      n_back=n_back+1
2     if (p_back(n_back).lt.p_front_high*(1.0_RealK+tol_p)) then
        IF (n_back < n_level_back) THEN
          n_back=n_back+1
          goto 2
        ELSE
!         No more points are needed.
          RETURN
        ENDIF
      ENDIF
      n_new=n_level_composite+n_level_back-n_back+1
      IF (n_new > npd_layer+1) THEN
        WRITE(iu_err, '(/a)')
     &    '*** Error: Merged array exceeds available size: '
     &    //'increase npd_layer and recompile.'
        ierr=i_err_fatal
        RETURN
      ENDIF
      DO i=1, n_level_back-n_back+1
        p_composite(n_level_composite+i)=p_back(n_back-1+i)
      ENDDO
      n_level_composite=n_new
!
!
!
      RETURN
      END
