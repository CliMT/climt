! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to split cloud into maximally overlapped C/S.
!
! Method:
!   Convective cloud is left-justified in the grid-box while
!   stratiform cloud is right-justified.
!
!- ---------------------------------------------------------------------
MODULE cloud_maxcs_split_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CLOUD_MAXCS_SPLIT_MOD'
CONTAINS
SUBROUTINE cloud_maxcs_split(ierr, n_profile, n_layer, n_cloud_top      &
    , w_cloud, frac_cloud                                               &
    , n_cloud_type                                                      &
    , n_column_cld, n_column_slv, list_column_slv                       &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
    , nd_profile, nd_layer, id_ct, nd_column, nd_cloud_type)


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, i_normal, ip_cloud_type_conv,         &
                     ip_cloud_type_strat
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength
  USE shell_sort_mod, ONLY: shell_sort

  IMPLICIT NONE


! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_column                                                         &
!       Size allocated for columns at a point
    , nd_cloud_type                                                     &
!       Size allocated for columns at a point
    , id_ct
!       Topmost allocated cloudy layer

! Dummy variables
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_cloud_type
!       Number of types of cloud
  REAL (RealK), INTENT(IN) ::                                           &
      w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Amount of cloud
    , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)
!       Fractions of different types of cloud

  INTEGER, INTENT(OUT) ::                                               &
      n_column_cld(nd_profile)                                          &
!       Number of columns in each profile (including those of
!       zero width)
    , n_column_slv(nd_profile)                                          &
!       Number of columns to be solved in each profile
    , list_column_slv(nd_profile, nd_column)                            &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(nd_profile, nd_column)                              &
!       Layer in the current column to change
    , i_clm_cld_typ(nd_profile, nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK), INTENT(OUT) ::                                          &
      area_column(nd_profile, nd_column)
!       Area of each column


! Local variables
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , ii                                                                &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , n_cld_layer                                                       &
!       Number of cloudy layers
    , ptr_st                                                            &
!       Pointer to stratiform cloud in arrays
    , ptr_cnv                                                           &
!       Pointer to convective cloud in arrays
    , key_st(nd_layer+1-id_ct)                                          &
!       Pointers to layers listing left edge of stratiform cloud
!       in increasing order
    , key_cnv(nd_layer+1-id_ct)                                         &
!       Pointers to layers listing right edge of convective cloud
!       in increasing order
    , i_key_cnv                                                         &
!       Current pointer to convective cloud
    , i_key_st
!       Current pointer to stratiform cloud
  REAL (RealK) ::                                                       &
      cnv_right(nd_layer+1-id_ct)                                       &
!       Right edges of convective cloud
    , strat_left(nd_layer+1-id_ct)                                      &
!       Left edges of stratiform cloud
    , x_cnv                                                             &
!       Right edge of current convective cloud
    , x_st                                                              &
!       Left edge of current stratiform cloud
    , x_done                                                            &
!       Fraction of the column treated
    , x_new_done                                                        &
!       Fraction of the column treated after adding new column
    , dx_col
!       Width of the current column
  REAL (RealK) ::                                                       &
      tol_cloud
!       Tolerance used in neglecting cloudy columns

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CLOUD_MAXCS_SPLIT'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerance used for clouds.
  tol_cloud=1.0e+04_RealK*EPSILON(tol_cloud)

  IF (n_cloud_type == 2) THEN
    ptr_st=0
    ptr_cnv=0
    DO k=1, n_cloud_type
      IF (k == ip_cloud_type_strat) ptr_st=k
      IF (k == ip_cloud_type_conv) ptr_cnv=k
    END DO
    IF ( (ptr_st == 0).OR.(ptr_cnv == 0) ) THEN
      cmessage = '*** Error: A type of cloud is missing.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  ELSE IF (n_cloud_type == 1) THEN
!   Only stratiform cloud is present.
    ptr_cnv=0
    ptr_st=1
  ELSE
    cmessage = '*** Error: There are too many types of cloud for '      &
      //'the type of overlap.'
    ierr=i_err_fatal
    GO TO 9999
  END IF


  n_cld_layer=n_layer+1-n_cloud_top
! We decompose a column at a time, as this is algorithmically
! easier, even if not compatible with vectorization.
  DO l=1, n_profile

!   Cloud is aligned with convective cloud against the left-hand
!   edge of the grid-box and stratiform cloud to the right. We
!   therefore need to find the right-hand edge of convective
!   cloud and the left-hand edge of stratiform cloud to partition.
    DO i=n_cloud_top, n_layer
      ii=i+1-n_cloud_top

      IF (n_cloud_type == 2) THEN
!       Calculate an explicit position for convective cloud if
!       included.
        cnv_right(ii)=w_cloud(l, i)*frac_cloud(l, i, ptr_cnv)
      ELSE
!       In the absence of convective cloud set its width to 0.
        cnv_right(ii)=0.0e+00_RealK
      END IF

      strat_left(ii)=1.0e+00_RealK-w_cloud(l, i)*frac_cloud(l, i        &
        , ptr_st)
!     Initialize the sorting key.
      key_st(ii)=ii
      key_cnv(ii)=ii
    END DO

!   Find the key ranking these edges in increasing order.
    CALL shell_sort(n_cld_layer, key_cnv, cnv_right)
    CALL shell_sort(n_cld_layer, key_st, strat_left)


!   Build up the list of notional columns and the list of those
!   where a solution is required.
    n_column_cld(l)=0
    n_column_slv(l)=0
!   Set the changes from a totally clear column to the first
!   actually used.
    DO i=n_cloud_top, n_layer
      ii=i+1-n_cloud_top
      IF (cnv_right(ii) >  tol_cloud) THEN
        IF (n_column_cld(l) <  nd_column) THEN
          n_column_cld(l)=n_column_cld(l)+1
        ELSE
          cmessage = '*** Error: ND_COLUMN is too small for the cloud ' &
            //'geometry selected.'
          ierr=i_err_fatal
          GO TO 9999
        END IF
        i_clm_lyr_chn(l, n_column_cld(l))=i
        i_clm_cld_typ(l, n_column_cld(l))=ptr_cnv
        area_column(l, n_column_cld(l))=0.0e+00_RealK
      END IF
      IF (strat_left(ii) <= tol_cloud) THEN
        IF (n_column_cld(l) <  nd_column) THEN
          n_column_cld(l)=n_column_cld(l)+1
        ELSE
          cmessage = '*** Error: ND_COLUMN is too small for the cloud ' &
            //'geometry selected.'
          ierr=i_err_fatal
          GO TO 9999
        END IF
        i_clm_lyr_chn(l, n_column_cld(l))=i
        i_clm_cld_typ(l, n_column_cld(l))=ptr_st
        area_column(l, n_column_cld(l))=0.0e+00_RealK
      END IF
    END DO

!   Now set up the mapping changing the contents of each layer
!   proceeding to the right.
    x_done=0.0e+00_RealK
!   Set the positions of the next convective and stratiform
!   changes, together with their corresponding indices.
    i_key_cnv=1
    x_cnv=cnv_right(key_cnv(1))
    DO WHILE ( (i_key_cnv <  n_cld_layer).AND.(x_cnv <  tol_cloud) )
      i_key_cnv=i_key_cnv+1
      x_cnv=cnv_right(key_cnv(i_key_cnv))
    END DO
    i_key_st=1
    x_st=strat_left(key_st(1))
    DO WHILE ( (i_key_st <  n_cld_layer).AND.(x_st <  tol_cloud) )
      i_key_st=i_key_st+1
      x_st=strat_left(key_st(i_key_st))
    END DO

!   Proceed throught the grid-box making the changes.
    DO WHILE (x_done <  1.0e+00_RealK-tol_cloud)

      IF (x_cnv <= x_st) THEN
!       The next change involves clearing convective cloud.
        IF (n_column_cld(l) <  nd_column) THEN
          n_column_cld(l)=n_column_cld(l)+1
        ELSE
          cmessage = '*** Error: ND_COLUMN is too small for the cloud ' &
            //'geometry selected.'
          ierr=i_err_fatal
          GO TO 9999
        END IF
        i_clm_lyr_chn(l, n_column_cld(l))=key_cnv(i_key_cnv)            &
          +n_cloud_top-1
        i_clm_cld_typ(l, n_column_cld(l))=0
        x_new_done=x_cnv
        i_key_cnv=i_key_cnv+1
        IF (i_key_cnv <= n_cld_layer) THEN
          x_cnv=cnv_right(key_cnv(i_key_cnv))
        ELSE
!         There are no further changes to convective cloud
!         right of this.
          x_cnv=1.0e+00_RealK
        END IF
      ELSE IF (x_cnv >  x_st) THEN
!       The next change involves introducing stratiform cloud.
        IF (n_column_cld(l) <  nd_column) THEN
          n_column_cld(l)=n_column_cld(l)+1
        ELSE
          cmessage = '*** Error: ND_COLUMN is too small for the cloud ' &
            //'geometry selected.'
          ierr=i_err_fatal
          GO TO 9999
        END IF
        i_clm_lyr_chn(l, n_column_cld(l))=key_st(i_key_st)              &
          +n_cloud_top-1
        i_clm_cld_typ(l, n_column_cld(l))=ptr_st
        x_new_done=x_st
        i_key_st=i_key_st+1
        IF (i_key_st <= n_cld_layer) THEN
          x_st=strat_left(key_st(i_key_st))
        ELSE
!         There are no further changes to stratiform cloud
!         right of this.
          x_st=1.0e+00_RealK
        END IF
      END IF

!     If both convective and stratiform right markers have
!     reached 1 we have a closing column.
      IF ( (x_st >  1.0e+00_RealK-tol_cloud).AND.                       &
           (x_cnv >  1.0e+00_RealK-tol_cloud) )                         &
        x_new_done=1.0e+00


!     If this new column is wide enough we solve within it.
      dx_col=x_new_done-x_done
      IF (dx_col >  tol_cloud) THEN
        n_column_slv(l)=n_column_slv(l)+1
        list_column_slv(l, n_column_slv(l))=n_column_cld(l)
        area_column(l, n_column_cld(l))=dx_col
        x_done=x_new_done
      ELSE
        area_column(l, n_column_cld(l))=0.0e+00_RealK
      END IF

    END DO

  END DO


  9999 CONTINUE
! Check error condition
  IF (ierr /= i_normal) THEN
    CALL ereport(RoutineName, ierr, cmessage)
  END IF

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cloud_maxcs_split
END MODULE cloud_maxcs_split_mod
