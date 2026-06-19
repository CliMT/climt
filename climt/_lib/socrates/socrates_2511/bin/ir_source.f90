! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calcaulate IR source function for differential flux.
!
! Method:
!   The linear contribution to the source function is proportional
!   to the absorption divided by the optical depth. A tolerance is
!   added to the optical depth to allow for the depth's being 0.
!   A correction may also be made for a quadratic variation in the
!   temperature across the layer.
!
!- ---------------------------------------------------------------------
MODULE ir_source_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'IR_SOURCE_MOD'
CONTAINS
SUBROUTINE ir_source(n_profile, i_layer_first, i_layer_last             &
     , source_coeff, del_planck, l_ir_source_quad, diff_planck_2        &
     , s_down, s_up                                                     &
     , nd_profile, nd_layer, nd_source_coeff                            &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_scf_ir_1d, ip_scf_ir_2d
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_source_coeff
!       Size allocated for source coefficients


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to consider
    , i_layer_last
!       Last layer to consider

  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Use a quadratic representation

  REAL (RealK), INTENT(IN) ::                                           &
      source_coeff(nd_profile, nd_layer, nd_source_coeff)               &
!       Coefficients for source terms
    , del_planck(nd_profile, nd_layer)                                  &
!       Difference in Planckian function across the layer
    , diff_planck_2(nd_profile, nd_layer)
!         2x2nd difference of Planckian

  REAL (RealK), INTENT(OUT) ::                                          &
      s_down(nd_profile, nd_layer)                                      &
!       Upward source function
    , s_up(nd_profile, nd_layer)
!       Upward source function


! Local variables.

  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='IR_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Multiply the source coefficients by the Planckian differences
! to the order required.
  IF (l_ir_source_quad) THEN

    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        s_up(l, i)=source_coeff(l, i, ip_scf_ir_1d)                     &
          *del_planck(l, i)                                             &
          +source_coeff(l, i, ip_scf_ir_2d)                             &
          *diff_planck_2(l, i)
        s_down(l, i)=-source_coeff(l, i, ip_scf_ir_1d)                  &
          *del_planck(l, i)                                             &
          +source_coeff(l, i, ip_scf_ir_2d)                             &
          *diff_planck_2(l, i)
      END DO

    END DO

  ELSE

    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        s_up(l, i)=source_coeff(l, i, ip_scf_ir_1d)                     &
          *del_planck(l, i)
        s_down(l, i)=-s_up(l, i)
      END DO
    END DO

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ir_source
END MODULE ir_source_mod
