! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set pointers to types of clouds
!
! Method:
!   The types of condensate included are examined. Their phases
!   are set and depending on the representation of clouds adopted
!   it is determined to which type of cloud they contribute.
!
!- ---------------------------------------------------------------------
MODULE set_cloud_pointer_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_CLOUD_POINTER_MOD'
CONTAINS
SUBROUTINE set_cloud_pointer(ierr                                       &
     , n_condensed, type_condensed, i_cloud_representation              &
     , l_drop, l_ice                                                    &
     , i_phase_cmp, i_cloud_type, l_cloud_cmp                           &
     , nd_cloud_component                                               &
     )


  USE rad_pcf, ONLY: i_err_fatal, &
   ip_cloud_homogen, ip_cloud_combine_homogen, &
   ip_cloud_ice_water, ip_cloud_combine_ice_water, &
   ip_cloud_conv_strat, ip_cloud_split_homogen, &
   ip_cloud_csiw, ip_cloud_split_ice_water, &
   ip_clcmp_st_water, ip_clcmp_st_ice, &
   ip_clcmp_cnv_water, ip_clcmp_cnv_ice, &
   ip_cloud_type_homogen, ip_cloud_type_water, ip_cloud_type_ice, &
   ip_cloud_type_strat, ip_cloud_type_conv, &
   ip_cloud_type_sw, ip_cloud_type_si, &
   ip_cloud_type_cw, ip_cloud_type_ci, &
   ip_phase_water, ip_phase_ice
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Dimensions of arrays
  INTEGER, INTENT(IN) ::                                                &
    nd_cloud_component
!     Maximum number of condensed components allowed in clouds

! Dummy variables.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_condensed                                                       &
!       Number of condensed components
    , type_condensed(nd_cloud_component)                                &
!       Types of components
    , i_cloud_representation
!       Representation of clouds used
  LOGICAL, INTENT(IN) ::                                                &
      l_drop                                                            &
!       Flag for inclusion of droplets
    , l_ice
!       Flag for inclusion of ice crystals

  INTEGER, INTENT(OUT) ::                                               &
      i_phase_cmp(nd_cloud_component)                                   &
!       Phases of components
    , i_cloud_type(nd_cloud_component)
!       Types of cloud to which each component contributes
  LOGICAL, INTENT(OUT) ::                                               &
      l_cloud_cmp(nd_cloud_component)
!       Logical switches to `include' components


! Local variables
  INTEGER                                                               &
      k
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SET_CLOUD_POINTER'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO k=1, n_condensed

    i_cloud_type(k)=0
!   Set pointers for valid condensed components.
    SELECT CASE (i_cloud_representation)
    CASE (ip_cloud_homogen, ip_cloud_combine_homogen)

      IF (type_condensed(k) == ip_clcmp_st_water) THEN
        i_cloud_type(k)=ip_cloud_type_homogen
      ELSE IF (type_condensed(k) == ip_clcmp_st_ice) THEN
        i_cloud_type(k)=ip_cloud_type_homogen
      END IF

    CASE (ip_cloud_ice_water, ip_cloud_combine_ice_water)

      IF (type_condensed(k) == ip_clcmp_st_water) THEN
        i_cloud_type(k)=ip_cloud_type_water
      ELSE IF (type_condensed(k) == ip_clcmp_st_ice) THEN
        i_cloud_type(k)=ip_cloud_type_ice
      END IF

    CASE (ip_cloud_conv_strat, ip_cloud_split_homogen)

      IF (type_condensed(k) == ip_clcmp_st_water) THEN
        i_cloud_type(k)=ip_cloud_type_strat
      ELSE IF (type_condensed(k) == ip_clcmp_st_ice) THEN
        i_cloud_type(k)=ip_cloud_type_strat
      ELSE IF (type_condensed(k) == ip_clcmp_cnv_water) THEN
        i_cloud_type(k)=ip_cloud_type_conv
      ELSE IF (type_condensed(k) == ip_clcmp_cnv_ice) THEN
        i_cloud_type(k)=ip_cloud_type_conv
      END IF

    CASE (ip_cloud_csiw, ip_cloud_split_ice_water)

      IF (type_condensed(k) == ip_clcmp_st_water) THEN
        i_cloud_type(k)=ip_cloud_type_sw
      ELSE IF (type_condensed(k) == ip_clcmp_st_ice) THEN
        i_cloud_type(k)=ip_cloud_type_si
      ELSE IF (type_condensed(k) == ip_clcmp_cnv_water) THEN
        i_cloud_type(k)=ip_cloud_type_cw
      ELSE IF (type_condensed(k) == ip_clcmp_cnv_ice) THEN
        i_cloud_type(k)=ip_cloud_type_ci
      END IF

    CASE DEFAULT

      cmessage = 'Invalid cloud representation'
      ierr = i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)

    END SELECT

!   Check for 0 flagging illegal types.
    IF (i_cloud_type(k) == 0) THEN
      cmessage = 'A component is not compatible with the ' &
        //'representation of clouds selected.'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END IF

    IF (type_condensed(k) == ip_clcmp_st_water) THEN

      i_phase_cmp(k)=ip_phase_water
      l_cloud_cmp(k)=l_drop

    ELSE IF (type_condensed(k) == ip_clcmp_st_ice) THEN

      i_phase_cmp(k)=ip_phase_ice
      l_cloud_cmp(k)=l_ice

    ELSE IF (type_condensed(k) == ip_clcmp_cnv_water) THEN

      i_phase_cmp(k)=ip_phase_water
      l_cloud_cmp(k)=l_drop

    ELSE IF (type_condensed(k) == ip_clcmp_cnv_ice) THEN

      i_phase_cmp(k)=ip_phase_ice
      l_cloud_cmp(k)=l_ice

    END IF

  END DO


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_cloud_pointer
END MODULE set_cloud_pointer_mod
