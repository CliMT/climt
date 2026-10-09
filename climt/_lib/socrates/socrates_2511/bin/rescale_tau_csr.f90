! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to rescale optical depth for direct flux
!
! Method:
!   The phase function scaling that includes circumsolar contribution
!   within FOV of an instrument is used to replace 
!   the Delta-Eddington scaling      
!
!----------------------------------------------------------------------
MODULE rescale_tau_csr_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE rescale_tau_csr(n_profile                                    &
    , i_layer_first, i_layer_last                                       &
    , forward_scatter                                                   &
    , tau, tau_dir, omega                                               &
    , nd_profile, nd_layer, id_1                                        &
    )

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE


! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , id_1
!       Topmost declared layer for optical properties
!       Size allocated for layers
! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to rescale
    , i_layer_last
!       First layer to rescale
  REAL (RealK), INTENT(IN) ::                                           &
   &  forward_scatter(nd_profile, id_1: nd_layer)                       &                
!       Circumsolar forward scatter fraction                             
    , omega(nd_profile, id_1: nd_layer)
!       Single scattering albedo
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, id_1: nd_layer)                                          
!       Optical depth
 
  REAL (RealK), INTENT(OUT) ::                                          &
      tau_dir(nd_profile, id_1: nd_layer)                                      
!       scaled Optical depth by circumsolar fraction
!
! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l                                                                 
!       Loop variable
!
  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile
       tau_dir(l, i) = tau(L, I)*(1.0e+00_RealK                         &
         - omega(l,i)*forward_scatter(l,i))
    END DO
  END DO
! 
END SUBROUTINE rescale_tau_csr     
END MODULE rescale_tau_csr_mod
