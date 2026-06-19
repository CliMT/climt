! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to rebin ESFT terms using the random overlap method with
! resorting and rebinning.
!
!- ---------------------------------------------------------------------
MODULE rebin_esft_terms_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE rebin_esft_terms(i_band_esft_mix, n_esft_red                 &
    , i_profile, i_layer                                                &
    , w_esft_target, glim_target                                        &
    , k_esft_layer_mix, w_esft_mix                                      &
    , k_esft_layer_mix_red, glim_mix_red                                &
    , nd_profile, nd_layer, nd_esft_max)
    
  USE realtype_rd, ONLY: RealK
  
  IMPLICIT NONE
  
! Input variables
  INTEGER ::                                                            &
      i_band_esft_mix                                                   &
!       Number of unbinned ESFT terms
    , n_esft_red                                                        &
!       Number of binned (reduced) ESFT terms
    , i_profile                                                         &
!       Index of current profile
    , i_layer                                                           &
!       Index of current layer
    , nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_esft_max
!       Maximum number of ESFT terms
  REAL (RealK), INTENT(IN) ::                                           &
      k_esft_layer_mix(nd_profile,nd_layer,nd_esft_max*nd_esft_max)     &
!       ESFT monochromatic exponents for the mixture of two gases
    , w_esft_mix(nd_esft_max*nd_esft_max)                               &
!       ESFT weights for the mixture of two gases
    , w_esft_target(nd_esft_max)                                        &
!       Target weights when performing rebinning
    , glim_target(nd_esft_max+1)
!       Target g-coordinate limits for the ESFT terms when performing
!       rebinning

! Output variables
  REAL (RealK), INTENT(OUT) ::                                          &
      k_esft_layer_mix_red(nd_profile,nd_layer,nd_esft_max)             &
!       Reduced (resorted and rebinned) ESFT monochromatic exponents
!       for the mixture of two gases
    , glim_mix_red(nd_esft_max+1)
!       g-coordinate limits for the reduced ESFT terms

! Local variables
  REAL (RealK) ::                                                       &
      w_esft_mix_current                                                &
!       Current weight for the gas mixture
    , w_esft_part                                                       &
!       Weight used when splitting mixed weight over several bins
    , eps = 1.0E-08_RealK
!       Tolerance used to take into account loss of precision due to
!       number of decimal points on weights in spectral files
  INTEGER ::                                                            &
      iex_mix                                                           &
!       Loop variable for unbinned ESFT terms
    , iex_mix_red
!       Loop variable for binned ESFT terms
  
  glim_mix_red=0.0_RealK
  iex_mix_red=1
  iex_mix=1
  DO WHILE (iex_mix <= i_band_esft_mix)
    w_esft_mix_current = w_esft_mix(iex_mix)
    
    DO WHILE (glim_mix_red(iex_mix_red+1) + w_esft_mix_current >        &
      glim_target(iex_mix_red+1)*(1.0_RealK + eps))
!     End of this reduced ESFT term has been reached
      
!     Weights for the bins the ESFT term will be split over
      w_esft_part = glim_target(iex_mix_red+1) -                        &
        glim_mix_red(iex_mix_red+1)
      
!     Include ESFT term with corresponding weights in current bin
      glim_mix_red(iex_mix_red+1) = glim_target(iex_mix_red+1)
      k_esft_layer_mix_red(i_profile,i_layer,iex_mix_red) =             &
        k_esft_layer_mix_red(i_profile,i_layer,iex_mix_red) +           &
        w_esft_part*k_esft_layer_mix(i_profile,i_layer,iex_mix)
      
!     Adjust current weight
      w_esft_mix_current = w_esft_mix_current - w_esft_part
      
!     Move on to next reduced term
      iex_mix_red = iex_mix_red + 1
      glim_mix_red(iex_mix_red+1) = glim_mix_red(iex_mix_red)
    END DO
    
!   Add this ESFT term to the reduced ESFT term
    glim_mix_red(iex_mix_red+1) = glim_mix_red(iex_mix_red+1) +         &
      w_esft_mix_current
    k_esft_layer_mix_red(i_profile,i_layer,iex_mix_red) =             &
      k_esft_layer_mix_red(i_profile,i_layer,iex_mix_red) +           &
      w_esft_mix_current*k_esft_layer_mix(i_profile,i_layer,iex_mix)
    iex_mix=iex_mix+1
  END DO
  
! Divide by norm in each bin to get the reduced ESFT terms
  k_esft_layer_mix_red(i_profile,i_layer,1:n_esft_red) =                &
    k_esft_layer_mix_red(i_profile,i_layer,                             &
    1:n_esft_red)/w_esft_target(1:n_esft_red)
  
END SUBROUTINE rebin_esft_terms
END MODULE rebin_esft_terms_mod
