! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate fractional contribution to the direct solar 
!  radiation due to scattering in the circumsolar region. 
!  Delta-scaling approxmation is tau1=(1-omega*g^2) tau. Here we 
!  replace g^2 with a fraction of circumsolar area PQ such that
!  tau1=(1-omega * PQ) tau. PQ is a phase function integrated 
!  from 0 to 2.5 of a solid angle of FOV. H-G phase function is used 
!  in the integration so the momont is just g^k. 
!  Note: forward_scatter_csr is determined using the weighted g 
!  and omega. This is consistent with the treatment of Delta-scaling.
!
!----------------------------------------------------------------------
MODULE circumsolar_fraction_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE circumsolar_fraction(n_profile                               &
    , indx, half_angle                                                  &
    , asymmetry_factor                                                  &
    , forward_scatter_csr                                               &
    , nd_profile                                                        &
    )

  USE realtype_rd, ONLY: RealK
  USE legendre_mod, ONLY: nstream, cos_half_angle, p_legendre

  IMPLICIT NONE

! Sizes of arrays
  INTEGER, INTENT(IN) ::                                                &
       nd_profile                                                       &
!       Size allocated for profiles
    , indx(nd_profile)
!       Indices satifying the test
! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         
!       Number of profiles
  REAL (RealK), INTENT(IN) ::                                           &
      asymmetry_factor(nd_profile)                                      &
!       Asymmetry factor (first moment of phase function)
    , half_angle
!       Half angle of pyrheliometer instrument FOV  
!
  REAL  (RealK), Intent(OUT) ::                                         &
      forward_scatter_csr(nd_profile)
!       Forward scattering in circumsolar region
!
! Local variables.
  INTEGER                                                               &
      l                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , ind_ang
!      index of half angle and cos_half_angle(angle_ind)=cos(half_angle)
     
!       
  REAL (RealK) ::                                                       &
      p_legendre_2d(32,1:20)                                            &
!      Legendre polynomial 32 series for 20 scattering angles
    , one_minus_cos
!      First term of integration
!
!
! Calculate index of Legedre polynomial based corresponding to
! the input half viewing angle of pyrheliometer 
  ind_ang=INT(half_angle/0.25)
  one_minus_cos = 1.0 - cos_half_angle(ind_ang)
  p_legendre_2d = reshape(p_legendre, (/32,20/))
  forward_scatter_csr=0.0
  DO k=2, nstream-1
     DO l=1, n_profile
        forward_scatter_csr(indx(l))=forward_scatter_csr(indx(l))        &
          +asymmetry_factor(indx(l))**(k-1)                              &
          *(p_legendre_2d(k-1, ind_ang) - p_legendre_2d(k+1, ind_ang))
     END DO
  END DO
  DO l=1, n_profile
    forward_scatter_csr(indx(l))=0.5*(one_minus_cos                      &
      +forward_scatter_csr(indx(l)))
  END DO
! 
END SUBROUTINE circumsolar_fraction
END MODULE circumsolar_fraction_mod
