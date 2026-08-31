! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module defining the size distributions.
!
MODULE def_size_dist
!
! Description:
!   This module defines the available size distributions and
!   associated information on particle shapes.
!
!- End of header
!
!
! Modules used:
  USE realtype_rd
!
!
!
  TYPE STR_size_dist
!
!
    INTEGER :: i_shape_particle
!                Shape of the particles

    INTEGER :: i_measure
!                The way in which the dimension of the particle
!                is measured
!
    INTEGER :: i_distribution
!                Type of size distribution
!
!   Externally specified distribution:
    INTEGER :: n_size
!                Number of sizes in the distribution
!
    REAL (RealK), Pointer :: dimen(:)
!                     Dimensions of size distribution
    REAL (RealK), Pointer :: number(:)
!                     Populations in distribution
    REAL (RealK), Pointer :: aspect(:)
!                     Aspect ratio of particle
    REAL (RealK), Pointer :: d2_number(:)
!                     2nd derivative for spline
    REAL (RealK), Pointer :: d2_aspect(:)
!                     2nd derivative for spline
!
!
!   Analytic Distributions:
    INTEGER :: n_mode
!                Number of modes in the distribution
    REAL (RealK), Pointer :: weight_mode(:)
!                     Weights for each mode
!
!   Log-normal Distribution:
    REAL (RealK) :: n0_ln
!                     Scaling factor for log normal distribution
    REAL (RealK), Pointer :: ln_r0_ln(:)
!                     Reference radius for log normal distribution
    REAL (RealK), Pointer :: ln_sigma_ln(:)
!                     Standard deviation of log normal distribution
!
!   Modified Gamma Distribution:
    REAL (RealK), Pointer :: coeff_mg(:)
!                     Parameter for modified gamma distribution
    REAL (RealK), Pointer :: alpha_mg(:)
!                     Parameter for modified gamma distribution
    REAL (RealK), Pointer :: rm_mg(:)
!                     Parameter for modified gamma distribution
    REAL (RealK), Pointer :: beta_mg(:)
!                     Parameter for modified gamma distribution
!
!   Heymsfield and Platt's Distributions:
    REAL (RealK) :: d_hp_min
!                     Least maximum dimension
    REAL (RealK) :: a_hp_1
!                       1st intercept
    REAL (RealK) :: b_hp_1
!                       1st slope
    REAL (RealK) :: a_hp_2
!                       2nd intercept
    REAL (RealK) :: b_hp_2
!                       2nd slope
    REAL (RealK) :: d_hp_0
!                     Transitional dimension
!
!   Parameters for Mitchell and Ivanova's distributions
!   (These are bimodal modified-gamma distributions, but
!   the mean maximum dimension which characterizes the
!   distribution is stored for reference)
    REAL  (RealK) :: d_mean_max_l
!                     Mean maximum dimension of the large mode
    REAL  (RealK) :: lambda_l
!                     Slope parameter of the large mode
    REAL  (RealK) :: nu_l
!                     Nu for the large mode
    REAL  (RealK) :: n0_l
!                     N0 for the large mode
    REAL  (RealK) :: d_mean_max_s
!                     Mean maximum dimension of the small mode
    REAL  (RealK) :: lambda_s
!                     Slope parameter of the large mode
    REAL  (RealK) :: nu_s
!                     Nu for the small mode
    REAL  (RealK) :: n0_s
!                     N0 for the small mode
!
!
    REAL (RealK) :: dimen_0
!                   Typical dimension used in estimating integrals
!
  END TYPE STR_size_dist
!
!
!
END MODULE
