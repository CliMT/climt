! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to calculate number denisty of particles.
!
FUNCTION number_particle_90 &
!
(dimen_point, SizeDist, ierr)
!
! Description:
!   This routine returns the number of particles per unit size range
!   as evaluated from the distribution.

! Modules used.
  USE realtype_rd
  USE distribution_pcf
  USE def_size_dist
  USE error_pcf
  USE rad_ccf, ONLY: pi
  USE max_size_acf
  USE spline_evaluate_mod, ONLY: spline_evaluate
!
!
  IMPLICIT NONE
!
!
! Dummy arguments.
  REAL  (RealK), Intent(IN) :: dimen_point
!           Dimension for evaluation
  TYPE (STR_size_dist), Intent(IN) :: SizeDist
!           Size distribution
!
  REAL  (RealK) :: number_particle_90
!           Differential population
  INTEGER, Intent(InOut) :: ierr
!           Error flag
!
! Local variables.
  INTEGER :: n_knot
!           Number of knot points
!! ----------------------------------------------------
!! Variables not used at present: see below.
!  INTEGER :: i_low
!!           Bisection variable
!  INTEGER :: i_high
!!           Bisection variable
!  INTEGER :: i_split
!!           Bisection variable
!! ----------------------------------------------------
  INTEGER :: i
!           Loop variable
!
  REAL  (RealK) :: log_dimen
!           Log of dimension
  REAL  (RealK) :: max_dimension
!           Maximum dimension of particle
  REAL  (RealK) :: length
!           Length of crystal
  REAL  (RealK) :: discriminant
!           Discriminant
  REAL  (RealK) :: conversion
!           Conversion factor
!
  REAL  (RealK), Allocatable, Dimension(:) :: zeros
!   Dummy array of zeros to make splining routine behave like a 
!   linear fit.
!
!
! A cutoff is placed at a prescribed limit.
  IF (dimen_point > max_dimen_scatterer) THEN
    number_particle_90 = 0.0_RealK
    RETURN
  ENDIF
!
! Evaluate the differential population as directed by the distribution.
  IF (SizeDist%i_distribution == ip_external) THEN
    IF (SizeDist%n_size > 1) THEN
!
!     The full splined fit is not used at the moment because observed
!     distributions are very irregular, so a dummy array of zeros is
!     passed to the splining routine.
      n_knot = SizeDist%n_size+1
!     N.B. Check that d2_number has been allocated if you ever make
!     use of this.
      ALLOCATE(zeros(SIZE(SizeDist%number)))
      zeros = 0.0_RealK
      ierr=i_normal
      CALL spline_evaluate(ierr, n_knot, &
             SizeDist%dimen, SizeDist%number, zeros, &
             dimen_point, number_particle_90)
      DEALLOCATE(zeros)
      IF (ierr /= i_normal) THEN
        IF (ierr == i_err_range) THEN
          number_particle_90 = 0.0_RealK
          ierr = i_normal
        ELSE
          ierr = i_err_fatal
          RETURN
        ENDIF
      ELSE
!       Deal with the possibility of small negative numbers.
        number_particle_90 = MAX(0.0_RealK, number_particle_90)
      ENDIF
!
!!     ----------------------------------------------------------------
!!     The following checks are unnecessary unless full splinig is used.
!!
!!     The spline fit could give negative results, and data may
!!     have stretches over which the population is zero: the fit
!!     must therefore be adjusted to give zero in such cases.
!      IF (number_particle_90 < 0.0_RealK) number_particle_90 = 0.0_RealK
!!     Chop to find radius interval for check.
!      i_low = 0
!      i_high = n_knot
!      DO
!        i_split = (i_low + i_high) / 2
!        IF (dimen_point < SizeDist%dimen(i_split)) THEN
!          i_high = i_split
!        ELSE
!          i_low = i_split
!        ENDIF
!        IF ( (i_high - i_low) <= 1) EXIT
!      ENDDO
!      IF ( (SizeDist%number(i_low)  < 0.0_RealK) .AND. &
! &         (SizeDist%number(i_high) < 0.0_RealK) ) THEN
!        number_particle_90 = 0.0_RealK
!      ENDIF
!     -------------------------------------------------------------------
!
    ELSE
      number_particle_90 = SizeDist%number(SizeDist%n_size)
    ENDIF
! 
  ELSE IF (SizeDist%i_distribution.EQ.IP_LOG_NORMAL) THEN
!
    log_dimen = LOG(dimen_point)
    number_particle_90 = 0.0_RealK
    DO i=1, SizeDist%n_mode
      number_particle_90 = number_particle_90 + SizeDist%weight_mode(i) * &
        (1.0_RealK/dimen_point) * &
        (SizeDist%n0_ln / (SizeDist%ln_sigma_ln(i) * &
        SQRT(2.0_RealK * PI))) * &
        EXP(-(log_dimen-SizeDist%ln_r0_ln(i))**2 / &
        (2.0_RealK * SizeDist%ln_sigma_ln(i)**2))
    ENDDO
!
  ELSE IF (SizeDist%i_distribution == ip_modified_gamma) THEN
!
    number_particle_90 = 0.0_RealK
    DO i=1, SizeDist%n_mode
      log_dimen=LOG(dimen_point/SizeDist%rm_mg(i))
      number_particle_90 = number_particle_90 + SizeDist%weight_mode(i) * &
        SizeDist%coeff_mg(i) * &
        EXP((SizeDist%alpha_mg(i)-1.0_RealK) * log_dimen - &
        EXP(SizeDist%beta_mg(i) * log_dimen))
    ENDDO
!
  ELSE IF (SizeDist%i_distribution == ip_heymsfield_platt) THEN
!
!   Find the maximum dimension of particle to apply in the
!   distribution.
!   The formulae for the aspect ratio are taken from Auer and Veal
!   J. Atm. Sci. [27] p. 919 (1970).
    discriminant = 1.004004_RealK - 9.36E+03_RealK * &
      (2.0_RealK * dimen_point + 8.479E-06_RealK)
    IF (discriminant > 0.0_RealK) THEN
      length = 2.1368E-04_RealK * (1.002_RealK - SQRT(discriminant))
      IF (length > 2.0E-04_RealK) THEN
        length = 4.7185E+06_RealK * &
          EXP(2.415_RealK * LOG(dimen_point))
        conversion = 1.13952E+07_RealK * &
          EXP(1.415_RealK * LOG(dimen_point))
      ELSE
        conversion = 2.0_RealK / &
          (1.002_RealK - 4.68E+03_RealK * length)
      ENDIF
    ELSE
      length = 4.7185E+06_RealK * &
        EXP(2.415_RealK * LOG(dimen_point))
      conversion = 1.13952E+07_RealK * &
 &      EXP(1.415_RealK * LOG(dimen_point))
    ENDIF
    max_dimension = MAX(2.0_RealK * dimen_point, length)
    IF (max_dimension < SizeDist%d_hp_min) THEN
      number_particle_90 = 0.0_RealK
    ELSE IF (max_dimension <= SizeDist%d_hp_0) THEN
      number_particle_90 = SizeDist%a_hp_1 * &
        EXP(SizeDist%b_hp_1 * LOG(max_dimension))
    ELSE IF (dimen_point < max_dimen_scatterer) THEN
      number_particle_90 = SizeDist%a_hp_2 * &
        EXP(SizeDist%b_hp_2 * LOG(max_dimension))
    ELSE
      number_particle_90 = 0.0_RealK
    ENDIF
!   Since this distribution is given in terms of the maximum 
!   dimension of the particle, we now have the number density in 
!   a unit range of maximum dimensions, and this must be scaled 
!   for a unit range of radii.
    number_particle_90 = number_particle_90 * conversion
!
  
  ELSE IF ( (SizeDist%i_distribution == ip_mitchell_97)      .OR. &
            (SizeDist%i_distribution == ip_mitchell_trop_00) .OR. &
            (SizeDist%i_distribution == ip_ivanova_mlat_00) ) THEN
!
    number_particle_90 = SizeDist%n0_s * dimen_point**SizeDist%nu_s * &
      EXP(-SizeDist%lambda_s * dimen_point) + &
      SizeDist%n0_l * dimen_point**SizeDist%nu_l * &
      EXP(-SizeDist%lambda_l * dimen_point)
!
  ENDIF
!
!
!
  RETURN
END
