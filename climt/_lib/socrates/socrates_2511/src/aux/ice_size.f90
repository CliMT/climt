! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Program to set ice sizes in a particular atmosphere
!
PROGRAM ice_size
!
! Description:
!    This code deteermines the sizes of ice crystals in a given
! atmosphere from the macrophysical fields.
!
! Method:
!   A character string is read to determine the parametrization to be
! used.
!
! NOTE: This program antedates the conversion of the CDL-input system 
! to F90. File suffixes are hard-wired.
!
! Modules used:
  USE realtype_rd
  USE dimensions_field_ucf
  USE dimensions_cdl_ucf
  USE def_atm
  USE def_cld
  USE def_std_io_icf
  USE error_pcf
!
!
!
  IMPLICIT NONE
!
!
!
  INTEGER :: ierr
!   Error flag
  CHARACTER (LEN=80) :: base_name
!   Base name of CDL input files
  CHARACTER (LEN=80) :: cdl_file_name
!   Base name of CDL input files
  CHARACTER (LEN=12) :: suffix
!   Suffix appended to CDL files indicating the type
  CHARACTER (LEN=24) :: name_vert_coord
!   Name of the vertical coordinate
  LOGICAL :: l_vert_coord = .FALSE.
!   Flag indicating that the vertical coordinate has been set
  LOGICAL :: l_vert_coord_level = .FALSE.
!   Flag indicating that the vertical coordinate for fields at the
!   edges of layers has been set
!
  TYPE  (StrAtm) :: Atm
!   Atmospheric profile
  TYPE  (StrCld) :: Cld
!   Cloud profile

  CHARACTER (LEN=24) :: param
!   Name of the vertical coordinate
!
  INTEGER :: n_lat = 0
!   Number of latitudes
  INTEGER :: n_lon = 0
!   Number of longitudes
  INTEGER :: n_level
!   Dummy number of levels
  INTEGER :: i
!   Loop variable
  INTEGER :: l
!   Loop variable
  INTEGER :: seed
!   Seed for random number generator
!
  REAL  (RealK) :: mean_size
!   Mean size of crystal
  REAL  (RealK) :: stdev_size
!   Standard deviation of crystal size
  REAL  (RealK) :: rand_gauss
!   Random variable from a Gaussian distribution with zero mean
!   and unit variance
!
!
!
! Provisional code: once the CDL routines have been converted to F90,
! this can be removed, as allocation will be done in the CDL routines.
  ALLOCATE(Atm%lat(npd_profile))
  ALLOCATE(Atm%lon(npd_profile))
  ALLOCATE(Atm%p(npd_profile, npd_layer))
  ALLOCATE(Atm%t(npd_profile, npd_layer))
  ALLOCATE(Atm%p_level(npd_profile, 0:npd_layer))
  ALLOCATE(Atm%t_level(npd_profile, 0:npd_layer))
  ALLOCATE(Cld%condensed_mix_ratio(npd_profile, npd_layer, 1))
  ALLOCATE(Cld%condensed_dim_char(npd_profile, npd_layer, 1))
!
! Read in those properties of the atmosphere which may affect the size.
  WRITE(iu_stdout, "(a)") "Enter the basename for the CDL files."
  READ(iu_stdin, "(a)") base_name
!
! Read the temperatures at the centres of layers
  cdl_file_name = TRIM(base_name)//'.t'
  CALL assign_input_vert_cdl(ierr, &
    cdl_file_name, 'temperatures', &
    l_vert_coord, name_vert_coord, &
    .TRUE., Atm%n_layer, .NOT.l_vert_coord, &
    n_lat, Atm%lat, n_lon, Atm%lon, 1, &
    Atm%n_profile, Atm%n_layer, &
    Atm%p, Atm%t, &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP
!
  cdl_file_name = TRIM(base_name)//".iwm"
  CALL assign_input_vert_cdl(ierr, &
    cdl_file_name, 'ice water mixing ratios', &
    l_vert_coord, name_vert_coord, &
    .TRUE., Atm%n_layer, .FALSE., &
    n_lat, Atm%lat, n_lon, Atm%lon, 1, &
    Atm%n_profile, Atm%n_layer, &
    Atm%p, Cld%condensed_mix_ratio, &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP
!
!
!
!
! Get the type of parametrization
  WRITE(iu_stdout, '(a)') 'Enter parametrization scheme'
  READ(iu_stdin, '(a)') param
!
  SELECT CASE(TRIM(param))
    CASE('Kristj_dl')
      Cld%condensed_dim_char(:, :, 1) = &
              Exp(5.522E-02 * (Atm%t - 2.7965E+02)) / 9.702E+02
      WHERE(Cld%condensed_dim_char(:, :, 1) > 7.198755E-04)
        Cld%condensed_dim_char(:, :, 1) = 7.198755E-04
      ENDWHERE
    CASE('Havem_dm')
      Cld%condensed_dim_char(:, :, 1) = 1.39265E-04 + &
        1.4431E-06 * (Atm%t - 2.7315E+02)
      WHERE(Cld%condensed_dim_char(:, :, 1) < 7.0E-06)
        Cld%condensed_dim_char(:, :, 1) = 7.0E-06
      ENDWHERE
    CASE('Havem_de')
      Cld%condensed_dim_char(:, :, 1) = 2.436614E-04 + &
        3.7199E-06 * (Atm%t - 2.7315E+02)
      WHERE(Cld%condensed_dim_char(:, :, 1) < 7.0E-06)
        Cld%condensed_dim_char(:, :, 1) = 7.0E-06
      ENDWHERE
      WHERE(Cld%condensed_dim_char(:, :, 1) > 1.70E-04)
        Cld%condensed_dim_char(:, :, 1) = 1.70E-04
      ENDWHERE
    CASE('Havem_dge')
      Cld%condensed_dim_char(:, :, 1) = 1.875706E-04 + &
        2.8636E-06 * (Atm%t - 2.7315E+02)
      WHERE(Cld%condensed_dim_char(:, :, 1) < 5.3886E-06)
        Cld%condensed_dim_char(:, :, 1) = 5.3886E-06
      ENDWHERE
    CASE('Mitch_dge')
!     Stephan's Fit to the Kristjansson/Mitchell scheme
      WHERE(Atm%t <  216.208)
        Cld%condensed_dim_char(:, :, 1) = 7.5094588E-04 * &
          EXP(0.05*(Atm%t(:, :) - 279.5)) + 5.0830326E-07
      ELSEWHERE
        Cld%condensed_dim_char(:, :, 1) = 1.3505403E-04 * &
          EXP(0.05*(Atm%t(:, :) - 279.5)) + 2.6517429E-05
      ENDWHERE
!     Limit the particle sizes to sensible ranges.
      WHERE(Cld%condensed_dim_char(:, :, 1) < 8.0E-06)
        Cld%condensed_dim_char(:, :, 1) = 8.0E-06
      ENDWHERE
      WHERE(Cld%condensed_dim_char(:, :, 1) > 1.24E-04)
        Cld%condensed_dim_char(:, :, 1) = 1.24E-04
      ENDWHERE
    CASE('Mitch_de')
!     Stephan's Fit to the Kristjansson/Mitchell scheme
!     We get dge first and convert.
      WHERE(Atm%t <  216.208)
        Cld%condensed_dim_char(:, :, 1) = 7.5094588E-04 * &
          EXP(0.05*(Atm%t(:, :) - 279.5)) + 5.0830326E-07
      ELSEWHERE
        Cld%condensed_dim_char(:, :, 1) = 1.3505403E-04 * &
          EXP(0.05*(Atm%t(:, :) - 279.5)) + 2.6517429E-05
      ENDWHERE
!     Limit the particle sizes to sensible ranges.
      WHERE(Cld%condensed_dim_char(:, :, 1) < 8.0E-06)
        Cld%condensed_dim_char(:, :, 1) = 8.0E-06
      ENDWHERE
      WHERE(Cld%condensed_dim_char(:, :, 1) > 1.24E-04)
        Cld%condensed_dim_char(:, :, 1) = 1.24E-04
      ENDWHERE
      Cld%condensed_dim_char(:, :, 1) = Cld%condensed_dim_char(:, :, 1) * &
        (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))
    CASE('Random_G')
      WRITE(iu_stdout, '(a)') &
        'Enter the mean size and standard deviation'
      READ(iu_stdin, *) mean_size, stdev_size
      DO i=1, Atm%n_layer
        DO l=1, Atm%n_profile
          Cld%condensed_dim_char(l, i, 1) = mean_size + stdev_size * &
            MAX( -2.0_RealK, MIN( 2.0_RealK, rand_gauss(seed)))
        ENDDO
      ENDDO
  END SELECT
!
  cdl_file_name = TRIM(base_name)//'.ire'
  CALL output_vert_cdl(ierr, &
    cdl_file_name, &
    n_lat, Atm%lat, n_lon, Atm%lon, &
    Atm%n_profile, Atm%n_layer, &
    TRIM(name_vert_coord), LEN(TRIM(name_vert_coord)), &
    Atm%p, 'ire', 'float', 'kg.kg-1', 'Characteristic Dimension', &
    Cld%condensed_dim_char(:, :, 1), &
    npd_profile, npd_latitude, npd_longitude, 1, npd_layer, &
    npd_cdl_dimen, npd_cdl_dimen_size, npd_cdl_data, npd_cdl_var)
  IF (ierr /= i_normal) STOP
!
!
!
END PROGRAM ice_size
