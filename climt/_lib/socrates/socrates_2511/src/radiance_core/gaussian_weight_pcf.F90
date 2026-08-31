! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module to define arrays for Gaussian integration.
!
! Description:
!   This module defines arrays for Gaussian integration at low
!   orders (for non-scattering IR flux calculations). For
!   higher orders (used in preprocessing) the weights and points
!   are calculated directly and iteratively.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
MODULE gaussian_weight_pcf


  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE


  INTEGER, Parameter :: npd_gauss_ord  = 10
!   Maximum order of Gaussian quadrature

  INTEGER :: i

  REAL (RealK), PARAMETER ::                                            &
!   Points of Gaussian integration
    gauss_point(npd_gauss_ord, npd_gauss_ord) = RESHAPE( (/             &
!     First order:
      0.0_RealK, (0.0_RealK, i=1, NPD_gauss_ord-1) ,                    &
!     Second order:
      -5.77350269189626E-01_RealK, 5.77350269189626E-01_RealK,          &
        (0.0_RealK, i=1, NPD_gauss_ord-2),                              &
!     Third order:
      -7.74596669241484E-01_RealK, 0.0_RealK,                           &
        7.74596669241484E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-3),  &
!     Fourth order:
      -8.61136311594053E-01_RealK, -3.39981043584856E-01_RealK,         &
        3.39981043584856E-01_RealK, 8.61136311594053E-01_RealK,         &
        (0.0_RealK, i=1, NPD_gauss_ord-4),                              &
!     Fifth order:
      -9.06179845938664E-01_RealK, -5.38469310105683E-01_RealK,         &
        0.0_RealK, 5.38469310105683E-01_RealK,                          &
        9.06179845938664E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-5),  &
!     Sixth order:
      -9.32469514203152E-01_RealK, -6.61209386466265E-01_RealK,         &
        -2.38619186083197E-01_RealK, 2.38619186083197E-01_RealK,        &
        6.61209386466265E-01_RealK, 9.32469514203152E-01_RealK,         &
        (0.0_RealK, i=1, NPD_gauss_ord-6),                              &
!     Seventh order:
      -9.49107912342759E-01_RealK, -7.41531185599394E-01_RealK,         &
        -4.05845151377397E-01_RealK, 0.0_RealK,                         &
        4.05845151377397E-01_RealK, 7.41531185599394E-01_RealK,         &
        9.49107912342759E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-7),  &
!     Eighth order:
      -9.60289856497536E-01_RealK, -7.96666477413627E-01_RealK,         &
        -5.25532409916329E-01_RealK, -1.83434642495650E-01_RealK,       &
        1.83434642495650E-01_RealK, 5.25532409916329E-01_RealK,         &
        7.96666477413627E-01_RealK, 9.60289856497536E-01_RealK,         &
        (0.0_RealK, i=1, NPD_gauss_ord-8),                              &
!     Ninth order:
      -9.68160239507626E-01_RealK, -8.36031107326636E-01_RealK,         &
        -6.13371432700590E-01_RealK, -3.24253423403809E-01_RealK,       &
        0.0_RealK, 3.24253423403809E-01_RealK,                          &
        6.13371432700590E-01_RealK, 8.36031107326636E-01_RealK,         &
        9.68160239507626E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-9),  &
!     Tenth order:
      -9.73906528517172E-01_RealK, -8.65063366688985E-01_RealK,         &
        -6.79409568299024E-01_RealK, -4.33395394129427E-01_RealK,       &
        -1.48874338981631E-01_RealK, 1.48874338981631E-01_RealK,        &
        4.33395394129427E-01_RealK, 6.79409568299024E-01_RealK,         &
        8.65063366688985E-01_RealK, 9.73906528517172E-01_RealK          &
      /), (/ NPD_gauss_ord, NPD_gauss_ord /) )

  REAL (RealK), PARAMETER ::                                            &
!   Weights for Gaussian integration
    gauss_weight(npd_gauss_ord, npd_gauss_ord) = RESHAPE( (/            &
!     First order:
      2.0_RealK, (0.0_RealK, i=1, NPD_gauss_ord-1),                     &
!     Second order:
      1.0_RealK, 1.0_RealK, (0.0_RealK, i=1, NPD_gauss_ord-2),          &
!     Third order:
      5.55555555555556E-01_RealK, 8.88888888888889E-01_RealK,           &
        5.55555555555556E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-3),  &
!     Fourth order:
      3.47854845137454E-01_RealK, 6.52145154862546E-01_RealK,           &
        6.52145154862546E-01_RealK, 3.47854845137454E-01_RealK,         &
        (0.0_RealK, i=1, NPD_gauss_ord-4),                              &
!     Fifth order:
      2.36926885056189E-01_RealK, 4.78628670499366E-01_RealK,           &
        4.67913934572691E-01_RealK, 4.78628670499366E-01_RealK,         &
        2.36926885056189E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-5),  &
!     Sixth order:
      1.71324492379170E-01_RealK, 3.60761573048139E-01_RealK,           &
        4.67913934572691E-01_RealK, 4.67913934572691E-01_RealK,         &
        3.60761573048139E-01_RealK, 1.71324492379170E-01_RealK,         &
        (0.0_RealK, i=1, NPD_gauss_ord-6),                              &
!     Seventh order:
      1.29484966168870E-01_RealK, 2.79705391489277E-01_RealK,           &
        3.81830050505119E-01_RealK, 4.17959183673469E-01_RealK,         &
        3.81830050505119E-01_RealK, 2.79705391489277E-01_RealK,         &
        1.29484966168870E-01_RealK, (0.0_RealK, i=1, NPD_gauss_ord-7),  &
!     Eighth order:
      1.01228536290376E-01_RealK, 2.22381034453374E-01_RealK,           &
        3.13706645877887E-01_RealK, 3.62683783378362E-01_RealK,         &
        3.62683783378362E-01_RealK, 3.13706645877887E-01_RealK,         &
        2.22381034453374E-01_RealK, 1.01228536290376E-01_RealK,         &
        (0.0_RealK, i=1, NPD_gauss_ord-8),                              &
!     Ninth order:
      8.1274388361574E-02_RealK, 1.80648160694857E-01_RealK,            &
        2.60610696402935E-01_RealK, 3.12347077040003E-01_RealK,         &
        3.30239355001260E-01_RealK, 3.12347077040003E-01_RealK,         &
        2.60610696402935E-01_RealK, 1.80648160694857E-01_RealK,         &
        8.1274388361574E-02_RealK, (0.0_RealK, i=1, NPD_gauss_ord-9),   &
!     Tenth order:
      6.6671344308688E-02_RealK, 1.49451349150581E-01_RealK,            &
        2.19086362515982E-01_RealK, 2.69266719309996E-01_RealK,         &
        2.95524224714753E-01_RealK, 2.95524224714753E-01_RealK,         &
        2.69266719309996E-01_RealK, 2.19086362515982E-01_RealK,         &
        1.49451349150581E-01_RealK, 6.6671344308688E-02_RealK           &
      /), (/ NPD_gauss_ord, NPD_gauss_ord /) )

END MODULE gaussian_weight_pcf
