! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate pressure and temperature interpolation factor for ESFT terms
!
! Method:
!   Input variables: P, T, N_PROFILE, N_LAYER
!   Output variables: FAC00, FAC01, FAC10, FAC11, JP, JT
!   Zhian Sun   JUN.  1997
!
!-----------------------------------------------------------------------

MODULE inter_pt_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTER_PT_MOD'
CONTAINS
SUBROUTINE inter_pt(nd_profile, nd_layer                                &
     , n_profile, n_layer, gas_mix_ratio                                &
     , p, t, fac00, fac01, fac10, fac11                                 &
     , fac00c, fac01c, fac10c, fac11c                                   &
     , facc00, facc01                                                   &
     , jp, jph2oc, jt, jtt, jto2c, jtswo3 )

  USE realtype_rd, ONLY: RealK
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


  INTEGER ::                                                            &
       nd_profile                                                       &
!       Max number of profile
     , nd_layer                                                         &
!       Max number of layer
     , n_layer                                                          &
!       Number of layers
     , n_profile                                                        &
!       Number of profiles
     , jp(nd_profile, nd_layer), jp1                                    &
!       Index of reference pressure level such that the actual
!       pressure is between JP and JP1
     , jph2oc(nd_profile, nd_layer)                                     &
!       Index of Ho2 continuum reference water vapour pressure at level I
!       such that the actual wvp is between JPHO2C and JPHO2C+1
     , jt(nd_profile, nd_layer)                                         &
!       Index of reference temperature at level I such that the actual
!       temperature is between JT and JT+1
     , jtt(nd_profile, nd_layer)                                        &
!       Index of reference temperature at level I+1 such that the actual
!       temperature is between JT and JT+1
     , jto2c(nd_profile, nd_layer)                                      &
!       Index of o2 continuum  reference temperature at level I
!       such that the actual temperature is between JTO2C and JTO2C+1
     , jtswo3(nd_profile, nd_layer)
!       Index of sw o3 reference temp

  REAL (RealK) ::                                                       &
       p(nd_profile, nd_layer)                                          &
!          Actual pressure
     , t(nd_profile, nd_layer)                                          &
!          Layer temperature
     , gas_mix_ratio(nd_profile, nd_layer)                              &
!          Mixing ratio of water vapour
     , fac00(nd_profile, nd_layer)                                      &
     , fac01(nd_profile, nd_layer)                                      &
     , fac10(nd_profile, nd_layer)                                      &
     , fac11(nd_profile, nd_layer)                                      &
!          Multiplication factors for P & T interpolation
     , fac00c(nd_profile, nd_layer)                                     &
     , fac01c(nd_profile, nd_layer)                                     &
     , fac10c(nd_profile, nd_layer)                                     &
     , fac11c(nd_profile, nd_layer)                                     &
!          Multiplication factors for P & T interpolation of h2o
     , facc00(nd_profile, nd_layer)                                     &
     , facc01(nd_profile, nd_layer)                                     &
!          Multiplication factors for continuum T interpolation
     , facswoz00(nd_profile, nd_layer)                                  &
     , facswoz01(nd_profile, nd_layer)
!          Multiplication factors for sw o3 T interpolation

! Local variables
  REAL ::                                                               &
       p_reflog(59)                                                     &
!          Log of p_ref
     , t_refc(3)                                                        &
!          Reference temperature for continuum
     , t_ref_swo3(5)                                                    &
!          Reference temperature for sw o3 & h2o continuum
!          interpolation (abs k)
     , wvp_ref(21)                                                      &
!          Reference water vapour pressure
     , fp, fph2oc, ft, ft1, compfp, compft                              &
!          Fraction factor for P & T interpolation
     , plog                                                             &
!          Inter-medium variable
     , dt, dtc, dtswo3                                                  &
!          Fraction of reference temperature difference
     , dwvp, wvp                                                        &
!          Fraction of reference water vapour pressure
     , ratio_air_h2o
!          ratio of molecular air weight to h2o weight
!          for converting mass mixing ratio to volumn mixing ratio

REAL, PARAMETER ::       &
!          Reference temperature
    t_ref(59, 5) =  RESHAPE( [     &

    ! These are temperature profile corresponding to MLS values
    ! (middle column) +/- 15 and +/- 30 K

      1.8434e+02, 1.9934e+02, 2.1434e+02, 2.2934e+02, 2.4434e+02,           &
      1.8466e+02, 1.9966e+02, 2.1466e+02, 2.2966e+02, 2.4466e+02,           &
      1.8498e+02, 1.9998e+02, 2.1498e+02, 2.2998e+02, 2.4498e+02,           &
      1.8530e+02, 2.0030e+02, 2.1530e+02, 2.3030e+02, 2.4530e+02,           &
      1.8563e+02, 2.0063e+02, 2.1563e+02, 2.3063e+02, 2.4563e+02,           &
      1.8595e+02, 2.0095e+02, 2.1595e+02, 2.3095e+02, 2.4595e+02,           &
      1.8627e+02, 2.0127e+02, 2.1627e+02, 2.3127e+02, 2.4627e+02,           &
      1.8661e+02, 2.0161e+02, 2.1661e+02, 2.3161e+02, 2.4661e+02,           &
      1.8703e+02, 2.0203e+02, 2.1703e+02, 2.3203e+02, 2.4703e+02,           &
      1.8805e+02, 2.0305e+02, 2.1805e+02, 2.3305e+02, 2.4805e+02,           &
      1.9045e+02, 2.0545e+02, 2.2045e+02, 2.3545e+02, 2.5045e+02,           &
      1.9409e+02, 2.0909e+02, 2.2409e+02, 2.3909e+02, 2.5409e+02,           &
      1.9809e+02, 2.1309e+02, 2.2809e+02, 2.4309e+02, 2.5809e+02,           &
      2.0222e+02, 2.1722e+02, 2.3222e+02, 2.4722e+02, 2.6222e+02,           &
      2.0646e+02, 2.2146e+02, 2.3646e+02, 2.5146e+02, 2.6646e+02,           &
      2.1068e+02, 2.2568e+02, 2.4068e+02, 2.5568e+02, 2.7068e+02,           &
      2.1503e+02, 2.3003e+02, 2.4503e+02, 2.6003e+02, 2.7503e+02,           &
      2.1949e+02, 2.3449e+02, 2.4949e+02, 2.6449e+02, 2.7949e+02,           &
      2.2397e+02, 2.3897e+02, 2.5397e+02, 2.6897e+02, 2.8397e+02,           &
      2.2852e+02, 2.4352e+02, 2.5852e+02, 2.7352e+02, 2.8852e+02,           &
      2.3320e+02, 2.4820e+02, 2.6320e+02, 2.7820e+02, 2.9320e+02,           &
      2.3790e+02, 2.5290e+02, 2.6790e+02, 2.8290e+02, 2.9790e+02,           &
      2.4242e+02, 2.5742e+02, 2.7242e+02, 2.8742e+02, 3.0242e+02,           &
      2.4552e+02, 2.6052e+02, 2.7552e+02, 2.9052e+02, 3.0552e+02,           &
      2.4532e+02, 2.6032e+02, 2.7532e+02, 2.9032e+02, 3.0532e+02,           &
      2.4288e+02, 2.5788e+02, 2.7288e+02, 2.8788e+02, 3.0288e+02,           &
      2.3956e+02, 2.5456e+02, 2.6956e+02, 2.8456e+02, 2.9956e+02,           &
      2.3580e+02, 2.5080e+02, 2.6580e+02, 2.8080e+02, 2.9580e+02,           &
      2.3198e+02, 2.4698e+02, 2.6198e+02, 2.7698e+02, 2.9198e+02,           &
      2.2816e+02, 2.4316e+02, 2.5816e+02, 2.7316e+02, 2.8816e+02,           &
      2.2434e+02, 2.3934e+02, 2.5434e+02, 2.6934e+02, 2.8434e+02,           &
      2.2063e+02, 2.3563e+02, 2.5063e+02, 2.6563e+02, 2.8063e+02,           &
      2.1712e+02, 2.3212e+02, 2.4712e+02, 2.6212e+02, 2.7712e+02,           &
      2.1381e+02, 2.2881e+02, 2.4381e+02, 2.5881e+02, 2.7381e+02,           &
      2.1062e+02, 2.2562e+02, 2.4062e+02, 2.5562e+02, 2.7062e+02,           &
      2.0754e+02, 2.2254e+02, 2.3754e+02, 2.5254e+02, 2.6754e+02,           &
      2.0447e+02, 2.1947e+02, 2.3447e+02, 2.4947e+02, 2.6447e+02,           &
      2.0141e+02, 2.1641e+02, 2.3141e+02, 2.4641e+02, 2.6141e+02,           &
      1.9852e+02, 2.1352e+02, 2.2852e+02, 2.4352e+02, 2.5852e+02,           &
      1.9594e+02, 2.1094e+02, 2.2594e+02, 2.4094e+02, 2.5594e+02,           &
      1.9391e+02, 2.0891e+02, 2.2391e+02, 2.3891e+02, 2.5391e+02,           &
      1.9226e+02, 2.0726e+02, 2.2226e+02, 2.3726e+02, 2.5226e+02,           &
      1.9070e+02, 2.0570e+02, 2.2070e+02, 2.3570e+02, 2.5070e+02,           &
      1.8913e+02, 2.0413e+02, 2.1913e+02, 2.3413e+02, 2.4913e+02,           &
      1.8760e+02, 2.0260e+02, 2.1760e+02, 2.3260e+02, 2.4760e+02,           &
      1.8635e+02, 2.0135e+02, 2.1635e+02, 2.3135e+02, 2.4635e+02,           &
      1.8581e+02, 2.0081e+02, 2.1581e+02, 2.3081e+02, 2.4581e+02,           &
      1.8575e+02, 2.0075e+02, 2.1575e+02, 2.3075e+02, 2.4575e+02,           &
      1.8580e+02, 2.0080e+02, 2.1580e+02, 2.3080e+02, 2.4580e+02,           &
      1.8710e+02, 2.0210e+02, 2.1710e+02, 2.3210e+02, 2.4710e+02,           &
      1.9321e+02, 2.0821e+02, 2.2321e+02, 2.3821e+02, 2.5321e+02,           &
      2.0187e+02, 2.1687e+02, 2.3187e+02, 2.4687e+02, 2.6187e+02,           &
      2.1090e+02, 2.2590e+02, 2.4090e+02, 2.5590e+02, 2.7090e+02,           &
      2.2025e+02, 2.3525e+02, 2.5025e+02, 2.6525e+02, 2.8025e+02,           &
      2.2974e+02, 2.4474e+02, 2.5974e+02, 2.7474e+02, 2.8974e+02,           &
      2.3916e+02, 2.5416e+02, 2.6916e+02, 2.8416e+02, 2.9916e+02,           &
      2.4880e+02, 2.6380e+02, 2.7880e+02, 2.9380e+02, 3.0880e+02,           &
      2.5827e+02, 2.7327e+02, 2.8827e+02, 3.0327e+02, 3.1827e+02,           &
      2.6400e+02, 2.7900e+02, 2.9400e+02, 3.0900e+02, 3.2400e+02 &
    ], [59, 5], order=[2,1] )

  INTEGER ::                                                            &
       i, l
!       Vertical and horizontal loop index

! These three temperature are reference for o2 continuum
  DATA t_refc /200., 250., 300./

! These are reference T for sw o3 & h2o continuum k interpolation
  DATA t_ref_swo3/180., 215., 250., 285., 320/

! Reference water vapour pressures
  DATA wvp_ref/0.005,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0     &
              ,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,47.5,50.0/

  PARAMETER(dt=1./15, dtc=0.02, dtswo3=1./35, dwvp=0.4                  &
           , ratio_air_h2o=1.60777e-02)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='INTER_PT'

! ratio_air_h2o=0.01*28.964/18.015 e-2 for changing
! partial pressure to hPa

! DTC=0.02              ! 1/50
! DTSWO3=1.0/35.        ! 1/35
! DWVP=1/2.5

! These ln(pressures) are chosen such that ln(PRE(1)) = -4.64
! ln(PRE(59)) = 6.96000) and
! each subsequent ln(pressure) differs from the previous one by 0.2.
! The order of pressure is reversed
!     1.05363E+03,8.62642E+02,7.06272E+02,5.78246E+02,4.73428E+02,
!     3.87610E+02,3.17348E+02,2.59823E+02,2.12725E+02,1.74164E+02,
!     1.42594E+02,1.16746E+02,9.55835E+01,7.82571E+01,6.40715E+01,
!     5.24573E+01,4.29484E+01,3.51632E+01,2.87892E+01,2.35706E+01,
!     1.92980E+01,1.57998E+01,1.29358E+01,1.05910E+01,8.67114E+00,
!     7.09933E+00,5.81244E+00,4.75882E+00,3.89619E+00,3.18993E+00,
!     2.61170E+00,2.13828E+00,1.75067E+00,1.43333E+00,1.17351E+00,
!     9.60789E-01,7.86628E-01,6.44036E-01,5.27292E-01,4.31710E-01,
!     3.53455E-01,2.89384E-01,2.36928E-01,1.93980E-01,1.58817E-01,
!     1.30029E-01,1.06458E-01,8.71608E-02,7.13612E-02,5.84256E-02,
!     4.78349E-02,3.91639E-02,3.20647E-02,2.62523E-02,2.14936E-02,
!     1.75975E-02,1.44076E-02,1.17959E-02,9.65769E-03
  DATA p_reflog/                                                        &
  -4.64000e+00,-4.44000e+00,-4.24000e+00,-4.04000e+00,-3.84000e+00,     &
  -3.64000e+00,-3.44000e+00,-3.24000e+00,-3.04000e+00,-2.84000e+00,     &
  -2.64000e+00,-2.44000e+00,-2.24000e+00,-2.04000e+00,-1.84000e+00,     &
  -1.64000e+00,-1.44000e+00,-1.24000e+00,-1.04000e+00,-8.40001e-01,     &
  -6.40001e-01,-4.40001e-01,-2.40000e-01,-4.00004e-02, 1.59999e-01,     &
   3.60000e-01, 5.59999e-01, 7.60002e-01, 9.60001e-01, 1.16000e+00,     &
   1.36000e+00, 1.56000e+00, 1.76000e+00, 1.96000e+00, 2.16000e+00,     &
   2.36000e+00, 2.56000e+00, 2.76000e+00, 2.96000e+00, 3.16000e+00,     &
   3.36000e+00, 3.56000e+00, 3.76000e+00, 3.96000e+00, 4.16000e+00,     &
   4.36000e+00, 4.56000e+00, 4.76000e+00, 4.96000e+00, 5.16000e+00,     &
   5.36000e+00, 5.56000e+00, 5.76000e+00, 5.96000e+00, 6.16000e+00,     &
   6.36000e+00, 6.56000e+00, 6.76000e+00, 6.96000e+00/

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  DO i=1, n_layer
    DO l=1, n_profile

! Find out two reference temperature on either side of the
! layer temperature. Store these indices in JTO2C for o2 continuum
! interpolation

      jto2c(l,i)=INT(1.+dtc*(t(l,i)-200.0))
      jto2c(l,i) = MIN(2, MAX(1, jto2c(l,i)) )
      facc00(l,i)=dtc*(t_refc(jto2c(l,i)+1)-t(l,i))
      facc01(l,i)=1.-facc00(l,i)

! Find out two reference temperature on either side of the
! layer temperature. Store these indices in JTSWO3 for
! sw ozone absorption and water vapor continuum interpolation

      jtswo3(l,i)=INT(1.+dtswo3*(t(l,i)-180.0))
      jtswo3(l,i) = MIN( 4, MAX(1, jtswo3(l,i)) )
      facswoz00(l,i)=dtswo3*(t_ref_swo3(jtswo3(l,i)+1)-t(l,i))
      facswoz01(l,i)=1.-facswoz00(l,i)
    END DO

! Find out two reference pressures on either side of the layer
! pressure. Store them in JP and JP1. Store in FP the fraction
! of the difference (in ln(pressure)) between these
! two values that the layer pressure lies for gas absorption
! coefficient interpolation.

    DO l=1, n_profile
      plog=LOG(p(l,i)/100.0)
      jp(l,i)=INT(24+5*(plog+0.04))
      jp(l,i)=MAX(1, MIN(58, jp(l,i) ) )
      jp1=jp(l,i)+1
      fp=5.0*(p_reflog(jp1)-plog)

! For water vapour, pressure is replaced by vapour pressure wvp

      wvp = ratio_air_h2o * p(l,i) * gas_mix_ratio(l,i)
      jph2oc(l,i) = INT(1 + dwvp * wvp )
      jph2oc(l,i) = MAX(1, MIN(20, jph2oc(l,i) ) )
      fph2oc = dwvp * (wvp_ref(jph2oc(l,i)+1) - wvp)

! Find out two reference temperature on either side of the
! layer temperature. Store these indices in JT and JTT, resp.
! Store in FT (resp. FT1) the fraction of the way between JT
! (JTT) and the next highest reference temperature that the
! layer temperature falls.

      jt(l,i)=INT(3. + dt*(t(l,i)-t_ref(jp(l,i),3)))

      jt(l,i) = MAX(1, MIN(4, jt(l,i)) )


      ft=dt*(t_ref(jp(l,i),jt(l,i)+1)-t(l,i))
      jtt(l,i)=INT(3. + dt*(t(l,i)-t_ref(jp1, 3)))
      jtt(l,i) = MAX(1, MIN(4, jtt(l,i)) )

      ft1=dt*(t_ref(jp1, jtt(l,i)+1)-t(l,i))

! Multiply the pressure fraction with the approperiate temperature
! fraction  for use in the interpolations

      compfp=1. - fp
      compft=1. - ft
      fac00(l,i)=fp*ft
      fac10(l,i)=compfp*ft1
      fac01(l,i)=fp*compft
      fac11(l,i)=compfp*(1.-ft1)

! Water vapour continuum

      fac00c(l,i)=fph2oc*facswoz00(l,i)
      fac10c(l,i)=(1.-fph2oc)*facswoz00(l,i)
      fac01c(l,i)=fph2oc*facswoz01(l,i)
      fac11c(l,i)=(1.-fph2oc)*(1.-facswoz00(l,i))
    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE inter_pt
END MODULE inter_pt_mod
