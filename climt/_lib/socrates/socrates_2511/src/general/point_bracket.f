! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to bracket a point in an array of wavelengths.
!
! Method:
!   An array of points is searched to find those bracketing
!   a given wavelength.
!
!- ---------------------------------------------------------------------
      SUBROUTINE point_bracket(lambda
     &  , n_points, spectrum
     &  , i_short, i_long)
!
!
!     Modules to set types of variables:
      USE realtype_rd
!
!
      IMPLICIT NONE
!
!
!     Include headers.
!
!     Dummy arguments
      INTEGER, Intent(IN) ::
     &    n_points
!           Number of points
      INTEGER, Intent(OUT) ::
     &    i_short
!           Lower bracketing point
     &  , i_long
!           Upper bracketing point
      REAL  (RealK), Intent(IN) ::
     &    lambda
!           Point to search for
     &  , spectrum(n_points)
!           Spectrum to search
!
!     Local variables.
      INTEGER
     &    i_mid
!           Mid-point of interval
!
!
!     Initialize the upper and lower points and then 
!     begin a binary search.
      i_short=0
      i_long=n_points+1
1     if (i_long-i_short.gt.1) then
!       Interval requires chopping.
        i_mid=(i_long+i_short)/2
        IF (spectrum(i_mid) <= lambda) THEN
          i_short=i_mid
        ENDIF
        IF (spectrum(i_mid) > lambda) THEN
          i_long=i_mid
        ENDIF
        goto 1
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE point_bracket
