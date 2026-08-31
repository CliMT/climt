MODULE coord_transforms

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE CartCoord
  REAL (RealK) :: x
  REAL (RealK) :: y
  REAL (RealK) :: z
END TYPE CartCoord

TYPE SphCoord
  REAL (RealK) :: r
  REAL (RealK) :: theta
  REAL (RealK) :: phi
END TYPE SphCoord


CONTAINS


SUBROUTINE cart2sph(cart, sph)

IMPLICIT NONE

TYPE (CartCoord), INTENT(IN)  :: cart
TYPE (SphCoord),  INTENT(OUT) :: sph

sph%r = SQRT(cart%x**2 + cart%y**2 + cart%z**2)
sph%theta = ACOS(cart%z / sph%r)
sph%phi = ATAN2(cart%y, cart%x)

END SUBROUTINE cart2sph


SUBROUTINE sph2cart(sph, cart)

IMPLICIT NONE

TYPE (SphCoord),  INTENT(IN)  :: sph
TYPE (CartCoord), INTENT(OUT) :: cart

cart%x = sph%r*SIN(sph%theta)*COS(sph%phi)
cart%y = sph%r*SIN(sph%theta)*SIN(sph%phi)
cart%z = sph%r*COS(sph%theta)

END SUBROUTINE sph2cart


SUBROUTINE rotate_x(ang, co_in, co_out)

! Rotate about the x axis

IMPLICIT NONE

REAL (RealK),     INTENT(IN)  :: ang
TYPE (CartCoord), INTENT(IN)  :: co_in
TYPE (CartCoord), INTENT(OUT) :: co_out

co_out%x = co_in%x
co_out%y = COS(ang)*co_in%y + SIN(ang)*co_in%z
co_out%z = -SIN(ang)*co_in%y + COS(ang)*co_in%z

END SUBROUTINE rotate_x


SUBROUTINE rotate_z(ang, co_in, co_out)

! Rotate about the z axis

IMPLICIT NONE

REAL (RealK),     INTENT(IN)  :: ang
TYPE (CartCoord), INTENT(IN)  :: co_in
TYPE (CartCoord), INTENT(OUT) :: co_out

co_out%x = COS(ang)*co_in%x + SIN(ang)*co_in%y
co_out%y = -SIN(ang)*co_in%x + COS(ang)*co_in%y
co_out%z = co_in%z

END SUBROUTINE rotate_z


SUBROUTINE icrf2ecliptic(obl, icrf, ecliptic)

! Convert from the Internation Celestial Reference frame (Earth Equatorial)
! to the Ecliptic coordinate system (Earth orbit)

IMPLICIT NONE

REAL (RealK),     INTENT(IN)  :: obl   ! Obliquity of the ecliptic
TYPE (CartCoord), INTENT(IN)  :: icrf
TYPE (CartCoord), INTENT(OUT) :: ecliptic

ecliptic%x = icrf%x
ecliptic%y = COS(obl)*icrf%y + SIN(obl)*icrf%z
ecliptic%z = -SIN(obl)*icrf%y + COS(obl)*icrf%z

END SUBROUTINE icrf2ecliptic


END MODULE coord_transforms
