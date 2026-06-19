! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module flexchem_test_mod

use flexchem_kinds_mod, only: wp

contains

subroutine compare_number_density(n_calc, n_truth)

use flexchem_mod, only: network_name

implicit none

real(wp), intent(in) :: n_calc(:)
real(wp), intent(in) :: n_truth(:)

integer :: i
real(wp), parameter :: abs_tol = 1.0e-8_wp
real(wp) :: max_err

max_err = 0.0_wp
do i = 1, size(n_calc)
  max_err = max(max_err, abs(n_calc(i) - n_truth(i)))
end do

if (max_err < abs_tol) then
  write(*, "(a40, e10.3)") "Test "//trim(network_name)//": passed at abs_tol=", abs_tol
else
  write(*, "(a40, e10.3, a14, e10.3)") "Test "//trim(network_name)// &
  ": failed at abs_tol=", abs_tol, " with max_err=", max_err
end if

end subroutine compare_number_density

subroutine test_heng2017(t, n_init, n_final)

implicit none

integer, intent(in) :: t
real(wp), intent(in) :: n_init(:)
real(wp), intent(out) :: n_final(:)
! Local
real(wp) :: n01, n02
real(wp) :: C1, C2
real(wp), parameter :: k1 = 100.0_wp
real(wp), parameter :: k2 = 1.0_wp
real(wp), parameter :: J = 1.0_wp

n01 = n_init(1) + n_init(2)
n02 = n_init(1) + 2 * n_init(2)
C1 = 2*k1 + k2
C2 = -(k1+k2+J)

! n_A
n_final(1) = 2*n01*exp(-J*t) - (1/exp(-C2*t)) * ( n02 - (C1*n01*exp(-(C2+J)*t))/(C2+J) + C1*n01/(C2+J) )

! n_B
n_final(2) = (1/exp(-C2*t)) * ( n02 - C1*n01*exp(-(C2+J)*t)/(C2+J) + C1*n01/(C2+J) ) - n01*exp(-J*t)

! n_C
n_final(3) = ( 2*n01 + n01*C1/(C2+J) ) * ( 1 - exp(-J*t) ) + ( J*n02/C2 + n01*C1*J/((C2+J)*C2) ) * ( 1 - exp(C2*t) )

! n_D
n_final(4) = ( n01 + n01*C1/(C2+J)) * (exp(-J*t) - 1) + ( J*n02/C2 + n01*C1*J/((C2+J)*C2) ) * (exp(C2*t) - 1)

end subroutine test_heng2017

end module flexchem_test_mod
