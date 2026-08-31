! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convert water vapour self-continuum data from Jon Elsey in cm2 molec-1 atm-1
! to Socrates format files in cm3 molec-1 (atm-1).
program convert_elsey_shine

use realtype_rd, only: RealK
use hitran_cnst, only: c2

implicit none

integer, parameter :: n=2001
real(RealK), parameter :: cont_inc = 10.0_RealK ! cm-1
character(len=*), parameter :: file_in(2) &
  = (/'SOC_new_296K.dat','SOC_new_260K.dat'/)
character(len=*), parameter :: file_out(2) &
  = (/'elsey_shine_s296','elsey_shine_s260'/)
character(len=*), parameter :: file_plot(2) &
  = (/'elsey_shine_s296_plot', 'elsey_shine_s260_plot'/)
character(len=*), parameter :: caviar_in(2) &
  = (/'caviar_s296', 'caviar_s260'/)
character(len=*), parameter :: caviar_plot(2) &
  = (/'caviar_s296_plot', 'caviar_s260_plot'/)
character(len=*), parameter :: mt_ckd_in(2) &
  = (/'mt_ckd3p2_s296', 'mt_ckd3p2_s260'/)
character(len=*), parameter :: mt_ckd_plot(2) &
  = (/'mt_ckd3p2_s296_plot', 'mt_ckd3p2_s260_plot'/)
real(RealK), parameter :: t(2) = (/ 296.0_RealK, 260.0_RealK /)

integer :: i, l, ios, iu_cont=30
real(RealK) :: data(2, n)
real(RealK) :: c_self(n+2)
character(len=80) :: line

do l=1, 2
  ! Read data file from John Elsey / Keith Shine
  open(unit=iu_cont, file=file_in(l), status='old')
  read(iu_cont, *) line ! header line
  read(iu_cont, *) data
  close(iu_cont)

  ! Write out file without header (for quick read in python)
  open(unit=iu_cont, file=file_plot(l), status='unknown', iostat=ios)
  do i=1, n
    write(iu_cont, '(2(1pe11.4))') data(:,i)
  end do
  close(iu_cont)

  ! Convert from standard format for continua: cm2 molec-1 atm-1 to
  ! format used with MT_CKD and Socrates: cm3 molec-1 (atm-1) * 1.0e-20
  ! Conversion involves dividing by the "radiation field term" which
  ! Socrates will include in the routine h2o_continuum.f90.
  c_self(4:) = data(2,2:)*1.0E+20_RealK / &
    (data(1,2:) * TANH(data(1,2:) * 100.0_RealK*c2/(2.0_RealK*t(l))))
  c_self(3) = 0.0_RealK
  c_self(2) = c_self(4)
  c_self(1) = c_self(5)

  ! Write out in Socrates continuum format 
  open(unit=iu_cont, file=file_out(l), status='unknown')
  write(iu_cont, '(A/)') '*BEGIN_DATA'
  write(iu_cont, '(A, 1PE12.5, A)') '     Start of table                = ', &
    data(1,1)-cont_inc*2.0_RealK, '     cm-1'
  write(iu_cont, '(A, 1PE12.5, A)') '     End of table                  = ', &
    data(1,n), '     cm-1'
  write(iu_cont, '(A, 1PE12.5, A)') '     Increment in table            = ', &
    cont_inc, '     cm-1'
  write(iu_cont, '(A, I5/)') '     Number of points in table     = ', n+2
  write(iu_cont, '(A)') '     Continuum coefficients (cm^3.mol-1 * 1.0E-20)'
  write(iu_cont, '(4(3X, 1PE12.5))') c_self
  write(iu_cont, '(A)') '*END'
  close(iu_cont)
end do

! Convert Socrates CAVIAR files into cm2 molec-1 atm-1 for plotting
do l=1, 2
  open(unit=iu_cont, file=caviar_in(l), status='old')
  read(iu_cont, '(////////, 4(3X, 1E12.5))') c_self
  close(iu_cont)
  c_self(4:) = c_self(4:) * &
   (data(1,2:) * TANH(data(1,2:) * 100.0_RealK*c2/(2.0_RealK*t(l)))) &
   / 1.0E+20_RealK
  
  open(unit=iu_cont, file=caviar_plot(l), status='unknown')
  do i=1,n
    write(iu_cont, '(2(1pe11.4))') data(1,i), c_self(i+2)
  end do
  close(iu_cont)
end do

! Convert Socrates MT_CKD files into cm2 molec-1 atm-1 for plotting
do l=1, 2
  open(unit=iu_cont, file=mt_ckd_in(l), status='old')
  read(iu_cont, '(////////, 4(3X, 1E12.5))') c_self
  close(iu_cont)
  c_self(4:) = c_self(4:) * &
   (data(1,2:) * TANH(data(1,2:) * 100.0_RealK*c2/(2.0_RealK*t(l)))) &
   / 1.0E+20_RealK
  
  open(unit=iu_cont, file=mt_ckd_plot(l), status='unknown')
  do i=1,n
    write(iu_cont, '(2(1pe11.4))') data(1,i), c_self(i+2)
  end do
  close(iu_cont)
end do

end program
