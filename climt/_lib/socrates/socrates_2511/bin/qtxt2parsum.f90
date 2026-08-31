! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to convert HITRAN isotope partition sum files to parsum.dat
program qtxt2parsum

use realtype_rd, only: RealK
use hitran_cnst

implicit none

integer :: ierr, ios, i, j, k, iu_qtxt, iu_parsum, temp
real(RealK) :: qtemp
character(len=256) :: filename
character(len=256) :: parsum_file='parsum_2021.dat'
character(len=10) :: isotope_names(number_species)
logical :: l_exist

! Read q???.txt files
call get_free_unit(ierr, iu_qtxt)
do i = 1, number_species
  write (filename, '(a,i0,a)') 'q', global_id(i), '.txt'
  inquire(FILE=trim(filename), EXIST=l_exist)
  if (l_exist) then
    open(UNIT=iu_qtxt, FILE=trim(filename), iostat=ios, status='old')
    do
      read(iu_qtxt, '(i4,f20.5)', iostat=ios) temp, qtemp
      if (ios /= 0) exit
      qcoeff(temp:, i) = qtemp
    end do
    close(iu_qtxt)
  else
    qcoeff(:,i) = q296(i)
  end if
end do

i=0
do j = 1, number_molecules
  do k = 1, number_isotopes(j)
    i = i + 1
    ! Constuct isotope names
    write(isotope_names(i), '(a,a,i0)') molecule_names(j), '_', global_id(i)
  end do
end do


! Write parsum.dat file
call get_free_unit(ierr, iu_parsum)
open(UNIT=iu_parsum, FILE=trim(parsum_file), iostat=ios, status='new')
write(iu_parsum, '(a7,9x,a10,999(17x,a10))') 'Temp(K)', isotope_names
do i=1, number_temperatures
  write(iu_parsum, '(f6.1,f20.6,999f27.6)') real(i), qcoeff(i, :)
end do
close(iu_parsum)

end program qtxt2parsum
