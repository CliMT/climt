! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to concatenate LbL absorption coefficient netCDF files
program lblcat

use netcdf
use realtype_rd, only: RealK

implicit none

integer :: dimid1, dimid2, dimid3
integer :: varid

integer :: i, j, ios, n_file_in
character(len=256), allocatable :: file_in(:)
character(len=256) :: file_out
character(len=256) :: kabs_units

integer :: ncidin_lbl, ncidout_lbl
integer :: n_pt_pair, n_nu_tot, n_gas_frac, n_nu_written
integer, allocatable :: n_nu(:), nu_order(:)
real(RealK), allocatable :: nu(:), kabs(:, :), kabs_self(:, :, :)
real(RealK), allocatable :: p_calc(:), t_calc(:)
real(RealK), allocatable :: gas_frac(:)
real(RealK), allocatable :: first_nu(:), last_nu(:)
real(RealK) :: nu_inc, nu_inc_next, previous_nu
logical :: l_self_broadening

real(RealK), parameter :: eps = sqrt(epsilon(1.0_RealK))

interface
  subroutine map_heap_func(a, map)
    use realtype_rd, only: RealK
    real(RealK), intent(in), dimension(:) :: a
    integer, intent(out), dimension(:) :: map
  end subroutine map_heap_func
end interface

n_file_in = command_argument_count() - 1
if (n_file_in < 1) then
  stop 'Usage: lblcat <file1> <file2> ... <file_out>'
end if

allocate(file_in(n_file_in))
allocate(n_nu(n_file_in))
allocate(first_nu(n_file_in))
allocate(last_nu(n_file_in))
allocate(nu_order(n_file_in))

do i=1, n_file_in
  call get_command_argument(i, file_in(i), status=ios)
end do
call get_command_argument(n_file_in+1, file_out, status=ios)

! Open the first file for reading
call nf(nf90_open(trim(file_in(1)),NF90_NOWRITE,ncidin_lbl))

! Get length of dimensions
call nf(nf90_inq_dimid(ncidin_lbl, 'nu', dimid1))
call nf(nf90_inquire_dimension(ncidin_lbl, dimid1, len=n_nu(1)))
call nf(nf90_inq_dimid(ncidin_lbl, 'pt_pair', dimid2))
call nf(nf90_inquire_dimension(ncidin_lbl, dimid2, len=n_pt_pair))
ios = nf90_inq_dimid(ncidin_lbl, 'gas_frac', dimid3)
if (ios == NF90_NOERR) then
  call nf(nf90_inquire_dimension(ncidin_lbl, dimid3, len=n_gas_frac))
  l_self_broadening = .true.
else
  l_self_broadening = .false.
end if

! Read and get step in wavenumber array 
call nf(nf90_inq_varid(ncidin_lbl,'nu',varid))
call nf(nf90_get_att(ncidin_lbl,varid,'step',nu_inc))
allocate(nu(n_nu(1)))
call nf(nf90_get_var(ncidin_lbl,varid,nu))
first_nu(1)=nu(1)
last_nu(1)=nu(n_nu(1))
deallocate(nu)

! Read pressures and temperatures
allocate(p_calc(n_pt_pair))
allocate(t_calc(n_pt_pair))
call nf(nf90_inq_varid(ncidin_lbl,'p_calc',varid)) 
call nf(nf90_get_var(ncidin_lbl,varid,p_calc))
call nf(nf90_inq_varid(ncidin_lbl,'t_calc',varid))
call nf(nf90_get_var(ncidin_lbl,varid,t_calc))

! Read gas fractions
if (l_self_broadening) then
  allocate(gas_frac(n_gas_frac))
  call nf(nf90_inq_varid(ncidin_lbl,'gas_frac',varid)) 
  call nf(nf90_get_var(ncidin_lbl,varid,gas_frac))
end if

! Read units of k_abs
call nf(nf90_inq_varid(ncidin_lbl,'kabs',varid))
call nf(nf90_get_att(ncidin_lbl,varid,'units',kabs_units))

! Close file
call nf(nf90_close(ncidin_lbl))

! Find total number of wavenumbers in all files
do i=2, n_file_in
  call nf(nf90_open(trim(file_in(i)),NF90_NOWRITE,ncidin_lbl))
  call nf(nf90_inq_dimid(ncidin_lbl, 'nu', dimid1))
  call nf(nf90_inquire_dimension(ncidin_lbl, dimid1, len=n_nu(i)))
  call nf(nf90_inq_varid(ncidin_lbl,'nu',varid))
  call nf(nf90_get_att(ncidin_lbl,varid,'step',nu_inc_next))
  if (abs(nu_inc_next-nu_inc) > eps) then
    stop 'Error: files have different wavenumber increments.'
  end if
  allocate(nu(n_nu(i)))
  call nf(nf90_get_var(ncidin_lbl,varid,nu))
  first_nu(i)=nu(1)
  last_nu(i)=nu(n_nu(i))
  deallocate(nu)
  call nf(nf90_close(ncidin_lbl))
end do
n_nu_tot = sum(n_nu(1:n_file_in))

! Determine wavenumber order of files and check they are sequential
call map_heap_func(first_nu, nu_order)
previous_nu=last_nu(nu_order(1))
do j=2, n_file_in
  i=nu_order(j)
write(*,'(3f16.4)') previous_nu, first_nu(i), nu_inc
  if (abs(previous_nu+nu_inc-first_nu(i)) > eps) then
    stop 'Error: cannot determine sequential file order.'
  end if
  previous_nu=last_nu(i)
end do

! Open output file
call nf(nf90_create(trim(file_out),NF90_NOCLOBBER,ncidout_lbl))

! Create dimensions
call nf(nf90_def_dim(ncidout_lbl, 'nu', n_nu_tot, dimid1))
call nf(nf90_def_dim(ncidout_lbl, 'pt_pair', n_pt_pair, dimid2))
if (l_self_broadening) &
  call nf(nf90_def_dim(ncidout_lbl, 'gas_frac', n_gas_frac, dimid3))

! Create variables and write p_calc and t_calc
call nf(nf90_def_var(ncidout_lbl, 'p_calc', NF90_FLOAT, dimid2, varid))
call nf(nf90_put_att(ncidout_lbl, varid, 'title', 'pressure' ))
call nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'pressure' ))
call nf(nf90_put_att(ncidout_lbl, varid, 'units', 'Pa' ))
call nf(nf90_enddef(ncidout_lbl))
call nf(nf90_put_var(ncidout_lbl, varid, p_calc))
call nf(nf90_redef(ncidout_lbl))

call nf(nf90_def_var(ncidout_lbl, 't_calc', NF90_FLOAT, dimid2, varid))
call nf(nf90_put_att(ncidout_lbl, varid, 'title', 'temperature' ))
call nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'temperature' ))
call nf(nf90_put_att(ncidout_lbl, varid, 'units', 'K' ))
call nf(nf90_enddef(ncidout_lbl))
call nf(nf90_put_var(ncidout_lbl, varid, t_calc))
call nf(nf90_redef(ncidout_lbl))

call nf(nf90_def_var(ncidout_lbl, 'nu', NF90_DOUBLE, dimid1, varid))
call nf(nf90_put_att(ncidout_lbl, varid, 'title', 'wavenumber'))
call nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'wavenumber'))
call nf(nf90_put_att(ncidout_lbl, varid, 'units', 'm-1'))
call nf(nf90_put_att(ncidout_lbl, varid, 'step', nu_inc))

if (l_self_broadening) then
  call nf(nf90_def_var(ncidout_lbl, 'gas_frac', NF90_FLOAT, dimid3, varid))
  call nf(nf90_put_att(ncidout_lbl, varid, 'title', 'gas fraction'))
  call nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'gas fraction'))
  call nf(nf90_enddef(ncidout_lbl))
  call nf(nf90_put_var(ncidout_lbl, varid, gas_frac))
  call nf(nf90_redef(ncidout_lbl))

  call nf(nf90_def_var(ncidout_lbl, 'kabs', NF90_FLOAT, &
    (/dimid1,dimid2,dimid3/), varid))
else
  call nf(nf90_def_var(ncidout_lbl, 'kabs', NF90_FLOAT, &
    (/dimid1,dimid2/), varid))
end if
call nf(nf90_put_att(ncidout_lbl, varid, 'title', 'absorption' ))
call nf(nf90_put_att(ncidout_lbl, varid, 'long_name', 'absorption' ))
call nf(nf90_put_att(ncidout_lbl, varid, 'units', trim(kabs_units) ))

CALL nf(nf90_enddef(ncidout_lbl))

! Initialise number of wavenumbers written to file
n_nu_written=0

! Read each file and write data to output file sequentially
do j=1, n_file_in
  i = nu_order(j)
  write(*,'(a)') trim(file_in(i))
  call nf(nf90_open(trim(file_in(i)),NF90_NOWRITE,ncidin_lbl))

  call nf(nf90_inq_varid(ncidin_lbl,'nu',varid))
  allocate(nu(n_nu(i)))
  call nf(nf90_get_var(ncidin_lbl,varid,nu))

  call nf(nf90_inq_varid(ncidin_lbl,'kabs',varid))
  if (l_self_broadening) then
    allocate(kabs_self(n_nu(i), n_pt_pair, n_gas_frac))
    call nf(nf90_get_var(ncidin_lbl,varid,kabs_self))
  else
    allocate(kabs(n_nu(i), n_pt_pair))
    call nf(nf90_get_var(ncidin_lbl,varid,kabs))
  end if

  call nf(nf90_close(ncidin_lbl))

  call nf(nf90_inq_varid(ncidout_lbl,'nu',varid))
  call nf(nf90_put_var(ncidout_lbl,varid,nu, &
    start=(/n_nu_written+1/), count=(/n_nu(i)/)))

  call nf(nf90_inq_varid(ncidout_lbl,'kabs',varid))
  if (l_self_broadening) then
    call nf(nf90_put_var(ncidout_lbl,varid,kabs_self, &
      start=(/n_nu_written+1, 1, 1/), &
      count=(/n_nu(i) ,n_pt_pair, n_gas_frac/)))
  else
    call nf(nf90_put_var(ncidout_lbl,varid,kabs, &
      start=(/n_nu_written+1, 1/), count=(/n_nu(i), n_pt_pair/)))
  end if

  if (l_self_broadening) then
    deallocate(kabs_self)
  else
    deallocate(kabs)
  end if
  deallocate(nu)
  n_nu_written=n_nu_written+n_nu(i)
end do

call nf(nf90_close(ncidout_lbl))
write(*,'(a)') 'Written concatenated file: '//trim(file_out)

contains

subroutine nf(status)
  use netcdf
  integer, intent(in):: status
  if (status /= NF90_NOERR) then
     write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
     stop 'STOPPED!'
  end if
end subroutine nf

end program lblcat
