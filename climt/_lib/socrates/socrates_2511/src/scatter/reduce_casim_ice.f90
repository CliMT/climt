! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program to reduce the CASIM ice database to a manageable size
!
! Method:
!   Takes as input the large database of ice optical properties calculated
!   for the CASIM cloud ice and snow size distributions within Anthony Baran's
!   ensemble ice optics model.
!   Repeatedly bisects the data into parts containing equal numbers of data
!   points in order of mass mixing ratio and mean particle mass.
!   Outputs a reduced database of 1026 blocks of meaned properties in a format
!   readable by the scatter_average routine.
!
!------------------------------------------------------------------------------
program reduce_casim_ice

use realtype_rd, only: RealK
use file_type_pcf, only: it_file_scat_mass
use scatter_algorithm_pcf, only: ip_scat_database, name_scatter_algorithm
use scatter_pp_pcf, only: ip_type_ice, name_scatter_type
use rad_ccf, only: pi

implicit none

character (len=256) :: filename(2)
integer :: filename_len(2)

integer, parameter :: n_wavelength = 169
integer, parameter :: n_bisects = 5
integer, parameter :: n_block_edge = 2
integer, parameter :: n_block_out = 4**n_bisects + n_block_edge

integer :: n_block_in
integer :: n_edge      ! n_block_in/(2*n_block_out)
integer :: n_to_bisect ! n_block_in - n_block_edge*n_edge

real (RealK), allocatable :: mass_mixing_ratio(:) ! (n_block_in)
real (RealK), allocatable :: mean_mass(:)         ! (n_block_in)
real (RealK), allocatable :: air_density(:)       ! (n_block_in)
real (RealK), allocatable :: absorption(:, :) ! (n_wavelength, n_block_in)
real (RealK), allocatable :: scattering(:, :) ! (n_wavelength, n_block_in)
real (RealK), allocatable :: asymmetry(:, :)  ! (n_wavelength, n_block_in)

real (RealK) :: number_by_mass
real (RealK) :: wavelength(n_wavelength)

real (RealK) :: mmr_out(n_block_out)
real (RealK) :: mass_out(n_block_out)
real (RealK) :: r_eqv(n_block_out)
real (RealK) :: air_density_out(n_block_out)
real (RealK) :: ice_density(n_block_out)
real (RealK) :: absorption_out(n_wavelength, n_block_out)
real (RealK) :: scattering_out(n_wavelength, n_block_out)
real (RealK) :: asymmetry_out(n_wavelength, n_block_out)
real (RealK) :: n_data(n_block_out)

integer, allocatable :: block_map(:)            ! (n_block_in)
integer, allocatable :: mmr_map(:), mass_map(:) ! (n_block_in)
integer :: ierr, ios, iu_db, iu_scatter
integer :: i, j, l, n, a, b

interface
  subroutine map_heap_func(a, map)
    use realtype_rd
    real (RealK), intent(in) :: a(:)
    integer, intent(out) :: map(:)
  end subroutine map_heap_func
end interface


! Get filenames from command-line arguments
do i=1, 2
  call get_command_argument(i, filename(i), filename_len(i), ios)
  if (ios > 0) then
    write(*, *) 'Usage: reduce_casim_ice <input database> <output database>'
    stop
  else if (ios < 0) then
    write(*,*) 'Error: File name too long (>256 characters)'
    stop
  end if
end do

! Count number of blocks in database
call get_free_unit(ierr, iu_db)
open(unit=iu_db, file=filename(1)(1:filename_len(1)), iostat=ios, status='old')
n_block_in = 0
outer: do
  do i = 1, n_wavelength + 4
    read(iu_db, *, iostat=ios)
    if (ios /= 0) exit outer
  end do
  n_block_in = n_block_in + 1
end do outer
close(iu_db)
write(*, '(a,i0, a)') 'Reading ', n_block_in, ' blocks'

allocate(mass_mixing_ratio(n_block_in))
allocate(mean_mass(n_block_in))
allocate(air_density(n_block_in))
allocate(absorption(n_wavelength, n_block_in))
allocate(scattering(n_wavelength, n_block_in))
allocate(asymmetry(n_wavelength, n_block_in))

! Read database
call get_free_unit(ierr, iu_db)
open(unit=iu_db, file=filename(1)(1:filename_len(1)), iostat=ios, status='old')
do i = 1, n_block_in
  read(iu_db, '(16x,f20.8)') mass_mixing_ratio(i)
  read(iu_db, '(16x,f20.8)') number_by_mass
  read(iu_db, '(28x,f20.8)') air_density(i)
  mean_mass(i) = mass_mixing_ratio(i)/number_by_mass
! Uncomment to debug:
!  write(101, '(i5, 3(4x, 1pe16.9))') i, mass_mixing_ratio(i), mean_mass(i), &
!    air_density(i)
  read(iu_db, *)
  do j = 1, n_wavelength
    read(iu_db, '(3E15.6, F15.4)')  wavelength(j), &
      scattering(j, i), absorption(j, i), asymmetry(j, i)
  end do
end do
close(iu_db)

allocate(block_map(n_block_in))
allocate(mmr_map(n_block_in))
allocate(mass_map(n_block_in))

! Determine the order of the blocks in terms of mass mixing ratio and mass.
call map_heap_func(mass_mixing_ratio, mmr_map)
call map_heap_func(mean_mass, mass_map)

! First remove the smallest mass particles from the data and place them into
! their own output blocks in order to ensure a good fit where the optical
! properties are changing most rapidly.
n_edge = n_block_in/(2*n_block_out)
block_map(:) = 0
do i = 1, n_block_edge
  block_map(mass_map(n_edge*(i-1)+1:n_edge*i)) = i - 1 - n_block_edge
end do

! Reduce the remaining data by bisecting into parts containing equal numbers
! of data points. This proceeds by alternately splitting each region into two
! parts by order of mass mixing ratio and then particle mass.
n_to_bisect = n_block_in - n_block_edge*n_edge
do i = 1, n_bisects
  ! Each i loop bisects the data by MMR and then by particle mass.
  do j = 1, 4**(i-1)
    ! Each j loop runs through the blocks to be bisected
    a = 0
    do n = 1, n_block_in
      ! Find the data currently in this j block and loop
      ! through it in order of increasing MMR
      l = mmr_map(n)
      if ( block_map(l) == j-1 ) then
        ! Move the first half of the data in this block into a new block
        block_map(l) = block_map(l) + 4**(i-1)
        a = a + 1
      end if
      if (a > n_to_bisect/(2*4**(i-1))) exit
    end do
    a = 0
    b = 0
    do n = 1, n_block_in
      ! Now loop through the data in the block that was just split
      ! in order of increasing mass.
      l = mass_map(n)
      ! For both parts of the split block move the first half
      ! of the data into a new block
      if ( block_map(l) == j-1 .and. a <= n_to_bisect/(4**i)) then
        block_map(l) = block_map(l) + 2*4**(i-1)
        a = a + 1
      else if ( block_map(l) == j-1+4**(i-1) .and. b <= n_to_bisect/(4**i)) then
        block_map(l) = block_map(l) + 2*4**(i-1)
        b = b + 1
      end if
      if (a > n_to_bisect/(4**i) .and. b > n_to_bisect/(4**i)) exit
    end do
    ! The j block has now been split into four parts
  end do
  ! All the j blocks have now been split into four parts
end do
! Start the block indexing from 1 rather than -n_block_edge
block_map(:) = block_map(:) + n_block_edge + 1

deallocate(mass_map)
deallocate(mmr_map)

! Average the data in each of the output block regions
ice_density(:) = 917.0_RealK
air_density_out(:) = 0.0_RealK
n_data(:) = 0.0_RealK
mmr_out(:) = 0.0_RealK
mass_out(:) = 0.0_RealK
absorption_out(:, :) = 0.0_RealK
scattering_out(:, :) = 0.0_RealK
asymmetry_out(:, :) = 0.0_RealK
do i = 1, n_block_in
  n_data(block_map(i)) = n_data(block_map(i)) + 1.0_RealK
  mmr_out(block_map(i)) = mmr_out(block_map(i)) + mass_mixing_ratio(i)
  mass_out(block_map(i)) = mass_out(block_map(i)) + mean_mass(i)
  air_density_out(block_map(i)) = air_density_out(block_map(i)) + air_density(i)
  ! Average data in units of m2/kg
  absorption_out(:, block_map(i)) = absorption_out(:, block_map(i)) &
    + absorption(1:n_wavelength, i) / (mass_mixing_ratio(i)*air_density(i))
  scattering_out(:, block_map(i)) = scattering_out(:, block_map(i)) &
    + scattering(:, i) / (mass_mixing_ratio(i)*air_density(i))
  asymmetry_out(:, block_map(i)) = asymmetry_out(:, block_map(i)) &
    + asymmetry(:, i) * scattering(:, i) &
    / (mass_mixing_ratio(i)*air_density(i))
end do

deallocate(block_map)

deallocate(asymmetry)
deallocate(scattering)
deallocate(absorption)
deallocate(air_density)
deallocate(mean_mass)
deallocate(mass_mixing_ratio)

! Write output
call get_free_unit(ierr, iu_scatter)
open(unit=iu_scatter, file=filename(2)(1:filename_len(2)), &
     iostat=ios, status='unknown')
do i = 1, n_block_out
  mass_out(i) = mass_out(i)/n_data(i)
  mmr_out(i) = mmr_out(i)/n_data(i)
  r_eqv(i) = (3.0*mass_out(i)/(4.0*pi*ice_density(i)))**(1.0/3.0)
  air_density_out(i) = air_density_out(i)/n_data(i)
! Uncomment to debug:
!  write(102, '(2i5, 3(4x, 1pe16.9))') &
!    i, int(n_data(i)), mmr_out(i), mass_out(i), r_eqv(i)
  absorption_out(:, i) = absorption_out(:, i)/n_data(i)
  asymmetry_out(:, i) = asymmetry_out(:, i)/scattering_out(:, i)
  scattering_out(:, i) = scattering_out(:, i)/n_data(i)
  ! Convert back to m-1
  absorption_out(:, i) = absorption_out(:, i)*mmr_out(i)*air_density_out(i)
  scattering_out(:, i) = scattering_out(:, i)*mmr_out(i)*air_density_out(i)
  write(iu_scatter, '(/, a13, i5)') '*FILE TYPE = ', it_file_scat_mass
  write(iu_scatter, '(/, 16x, a, //)') 'Single scattering parameters.'
  write(iu_scatter, '(a23, i5, 3x, a1, a40, a1)') &
    'Scattering algorithm = ', ip_scat_database, &
    '(', name_scatter_algorithm(ip_scat_database), ')'
  write(iu_scatter, '(a23, i5, 3x, a1, a20, a1)') &
    'Scattering type      = ', ip_type_ice, &
    '(', name_scatter_type(ip_type_ice), ')'
  write(iu_scatter, '(/, a30, 1pe12.5, 1x, a5)') &
    'Mass mixing ratio           = ', mmr_out(i), 'kg/kg'
  write(iu_scatter, '(a30, 1pe12.5, 1x, a2)') &
    'Mean particle mass          = ', mass_out(i), 'kg'
  write(iu_scatter, '(a30, 1pe12.5, 1x, a5)') &
    'Air density                 = ', air_density_out(i), 'kg/m3'
  write(iu_scatter, '(a30, 1pe12.5, 1x, a5)') &
    'Ice density                 = ', ice_density(i), 'kg/m3'
  write(iu_scatter, '(a30, 1pe12.5, 1x, a1)') &
    'Equivalent spherical radius = ', r_eqv(i), 'm'
  write(iu_scatter, '(a40, i3)') &
    'Number of terms in the phase function = ', 1
  write(iu_scatter, '(/, 4x, a14, 6x, a16, 4x, a16, 4x, a9)') &
    'Wavelength (m)','Absorption (m-1)','Scattering (m-1)','Asymmetry'
  do j = 1, n_wavelength
    write(iu_scatter, '(4(4x, 1pe16.9))') wavelength(j), &
      absorption_out(j, i), scattering_out(j, i), asymmetry_out(j, i)
  end do
  write(iu_scatter, '(1x)')
end do
close(iu_scatter)

end program reduce_casim_ice
