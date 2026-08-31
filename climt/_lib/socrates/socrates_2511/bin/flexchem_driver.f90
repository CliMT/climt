! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
program flexchem_driver

  use realtype_rd, only: realext

  use flexchem_kinds_mod, only: wp
  use flexchem_kinetics_mod, only: kinetics_solver
  use flexchem_mod, only: flexchem_alloc, flexchem_dealloc, &
    network_name, molnames, tmax_kinetics, &
    n_init, n_final, nmol_neq
  use flexchem_input_mod, only: read_molecule_names, read_nasa_coeffs, &
    read_network_params, read_reactions, read_number_density
  use flexchem_test_mod, only: compare_number_density, test_heng2017

  implicit none

  character(len=132) :: dir_common
  character(len=132) :: dir_network
  character(len=132) :: fnasa7pol, fnasa9pol

  integer :: i

  real(realext) :: pressure, temperature
  real(wp) :: pres_wp, temp_wp
  real(wp), allocatable :: number_density(:)

  dir_common = "."
  pressure = 1.0e5_realext ! [Pa]
  temperature = 2500.0_realext ! [K]

  fnasa7pol = trim(adjustl(dir_common))//"/coeff7_NASA_sc.dat"
  fnasa9pol = trim(adjustl(dir_common))//"/coeff9_NASA_sc.dat"
  call get_command_argument(1, network_name)
  network_name = trim(network_name)
  dir_network = trim(adjustl(dir_common))//"/"//trim(adjustl(network_name))//"/"

  call read_network_params(dir_network)

  call flexchem_alloc
  allocate(number_density(nmol_neq))

  call read_molecule_names(dir_network)
  call read_reactions(dir_network)
  call read_nasa_coeffs(fnasa7pol, fnasa9pol)

  call read_number_density(dir_network)

  number_density = real(n_init, kind=wp)
  pres_wp = real(pressure, kind=wp)
  temp_wp = real(temperature, kind=wp)

  if ( network_name == "heng2017" ) then
    ! Heng2017 analytical test
    tmax_kinetics = 6 ! [s]
    call test_heng2017(6, n_init, n_final)
    ! call compare_number_density(number_density, n_final)
    write(*, "(a35)") "Test: heng2017 (analytic solution)"
    write(*, "(a5, a16)") "i", "number_density"
    do i = 1, nmol_neq
      write(*, "(i5, e16.7)") i, n_final(i)
    end do
  else if ( network_name == "venot2019" ) then
    ! Venot2019 CHON reduced chemical network
    tmax_kinetics = 1e15 ! [s]
    ! call compare_number_density(number_density, n_final)
  else
    stop "ERROR: invalid network_name"
  end if

  call kinetics_solver(pres_wp, temp_wp, number_density)
  ! Print out final values
  write(*, "(a20)") "Test: "//network_name
  write(*, "(a5, a16)") "i", "number_density"
  do i = 1, nmol_neq
    write(*, "(i5, e16.7)") i, number_density(i)
  end do

  call flexchem_dealloc
  deallocate(number_density)

end program flexchem_driver
