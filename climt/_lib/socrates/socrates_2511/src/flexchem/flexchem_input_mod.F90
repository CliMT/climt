! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module flexchem_input_mod

use flexchem_kinds_mod, only: wp
use flexchem_mod, only: nmol, molnames, &
  mol_param, nmol_neq, neff_network, nreac, ko_param, kinf_param, troe_param, &
  eff_param, sri_param, imol_reac, imol_prod, nb_reac, nb_prod, mod_param, &
  eff_list, eff_list_mol, n_init

implicit none

contains


subroutine read_molecule_names(dir_network)
  ! Read in the list of molecules from a text file

  character(len=132), intent(in) :: dir_network
  integer :: i, iunit

  open(newunit=iunit, file=trim(dir_network)//"composes.dat")
  read(iunit, *)
  do i = 1, nmol
    read(iunit, '(a10)') molnames(i)
  enddo
  close(iunit)

end subroutine read_molecule_names


subroutine read_number_density(dir_network)
  ! Read in number density from a text file for testing

  character(len=132), intent(in) :: dir_network
  integer :: i, iunit

  open(newunit=iunit, file=trim(dir_network)//"number_density.dat")
  read(iunit, *)
  do i = 1, nmol_neq
    read(iunit, '(e16.7)') n_init(i)
  enddo
  close(iunit)

end subroutine read_number_density


subroutine read_network_params(dir_network)

  use flexchem_mod, only: nfile, reac_method, is_reversible

  implicit none

  character(len=132), intent(in) :: dir_network

  ! local variables
  integer :: stat
  integer :: i, ireac, iunit
  character(2) :: dummy, ifile

  ! Read in network parameters
  open(newunit=iunit, file=trim(dir_network)//'network_params.dat')
  read(iunit, *) ! Top header
  read(iunit, *) ! Next header
  read(iunit, '(i10)') nfile
  read(iunit, *) ! Next header
  read(iunit, '(i10)') nmol
  read(iunit, *) ! Next header
  read(iunit, '(i10)') nmol_neq
  read(iunit, *) ! Next header

  allocate(reac_method(nfile))
  allocate(is_reversible(nfile))

  do i = 1, nfile
    read(iunit, '(11x, i10, i10)') reac_method(i), is_reversible(i)
  end do
  close(iunit)

  ! count number of reactions
  nreac = 0
  do i = 1, nfile
    write(ifile, '(i2)') i
    open(newunit=iunit, &
      file=trim(dir_network)//'reactions_'//trim(adjustl(ifile))//'.dat')
    do
      read(iunit, '(a)', iostat=stat) dummy
      if (stat < 0) exit
      nreac = nreac + 1
    end do
    close(iunit)
  end do

end subroutine read_network_params


subroutine read_reactions(dir_network)

use flexchem_mod, only: nfile, reac_method, is_reversible, &
  nreactants, nproducts, nko, nkinf, ntroe, neff, nsri, nmod, &
  reactants, products, effmol

implicit none

character(len=132), intent(in) :: dir_network

! local variables
integer :: stat
integer :: i, j, ireac, iunit
character(2) :: ifile

! initialize to 0
ko_param  = 0. ; kinf_param = 0. ; troe_param = 0. ; sri_param = 0.
eff_param = 0. ; mod_param  = 0
imol_reac = 0 ; imol_prod = 0
nb_reac   = 0 ; nb_prod   = 0
neff_network = 0 ; eff_list = 0 ; eff_list_mol = 0

! read in thermochemical reactions
ireac = 1
do i = 1, nfile
  write(ifile, '(i2)') i
  open(newunit=iunit, &
    file=trim(dir_network)//'reactions_'//trim(adjustl(ifile))//'.dat')
  do
    ! read in format depends on the method of calculation for reactions in file
    if (reac_method(i) == 1) then
      read(iunit, '(5(1x, a10),5(1x, a10),5(1x, e10.3))', end=1) &
        reactants(ireac, :), &
        products(ireac, :), &
        ko_param(ireac, :)
    else if (reac_method(i) == 2) then
      read(iunit, '(5(1x, a10),5(1x, a10),5(1x, e10.3),22(1x, e10.3))', end=1) &
        reactants(ireac, :), &
        products(ireac, :), &
        ko_param(ireac, :), &
        eff_param(ireac, :)
    else if (reac_method(i) == 3) then
      read(iunit, '(5(1x, a10),5(1x, a10),10(1x, e10.3),4(1x, e10.3),22(1x, e10.3))', end=1) &
        reactants(ireac, :), &
        products(ireac, :), &
        ko_param(ireac, :), &
        kinf_param(ireac, :), &
        troe_param(ireac, :), &
        eff_param(ireac, :)
    else if (reac_method(i) == 4) then
      read(iunit, '(5(1x, a10),5(1x, a10),10(1x, e10.3),5(1x, e10.3),22(1x, e10.3))', end=1) &
        reactants(ireac, :), &
        products(ireac, :), &
        ko_param(ireac, :), &
        kinf_param(ireac, :), &
        sri_param(ireac, :), &
        eff_param(ireac, :)
    end if

    mod_param(ireac, 1) = i
    mod_param(ireac, 2) = is_reversible(i)

    ireac = ireac + 1
  end do
1 close(iunit)
end do

! Molecules for efficiencies in three-body reactions
effmol = (/'O2   ', 'CO   ', 'CO2  ', 'H2O  ', 'CH4  ', 'H2   ', 'C2H6 ', 'Ar   ', 'N2   ', &
'He   ', 'C2H4 ', 'cC6H6', 'C7H8 ', 'H2O2 ', 'N2O  ', 'O-3P ', 'NH3  ', 'N2H4 ', 'N2O4 ', 'NO2  ', 'NO   ', 'N-4S '/)

! Build list of indices for efficiency molecules
do i = 1, neff
  if ( check_mol_name( effmol(i) ) ) then
    eff_list(i) = i
    eff_list_mol(i) = get_mol_idx( effmol(i) )
    neff_network = neff_network + 1
  end if
end do

! Map reactants and products
! and count number of products and reactants in each reaction
do ireac = 1, nreac
  do j = 1, nreactants
    i = get_mol_idx( reactants(ireac, j) )
    if ( i > 0 ) then
      nb_reac(ireac) = nb_reac(ireac) + 1
      imol_reac(ireac, j) = i
    end if
  end do

  do j = 1, nproducts
    i = get_mol_idx( products(ireac, j) )
    if ( i > 0 ) then
      nb_prod(ireac) = nb_prod(ireac) + 1
      imol_prod(ireac, j) = i
    end if
  end do
end do

end subroutine read_reactions


! READ in NASA polynomial coefficients
subroutine read_nasa_coeffs(fnasa7pol,fnasa9pol)

  implicit none

  character(len=132), intent(in)   :: fnasa7pol, fnasa9pol

  ! Local variables
  integer  :: i, iunit
  real(wp) :: tmin, tmax, tmean
  logical :: mol_inc
  character(len=10) :: name_loc

  mol_param = 0.

  do i =1,nmol
    mol_inc = .false.
    open(newunit=iunit, file=fnasa7pol)
    fileloop: do
    read(iunit, *, end=1) name_loc,tmin,tmax,tmean
      if (molnames(i) == name_loc) then
        mol_inc = .true.
        mol_param(1,i) = tmin
        mol_param(2,i) = tmax
        mol_param(3,i) = tmean
        read(iunit, *) mol_param(6:12,i)
        read(iunit, *) mol_param(15:21,i)
        exit fileloop
        endif
    enddo fileloop
  1 close(iunit)

  ! if not found in 7 coefficient file, check 9 coefficient file
  if (.not. mol_inc) then
    open(newunit=iunit, file=fnasa9pol)
    fileloop2: do
    read(iunit, *, end=2) name_loc,tmin,tmax,tmean
    if (molnames(i) == name_loc) then
      mol_inc = .true.
      mol_param(1,i) = tmin
      mol_param(2,i) = tmax
      mol_param(3,i) = tmean

      read(iunit, *) mol_param(4:12,i)
      read(iunit, *) mol_param(13:21,i)

      exit fileloop2
      endif

    enddo fileloop2
  2 close(iunit)
  end if

  enddo

end subroutine read_nasa_coeffs


function get_mol_idx(molrequest) result(idx)

  implicit none

  integer :: i
  integer :: idx
  character(*) :: molrequest

  idx = 0

  ! First, check for special cases before looping
  ! HV = photochemical reaction
  if (molrequest == 'HV        ' .or. molrequest == '          ') then
     idx = -1
     return
  end if

  ! Loop through molnames to find the matching molecule
  do i = 1, nmol
     if (molrequest == molnames(i)) then
        idx = i
        return  ! Exit immediately once a match is found
     end if
  end do

  ! Error handling if not found
  if (idx == 0) then
     write(*, *) 'Error in get_mol_idx: '//trim(adjustl(molrequest))//' not found'
     stop
  end if

end function get_mol_idx


! Check if requested molecule is in molecules list
function check_mol_name(molrequest) result(check)

implicit none

integer :: i
logical :: check
character(*) :: molrequest

check = .false.

do i = 1, nmol
  if (molrequest == molnames(i)) then
    check = .true.
    return  ! Exit immediately once a match is found
  end if
end do

end function check_mol_name


end module flexchem_input_mod
