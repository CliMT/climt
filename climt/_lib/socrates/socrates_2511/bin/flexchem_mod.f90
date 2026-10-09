! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module flexchem_mod

use flexchem_kinds_mod, only: wp

implicit none

! Network name
character(len=10) :: network_name

! Names of molecules
character(len=10),allocatable :: molnames(:)
integer :: nmol

! Network parameters
integer :: nfile                         ! number of reaction files
integer :: nreac                         ! number of reactions
integer, parameter :: max_len_name = 10  ! max length of molecule name
integer, allocatable :: reac_method(:)   ! method to calculate reaction rate
integer, allocatable :: is_reversible(:) ! reaction reversible (-1) or not (1)

! Thermochemistry
integer, parameter :: nreactants = 5 ! max number of reactants
integer, parameter :: nproducts = 5  ! max number of products
integer, parameter :: nko = 5        ! number of parameters for ko
integer, parameter :: nkinf = 5      ! number of parameters for kinf
integer, parameter :: ntroe = 4      ! number of parameters for troe
integer, parameter :: nsri = 5       ! number of parameters for sri
integer, parameter :: neff = 22      ! number of efficiency species
integer, parameter :: nmod = 2       ! number of parameters for mod_param

integer, allocatable :: imol_reac(:,:) ! ids of the reactants
integer, allocatable :: imol_prod(:,:) ! ids of the products
integer, allocatable :: nb_reac(:)     ! number of reactants in each reaction
integer, allocatable :: nb_prod(:)     ! number of products in each reaction
integer, allocatable :: mod_param(:,:) ! contains file id of reactions

integer :: nmol_neq
integer :: neff_network ! actual number of efficiency molecules in network

integer, allocatable :: eff_list(:)   ! indices of molecules included as efficienies
integer, allocatable :: eff_list_mol(:)

real(wp), allocatable :: ko_param(:,:)   ! low pressure rate parameter
real(wp), allocatable :: kinf_param(:,:) ! high pressure rate parameter
real(wp), allocatable :: troe_param(:,:) ! troe parameters
real(wp), allocatable :: eff_param(:,:)  ! efficiency parameters, third body reactions
real(wp), allocatable :: sri_param(:,:)  ! SRI formalism parameters

character(max_len_name), allocatable :: reactants(:,:) ! reactant names
character(max_len_name), allocatable :: products(:,:)  ! product names
character(max_len_name) :: effmol(neff) ! names of molecules which act as third bodies

! Number of NASA polynomial coefficients
integer, parameter :: nparam  = 21
! NASA polynomial coefficients
real(wp), allocatable :: mol_param(:, :) ! nparam, nmol

! Number densities
real(wp), allocatable :: n_init(:)
real(wp), allocatable :: n_final(:)

real(wp) :: tmax_kinetics


contains

  subroutine flexchem_alloc

    implicit none

    ! Number densities
    allocate(n_init(nmol_neq))
    allocate(n_final(nmol_neq))

    ! Network parameters
    allocate(mol_param(nparam, nmol))
    allocate(molnames(nmol))

    ! Parameters for thermochemical reactions
    allocate( reactants( nreac, nreactants ) )
    allocate( products( nreac, nproducts ) )
    allocate( ko_param( nreac, nko ) )
    allocate( kinf_param( nreac, nkinf ) )
    allocate( troe_param( nreac, ntroe ) )
    allocate( eff_param( nreac, neff ) )
    allocate( sri_param( nreac, nsri ) )
    allocate( mod_param( nreac, nmod ) )
    allocate( imol_reac( nreac, nreactants ) )
    allocate( imol_prod( nreac, nproducts ) )
    allocate( nb_reac( nreac ) )
    allocate( nb_prod( nreac ) )
    allocate( eff_list( neff ) )
    allocate( eff_list_mol( neff ) )

  end subroutine flexchem_alloc


  subroutine flexchem_dealloc

    implicit none

    ! Number densities
    deallocate(n_init)
    deallocate(n_final)

    ! Network parameters
    deallocate(mol_param)
    deallocate(molnames)
    deallocate(reac_method)
    deallocate(is_reversible)

    ! Deallocate arrays for thermochemical reactions
    deallocate(reactants)
    deallocate(products)
    deallocate(ko_param)
    deallocate(kinf_param)
    deallocate(troe_param)
    deallocate(eff_param)
    deallocate(sri_param)
    deallocate(mod_param)
    deallocate(imol_reac)
    deallocate(imol_prod)
    deallocate(nb_reac)
    deallocate(nb_prod)
    deallocate(eff_list)
    deallocate(eff_list_mol)

  end subroutine flexchem_dealloc
end module flexchem_mod
