! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module flexchem_kinetics_mod

use flexchem_kinds_mod, only: wp, qp
use flexchem_mod, only: nmol, nreac, reac_method, imol_reac, imol_prod, nb_prod, nb_reac, &
  mod_param, eff_list, eff_list_mol, ko_param, kinf_param, troe_param, eff_param, sri_param, &
  tmax_kinetics, mol_param, neff_network, nmol_neq

implicit none

real(wp), allocatable :: kf(:)           ! forward reaction rate
real(wp), allocatable :: kr(:)           ! reverse reaction rate

contains

  subroutine kinetics_solver(p, t, n)

    implicit none

    real(wp), intent(in)    :: p    ! pressure
    real(wp), intent(in)    :: t    ! temperature

    ! DLSODES required inputs (see opkdmain.f for details)
    real(wp), intent(inout) :: n(:) ! array of (initial) values of length node
    integer  :: node        ! number of first order ODEs
    real(wp) :: time        ! (initial) value of time
    real(wp) :: dt          ! first point where output is desired (timestep)
    integer  :: itol        ! 1 for atol as scalar, 2 for atol as array
    real(wp) :: rtol        ! relative tolerance parameter (scalar)
    real(wp) :: atol        ! absolute tolerance parameter (scalar or array)
    integer  :: itask       ! 1 for normal computation
    integer  :: istate      ! flag for the state of the calculation
    integer  :: iopt        ! flag for optional inputs
    real(wp), allocatable :: rwork(:) ! real work array
    integer  :: lrw         ! declared length of rwork
    integer, allocatable :: iwork(:) ! integer work array
    integer  :: liw         ! declared length of iwork
    integer  :: mf          ! method flag
    ! DLSODES optional inputs
    integer :: maxord ! max ODEs order to be allowed
    integer :: mxstep ! max number of internal steps allowed per call
    integer :: mxhnil ! max number of warnings printed to screen
    integer :: mferr  ! error messages flag

    ! Non-DLSODES parameters
    real(wp), parameter :: dtmin = 1.0e-10_wp ! minimum allowed timestep [s]
    real(wp), parameter :: dtmax = 1.0e10_wp  ! maximum allowed timestep [s]
    real(wp), parameter :: ttest = 1.1_wp     ! time factor
    real(wp), parameter :: nmin = 1.0e-100_wp ! minimum allowed number density
    integer,  parameter :: mod_jac = 50       ! Jacobian recalculation frequency
    integer,  parameter :: nit_max = 5000      ! maximum number of iterations
    real(wp), parameter :: dnn_tol   = 1e-2_wp ! tolerance for the change in number density
    real(wp), parameter :: dnndt_tol = 1e-4_wp ! tolerance for the change in number density with time
    real(wp), parameter :: mf_tol = 1e-30_wp   ! ?
    real(wp), parameter :: tmin_kinetics = 1.0_wp ! ? 1 second

    ! Non-DLSODES variables
    integer :: i
    integer :: nit_neq          ! number of iterations
    real(wp) :: told            ! time at a previous iteration
    real(wp) :: tvold           ! time at two previous iterations
    real(wp) :: nold(nmol_neq)  ! number density at a previous iteration
    real(wp) :: nvold(nmol_neq) ! number density at two previous iterations
    integer :: nbstep           ! number of steps taken for the problem so far
    integer :: nbstep_old       ! number of steps on the previous iteration
    real(wp) :: dnn             ! change in number density
    real(wp) :: dnndt           ! change in number density divided by change in time
    integer :: stop_flag        ! flag for stopping the computation

    logical :: l_continue

    ! DLSODES settings (see opkdmain.f for details)
    node   = nmol_neq
    itol   = 1            ! atol is a scalar
    rtol   = 1.0e-3_wp
    atol   = 1.0e-10_wp
    itask  = 1            ! normal computation
    iopt   = 1            ! optional inputs are used
    if (node < 5) then
        lrw = 20 + (2 + 1/2) * node*nreac*5 + (11 + 9/2) * node
        liw = 30
    else
        lrw    = 20 + 11*node + 2*node**2
        liw    = 30 + 2*node  ! must be at least 30
    end if
    mf     = 222          ! stiff problem with Jacobian internally generated
    maxord = 2            ! opt input: default is 12 or 5
    mxstep = 200          ! opt input: default is 500
    mxhnil = 0            ! opt input: default is 10
    mferr  = 1            ! opt input: 0 error messages off, 1 on

    ! Allocate arrays
    allocate(iwork(liw))
    allocate(rwork(lrw))
    allocate(kf(nreac))
    allocate(kr(nreac))

    ! Initialise variables
    ! Initialise DLSODES
    iwork = 0
    rwork = 0.0_wp
    iwork(5) = maxord ! apply optional setting to DLSODES call
    iwork(6) = mxstep ! apply optional setting to DLSODES call
    iwork(7) = mxhnil ! apply optional setting to DLSODES call
    istate   = 1      ! flag for first call
    ! Initialise other variables
    time     = 0.0_wp
    dt       = dtmin
    nit_neq  = 0
    told     = 0.0_wp
    tvold    = 0.0_wp
    nold     = 0.0_wp
    nvold    = 0.0_wp
    kf       = 0.0_wp
    kr       = 0.0_wp

    l_continue = .true.

    nbstep = 0
    nbstep_old = 0

    ! Calculate initial reaction rates
    call get_krkf(p, t, n)

    ! Start iterations
    do while (l_continue)

      ! Save previous number densities
      if (time >= told*ttest) then
        tvold = told
        told  = time
        nvold = nold
        nold  = n
      end if

      call xsetf(mferr)
      call dlsodes(dndt, node, n, time, time+dt, itol, rtol, atol, itask, &
                   istate, iopt, rwork, lrw, iwork, liw, jac, mf)

      if (nmin /= 0.0_wp) then
        where(n<nmin) n(:) = 0.0_wp
      end if

      nbstep_old = nbstep
      nbstep = iwork(11)

      nit_neq  = nit_neq  + 1

      ! If solver error or requested, recalculate on next iteration
      if (istate<0 .or. nbstep-nbstep_old>mxstep .or. mod(nit_neq, mod_jac)==0) then
        istate = 1
        call get_krkf(p, t, n)
        nbstep_old = 0
        nbstep = 0
      end if

      !calculate new timestep from lsode
      dt = 0.9_wp * rwork(11) ! rwork(11) is step size in time last used (successfully)

      ! apply limits to timestep
      dt = min(dtmax, dt)
      dt = max(dtmin, dt)

      ! ensure that we do not overshoot tmax
      if ((time+dt)>tmax_kinetics) then
        dt = tmax_kinetics - time
      end if

      ! find maximum variations in number densities
      ! loop over species and levels
      dnn   = 0.0_wp
      dnndt = 0.0_wp
      do i = 1, node
        if (n(i)>0.0_wp .and. nvold(i)>0.0_wp .and. n(i)/sum(n)>mf_tol) then
          dnn   = max( abs(n(i)-nvold(i))/n(i), dnn )
          dnndt = max( abs(n(i)-nvold(i))/n(i)/(time-tvold), dnndt )
        end if
      end do

      ! Stop check
      !   reset counter if criteria not met on this timestep
      if (dnn>dnn_tol .or. dnndt>dnndt_tol)  then
        stop_flag = 0
      end if
      !   if criteria met add one to counter
      if (time>tmin_kinetics .and. dnn<dnn_tol .and. dnndt<dnndt_tol) then
        stop_flag = stop_flag + 1
      end if
      !   if convergence criteria met for three consecutive time steps stop iterations
      !   or maximum time reached or maximum allowed iterations reached
      if (stop_flag==3 .or. time>=tmax_kinetics .or. nit_neq>nit_max) then
        l_continue = .false.
      end if

    end do

    if (nit_neq > nit_max) then
     write(*,*) 'ERROR: nit exceeded ',nit_max,' with t = ',time
    end if

    deallocate(iwork, rwork, kf, kr)

  end subroutine kinetics_solver


  subroutine get_krkf(p, t, n)
  ! ------------------------------------------------------------------------------
  ! Subroutine: get_krkf
  ! Purpose:
  !   This subroutine calculates the forward and reverse rate coefficients (kf and kr)
  !   for a set of chemical reactions based on various reaction rate formalisms.
  !   The method used depends on the `reac_method` and includes Kooij, Troe, SRI,
  !   and third-body effects. The subroutine also computes thermodynamic properties
  !   (chemical potentials, enthalpies, and heat capacities) for the reactions.
  !
  ! Arguments:
  !   - p (input): Real, pressure (in appropriate units such as Pa or dyn/cm^2).
  !   - t (input): Real, temperature (in Kelvin).
  !   - n (input): Real array, concentrations or number densities of species involved
  !     in the reactions.
  !
  ! Local Variables:
  !   - Various local variables to store intermediate calculation results, including
  !     rate coefficients, chemical potentials, thermodynamic properties, and reaction-specific
  !     parameters for Kooij, Troe, SRI, and third-body effects.
  !
  ! The subroutine computes:
  !   - Forward reaction rate coefficients (kf) using the appropriate formalism.
  !   - Reverse reaction rate coefficients (kr) by calculating the equilibrium constant
  !     and using the relation `kr = kf / K_eq`.
  !
  ! Reaction rate formalisms:
  !   - Kooij formalism
  !   - Kooij with third-body effects
  !   - Troe formalism
  !   - SRI formalism
  !
  ! Thermodynamic quantities (chemical potentials, enthalpies, heat capacities) are
  ! computed using the `get_thermodynamics` subroutine.
  ! ------------------------------------------------------------------------------

    implicit none

    ! Input arguments
    real(wp), intent(in) :: p        ! Pressure (Pa or appropriate units)
    real(wp), intent(in) :: t        ! Temperature (K)
    real(wp), intent(in) :: n(:)     ! Array of species concentrations (or number densities)

    ! Local variables
    integer :: i, ireac, idepth, ipr, exp_pt, iloc, istart
    real(wp) :: ko_loc, kinf_loc, tloc, ceff_loc, zfc, exp_zfc, c_zfc, n_zfc, d_zfc, &
               k_loc, nloc, nref, ntot
    real(qp) :: k_loc16
    real(wp) :: mu(nmol)            ! Chemical potential for species
    real(wp) :: h(nmol)             ! Enthalpy for species
    real(wp) :: cp(nmol)            ! Specific heat capacity for species
    real(wp) :: one_over_t          ! Inverse of temperature
    real(wp) :: t_over_300          ! Temperature ratio to 300 K

    ! Constants
    real(wp), parameter :: p0 = 1.013250e+6_wp  ! Reference pressure (dyn/cm^2)
    real(wp), parameter :: kb = 1.380620e-16_wp  ! Boltzmann constant (cgs units)

    ! Initialize variables
    one_over_t = 1.0_wp / t
    t_over_300 = t / 300.0_wp
    ntot = sum(n)
    nref = p0 * one_over_t / kb

    ! Calculate thermodynamic quantities
    call get_thermodynamics(t, mu, h, cp)

    ! Loop over all reactions
    do ireac = 1, nreac

      ! Calculate rate coefficients for both low and high pressure limits
      ko_loc = ko_param(ireac, 1) * t_over_300**ko_param(ireac, 2) * &
               exp(-ko_param(ireac, 3) * one_over_t)
      kinf_loc = kinf_param(ireac, 1) * (t_over_300)**kinf_param(ireac, 2) * &
                 exp(-kinf_param(ireac, 3) * one_over_t)

      ! Reaction method: Kooij formalism
      if (reac_method(mod_param(ireac, 1)) == 1) then
        kf(ireac) = ko_loc

      ! Kooij formalism with third body
      else if (reac_method(mod_param(ireac, 1)) == 2) then
        ! Calculate efficiencies for third-body effect
        ceff_loc = ntot
        do i = 1, neff_network
          if (eff_list(i) > 0) then
            ceff_loc = ceff_loc - n(eff_list_mol(i)) + eff_param(ireac, eff_list(i)) * n(eff_list_mol(i))
          end if
        end do
        kf(ireac) = ko_loc * ceff_loc

      ! Troe formalism
      else if (reac_method(mod_param(ireac, 1)) == 3) then
        ! Calculate efficiencies for Troe formalism
        ceff_loc = ntot
        do i = 1, neff_network
          if (eff_list(i) > 0) then
            ceff_loc = ceff_loc - n(eff_list_mol(i)) + eff_param(ireac, eff_list(i)) * n(eff_list_mol(i))
          end if
        end do
        ko_loc = ko_loc * ceff_loc
        if (ko_loc == 0.0_wp) then
          kf(ireac) = kinf_loc
        else if (kinf_loc == 0.0_wp) then
          kf(ireac) = ko_loc
        else
          zfc = exp(-troe_param(ireac, 4) / t)
          if (.not. (troe_param(ireac, 2) == 0.0_wp)) then
            zfc = zfc + (1.0_wp - troe_param(ireac, 1)) * exp(-t / troe_param(ireac, 2))
          end if
          if (.not. (troe_param(ireac, 3) == 0.0_wp)) then
            zfc = zfc + troe_param(ireac, 1) * exp(-t / troe_param(ireac, 3))
          end if

          c_zfc = -0.40_wp - 0.67_wp * log10(zfc)
          n_zfc = 0.75_wp - 1.27_wp * log10(zfc)
          d_zfc = 0.14_wp

          exp_zfc = 1.0_wp / (1.0_wp + ((log10(ko_loc / kinf_loc) + c_zfc) / &
            (n_zfc - d_zfc * (log10(ko_loc / kinf_loc) + c_zfc)))**2)

          kf(ireac) = (ko_loc / (1.0_wp + ko_loc / kinf_loc)) * zfc**exp_zfc
        end if

      ! SRI formalism
      else if (reac_method(mod_param(ireac, 1)) == 4) then
        ! Calculate efficiencies for SRI formalism
        ceff_loc = ntot
        do i = 1, neff_network
          if (eff_list(i) > 0) then
            ceff_loc = ceff_loc - n(eff_list_mol(i)) + eff_param(ireac, eff_list(i)) * n(eff_list_mol(i))
          end if
        end do
        ko_loc = ko_loc * ceff_loc
        if (ko_loc == 0.0_wp) then
          kf(ireac) = kinf_loc
        else if (kinf_loc == 0.0_wp) then
          kf(ireac) = ko_loc
        else
          exp_zfc = 1.0_wp / (1.0_wp + log10(ko_loc / kinf_loc)**2)
          if (sri_param(ireac, 3) == 0.0_wp) then
            zfc = sri_param(ireac, 4) * (sri_param(ireac, 1) * exp(-sri_param(ireac, 2) * one_over_t))**exp_zfc * &
                  tloc**sri_param(ireac, 5)
          else
            zfc = sri_param(ireac, 4) * (sri_param(ireac, 1) * exp(-sri_param(ireac, 2) * one_over_t) + &
                  exp(-tloc / sri_param(ireac, 3)))**exp_zfc * tloc**sri_param(ireac, 5)
          end if
          kf(ireac) = (ko_loc / (1.0_wp + ko_loc / kinf_loc)) * zfc
        end if
      end if

      ! Calculate reverse reaction rate coefficients (kr)
      if (mod_param(ireac, 2) == -1) then
        k_loc = 0.0_wp
        exp_pt = 0
        do i = 1, nb_prod(ireac)
          iloc = imol_prod(ireac, i)
          ipr = 1
          exp_pt = exp_pt + 1
          k_loc = k_loc + ipr * mu(iloc)
        end do
        do i = 1, nb_reac(ireac)
          iloc = imol_reac(ireac, i)
          ipr = -1
          exp_pt = exp_pt - 1
          k_loc = k_loc + ipr * mu(iloc)
        end do

        ! Calculate equilibrium constant
        k_loc16 = k_loc
        k_loc16 = exp(-k_loc16) * nref**exp_pt

        ! Calculate reverse rate (kr)
        kr(ireac) = real(kf(ireac) / k_loc16, kind=wp)
      end if

    end do ! Loop over reactions

  end subroutine get_krkf


  subroutine get_thermodynamics(tloc, mu, h, cp)

    implicit none

    ! Input/output arguments
    real(wp), intent(in)  :: tloc
    real(wp), intent(out) :: mu(nmol)
    real(wp), intent(out) :: h(nmol)
    real(wp), intent(out) :: cp(nmol)

    ! Local variables
    integer :: i, istart

    ! Calculate chemical potential, enthalpy, and specific heat for all species
    ! using NASA polynomial coefficients (McBride, Gordon, and Reno, 1993)
    ! See Venot et al. (2012) for details

    do i = 1, nmol
      ! Determine the start index based on temperature
      if (tloc > mol_param(3, i)) then
        istart = 4
      else
        istart = 13
      end if

      ! Chemical potential (mu)
      mu(i) =  -(mol_param(istart, i) / 2.0_wp) * tloc**(-2.0_wp) + &
               (mol_param(istart+1, i) / tloc) * (1.0_wp + log(tloc)) + &
               mol_param(istart+2, i) * (1.0_wp - log(tloc)) - &
               mol_param(istart+3, i) * tloc / 2.0_wp - &
               mol_param(istart+4, i) * tloc**2 / 6.0_wp - &
               mol_param(istart+5, i) * tloc**3 / 12.0_wp - &
               mol_param(istart+6, i) * tloc**4 / 20.0_wp + &
               mol_param(istart+7, i) / tloc - &
               mol_param(istart+8, i)

      ! Specific heat (cp)
      cp(i) =   mol_param(istart, i) * tloc**(-2.0_wp) + &
                mol_param(istart+1, i) * tloc**(-1.0_wp) + &
                mol_param(istart+2, i) + &
                mol_param(istart+3, i) * tloc + &
                mol_param(istart+4, i) * tloc**2 + &
                mol_param(istart+5, i) * tloc**3 + &
                mol_param(istart+6, i) * tloc**4

      ! Enthalpy (h)
      h(i)  =  -mol_param(istart, i) * tloc**(-2.0_wp) + &
               (mol_param(istart+1, i) / tloc) * log(tloc) + &
               mol_param(istart+2, i) + &
               mol_param(istart+3, i) * tloc / 2.0_wp + &
               mol_param(istart+4, i) * tloc**2 / 3.0_wp + &
               mol_param(istart+5, i) * tloc**3 / 4.0_wp + &
               mol_param(istart+6, i) * tloc**4 / 5.0_wp + &
               mol_param(istart+7, i) / tloc
    end do

  end subroutine get_thermodynamics

  subroutine dndt(node, time, n, dndt_loc)
  ! ------------------------------------------------------------------------------
  ! Subroutine: dndt
  ! Purpose:
  !   This subroutine calculates the rate of change of species concentrations
  !   (`dndt_loc`) for a given node in a system of chemical reactions. It
  !   handles both forward and reversed reactions and updates the rate of change
  !   of each species based on the reaction parameters, concentrations, and
  !   reaction rates.
  !
  ! Arguments:
  !   - node (input): Integer, index of the node for which the rate of change is calculated.
  !   - time (input): Real, current time (not used in calculations, but may be part of
  !     a larger system).
  !   - n (input): Real array of size `node`, containing the concentrations or quantities
  !     of species involved in the reactions.
  !   - dndt_loc (output): Real array of size `node`, containing the calculated rates
  !     of change of the concentrations for the species.
  !
  ! Local Variables:
  !   - ireac: Integer, loop counter for reactions.
  !   - ip: Integer, loop counter for products.
  !   - ir: Integer, loop counter for reactants.
  !   - prls_loc: Real variable used to store the production/loss term for a reaction.
  !
  ! This subroutine assumes that the reaction rate constants (kf, kr) and the
  ! stoichiometric coefficients (nb_reac, nb_prod, imol_reac, imol_prod) are
  ! already defined and available in the program context.
  ! ------------------------------------------------------------------------------

    implicit none

    ! Input arguments
    integer, intent(in) :: node         ! Index of the node
    real(wp), intent(in) :: time        ! Current time
    real(wp), intent(in) :: n(node)     ! Array of concentrations or quantities (input)

    ! Output argument
    real(wp), intent(out) :: dndt_loc(node)  ! Array of rate of change (output)

    ! Local variables
    integer :: ireac, ip, ir           ! Loop counters for reactions, products, and reactants
    real(wp) :: prls_loc               ! Local variable to store production/loss term

    ! Initialize the rate of change array to zero
    dndt_loc = 0.0_wp

    ! Loop over all reactions
    do ireac = 1, nreac

      ! Check if the reaction is not a reversible reaction
      if (.not. mod_param(ireac, 2) == 0) then
        ! Compute the production/loss term prls for forward reaction
        prls_loc = kf(ireac)

        ! Loop over the reactants and calculate the production/loss term
        do ir = 1, nb_reac(ireac)
          prls_loc = prls_loc * n(imol_reac(ireac, ir))
        end do

        ! Update the derivatives for the product species
        do ip = 1, nb_prod(ireac)
          dndt_loc(imol_prod(ireac, ip)) = dndt_loc(imol_prod(ireac, ip)) + prls_loc
        end do

        ! Update the derivatives for the reactant species (loss terms)
        do ir = 1, nb_reac(ireac)
          dndt_loc(imol_reac(ireac, ir)) = dndt_loc(imol_reac(ireac, ir)) - prls_loc
        end do
      end if

      ! If the reaction is reversed, compute the production/loss term for the reversed reaction
      if (mod_param(ireac, 2) == -1) then
        ! Compute the production/loss term prls for the reversed reaction
        prls_loc = kr(ireac)

        ! Loop over the products and calculate the production/loss term
        do ip = 1, nb_prod(ireac)
          prls_loc = prls_loc * n(imol_prod(ireac, ip))
        end do

        ! Update the derivatives for the reactant species in the reversed reaction
        do ir = 1, nb_reac(ireac)
          dndt_loc(imol_reac(ireac, ir)) = dndt_loc(imol_reac(ireac, ir)) + prls_loc
        end do

        ! Update the derivatives for the product species (loss terms) in the reversed reaction
        do ip = 1, nb_prod(ireac)
          dndt_loc(imol_prod(ireac, ip)) = dndt_loc(imol_prod(ireac, ip)) - prls_loc
        end do
      end if

    end do

  end subroutine dndt


  subroutine jac(ndim_loc, T, n, j, ian, jan, pd)

    implicit none

    integer, intent(in) :: ndim_loc, j
    real(wp), intent(in) :: T
    real(wp), intent(in) :: n(ndim_loc)
    real(wp), intent(in) :: ian(ndim_loc)
    real(wp), intent(in) :: jan(ndim_loc)
    real(wp), intent(out) :: pd(ndim_loc)

    ! Local variables
    integer :: ireac, ip, ir, ir2, ip2, i1
    real(wp) :: prls_loc, ord

    logical :: incl

!    pd = 0.0_wp
!    DO ireac = 1,nreac
!       IF (.not. mod_param(ireac,2) == 0) THEN
!          ! compute the production/loss term prls
!          incl = .false.
!          ord = 0.0_wp
!          DO ir = 1,nb_reac(ireac)
!             if (imol_reac(ireac,ir) .eq. j) THEN
!                incl = .true.
!                ord = ord + 1.0_wp
!             ENDIF
!          ENDDO
!
!          if (incl) THEN
!             prls_loc = kf(ireac)
!             DO ir = 1,nb_reac(ireac)
!                prls_loc = prls_loc*n(imol_reac(ireac,ir))
!             END DO
!             prls_loc = ord*prls_loc/n(j)
!             ! complete the derivatives with the production terms
!             DO ip = 1,nb_prod(ireac)
!                i1 = imol_prod(ireac,ip)
!                pd(i1) = pd(i1) + prls_loc
!             END DO
!             ! complete the derivatives with the loss terms
!             DO ir = 1,nb_reac(ireac)
!                i1 = imol_reac(ireac,ir)
!                pd(i1) = pd(i1) - prls_loc
!             END DO
!          ENDIF
!       END IF
!
!       ! if reversed compute the production/loss term of the reversed reactions
!       IF (mod_param(ireac,2) == -1) THEN
!
!          incl = .false.
!          ord = 0.0_wp
!          DO ip = 1,nb_prod(ireac)
!             if (imol_prod(ireac,ip) .eq. j) THEN
!                incl = .true.
!                ord = ord + 1.0_wp
!             ENDIF
!          ENDDO
!
!          if (incl) then
!             !compute the production/loss term prls
!             prls_loc = kr(ireac)
!             DO ip = 1,nb_prod(ireac)
!                prls_loc = prls_loc*n(imol_prod(ireac,ip))
!             END DO
!             prls_loc = ord*prls_loc/n(j)
!             ! complete the derivatives with the production terms
!             DO ir = 1,nb_reac(ireac)
!                i1 = imol_reac(ireac,ir)
!                pd(i1) = pd(i1) + prls_loc
!             END DO
!             ! complete the derivatives with the loss terms
!             DO ip = 1,nb_prod(ireac)
!                i1 = imol_prod(ireac,ip)
!                pd(i1) = pd(i1) - prls_loc
!             END DO
!          endif
!       END IF
!    END DO

  end subroutine

end module flexchem_kinetics_mod
