! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Centre National de la Recherche Scientifique
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2007: ActSim code of M. Chiriaco and H. Chepfer rewritten by S. Bony
!
! May 2008, H. Chepfer:
! - Units of pressure inputs: Pa 
! - Non Spherical particles : LS Ice NS coefficients, CONV Ice NS coefficients
! - New input: ice_type (0=ice-spheres ; 1=ice-non-spherical)
!
! June 2008, A. Bodas-Salcedo:
! - Ported to Fortran 90 and optimisation changes
!
! August 2008, J-L Dufresne:
! - Optimisation changes (sum instructions suppressed)
!
! October 2008, S. Bony,  H. Chepfer and J-L. Dufresne :  
! - Interface with COSP v2.0:
!      cloud fraction removed from inputs
!      in-cloud condensed water now in input (instead of grid-averaged value)
!      depolarisation diagnostic removed
!      parasol (polder) reflectances (for 5 different solar zenith angles) added
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - Modification of the integration of the lidar equation.
! - change the cloud detection threshold
!
! April 2008, A. Bodas-Salcedo:
! - Bug fix in computation of pmol and pnorm of upper layer
!
! April 2008, J-L. Dufresne
! - Bug fix in computation of pmol and pnorm, thanks to Masaki Satoh: a factor 2 
! was missing. This affects the ATB values but not the cloud fraction. 
!
! January 2013, G. Cesana and H. Chepfer:
! - Add the perpendicular component of the backscattered signal (pnorm_perp_tot) in the arguments
! - Add the temperature for each levels (temp) in the arguments
! - Add the computation of the perpendicular component of the backscattered lidar signal 
! Reference: Cesana G. and H. Chepfer (2013): Evaluation of the cloud water phase
! in a climate model using CALIPSO-GOCCP, J. Geophys. Res., doi: 10.1002/jgrd.50376
!
! May 2015 - D. Swales - Modified for COSPv2.0
!
! Mar 2018 - R. Guzman - Added OPAQ subroutines
! References OPAQ:
!
!       Guzman et al. (2017): Direct atmosphere opacity observations from CALIPSO provide 
! new constraints on cloud-radiation interactions. JGR-Atmospheres, DOI: 10.1002/2016JD025946
!       Vaillant de Guelis et al. (2017a): The link between outgoing longwave radiation and
! the altitude at which a spaceborne lidar beam is fully attenuated. AMT, 10, 4659-4685,
! https://doi.org/10.5194/amt-10-4659-2017
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_lidar_simulator
  USE COSP_KINDS,         ONLY: wp
  USE MOD_COSP_CONFIG,    ONLY: SR_BINS,S_CLD,S_ATT,S_CLD_ATT,R_UNDEF,calipso_histBsct,  &
                                use_vgrid,vgrid_zl,vgrid_zu,vgrid_z,atlid_histBsct,      &
                                grLidar532_histBsct,S_CLD_ATLID,S_ATT_ATLID,S_CLD_ATT_ATLID
  USE MOD_COSP_STATS,     ONLY: COSP_CHANGE_VERTICAL_GRID,hist1d
  implicit none
  
  ! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the 
  ! ATBperpendicular as a function of the ATB for ice or liquid cloud particles 
  ! derived from CALIPSO-GOCCP observations at 120m vertical grid 
  ! (Cesana and Chepfer, JGR, 2013).
  !
  ! Relationship between ATBice and ATBperp,ice for ice particles:
  !                ATBperp,ice = Alpha*ATBice 
  ! Relationship between ATBice and ATBperp,ice for liquid particles:
  !          ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
  real(wp) :: &
       alpha,beta,gamma    

contains
  ! ######################################################################################
  ! SUBROUTINE lidar_subcolumn
  ! Inputs with a vertical dimensions (nlev) should ordered in along the vertical 
  ! dimension from TOA-2-SFC, for example: varIN(nlev) is varIN @ SFC. 
  ! ######################################################################################
  subroutine lidar_subcolumn(npoints, ncolumns, nlev, lground, beta_mol, tau_mol,        &
       betatot, tautot, pmol, pnorm, betatot_ice, tautot_ice, betatot_liq, tautot_liq,   &
       pnorm_perp_tot)

    ! INPUTS
    INTEGER,intent(in) :: & 
         npoints,      & ! Number of gridpoints
         ncolumns,     & ! Number of subcolumns
         nlev            ! Number of levels
    logical,intent(in) :: &
         lground         ! True for ground-based lidar simulator
    REAL(WP),intent(in),dimension(npoints,nlev) :: &
         beta_mol,     & ! Molecular backscatter coefficient
         tau_mol         ! Molecular optical depth
    REAL(WP),intent(in),dimension(npoints,ncolumns,nlev)       :: &
         betatot,      & ! 
         tautot          ! Optical thickess integrated from top
    ! Optional Inputs
    REAL(WP),intent(in),dimension(npoints,ncolumns,nlev),optional       :: &         
         betatot_ice,  & ! Backscatter coefficient for ice particles
         betatot_liq,  & ! Backscatter coefficient for liquid particles
         tautot_ice,   & ! Total optical thickness of ice
         tautot_liq      ! Total optical thickness of liq

    ! OUTPUTS
    REAL(WP),intent(out),dimension(npoints,nlev) :: &
         pmol            ! Molecular attenuated backscatter lidar signal power(m^-1.sr^-1)
    REAL(WP),intent(out),dimension(npoints,ncolumns,nlev) :: &
         pnorm           ! Molecular backscatter signal power (m^-1.sr^-1)
    ! Optional outputs
    REAL(WP),intent(out),dimension(npoints,ncolumns,nlev),optional :: &
         pnorm_perp_tot  ! Perpendicular lidar backscattered signal power

    ! LOCAL VARIABLES
    INTEGER :: k,icol,zi,zf,zinc
    logical :: lphaseoptics
    REAL(WP),dimension(npoints) :: &
         tautot_lay        !
    REAL(WP),dimension(npoints,ncolumns,nlev) :: &
         pnorm_liq,      & ! Lidar backscattered signal power for liquid
         pnorm_ice,      & ! Lidar backscattered signal power for ice
         pnorm_perp_ice, & ! Perpendicular lidar backscattered signal power for ice
         pnorm_perp_liq, & ! Perpendicular lidar backscattered signal power for liq
         beta_perp_ice,  & ! Perpendicular backscatter coefficient for ice
         beta_perp_liq     ! Perpendicular backscatter coefficient for liquid
    REAL(WP),dimension(npoints,nlev) :: &
         tautot_ice_tmp,     &
         tautot_liq_tmp,     &
         pnorm_perp_ice_tmp, &
         pnorm_perp_liq_tmp, &
         beta_perp_ice_tmp,  &
         beta_perp_liq_tmp,  &
         beta_mol_tmp,       &
         tau_mol_tmp,        &
         betatot_tmp,        &
         tautot_tmp,         &
         betatot_ice_tmp,    &
         betatot_liq_tmp,    &
         pnorm_tmp,          &
         pmol_tmp,           &
         pnorm_ice_tmp,      &
         pnorm_liq_tmp

    ! Phase optics?
    lphaseoptics=.false.
    if (present(betatot_ice) .and. present(betatot_liq) .and. present(tautot_liq) .and. &
         present(tautot_ice)) lphaseoptics=.true.
    
    ! Is this lidar spaceborne (default) or ground-based?
    if (lground) then
       zi   = nlev
       zf   = 1
       zinc = -1
    else
       zi   = 1
       zf   = nlev
       zinc = 1
    endif    

    ! ####################################################################################
    ! *) Molecular signal
    ! ####################################################################################
    beta_mol_tmp(:,:) = beta_mol(1:npoints,zi:zf:zinc)
    tau_mol_tmp(:,:) = tau_mol(1:npoints,zi:zf:zinc)
    call cmp_backsignal(nlev,npoints,beta_mol_tmp,tau_mol_tmp,pmol_tmp)
    pmol(1:npoints,zi:zf:zinc) = pmol_tmp(:,:)
    ! ####################################################################################
    ! PLANE PARALLEL FIELDS
    ! ####################################################################################
    do icol=1,ncolumns
       ! #################################################################################
       ! *) Total Backscatter signal
       ! #################################################################################
       betatot_tmp(:,:) = betatot(1:npoints,icol,zi:zf:zinc)
       tautot_tmp(:,:) = tautot(1:npoints,icol,zi:zf:zinc)
       call cmp_backsignal(nlev,npoints,betatot_tmp,tautot_tmp,pnorm_tmp)
       pnorm(1:npoints,icol,zi:zf:zinc) = pnorm_tmp(:,:)
       ! #################################################################################
       ! *) Ice/Liq Backscatter signal
       ! #################################################################################
       if (lphaseoptics) then
          ! Computation of the ice and liquid lidar backscattered signal (ATBice and ATBliq)
          ! Ice only
          betatot_ice_tmp(:,:) = betatot_ice(1:npoints,icol,zi:zf:zinc)
          tautot_ice_tmp(:,:) = tautot_ice(1:npoints,icol,zi:zf:zinc)
          call cmp_backsignal(nlev,npoints,betatot_ice_tmp,tautot_ice_tmp, pnorm_ice_tmp)
          pnorm_ice(1:npoints,icol,zi:zf:zinc) = pnorm_ice_tmp(:,:)

          ! Liquid only
          betatot_liq_tmp(:,:) = betatot_liq(1:npoints,icol,zi:zf:zinc)
          tautot_liq_tmp(:,:) = tautot_liq(1:npoints,icol,zi:zf:zinc)
          call cmp_backsignal(nlev,npoints,betatot_liq_tmp,tautot_liq_tmp, pnorm_liq_tmp)
          pnorm_liq(1:npoints,icol,zi:zf:zinc) = pnorm_liq_tmp(:,:)
       endif
    enddo

    ! ####################################################################################
    ! PERDENDICULAR FIELDS (Only needed if distinguishing by phase (ice/liquid))
    ! ####################################################################################
    if (lphaseoptics) then
       do icol=1,ncolumns
          ! #################################################################################
          ! *) Ice/Liq Perpendicular Backscatter signal
          ! #################################################################################
          ! Computation of ATBperp,ice/liq from ATBice/liq including the multiple scattering 
          ! contribution (Cesana and Chepfer 2013, JGR)
          do k=1,nlev
             ! Ice particles
             pnorm_perp_ice(1:npoints,icol,k) = Alpha * pnorm_ice(1:npoints,icol,k)
             
             ! Liquid particles
             pnorm_perp_liq(1:npoints,icol,k) = 1000._wp*Beta*pnorm_liq(1:npoints,icol,k)**2+&
                  Gamma*pnorm_liq(1:npoints,icol,k) 
          enddo
          
          ! #################################################################################
          ! *) Computation of beta_perp_ice/liq using the lidar equation
          ! #################################################################################
          ! Ice only
          pnorm_perp_ice_tmp(:,:) = pnorm_perp_ice(1:npoints,icol,zi:zf:zinc)
          tautot_ice_tmp(:,:) = tautot_ice(1:npoints,icol,zi:zf:zinc)
          call cmp_beta(nlev,npoints,pnorm_perp_ice_tmp, tautot_ice_tmp,beta_perp_ice_tmp)
          beta_perp_ice(1:npoints,icol,zi:zf:zinc) = beta_perp_ice_tmp(:,:)

          ! Liquid only
          pnorm_perp_liq_tmp(:,:) = pnorm_perp_liq(1:npoints,icol,zi:zf:zinc)
          tautot_liq_tmp(:,:) = tautot_liq(1:npoints,icol,zi:zf:zinc)
          call cmp_beta(nlev,npoints,pnorm_perp_liq_tmp,tautot_liq_tmp,beta_perp_liq_tmp)
          beta_perp_liq(1:npoints,icol,zi:zf:zinc) = beta_perp_liq_tmp(:,:)
          ! #################################################################################
          ! *) Perpendicular Backscatter signal
          ! #################################################################################
          ! Computation of the total perpendicular lidar signal (ATBperp for liq+ice)
          ! Upper layer
          WHERE(tautot(1:npoints,icol,1) .gt. 0)
             pnorm_perp_tot(1:npoints,icol,1) = (beta_perp_ice(1:npoints,icol,1)+           &
                  beta_perp_liq(1:npoints,icol,1)-                                          &
                  (beta_mol(1:npoints,1)/(1._wp+1._wp/0.0284_wp))) /                        &
                  (2._wp*tautot(1:npoints,icol,1))*                                         &
                  (1._wp-exp(-2._wp*tautot(1:npoints,icol,1)))
          ELSEWHERE
             pnorm_perp_tot(1:npoints,icol,1) = 0._wp
          ENDWHERE
          
          ! Other layers
          do k=2,nlev
             ! Optical thickness of layer k
             tautot_lay(1:npoints) = tautot(1:npoints,icol,k)-tautot(1:npoints,icol,k-1) 
             
             ! The perpendicular component of the molecular backscattered signal (Betaperp) 
             ! has been taken into account two times (once for liquid and once for ice). 
             ! We remove one contribution using 
             ! Betaperp=beta_mol(:,k)/(1+1/0.0284)) [bodhaine et al. 1999] in the following 
             ! equations:
             WHERE (pnorm(1:npoints,icol,k) .eq. 0)
                pnorm_perp_tot(1:npoints,icol,k)=0._wp
             ELSEWHERE
                where(tautot_lay(1:npoints) .gt. 0.)
                   pnorm_perp_tot(1:npoints,icol,k) = (beta_perp_ice(1:npoints,icol,k)+     &
                        beta_perp_liq(1:npoints,icol,k)-(beta_mol(1:npoints,k)/(1._wp+1._wp/  &
                        0.0284_wp)))*EXP(-2._wp*tautot(1:npoints,icol,k-1))/                  &
                        (2._wp*tautot_lay(1:npoints))* (1._wp-EXP(-2._wp*tautot_lay(1:npoints)))
                elsewhere
                   pnorm_perp_tot(1:npoints,icol,k) = (beta_perp_ice(1:npoints,icol,k)+     &
                        beta_perp_liq(1:npoints,icol,k)-(beta_mol(1:npoints,k)/(1._wp+1._wp/  &
                        0.0284_wp)))*EXP(-2._wp*tautot(1:npoints,icol,k-1))
                endwhere
             ENDWHERE
          END DO
       enddo
    end if
  end subroutine lidar_subcolumn

  ! ######################################################################################
  ! SUBROUTINE lidar_column
  ! ######################################################################################
  subroutine lidar_column(npoints, ncol, nlevels, llm, max_bin, ntype, platform, pnorm, pmol,             &
       pplay, zlev, zlev_half, vgrid_z, ok_lidar_cfad, ncat, cfad2, lidarcld, cldlayer,  &
       ! Optional stuff below
       tmp, pnorm_perp, surfelev,  lidarcldphase, lidarcldtype, cldtype, cldtypetemp, &
       cldtypemeanz, cldtypemeanzse, cldthinemis, cldlayerphase, lidarcldtmp)

    integer,parameter :: &
         nphase = 6 ! Number of cloud layer phase types

    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal grid points
         ncol,    & ! Number of subcolumns
         nlevels, & ! Number of vertical layers (OLD grid)
         llm,     & ! Number of vertical layers (NEW grid)
         max_bin, & ! Number of bins for SR CFADs
         ncat,    & ! Number of cloud layer types (low,mid,high,total) 
         ntype      ! Number of OPAQ products (opaque/thin cloud + z_opaque)
    character(len=*),intent(in) :: &
         platform   ! Name of platform (e.g. calipso,atlid,grLidar532)
    real(wp),intent(in),dimension(npoints,ncol,Nlevels) :: &
         pnorm      ! Lidar ATB
    real(wp),intent(in),dimension(npoints,Nlevels) :: &
         pmol,    & ! Molecular ATB
         pplay      ! Pressure on model levels (Pa)
    logical,intent(in) :: &
         ok_lidar_cfad ! True if lidar CFAD diagnostics need to be computed
    real(wp),intent(in),dimension(npoints,nlevels) :: &
         zlev        ! Model full levels
    real(wp),intent(in),dimension(npoints,nlevels) :: &
         zlev_half   ! Model half levels
    real(wp),intent(in),dimension(llm) :: & 
         vgrid_z     ! mid-level altitude of the output vertical grid
    ! Optional Inputs
    real(wp),intent(in),dimension(npoints,ncol,Nlevels),optional :: &
         pnorm_perp ! Lidar perpendicular ATB
    real(wp),intent(in),dimension(npoints),optional :: &
         surfelev   ! Surface Elevation (m)
    real(wp),intent(in),dimension(npoints,Nlevels),optional :: &
         tmp        ! Temperature at each levels
    
    ! Outputs
    real(wp),intent(inout),dimension(npoints,llm) :: &
         lidarcld      ! 3D "lidar" cloud fraction
    real(wp),intent(inout),dimension(npoints,ncat) :: &
         cldlayer      ! "lidar" cloud layer fraction (low, mid, high, total)
    real(wp),intent(inout),dimension(npoints,max_bin,llm) :: &
         cfad2         ! CFADs of SR
    ! Optional Outputs
    real(wp),intent(out),dimension(npoints,ntype),optional :: & 
         cldtype,    & ! "lidar" OPAQ type covers (opaque/thin cloud + z_opaque)
         cldtypetemp   ! Opaque and thin clouds + z_opaque temperature
    real(wp),intent(out),dimension(npoints,2),optional :: &  
         cldtypemeanz  ! Opaque and thin clouds altitude 
    real(wp),intent(out),dimension(npoints,3),optional :: &  
         cldtypemeanzse ! Opaque, thin clouds and z_opaque altitude with respect to SE
    real(wp),intent(out),dimension(npoints),optional :: &   
         cldthinemis   ! Thin clouds emissivity computed from SR
    real(wp),intent(out),dimension(npoints,llm,nphase),optional :: &
         lidarcldphase ! 3D "lidar" phase cloud fraction
    real(wp),intent(out),dimension(npoints,llm,ntype+1),optional :: & 
         lidarcldtype ! 3D "lidar" OPAQ type fraction 
    real(wp),intent(out),dimension(npoints,40,5),optional :: &
         lidarcldtmp   ! 3D "lidar" phase cloud fraction as a function of temp
    real(wp),intent(out),dimension(npoints,ncat,nphase),optional :: &
         cldlayerphase ! "lidar" phase low mid high cloud fraction

    ! Local Variables
    integer :: ic,i,j
    logical :: lcalipso,latlid,lgrlidar532
    real(wp),dimension(npoints,ncol,llm) :: &
         x3d
    real(wp),dimension(npoints,llm) :: &
         x3d_c,pnorm_c
    real(wp)  :: &
         xmax
    real(wp),dimension(npoints,1,Nlevels) :: t_in,ph_in,betamol_in
    real(wp),dimension(npoints,ncol,llm)  :: pnormFlip,pnorm_perpFlip,pnorm_perpFlip_tmp
    real(wp),dimension(npoints,1,llm)     :: tmpFlip,pplayFlip,betamolFlip
    real(wp),dimension(SR_BINS+1)         :: histBsct
    real(wp),dimension(npoints,Nlevels)   :: zlev_tmp,zlev_half_tmp
    real(wp),dimension(llm)               :: vgrid_zl_tmp,vgrid_zu_tmp
    real(wp),dimension(npoints,1,llm)     :: pplayFlip_tmp,betamolFlip_tmp,tmpFlip_tmp
    real(wp),dimension(npoints,ncol,llm)  :: pnormFlip_tmp
    real(wp),dimension(npoints,ncol,Nlevels) :: pnorm_in, pnorm_perp_in
    real(wp),dimension(ncol) :: x3d_tmp


    ! Which lidar platform?
    lcalipso = .false.
    latlid = .false.
    lgrlidar532 = .false.
    if (platform .eq. 'calipso') lcalipso=.true.
    if (platform .eq. 'atlid') latlid=.true.
    if (platform .eq. 'grlidar532') lgrlidar532=.true.
        
    ! Vertically regrid input data
    if (use_vgrid) then 
       ph_in(:,1,:) = pplay(:,nlevels:1:-1)
       zlev_tmp(:,:) = zlev(:,nlevels:1:-1)
       zlev_half_tmp(:,:) = zlev_half(:,nlevels:1:-1)
       vgrid_zl_tmp(:) = vgrid_zl(llm:1:-1)
       vgrid_zu_tmp(:) = vgrid_zu(llm:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev_tmp,zlev_half_tmp,&
            ph_in,llm,vgrid_zl_tmp,vgrid_zu_tmp,pplayFlip_tmp)
       pplayFlip(:,:,llm:1:-1) = pplayFlip_tmp(:,:,:)

       betamol_in(:,1,:) = pmol(:,nlevels:1:-1)
       pnorm_in(:,:,:) = pnorm(:,:,nlevels:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev_tmp,zlev_half_tmp,&
            betamol_in,llm,vgrid_zl_tmp,vgrid_zu_tmp,betamolFlip_tmp)
       call cosp_change_vertical_grid(Npoints,Ncol,Nlevels,zlev_tmp,zlev_half_tmp,&
            pnorm_in,llm,vgrid_zl_tmp,vgrid_zu_tmp,pnormFlip_tmp)
       betamolFlip(:,:,llm:1:-1) = betamolFlip_tmp(:,:,:)
       pnormFlip(:,:,llm:1:-1) = pnormFlip_tmp(:,:,:)

       if (lcalipso) then
          t_in(:,1,:) = tmp(:,nlevels:1:-1)
          pnorm_perp_in = pnorm_perp(:,:,nlevels:1:-1)
          call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev_tmp,zlev_half_tmp,&
               t_in,llm,vgrid_zl_tmp,vgrid_zu_tmp,tmpFlip_tmp)
          call cosp_change_vertical_grid(Npoints,Ncol,Nlevels,zlev_tmp,zlev_half_tmp,&
               pnorm_perp_in,llm,vgrid_zl_tmp,vgrid_zu_tmp,pnorm_perpFlip_tmp)
          tmpFlip(:,:,llm:1:-1) = tmpFlip_tmp(:,:,:)
          pnorm_perpFlip(:,:,llm:1:-1) = pnorm_perpFlip_tmp(:,:,:)
       endif
    endif

    ! Initialization (The histogram bins, are set up during initialization and the
    ! maximum value is used as the upper bounds.)
    if (lcalipso)    then
       xmax = maxval(calipso_histBsct)
       histBsct = calipso_histBsct
    endif
    if (latlid)      then
       xmax = maxval(atlid_histBsct)
       histBsct = atlid_histBsct
    endif
    if (lgrlidar532) then
       xmax = maxval(grLidar532_histBsct)
       histBsct = grLidar532_histBsct
    endif
       
    ! Compute LIDAR scattering ratio
    if (use_vgrid) then
       do ic = 1, ncol
          pnorm_c = pnormFlip(:,ic,:)
          where ((pnorm_c .lt. xmax) .and. (betamolFlip(:,1,:) .lt. xmax) .and.          &
                (betamolFlip(:,1,:) .gt. 0.0 ))
             x3d_c = pnorm_c/betamolFlip(:,1,:)
          elsewhere
             x3d_c = R_UNDEF
          end where
          x3d(:,ic,:) = x3d_c
       enddo
       if (lcalipso) then
          ! Diagnose cloud fractions for subcolumn lidar scattering ratios
          CALL COSP_CLDFRAC(npoints,ncol,llm,ncat,nphase,tmpFlip,x3d,pnormFlip,pnorm_perpFlip,&
               pplayFlip,S_att,S_cld,S_cld_att,R_UNDEF,lidarcld,cldlayer,lidarcldphase,&
               cldlayerphase,lidarcldtmp)                         

          ! Calipso opaque cloud diagnostics
          CALL COSP_OPAQ(npoints,ncol,llm,ntype,tmpFlip,x3d,S_att,S_cld,R_UNDEF,lidarcldtype, &
               cldtype,cldtypetemp,cldtypemeanz,cldtypemeanzse,cldthinemis,vgrid_z,surfelev)
       endif
       if (latlid) then
          CALL COSP_CLDFRAC_NOPHASE(npoints,ncol,llm,ncat,x3d,pnormFlip,pplayFlip,  &
               S_att_atlid,S_cld_atlid,S_cld_att_atlid,R_UNDEF,lidarcld,cldlayer)   
       endif
       if (lgrLidar532) then
          CALL COSP_CLDFRAC_NOPHASE(npoints,ncol,llm,ncat,x3d,pnormFlip,pplayFlip,  &
               S_att,S_cld,S_cld_att,R_UNDEF,lidarcld,cldlayer)
       endif
    else
       do ic = 1, ncol
          pnorm_c = pnorm(:,ic,:)
          where ((pnorm_c.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0 ))
             x3d_c = pnorm_c/pmol
          elsewhere
             x3d_c = R_UNDEF
          end where
          x3d(:,ic,:) = x3d_c
       enddo
       if (lcalipso) then
          ! Diagnose cloud fractions for subcolumn lidar scattering ratios
          CALL COSP_CLDFRAC(npoints,ncol,nlevels,ncat,nphase,tmp,x3d,pnorm,pnorm_perp,pplay,&
               S_att,S_cld,S_cld_att,R_UNDEF,lidarcld,cldlayer,lidarcldphase,  &
               cldlayerphase,lidarcldtmp)
          ! Calipso opaque cloud diagnostics
          CALL COSP_OPAQ(npoints,ncol,nlevels,ntype,tmp,x3d,S_att,S_cld,R_UNDEF,lidarcldtype, &
               cldtype,cldtypetemp,cldtypemeanz,cldtypemeanzse,cldthinemis,vgrid_z,surfelev)
       endif
       if (latlid) then
          CALL COSP_CLDFRAC_NOPHASE(npoints,ncol,nlevels,ncat,x3d,pnorm,pplay,  &
               S_att_atlid,S_cld_atlid,S_cld_att_atlid, R_UNDEF,lidarcld,cldlayer)
       endif
       if (lgrlidar532) then
          CALL COSP_CLDFRAC_NOPHASE(npoints,ncol,nlevels,ncat,x3d,pnorm,pplay,      &
               S_att,S_cld,S_cld_att,R_UNDEF,lidarcld,cldlayer)
       endif
    endif

    ! CFADs
    if (ok_lidar_cfad) then
       ! CFADs of subgrid-scale lidar scattering ratios
       do i=1,Npoints
          do j=1,llm
             x3d_tmp(:) = x3d(i,:,j)
             cfad2(i,:,j) = hist1D(ncol,x3d_tmp,SR_BINS,histBsct)
          enddo
       enddo
       where(cfad2 .ne. R_UNDEF) cfad2=cfad2/ncol
    endif 
    
    ! Unit conversions
    where(lidarcld /= R_UNDEF)      lidarcld      = lidarcld*100._wp
    where(cldlayer /= R_UNDEF)      cldlayer      = cldlayer*100._wp
    if (lcalipso) then
       where(cldtype(:,1) /= R_UNDEF)  cldtype(:,1)  = cldtype(:,1)*100._wp
       where(cldtype(:,2) /= R_UNDEF)  cldtype(:,2)  = cldtype(:,2)*100._wp
       where(cldlayerphase /= R_UNDEF) cldlayerphase = cldlayerphase*100._wp
       where(lidarcldphase /= R_UNDEF) lidarcldphase = lidarcldphase*100._wp
       where(lidarcldtype /= R_UNDEF)  lidarcldtype  = lidarcldtype*100._wp
       where(lidarcldtmp /= R_UNDEF)   lidarcldtmp   = lidarcldtmp*100._wp
    endif
  end subroutine lidar_column

  ! ######################################################################################
  ! The subroutines below compute the attenuated backscatter signal and the lidar 
  ! backscatter coefficients using eq (1) from doi:0094-8276/08/2008GL034207
  ! ######################################################################################
  subroutine cmp_backsignal(nlev,npoints,beta,tau,pnorm)
    ! INPUTS
    integer, intent(in) :: nlev,npoints
    real(wp),intent(in),dimension(npoints,nlev) :: beta,tau

    ! OUTPUTS
    real(wp),intent(out),dimension(npoints,nlev) :: pnorm

    ! Internal Variables
    real(wp), dimension(npoints) :: tautot_lay
    integer :: k

    ! Uppermost layer 
    pnorm(:,1) = beta(:,1) / (2._wp*tau(:,1)) * (1._wp-exp(-2._wp*tau(:,1)))

    ! Other layers
    do k=2,nlev
       tautot_lay(:) = tau(:,k)-tau(:,k-1) 
       WHERE (tautot_lay(:) .gt. 0.)
          pnorm(:,k) = beta(:,k)*EXP(-2._wp*tau(:,k-1)) /&
               (2._wp*tautot_lay(:))*(1._wp-EXP(-2._wp*tautot_lay(:)))
       ELSEWHERE
          ! This must never happen, but just in case, to avoid div. by 0
          pnorm(:,k) = beta(:,k) * EXP(-2._wp*tau(:,k-1))
       END WHERE

    END DO
  end subroutine cmp_backsignal

  subroutine cmp_beta(nlev,npoints,pnorm,tau,beta)
    ! INPUTS
    integer, intent(in) :: nlev,npoints
    real(wp),intent(in),dimension(npoints,nlev) :: pnorm,tau

    ! OUTPUTS
    real(wp),intent(out),dimension(npoints,nlev) :: beta

    ! Internal Variables
    real(wp), dimension(npoints) :: tautot_lay
    integer :: k
    real(wp) :: epsrealwp

    epsrealwp = epsilon(1._wp)
    beta(:,1) = pnorm(:,1) * (2._wp*tau(:,1))/(1._wp-exp(-2._wp*tau(:,1)))
    do k=2,nlev
       tautot_lay(:) = tau(:,k)-tau(:,k-1)       
       WHERE ( EXP(-2._wp*tau(:,k-1)) .gt. epsrealwp )
          WHERE (tautot_lay(:) .gt. 0.)
             beta(:,k) = pnorm(:,k)/ EXP(-2._wp*tau(:,k-1))* &
                  (2._wp*tautot_lay(:))/(1._wp-exp(-2._wp*tautot_lay(:)))
          ELSEWHERE
             beta(:,k)=pnorm(:,k)/EXP(-2._wp*tau(:,k-1))
          END WHERE
       ELSEWHERE
          beta(:,k)=pnorm(:,k)/epsrealwp
       END WHERE
    ENDDO

  end subroutine cmp_beta
    ! ####################################################################################
    ! SUBROUTINE cosp_cldfrac
    ! Conventions: Ncat must be equal to 4
    ! ####################################################################################
    SUBROUTINE COSP_CLDFRAC(Npoints,Ncolumns,Nlevels,Ncat,Nphase,tmp,x,ATB,ATBperp,      &
                               pplay,S_att,S_cld,S_cld_att,undef,lidarcld,cldlayer,      &
                               lidarcldphase,cldlayerphase,lidarcldtemp)
    ! Parameters
    integer,parameter :: Ntemp=40 ! indice of the temperature vector
    real(wp),parameter,dimension(Ntemp+1) :: &
       tempmod = [0.0,   183.15,186.15,189.15,192.15,195.15,198.15,201.15,204.15,207.15, &
                  210.15,213.15,216.15,219.15,222.15,225.15,228.15,231.15,234.15,237.15, &
                  240.15,243.15,246.15,249.15,252.15,255.15,258.15,261.15,264.15,267.15, &
                  270.15,273.15,276.15,279.15,282.15,285.15,288.15,291.15,294.15,297.15, &
                  473.15]
         
    ! Polynomial coefficient of the phase discrimination line used to separate liquid from ice
    ! (Cesana and Chepfer, JGR, 2013)
    ! ATBperp = ATB^5*alpha50 + ATB^4*beta50 + ATB^3*gamma50 + ATB^2*delta50 + ATB*epsilon50 + zeta50
    real(wp),parameter :: &
       alpha50   = 9.0322e+15_wp,  & !
       beta50    = -2.1358e+12_wp, & !
       gamma50   = 173.3963e06_wp, & !
       delta50   = -3.9514e03_wp,  & !
       epsilon50 = 0.2559_wp,      & !
       zeta50    = -9.4776e-07_wp    ! 
       
    ! Inputs
    integer,intent(in) :: &
       Npoints,  & ! Number of gridpoints
       Ncolumns, & ! Number of subcolumns
       Nlevels,  & ! Number of vertical levels
       Ncat,     & ! Number of cloud layer types
       Nphase      ! Number of cloud layer phase types
                   ! [ice,liquid,undefined,false ice,false liquid,Percent of ice]
    real(wp),intent(in) :: &
       S_att,    & !
       S_cld,    & !
       S_cld_att,& ! New threshold for undefine cloud phase detection
       undef       ! Undefined value
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: &
       x,        & ! 
       ATB,      & ! 3D attenuated backscatter
       ATBperp     ! 3D attenuated backscatter (perpendicular)
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
       tmp,      & ! Temperature   
       pplay       ! Pressure

    ! Outputs
    real(wp),intent(out),dimension(Npoints,Ntemp,5) :: &
       lidarcldtemp  ! 3D Temperature 1=tot,2=ice,3=liq,4=undef,5=ice/ice+liq
    real(wp),intent(out),dimension(Npoints,Nlevels,Nphase) :: &
       lidarcldphase ! 3D cloud phase fraction
    real(wp),intent(out),dimension(Npoints,Nlevels) :: &
       lidarcld      ! 3D cloud fraction
    real(wp),intent(out),dimension(Npoints,Ncat) :: &
       cldlayer      ! Low, middle, high, total cloud fractions
    real(wp),intent(out),dimension(Npoints,Ncat,Nphase) :: &
       cldlayerphase ! Low, middle, high, total cloud fractions for ice liquid and undefine phase    
    
    ! Local variables
    integer  :: &
       ip, k, iz, ic, ncol, nlev, i, itemp, toplvlsat 
    real(wp) :: &
       p1,checktemp, ATBperp_tmp,checkcldlayerphase, checkcldlayerphase2
    real(wp),dimension(Npoints,Nlevels) :: &
       nsub,lidarcldphasetmp   
    real(wp),dimension(Npoints,Ntemp) :: &
       sumlidarcldtemp,lidarcldtempind
    real(wp),dimension(Npoints,Ncolumns,Ncat) :: &
       cldlay,nsublay   
    real(wp),dimension(Npoints,Ncat) :: &
       nsublayer,cldlayerphasetmp,cldlayerphasesum
    real(wp),dimension(Npoints,Ncolumns,Nlevels) :: &   
       tmpi, & ! Temperature of ice cld
       tmpl, & ! Temperature of liquid cld
       tmpu, & ! Temperature of undef cld
       cldy, & ! 
       srok    !
    real(wp),dimension(Npoints,Ncolumns,Ncat,Nphase) :: &
       cldlayphase ! subgrided low mid high phase cloud fraction
             
    ! ####################################################################################
    ! 1) Initialize    
    ! ####################################################################################
    lidarcld              = 0._wp
    nsub                  = 0._wp
    cldlay                = 0._wp
    nsublay               = 0._wp
    ATBperp_tmp           = 0._wp
    lidarcldphase(:,:,:)  = 0._wp
    cldlayphase(:,:,:,:)  = 0._wp
    cldlayerphase(:,:,:)  = 0._wp
    tmpi(:,:,:)           = 0._wp
    tmpl(:,:,:)           = 0._wp
    tmpu(:,:,:)           = 0._wp
    cldlayerphasesum(:,:) = 0._wp
    lidarcldtemp(:,:,:)   = 0._wp
    lidarcldtempind(:,:)  = 0._wp
    sumlidarcldtemp(:,:)  = 0._wp
    lidarcldphasetmp(:,:) = 0._wp
    toplvlsat             = 0

    ! ####################################################################################
    ! 2) Cloud detection
    ! ####################################################################################
    do k=1,Nlevels
       ! Cloud detection at subgrid-scale:
       where ((x(:,:,k) .gt. S_cld) .and. (x(:,:,k) .ne. undef) )
          cldy(:,:,k)=1._wp
       elsewhere
          cldy(:,:,k)=0._wp
       endwhere
       
       ! Number of usefull sub-columns:
       where ((x(:,:,k) .gt. S_att) .and. (x(:,:,k) .ne. undef) )
          srok(:,:,k)=1._wp
       elsewhere
          srok(:,:,k)=0._wp
       endwhere
    enddo    
    
    ! ####################################################################################
    ! 3) Grid-box 3D cloud fraction and layered cloud fractions(ISCCP pressure categories)
    ! ####################################################################################
    lidarcld = 0._wp
    nsub     = 0._wp
    cldlay   = 0._wp
    nsublay  = 0._wp
    do k=1,Nlevels
       do ic = 1, Ncolumns
          do ip = 1, Npoints
          
             ! Computation of the cloud fraction as a function of the temperature instead
             ! of height, for ice,liquid and all clouds
             if(srok(ip,ic,k).gt.0.)then
                do itemp=1,Ntemp
                   if( (tmp(ip,k).ge.tempmod(itemp)).and.(tmp(ip,k).lt.tempmod(itemp+1)) )then
                      lidarcldtempind(ip,itemp)=lidarcldtempind(ip,itemp)+1._wp
                   endif
                enddo
             endif
             
             if(cldy(ip,ic,k).eq.1.)then
                do itemp=1,Ntemp 
                   if( (tmp(ip,k) .ge. tempmod(itemp)).and.(tmp(ip,k) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(ip,itemp,1)=lidarcldtemp(ip,itemp,1)+1._wp
                   endif
                enddo
             endif

             iz=1
             p1 = pplay(ip,k)
             if ( p1.gt.0. .and. p1.lt.(440._wp*100._wp)) then ! high clouds
                iz=3
             else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid clouds
                iz=2
             endif
             
             cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
             cldlay(ip,ic,4)  = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
             lidarcld(ip,k)   = lidarcld(ip,k) + cldy(ip,ic,k)
             
             nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
             nsublay(ip,ic,4)  = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
             nsub(ip,k)        = nsub(ip,k) + srok(ip,ic,k)
             
          enddo
       enddo
    enddo   
    
    ! Grid-box 3D cloud fraction
    where ( nsub(:,:).gt.0.0 )
       lidarcld(:,:) = lidarcld(:,:)/nsub(:,:)
    elsewhere
       lidarcld(:,:) = undef
    endwhere
    
    ! Layered cloud fractions
    cldlayer  = 0._wp
    nsublayer = 0._wp
    do iz = 1, Ncat
       do ic = 1, Ncolumns
          cldlayer(:,iz)  = cldlayer(:,iz)  + cldlay(:,ic,iz)
          nsublayer(:,iz) = nsublayer(:,iz) + nsublay(:,ic,iz)
       enddo
    enddo
    where (nsublayer(:,:) .gt. 0.0)
       cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
    elsewhere
       cldlayer(:,:) = undef
    endwhere
              
    ! ####################################################################################
    ! 4) Grid-box 3D cloud Phase
    ! ####################################################################################
    
    ! ####################################################################################
    ! 4.1) For Cloudy pixels with 8.16km < z < 19.2km
    ! ####################################################################################
    do ncol=1,Ncolumns
       do i=1,Npoints          
          do nlev=1,23 ! from 19.2km until 8.16km
               p1 = pplay(1,nlev)

             ! Avoid zero values
             if( (cldy(i,ncol,nlev).eq.1.) .and. (ATBperp(i,ncol,nlev).gt.0.) )then
                ! Computation of the ATBperp along the phase discrimination line
                ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 + (ATB(i,ncol,nlev)**4)*beta50 + &
                     (ATB(i,ncol,nlev)**3)*gamma50 + (ATB(i,ncol,nlev)**2)*delta50 + &
                     ATB(i,ncol,nlev)*epsilon50 + zeta50     
                ! ########################################################################
                ! 4.1.a) Ice: ATBperp above the phase discrimination line
                ! ########################################################################
                if((ATBperp(i,ncol,nlev)-ATBperp_tmp) .ge. 0.)then ! Ice clouds

                   ! ICE with temperature above 273,15°K = Liquid (false ice)
                   if(tmp(i,nlev) .gt. 273.15) then ! Temperature above 273,15 K
                     ! Liquid: False ice corrected by the temperature to Liquid
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp ! False ice detection ==> added to Liquid
                                    
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,5) = lidarcldphase(i,nlev,5)+1._wp ! Keep the information "temperature criterium used"                      
                                                                              ! to classify the phase cloud
                      cldlayphase(i,ncol,4,2) = 1. ! tot cloud
                      if (p1 .gt. 0. .and. p1.lt.(440._wp*100._wp)) then ! high cloud
                         cldlayphase(i,ncol,3,2) = 1._wp
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then ! mid cloud
                         cldlayphase(i,ncol,2,2) = 1._wp
                      else ! low cloud
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                      cldlayphase(i,ncol,4,5) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then 
                         cldlayphase(i,ncol,3,5) = 1._wp
                      ! Middle cloud
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,5) = 1._wp
                      ! Low cloud
                      else 
                         cldlayphase(i,ncol,1,5) = 1._wp
                      endif
                   else
                      ! ICE with temperature below 273,15°K
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp
                      tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud 
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then 
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then 
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                ! ########################################################################
                ! 4.1.b) Liquid: ATBperp below the phase discrimination line
                ! ########################################################################
                else
                   ! Liquid with temperature above 231,15°K
                   if(tmp(i,nlev) .gt. 231.15_wp) then
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,2) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,2) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,2) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                   else
                      ! Liquid with temperature below 231,15°K = Ice (false liquid)
                      tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp ! false liquid detection ==> added to ice
                      lidarcldphase(i,nlev,4) = lidarcldphase(i,nlev,4)+1._wp
                      cldlayphase(i,ncol,4,4) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,4) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,4) = 1._wp
                      ! Low cloud
                      else
                         cldlayphase(i,ncol,1,4) = 1._wp
                      endif
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                endif ! end of discrimination condition
             endif ! end of cloud condition
          enddo ! end of altitude loop

          ! ##############################################################################
          ! 4.2) For Cloudy pixels with 0km < z < 8.16km
          ! ##############################################################################
          toplvlsat = 0
          do nlev=24,Nlevels! from 8.16km until 0km
             p1 = pplay(i,nlev)

             if((cldy(i,ncol,nlev) .eq. 1.) .and. (ATBperp(i,ncol,nlev) .gt. 0.) )then
                ! Computation of the ATBperp of the phase discrimination line
                ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 + (ATB(i,ncol,nlev)**4)*beta50 + &
                     (ATB(i,ncol,nlev)**3)*gamma50 + (ATB(i,ncol,nlev)**2)*delta50 + &
                     ATB(i,ncol,nlev)*epsilon50 + zeta50
                ! ########################################################################
                ! 4.2.a) Ice: ATBperp above the phase discrimination line
                ! ########################################################################
                ! ICE with temperature above 273,15°K = Liquid (false ice)
                if((ATBperp(i,ncol,nlev)-ATBperp_tmp) .ge. 0.)then ! Ice clouds
                   if(tmp(i,nlev) .gt. 273.15)then
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp ! false ice ==> liq
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,5) = lidarcldphase(i,nlev,5)+1._wp
                      cldlayphase(i,ncol,4,2) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then 
                         cldlayphase(i,ncol,3,2) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then 
                         cldlayphase(i,ncol,2,2) = 1._wp
                      ! Low cloud
                      else 
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                      
                      cldlayphase(i,ncol,4,5) = 1. ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,5) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,5) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,5) = 1._wp
                      endif
                   else
                      ! ICE with temperature below 273,15°K
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp
                     tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt.(680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                   
                ! ########################################################################
                ! 4.2.b) Liquid: ATBperp below the phase discrimination line
                ! ########################################################################
                else
                   ! Liquid with temperature above 231,15°K
                   if(tmp(i,nlev) .gt. 231.15)then
                      lidarcldphase(i,nlev,2) = lidarcldphase(i,nlev,2)+1._wp
                      tmpl(i,ncol,nlev)       = tmp(i,nlev)
                      cldlayphase(i,ncol,4,2) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,2) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,2) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,2) = 1._wp
                      endif
                   else
                      ! Liquid with temperature below 231,15°K = Ice (false liquid)
                      tmpi(i,ncol,nlev)       = tmp(i,nlev)
                      lidarcldphase(i,nlev,1) = lidarcldphase(i,nlev,1)+1._wp ! false liq ==> ice
                      lidarcldphase(i,nlev,4) = lidarcldphase(i,nlev,4)+1._wp ! false liq ==> ice
                      cldlayphase(i,ncol,4,4) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,4) = 1._wp
                      ! Middle   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,4) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,4) = 1._wp
                      endif
                      
                      cldlayphase(i,ncol,4,1) = 1._wp ! tot cloud
                      ! High cloud
                      if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                         cldlayphase(i,ncol,3,1) = 1._wp
                      ! Middle cloud   
                      else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                         cldlayphase(i,ncol,2,1) = 1._wp
                      ! Low cloud   
                      else
                         cldlayphase(i,ncol,1,1) = 1._wp
                      endif
                   endif
                endif ! end of discrimination condition
                
                toplvlsat=0
                
                ! Find the level of the highest cloud with SR>30
                if(x(i,ncol,nlev) .gt. S_cld_att) then ! SR > 30.
                    toplvlsat = nlev+1
                    exit
                endif
             endif ! end of cloud condition
          enddo ! end of altitude loop
          
          ! ##############################################################################
          ! Undefined phase: For a cloud located below another cloud with SR>30
          ! see Cesana and Chepfer 2013 Sect.III.2
          ! ##############################################################################
          if(toplvlsat.ne.0) then
             do nlev = toplvlsat,Nlevels
                p1 = pplay(i,nlev)
                if(cldy(i,ncol,nlev).eq.1.)then
                   lidarcldphase(i,nlev,3) = lidarcldphase(i,nlev,3)+1._wp
                   tmpu(i,ncol,nlev)       = tmp(i,nlev)
                   cldlayphase(i,ncol,4,3) = 1._wp ! tot cloud
                   ! High cloud
                   if (p1 .gt. 0. .and. p1 .lt. (440._wp*100._wp)) then
                      cldlayphase(i,ncol,3,3) = 1._wp
                   ! Middle cloud   
                   else if(p1 .ge. (440._wp*100._wp) .and. p1 .lt. (680._wp*100._wp)) then
                      cldlayphase(i,ncol,2,3) = 1._wp
                   ! Low cloud   
                   else
                      cldlayphase(i,ncol,1,3) = 1._wp
                   endif
                endif
             enddo
          endif
          toplvlsat=0
       enddo
    enddo
     
    ! ####################################################################################
    ! Computation of final cloud phase diagnosis
    ! ####################################################################################

    ! Compute the Ice percentage in cloud = ice/(ice+liq) as a function of the occurrences
    lidarcldphasetmp(:,:) = lidarcldphase(:,:,1)+lidarcldphase(:,:,2);
    WHERE (lidarcldphasetmp(:,:) .gt. 0.)
       lidarcldphase(:,:,6)=lidarcldphase(:,:,1)/lidarcldphasetmp(:,:)
    ELSEWHERE
       lidarcldphase(:,:,6) = undef
    ENDWHERE
    
    ! Compute Phase 3D Cloud Fraction
    !WHERE (nsub(:,Nlevels:1:-1) .gt. 0.0 )
    WHERE (nsub(:,:) .gt. 0.0 )  
       lidarcldphase(:,:,1)=lidarcldphase(:,:,1)/nsub(:,:)
       lidarcldphase(:,:,2)=lidarcldphase(:,:,2)/nsub(:,:)
       lidarcldphase(:,:,3)=lidarcldphase(:,:,3)/nsub(:,:)
       lidarcldphase(:,:,4)=lidarcldphase(:,:,4)/nsub(:,:)
       lidarcldphase(:,:,5)=lidarcldphase(:,:,5)/nsub(:,:)
    ELSEWHERE
       lidarcldphase(:,:,1) = undef
       lidarcldphase(:,:,2) = undef
       lidarcldphase(:,:,3) = undef
       lidarcldphase(:,:,4) = undef
       lidarcldphase(:,:,5) = undef
    ENDWHERE

    ! Compute Phase low mid high cloud fractions
    do iz = 1, Ncat
       do i=1,Nphase-3
          do ic = 1, Ncolumns
             cldlayerphase(:,iz,i)  = cldlayerphase(:,iz,i)  + cldlayphase(:,ic,iz,i)
             cldlayerphasesum(:,iz) = cldlayerphasesum(:,iz) + cldlayphase(:,ic,iz,i)
          enddo
       enddo
    enddo
    do iz = 1, Ncat
       do i=4,5
          do ic = 1, Ncolumns
             cldlayerphase(:,iz,i) = cldlayerphase(:,iz,i) + cldlayphase(:,ic,iz,i)
          enddo
       enddo
    enddo
    
    ! Compute the Ice percentage in cloud = ice/(ice+liq)
    cldlayerphasetmp(:,:)=cldlayerphase(:,:,1)+cldlayerphase(:,:,2)
    WHERE (cldlayerphasetmp(:,:).gt. 0.)
       cldlayerphase(:,:,6)=cldlayerphase(:,:,1)/cldlayerphasetmp(:,:)
    ELSEWHERE
       cldlayerphase(:,:,6) = undef
    ENDWHERE
    
    do i=1,Nphase-1
       WHERE ( cldlayerphasesum(:,:).gt.0.0 )
          cldlayerphase(:,:,i) = (cldlayerphase(:,:,i)/cldlayerphasesum(:,:)) * cldlayer(:,:)
       ENDWHERE
    enddo
    
    do i=1,Npoints
       do iz=1,Ncat
          checkcldlayerphase=0.
          checkcldlayerphase2=0.
          if (cldlayerphasesum(i,iz) .gt. 0.0 )then
             do ic=1,Nphase-3
                checkcldlayerphase = checkcldlayerphase+cldlayerphase(i,iz,ic)
             enddo
             checkcldlayerphase2 = cldlayer(i,iz)-checkcldlayerphase
             if((checkcldlayerphase2 .gt. 0.01) .or. (checkcldlayerphase2 .lt. -0.01) ) print *, checkcldlayerphase,cldlayer(i,iz)
          endif
       enddo
    enddo
    
    do i=1,Nphase-1
       WHERE (nsublayer(:,:) .eq. 0.0)
          cldlayerphase(:,:,i) = undef
       ENDWHERE
    enddo
 
    ! Compute Phase 3D as a function of temperature
    do nlev=1,Nlevels
       do ncol=1,Ncolumns
          do i=1,Npoints
             do itemp=1,Ntemp
                if(tmpi(i,ncol,nlev).gt.0.)then
                   if((tmpi(i,ncol,nlev) .ge. tempmod(itemp)) .and. (tmpi(i,ncol,nlev) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(i,itemp,2)=lidarcldtemp(i,itemp,2)+1._wp
                   endif
                elseif(tmpl(i,ncol,nlev) .gt. 0.)then
                   if((tmpl(i,ncol,nlev) .ge. tempmod(itemp)) .and. (tmpl(i,ncol,nlev) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(i,itemp,3)=lidarcldtemp(i,itemp,3)+1._wp
                   endif
                elseif(tmpu(i,ncol,nlev) .gt. 0.)then
                   if((tmpu(i,ncol,nlev) .ge. tempmod(itemp)) .and. (tmpu(i,ncol,nlev) .lt. tempmod(itemp+1)) )then
                      lidarcldtemp(i,itemp,4)=lidarcldtemp(i,itemp,4)+1._wp
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
    
    ! Check temperature cloud fraction
    do i=1,Npoints
       do itemp=1,Ntemp
          checktemp=lidarcldtemp(i,itemp,2)+lidarcldtemp(i,itemp,3)+lidarcldtemp(i,itemp,4)
          !if(checktemp .NE. lidarcldtemp(i,itemp,1))then
          !   print *, i,itemp
          !   print *, lidarcldtemp(i,itemp,1:4)
          !endif
          
       enddo
    enddo
    
    ! Compute the Ice percentage in cloud = ice/(ice+liq)
    sumlidarcldtemp(:,:)=lidarcldtemp(:,:,2)+lidarcldtemp(:,:,3)    
    WHERE(sumlidarcldtemp(:,:) .gt. 0.)
       lidarcldtemp(:,:,5)=lidarcldtemp(:,:,2)/sumlidarcldtemp(:,:)
    ELSEWHERE
       lidarcldtemp(:,:,5)=undef
    ENDWHERE
    
    do i=1,4
       WHERE(lidarcldtempind(:,:) .gt. 0.)
          lidarcldtemp(:,:,i) = lidarcldtemp(:,:,i)/lidarcldtempind(:,:)
       ELSEWHERE
          lidarcldtemp(:,:,i) = undef
       ENDWHERE
    enddo
    
    RETURN
  END SUBROUTINE COSP_CLDFRAC
  
  ! ####################################################################################
  ! SUBROUTINE cosp_cldfrac_nophase
  ! Conventions: Ncat must be equal to 4
  ! ####################################################################################
  SUBROUTINE COSP_CLDFRAC_NOPHASE(Npoints,Ncolumns,Nlevels,Ncat,x,ATB,pplay,      &
       S_att,S_cld,S_cld_att,undef,lidarcld,cldlayer)
    
    ! Inputs
    integer,intent(in) :: &
         Npoints,  & ! Number of gridpoints
         Ncolumns, & ! Number of subcolumns
         Nlevels,  & ! Number of vertical levels
         Ncat        ! Number of cloud layer types
    real(wp),intent(in) :: &
         S_att,    & !
         S_cld,    & !
         S_cld_att,& ! New threshold for undefine cloud phase detection
         undef       ! Undefined value
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: &
         x,        & ! 
         ATB         ! 3D attenuated backscatter
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
         pplay       ! Pressure
    
    ! Outputs
    real(wp),intent(out),dimension(Npoints,Nlevels) :: &
         lidarcld      ! 3D cloud fraction
    real(wp),intent(out),dimension(Npoints,Ncat) :: &
         cldlayer      ! Low, middle, high, total cloud fractions
    
    ! Local variables
    integer  :: &
         ip, k, iz, ic
    real(wp) :: &
         p1
    real(wp),dimension(Npoints,Nlevels) :: &
         nsub
    real(wp),dimension(Npoints,Ncolumns,Ncat) :: &
         cldlay,nsublay   
    real(wp),dimension(Npoints,Ncat) :: &
         nsublayer
    real(wp),dimension(Npoints,Ncolumns,Nlevels) :: &   
         cldy, & ! 
         srok    !
    
    ! ####################################################################################
    ! 1) Initialize    
    ! ####################################################################################
    lidarcld           = 0._wp
    nsub               = 0._wp
    cldlay             = 0._wp
    nsublay            = 0._wp
    
    ! ####################################################################################
    ! 2) Cloud detection
    ! ####################################################################################
    do k=1,Nlevels
       ! Cloud detection at subgrid-scale:
       where ((x(:,:,k) .gt. S_cld) .and. (x(:,:,k) .ne. undef) )
          cldy(:,:,k)=1._wp
       elsewhere
          cldy(:,:,k)=0._wp
       endwhere
       
       ! Number of usefull sub-columns:
       where ((x(:,:,k) .gt. S_att) .and. (x(:,:,k) .ne. undef) )
          srok(:,:,k)=1._wp
       elsewhere
          srok(:,:,k)=0._wp
       endwhere
    enddo    

    ! ####################################################################################
    ! 3) Grid-box 3D cloud fraction and layered cloud fractions(ISCCP pressure categories)
    ! ####################################################################################
    do k=1,Nlevels
       do ic = 1, Ncolumns
          do ip = 1, Npoints

             iz=1
             p1 = pplay(ip,k)
             if ( p1.gt.0. .and. p1.lt.(440._wp*100._wp)) then ! high clouds
                iz=3
             else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid clouds
                iz=2
             endif
             
             cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
             cldlay(ip,ic,4)  = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
             lidarcld(ip,k)   = lidarcld(ip,k) + cldy(ip,ic,k)
             
             nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
             nsublay(ip,ic,4)  = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
             nsub(ip,k)        = nsub(ip,k) + srok(ip,ic,k)
             
          enddo
       enddo
    enddo   
    
    ! Grid-box 3D cloud fraction
    where ( nsub(:,:).gt.0.0 )
       lidarcld(:,:) = lidarcld(:,:)/nsub(:,:)
    elsewhere
       lidarcld(:,:) = undef
    endwhere
    
    ! Layered cloud fractions
    cldlayer  = 0._wp
    nsublayer = 0._wp
    do iz = 1, Ncat
       do ic = 1, Ncolumns
          cldlayer(:,iz)  = cldlayer(:,iz)  + cldlay(:,ic,iz)
          nsublayer(:,iz) = nsublayer(:,iz) + nsublay(:,ic,iz)
       enddo
    enddo
    where (nsublayer(:,:) .gt. 0.0)
       cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
    elsewhere
       cldlayer(:,:) = undef
    endwhere

    RETURN
  END SUBROUTINE COSP_CLDFRAC_NOPHASE

    ! ####################################################################################
    ! SUBROUTINE cosp_opaq
    ! Conventions: Ntype must be equal to 3
    ! ####################################################################################
    SUBROUTINE COSP_OPAQ(Npoints,Ncolumns,Nlevels,Ntype,tmp,x,S_att,S_cld,undef,lidarcldtype,   &
                         cldtype,cldtypetemp,cldtypemeanz,cldtypemeanzse,cldthinemis,vgrid_z,   &
                         surfelev)

    ! Local parameter
    real(wp),parameter  :: &
       S_att_opaq = 0.06_wp, & ! Fully Attenuated threshold (Guzman et al. 2017, JGR-Atmospheres)
       eta = 0.6_wp            ! Multiple-scattering factor (Vaillant de Guelis et al. 2017a, AMT)

    ! Inputs
    integer,intent(in) :: &
       Npoints,  & ! Number of gridpoints
       Ncolumns, & ! Number of subcolumns
       Nlevels,  & ! Number of vertical levels
       Ntype       ! Number of OPAQ cloud types (opaque, thin clouds and z_opaque)
    real(wp),intent(in) :: &
       S_att,    & ! Fully Attenuated legacy threshold
       S_cld,    & ! Cloud detection threshold
       undef       ! Undefined value
    real(wp),intent(in),dimension(Nlevels) :: &
       vgrid_z     ! mid-level vertical profile altitude (subcolumns)
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: &
       x           ! SR profiles (subcolumns)
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
       tmp         ! Temperature profiles
    real(wp),intent(in),dimension(Npoints) :: &
       surfelev    ! Surface Elevation (SE)

    ! Outputs
    real(wp),intent(out),dimension(Npoints,Nlevels,Ntype+1) :: &
       lidarcldtype   ! 3D OPAQ product fraction (opaque clouds, thin clouds, z_opaque, opacity)
    real(wp),intent(out),dimension(Npoints,Ntype) :: &
       cldtype,     & ! Opaque/thin cloud covers + z_opaque altitude
       cldtypetemp    ! Opaque and thin clouds + z_opaque temperature
    real(wp),intent(out),dimension(Npoints,2) :: &
       cldtypemeanz   ! Opaque and thin clouds altitude
    real(wp),intent(out),dimension(Npoints,3) :: & 
       cldtypemeanzse ! Opaque, thin clouds and z_opaque altitude with respect to SE
    real(wp),intent(out),dimension(Npoints) :: &
       cldthinemis    ! Thin clouds emissivity
   
    ! Local variables
    integer  :: &
       ip, k, zopac, ic, iz, z_top, z_base, topcloud
    real(wp)  :: &
       srmean, srcount, trans2, tau_app, tau_vis, tau_ir, cloudemis
    real(wp),dimension(Npoints) :: &
       count_emis
    real(wp),dimension(Npoints,Nlevels) :: &
       nsub, nsubopaq 
    real(wp),dimension(Npoints,Ncolumns,Ntype+1) :: & ! Opaque, thin, z_opaque and all cloud cover
       cldlay, nsublay   
    real(wp),dimension(Npoints,Ntype) :: &
       nsublayer
    real(wp),dimension(Npoints,Ncolumns,Nlevels) :: &   
       cldy,     & ! 
       cldyopaq, & ! 
       srok,     & !
       srokopaq    !

    ! ####################################################################################
    ! 1) Initialize    
    ! ####################################################################################
    cldtype(:,:)          = 0._wp
    cldtypetemp(:,:)      = 0._wp
    cldtypemeanz(:,:)     = 0._wp
    cldtypemeanzse(:,:)   = 0._wp
    cldthinemis(:)        = 0._wp
    count_emis(:)         = 0._wp
    lidarcldtype(:,:,:)   = 0._wp
    nsub                  = 0._wp
    nsubopaq              = 0._wp
    cldlay                = 0._wp
    nsublay               = 0._wp
    nsublayer             = 0._wp

    ! ####################################################################################
    ! 2) Cloud detection and Fully attenuated layer detection
    ! ####################################################################################
    do k=1,Nlevels
       ! Cloud detection at subgrid-scale:
       where ( (x(:,:,k) .gt. S_cld) .and. (x(:,:,k) .ne. undef) )
          cldy(:,:,k)=1._wp
       elsewhere
          cldy(:,:,k)=0._wp
       endwhere
       ! Fully attenuated layer detection at subgrid-scale:
       where ( (x(:,:,k) .lt. S_att_opaq) .and. (x(:,:,k) .ge. 0.) .and. (x(:,:,k) .ne. undef) )
          cldyopaq(:,:,k)=1._wp
       elsewhere
          cldyopaq(:,:,k)=0._wp
       endwhere


       ! Number of usefull sub-column layers:
       where ( (x(:,:,k) .gt. S_att) .and. (x(:,:,k) .ne. undef) )
          srok(:,:,k)=1._wp
       elsewhere
          srok(:,:,k)=0._wp
       endwhere
       ! Number of usefull sub-columns layers for z_opaque 3D fraction:
       where ( (x(:,:,k) .ge. 0.) .and. (x(:,:,k) .ne. undef) )
          srokopaq(:,:,k)=1._wp
       elsewhere
          srokopaq(:,:,k)=0._wp
       endwhere
    enddo

    ! ####################################################################################
    ! 3) Grid-box 3D OPAQ product fraction and cloud type cover (opaque/thin) + mean z_opaque
    ! ####################################################################################

    do k=1,Nlevels
       do ic = 1, Ncolumns
          do ip = 1, Npoints

             cldlay(ip,ic,1)   = MAX(cldlay(ip,ic,1),cldyopaq(ip,ic,k)) ! Opaque cloud
             cldlay(ip,ic,4)   = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))     ! All cloud

             nsublay(ip,ic,1)  = MAX(nsublay(ip,ic,1),srok(ip,ic,k))
             nsublay(ip,ic,2)  = MAX(nsublay(ip,ic,2),srok(ip,ic,k))
!             nsublay(ip,ic,4)  = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
             nsub(ip,k)        = nsub(ip,k) + srok(ip,ic,k)
             nsubopaq(ip,k)    = nsubopaq(ip,k) + srokopaq(ip,ic,k)

          enddo
       enddo
    enddo   

    ! OPAQ variables
    do ic = 1, Ncolumns
       do ip = 1, Npoints

          ! Declaring non-opaque cloudy profiles as thin cloud profiles
          if ( cldlay(ip,ic,4).gt. 0. .and. cldlay(ip,ic,1) .eq. 0. ) then
             cldlay(ip,ic,2)  =  1._wp
          endif

          ! Filling in 3D and 2D variables

          ! Opaque cloud profiles
          if ( cldlay(ip,ic,1) .eq. 1. ) then
             zopac = 0._wp
             z_top = 0._wp
             do k=1,Nlevels-1
                ! Declaring z_opaque altitude and opaque cloud fraction for 3D and 2D variables
                ! From SFC-2-TOA ( actually from vgrid_z(SFC+1) = vgrid_z(Nlevels-1) )
                if ( cldy(ip,ic,Nlevels-k) .eq. 1. .and. zopac .eq. 0. ) then
                   lidarcldtype(ip,Nlevels-k + 1,3) = lidarcldtype(ip,Nlevels-k + 1,3) + 1._wp
                   cldlay(ip,ic,3)                  = vgrid_z(Nlevels-k+1)      ! z_opaque altitude
                   nsublay(ip,ic,3)                 = 1._wp
                   zopac = Nlevels-k+1                        ! z_opaque vertical index on vgrid_z
                endif
                if ( cldy(ip,ic,Nlevels-k) .eq. 1. ) then
                   lidarcldtype(ip,Nlevels-k ,1)    = lidarcldtype(ip,Nlevels-k ,1) + 1._wp
                   z_top = Nlevels-k    ! top cloud layer vertical index on vgrid_z
                endif
             enddo
             ! Summing opaque cloud mean temperatures and altitudes
             ! as defined in Vaillant de Guelis et al. 2017a, AMT
             if (zopac .ne. 0) then 
                cldtypetemp(ip,1) = cldtypetemp(ip,1) + ( tmp(ip,zopac) + tmp(ip,z_top) )/2.
                cldtypetemp(ip,3) = cldtypetemp(ip,3) + tmp(ip,zopac)                 ! z_opaque
                cldtypemeanz(ip,1) = cldtypemeanz(ip,1) + ( vgrid_z(zopac) + vgrid_z(z_top) )/2.
                cldtypemeanzse(ip,1) = cldtypemeanzse(ip,1) + (( vgrid_z(zopac) + vgrid_z(z_top) )/2.) - surfelev(ip)
                cldtypemeanzse(ip,3) = cldtypemeanzse(ip,3) + ( vgrid_z(zopac) - surfelev(ip) )
             else
                cldlay(ip,ic,1) = 0
             endif
          endif

          ! Thin cloud profiles
          if ( cldlay(ip,ic,2) .eq. 1. ) then
             topcloud = 0._wp
             z_top = 0._wp
             z_base = 0._wp
             do k=1,Nlevels
                ! Declaring thin cloud fraction for 3D variable
                ! From TOA-2-SFC
                if ( cldy(ip,ic,k) .eq. 1. .and. topcloud .eq. 1. ) then
                   lidarcldtype(ip,k,2) = lidarcldtype(ip,k,2) + 1._wp
                   z_base = k ! bottom cloud layer
                endif
                if ( cldy(ip,ic,k) .eq. 1. .and. topcloud .eq. 0. ) then
                   lidarcldtype(ip,k,2) = lidarcldtype(ip,k,2) + 1._wp
                   z_top = k  ! top cloud layer
                   z_base = k ! bottom cloud layer
                   topcloud = 1._wp
                endif
             enddo
             ! Computing mean emissivity using layers below the bottom cloud layer to the surface
             srmean = 0._wp
             srcount = 0._wp
             cloudemis = 0._wp
             do k=z_base+1,Nlevels
                if (  (x(ip,ic,k) .gt. S_att_opaq) .and. (x(ip,ic,k) .lt. 1.0) .and. (x(ip,ic,k) .ne. undef)  ) then
                   srmean = srmean + x(ip,ic,k)
                   srcount = srcount + 1.
                endif
             enddo
             ! If clear sky layers exist below bottom cloud layer
             if ( srcount .gt. 0. ) then
                trans2 = srmean/srcount              ! thin cloud transmittance**2
                tau_app = -(log(trans2))/2.          ! apparent cloud optical depth
                tau_vis = tau_app/eta                ! cloud visible optical depth (multiple scat.)
                tau_ir = tau_vis/2.                  ! approx. relation between visible and IR ODs
                cloudemis = 1. - exp(-tau_ir)        ! no diffusion in IR considered : emis = 1-T
                count_emis(ip) = count_emis(ip) + 1.
             endif
             ! Summing thin cloud mean temperatures and altitudes
             ! as defined in Vaillant de Guelis et al. 2017a, AMT
             cldtypetemp(ip,2) = cldtypetemp(ip,2) + ( tmp(ip,z_base) + tmp(ip,z_top) )/2.
             cldtypemeanz(ip,2) = cldtypemeanz(ip,2) + ( vgrid_z(z_base) + vgrid_z(z_top) )/2.
             cldtypemeanzse(ip,2) = cldtypemeanzse(ip,2) + (( vgrid_z(z_base) + vgrid_z(z_top) )/2.) - surfelev(ip)
             cldthinemis(ip) = cldthinemis(ip) + cloudemis
          endif

       enddo
    enddo   

    ! 3D cloud types fraction (opaque=1 and thin=2 clouds)
    where ( nsub(:,:) .gt. 0. )
       lidarcldtype(:,:,1) = lidarcldtype(:,:,1)/nsub(:,:)
       lidarcldtype(:,:,2) = lidarcldtype(:,:,2)/nsub(:,:)
    elsewhere
       lidarcldtype(:,:,1) = undef
       lidarcldtype(:,:,2) = undef
    endwhere
    ! 3D z_opaque fraction (=3)
    where ( nsubopaq(:,:) .gt. 0. )
       lidarcldtype(:,:,3) = lidarcldtype(:,:,3)/nsubopaq(:,:)
    elsewhere
       lidarcldtype(:,:,3) = undef
       lidarcldtype(:,:,4) = undef !declaring undef for opacity as well
    endwhere
    ! 3D opacity fraction (=4) !Summing z_opaque fraction from TOA(k=1) to SFC(k=Nlevels)
       lidarcldtype(:,1,4) = lidarcldtype(:,1,3) !top layer equal to 3D z_opaque fraction
    do ip = 1, Npoints
       do k = 2, Nlevels
          if ( (lidarcldtype(ip,k,3) .ne. undef) .and. (lidarcldtype(ip,k-1,4) .ne. undef) ) then
             lidarcldtype(ip,k,4) = lidarcldtype(ip,k,3) + lidarcldtype(ip,k-1,4)
          else
             lidarcldtype(ip,k,4) = undef
          endif
       enddo
    enddo

    ! Layered cloud types (opaque, thin and z_opaque 2D variables)

    do iz = 1, Ntype
       do ic = 1, Ncolumns
          cldtype(:,iz)  = cldtype(:,iz)  + cldlay(:,ic,iz)
          nsublayer(:,iz) = nsublayer(:,iz) + nsublay(:,ic,iz)
       enddo
    enddo

    ! Mean temperature and altitude
    where (cldtype(:,1) .gt. 0.)
       cldtypetemp(:,1) = cldtypetemp(:,1)/cldtype(:,1) ! opaque cloud temp
       cldtypetemp(:,3) = cldtypetemp(:,3)/cldtype(:,1) ! z_opaque
       cldtypemeanz(:,1) = cldtypemeanz(:,1)/cldtype(:,1) ! opaque cloud alt
       cldtypemeanzse(:,1) = cldtypemeanzse(:,1)/cldtype(:,1) ! opaque cloud alt - SE 
       cldtypemeanzse(:,3) = cldtypemeanzse(:,3)/cldtype(:,1) ! z_opaque - SE 
    elsewhere
       cldtypetemp(:,1) = undef
       cldtypetemp(:,3) = undef
       cldtypemeanz(:,1) = undef
       cldtypemeanzse(:,1) = undef
       cldtypemeanzse(:,3) = undef
    endwhere

    where (cldtype(:,2) .gt. 0.) ! thin cloud
       cldtypetemp(:,2) = cldtypetemp(:,2)/cldtype(:,2)
       cldtypemeanz(:,2) = cldtypemeanz(:,2)/cldtype(:,2)
       cldtypemeanzse(:,2) = cldtypemeanzse(:,2)/cldtype(:,2)
    elsewhere
       cldtypetemp(:,2) = undef
       cldtypemeanz(:,2) = undef
       cldtypemeanzse(:,2) = undef
    endwhere

    ! Mean thin cloud emissivity
    where (count_emis(:) .gt. 0.) ! thin cloud
       cldthinemis(:) = cldthinemis(:)/count_emis(:)
    elsewhere
       cldthinemis(:) = undef
    endwhere

    where (nsublayer(:,:) .gt. 0.)
       cldtype(:,:) = cldtype(:,:)/nsublayer(:,:)
    elsewhere
       cldtype(:,:) = undef
    endwhere

  END SUBROUTINE COSP_OPAQ

end module mod_lidar_simulator
