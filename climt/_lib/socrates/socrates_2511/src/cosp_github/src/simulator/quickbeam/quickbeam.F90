! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
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
! 11/2005: John Haynes - Created
! 09/2006  placed into subroutine form (Roger Marchand,JMH)
! 08/2007  added equivalent volume spheres, Z and N scalling most distrubtion types (Roger Marchand)
! 01/2008  'Do while' to determine if hydrometeor(s) present in volume
!           changed for vectorization purposes (A. Bodas-Salcedo)
!
! 07/2010  V3.0 ... Modified to load or save scale factors to disk as a Look-Up Table (LUT)
!  ... All hydrometeor and radar simulator properties now included in hp structure
!  ... hp structure should be initialized by call to radar_simulator_init prior 
!  ... to calling this subroutine.  
!     Also ... Support of Morrison 2-moment style microphyscis (Np_matrix) added 
!  ... Changes implement by Roj Marchand following work by Laura Fowler
!
!   10/2011  Modified ngate loop to go in either direction depending on flag 
!     hp%radar_at_layer_one.  This affects the direction in which attenuation is summed.
!
!     Also removed called to AVINT for gas and hydrometeor attenuation and replaced with simple
!     summation. (Roger Marchand)
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module quickbeam
  USE COSP_KINDS,           ONLY: wp
  USE MOD_COSP_CONFIG,      ONLY: R_UNDEF,cloudsat_histRef,use_vgrid,vgrid_zl,vgrid_zu,&
                                  pClass_noPrecip, pClass_Rain1, pClass_Rain2, pClass_Rain3,&
                                  pClass_Snow1, pClass_Snow2, pClass_Mixed1, pClass_Mixed2, &
                                  pClass_Rain4, pClass_default, Zenonbinval, Zbinvallnd,    &
                                  N_HYDRO,nCloudsatPrecipClass,cloudsat_preclvl

  USE MOD_COSP_STATS,       ONLY: COSP_LIDAR_ONLY_CLOUD,hist1D,COSP_CHANGE_VERTICAL_GRID
  implicit none

  integer,parameter :: &
       maxhclass     = 20,  & ! Qucikbeam maximum number of hydrometeor classes.
       nRe_types     = 550, & ! Quickbeam maximum number or Re size bins allowed in N and Z_scaled look up table.
       nd            = 85,  & ! Qucikbeam number of discrete particles used in construction DSDs.
       mt_ntt        = 39,  & ! Quickbeam number of temperatures in mie LUT.
       Re_BIN_LENGTH = 10,  & ! Quickbeam minimum Re interval in scale LUTs  
       Re_MAX_BIN    = 250    ! Quickbeam maximum Re interval in scale LUTs
  real(wp),parameter :: &
       dmin          = 0.1, & ! Quickbeam minimum size of discrete particle
       dmax          = 10000. ! Quickbeam maximum size of discrete particle
  
  !djs logical :: radar_at_layer_one   ! If true radar is assume to be at the edge 
                                  ! of the first layer, if the first layer is the
                                  ! surface than a ground-based radar.   If the
                                  ! first layer is the top-of-atmosphere, then
                                  ! a space borne radar.

  ! ##############################################################################################
  type radar_cfg
     ! Radar properties
     real(wp) :: freq,k2
     integer  :: nhclass               ! Number of hydrometeor classes in use
     integer  :: use_gas_abs, do_ray
     logical  :: radar_at_layer_one    ! If true radar is assume to be at the edge 
                                       ! of the first layer, if the first layer is the
                                       ! surface than a ground-based radar.   If the
                                       ! first layer is the top-of-atmosphere, then
                                       ! a space borne radar.
     
     ! Variables used to store Z scale factors
     character(len=240)                             :: scale_LUT_file_name
     logical                                        :: load_scale_LUTs, update_scale_LUTs
     logical,  allocatable, dimension(:,:)   :: N_scale_flag
     logical,  allocatable, dimension(:,:,:) :: Z_scale_flag, Z_scale_added_flag
     real(wp), allocatable, dimension(:,:,:) :: Ze_scaled, Zr_scaled, kr_scaled
     real(wp), allocatable, dimension(:,:,:) :: fc, rho_eff
     real(wp), allocatable, dimension(:)     :: base_list, step_list
  end type radar_cfg

contains
  ! ######################################################################################
  ! SUBROUTINE quickbeam_subcolumn
  ! ######################################################################################
  subroutine quickbeam_subcolumn(rcfg,nprof,ngate,hgt_matrix,z_vol,kr_vol,g_vol,dBZe,Ze_non)

    ! INPUTS
    type(radar_cfg),intent(inout) :: &
         rcfg             ! Derived type for radar simulator setup
    integer,intent(in) :: &
         nprof,         & ! Number of hydrometeor profiles
         ngate            ! Number of vertical layers
    real(wp),intent(in),dimension(nprof,ngate) :: &
         hgt_matrix,    & ! Height of hydrometeors (km)
         z_vol,         & ! Effective reflectivity factor (mm^6/m^3)
         kr_vol,        & ! Attenuation coefficient hydro (dB/km)
         g_vol            ! Attenuation coefficient gases (dB/km)
    
    ! OUTPUTS
    real(wp), intent(out),dimension(nprof,ngate) :: &
         Ze_non,        & ! Radar reflectivity without attenuation (dBZ)
         dBZe             ! Effective radar reflectivity factor (dBZ)

    ! LOCAL VARIABLES
    integer :: k,pr,start_gate,end_gate,d_gate
    real(wp),dimension(nprof,ngate) :: &
         g_to_vol,      & ! Gaseous atteunation, radar to vol (dB)
         a_to_vol         ! Hydromets attenuation, radar to vol (dB) 

    ! Load scaling matricies from disk -- but only the first time this subroutine is called
    if(rcfg%load_scale_LUTs) then
       call load_scale_LUTs(rcfg)
       rcfg%load_scale_LUTs=.false.
       rcfg%Z_scale_added_flag = .false. ! will be set true if scaling Look Up Tables are modified during run
    endif

    ! Initialization
    g_to_vol = 0._wp
    a_to_vol = 0._wp

    ! Loop over each range gate (ngate) ... starting with layer closest to the radar !
    if(rcfg%radar_at_layer_one) then
       start_gate = 1
       end_gate   = ngate
       d_gate     = 1
    else
       start_gate = ngate
       end_gate   = 1
       d_gate     = -1
    endif
    do k=start_gate,end_gate,d_gate
       ! Loop over each profile (nprof)
       do pr=1,nprof
          ! Attenuation due to hydrometeors between radar and volume
          
          ! NOTE old scheme integrates attenuation only for the layers ABOVE
          ! the current layer ... i.e. 1 to k-1 rather than 1 to k ...
          ! which may be a problem.   ROJ
          ! in the new scheme I assign half the attenuation to the current layer
          if(d_gate==1) then
             ! dheight calcuations assumes hgt_matrix points are the cell mid-points.
             if (k>2) then
                ! add to previous value to half of above layer + half of current layer
                a_to_vol(pr,k)=  a_to_vol(pr,k-1) + &
                     (kr_vol(pr,k-1)+kr_vol(pr,k))*(hgt_matrix(pr,k-1)-hgt_matrix(pr,k))
             else
                a_to_vol(pr,k)=  kr_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k+1))
             endif
          else   ! d_gate==-1
             if(k<ngate) then
                ! Add to previous value half of above layer + half of current layer
                a_to_vol(pr,k) = a_to_vol(pr,k+1) + &
                     (kr_vol(pr,k+1)+kr_vol(pr,k))*(hgt_matrix(pr,k+1)-hgt_matrix(pr,k))
             else
                a_to_vol(pr,k)= kr_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k-1))
             endif
          endif
          
          ! Attenuation due to gaseous absorption between radar and volume
          if ((rcfg%use_gas_abs == 1) .or. (rcfg%use_gas_abs == 2 .and. pr .eq. 1)) then
             if (d_gate==1) then
                if (k>1) then
                   ! Add to previous value to half of above layer + half of current layer
                   g_to_vol(pr,k) =  g_to_vol(pr,k-1) + &
                        0.5*(g_vol(pr,k-1)+g_vol(pr,k))*(hgt_matrix(pr,k-1)-hgt_matrix(pr,k))
                else
                   g_to_vol(pr,k)=  0.5_wp*g_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k+1))
                endif
             else   ! d_gate==-1
                if (k<ngate) then
                   ! Add to previous value to half of above layer + half of current layer
                   g_to_vol(pr,k) = g_to_vol(pr,k+1) + &
                        0.5_wp*(g_vol(pr,k+1)+g_vol(pr,k))*(hgt_matrix(pr,k+1)-hgt_matrix(pr,k))
                else
                   g_to_vol(pr,k)= 0.5_wp*g_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k-1))
                endif
             endif
          elseif(rcfg%use_gas_abs == 2) then
             ! Using value calculated for the first column
             g_to_vol(pr,k) = g_to_vol(1,k)
          elseif (rcfg%use_gas_abs == 0) then
             g_to_vol(pr,k) = 0._wp
          endif
       enddo   ! End loop over pr (profile)
    enddo ! End loop of k (range gate)
    
    ! Compute full and attenuated reflectivity
    where(z_vol(1:nprof,1:ngate) > 0._wp) 
      Ze_non(1:nprof,1:ngate) = 10._wp*log10(z_vol(1:nprof,1:ngate))
      dBZe(1:nprof,1:ngate) = Ze_non(1:nprof,1:ngate)-a_to_vol(1:nprof,1:ngate)-g_to_vol(1:nprof,1:ngate)
    elsewhere
      dBZe(1:nprof,1:ngate) = R_UNDEF
      Ze_non(1:nprof,1:ngate) = R_UNDEF
    end where 

    ! Save any updates made 
    if (rcfg%update_scale_LUTs) call save_scale_LUTs(rcfg)
 
  end subroutine quickbeam_subcolumn
  ! ######################################################################################
  ! SUBROUTINE quickbeam_column
  ! ######################################################################################
  subroutine quickbeam_column(npoints, ncolumns, nlevels, llm, DBZE_BINS, platform,      &
       Ze_tot, Ze_tot_non, land, surfelev, t2m, fracPrecipIce, zlev, zlev_half, cfad_ze, &
       cloudsat_precip_cover, cloudsat_pia)
    ! Inputs
    integer,intent(in) :: &
         npoints,    & ! Number of horizontal grid points
         ncolumns,   & ! Number of subcolumns
         nlevels,    & ! Number of vertical layers in OLD grid
         llm,        & ! Number of vertical layers in NEW grid
         DBZE_BINS     ! Number of bins for cfad.
    character(len=*),intent(in) :: &
         platform      ! Name of platform (e.g. cloudsat)
    real(wp),dimension(Npoints),intent(in) :: &
         land,               & ! Land/Sea mask. (1/0)
         surfelev,           & ! Surface Elevation (m)
         t2m                   ! Near-surface temperature
    real(wp),dimension(Npoints,Ncolumns),intent(in) :: &
         fracPrecipIce         ! Fraction of precipitation which is frozen.     (1)
    real(wp),intent(in),dimension(npoints,ncolumns,Nlevels) :: &
         Ze_tot,     & ! Effective reflectivity factor                 (dBZ)
         Ze_tot_non    ! Effective reflectivity factor w/o attenuation (dBZ)
    real(wp),intent(in),dimension(npoints,Nlevels) :: &
         zlev          ! Model full levels
    real(wp),intent(in),dimension(npoints,Nlevels) :: &
         zlev_half     ! Model half levels
         
    ! Outputs
    real(wp),intent(inout),dimension(npoints,DBZE_BINS,llm) :: &
         cfad_ze    !
    real(wp),dimension(Npoints,nCloudsatPrecipClass),intent(out) :: &
         cloudsat_precip_cover ! Model precip rate in by CloudSat precip flag
    real(wp),dimension(Npoints),intent(out) :: &
         cloudsat_pia          ! Cloudsat path integrated attenuation
    
    ! Local variables
    integer :: i,j
    real(wp) :: zstep
    real(wp),dimension(npoints,ncolumns,llm) :: ze_toti,ze_noni
    logical :: lcloudsat = .false.

    ! Which platforms to create diagnostics for?
    if (platform .eq. 'cloudsat') lcloudsat=.true.

    ! Create Cloudsat diagnostics.
    if (lcloudsat) then
       if (use_vgrid) then
          ! Regrid in the vertical (*NOTE* This routine requires SFC-2-TOA ordering, so flip
          ! inputs and outputs to maintain TOA-2-SFC ordering convention in COSP2.)
          call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,zlev(:,nlevels:1:-1),&
               zlev_half(:,nlevels:1:-1),Ze_tot(:,:,nlevels:1:-1),llm,vgrid_zl(llm:1:-1),&
               vgrid_zu(llm:1:-1),Ze_toti(:,:,llm:1:-1),log_units=.true.)
          
          ! Effective reflectivity histogram
          do i=1,Npoints
             do j=1,llm
                cfad_ze(i,:,j) = hist1D(Ncolumns,Ze_toti(i,:,j),DBZE_BINS,cloudsat_histRef)
             enddo
          enddo
          where(cfad_ze .ne. R_UNDEF) cfad_ze = cfad_ze/Ncolumns

          ! Compute cloudsat near-surface precipitation diagnostics
          ! First, regrid in the vertical Ze_tot_non.
          call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,zlev(:,nlevels:1:-1),&
               zlev_half(:,nlevels:1:-1),Ze_tot_non(:,:,nlevels:1:-1),llm,vgrid_zl(llm:1:-1),&
               vgrid_zu(llm:1:-1),Ze_noni(:,:,llm:1:-1),log_units=.true.)
          ! Compute the zstep distance between two atmopsheric layers
          zstep = vgrid_zl(1)-vgrid_zl(2)
          ! Now call routine to generate diagnostics.
          call cloudsat_precipOccurence(Npoints, Ncolumns, llm, N_HYDRO, Ze_toti, Ze_noni, &
               land, surfelev, t2m, fracPrecipIce, cloudsat_precip_cover, cloudsat_pia, zstep)
       else
          ! Effective reflectivity histogram
          do i=1,Npoints
             do j=1,llm
                cfad_ze(i,:,j) = hist1D(Ncolumns,Ze_tot(i,:,j),DBZE_BINS,cloudsat_histRef)
             enddo
          enddo
          where(cfad_ze .ne. R_UNDEF) cfad_ze = cfad_ze/Ncolumns
       endif
    endif

  end subroutine quickbeam_column
  ! ##############################################################################################
  ! ##############################################################################################
  ! ######################################################################################
  ! SUBROUTINE cloudsat_precipOccurence
  !
  ! Notes from July 2016: Add precip flag also looped over subcolumns
  ! Modified by Tristan L'Ecuyer (TSL) to add precipitation flagging
  ! based on CloudSat's 2C-PRECIP-COLUMN algorithm (Haynes et al, JGR, 2009).
  ! To mimic the satellite algorithm, this code applies thresholds to non-attenuated
  ! reflectivities, Ze_non, consistent with those outlined in Haynes et al, JGR (2009).
  !
  ! Procedures/Notes:
  !
  ! (1) If the 2-way attenuation exceeds 40 dB, the pixel will be flagged as 'heavy rain'
  ! consistent with the multiple-scattering analysis of Battaglia et al, JGR (2008).
  ! (2) Rain, snow, and mixed precipitation scenes are separated according to the fraction
  ! of the total precipitation hydrometeor mixing ratio that exists as ice.
  ! (3) The entire analysis is applied to the range gate from 480-960 m to be consistent with
  ! CloudSat's 720 m ground-clutter.
  ! (4) Only a single flag is assigned to each model grid since there is no variation in
  ! hydrometeor contents across a single model level.  Unlike CFADs, whose variation enters
  ! due to differening attenuation corrections from hydrometeors aloft, the non-attenuated
  ! reflectivities used in the computation of this flag cannot vary across sub-columns.
  !
  ! radar_prec_flag = 1-Rain possible 2-Rain probable 3-Rain certain 
  !                   4-Snow possible 5-Snow certain 
  !                   6-Mixed possible 7-Mixed certain 
  !                   8-Heavy Rain
  !                   9- default value
  !
  ! Modified by Dustin Swales (University of Colorado) for use with COSP2.
  ! *NOTE* All inputs (Ze_out, Ze_non_out, fracPrecipIce) are at a single level from the
  !        statistical output grid used by Cloudsat. This level is controlled by the
  !        parameter cloudsat_preclvl, defined in src/cosp_config.F90
  ! ######################################################################################
  subroutine cloudsat_precipOccurence(Npoints, Ncolumns, llm, Nhydro, Ze_out, Ze_non_out, &
       land, surfelev, t2m, fracPrecipIce,  cloudsat_precip_cover, cloudsat_pia, zstep)
 
    ! Inputs
    integer,intent(in) :: &
         Npoints,            & ! Number of columns
         Ncolumns,           & ! Numner of subcolumns
         Nhydro,             & ! Number of hydrometeor types
         llm                   ! Number of levels
    real(wp),dimension(Npoints),intent(in) :: &
         land,               & ! Land/Sea mask. (1/0)
         surfelev,           & ! Surface Elevation (m)
         t2m                   ! Near-surface temperature
    real(wp),dimension(Npoints,Ncolumns,llm),intent(in) :: &
         Ze_out,             & ! Effective reflectivity factor                  (dBZ)
         Ze_non_out            ! Effective reflectivity factor, w/o attenuation (dBZ)
    real(wp),dimension(Npoints,Ncolumns),intent(in) :: &
         fracPrecipIce         ! Fraction of precipitation which is frozen.     (1)
    real(wp),intent(in) :: &
         zstep                 ! Distance between two atmopsheric layers (m)

    ! Outputs 
    real(wp),dimension(Npoints,nCloudsatPrecipClass),intent(out) :: &
         cloudsat_precip_cover ! Model precip rate in by CloudSat precip flag
    real(wp),dimension(Npoints),intent(out) :: &
         cloudsat_pia          ! Cloudsat path integrated attenuation

    ! Local variables 
    integer,dimension(Npoints,Ncolumns) :: &
         cloudsat_pflag,      & ! Subcolumn precipitation flag
         cloudsat_precip_pia    ! Subcolumn path integrated attenutation.
    integer,dimension(Npoints) :: &
         cloudsat_preclvl_index ! Altitude index for precip flags calculation
                                ! in 40-level grid (one layer above surfelev) 
    integer :: pr,i,k
    real(wp) :: Zmax
    
    ! Initialize 
    cloudsat_pflag(:,:)        = pClass_default
    cloudsat_precip_pia(:,:)   = 0._wp
    cloudsat_precip_cover(:,:) = 0._wp
    cloudsat_pia(:)            = 0._wp
    cloudsat_preclvl_index(:)  = 0._wp

    ! Computing altitude index for precip flags calculation
    cloudsat_preclvl_index(:) = cloudsat_preclvl - floor( surfelev(:)/zstep )

    ! ######################################################################################
    ! SUBCOLUMN processing
    ! ######################################################################################
    do i=1, Npoints
       do pr=1,Ncolumns
          ! Compute precipitation flag
          ! ################################################################################
          ! 1) Oceanic points.
          ! ################################################################################
          if (land(i) .eq. 0) then

             ! 1a) Compute the PIA in all profiles containing hydrometeors
             if ( (Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.-100) .and. &
                  (Ze_out(i,pr,cloudsat_preclvl_index(i)).gt.-100) ) then
                if ( (Ze_non_out(i,pr,cloudsat_preclvl_index(i)).lt.100) .and. &
                     (Ze_out(i,pr,cloudsat_preclvl_index(i)).lt.100) ) then
                   cloudsat_precip_pia(i,pr) = &
                   Ze_non_out(i,pr,cloudsat_preclvl_index(i)) - Ze_out(i,pr,cloudsat_preclvl_index(i))
                endif
             endif

             ! Snow
             if(fracPrecipIce(i,pr).gt.0.9) then
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(2)) then
                   cloudsat_pflag(i,pr) = pClass_Snow2                   ! TSL: Snow certain
                endif
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(4).and. &
                     Ze_non_out(i,pr,cloudsat_preclvl_index(i)).le.Zenonbinval(2)) then
                   cloudsat_pflag(i,pr) = pClass_Snow1                   ! TSL: Snow possible
                endif
             endif
             
             ! Mixed
             if(fracPrecipIce(i,pr).gt.0.1.and.fracPrecipIce(i,pr).le.0.9) then
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(2)) then
                   cloudsat_pflag(i,pr) = pClass_Mixed2                  ! TSL: Mixed certain
                endif
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(4).and. &
                     Ze_non_out(i,pr,cloudsat_preclvl_index(i)).le.Zenonbinval(2)) then
                   cloudsat_pflag(i,pr) = pClass_Mixed1                  ! TSL: Mixed possible
                endif
             endif
             
             ! Rain
             if(fracPrecipIce(i,pr).le.0.1) then
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(1)) then
                   cloudsat_pflag(i,pr) = pClass_Rain3                   ! TSL: Rain certain
                endif
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(3).and. &
                     Ze_non_out(i,pr,cloudsat_preclvl_index(i)).le.Zenonbinval(1)) then
                   cloudsat_pflag(i,pr) = pClass_Rain2                   ! TSL: Rain probable
                endif
                if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).gt.Zenonbinval(4).and. &
                     Ze_non_out(i,pr,cloudsat_preclvl_index(i)).le.Zenonbinval(3)) then
                   cloudsat_pflag(i,pr) = pClass_Rain1                   ! TSL: Rain possible
                endif
                if(cloudsat_precip_pia(i,pr).gt.40) then
                   cloudsat_pflag(i,pr) = pClass_Rain4                   ! TSL: Heavy Rain
                endif
             endif
             
             ! No precipitation
             if(Ze_non_out(i,pr,cloudsat_preclvl_index(i)).le.-15) then
                cloudsat_pflag(i,pr) = pClass_noPrecip                   ! TSL: Not Raining
             endif
          endif ! Ocean points
          
          ! ################################################################################
          ! 2) Land points.
          !    *NOTE* For land points we go up a layer higher, so cloudsat_preclvl_index(i)-1
          !                  
          ! ################################################################################
          if (land(i) .eq. 1) then             
             ! 2a) Compute the PIA in all profiles containing hydrometeors
             if ( (Ze_non_out(i,pr,cloudsat_preclvl_index(i)-1).gt.-100) .and. &
                  (Ze_out(i,pr,cloudsat_preclvl_index(i)-1).gt.-100) ) then
                if ( (Ze_non_out(i,pr,cloudsat_preclvl_index(i)-1).lt.100) .and. &
                  (Ze_out(i,pr,cloudsat_preclvl_index(i)-1).lt.100) ) then
                   cloudsat_precip_pia(i,pr) = &
                   Ze_non_out(i,pr,cloudsat_preclvl_index(i)-1) - Ze_out(i,pr,cloudsat_preclvl_index(i)-1)
                endif
             endif

             ! Find Zmax, the maximum reflectivity value in the attenuated profile (Ze_out);
             Zmax=maxval(Ze_out(i,pr,:))

             ! Snow (T<273)
             if(t2m(i) .lt. 273._wp) then
                if(Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(5)) then
                   cloudsat_pflag(i,pr) = pClass_Snow2                      ! JEK: Snow certain
                endif
                if(Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(6) .and. &
                     Ze_out(i,pr,cloudsat_preclvl_index(i)-1).le.Zbinvallnd(5)) then
                   cloudsat_pflag(i,pr) = pClass_Snow1                      ! JEK: Snow possible
                endif
             endif
             
             ! Mized phase (273<T<275)
             if(t2m(i) .ge. 273._wp .and. t2m(i) .le. 275._wp) then
                if ((Zmax .gt. Zbinvallnd(1) .and. cloudsat_precip_pia(i,pr).gt.30) .or. &
                     (Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(4))) then
                   cloudsat_pflag(i,pr) = pClass_Mixed2                     ! JEK: Mixed certain
                endif
                if ((Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(6)  .and. &
                     Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .le. Zbinvallnd(4)) .and. &
                     (Zmax .gt. Zbinvallnd(5)) ) then
                   cloudsat_pflag(i,pr) = pClass_Mixed1                     ! JEK: Mixed possible
                endif
             endif

             ! Rain (T>275)
             if(t2m(i) .gt. 275) then
                if ((Zmax .gt. Zbinvallnd(1) .and. cloudsat_precip_pia(i,pr).gt.30) .or. &
                     (Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(2))) then
                   cloudsat_pflag(i,pr) = pClass_Rain3                      ! JEK: Rain certain
                endif
                if((Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(6)) .and. &
                     (Zmax .gt. Zbinvallnd(3))) then
                   cloudsat_pflag(i,pr) = pClass_Rain2                      ! JEK: Rain probable
                endif
                if((Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .gt. Zbinvallnd(6)) .and. &
                     (Zmax.lt.Zbinvallnd(3))) then
                   cloudsat_pflag(i,pr) = pClass_Rain1                      ! JEK: Rain possible
                endif
                if(cloudsat_precip_pia(i,pr).gt.40) then
                   cloudsat_pflag(i,pr) = pClass_Rain4                      ! JEK: Heavy Rain
                endif
             endif
             
             ! No precipitation
             if(Ze_out(i,pr,cloudsat_preclvl_index(i)-1) .le. -15) then
                cloudsat_pflag(i,pr) =  pClass_noPrecip                     ! JEK: Not Precipitating
             endif         
          endif ! Land points
       enddo ! Sub-columns
    enddo    ! Gridpoints

    ! ######################################################################################
    ! COLUMN processing
    ! ######################################################################################

    ! Aggregate subcolumns
    do i=1,Npoints
       ! Gridmean precipitation fraction for each precipitation type
       do k=1,nCloudsatPrecipClass
          if (any(cloudsat_pflag(i,:) .eq. k-1)) then
             cloudsat_precip_cover(i,k) = count(cloudsat_pflag(i,:) .eq. k-1)
          endif
       enddo
 
       ! Gridmean path integrated attenuation (pia)
       cloudsat_pia(i)=sum(cloudsat_precip_pia(i,:))
    enddo
  
    ! Normalize by number of subcolumns
    where ((cloudsat_precip_cover /= R_UNDEF).and.(cloudsat_precip_cover /= 0.0)) &
         cloudsat_precip_cover = cloudsat_precip_cover / Ncolumns
    where ((cloudsat_pia/= R_UNDEF).and.(cloudsat_pia/= 0.0)) &
         cloudsat_pia = cloudsat_pia / Ncolumns 

  end subroutine cloudsat_precipOccurence
  
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine load_scale_LUTs(rcfg)
    
    type(radar_cfg), intent(inout) :: rcfg
    logical                        :: LUT_file_exists
    integer                        :: i,j,k,ind
    
    ! Load scale LUT from file 
    inquire(file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
         exist=LUT_file_exists)
    
    if(.not.LUT_file_exists) then  
       write(*,*) '*************************************************'
       write(*,*) 'Warning: Could NOT FIND radar LUT file: ', &
            trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'        
       write(*,*) 'Will calculated LUT values as needed'
       write(*,*) '*************************************************'
       return
    else
       OPEN(unit=12,file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
            form='unformatted', &
            err= 89, &
            access='DIRECT',&
            recl=28)
       write(*,*) 'Loading radar LUT file: ', &
            trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
       
       do i=1,maxhclass
          do j=1,mt_ntt
             do k=1,nRe_types
                ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
                read(12,rec=ind) rcfg%Z_scale_flag(i,j,k), &
                     rcfg%Ze_scaled(i,j,k), &
                     rcfg%Zr_scaled(i,j,k), &
                     rcfg%kr_scaled(i,j,k)
             enddo
          enddo
       enddo
       close(unit=12)
       return 
    endif
    
89  write(*,*) 'Error: Found but could NOT READ radar LUT file: ', &
         trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
  end subroutine load_scale_LUTs
  
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine save_scale_LUTs(rcfg)
    type(radar_cfg), intent(inout) :: rcfg
    logical                        :: LUT_file_exists
    integer                        :: i,j,k,ind
    
    inquire(file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
         exist=LUT_file_exists)
    
    OPEN(unit=12,file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
         form='unformatted',err= 99,access='DIRECT',recl=28)
    
    write(*,*) 'Creating or Updating radar LUT file: ', &
         trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
    do i=1,maxhclass
       do j=1,mt_ntt
          do k=1,nRe_types
             ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
             if(.not.LUT_file_exists .or. rcfg%Z_scale_added_flag(i,j,k)) then
                rcfg%Z_scale_added_flag(i,j,k)=.false.
                write(12,rec=ind) rcfg%Z_scale_flag(i,j,k), &
                     rcfg%Ze_scaled(i,j,k), &
                     rcfg%Zr_scaled(i,j,k), &
                     rcfg%kr_scaled(i,j,k)
             endif
          enddo
       enddo
    enddo
    close(unit=12)
    return 
    
99  write(*,*) 'Error: Unable to create/update radar LUT file: ', &
         trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    return  
    
  end subroutine save_scale_LUTs
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine quickbeam_init()

    
  end subroutine quickBeam_init
  ! ##############################################################################################
  ! ##############################################################################################


end module quickbeam


