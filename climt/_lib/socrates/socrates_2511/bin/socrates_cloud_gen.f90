! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Generate sub-grid cloud columns for MCICA
!
! Method:
!   Significantly adapted from a routine supplied by Environment Canada
!   See Raisanen et al, 2004, Stochastic generation of subgrid-scale 
!   cloudy columns for large-scale models,
!   Quart. J. Roy. Meteorol. Soc., 2004 , 130 , 2047-2068
!
!   Generate a nd_profile GCM column by n_layer cloud field with n_subcol
!   subcolumns in each GCM column. Method and model evaluation is 
!   described in: Raisanen et al, 2004, Stochastic generation of 
!   subgrid-scale cloudy columns for large-scale models,
!   Quart. J. Roy. Meteorol. Soc., 2004 , 130 , 2047-2068
!   Input profiles must be from top of the model to the bottom as is the 
!   output. Note that all the cloudy subcolumns are placed at the "front".
!   For example, if there are N cloudy subcolumns for a GCM column then 
!   subcolumns 1 through N will contain the cloudy subcolumns while the 
!   remaining subcolumns will not.
!   The vertical coordinate used here is pressure, rather than height
!   (as used in the original Raisanen 2004 version).
!   -------------------------------------------------------------------- 
!   This version generates the cloud overlap using 4 methods (ioverlap):
!   0. Maximum-Random overlap as described in Geleyn and Hollingsworth, 1979
!   1. Random overlap
!   2. "Generalized overlap" approach of Hogan and Illingsworth, 2000
!   3. "Generalized overlap" approach of Hogan and Illingsworth, 2000
!      with non-contiguous cloudy layers randomly overlapped
!
!   The cloud can be either horizontally homogeneous (ipph=0) or 
!   horizontally inhomogeneous (ipph=1)
!   The inhomogeneity can be described by either a Gamma distribution 
!   or a Beta distribution
!   The values describing the distribution are read in from a data file
!   by the subroutine read_mcica_data, which is contained in the module
!   def_mcica. To switch between distributions, a different data file 
!   is required.
!------------------------------------------------------------------------------
module socrates_cloud_gen
implicit none
character(len=*), parameter, private :: ModuleName = 'SOCRATES_CLOUD_GEN'

integer, parameter :: m = 86436
integer, parameter :: a = 1093
integer, parameter :: c = 18257

contains

subroutine cloud_gen(nd_layer, cloud_top, n_layer, nd_profile, il1, il2, &
                     n_subcol, n1, n2, &
                     ipph, ioverlap, rand_seed, &
                     rlc_cf, rlc_cw, sigma_qcw, avg_cf, &
                     c_cloud, c_ratio, zf, xcw, &
                     n_subcol_cld, c_sub)

use realtype_rd, only: RealK
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use ereport_mod, only: ereport

implicit none

! Input
  integer, intent(in) :: &
    nd_layer,            & ! Size of layer dimension
    cloud_top,           & ! Top cloudy layer
    n_layer,             & ! Number of layers in GCM
    nd_profile,          & ! Dimension size for GCM columns
    il1,                 & ! Start GCM column
    il2,                 & ! End GCM column
    n_subcol,            & ! Number of sub-columns to generate
    n1, n2,              & ! Dimensions of xcw array
    ipph,                & ! Horizontal homogeneity
    ioverlap,            & ! Vertical overlap
    rand_seed(nd_profile)  ! Seed for generating random numbers

  real(RealK), intent(in) :: &
    zf(nd_profile, nd_layer), &
!     Full-level (layer midpoint) pressure (Pa)
    avg_cf(nd_profile, cloud_top:nd_layer), &
!     Cloud amount for each layer
    c_cloud(nd_profile, cloud_top:nd_layer), &
!     Convective cloud amount for each layer
    c_ratio(nd_profile, cloud_top:nd_layer), &
!     Convective condensate ratio
    sigma_qcw(nd_profile, cloud_top:nd_layer), &
!     Normalized cloud condensate std. dev. (Std. dev. over mean)
    rlc_cf(nd_profile, cloud_top:nd_layer), &
!     Cloud fraction decorrelation scale (Pa)
    rlc_cw(nd_profile, cloud_top:nd_layer), &
!     Cloud condensate decorrelation scale (Pa)
    xcw(n1, n2)
!     Distribution of normalised condensate amount as a function of
!     cumulative probability (n1) and relative standard deviation (n2)

! Output
  integer, intent(out) :: n_subcol_cld(nd_profile)
!     Number of cloudy sub-columns

  real(RealK), intent(inout) :: c_sub(nd_profile, cloud_top:nd_layer, n_subcol)
!     Sub-grid cloud water content

! Local variables
  real(RealK) :: ls_ratio(nd_profile, cloud_top:nd_layer)
!     Large-scale condensate ratio
  real(RealK) :: sigma_ccw(nd_profile, cloud_top:nd_layer)
!     Normalized cloud condensate std. dev. for convective cloud

  integer :: &
    int_x(nd_profile), &
    int_y(nd_profile)
  real(RealK) ::       &
    x(nd_profile),     &
    y(nd_profile),     &
    x1(nd_profile),    &
    y1(nd_profile)
!     Random number vectors

  real(RealK) :: alpha(nd_profile, 0:nd_layer)
!     Fraction of maximum/random cloud overlap
  real(RealK) :: rcorr(nd_profile, 0:nd_layer)
!     Fraction of maximum/random cloud condensate overlap

  integer :: &
    rand_seed_x(il1:il2, n_subcol), &
    rand_seed_y(il1:il2, n_subcol)
!     Random numbers for first level in cloud gen

  integer :: random_dummy(il1:il2)

  integer :: i_loc(nd_profile)
!     Index to place the new subcolumns into the arrays starting from the front

  logical :: l_cld(nd_profile)
!     Flag that cloud has been generated in subcolumn


  integer ::   &
    il,        & ! Counter over GCM columns
    i,         & ! Counter
    k,         & ! Counter over vertical layers
    k_top,     & ! Index of top most cloud layer
    k_base,    & ! Index of lowest cloud layer
    ind1,      & ! Index in variability calculation
    ind2,      & ! Index in variability calculation
    ErrorStatus  ! ErrorStatus for Ereport

  real(RealK) :: &
    rind1,       & ! Real index in variability calculation
    rind2,       & ! Real index in variability calculation
    zcw            ! Ratio of cloud condensate mixing ratio for
                   ! this cell to its layer cloud-mean value

  logical :: &
    l_top,   & ! Flag for cloud top
    l_bot      ! Flag for cloud bottom

  character(len=*), parameter :: RoutineName = 'CLOUD_GEN'
  character(len=80)           :: cmessage

  integer(kind=jpim), parameter :: zhook_in  = 0
  integer(kind=jpim), parameter :: zhook_out = 1
  real(kind=jprb)               :: zhook_handle

  integer :: temp
  integer, parameter :: rand_seed_max = (huge(m)-c)/a
  real(RealK), parameter :: rm = 1.0_RealK/real(m, RealK)
  real(RealK) :: cut = epsilon(cut)
!     Minimum cloud fraction


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Generate random seeds for each sub-column
  do il=il1, il2
    ! Ensure the integers do not overflow
    random_dummy(il) = modulo(abs(rand_seed(il)), rand_seed_max)
  end do
  do i=1, n_subcol
    do il=il1, il2
      temp = random_dummy(il)*a+c
      random_dummy(il) = temp-int(real(temp, RealK)*rm)*m
    end do
  end do
  do i=1, n_subcol
    do il=il1, il2
      temp = random_dummy(il)*a+c
      random_dummy(il) = temp-int(real(temp, RealK)*rm)*m
      rand_seed_x(il, i) = random_dummy(il)
    end do
  end do
  do i=1, n_subcol
    do il=il1, il2
      temp = random_dummy(il)*a+c
      random_dummy(il) = temp-int(real(temp, RealK)*rm)*m
      rand_seed_y(il, i) = random_dummy(il)
    end do
  end do

  ! Set the large-scale condensate ratio to be consistent with c_ratio
  where (c_cloud > 0.0_RealK .and. c_cloud < avg_cf)
    ls_ratio = max(0.0_RealK, (avg_cf - c_ratio*c_cloud)) / (avg_cf - c_cloud)
  elsewhere
    ls_ratio = 1.0_RealK
  end where

  ! Set convective cloud RSD equal to large-scale value for now.
  sigma_ccw=sigma_qcw

  ! Initialize the arrays
  do il = il1, il2
     i_loc(il) = 1
     l_cld(il) = .false.
  end do
  c_sub(il1:il2, :, :) = 0.0_RealK

  ! Find uppermost cloudy layer
  l_top = .false.
  k_top = 0
  do k=cloud_top, n_layer
     do il = il1, il2
        k_top = k
        if (avg_cf(il,k) > cut) l_top = .true.
     end do ! il
     if (l_top) exit
  end do ! k

  ! If no cloudy layers in any GCM column then exit
  if (k_top == 0) then
    if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                            zhook_handle)
    return
  end if

  ! Find lowermost cloudy layer
  l_bot = .false.
  k_base = 0
  do k=n_layer, cloud_top, -1
     do il = il1, il2
        k_base = k
        if (avg_cf(il,k) >= cut) l_bot = .true.
     end do ! il
     if (l_bot) exit
  end do ! k

  select case (ioverlap)
  case (0)
    ! Maximum-random overlap following the definition 
    ! given in Geleyn and Hollingsworth
    do i=1, n_subcol

      ! Generate all subcolumns for latitude chain
      do il = il1, il2
         int_x(il) = rand_seed_x(il,i)
         int_y(il) = rand_seed_y(il,i)
      end do

      do k = k_top, k_base
        do il = il1, il2
          if (avg_cf(il,k) > 0.0_RealK) then

            temp=int_x(il)*a+c
            int_x(il)=temp-int(real(temp, RealK)*rm)*m
            x(il)=real(int_x(il), RealK)*rm

            temp=int_y(il)*a+c
            int_y(il)=temp-int(real(temp, RealK)*rm)*m
            y(il)=real(int_y(il), RealK)*rm

            if (k == cloud_top) then
              temp=int_x(il)*a+c
              int_x(il)=temp-int(real(temp, RealK)*rm)*m
              x(il)=real(int_x(il), RealK)*rm

              temp=int_y(il)*a+c
              int_y(il)=temp-int(real(temp, RealK)*rm)*m
              y(il)=real(int_y(il), RealK)*rm

            else if (x(il) < (1.0_RealK-avg_cf(il,k-1))) then
              ! It is clear above
              temp=int_x(il)*a+c
              int_x(il)=temp-int(real(temp, RealK)*rm)*m
              x1(il)=real(int_x(il), RealK)*rm

              temp=int_y(il)*a+c
              int_y(il)=temp-int(real(temp, RealK)*rm)*m
              y1(il)=real(int_y(il), RealK)*rm

              x(il)=x1(il) * (1.0_RealK - avg_cf(il,k-1))
              y(il)=y1(il)
            end if

            ! Treatment of cloudy cells
            if (x(il) < avg_cf(il,k)) then
              ! Generate cloud here
              if (ipph == 0) then
                ! Homogeneous clouds
                zcw = 1.0_RealK
              else
                ! Horizontally variable clouds:
                ! Determine ZCW = ratio of cloud condensate mixing ratio
                ! QC for this cell to its mean value for all cloudy cells
                ! in this layer.
                ! Use bilinear interpolation of ZCW tabulated in array
                ! XCW as a function of
                !    * cumulative probability Y
                !    * relative standard deviation SIGMA
                rind1 = y(il) * (n1 - 1) + 1.0_RealK
                ind1  = max(1, min(int(rind1), n1-1))
                rind1 = max(0.0_RealK, min(rind1-real(ind1, RealK), 1.0_RealK))
                rind2 = 40.0_RealK * sigma_qcw(il,k) - 3.0_RealK
                ind2  = max(1, min(int(rind2), n2-1))
                rind2 = max(0.0_RealK, min(rind2-real(ind2, RealK), 1.0_RealK))

                zcw = (1.0_RealK-rind1) * (1.0_RealK-rind2)* xcw(ind1,ind2) &
                    + (1.0_RealK-rind1) * rind2 * xcw(ind1,ind2+1) &
                    + rind1 * (1.0_RealK-rind2) * xcw(ind1+1,ind2) &
                    + rind1 * rind2 * xcw(ind1+1,ind2+1)
              end if

              ! A horizontally constant IWC/LWC ratio is assumed
              ! for each layer so far
              l_cld(il) = .true.
              c_sub(il,k,i_loc(il)) = zcw
            end if

          end if
        end do ! il
      end do ! k

      ! Need to check if a cloudy subcolumn was generated
      do il = il1, il2
         if (l_cld(il)) then
            i_loc(il) = i_loc(il) + 1
            l_cld(il) = .false.
         end if
      end do
    end do ! i

  case (1, 2, 3)
    ! ioverlap = 1, random overlap
    ! ioverlap = 2, exponential overlap
    ! ioverlap = 3, exponential-random overlap
    do k = k_top-1, k_base-1
      ! Calculate overlap factors ALPHA for cloud fraction and RCORR for cloud
      ! condensate based on layer midpoint distances and decorrelation depths
      if (ioverlap == 1 .or. k < cloud_top) then
        alpha(:,k) = 0.0_RealK
        rcorr(:,k) = 0.0_RealK
      else
        do il = il1, il2
          if (ioverlap == 2 .or. avg_cf(il,k) > 0.0_RealK) then
            if (rlc_cf(il,k) > 0.0_RealK) then
              alpha(il,k) = exp((zf(il,k) - zf(il,k+1)) / rlc_cf(il,k))
            else
              alpha(il,k) = 0.0_RealK
            end if
            if (rlc_cw(il,k) > 0.0_RealK) then
              rcorr(il,k) = exp((zf(il,k) - zf(il,k+1)) / rlc_cw(il,k))
            else
              rcorr(il,k) = 0.0_RealK
            end if
          else
            ! If ioverlap == 3, non-contiguous cloud layers are
            ! randomly overlapped.
            alpha(il,k) = 0.0_RealK
            rcorr(il,k) = 0.0_RealK
          end if
        end do ! il
      end if
    end do ! k

    do i=1, n_subcol
      ! Generate all subcolumns for latitude chain
      do il = il1, il2
        int_x(il)=rand_seed_x(il,i)
        int_y(il)=rand_seed_y(il,i)
      end do

      do k = k_top, k_base
        do il = il1, il2
          if (avg_cf(il,k) > 0.0_RealK) then

            temp=int_x(il)*a+c
            int_x(il)=temp-int(real(temp, RealK)*rm)*m
            x1(il)=real(int_x(il), RealK)*rm

            temp=int_y(il)*a+c
            int_y(il)=temp-int(real(temp, RealK)*rm)*m
            y1(il)=real(int_y(il), RealK)*rm

            if (x1(il) >= alpha(il,k-1)) then
              temp=int_x(il)*a+c
              int_x(il)=temp-int(real(temp, RealK)*rm)*m
              x(il)=real(int_x(il), RealK)*rm
            end if
            if (y1(il) >= rcorr(il,k-1)) then
              temp=int_y(il)*a+c
              int_y(il)=temp-int(real(temp, RealK)*rm)*m
              y(il)=real(int_y(il), RealK)*rm
            end if

            ! Treatment of cloudy cells
            if (x(il) < avg_cf(il,k)) then
              ! Generate cloud here
              if (ipph == 0) then
                ! Homogeneous clouds
                zcw = 1.0_RealK
              else if (x(il) < c_cloud(il,k)) then
                ! Horizontally variable convective cloud
                rind1 = y(il) * (n1 - 1) + 1.0_RealK
                ind1  = max(1, min(int(rind1), n1-1))
                rind1 = max(0.0_RealK, min(rind1-real(ind1, RealK), 1.0_RealK))
                rind2 = 40.0_RealK * sigma_ccw(il,k) - 3.0_RealK
                ind2  = max(1, min(int(rind2), n2-1))
                rind2 = max(0.0_RealK, min(rind2-real(ind2, RealK), 1.0_RealK))

                zcw = (1.0_RealK-rind1) * (1.0_RealK-rind2)* xcw(ind1,ind2) &
                    + (1.0_RealK-rind1) * rind2 * xcw(ind1,ind2+1) &
                    + rind1 * (1.0_RealK-rind2) * xcw(ind1+1,ind2) &
                    + rind1 * rind2 * xcw(ind1+1,ind2+1)

                zcw = zcw*c_ratio(il,k)
              else
                ! Horizontally variable cloud
                rind1 = y(il) * (n1 - 1) + 1.0_RealK
                ind1  = max(1, min(int(rind1), n1-1))
                rind1 = max(0.0_RealK, min(rind1-real(ind1, RealK), 1.0_RealK))
                rind2 = 40.0_RealK * sigma_qcw(il,k) - 3.0_RealK
                ind2  = max(1, min(int(rind2), n2-1))
                rind2 = max(0.0_RealK, min(rind2-real(ind2, RealK), 1.0_RealK))

                zcw = (1.0_RealK-rind1) * (1.0_RealK-rind2)* xcw(ind1,ind2) &
                    + (1.0_RealK-rind1) * rind2 * xcw(ind1,ind2+1) &
                    + rind1 * (1.0_RealK-rind2) * xcw(ind1+1,ind2) &
                    + rind1 * rind2 * xcw(ind1+1,ind2+1)

                zcw = zcw*ls_ratio(il,k)
              end if

              ! A horizontally constant IWC/LWC ratio is assumed
              ! for each layer so far
              l_cld(il) = .true.
              c_sub(il,k,i_loc(il)) = zcw
            end if

          end if
        end do ! il
      end do ! k

      ! Need to check if a cloudy subcolumn was generated
      do il = il1, il2
        if (l_cld(il)) then
          i_loc(il) = i_loc(il) + 1
          l_cld(il) = .false.
        end if
      end do
    end do ! i

  case default
    cmessage = 'Choice of overlap not compatible with cloud generator'
    ErrorStatus = 1
    call ereport(RoutineName, ErrorStatus, cmessage)

  end select ! ioverlap

  ! Record the number of cloudy subcolumns generated
  do il = il1, il2
    n_subcol_cld(il) = i_loc(il)-1
  end do

  if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine cloud_gen
end module socrates_cloud_gen
