!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Grell-Freitas deep convection scheme
!
! Reference: Grell and Freitas (2014, ACP)
!
! Author: Jihyeon Jang (NCAR/MMM)

! December 2017 implemented and modified by J. Jang from WRF V3.9 to CAM physics
!               only available with MPAS dynamical core
!---------------------------------------------------------------------------------

module gf_conv_intr

   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use gf_conv,         only: gf_convr=>cup_gf,neg_check,autoconv,aeroevap
   use perf_mod
   use cam_logfile,  only: iulog
   use constituents, only: cnst_add
   use phys_control,      only : phys_getopts
   use pmgrid, only : MPAS_SPHERE_RAD  !MPAS

   implicit none
   private
   save

   ! Public methods

   public ::&
      gf_conv_register,           &! register fields in physics buffer
      gf_conv_readnl,             &! read namelist
      gf_conv_init,               &! initialize
      gf_conv_tend                 ! return tendencies

   character(len=16) :: shallow_scheme

   integer ::& ! indices for fields in the physics buffer
      zm_jt_idx,      &
      zm_maxg_idx,    &
      zm_ideep_idx,   &
      rthften_idx,    &
      rthblten_idx,   &
      rthraten_idx,   &
      rqvften_idx,    &
      rqvblten_idx,   &
      rthcuten_idx,   &
      rqvcuten_idx,   &
      rqccuten_idx,   &
      rqicuten_idx,   &
      rucuten_idx,    &
      rvcuten_idx,    &
      prec_dp_idx,    &
      snow_dp_idx,    &
      prec_sh_idx

   real(r8), parameter :: unset_r8 = huge(1.0_r8)

!  indices for fields in the physics buffer
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    nevapr_dpcu_idx  = 0    
   integer  ::    cmfmc_sh_idx     = 0     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This convective parameterization is build to attempt     !
!     a smooth transition to cloud resolving scales as proposed!
!     by Arakawa et al (2011, ACP). It currently does not use  !
!     subsidencespreading as in G3. Difference and details     !
!     will be described in a paper by                          !
!     Grell and Freitas (2014). The parameterization also      !
!     offers options to couple with aerosols. Both, the smooth !
!     transition part as well as the aerosol coupling are      !
!     experimental. While the smooth transition part is turned !
!     on, nd has been tested dow to a resolution of about 3km  !
!     the aerosol coupling is turned off.                      !
!     More clean-up as well as a direct coupling to chemistry  !
!     will follow for V3.5.1                                   !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS

subroutine gf_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8, dtype_i4

  implicit none

  integer idx

   ! wg top level index of deep cumulus convection.
   call pbuf_add_field('ZM_JT', 'physpkg', dtype_i4, (/pcols/), zm_jt_idx)

   ! wg gathered values of maxi.
   call pbuf_add_field('ZM_MAXG', 'physpkg', dtype_i4, (/pcols/), zm_maxg_idx)

   ! map gathered points to chunk index
   call pbuf_add_field('ZM_IDEEP', 'physpkg', dtype_i4, (/pcols/), zm_ideep_idx)

! deep gbm cloud liquid water (kg/kg)
!   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)

! deep gbm cloud liquid water (kg/kg)
!   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)

   call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
   call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
   call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)

   call pbuf_add_field('PREC_DP',   'physpkg',dtype_r8,(/pcols/),      prec_dp_idx)
   call pbuf_add_field('SNOW_DP',   'physpkg',dtype_r8,(/pcols/),      snow_dp_idx)

   call pbuf_add_field('rthften', 'physpkg', dtype_r8, (/pcols,pver/), rthften_idx)
   call pbuf_add_field('rqvften', 'physpkg', dtype_r8, (/pcols,pver/), rqvften_idx)

   call pbuf_add_field('rucuten', 'physpkg', dtype_r8, (/pcols,pver/), rucuten_idx)
   call pbuf_add_field('rvcuten', 'physpkg', dtype_r8, (/pcols,pver/), rvcuten_idx)
   call pbuf_add_field('rthcuten', 'physpkg', dtype_r8, (/pcols,pver/), rthcuten_idx)
   call pbuf_add_field('rqvcuten', 'physpkg', dtype_r8, (/pcols,pver/), rqvcuten_idx)
   call pbuf_add_field('rqccuten', 'physpkg', dtype_r8, (/pcols,pver/), rqccuten_idx)
   call pbuf_add_field('rqicuten', 'physpkg', dtype_r8, (/pcols,pver/), rqicuten_idx)

end subroutine gf_conv_register

!=========================================================================================

subroutine gf_conv_readnl(nlfile)

   use cam_abortutils,  only: endrun
   use spmd_utils,      only: masterproc
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input


end subroutine gf_conv_readnl

!=========================================================================================

subroutine gf_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: addfld, add_default, horiz_only
  use ppgrid,         only: pcols, pver
  use gf_conv,        only: gfinit
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces


  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields


!
! Register fields with the output buffer
!

    call addfld ('PRECZ',    horiz_only,   'A', 'm/s','total precipitation from ZM convection')
    call addfld ('ZMDT',     (/ 'lev' /),  'A', 'K/s','T tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDQ',     (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDICE',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud ice tendency - Zhang-McFarlane convection')
    call addfld ('ZMDLIQ',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud liq tendency - Zhang-McFarlane convection')
    call addfld ('EVAPTZM',  (/ 'lev' /),  'A', 'K/s','T tendency - Evaporation/snow prod from Zhang convection')
    call addfld ('EVAPQZM',  (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Evaporation from Zhang-McFarlane moist convection')

    call addfld ('CMFMCDZM', (/ 'ilev' /), 'A', 'kg/m2/s','Convection mass flux from ZM deep ')
    call addfld ('PRECCDZM', horiz_only,   'A', 'm/s','Convective precipitation rate from ZM deep')
    call addfld ('CBMF_DP',  horiz_only,   'A', 'kg/m2/s', 'Cloud base mass flux from deep convection')

    call addfld ('PCONVB',   horiz_only ,  'A', 'Pa'    ,'convection base pressure')
    call addfld ('PCONVT',   horiz_only ,  'A', 'Pa'    ,'convection top  pressure')

    call addfld ('CAPE',     horiz_only,   'A', 'J/kg', 'Convectively available potential energy')
    call addfld ('CLDWRK',   horiz_only,   'A', 'J/kg', 'Cloud Work Function')
    call addfld ('TAUZM',      horiz_only,   'A', 'J/kg', 'CAPE consumption time scale')
    call addfld ('rthften', (/ 'lev' /) , 'A', 'K/s'  ,'Heating by Advection')
    call addfld ('rqvften', (/ 'lev' /) , 'A', 'kg/kg/s'  ,'QV tendency by Advection')
    call addfld ('rucuten', (/ 'lev' /) , 'A', 'm/s/s'  ,'U tendency by GF convection')
    call addfld ('rvcuten', (/ 'lev' /) , 'A', 'm/s/s'  ,'V tendency by GF convection')
    call addfld ('rthcuten', (/ 'lev' /) , 'A', 'K/s'   ,'Heating by GF convection')
    call addfld ('rqvcuten', (/ 'lev' /) , 'A', 'kg/kg/s'   ,'QV tendency by GF convection')
    call addfld ('rqccuten', (/ 'lev' /) , 'A', 'kg/kg/s'   ,'QC tendency by GF convection')
    call addfld ('rqicuten', (/ 'lev' /) , 'A', 'kg/kg/s'   ,'QI tendency by GF convection')

    call addfld ('FREQZM',   horiz_only  , 'A', 'fraction', 'Fractional occurance of ZM convection')

    call addfld ('ZMMU',     (/ 'lev' /),  'A', 'kg/m2/s', 'ZM convection updraft mass flux')
    call addfld ('ZMMD',     (/ 'lev' /),  'A', 'kg/m2/s', 'ZM convection downdraft mass flux')
    call addfld ('CMF_DP', (/ 'ilev' /), 'A', 'kg/m2/s','Convection mass flux from deep convection ')
    call addfld ('DLF_DP', (/ 'lev' /), 'A', 'kg/m2/s',' Detrained liquid water from deep convection ')
    call addfld ('RPRD_DP', (/ 'lev' /), 'A', 'kg/kg/s',' Precipitation production from deep convection ')
    call addfld ('EVAPQ_DP',  (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Evaporation from deep convection')
    call addfld ('XF_ENS1_DP',  (/ 'lev' /),  'A', 'kg/m2/s','mass flux ensemble from deep convection')
    call addfld ('XF_ENS2_DP',  (/ 'lev' /),  'A', 'kg/m2/s','mass flux ensemble from deep convection')
    call addfld ('SIG_DP',  horiz_only,  'A', 'no','sigma from deep convection')

! mid convection
    call addfld ('CMF_MID', (/ 'ilev' /), 'A', 'kg/m2/s','Convection mass flux from mid convection ')
    call addfld ('PREC_MID', horiz_only,   'A', 'm/s','Convective precipitation rate from mid convection')
    call addfld ('CBMF_MID',  horiz_only,   'A', 'kg/m2/s', 'Cloud base mass flux from mid convection')

    call addfld ('PCONVB_MID',   horiz_only ,  'A', 'Pa'    ,'mid convection base pressure')
    call addfld ('PCONVT_MID',   horiz_only ,  'A', 'Pa'    ,'mid convection top  pressure')
    call addfld ('CLDWRK_MID',   horiz_only,   'A', 'J/kg', 'Cloud Work Function from mid convection')
    call addfld ('TAUZM_MID',      horiz_only,   'A', 'J/kg', 'CAPE consumption time scale from mid convection')

    call addfld ('MU_MID',     (/ 'lev' /),  'A', 'kg/m2/s', 'mid convection updraft mass flux')
    call addfld ('MD_MID',     (/ 'lev' /),  'A', 'kg/m2/s', 'mid convection downdraft mass flux')
    call addfld ('DLF_MID',     (/ 'lev' /),  'A', 'kg/m2/s', 'Detrained liquid water from mid convection')
    call addfld ('FREQ_MID',   horiz_only  , 'A', 'fraction', 'Fractional occurance of mid convection')
    call addfld ('RPRD_MID', (/ 'lev' /), 'A', 'kg/kg/s',' Precipitation production from mid convection ')
    call addfld ('EVAPQ_MID',  (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Evaporation from mid convection')
    call addfld ('XF_ENS1_MID',  (/ 'lev' /),  'A', 'kg/m2/s','mass flux ensemble from mid convection')
    call addfld ('XF_ENS2_MID',  (/ 'lev' /),  'A', 'kg/m2/s','mass flux ensemble from mid convection')
    call addfld ('SIG_MID',  horiz_only,  'A', 'no','sigma from mid convection')

    call addfld ('PREC_SH', horiz_only,   'A', 'm/s','Convective precipitation rate from GF shallow convection')
    call addfld ('DLF_SH',     (/ 'lev' /),  'A', 'kg/m2/s', 'Detrained liquid water from GF shallow convection')

    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num)

    if ( history_budget ) then
       call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
       call add_default('ZMDT     ', history_budget_histfile_num, ' ')
       call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
       call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('ZMDICE   ', history_budget_histfile_num, ' ')
    end if
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
!    limcnv = 0   ! null value to check against below
!    if (pref_edge(1) >= 4.e3_r8) then
!       limcnv = 1
!    else
!       do k=1,plev
!          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
!             limcnv = k
!             exit
!          end if
!       end do
!       if ( limcnv == 0 ) limcnv = plevp
!    end if
!
!    if (masterproc) then
!       write(iulog,*)'GF_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
!            ' which is ',pref_edge(limcnv),' pascals'
!    end if

    rthften_idx   = pbuf_get_index('rthften')
    rthblten_idx  = pbuf_get_index('rthblten')
    rthraten_idx  = pbuf_get_index('rthraten')
    rqvften_idx   = pbuf_get_index('rqvften')
    rqvblten_idx  = pbuf_get_index('rqvblten')

    rucuten_idx  = pbuf_get_index('rucuten')
    rvcuten_idx  = pbuf_get_index('rvcuten')
    rthcuten_idx  = pbuf_get_index('rthcuten')
    rqvcuten_idx  = pbuf_get_index('rqvcuten')
    rqccuten_idx  = pbuf_get_index('rqccuten')
    rqicuten_idx  = pbuf_get_index('rqicuten')

!    icwmrdp_idx     = pbuf_get_index('ICWMRDP')
!    rprddp_idx      = pbuf_get_index('RPRDDP')
!    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')

!    prec_dp_idx     = pbuf_get_index('PREC_DP')
!    snow_dp_idx     = pbuf_get_index('SNOW_DP')

    prec_sh_idx     = pbuf_get_index('PREC_SH')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')

end subroutine gf_conv_init
!=========================================================================================

subroutine gf_conv_tend(                                        &
              kpbl_in ,                                         &
              ztodt   ,                                         &
              shf     ,cflx    ,                                &
              mcon    ,                                         &
              dlf     ,rliq    ,                                &
              jctop   ,jcbot   ,                                &
              state   ,ptend_all   ,landfrac,  pbuf )
              !dlf2    ,rliq2 )

! - for CAM, the GF shallow convection is not ready yet (need more tests)
! - removed WRF_DFI_RADAR flag dependency in the WRF

!-------------------------------------------------------------
   use spmd_utils,    only: masterproc

   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init, physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc

   use phys_grid,     only: get_area_all_p
   use time_manager,  only: get_nstep, is_first_step
   use physics_buffer,only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,  only: pcnst, cnst_get_ind
   use physconst,    only: cp=>cpair, g=>gravit, xlv=>latvap, cappa

   use gf_conv_sh, only: cup_gf_sh
!++JHJ
   use physics_buffer, only: pbuf_get_index
!--JHJ

!-------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------
      integer, parameter :: ideep=1
      !integer, parameter :: imid_gf=0 !org in WRF V3.9
      integer, parameter :: imid_gf=1  !gfmid
      integer, parameter :: ichoicem=0  ! 0 1 2 8 11 GG
      integer, parameter :: ichoice_s=0 ! 0 1 2 3
      integer, parameter :: dicycle=1 !- diurnal cycle flag
      integer, parameter :: dicycle_m=0 !- diurnal cycle flag
      real(r8), parameter :: aodccn=0.1_r8
      real(r8), parameter :: p1000mb = 100000._r8

      integer, parameter :: ichoice=0
      !integer, parameter :: ishallow_g3=0 !org => not parameter, set from namelist

      integer, parameter :: irqcten=1
      integer, parameter :: irqiten=1

!     stochastic
      integer, parameter :: spp_conv=0
      real(r8), parameter, dimension(pcols,4) :: pattern_spp_conv=0._r8

!-------------------------------------------------------------

! inout variables

   type(physics_state), intent(in),target   :: state          ! Physics state variables
   type(physics_ptend), intent(out)         :: ptend_all      ! individual parameterization tendencies
   type(physics_buffer_desc), pointer       :: pbuf(:)

   real(r8), intent(in) :: ztodt
   real(r8), intent(in) :: landfrac(pcols)
   real(r8), intent(in) :: shf(pcols)
   real(r8), intent(in) :: cflx(pcols)
   real(r8), intent(in) :: kpbl_in(pcols)

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
!   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
!   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
!   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux
   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out
!
!   real(r8), intent(out) :: dlf2(pcols,pver)   ! for shallow conv
!   real(r8), intent(out) :: rliq2(pcols)       ! for shallow conv

! local variables

   ! shallow convection 
   integer :: ishallow_g3

   ! physics types
   type(physics_ptend),target :: ptend_loc     ! package tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection
   real(r8), pointer, dimension(:,:) :: ql  
   real(r8), pointer, dimension(:,:) :: rprd
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation

   integer,  pointer :: jt(:)      ! (pcols)
   integer,  pointer :: maxg(:)    ! (pcols)
   integer,  pointer :: ideepcon(:)   ! (pcols)
   integer           :: lengath

   real(r8), pointer, dimension(:,:)   :: rthften
   real(r8), pointer, dimension(:,:)   :: rqvften
   real(r8), pointer, dimension(:,:)   :: rthblten
   real(r8), pointer, dimension(:,:)   :: rqvblten
   real(r8), pointer, dimension(:,:)   :: rthraten

   real(r8), pointer, dimension(:,:)   :: rucuten
   real(r8), pointer, dimension(:,:)   :: rvcuten
   real(r8), pointer, dimension(:,:)   :: rthcuten
   real(r8), pointer, dimension(:,:)   :: rqvcuten
   real(r8), pointer, dimension(:,:)   :: rqccuten
   real(r8), pointer, dimension(:,:)   :: rqicuten

   real(r8), pointer, dimension(:)     :: prec_sh      ! total precipitation
   real(r8), pointer, dimension(:,:)   :: cmfmc_sh ! Shallow convective mass flux (pcols,pverp) [ kg/s/m^2 ]

   ! for cam
   integer :: kk,kku
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns

   ! history output fields
   real(r8) :: ftem(pcols,pver)       ! Temporary workspace for outfld variables
   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)
   real(r8) :: cldwrk(pcols)          ! cloud work function
   real(r8) :: tau(pcols)             ! convective adjustment time scale
   real(r8) :: mu_out(pcols,pver)     ! updraft mass flux
   real(r8) :: md_out(pcols,pver)     ! downdraft mass flux
   real(r8) :: mb_out(pcols)          ! cloud base mass flux
   real(r8) :: mcond(pcols,pverp)     ! convective mass flux
   real(r8) :: dlfd(pcols,pver)       ! 
   real(r8) :: dlfd_out(pcols,pver)   ! 
   real(r8) :: rliqd(pcols)           ! 
   real(r8) :: rprdd(pcols,pver)      ! 
   real(r8) :: rprdd_out(pcols,pver)  !
   real(r8) :: evapd(pcols,pver)      ! evap q tendency
   real(r8) :: evapd_out(pcols,pver)  ! evap q tendency
!   real(r8) :: evapt(pcols,pver)     ! evap t tendency

   real(r8) :: xf_ens1_out(pcols,pver)   ! cloud base mass flux for debugging
   real(r8) :: xf_ens2_out(pcols,pver)   ! cloud base mass flux for debugging
   real(r8) :: sig_out(pcols)
   real(r8) :: xfm_ens1_out(pcols,pver)  ! cloud base mass flux for debugging
   real(r8) :: xfm_ens2_out(pcols,pver)  ! cloud base mass flux for debugging
   real(r8) :: sigm_out(pcols)

   ! history output fields - mid-level convection
   integer  :: jtm(pcols), maxgm(pcols)
   real(r8) :: pcontm(pcols), pconbm(pcols), freqzmm(pcols)
   real(r8) :: cldwrkm(pcols)      ! cloud work function
   real(r8) :: taum(pcols)         ! convective adjustment time scale
   real(r8) :: mum_out(pcols,pver) ! updraft mass flux
   real(r8) :: mdm_out(pcols,pver) ! downdraft mass flux
   real(r8) :: mbm_out(pcols)      ! cloud base mass flux
   real(r8) :: mconm(pcols,pverp)  ! convective mass flux
   real(r8) :: dlfm(pcols,pver)      !
   real(r8) :: dlfm_out(pcols,pver)  !
   real(r8) :: rliqm(pcols)          !
   real(r8) :: rprdm(pcols,pver)     !
   real(r8) :: rprdm_out(pcols,pver) !
   real(r8) :: evapm(pcols,pver)     ! evap q tendency - mid conv
   real(r8) :: evapm_out(pcols,pver) ! evap q tendency
!   real(r8) :: evaptm(pcols,pver)   ! evap t tendency - mid conv
   real(r8) :: precm(pcols)          !

   ! history output fields - shallow convection
   real(r8) :: mbs_out(pcols)      ! cloud base mass flux - shallow convection
   real(r8) :: mcons(pcols,pverp)  ! Convective mass flux
   !real(r8) :: dlfs(pcols,pver)   !

   ! for state update
   logical  :: lq(pcnst)

   ! for grid-size dependency
   real(r8) :: area(pcols)

   ! local variables for interface between here and gf_convr
   real(r8) :: ht(pcols), xland(pcols)
   integer  :: kpbl(pcols)
   integer  :: k22_shallow(pcols)
   integer  :: kbcon_shallow(pcols), ktop_shallow(pcols)
   real(r8) :: xmb_shallow(pcols)
   integer  :: ktop_deep(pcols)

   ! stochastic
   real(r8) :: rstochcol(pcols,4)
   real(r8) :: rand_mom(pcols), rand_vmas(pcols)
   real(r8) :: rand_clos(pcols,4)

   ! indexes
   integer  :: its,ite,itf
   integer  :: kts,kte,ktf

! LOCAL VARS in WRF

     real(r8)    :: dt
     real(r8),  dimension( pcols ) :: dx
     real(r8),  dimension( pcols ) :: phis
!     real(r8),  dimension( ims:ime , kms:kme , jms:jme )  ::    & !org
     real(r8),  dimension( pcols, pver )                  ::    &
                                                          u,    &
                                                          v,    &
                                                          w,    &
                                                         pi,    &
                                                          t,    &
                                                          q,    &
                                                          p,    &
                                                         zl,    &
                                                       dz8w,    &
!                                                        p8w,    &
                                                        rho

     real(r8),  dimension( pcols, pverp )                 ::    &
                                                        p8w,    &
                                                        zi

     real(r8),  dimension( pcols, pver )                  ::    &
                                                rthften_b2t,    & !bottom to top
                                                rqvften_b2t,    &
                                               rthblten_b2t,    &
                                               rqvblten_b2t,    &
                                               rthraten_b2t

!     real(r8),  dimension( ims:ime , kms:kme , jms:jme ),       &
!                optional                                  ::    &
!               GDC,GDC2

     real(r8),    dimension (pcols,pver) :: &
        dhdt
     real(r8),    dimension (pcols,pver) :: &
        OUTT,OUTQ,OUTQC,cupclw,outu,outv,cnvwt
     real(r8),    dimension (pcols,pver) :: &
        OUTTs,OUTQs,OUTQCs,cupclws,outus,outvs,cnvwts
     real(r8),    dimension (pcols,pver) :: &
        OUTTm,OUTQm,OUTQCm,cupclwm,outum,outvm,cnvwtm
     real(r8),    dimension (pcols)         ::                   &
        pret, prets,pretm,ter11, aa0, xlandi
     real(r8),    dimension (pcols)         ::                   &
        aa1
     real(r8),    dimension (pcols)         ::                   &
        hfxi,qfxi,dxi
     real(r8),    dimension (pcols)         ::                   &
        raincv,pratec
!+lxz
     integer, dimension (pcols) ::                               &
        ierr,ierrs,ierrm
     integer, dimension (pcols) ::                               &
        kbcon, kbcons, kbconm,                                   &
        ktop, ktops, ktopm,                                      &
        kpbli, k22, k22s, k22m
!.lxz
     integer :: ibegc,iendc,jbegc,jendc

     integer, dimension (pcols)         :: jmin,jminm

!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs)
!
     real(r8),    dimension (pcols,pver) :: &
        zo,T2d,q2d,PO,P2d,US,VS,rhoi,tn,qo,tshall,qshall
! output from cup routines, can be used for diagnostics
     real(r8),    dimension (pcols,pver) :: &
        zus,zum,zu,zdm,zd
     real(r8),    dimension (pcols,pver) :: &
        omeg
     real(r8), dimension (pcols)            ::                    &
        ccn,Z1,PSUR,cuten,cutens,cutenm,                        &
        umean,vmean,pmean,xmb,xmbs,                             &
        xmbm,xmb_out,tau_ecmwf_out,xmb_dumm
     real(r8), dimension (pcols)     ::                    &
        edt,edtm,mconv
     real(r8), dimension (pcols)            ::                    &
        tau_ecmwf_outm,aa1m

     integer :: i,j,k,ICLDCK,ipr,jpr,n
     real(r8) :: tcrit,dp,dq
     integer :: iss,jss,nbegin,nend
     real(r8)    :: rkbcon,rktop        !-lxz
     character*50 :: ierrc(pcols)
     character*50 :: ierrcs(pcols)
     character*50 :: ierrcm(pcols)

     real(r8),    dimension (pcols,pver) :: hco,hcdo,zdo
     real(r8),    dimension (pcols,10)   :: forcing,forcing2

     integer, dimension (pcols) :: cactiv
     real(r8),    dimension (pcols,pver) ::  qcheck

!++JHJ
   !all ten - v2: dyn tendency computed for one phy time step
     real(r8), pointer, dimension(:,:) :: dtcore
     real(r8), pointer, dimension(:,:) :: dqcore
     integer itim_old, ifld
!--JHJ

!-------------------------------------------------------------------------------

! check shallow scheme

   call phys_getopts( shallow_scheme_out = shallow_scheme )
   if (shallow_scheme .eq. 'GF') then
      ishallow_g3=1
   else 
      ishallow_g3=0
   end if

! initialize

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

!   call physics_state_copy(state,state1)  ! copy state to local state1.

   ftem = 0._r8

   lq(:) = .FALSE.
   lq(1:3) = .TRUE.
   call physics_ptend_init(ptend_loc, state%psetcols, 'gf_convr', ls=.true., lu=.true., lv=.true., lq=lq)
        ! initialize local ptend type
        ! ls:  if true, then fields to support dsdt are allocated
        ! lu:  if true, then fields to support dudt are allocated
        ! lv:  if true, then fields to support dvdt are allocated
        ! lq:  if true, then fields to support dqdt are allocated
!
! Associate pointers with physics buffer fields
!
   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

   call pbuf_get_field(pbuf, snow_dp_idx,     snow )  ! in gf?? no..
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )  ! [kg/m2/s]

!++KSA check the values of these indices have valid values
    if (masterproc) then
    !   write(iulog,*)'GF_CONV_TEND: prec_dp_idx, snow_dp_idx ',snow_dp_idx, prec_dp_idx
    end if
!--KSA

   call pbuf_get_field(pbuf, rthften_idx,     rthften )
   call pbuf_get_field(pbuf, rqvften_idx,     rqvften )
   call pbuf_get_field(pbuf, rthblten_idx,    rthblten )
   call pbuf_get_field(pbuf, rqvblten_idx,    rqvblten )
   call pbuf_get_field(pbuf, rthraten_idx,    rthraten )

   call pbuf_get_field(pbuf, rucuten_idx,     rucuten )
   call pbuf_get_field(pbuf, rvcuten_idx,     rvcuten )
   call pbuf_get_field(pbuf, rthcuten_idx,    rthcuten )
   call pbuf_get_field(pbuf, rqvcuten_idx,    rqvcuten )
   call pbuf_get_field(pbuf, rqccuten_idx,    rqccuten )
   call pbuf_get_field(pbuf, rqicuten_idx,    rqicuten )

   call pbuf_get_field(pbuf, zm_jt_idx,      jt)
   call pbuf_get_field(pbuf, zm_maxg_idx,    maxg)
   call pbuf_get_field(pbuf, zm_ideep_idx,   ideepcon)

   if(ishallow_g3 .eq. 1 )then
   call pbuf_get_field(pbuf, prec_sh_idx,  prec_sh)
   call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)
   end if

   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )

!++JHJ
   !all ten - v2: dyn tendency computed for one phy time step
    itim_old = pbuf_old_tim_idx()
    ifld = pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    ifld = pbuf_get_index('DQCORE')
    call pbuf_get_field(pbuf, ifld, dqcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
!--JHJ

   call t_startf ('gf_convr')

   its=1 ; ite=ncol ; itf=ite
   kts=1 ; kte=pver ; ktf=kte
   ibegc=its  ;  iendc=itf

   tcrit=258._r8
   ipr=0  !639
   !jpr=0 !141
   rand_mom(:ncol)    = 0._r8
   rand_vmas(:ncol)   = 0._r8
   rand_clos(:ncol,:) = 0._r8

   dt=ztodt

!  CAM: compute delx from area

   call get_area_all_p(lchnk, ncol, area)
   do i = its,ite
     dx(i)=sqrt((area(i)*MPAS_SPHERE_RAD*MPAS_SPHERE_RAD))
   enddo

!  variables (top-bottom ==> bottom-top)

   do k = kts,kte
     kk = kte-k+1
     do i = its,ite
       u(i,kk)=state%u(i,k)
       v(i,kk)=state%v(i,k)
       p(i,kk)=state%pmid(i,k) !p  =grid%p_hyd = 0.5*(p_hyd_w(k)+p_hyd_w(k+1))
       t(i,kk)=state%t(i,k)
       q(i,kk)=state%q(i,k,1) !specific humidity
!      q1(i,k,j)=qv3d(i,k,j)/(1.+qv3d(i,k,j))
       zl(i,kk)=state%zm(i,k) 
       pi(i,kk)=( state%pmid(i,k)/p1000mb )**cappa
       !pi(i,kk)=( state%pressure(i,k)/p1000mb )**cappa
       !pi(i,kk)=1._r8/state%exner(i,k)  !state%exner: inverse of exner function
     enddo
   enddo
!
   do k = kts,kte+1
     kk = kte+1-k+1
     do i = its,ite
       p8w(i,kk)=state%pint(i,k)
       zi(i,kk)=state%zi(i,k) 
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       dz8w(i,k)=zi(i,k+1)-zi(i,k)
       rho(i,k)= (p8w(i,k)-p8w(i,k+1))/g/dz8w(i,k)
       if ( dz8w(i,k) .lt. 0._r8 ) then
         write(iulog,*) 'i k dz8w',i,k,dz8w(i,k)
       endif
     enddo
   enddo
!
   do i = its,ite
     ht(i)=state%phis(i)/g
     kpbl(i)=pver-nint(kpbl_in(i))+1
     !write(iulog,*) 'gf b4 i kpblin kpbl',i,kpbl_in(i),kpbl(i)
   enddo
!
   do i = its,ite
     if ( landfrac(i) .gt. 0.5_r8 ) then
       xland(i)=1.0_r8 ! 1 land, 2 water
     else
       xland(i)=2.0_r8 ! 1 land, 2 water
     endif
   enddo
!
   do k = kts,kte
     kk = kte-k+1
     do i = its,ite
       !all ten - v1: original dyn tendency for one dyn time step
       rthften_b2t(i,kk)=state%rthften(i,k) !dyn
       rqvften_b2t(i,kk)=state%rqvften(i,k) !dyn

!++JHJ
       !all ten - v2: dyn tendency computed for one phy time step
!       write(iulog,*) 'JHJ: ik rthf rqvf',i,k,state%rthften(i,k),state%rthften(i,k)
!       write(iulog,*) 'JHJ: ik dtco dqco',i,k,dtcore(i,k),dqcore(i,k)
!       rthften_b2t(i,kk)=dtcore(i,k)/pi(i,kk)
!       rqvften_b2t(i,kk)=dqcore(i,k)
!--JHJ

       rthblten_b2t(i,kk)=rthblten(i,k)      !pbl
       rqvblten_b2t(i,kk)=rqvblten(i,k)      !pbl
       rthraten_b2t(i,kk)=rthraten(i,k)      !rad

       ! dyn ten
       !rthften_b2t(i,kk)=state%rthften(i,k) !dyn
       !rqvften_b2t(i,kk)=state%rqvften(i,k) !dyn
       !rthblten_b2t(i,kk)=0._r8
       !rqvblten_b2t(i,kk)=0._r8
       !rthraten_b2t(i,kk)=0._r8

       ! phy ten
       !rthften_b2t(i,kk)=0._r8
       !rqvften_b2t(i,kk)=0._r8
       !rthblten_b2t(i,kk)=rthblten(i,k)      !pbl
       !rqvblten_b2t(i,kk)=rqvblten(i,k)      !pbl
       !rthraten_b2t(i,kk)=rthraten(i,k)      !rad

       ! no ls ten
       !rthften_b2t(i,kk)=0._r8
       !rqvften_b2t(i,kk)=0._r8
       !rthblten_b2t(i,kk)=0._r8
       !rqvblten_b2t(i,kk)=0._r8
       !rthraten_b2t(i,kk)=0._r8
     enddo
   enddo

   !IF(PRESENT(k22_shallow)) THEN
   if(ishallow_g3 == 1 )then
   do i=its,ite
     k22_shallow(i)=0
     kbcon_shallow(i)=0
     ktop_shallow(i)=0
     xmb_shallow(i)=0._r8
   enddo
   endif
   do i=its,ite
     ktop_deep(i)=0
   enddo
   rstochcol=0.0_r8

     DO I= its,ite
     do k=kts,kte
       rucuten(i,k)=0._r8
       rvcuten(i,k)=0._r8
       rthcuten(i,k)=0._r8
       rqvcuten(i,k)=0._r8
       IF(irqcten.eq.1)rqccuten(i,k)=0._r8
       IF(irqiten.eq.1)rqicuten(i,k)=0._r8
     enddo
     enddo

     DO I= its,itf
! Stochastic
        if (spp_conv==1) then
        do n=1,4
        rstochcol(i,n)= pattern_spp_conv(i,n)
        if (pattern_spp_conv(i,n) .le. -1.0_r8) then
          rstochcol(i,n)= -1.0_r8
        endif
        if (pattern_spp_conv(i,n) .ge.  1.0_r8) then
          rstochcol(i,n)=  1.0_r8
        endif
        enddo
        endif
        ierrc(i)=" "
        ierrcs(i)=" "
        ierrcm(i)=" "
        ierr(i)=0
        ierrs(i)=0
        ierrm(i)=0

        cuten(i)=0._r8
        cutenm(i)=0._r8
        cutens(i)=1._r8
        if(ishallow_g3.eq.0)cutens(i)=0._r8

        kbcon(i)=0
        kbcons(i)=0
        kbconm(i)=0
        ktop(i)=0
        ktops(i)=0
        ktopm(i)=0

        xmb(i)=0._r8
        xmbs(i)=0._r8
        xmbm(i)=0._r8
        xmb_out(i)=0._r8
        xmb_dumm(i)=0._r8

        k22(i)=0
        k22s(i)=0
        k22m(i)=0

        raincv(i)=0._r8
        pratec (i)=0._r8
        xlandi(i)=xland(i)
        hfxi(i)=shf(i)
        qfxi(i)=cflx(i)

        cactiv(i) = 0
        jmin(i) = 0
        jminm(i) = 0
        forcing(i,:)=0._r8
        forcing2(i,:)=0._r8
        tau_ecmwf_out(i) = 0._r8
        tau_ecmwf_outm(i) = 0._r8

        pret(i)=0._r8
        prets(i) = 0._r8
        pretm(i) = 0._r8

        mconv(i)=0._r8
        ccn(i)=150._r8

        aa1(i) = 0._r8
        aa1m(i) = 0._r8
     ENDDO

     do k=kts,kte
     DO I= its,itf
         omeg(i,k)=0._r8
     ENDDO
     ENDDO

     do i=its,itf
       prec(i)  = 0._r8
       snow(i)  = 0._r8  !KSA, initialize snow 
       rliq(i)  = 0._r8
       !rliq2(i)  = 0._r8
       rliqd(i)  = 0._r8
       mb_out(i) = 0._r8
       tau(i) = 0._r8
       cldwrk(i) = 0._r8
       freqzm(i)=0._r8

       jtm(i)   = 0
       maxgm(i) = 0
       precm(i) = 0._r8
       pcontm(i) = 0._r8
       pconbm(i) = 0._r8
       taum(i) = 0._r8
       cldwrkm(i) = 0._r8
       freqzmm(i)=0._r8
       rliqm(i)  = 0._r8
       mbm_out(i) = 0._r8

       mbs_out(i) = 0._r8
     enddo

     do i=its,itf
       sig_out(i) = 0._r8
       sigm_out(i) = 0._r8
     do k=kts,kte
       xf_ens1_out(i,k) = 0._r8
       xf_ens2_out(i,k) = 0._r8
       xfm_ens1_out(i,k) = 0._r8
       xfm_ens2_out(i,k) = 0._r8
     enddo
     enddo

     do k=kts,kte
     do i=its,itf
       dlf(i,k) = 0._r8
       !dlf2(i,k) = 0._r8
       ql(i,k) = 0._r8
       rprd(i,k) = 0._r8
       evapcdp(i,k) = 0._r8

       mu_out(i,k) = 0._r8
       md_out(i,k) = 0._r8
       dlfd(i,k) = 0._r8
       dlfd_out(i,k) = 0._r8
       dlfm(i,k) = 0._r8
       dlfm_out(i,k) = 0._r8
       rprdd(i,k) = 0._r8
       rprdd_out(i,k) = 0._r8
       rprdm(i,k) = 0._r8
       rprdm_out(i,k) = 0._r8
       evapm(i,k) = 0._r8
       evapd(i,k) = 0._r8
       evapm_out(i,k) = 0._r8
       evapd_out(i,k) = 0._r8
!      evapt(i,k) = 0._r8
!      evaptm(i,k) = 0._r8
       mum_out(i,k) = 0._r8
       mdm_out(i,k) = 0._r8
     enddo
     enddo

     do k=kts,kte+1
     do i=its,itf
        mcon(i,k) = 0._r8
        mcons(i,k) = 0._r8
        mcond(i,k) = 0._r8
        mconm(i,k) = 0._r8
     enddo
     enddo

!ipr= 33 !78
!jpr= 17 !110
     DO I=ITS,ITF
         dxi(i)=dx(i)
         PSUR(I)=p8w(I,1)*.01_r8
!        PSUR(I)=p(I,1,J)*.01
         TER11(I)=max(0._r8,HT(i))
         umean(i)=0._r8
         vmean(i)=0._r8
         pmean(i)=0._r8
         kpbli(i)=kpbl(i)
         zo(i,kts)=ter11(i)+.5_r8*dz8w(i,1)
         DO K=kts+1,ktf
         zo(i,k)=zo(i,k-1)+.5_r8*(dz8w(i,k-1)+dz8w(i,k))
         enddo
     ENDDO

     DO K=kts,ktf
     DO I=ITS,ITF
         po(i,k)=p(i,k)*.01_r8
         P2d(I,K)=PO(i,k)
         rhoi(i,k)=rho(i,k)
         US(I,K) =u(i,k)
         VS(I,K) =v(i,k)

!++JHJ
! ++ org (same as wrf)
!         T2d(I,K)=t(i,k)
!         q2d(I,K)=q(i,k)
!         IF(Q2d(I,K).LT.1.E-08_r8)Q2d(I,K)=1.E-08_r8
!         TN(I,K)=t2d(i,k)+(rthften_b2t(i,k)+rthraten_b2t(i,k)+rthblten_b2t(i,k)) &
!                          *pi(i,k)*dt
!         QO(I,K)=q2d(i,k)+(rqvften_b2t(i,k)+rqvblten_b2t(i,k))*dt
!         TSHALL(I,K)=t2d(i,k)+rthblten_b2t(i,k)*pi(i,k)*dt
!         DHDT(I,K)=cp*rthblten_b2t(i,k)*pi(i,k)+ XLV*rqvblten_b2t(i,k)
!         QSHALL(I,K)=q2d(i,k)+rqvblten_b2t(i,k)*dt
! -- org
!
         TN(I,K)=t(i,k)
         QO(I,K)=q(i,k)
         IF(QO(I,K).LT.1.E-08_r8)QO(I,K)=1.E-08_r8
!         t2d(I,K)=t(i,k)-(rthften_b2t(i,k)+rthraten_b2t(i,k)+rthblten_b2t(i,k)) &
!                          *pi(i,k)*dt
!         q2d(I,K)=q(i,k)-(rqvften_b2t(i,k)+rqvblten_b2t(i,k))*dt

! for midlevel convection (use only BL forcing)

         if( imid_gf == 1 .or. ishallow_g3 == 1 )then
         t2d(I,K)=tn(i,k)-(rthblten_b2t(i,k))*pi(i,k)*dt
         q2d(I,K)=qo(i,k)-(rqvblten_b2t(i,k))*dt
         IF(Q2d(I,K).LT.1.E-08_r8)Q2d(I,K)=1.E-08_r8
         TSHALL(I,K)=t(i,k)
         QSHALL(I,K)=q(i,k)
         end if
         DHDT(I,K)=cp*rthblten_b2t(i,k)*pi(i,k)+ XLV*rqvblten_b2t(i,k)
!         IF(TN(I,K).LT.200._r8)TN(I,K)=T2d(I,K)
!         IF(QO(I,K).LT.1.E-08_r8)QO(I,K)=1.E-08_r8
!--JHJ
         OUTT(I,K)=0._r8
         OUTu(I,K)=0._r8
         OUTv(I,K)=0._r8
         OUTQ(I,K)=0._r8
         OUTQC(I,K)=0._r8
         OUTTm(I,K)=0._r8
         OUTum(I,K)=0._r8
         OUTvm(I,K)=0._r8
         OUTQm(I,K)=0._r8
         OUTQCm(I,K)=0._r8
         OUTTs(I,K)=0._r8
         OUTus(I,K)=0._r8
         OUTvs(I,K)=0._r8
         OUTQs(I,K)=0._r8
         OUTQCs(I,K)=0._r8
         cupclws(i,k) = 0._r8
         cupclw(i,k) = 0._r8
         cupclwm(i,k) = 0._r8
         qcheck(i,k) = 0._r8
     ENDDO
     ENDDO

! for WRf NMM CORE
!#if (NMM_CORE==1)
!! for NMM, tendencies have already been added to T,Q, and total tendencies
!! are stored in *FTEN variables
!     DO K=kts,ktf
!     DO I=ITS,ITF
!         TN(I,K)=t2d(i,k) + RTHFTEN(i,k,j)*pi(i,k,j)*dt
!         QO(I,K)=q2d(i,k) + RQVFTEN(i,k,j)*dt
!         IF(TN(I,K).LT.200.)TN(I,K)=T2d(I,K)
!         IF(QO(I,K).LT.1.E-08)QO(I,K)=1.E-08
!     ENDDO
!     ENDDO
!#endif

     DO K=kts,ktf
     DO I=ITS,ITF
       kk = ktf-k+1
       omeg(i,kk)=state%omega(i,k)
     enddo
     enddo
     do k=  kts+1,ktf-1
     DO I = its,itf
         if((p2d(i,1)-p2d(i,k)).gt.150._r8 .and. p2d(i,k).gt.300._r8)then
            dp=-.5_r8*(p2d(i,k+1)-p2d(i,k-1))
            umean(i)=umean(i)+us(i,k)*dp
            vmean(i)=vmean(i)+vs(i,k)*dp
            pmean(i)=pmean(i)+dp
         endif
     enddo
     enddo
      DO K=kts,ktf-1
      DO I = its,itf
        dq=(q2d(i,k+1)-q2d(i,k))
        mconv(i)=mconv(i)+omeg(i,k)*dq/g
      enddo
      ENDDO
      DO I = its,itf
        if(mconv(i).lt.0._r8)mconv(i)=0._r8
      ENDDO

!-------------------------------------------------------------------------------
!
!---- SHALLOW CONVECTION PARAMETERIZATION
!
       if(ishallow_g3 == 1 )then

          call CUP_gf_sh (                                              &
! input variables, must be supplied
              zo,t2d,q2d,ter11,tshall,qshall,p2d,psur,dhdt,kpbli,      &
              rhoi,hfxi,qfxi,xlandi,ichoice_s,tcrit,dt,                  &
! input variables. Ierr should be initialized to zero or larger than zero for
! turning off shallow convection for grid points
              zus,xmbs,kbcons,ktops,k22s,ierrs,ierrcs,    &
! output tendencies
              outts,outqs,outqcs,cnvwt,prets,cupclws,             &
! shallow inout variables
              !dlfs, &
! dimesnional variables
              itf,ktf,its,ite, kts,kte,ipr)
          do i=its,itf
           if(xmbs(i).le.0._r8)cutens(i)=0._r8
          enddo
          CALL neg_check('shallow',ipr,dt,q2d,outqs,outts,outus,outvs,   &
                                 outqcs,prets,its,ite,kts,kte,itf,ktf)

        endif

!-------------------------------------------------------------------------------
!
!---- MID-LEVEL CONVECTION PARAMETERIZATION
!
   if(imid_gf == 1)then

      call gf_convr(        &
               itf,ktf,its,ite, kts,kte  &

              ,dicycle_m     &
              ,ichoicem      &
              ,ipr           &
              ,ccn           &
              ,dt            &
              ,imid_gf       &

              ,kpbli         &
              ,dhdt          &
              ,xlandi        &

              ,zo            &
              ,forcing2      &
              ,t2d           &
              ,q2d           &
              ,ter11         &
              ,tshall        &
              ,qshall        &
              ,p2d           &
              ,psur          &
              ,us            &
              ,vs            &
              ,rhoi          &
              ,hfxi          &
              ,qfxi          &
              ,dxi           &
              ,mconv         &
              ,omeg          &

              ,cactiv        &
              ,cnvwtm        &
              ,zum           &
              ,zdm           &
              ,edtm          &
              ,xmbm          &
              ,xmb_dumm      &
              ,xmbs          &
              ,pretm         &
              ,outum         &
              ,outvm         &
              ,outtm         &
              ,outqm         &
              ,outqcm        &
              ,kbconm        &
              ,ktopm         &
              ,cupclwm       &
              ,ierrm         &
              ,ierrcm        &
!    the following should be set to zero if not available
              ,rand_mom      & ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     & ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     & ! for stochastics closures, if temporal and spatial patterns exist
              ,0             & ! flag to what you want perturbed
                               ! 1 = momentum transport 
                               ! 2 = normalized vertical mass flux profile
                               ! 3 = closures
                               ! more is possible, talk to developer or
                               ! implement yourself. pattern is expected to be
                               ! betwee -1 and +1
              ,k22m          &
              ,jminm         &
              ,tau_ecmwf_outm &
              ,aa1m   &
              ,evapm  &
              ,dlfm   &
              ,rprdm  &
              ,xfm_ens1_out(:,1:16) &
              ,xfm_ens2_out(:,1:16) &
              ,sigm_out )

            DO I=its,itf
            DO K=kts,ktf
              qcheck(i,k)=q2d(i,k) +outqs(i,k)*dt
            enddo
            enddo
      CALL neg_check('mid',ipr,dt,qcheck,outqm,outtm,outum,outvm,   &
                     outqcm,pretm,its,ite,kts,kte,itf,ktf)
    endif

!-------------------------------------------------------------------------------
!
!---- DEEP CONVECTION PARAMETERIZATION
!

     DO K=kts,ktf
     DO I=ITS,ITF
         t2d(I,K)=tn(i,k)-(rthften_b2t(i,k)+rthraten_b2t(i,k)+rthblten_b2t(i,k)) &
                          *pi(i,k)*dt
         q2d(I,K)=qo(i,k)-(rqvften_b2t(i,k)+rqvblten_b2t(i,k))*dt
         IF(Q2d(I,K).LT.1.E-08_r8)Q2d(I,K)=1.E-08_r8
      ENDDO
      ENDDO

   if(ideep.eq.1)then

      call gf_convr(        &
               itf,ktf,its,ite, kts,kte  &

              ,dicycle       &
              ,ichoice       &
              ,ipr           &
              ,ccn           &
              ,dt            &
              ,0             &

              ,kpbli         &
              ,dhdt          &
              ,xlandi        &

              ,zo            &
              ,forcing       &
              ,t2d           &
              ,q2d           &
              ,ter11         &
              ,tn            &
              ,qo            &
              ,p2d           &
              ,psur          &
              ,us            &
              ,vs            &
              ,rhoi          &
              ,hfxi          &
              ,qfxi          &
              ,dxi           &
              ,mconv         &
              ,omeg          &

              ,cactiv       &
              ,cnvwt        &
              ,zu           &
              ,zd           &
              ,edt          &
              ,xmb          &
              ,xmbm         &
              ,xmbs         &
              ,pret         &
              ,outu         &
              ,outv         &
              ,outt         &
              ,outq         &
              ,outqc        &
              ,kbcon        &
              ,ktop         &
              ,cupclw       &
              ,ierr         &
              ,ierrc        &
!    the following should be set to zero if not available
              ,rand_mom      & ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     & ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     & ! for stochastics closures, if temporal and spatial patterns exist
              ,0             & ! flag to what you want perturbed
                               ! 1 = momentum transport 
                               ! 2 = normalized vertical mass flux profile
                               ! 3 = closures
                               ! more is possible, talk to developer or
                               ! implement yourself. pattern is expected to be
                               ! betwee -1 and +1
              ,k22           &
              ,jmin          &
              ,tau_ecmwf_out &
              ,aa1    &
              ,evapd  &
              ,dlfd   &
              ,rprdd  & 
              ,xf_ens1_out(:,1:16) &
              ,xf_ens2_out(:,1:16) &
              ,sig_out )


        !jpr=0
        ipr=0
            DO I=its,itf
            DO K=kts,ktf
              qcheck(i,k)=q2d(i,k) +(outqs(i,k)+outqm(i,k))*dt
            enddo
            enddo
      CALL neg_check('deep',ipr,dt,qcheck,outq,outt,outu,outv,   &
                                         outqc,pret,its,ite,kts,kte,itf,ktf)
!
   endif

!-------------------------------------------------------------------------------

   if(ishallow_g3.eq.1)then
        DO I=ibegc,iendc
          xmb_shallow(i)=xmbs(i)
          k22_shallow(i)=k22s(i)
          kbcon_shallow(i)=kbcons(i)
          ktop_shallow(i)=ktops(i)
          if (k22s(i).eq.0)   k22_shallow(i)=pver         !JHJ (add if)
          if (kbcons(i).eq.0) kbcon_shallow(i)=pver       !JHJ (add if)
          if (ktops(i).eq.0)  ktop_shallow(i)=1           !JHJ (add if)
        ENDDO
   endif

   DO I=ibegc,iendc
              cuten(i)=0._r8
              ktop_deep(i) = ktop(i)
              if(pret(i).gt.0._r8)then
                 cuten(i)=1._r8
                 freqzm(i)=1.0_r8
              else
                 cuten(i)=0._r8
                 kbcon(i)=0
                 ktop(i)=0
              endif
              if(pretm(i).gt.0._r8)then
                 cutenm(i)=1._r8
                 freqzmm(i)=1.0_r8
              else
                 cutenm(i)=0._r8
                 kbconm(i)=0
                 ktopm(i)=0
              endif

   ENDDO

!++JHJ (make layer variables not zero)
   DO I=ibegc,iendc
              if (k22m(i).eq.0)   k22m(i)  =pver
              if (kbconm(i).eq.0) kbconm(i)=pver
              if (ktopm(i).eq.0)  ktopm(i) =1
              if (k22(i).eq.0)    k22(i)=pver
              if (kbcon(i).eq.0)  kbcon(i) =pver
              if (ktop(i).eq.0)   ktop(i)  =1
              if (ktop_deep(i).eq.0) ktop_deep(i) = 1
   ENDDO
!--JHJ

   DO I=ibegc,iendc
     DO K=kts,ktf
       kk=ktf-k+1
       RTHCUTEN(I,kk)= (cutens(i)*outts(i,k)+ &
                                 cutenm(i)*outtm(i,k)+ &
                                 cuten(i)* outt(i,k)  )/pi(i,k)
       RQVCUTEN(I,kk)= cuten(i)*outq(i,k)   + &
                                cutens(i)*outqs(i,k)+  &
                                cutenm(i)*outqm(i,k)
       RUCUTEN(I,kk)=outum(i,k)*cutenm(i)+outu(i,k)*cuten(i)
       RVCUTEN(I,kk)=outvm(i,k)*cutenm(i)+outv(i,k)*cuten(i)
     ENDDO
   ENDDO

   DO I=ibegc,iendc
     if(pret(i).gt.0. .or. pretm(i).gt.0. .or. prets(i).gt.0.)then
       pratec(i)=cuten(i)*pret(i)+cutenm(i)*pretm(i)+cutens(i)*prets(i)
       raincv(i)=pratec(i)*dt
       prec(i)=pratec(i)*0.001_r8
       precm(i)=cutenm(i)*pretm(i)*0.001_r8
       rkbcon = kte - kbcon(i) + 1
       rktop  = kte -  ktop(i) + 1
     endif
   ENDDO

   DO I=ibegc,iendc
     !tau(i)=tau_ecmwf_out(i)*cuten(i)
     !cldwrk(i)=aa1(i)*cuten(i)
     !taum(i)=tau_ecmwf_outm(i)*cutenm(i)
     !cldwrkm(i)=aa1m(i)*cutenm(i)
     tau(i)=tau_ecmwf_out(i)
     cldwrk(i)=aa1(i)
     taum(i)=tau_ecmwf_outm(i)
     cldwrkm(i)=aa1m(i)
   ENDDO

   IF(irqcten.eq.1) THEN
     DO K=kts,ktf
       kk=ktf-k+1
       DO I=ibegc,iendc
         !rqccuten(I,kk)=outqcm(i,k)+outqcs(i,k)+outqc(I,K)*cuten(i)
         rqccuten(I,kk)=outqcm(i,k)*cutenm(i)+outqcs(i,k)*cutens(i)+outqc(I,K)*cuten(i)
         ql(i,kk)=cupclwm(i,k)*cutenm(i)+cupclws(i,k)*cutens(i)+CUPCLW(I,K)*cuten(i)
         !IF ( PRESENT( GDC ) ) GDC(I,K)=cupclwm(i,k)+cupclws(i,k)+CUPCLW(I,K)*cuten(i)
         !IF ( PRESENT( GDC2 ) ) GDC2(I,K)=0.
       ENDDO
     ENDDO
   ENDIF

!++JHJ - rqicuten will be updated in macrop_driver_tend
            !IF(PRESENT(rqicuten).AND.PRESENT(rqccuten))THEN
            !IF(irqiten.eq.1 .AND. irqcten.eq.1)THEN
            !    DO K=kts,ktf
            !       kk=ktf-k+1
            !      DO I=ibegc,iendc
            !       if(t2d(i,k).lt.258._r8)then
            !          rqicuten(I,kk)=outqcm(i,k)+outqcs(i,k)+outqc(I,K)*cuten(i)
            !          rqccuten(I,kk)=0._r8
            !          !IF ( PRESENT( GDC2 ) ) GDC2(I,K)=cupclwm(i,k)+cupclws(i,k)+CUPCLW(I,K)*cuten(i)
            !       else
            !          rqicuten(I,kk)=0._r8
            !          rqccuten(I,kk)=outqcm(i,k)+outqcs(i,k)+outqc(I,K)*cuten(i)
            !          !IF ( PRESENT( GDC ) ) GDC(I,K)=cupclwm(i,k)+cupclws(i,k)+CUPCLW(I,K)*cuten(i)
            !       endif
            !    ENDDO
            !    ENDDO
            !ENDIF
!--JHJ

   DO I=ibegc,iendc

     mb_out(i)   =  xmb(i)
     mbm_out(i)  =  xmbm(i)

     if(pret(i).gt.0. .or. pretm(i).gt.0.) then

       DO K=kts,ktf
         kk=ktf-k+1
         mu_out(i,kk)=   zu(i,k)*xmb(i)*cuten(i)         ! mu (+ zu updraft)
         md_out(i,kk)= - zd(i,k)*edt(i)*xmb(i)*cuten(i)  ! md (- zd downdraft)
         mum_out(i,kk)=   zum(i,k)*xmbm(i)*cutenm(i)         ! mu (+ zu updraft)
         mdm_out(i,kk)= - zdm(i,k)*edtm(i)*xmbm(i)*cutenm(i)  ! md (- zd downdraft)
         dlf(i,kk) = dlfd(i,k)*cuten(i)+dlfm(i,k)*cutenm(i)
         dlfd_out(i,kk)= dlfd(i,k)*cuten(i)
         dlfm_out(i,kk)= dlfm(i,k)*cutenm(i)

         rliqd(i) = rliqd(i) + dlfd(i,kk)*state%pdel(i,k)/g*cuten(i)
         rliqm(i) = rliqm(i) + dlfm(i,kk)*state%pdel(i,k)/g*cutenm(i)

         rprd(i,kk) = rprdd(i,k)*cuten(i)+rprdm(i,k)*cutenm(i)
         rprdd_out(i,kk)= rprdd(i,k)*cuten(i)
         rprdm_out(i,kk)= rprdm(i,k)*cutenm(i)
         evapcdp(i,kk) = evapd(i,k)*cuten(i)+evapm(i,k)*cutenm(i)
         evapd_out(i,kk)= evapd(i,k)*cuten(i)
         evapm_out(i,kk)= evapm(i,k)*cutenm(i)
       ENDDO

       DO K=kts,ktf
         kku=ktf+1-k+1
         mcond(i,kku) = (zu(i,k) - zd(i,k)*edt(i))*xmb(i)*cuten(i)
         mconm(i,kku) = (zum(i,k) - zdm(i,k)*edtm(i))*xmbm(i)*cutenm(i)
         mcon(i,kku)  = mcond(i,kku) + mconm(i,kku)
       ENDDO

     endif

   ENDDO

   rliq(:itf) = (rliqd(:itf)+rliqm(:itf)) /1000._r8

   if(ishallow_g3.eq.1)then
     DO I=ibegc,iendc
       !DO K=kts,ktf
       !   kk=ktf-k+1
       !   dlf2(i,kk) = dlfs(i,k)*cutens(i)
       !   rliq2(i) = rliq2(i) + dlf2(i,kk)*state%pdel(i,k)/g*cutens(i)
       !ENDDO
       DO K=kts,ktf
         !kk=ktf-k+1
         !mcons(i,kk)= zus(i,k)*xmbs(i)
         kku=ktf+1-k+1
         mcons(i,kku)= zus(i,k)*xmbs(i)*cutens(i)
       ENDDO
       mbs_out(i)   =  xmbs(i)*cutens(i)
       prec_sh(i)=cutens(i)*prets(i)*0.001_r8 !JHJ
     ENDDO
     cmfmc_sh(:ncol,:) = mcons(:ncol,:)
     !rliq2(:itf) = rliq2(:itf)/1000._r8

     call outfld('CMFMC_SH',cmfmc_sh     , pcols   ,lchnk   ) !JHJ
     call outfld( 'CBMF'   , mbs_out     , pcols   , lchnk )
     call outfld('PREC_SH   ', prec_sh   , pcols, lchnk)
     !call outfld('DLF_SH   ',dlf2 ,  pcols   ,lchnk   )
   endif

   call outfld('ZMMU', mu_out, pcols, lchnk)
   call outfld('ZMMD', md_out, pcols, lchnk)
   call outfld('CMFMCDZM   ',mcon ,  pcols   ,lchnk   )
   call outfld('CBMF_DP'    ,mb_out  ,  pcols   ,lchnk   )
   call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )
   call outfld('CMF_DP   ',mcond ,  pcols   ,lchnk   )
   call outfld('DLF_DP   ',dlfd_out ,  pcols   ,lchnk   )
   call outfld('RPRD_DP   ',rprdd_out ,  pcols   ,lchnk   )
   call outfld('XF_ENS1_DP', xf_ens1_out, pcols   ,lchnk   )
   call outfld('XF_ENS2_DP', xf_ens2_out, pcols   ,lchnk   )
   call outfld('SIG_DP',  sig_out, pcols   ,lchnk   )

   call outfld('MU_MID', mum_out, pcols, lchnk)
   call outfld('MD_MID', mdm_out, pcols, lchnk)
   call outfld('CBMF_MID'  ,mbm_out  ,  pcols   ,lchnk   )
   call outfld('FREQ_MID  ',freqzmm          ,pcols   ,lchnk   )
   call outfld('CMF_MID   ',mconm ,  pcols   ,lchnk   )
   call outfld('DLF_MID   ',dlfm_out ,  pcols   ,lchnk   )
   call outfld('RPRD_MID   ',rprdm_out ,  pcols   ,lchnk   )
   call outfld('XF_ENS1_MID', xfm_ens1_out, pcols   ,lchnk   )
   call outfld('XF_ENS2_MID', xfm_ens2_out, pcols   ,lchnk   )
   call outfld('SIG_MID',  sigm_out, pcols   ,lchnk   )

   mcon(:ncol,:) = mcon(:ncol,:) + mcons(:ncol,:)
!
!  tendency computed from updated variables (bottom-top ==> top-bottom)
!
   do k = kts,kte
     kk = kte-k+1
     do i = its,ite
       ptend_loc%u(i,k) = rucuten(i,k)
       ptend_loc%v(i,k) = rvcuten(i,k)
       ptend_loc%s(i,k) = rthcuten(i,k)*pi(i,kk)*cp  !dry static energy
       ptend_loc%q(i,k,1) = rqvcuten(i,k)
       ptend_loc%q(i,k,ixcldliq) = rqccuten(i,k)
       !ptend_loc%q(i,k,ixcldice) = rqicuten(i,k) !JHJ not update ice (will be updated in macrop)
       !ptend_loc%q(i,k,1) = (q1(i,kk)-state%q(i,k,1))*rdelt
     enddo
   enddo
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cp
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   pcont(:ncol) = state%ps(:ncol)
   pconb(:ncol) = state%ps(:ncol)
   pcontm(:ncol) = state%ps(:ncol)
   pconbm(:ncol) = state%ps(:ncol)

   do i = its,ite
     if(pret(i).gt.0.) then
     jt(i)   = kte-ktop(i)+1
     maxg(i) = kte-k22(i)+1
     jctop(i) = jt(i)
     jcbot(i) = kte-kbcon(i)+1
     if (maxg(i).gt.jt(i)) then
       pcont(i) = state%pmid(i,jt(i))  ! gathered array (or jctop ungathered)
       pconb(i) = state%pmid(i,maxg(i))! gathered array
     endif
     endif

     if(pretm(i).gt.0.) then
     jtm(i)   = kte-ktopm(i)+1
     maxgm(i) = kte-k22m(i)+1
     if (maxgm(i).gt.jtm(i)) then
       pcontm(i) = state%pmid(i,jtm(i))  ! gathered array (or jctop ungathered)
       pconbm(i) = state%pmid(i,maxgm(i))! gathered array
     endif
     endif
   end do

   call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
   call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )
   call outfld('PCONVT_MID  ',pcontm     ,pcols   ,lchnk   )
   call outfld('PCONVB_MID  ',pconbm     ,pcols   ,lchnk   )
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )

   call outfld('EVAPQZM ',evapcdp ,pcols   ,lchnk   )
   call outfld('EVAPQ_DP ',evapd_out ,pcols   ,lchnk   )
   call outfld('EVAPQ_MID ',evapm_out ,pcols   ,lchnk   )

   call t_stopf ('gf_convr')

   call outfld('PRECZ   ', prec   , pcols, lchnk)

!  check for closure.. stability and adjustment time scale
   call outfld('CLDWRK', cldwrk, pcols, lchnk)
   call outfld('TAUZM', tau, pcols, lchnk)        

   call outfld('CLDWRK_MID', cldwrkm, pcols, lchnk)
   call outfld('TAUZM_MID', taum, pcols, lchnk)        
   call outfld('PREC_MID', precm, pcols, lchnk)        

!  check tendency
   call outfld('rthraten',rthraten       ,pcols   ,lchnk   )
   call outfld('rthblten',rthblten       ,pcols   ,lchnk   )
   call outfld('rqvblten',rqvblten       ,pcols   ,lchnk   )
   call outfld('rthcuten',rthcuten       ,pcols   ,lchnk   )
   call outfld('rqvcuten',rqvcuten       ,pcols   ,lchnk   )
   call outfld('rthften',state%rthften       ,pcols   ,lchnk   )
   call outfld('rqvften',state%rqvften       ,pcols   ,lchnk   )

   call physics_ptend_init(ptend_all, state%psetcols, 'gf_conv_tend')

   ! add tendency from this process to tendencies from other processes
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   call physics_ptend_dealloc(ptend_loc)

end subroutine gf_conv_tend

end module gf_conv_intr
