!BOP
!
! !MODULE: stepon -- FV Dynamics specific time-stepping
!
! !INTERFACE:
module stepon

! !USES:
! from cam
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_sys_mod,    only: shr_sys_flush
   use spmd_utils,     only: iam, mpicom
   use constituents,   only: pcnst, cnst_name, cnst_longname
   use cam_abortutils, only: endrun
   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use dyn_comp,       only: dyn_import_t, dyn_export_t
   use perf_mod,       only: t_startf, t_stopf, t_barrierf
   use time_manager,   only: get_step_size
! from SE
   use derivative_mod, only: derivative_t
   use quadrature_mod, only: quadrature_t
   use edge_mod,       only: initEdgeBuffer, edgeVpack, edgeVunpack
   use edgetype_mod,   only: EdgeBuffer_t
   use parallel_mod,   only: par
#ifdef debug_coupling
   use dyn_inic_baroclinic, only: test_func
#endif

   implicit none
   private
   save

!
! !PUBLIC MEMBER FUNCTIONS: 
!
  public stepon_init   ! Initialization
  public stepon_run1    ! run method phase 1
  public stepon_run2    ! run method phase 2
  public stepon_run3    ! run method phase 3
  public stepon_final  ! Finalization

!----------------------------------------------------------------------
!
! !DESCRIPTION: Module for dynamics time-stepping.
!
! !REVISION HISTORY:
!
! 2006.05.31  JPE    Created
!
!EOP
!----------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
  type (derivative_t)   :: deriv           ! derivative struct
  type (quadrature_t)   :: gv,gp           ! quadratures on velocity and pressure grids
  type (EdgeBuffer_t) :: edgebuf              ! edge buffer
!-----------------------------------------------------------------------


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_init --- Time stepping initialization
!
! !INTERFACE:
subroutine stepon_init(dyn_in, dyn_out )
! !USES:
  use dimensions_mod, only: nlev, nc, qsize, ntrac, np, qsize_condensate_loading,ldry_mass_vertical_coordinates
  use cam_history,    only: addfld, add_default, horiz_only
  use cam_history,    only: register_vector_field
  use control_mod,    only: smooth_phis_numcycle,rsplit
  use control_mod,    only: tracer_grid_type, TRACER_GRIDTYPE_FVM, ftype
  use gravity_waves_sources, only: gws_init
  use phys_control,   only: use_gw_front, use_gw_front_igw
  use dimensions_mod, only: fv_nphys,ntrac

!!XXgoldyXX:
#ifdef TEST_GRID_COUPLING
  use phys_grid,              only: get_ncols_p, get_rlon_all_p, get_rlat_all_p
  use ppgrid,                 only: pcols, begchunk, endchunk, pver
  use test_dp_grids,          only: dptest_init
#endif
!!XXgoldyXX:

! !OUTPUT PARAMETERS
!
  type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

  integer  :: m, ie,k,tl_f,i,j,tl_fQdp
  real(r8) :: dp_wet, dp_dry(np,np,nlev), factor, factor_array(np,np,nlev),factor_array_nc(nc,nc,nlev)
!!XXgoldyXX:
#ifdef TEST_GRID_COUPLING
  integer  :: lchnk
  real(r8) :: plats(pcols,begchunk:endchunk)
  real(r8) :: plons(pcols,begchunk:endchunk)
  real(r8) :: cols(begchunk:endchunk)
#endif
!!XXgoldyXX
! !DESCRIPTION:
!
! Allocate data, initialize values, setup grid locations and other
! work needed to prepare the dynamics to run. Return weights and 
! vertical coords to atmosphere calling program.
!
!EOP
!-----------------------------------------------------------------------
!BOC

  if (rsplit==0) call endrun('stepon_init: rsplit==0 not supported in CAM')
  if (qsize_condensate_loading>1.and..not.ldry_mass_vertical_coordinates) &
       call endrun('qsize_condensate_loading>1 requires ldry_mass_vertical_coordinates=.true.')

  ! If using FVM transport, check to make sure we are set up correctly
  ! NB: If generalized mapping is implemented between physics grid and
  !     FVM, this restriction can be lifted --goldy
!  if (tracer_grid_type == TRACER_GRIDTYPE_FVM) then
!    if (fv_nphys /= nc) then
!      call endrun('stepon_init: fv_nphys /= nc')
!    end if
!  end if
  ! This is not done in dyn_init due to a circular dependency issue.
  if(iam < par%nprocs) then
     call initEdgeBuffer(par, edgebuf, dyn_in%elem, (3+pcnst)*nlev,nthreads=1)
     if (use_gw_front .or. use_gw_front_igw) call gws_init(dyn_in%elem)
  end if


  ! fields that are written by the dycore
  ! these calls cant be done in dyn_init() because physics grid
  ! is not initialized at that point if making a restart runs
  !
  call addfld ('VOR', (/ 'lev' /), 'A', '1/s',  'Vorticity',                   gridname='GLL')
  call addfld ('DIV', (/ 'lev' /), 'A', '1/s',  'Divergence',                  gridname='GLL')

  if (smooth_phis_numcycle>0) then
     call addfld ('PHIS_SM',  horiz_only, 'I', 'm2/s2', 'Surface geopotential (smoothed)',                gridname='GLL')
     call addfld ('SGH_SM',   horiz_only, 'I', 'm',     'Standard deviation of orography (smoothed)',     gridname='GLL')
     call addfld ('SGH30_SM', horiz_only, 'I', 'm',     'Standard deviation of 30s orography (smoothed)', gridname='GLL')
  endif

  call addfld ('CONVU   ', (/ 'ilev' /),'A', 'm/s2    ','Zonal component IE->KE conversion term',      gridname='physgrid')
  call addfld ('CONVV   ', (/ 'ilev' /),'A', 'm/s2    ','Meridional component IE->KE conversion term', gridname='physgrid')
  call register_vector_field('CONVU', 'CONVV')
  call addfld ('DIFFU   ', (/ 'ilev' /),'A', 'm/s2    ','U horizontal diffusion',                      gridname='physgrid')
  call addfld ('DIFFV   ', (/ 'ilev' /),'A', 'm/s2    ','V horizontal diffusion',                      gridname='physgrid')
  call register_vector_field('DIFFU', 'DIFFV')
  
  call addfld ('ETADOT', (/ 'ilev' /), 'A', '1/s', 'Vertical (eta) velocity', gridname='physgrid')
  ! 
  ! State fields written from the dynamics for analysis
  ! 
 if (fv_nphys>0) then
    do m = 1, qsize
      call addfld (trim(cnst_name(m))//'_gll',  (/ 'lev' /), 'I', 'kg/kg',   &
           cnst_longname(m), gridname='GLL')
    end do
    do m = 1, ntrac
      call addfld (trim(cnst_name(m))//'_fvm',  (/ 'lev' /), 'I', 'kg/kg',   &
           cnst_longname(m), gridname='physgrid')
    end do
    call addfld ('U_gll'     ,(/ 'lev' /), 'I', 'm/s ','U wind on gll grid',gridname='GLL')
    call addfld ('V_gll'     ,(/ 'lev' /), 'I', 'm/s ','V wind on gll grid',gridname='GLL')
    call addfld ('T_gll'     ,(/ 'lev' /), 'I', 'K '  ,'T on gll grid'     ,gridname='GLL')
    call addfld ('PSDRY_gll' ,horiz_only , 'I', 'Pa ' ,'psdry on gll grid' ,gridname='GLL')
    call addfld ('PS_gll'    ,horiz_only , 'I', 'Pa ' ,'ps on gll grid'    ,gridname='GLL')
  end if

  ! Forcing from physics
  ! FU, FV, other dycores, doc, says "m/s" but I think that is m/s^2
  call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term on GLL grid',     gridname='GLL')
  call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term on GLL grid',gridname='GLL')
  call register_vector_field('FU', 'FV')
  call addfld ('FT',  (/ 'lev' /), 'A', 'K/s', 'Temperature forcing term on GLL grid',gridname='GLL')
  do m = 1, qsize
    call addfld ('F'//trim(cnst_name(m)),  (/ 'lev' /), 'I', 'kg/kg/s',   &
         cnst_longname(m)//' mixing ratio forcing term (q_new-q_old) on GLL grid', gridname='GLL')
  end do
  if (ntrac>0) then
    do m = 1, max(ntrac,qsize)
      call addfld ('F'//trim(cnst_name(m))//'_fvm',  (/ 'lev' /), 'I', 'kg/kg/s',   &
           cnst_longname(m)//' mixing ratio forcing term (q_new-q_old) on fvm grid', gridname='physgrid')
    end do
    call addfld ('FU_fvm',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term on fvm grid'     ,gridname='physgrid')
    call addfld ('FV_fvm',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term on fvm grid',gridname='physgrid')
    call register_vector_field('FU_fvm', 'FV_fvm')
    call addfld ('FT_fvm',  (/ 'lev' /), 'A', 'K/s' , 'Temperature forcing term on fvm grid'    ,gridname='physgrid')
  end if





  ! Fields for generating new initial condition files
  call addfld ('U&IC',   (/ 'lev' /),  'I', 'm/s', 'Zonal wind',              gridname='physgrid' )
  call addfld ('V&IC',   (/ 'lev' /),  'I', 'm/s', 'Meridional wind',         gridname='physgrid' )
  ! Don't need to register U&IC V&IC since we don't interpolate IC files
  call add_default ('U&IC',0, 'I')
  call add_default ('V&IC',0, 'I')

  call addfld ('PS&IC', horiz_only,  'I', 'Pa', 'Surface pressure',gridname='physgrid')
  call addfld ('T&IC',  (/ 'lev' /), 'I', 'K',  'Temperature',     gridname='physgrid')

  call add_default ('PS&IC      ',0, 'I')
  call add_default ('T&IC       ',0, 'I')
  if (fv_nphys>0) then
    call addfld ('Ugll&IC ',(/ 'lev' /), 'I', 'm/s ','Zonal wind', gridname='GLL')
    call addfld ('Vgll&IC ',(/ 'lev' /), 'I', 'm/s ','Meridional wind', gridname='GLL')
    call add_default ('Ugll&IC    ',0, 'I')
    call add_default ('Vgll&IC    ',0, 'I')

    call addfld ('PSgll&IC   ', horiz_only,  'I', 'Pa','Surface pressure', gridname='GLL')
    call addfld ('Tgll&IC    ', (/ 'lev' /), 'I', 'K', 'Temperature', gridname='GLL')

    call add_default ('PSgll&IC      ',0, 'I')
    call add_default ('Tgll&IC       ',0, 'I')
  end if
  do m = 1,pcnst
     call addfld (trim(cnst_name(m))//'&IC', (/ 'lev' /), 'I', 'kg/kg', cnst_longname(m), gridname='physgrid')
  end do
  do m = 1,pcnst
     call add_default(trim(cnst_name(m))//'&IC',0, 'I')
  end do
!!XXgoldyXX:
#ifdef TEST_GRID_COUPLING
  do lchnk = begchunk, endchunk
    cols(lchnk) = get_ncols_p(lchnk)
    call get_rlon_all_p(lchnk, pcols, plons(:, lchnk))
    call get_rlat_all_p(lchnk, pcols, plats(:, lchnk))
  end do
  call dptest_init(dyn_in%elem, pcols, begchunk, endchunk, pver, plats, plons, cols)
#endif
!!XXgoldyXX:
  call addfld ('FQ_dBD',  horiz_only, 'I', 'kg/m2', 'Column integrated water vapor mass change passed to dynamics',gridname='GLL')

  call addfld ('ABS_dPSdt',  horiz_only, 'A', 'Pa/s', 'Absolute surface pressure tendency',gridname='GLL')
!  call addfld ('WV_dBD',  horiz_only, 'I', 'kg/m2', 'Column integrated water vapor mass before dynamics',gridname='GLL')
!  call addfld ('WV_dAR',  horiz_only, 'I', 'kg/m2', 'Column integrated water vapor mass after dynamics',gridname='GLL')
!  call addfld ('TT_dBD',  horiz_only, 'I', 'kg/m2', 'Column integrated test tracer TT_LW mass before dynamics',gridname='GLL')
!  call addfld ('TT_dAR',  horiz_only, 'I', 'kg/m2', 'Column integrated test tracer TT_LW mass after dynamics',gridname='GLL')

  !
  ! energy diagnostics
  !
   call addfld ('WV_PDC',   horiz_only, 'A', 'kg/m2','Total column water vapor lost in physics-dynamics coupling',gridname='GLL')
   call addfld ('WL_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud water lost in physics-dynamics coupling',gridname='GLL')
   call addfld ('WI_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud ice lost in physics-dynamics coupling'  ,gridname='GLL')
   call addfld ('TT_PDC',   horiz_only, 'A', 'kg/m2','Total column test tracer lost in physics-dynamics coupling'  ,gridname='GLL')


   call addfld ('WV_dED',   horiz_only, 'A', 'kg/m2','Total column water vapor end of previous dynamics',gridname='GLL')
   call addfld ('WL_dED',   horiz_only, 'A', 'kg/m2','Total column cloud water end of previous dynamics',gridname='GLL')
   call addfld ('WI_dED',   horiz_only, 'A', 'kg/m2','Total column cloud ice end of previous dynamics',gridname='GLL')
   call addfld ('TT_dED',   horiz_only, 'A', 'kg/m2','Total column test tracer ice end of previous dynamics',gridname='GLL')
   call addfld ('SE_dED',   horiz_only, 'A', 'J/m2','Total column Dry Static Energy end of previous dynamics',gridname='GLL')
   call addfld ('KE_dED',   horiz_only, 'A', 'J/m2','Total column Kinetic Energy end of previous dynamics',gridname='GLL')
   !
   ! state in beginning of nsplit loop
   !
   call addfld ('WV_dAF',   horiz_only, 'A', 'kg/m2', &
        'Total column water vapor from previous remapping or state passed to dynamics',gridname='GLL')
   call addfld ('WL_dAF',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud water from previous remapping or state passed to dynamics',gridname='GLL')
   call addfld ('WI_dAF',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud ice from previous remapping or state passed to dynamics',gridname='GLL')
   call addfld ('TT_dAF',   horiz_only, 'A', 'kg/m2', &
        'Total column test tracer ice from previous remapping or state passed to dynamics',gridname='GLL')
   call addfld ('SE_dAF',   horiz_only, 'A', 'J/m2', &
        'Total column Dry Static Energy from previous remapping or state passed to dynamics',gridname='GLL')
   call addfld ('KE_dAF',   horiz_only, 'A', 'J/m2', &
        'Total column Kinetic Energy from previous remapping or state passed to dynamics',gridname='GLL')
   !
   ! state after applyCAMforcing but without the dry mass correction (must be diagnosed)
   !
   call addfld ('WV_dBM',   horiz_only, 'A', 'kg/m2', &
        'Total column water vapor passed to dynamics before dry mass correction',gridname='GLL')
   call addfld ('WL_dBM',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud water passed to dynamics before dry mass correction',gridname='GLL')
   call addfld ('WI_dBM',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud ice passed to dynamics before dry mass correction',gridname='GLL')
   call addfld ('TT_dBM',   horiz_only, 'A', 'kg/m2',  &
        'Total column test tracer ice passed to dynamics before dry mass correction',gridname='GLL')
   call addfld ('SE_dBM',   horiz_only, 'A', 'J/m2', &
        'Total column Dry Static Energy passed to dynamics before dry mass correction',gridname='GLL')
   call addfld ('KE_dBM',   horiz_only, 'A', 'J/m2', &
        'Total column Kinetic Energy passed to dynamics before dry mass correction',gridname='GLL')

   !
   ! state after applyCAMforcing
   !
   call addfld ('WV_dBD',   horiz_only, 'A', 'kg/m2', &
        'Total column water vapor state before dynamics, after dry mass adjustment',gridname='GLL')
   call addfld ('WL_dBD',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud water state before dynamics, after dry mass adjustment',gridname='GLL')
   call addfld ('WI_dBD',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud ice state before dynamics, after dry mass adjustment',gridname='GLL')
   call addfld ('TT_dBD',   horiz_only, 'A', 'kg/m2', &
        'Total column test tracer ice state before dynamics, after dry mass adjustment',gridname='GLL')
   call addfld ('SE_dBD',   horiz_only, 'A', 'J/m2', &
        'Total column Dry Static Energy state before dynamics, after dry mass adjustment',gridname='GLL')
   call addfld ('KE_dBD',   horiz_only, 'A', 'J/m2', &
        'Total column Kinetic Energy state before dynamics, after dry mass adjustment',gridname='GLL')
   !
   ! state before vertical remapping
   !
   call addfld ('WV_dAD',   horiz_only, 'A', 'kg/m2','Total column water vapor before vertical remapping',gridname='GLL')
   call addfld ('WL_dAD',   horiz_only, 'A', 'kg/m2','Total column cloud water before vertical remapping',gridname='GLL')
   call addfld ('WI_dAD',   horiz_only, 'A', 'kg/m2','Total column cloud ice before vertical remapping',gridname='GLL')
   call addfld ('TT_dAD',   horiz_only, 'A', 'kg/m2','Total column test tracer ice before vertical remapping',gridname='GLL')
   call addfld ('SE_dAD',   horiz_only, 'A', 'J/m2','Total column Dry Static Energy before vertical remapping',gridname='GLL')
   call addfld ('KE_dAD',   horiz_only, 'A', 'J/m2','Total column Kinetic Energy before vertical remapping',gridname='GLL')
   !
   ! state at end of nsplit loop
   !
   call addfld ('WV_dAR',   horiz_only, 'A', 'kg/m2','Total column water vapor after vertical remapping',gridname='GLL')
   call addfld ('WL_dAR',   horiz_only, 'A', 'kg/m2','Total column cloud water after vertical remapping',gridname='GLL')
   call addfld ('WI_dAR',   horiz_only, 'A', 'kg/m2','Total column cloud ice after vertical remapping',gridname='GLL')
   call addfld ('TT_dAR',   horiz_only, 'A', 'kg/m2','Total column test tracer ice after vertical remapping',gridname='GLL')
   call addfld ('SE_dAR',   horiz_only, 'A', 'J/m2','Total column Dry Static Energy after vertical remapping',gridname='GLL')
   call addfld ('KE_dAR',   horiz_only, 'A', 'J/m2','Total column Kinetic Energy after vertical remapping',gridname='GLL')

   call addfld ('WV_dBF',   horiz_only, 'A', 'kg/m2','Total column water vapor state passed to parameterizations',gridname='GLL')
   call addfld ('WL_dBF',   horiz_only, 'A', 'kg/m2','Total column cloud water state passed to parameterizations',gridname='GLL')
   call addfld ('WI_dBF',   horiz_only, 'A', 'kg/m2','Total column cloud ice state passed to parameterizations',gridname='GLL')
   call addfld ('TT_dBF',   horiz_only, 'A', 'kg/m2', &
        'Total column test tracer ice state passed to parameterizations',gridname='GLL')
   call addfld ('SE_dBF',   horiz_only, 'A', &
        'J/m2','Total column Dry Static Energy state passed to parameterizations',gridname='GLL')
   call addfld ('KE_dBF',   horiz_only, 'A', 'J/m2','Total column Kinetic Energy state passed to parameterizations',gridname='GLL')


   call addfld ('WV_dBH',   horiz_only, 'A', 'kg/m2','Total column water vapor state before hypervis',gridname='GLL')
   call addfld ('WL_dBH',   horiz_only, 'A', 'kg/m2','Total column cloud water state before hypervis',gridname='GLL')
   call addfld ('WI_dBH',   horiz_only, 'A', 'kg/m2','Total column cloud icestate before hypervis',gridname='GLL')
   call addfld ('TT_dBH',   horiz_only, 'A', 'kg/m2','Total column test tracer icestate before hypervis',gridname='GLL')
   call addfld ('SE_dBH',   horiz_only, 'A', 'J/m2','Total column Dry Static Energystate before hypervis',gridname='GLL')
   call addfld ('KE_dBH',   horiz_only, 'A', 'J/m2','Total column Kinetic Energystate before hypervis',gridname='GLL')

   call addfld ('WV_dCH',   horiz_only, 'A', 'kg/m2', &
        'Total column water vapor state after hypervis but before adding heating term',gridname='GLL')
   call addfld ('WL_dCH',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud water state after hypervis but before adding heating term',gridname='GLL')
   call addfld ('WI_dCH',   horiz_only, 'A', 'kg/m2', &
        'Total column cloud icestate after hypervis but before adding heating term',gridname='GLL')
   call addfld ('TT_dCH',   horiz_only, 'A', 'kg/m2', &
        'Total column test tracer ice state after hypervis but before adding heating term',gridname='GLL')
   call addfld ('SE_dCH',   horiz_only, 'A', 'J/m2', &
        'Total column Dry Static Energy state after hypervis but before adding heating term',gridname='GLL')
   call addfld ('KE_dCH',   horiz_only, 'A', 'J/m2', &
        'Total column Kinetic Energy state after hypervis but before adding heating term',gridname='GLL')


   call addfld ('WV_dAH',   horiz_only, 'A', 'kg/m2','Total column water vapor state after hypervis',gridname='GLL')
   call addfld ('WL_dAH',   horiz_only, 'A', 'kg/m2','Total column cloud water state after hypervis',gridname='GLL')
   call addfld ('WI_dAH',   horiz_only, 'A', 'kg/m2','Total column cloud icestate after hypervis',gridname='GLL')
   call addfld ('TT_dAH',   horiz_only, 'A', 'kg/m2','Total column test tracer icestate after hypervis',gridname='GLL')
   call addfld ('SE_dAH',   horiz_only, 'A', 'J/m2','Total column Dry Static Energy state after hypervis',gridname='GLL')
   call addfld ('KE_dAH',   horiz_only, 'A', 'J/m2','Total column Kinetic Energy state after hypervis',gridname='GLL')

  if (ntrac>0) then
    call addfld ('dp_fvm' ,(/ 'lev' /), 'I', 'Pa','CSLAM Pressure level thickness', gridname='physgrid')
    call addfld ('dp3d_fvm' ,(/ 'lev' /), 'I', 'Pa','SE pressure level thickness integrated over CSLAM control volumes', &
         gridname='physgrid')
    call addfld ('ps_se_fvm',(/ 'lev' /) , 'I', 'Pa','SE surface pressure integrated over CSLAM control volumes', &
         gridname='physgrid')
    call addfld ('PS_fvm',(/ 'lev' /) , 'I', 'Pa','CSLAM surface pressure', gridname='physgrid')
  end if

#ifdef debug_coupling
    call addfld ('xpdel',(/ 'lev' /), 'I', 'm/s ','DP'    )
    call addfld ('xpdeldry',(/ 'lev' /), 'I', 'm/s ','DPdry')
    call addfld ('xpdeldrygll',(/ 'lev' /), 'I', 'm/s ','DPdry gll',gridname='GLL')
    call addfld ('xrpdel',(/ 'lev' /), 'I', 'm/s ','rDP'    )
    call addfld ('xrpdeldry',(/ 'lev' /), 'I', 'm/s ','rDPdry'    )
    call addfld ('xlnpmiddry',(/ 'lev' /), 'I', 'm/s ','lnDPmiddry'    )
    call addfld ('xlnpintdry',(/ 'lev' /), 'I', 'm/s ','lnDPintdry'    )
    call addfld ('xexner',(/ 'lev' /), 'I', 'm/s ','exner'    )
    call addfld ('xDPdry',(/ 'lev' /), 'I', 'Pa ','DPdry'    )
    call addfld ('xQ',(/ 'lev' /), 'I', 'kg/kg ','Q'    )
    call addfld ('xQI',(/ 'lev' /), 'I', 'kg/kg ','QL'    )
    call addfld ('xQL',(/ 'lev' /), 'I', 'kg/kg ','QI'    )
    call addfld ('xTT',(/ 'lev' /), 'I', 'kg/kg ','tracer TT_UN'    )
    call addfld ('xTTgll',(/ 'lev' /), 'I', 'kg/kg ','tracer TT_UN',gridname='GLL')
    call addfld ('xTT_reverse',(/ 'lev' /), 'I', 'kg/kg ','tracer TT_UN'    )
    call addfld ('xTTgll_reverse',(/ 'lev' /), 'I', 'kg/kg ','tracer TT_UN',gridname='GLL')
    call addfld ('xomega',(/ 'lev' /), 'I', 'Pa/s ','omega'    )
    call addfld ('xpint',(/ 'lev' /), 'I', 'Pa ','pint'    )
    call addfld ('xlnpint',(/ 'lev' /), 'I', 'Pa ','lnpint'    )
    call addfld ('xpmid',(/ 'lev' /), 'I', 'Pa ','pmid'    )
    call addfld ('xlnpmid',(/ 'lev' /), 'I', 'Pa ','lnpmid'    )
    call addfld ('xzi',(/ 'lev' /), 'I', 'Pa ','zi'    )
    call addfld ('xzm',(/ 'lev' /), 'I', 'Pa ','zm'    )
    call addfld ('xs',(/ 'lev' /), 'I', 'Pa ','s'    )
    call addfld ('xphis',horiz_only, 'I', 'Pa ','phis'    )
    call addfld ('xpsdry',horiz_only, 'I', 'Pa ','psdry'    )
    call addfld ('xpsdrygll',horiz_only, 'I', 'Pa ','psdry gll',gridname='GLL')
    call addfld ('xps',horiz_only, 'I', 'Pa ','ps'    )
    call addfld ('xpswetgll',horiz_only, 'I', 'Pa ','psdry gll',gridname='GLL')
    call addfld ('xpintdry',(/ 'lev' /), 'I', 'Pa ','pintdry'    )
    call addfld ('xpmiddry',(/ 'lev' /), 'I', 'Pa ','pmiddry'    )
#endif


end subroutine stepon_init

!-----------------------------------------------------------------------
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run1 -- Phase 1 of dynamics run method.
!
! !INTERFACE:
subroutine stepon_run1( dtime_out, phys_state, phys_tend,               &
                        pbuf2d, dyn_in, dyn_out )
  
  use dp_coupling, only: d_p_coupling
  use time_mod,    only: tstep, phys_tscale       ! dynamics timestep
  use physics_buffer, only : physics_buffer_desc
  use dimensions_mod, only: nlev, nc, qsize, ntrac, np, qsize_condensate_loading,ldry_mass_vertical_coordinates!xxx
  implicit none
!
! !OUTPUT PARAMETERS:
!

   real(r8), intent(out) :: dtime_out   ! Time-step
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   type (physics_buffer_desc), pointer :: pbuf2d(:,:)

!-----------------------------------------------------------------------

   ! NOTE: dtime_out computed here must match formula below
   dtime_out = get_step_size()
   if (phys_tscale/=0) then
      dtime_out=phys_tscale  ! set by user in namelist
   endif

   if(iam < par%nprocs) then
      if(tstep <= 0)  call endrun( 'bad tstep')
      if(dtime_out <= 0)  call endrun( 'bad dtime')
   end if
   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------

   call t_barrierf('sync_d_p_coupling', mpicom)
#ifdef debug_coupling   
   call overwrite_dynamics_state(dyn_out)
#endif
   call write_gll_vars(dyn_out)
   call t_startf('d_p_coupling')
   call d_p_coupling (phys_state, phys_tend,  pbuf2d, dyn_out )
   call t_stopf('d_p_coupling')
   call write_phys_state_fvm(phys_state,dyn_out)
end subroutine stepon_run1

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out )
   use bndry_mod,      only: bndry_exchangeV
   use dimensions_mod, only: nlev, nelemd, np, npsq, qsize, ntrac, nc, fv_nphys
   use dimensions_mod, only: qsize_condensate_loading
   use dp_coupling,    only: p_d_coupling
   use parallel_mod,   only: par
   use dyn_grid,       only: TimeLevel
   
   use time_mod,        only: phys_tscale, TimeLevel_Qdp   !  dynamics typestep
   use control_mod,     only: qsplit, smooth_phis_numcycle
   use hycoef,          only: hyai, hybi, ps0
   use cam_history,     only: outfld, hist_fld_active
   use nctopo_util_mod, only: phisdyn,sghdyn,sgh30dyn

   use prim_advance_mod,  only: calc_tot_energy_dynamics
   use fvm_control_volume_mod, only: n0_fvm
   use cam_history,    only: outfld, hist_fld_active, fieldname_len
#ifdef debug_coupling
   use constituents,        only: cnst_get_ind
#endif
   integer :: ixtt

   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   integer :: kptr, ie, ic, i, j, k, m, tl_f, tl_fQdp, qsize_local
   real(r8) :: rec2dt, dyn_ps0
   real(r8) :: dp,dp_tmp,fq,fc,fq0,qn0, ftmp(npsq,nlev,3), fdpq,fdq
   real(r8) :: dtime
   character(len=fieldname_len) :: tfname

   dtime = get_step_size()

   tl_f = TimeLevel%n0   ! timelevel which was adjusted by physics     
   call TimeLevel_Qdp(TimeLevel, qsplit, tl_fQdp)

   ! copy from phys structures -> dynamics structures
   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   call p_d_coupling(phys_state, phys_tend,  dyn_in,tl_fQdp)
   call t_stopf('p_d_coupling')




   call calc_tot_energy_dynamics(dyn_in%elem,1,nelemd,tl_f,tl_fQdp,'dED')

   if(iam >= par%nprocs) return

   ! NOTE: rec2dt MUST be 1/dtime_out as computed above
   
   rec2dt = 1._r8/dtime
   if (phys_tscale/=0) then
     rec2dt = 1._r8/phys_tscale
   endif

   if (ntrac>0) then
     qsize_local = qsize_condensate_loading
   else
     qsize_local = qsize
   end if
#ifdef debug_coupling
   qsize_local = qsize
#endif
   
   dyn_ps0=ps0
   
   call t_startf('bndry_exchange')
   ! do boundary exchange
   ! for physics on GLL grid, for points with duplicate degrees of freedom, 
   ! putuniquepoints() set one of the element values and set the others to zero,
   ! so do a simple sum (boundary exchange with no weights)
   ! for physics grid, we interpolated into all points, so do weighted average:
   do ie=1,nelemd     
      if (fv_nphys>0) then
         do k=1,nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k,1) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k,1) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k,1) =                          &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k,1) *                     &
                 dyn_in%elem(ie)%spheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k,1) =                            &
                 dyn_in%elem(ie)%derived%FT(:,:,k,1) *                       &
                 dyn_in%elem(ie)%spheremp(:,:)
            do m=1, qsize_local
               dyn_in%elem(ie)%derived%FQ(:,:,k,m,1) =                       &
                    dyn_in%elem(ie)%derived%FQ(:,:,k,m,1) *                  &
                    dyn_in%elem(ie)%spheremp(:,:)             
            enddo
         enddo
      endif
      kptr=0
      call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FM(:,:,:,:,1),2*nlev,kptr,ie)
      kptr=kptr+2*nlev
      
      call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FT(:,:,:,1),nlev,kptr,ie)
      kptr=kptr+nlev
      call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FQ(:,:,:,:,1),nlev*qsize_local,kptr,ie)
   end do
   
   call bndry_exchangeV(par, edgebuf,location='stepon_run2')
   
   do ie=1,nelemd
      kptr=0
      call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FM(:,:,:,:,1),2*nlev,kptr,ie)
      kptr=kptr+2*nlev
      call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FT(:,:,:,1),nlev,kptr,ie)
      kptr=kptr+nlev
      call edgeVunpack(edgebuf, dyn_in%elem(ie)%derived%FQ(:,:,:,:,1),        &
           nlev*qsize_local, kptr, ie)
      if (fv_nphys>0) then
         do k=1,nlev
            dyn_in%elem(ie)%derived%FM(:,:,1,k,1) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,1,k,1) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FM(:,:,2,k,1) =                             &
                 dyn_in%elem(ie)%derived%FM(:,:,2,k,1) *                        &
                 dyn_in%elem(ie)%rspheremp(:,:)
            dyn_in%elem(ie)%derived%FT(:,:,k,1) =                               &
                 dyn_in%elem(ie)%derived%FT(:,:,k,1) *                          &
                 dyn_in%elem(ie)%rspheremp(:,:)
            do m=1, qsize_local
               dyn_in%elem(ie)%derived%FQ(:,:,k,m,1) =                           &
                    dyn_in%elem(ie)%derived%FQ(:,:,k,m,1) *                      &
                    dyn_in%elem(ie)%rspheremp(:,:)
            end do
         end do
      end if
   end do

   
   if (hist_fld_active('FU') .or. hist_fld_active('FV') .or.hist_fld_active('FT')) then
     do ie=1,nelemd
       do k=1,nlev
         do j=1,np
           do i=1,np
#ifdef debug_coupling
             dyn_in%elem(ie)%derived%FM(i,j,1,k,1) = dyn_in%elem(ie)%derived%FM(i,j,1,k,1)-test_func(dyn_in%elem(ie)%spherep(i,j)%lat,dyn_in%elem(ie)%spherep(i,j)%lon,k,12)
             dyn_in%elem(ie)%derived%FM(i,j,2,k,1) = dyn_in%elem(ie)%derived%FM(i,j,2,k,1)-test_func(dyn_in%elem(ie)%spherep(i,j)%lat,dyn_in%elem(ie)%spherep(i,j)%lon,k,13)
             dyn_in%elem(ie)%derived%FT(i,j,k,1) = dyn_in%elem(ie)%derived%FT(i,j,k,1)-test_func(dyn_in%elem(ie)%spherep(i,j)%lat,dyn_in%elem(ie)%spherep(i,j)%lon,k,9)
#endif
             ftmp(i+(j-1)*np,k,1) = dyn_in%elem(ie)%derived%FM(i,j,1,k,1)
             ftmp(i+(j-1)*np,k,2) = dyn_in%elem(ie)%derived%FM(i,j,2,k,1)
             ftmp(i+(j-1)*np,k,3) = dyn_in%elem(ie)%derived%FT(i,j,k,1)
           end do
         end do
       end do
       
       call outfld('FU',ftmp(:,:,1),npsq,ie)
       call outfld('FV',ftmp(:,:,2),npsq,ie)
       call outfld('FT',ftmp(:,:,3),npsq,ie)
     end do
   endif

   do m = 1, qsize_local
     tfname = 'F'//trim(cnst_name(m))
     if (hist_fld_active(tfname)) then
       do ie = 1, nelemd
         do j = 1, np
           do i = 1, np
#ifdef debug_coupling
             if (m>1) then
               do k=1,nlev
                 dyn_in%elem(ie)%derived%FQ(i,j,k,m,1) = dyn_in%elem(ie)%derived%FQ(i,j,k,m,1)/dyn_in%elem(ie)%state%dp3d(i,j,k,tl_f)&
                      -test_func(dyn_in%elem(ie)%spherep(i,j)%lat,dyn_in%elem(ie)%spherep(i,j)%lon,k,m)
               end do
             end if
             ftmp(i+(j-1)*np,:,1) = dyn_in%elem(ie)%derived%FQ(i,j,:,m,1)
#else
             ftmp(i+(j-1)*np,:,1) = dyn_in%elem(ie)%derived%FQ(i,j,:,m,1)!/dyn_in%elem(ie)%state%dp3d(i,j,:,tl_f)
#endif
           end do
         end do
         call outfld(tfname, ftmp(:,:,1), npsq, ie)
       end do
     end if
   end do
   
#ifdef debug_coupling
   do ie = 1, nelemd
     dyn_in%elem(ie)%derived%FQ(:,:,:,:,1) = 0.0D0
     dyn_in%elem(ie)%derived%FM(:,:,:,:,1) = 0.0D0
     dyn_in%elem(ie)%derived%FT(:,:,:,1)   = 0.0D0
   end do
#endif
   !
   ! convert elem(ie)%derived%fq to tendency
   !
   do ie=1,nelemd
     do ic = 1, qsize_local
       !$omp parallel do private(k, j, i)
       do k=1,nlev
         do j=1,np
           do i=1,np  
             dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1)=dyn_in%elem(ie)%derived%FQ(i,j,k,ic,1)*&
                  rec2dt*dyn_in%elem(ie)%state%dp3d(i,j,k,tl_f)
           end do
         end do
       end do
     end do
   end do
   if (ntrac>0) then
     do ie=1,nelemd
       do ic=1,ntrac              
         do k=1,nlev
           do j=1,nc
             do i=1,nc
               dyn_in%fvm(ie)%fc(i,j,k,ic) = dyn_in%fvm(ie)%fc(i,j,k,ic)*&
                    rec2dt*dyn_in%fvm(ie)%dp_fvm(i,j,k,n0_fvm)
             end do
           end do
         end do
       end do
     end do
   end if
   call t_stopf('bndry_exchange')

   ! Most output is done by physics.  We pass to the physics state variables
   ! at timelevel "tl_f".  
   ! we will output dycore variables here to ensure they are always at the same
   ! time as what the physics is writing.  
#if 0
   if (hist_fld_active('VOR')) then
      call compute_zeta_C0(tmp_dyn,dyn_in%elem,hybrid,1,nelemd,tl_f,k)
      do ie=1,nelemd
         do j=1,np
            do i=1,np
               ftmp(i+(j-1)*np,1:pver,1) = tmp_dyn(i,j,1:pver)
            end do
         end do
         call outfld('VOR',ftmp(:,:,1),npsq,ie)
       enddo
   endif
   if (hist_fld_active('DIV')) then
      call compute_div_C0(tmp_dyn,dyn_in%elem,hybrid,1,nelemd,tl_f,k)
      do ie=1,nelemd
         do j=1,np
            do i=1,np
               ftmp(i+(j-1)*np,1:pver,1) = tmp_dyn(i,j,1:pver)
            end do
         end do
         call outfld('DIV',ftmp(:,:,1),npsq,ie)
      enddo
   endif
#endif
   if (smooth_phis_numcycle>0) then
      if (hist_fld_active('PHIS_SM')) then
         do ie=1,nelemd
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,1,1) = phisdyn(i,j,ie)
               end do
            end do
            call outfld('PHIS_SM',ftmp(:,1,1),npsq,ie)
         enddo
      endif
      if (hist_fld_active('SGH_SM')) then
         do ie=1,nelemd
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,1,1) = sghdyn(i,j,ie)
               end do
            end do
            call outfld('SGH_SM',ftmp(:,1,1),npsq,ie)
         enddo
      endif
      if (hist_fld_active('SGH30_SM')) then
         do ie=1,nelemd
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,1,1) = sgh30dyn(i,j,ie)
               end do
            end do
            call outfld('SGH30_SM',ftmp(:,1,1),npsq,ie)
         enddo
      endif
   end if
   
   
#ifdef debug_coupling   
   call cnst_get_ind('TT_UN'   , ixtt    , abort=.false.)
   if (ixtt>0) then
     do ie=1,nelemd
       do k=1,nlev
         do j=1,np
           do i=1,np
             ftmp(i+(j-1)*np,k,1) = dyn_in%elem(ie)%derived%fQ(i,j,k,ixtt,1)/rec2dt!-test_func(dyn_in%elem(ie)%spherep(i,j)%lon,dyn_in%elem(ie)%spherep(i,j)%lat, MIN(k-1,6))
           end do
         end do
       end do
       call outfld('xTTgll_reverse',ftmp(:,:,1),npsq,ie)
     end do
   end if
#endif



   end subroutine stepon_run2
   

subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)
   use fvm_control_volume_mod, only: n0_fvm
   use camsrfexch,     only: cam_out_t     
   use dyn_comp,       only: dyn_run
   use cam_history,    only: outfld, hist_fld_active, fieldname_len
   use dimensions_mod, only: nelemd, npsq, nc, np, qsize, ntrac, nlev,fv_nphys
   use dyn_grid,       only: get_gcol_block_d
   use phys_grid,      only: get_ncols_p, get_gcol_all_p, transpose_block_to_chunk
   use phys_grid,      only: block_to_chunk_send_pters, block_to_chunk_recv_pters
   use ppgrid,         only: pcols, pver
   use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
   use derivative_mod, only: subcell_integration
   use dyn_grid,       only: TimeLevel
   use thread_mod,     only: horz_num_threads

   real(r8), intent(in) :: dtime   ! Time-step
   type(cam_out_t),     intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

   character(len=fieldname_len) :: tfname
   integer                      :: ie, i, j, k, m, ioff, icol, pgcols(pcols)
   integer                      :: lchnk, ncols, idmb1(1), idmb2(1), idmb3(1)
   integer                      :: bpter(npsq,0:pver)
   integer                      :: cpter(pcols,0:pver)
   integer                      :: tl_f
   real(r8), allocatable        :: ftmp(:,:,:)
   real(r8), allocatable        :: ptmp(:,:)
   real(r8), allocatable        :: bbuffer(:), cbuffer(:)
   real(r8), allocatable        :: dp_se_fvm_tmp(:,:,:)
   integer :: rc, num_trac

   tl_f = TimeLevel%n0


   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
   call dyn_run(dyn_out,rc)
   call t_stopf  ('dyn_run')


   if (allocated(ftmp)) then
     deallocate(ftmp)
   end if
!#ifdef xxx
   if (fv_nphys>0) then
     !
     ! output FVM variables
     !     
     allocate(ftmp(nc*nc, nlev, nelemd))
     allocate(ptmp(pcols, pver))
     num_trac = MAX(ntrac,qsize)
     do m = 1, 2*num_trac+7
       if (m.le.num_trac.and.m.le.ntrac) then
         tfname = trim(cnst_name(m))//'_fvm'
       else if (ntrac>0.and.m.le.2*num_trac) then
         tfname = 'F'//trim(cnst_name(m-num_trac))//'_fvm'
       else
         if (m==2*num_trac+1.and.ntrac>0) tfname = 'dp_fvm'
         if (m==2*num_trac+2.and.ntrac>0) tfname = 'dp3d_fvm'
         if (m==2*num_trac+3.and.ntrac>0) tfname = 'PS_fvm'
         if (m==2*num_trac+4.and.ntrac>0) tfname = 'ps_se_fvm'
         if (m==2*num_trac+5) tfname = 'FU_fvm'
         if (m==2*num_trac+6) tfname = 'FV_fvm'
         if (m==2*num_trac+7) tfname = 'FT_fvm'
       end if
       if (hist_fld_active(tfname)) then
         !
         ! prepare output field
         !
         if (m.le.num_trac) then
           !$omp parallel do num_threads(horz_num_threads) private(ie, i, j, ioff)
           do ie = 1, nelemd
             ioff = 1
             do j = 1, nc
               do i = 1, nc
                 ftmp(ioff,:,ie) = dyn_out%fvm(ie)%c(i,j,:,m,n0_fvm)
                 ioff = ioff + 1
               end do
             end do
           end do
         else if (m.le.2*num_trac) then
           !$omp parallel do num_threads(horz_num_threads) private(ie, i, j, ioff)
           do ie = 1, nelemd
             ioff = 1
             do j = 1, nc
               do i = 1, nc
                 ftmp(ioff,:,ie) = dyn_out%fvm(ie)%fc(i,j,:,m)
                 ioff = ioff + 1
               end do
             end do
           end do
         else
           if (m==2*num_trac+1) then
             !$omp parallel do num_threads(horz_num_threads) private(ie, i, j, ioff)
             do ie = 1, nelemd
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dyn_out%fvm(ie)%dp_fvm(i,j,:,n0_fvm)
                   ioff = ioff + 1
                 end do
               end do
             end do
           end if
           if (m==2*num_trac+2) then
             allocate(dp_se_fvm_tmp(nc,nc,nlev))
             do ie = 1, nelemd
               do k = 1, nlev
                 dp_se_fvm_tmp(1:nc,1:nc,k) = &
                      subcell_integration(dyn_out%elem(ie)%state%dp3d(:,:,k,tl_f), np, nc, dyn_out%elem(ie)%metdet) / &
                      dyn_out%fvm(ie)%area_sphere(1:nc,1:nc)              
               end do
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dp_se_fvm_tmp(i,j,:)- dyn_out%fvm(ie)%dp_fvm(i,j,:,n0_fvm)
                   ioff = ioff + 1
                 end do
               end do
             end do
             deallocate(dp_se_fvm_tmp)
           end if
           if (m==2*num_trac+3) then
             do ie = 1, nelemd
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dyn_out%fvm(ie)%psC(i,j)
                   ioff = ioff + 1
                 end do
               end do
             end do
           end if
           if (m==2*num_trac+4) then
             allocate(dp_se_fvm_tmp(nc,nc,1))
             do ie = 1, nelemd
               dp_se_fvm_tmp(1:nc,1:nc,1) = &
                    subcell_integration(dyn_out%elem(ie)%state%ps(:,:,tl_f), np, nc, dyn_out%elem(ie)%metdet) / &
                    dyn_out%fvm(ie)%area_sphere(1:nc,1:nc)              
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dp_se_fvm_tmp(i,j,1)-dyn_out%fvm(ie)%psC(i,j)
                   ioff = ioff + 1
                 end do
               end do
             end do
             deallocate(dp_se_fvm_tmp)
           end if
           if (m==2*num_trac+5) then
             do ie = 1, nelemd
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dyn_out%fvm(ie)%fm(i,j,1,:)
                   ioff = ioff + 1
                 end do
               end do
             end do
           end if
           if (m==2*num_trac+6) then
             do ie = 1, nelemd
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dyn_out%fvm(ie)%fm(i,j,2,:)
                   ioff = ioff + 1
                 end do
               end do
             end do
           end if
           if (m==2*num_trac+7) then
             do ie = 1, nelemd
               ioff = 1
               do j = 1, nc
                 do i = 1, nc
                   ftmp(ioff,:,ie) = dyn_out%fvm(ie)%ft(i,j,:)
                   ioff = ioff + 1
                 end do
               end do
             end do
           end if
           
         end if
         !
         ! done preparing output field - change layout
         !
         if (local_dp_map) then
           do lchnk = begchunk, endchunk
             ncols = get_ncols_p(lchnk)
             call get_gcol_all_p(lchnk, pcols, pgcols)
             ptmp = 0._r8
             do icol = 1, ncols
               call get_gcol_block_d(pgcols(icol), 1, idmb1, idmb2, idmb3)
               ie = idmb3(1)
               ioff = idmb2(1)
               ptmp(icol, :) = ftmp(ioff, :, ie)
             end do
             call outfld(tfname, ptmp, pcols, lchnk)
           end do
         else
           allocate(bbuffer(block_buf_nrecs))
           allocate(cbuffer(chunk_buf_nrecs))
           bbuffer = 0.0_r8
           if(iam .lt. par%nprocs) then
             do ie = 1, nelemd
               if (fv_nphys>0) then
                 call block_to_chunk_send_pters(dyn_out%elem(ie)%GlobalID,&
                      fv_nphys*fv_nphys, pver+1, 1, bpter)
                 ncols = fv_nphys*fv_nphys
               else
                 call block_to_chunk_send_pters(dyn_out%elem(ie)%GlobalID,&
                      npsq, pver+1, 1, bpter)
                 ncols = dyn_out%elem(ie)%idxP%NumUniquePts
               end if
               do icol = 1, ncols
                 do i = 0, pver - 1
                   bbuffer(bpter(icol,i)) = ftmp(icol, i, ie)
                 end do
               end do
             end do
           end if
           call t_barrierf ('sync_blk_to_chk', mpicom)
           call transpose_block_to_chunk(1, bbuffer, cbuffer)
           
           !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, icol, i)
           do lchnk = begchunk, endchunk
             ncols = phys_state(lchnk)%ncol
             call block_to_chunk_recv_pters(lchnk, pcols, pver+1, 1, cpter)
             do icol = 1, ncols
               do i = 1, pver
                 ptmp(icol, i) = cbuffer(cpter(icol,i))
               end do
             end do
             call outfld(tfname, ptmp(:,:), pcols, lchnk)
           end do
           
           deallocate( bbuffer )
           deallocate( cbuffer )
         end if
       end if
     end do
   end if
!#endif
   if (allocated(ftmp)) then
     deallocate(ftmp)
   end if
   if (allocated(ptmp)) then
     deallocate(ptmp)
   end if
   if (allocated(bbuffer)) then
     deallocate(bbuffer)
   end if
   if (allocated(cbuffer)) then
     deallocate(cbuffer)
   end if
   
 end subroutine stepon_run3


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)

! !PARAMETERS:
  ! WARNING: intent(out) here means that pointers in dyn_in and dyn_out
  ! are nullified. Unless this memory is released in some other routine,
  ! this is a memory leak.
  ! These currently seem to point to dyn_grid global data.
  type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
!
! !DESCRIPTION:
!
! Deallocate data needed for dynamics. Finalize any dynamics specific
! files or subroutines.
!
!EOP
!-----------------------------------------------------------------------
!BOC


!EOC
end subroutine stepon_final

  subroutine write_phys_state_fvm(state,dyn_out)
    use cam_history,    only: outfld, hist_fld_active, fieldname_len
    use dimensions_mod, only: qsize,ntrac,nelemd, npsq, np, qsize, nlev,fv_nphys
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use constituents,   only: cnst_get_ind

!    use dyn_inic_baroclinic, only: test_func
!------------------------------Arguments--------------------------------
    type(physics_state), intent(inout) :: state(begchunk:endchunk)
    type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
!---------------------------Local storage-------------------------------
    character(len=fieldname_len) :: tfname
    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k,m                                 ! column, level, tracer indices
    integer :: ixcldice, ixcldliq,ixtt             ! CLDICE and CLDLIQ indices
    integer :: ie
!-----------------------------------------------------------------------
    if (fv_nphys>0) then
#ifdef debug_coupling
      do lchnk=begchunk, endchunk
        do m = 2, max(qsize,ntrac)
          do k = 1, pver
            do i = 1, state(lchnk)%ncol 
              state(lchnk)%Q(i,k,m) = state(lchnk)%Q(i,k,m)-test_func(state(lchnk)%lat(i),state(lchnk)%lon(i),k,m)
            end do
          end do
        end do
      end do
      do lchnk=begchunk, endchunk
        do k = 1, pver
          do i = 1, state(lchnk)%ncol 
            state(lchnk)%U(i,k) = state(lchnk)%U(i,k)-test_func(state(lchnk)%lat(i),state(lchnk)%lon(i),k,12)
            state(lchnk)%V(i,k) = state(lchnk)%V(i,k)-test_func(state(lchnk)%lat(i),state(lchnk)%lon(i),k,13)
          end do
        end do
      end do
#endif
#ifdef xxx
      if (hist_fld_active('U_fvm')) then
        do lchnk=begchunk, endchunk
          call outfld('U_fvm'        ,state(lchnk)%U     , pcols   ,lchnk   )
        end do
      end if
      if (hist_fld_active('V_fvm')) then
        do lchnk=begchunk, endchunk
          call outfld('V_fvm'        ,state(lchnk)%V     , pcols   ,lchnk   )
        end do
      end if
      if (hist_fld_active('T_fvm')) then
        do lchnk=begchunk, endchunk
          call outfld('T_fvm'        ,state(lchnk)%T     , pcols   ,lchnk   )
        end do
      end if
      if (hist_fld_active('PSDRY_fvm')) then
        do lchnk=begchunk, endchunk
          call outfld('PSDRY_fvm'    ,state(lchnk)%PSDRY , pcols   ,lchnk   )
        end do
      end if
      if (hist_fld_active('PS_fvm')) then
        do lchnk=begchunk, endchunk
          call outfld('PS_fvm'    ,state(lchnk)%PS       , pcols   ,lchnk   )
        end do
      end if
      do m = 1, max(qsize,ntrac)
        tfname = trim(cnst_name(m))//'_fvm'
        if (hist_fld_active(tfname)) then
          do lchnk=begchunk, endchunk
            call outfld(tfname  ,state(lchnk)%Q(:,:,m), pcols   ,lchnk   )
          end do
        end if
      end do
#endif

    end if
  end subroutine write_phys_state_fvm

  subroutine write_gll_vars(dyn_out)
    use cam_history,    only: outfld, hist_fld_active, fieldname_len
    use dimensions_mod, only: nelemd, npsq, np, qsize, nlev, fv_nphys, ntrac
    use dyn_grid,       only: TimeLevel

    type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
    
    character(len=fieldname_len)       :: tfname
    integer                            :: ie, i, j, k, m, tl_f,ixtt
    real(r8), allocatable              :: ftmp(:,:,:)

    tl_f = TimeLevel%n0


    if (fv_nphys>0) then
      !
      ! Output tracer fields for analysis of advection schemes
      !
      allocate(ftmp(npsq, nlev, 3))
      do m = 1, qsize
        tfname = trim(cnst_name(m))//'_gll'
        if (hist_fld_active(tfname)) then
          do ie = 1, nelemd
            do j = 1, np
              do i = 1, np
                ftmp(i+(j-1)*np,:,1) = dyn_out%elem(ie)%state%Q(i,j,:,m)
              end do
            end do
            call outfld(tfname, ftmp(:,:,1), npsq, ie)
          end do
        end if
      end do
      
      if (hist_fld_active('U_gll').or.hist_fld_active('V_gll')) then
        do ie = 1, nelemd
          do j = 1, np
            do i = 1, np
              ftmp(i+(j-1)*np,:,1) = dyn_out%elem(ie)%state%v(i,j,1,:,tl_f)
              ftmp(i+(j-1)*np,:,2) = dyn_out%elem(ie)%state%v(i,j,2,:,tl_f)
            end do
          end do
          call outfld('U_gll', ftmp(:,:,1), npsq, ie)
          call outfld('V_gll', ftmp(:,:,2), npsq, ie)
        end do
      end if
      
      if (hist_fld_active('T_gll')) then
        do ie = 1, nelemd
          do j = 1, np
            do i = 1, np
              ftmp(i+(j-1)*np,:,1) = dyn_out%elem(ie)%state%T(i,j,:,tl_f)
            end do
          end do
          call outfld('T_gll', ftmp(:,:,1), npsq, ie)
        end do
      end if
      
      if (hist_fld_active('PSDRY_gll')) then
        do ie = 1, nelemd
          do j = 1, np
            do i = 1, np
              ftmp(i+(j-1)*np,1,1) = dyn_out%elem(ie)%state%ps(i,j,tl_f)
            end do
          end do
          call outfld('PSDRY_gll', ftmp(:,1,1), npsq, ie)
        end do
      end if
      
      if (hist_fld_active('PS_gll')) then
        do ie = 1, nelemd
          do j = 1, np
            do i = 1, np
              ftmp(i+(j-1)*np,1,1) = dyn_out%elem(ie)%state%pswet(i,j)
            end do
          end do
          call outfld('PS_gll', ftmp(:,1,1), npsq, ie)
        end do
      end if
    end if
  end subroutine write_gll_vars
  
#ifdef debug_coupling
  subroutine overwrite_dynamics_state(dyn)
    use dimensions_mod, only: nlev, nelemd, np, npsq, qsize, ntrac, nc, fv_nphys
    use dyn_comp,       only: TimeLevel
    use dyn_inic_baroclinic, only: test_func
    use fvm_control_volume_mod, only: n0_fvm
    implicit none
    type (dyn_export_t), intent(inout) :: dyn! Dynamics export container
    integer                            :: ie, i, j, k, m, tl_f
    
    tl_f = TimeLevel%n0    
    
    do m = 2, qsize
       do ie = 1, nelemd
          do k=1,nlev
             do j = 1, np
                do i = 1, np
                   dyn%elem(ie)%state%Q(i,j,k,m)        = &
                        test_func(dyn%elem(ie)%spherep(i,j)%lat,dyn%elem(ie)%spherep(i,j)%lon,k,m)
                   dyn%elem(ie)%state%Qdp(i,j,k,m,tl_f) = &
                        dyn%elem(ie)%state%Q(i,j,k,m)*dyn%elem(ie)%state%dp3d(i,j,k,tl_f)
                end do
             end do
          end do
       end do
    end do
    if (ntrac>0) then
       do m = 2, qsize
          do ie = 1, nelemd             
             do k=1,nlev
                do j=1,nc
                   do i=1,nc
                      dyn%fvm(ie)%c(i,j,k,m,n0_fvm) = test_func(dyn%fvm(ie)%center_cart(i,j)%lat,dyn%fvm(ie)%center_cart(i,j)%lon,k,m)
                   end do
                end do
             end do
          end do
       end do
    end if
    do ie = 1, nelemd
       do k=1,nlev
          do j = 1, np
             do i = 1, np
                dyn%elem(ie)%state%V(i,j,1,k,tl_f)   = test_func(dyn%elem(ie)%spherep(i,j)%lat,dyn%elem(ie)%spherep(i,j)%lon,k,12)
                dyn%elem(ie)%state%V(i,j,2,k,tl_f)   = test_func(dyn%elem(ie)%spherep(i,j)%lat,dyn%elem(ie)%spherep(i,j)%lon,k,13)
             end do
          end do
       end do
    end do
  end subroutine overwrite_dynamics_state
#endif

  

end module stepon
