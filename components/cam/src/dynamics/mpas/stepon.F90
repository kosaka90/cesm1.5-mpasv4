!
!MODULE: stepon -- MPAS Dynamics specific time-stepping
!
module stepon

! from cam
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_sys_mod,    only: shr_sys_flush
   use pmgrid,         only: plev, plevp, plat
   use spmd_utils,     only: iam, masterproc
   use constituents,   only: pcnst, cnst_name, cnst_longname
   use cam_abortutils, only: endrun
   use ppgrid,         only: begchunk, endchunk
   use physics_types,  only: physics_state, physics_tend
   use dyn_comp,       only: dyn_import_t, dyn_export_t
#if defined ( SPMD )
   use mpishorthand,   only: mpicom
#endif
   use perf_mod
 
   implicit none

   private

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
! 2010.01.17  AAM    Created from Homme version
!
!EOP
!----------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
!-----------------------------------------------------------------------


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_init --- Time stepping initialization
!
! !INTERFACE:
subroutine stepon_init( dyn_in, dyn_out )
! !USES:
  use cam_history,    only: addfld, add_default, horiz_only
  use cam_history,    only: register_vector_field


! !OUTPUT PARAMETERS
!
  type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container (unused)
  type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container (unused)

  integer :: m
! !DESCRIPTION:
!
! Allocate data, initialize values, setup grid locations and other
! work needed to prepare the dynamics to run. Return weights and 
! vertical coords to atmosphere calling program.
!
!EOP
!-----------------------------------------------------------------------
!BOC


  ! fields that are written by the dycore
  ! these calls cant be done in dyn_init() because physics grid
  ! is not initialized at that point if making a restart runs
  !
  ! Forcing from physics
  ! FU, FV, other dycores, doc, says "m/s" but I think that is m/s^2
  call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term', gridname='physgrid')
  call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term', gridname='physgrid')
  call register_vector_field('FU', 'FV')

  call addfld ('U&IC',   (/ 'lev' /),  'I', 'm/s', 'Zonal wind', gridname='physgrid' )
  call addfld ('V&IC',   (/ 'lev' /),  'I', 'm/s', 'Meridional wind', gridname='physgrid' )
  call add_default ('U&IC       ',0, 'I')
  call add_default ('V&IC       ',0, 'I')

  call addfld ('PS&IC', horiz_only,  'I', 'Pa', 'Surface pressure',gridname='physgrid')
  call addfld ('T&IC',  (/ 'lev' /), 'I', 'K',  'Temperature', gridname='physgrid')
  call add_default ('PS&IC      ',0, 'I')
  call add_default ('T&IC       ',0, 'I')
  do m = 1,pcnst
     call addfld (trim(cnst_name(m))//'&IC', (/ 'lev' /), 'I', 'kg/kg', cnst_longname(m), gridname='physgrid')
  end do
  do m = 1,pcnst
     call add_default(trim(cnst_name(m))//'&IC',0, 'I')
  end do
  call addfld ('PHIS_SM',  horiz_only, 'I', 'm2/s2', 'Surface geopotential (smoothed)',                gridname='physgrid')
  call addfld ('SGH_SM',   horiz_only, 'I', 'm',     'Standard deviation of orography (smoothed)',     gridname='physgrid')
  call addfld ('SGH30_SM', horiz_only, 'I', 'm',     'Standard deviation of 30s orography (smoothed)', gridname='physgrid')


end subroutine stepon_init

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run1 -- Phase 1 of dynamics run method.
!
! !INTERFACE:
subroutine stepon_run1( dtime_out, phys_state, phys_tend,               &
                        pbuf2d, dyn_in, dyn_out )

  use physics_buffer, only : physics_buffer_desc
  use dyn_comp,       only: dyn_run
  use dp_coupling,    only: d_p_coupling
  use time_manager,   only: get_step_size

  implicit none

!
! !OUTPUT PARAMETERS:
!

   real(r8), intent(out) :: dtime_out   ! Time-step
   type(physics_state), intent(inout):: phys_state(begchunk:endchunk)
   !type(physics_tend), intent(out):: phys_tend(begchunk:endchunk)
   !SHP-DEBUG
   type(physics_tend), intent(inout):: phys_tend(begchunk:endchunk)
   type (physics_buffer_desc), pointer :: pbuf2d(:,:)
   type (dyn_import_t), intent(in) :: dyn_in  ! Dynamics import container
   !type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container

#if (! defined SPMD)
   integer  :: mpicom = 0
#endif
!-----------------------------------------------------------------------

!   call t_barrierf('sync_dyn_run', mpicom)
!   call t_startf ('dyn_run')
!   call dyn_run(dyn_in, dyn_out)
!   call t_stopf  ('dyn_run')

   dtime_out = get_step_size()
   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------
   
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling (phys_state, phys_tend, pbuf2d, dyn_out )
   call t_stopf('d_p_coupling')
   
end subroutine stepon_run1

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out )

   use dp_coupling,    only: p_d_coupling

   type(physics_state), intent(inout):: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout):: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif
 
   ! copy from phys structures -> dynamics structures
   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   call p_d_coupling(phys_state, phys_tend, dyn_in)
   call t_stopf('p_d_coupling')


end subroutine stepon_run2


subroutine stepon_run3( dtime, cam_out, phys_state, dyn_in, dyn_out )
   use camsrfexch,  only: cam_out_t     
   use dyn_comp,    only: dyn_run
   use time_manager, only: get_curr_date, get_nstep, is_last_step
   real(r8), intent(in) :: dtime   ! Time-step
   type(cam_out_t), intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout):: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif
   integer :: ncdate            ! current date in integer format [yyyymmdd]
   integer :: ncsec             ! time of day relative to current date [seconds]
   integer :: yr, mon, day      ! year, month, day components of a date

!   call get_curr_date(yr, mon, day, ncsec)
!   ncdate = yr*10000 + mon*100 + day

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
!   write (0,*) 'SHP-MPAS :: dyn nstep, ', get_nstep()
   call dyn_run(dyn_in, dyn_out)	
   call t_stopf  ('dyn_run')


end subroutine stepon_run3


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)
  use dyn_comp, only: dyn_final

! !PARAMETERS:
  type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
!
! !DESCRIPTION:
!
! Deallocate data needed for dynamics. Finalize any dynamics specific
! files or subroutines.
!
!EOP
!-----------------------------------------------------------------------
!BOC

  call dyn_final(dyn_in, dyn_out)

!EOC
end subroutine stepon_final

!-----------------------------------------------------------------------


end module stepon
