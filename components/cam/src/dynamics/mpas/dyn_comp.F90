!
! MODULE: dyn_comp --- Dynamical Core Component

Module dyn_comp

   use shr_kind_mod,   only: r8=>shr_kind_r8
   implicit none
   private
   save

   ! PUBLIC MEMBER FUNCTIONS:
   public dyn_init, dyn_run, dyn_final, dyn_register, dyn_readnl 

   ! PUBLIC DATA MEMBERS:
   public dyn_import_t, dyn_export_t

   type dyn_import_t
       real(r8), dimension(: ),     pointer     :: phis   ! Surface geopotential      (ncol)
       real(r8), dimension(: ),     pointer     :: psd    ! Dry surface pressure      (ncol)
       real(r8), dimension(:,:,: ), pointer     :: uvperp ! Normal velocity at edges  (ncol,nedge,nver)
       real(r8), dimension(:,:   ), pointer     :: ux     ! Lon veloc at center       (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: uy     ! Lat veloc at center       (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: t      ! Temperature               (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: omega  ! Omega                     (ncol,nver+1)
       real(r8), dimension(:,:,: ), pointer     :: tracer ! Tracers                   (ncol,nver,nq)
       real(r8), dimension(:,:   ), pointer     :: ux_tend! Lon veloc tend at center  (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: uy_tend! Lat veloc tend at center  (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: t_tend ! Temperature tendency      (ncol,nver)
   end type dyn_import_t

   type dyn_export_t
       real(r8), dimension(: ),     pointer     :: phis   ! Surface geopotential      (ncol)
       real(r8), dimension(: ),     pointer     :: psd    ! Dry surface pressure      (ncol)
       real(r8), dimension(:,: ),   pointer     :: pint   ! Dry pressure at layer interfaces (ncol,nver+1)
       real(r8), dimension(:,: ),   pointer     :: pmid   ! Dry pressure at layer mid-points (ncol,nver)
       real(r8), dimension(:,: ),   pointer     :: zint   ! Geopotential height 
                                                          !              at layer interfaces (ncol,nver+1)
       real(r8), dimension(:,: ),   pointer     :: zmid   ! Geopotential height 
                                                          !              at layer mid-points (ncol,nver)
       real(r8), dimension(:,:,: ), pointer     :: uvperp ! Normal velocity at edges  (ncol,nedge,nver)
       real(r8), dimension(:,:   ), pointer     :: ux     ! Lon veloc at center       (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: uy     ! Lat veloc at center       (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: t      ! Temperature               (ncol,nver)
       real(r8), dimension(:,:   ), pointer     :: omega  ! Omega                     (ncol,nver+1)
       real(r8), dimension(:,:,: ), pointer     :: tracer ! Tracers                   (ncol,nver,nq)
       real(r8), dimension(:,:   ), pointer     :: pressure! Pressure                 (ncol,nver)
   end type dyn_export_t

! Note: when operating without dynamics, psd will be the total surface pressure
! Note: when operating with dynamics, during initializaiton and restarts t will contain
!       the potential temperature

   character(*), parameter, public :: MODULE_NAME = "dyn_comp"
   character(*), parameter, public :: VERSION     = "$Id$" 

CONTAINS

!========================================================================

subroutine dyn_readnl(NLFileName)
!  use namelist_utils, only: find_group_name
!  use namelist_mod,   only: homme_set_defaults, homme_postprocess_namelist
!  use units,          only: getunit, freeunit
  use cam_logfile,            only: iulog
  use spmd_utils,     only: masterproc, masterprocid, mpicom, npes
!  use spmd_utils,     only: mpi_real8, mpi_integer, mpi_character, mpi_logical
!  use control_mod,    only: TRACERTRANSPORT_SE_GLL, TRACERTRANSPORT_LAGRANGIAN_FVM
!  use control_mod,    only: TRACERTRANSPORT_FLUXFORM_FVM, tracer_transport_type
!  use control_mod,    only: TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM, tracer_grid_type
!  use control_mod,    only: energy_fixer, hypervis_order, hypervis_subcycle
!  use control_mod,    only: hypervis_subcycle_q, integration, statefreq, runtype
!  use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_top, qsplit, rsplit
!  use control_mod,    only: vert_remap_q_alg, tstep_type, rk_stage_user
!  use control_mod,    only: ftype, limiter_option, partmethod
!  use control_mod,    only: topology, sub_case, numnodes, tasknum, moisture
!  use control_mod,    only: columnpackage, remapfreq, remap_type
!  use control_mod,    only: initial_total_mass, use_semi_lagrange_transport
!  use control_mod,    only: disable_diagnostics
!  use control_mod,    only: fine_ne, hypervis_power, hypervis_scaling
!  use control_mod,    only: max_hypervis_courant
!  use fvm_mod,        only: fvm_ideal_test, fvm_test_type
!  use fvm_mod,        only: fvm_get_test_type
!  use dimensions_mod, only: qsize, qsize_d, ntrac, ntrac_d, npsq, ne, npart
!  use constituents,   only: pcnst
!  use params_mod,     only: SFCURVE
!  use parallel_mod,   only: par, initmp
!  use native_mapping, only: native_mapping_readnl
!!XXgoldyXX: v For future CSLAM / physgrid commit
!    use dp_grids,       only: fv_nphys, fv_nphys2, nphys_pts, write_phys_grid, phys_grid_file
!!XXgoldyXX: ^ For future CSLAM / physgrid commit

  ! Dummy argument
  character(len=*), intent(in) :: NLFileName

  ! Local variables
!  integer                      :: unitn, ierr

    if (masterproc) then
      write(iulog, *) "dyn_readnl: reading dyn_mpas_inparm namelist..."
    end if

end subroutine dyn_readnl

  !-------------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  dyn_init --- Initialize the dynamical core
  !
  ! !INTERFACE:
  subroutine dyn_init(dyn_in, dyn_out)

    ! USES:
    use pmgrid,           only: plev, plevp, numcols, maxEdges
    use constituents,     only: pcnst
    use hycoef,           only: hycoef_init
    use pio,              only: file_desc_t
    use cam_history_support, only: add_vert_coord
    use time_manager,     only: get_step_size
    use mpas_cam_interface, only: mpas_init2
#if defined ( SPMD )
    use mpishorthand,   only: mpicom
#else
    integer :: mpicom = 0
#endif

    implicit none

    ! ARGUMENTS:
    type (dyn_import_t), intent(out)  :: dyn_in
    type (dyn_export_t), intent(out)  :: dyn_out

    integer :: k
    real(r8) :: alev(plev)
    real(r8) :: ailev(plevp)

!    real(kind=r8) :: dtime

!    dtime = get_step_size()
!    call mpas_init1(mpicom, real(dtime,8) )

    do k=1,plev
        alev(k) = real(k,r8)
    end do
    do k=1,plevp
        ailev(k) = real(k,r8)
    end do

    call add_vert_coord('lev', plev,                                         &
         'zeta level index at vertical midpoints', '-', alev)

    call add_vert_coord('ilev', plevp,                                       &
         'zeta level index at vertical interfaces', '-', ailev)

! Create dynamics interface arrays

    allocate ( dyn_in%phis   (numcols) )
    allocate ( dyn_in%psd    (numcols) )
    allocate ( dyn_in%uvperp (numcols,maxEdges,plev) )
    allocate ( dyn_in%ux     (numcols,plev) )
    allocate ( dyn_in%uy     (numcols,plev) )
    allocate ( dyn_in%t      (numcols,plev) )
    allocate ( dyn_in%omega  (numcols,plev+1) )
    allocate ( dyn_in%tracer (numcols,plev,pcnst) )
    allocate ( dyn_in%ux_tend(numcols,plev) )
    allocate ( dyn_in%uy_tend(numcols,plev) )
    allocate ( dyn_in%t_tend (numcols,plev) )

    dyn_in%phis(:)       = 0._r8
    dyn_in%psd(:)        = 0._r8
    dyn_in%uvperp(:,:,:) = 0._r8
    dyn_in%ux(:,:)       = 0._r8
    dyn_in%uy(:,:)       = 0._r8
    dyn_in%t(:,:)        = 0._r8
    dyn_in%omega(:,:)    = 0._r8
    dyn_in%tracer(:,:,:) = 0._r8
    dyn_in%ux_tend(:,:)  = 0._r8
    dyn_in%uy_tend(:,:)  = 0._r8
    dyn_in%t_tend(:,:)   = 0._r8

    allocate ( dyn_out%phis   (numcols) )
    allocate ( dyn_out%psd    (numcols) )
    allocate ( dyn_out%pint   (numcols,plev+1) )
    allocate ( dyn_out%pmid   (numcols,plev) )
    allocate ( dyn_out%zint   (numcols,plev+1) )
    allocate ( dyn_out%zmid   (numcols,plev) )
    allocate ( dyn_out%uvperp (numcols,maxEdges,plev) )
    allocate ( dyn_out%ux     (numcols,plev) )
    allocate ( dyn_out%uy     (numcols,plev) )
    allocate ( dyn_out%t      (numcols,plev) )
    allocate ( dyn_out%omega  (numcols,plev+1) )
    allocate ( dyn_out%tracer (numcols,plev,pcnst) )
    allocate ( dyn_out%pressure(numcols,plev) )

    dyn_out%phis(:)       = 0._r8
    dyn_out%psd(:)        = 0._r8
    dyn_out%pint(:,:)     = 0._r8
    dyn_out%pmid(:,:)     = 0._r8
    dyn_out%zint(:,:)     = 0._r8
    dyn_out%zmid(:,:)     = 0._r8
    dyn_out%uvperp(:,:,:) = 0._r8
    dyn_out%ux(:,:)       = 0._r8
    dyn_out%uy(:,:)       = 0._r8
    dyn_out%t(:,:)        = 0._r8
    dyn_out%omega(:,:)    = 0._r8
    dyn_out%tracer(:,:,:) = 0._r8
    dyn_out%pressure(:,:) = 0._r8

    call mpas_init2

  end subroutine dyn_init
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  Driver for the MPAS dynamical core 
  !
  ! !INTERFACE:
  subroutine dyn_run (dyn_in, dyn_out)

    use mpas_cam_interface,  only: mpas_dyn_run

    implicit none

    type (dyn_import_t), intent(inout)  :: dyn_in
    type (dyn_export_t), intent(inout)  :: dyn_out

!    if (skip_dynamics) then
! 'Advance' dynamics by copying dyn_in into dyn_out
!       dyn_out%phis(:)       = dyn_in%phis(:)
!       dyn_out%psd(:)        = dyn_in%psd(:)
!       dyn_out%uvperp(:,:,:) = dyn_in%uvperp(:,:,:)
!       dyn_out%ux(:,:)       = dyn_in%ux(:,:)
!       dyn_out%uy(:,:)       = dyn_in%uy(:,:)
!       dyn_out%t(:,:)        = dyn_in%t(:,:)
!       dyn_out%omega(:,:)    = dyn_in%omega(:,:)
!       dyn_out%tracer(:,:,:) = dyn_in%tracer(:,:,:)
!    else
       call mpas_dyn_run()
!    endif

  end subroutine dyn_run
  !-----------------------------------------------------------------------
  !BOP
  ! !ROUTINE:  Finalization for the MPAS dynamical core 
  !
  ! !INTERFACE:
  subroutine dyn_final (dyn_in, dyn_out)

    use mpas_cam_interface,  only: mpas_final

    implicit none

    type (dyn_import_t), intent(inout)  :: dyn_in
    type (dyn_export_t), intent(inout)  :: dyn_out

    deallocate ( dyn_in%phis )
    deallocate ( dyn_in%psd )
    deallocate ( dyn_in%uvperp )
    deallocate ( dyn_in%ux )
    deallocate ( dyn_in%uy )
    deallocate ( dyn_in%t )
    deallocate ( dyn_in%omega )
    deallocate ( dyn_in%tracer )
    deallocate ( dyn_in%ux_tend )
    deallocate ( dyn_in%uy_tend )
    deallocate ( dyn_in%t_tend )

    deallocate ( dyn_out%phis )
    deallocate ( dyn_out%psd )
    deallocate ( dyn_out%pint )
    deallocate ( dyn_out%pmid )
    deallocate ( dyn_out%zint )
    deallocate ( dyn_out%zmid )
    deallocate ( dyn_out%uvperp )
    deallocate ( dyn_out%ux )
    deallocate ( dyn_out%uy )
    deallocate ( dyn_out%t )
    deallocate ( dyn_out%omega )
    deallocate ( dyn_out%tracer )

    call mpas_final()

  end subroutine dyn_final

!========================================================================
subroutine read_inidat(dyn_in )

   use pio, only : file_desc_t
!   use dyn_comp, only : dyn_import_t
   use pmgrid, only : plev, numcols, maxEdges
   use constituents, only : cnst_name, cnst_read_iv, pcnst
   use shr_kind_mod, only : r8 => shr_kind_r8
   use spmd_utils, only : iam, npes
   use dyn_grid, only : nCells, ncells_per_proc
   use cam_initfiles, only: initial_file_get_id, topo_file_get_id
   use physconst, only : cappa
   
   implicit none

   type(dyn_import_t), target, intent(out) :: dyn_in

   type(file_desc_t), pointer :: fh_ini
   type(file_desc_t), pointer :: fh_topo

   integer :: p, nCellsLocal, io_master_id, m, i, k
   integer, allocatable, dimension(:) :: proc_start_i
   real(r8), allocatable, dimension(:,:) :: tmp_in, tmp_out

! primary output quantities are psd, phis, t, uvperp
! tracers are currently initialized to zero if not on IC file
! for cases without dynamics, psd will contain total surface pressure
!   and t will contain temperature

   real(r8), pointer :: psd(:)        !  psd(numcols)               ! dry surface pressure
   real(r8), pointer :: phis(:)       !  phis(numcols)              ! surface geopotential
   real(r8), pointer :: pt(:,:)       !  pt(numcols,plev)           ! potential temperature
   real(r8), pointer :: ux(:,:)       !  ux(numcols,plev)           ! cell centered velocity
   real(r8), pointer :: uy(:,:)       !  uy(numcols,plev)           ! cell centered velocity
   real(r8), pointer :: uvperp(:,:,:) !  uvperp(numcols,maxEdges,plev) ! edge normal velocity
   real(r8), pointer :: tracer(:,:,:) !  tracer(numcols,plev,pcnst) ! tracers

   real(r8), allocatable :: pdry(:,:)
   real(r8), allocatable :: delpp(:,:)
   real(r8), allocatable :: ptot(:,:)
   real(r8), allocatable :: pcap(:,:)

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()

   psd    => dyn_in%psd
   phis   => dyn_in%phis
   pt     => dyn_in%t
   uvperp => dyn_in%uvperp
   ux     => dyn_in%ux
   uy     => dyn_in%uy
   tracer => dyn_in%tracer

   allocate( proc_start_i( npes ) )

   proc_start_i(1) = 1
   do p = 2, npes
      proc_start_i(p) = proc_start_i(p-1) + ncells_per_proc(p-1)
   enddo

! For CAM/MPAS simulations, PS contains the dry surface pressure
! For CAM with fixed dynamics, PS contains the total pressure

   !MGD IO -- call MPAS_read_initial(ncid_ini, ncid_topo?) 
   !MGD IO -- fill in dyn_in at this point?

end subroutine read_inidat
!========================================================================

subroutine dyn_register()

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use ppgrid,          only: pcols, pver
   use phys_control,    only: use_gw_front, use_gw_front_igw

   ! These fields are computed by the dycore and passed to the physics via the
   ! physics buffer.

!   if (use_gw_front .or. use_gw_front_igw) then
!      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
!         frontgf_idx)
!      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
!         frontga_idx)
!   end if

end subroutine dyn_register


end module dyn_comp
