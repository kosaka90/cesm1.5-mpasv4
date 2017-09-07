module restart_dynamics

use shr_kind_mod,   only: r8 => shr_kind_r8, r4=> shr_kind_r4
use pmgrid,         only: plev, nEdges, nVertices, numcols, maxEdges, plevp
use dyn_grid,       only: gcol, ncells

use pio, only: var_desc_t
use pio, only: pio_double, pio_def_dim, pio_def_var
use perf_mod
use time_manager,   only: get_step_size

implicit none

public read_restart_dynamics, write_restart_dynamics, init_restart_dynamics

type(var_desc_t) :: PHISdesc, PSdesc, timedesc, Tdesc, Tracerdesc, UVPERPdesc
type(var_desc_t) :: UXdesc, UYdesc, OMEGAdesc
integer :: ncol_dimid, nlev_dimid, nlevp_dimid
real(r8), parameter ::  SECONDS_IN_DAY          =  86400_r8
logical, save :: initialized=.false.

!=========================================================================================
CONTAINS
!=========================================================================================

  subroutine init_restart_dynamics( File, dyn_out )

    use dyn_comp, only : dyn_export_t
    use pio, only : file_desc_t, pio_def_dim, pio_put_att, pio_global, &
         pio_setdebuglevel, pio_unlimited
    use dyn_grid, only : get_horiz_grid_dim_d
    ! ... nCells,nEdges,nVertices are global numbers, numcols local (to the processor)
    ! where column and cell are synonyms.
    ! maxEdges is the maximum number of edges per cell.
    use constituents,     only: pcnst
    use mpas_cam_interface, only : MPAS_init_restart

    type(file_desc_t) :: File
    type(dyn_export_t), intent(in)  :: dyn_out

    integer :: ncols
    integer :: pcnst_dimid, time_dimid
    integer :: vdimids(2)
    integer :: ierr

    call PIO_Setdebuglevel(0) ! 1&3 set debug, 2&3 set debugio, 0 for none

    ierr = PIO_Put_Att( File, PIO_GLOBAL, 'pcnst', pcnst )
    ierr = PIO_Put_Att( File, PIO_GLOBAL, 'plev', plev )

    call get_horiz_grid_dim_d(ncols)
    ierr = PIO_Def_Dim(File, 'ncol', ncols, ncol_dimid)

    ierr = PIO_Def_Dim(File, 'pcnst', pcnst, pcnst_dimid )

    nlev_dimid  = vdimids(1) ! 
    nlevp_dimid = vdimids(2) ! nlev_dimid+1
    !ierr = PIO_Def_Dim(File, 'nlev', nlev, nlev_dimid)
    !ierr = PIO_Def_Dim(File, 'nlevp', nlevp, nlevp_dimid)

    !ierr = PIO_Def_Dim(File, 'time', PIO_UNLIMITED,time_dimid)
    !ierr = PIO_Def_Var(File,'time', pio_double, (/time_dimid/), timedesc)

    ierr = PIO_Def_Var(File,'time', pio_double, timedesc)

    call MPAS_init_restart(File)

  end subroutine init_restart_dynamics

  !========================================================================

  subroutine write_restart_dynamics( File, dyn_out )

    use dyn_comp,         only: dyn_import_t, dyn_export_t
    use pio, only : file_desc_t, io_desc_t, pio_double, pio_put_var, &
         pio_setdebuglevel, pio_write_darray
    use time_manager,       only: get_curr_time
    use spmd_utils,         only: iam
    use constituents,       only: pcnst
    use shr_sys_mod,        only: shr_sys_flush
#ifdef SPMD
    use mpishorthand,   only: mpicom
#endif
    use mpas_cam_interface, only : MPAS_write_restart
   

    ! ARGUMENTS
    type(File_desc_t), intent(inout) :: File     ! Unit number
    type(dyn_export_t), intent(in)  :: dyn_out

    ! LOCALS
    type(io_desc_t)       :: iodesc2d, iodesc3d, iodesctr, iodescu, iodescw  ! I/O descriptors
    integer               :: ndcur, nscur        ! Time (days/seconds)
    real(r8)              :: time                ! current time 
    integer               :: ierr                ! error flag
    real(kind=r8),pointer :: var3d(:), var2d(:), vartr(:), varu(:)
    integer :: i, j, k, icnt
    type(dyn_import_t)    :: dyn_in ! for testing ounly
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif

    call MPAS_write_restart(File)

  end subroutine write_restart_dynamics

  !========================================================================

  subroutine read_restart_dynamics( File, dyn_in, dyn_out )

    use dyn_comp, only : dyn_init, dyn_import_t, dyn_export_t
    use pio,      only : file_desc_t, pio_inq_varid, pio_global, pio_get_att, &
         pio_read_darray, io_desc_t
    use cam_logfile, only : iulog
    use spmd_utils,  only : iam
    use dyn_grid,    only : dyn_grid_init
    use constituents, only :  pcnst
    use cam_control_mod,only: initial_run, restart_run, branch_run
#if defined ( SPMD )
    use mpishorthand, only : mpicom
#else
    integer :: mpicom = 0
#endif
    use mpas_cam_interface, only : mpas_init1, mpas_init2
    !use mpas_configure, only : config_input_name
   

    implicit none

    ! ARGUMENTS
    type(File_desc_t),  intent(inout)    :: File
    type(dyn_import_t), intent(inout)    :: dyn_in
    type(dyn_export_t), intent(inout)    :: dyn_out

    ! LOCALS
    type(io_desc_t)       :: iodesc2d, iodesc3d, iodesctr, iodescu, iodescw  ! I/O descriptors
    real(kind=r8),pointer :: var3d(:), var2d(:), vartr(:), varu(:)
    integer               :: i, j, icnt, ierr, k
    integer               :: fplev, fnCells, fnEdges, fnVertices, fnumcols, fpcnst
    real(r8) :: dtime
    logical  :: initialization

    if (restart_run .or. branch_run) initialization=.true.

    dtime = get_step_size() 
    !call mpas_init1(mpicom, real(dtime,8), initial_run, File )
    call mpas_init1(mpicom, real(dtime,8), File )

    ! Initialize the dynamics
    !
    ! MGD TODO: The filename here should actually be the restart filename, but this
    ! can be difficult to construct without accessing quite a few variables from
    ! MPAS modules; for now, just use the startup filename, since grid
    ! information won't change over time...
    !
    !call dyn_grid_init( trim(config_input_name) )
    !call dyn_grid_init()

    !call dyn_init( file, nlfilename, dyn_in, dyn_out )
    call dyn_init( dyn_in, dyn_out )

    ierr = PIO_Get_Att( File, PIO_GLOBAL, 'pcnst', fpcnst )
    ierr = PIO_Get_Att( File, PIO_GLOBAL, 'plev', fplev )
    ierr = PIO_Get_Att( File, PIO_GLOBAL, 'nCells', fnCells )
    ierr = PIO_Get_Att( File, PIO_GLOBAL, 'nEdges', fnEdges )
    ierr = PIO_Get_Att( File, PIO_GLOBAL, 'nVertices', fnVertices )
    if ( plev/=fplev .or. pcnst/=fpcnst .or. nCells/=fnCells .or. &
         nEdges/=fnEdges .or. nVertices/=fnVertices ) then
       write(iulog,*) 'Restart file grid attributes do not match current model.', &
            ' (File,model) values of:  plev ', fplev,plev, ' pcnst ', fpcnst,pcnst, &
            '  nCells ', fnCells,nCells, '  nEdges ', fnEdges,nEdges, &
            '  nVertices ', fnVertices,nVertices
    end if

!    if (skip_dynamics) then
! In absence of actual dycore, copy initial conditions into dyn_out
!   structure in preparation for d_p_coupling in stepon_run1
! In this case,
!    dyn_in%t contains the temperature
!    dyn_in%psd contains the total surface pressure
!
!       dyn_out%phis(:)       = dyn_in%phis(:)
!       dyn_out%psd(:)        = dyn_in%psd(:)
!       dyn_out%uvperp(:,:,:) = dyn_in%uvperp(:,:,:)
!       dyn_out%ux(:,:)       = dyn_in%ux(:,:)
!       dyn_out%uy(:,:)       = dyn_in%uy(:,:)
!       dyn_out%t(:,:)        = dyn_in%t(:,:)
!       dyn_out%omega(:,:)    = dyn_in%omega(:,:)
!       dyn_out%tracer(:,:,:) = dyn_in%tracer(:,:,:)
!    else
! dyn_in%t contains the potential temperature
       !MGD IO -- call MPAS_read_restart(File)
       !call mpas_init2
!    endif

    return

  end subroutine read_restart_dynamics

  !========================================================================

end module restart_dynamics
