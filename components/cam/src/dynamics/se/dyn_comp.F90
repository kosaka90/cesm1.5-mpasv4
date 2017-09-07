module dyn_comp

  ! This module implements the CAM interfaces to the SE Dynamical Core

  use shr_kind_mod,           only: r8=>shr_kind_r8, shr_kind_cl
  use spmd_utils,             only: iam, masterproc
  use dyn_grid,               only: timelevel, hvcoord
  use cam_control_mod,        only: initial_run

  use cam_logfile,            only: iulog
  use cam_abortutils,         only: endrun

  use fvm_control_volume_mod, only: fvm_struct
  use element_mod,            only: element_t, elem_state_t
  use time_mod,               only: nsplit
  use hybrid_mod,             only: hybrid_t

  implicit none

  private
  save

  public ::          &
       dyn_import_t, &
       dyn_export_t, &
       dyn_readnl,   &
       dyn_register, &
       dyn_init,     &
       dyn_run,      &
       dyn_final

  type dyn_import_t
    type (element_t),  pointer :: elem(:) => null()
    type (fvm_struct), pointer :: fvm(:) => null()
  end type dyn_import_t

  type dyn_export_t
    type (element_t),  pointer :: elem(:) => null()
    type (fvm_struct), pointer :: fvm(:) => null()
  end type dyn_export_t

  ! Parameters for namelist string lengths (from namelist_definition.xml)
  integer, parameter            :: METHOD_LEN             = 32

  integer, parameter  ::  DYN_RUN_SUCCESS                   = 0
  integer, parameter  ::  DYN_RUN_FAILURE                   = -1

  ! Frontogenesis indices
  integer, public :: frontgf_idx = -1
  integer, public :: frontga_idx = -1

  ! Private interfaces
  interface read_dyn_var
    module procedure read_dyn_field_2d
    module procedure read_dyn_field_3d
  end interface read_dyn_var

!===============================================================================
contains
!===============================================================================

  subroutine dyn_readnl(NLFileName)
    use namelist_utils, only: find_group_name
    use namelist_mod,   only: homme_set_defaults, homme_postprocess_namelist
    use units,          only: getunit, freeunit
    use spmd_utils,     only: masterproc, masterprocid, mpicom, npes
    use spmd_utils,     only: mpi_real8, mpi_integer, mpi_character, mpi_logical
    use control_mod,    only: TRACERTRANSPORT_SE_GLL
    use control_mod,    only: TRACERTRANSPORT_CONSISTENT_SE_FVM
    use control_mod,    only: tracer_transport_type
    use control_mod,    only: TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM, tracer_grid_type
    use control_mod,    only: energy_fixer, hypervis_order, hypervis_subcycle
    use control_mod,    only: hypervis_subcycle_q, integration, statefreq, runtype
    use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_top, qsplit, rsplit
    use control_mod,    only: vert_remap_q_alg, tstep_type, rk_stage_user
    use control_mod,    only: ftype, limiter_option, partmethod
    use control_mod,    only: topology, sub_case, numnodes, tasknum
    use control_mod,    only: columnpackage, remapfreq, remap_type
    use control_mod,    only: initial_total_mass
    use control_mod,    only: disable_diagnostics
    use control_mod,    only: fine_ne, hypervis_power, hypervis_scaling
    use control_mod,    only: max_hypervis_courant
    use control_mod,    only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
    use dimensions_mod, only: qsize, qsize_d, ntrac, npsq, ne, npart, ldry_mass_vertical_coordinates
    use dimensions_mod, only: fv_nphys, qsize_condensate_loading, lcp_moist
    use constituents,   only: pcnst
    use params_mod,     only: SFCURVE
    use parallel_mod,   only: par, initmpi
    use thread_mod,     only: initomp, max_num_threads
    use thread_mod,     only: horz_num_threads, vert_num_threads, tracer_num_threads
    use dp_grids,       only: nphys_pts, se_write_phys_grid, se_phys_grid_file
    use native_mapping, only: native_mapping_readnl

    ! Dummy argument
    character(len=*), intent(in) :: NLFileName

    ! Local variables
    integer                      :: unitn, ierr
   ! SE Namelist variables
    integer                      :: se_qsize_condensate_loading
    integer                      :: se_fine_ne
    integer                      :: se_ftype
    integer                      :: se_fv_nphys
    integer                      :: se_hypervis_order
    real(r8)                     :: se_hypervis_power
    real(r8)                     :: se_hypervis_scaling
    integer                      :: se_hypervis_subcycle
    integer                      :: se_hypervis_subcycle_q
    integer                      :: se_limiter_option
    real(r8)                     :: se_max_hypervis_courant
    character(len=SHR_KIND_CL)   :: se_mesh_file
    integer                      :: se_ne
    integer                      :: se_npes
    integer                      :: se_nsplit
    real(r8)                     :: se_nu
    real(r8)                     :: se_nu_div
    real(r8)                     :: se_nu_p
    real(r8)                     :: se_nu_q
    real(r8)                     :: se_nu_top
    integer                      :: se_qsplit
    logical                      :: se_refined_mesh
    integer                      :: se_rsplit
    integer                      :: se_statefreq
    integer                      :: se_tstep_type
    integer                      :: se_vert_remap_q_alg
    character(len=METHOD_LEN)    :: se_tracer_transport_method
    integer                      :: se_horz_num_threads
    integer                      :: se_vert_num_threads
    integer                      :: se_tracer_num_threads

    namelist /dyn_se_inparm/          &
         se_qsize_condensate_loading, &
         se_fine_ne,                  & ! For refined meshes
         se_ftype,                    & ! forcing type
         se_fv_nphys,                 &
         se_hypervis_order,           &
         se_hypervis_power,           &
         se_hypervis_scaling,         &
         se_hypervis_subcycle,        &
         se_hypervis_subcycle_q,      &
         se_limiter_option,           &
         se_max_hypervis_courant,     &
         se_mesh_file,                & ! Refined mesh definition file
         se_ne,                       &
         se_npes,                     &
         se_nsplit,                   & ! # of dyn steps per physics timestep
         se_nu,                       &
         se_nu_div,                   &
         se_nu_p,                     &
         se_nu_q,                     &
         se_nu_top,                   &
         se_qsplit,                   &
         se_refined_mesh,             &
         se_rsplit,                   &
         se_statefreq,                & ! number of steps per printstate call
         se_tstep_type,               &
         se_vert_remap_q_alg,         &
         se_met_nudge_u,              &
         se_met_nudge_p,              &
         se_met_nudge_t,              &
         se_met_tevolve,              &
         se_tracer_transport_method,  & ! se_gll or cslam_fvm
         se_write_phys_grid,          &
         se_phys_grid_file,           &
         se_horz_num_threads,         &
         se_vert_num_threads,         &
         se_tracer_num_threads

    !--------------------------------------------------------------------------

   ! namelist default values should be in namelist (but you know users . . .)
    ! NB: Of course, these should keep up with what is in namelist_defaults ...
    se_qsize_condensate_loading = 1
    se_fine_ne                  = -1
    se_ftype                    = 0
    se_fv_nphys                 = 0
    se_hypervis_order           = 2
    se_hypervis_power           = 0
    se_hypervis_scaling         = 0
    se_hypervis_subcycle        = 3
    se_hypervis_subcycle_q      = 1
    se_limiter_option           = 8
    se_max_hypervis_courant     = 1.0e99_r8
    se_mesh_file                = ''
    se_ne                       = -1
    se_npes                     = npes
    se_nsplit                   = 2
    se_nu                       = 1.0e15_r8
    se_nu_div                   = 2.5e15_r8
    se_nu_p                     = 1.0e15_r8
    se_nu_q                     = -1.0_r8
    se_nu_top                   = 2.5e5_r8
    se_qsplit                   = 1
    se_refined_mesh             = .false.
    se_rsplit                   = 3
    se_statefreq                = 480
    se_tracer_transport_method  = 'se_gll'
    se_tstep_type               = 5
    se_vert_remap_q_alg         = 1
    ! Physgrid
    se_write_phys_grid          = 'no'
    se_phys_grid_file           = 'phys_grid.nc' ! Initialized in dp_grids
    se_horz_num_threads         = 0
    se_vert_num_threads         = 0
    se_tracer_num_threads       = 0

    ! Read the namelist (dyn_se_inparm)
    call MPI_barrier(mpicom, ierr)
    if (masterproc) then
      write(iulog, *) "dyn_readnl: reading dyn_se_inparm namelist..."
      unitn = getunit()
      open( unitn, file=trim(NLFileName), status='old' )
      call find_group_name(unitn, 'dyn_se_inparm', status=ierr)
      if (ierr == 0) then
        read(unitn, dyn_se_inparm, iostat=ierr)
        if (ierr /= 0) then
          call endrun('dyn_readnl: ERROR reading dyn_se_inparm namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if

    ! Broadcast namelist values to all PEs
    call MPI_bcast(se_qsize_condensate_loading, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_fine_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_ftype, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_order, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_power, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_scaling, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_subcycle, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_hypervis_subcycle_q, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_limiter_option, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_max_hypervis_courant, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_mesh_file, SHR_KIND_CL,  mpi_character, masterprocid, mpicom, ierr)
    call MPI_bcast(se_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_npes, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_div, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_p, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_q, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_nu_top, 1, mpi_real8, masterprocid, mpicom, ierr)
    call MPI_bcast(se_qsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_refined_mesh, 1, mpi_logical, masterprocid, mpicom, ierr)
    call MPI_bcast(se_rsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_statefreq, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_tracer_transport_method, METHOD_LEN,  mpi_character, masterprocid, mpicom, ierr)
    call MPI_bcast(se_tstep_type, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_vert_remap_q_alg, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_met_nudge_u, 1, MPI_real8, masterprocid, mpicom,ierr)
    call MPI_bcast(se_met_nudge_p, 1, MPI_real8, masterprocid, mpicom,ierr)
    call MPI_bcast(se_met_nudge_t, 1, MPI_real8, masterprocid, mpicom,ierr)
    call MPI_bcast(se_met_tevolve, 1, MPI_integer, masterprocid, mpicom,ierr)
    call MPI_bcast(se_fv_nphys, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(se_write_phys_grid, 80,  mpi_character, masterprocid, mpicom, ierr)
    call MPI_bcast(se_phys_grid_file, shr_kind_cl, mpi_character, masterprocid, mpicom, ierr)
    call MPI_bcast(se_horz_num_threads, 1, MPI_integer, masterprocid, mpicom,ierr)
    call MPI_bcast(se_vert_num_threads, 1, MPI_integer, masterprocid, mpicom,ierr)
    call MPI_bcast(se_tracer_num_threads, 1, MPI_integer, masterprocid, mpicom,ierr)

    ! Initialize the SE structure that holds the MPI decomposition information
    if (se_npes <= 0) then
      se_npes = npes
    end if
    par = initmpi(se_npes)
    call initomp()

    ! Fix up unresolved default values
    ! default diffusion coefficiets
    if (se_nu_q < 0) then
      se_nu_q = se_nu
    end if
    if (se_nu_div < 0) then
      se_nu_div = se_nu
    end if
    ! Go ahead and enforce ne = 0 for refined mesh runs
    if (se_refined_mesh) then
      se_ne = 0
    end if

    ! Set HOMME defaults
    call homme_set_defaults()
    ! Set HOMME variables not in CAM's namelist but with different CAM defaults
    partmethod               = SFCURVE
    npart                    = se_npes
    energy_fixer             = -1    ! no fixer, non-staggered-in-time formulas
    ! CAM requires forward-in-time, subcycled dynamics
    ! RK2 3 stage tracers, sign-preserving conservative
    rk_stage_user            = 3
    integration              = 'explicit'
    qsize                    = qsize_d
    topology                 = "cube"
    ! Finally, set the HOMME variables which have different names
    qsize_condensate_loading = se_qsize_condensate_loading
    fine_ne                  = se_fine_ne
    ftype                    = se_ftype
    hypervis_order           = se_hypervis_order
    hypervis_power           = se_hypervis_power
    hypervis_scaling         = se_hypervis_scaling
    hypervis_subcycle        = se_hypervis_subcycle
    hypervis_subcycle_q      = se_hypervis_subcycle_q
    limiter_option           = se_limiter_option
    max_hypervis_courant     = se_max_hypervis_courant
    ne                       = se_ne
    nsplit                   = se_nsplit
    nu                       = se_nu
    nu_div                   = se_nu_div
    nu_p                     = se_nu_p
    nu_q                     = se_nu_q
    nu_top                   = se_nu_top
    qsplit                   = se_qsplit
    rsplit                   = se_rsplit
    statefreq                = se_statefreq
    tstep_type               = se_tstep_type
    vert_remap_q_alg         = se_vert_remap_q_alg
    fv_nphys                 = se_fv_nphys
    if (fv_nphys>0) then
      nphys_pts = fv_nphys*fv_nphys
    else
      nphys_pts = npsq
      ! Cannot write phys_grid file if the physics grid is not active
      se_write_phys_grid = 'no'
      se_phys_grid_file = ''
    end if

    ! Set tracer transport type and tracer numbers
    if (trim(se_tracer_transport_method) == 'se_gll') then
      tracer_transport_type = TRACERTRANSPORT_SE_GLL
      tracer_grid_type = TRACER_GRIDTYPE_GLL
      qsize = pcnst
      if (fv_nphys>0) then!shouldn't this be ntrac>0xxxxxxxx
        ntrac = pcnst !needed for mapping to/from physgrid
      else
        ntrac = 0
      end if
    else if (trim(se_tracer_transport_method) == 'cslam_fvm') then
      tracer_transport_type = TRACERTRANSPORT_CONSISTENT_SE_FVM
      tracer_grid_type = TRACER_GRIDTYPE_FVM
      qsize = 1
      qsize = pcnst !add phl - xxx
      ntrac = pcnst
      if (fv_nphys<0) &
           call endrun('dyn_comp: if running cslam_fvm then fv_nphys>0')
    else
      call endrun('dyn_comp: Unknown tracer transport method: '//trim(se_tracer_transport_method))
    end if

    if ((ntrac>0.or.fv_nphys>0).and..not.ldry_mass_vertical_coordinates) then
      call endrun('fvm transport and/or physics grid requires ldry_mass_vertical_coordinates-.true.')
    end if

    ! if restart or branch run
    if (.not. initial_run) runtype = 1

    ! HOMME wants 'none' to indicate no mesh file
    if (len_trim(se_mesh_file) == 0) then
      se_mesh_file = 'none'
      if (se_refined_mesh) then
        call endrun('dyn_readnl ERROR: se_refined_mesh=.true. but no se_mesh_file')
      end if
    end if
    call homme_postprocess_namelist(se_mesh_file, par)

    ! Set threading numbers to reasonable values
    if ((se_horz_num_threads == 0) .and. (se_vert_num_threads == 0) .and. (se_tracer_num_threads == 0)) then
      ! The user has not set any threading values, choose defaults
      se_horz_num_threads = max_num_threads
      se_vert_num_threads = 1
      se_tracer_num_threads = se_vert_num_threads
    end if
    if (se_horz_num_threads < 1) then
      se_horz_num_threads = 1
    end if
    if (se_vert_num_threads < 1) then
      se_vert_num_threads = 1
    end if
    if (se_tracer_num_threads < 1) then
      se_tracer_num_threads = 1
    end if
    horz_num_threads = se_horz_num_threads
    vert_num_threads = se_vert_num_threads
    tracer_num_threads = se_tracer_num_threads

    if (qsize_condensate_loading==1) then
      lcp_moist=.false.
    else if (qsize_condensate_loading==3) then
      lcp_moist=.true.
    end if

    if (masterproc) then
      write(iulog, *       ) '            lcp_moist= ',lcp_moist
      write(iulog, '(a,i0)') 'dyn_readnl: se_ftype = ',ftype
      write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_order = ',se_hypervis_order
      write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_subcycle = ',se_hypervis_subcycle
      write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_subcycle_q = ',se_hypervis_subcycle_q
      write(iulog, '(a,i0)') 'dyn_readnl: se_limiter_option = ',se_limiter_option
      if (.not. se_refined_mesh) then
        write(iulog, '(a,i0)') 'dyn_readnl: se_ne = ',se_ne
      end if
      write(iulog, '(a,i0)') 'dyn_readnl: se_npes = ',se_npes
      write(iulog, '(a,i0)') 'dyn_readnl: se_nsplit = ',se_nsplit
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu = ',se_nu
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_div = ',se_nu_div
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_p = ',se_nu_p
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_q = ',se_nu_q
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_top = ',se_nu_top
      write(iulog, '(a,i0)') 'dyn_readnl: se_qsplit = ',se_qsplit
      write(iulog, '(a,i0)') 'dyn_readnl: se_rsplit = ',se_rsplit
      write(iulog, '(a,i0)') 'dyn_readnl: se_statefreq = ',se_statefreq
      write(iulog, '(a,i0)') 'dyn_readnl: se_tstep_type = ',se_tstep_type
      write(iulog, '(a,i0)') 'dyn_readnl: se_vert_remap_q_alg = ',se_vert_remap_q_alg
      write(iulog, '(a,i0)') 'dyn_readnl: se_qsize_condensate_loading = ',se_qsize_condensate_loading
      if (se_refined_mesh) then
        write(iulog, '(a)') 'dyn_readnl: Refined mesh simulation'
        write(iulog, '(a)') 'dyn_readnl: se_mesh_file = ',trim(se_mesh_file)
        if (abs(se_hypervis_power) < 1.0e-12_r8) then
          write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power, ', (tensor hyperviscosity)'
          write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_scaling = ',se_hypervis_scaling
        else if (abs(se_hypervis_power - 3.322_r8) < 1.0e-12_r8) then
          write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power, ', (scalar hyperviscosity)'
          write(iulog, '(a,i0)') 'dyn_readnl: se_fine_ne = ',se_fine_ne
        else
          write(iulog, '(a,i0)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power
          write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_scaling = ',se_hypervis_scaling
          write(iulog, '(a,e11.4)') 'dyn_readnl: se_fine_ne = ',se_fine_ne
        end if
        write(iulog, '(a,e11.4)') 'dyn_readnl: se_max_hypervis_courant = ',se_max_hypervis_courant
      end if
      if ((se_met_nudge_u /= 0._r8) .or. (se_met_nudge_p /= 0._r8) .or.       &
           (se_met_nudge_t /= 0._r8) .or. (se_met_tevolve /= 0)) then
        write(iulog, '(a)') 'dyn_readnl: Nudging:'
        write(iulog,'(a,e14.6)') "          :  se_met_nudge_u = ", se_met_nudge_u
        write(iulog,'(a,e14.6)') "          :  se_met_nudge_p = ", se_met_nudge_p
        write(iulog,'(a,e14.6)') "          :  se_met_nudge_t = ", se_met_nudge_t
        write(iulog,'(a,i0)')    "          :  se_met_tevolve = ", se_met_tevolve
      else
        write(iulog, '(a)') 'dyn_readnl: Nudging off'
      end if

      write(iulog,'(a,i0)') 'dyn_readnl: nphys_pts = ', nphys_pts
      if (fv_nphys>0) then
        write(iulog,'(a,i0)') 'dyn_readnl: nphys_pts = ', fv_nphys
        if (trim(se_write_phys_grid) == 'SCRIP') then
          write(iulog,'(2a)') "dyn_readnl: write physics grid file = ", trim(se_phys_grid_file)
        else
          write(iulog,'(a)') "dyn_readnl: do not write physics grid file"
        end if
      else
        write(iulog, '(a)') 'dyn_readnl: physics will run on SE GLL points'
      end if
      write(iulog, '(2a)') 'dyn_readnl: tracer transport method: ', trim(se_tracer_transport_method)
      write(iulog, '(a,i0)') 'dyn_readnl: se_horz_num_threads = ',horz_num_threads
      write(iulog, '(a,i0)') 'dyn_readnl: se_vert_num_threads = ',vert_num_threads
      write(iulog, '(a,i0)') 'dyn_readnl: se_tracer_num_threads = ',tracer_num_threads
    end if

    call native_mapping_readnl(NLFileName)

  end subroutine dyn_readnl

  !=============================================================================
  subroutine dyn_register()

    use physics_buffer,  only: pbuf_add_field, dtype_r8
    use ppgrid,          only: pcols, pver
    use phys_control,    only: use_gw_front, use_gw_front_igw

    ! These fields are computed by the dycore and passed to the physics via the
    ! physics buffer.

    if (use_gw_front .or. use_gw_front_igw) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/),       &
           frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/),       &
           frontga_idx)
    end if

  end subroutine dyn_register

  !=============================================================================

  subroutine dyn_init(dyn_in, dyn_out)

    use dyn_grid,         only: elem, fvm
    use cam_control_mod,  only: aqua_planet, ideal_phys, adiabatic
    use constituents,     only: cnst_get_ind
    use cam_instance,     only: inst_index
    use native_mapping,   only: create_native_mapping_files
    use cam_pio_utils,    only: clean_iodesc_list

    use dimensions_mod,   only: nlev, nelemd, qsize_condensate_loading
    use prim_driver_mod,  only: prim_init2, prim_set_mass
    use parallel_mod,     only: par
    use time_mod,         only: time_at
    use control_mod,      only: runtype
    use thread_mod,       only: horz_num_threads
    use nctopo_util_mod,  only: nctopo_util_driver
    use hybrid_mod,       only: get_loop_ranges, config_thread_region
    use hycoef,           only: ps0
    use inic_analytic,    only: analytic_ic_active

    ! Dummy arguments:
    type(dyn_import_t), intent(out) :: dyn_in
    type(dyn_export_t), intent(out) :: dyn_out

    ! Local variables
    integer             :: ithr, nets, nete, ie, k
    real(r8), parameter :: Tinit = 300.0_r8
    type(hybrid_t)      :: hybrid
    real(r8) :: dyn_ps0
    !--------------------------------------------------------------------------

    ! Check condensate loading constituents
    if (qsize_condensate_loading > 1) then
      if (qsize_condensate_loading /= 3) then
        call endrun('dyn_init: se_qsize_condensate_loading must be 1 or 3')
      end if
      call cnst_get_ind('CLDLIQ', k,  abort=.false.)
      call cnst_get_ind('CLDICE', ie, abort=.false.)
      if (k /= 2) then
        call endrun("dyn_init: tracer, CLDLIQ, must have index 2")
      end if
      if (ie /= 3) then
        call endrun("dyn_init: tracer, CLDICE, must have index 3")
      end if
    end if

    ! Initialize the import/export objects
    dyn_in%elem  => elem
    dyn_in%fvm   => fvm

    dyn_out%elem => elem
    dyn_out%fvm  => fvm

    dyn_ps0=ps0
    !
    ! This subroutine creates mapping files using SE basis functions if requested
    !
    call create_native_mapping_files(par, elem, 'native')
    call create_native_mapping_files(par, elem, 'bilin')

    if (initial_run) then
      call read_inidat(dyn_in)
      call clean_iodesc_list()
    end if

    if(iam < par%nprocs) then

!$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,ie)
      hybrid = config_thread_region(par,'horizontal')
!      hybrid = config_thread_region(par,'serial')
      call get_loop_ranges(hybrid, ibeg=nets, iend=nete)

      if(adiabatic) then
        if((.not. analytic_ic_active()) .and. (runtype == 0)) then
          do ie=nets,nete
            elem(ie)%state%q(:,:,:,1)=0.0_r8
          end do
        end if
      else if(ideal_phys) then
        if(runtype == 0) then
          do ie=nets,nete
            elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)
            elem(ie)%state%ps(:,:,:) =dyn_ps0
            elem(ie)%state%phis(:,:)=0.0_r8
            elem(ie)%state%T(:,:,:,:) =Tinit
            elem(ie)%state%v(:,:,:,:,:) =0.0_r8
            elem(ie)%state%q(:,:,:,1)=0.0_r8
          end do
        end if
      else if(aqua_planet .and. runtype==0)  then
        do ie=nets,nete
          !          elem(ie)%state%lnps(:,:,:) =LOG(dyn_ps0)
          !          elem(ie)%state%ps(:,:,:) =dyn_ps0
          elem(ie)%state%phis(:,:)=0.0_r8
        end do
      end if

      do ie=nets,nete
        elem(ie)%derived%FM=0.0_r8
        elem(ie)%derived%FT=0.0_r8
        elem(ie)%derived%FQ=0.0_r8
        elem(ie)%derived%omega=0.0_r8
      end do

      ! scale PS to achieve prescribed dry mass
      if (runtype == 0) then
        ! new run, scale mass to value given in namelist, if needed
        call prim_set_mass(elem, TimeLevel,hybrid,hvcoord,nets,nete)
      endif

      call prim_init2(elem,fvm,hybrid,nets,nete, TimeLevel, hvcoord)

      ! This subroutine is used to create nc_topo files, if requested
      call nctopo_util_driver(elem,hybrid,nets,nete)
!$OMP END PARALLEL
    end if

    if (inst_index == 1) then
      call write_grid_mapping(par, elem)
    end if

  end subroutine dyn_init

!===============================================================================

  subroutine dyn_run( dyn_state, rc )

    use parallel_mod,     only: par
    use prim_advance_mod, only: calc_tot_energy_dynamics
    use prim_driver_mod,  only: prim_run_subcycle
    use dimensions_mod,   only: nlev,npsq,nelemd,np
    use time_mod,         only: tstep, nsplit, timelevel_qdp
    use hybrid_mod,       only: config_thread_region, get_loop_ranges
    use control_mod,      only: qsplit ,rsplit
    use cam_history,      only: outfld
    use thread_mod,       only: horz_num_threads

    type (dyn_export_t), intent(inout) :: dyn_state   !  container
    integer,             intent(out)   :: rc      ! Return code

    type(hybrid_t) :: hybrid
    integer        :: n
    integer        :: nets, nete, ithr
    integer        :: ie, i, j
    integer        :: n0_qdp

#if 1
    real(r8) :: ps_before(np,np,nelemd),abs_ps_tend(np,np,nelemd),ftmp(npsq)!temporary diagnostic
#endif

    if(iam < par%nprocs) then
       !$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n,ie)
       hybrid = config_thread_region(par,'horizontal')
!JMD       hybrid = config_thread_region(par,'serial')
       call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

       do n=1, nsplit
#if 1
         do ie=nets,nete
           if (n==1) abs_ps_tend(:,:,ie) = 0.0_r8
           ps_before(:,:,ie) = dyn_state%elem(ie)%state%ps(:,:,TimeLevel%n0)
         end do
#endif
         ! forward-in-time RK, with subcycling
         call prim_run_subcycle(dyn_state%elem,dyn_state%fvm,hybrid,nets,nete,&
              tstep, TimeLevel, hvcoord, n)
#if 1
         do ie=nets,nete
           abs_ps_tend(:,:,ie) = abs_ps_tend(:,:,ie)+&
                ABS(ps_before(:,:,ie)-dyn_state%elem(ie)%state%ps(:,:,TimeLevel%n0))/(tstep*qsplit*rsplit)
         end do
         if (n==nsplit) then
           !$OMP MASTER
           do ie=1,nelemd

             abs_ps_tend(:,:,ie)=abs_ps_tend(:,:,ie)/DBLE(nsplit)
             do j=1,np
               do i=1,np
                 ftmp(i+(j-1)*np) = abs_ps_tend(i,j,ie)
               end do
             end do
             call outfld('ABS_dPSdt',ftmp(:),npsq,ie)
           end do
           !$OMP END MASTER
         end if
        !$OMP BARRIER
#endif
       end do

!JMD       hybrid_serial = config_thread_region(par,'serial')
!JMD       call get_loop_ranges(hybrid_serial,ibeg=nets,iend=nete)

       call TimeLevel_Qdp(TimeLevel, qsplit, n0_qdp)!get n0_qdp for diagnostics call
       call calc_tot_energy_dynamics(dyn_state%elem,nets,nete,TimeLevel%n0,n0_qdp,'dBF')
       !$OMP END PARALLEL

     end if
     rc = DYN_RUN_SUCCESS

  end subroutine dyn_run

!===============================================================================

  subroutine dyn_final(DYN_STATE, RESTART_FILE)

    type (elem_state_t), target     :: DYN_STATE
    character(LEN=*)   , intent(IN) :: RESTART_FILE

  end subroutine dyn_final

!===============================================================================

  subroutine read_inidat(dyn_in)
    use shr_sys_mod,         only: shr_sys_flush
    use shr_vmath_mod,       only: shr_vmath_log
    use physconst,           only: pi
    use hycoef,              only: hyai, hybi, ps0
    use dyn_grid,            only: dyn_decomp, pelat_deg, pelon_deg
    use constituents,        only: cnst_name, cnst_read_iv, qmin, pcnst, cnst_type
    use const_init,          only: cnst_init_default
    use cam_control_mod,     only: ideal_phys, aqua_planet
    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id, pertlim
    use cam_history_support, only: max_fieldname_len
    use cam_grid_support,    only: cam_grid_get_local_size, cam_grid_get_gcid
    use cam_grid_support,    only: cam_grid_dimensions, cam_grid_get_dim_names
    use cam_map_utils,       only: iMap

    use parallel_mod,        only: par
    use bndry_mod,           only: bndry_exchangev
    use element_mod,         only: timelevels
    use edgetype_mod,        only: EdgeBuffer_t
    use edge_mod,            only: edgevpack, edgevunpack, InitEdgeBuffer, FreeEdgeBuffer
    use dimensions_mod,      only: nelemd, nlev, nc, np, npsq
    use dimensions_mod,      only: ntrac, qsize, qsize_d, qsize_condensate_loading
    use dimensions_mod,      only: ldry_mass_vertical_coordinates
    use fvm_mapping,         only: dyn2fvm_mass_vars
    use inic_analytic,       only: analytic_ic_active, analytic_ic_set_ic
    use cam_grid_support,    only: cam_grid_id, cam_grid_get_latvals, cam_grid_get_lonvals
    use dyn_tests_utils,     only: vc_moist_pressure, vc_dry_pressure
    use pio,                 only: file_desc_t, pio_seterrorhandling, PIO_BCAST_ERROR

    type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

    type(file_desc_t), pointer                 :: fh_ini, fh_topo

    type(element_t), pointer         :: elem(:)
    real(r8), allocatable            :: dbuf2(:,:)     ! (npsq,nelemd)
    real(r8), allocatable            :: dbuf3(:,:,:)   ! (npsq,nlev,nelemd)
    real(r8), allocatable            :: qtmp(:,:,:,:,:) ! (np,np,nlev,nelemd,n)
    real(r8), allocatable            :: factor_array(:,:,:,:) ! (np,np,nlev,nelemd)
    real(r8),                pointer :: dp_phys(:,:,:)
    real(r8),                pointer :: dp_gll(:,:,:,:)
    logical,  allocatable            :: pmask(:) ! (npsq*nelemd) unique grid vals
    real(r8), allocatable            :: latvals(:)
    real(r8), allocatable            :: lonvals(:)
    integer                          :: ie, k, t
    character(len=max_fieldname_len) :: fieldname
    logical                          :: found, required
    integer                          :: kptr, m_cnst
    type(EdgeBuffer_t)               :: edge
    integer                          :: lsize

    integer(iMap), pointer           :: ldof(:) ! Basic (2D) grid dof

    ! Information on dimension names for this run
    character(len=max_fieldname_len) :: dyn_dimname
    character(len=4)                 :: field_suffix
    integer                          :: ncol_d_dimid

    integer                          :: rndm_seed_sz
    integer, allocatable             :: rndm_seed(:)
    integer                          :: dims(2)
    integer                          :: pio_errtype
    real(r8)                         :: pertval
    integer                          :: i, j, indx
    integer                          :: dyn_cols
    real(r8), parameter              :: rad2deg = 180.0_r8 / pi
    real(r8), parameter              :: deg2rad = pi / 180.0_r8
    character(len=128)               :: errmsg
    character(len=*), parameter      :: subname='READ_INIDAT'

    integer                          :: ioff
    !
    ! fvm vars
    !
    real(r8), allocatable            :: phis_tmp(:,:) ! (npsp,nelemd)
    real(r8), allocatable            :: inv_dp_darea_fvm(:,:,:)
    real(r8)                         :: min_val, max_val

    real(r8)                         :: dp_tmp, pstmp(np,np)

    ! Variables for analytic initial conditions
    integer,  allocatable            :: glob_ind(:)
    integer,  allocatable            :: m_ind(:)
    real(r8), allocatable            :: dbuf4(:,:,:,:)
    integer                          :: vcoord


    nullify(dp_phys)
    nullify(dp_gll)
    nullify(ldof)
    nullify(fh_ini)
    nullify(fh_topo)

    if(iam < par%nprocs) then
      elem=> dyn_in%elem
    else
      nullify(elem)
    end if

    lsize = cam_grid_get_local_size(dyn_decomp)

    if (lsize /= (np*np*nelemd)) then
      call endrun(trim(subname)//': mismatch in local input array size')
    end if

    if (iam < par%nprocs) then
      if(elem(1)%idxP%NumUniquePts <=0 .or. elem(1)%idxP%NumUniquePts > np*np) then
        write(iulog,*)  elem(1)%idxP%NumUniquePts
        call endrun(trim(subname)//': invalid idxP%NumUniquePts')
      end if
    end if

    ! What are the dimension names for the input grid?
    call cam_grid_get_dim_names('GLL', dyn_dimname, fieldname)
    if (trim(dyn_dimname) /= trim(fieldname)) then
      call endrun(trim(subname)//': dynamics dimensions not equal')
    end if
    call cam_grid_dimensions('GLL', dims, indx)
    if (indx /= 1) then
      call endrun(trim(subname)//': dynamics grid has incorrect dims')
    end if
    dyn_cols = dims(1)
    ! Allocate dynamics input buffers
    allocate(phis_tmp(npsq,nelemd))

    ! We need to know which elements are active
    if (associated(ldof)) then
      call endrun(trim(subname)//': ldof should not be associated')
    end if
    call cam_grid_get_gcid(dyn_decomp, ldof)

    allocate(latvals(size(pelat_deg)))
    allocate(lonvals(size(pelon_deg)))
    latvals(:) = pelat_deg(:)*deg2rad
    lonvals(:) = pelon_deg(:)*deg2rad

    allocate(pmask(npsq*nelemd))
    pmask(:) = (ldof /= 0)

    if (analytic_ic_active()) then
      if (ldry_mass_vertical_coordinates) then
        vcoord = vc_dry_pressure
      else
        vcoord = vc_moist_pressure
      endif
      allocate(glob_ind(npsq * nelemd))
      j = 1
      do ie = 1, nelemd
        do i = 1, npsq
          ! Create a global(ish) column index
          glob_ind(j) = elem(ie)%GlobalId
          j = j + 1
        end do
      end do

      ! First, initialize all the variables, then assign
      allocate(dbuf4(npsq, nlev, nelemd, (qsize + 4)))
      dbuf4 = 0.0_r8
      allocate(m_ind(qsize))
      do m_cnst = 1, qsize
        m_ind(m_cnst) = m_cnst
      end do
      call analytic_ic_set_ic(vcoord, latvals, lonvals, glob_ind,              &
           PS=dbuf4(:,1,:,(qsize+1)), PHIS=phis_tmp, U=dbuf4(:,:,:,(qsize+2)), &
           V=dbuf4(:,:,:,(qsize+3)), T=dbuf4(:,:,:,(qsize+4)),                 &
           Q=dbuf4(:,:,:,1:qsize), m_cnst=m_ind, mask=pmask(:))
      deallocate(m_ind)
      deallocate(glob_ind)
      do ie=1,nelemd
        indx = 1
        do j = 1, np
          do i = 1, np
            ! PS
            elem(ie)%state%ps(i,j,1) = dbuf4(indx, 1, ie, (qsize+1))
            ! PHIS is done assigned later
            ! U
            elem(ie)%state%v(i,j,1,:,1) = dbuf4(indx, :, ie, (qsize+2))
            ! V
            elem(ie)%state%v(i,j,2,:,1) = dbuf4(indx, :, ie, (qsize+3))
            ! T
            elem(ie)%state%T(i,j,:,1) = dbuf4(indx, :, ie, (qsize+4))
            indx = indx + 1
          end do
        end do
      end do
      !
      ! tracers
      !
      do m_cnst = 1, qsize
        do ie = 1, nelemd
          elem(ie)%state%Q(:,:,:,m_cnst) = 0.0_r8
          indx = 1
          do j = 1, np
            do i = 1, np
              elem(ie)%state%Q(i, j, :, m_cnst) = dbuf4(indx, :, ie, m_cnst)
              indx = indx + 1
            end do
          end do
        end do
      end do
      deallocate(dbuf4)
      ! Note that fvm tracers are initialized below
    else
      ! Find out what is in the input file
      fh_ini  => initial_file_get_id()
      fh_topo => topo_file_get_id()

      if (associated(fh_ini)) then
        ! We need to be able to see the PIO return values
        call pio_seterrorhandling(fh_ini, PIO_BCAST_ERROR, pio_errtype)
      end if

      call init_ncid_data(fh_ini, dyn_dimname, dyn_cols, ncol_d_dimid, field_suffix)

      allocate(dbuf2(npsq, nelemd))
      allocate(dbuf3(npsq, nlev, nelemd))

      !! Read 2-D field
      fieldname = 'PS'
      if (dyn_field_exists(fh_ini, trim(fieldname), field_suffix)) then
        call read_dyn_var(trim(fieldname), fh_ini, dyn_dimname, field_suffix, dbuf2)
      else
        call endrun(trim(subname)//': PS must be on GLL grid')
      end if

      if(minval(dbuf2, mask=reshape(pmask, (/npsq,nelemd/))) < 10000._r8) then
        call endrun(trim(subname)//': Problem reading ps field -- bad values')
      end if

      do ie=1,nelemd
        indx = 1
        do j = 1, np
          do i = 1, np
            elem(ie)%state%ps(i,j,1) = dbuf2(indx,ie)
            indx = indx + 1
          end do
        end do
      end do

      !! Read in 3-D fields

      if (dyn_field_exists(fh_ini, 'U', field_suffix)) then
        call read_dyn_var('U', fh_ini, dyn_dimname, field_suffix, dbuf3)
      else
        call endrun(trim(subname)//': U'//trim(field_suffix)//' not found')
      end if
      do ie = 1, nelemd
        elem(ie)%state%v = 0.0_r8
        indx = 1
        do j = 1, np
          do i = 1, np
            elem(ie)%state%v(i,j,1,:,1) = dbuf3(indx,:,ie)
            indx = indx + 1
          end do
        end do
      end do

      if (dyn_field_exists(fh_ini, 'V', field_suffix)) then
        call read_dyn_var('V', fh_ini, dyn_dimname, field_suffix, dbuf3)
      else
        call endrun(trim(subname)//': V'//trim(field_suffix)//' not found')
      end if
      do ie = 1, nelemd
        indx = 1
        do j = 1, np
          do i = 1, np
            elem(ie)%state%v(i,j,2,:,1) = dbuf3(indx,:,ie)
            indx = indx + 1
          end do
        end do
      end do

      if (dyn_field_exists(fh_ini, 'T', field_suffix)) then
        call read_dyn_var('T', fh_ini, dyn_dimname, field_suffix, dbuf3)
      else
        call endrun(trim(subname)//': T'//trim(field_suffix)//' not found')
      end if
      do ie=1,nelemd
        elem(ie)%state%T=0.0_r8
        indx = 1
        do j = 1, np
          do i = 1, np
            elem(ie)%state%T(i,j,:,1) = dbuf3(indx,:,ie)
            indx = indx + 1
          end do
        end do
      end do

      if (pertlim .ne. 0.0_r8) then
        if(masterproc) then
          write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
               'by +/- ', pertlim, ' to initial temperature field'
        end if

        call random_seed(size=rndm_seed_sz)
        allocate(rndm_seed(rndm_seed_sz))

        do ie=1,nelemd
          ! seed random number generator based on element ID
          ! (possibly include a flag to allow clock-based random seeding)
          rndm_seed = elem(ie)%GlobalId
          call random_seed(put=rndm_seed)
          do i=1,np
            do j=1,np
              do k=1,nlev
                call random_number(pertval)
                pertval = 2.0_r8*pertlim*(0.5_r8 - pertval)
                elem(ie)%state%T(i,j,k,1) = elem(ie)%state%T(i,j,k,1)*(1.0_r8 + pertval)
              end do
            end do
          end do
        end do

        deallocate(rndm_seed)
      end if

      ! qmin = 1e-12,0,0

      ! Read in or cold-initialize all the tracer fields
      ! Data is read in on the GLL grid
      ! Both GLL and FVM tracer fields are initialized based on the
      ! dimension qsize or ntrac for GLL or FVM tracers respectively.
      ! Data is only read in on GLL so if FVM tracers are active,
      ! interpolation is performed.
      if (ntrac > qsize) then
        if (ntrac < pcnst) then
          write(errmsg, '(a,3(i0,a))') ': ntrac (',ntrac,') > qsize (',qsize, &
               ') but < pcnst (',pcnst,')'
          call endrun(trim(subname)//errmsg)
        end if
        allocate(qtmp(np,np,nlev,nelemd,ntrac-qsize))
      else if (qsize < pcnst) then
        write(errmsg, '(a,2(i0,a))') ': qsize (',qsize,') < pcnst (',pcnst,')'
        call endrun(trim(subname)//errmsg)
      end if
      do m_cnst = 1, pcnst

        found = .false.

        if(cnst_read_iv(m_cnst)) then
          required = .false.
          found = dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)),            &
               field_suffix, field_required_in=required)
        end if

        if(found) then
          call read_dyn_var(trim(cnst_name(m_cnst)), fh_ini, dyn_dimname, field_suffix, dbuf3)
        else
          call cnst_init_default(m_cnst, latvals, lonvals, dbuf3, pmask)
        end if

        ! Here, we initialize the Eulerian tracers
        ! FVM-only tracers are stored (in qtmp) for later use.
        do ie = 1, nelemd
          ! Copy tracers defined on GLL grid into Eulerian array
          ! Make sure tracers have at least minimum value
          do k=1, nlev
            indx = 1
            do j = 1, np
              do i = 1, np
                ! Zero out the dbuf3 values which might have been set
                ! erroneously by <param>_init_const
                if (pmask(((ie - 1) * npsq) + indx)) then
                  if (qsize >= m_cnst) then
                    elem(ie)%state%Q(i,j,k,m_cnst) = max(qmin(m_cnst),dbuf3(indx,k,ie))
                  else
                    qtmp(i,j, k, nelemd, m_cnst-qsize) = max(qmin(m_cnst),dbuf3(indx,k,ie))
                  end if
                else
                  if (qsize >= m_cnst) then
                    elem(ie)%state%Q(i,j,k,m_cnst) = 0.0_r8
                  else
                    qtmp(i,j, k, nelemd, m_cnst-qsize) = 0.0_r8
                  end if
                end if
                indx = indx + 1
              end do
            end do
          end do
        end do
      end do ! pcnst

      ! Process topography
      fieldname = 'PHIS'
      if (ideal_phys .or. aqua_planet .or. (.not. associated(fh_topo))) then
        phis_tmp = 0._r8
      else if (dyn_field_exists(fh_topo, trim(fieldname), field_suffix)) then
        call read_dyn_var(trim(fieldname), fh_topo, dyn_dimname, field_suffix, phis_tmp)
      else
        call endrun('Could not find PHIS field on input datafile')
      end if

      ! Cleanup
      deallocate(dbuf2)
      deallocate(dbuf3)
    end if ! analytic_ic_active

    ! Cleanup
    if (allocated(pmask)) then
      deallocate(pmask)
    end if ! analytic_ic_active
    if (allocated(latvals)) then
      deallocate(latvals)
    end if
    if (allocated(lonvals)) then
      deallocate(lonvals)
    end if
    if (associated(ldof)) then
      deallocate(ldof)
      nullify(ldof)
    end if

    ! Process phis_tmp
    do ie = 1, nelemd
      elem(ie)%state%phis = 0.0_r8
      indx = 1
      do j = 1, np
        do i = 1, np
          elem(ie)%state%phis(i,j) = phis_tmp(indx, ie)
          indx = indx + 1
        end do
      end do
    end do

    ! once we've read or initialized all the fields we do a boundary exchange to
    ! update the redundent columns in the dynamics
    if(iam < par%nprocs) then
      call initEdgeBuffer(par, edge, elem, ((3+max(qsize,ntrac)) * nlev) + 2)
    end if
    do ie = 1, nelemd
      kptr = 0
      call edgeVpack(edge, elem(ie)%state%ps,1,kptr,ie)
      kptr = kptr + 1
      call edgeVpack(edge, elem(ie)%state%phis,1,kptr,ie)
      kptr = kptr + 1
      call edgeVpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,ie)
      kptr = kptr + (2 * nlev)
      call edgeVpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,ie)
      kptr = kptr + nlev
      call edgeVpack(edge, elem(ie)%state%Q(:,:,:,1:qsize),nlev*qsize,kptr,ie)
      if (ntrac > qsize) then
        kptr = kptr + (qsize * nlev)
        call edgeVpack(edge, qtmp,nlev*(ntrac-qsize),kptr,ie)
      end if
    end do
    if(iam < par%nprocs) then
      call bndry_exchangeV(par,edge,location='read_inidat')
    end if
    do ie = 1, nelemd
      kptr = 0
      call edgeVunpack(edge, elem(ie)%state%ps,1,kptr,ie)
      kptr = kptr + 1
      call edgeVunpack(edge, elem(ie)%state%phis,1,kptr,ie)
      kptr = kptr + 1
      call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,ie)
      kptr = kptr + (2 * nlev)
      call edgeVunpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,ie)
      kptr = kptr + nlev
      call edgeVunpack(edge, elem(ie)%state%Q(:,:,:,1:qsize),nlev*qsize,kptr,ie)
      if (ntrac > qsize) then
        kptr = kptr + (qsize * nlev)
        call edgeVunpack(edge, qtmp,nlev*(ntrac-qsize),kptr,ie)
      end if
    end do

    if (ldry_mass_vertical_coordinates) then
      !
      ! convert to dry
      !
      ! (this has to be done after edge-exchange since shared points between elements are only
      !  initialized in one element and not the other!)
      !
      if(par%masterproc) then
        write(iulog,*) 'Convert specific/wet mixing ratios to dry'
      end if
      allocate(factor_array(np,np,nlev, nelemd))
      do ie=1,nelemd
        do k=1,nlev
          do j = 1, np
            do i = 1, np
              factor_array(i,j,k,ie) = 1.0_r8/(1.0_r8-SUM(elem(ie)%state%Q(i,j,k,1:qsize_condensate_loading)))
            end do
          end do
        end do
      end do
      do m_cnst = 1,pcnst
        if (cnst_type(m_cnst).eq.'wet') then
          do ie=1,nelemd
            !
            ! convert wet mixing ratio to dry
            !
            do k=1,nlev
              do j = 1, np
                do i = 1, np
                  if (qsize >= m_cnst) then
                    elem(ie)%state%Q(i,j,k,m_cnst)=elem(ie)%state%Q(i,j,k,m_cnst)*factor_array(i,j,k,ie)
                  else
                    qtmp(i,j,k,ie,m_cnst-qsize) = qtmp(i,j,k,ie,m_cnst-qsize) * factor_array(i,j,k,ie)
                  end if
                  if (.not. analytic_ic_active()) then
                    if (qsize >= m_cnst) then
                      elem(ie)%state%Q(i,j,k,m_cnst) = max(qmin(m_cnst),elem(ie)%state%Q(i,j,k,m_cnst))
                    else
                      qtmp(i,j,k,ie,m_cnst-qsize) = max(qmin(m_cnst),qtmp(i,j,k,ie,m_cnst-qsize))
                    end if
                  end if
                end do
              end do
            end do
          end do
        end if
      end do
      !
      ! initialize dp3d and qdp
      !
      do ie=1,nelemd
        do k=1,nlev
          do j = 1, np
            do i = 1, np
              factor_array(i,j,k,ie) = 1.0_r8/(1.0_r8+SUM(elem(ie)%state%Q(i,j,k,1:qsize_condensate_loading)))
            end do
          end do
        end do
      end do
      do ie=1,nelemd
        pstmp = elem(ie)%state%ps(:,:,1)
        elem(ie)%state%ps(:,:,1) = hyai(1)*ps0
        do k=1,nlev
          do j = 1,np
            do i = 1,np
              dp_tmp = ((hyai(k+1) - hyai(k))*ps0)         +                &
                   ((hybi(k+1) - hybi(k))*pstmp(i,j))
              if (.not. analytic_ic_active()) then
                !
                ! if analytic_ic then the surface pressure is already dry
                ! (note that it is not correct to convert to moist pressure
                ! in analytic_ic and not have the #ifndef statement here
                ! since the dry levels are in a different location than
                ! what is obtained from algorithm below)
                !
                dp_tmp = dp_tmp*factor_array(i,j,k,ie)
              end if
              elem(ie)%state%dp3d(i,j,k,:) = dp_tmp
              !
              ! compute dry surface pressure; note that at this point
              !
              !     dp3d .NE. ((hyai(k+1) - hyai(k))*ps0)+((hybi(k+1) - hybi(k))*ps(i,j))
              !
              !
              elem(ie)%state%ps(i,j,1) = elem(ie)%state%ps(i,j,1)+elem(ie)%state%dp3d(i,j,k,1)
              do m_cnst = 1, qsize
                elem(ie)%state%Qdp(i,j,k,m_cnst,:)=dp_tmp*elem(ie)%state%q(i,j,k,m_cnst)
              end do
            end do
          end do
        end do
      end do
      deallocate(factor_array)
    else
      if (qsize_condensate_loading /= 1) then
        call endrun('qsize_condensate_loading must be 1 for wet vertical coordinates')
      end if
      do m_cnst = 1, pcnst
        if (cnst_type(m_cnst).eq.'dry') then
          do ie = 1, nelemd
            !
            ! convert dry mixing ratio to wet - note that qv is always wet
            !
            do k = 1, nlev
              do j = 1, np
                do i = 1, np
                  if (qsize >= m_cnst) then
                    elem(ie)%state%Q(i,j,k,m_cnst)=elem(ie)%state%Q(i,j,k,m_cnst)*&
                         (1.0_r8-elem(ie)%state%Q(i,j,k,1))
                    elem(ie)%state%Q(i,j,k,m_cnst) = max(qmin(m_cnst),elem(ie)%state%Q(i,j,k,m_cnst))
                  else
                    qtmp(i,j,k,ie,m_cnst-qsize) = qtmp(i,j,k,ie,m_cnst-qsize)*&
                         (1.0_r8-elem(ie)%state%Q(i,j,k,1))
                    qtmp(i,j,k,ie,m_cnst-qsize) = max(qmin(m_cnst),qtmp(i,j,k,ie,m_cnst-qsize))
                  end if
                end do
              end do
            end do
          end do
        end if
      end do
      !
      ! initialize dp3d and qdp
      !
      do ie = 1, nelemd
        do k = 1, nlev
          do j = 1, np
            do i = 1, np
              dp_tmp = ((hyai(k+1) - hyai(k))*ps0)         +                &
                   ((hybi(k+1) - hybi(k))*elem(ie)%state%ps(i,j,1))
              elem(ie)%state%dp3d(i,j,k,:) = dp_tmp
              do m_cnst = 1, qsize
                elem(ie)%state%Qdp(i,j,k,m_cnst,:)=dp_tmp*elem(ie)%state%q(i,j,k,m_cnst)
              end do
            end do
          end do
        end do
      end do
    end if ! ldry_mass_vertical_coordinates
    !
    ! interpolate fvm tracers and fvm pressure variables
    !
    if (ntrac>0) then
      if(par%masterproc) then
        write(iulog,*) 'Initializing dp_fvm from spectral element dp'
      end if
      do ie = 1, nelemd
        !
        ! note that the area over fvm cells as computed from subcell_integration is up to 1.0E-6
        ! different than the areas (exact) computed by CSLAM
        !
        do k = 1, nlev
          dyn_in%fvm(ie)%dp_ref(k)         = ( hyai(k+1) - hyai(k) )*ps0 + ( hybi(k+1) - hybi(k) )*ps0
          dyn_in%fvm(ie)%dp_ref_inverse(k) = 1.0_r8/dyn_in%fvm(ie)%dp_ref(k)
        end do
        ! Map the constituents which are also to be transported by dycore
        if (analytic_ic_active()) then
          lsize = 1
        else
          lsize = qsize
        end if
        call dyn2fvm_mass_vars(elem(ie)%state%dp3d(:,:,:,1),elem(ie)%state%ps(:,:,1),&
             elem(ie)%state%q(:,:,:,1:lsize),&
             dyn_in%fvm(ie)%dp_fvm(1:nc,1:nc,:,1),dyn_in%fvm(ie)%psC(1:nc,1:nc),&
             dyn_in%fvm(ie)%c(1:nc,1:nc,:,1:lsize,1),&
             lsize,elem(ie)%metdet,dyn_in%fvm(ie)%inv_se_area_sphere(1:nc,1:nc),.false.)
        ! Map the fvm-only constituents (when not using analytic initial cond.)
        if ((.not. analytic_ic_active()) .and. (ntrac > lsize)) then
          call dyn2fvm_mass_vars(elem(ie)%state%dp3d(:,:,:,1),elem(ie)%state%ps(:,:,1),&
               qtmp(:,:,:,ie,1:ntrac-qsize),&
               dyn_in%fvm(ie)%dp_fvm(1:nc,1:nc,:,1),dyn_in%fvm(ie)%psC(1:nc,1:nc),&
               dyn_in%fvm(ie)%c(1:nc,1:nc,:,qsize+1:ntrac,1),&
               ntrac-qsize,elem(ie)%metdet,dyn_in%fvm(ie)%inv_se_area_sphere(1:nc,1:nc),.false.)
        end if

        dyn_in%fvm(ie)%dp_fvm(1:nc,1:nc,:,2)    = dyn_in%fvm(ie)%dp_fvm(1:nc,1:nc,:,1)
        dyn_in%fvm(ie)%c(1:nc,1:nc,:,1:ntrac,2) = dyn_in%fvm(ie)%c(1:nc,1:nc,:,1:ntrac,1)
      end do

      if (analytic_ic_active()) then
        !
        ! initialize tracers
        !
        allocate(latvals(nc*nc*nelemd))
        allocate(lonvals(nc*nc*nelemd))
        indx = 1
        do ie = 1, nelemd
          do j = 1, nc
            do i = 1, nc
              latvals(indx) = dyn_in%elem(ie)%spherep(i,j)%lat
              lonvals(indx) = dyn_in%elem(ie)%spherep(i,j)%lon
              indx = indx + 1
            end do
          end do
        end do

        allocate(pmask(nc*nc*nelemd))
        pmask(:) = .true.

        allocate(dbuf4(nc*nc, nlev, nelemd, ntrac))
        allocate(m_ind(ntrac))
        allocate(glob_ind(nc*nc*nelemd))
        j = 1
        do ie = 1, nelemd
          do i = 1, nc*nc
            ! Create a global(ish) column index
            glob_ind(j) = elem(ie)%GlobalId
            j = j + 1
          end do
        end do

        dbuf4 = 0.0_r8
        do m_cnst = 1, ntrac
          m_ind(m_cnst) = m_cnst
        end do
        call analytic_ic_set_ic(vcoord, latvals, lonvals, glob_ind, Q=dbuf4, m_cnst=m_ind, mask=pmask)
        !
        ! it is more balanced to use dyn2fvm for Q than to use the "analytical" value
        ! on the fvm grid
        !
        do m_cnst = 2, ntrac
          do ie = 1, nelemd
            indx = 1
            do j = 1, nc
              do i = 1, nc
                dyn_in%fvm(ie)%c(i,j,:,m_cnst,1) = dbuf4(indx, :, ie, m_cnst)
                dyn_in%fvm(ie)%c(i,j,:,m_cnst,2) = dbuf4(indx, :, ie, m_cnst)
                indx = indx + 1
              end do
            end do
          end do
        end do
        deallocate(dbuf4)
        deallocate(m_ind)
        deallocate(latvals)
        deallocate(lonvals)
        deallocate(glob_ind)
        deallocate(pmask)
      end if

      if(par%masterproc) then
        write(iulog,*) 'FVM tracers, FVM pressure variables and se_area_sphere initialized.'
      end if
    end if !ntrac>0
    ! Cleanup
    if (allocated(qtmp)) then
      deallocate(qtmp)
    end if

    do ie = 1, nelemd
      do t = 2, timelevels
        elem(ie)%state%ps(:,:,t)=elem(ie)%state%ps(:,:,1)
        elem(ie)%state%v(:,:,:,:,t)=elem(ie)%state%v(:,:,:,:,1)
        elem(ie)%state%T(:,:,:,t)=elem(ie)%state%T(:,:,:,1)
      end do
      call shr_vmath_log(elem(ie)%state%ps,elem(ie)%state%lnps,size(elem(ie)%state%lnps))
    end do

    if(iam < par%nprocs) then
      call FreeEdgeBuffer(edge)
    end if

    ! Cleanup
    if (associated(fh_ini)) then
      ! Put the error handling back the way it was
      call pio_seterrorhandling(fh_ini, pio_errtype)
    end if
  end subroutine read_inidat

  subroutine init_ncid_data(ncid, dyn_dimname, dyn_cols, ncol_d_dimid, field_suffix)
    use shr_sys_mod,      only: shr_sys_flush
    use pio,              only: pio_inq_dimid, pio_inq_dimlen, PIO_NOERR
    use pio,              only: file_desc_t
    use dimensions_mod,   only: np, nelemd
    use dyn_grid,         only: dyn_decomp
    use cam_grid_support, only: cam_grid_get_local_size

    type(file_desc_t), intent(in)    :: ncid
    character(len=*),  intent(in)    :: dyn_dimname
    integer,           intent(in)    :: dyn_cols     ! Global # of dyn cols
    integer,           intent(out)   :: ncol_d_dimid ! Dynamics dimension ID
    character(len=*),  intent(inout) :: field_suffix

    ! Local variables
    integer                          :: ncol_file, ierr
    integer                          :: col_size
    integer                          :: ncol_dimid
    character(len=*), parameter      :: subname = 'INIT_NCID_DATA'

    field_suffix = ''
    ncol_dimid = -1

    ! Check to see if we have the dyn_decomp dimension
    ierr = pio_inq_dimid(ncid, trim(dyn_dimname), ncol_d_dimid)
    if (ierr /= PIO_NOERR) then
      ncol_d_dimid = -1
    end if
    ! As a backup, see if we have an 'ncol' dimension (for legacy files)
    if ((ncol_d_dimid == -1) .and. (trim(dyn_dimname) /= 'ncol')) then
      ierr = pio_inq_dimid(ncid, 'ncol', ncol_dimid)
      if (ierr /= PIO_NOERR) then
        ncol_dimid = -1
      end if
    end if

    if (ncol_d_dimid /= -1) then
      ! We have the dyn_decomp dimension. Check its length
      ierr = pio_inq_dimlen(ncid, ncol_d_dimid, ncol_file)  ! value in file
      if (ncol_file /= dyn_cols) then
        if (masterproc) then
          write(iulog,*) trim(dyn_dimname),' dimension in IC file =',ncol_file
          write(iulog,*) 'Expected dimension size     =', dyn_cols
          call shr_sys_flush(iulog)
        end if
        call endrun(trim(subname)//": IC file dynamics '"//trim(dyn_dimname)//"' dimension mismatch")
      end if
      col_size = cam_grid_get_local_size(dyn_decomp)

      if (col_size /= (np*np*nelemd)) then
        call endrun(trim(subname)//': mismatch in local GLL input array size')
      end if
      if (trim(dyn_dimname) == 'ncol') then
        ! This is (should be) a legacy IC file
        field_suffix = ''
      else
        ! We have a special ncol_d dimension so variables will have gll suffix
        field_suffix = 'gll'
      end if
    else if (ncol_dimid /= -1) then
      ! No specific dynamics dimension, try 'ncol'
      ierr = pio_inq_dimlen(ncid, ncol_dimid, ncol_file)  ! value in file
      ! If we are here, there was no 'ncol_d' so 'ncol' should be the GLL grid
      if (ncol_file /= dyn_cols) then
        ! This is bad. None of the dimensions match up
        if (masterproc) then
          write(iulog,*) trim(dyn_dimname), ' dimension does not exist in IC but ncol does not have an appropriate size'
          call shr_sys_flush(iulog)
        end if
        call endrun(trim(subname)//": IC file 'ncol' dimension mismatch")
      end if
      ! ncol is a valid dynamics dimension, this is a legacy file
      ncol_d_dimid = ncol_dimid
      field_suffix = '' ! redundant but just making it clear
    else
      ! Mayday! We don't have any appropriate dimension
      call endrun(trim(subname)//": No dynamics dimension found in IC file")
    end if

  end subroutine init_ncid_data

  logical function dyn_field_exists(ncid, fieldname_in, field_suffix, field_required_in)
    use pio,            only: file_desc_t, var_desc_t, PIO_inq_varid
    use pio,            only: PIO_NOERR

    type(file_desc_t), intent(in)    :: ncid
    character(len=*),  intent(in)    :: fieldname_in
    character(len=*),  intent(in)    :: field_suffix
    logical, optional, intent(inout) :: field_required_in

    ! Local variables
    character(len=40)        :: fieldname
    logical                  :: found
    logical                  :: field_required
    integer                  :: ret
    type(var_desc_t)         :: varid
    character(len=128)       :: errormsg

    if (present(field_required_in)) then
      field_required = field_required_in
    else
      field_required = .true.
    end if

    fieldname = trim(fieldname_in)//trim(field_suffix)
    ret = PIO_inq_varid(ncid, trim(fieldname), varid)
    found = (ret == PIO_NOERR)
    if (.not. found) then
      if (field_required) then
        write(errormsg, *) trim(fieldname),' was not present in the input file.'
        call endrun('DYN_FIELD_EXISTS: '//errormsg)
      end if
    end if

    ! In case they passed this in as false, we set it to found
    if (present(field_required_in)) then
      field_required_in = found
    end if

    dyn_field_exists = found

  end function dyn_field_exists

  subroutine read_dyn_field_2d(fieldname_in, ncid, dimname, field_suffix, buffer)
    use pio,                 only: file_desc_t
    use dimensions_mod,      only: npsq, nelemd
    use ncdio_atm,           only: infld

    ! Dummy arguments
    character(len=*),  intent(in)    :: fieldname_in
    type(file_desc_t), intent(inout) :: ncid
    character(len=*),  intent(in)    :: dimname
    character(len=*),  intent(in)    :: field_suffix
    real(r8),          intent(inout) :: buffer(:, :)

    ! Local variables
    character(len=40)        :: fieldname
    logical                  :: found

    fieldname = trim(fieldname_in)//trim(field_suffix)

    buffer = 0.0_r8
    call infld(trim(fieldname), ncid, dimname, 1, npsq, 1, nelemd, buffer,    &
         found, gridname='GLL')
    if(.not. found) then
      call endrun('READ_DYN_FIELD_2D: Could not find '//trim(fieldname)//' field on input datafile')
    end if
    if (any(buffer(:,:) /= buffer(:,:))) then
      call endrun('READ_DYN_FIELD_2D: NaN found in field '//trim(fieldname))
    end if

  end subroutine read_dyn_field_2d

  subroutine read_dyn_field_3d(fieldname_in, ncid, dimname, field_suffix, buffer)
    use pio,                 only: file_desc_t
    use dimensions_mod,      only: npsq, nelemd, nlev
    use ncdio_atm,           only: infld

    ! Dummy arguments
    character(len=*),  intent(in)    :: fieldname_in
    type(file_desc_t), intent(inout) :: ncid
    character(len=*),  intent(in)    :: dimname
    character(len=*),  intent(in)    :: field_suffix
    real(r8),          intent(inout) :: buffer(:,:,:)

    ! Local variables
    character(len=40)        :: fieldname
    logical                  :: found

    fieldname = trim(fieldname_in)//trim(field_suffix)

    buffer = 0.0_r8
    call infld(trim(fieldname), ncid, dimname, 'lev',  1, npsq, 1, nlev,      &
         1, nelemd, buffer, found, gridname='GLL')
    if(.not. found) then
      call endrun('READ_DYN_FIELD_3D: Could not find '//trim(fieldname)//' field on input datafile')
    end if
    if (any(buffer(:,:,:) /= buffer(:,:,:))) then
      call endrun('READ_DYN_FIELD_3D: NaN found in field '//trim(fieldname))
    end if

  end subroutine read_dyn_field_3d

!===============================================================================

  subroutine write_grid_mapping(par, elem)
    use parallel_mod,     only: parallel_t
    use element_mod, only: element_t
    use cam_pio_utils, only: cam_pio_createfile, pio_subsystem
    use pio, only: file_desc_t, pio_def_dim, var_desc_t, pio_int, pio_def_var, &
         pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, pio_write_darray, &
         pio_freedecomp
    use dimensions_mod, only: np, nelem, nelemd
    use dof_mod, only: createmetadata

    type(parallel_t) :: par
    type(element_t) :: elem(:)
    type(file_desc_t) :: nc
    type(var_desc_t) :: vid
    type(io_desc_t) :: iodesc
    integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
    integer, parameter :: npm12 = (np-1)*(np-1)
    integer :: subelement_corners(npm12*nelemd,4)
    integer :: dof(npm12*nelemd*4)


    ! Create a CS grid mapping file for postprocessing tools

       ! write meta data for physics on GLL nodes
       call cam_pio_createfile(nc, 'SEMapping.nc', 0)

       ierr = pio_def_dim(nc, 'ncenters', npm12*nelem, dim1)
       ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
       ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/),vid)

       ierr = pio_enddef(nc)
       call createmetadata(par, elem, subelement_corners)

       jj=0
       do cc = 0, 3
          do ie = 1, nelemd
             base = ((elem(ie)%globalid-1)+cc*nelem)*npm12
             ii=0
             do j = 1, np-1
                do i = 1, np-1
                   ii=ii+1
                   jj=jj+1
                   dof(jj) = base+ii
                end do
             end do
          end do
       end do

       call pio_initdecomp(pio_subsystem, pio_int, (/nelem*npm12,4/), dof, iodesc)

       call pio_write_darray(nc, vid, iodesc, reshape(subelement_corners,(/nelemd*npm12*4/)), ierr)

       call pio_freedecomp(nc, iodesc)

       call pio_closefile(nc)

  end subroutine write_grid_mapping

!===============================================================================

end module dyn_comp
