#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module namelist_mod
  !-----------------
  use kinds, only : real_kind, iulog
  !-----------------
  use params_mod, only : recursive, sfcurve
  !-----------------
  USE shr_string_mod, ONLY : shr_string_toUpper
  !-----------------
  use cube_mod, only : rotate_grid
  use control_mod, only : &
       MAX_STRING_LEN,&
       MAX_FILE_LEN,  &
       partmethod,    &       ! Mesh partitioning method (METIS)
       topology,      &       ! Mesh topology
       test_case,     &       ! test case
       uselapi,       &
       multilevel,    &
       numnodes,      &
       sub_case,      &
       tasknum,       &       ! used dg model in AIX machine
       remapfreq,     &       ! number of steps per remapping call
       remap_type,    &       ! selected remapping option
       statefreq,     &       ! number of steps per printstate call
       restartfreq,   &
       restartfile,   &       ! name of the restart file for INPUT
       restartdir,    &       ! name of the restart directory for OUTPUT
       runtype,       &
       integration,   &       ! integration method
       tracer_advection_formulation, &   ! conservation or non-conservation formulaton
       tstep_type, &
       cubed_sphere_map, &
       qsplit, &
       rsplit, &
       physics, &
       rk_stage_user, &
       LFTfreq,       &
       TRACERADV_TOTAL_DIVERGENCE, &
       TRACERADV_UGRADQ, &
       prescribed_wind, &
       ftype,        &
       energy_fixer,        &
       limiter_option, &
       fine_ne,       &
       max_hypervis_courant, &
       nu,            &
       nu_s,          &
       nu_q,          &
       nu_div,          &
       nu_p,          &
       nu_top,        &
       hypervis_scaling,   & ! use tensor HV instead of scalar coefficient
       disable_diagnostics, & ! Use to disable diagnostics for timing reasons
       psurf_vis,    &
       hypervis_order,    &
       hypervis_power,    &
       hypervis_subcycle, &
       hypervis_subcycle_q, &
       smooth_phis_numcycle, &
       smooth_sgh_numcycle, &
       smooth_phis_nudt, &
       initial_total_mass, &  ! set > 0 to set the initial_total_mass
       u_perturb,     &          ! J&W bareclinic test perturbation size
       columnpackage, &
       filter_type,   &
       transfer_type, &
       filter_freq,   &
       filter_mu,     &
       filter_freq_advection,   &
       filter_mu_advection,     &
       p_bv,          &
       s_bv,          &
       wght_fm,       &
       kcut_fm,       &
       vform,           &
       vfile_mid,       &
       vfile_int,       &
       precon_method, &
       maxits,        &
       tol,           &
       debug_level,   &
       vert_remap_q_alg, &
       tracer_transport_type,            &
       TRACERTRANSPORT_SE_GLL,           &
       TRACERTRANSPORT_CONSISTENT_SE_FVM,&
       tracer_grid_type,                 &
       TRACER_GRIDTYPE_GLL,              &
       TRACER_GRIDTYPE_FVM,              &
       test_cfldep, &
       se_prescribed_wind_2d


  !-----------------
  use thread_mod, only : omp_get_max_threads, max_num_threads, horz_num_threads, vert_num_threads, tracer_num_threads
  !-----------------
  use dimensions_mod, only : ne, np, npdg, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
  !-----------------
  !-----------------
  use parallel_mod, only : parallel_t,  iam, abortmp, &
       partitionfornodes, useframes, mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t, &
       boundaryCommMethod, HME_BNDRY_P2P, HME_BNDRY_A2A, HME_BNDRY_MASHM
  !-----------------
  use cg_mod, only : cg_no_debug
  !-----------------


  use interpolate_mod, only : set_interp_parameter, get_interp_parameter


!=======================================================================================================!
!    Adding for SW DG                                                                                   !
!=======================================================================================================!
#ifdef _SWDG
  ! ------------------------
  use dg_flux_mod, only: riemanntype
  ! ------------------------
  use dg_tests_mod, only : alpha_dg, alphatype
  ! ------------------------
  use dg_sweq_mod, only: stage_rk
  ! ------------------------
  use physical_constants, only: dd_pi
  ! ------------------------
#endif

!=======================================================================================================!
  implicit none
  private
!
! This module should contain no global data and should only be 'use'd to
!    call one of the public interfaces below
!
  public :: homme_set_defaults
  public :: homme_postprocess_namelist

 contains

  ! ============================================
  ! homme_set_defaults:
  !
  !  Set default values for namelist variables
  !
  ! ============================================
  subroutine homme_set_defaults()

    PARTMETHOD                                     = RECURSIVE
    npart                                          = 1
    useframes                                      = 0
    multilevel                                     = 1
    uselapi                                        = .TRUE.
    sub_case                                       = 1
    numnodes                                       = -1
    restartfreq                                    = -100
    restartdir                                     = "./restart/"
    runtype                                        = 0
    statefreq                                      = 1
    remapfreq                                      = 240
    remap_type                                     = "parabolic"
    tasknum                                        =-1
    integration                                    = "explicit"
    columnpackage                                  = "none"
    nu_top                                         = 0
    initial_total_mass                             = 0
    ne                                             = 0
    disable_diagnostics                            = .false.

  end subroutine homme_set_defaults

  subroutine homme_postprocess_namelist(mesh_file, par)
    use mesh_mod,        only: MeshOpen

    ! Dummy arguments
    character(len=*),  intent(in) :: mesh_file
    type (parallel_t), intent(in) :: par

    ! Local variable
    real(kind=real_kind) :: dt_max

    !!XXgoldyXX: Why is npart in HOMME namelist if it is just reset here?
    npart  = par%nprocs


    if(par%masterproc) print *,'omp_get_max_threads() = ',max_num_threads

    if((vert_num_threads > 1) .and. (limiter_option .ne. 8)) then
       if(par%masterproc) print *,'WARNING: vertical threading on supported for limiter_option != 8 '
       vert_num_threads = 1
    endif

    if (ne /= 0) then
      if (mesh_file /= "none" .and. mesh_file /= "/dev/null") then
        write (*,*) "namelist_mod: mesh_file:",trim(mesh_file), &
             " and ne:",ne," are both sepcified in the input file."
        write (*,*) "Specify one or the other, but not both."
        call abortmp("Do not specify ne if using a mesh file input.")
      end if
    end if
    if (par%masterproc) then
      write(iulog,*) "Mesh File:", trim(mesh_file)
    end if
    if (ne == 0) then
      if (par%masterproc) then
        write (iulog,*) "Opening Mesh File:", trim(mesh_file)
      end if
      call set_mesh_dimensions()
      call MeshOpen(mesh_file, par)
    end if

    ! set map
    if (cubed_sphere_map < 0) then
      if (ne == 0) then
        cubed_sphere_map = 2  ! element_local for var-res grids
      else
        cubed_sphere_map = 0  ! default is equi-angle gnomonic
      end if
    end if

    if ((cubed_sphere_map /= 0) .AND.                                         &
        tracer_transport_type .eq. TRACERTRANSPORT_CONSISTENT_SE_FVM) then
      print *,' fvm transport and require equi-angle gnomonic cube sphere mapping.'
      print *,' Set cubed_sphere_map = 0 or comment it out all together.                          '
        call abortmp("Error: fvm transport and cubed_sphere_map>0")
    end if
    if (par%masterproc) write (iulog,*) "Reference element projection: cubed_sphere_map=",cubed_sphere_map

    !logic around different hyperviscosity options
    if (hypervis_power /= 0) then
      if (hypervis_scaling /= 0) then
        print *,'Both hypervis_power and hypervis_scaling are nonzero.'
        print *,'(1) Set hypervis_power=1, hypervis_scaling=0 for HV based on an element area.'
        print *,'(2) Set hypervis_power=0 and hypervis_scaling=1 for HV based on a tensor.'
        print *,'(3) Set hypervis_power=0 and hypervis_scaling=0 for constant HV.'
        call abortmp("Error: hypervis_power>0 and hypervis_scaling>0")
      end if
    end if

    if((prescribed_wind /= 0) .and. (prescribed_wind /= 1))then
      call abortmp('prescribed_wind should be either 0 or 1')
    end if

    ! some default diffusion coefficiets
    if (nu_s < 0) then
      nu_s = nu
    end if
    if (nu_q < 0) then
      nu_q = nu
    end if
    if (nu_div < 0) then
      nu_div = nu
    end if

    if (multilevel <= 0) then
      nmpi_per_node = 1
    end if

    nnodes = npart/nmpi_per_node

    if((numnodes > 0) .and. (multilevel == 1)) then
      nnodes = numnodes
      nmpi_per_node = npart/nnodes
    end if

    ! ====================================================================
    !  Do not perform node level partitioning if you are only on one node
    ! ====================================================================
    ! PartitionForNodes=.FALSE.
    if((nnodes .eq. 1) .and. PartitionForNodes) then
      PartitionForNodes = .FALSE.
    end if

  end subroutine homme_postprocess_namelist
end module namelist_mod
