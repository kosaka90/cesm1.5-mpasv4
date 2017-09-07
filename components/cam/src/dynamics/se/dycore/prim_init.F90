#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_init

  use kinds, only : real_kind, longdouble_kind
  use dimensions_mod, only : nc
  use reduction_mod, only : reductionbuffer_ordered_1d_t
  use quadrature_mod, only : quadrature_t, gausslobatto

  implicit none
  private
  save

  public :: prim_init1

  real(kind=longdouble_kind), public :: fvm_corners(nc+1) ! fvm cell corners on reference element
  real(kind=longdouble_kind), public :: fvm_points(nc)    ! fvm cell centers on reference element

  type (quadrature_t), public               :: gp    ! element GLL points
  type (ReductionBuffer_ordered_1d_t)       :: red   ! reduction buffer (shared)

!=============================================================================!
contains
!=============================================================================!

  subroutine prim_init1(elem, fvm, par, Tl)

    ! --------------------------------
    use cam_logfile, only : iulog
    ! --------------------------------
    use thread_mod, only : max_num_threads, omp_get_thread_num, omp_set_num_threads
    ! --------------------------------
    use dimensions_mod, only : np, nlev, nelem, nelemd, nelemdmax
    use dimensions_mod, only : GlobalUniqueCols, fv_nphys,irecons_tracer
    ! --------------------------------
    use control_mod, only : runtype, restartfreq, filter_counter, integration, topology, &
         partmethod, while_iter
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use element_mod, only : element_t, allocate_element_desc
    ! --------------------------------
    use fvm_control_volume_mod, only : fvm_struct, allocate_physgrid_vars
    use fvm_mod               , only : fvm_init1
    ! --------------------------------
    use mesh_mod, only : MeshUseMeshFile
    ! --------------------------------
    use time_mod, only : nmax, time_at, timelevel_init, timelevel_t
    ! --------------------------------
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use derivative_mod, only : allocate_subcell_integration_matrix_cslam,allocate_subcell_integration_matrix_physgrid
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
    ! --------------------------------
    use mesh_mod, only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, &
         MeshCubeElemCount, MeshCubeEdgeCount
    use cube_mod, only : cube_init_atomic, rotation_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, localelemcount, initmetagraph, printmetavertex
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t, allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    ! --------------------------------
    use schedtype_mod, only : schedule
    ! --------------------------------
    use schedule_mod, only : genEdgeSched,  PrintSchedule
    ! --------------------------------
    use prim_advection_mod, only: prim_advec_init1
    ! --------------------------------
    use diffusion_mod, only      : diffusion_init
    ! --------------------------------
    use parallel_mod, only : iam, parallel_t, syncmp, abortmp, global_shared_buf, nrepro_vars
    use parallel_mod, only : mpiinteger_t, mpireal_t, mpi_max, mpi_sum, haltmp
    ! --------------------------------
    use metis_mod, only : genmetispart
    ! --------------------------------
    use spacecurve_mod, only : genspacepart
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! --------------------------------
    use physical_constants, only : dd_pi
    ! --------------------------------
    use bndry_mod, only : sort_neighbor_buffer_mapping
    ! --------------------------------
    use reduction_mod, only : red_min, red_max, red_max_int, red_flops
    use reduction_mod, only : red_sum, red_sum_int, initreductionbuffer
    ! --------------------------------
    use infnan,           only: nan, assignment(=)
    ! --------------------------------
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc

    use fvm_analytic_mod , only : compute_basic_coordinate_vars

    implicit none
    type (element_t), pointer :: elem(:)
    type (fvm_struct), pointer   :: fvm(:)
    type (parallel_t), intent(inout) :: par
    type (timelevel_t), intent(out) :: Tl
    ! Local Variables

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)

    integer :: ii,ie, ith
    integer :: nets, nete
    integer :: nelem_edge,nedge
    integer :: nstep
    integer :: nlyr
    integer :: iMv
    integer :: err, ierr, l, j
    logical, parameter :: Debug = .FALSE.

    real(kind=real_kind), allocatable :: aratio(:,:)
    real(kind=real_kind) :: area(1),xtmp
    character(len=80) rot_type   ! cube edge rotation type

    integer  :: i
    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)

    integer total_nelem
    real(kind=real_kind) :: approx_elements_per_task


    ! ====================================
    ! Set cube edge rotation type for model
    ! unnecessary complication here: all should
    ! be on the same footing. RDL
    ! =====================================
    rot_type="contravariant"

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================

    if (topology=="cube") then

       if (par%masterproc) then
          write(iulog,*)"creating cube topology..."
       end if

       if (MeshUseMeshFile) then
           nelem = MeshCubeElemCount()
           nelem_edge = MeshCubeEdgeCount()
       else
           nelem      = CubeElemCount()
           nelem_edge = CubeEdgeCount()
       end if

       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))

       do j =1,nelem
          call allocate_gridvertex_nbrs(GridVertex(j))
       end do

       if (MeshUseMeshFile) then
           if (par%masterproc) then
               write(iulog,*) "Set up grid vertex from mesh..."
           end if
           call MeshCubeTopology(GridEdge, GridVertex)
       else
           call CubeTopology(GridEdge,GridVertex)
        end if

       if(par%masterproc)       write(iulog,*)"...done."
    end if
    if(par%masterproc) write(iulog,*)"total number of elements nelem = ",nelem

    !DBG if(par%masterproc) call PrintGridVertex(GridVertex)


    if(partmethod .eq. SFCURVE) then
       if(par%masterproc) write(iulog,*)"partitioning graph using SF Curve..."
       call genspacepart(GridEdge,GridVertex)
    else
        if(par%masterproc) write(iulog,*)"partitioning graph using Metis..."
       call genmetispart(GridEdge,GridVertex)
    endif

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
    allocate(MetaVertex(1))
    allocate(Schedule(1))


    nelem_edge=SIZE(GridEdge)

    allocate(TailPartition(nelem_edge))
    allocate(HeadPartition(nelem_edge))
    do i=1,nelem_edge
       TailPartition(i)=GridEdge(i)%tail%processor_number
       HeadPartition(i)=GridEdge(i)%head%processor_number
    enddo

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)


    nelemd = LocalElemCount(MetaVertex(1))
    if(par%masterproc .and. Debug) then 
        call PrintMetaVertex(MetaVertex(1))
    endif

    if(nelemd .le. 0) then
       call abortmp('Not yet ready to handle nelemd = 0 yet' )
       stop
    endif
#ifdef _MPI
    call mpi_allreduce(nelemd,nelemdmax,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
#else
    nelemdmax=nelemd
#endif


    if (nelemd>0) then
       allocate(elem(nelemd))
       call allocate_element_desc(elem)
    endif

    if (fv_nphys>0) then
       allocate(fvm(nelemd))
       call allocate_physgrid_vars(fvm,par)
    else
       ! Even if fvm not needed, still desirable to allocate it as empty
       ! so it can be passed as a (size zero) array rather than pointer.
       allocate(fvm(0))
    end if

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(par, elem,iam,Schedule(1),MetaVertex(1))

    allocate(global_shared_buf(nelemd,nrepro_vars))
    global_shared_buf=0.0_real_kind
    !  nlyr=edge3p1%nlyr
    !  call MessageStats(nlyr)
    !  call testchecksum(par,GridEdge)

    ! ========================================================
    ! load graph information into local element descriptors
    ! ========================================================

    !  do ii=1,nelemd
    !     elem(ii)%vertex = MetaVertex(iam)%members(ii)
    !  enddo

    call syncmp(par)

    ! =================================================================
    ! Set number of domains (for 'decompose') equal to number of threads
    !  for OpenMP across elements, equal to 1 for OpenMP within element
    ! =================================================================

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'init shared boundary_exchange buffers'
    call InitReductionBuffer(red,3*nlev,max_num_threads)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_max_int,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)


    gp=gausslobatto(np)  ! GLL points

    ! fvm nodes are equally spaced in alpha/beta
    ! HOMME with equ-angular gnomonic projection maps alpha/beta space
    ! to the reference element via simple scale + translation
    ! thus, fvm nodes in reference element [-1,1] are a tensor product of
    ! array 'fvm_corners(:)' computed below:
    xtmp=nc
    do i=1,nc+1
       fvm_corners(i)= 2*(i-1)/xtmp - 1  ! [-1,1] including end points
    end do
    do i=1,nc
       fvm_points(i)= ( fvm_corners(i)+fvm_corners(i+1) ) /2
    end do

    if (topology=="cube") then
       if(par%masterproc) write(iulog,*) "initializing cube elements..."
       if (MeshUseMeshFile) then
           call MeshSetCoordinates(elem)
       else
           do ie=1,nelemd
               call set_corner_coordinates(elem(ie))
           end do
           call assign_node_numbers_to_elem(elem, GridVertex)
       end if
       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points)
       enddo
    end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running mass_matrix'
    call mass_matrix(par,elem)
    allocate(aratio(nelemd,1))

    if (topology=="cube") then
       area = 0
       do ie=1,nelemd
          aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
       enddo
       call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
       area(1) = 4*dd_pi/area(1)  ! ratio correction
       deallocate(aratio)
       if (par%masterproc) &
            write(iulog,'(a,f20.17)') " re-initializing cube elements: area correction=",area(1)

       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points,area(1))
          call rotation_init_atomic(elem(ie),rot_type)
       enddo
    end if


    if(par%masterproc) write(iulog,*) 're-running mass_matrix'
    call mass_matrix(par,elem)


    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running global_dof'
    call global_dof(par,elem)

    ! =================================================================
    ! Create Unique Indices
    ! =================================================================

    do ie=1,nelemd
       call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    enddo

    call SetElemOffset(par,elem, GlobalUniqueCols)

    do ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    end do

    !JMD call PrintDofP(elem)
    !JMD call PrintDofV(elem)



    call prim_printstate_init(par)
    ! Initialize output fields for plotting...


    while_iter = 0
    filter_counter = 0

    ! initialize flux terms to 0

    do ie=1,nelemd
       elem(ie)%derived%FM=0.0_real_kind
       elem(ie)%derived%FQ=0.0_real_kind
       elem(ie)%derived%FT=0.0_real_kind
       elem(ie)%derived%pecnd=0.0_real_kind

       elem(ie)%derived%Omega=0
       elem(ie)%state%dp3d=0

       elem(ie)%derived%etadot_prescribed = nan
       elem(ie)%derived%u_met = nan
       elem(ie)%derived%v_met = nan
       elem(ie)%derived%dudt_met = nan
       elem(ie)%derived%dvdt_met = nan
       elem(ie)%derived%T_met = nan
       elem(ie)%derived%dTdt_met = nan
       elem(ie)%derived%ps_met = nan
       elem(ie)%derived%dpsdt_met = nan
       elem(ie)%derived%nudge_factor = nan

       elem(ie)%derived%Utnd=0._real_kind
       elem(ie)%derived%Vtnd=0._real_kind
       elem(ie)%derived%Ttnd=0._real_kind
    enddo


    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================
    !DBG  write(iulog,*) 'prim_init: after call to initRestartFile'

    deallocate(GridEdge)
    do j =1,nelem
       call deallocate_gridvertex_nbrs(GridVertex(j))
    end do
    deallocate(GridVertex)

    do j = 1, MetaVertex(1)%nmembers
       call deallocate_gridvertex_nbrs(MetaVertex(1)%members(j))
    end do
    deallocate(MetaVertex)
    deallocate(TailPartition)
    deallocate(HeadPartition)

    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) then
       write(iulog,*) "Main:max_num_threads=",max_num_threads
    endif


    ith=0
    nets=1
    nete=nelemd
    call Prim_Advec_Init1(par, elem,max_num_threads)
    call diffusion_init(par,elem)
    if (fv_nphys>0) then
      call fvm_init1(par,elem)
    endif

    ! =======================================================
    ! Allocate memory for subcell flux calculations.
    ! =======================================================
    call allocate_subcell_integration_matrix_cslam(np, nc)
    if (fv_nphys>0) &
         call allocate_subcell_integration_matrix_physgrid(np, fv_nphys)

    call TimeLevel_init(tl)

    if (fv_nphys>0) then
      if(par%masterproc) write(iulog,*) 'initialize basic fvm coordinate variables'
      do ie=1,nelemd
        !xxxxx is this needed if ntrac==0??????
        call compute_basic_coordinate_vars(elem(ie),&
             nc,irecons_tracer,fvm(ie)%dalpha,fvm(ie)%dbeta,fvm(ie)%vtx_cart(1:nc,1:nc,:,:),&
             fvm(ie)%center_cart(1:nc,1:nc),fvm(ie)%area_sphere(1:nc,1:nc),&
             fvm(ie)%spherecentroid(1:nc,1:nc,:))
        call compute_basic_coordinate_vars(elem(ie),&
             fv_nphys,irecons_tracer,fvm(ie)%dalpha_physgrid,fvm(ie)%dbeta_physgrid,&
             fvm(ie)%vtx_cart_physgrid   (1:fv_nphys,1:fv_nphys,:,:),&
             fvm(ie)%center_cart_physgrid(1:fv_nphys,1:fv_nphys),&
             fvm(ie)%area_sphere_physgrid(1:fv_nphys,1:fv_nphys),&
             fvm(ie)%spherecentroid_physgrid(1:fv_nphys,1:fv_nphys,:))
      end do
    end if

    if(par%masterproc) write(iulog,*) 'end of prim_init'
  end subroutine prim_init1

!=============================================================================!

end module prim_init
