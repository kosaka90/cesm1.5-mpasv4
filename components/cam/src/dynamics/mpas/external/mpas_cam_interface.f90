module mpas_cam_interface

   use mpas_configure
   use mpas_io
   use mpas_dmpar
   use mpas_constants
   use atm_core
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units, only : stdoutUnit, stderrUnit
   use mpas_framework
   use mpas_kind_types
   use atm_core_interface

   type (core_type), pointer :: corelist => null()
   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain_ptr
   real (kind=RKIND) :: dt_dynamics, dt_physics, p0
   integer :: itimestep, n_subcycle_steps
   type(iosystem_desc_t), pointer :: pio_subsystem
   integer :: output_frame

   type (io_output_object), save :: cam_restart_obj

   contains

   subroutine mpas_init1(mpi_comm, phys_dt, PIOFileDesc)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_INIT1
   !
   ! Initializes the MPAS software infrastructure, reads run-time configuration
   !    information, reads grid information from an MPAS grid file, and allocates
   !    storage for fields to be provided by CAM through either 
   !    cam_inidat_to_mpas() or cam_restart_to_mpas().
   !
   ! Input: mpi_comm - an MPI communicator supplied by CAM to be used by MPAS
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use shr_pio_mod, only : shr_pio_getiosys
      use mpas_stream_manager, only : MPAS_stream_mgr_init, MPAS_build_stream_filename
      use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
      use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string
      use mpas_bootstrapping, only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2
      use mpas_timekeeping

#ifdef 1
      use cam_initfiles,    only: ncdata, rest_pfile
      use cam_control_mod,  only: initial_run
      use ioFileMod,        only: opnfil
      use spmd_utils,       only: masterproc !added by KSA
#endif 
   
      implicit none

      integer :: iArg, nArgs
      logical :: readNamelistArg, readStreamsArg
      character(len=StrKIND) :: argument, namelistFile, streamsFile
      character(len=StrKIND) :: timeStamp

      character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
      integer(kind=c_int) :: c_comm
      integer(kind=c_int) :: c_ierr
      type (c_ptr) :: mgr_p
      character(len=StrKIND) :: mesh_stream
      character(len=StrKIND) :: mesh_filename
      character(len=StrKIND) :: mesh_filename_temp
      character(len=StrKIND) :: ref_time_temp
      character(len=StrKIND) :: filename_interval_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_ref_time_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_filename_interval_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_iotype

      type (MPAS_Time_type) :: start_time, ref_time, currTime
      type (MPAS_TimeInterval_type) :: filename_interval
      character(len=StrKIND) :: start_timestamp
      character(len=StrKIND) :: iotype
      logical :: streamsExists
      integer :: mesh_iotype

      integer, intent(in) :: mpi_comm
      real (kind=RKIND), pointer :: dt
      real (kind=RKIND), intent(in) :: phys_dt
      type (file_desc_t), intent(inout), optional :: PIOFileDesc

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh

      integer,pointer:: nCells
      integer :: iCell, iEdge, j, ierr
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), pointer :: east, north

      logical, pointer :: config_do_restart

      interface
         subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
            use iso_c_binding, only : c_char, c_ptr, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            type (c_ptr), intent(inout) :: mgr_p
            integer(kind=c_int), intent(inout) :: comm
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_parser

         subroutine xml_stream_get_attributes(xmlname, streamname, comm, filename, ref_time, filename_interval, io_type, ierr) bind(c)
            use iso_c_binding, only : c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            character(kind=c_char), dimension(*), intent(in) :: streamname
            integer(kind=c_int), intent(inout) :: comm
            character(kind=c_char), dimension(*), intent(out) :: filename
            character(kind=c_char), dimension(*), intent(out) :: ref_time
            character(kind=c_char), dimension(*), intent(out) :: filename_interval
            character(kind=c_char), dimension(*), intent(out) :: io_type
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_get_attributes
      end interface

      if(masterproc) then
          write(stderrUnit,*) 'Called MPAS_INIT1'
      end if

      readNamelistArg = .false.
      readStreamsArg = .false.

      nArgs = command_argument_count()
      iArg = 1
      do while (iArg < nArgs)
         call get_command_argument(iArg, argument)
         if (len_trim(argument) == 0) exit

         if ( trim(argument) == '-n' ) then
            iArg = iArg + 1
            readNamelistArg = .true.
            call get_command_argument(iArg, namelistFile)
            if ( len_trim(namelistFile) == 0 ) then
                write(0,*) 'ERROR: The -n argument requires a namelist file argument.'
                stop
            else if ( trim(namelistFile) == '-s' ) then
                write(0,*) 'ERROR: The -n argument requires a namelist file argument.'
                stop
            end if
         else if ( trim(argument) == '-s' ) then
            iArg = iArg + 1
            readStreamsArg = .true.
            call get_command_argument(iArg, streamsFile)
            if ( len_trim(streamsFile) == 0 ) then
                write(0,*) 'ERROR: The -s argument requires a streams file argument.'
                stop
            else if ( trim(streamsFile) == '-n' ) then
                write(0,*) 'ERROR: The -s argument requires a streams file argument.'
                stop
            end if
         end if

         iArg = iArg + 1
      end do

      pio_subsystem => shr_pio_getiosys('ATM')

      allocate(corelist)
      nullify(corelist % next)

      allocate(corelist % domainlist)
      nullify(corelist % domainlist % next)

      domain_ptr => corelist % domainlist
      domain_ptr % core => corelist

      call mpas_allocate_domain(domain_ptr)

      !
      ! Initialize infrastructure
      !
      call mpas_framework_init_phase1(domain_ptr % dminfo, mpi_comm)

      call atm_setup_core(corelist)
      call atm_setup_domain(domain_ptr)

      if ( readNamelistArg ) then
         domain_ptr % namelist_filename = namelistFile
      end if

      if ( readStreamsArg ) then
         domain_ptr % streams_filename = streamsFile
      end if

      ierr = domain_ptr % core % setup_namelist(domain_ptr % configs, domain_ptr % namelist_filename, domain_ptr % dminfo)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Namelist setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      !call mpas_framework_init_phase2(domain_ptr)
      call mpas_framework_init_phase2(domain_ptr, io_system=pio_subsystem)

      ierr = domain_ptr % core % define_packages(domain_ptr % packages)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Package definition failed for core ' // trim(domain_ptr % core % coreName))
      end if

      ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Package setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Decomposition setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Clock setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      write(stderrUnit,*) 'Reading streams configuration from file '//trim(domain_ptr % streams_filename)
      inquire(file=trim(domain_ptr % streams_filename), exist=streamsExists)

      if ( .not. streamsExists ) then
         call mpas_dmpar_global_abort('ERROR: Streams file '//trim(domain_ptr % streams_filename)//' does not exist.')
      end if

      !
      ! Using information from the namelist, a graph.info file, and a file containing
      !    mesh fields, build halos and allocate blocks in the domain
      !

#ifdef 1
      !SHP-CAM:: re-setup for initial/restart information 
      write (0,*) 'SHP:: initial_run? ->', initial_run
      call mpas_pool_get_config(domain_ptr % configs, 'config_do_restart', config_do_restart)
      write (0,*) 'SHP:: do_restart? ->', config_do_restart

      if (initial_run) then
         config_do_restart = .false.
         write (0,*) 'SHP:: do_restart? ', config_do_restart
      else
         config_do_restart = .true.
         write (0,*) 'SHP:: do_restart? ', config_do_restart
      end if
#endif

      ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Failed to find mesh stream for core ' // trim(domain_ptr % core % coreName))
      end if

      call mpas_f_to_c_string(domain_ptr % streams_filename, c_filename)
      call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
      c_comm = domain_ptr % dminfo % comm
      call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                     c_mesh_filename_temp, c_ref_time_temp, &
                                     c_filename_interval_temp, c_iotype, c_ierr)

      if (c_ierr /= 0) then
         call mpas_dmpar_abort(domain_ptr % dminfo)
      end if
      call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
      call mpas_c_to_f_string(c_ref_time_temp, ref_time_temp)
      call mpas_c_to_f_string(c_filename_interval_temp, filename_interval_temp)
      call mpas_c_to_f_string(c_iotype, iotype)

      if (trim(iotype) == 'pnetcdf') then
         mesh_iotype = MPAS_IO_PNETCDF
      else if (trim(iotype) == 'pnetcdf,cdf5') then
         mesh_iotype = MPAS_IO_PNETCDF5
      else if (trim(iotype) == 'netcdf') then
         mesh_iotype = MPAS_IO_NETCDF
      else if (trim(iotype) == 'netcdf4') then
         mesh_iotype = MPAS_IO_NETCDF4
      else
         mesh_iotype = MPAS_IO_PNETCDF
      end if

#ifdef 1
      !SHP-CAM:: re-setup for mesh stream (from CAM)
      if (initial_run) then
         mesh_filename_temp(1:len(ncdata)) = ncdata(:)
         write (0,*) 'SHP:: initial_run ', mesh_filename_temp
      else
         call opnfil(rest_pfile, 1000, 'f', status="old")
         read (1000, '(a)', iostat=ierr) mesh_filename_temp
         if (ierr /= 0) then
           write(stderrUnit,*) 'Error:opening file for ', rest_pfile 
         end if
         close(1000)
         write (0,*) 'SHP:: restart_run', mesh_filename_temp  
      end if
#endif

      start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
      if ( trim(ref_time_temp) == 'initial_time' ) then
          call mpas_get_time(start_time, dateTimeString=ref_time_temp, ierr=ierr)
      end if

      if ( trim(filename_interval_temp) == 'none' ) then
          call mpas_expand_string(ref_time_temp, mesh_filename_temp, mesh_filename)
      else
          call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
          call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
          call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, mesh_filename, ierr)
      end if

      itimestep = 1

      !
      ! Set physics timestep and verify that it is evenly divided by dynamics timestep
      !
      call mpas_pool_get_config(domain_ptr % configs, 'config_dt', dt)

      dt_physics = phys_dt
      dt_dynamics = dt_physics / ceiling(dt_physics / dt)
      n_subcycle_steps = ceiling(dt_physics / dt)

#ifdef 1
      !SHP-CAM:: re-setup for clock information for restart_run 
      !          MPAS restart_run start on (+ phys_dt) step.
      !          So, we need some trick to create MPAS output on regular clock time.
      !          MGD & SHP quite don't understand why CAM's time for output and restart files are different. 
      if (.not. initial_run) then
        do j=1,n_subcycle_steps
          call mpas_advance_clock(domain_ptr % clock)
        end do
        currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
      end if
#endif

      write(stderrUnit, *) ' ** Attempting to bootstrap MPAS framework using stream: ', trim(mesh_stream)
      call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

      !
      ! Set up run-time streams
      !
      call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % clock, domain_ptr % blocklist % allFields, domain_ptr % packages, domain_ptr % blocklist % allStructs)

      call add_stream_attributes(domain_ptr)

      ierr = domain_ptr % core % setup_immutable_streams(domain_ptr % streamManager)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Immutable streams setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      mgr_p = c_loc(domain_ptr % streamManager)
      call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
      if (c_ierr /= 0) then
         call mpas_dmpar_abort(domain_ptr % dminfo)
      end if

      !
      ! Finalize the setup of blocks and fields
      !
      call mpas_bootstrap_framework_phase2(domain_ptr)

      !
      ! Initialize core
      !
      iErr = domain_ptr % core % core_init(domain_ptr, timeStamp)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Core init failed for core ' // trim(domain_ptr % core % coreName))
      end if
      !call mpas_core_setup_packages(ierr)

      block => domain_ptr % blocklist
      do while (associated(block))

        call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
        call mpas_pool_get_dimension(mesh,'nCells',nCells)
  
        call mpas_pool_get_array(mesh, 'latCell', latCell)
        call mpas_pool_get_array(mesh, 'lonCell', lonCell)
        call mpas_pool_get_array(mesh, 'east'   , east   )
        call mpas_pool_get_array(mesh, 'north'  , north  )

        block => block % next
      end do

      !
      ! Compute unit vectors in east and north directions for each cell
      !
      do iCell = 1,nCells

         east(1,iCell) = -sin(lonCell(iCell))
         east(2,iCell) =  cos(lonCell(iCell))
         east(3,iCell) =  0.0
         call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

         north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
         north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
         north(3,iCell) =  cos(latCell(iCell))
         call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

      end do

      p0 = 1.0E5    ! This should be consistent with the value used in MPAS to
                    !    originally compute theta from temperature
   
   end subroutine mpas_init1

   
   subroutine mpas_global_dimensions(nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdge)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_GLOBAL_DIMENSIONS
   !
   ! Returns global grid dimensions
   !
   ! Output: nCellsGlobal - global number of horizontal columns
   !         nEdgesGlobal - global number of horizontal cell faces
   !         nVerticesGlobal - global number of horizontal cell corners
   !         maxEdge - maximum number of horizontal faces for any cell in the mesh
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh

      integer, pointer     :: maxEdges, nCellsSolve, nEdgesSolve, nVerticesSolve
      integer, intent(out) :: nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdge
      integer              :: nCellsSolveSum, nEdgesSolveSum, nVerticesSolveSum

      write(stderrUnit,*) 'Called MPAS_GLOBAL_DIMENSIONS()'

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh,'maxEdges',maxEdges)
      maxEdge = maxEdges

      nCellsSolveSum = 0
      nEdgesSolveSum = 0
      nVerticesSolveSum = 0

      block => domain_ptr % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
         call mpas_pool_get_dimension(mesh,'nEdgesSolve',nEdgesSolve)
         call mpas_pool_get_dimension(mesh,'nVerticesSolve',nVerticesSolve)

         nCellsSolveSum = nCellsSolveSum + nCellsSolve
         nEdgesSolveSum = nEdgesSolveSum + nEdgesSolve
         nVerticesSolveSum = nVerticesSolveSum + nVerticesSolve

         block => block % next
      end do

      call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolveSum, nCellsGlobal)
      call mpas_dmpar_sum_int(domain_ptr % dminfo, nEdgesSolveSum, nEdgesGlobal)
      call mpas_dmpar_sum_int(domain_ptr % dminfo, nVerticesSolveSum, nVerticesGlobal)

   end subroutine mpas_global_dimensions
   

   subroutine mpas_global_coordinates(nCellsGlobal, nVerticesGlobal, maxEdges,      &
                                      latCellGlobal, lonCellGlobal, areaCellGlobal, &
                                      latVertexGlobal, lonVertexGlobal, nEdgesOnCellGlobal, verticesOnCellGlobal)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_GLOBAL_COORDINATES
   !
   ! Returns global coordinate information for the mesh
   !
   ! Input: nCellsGlobal - global number of horizontal columns
   !        nVerticesGlobal - global number of horizontal cell corners
   !        maxEdges - maximum number of horizontal faces for any cell in the mesh
   !
   ! Output: latCell - latitude at cell centers
   !         lonCell - longitude at cell centers
   !         areaCell - area of cells (m^2)
   !         latVertex - latitude of cell corners
   !         lonVertex - longitude ov cell corners
   !         nEdgesOnCell - number of corners/faces for each cell
   !         verticesOnCell - global indices in vertex-based arrays of cell corners
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: mesh
      integer, pointer               :: nCellsSolve, nVerticesSolve 
      integer, dimension(:), pointer :: indexToCellID, indexToVertexID
      integer, dimension(:), pointer :: nEdgesOnCell, itemp, itemp2
      integer, dimension(:,:), pointer :: verticesOnCell
      real (kind=RKIND), dimension(:), pointer :: rtemp, latCell, lonCell, areaCell, latVertex, lonVertex 

      integer, intent(in) :: nCellsGlobal, nVerticesGlobal, maxEdges
      real (kind=RKIND), dimension(nCellsGlobal), intent(out) :: latCellGlobal, lonCellGlobal, areaCellGlobal
      real (kind=RKIND), dimension(nVerticesGlobal), intent(out) :: latVertexGlobal, lonVertexGlobal
      integer, dimension(nCellsGlobal), intent(out) :: nEdgesOnCellGlobal
      integer, dimension(maxEdges, nCellsGlobal), intent(out) :: verticesOnCellGlobal

      integer :: i, j

      write(stderrUnit,*) 'Called MPAS_GLOBAL_COORDINATES()'

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)

      call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
      call mpas_pool_get_dimension(mesh,'nVerticesSolve',nVerticesSolve)

      call mpas_pool_get_array(mesh,'indexToCellID',indexToCellID)
      call mpas_pool_get_array(mesh,'indexToVertexID',indexToVertexID)

      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'latVertex', latVertex)
      call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)

      !
      ! Handle 1-d real cell-based fields
      !
      allocate(rtemp(nCellsGlobal))

      rtemp(:) = -999.0
      do i=1,nCellsSolve
         rtemp(indexToCellID(i)) = latCell(i)
      end do
      call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, rtemp, latCellGlobal) 

      rtemp(:) = -999.0
      do i=1,nCellsSolve 
         rtemp(indexToCellID(i)) = lonCell(i)
      end do
      call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, rtemp, lonCellGlobal) 

      rtemp(:) = -999.0
      do i=1,nCellsSolve
         rtemp(indexToCellID(i)) = areaCell(i)
      end do
      call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, rtemp, areaCellGlobal) 

      deallocate(rtemp)


      !
      ! Handle 1-d integer cell-based fields
      !
      allocate(itemp(nCellsGlobal))

      itemp(:) = -1
      do i=1,nCellsSolve
         itemp(indexToCellID(i)) = nEdgesOnCell(i)
      end do
      call mpas_dmpar_max_int_array(domain_ptr % dminfo, nCellsGlobal, itemp, nEdgesOnCellGlobal) 

      deallocate(itemp)


      !
      ! Handle 2-d integer cell-based fields
      !
      allocate(itemp(maxEdges*nCellsGlobal))
      allocate(itemp2(maxEdges*nCellsGlobal))

      itemp(:) = -1
      do i=1,nCellsSolve
         do j=1,nEdgesOnCell(i)
            itemp(j+(indexToCellID(i)-1)*maxEdges) = indexToVertexID(verticesOnCell(j,i))
         end do
         do j=nEdgesOnCell(i)+1,maxEdges
            itemp(j+(indexToCellID(i)-1)*maxEdges) = 0
         end do
      end do
      call mpas_dmpar_max_int_array(domain_ptr % dminfo, nCellsGlobal*maxEdges, itemp, itemp2) 

      do i=1,nCellsGlobal
         do j=1,maxEdges
            verticesOnCellGlobal(j,i) = itemp2(j+(i-1)*maxEdges)
         end do
      end do

      deallocate(itemp)
      deallocate(itemp2)


      !
      ! Handle 1-d real vertex-based fields
      !
      allocate(rtemp(nVerticesGlobal))

      rtemp(:) = -999.0
      do i=1,nVerticesSolve
         rtemp(indexToVertexID(i)) = latVertex(i)
      end do
      call mpas_dmpar_max_real_array(domain_ptr % dminfo, nVerticesGlobal, rtemp, latVertexGlobal) 

      rtemp(:) = -999.0
      do i=1,nVerticesSolve
         rtemp(indexToVertexID(i)) = lonVertex(i)
      end do
      call mpas_dmpar_max_real_array(domain_ptr % dminfo, nVerticesGlobal, rtemp, lonVertexGlobal) 

      deallocate(rtemp)

   end subroutine mpas_global_coordinates


   subroutine mpas_decomp(nCellsGlobal, npes, ncells_per_proc, global_cell_owner, global_to_local_cell)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_DECOMP
   !
   ! Returns information on the parallel decomposition used by the dynamical core.
   !
   ! Input: nCellsGlobal - global (i.e., total) number of grid columns in the mesh
   !        npes - number of MPI tasks to be used by the dynamical core
   !
   ! Output: ncells_per_proc - number of columns owned by each task
   !         global_cell_owner - task ID that owns each column
   !         global_to_local_cell - maps global indices to local indices for columns
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: mesh
      integer, pointer               :: nCellsSolve
      integer, dimension(:), pointer :: indexToCellID

      integer, intent(in) :: nCellsGlobal
      integer, intent(in) :: npes
      integer, dimension(npes), intent(out) :: ncells_per_proc
      integer, dimension(nCellsGlobal), intent(out) :: global_cell_owner
      integer, dimension(nCellsGlobal), intent(out) :: global_to_local_cell

      integer :: i
      integer, dimension(:), allocatable :: temp

      write(stderrUnit,*) 'Called MPAS_DECOMP()'

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)

      call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

      call mpas_pool_get_array(mesh,'indexToCellID',indexToCellID)

      !
      !  Perform basic sanity check on expected and available field dimensions
      !
      if (npes /= domain_ptr % dminfo % nprocs ) then
         write(stderrUnit,*) 'Error: mismatch between npes and nprocs: ', npes, domain_ptr % dminfo % nprocs
         return
      end if

      allocate(temp(npes))
      temp(:) = -1
      temp(domain_ptr % dminfo % my_proc_id + 1) = nCellsSolve
      call mpas_dmpar_max_int_array(domain_ptr % dminfo, npes, temp, ncells_per_proc)
      deallocate(temp)

      allocate(temp(nCellsGlobal))
      temp(:) = -1
      do i=1,nCellsSolve
         temp(indexToCellID(i)) = domain_ptr % dminfo % my_proc_id + 1
      end do
      call mpas_dmpar_max_int_array(domain_ptr % dminfo, nCellsGlobal, temp, global_cell_owner)
      deallocate(temp)

      allocate(temp(nCellsGlobal))
      temp(:) = -1
      do i=1,nCellsSolve
         temp(indexToCellID(i)) = i
      end do
      call mpas_dmpar_max_int_array(domain_ptr % dminfo, nCellsGlobal, temp, global_to_local_cell)
      deallocate(temp)

   end subroutine mpas_decomp

   
   subroutine mpas_get_pref_profile(pref_edge, pref_mid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_GET_PREF_PROFILE
   !
   ! Returns reference profiles for pressure at layer interfaces and layer
   !    midpoints. Pressures should be in ascending order, i.e., from the top of
   !    the atmosphere downward.
   !
   ! Output: pref_edge - reference pressure at layer interfaces (vertical edges)
   !         pref_mid - reference pressure at layer midpoints
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), pointer :: mesh, diag
      integer, pointer               :: nCellsSolve, nVertLevels
      real (kind=RKIND), dimension(:,:), pointer :: pressure_b 

      real (kind=RKIND), dimension(:), intent(out) :: pref_edge
      real (kind=RKIND), dimension(:), intent(out) :: pref_mid


      integer :: k
      real (kind=RKIND), dimension(:), pointer :: temp
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag)

      call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
      call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

      call mpas_pool_get_array(mesh,'fzm',fzm)
      call mpas_pool_get_array(mesh,'fzp',fzp)

      call mpas_pool_get_array(diag, 'pressure_base', pressure_b)

      allocate(temp(nVertLevels))

      !
      ! Take layer midpoint pressures directly from pressure_base
      !
      do k=1,nVertLevels
         temp(nVertLevels-k+1) = maxval( pressure_b(k,1:nCellsSolve))
      end do

      call mpas_dmpar_max_real_array(domain_ptr % dminfo, nVertLevels, temp, pref_mid)

      deallocate(temp)

      !
      ! Interpolate interface pressure for interior edges, and extrapolate for
      !    top and bottom pressure
      !
      do k=2,nVertLevels
         pref_edge(k) = fzm(k)*pref_mid(k) + fzp(k)*pref_mid(k-1)
      end do
      pref_edge(1) = pref_mid(1) + pref_mid(1) - pref_edge(2)
      pref_edge(nVertLevels+1) = pref_mid(nVertLevels) + pref_mid(nVertLevels) - pref_edge(nVertLevels)

   end subroutine mpas_get_pref_profile

   
   subroutine mpas_init2()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_INIT2
   !
   ! Calls the core-specific (most likely the hydrostatic atmosphere core) 
   !    initialization routine mpas_core_init() after initial fields have been provided
   !    through a call to either cam_inidat_to_mpas() or cam_restart_to_mpas().
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      type (block_type), pointer :: block

      character (len=StrKIND) :: initialTimeStamp

      integer :: iErr
   
      write(stderrUnit,*) 'Called MPAS_INIT2' 

      ! Compute diagnostic fields needed in solve loop, and initialize
      !    simulation time to 0 for all blocks
!      block => domain % blocklist
!      do while (associated(block))
!         call mpas_core_init(domain, initialTimeStamp)
!         block => block % next
!      end do

!!!!! MGD IO -- maybe we don't need to do anything for output here, as we will be able to do this later?

!      output_frame = 1
!      if(config_frames_per_outfile > 0) then
!         call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(initialTimeStamp))
!      else
!         call mpas_output_state_init(output_obj, domain, "OUTPUT")
!      end if
!      call atm_write_output_frame(output_obj, output_frame, domain)

   end subroutine mpas_init2
   
   
   subroutine mpas_to_cam(Numcols, Plev, Pcnst, Psd, Phis, Pint, Pmid, Zint, Zmid, &
                          T, Ux, Uy, Omega, Tracer, Div, Vor) !KSA, Div and Vor added
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_TO_CAM
   !
   ! Input: Numcols - number of columns/cells
   !        Plev - number of vertical levels
   !        Pcnst - number of tracers
   !
   ! Output: Psd - dry surface pressure (Pa)
   !         Phis - surface geopotential (m^2/s^2)
   !         Pint - dry pressure at vertical layer interfaces
   !         Pmid - dry pressure at vertical layer mid-points
   !         Zint - geopotential height at vertical layer interfaces
   !         Zmid - geopotential height at vertical layer mid-points
   !         T - temperature (K)
   !         Ux - longitudinal velocity at cell centers (m/s)
   !         Uy - latitudinal velocity at cell centers (m/s)
   !         Omega - omega (Pa/s)
   !         Tracer - tracer mixing ratios (kg/kg)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      integer, intent(in) :: Numcols
      integer, intent(in) :: Plev
      integer, intent(in) :: Pcnst
      real (kind=RKIND), dimension(Numcols), intent(out) :: Psd
      real (kind=RKIND), dimension(Numcols), intent(out) :: Phis
      real (kind=RKIND), dimension(Numcols,Plev+1), intent(out) :: Pint
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Pmid
      real (kind=RKIND), dimension(Numcols,Plev+1), intent(out) :: Zint
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Zmid
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: T
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Ux
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Uy
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Omega
      real (kind=RKIND), dimension(Numcols,Plev,Pcnst), intent(out) :: Tracer
!++KSA
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Div
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Vor
!--KSA

      type (mpas_pool_type), pointer :: mesh, state, diag

      integer, pointer :: nCellsSolve, nVertLevels, index_qv, num_scalars
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), pointer :: theta, pressure_base, pressure_p, ww, uReconstZonal, uReconstMeridional, east, north
      real (kind=RKIND), dimension(:,:), pointer :: w, rho_p, rho_base, zz, zgrid, theta_m
!++KSA
      real (kind=RKIND), dimension(:,:), pointer :: divergence, vor_cell
!--KSA

      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      integer :: iCell, k, iScalar
      real (kind=RKIND) :: z0, z1, z2, w1, w2, rho_a

      logical :: DEBUG=.false.

      !write(stderrUnit,*) 'Called MPAS_TO_CAM'

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag)

      call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
      call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)

      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'w', w, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)

      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstMeridional)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_base', rho_base)

!++KSA
      call mpas_pool_get_array(diag, 'vor_cell', vor_cell)
      call mpas_pool_get_array(diag, 'divergence', divergence)
!--KSA
      !
      !  Perform basic sanity check on expected and available field dimensions
      !
      if (Numcols /= nCellsSolve) then
         write(stderrUnit,*) 'Error: mismatch between Numcols and nCellsSolve: ', Numcols, nCellsSolve
         return
      end if
      if (Plev /= nVertLevels) then
         write(stderrUnit,*) 'Error: mismatch between Plev and nVertLevels: ', Plev, nVertLevels
         return
      end if
      if (Pcnst /= num_scalars) then
         write(stderrUnit,*) 'Error: mismatch between Pcnst and num_scalars: ', Pcnst, num_scalars
         return
      end if

      !
      ! Fill in CAM arrays from block % state % time_levs(1) % state arrays
      !
      do iCell=1,nCellsSolve

         Phis(iCell) = zgrid(1,iCell) * gravity

         do k=1,nVertLevels
            theta(k,iCell) = theta_m(k,iCell)/(1.0 + 1.61 * scalars(index_qv,k,iCell)) 
            T(iCell,k) = theta(k,iCell) * ((pressure_base(k,iCell)+pressure_p(k,iCell)) / p0) ** (rgas/cp)

            !SHP :: Zint, Zmid are the height above the surface - RELATIVE HEIGHT
            !Zint(iCell,k) = zgrid(k,iCell)
            !Zmid(iCell,k) = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell))
            Zint(iCell,k) = zgrid(k,iCell) - zgrid(1,iCell)
            Zmid(iCell,k) = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell)) - zgrid(1,iCell)

            Omega(iCell,k) = -1.0 * gravity * (rho_p(k,iCell)+rho_base(k,iCell)) * zz(k,iCell) * (w(k,iCell)+w(k+1,iCell)) * 0.5

!++KSA
            Div(iCell,k) = divergence(k,iCell)
            Vor(iCell,k) = vor_cell(k,iCell)
!--KSA
            !
            ! NOTE: Eventually, we need to ensure that the moisture variables we return to CAM
            !       are what CAM expects, rather than simply what we have in our scalars array
            !
            do iScalar=1,num_scalars
               !SHP :: CAM-physics require wet mixing ratio for tracer #1~5 .. other than these, they require dry mixing ratio (MPAS needs dry mixing ratio for all tracers)
               !    :: In here, coupler should converts all to wet mixing ratio like FV .. 
               !    :: Another routine will be applied to change wet mixing ratio to dry mixing ratio before CAM-physics is called (for tracer #6~25) - "set_wet_to_dry".
               !    :: tracer 1 - Specific Humidity
               !    ::        2 - CLDLIQ    Grid box averaged cloud liquid amount
               !    ::        3 - CLDICE    Grid box averaged cloud ice amount  
               !    ::        4 - NUMLIQ    Grid box averaged cloud liquid number 
               !    ::        5 - NUMICE    Grid box averaged cloud ice number           
                 Tracer(iCell,k,iScalar) = scalars(iScalar,k,iCell)/(1.0 + scalars(index_qv,k,iCell))
            end do
         end do

         !SHP :: Zint, Zmid are the height above the surface - RELATIVE HEIGHT
         !Zint(iCell,Plev+1) = zgrid(Plev+1,iCell)
         Zint(iCell,Plev+1) = zgrid(Plev+1,iCell) - zgrid(1,iCell)

      end do

      k = nVertLevels + 1
      do iCell = 1, nCellsSolve
         z0 = Zint(iCell,k)
         z1 = 0.5*(Zint(iCell,k)   + Zint(iCell,k-1)) 
         z2 = 0.5*(Zint(iCell,k-1) + Zint(iCell,k-2))
         w1 = (z0-z2)/(z1-z2)
         w2 = 1.0 - w1

         !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
         Pint(iCell,k) = exp(w1*log(pressure_base(k-1,iCell)+pressure_p(k-1,iCell)) &
                            +w2*log(pressure_base(k-2,iCell)+pressure_p(k-2,iCell)))
      end do

      do iCell = 1, nCellsSolve

         !pressure at w-points:
         do k = nVertLevels,1,-1
             !SHP :: full hydrostatic pressure .. not dry pressure
             !rho_a = (rho_p(k,iCell) + rho_base(k,iCell)) * zz(k,iCell) * (1.+scalars(index_qv,k,iCell))
             !Pint(iCell,k) = Pint(iCell,k+1) + gravity*rho_a*(Zint(iCell,k+1) - Zint(iCell,k))

             !SHP :: full hydrostatic dry pressure 
             rho_a = (rho_p(k,iCell) + rho_base(k,iCell)) * zz(k,iCell)
             Pint(iCell,k) = Pint(iCell,k+1) + gravity*rho_a*(Zint(iCell,k+1) - Zint(iCell,k))
         end do                                                                                              

         !pressure at theta-points:
         do k = nVertLevels,1,-1
             Pmid(iCell,k) = 0.5*(Pint(iCell,k+1) + Pint(iCell,k))
             !Pressure(iCell,k) = Pmid(iCell,k)
             !Pressure(iCell,k) = pressure_base(k,iCell)+pressure_p(k,iCell)
         end do

         !surface pressure:
         Psd(iCell) = Pint(iCell,1)
      end do

      do iCell=1, nCellsSolve
         do k=1, nVertLevels
            Ux(iCell,k) =   uReconstZonal(k,iCell)
            Uy(iCell,k) =   uReconstMeridional(k,iCell)
         end do
      end do 

      if (DEBUG) then
        write (0,*) '********* CHECK VARIABLES TO CAM *********'
        write (0,*) 'zonal wind =', minval(Ux), maxval(Ux)
        write (0,*) 'meridional wind =', minval(Uy), maxval(Uy)
        write (0,*) 'temperature =', minval(T), maxval(T)
        write (0,*) 'Psd =', minval(Psd), maxval(Psd)
        write (0,*) 'Pmid =', minval(Pmid), maxval(Pmid)
        write (0,*) 'Pint =', minval(Pint), maxval(Pint)
        write (0,*) 'Zint =', minval(Zint), maxval(Zint)
        write (0,*) 'num_scalars =', num_scalars
        do iScalar=1,num_scalars
           write (0,*) 'scalars ->',iScalar,' :: ',minval(scalars(iScalar,:,:)),maxval(scalars(iScalar,:,:))
        end do
      end if

   end subroutine mpas_to_cam
   
   
   subroutine cam_to_mpas(Numcols, Plev, Pcnst, T_tend, Ux_tend, Uy_tend, Tracer)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CAM_TO_MPAS
   !
   ! Input: Numcols - number of columns/cells
   !        Plev - number of vertical levels
   !        Pcnst - number of tracers
   !        T_tend - temperature tendency (K/s)
   !        Ux_tend - longitudinal velocity tendency at cell centers (m/s^2)
   !        Uy_tend - latitudinal velocity tendency at cell centers (m/s^2)
   !        Tracer - tracer mixing ratios (kg/kg)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      integer, intent(in) :: Numcols
      integer, intent(in) :: Plev
      integer, intent(in) :: Pcnst
      real (kind=RKIND), dimension(Numcols,Plev), intent(in) :: T_tend
      real (kind=RKIND), dimension(Numcols,Plev), intent(in) :: Ux_tend
      real (kind=RKIND), dimension(Numcols,Plev), intent(in) :: Uy_tend
      real (kind=RKIND), dimension(Numcols,Plev,Pcnst), intent(in) :: Tracer

      type (mpas_pool_type), pointer :: mesh, state, diag, tend_physics

      integer, pointer :: nCells, nCellsSolve, nVertLevels, num_scalars
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: pressure_base, pressure_p
      real (kind=RKIND), dimension(:,:), pointer :: east, north, edge_normal, pressure, w, theta_m, zz, zgrid
      real (kind=RKIND), dimension(:,:), pointer :: theta_tend, u_tend
      real (kind=RKIND), dimension(:,:), pointer :: zonal_tend, meridional_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, scalars_tend
      type (field2DReal), pointer :: zonal_tend_field
      type (field2DReal), pointer :: meridional_tend_field
      type (field2DReal), pointer :: u_tend_field
      type (field2DReal), pointer :: theta_tend_field
      type (field3DReal), pointer :: scalars_tend_field
      integer, pointer :: index_qv

      integer :: iCell, iEdge, iScalar, k, j

      !write(stderrUnit,*) 'Called CAM_TO_MPAS'

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_physics)

      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
      call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'east'   ,east   )
      call mpas_pool_get_array(mesh, 'north'  ,north  )
      call mpas_pool_get_array(mesh, 'edgeNormalVectors', edge_normal)

      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'w', w, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)

      call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)

      call mpas_pool_get_array(tend_physics, 'scalars', scalars_tend)
      call mpas_pool_get_array(tend_physics, 'theta', theta_tend)
      call mpas_pool_get_array(tend_physics, 'u', u_tend)
      call mpas_pool_get_array(tend_physics, 'ux', zonal_tend)
      call mpas_pool_get_array(tend_physics, 'uy', meridional_tend)

      !
      !  Perform basic sanity check on expected and available field dimensions
      !
      if (Numcols /= nCellsSolve) then
         write(stderrUnit,*) 'Error: mismatch between Numcols and nCellsSolve: ', Numcols, nCellsSolve
         return
      end if
      if (Plev /= nVertLevels) then
         write(stderrUnit,*) 'Error: mismatch between Plev and nVertLevels: ', Plev, nVertLevels
         return
      end if
      if (Pcnst /= num_scalars) then
         write(stderrUnit,*) 'Error: mismatch between Pcnst and num_scalars: ', Pcnst, num_scalars
         return
      end if

      !
      ! Fill in MPAS tendency arrays from arguments
      !
      do iCell=1,nCellsSolve
         do k=1,nVertLevels
            zonal_tend(k,iCell) = Ux_tend(iCell,k)
            meridional_tend(k,iCell) = Uy_tend(iCell,k)
         end do
      end do

      call mpas_pool_get_field(tend_physics, 'ux', zonal_tend_field)
      call mpas_pool_get_field(tend_physics, 'uy', meridional_tend_field)

      call mpas_dmpar_exch_halo_field(zonal_tend_field)
      call mpas_dmpar_exch_halo_field(meridional_tend_field)

      u_tend(:,:) = 0.0
      do iCell=1,nCells
         do j=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(j,iCell)
            do k=1,nVertLevels
               u_tend(k,iEdge) = u_tend(k,iEdge) + 0.5 * zonal_tend(k,iCell) * (edge_normal(1,iEdge) * east(1,iCell) + &
                                                                                edge_normal(2,iEdge) * east(2,iCell) + &
                                                                                edge_normal(3,iEdge) * east(3,iCell))  &
                                                 + 0.5 * meridional_tend(k,iCell) * (edge_normal(1,iEdge) * north(1,iCell) + &
                                                                                     edge_normal(2,iEdge) * north(2,iCell) + &
                                                                                     edge_normal(3,iEdge) * north(3,iCell))
            end do
         end do

      end do

      do iCell=1,nCellsSolve
         do k=1,nVertLevels

            theta_tend(k,iCell) = T_tend(iCell,k) * (p0 / (pressure_base(k,iCell)+pressure_p(k,iCell))) ** (rgas/cp)
            
            !
            ! NOTE: Once we begin to use more than just qv in MPAS, we need to make sure to
            !       properly assign variables from Tracer array provided from CAM physics
            !
            do iScalar=1,num_scalars
               !SHP :: CAM-physics require wet mixing ratio for tracer #1~5 .. other than these, they require dry mixing ratio (MPAS needs dry mixing ratio for all tracers)
               !    :: In here, every Tracer is wet mixing ratio like FV .. thus, coupler should convert all to dry mixing ratio .. 
               !    :: Another routine have been applied to change dry mixing ratio to wet mixing ratio after CAM-physics is called (for tracer #6~25) - "set_dry_to_wet".
               !    :: tracer 1 - Specific Humidity
               !    ::        2 - CLDLIQ    Grid box averaged cloud liquid amount
               !    ::        3 - CLDICE    Grid box averaged cloud ice amount  
               !    ::        4 - NUMLIQ    Grid box averaged cloud liquid number 
               !    ::        5 - NUMICE    Grid box averaged cloud ice number           
                 scalars_tend(iScalar,k,iCell) = ( Tracer(iCell,k,iScalar)/(1.0 - Tracer(iCell,k,index_qv)) - scalars(iScalar,k,iCell) ) / dt_physics
            end do

         end do
      end do


      !
      ! Do a ghost-cell update on tendency arrays -- SHP:: is this right place?
      !
      call mpas_pool_get_field(tend_physics, 'scalars', scalars_tend_field)
      call mpas_pool_get_field(tend_physics, 'theta', theta_tend_field)
      call mpas_pool_get_field(tend_physics, 'u', u_tend_field)

      call mpas_dmpar_exch_halo_field(theta_tend_field)
      call mpas_dmpar_exch_halo_field(u_tend_field)
      call mpas_dmpar_exch_halo_field(scalars_tend_field)

   end subroutine cam_to_mpas
   
   
   subroutine mpas_dyn_run()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_DYN_RUN
   !
   ! Advance MPAS dynamics by one timestep.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_timekeeping
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_derived_types, only : MPAS_STREAM_LATEST_BEFORE, MPAS_STREAM_INPUT, MPAS_STREAM_INPUT_OUTPUT
      use mpas_timer

      implicit none

      !type (domain_type), intent(inout) :: domain
      integer :: ierr

      real (kind=RKIND) :: dt
      logical, pointer :: config_do_restart
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      integer :: itimestep

      integer :: stream_dir, idynstep
      character(len=StrKIND) :: input_stream, read_time

      type (mpas_pool_type), pointer :: state, diag, diag_physics, mesh

      ! For high-frequency diagnostics output

      ! Eventually, dt should be domain specific
      !call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)
      call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Avoid writing a restart file at the initial time
      call MPAS_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      ! Also, for restart runs, avoid writing the initial history or diagnostics fields to avoid overwriting those from the preceding run
      if (config_do_restart) then
         call MPAS_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID='output', direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         call MPAS_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID='diagnostics', direction=MPAS_STREAM_OUTPUT, ierr=ierr)
      end if

      if (MPAS_stream_mgr_ringing_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
         block_ptr => domain_ptr % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
            call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
            call atm_compute_output_diagnostics(state, 1, diag, mesh)

            block_ptr => block_ptr % next
         end do
      end if
      call mpas_stream_mgr_write(domain_ptr % streamManager, ierr=ierr)
      if (ierr /= MPAS_STREAM_MGR_NOERR .and. &
          ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_FILE .and. &
          ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_REC) then
         write(0,*) ' '
         write(0,*) '********************************************************************************'
         write(0,*) 'Error writing one or more output streams'
         call mpas_dmpar_global_abort('********************************************************************************')
      end if
      call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      !itimestep = 1
      !do while (.not. mpas_is_clock_stop_time(clock))

      dt = dt_dynamics
      write(stderrUnit,*) 'MPAS running dynamics for n_subcycle_steps=',n_subcycle_steps
      do idynstep=1,n_subcycle_steps

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

         write(0,*) ' '
         write(0,*) 'Begin timestep ', trim(timeStamp)

         !
         ! Read external field updates
         !
         call MPAS_stream_mgr_begin_iteration(domain_ptr % streamManager, ierr=ierr)
         do while (MPAS_stream_mgr_get_next_stream(domain_ptr % streamManager, streamID=input_stream, directionProperty=stream_dir))
            if (stream_dir == MPAS_STREAM_INPUT .or. stream_dir == MPAS_STREAM_INPUT_OUTPUT) then
               if (MPAS_stream_mgr_ringing_alarms(domain_ptr % streamManager, streamID=input_stream, &
                                                  direction=MPAS_STREAM_INPUT, ierr=ierr)) then
                  call MPAS_stream_mgr_read(domain_ptr % streamManager, streamID=input_stream, whence=MPAS_STREAM_LATEST_BEFORE, &
                                            actualWhen=read_time, ierr=ierr)
                  if (ierr /= MPAS_STREAM_MGR_NOERR) then
                     write(0,*) ' '
                     write(0,*) '********************************************************************************'
                     write(0,*) 'Error reading input stream '//trim(input_stream)
                     call mpas_dmpar_global_abort('********************************************************************************')
                  end if

                  write(0,*) '----------------------------------------------------------------------'
                  write(0,*) '  Read '''//trim(input_stream)//''' input stream valid at '//trim(read_time)
                  write(0,*) '----------------------------------------------------------------------'

                  call MPAS_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID=input_stream, direction=MPAS_STREAM_INPUT, ierr=ierr)
               end if
            end if
         end do

         call mpas_timer_start("time integration")
         call atm_do_timestep(domain_ptr, dt, itimestep)
         call mpas_timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
         call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
         call mpas_pool_shift_time_levels(state)

         ! Advance clock before writing output
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)

         !
         ! Write any output streams that have alarms ringing, after computing diagnostics fields
         !
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         if (MPAS_stream_mgr_ringing_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
            block_ptr => domain_ptr % blocklist
            do while (associated(block_ptr))

               call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
               call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
               call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
               call atm_compute_output_diagnostics(state, 1, diag, mesh)

               block_ptr => block_ptr % next
            end do
         end if

!         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
!            block_ptr => domain % blocklist
!            do while (associated(block_ptr))
!
!               call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
!               call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
!               call mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
!               call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
!               call atm_compute_restart_diagnostics(state, 1, diag, mesh)
!
!               block_ptr => block_ptr % next
!            end do
!         end if

         call mpas_stream_mgr_write(domain_ptr % streamManager, ierr=ierr)

!         if (ierr /= MPAS_STREAM_MGR_NOERR .and. &
!             ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_FILE .and. &
!             ierr /= MPAS_STREAM_MGR_ERR_CLOBBER_REC) then
!            write(0,*) ' '
!            write(0,*) '********************************************************************************'
!            write(0,*) 'Error writing one or more output streams'
!            call mpas_dmpar_global_abort('********************************************************************************')
!         end if
!
!         ! Only after we've successfully written the restart file should we we
!         !    write the restart_timestamp file
!         if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then
!            if (domain % dminfo % my_proc_id == 0) then
!               open(22,file=trim(config_restart_timestamp_name),form='formatted',status='replace')
!               write(22,*) trim(timeStamp)
!               close(22)
!            end if
!         end if

         call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      end do

   end subroutine mpas_dyn_run


   subroutine mpas_init_restart(PIOFileDesc)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_INIT_RESTART
   !
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_stream_manager
      use mpas_io_restart
   
      implicit none

      type (mpas_pool_type), pointer :: mesh
      type (file_desc_t), intent(inout) :: PIOFileDesc
      type (mpas_pool_field_info_type) :: fieldInfo

      character(len=StrKIND) :: fieldName
      logical :: fieldActive
      integer :: ierr

      real(kind=RKIND), dimension(:,:,:), pointer :: real3d
      real(kind=RKIND), dimension(:,:), pointer :: real2d
      real(kind=RKIND), dimension(:), pointer :: real1d
      real(kind=RKIND), pointer :: real0d

      integer, dimension(:,:), pointer :: int2d
      integer, dimension(:), pointer :: int1d
      integer, pointer :: int0d

      character(len=StrKIND), dimension(:), pointer :: char1d
      character(len=StrKIND), pointer :: char0d

      integer, parameter :: RESTART = 2

      write(stderrUnit,*) 'Called MPAS_INIT_RESTART'

      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)

      cam_restart_obj % stream = RESTART
      call mpas_io_output_init_MPAS2(domain_ptr, cam_restart_obj, domain_ptr % dminfo, mesh, PIOFileDesc)

   end subroutine mpas_init_restart
   
   
   subroutine mpas_write_restart(PIOFileDesc)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_WRITE_RESTART
   !
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use mpas_stream_manager
      use mpas_io_restart
   
      implicit none

      type (file_desc_t), intent(inout) :: PIOFileDesc
      integer :: ierr
   
      write(stderrUnit,*) 'Called MPAS_WRITE_RESTART'

      call mpas_output_state_for_domain_MPAS2(cam_restart_obj, domain_ptr, 1)

    !  call prewrite_reindex(domain_ptr % blocklist % allFields, stream % field_pool)
    !  call MPAS_writeStream_MPAS2(cam_restart_obj % io_stream, 1, ierr=ierr)

   end subroutine mpas_write_restart


   subroutine mpas_final()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_FINAL
   !
   ! Deallocates storage allocated by MPAS and shuts down the MPAS software
   !    infrastructure.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_stream_manager, only : MPAS_stream_mgr_finalize

      implicit none

      integer :: iErr
   
      write(stderrUnit,*) 'Called MPAS_FINAL'

      iErr = domain_ptr % core % core_finalize(domain_ptr)
      if ( iErr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Core finalize failed for core ' // trim(domain_ptr % core % coreName))
      end if

      !
      ! Finalize infrastructure
      !
      call MPAS_stream_mgr_finalize(domain_ptr % streamManager)

      call mpas_framework_finalize(domain_ptr % dminfo, domain_ptr, io_system=pio_subsystem)

      deallocate(corelist % domainlist)
      deallocate(corelist)
   
   end subroutine mpas_final


   subroutine add_stream_attributes(domain)

      use mpas_stream_manager, only : MPAS_stream_mgr_add_att

      implicit none

      type (domain_type), intent(inout) :: domain

      type (MPAS_Pool_iterator_type) :: itr
      integer, pointer :: intAtt
      logical, pointer :: logAtt
      character (len=StrKIND), pointer :: charAtt
      real (kind=RKIND), pointer :: realAtt
      character (len=StrKIND) :: histAtt

      integer :: local_ierr

      if (domain % dminfo % nProcs < 10) then
          write(histAtt, '(A,I1,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100) then
          write(histAtt, '(A,I2,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 1000) then
          write(histAtt, '(A,I3,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 10000) then
          write(histAtt, '(A,I4,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100000) then
          write(histAtt, '(A,I5,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else
          write(histAtt, '(A,I6,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      end if

      call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'source', domain % core % source)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', domain % core % Conventions)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', domain % core % git_version)

      call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
!     call MPAS_stream_mgr_add_att(domain % streamManager, 'history', domain % history)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', histAtt)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain %  parent_id)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)

      call mpas_pool_begin_iteration(domain % configs)

      do while (mpas_pool_get_next_member(domain % configs, itr))

         if ( itr % memberType == MPAS_POOL_CONFIG) then

            if ( itr % dataType == MPAS_POOL_REAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, realAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_INTEGER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, intAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, charAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
               if (logAtt) then
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'YES', ierr=local_ierr)
               else
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'NO', ierr=local_ierr)
               end if
            end if

          end if
      end do

   end subroutine add_stream_attributes


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ROTATE_ABOUT_VECTOR
   !
   ! Rotates the point (x,y,z) through an angle theta about the vector
   !   originating at (a, b, c) and having direction (u, v, w).
   !
   ! Reference: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine rotate_about_vector(x, y, z, theta, a, b, c, u, v, w, xp, yp, zp)
   
      implicit none
   
      real (kind=RKIND), intent(in) :: x, y, z, theta, a, b, c, u, v, w
      real (kind=RKIND), intent(out) :: xp, yp, zp
   
      real (kind=RKIND) :: vw2, uw2, uv2
      real (kind=RKIND) :: m
   
      vw2 = v**2.0 + w**2.0
      uw2 = u**2.0 + w**2.0
      uv2 = u**2.0 + v**2.0
      m = sqrt(u**2.0 + v**2.0 + w**2.0)
   
      xp = (a*vw2 + u*(-b*v-c*w+u*x+v*y+w*z) + ((x-a)*vw2+u*(b*v+c*w-v*y-w*z))*cos(theta) + m*(-c*v+b*w-w*y+v*z)*sin(theta))/m**2.0
      yp = (b*uw2 + v*(-a*u-c*w+u*x+v*y+w*z) + ((y-b)*uw2+v*(a*u+c*w-u*x-w*z))*cos(theta) + m*( c*u-a*w+w*x-u*z)*sin(theta))/m**2.0
      zp = (c*uv2 + w*(-a*u-b*v+u*x+v*y+w*z) + ((z-c)*uv2+w*(a*u+b*v-u*x-v*y))*cos(theta) + m*(-b*u+a*v-v*x+u*y)*sin(theta))/m**2.0
   
   end subroutine rotate_about_vector


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE R3_CROSS
   !
   ! Computes the cross product of (ax, ay, az) and (bx, by, bz).
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine r3_cross(ax, ay, az, bx, by, bz, cx, cy, cz)

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az
      real (kind=RKIND), intent(in) :: bx, by, bz
      real (kind=RKIND), intent(out) :: cx, cy, cz

      cx = ay * bz - az * by
      cy = az * bx - ax * bz
      cz = ax * by - ay * bx

   end subroutine r3_cross 


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE R3_NORMALIZE
   !
   ! Normalizes the vector (ax, ay, az)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine r3_normalize(ax, ay, az)

      implicit none

      real (kind=RKIND), intent(inout) :: ax, ay, az

      real (kind=RKIND) :: mi

      mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
      ax = ax * mi
      ay = ay * mi
      az = az * mi

   end subroutine r3_normalize 


end module mpas_cam_interface
