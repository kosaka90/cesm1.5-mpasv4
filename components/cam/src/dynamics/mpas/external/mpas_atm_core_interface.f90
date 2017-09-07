! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module atm_core_interface


   contains


   !***********************************************************************
   !
   !  routine atm_setup_core
   !
   !> \brief   Atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use atm_core, only : atm_core_init, atm_core_run, atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => atm_core_init
      core % core_run => atm_core_run
      core % core_finalize => atm_core_finalize
      core % define_packages => atm_define_packages
      core % setup_packages => atm_setup_packages
      core % setup_decompositions => atm_setup_decompositions
      core % setup_clock => atm_setup_clock
      core % get_mesh_stream => atm_get_mesh_stream
      core % setup_immutable_streams => atm_setup_immutable_streams
      core % setup_derived_dimensions => atm_setup_derived_dimensions
      core % setup_block => atm_setup_block
      core % setup_namelist => atm_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

       core % modelName = 'mpas'
       core % coreName = 'atmosphere'
       core % modelVersion = '4.0'
       core % executableName = 'atmosphere_model'
       core % git_version = 'v4.0-4-gcba1d7e-dirty'

   end subroutine atm_setup_core


   !***********************************************************************
   !
   !  routine atm_setup_domain
   !
   !> \brief   Atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

       domain % namelist_filename = 'namelist.atmosphere'
       domain % streams_filename = 'streams.atmosphere'

   end subroutine atm_setup_domain


   !***********************************************************************
   !
   !  function atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function atm_setup_packages(configs, packages) result(ierr)

      use mpas_derived_types, only : mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (mpas_pool_type), intent(inout) :: packages
      integer :: ierr

      ierr = 0

   end function atm_setup_packages


   !***********************************************************************
   !
   !  function atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_derived_types, only : MPAS_Clock_type, mpas_pool_type
      use atm_core, only : atm_simulation_clock_init

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      ierr = 0

      call atm_simulation_clock_init(core_clock, configs, ierr)

   end function atm_setup_clock


   !***********************************************************************
   !
   !  function atm_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function atm_get_mesh_stream(configs, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart

      ierr = 0

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         call mpas_dmpar_global_abort('ERROR: config_do_restart was not found when defining mesh stream.')
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end function atm_get_mesh_stream


   !***********************************************************************
   !
   !  function atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function atm_setup_decompositions(decompList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function atm_setup_decompositions


   !***********************************************************************
   !
   !  function atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be
   !>  an atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ierr = 0

      call atm_generate_structs(block, block % structs, block % dimensions, block % packages)

   end function atm_setup_block


function atm_setup_immutable_streams(manager) result(iErr)

   use MPAS_derived_types, only : MPAS_streamManager_type, &
                                  MPAS_STREAM_INPUT_OUTPUT, MPAS_STREAM_INPUT, &
                                  MPAS_STREAM_OUTPUT, MPAS_STREAM_NONE, MPAS_STREAM_PROPERTY_IMMUTABLE
   use MPAS_stream_manager, only : MPAS_stream_mgr_create_stream, MPAS_stream_mgr_set_property, &
                                   MPAS_stream_mgr_add_field, MPAS_stream_mgr_add_pool
   use mpas_io_units

   implicit none

   type (MPAS_streamManager_type), pointer :: manager
   character (len=StrKIND) :: packages
   integer :: iErr

   iErr = 0

   call MPAS_stream_mgr_create_stream(manager, 'input', MPAS_STREAM_INPUT, 'x1.40962.init.nc', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'scalars', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'latCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'lonCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'yCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'indexToCellID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'latEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'lonEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'yEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'indexToEdgeID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'latVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'lonVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'yVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'indexToVertexID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'nEdgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'nEdgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'weightsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'dvEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'dcEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'angleEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'areaCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'areaTriangle', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgeNormalVectors', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'localVerticalUnitVectors', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellTangentPlane', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellsOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'verticesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'verticesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'edgesOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cellsOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'kiteAreasOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'fEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'fVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'meshDensity', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cf1', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cf2', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'cf3', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zgrid', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'rdzw', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'dzu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'rdzu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'fzm', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'fzp', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zx', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zz', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zb', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'zb3', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'dss', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'u_init', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 't_init', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'qv_init', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'deriv_two', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'defc_a', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'defc_b', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'coeffs_reconstruct', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'xtime', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'u', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'w', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'rho', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'theta', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'relhum', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'exner_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'pressure_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'rho_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'theta_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'input', 'surface_pressure', ierr=ierr)
   call MPAS_stream_mgr_set_property(manager, 'input', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)

   call MPAS_stream_mgr_create_stream(manager, 'restart', MPAS_STREAM_INPUT_OUTPUT, 'restart.$Y-$M-$D_$h.$m.$s.nc', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'scalars', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'latCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'lonCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'xCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'yCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'indexToCellID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'latEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'lonEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'xEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'yEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'indexToEdgeID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'latVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'lonVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'xVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'yVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'indexToVertexID', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cellsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'nEdgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'nEdgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'edgesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'edgesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'weightsOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'dvEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'dcEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'angleEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'areaCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'areaTriangle', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'edgeNormalVectors', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'localVerticalUnitVectors', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cellTangentPlane', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cellsOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'verticesOnCell', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'verticesOnEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'edgesOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cellsOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'kiteAreasOnVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'fEdge', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'fVertex', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'meshDensity', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'meshScalingDel2', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'meshScalingDel4', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cf1', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cf2', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cf3', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cpr', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'cpl', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zgrid', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rdzw', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'dzu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rdzu', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'fzm', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'fzp', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zx', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zz', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zb', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'zb3', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'pzm', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'pzp', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'dss', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'u_init', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 't_init', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'qv_init', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'deriv_two', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'defc_a', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'defc_b', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'coeffs_reconstruct', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'east', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'north', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'xtime', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'u', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'w', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rho_zz', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'theta_m', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'pressure_p', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rho', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'theta', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'relhum', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'uReconstructZonal', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'uReconstructMeridional', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rv', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'circulation', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'exner', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'exner_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rtheta_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'pressure_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rho_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'theta_base', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'ru', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'ru_p', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rw', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rw_p', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rtheta_p', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rho_p', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'surface_pressure', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'tend_sfc_pressure', ierr=ierr)
   call MPAS_stream_mgr_add_field(manager, 'restart', 'rt_diabatic_tend', ierr=ierr)
   call MPAS_stream_mgr_set_property(manager, 'restart', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)

end function atm_setup_immutable_streams



   function atm_setup_derived_dimensions(readDimensions, dimensionPool, configPool) result(iErr)

      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units

      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: configPool !< Input: Pool containing namelist options with configs
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into

      integer :: iErr, errLevel

      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, pointer :: maxEdges
      integer, pointer :: maxEdges2
      integer, pointer :: nVertices
      integer, pointer :: TWO
      integer, pointer :: THREE
      integer, pointer :: vertexDegree
      integer, pointer :: FIFTEEN
      integer, pointer :: TWENTYONE
      integer, pointer :: R3
      integer, pointer :: nVertLevels
      integer, pointer :: nVertLevelsP1

      iErr = 0
      errLevel = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)


      nullify(nCells)
      call mpas_pool_get_dimension(dimensionPool, 'nCells', nCells)
      nullify(nEdges)
      call mpas_pool_get_dimension(dimensionPool, 'nEdges', nEdges)
      nullify(maxEdges)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges', maxEdges)
      nullify(maxEdges2)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      nullify(nVertices)
      call mpas_pool_get_dimension(dimensionPool, 'nVertices', nVertices)
      nullify(TWO)
      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      nullify(THREE)
      call mpas_pool_get_dimension(dimensionPool, 'THREE', THREE)
      nullify(vertexDegree)
      call mpas_pool_get_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      nullify(FIFTEEN)
      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      nullify(TWENTYONE)
      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      nullify(R3)
      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      nullify(nVertLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      nullify(nVertLevelsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)

write(stderrUnit,'(a)') 'Assigning remaining dimensions from definitions in Registry.xml ...'
      if ( .not. associated(nCells) ) then
         allocate(nCells)
         nCells = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nCells', nCells)
      end if

      if ( .not. associated(nEdges) ) then
         allocate(nEdges)
         nEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nEdges', nEdges)
      end if

      if ( .not. associated(maxEdges) ) then
         allocate(maxEdges)
         maxEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges', maxEdges)
      end if

      if ( .not. associated(maxEdges2) ) then
         allocate(maxEdges2)
         maxEdges2 = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      end if

      if ( .not. associated(nVertices) ) then
         allocate(nVertices)
         nVertices = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nVertices', nVertices)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      if ( .not. associated(TWO) ) then
         allocate(TWO)
         TWO = 2
write(stderrUnit,'(a,a20,a,i8)') '       ', 'TWO', ' =', 2
         call mpas_pool_add_dimension(dimensionPool, 'TWO', TWO)
          else if ( TWO == MPAS_MISSING_DIM ) then
         TWO = 2
          end if

      call mpas_pool_get_dimension(dimensionPool, 'THREE', THREE)
      if ( .not. associated(THREE) ) then
         allocate(THREE)
         THREE = 3
write(stderrUnit,'(a,a20,a,i8)') '       ', 'THREE', ' =', 3
         call mpas_pool_add_dimension(dimensionPool, 'THREE', THREE)
          else if ( THREE == MPAS_MISSING_DIM ) then
         THREE = 3
          end if

      if ( .not. associated(vertexDegree) ) then
         allocate(vertexDegree)
         vertexDegree = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      if ( .not. associated(FIFTEEN) ) then
         allocate(FIFTEEN)
         FIFTEEN = 15
write(stderrUnit,'(a,a20,a,i8)') '       ', 'FIFTEEN', ' =', 15
         call mpas_pool_add_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
          else if ( FIFTEEN == MPAS_MISSING_DIM ) then
         FIFTEEN = 15
          end if

      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      if ( .not. associated(TWENTYONE) ) then
         allocate(TWENTYONE)
         TWENTYONE = 21
write(stderrUnit,'(a,a20,a,i8)') '       ', 'TWENTYONE', ' =', 21
         call mpas_pool_add_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
          else if ( TWENTYONE == MPAS_MISSING_DIM ) then
         TWENTYONE = 21
          end if

      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      if ( .not. associated(R3) ) then
         allocate(R3)
         R3 = 3
write(stderrUnit,'(a,a20,a,i8)') '       ', 'R3', ' =', 3
         call mpas_pool_add_dimension(dimensionPool, 'R3', R3)
          else if ( R3 == MPAS_MISSING_DIM ) then
         R3 = 3
          end if

      if ( .not. associated(nVertLevels) ) then
         allocate(nVertLevels)
         nVertLevels = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      if ( .not. associated(nVertLevelsP1) ) then
         allocate(nVertLevelsP1)
         nVertLevelsP1 = nVertLevels+1
write(stderrUnit,'(a,a20,a,i8)') '       ', 'nVertLevelsP1', ' =', nVertLevels+1
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
          else if ( nVertLevelsP1 == MPAS_MISSING_DIM ) then
         nVertLevelsP1 = nVertLevels+1
          end if

write(stderrUnit,*) ' '
write(stderrUnit,'(a)') ' ----- done assigning dimensions from Registry.xml -----'
write(stderrUnit,*) ' '
write(stderrUnit,*) ' '
      call mpas_pool_set_error_level(errLevel)

   end function atm_setup_derived_dimensions


   subroutine mpas_setupatmdecomposed_dimensions(block, manager, readDimensions, dimensionPool, totalBlocks)

      use mpas_derived_types
      use mpas_decomp
      use mpas_pool_routines
      use mpas_io_units

      type (block_type), intent(inout) :: block !< Input: Pointer to block
      type (mpas_streamManager_type), intent(inout) :: manager !< Input: Stream manager
      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into
      integer, intent(in) :: totalBlocks !< Input: Number of blocks

      integer :: iErr
      type (field1DInteger), pointer :: ownedIndices
      procedure (mpas_decomp_function), pointer :: decompFunc


write(stderrUnit,'(a)') 'Processing decomposed dimensions ...'

write(stderrUnit,*) ' '
write(stderrUnit,'(a)') ' ----- done processing decomposed dimensions -----'
write(stderrUnit,*) ' '
write(stderrUnit,*) ' '

   end subroutine mpas_setupatmdecomposed_dimensions

   function atm_define_packages(packagePool) result(iErr)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: packagePool !< Input: MPAS Pool for containing package logicals.

      integer :: iErr

      iErr = 0
      call mpas_pool_add_package(packagePool, 'cam5Active', .true.)
   end function atm_define_packages

   subroutine atm_generate_pool_mesh(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: cam5Active


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'cam5Active', cam5Active)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'mesh', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'mesh', newSubPool)

! Define variable latCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'latCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latCell', r1Ptr)

! Define variable lonCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lonCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonCell', r1Ptr)

! Define variable xCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xCell', r1Ptr)

! Define variable yCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'yCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yCell', r1Ptr)

! Define variable zCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zCell', r1Ptr)

! Define variable indexToCellID
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'indexToCellID'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
     i1Ptr(1) % defaultValue = 0
     nullify(i1Ptr(1) % array)
      nullify(i1Ptr(1) % next)
      nullify(i1Ptr(1) % prev)
      nullify(i1Ptr(1) % sendList)
      nullify(i1Ptr(1) % recvList)
      nullify(i1Ptr(1) % copyList)
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToCellID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToCellID', i1Ptr)

! Define variable latEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'latEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latEdge', r1Ptr)

! Define variable lonEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lonEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonEdge', r1Ptr)

! Define variable xEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xEdge', r1Ptr)

! Define variable yEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'yEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yEdge', r1Ptr)

! Define variable zEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zEdge', r1Ptr)

! Define variable indexToEdgeID
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'indexToEdgeID'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
     i1Ptr(1) % defaultValue = 0
     nullify(i1Ptr(1) % array)
      nullify(i1Ptr(1) % next)
      nullify(i1Ptr(1) % prev)
      nullify(i1Ptr(1) % sendList)
      nullify(i1Ptr(1) % recvList)
      nullify(i1Ptr(1) % copyList)
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToEdgeID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToEdgeID', i1Ptr)

! Define variable latVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'latVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latVertex', r1Ptr)

! Define variable lonVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lonVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonVertex', r1Ptr)

! Define variable xVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xVertex', r1Ptr)

! Define variable yVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'yVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yVertex', r1Ptr)

! Define variable zVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zVertex', r1Ptr)

! Define variable indexToVertexID
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'indexToVertexID'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nVertices'
     i1Ptr(1) % defaultValue = 0
     nullify(i1Ptr(1) % array)
      nullify(i1Ptr(1) % next)
      nullify(i1Ptr(1) % prev)
      nullify(i1Ptr(1) % sendList)
      nullify(i1Ptr(1) % recvList)
      nullify(i1Ptr(1) % copyList)
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToVertexID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToVertexID', i1Ptr)

! Define variable cellsOnEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'cellsOnEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'TWO'
      i2Ptr(1) % dimNames(2) = 'nEdges'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnEdge', i2Ptr)

! Define variable nEdgesOnCell
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nEdgesOnCell'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
     i1Ptr(1) % defaultValue = 0
     nullify(i1Ptr(1) % array)
      nullify(i1Ptr(1) % next)
      nullify(i1Ptr(1) % prev)
      nullify(i1Ptr(1) % sendList)
      nullify(i1Ptr(1) % recvList)
      nullify(i1Ptr(1) % copyList)
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nEdgesOnCell', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nEdgesOnCell', i1Ptr)

! Define variable nEdgesOnEdge
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nEdgesOnEdge'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
     i1Ptr(1) % defaultValue = 0
     nullify(i1Ptr(1) % array)
      nullify(i1Ptr(1) % next)
      nullify(i1Ptr(1) % prev)
      nullify(i1Ptr(1) % sendList)
      nullify(i1Ptr(1) % recvList)
      nullify(i1Ptr(1) % copyList)
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nEdgesOnEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nEdgesOnEdge', i1Ptr)

! Define variable edgesOnCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'edgesOnCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnCell', i2Ptr)

! Define variable edgesOnEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'edgesOnEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges2'
      i2Ptr(1) % dimNames(2) = 'nEdges'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnEdge', i2Ptr)

! Define variable weightsOnEdge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'weightsOnEdge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges2'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'weightsOnEdge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'weightsOnEdge', r2Ptr)

! Define variable dvEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dvEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dvEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dvEdge', r1Ptr)

! Define variable dcEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dcEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dcEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dcEdge', r1Ptr)

! Define variable angleEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'angleEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'angleEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'angleEdge', r1Ptr)

! Define variable areaCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'areaCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'areaCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'areaCell', r1Ptr)

! Define variable invAreaCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'invAreaCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'invAreaCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'invAreaCell', r1Ptr)

! Define variable areaTriangle
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'areaTriangle'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'areaTriangle', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'areaTriangle', r1Ptr)

! Define variable edgeNormalVectors
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'edgeNormalVectors'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgeNormalVectors', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgeNormalVectors', r2Ptr)

! Define variable localVerticalUnitVectors
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'localVerticalUnitVectors'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'localVerticalUnitVectors', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'localVerticalUnitVectors', r2Ptr)

! Define variable cellTangentPlane
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'cellTangentPlane'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'R3'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nCells'
     r3Ptr(1) % defaultValue = 0.0
     nullify(r3Ptr(1) % array)
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellTangentPlane', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'cellTangentPlane', r3Ptr)

! Define variable cellsOnCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'cellsOnCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnCell', i2Ptr)

! Define variable verticesOnCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'verticesOnCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'verticesOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'verticesOnCell', i2Ptr)

! Define variable verticesOnEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'verticesOnEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'TWO'
      i2Ptr(1) % dimNames(2) = 'nEdges'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'verticesOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'verticesOnEdge', i2Ptr)

! Define variable edgesOnVertex
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'edgesOnVertex'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'vertexDegree'
      i2Ptr(1) % dimNames(2) = 'nVertices'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnVertex', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnVertex', i2Ptr)

! Define variable cellsOnVertex
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'cellsOnVertex'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'vertexDegree'
      i2Ptr(1) % dimNames(2) = 'nVertices'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnVertex', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnVertex', i2Ptr)

! Define variable kiteAreasOnVertex
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kiteAreasOnVertex'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'vertexDegree'
      r2Ptr(1) % dimNames(2) = 'nVertices'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'kiteAreasOnVertex', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kiteAreasOnVertex', r2Ptr)

! Define variable fEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fEdge', r1Ptr)

! Define variable fVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fVertex', r1Ptr)

! Define variable meshDensity
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshDensity'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshDensity', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshDensity', r1Ptr)

! Define variable meshScalingDel2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshScalingDel2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshScalingDel2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshScalingDel2', r1Ptr)

! Define variable meshScalingDel4
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshScalingDel4'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshScalingDel4', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshScalingDel4', r1Ptr)

! Define variable cf1
      allocate(r0Ptr(1))

! Setting up time level 1
      r0Ptr(1) % fieldName = 'cf1'
      r0Ptr(1) % isVarArray = .false.
      r0Ptr(1) % isDecomposed = .false.
      r0Ptr(1) % hasTimeDimension = .false.
     r0Ptr(1) % scalar = 0.0
      nullify(r0Ptr(1) % next)
      nullify(r0Ptr(1) % prev)
      nullify(r0Ptr(1) % sendList)
      nullify(r0Ptr(1) % recvList)
      nullify(r0Ptr(1) % copyList)
      r0Ptr(1) % block => block

      r0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf1', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf1', r0Ptr)

! Define variable cf2
      allocate(r0Ptr(1))

! Setting up time level 1
      r0Ptr(1) % fieldName = 'cf2'
      r0Ptr(1) % isVarArray = .false.
      r0Ptr(1) % isDecomposed = .false.
      r0Ptr(1) % hasTimeDimension = .false.
     r0Ptr(1) % scalar = 0.0
      nullify(r0Ptr(1) % next)
      nullify(r0Ptr(1) % prev)
      nullify(r0Ptr(1) % sendList)
      nullify(r0Ptr(1) % recvList)
      nullify(r0Ptr(1) % copyList)
      r0Ptr(1) % block => block

      r0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf2', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf2', r0Ptr)

! Define variable cf3
      allocate(r0Ptr(1))

! Setting up time level 1
      r0Ptr(1) % fieldName = 'cf3'
      r0Ptr(1) % isVarArray = .false.
      r0Ptr(1) % isDecomposed = .false.
      r0Ptr(1) % hasTimeDimension = .false.
     r0Ptr(1) % scalar = 0.0
      nullify(r0Ptr(1) % next)
      nullify(r0Ptr(1) % prev)
      nullify(r0Ptr(1) % sendList)
      nullify(r0Ptr(1) % recvList)
      nullify(r0Ptr(1) % copyList)
      r0Ptr(1) % block => block

      r0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf3', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf3', r0Ptr)

! Define variable cpr
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cpr'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'THREE'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cpr', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cpr', r2Ptr)

! Define variable cpl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cpl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'THREE'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cpl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cpl', r2Ptr)

! Define variable zgrid
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'zgrid'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zgrid', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zgrid', r2Ptr)

! Define variable rdzw
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rdzw'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rdzw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rdzw', r1Ptr)

! Define variable dzu
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dzu'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dzu', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dzu', r1Ptr)

! Define variable rdzu
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rdzu'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rdzu', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rdzu', r1Ptr)

! Define variable fzm
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fzm'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fzm', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fzm', r1Ptr)

! Define variable fzp
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fzp'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fzp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fzp', r1Ptr)

! Define variable zx
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'zx'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zx', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zx', r2Ptr)

! Define variable zz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'zz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zz', r2Ptr)

! Define variable zb
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'zb'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nEdges'
     r3Ptr(1) % defaultValue = 0.0
     nullify(r3Ptr(1) % array)
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zb', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb', r3Ptr)

! Define variable zb3
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'zb3'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nEdges'
     r3Ptr(1) % defaultValue = 0.0
     nullify(r3Ptr(1) % array)
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zb3', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb3', r3Ptr)

! Define variable pzm
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pzm'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pzm', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pzm', r2Ptr)

! Define variable pzp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pzp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pzp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pzp', r2Ptr)

! Define variable dss
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dss'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dss', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dss', r2Ptr)

! Define variable u_init
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'u_init'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'u_init', r1Ptr)

! Define variable t_init
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 't_init'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 't_init', r2Ptr)
      call mpas_pool_add_field(block % allFields, 't_init', r2Ptr)

! Define variable qv_init
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'qv_init'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qv_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qv_init', r1Ptr)

! Define variable deriv_two
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'deriv_two'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'FIFTEEN'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nEdges'
     r3Ptr(1) % defaultValue = 0.0
     nullify(r3Ptr(1) % array)
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'deriv_two', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'deriv_two', r3Ptr)

! Define variable adv_coefs
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'adv_coefs'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'FIFTEEN'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'adv_coefs', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'adv_coefs', r2Ptr)

! Define variable adv_coefs_3rd
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'adv_coefs_3rd'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'FIFTEEN'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'adv_coefs_3rd', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'adv_coefs_3rd', r2Ptr)

! Define variable advCellsForEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'advCellsForEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'FIFTEEN'
      i2Ptr(1) % dimNames(2) = 'nEdges'
     i2Ptr(1) % defaultValue = 0
     nullify(i2Ptr(1) % array)
      nullify(i2Ptr(1) % next)
      nullify(i2Ptr(1) % prev)
      nullify(i2Ptr(1) % sendList)
      nullify(i2Ptr(1) % recvList)
      nullify(i2Ptr(1) % copyList)
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'advCellsForEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'advCellsForEdge', i2Ptr)

! Define variable nAdvCellsForEdge
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nAdvCellsForEdge'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
     i1Ptr(1) % defaultValue = 0
     nullify(i1Ptr(1) % array)
      nullify(i1Ptr(1) % next)
      nullify(i1Ptr(1) % prev)
      nullify(i1Ptr(1) % sendList)
      nullify(i1Ptr(1) % recvList)
      nullify(i1Ptr(1) % copyList)
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nAdvCellsForEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nAdvCellsForEdge', i1Ptr)

! Define variable defc_a
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'defc_a'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'defc_a', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'defc_a', r2Ptr)

! Define variable defc_b
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'defc_b'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'defc_b', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'defc_b', r2Ptr)

! Define variable coeffs_reconstruct
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'coeffs_reconstruct'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'R3'
      r3Ptr(1) % dimNames(2) = 'maxEdges'
      r3Ptr(1) % dimNames(3) = 'nCells'
     r3Ptr(1) % defaultValue = 0.0
     nullify(r3Ptr(1) % array)
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'coeffs_reconstruct', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'coeffs_reconstruct', r3Ptr)

! Define variable east
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'east'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'east', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'east', r2Ptr)

! Define variable north
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'north'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'north', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'north', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
      end if

   end subroutine atm_generate_pool_mesh


   subroutine atm_generate_pool_state(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: cam5Active


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'cam5Active', cam5Active)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'state', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'state', newSubPool)

! Define var array scalars
      allocate(r3Ptr(2))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var qv
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var qc
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var qr
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
      if (.not. group_started) then
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
            call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
         end if
      else
         group_started = .false.
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
         end if
      end if
! End of group       
! Starting group conc
! Define constituent var qnr
! My Packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qnr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_qnr', -1)
           end if
      end if
! Define constituent var qni
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qni', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qni', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'conc_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_end', index_counter)
            end if
         end if
! End of group       
! Starting group aerosol
! Define constituent var aer1
! My Packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer1', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_aer1', -1)
           end if
      end if
! Define constituent var aer2
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer2', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer2', -1)
         end if
      end if
! Define constituent var aer3
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer3', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer3', -1)
         end if
      end if
! Define constituent var aer4
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer4', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer4', -1)
         end if
      end if
! Define constituent var aer5
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer5', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer5', -1)
         end if
      end if
! Define constituent var aer6
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer6', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer6', -1)
         end if
      end if
! Define constituent var aer7
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer7', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer7', -1)
         end if
      end if
! Define constituent var aer8
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer8', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer8', -1)
         end if
      end if
! Define constituent var aer9
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer9', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer9', -1)
         end if
      end if
! Define constituent var aer10
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer10', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer10', -1)
         end if
      end if
! Define constituent var aer11
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer11', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer11', -1)
         end if
      end if
! Define constituent var aer12
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer12', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer12', -1)
         end if
      end if
! Define constituent var aer13
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer13', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer13', -1)
         end if
      end if
! Define constituent var aer14
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer14', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer14', -1)
         end if
      end if
! Define constituent var aer15
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer15', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer15', -1)
         end if
      end if
! Define constituent var aer16
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer16', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer16', -1)
         end if
      end if
! Define constituent var aer17
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer17', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer17', -1)
         end if
      end if
! Define constituent var aer18
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer18', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer18', -1)
         end if
      end if
! Define constituent var aer19
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer19', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer19', -1)
         end if
      end if
!++CMZ
! Define constituent var aer21
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer21', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer21', -1)
         end if
      end if
! Define constituent var aer22
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer22', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer22', -1)
         end if
      end if
! Define constituent var aer23
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer23', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer23', -1)
         end if
      end if
! Define constituent var aer24
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer24', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer24', -1)
         end if
      end if
! Define constituent var aer25
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer25', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer25', -1)
         end if
      end if
! Define constituent var aer26
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer26', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer26', -1)
         end if
      end if
! Define constituent var aer27
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer27', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer27', -1)
         end if
      end if
! Define constituent var aer28
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer28', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer28', -1)
         end if
      end if
!--CMZ
! Define constituent var aer20
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer20', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer20', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'aerosol_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_end', index_counter)
            end if
         end if
! End of group       



      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars', numConstituents)
      end if

!++CMZ      
      write(stderrUnit,*) 'CMZ1 numConstituents: ', index_counter
!--CMZ

! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'scalars'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qnr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qnr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qni', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer1', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer2', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer3', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer4', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer5', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer5'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer6', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer6'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer7', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer7'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer8', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer8'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer9', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer9'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer10', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer10'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer11', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer11'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer12', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer12'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer13', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer13'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer14', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer14'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer15', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer15'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer16', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer16'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer17', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer17'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer18', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer18'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer19', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer19'
      end if


!++CMZ
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer21', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer21'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer22', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer22'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer23', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer23'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer24', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer24'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer25', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer25'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer26', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer26'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer27', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer27'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer28', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer28'
      end if
!--CMZ


      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer20', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer20'
      end if

!++CMZ      
      write(stderrUnit,*) 'CMZ2 numConstituents'
!--CMZ

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_scalars'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      nullify(r3Ptr(1) % array)
      r3Ptr(1) % defaultValue = 0.0
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block
! Defining time level 2
      allocate( r3Ptr(2) % constituentNames(numConstituents) )
      r3Ptr(2) % fieldName = 'scalars'
      r3Ptr(2) % isDecomposed = .true.
      r3Ptr(2) % hasTimeDimension = .true.
      r3Ptr(2) % isVarArray = .true.
      r3Ptr(2) % isPersistent = .true.
      r3Ptr(2) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qnr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qnr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qni', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer1', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer2', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer3', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer4', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer5', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer5'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer6', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer6'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer7', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer7'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer8', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer8'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer9', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer9'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer10', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer10'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer11', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer11'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer12', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer12'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer13', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer13'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer14', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer14'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer15', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer15'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer16', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer16'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer17', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer17'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer18', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer18'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer19', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer19'
      end if

!++CMZ
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer21', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer21'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer22', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer22'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer23', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer23'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer24', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer24'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer25', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer25'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer26', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer26'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer27', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer27'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer28', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer28'
      end if
!--CMZ

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer20', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'aer20'
      end if

! Setup dimensions for       
      r3Ptr(2) % dimNames(1) = 'num_scalars'
      r3Ptr(2) % dimNames(2) = 'nVertLevels'
      r3Ptr(2) % dimNames(3) = 'nCells'

      nullify(r3Ptr(2) % array)
      r3Ptr(2) % defaultValue = 0.0
      nullify(r3Ptr(2) % next)
      nullify(r3Ptr(2) % prev)
      nullify(r3Ptr(2) % sendList)
      nullify(r3Ptr(2) % recvList)
      nullify(r3Ptr(2) % copyList)
      r3Ptr(2) % block => block

            r3Ptr(1) % isActive = .true.
            r3Ptr(2) % isActive = .true.
            call mpas_pool_add_field(newSubPool, 'scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars', r3Ptr)

! Define variable xtime
      allocate(c0Ptr(2))

! Setting up time level 1
      c0Ptr(1) % fieldName = 'xtime'
      c0Ptr(1) % isVarArray = .false.
      c0Ptr(1) % isDecomposed = .false.
      c0Ptr(1) % hasTimeDimension = .true.
     c0Ptr(1) % scalar = ''
      nullify(c0Ptr(1) % next)
      nullify(c0Ptr(1) % prev)
      nullify(c0Ptr(1) % sendList)
      nullify(c0Ptr(1) % recvList)
      nullify(c0Ptr(1) % copyList)
      c0Ptr(1) % block => block

! Setting up time level 2
      c0Ptr(2) % fieldName = 'xtime'
      c0Ptr(2) % isVarArray = .false.
      c0Ptr(2) % isDecomposed = .false.
      c0Ptr(2) % hasTimeDimension = .true.
     c0Ptr(2) % scalar = ''
      nullify(c0Ptr(2) % next)
      nullify(c0Ptr(2) % prev)
      nullify(c0Ptr(2) % sendList)
      nullify(c0Ptr(2) % recvList)
      nullify(c0Ptr(2) % copyList)
      c0Ptr(2) % block => block

      c0Ptr(1) % isActive = .true.
      c0Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xtime', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'xtime', c0Ptr)

! Define variable u
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'u'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nEdges'
     r2Ptr(2) % defaultValue = 0.0
     nullify(r2Ptr(2) % array)
      nullify(r2Ptr(2) % next)
      nullify(r2Ptr(2) % prev)
      nullify(r2Ptr(2) % sendList)
      nullify(r2Ptr(2) % recvList)
      nullify(r2Ptr(2) % copyList)
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u', r2Ptr)

! Define variable w
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'w'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(2) % dimNames(2) = 'nCells'
     r2Ptr(2) % defaultValue = 0.0
     nullify(r2Ptr(2) % array)
      nullify(r2Ptr(2) % next)
      nullify(r2Ptr(2) % prev)
      nullify(r2Ptr(2) % sendList)
      nullify(r2Ptr(2) % recvList)
      nullify(r2Ptr(2) % copyList)
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w', r2Ptr)

! Define variable rho_zz
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_zz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'rho_zz'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
     r2Ptr(2) % defaultValue = 0.0
     nullify(r2Ptr(2) % array)
      nullify(r2Ptr(2) % next)
      nullify(r2Ptr(2) % prev)
      nullify(r2Ptr(2) % sendList)
      nullify(r2Ptr(2) % recvList)
      nullify(r2Ptr(2) % copyList)
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_zz', r2Ptr)

! Define variable theta_m
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_m'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'theta_m'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
     r2Ptr(2) % defaultValue = 0.0
     nullify(r2Ptr(2) % array)
      nullify(r2Ptr(2) % next)
      nullify(r2Ptr(2) % prev)
      nullify(r2Ptr(2) % sendList)
      nullify(r2Ptr(2) % recvList)
      nullify(r2Ptr(2) % copyList)
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_m', r2Ptr)

! Define variable m_ps
      allocate(r1Ptr(2))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'm_ps'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

! Setting up time level 2
      r1Ptr(2) % fieldName = 'm_ps'
      r1Ptr(2) % isVarArray = .false.
      r1Ptr(2) % isDecomposed = .true.
      r1Ptr(2) % hasTimeDimension = .true.
      r1Ptr(2) % isPersistent = .true.
      r1Ptr(2) % isActive = .false.
! Setting up dimensions
      r1Ptr(2) % dimNames(1) = 'nCells'
     r1Ptr(2) % defaultValue = 0.0
     nullify(r1Ptr(2) % array)
      nullify(r1Ptr(2) % next)
      nullify(r1Ptr(2) % prev)
      nullify(r1Ptr(2) % sendList)
      nullify(r1Ptr(2) % recvList)
      nullify(r1Ptr(2) % copyList)
      r1Ptr(2) % block => block

      r1Ptr(1) % isActive = .true.
      r1Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'm_ps', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'm_ps', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
      end if

   end subroutine atm_generate_pool_state


   subroutine atm_generate_pool_diag(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: cam5Active


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'cam5Active', cam5Active)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'diag', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'diag', newSubPool)

! Define variable cofrz
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cofrz'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofrz', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cofrz', r1Ptr)

! Define variable cofwr
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cofwr'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofwr', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cofwr', r2Ptr)

! Define variable cofwz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cofwz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofwz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cofwz', r2Ptr)

! Define variable coftz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'coftz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'coftz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'coftz', r2Ptr)

! Define variable cofwt
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cofwt'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofwt', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cofwt', r2Ptr)

! Define variable a_tri
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'a_tri'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'a_tri', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'a_tri', r2Ptr)

! Define variable alpha_tri
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'alpha_tri'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'alpha_tri', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'alpha_tri', r2Ptr)

! Define variable gamma_tri
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'gamma_tri'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gamma_tri', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gamma_tri', r2Ptr)

! Define variable pressure_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pressure_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pressure_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure_p', r2Ptr)

! Define variable rho
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho', r2Ptr)

! Define variable theta
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta', r2Ptr)

! Define variable relhum
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'relhum'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum', r2Ptr)

! Define variable v
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'v'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'v', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'v', r2Ptr)

! Define variable divergence
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'divergence'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'divergence', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'divergence', r2Ptr)

! Define variable vorticity
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'vorticity'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity', r2Ptr)

! Define variable pv_edge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_edge', r2Ptr)

! Define variable rho_edge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_edge', r2Ptr)

! Define variable ke
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ke'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ke', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ke', r2Ptr)

! Define variable pv_vertex
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_vertex'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_vertex', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_vertex', r2Ptr)

! Define variable pv_cell
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_cell'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_cell', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_cell', r2Ptr)

! Define variable uReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructX', r2Ptr)

! Define variable uReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructY', r2Ptr)

! Define variable uReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZ', r2Ptr)

! Define variable uReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZonal', r2Ptr)

! Define variable uReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructMeridional', r2Ptr)

! Define variable rv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rv', r2Ptr)

! Define variable circulation
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'circulation'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'circulation', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'circulation', r2Ptr)

! Define variable gradPVt
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'gradPVt'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gradPVt', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gradPVt', r2Ptr)

! Define variable gradPVn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'gradPVn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gradPVn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gradPVn', r2Ptr)

! Define variable h_divergence
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'h_divergence'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'h_divergence', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'h_divergence', r2Ptr)

! Define variable exner
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'exner'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'exner', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exner', r2Ptr)

! Define variable exner_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'exner_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'exner_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exner_base', r2Ptr)

! Define variable rtheta_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_base', r2Ptr)

! Define variable pressure
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pressure'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pressure', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure', r2Ptr)

! Define variable pressure_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pressure_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pressure_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure_base', r2Ptr)

! Define variable rho_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_base', r2Ptr)

! Define variable theta_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_base', r2Ptr)

! Define variable rho_zz_old_split
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_zz_old_split'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_zz_old_split', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_zz_old_split', r2Ptr)

! Define variable ruAvg
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ruAvg'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ruAvg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ruAvg', r2Ptr)

! Define variable wwAvg
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'wwAvg'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wwAvg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'wwAvg', r2Ptr)

! Define variable ruAvg_split
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ruAvg_split'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ruAvg_split', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ruAvg_split', r2Ptr)

! Define variable wwAvg_split
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'wwAvg_split'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wwAvg_split', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'wwAvg_split', r2Ptr)

! Define variable cqu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cqu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cqu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cqu', r2Ptr)

! Define variable cqw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cqw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cqw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cqw', r2Ptr)

! Define variable ru
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ru', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru', r2Ptr)

! Define variable ru_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ru_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_p', r2Ptr)

! Define variable ru_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ru_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_save', r2Ptr)

! Define variable rw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw', r2Ptr)

! Define variable rw_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rw_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_p', r2Ptr)

! Define variable rw_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rw_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_save', r2Ptr)

! Define variable rtheta_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_p', r2Ptr)

! Define variable rtheta_pp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_pp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_pp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_pp', r2Ptr)

! Define variable rtheta_p_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_p_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_p_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_p_save', r2Ptr)

! Define variable rtheta_pp_old
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_pp_old'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_pp_old', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_pp_old', r2Ptr)

! Define variable rho_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_p', r2Ptr)

! Define variable rho_pp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_pp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_pp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_pp', r2Ptr)

! Define variable rho_p_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_p_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_p_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_p_save', r2Ptr)

! Define variable kdiff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kdiff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'kdiff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kdiff', r2Ptr)

! Define variable surface_pressure
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'surface_pressure'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'surface_pressure', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'surface_pressure', r1Ptr)

! Define variable mslp
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'mslp'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'mslp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'mslp', r1Ptr)

! Define variable temperature_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_200hPa', r1Ptr)

! Define variable temperature_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_500hPa', r1Ptr)

! Define variable temperature_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_700hPa', r1Ptr)

! Define variable temperature_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_850hPa', r1Ptr)

! Define variable relhum_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_200hPa', r1Ptr)

! Define variable relhum_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_500hPa', r1Ptr)

! Define variable relhum_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_700hPa', r1Ptr)

! Define variable relhum_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_850hPa', r1Ptr)

! Define variable height_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_200hPa', r1Ptr)

! Define variable height_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_500hPa', r1Ptr)

! Define variable height_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_700hPa', r1Ptr)

! Define variable height_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_850hPa', r1Ptr)

! Define variable uzonal_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_200hPa', r1Ptr)

! Define variable uzonal_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_500hPa', r1Ptr)

! Define variable uzonal_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_700hPa', r1Ptr)

! Define variable uzonal_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_850hPa', r1Ptr)

! Define variable umeridional_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_200hPa', r1Ptr)

! Define variable umeridional_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_500hPa', r1Ptr)

! Define variable umeridional_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_700hPa', r1Ptr)

! Define variable umeridional_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_850hPa', r1Ptr)

! Define variable w_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_200hPa', r1Ptr)

! Define variable w_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_500hPa', r1Ptr)

! Define variable w_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_700hPa', r1Ptr)

! Define variable w_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_850hPa', r1Ptr)

! Define variable vorticity_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_200hPa', r1Ptr)

! Define variable vorticity_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_500hPa', r1Ptr)

! Define variable vorticity_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_700hPa', r1Ptr)

! Define variable vorticity_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_850hPa', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
      end if

   end subroutine atm_generate_pool_diag


   subroutine atm_generate_pool_tend(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: cam5Active


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'cam5Active', cam5Active)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'tend', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'tend', newSubPool)

! Define var array scalars_tend
      allocate(r3Ptr(1))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var tend_qv
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var tend_qc
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var tend_qr
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
      if (.not. group_started) then
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
            call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
         end if
      else
         group_started = .false.
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
         end if
      end if
! End of group       
! Starting group conc
! Define constituent var tend_qnr
! My Packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qnr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_qnr', -1)
           end if
      end if
! Define constituent var tend_qni
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qni', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qni', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'conc_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_end', index_counter)
            end if
         end if
! End of group       
! Starting group aerosol
! Define constituent var tend_aer1
! My Packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer1', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_aer1', -1)
           end if
      end if
! Define constituent var tend_aer2
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer2', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer2', -1)
         end if
      end if
! Define constituent var tend_aer3
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer3', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer3', -1)
         end if
      end if
! Define constituent var tend_aer4
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer4', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer4', -1)
         end if
      end if
! Define constituent var tend_aer5
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer5', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer5', -1)
         end if
      end if
! Define constituent var tend_aer6
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer6', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer6', -1)
         end if
      end if
! Define constituent var tend_aer7
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer7', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer7', -1)
         end if
      end if
! Define constituent var tend_aer8
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer8', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer8', -1)
         end if
      end if
! Define constituent var tend_aer9
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer9', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer9', -1)
         end if
      end if
! Define constituent var tend_aer10
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer10', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer10', -1)
         end if
      end if
! Define constituent var tend_aer11
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer11', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer11', -1)
         end if
      end if
! Define constituent var tend_aer12
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer12', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer12', -1)
         end if
      end if
! Define constituent var tend_aer13
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer13', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer13', -1)
         end if
      end if
! Define constituent var tend_aer14
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer14', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer14', -1)
         end if
      end if
! Define constituent var tend_aer15
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer15', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer15', -1)
         end if
      end if
! Define constituent var tend_aer16
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer16', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer16', -1)
         end if
      end if
! Define constituent var tend_aer17
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer17', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer17', -1)
         end if
      end if
! Define constituent var tend_aer18
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer18', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer18', -1)
         end if
      end if
! Define constituent var tend_aer19
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer19', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer19', -1)
         end if
      end if

!++CMZ

! Define constituent var tend_aer21
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer21', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer21', -1)
         end if
      end if

! Define constituent var tend_aer22
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer22', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer22', -1)
         end if
      end if

! Define constituent var tend_aer23
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer23', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer23', -1)
         end if
      end if

! Define constituent var tend_aer24
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer24', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer24', -1)
         end if
      end if

! Define constituent var tend_aer25
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer25', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer25', -1)
         end if
      end if

! Define constituent var tend_aer26
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer26', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer26', -1)
         end if
      end if

! Define constituent var tend_aer27
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer27', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer27', -1)
         end if
      end if

! Define constituent var tend_aer28
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer28', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer28', -1)
         end if
      end if


!--CMZ



! Define constituent var tend_aer20
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer20', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer20', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'aerosol_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars_tend', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'scalars_tend'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qnr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qnr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qni', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer1', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer2', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer3', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer4', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer5', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer5'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer6', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer6'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer7', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer7'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer8', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer8'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer9', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer9'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer10', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer10'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer11', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer11'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer12', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer12'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer13', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer13'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer14', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer14'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer15', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer15'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer16', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer16'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer17', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer17'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer18', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer18'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer19', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer19'
      end if

!++CMZ
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer21', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer21'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer22', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer22'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer23', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer23'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer24', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer24'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer25', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer25'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer26', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer26'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer27', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer27'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer28', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer28'
      end if
!--CMZ

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer20', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_aer20'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_scalars_tend'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      nullify(r3Ptr(1) % array)
      r3Ptr(1) % defaultValue = 0.0
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

            r3Ptr(1) % isActive = .true.
            call mpas_pool_add_field(newSubPool, 'scalars_tend', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars_tend', r3Ptr)

! Define variable tend_u
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_u', r2Ptr)

! Define variable tend_w
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_w', r2Ptr)

! Define variable tend_rho
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_rho'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_rho', r2Ptr)

! Define variable tend_theta
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta', r2Ptr)

! Define variable rt_diabatic_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rt_diabatic_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rt_diabatic_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rt_diabatic_tend', r2Ptr)

! Define variable euler_tend_u
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'euler_tend_u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u_euler', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'euler_tend_u', r2Ptr)

! Define variable euler_tend_w
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'euler_tend_w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_euler', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'euler_tend_w', r2Ptr)

! Define variable euler_tend_theta
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'euler_tend_theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_euler', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'euler_tend_theta', r2Ptr)

! Define variable tend_w_pgf
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_w_pgf'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_pgf', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_w_pgf', r2Ptr)

! Define variable tend_w_buoy
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_w_buoy'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_buoy', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_w_buoy', r2Ptr)

! Define variable tend_sfc_pressure
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'tend_sfc_pressure'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
     r1Ptr(1) % defaultValue = 0.0
     nullify(r1Ptr(1) % array)
      nullify(r1Ptr(1) % next)
      nullify(r1Ptr(1) % prev)
      nullify(r1Ptr(1) % sendList)
      nullify(r1Ptr(1) % recvList)
      nullify(r1Ptr(1) % copyList)
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tend_sfc_pressure', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_sfc_pressure', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
      end if

   end subroutine atm_generate_pool_tend


   subroutine atm_generate_pool_tend_physics(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: cam5Active


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'cam5Active', cam5Active)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'tend_physics', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'tend_physics', newSubPool)

! Define var array scalars_cam_tend
      allocate(r3Ptr(1))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var qv_cam_tend
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var qc_cam_tend
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var qr_cam_tend
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
      if (.not. group_started) then
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
            call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
         end if
      else
         group_started = .false.
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
         end if
      end if
! End of group       
! Starting group conc
! Define constituent var qnr_cam_tend
! My Packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qnr_c', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_qnr_c', -1)
           end if
      end if
! Define constituent var qni_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qni_c', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qni_c', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'conc_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'conc_end', index_counter)
            end if
         end if
! End of group       
! Starting group aerosol
! Define constituent var aer1_cam_tend
! My Packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer1', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_aer1', -1)
           end if
      end if
! Define constituent var aer2_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer2', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer2', -1)
         end if
      end if
! Define constituent var aer3_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer3', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer3', -1)
         end if
      end if
! Define constituent var aer4_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer4', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer4', -1)
         end if
      end if
! Define constituent var aer5_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer5', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer5', -1)
         end if
      end if
! Define constituent var aer6_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer6', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer6', -1)
         end if
      end if
! Define constituent var aer7_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer7', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer7', -1)
         end if
      end if
! Define constituent var aer8_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer8', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer8', -1)
         end if
      end if
! Define constituent var aer9_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer9', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer9', -1)
         end if
      end if
! Define constituent var aer10_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer10', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer10', -1)
         end if
      end if
! Define constituent var aer11_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer11', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer11', -1)
         end if
      end if
! Define constituent var aer12_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer12', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer12', -1)
         end if
      end if
! Define constituent var aer13_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer13', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer13', -1)
         end if
      end if
! Define constituent var aer14_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer14', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer14', -1)
         end if
      end if
! Define constituent var aer15_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer15', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer15', -1)
         end if
      end if
! Define constituent var aer16_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer16', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer16', -1)
         end if
      end if
! Define constituent var aer17_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer17', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer17', -1)
         end if
      end if
! Define constituent var aer18_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer18', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer18', -1)
         end if
      end if
! Define constituent var aer19_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer19', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer19', -1)
         end if
      end if

!++CMZ

! Define constituent var aer21_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer21', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer21', -1)
         end if
      end if

! Define constituent var aer22_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer22', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer22', -1)
         end if
      end if

! Define constituent var aer23_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer23', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer23', -1)
         end if
      end if

! Define constituent var aer24_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer24', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer24', -1)
         end if
      end if

! Define constituent var aer25_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer25', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer25', -1)
         end if
      end if

! Define constituent var aer26_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer26', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer26', -1)
         end if
      end if

! Define constituent var aer27_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer27', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer27', -1)
         end if
      end if

! Define constituent var aer28_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer28', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer28', -1)
         end if
      end if

!--CMZ


! Define constituent var aer20_cam_tend
! My packages are cam5
      if (cam5Active) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer20', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_aer20', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'aerosol_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'aerosol_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars_cam_tend', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'scalars_cam_tend'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qv_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qc_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qr_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qnr_c', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qnr_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qni_c', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qni_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer1', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer1_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer2', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer2_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer3', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer3_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer4', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer4_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer5', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer5_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer6', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer6_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer7', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer7_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer8', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer8_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer9', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer9_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer10', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer10_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer11', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer11_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer12', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer12_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer13', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer13_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer14', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer14_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer15', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer15_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer16', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer16_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer17', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer17_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer18', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer18_cam_tend'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer19', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer19_cam_tend'
      end if

!++CMZ
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer21', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer21_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer22', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer22_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer23', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer23_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer24', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer24_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer25', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer25_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer26', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer26_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer27', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer27_cam_tend'
      end if

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer28', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer28_cam_tend'
      end if
!--CMZ
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_aer20', const_index)
      end if
      if (index_counter > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'aer20_cam_tend'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_scalars_cam_tend'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      nullify(r3Ptr(1) % array)
      r3Ptr(1) % defaultValue = 0.0
      nullify(r3Ptr(1) % next)
      nullify(r3Ptr(1) % prev)
      nullify(r3Ptr(1) % sendList)
      nullify(r3Ptr(1) % recvList)
      nullify(r3Ptr(1) % copyList)
      r3Ptr(1) % block => block

            r3Ptr(1) % isActive = .true.
            call mpas_pool_add_field(newSubPool, 'scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars_cam_tend', r3Ptr)

! Define variable u_cam_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_cam_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_cam_tend', r2Ptr)

! Define variable ux_cam_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ux_cam_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ux', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ux_cam_tend', r2Ptr)

! Define variable uy_cam_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uy_cam_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uy', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uy_cam_tend', r2Ptr)

! Define variable theta_cam_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_cam_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
     r2Ptr(1) % defaultValue = 0.0
     nullify(r2Ptr(1) % array)
      nullify(r2Ptr(1) % next)
      nullify(r2Ptr(1) % prev)
      nullify(r2Ptr(1) % sendList)
      nullify(r2Ptr(1) % recvList)
      nullify(r2Ptr(1) % copyList)
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_cam_tend', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
      end if

   end subroutine atm_generate_pool_tend_physics


   subroutine atm_generate_structs(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_io_units
      type (block_type), pointer, intent(inout) :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      call atm_generate_pool_mesh(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_state(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_diag(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_tend(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_tend_physics(block, structPool, dimensionPool, packagePool)

   end subroutine atm_generate_structs

   function atm_setup_namelists(configPool, namelistFilename, dminfo) result(iErr)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      character (len=*), intent(in) :: namelistFilename
      type (dm_info), intent(in) :: dminfo
      integer :: iErr

      integer :: unitNumber
      logical :: nmlExists

      iErr = 0
      unitNumber = 21
      write(stderrUnit, *) 'Reading namelist from file '//trim(namelistFilename)
      inquire(file=trim(namelistFilename), exist=nmlExists)
      if ( .not. nmlExists ) then
         call mpas_dmpar_global_abort('ERROR: Namelist file '//trim(namelistFilename)//' does not exist.')
      end if
      open(unitNumber,file=trim(namelistFilename),status='old',form='formatted')

      call atm_setup_nmlrec_nhyd_model(configPool, unitNumber, dminfo)
      call atm_setup_nmlrec_damping(configPool, unitNumber, dminfo)
      call atm_setup_nmlrec_io(configPool, unitNumber, dminfo)
      call atm_setup_nmlrec_decomposition(configPool, unitNumber, dminfo)
      call atm_setup_nmlrec_restart(configPool, unitNumber, dminfo)
      call atm_setup_nmlrec_printout(configPool, unitNumber, dminfo)
      call atm_setup_nmlrec_physics(configPool, unitNumber, dminfo)

      close(unitNumber)
   end function atm_setup_namelists

   subroutine atm_setup_nmlrec_nhyd_model(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      character (len=StrKIND) :: config_time_integration = 'SRK3'
      real (kind=RKIND) :: config_dt = 720.000000
      character (len=StrKIND) :: config_calendar_type = 'gregorian'
      character (len=StrKIND) :: config_start_time = '2010-10-23_00:00:00'
      character (len=StrKIND) :: config_stop_time = 'none'
      character (len=StrKIND) :: config_run_duration = '5_00:00:00'
      logical :: config_split_dynamics_transport = .true.
      integer :: config_number_of_sub_steps = 2
      integer :: config_dynamics_split_steps = 3
      real (kind=RKIND) :: config_h_mom_eddy_visc2 = 0.000000
      real (kind=RKIND) :: config_h_mom_eddy_visc4 = 0.000000
      real (kind=RKIND) :: config_v_mom_eddy_visc2 = 0.000000
      real (kind=RKIND) :: config_h_theta_eddy_visc2 = 0.000000
      real (kind=RKIND) :: config_h_theta_eddy_visc4 = 0.000000
      real (kind=RKIND) :: config_v_theta_eddy_visc2 = 0.000000
      character (len=StrKIND) :: config_horiz_mixing = '2d_smagorinsky'
      real (kind=RKIND) :: config_len_disp = 120000.000000
      real (kind=RKIND) :: config_visc4_2dsmag = 0.050000
      real (kind=RKIND) :: config_del4u_div_factor = 1.000000
      integer :: config_w_adv_order = 3
      integer :: config_theta_adv_order = 3
      integer :: config_scalar_adv_order = 3
      integer :: config_u_vadv_order = 3
      integer :: config_w_vadv_order = 3
      integer :: config_theta_vadv_order = 3
      integer :: config_scalar_vadv_order = 3
      logical :: config_scalar_advection = .true.
      logical :: config_positive_definite = .false.
      logical :: config_monotonic = .true.
      real (kind=RKIND) :: config_coef_3rd_order = 0.250000
      logical :: config_mix_full = .true.
      real (kind=RKIND) :: config_epssm = 0.100000
      real (kind=RKIND) :: config_smdiv = 0.100000
      logical :: config_newpx = .false.
      real (kind=RKIND) :: config_apvm_upwinding = 0.500000
      logical :: config_h_ScaleWithMesh = .true.
      integer :: config_num_halos = 2

      namelist /nhyd_model/ &
         config_time_integration, &
         config_dt, &
         config_calendar_type, &
         config_start_time, &
         config_stop_time, &
         config_run_duration, &
         config_split_dynamics_transport, &
         config_number_of_sub_steps, &
         config_dynamics_split_steps, &
         config_h_mom_eddy_visc2, &
         config_h_mom_eddy_visc4, &
         config_v_mom_eddy_visc2, &
         config_h_theta_eddy_visc2, &
         config_h_theta_eddy_visc4, &
         config_v_theta_eddy_visc2, &
         config_horiz_mixing, &
         config_len_disp, &
         config_visc4_2dsmag, &
         config_del4u_div_factor, &
         config_w_adv_order, &
         config_theta_adv_order, &
         config_scalar_adv_order, &
         config_u_vadv_order, &
         config_w_vadv_order, &
         config_theta_vadv_order, &
         config_scalar_vadv_order, &
         config_scalar_advection, &
         config_positive_definite, &
         config_monotonic, &
         config_coef_3rd_order, &
         config_mix_full, &
         config_epssm, &
         config_smdiv, &
         config_newpx, &
         config_apvm_upwinding, &
         config_h_ScaleWithMesh, &
         config_num_halos
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, nhyd_model, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record nhyd_model.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record nhyd_model not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_char(dminfo, config_time_integration)
         call mpas_dmpar_bcast_real(dminfo, config_dt)
         call mpas_dmpar_bcast_char(dminfo, config_calendar_type)
         call mpas_dmpar_bcast_char(dminfo, config_start_time)
         call mpas_dmpar_bcast_char(dminfo, config_stop_time)
         call mpas_dmpar_bcast_char(dminfo, config_run_duration)
         call mpas_dmpar_bcast_logical(dminfo, config_split_dynamics_transport)
         call mpas_dmpar_bcast_int(dminfo, config_number_of_sub_steps)
         call mpas_dmpar_bcast_int(dminfo, config_dynamics_split_steps)
         call mpas_dmpar_bcast_real(dminfo, config_h_mom_eddy_visc2)
         call mpas_dmpar_bcast_real(dminfo, config_h_mom_eddy_visc4)
         call mpas_dmpar_bcast_real(dminfo, config_v_mom_eddy_visc2)
         call mpas_dmpar_bcast_real(dminfo, config_h_theta_eddy_visc2)
         call mpas_dmpar_bcast_real(dminfo, config_h_theta_eddy_visc4)
         call mpas_dmpar_bcast_real(dminfo, config_v_theta_eddy_visc2)
         call mpas_dmpar_bcast_char(dminfo, config_horiz_mixing)
         call mpas_dmpar_bcast_real(dminfo, config_len_disp)
         call mpas_dmpar_bcast_real(dminfo, config_visc4_2dsmag)
         call mpas_dmpar_bcast_real(dminfo, config_del4u_div_factor)
         call mpas_dmpar_bcast_int(dminfo, config_w_adv_order)
         call mpas_dmpar_bcast_int(dminfo, config_theta_adv_order)
         call mpas_dmpar_bcast_int(dminfo, config_scalar_adv_order)
         call mpas_dmpar_bcast_int(dminfo, config_u_vadv_order)
         call mpas_dmpar_bcast_int(dminfo, config_w_vadv_order)
         call mpas_dmpar_bcast_int(dminfo, config_theta_vadv_order)
         call mpas_dmpar_bcast_int(dminfo, config_scalar_vadv_order)
         call mpas_dmpar_bcast_logical(dminfo, config_scalar_advection)
         call mpas_dmpar_bcast_logical(dminfo, config_positive_definite)
         call mpas_dmpar_bcast_logical(dminfo, config_monotonic)
         call mpas_dmpar_bcast_real(dminfo, config_coef_3rd_order)
         call mpas_dmpar_bcast_logical(dminfo, config_mix_full)
         call mpas_dmpar_bcast_real(dminfo, config_epssm)
         call mpas_dmpar_bcast_real(dminfo, config_smdiv)
         call mpas_dmpar_bcast_logical(dminfo, config_newpx)
         call mpas_dmpar_bcast_real(dminfo, config_apvm_upwinding)
         call mpas_dmpar_bcast_logical(dminfo, config_h_ScaleWithMesh)
         call mpas_dmpar_bcast_int(dminfo, config_num_halos)
      end if

      call mpas_pool_add_config(configPool, 'config_time_integration', config_time_integration)
      call mpas_pool_add_config(configPool, 'config_dt', config_dt)
      call mpas_pool_add_config(configPool, 'config_calendar_type', config_calendar_type)
      call mpas_pool_add_config(configPool, 'config_start_time', config_start_time)
      call mpas_pool_add_config(configPool, 'config_stop_time', config_stop_time)
      call mpas_pool_add_config(configPool, 'config_run_duration', config_run_duration)
      call mpas_pool_add_config(configPool, 'config_split_dynamics_transport', config_split_dynamics_transport)
      call mpas_pool_add_config(configPool, 'config_number_of_sub_steps', config_number_of_sub_steps)
      call mpas_pool_add_config(configPool, 'config_dynamics_split_steps', config_dynamics_split_steps)
      call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc2', config_h_mom_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)
      call mpas_pool_add_config(configPool, 'config_v_mom_eddy_visc2', config_v_mom_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc2', config_h_theta_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc4', config_h_theta_eddy_visc4)
      call mpas_pool_add_config(configPool, 'config_v_theta_eddy_visc2', config_v_theta_eddy_visc2)
      call mpas_pool_add_config(configPool, 'config_horiz_mixing', config_horiz_mixing)
      call mpas_pool_add_config(configPool, 'config_len_disp', config_len_disp)
      call mpas_pool_add_config(configPool, 'config_visc4_2dsmag', config_visc4_2dsmag)
      call mpas_pool_add_config(configPool, 'config_del4u_div_factor', config_del4u_div_factor)
      call mpas_pool_add_config(configPool, 'config_w_adv_order', config_w_adv_order)
      call mpas_pool_add_config(configPool, 'config_theta_adv_order', config_theta_adv_order)
      call mpas_pool_add_config(configPool, 'config_scalar_adv_order', config_scalar_adv_order)
      call mpas_pool_add_config(configPool, 'config_u_vadv_order', config_u_vadv_order)
      call mpas_pool_add_config(configPool, 'config_w_vadv_order', config_w_vadv_order)
      call mpas_pool_add_config(configPool, 'config_theta_vadv_order', config_theta_vadv_order)
      call mpas_pool_add_config(configPool, 'config_scalar_vadv_order', config_scalar_vadv_order)
      call mpas_pool_add_config(configPool, 'config_scalar_advection', config_scalar_advection)
      call mpas_pool_add_config(configPool, 'config_positive_definite', config_positive_definite)
      call mpas_pool_add_config(configPool, 'config_monotonic', config_monotonic)
      call mpas_pool_add_config(configPool, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_add_config(configPool, 'config_mix_full', config_mix_full)
      call mpas_pool_add_config(configPool, 'config_epssm', config_epssm)
      call mpas_pool_add_config(configPool, 'config_smdiv', config_smdiv)
      call mpas_pool_add_config(configPool, 'config_newpx', config_newpx)
      call mpas_pool_add_config(configPool, 'config_apvm_upwinding', config_apvm_upwinding)
      call mpas_pool_add_config(configPool, 'config_h_ScaleWithMesh', config_h_ScaleWithMesh)
      call mpas_pool_add_config(configPool, 'config_num_halos', config_num_halos)

   end subroutine atm_setup_nmlrec_nhyd_model


   subroutine atm_setup_nmlrec_damping(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      real (kind=RKIND) :: config_zd = 22000.000000
      real (kind=RKIND) :: config_xnutr = 0.200000
      logical :: config_add_diffusive_damping = .true.

      namelist /damping/ &
         config_zd, &
         config_xnutr, &
         config_add_diffusive_damping
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, damping, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record damping.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record damping not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_real(dminfo, config_zd)
         call mpas_dmpar_bcast_real(dminfo, config_xnutr)
         call mpas_dmpar_bcast_logical(dminfo, config_add_diffusive_damping)
      end if

      call mpas_pool_add_config(configPool, 'config_zd', config_zd)
      call mpas_pool_add_config(configPool, 'config_xnutr', config_xnutr)
      call mpas_pool_add_config(configPool, 'config_add_diffusive_damping', config_add_diffusive_damping)

   end subroutine atm_setup_nmlrec_damping


   subroutine atm_setup_nmlrec_io(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      character (len=StrKIND) :: config_restart_timestamp_name = 'restart_timestamp'
      integer :: config_pio_num_iotasks = 0
      integer :: config_pio_stride = 1

      namelist /io/ &
         config_restart_timestamp_name, &
         config_pio_num_iotasks, &
         config_pio_stride
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, io, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record io.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record io not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_char(dminfo, config_restart_timestamp_name)
         call mpas_dmpar_bcast_int(dminfo, config_pio_num_iotasks)
         call mpas_dmpar_bcast_int(dminfo, config_pio_stride)
      end if

      call mpas_pool_add_config(configPool, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_add_config(configPool, 'config_pio_num_iotasks', config_pio_num_iotasks)
      call mpas_pool_add_config(configPool, 'config_pio_stride', config_pio_stride)

   end subroutine atm_setup_nmlrec_io


   subroutine atm_setup_nmlrec_decomposition(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      character (len=StrKIND) :: config_block_decomp_file_prefix = 'x1.40962.graph.info.part.'
      integer :: config_number_of_blocks = 0
      logical :: config_explicit_proc_decomp = .false.
      character (len=StrKIND) :: config_proc_decomp_file_prefix = 'graph.info.part.'

      namelist /decomposition/ &
         config_block_decomp_file_prefix, &
         config_number_of_blocks, &
         config_explicit_proc_decomp, &
         config_proc_decomp_file_prefix
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, decomposition, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record decomposition.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record decomposition not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_char(dminfo, config_block_decomp_file_prefix)
         call mpas_dmpar_bcast_int(dminfo, config_number_of_blocks)
         call mpas_dmpar_bcast_logical(dminfo, config_explicit_proc_decomp)
         call mpas_dmpar_bcast_char(dminfo, config_proc_decomp_file_prefix)
      end if

      call mpas_pool_add_config(configPool, 'config_block_decomp_file_prefix', config_block_decomp_file_prefix)
      call mpas_pool_add_config(configPool, 'config_number_of_blocks', config_number_of_blocks)
      call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
      call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)

   end subroutine atm_setup_nmlrec_decomposition


   subroutine atm_setup_nmlrec_restart(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      logical :: config_do_restart = .false.
      logical :: config_do_DAcycling = .false.

      namelist /restart/ &
         config_do_restart, &
         config_do_DAcycling
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, restart, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record restart.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record restart not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_do_restart)
         call mpas_dmpar_bcast_logical(dminfo, config_do_DAcycling)
      end if

      call mpas_pool_add_config(configPool, 'config_do_restart', config_do_restart)
      call mpas_pool_add_config(configPool, 'config_do_DAcycling', config_do_DAcycling)

   end subroutine atm_setup_nmlrec_restart


   subroutine atm_setup_nmlrec_printout(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      logical :: config_print_global_minmax_vel = .true.
      logical :: config_print_global_minmax_sca = .false.

      namelist /printout/ &
         config_print_global_minmax_vel, &
         config_print_global_minmax_sca
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, printout, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record printout.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record printout not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_logical(dminfo, config_print_global_minmax_vel)
         call mpas_dmpar_bcast_logical(dminfo, config_print_global_minmax_sca)
      end if

      call mpas_pool_add_config(configPool, 'config_print_global_minmax_vel', config_print_global_minmax_vel)
      call mpas_pool_add_config(configPool, 'config_print_global_minmax_sca', config_print_global_minmax_sca)

   end subroutine atm_setup_nmlrec_printout


   subroutine atm_setup_nmlrec_physics(configPool, unitNumber, dminfo)
      use mpas_derived_types
      use mpas_dmpar
      use mpas_pool_routines
      use mpas_io_units
      type (mpas_pool_type), intent(inout) :: configPool
      integer, intent(in) :: unitNumber
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), pointer :: recordPool

      character (len=StrKIND) :: config_cam_physics = 'cam5'

      namelist /physics/ &
         config_cam_physics
      if (dminfo % my_proc_id == IO_NODE) then
         rewind(unitNumber)
         read(unitNumber, physics, iostat=ierr)
         if (ierr > 0) then
            write(stderrUnit, *) 'Error while reading namelist record physics.'
            call mpas_dmpar_abort(dminfo)
         else if (ierr < 0) then
            write(stderrUnit,*) 'Namelist record physics not found; using default values for variables in this namelist'
         end if
      end if
      call mpas_dmpar_bcast_int(dminfo, ierr)

      if (ierr == 0) then
         call mpas_dmpar_bcast_char(dminfo, config_cam_physics)
      end if

      call mpas_pool_add_config(configPool, 'config_cam_physics', config_cam_physics)

   end subroutine atm_setup_nmlrec_physics



end module atm_core_interface
