module dyn_grid

!----------------------------------------------------------------------- 
! 
! Definition of dynamics computational grid.
!
! Entry points:
!      get_block_ldof_d
!      blocks_init               determine number of blocks
!      get_block_bounds_d        get first and last indices in global 
!                                block ordering
!      get_block_gcol_d          get column indices for given block
!      get_block_gcol_cnt_d      get number of columns in given block
!      get_block_lvl_cnt_d       get number of vertical levels in column
!      get_block_levels_d        get vertical levels in column
!      get_gcol_block_d          get global block indices and local columns 
!                                index for given global column index
!      get_gcol_block_cnt_d      get number of blocks containing data
!                                from a given global column index
!      get_block_owner_d         get process "owning" given block
!      get_horiz_grid_d          get horizontal grid coordinates
!      get_horiz_grid_dim_d      get horizontal dimensions of dynamics grid
!      dyn_grid_get_pref         get reference pressures for the dynamics grid
!      get_dyn_grid_parm_real2d
!      get_dyn_grid_parm_real1d
!      get_dyn_grid_parm
!      dyn_grid_find_gcol
!      dyn_grid_find_gcols
!      dyn_grid_get_elem_coords
!      dyn_grid_get_colndx       get global lat and lon coordinate and MPI process indices 
!                                corresponding to a specified global column index
!
! Originally from the hydrostatic MPAS version by Michael Wickett and Art Mirin
! 
!-----------------------------------------------------------------------

   use shr_kind_mod,      only : r8 => shr_kind_r8
   use physconst,         only : pi
   use pmgrid,            only : plev, plevp, nEdges, numcols, maxEdges, nVertices

   use spmd_utils,        only : iam, masterproc, npes
   use cam_logfile,       only : iulog
   use cam_abortutils,    only : endrun

   use ref_pres,          only : ref_pres_init

#if ( defined SPMD )
   use mpishorthand, only : mpiint, mpicom, mpir8
#endif

   implicit none
!   private
!   save
!
   integer, public, parameter :: ptimelevels = 2

   integer, public :: nblocks_per_pe, nblocks_tot, max_col_per_block

   integer, public, allocatable :: col_indices_in_block(:,:)  !  global column indices
   integer, public, allocatable :: num_col_per_block(:)
   integer, public, allocatable :: global_blockid(:)
   integer, public, allocatable :: local_col_index(:)  !  local to block

   ! The MPAS dynamics grids
   integer, parameter, public :: dyn_decomp      = 101

!   public :: &
!      dyn_grid_init, &
!      blocks_init, &
!      get_dyn_grid_parm, &
!      get_horiz_grid_d, &
!      get_horiz_grid_dim_d, &
!      get_dyn_grid_parm_real1d, &
!      get_dyn_grid_parm_real2d, &
!      get_block_gcol_d, &
!      get_block_gcol_cnt_d, &
!      get_block_lvl_cnt_d, &
!      get_block_levels_d, &
!      get_block_bounds_d, &
!      get_gcol_block_cnt_d, &
!      get_block_owner_d, &
!      get_gcol_block_d, &
!      get_block_ldof_d, &
!      dyn_grid_find_gcol, &
!      dyn_grid_find_gcols, &
!      dyn_grid_get_elem_coords, &
!      dyn_grid_get_pref, & 
!      dyn_grid_get_colndx,&
!      physgrid_copy_attributes_d

   real(r8), parameter :: rad2deg=180.0_r8/pi ! convert radians to degrees

   !
   ! moved from pmgrid
   !
   integer, public :: nCells

   !
   ! The following has been moved here from commap.  
   ! The 2D versions of clon and londeg are for backwards compatibility. 
   ! Needs to be removed from cam_history and the other dycores all at once.
   !
   real(r8), pointer:: clon(:,:)   ! model longitudes (radians)
   real(r8), pointer:: londeg(:,:) ! model longitudes (degrees)

   !
   ! need to figure out whether these global arrays are really necessary
   !
   real(r8), pointer :: w(:)           ! area weights
   real(r8), pointer:: latdeg(:)       ! model latitudes (degrees)

   !
   ! module data moved here from global_grid
   !
   real(r8), dimension(:), pointer :: lonCell               ! global cell longitudes
   real(r8), dimension(:), pointer :: latCell               ! global cell latitudes
   real(r8), dimension(:), pointer :: latVertex             ! global vertex latitude
   real(r8), dimension(:), pointer :: lonVertex             ! global vertex longitude
   real(r8), dimension(:,:), pointer, public :: bounds_lon  ! global bounds longitude
   real(r8), dimension(:,:), pointer, public :: bounds_lat  ! global bounds latitude
   real(r8), dimension(:), pointer :: areaCell              ! global cell areas

   integer, dimension(:),  pointer, public :: nCells_per_proc        ! no. cells per process
   integer, dimension(:),  pointer, public :: global_cell_owner      ! owner of global cell
   integer, dimension(:),  pointer, public :: global_to_local_cell   ! local index of global cell
   integer, dimension(:),  pointer :: nEdgesOnCell           ! number of edges of cell
   integer, dimension(:,:),pointer :: verticesOnCell         ! index of vertices on cell

   integer, dimension(:), pointer, public :: gcol            ! global index of local cell
   real(r8), dimension(:), pointer :: llon                   ! local cell longitudes
   real(r8), dimension(:), pointer :: llat                   ! local cell latitudes
   real(r8), dimension(:), pointer :: larea => null()        ! local cell area (m2) !jjang (debug3,debug)


contains

   !subroutine dyn_grid_init( gridfile )
   subroutine dyn_grid_init
   
      use mpas_cam_interface, only : mpas_decomp, mpas_global_dimensions, mpas_global_coordinates
      use pmgrid, only : MPAS_SPHERE_RAD
   
      use time_manager,     only: get_step_size
      use mpas_cam_interface, only: mpas_init1

      use cam_control_mod,   only : initial_run, restart_run, branch_run
      use cam_initfiles,     only : initial_file_get_id
      use pio,              only: file_desc_t
        
      ! Initialize mpas grid
      !character (len=*), intent(in) :: gridfile        ! name of global grid file
      !character (len=*) :: gridfile        ! name of global grid file
   
      type(file_desc_t), pointer :: fh_ini

      integer :: i, owner, iunit, istatus, ierr, j
      integer :: ncid, idlon, idlat, idarea, idcells, idedges, idvertices, idmaxedges
      integer :: idlatvertex, idlonvertex, idverticesonCell, idnedgesonCell
   
      real(r8) :: pi
   
      character(len=*), parameter :: subname = "dyn_grid::dyn_grid_init"
      character(len=256) :: partition_filename

      ! arguments
      real(r8) :: pref_edge(plev+1) ! reference pressure at layer edges (Pa)
      real(r8) :: pref_mid(plev)    ! reference pressure at layer midpoints (Pa)
      integer  :: num_pr_lev        ! number of top levels using pure pressure representation

      real(kind=r8) :: dtime
   
      pi = 4.0_r8*atan(1.0_r8)
   
      !
      ! open initial file
      !

      ! get filehandle pointer to primary restart file
      fh_ini => initial_file_get_id()

      dtime = get_step_size()
      !call mpas_init1(mpicom, real(dtime,8), initial_run, fh_ini)
      call mpas_init1(mpicom, real(dtime,8), fh_ini)

      !
      ! Query global grid dimensions from MPAS
      !
      call MPAS_GLOBAL_DIMENSIONS(nCells, nEdges, nVertices, maxEdges)
   
      allocate( lonCell(nCells) )
      allocate( latCell(nCells) )
      allocate( areaCell(nCells) )
   
      allocate (lonVertex(nVertices))
      allocate (latVertex(nVertices))
   
      allocate (nEdgesOnCell(nCells))
      allocate (verticesOnCell(maxEdges, nCells))
   
      allocate( bounds_lon(maxEdges, nCells))
      allocate( bounds_lat(maxEdges, nCells))
   
      allocate( global_cell_owner( nCells ) )
      allocate( global_to_local_cell( nCells ) )
      allocate( nCells_per_proc( npes ) )
   
      !
      ! Query location information from MPAS
      !
      call MPAS_GLOBAL_COORDINATES(nCells, nVertices, maxEdges, latCell, lonCell, areaCell, latVertex, lonVertex, nEdgesOnCell, verticesOnCell)
   
      do i = 1, nCells
         if (lonCell(i) < 0.0) lonCell(i) = lonCell(i) + 2.0_r8*pi
      enddo
   
      ! Calculate bounds. These will be written to history file to enable visualization.
      ! Bounds are the latitude and longitude of each vertex surrounding a cell.
   
      do i = 1, nCells
         do j = 1, maxEdges
   
            if (j <= nEdgesOnCell(i)) then
               bounds_lon(j,i) = lonVertex(verticesOnCell(j,i)) * 180._r8 / pi
               if (bounds_lon(j,i) < 0._r8) bounds_lon(j,i) = 360._r8 + bounds_lon(j,i)
               if ((lonCell(i)*180._r8/pi > 300._r8) .and. (bounds_lon(j,i) < 100._r8)) &
                  bounds_lon(j,i) = bounds_lon(j,i) + 360._r8
               if ((lonCell(i)*180._r8/pi < 100._r8) .and. (bounds_lon(j,i) > 300._r8)) &
                  bounds_lon(j,i) = bounds_lon(j,i) - 360._r8
   
               bounds_lat(j,i) = latVertex(verticesOnCell(j,i)) * 180._r8 / pi
            else
               bounds_lon(j,i) = 1.0e20_r8
               bounds_lat(j,i) = 1.0e20_r8
            end if
   
         end do
      end do
   
      deallocate( lonVertex )
      deallocate( latVertex )
      deallocate( nEdgesOnCell )
      deallocate( verticesOnCell )
   
      !
      ! Query mesh decomposition information from MPAS
      !
      call MPAS_DECOMP(nCells, npes, nCells_per_proc, global_cell_owner, global_to_local_cell)
   
      !
      ! Compute local arrays
      !
      numcols = nCells_per_proc(iam+1)
      allocate( gcol(numcols) )
      allocate( llon(numcols) )
      allocate( llat(numcols) )
      allocate( larea(numcols) )
      do i = 1, nCells
         if (global_cell_owner(i) == iam+1) then
            j = global_to_local_cell(i)
            gcol(j) = i
            llon(j) = lonCell(i)
            llat(j) = latCell(i)
            larea(j) = areaCell(i)/(MPAS_SPHERE_RAD**2.0_r8)
         end if
      end do

      ! Get reference pressures from the dynamical core.
      call dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)
      call ref_pres_init (pref_edge, pref_mid, num_pr_lev)

      call blocks_init

      call define_cam_grids
   
   end subroutine dyn_grid_init
   
   
   subroutine get_block_ldof_d(nlev, ldof)
   
       use pio, only : pio_offset_kind
   
       integer, intent(in) ::  nlev
       integer(kind=pio_offset_kind), intent(out) :: ldof(:)
   
   end subroutine get_block_ldof_d
   
   
   subroutine blocks_init
   
      ! Determine number of blocks and column indices assigned to each block
      ! 
      ! Method: Assume uniform number of blocks per process
      !         Assume blocks are assigned to processes in order
      !         Assume that the columns assigned to a given process are assigned to
      !            the blocks in order; that is, the first group of columns is assigned
      !            to the first block, etc.
      !         Note: 'local' sometimes refers to process, and other times to block
   
   
      integer :: pcols, nblck
      integer :: n, ncpb_quot, ncpb_rem
      integer :: nc, i, iloc, ilocb, found, nb_loc, nblk
   
      integer :: nbeg(npes), iofst(npes)
   
      character(len=*), parameter :: subname = "dyn_grid::blocks_init"

      ! First determine number of blocks
   
      pcols = PCOLS
      nblck = (real(nCells) / real(pcols)) + 0.5_r8         ! Total number of blocks (approx)
      nblocks_per_pe = (real(nblck) / real(npes)) + 0.5_r8  ! Number of blocks per process (round)
      nblocks_per_pe = max(nblocks_per_pe, 1)               ! Guarantee at least one block per process
   
      !  nblocks_per_pe = 1                                 ! Uncomment to force one block per process
   
      nblocks_tot = npes * nblocks_per_pe                   ! Total number of blocks
   
      ! Assign cells (columns) to blocks
   
      !    nCells_per_proc(1:npes) is the number of cells assigned to each process
      !    global_cell_owner(1:nCells) is the process that owns a given cell (1-based)
      !    global_to_local_cell(1:nCells) gives the local (relative to process) cell index
      !      corresponding to the global cell index
   
      ! Compute number of columns in each block
      allocate( num_col_per_block(nblocks_tot) )
      max_col_per_block = 0
      do n = 1, npes
   
         ncpb_quot = nCells_per_proc(n) / nblocks_per_pe
         ncpb_rem = nCells_per_proc(n) - ncpb_quot * nblocks_per_pe
   
         do nc = 1, nblocks_per_pe
            nblk = (n-1) * nblocks_per_pe + nc
            num_col_per_block(nblk) = ncpb_quot
            if (nc <= ncpb_rem) num_col_per_block(nblk) = num_col_per_block(nblk) + 1
         enddo
   
         nblk = (n-1) * nblocks_per_pe + 1
         max_col_per_block = max(num_col_per_block(nblk), max_col_per_block)
   
      enddo
   
      allocate( global_blockid(nCells) )
      allocate( local_col_index(nCells) )  !  local is relative to block
      allocate( col_indices_in_block(max_col_per_block, nblocks_tot) )
      col_indices_in_block(:,:) = -1
   
      ! Perform column assignment
      nbeg(:) = 1
      iofst(:) = 0
      do i = 1, nCells
   
         n = global_cell_owner(i)
         iloc = global_to_local_cell(i)      ! local (to process) column index
   
         ! Search for owning block (among those blocks assigned to the global_cell_owner)
         ilocb = iloc - iofst(n)
         found = 0
         do nc = nbeg(n), nblocks_per_pe
            nblk = (n-1) * nblocks_per_pe + nc               ! global block index
   
            !  See if local (to process) column index falls within current block
            if (ilocb <= num_col_per_block(nblk)) then

               !  Found block; make sure column index is within extent of array; then exit loop
               if (ilocb > max_col_per_block) then
                  call endrun( subname // ':: column to block mapping failed 1' )
               endif

               col_indices_in_block(ilocb,nblk) = i  !  global column index
               found = nc
               global_blockid(i) = nblk
               local_col_index(i) = ilocb  !  local to block
               exit
            endif
   
            ! Offset column index by number of columns in just tested block, and test next block
            ilocb = ilocb - num_col_per_block(nblk)
         enddo
   
         if (found == 0) then
            call endrun( subname // ':: column to block mapping failed 2' )
         endif
      enddo
      nbeg(n) = found
      iofst(n) = iloc - ilocb
   
      !  Check that proper indices of col_indices_in_block have been set
      do nblk = 1, nblocks_tot
         do ilocb = 1, num_col_per_block(nblk)
            if (col_indices_in_block(ilocb,nblk) < 1) then
               call endrun( subname // ':: column to block mapping failed 3' )
            endif
         enddo
         do ilocb = num_col_per_block(nblk)+1, max_col_per_block 
            if (col_indices_in_block(ilocb,nblk) .ne. -1) then
               call endrun( subname // ':: column to block mapping failed 4' )
            endif
         enddo
      enddo
   
   end subroutine blocks_init
   
   
   subroutine get_block_bounds_d( block_first, block_last )
   
      ! Return first and last indices used in global block ordering
   
      integer, intent(out) :: block_first  ! first (global) index used for blocks
      integer, intent(out) :: block_last   ! last (global) index used for blocks
   
      block_first = 1
      block_last = nblocks_tot
   
   end subroutine get_block_bounds_d
   
   
   subroutine get_block_gcol_d( blockid, size, cdex )
   
      ! Return list of dynamics column indices in given block
      ! 
      ! Method: Assuming uniform number of blocks per process, also that
      !         blockid is one-based.
   
      integer, intent(in) :: blockid      ! global block id
      integer, intent(in) :: size         ! array size
   
      integer, intent(out):: cdex(size)   ! global column indices
   
      integer :: count
      
      do count = 1, num_col_per_block(blockid)
            cdex(count) = col_indices_in_block(count, blockid)
      enddo
   
   end subroutine get_block_gcol_d
   
   
   integer function get_block_gcol_cnt_d( blockid )
   
      ! Return number of dynamics columns in indicated block
   
      integer, intent(in) :: blockid
   
      get_block_gcol_cnt_d = num_col_per_block(blockid)
   
   end function get_block_gcol_cnt_d
   
   
   integer function get_block_lvl_cnt_d( blockid, bcid )
   
      ! Return number of levels in indicated column. If column
      ! includes surface fields, then it is defined to also
      ! include level 0.
   
      integer, intent(in) :: blockid  ! global block id
      integer, intent(in) :: bcid    ! column index within block
   
      get_block_lvl_cnt_d = plevp
   
   end function get_block_lvl_cnt_d
   
   
   subroutine get_block_levels_d( blockid, bcid, lvlsiz, levels )
   
      ! Return level indices in indicated column. If column
      ! includes surface fields, then it is defined to also
      ! include level 0.
   
      integer, intent(in) :: blockid  ! global block id
      integer, intent(in) :: bcid    ! column index within block
      integer, intent(in) :: lvlsiz   ! dimension of levels array
   
      integer, intent(out) :: levels(lvlsiz) ! levels indices for block
   
      integer :: k
   
      character(len=*), parameter :: subname = "dyn_grid::get_block_levels_d"
   
      !  latitude slice block
      if ( lvlsiz < plev + 1 ) then
         call endrun( subname // ':: level arrays not large enough' )
      else
         do k = 0, plev
            levels(k+1) = k
         enddo
         do k = plev+2, lvlsiz
            levels(k) = -1
         enddo
      endif
   
   end subroutine get_block_levels_d
   
   
   subroutine get_gcol_block_d( gcol, cnt, blockid, bcid, localblockid )
   
      ! Return global block index and local column index
      ! for global column index
   
      integer, intent(in) :: gcol     ! global column index
      integer, intent(in) :: cnt      ! size of blockid and bcid arrays
   
      integer, intent(out) :: blockid(cnt) ! block index
      integer, intent(out) :: bcid(cnt)    ! column index within block
      integer, intent(out), optional :: localblockid(cnt)
   
      integer :: j
   
      character(len=*), parameter :: subname = "dyn_grid::get_gcol_block_d"
   
      if ( cnt < 1 ) then
         call endrun( subname // ':: arrays not large enough' )
      endif
   
      blockid(1) = global_blockid(gcol)
      bcid(1) = local_col_index(gcol)
   
      do j=2,cnt
         blockid(j) = -1
         bcid(j)    = -1
      enddo
   
   end subroutine get_gcol_block_d
   
   
   integer function get_gcol_block_cnt_d( gcol )
   
      ! Return number of blocks containing data for the vertical column
      ! with the given global column index
   
      integer, intent(in) :: gcol     ! global column index
   
      get_gcol_block_cnt_d = 1
   
   end function get_gcol_block_cnt_d
   
   
   integer function get_block_owner_d( blockid )
   
      ! Return id of process that "owns" the indicated block
   
      integer, intent(in) :: blockid  ! global block id
   
      get_block_owner_d = (blockid - 1) / nblocks_per_pe
   
   end function get_block_owner_d
   
   
   subroutine get_horiz_grid_dim_d( hdim1_d, hdim2_d )
   
      ! Return declared horizontal dimensions of computational grid.
      ! For non-lon/lat grids, declare grid to be one-dimensional,
      ! i.e., (ncols x 1)
   
      integer, intent(out) :: hdim1_d             ! first horizontal dimension
      integer, intent(out), optional :: hdim2_d   ! second horizontal dimension
   
      hdim1_d = nCells
   
      if( present(hdim2_d) ) hdim2_d = 1
   
   end subroutine get_horiz_grid_dim_d
   
   
   subroutine get_horiz_grid_d( nxy, clat_d_out, clon_d_out, area_d_out, &
       wght_d_out, lat_d_out, lon_d_out )
   
      ! Return latitude and longitude (in radians), column surface
      ! area (in radians squared) and surface integration weights
      ! for global column indices that will be passed to/from physics
   
      use pmgrid, only : MPAS_SPHERE_RAD
      
      integer, intent(in) :: nxy                     ! array sizes
   
      real(r8), intent(out), optional :: clat_d_out(:) ! column latitudes
      real(r8), intent(out), optional :: clon_d_out(:) ! column longitudes
      real(r8), intent(out), target, optional :: area_d_out(:)
      real(r8), intent(out), target, optional :: wght_d_out(:)
      !  weight
      real(r8), intent(out), optional :: lat_d_out(:)  ! column degree latitudes
      real(r8), intent(out), optional :: lon_d_out(:)  ! column degree longitudes

      character(len=*), parameter :: subname = "dyn_grid::get_horiz_grid_d"

      if ( nxy /= nCells ) then
         call endrun( subname // ':: incorrect number of cells' )
      endif
   
      ! this is initialization of old commap vars that doesn't belong here
      if ( .not. associated(clon) ) then
         allocate( clon(nxy,1) )
         allocate( londeg(nxy,1) )
         allocate( latdeg(nxy) )
         allocate( w(nxy) )
      end if
      clon(:,1) = lonCell(:)
      londeg    = clon * rad2deg
      latdeg    = latcell * rad2deg
      w         = areacell
   
      if ( present( clat_d_out ) ) then
         clat_d_out(:) = latCell(:)
      endif
   
      if ( present( clon_d_out ) ) then
         clon_d_out(:) = lonCell(:)
      endif
   
      if ( present( area_d_out ) ) then
         area_d_out(:) = areaCell(:) / (MPAS_SPHERE_RAD**2.0_r8)
      endif
   
      if ( present( wght_d_out ) ) then
         wght_d_out(:) = area_d_out(:)
      endif

      if ( present( lat_d_out ) ) then
         lat_d_out(:) = latCell(:) * rad2deg
      endif

      if ( present( lon_d_out ) ) then
         lon_d_out(:) = lonCell(:) * rad2deg
      endif

   end subroutine get_horiz_grid_d
   
   
!   subroutine get_gcol_lat( gcid, lat )
!   
!      ! Return latitudes for global column indices
!   
!      integer, intent(in) :: gcid(:)
!      real(r8), intent(out) :: lat(:)
!      integer :: j, glen
!   
!      glen = size(gcid)
!   
!      do j = 1, glen
!         lat(j) = latCell(gcid(j))
!      end do
!   
!   end subroutine get_gcol_lat
!   
!   
!   subroutine get_gcol_lon( gcid, lon )
!   
!      ! Return longitudes for global column indices
!      
!      integer, intent(in) :: gcid(:)
!      real(r8), intent(out) :: lon(:)
!      integer :: i, glen
!   
!      glen = size(gcid)
!   
!      do i = 1, glen
!         lon(i) = lonCell(gcid(i))
!      end do
!   
!   end subroutine get_gcol_lon
   
   
   function get_dyn_grid_parm_real2d(name) result(rval)
   
      character(len=*), intent(in) :: name
      real(r8), pointer :: rval(:,:)
   
      if(name == 'clon') then
         rval => clon
      else if(name == 'londeg') then
         rval => londeg
      else if(name == 'bounds_lon') then
         rval => bounds_lon
      else if(name == 'bounds_lat') then
         rval => bounds_lat
      else
         nullify(rval)
      end if
   
   end function get_dyn_grid_parm_real2d
   
   
   function get_dyn_grid_parm_real1d(name) result(rval)
   
      character(len=*), intent(in) :: name
      real(r8), pointer :: rval(:)
   
      if(name == 'clat') then
         rval => latcell
      else if(name == 'latdeg') then
         rval => latdeg
      else if(name == 'w') then
         rval => w
      else
         nullify(rval)
      end if
   
   end function get_dyn_grid_parm_real1d
   
   
   integer function get_dyn_grid_parm(name) result(ival)
   
      character(len=*), intent(in) :: name
   
      if (name == 'plon') then
         ival = nCells
      else if(name == 'plev') then
         ival = plev
      else if(name == 'plevp') then
         ival = plevp
      else if(name == 'maxEdges') then
         ival = maxEdges
      else	
         ival = -1
      end if
   
   
   end function get_dyn_grid_parm

   
   subroutine dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)
   
      ! return reference pressures for the dynamics grid
   
      use pmgrid, only: plev
      use mpas_cam_interface, only: mpas_get_pref_profile
   
      ! arguments
      real(r8), intent(out) :: pref_edge(:) ! reference pressure at layer edges (Pa)
      real(r8), intent(out) :: pref_mid(:)  ! reference pressure at layer midpoints (Pa)
      integer,  intent(out) :: num_pr_lev   ! number of top levels using pure pressure representation
   
      integer :: k
      !-----------------------------------------------------------------------
   
      call mpas_get_pref_profile(pref_edge, pref_mid)
   
   !   write(0,*) 'DEBUG REF_PRES:'
   !   do k=1,plev
   !      write(0,*) k, pref_edge(k), pref_mid(k)
   !   end do
   !   write(0,*) plev+1, pref_edge(plev+1)
   
      num_pr_lev = 0
   
   end subroutine dyn_grid_get_pref
   
   
   subroutine dyn_grid_find_gcol( lat, lon, owner, col, lbk, rlat, rlon ) 
   
      !-------------------------------------------------------------------------------
      ! Return the lat/lon information (and corresponding MPI task number (owner)) 
      ! of the global model grid column nearest to the input lat, lon coordinate.
      !-------------------------------------------------------------------------------
   
      real(r8), intent(in) :: lat
      real(r8), intent(in) :: lon
      integer, intent(out) :: owner
      integer, intent(out) :: col
      integer, intent(out) :: lbk
   
      real(r8),optional, intent(out) :: rlon
      real(r8),optional, intent(out) :: rlat
   
      call endrun('dyn_grid_find_gcol not supported for mpas dycore')
   
   end subroutine dyn_grid_find_gcol

   
   subroutine dyn_grid_find_gcols( lat, lon, nclosest, owners, col, lbk, rlat, rlon, idyn_dists )
   
      !-------------------------------------------------------------------------------
      ! This returns the lat/lon information (and corresponding MPI task numbers
      ! (owners))
      ! of the global model grid columns nearest to the input satellite coordinate
      ! (lat,lon)
      !-------------------------------------------------------------------------------
   
      real(r8), intent(in) :: lat
      real(r8), intent(in) :: lon
      integer, intent(in)  :: nclosest
      integer, intent(out) :: owners(nclosest)
      integer, intent(out) :: col(nclosest)
      integer, intent(out) :: lbk(nclosest)
   
      real(r8),optional, intent(out) :: rlon(nclosest)
      real(r8),optional, intent(out) :: rlat(nclosest)
      real(r8),optional, intent(out) :: idyn_dists(nclosest)
   
      call endrun('dyn_grid_find_gcol not supported for mpas dycore')
   
   end subroutine dyn_grid_find_gcols


   subroutine dyn_grid_get_colndx( igcol, nclosest, owners, col, lbk )
     use spmd_utils, only: iam
     use pmgrid,     only: plon
   
     integer, intent(in)  :: nclosest
     integer, intent(in)  :: igcol(nclosest)
     integer, intent(out) :: owners(nclosest)
     integer, intent(out) :: col(nclosest)
     integer, intent(out) :: lbk(nclosest)
   
     integer  :: i
     integer :: blockid(1), bcid(1), lclblockid(1)
   
!     do i = 1,nclosest
!  
!       call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
!       owners(i) = get_block_owner_d(blockid(1))
!  
!       if ( iam==owners(i) ) then
!          lbk(i) = lclblockid(1)
!          col(i) = bcid(1)
!       else
!          lbk(i) = -1
!          col(i) = -1
!       endif
!  
!    end do

      call endrun('dyn_grid_find_gcol not supported for mpas dycore')

   end subroutine dyn_grid_get_colndx

   subroutine dyn_grid_get_elem_coords( icol, rlon, rlat, cdex )
   
      ! Return coordinates of a specified block element of the dyn grid
   
     integer, intent(in) :: icol ! block column index
   
     real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the block
     real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the block
     integer, optional, intent(out) :: cdex(:) ! global column index
   
     call endrun('dyn_grid_get_elem_coords not supported for mpas dycore')
   
   end subroutine dyn_grid_get_elem_coords


   subroutine define_cam_grids()
     use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
     use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
     use pmgrid, only : MPAS_SPHERE_RAD
 
     type(horiz_coord_t), pointer :: lat_coord
     type(horiz_coord_t), pointer :: lon_coord
     integer(iMap),       pointer :: coord_map(:)
     integer(iMap),       pointer :: grid_map(:,:)
     integer                :: i, j   ! Loop variables
 
     ! Map for dynamics grid
     allocate(grid_map(3, numcols))
     do i = 1, numcols
        grid_map(1, i) = i
        grid_map(2, i) = gcol(i) 
        grid_map(3, i) = gcol(i)
     end do

     ! Note: not using get_horiz_grid_dim_d or get_horiz_grid_d since those
     !       are deprecated ('cause I said so' -- goldy)
     ! Coordinates for MPAS points
     allocate(coord_map(numcols))
     coord_map = gcol
     lat_coord => horiz_coord_create('lat', 'ncol', nCells, 'latitude',      &
          'degrees_north', 1, size(llat), llat*rad2deg, map=coord_map)
     lon_coord => horiz_coord_create('lon', 'ncol', nCells, 'longitude',     &
          'degrees_east', 1, size(llon), llon*rad2deg, map=coord_map)
 
     ! MPAS grid
     call cam_grid_register('mpas_cell', dyn_decomp, lat_coord, lon_coord,     &
          grid_map, block_indexed=.false., unstruct=.true.)
     call cam_grid_attribute_register('mpas_cell', 'area', 'MPAS cell areas',  &
          'ncol', larea, coord_map)
     call cam_grid_attribute_register('mpas_cell', 'nCells', '', nCells)
     nullify(grid_map) ! Map belongs to grid now
 
   end subroutine define_cam_grids


   subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
     use cam_grid_support, only: max_hcoordname_len
 
     ! Dummy arguments
     character(len=max_hcoordname_len),          intent(out) :: gridname
     character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)
 
     gridname = 'mpas_cell'
     allocate(grid_attribute_names(2))
     grid_attribute_names(1) = 'area'
     grid_attribute_names(2) = 'nCells'
 
   end subroutine physgrid_copy_attributes_d


end module dyn_grid
