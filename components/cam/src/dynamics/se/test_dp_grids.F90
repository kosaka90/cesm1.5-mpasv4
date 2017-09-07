module test_dp_grids
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use GridGraph_mod,          only: num_neighbors
  use pio,                    only: io_desc_t, file_desc_t, var_desc_t

  implicit none
  private
  save

#ifdef TEST_GRID_COUPLING
  public dptest_print_dyn_elem
  public dptest_print_phys_elem
  public dptest_init
  public dptest_set_dyn_state_tvals
  public dptest_test_phys_state_tvals
  public dptest_test_dyn_state_ttends
  public dptest_set_phys_state_ttends
  public dptest_write_dyn_state_history
  public dptest_write_phys_state_history
  public dptest_write_dyn_tend_history
  public dptest_write_phys_tend_history
  public dptest_compute_error
  public dptest_print_array_stats
#if 0
  !! Obsolete?
  public dptest_init_coords
  public dptest_bndry_exchange
  public dptest_init_grid_file
  public dptest_write_grid_file
  public dptest_close_grid_file
#endif

  interface dptest_compute_error
    module procedure compute_error_1d
    module procedure compute_error_2d
  end interface

  interface dptest_print_array_stats
    module procedure print_array_stats_phys
    module procedure print_array_stats_dyn
  end interface

  interface dptest_write_phys_state_history
    module procedure dptest_write_phys_state_history_all
    module procedure dptest_write_phys_state_history_chunk
  end interface

  interface dptest_write_phys_tend_history
    module procedure dptest_write_phys_tend_history_chunk
    module procedure dptest_write_phys_tend_history_all
    module procedure dptest_write_phys_ptend_history_chunk
    module procedure dptest_write_phys_ptend_history_all
    module procedure dptest_write_tend_fields_history_chunk
  end interface

  interface fatal_error
    module procedure fatal_errori
  end interface

#if 0
  !! Obsolete?
  interface dptest_write_grid_file
     module procedure write_grid_file_2d, write_grid_file_3d
  end interface

  integer                         :: np_dim, nphys2_dim      ! diagfile dims
  integer                         :: nelem_dim               ! diagfile dims
  integer                         :: time_dim                ! diagfile time
  type(io_desc_t)                 :: phys_iodesc, dyn_iodesc ! PIO decomps
  type(file_desc_t)               :: diagfile    ! PIO file descriptor for diags
  type(var_desc_t), target        :: dyn_grid_id             ! diagfile vars
  type(var_desc_t), target        :: phys_grid_id            ! diagfile var
  type(var_desc_t), target        :: dyn_grid_error_id       ! diagfile var
#endif
  integer                         :: begchunk, endchunk, pver, pcols
  real(r8), allocatable           :: phys_lons(:)      ! (pcols*chunks)
  real(r8), allocatable           :: phys_lats(:)      ! (pcols*chunks)
  integer,  allocatable           :: ncols(:)          ! (chunks)
  integer                         :: total_cols        ! sum(ncols)
  real(r8), allocatable           :: dyn_lons(:)       ! (np*np*nelemd)
  real(r8), allocatable           :: dyn_lats(:)       ! (np*np*nelemd)
  logical                         :: coords_initialized = .false.
  integer,  public                :: test_cells(num_neighbors + 1) = 0

  real(r8), parameter             :: min_scale = 1.0E-18_r8  ! Min error scale
  real(r8), parameter             :: err_lev = 2.0E-2_r8     ! Error report lev
  real(r8), parameter             :: state_factor = 100.0_r8
  real(r8), parameter             :: tend_factor = 1.0e-9_r8
#endif

!!
contains
!!

#ifdef TEST_GRID_COUPLING
  subroutine fatal_errori(errnum, msg)
    use cam_abortutils,      only: endrun


    integer,                 intent(in)  :: errnum
    character(len=*),        intent(in)  :: msg

    character(len=80)                          :: errmsg

    write(errmsg, '(a,i6," ",a)') ' ERROR ', errnum, msg
    call endrun(errmsg)

  end subroutine fatal_errori

  ! A simple analytic function to test grid interpolation
  ! Test functions
  ! (0) fout = 2
  ! (1) fout = 2 + theta
  ! (2) fout = 2 + cos(theta)
  ! (3) fout = 2 + cos(lambda)
  ! (4) fout = 2 + cos^2(theta)*cos(8*lambda)
  ! (5) fout = 2 + cos^2(theta)*cos(16*lambda)
  ! (6) fout = 2 + cos^2(theta)*cos(32*lambda)
  ! All functions multiplied by factor (default 1.0)
  function test_func(lon, lat, funcnum, factor) result(fout)
    use spmd_utils,      only: masterproc
    use cam_logfile,     only: iulog
    use shr_sys_mod,     only: shr_sys_flush
    use spmd_utils,      only: masterproc, npes, iam, mpicom
    use physical_constants, only : DD_PI

    real(r8),           intent(in)  :: lon
    real(r8),           intent(in)  :: lat
    integer,            intent(in)  :: funcnum
    real(r8), optional, intent(in)  :: factor
    real(r8)                        :: fout
    real(r8)                        :: lon1,lat1,R0,Rg1,Rg2,lon2,lat2

    select case(funcnum)
    case(0)
      fout = 2.0_r8
    case(1)
      fout = 2.0_r8 + lat
    case(2)
      fout = 2.0_r8 + cos(lat)
    case(3)
      fout = 2.0_r8 + cos(lon)
    case(4, 5)
      fout = 2.0_r8 + (cos(lat)**2 * cos(ISHFT(1, funcnum-1) * lon))
    case(6)       !   Non-smooth scalar field (slotted cylinder)
       R0=0.5_r8
       lon1=4.0_r8*DD_PI/5.0_r8
       lat1=0.0_r8
       Rg1 = acos(sin(lat1)*sin(lat)+cos(lat1)*cos(lat)*cos(lon-lon1))
       lon2=6.0_r8*DD_PI/5.0_r8
       lat2=0.0_r8
       Rg2 = acos(sin(lat2)*sin(lat)+cos(lat2)*cos(lat)*cos(lon-lon2))

       if ((Rg1 .le. R0) .AND. (abs(lon-lon1).ge. R0/6)) then
          fout = 2.0_r8
       elseif ((Rg2 .le. R0) .AND. (abs(lon-lon2).ge. R0/6)) then
          fout = 2.0_r8
       elseif ((Rg1 .le. R0) .AND. (abs(lon-lon1) < R0/6) &
            .AND. (lat-lat1 < -5.0_r8*R0/12.0_r8)) then
          fout = 2.0_r8
       elseif ((Rg2 .le. R0) .AND. (abs(lon-lon2) < R0/6) &
            .AND. (lat-lat2 > 5.0_r8*R0/12.0_r8)) then
          fout = 2.0_r8
       else
          fout = 1.0_r8
       endif
!       fout=1.0_r8
    case default
       call fatal_error(funcnum, " is an illegal funcnum_arg")
    end select
    if (present(factor)) then
       fout = fout * factor
    end if
  end function test_func

  function dptest_test_func_string(funcnum) result(func_string)
    use cam_abortutils, only: endrun
    use spmd_utils,     only: masterproc
    use cam_logfile,    only: iulog

    integer,            intent(in)  :: funcnum
    character(len=40)               :: func_string

    select case(funcnum)
    case (0)
      func_string = '2'
    case (1)
      func_string = '2 + theta'
    case(2)
      func_string = '2 + cos(theta)'
    case(3)
      func_string = '2 + cos(lambda)'
    case (4)
      func_string = '2 + cos^2(theta)*cos(8*lambda)'
    case (5)
      func_string = '2 + cos^2(theta)*cos(16*lambda)'
    case (6)
      func_string = '2 + cos^2(theta)*cos(32*lambda)'
    case default
      if (masterproc) then
        write(iulog, '(a, i4)') "ERROR: Illegal funcnum_arg, ", funcnum
      end if
      call endrun('ERROR: Illegal funcnum_arg')
    end select
  end function dptest_test_func_string

  subroutine set_test_func(ncol, clon, clat, array, factor, funcnum_arg)

    integer,                 intent(in)    :: ncol
    real(r8),                intent(in)    :: clon(ncol)
    real(r8),                intent(in)    :: clat(ncol)
    real(r8),                intent(inout) :: array(ncol)
    real(r8),      optional, intent(in)    :: factor
    integer,       optional, intent(in)    :: funcnum_arg

    integer                                :: i
    integer                                :: funcnum

    if (present(funcnum_arg)) then
      funcnum = funcnum_arg
    else
      funcnum = 1
    end if


    do i = 1, ncol
      array(i) = test_func(clon(i), clat(i), funcnum, factor)
    end do

  end subroutine set_test_func

  subroutine test_func_error(ncol, clon, clat, testval, factor, fnum_arg)
    use cam_logfile,   only: iulog
    use spmd_utils,    only: iam, mpicom, npes, masterproc, masterprocid
    use spmd_utils,    only: mpi_integer, mpi_real8, mpi_sum, mpi_max
    use shr_const_mod, only: SHR_CONST_PI

    integer,                 intent(in)    :: ncol
    real(r8),                intent(in)    :: clon(ncol)
    real(r8),                intent(in)    :: clat(ncol)
    real(r8),                intent(in)    :: testval(ncol)
    real(r8),      optional, intent(in)    :: factor
    integer,       optional, intent(in)    :: fnum_arg

    integer                                :: n, i, ierr
    integer                                :: fnum
    integer                                :: numerrs
    real(r8)                               :: maxerr(3)    ! (err, lon, lat)
    real(r8)                               :: lonerr, laterr
    real(r8), allocatable                  :: checkval(:)
    real(r8), allocatable                  :: error(:)
    logical                                :: header_printed
    real(r8), parameter                    :: rad2deg = 180.0_r8 / SHR_CONST_PI

    if (present(fnum_arg)) then
      fnum = fnum_arg
    else
      fnum = 1
    end if

    header_printed = .false.

    allocate(checkval(ncol))
    allocate(error(ncol))
    call set_test_func(ncol, clon, clat, checkval, factor, fnum)
    call dptest_compute_error(testval, checkval, error)
    maxerr = -1.0_r8
    numerrs = 0
    do n = 0, npes - 1
      if (iam == n) then
        do i = 1, ncol
          if (error(i) > maxerr(1)) then
            maxerr(1) = error(i)
            maxerr(2) = clon(i)
            maxerr(3) = clat(i)
          end if
          if (error(i) > err_lev) then
            if (.not. header_printed) then
              write(iulog, *) 'Test function: ', dptest_test_func_string(fnum)
              write(iulog, *) 'PE     pt   Lon      Lat     Error       Value      Expected'
              header_printed = .true.
            end if
            numerrs = numerrs + 1
402         format(i3,'  ',i5,'  ',f6.2,'  ',f6.2,3('  ',es10.3))
            if (numerrs < 50) then
              lonerr = clon(i) * rad2deg
              laterr = clat(i) * rad2deg
              write(iulog, 402) n, i, lonerr, laterr, error(i), testval(i), checkval(i)
            end if
          end if
        end do
      end if
      call mpi_barrier(mpicom, ierr)
    end do
    i = numerrs
    call MPI_Reduce(i, numerrs, 1, mpi_integer, mpi_sum, masterprocid, mpicom, ierr)
    checkval(1:3) = maxerr(1:3)
    call MPI_Reduce(checkval, maxerr, 3, mpi_real8, mpi_max, masterprocid, mpicom, ierr)
    if (masterproc) then
      if (numerrs > 0) then
        lonerr = maxerr(2) * rad2deg
        laterr = maxerr(3) * rad2deg
        write(iulog, '(i5,a,es12.5," at (",f6.2,", ",f6.2,")")') numerrs,     &
             ' errors, max error = ', maxerr(1), maxerr(2), maxerr(3)
      end if
    end if
    deallocate(checkval)
    deallocate(error)

  end subroutine test_func_error

  function find_element(elem, number) result(elem_out)
    use element_mod,     only: element_t
    use dimensions_mod,  only: nelemd

    type(element_t), dimension(:), intent(in)  :: elem
    integer,                       intent(in)  :: number
    integer                                    :: elem_out
    integer                                    :: i

    elem_out = -1
    do i = 1, nelemd
      if ( elem(i)%vertex%number == number) then
        elem_out = i
        exit
      end if
    end do

  end function find_element

  subroutine dptest_print_dyn_elem(elem, array, ie, message)
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog
    use element_mod,    only: element_t
    use dimensions_mod, only: np
    use control_mod,    only: north,south,east,west,neast,nwest,seast,swest

    type(element_t),               intent(in)  :: elem(:)
    real(r8),                      intent(in)  :: array(:, :, :)
    integer,                       intent(in)  :: ie
    character(len=*), optional,    intent(in)  :: message

    character(len=32)            :: blkstr
    character(len=128)           :: fmtstr               ! Format string
    integer                      :: row, col, i, j, k, n
    real(r8)                     :: area
    real(r8)                     :: numbers(np*3)
    integer                      :: elementNums(3,3)
    integer                      :: nbr_ind(9) = (/swest, south, seast,       &
         west, 0, east, nwest, north, neast/)

    ! First, write a line for one element
    write(blkstr, '(i1,"(f8.4,'', ''),f8.4")') (np - 1)
    ! Expand it to three element lines
    write(fmtstr, '("(",a,",''      '',",a,",''      '',",a,")")')            &
         trim(blkstr), trim(blkstr), trim(blkstr)
    if (present(message)) then
      write(iulog, *) trim(message)
    else
      write(iulog,'("gtmp_dyn = ")')
    end if
    ! Gather up the numbers
    do row = 1, 3
      do i = 1, np
        do col = 1, 3
          j = ((row - 1) * 3) + col
          if (nbr_ind(j) > 0) then
            n = find_element(elem, elem(ie)%vertex%nbrs(nbr_ind(j)))
            if (n < 0) then
              write(fmtstr, *) 'unable to find neighbor element,', j,         &
                   ', to element, ',ie
              call endrun(fmtstr)
            endif
          else
            n = ie
          endif
          elementNums(col, row) = n
          do k = 1, np
            numbers(((col - 1) * np) + k) = array(k, i, n)
          end do
        end do
        write(iulog, fmtstr) (numbers(k), k = 1, (np * 3))
      end do
      write(iulog, *) ' '
    end do
!    area = dptest_get_dyn_area(array, np, ie)
!    write(iulog, '("dyn_area          = ",f10.6)') area
    write(iulog, '(3i6)') ((elementNums(j, k), j = 1, 3), k = 1, 3)
    write(iulog, *) ' '
    if (present(message)) then
      write(iulog, *) trim(message)
    end if
  end subroutine dptest_print_dyn_elem

  subroutine dptest_print_phys_elem(elem, array, ie, nphys, message)
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog
    use element_mod,    only: element_t
    use control_mod,    only: north,south,east,west,neast,nwest,seast,swest

    type(element_t), dimension(:), intent(in)  :: elem
    real(r8),                      intent(in)  :: array(:, :)
    integer,                       intent(in)  :: ie
    integer,                       intent(in)  :: nphys
    character(len=*), optional,    intent(in)  :: message

    character(len=32)            :: blkstr
    character(len=128)           :: fmtstr               ! Format string
    integer                      :: row, col, i, j, k, n
    real(r8)                     :: area
    real(r8)                     :: numbers(nphys*3)
    integer                      :: nbr_ind(9) = (/swest, south, seast,       &
         west, 0, east, nwest, north, neast/)

    ! First, write a line for one element
    write(blkstr, '(i1,"(f8.4,'', ''),f8.4")') (nphys - 1)
    ! Expand it to three element lines
    write(fmtstr, '("(",a,",''      '',",a,",''      '',",a,")")')            &
         trim(blkstr), trim(blkstr), trim(blkstr)
    if (present(message)) then
      write(iulog, *) trim(message)
    else
      write(iulog,'("gtmp_phys = ")')
    end if
    ! Gather up the numbers
    do row = 1, 3
      do i = 1, nphys
        do col = 1, 3
          j = ((row - 1) * 3) + col
          if (nbr_ind(j) > 0) then
            n = find_element(elem, elem(ie)%vertex%nbrs(nbr_ind(j)))
            if (n < 0) then
              call endrun('unable to find neighbor element')
            endif
          else
            n = ie
          endif
          do k = 1, nphys
            numbers(((col - 1) * nphys) + k) = array(k + ((i - 1) * nphys), n)
          end do
        end do
        write(iulog, fmtstr) (numbers(k), k = 1, (nphys * 3))
      end do
      write(iulog, *) ' '
    end do
!    area = dptest_get_phys_area(array, nphys, ie)
!    write(iulog, '("phys_area          = ",f10.6)') area
    write(iulog, *) ' '
    if (present(message)) then
      write(iulog, *) trim(message)
    end if
  end subroutine dptest_print_phys_elem

  subroutine dptest_init(elem, pcarg, bchunk, echunk, nver, plats, plons, cols)
    use cam_history,     only:  addfld, horiz_only, add_default
    use element_mod,     only: element_t
    use dimensions_mod,  only: nlev, nelemd, npsq, np
    use cam_logfile,     only: iulog
    use shr_sys_mod,     only: shr_sys_flush
    use spmd_utils,      only: mpi_integer, mpi_sum, iam, mpicom

    type(element_t), intent(in)  :: elem(:)
    integer,         intent(in)  :: pcarg
    integer,         intent(in)  :: bchunk
    integer,         intent(in)  :: echunk
    integer,         intent(in)  :: nver
    real(r8),        intent(in)  :: plats(pcarg, bchunk:echunk)
    real(r8),        intent(in)  :: plons(pcarg, bchunk:echunk)
    real(r8),        intent(in)  :: cols(bchunk:echunk)

    integer                      :: i, j, r, ie
    integer                      :: temp_cells(num_neighbors) = 0

    ! Initialize Physics module variables
    begchunk = bchunk
    endchunk = echunk
    pver = nver
    pcols = pcarg
    allocate(ncols(begchunk:endchunk))
    ncols = cols
    total_cols = sum(ncols)
    allocate(phys_lons(total_cols))
    allocate(phys_lats(total_cols))
    i = 1
    do r = begchunk, endchunk
      j = i + ncols(r) - 1
      phys_lons(i:j) = plons(1:ncols(r), r)
      phys_lats(i:j) = plats(1:ncols(r), r)
      i = j + 1
    end do
    ! Initialize Dynamics module variables
    allocate(dyn_lons(np*np*nelemd))
    allocate(dyn_lats(np*np*nelemd))
    r = 1
    do ie = 1, nelemd
      do j = 1, np
        do i = 1, np
          dyn_lons(r) = elem(ie)%spherep(i,j)%lon
          dyn_lats(r) = elem(ie)%spherep(i,j)%lat
          r = r + 1
        end do
      end do
    end do

    ! Dynamics fields added for debugging interpolation
    call addfld('dUglldbg',(/ 'lev' /),'I','m/s2',                            &
         'Zonal wind tendency (GLL)',        gridname='GLL')
    call add_default('dUglldbg', 3, 'I')
    call addfld('dVglldbg',(/ 'lev' /),'I','m/s2',                            &
         'Meridional wind tendency (GLL)',   gridname='GLL')
    call add_default('dVglldbg', 3, 'I')
    call addfld('dTglldbg',(/ 'lev' /),'I','K/s',                             &
         'Temperature tendency (GLL)',       gridname='GLL')
    call add_default('dTglldbg', 3, 'I')
    call addfld('dUdbg',(/ 'lev' /),'I','m/s2','Zonal wind tendency')
    call add_default('dUdbg   ', 3, 'I')
    call addfld('dVdbg',(/ 'lev' /),'I','m/s2','Meridional wind tendency')
    call add_default('dVdbg   ', 3, 'I')
    call addfld('dTdbg',(/ 'lev' /),'I','K/s','Temperature tendency')
    call add_default('dTdbg   ', 3, 'I')
    call addfld('Ugll_dbg',(/ 'lev' /),'I','m/s',                             &
         'Zonal wind (GLL)',                 gridname='GLL')
    call add_default('Ugll_dbg', 3, 'I')
    call addfld('Vgll_dbg',(/ 'lev' /),'I','m/s',                             &
         'Meridional wind (GLL)',            gridname='GLL')
    call add_default('Vgll_dbg', 3, 'I')
    call addfld('Tgll_dbg',(/ 'lev' /),'I','K',                               &
         'Temperature (GLL)',                gridname='GLL')
    call addfld('WVgll_dbg',(/ 'lev' /),'I','kg/kg',                               &
         'WV (GLL)',                gridname='GLL')
    call addfld('WIgll_dbg',(/ 'lev' /),'I','kg/kg',                               &
         'WI (GLL)',                gridname='GLL')
    call addfld('WLgll_dbg',(/ 'lev' /),'I','kg/kg',                               &
         'WL (GLL)',                gridname='GLL')
    call addfld('TT_UNgll_dbg',(/ 'lev' /),'I','kg/kg',                               &
         'TT_UN (GLL)',                gridname='GLL')
    call addfld('TT_LWgll_dbg',(/ 'lev' /),'I','kg/kg',                               &
         'TT_LW (GLL)',                gridname='GLL')
    call add_default('Tgll_dbg', 3, 'I')
    call addfld('PHISglld',horiz_only,'I','m2/s2',                            &
         'Surface geopotential (GLL)',       gridname='GLL')
    call add_default('PHISglld', 3, 'I')
    call addfld('PSglldbg',horiz_only,'I','Pa',                               &
         'Surface pressure (GLL)',           gridname='GLL')
    call add_default('PSglldbg', 3, 'I')
    call addfld('OMEGAgld',(/ 'lev' /),'I','Pa/s',                            &
         'Vertical velocity (press.) (GLL)', gridname='GLL')
    call add_default('OMEGAgld', 3, 'I')
    call addfld('U_dbg',(/ 'lev' /),'I','m/s','Zonal wind')
    call add_default('U_dbg   ', 3, 'I')
    call addfld('V_dbg',(/ 'lev' /),'I','m/s','Meridional wind')
    call add_default('V_dbg   ', 3, 'I')
    call addfld('T_dbg',(/ 'lev' /),'I','K','Temperature')
    call add_default('T_dbg   ', 3, 'I')
    call addfld('WV_dbg',(/ 'lev' /),'I','kg/kg','WV')
    call addfld('WI_dbg',(/ 'lev' /),'I','kg/kg','WI')
    call addfld('WL_dbg',(/ 'lev' /),'I','kg/kg','WL')
    call addfld('TT_UN_dbg',(/ 'lev' /),'I','kg/kg','TT_UN')
    call addfld('TT_LW_dbg',(/ 'lev' /),'I','kg/kg','TT_LW')
    call addfld('PHISdbg',horiz_only,'I','m2/s2','Surface geopotential')
    call add_default('PHISdbg', 3, 'I')
    call addfld('PSdbg',horiz_only,'A','Pa','Surface pressure')
    call add_default('PSdbg', 3, 'I')
    call addfld('OMEGAdbg',(/ 'lev' /),'I','Pa/s',                            &
         'Vertical velocity (pressure)')
    call add_default('OMEGAdbg', 3, 'I')

    ! First, set center cell
    test_cells(num_neighbors + 1) = 323
    do ie = 1, nelemd
      if (elem(ie)%GlobalId == test_cells(num_neighbors + 1)) then
        write(iulog, *) 'XXgoldyXX ====================================='
        do i = 1, num_neighbors
          j = elem(ie)%vertex%nbrs_ptr(i)
          if (j > 0) then
            write(iulog, *) i, elem(ie)%vertex%nbrs(j)
            temp_cells(i) = elem(ie)%vertex%nbrs(j)
          end if
        end do
        r = iam
        write(iulog, *) 'Found',test_cells(num_neighbors + 1),' on task',r
        call shr_sys_flush(iulog)
        write(iulog, *) '             lat           lon'
        do j = 1, np
          do i = 1, np
            write(iulog, '(2(a,i1),2(a,e12.6))') '(',i,', ',j,') = ', &
                 elem(ie)%spherep(i,j)%lat, '  ', elem(ie)%spherep(i,j)%lon
          end do
        end do
        write(iulog, *) 'XXgoldyXX ====================================='
        call shr_sys_flush(iulog)
      end if
    end do
    call mpi_allreduce(temp_cells, test_cells, num_neighbors, mpi_integer,    &
         mpi_sum, mpicom, j)

  end subroutine dptest_init

  subroutine dptest_set_dyn_state_tvals(elem)
    use element_mod,            only: element_t
    use dimensions_mod,         only: nelemd, np
    use dyn_comp,               only: TimeLevel

    type(element_t), intent(inout) :: elem(:)

    integer                        :: ie, i, j, k, ind
    real(r8), allocatable          :: test_vals(:)
    integer :: tl_f

    tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)

    ind = np*np*nelemd
    allocate(test_vals(ind))
    ! Change the upper levels of T with test function data
    do k = 1, 6
      call set_test_func(ind, dyn_lons, dyn_lats, test_vals, state_factor, k)
      ind = 1
      do ie = 1, nelemd
        do j = 1, np
          do i = 1, np
            elem(ie)%state%T(i,j,k,tl_f) = test_vals(ind)
            ind = ind + 1
          end do
        end do
      end do
    end do
    deallocate(test_vals)

  end subroutine dptest_set_dyn_state_tvals

  subroutine dptest_write_dyn_state_history(elem)
    use constituents,   only: cnst_get_ind
    use cam_history,            only: outfld, hist_fld_active
    use element_mod,            only: element_t
    use dimensions_mod,         only: nlev, nelemd, np, npsq
    use cam_logfile,            only: iulog
    use dyn_comp,               only: TimeLevel
    integer :: ixcldice, ixcldliq,ixtt_un,ixtt_lw

    type(element_t), intent(inout) :: elem(:)

    real(r8)                       :: ftmp(npsq, nlev, 2)
    integer                        :: ie, i, j, k
    integer :: tl_f

    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('TT_UN' , ixtt_un , abort=.false.)
    call cnst_get_ind('TT_LW' , ixtt_lw , abort=.false.)

    tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)


!    if ( hist_fld_active('Ugll_dbg') .or. hist_fld_active('Vgll_dbg') .or.    &
!         hist_fld_active('Tgll_dbg') .or. hist_fld_active('PHISglld') .or.    &
!         hist_fld_active('PSglldbg') .or. hist_fld_active('Qgll_dbg') .or.    &
!         hist_fld_active('CLDICEgll_dbg')) then
      do ie = 1, nelemd
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%V(i,j,1,k,tl_f)
              ftmp(i+(j-1)*np,k,2) = elem(ie)%state%V(i,j,2,k,tl_f)
            end do
          end do
        end do
        call outfld('Ugll_dbg', ftmp(:,:,1), npsq, ie)
        call outfld('Vgll_dbg', ftmp(:,:,2), npsq, ie)
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%T(i,j,k,tl_f)
              ftmp(i+(j-1)*np,k,2) = elem(ie)%derived%omega(i,j,k)
            end do
          end do
        end do
        call outfld('Tgll_dbg', ftmp(:,:,1), npsq, ie)
        call outfld('OMEGAgld', ftmp(:,:,2), npsq, ie)
        do j = 1, np
          do i = 1, np
            ftmp(i+(j-1)*np,1,1) = elem(ie)%state%phis(i,j)
            ftmp(i+(j-1)*np,1,2) = elem(ie)%state%ps(i,j,tl_f)
          end do
        end do
        call outfld('PHISglld', ftmp(:,1,1), npsq, ie)
        call outfld('PSglldbg', ftmp(:,1,2), npsq, ie)
        !    end if

        if (hist_fld_active('TT_LW')) then
          do j = 1, np
            do i = 1, np
              do k = 1, nlev
                ftmp(i+(j-1)*np,k,1) = elem(ie)%state%q(i,j,k,ixtt_lw)
              end do
            end do
          end do
          call outfld('TT_LWgll_dbg', ftmp(:,:,1), npsq, ie)
        end if
      if (hist_fld_active('TT_UN')) then
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%q(i,j,k,ixtt_un)
            end do
          end do
        end do
        call outfld('TT_UNgll_dbg', ftmp(:,:,1), npsq, ie)
      end if
      if (hist_fld_active('WVgll_dbg')) then
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%q(i,j,k,1)
            end do
          end do
        end do
        call outfld('WVgll_dbg', ftmp(:,:,1), npsq, ie)
      end if
      if (hist_fld_active('WIgll_dbg')) then
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%q(i,j,k,ixcldice)
            end do
          end do
        end do
        call outfld('WIgll_dbg', ftmp(:,:,1), npsq, ie)
      end if
      if (hist_fld_active('WLgll_dbg')) then
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%q(i,j,k,ixcldliq)
            end do
          end do
        end do
        call outfld('WLgll_dbg', ftmp(:,:,1), npsq, ie)
      end if
    end do

!
!, 'TT_MD', 'TT_HI', 'TTRMD' , 'TT_UN'

  end subroutine dptest_write_dyn_state_history

  subroutine dptest_test_phys_state_tvals(msg, phys_state)
    use physics_types,          only: physics_state
    use cam_logfile,            only: iulog
    use spmd_utils,             only: masterproc

    character(len=*),    intent(in)  :: msg
    type(physics_state), intent(in)  :: phys_state(begchunk:endchunk)

    integer                          :: lchnk, k, b, e
    real(r8)                         :: test_values(total_cols)

    if (masterproc) then
      write(iulog, *) trim(msg)
    end if
    ! Check the upper levels of T against test function data
    do k = 1, 6
      b = 1
      do lchnk = begchunk, endchunk
        e = b + ncols(lchnk) - 1
        test_values(b:e) = phys_state(lchnk)%t(1:ncols(lchnk),k)
        b = e + 1
      end do
      call test_func_error(total_cols, phys_lons(:), phys_lats(:),  &
           test_values, state_factor, k)
    end do

  end subroutine dptest_test_phys_state_tvals

  subroutine dptest_write_phys_state_history_all(phys_state)
    use physics_types,          only: physics_state
    use cam_history,            only: outfld
    use constituents,   only: cnst_get_ind
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    integer :: ixcldice, ixcldliq,ixtt_un,ixtt_lw
    integer                         :: lchnk

    do lchnk = begchunk, endchunk
      call outfld('U_dbg   ', phys_state(lchnk)%u(:,:), pcols, lchnk)
      call outfld('V_dbg   ', phys_state(lchnk)%v(:,:), pcols, lchnk)
      call outfld('T_dbg   ', phys_state(lchnk)%t(:,:), pcols, lchnk)
      call outfld('PHISdbg ', phys_state(lchnk)%phis(:), pcols, lchnk)
      call outfld('PSdbg   ', phys_state(lchnk)%ps(:), pcols, lchnk)
      call outfld('OMEGAdbg', phys_state(lchnk)%omega(:,:), pcols, lchnk)

      call cnst_get_ind('WV_dbg', ixcldice, abort=.false.)
      call outfld('WV_dbg   ', phys_state(lchnk)%q(:,:,1), pcols, lchnk)
      call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
      call outfld('WI_dbg   ', phys_state(lchnk)%q(:,:,ixcldice), pcols, lchnk)
      call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
      call outfld('WL_dbg   ', phys_state(lchnk)%q(:,:,ixcldliq), pcols, lchnk)
      call cnst_get_ind('TT_UN' , ixtt_un , abort=.false.)
      call outfld('TT_UN_dbg   ', phys_state(lchnk)%q(:,:,ixtt_un), pcols, lchnk)
      call cnst_get_ind('TT_LW' , ixtt_lw , abort=.false.)
      call outfld('TT_LW_dbg   ', phys_state(lchnk)%q(:,:,ixtt_lw), pcols, lchnk)
    end do

  end subroutine dptest_write_phys_state_history_all

  subroutine dptest_write_phys_state_history_chunk(phys_state, lchnk)
    use physics_types,          only: physics_state
    use cam_history,            only: outfld
    use constituents,   only: cnst_get_ind
    type(physics_state), intent(in) :: phys_state
    integer,             intent(in) :: lchnk
    integer :: ixcldice, ixcldliq,ixtt_un,ixtt_lw

    call outfld('U_dbg   ', phys_state%u(:,:), pcols, lchnk)
    call outfld('V_dbg   ', phys_state%v(:,:), pcols, lchnk)
    call outfld('T_dbg   ', phys_state%t(:,:), pcols, lchnk)
    call outfld('PHISdbg ', phys_state%phis(:), pcols, lchnk)
    call outfld('PSdbg   ', phys_state%ps(:), pcols, lchnk)
    call outfld('OMEGAdbg', phys_state%omega(:,:), pcols, lchnk)
    call outfld('WV_dbg   ', phys_state%q(:,:,1), pcols, lchnk)


    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call outfld('WI_dbg   ', phys_state%q(:,:,ixcldice), pcols, lchnk)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call outfld('WL_dbg   ', phys_state%q(:,:,ixcldliq), pcols, lchnk)
    call cnst_get_ind('TT_UN' , ixtt_un , abort=.false.)
    call outfld('TT_UN_dbg   ', phys_state%q(:,:,ixtt_un), pcols, lchnk)
    call cnst_get_ind('TT_LW' , ixtt_lw , abort=.false.)
    call outfld('TT_LW_dbg   ', phys_state%q(:,:,ixtt_lw), pcols, lchnk)

  end subroutine dptest_write_phys_state_history_chunk

  subroutine dptest_test_dyn_state_ttends(msg, elem)
    use element_mod,            only: element_t
    use dimensions_mod,         only: nelemd, np
    use cam_logfile,            only: iulog
    use spmd_utils,             only: masterproc

    character(len=*), intent(in)    :: msg
    type(element_t),  intent(inout) :: elem(:)

    integer                         :: ie, i, j, k, ind
    real(r8), allocatable           :: test_vals(:)

    allocate(test_vals(np*np*nelemd))
    if (masterproc) then
      write(iulog, *) trim(msg)
    end if
    ! Check the upper levels of T against test function data
    do k = 1, 6
      ind = 0
      do ie = 1, nelemd
        do j = 1, np
          do i = 1, np
            ind = ind + 1
            test_vals(ind) = elem(ie)%derived%fT(i,j,k,1)
          end do
        end do
      end do
      call test_func_error(ind, dyn_lons, dyn_lats, test_vals, tend_factor, k)
    end do
    deallocate(test_vals)

  end subroutine dptest_test_dyn_state_ttends

  subroutine dptest_write_dyn_tend_history(elem)
    use cam_history,            only: outfld, hist_fld_active
    use element_mod,            only: element_t
    use dimensions_mod,         only: nlev, nelemd, np, npsq

    type(element_t), intent(inout) :: elem(:)

    real(r8)                       :: ftmp(npsq, nlev, 2)
    integer                        :: ie, i, j, k

    do ie = 1, nelemd
      if (hist_fld_active('dUglldbg') .or. hist_fld_active('dVglldbg')) then
        do k = 1, nlev
          do j = 1, np
            do i = 1, np
              ftmp(i+(j-1)*np,k,1) = elem(ie)%derived%FM(i,j,1,k,1)
              ftmp(i+(j-1)*np,k,2) = elem(ie)%derived%FM(i,j,2,k,1)
            end do
          end do
        end do
        call outfld('dUglldbg', ftmp(:,:,1), npsq, ie)
        call outfld('dVglldbg', ftmp(:,:,2), npsq, ie)
      end if

      if (hist_fld_active('dTglldbg')) then
        do k = 1, nlev
          do j = 1, np
            do i = 1, np
              ftmp(i+(j-1)*np,k,1) = elem(ie)%derived%FT(i,j,k,1)
            end do
          end do
        end do
        call outfld('dTglldbg', ftmp(:,:,1), npsq, ie)
      end if
    end do

  end subroutine dptest_write_dyn_tend_history

  subroutine dptest_set_phys_state_ttends(phys_tend)
    use physics_types,          only: physics_tend
    use cam_logfile,            only: iulog
    use spmd_utils,             only: masterproc

    type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)

    integer                           :: lchnk, k, b, e
    real(r8), allocatable             :: test_values(:)

    allocate(test_values(total_cols))
    ! Set the upper levels of T with test function data
    do k = 1, 6
      call set_test_func(total_cols, phys_lons(:), phys_lats(:),              &
           test_values, tend_factor, k)
      b = 1
      do lchnk = begchunk, endchunk
        e = b + ncols(lchnk) - 1
        phys_tend(lchnk)%dtdt(1:ncols(lchnk),k) = test_values(b:e)
        b = e + 1
      end do
    end do
    deallocate(test_values)

  end subroutine dptest_set_phys_state_ttends

  subroutine dptest_write_phys_tend_history_chunk(phys_tend, lchnk)
    use physics_types,          only: physics_tend
    use cam_history,            only: outfld

    type(physics_tend), intent(in) :: phys_tend
    integer,            intent(in) :: lchnk

    call outfld('dUdbg   ', phys_tend%dudt(:,:), pcols, lchnk)
    call outfld('dVdbg   ', phys_tend%dvdt(:,:), pcols, lchnk)
    call outfld('dTdbg   ', phys_tend%dtdt(:,:), pcols, lchnk)

  end subroutine dptest_write_phys_tend_history_chunk

  subroutine dptest_write_phys_tend_history_all(phys_tend)
    use physics_types,          only: physics_tend
    use cam_history,            only: outfld

    type(physics_tend), intent(in) :: phys_tend(begchunk:endchunk)

    integer                        :: lchnk

    do lchnk = begchunk, endchunk
      call dptest_write_phys_tend_history(phys_tend(lchnk), lchnk)
    end do

  end subroutine dptest_write_phys_tend_history_all

  subroutine dptest_write_phys_ptend_history_chunk(phys_ptend, lchnk)
    use physics_types,          only: physics_ptend
    use cam_history,            only: outfld

    type(physics_ptend), intent(in) :: phys_ptend
    integer,             intent(in) :: lchnk

    if (allocated(phys_ptend%u)) then
      call outfld('dUdbg   ', phys_ptend%u(:,:), pcols, lchnk)
    end if
    if (allocated(phys_ptend%v)) then
      call outfld('dVdbg   ', phys_ptend%v(:,:), pcols, lchnk)
    end if
    if (allocated(phys_ptend%s)) then
      call outfld('dTdbg   ', phys_ptend%s(:,:), pcols, lchnk)
    end if

  end subroutine dptest_write_phys_ptend_history_chunk

  subroutine dptest_write_phys_ptend_history_all(phys_ptend)
    use physics_types,          only: physics_ptend
    use cam_history,            only: outfld

    type(physics_ptend), intent(in) :: phys_ptend(begchunk:endchunk)

    integer                         :: lchnk

    do lchnk = begchunk, endchunk
      call dptest_write_phys_tend_history(phys_ptend(lchnk), lchnk)
    end do

  end subroutine dptest_write_phys_ptend_history_all

  subroutine dptest_write_tend_fields_history_chunk(fldu, fldv, flds, lchnk)
    use cam_history,            only: outfld

    real(r8),           intent(in) :: fldu(:,:)
    real(r8),           intent(in) :: fldv(:,:)
    real(r8),           intent(in) :: flds(:,:)
    integer,            intent(in) :: lchnk

    call outfld('dUdbg   ', fldu(:,:), pcols, lchnk)
    call outfld('dVdbg   ', fldv(:,:), pcols, lchnk)
    call outfld('dTdbg   ', flds(:,:), pcols, lchnk)
  end subroutine dptest_write_tend_fields_history_chunk

  subroutine compute_error_1d(array_val, array_check, array_error)
    real(r8),  intent(in)     :: array_val(:)
    real(r8),  intent(in)     :: array_check(:)
    real(r8),  intent(inout)  :: array_error(:)

    array_error = (array_val - array_check) / MAX(ABS(array_check), min_scale)
  end subroutine compute_error_1d

  subroutine compute_error_2d(array_val, array_check, array_error)
    real(r8),  intent(in)     :: array_val(:,:)
    real(r8),  intent(in)     :: array_check(:,:)
    real(r8),  intent(inout)  :: array_error(:,:)

    array_error = (array_val - array_check) / MAX(ABS(array_check), min_scale)
  end subroutine compute_error_2d

  subroutine print_array_stats_phys(array, message)
    use spmd_utils,             only: masterproc, masterprocid,               &
                                      mpicom, mpi_real8, npes,                &
                                      mpi_sum, mpi_min, mpi_max, MPI_SUCCESS
    use cam_logfile,            only: iulog

    real(r8),                   intent(in)  :: array(:,:)
    character(len=*), optional, intent(in)  :: message

    real(r8)                                :: temp
    real(r8)                                :: average, amin, amax
    integer                                 :: ierr

    ! Average
    temp = sum(abs(array)) / max(size(array),1)
    call MPI_Allreduce(temp, average, 1, mpi_real8, mpi_sum, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) call fatal_error(ierr, 'finding array average')
    average = average / npes

    ! Max
    temp = maxval(array)
    call MPI_Reduce(temp, amax, 1, mpi_real8, mpi_max, masterprocid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) call fatal_error(ierr, 'finding array maximum')

    ! Min
    temp = minval(array)
    call MPI_Reduce(temp, amin, 1, mpi_real8, mpi_min, masterprocid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) call fatal_error(ierr, 'finding array minimum')

    if (masterproc) then
      if(present(message)) then
        write(iulog, '(a)') trim(message)
      end if
      write(iulog, '(a, f10.6)') 'Array absolute value average = ', average
      write(iulog, '(a, f10.6)') 'Array maximum value = ', amax
      write(iulog, '(a, f10.6)') 'Array minimum value = ', amin
    end if

  end subroutine print_array_stats_phys

  subroutine print_array_stats_dyn(array, message)
    use spmd_utils,             only: masterproc, masterprocid,               &
                                      mpicom, mpi_real8, npes,                &
                                      mpi_sum, mpi_min, mpi_max, MPI_SUCCESS
    use cam_logfile,            only: iulog

    real(r8),                   intent(in)  :: array(:,:,:)
    character(len=*), optional, intent(in)  :: message

    real(r8)                                :: temp
    real(r8)                                :: average, amin, amax
    integer                                 :: ierr

    ! Average
    temp = sum(abs(array)) / max(size(array),1)
    call MPI_Allreduce(temp, average, 1, mpi_real8, mpi_sum, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) call fatal_error(ierr, 'finding array average')
    average = average / npes

    ! Max
    temp = maxval(array)
    call MPI_Reduce(temp, amax, 1, mpi_real8, mpi_max, masterprocid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) call fatal_error(ierr, 'finding array maximum')

    ! Min
    temp = minval(array)
    call MPI_Reduce(temp, amin, 1, mpi_real8, mpi_min, masterprocid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) call fatal_error(ierr, 'finding array minimum')

    if (masterproc) then
      if(present(message)) then
        write(iulog, '(a)') trim(message)
      end if
      write(iulog, '(a, f10.6)') 'Array average absolute value = ', average
      write(iulog, '(a, f10.6)') 'Array average maximum value = ', amax
      write(iulog, '(a, f10.6)') 'Array average minimum value = ', amin
    end if

  end subroutine print_array_stats_dyn

#if 0
  !! Obsolete?
  subroutine dptest_init_grid_file(elem, diag_fname, nphys, test_num, coords)
    use element_mod,     only: element_t
    use dimensions_mod,  only: nelemd, np, nelem
    use shr_pio_mod,     only: shr_pio_getiosys
    use cam_instance,    only: atm_id
    use pio_kinds,       only: i4
    use coordinate_systems_mod, only: spherical_polar_t
!    use spmd_utils,      only: masterproc
!    use cam_logfile,     only: iulog
    use pio_types,       only: iosystem_desc_t
    use pio,             only: pio_seterrorhandling, pio_createfile,          &
                               pio_def_dim, pio_def_var, pio_initdecomp,      &
                               PIO_BCAST_ERROR, PIO_iotype_netcdf,            &
                               PIO_clobber, PIO_unlimited, PIO_DOUBLE,        &
                               PIO_global, pio_enddef, pio_put_att

    type(element_t),       intent(in)  :: elem(:)
    character(len=*),      intent(in)  :: diag_fname
    integer,               intent(in)  :: nphys
    integer,     optional, intent(in)  :: test_num
    logical,     optional, intent(in)  :: coords

    type(iosystem_desc_t), pointer     :: pio_subsystem
    integer(i4),           allocatable :: gdof(:)   ! Mapping to file
    integer(i4)                        :: pdims(2), ddims(3) ! For PIO
    type(var_desc_t)                   :: dyngrid_lat_id, dyngrid_lon_id
    type(var_desc_t)                   :: physgrid_lat_id, physgrid_lon_id
    type(spherical_polar_t), pointer   :: dyn_sphere(:,:,:)
    type(spherical_polar_t), pointer   :: phys_sphere(:,:)
    integer                            :: nphys2
    integer                            :: curr_peh
    integer                            :: ierr, i, j, ie

    nphys2 = nphys * nphys
    pio_subsystem => shr_pio_getiosys(atm_id)
    ! Hack to get around no pio_geterrorhandling
    curr_peh = pio_subsystem%error_handling
    call pio_seterrorhandling(pio_subsystem, PIO_BCAST_ERROR)
    ierr = pio_createfile(pio_subsystem, diagfile, PIO_iotype_netcdf,         &
         diag_fname, PIO_clobber)
    if (ierr /= 0) call fatal_error(ierr, 'unable to create diag file')
    ierr = pio_def_dim(diagfile, 'np', np, np_dim)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define np dimension')
    ierr = pio_def_dim(diagfile, 'nphys2', nphys2, nphys2_dim)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define nphys2 dimension')
    ierr = pio_def_dim(diagfile, 'nelem', nelem, nelem_dim)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define nelemd dimension')
    ierr = pio_def_dim(diagfile, 'time', PIO_unlimited, time_dim)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define time dimension')
    ierr = pio_def_var(diagfile, 'dyn_grid', pio_double,                      &
         (/np_dim, np_dim, nelem_dim, time_dim/), dyn_grid_id)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define dyn_grid variable')
    ierr = pio_def_var(diagfile, 'phys_grid', pio_double,                     &
         (/nphys2_dim, nelem_dim, time_dim/), phys_grid_id)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define phys_grid variable')
    ierr = pio_def_var(diagfile, 'dyn_grid_error', pio_double,                &
         (/np_dim, np_dim, nelem_dim, time_dim/), dyn_grid_error_id)
    if (ierr /= 0) call fatal_error(ierr, 'unable to define dyn_grid variable')
    if (present(coords) .and. coords) then
      ierr = pio_def_var(diagfile, 'dyn_grid_lat', pio_double,                &
           (/np_dim, np_dim, nelem_dim/), dyngrid_lat_id)
      if (ierr /= 0) then
        call fatal_error(ierr, 'unable to define dyn_grid_lat variable')
      end if
      ierr = pio_put_att(diagfile, dyngrid_lat_id, 'units', 'degrees')
      if (ierr /= 0) call fatal_error(ierr, 'unable to put "units" attribute')
      ierr = pio_def_var(diagfile, 'dyn_grid_lon', pio_double,                &
           (/np_dim, np_dim, nelem_dim/), dyngrid_lon_id)
      if (ierr /= 0) then
        call fatal_error(ierr, 'unable to define dyn_grid_lon variable')
      end if
      ierr = pio_put_att(diagfile, dyngrid_lon_id, 'units', 'degrees')
      if (ierr /= 0) call fatal_error(ierr, 'unable to put "units" attribute')
      ierr = pio_def_var(diagfile, 'phys_grid_lat', pio_double,               &
           (/nphys2_dim, nelem_dim/), physgrid_lat_id)
      if (ierr /= 0) then
        call fatal_error(ierr, 'unable to define phys_grid_lat variable')
      end if
      ierr = pio_put_att(diagfile, physgrid_lat_id, 'units', 'degrees')
      if (ierr /= 0) call fatal_error(ierr, 'unable to put "units" attribute')
      ierr = pio_def_var(diagfile, 'phys_grid_lon', pio_double,               &
           (/nphys2_dim, nelem_dim/), physgrid_lon_id)
      if (ierr /= 0) then
        call fatal_error(ierr, 'unable to define phys_grid_lon variable')
      end if
      ierr = pio_put_att(diagfile, physgrid_lon_id, 'units', 'degrees')
      if (ierr /= 0) call fatal_error(ierr, 'unable to put "units" attribute')
    end if
    if (present(test_num)) then
      ierr = pio_put_att(diagfile, PIO_global, 'Test',                        &
           trim(dptest_test_func_string(test_num)))
      if (ierr /= 0) call fatal_error(ierr, 'unable to put "Test" attribute')
    end if
    ierr = pio_enddef(diagfile)
    if (ierr /= 0) call fatal_error(ierr, 'unable to exit NetCDF define mode')
    ! Get physics decomposition
    pdims(1) = nphys2
    pdims(2) = nelem
    allocate(gdof(nphys2 * nelemd))
    do ie = 1, nelemd
       do i = 1, nphys2
          gdof(i + ((ie - 1) * nphys2)) =                                     &
               i + ((elem(ie)%GlobalId - 1) * nphys2)
       enddo
    enddo
!    if (masterproc) then
!      write(iulog, '(a,i2,a,i5,a,i5,a)') 'Calling pio_initdecomp with pdims(',&
!           pdims(1), ',', pdims(2), '), gdof(', size(gdof), ')'
!    end if
    call pio_initdecomp(pio_subsystem, PIO_DOUBLE, pdims, gdof, phys_iodesc)
    deallocate(gdof)

    ! Get dynamics decomposition
    ddims(1) = np
    ddims(2) = np
    ddims(3) = nelem
    allocate(gdof(np * np * nelemd))
    do ie = 1, nelemd
       do i = 1, np
          do j = 1, np
             gdof(j + (((i - 1) + (ie - 1) * np) * np)) =                  &
                  j + (((i - 1) + (elem(ie)%GlobalId - 1) * np) * np)
          enddo
       enddo
    enddo
    call pio_initdecomp(pio_subsystem, PIO_DOUBLE, ddims, gdof, dyn_iodesc)
    deallocate(gdof)

    ! Dump coordinate variables if instructed
    if (present(coords) .and. coords) then
      call dptest_init_coords(elem, nphys)
      call dptest_get_dyn_coords(dyn_sphere)
      call dptest_get_phys_coords(phys_sphere)
      call write_grid_file_coord_3d(dyn_sphere, dyngrid_lat_id, dyngrid_lon_id)
      call write_grid_file_coord_2d(phys_sphere, nphys2, physgrid_lat_id,     &
           physgrid_lon_id)
    end if

    ! Reset PIO error handling
    call pio_seterrorhandling(pio_subsystem, curr_peh)
  end subroutine dptest_init_grid_file

  subroutine write_grid_file_2d(array, frame, id)
    use cam_abortutils, only: endrun
    use shr_pio_mod,    only: shr_pio_getiosys
    use cam_instance,   only: atm_id
    use pio_types,      only: iosystem_desc_t
    use pio,            only: pio_seterrorhandling, pio_write_darray,         &
                              pio_setframe, PIO_BCAST_ERROR, pio_offset
!    use pio,             only: pio_setdebuglevel

    real(r8),      dimension(:,:), intent(in)  :: array
    integer(pio_offset),           intent(in)  :: frame
    character(len=*),              intent(in)  :: id

    type(iosystem_desc_t), pointer             :: pio_subsystem
    integer                                    :: ierr
    type(var_desc_t), pointer                  :: file_id
    integer                                    :: curr_peh

    if (trim(id) .eq. 'phys_grid') then
      file_id => phys_grid_id
    else
      call endrun('Unknown field id, '//trim(id))
    end if
    pio_subsystem => shr_pio_getiosys(atm_id)
    ! Hack to get around no pio_geterrorhandling
    curr_peh = pio_subsystem%error_handling
    call pio_seterrorhandling(pio_subsystem, PIO_BCAST_ERROR)
    call pio_setframe(file_id, frame)
    call pio_write_darray(diagfile, file_id, phys_iodesc, array, ierr)
    if (ierr /= 0) call fatal_error(ierr, 'writing 2d array')
    call pio_seterrorhandling(pio_subsystem, curr_peh)
  end subroutine write_grid_file_2d

  subroutine write_grid_file_3d(array, frame, id)
    use cam_abortutils,  only: endrun
    use shr_pio_mod,     only: shr_pio_getiosys
    use cam_instance,    only: atm_id
    use pio_types,       only: iosystem_desc_t
    use pio,             only: pio_seterrorhandling, pio_write_darray,        &
                               pio_setframe, PIO_BCAST_ERROR, pio_offset

    real(r8),    dimension(:,:,:), intent(in)  :: array
    integer(pio_offset),           intent(in)  :: frame
    character(len=*),              intent(in)  :: id

    type(iosystem_desc_t), pointer             :: pio_subsystem
    integer                                    :: ierr
    type(var_desc_t), pointer                  :: file_id
    integer                                    :: curr_peh

    if (trim(id) .eq. 'dyn_grid') then
      file_id => dyn_grid_id
    else if(trim(id) .eq. 'dyn_grid_error') then
      file_id => dyn_grid_error_id
    else
      call endrun('Unknown field id, '//trim(id))
    end if

    pio_subsystem => shr_pio_getiosys(atm_id)
    ! Hack to get around no pio_geterrorhandling
    curr_peh = pio_subsystem%error_handling
    call pio_seterrorhandling(pio_subsystem, PIO_BCAST_ERROR)
    call pio_setframe(file_id, frame)
    call pio_write_darray(diagfile, file_id, dyn_iodesc, array, ierr)
    if (ierr /= 0) call fatal_error(ierr, 'writing 3d array')
    call pio_seterrorhandling(pio_subsystem, curr_peh)
  end subroutine write_grid_file_3d

  subroutine write_grid_file_coord_2d(sarray, nphys2, lat_id, lon_id)
    use coordinate_systems_mod, only: spherical_polar_t
    use dimensions_mod,  only: nelemd
    use shr_pio_mod,     only: shr_pio_getiosys
    use cam_instance,    only: atm_id
    use pio_types,       only: iosystem_desc_t
    use pio,             only: pio_seterrorhandling, pio_write_darray,        &
                               PIO_BCAST_ERROR

    type(spherical_polar_t), pointer, intent(in)    :: sarray(:,:)
    integer,                          intent(in)    :: nphys2
    type(var_desc_t),                 intent(inout) :: lat_id
    type(var_desc_t),                 intent(inout) :: lon_id

    type(iosystem_desc_t),   pointer                :: pio_subsystem
    real(r8)                                        :: array(nphys2, nelemd)
    real(r8)                                        :: conv
    integer                                         :: ierr, i, ie
    integer                                         :: curr_peh

    ! Conversion between radians and degrees
    conv = 90.0_r8 / ASIN(1.0_r8)
    ! Retrieve the PIO iosys
    pio_subsystem => shr_pio_getiosys(atm_id)
    ! Hack to get around no pio_geterrorhandling
    curr_peh = pio_subsystem%error_handling
    call pio_seterrorhandling(pio_subsystem, PIO_BCAST_ERROR)
    do ie = 1, nelemd
      do i = 1, nphys2
        array(i, ie) = sarray(i, ie)%lat * conv
      end do
    end do
    call pio_write_darray(diagfile, lat_id, phys_iodesc, array, ierr)
    if (ierr /= 0) call fatal_error(ierr, 'writing 2d lat coordinate variable')
    do ie = 1, nelemd
      do i = 1, nphys2
        array(i, ie) = sarray(i, ie)%lon * conv
      end do
    end do
    call pio_write_darray(diagfile, lon_id, phys_iodesc, array, ierr)
    if (ierr /= 0) call fatal_error(ierr, 'writing 2d lon coordinate variable')
    call pio_seterrorhandling(pio_subsystem, curr_peh)
  end subroutine write_grid_file_coord_2d

  subroutine write_grid_file_coord_3d(sarray, lat_id, lon_id)
    use coordinate_systems_mod, only: spherical_polar_t
    use dimensions_mod,  only: nelemd, np
    use shr_pio_mod,     only: shr_pio_getiosys
    use cam_instance,    only: atm_id
    use pio_types,       only: iosystem_desc_t
    use pio,             only: pio_seterrorhandling, pio_write_darray,        &
                               PIO_BCAST_ERROR

    type(spherical_polar_t), pointer, intent(in)    :: sarray(:,:,:)
    type(var_desc_t),                 intent(inout) :: lat_id
    type(var_desc_t),                 intent(inout) :: lon_id

    type(iosystem_desc_t),   pointer                :: pio_subsystem
    real(r8)                                        :: array(np, np, nelemd)
    real(r8)                                        :: conv
    integer                                         :: ierr, i, j, ie
    integer                                         :: curr_peh

    ! Conversion between radians and degrees
    conv = 90.0_r8 / ASIN(1.0_r8)
    ! Retrieve the PIO iosys
    pio_subsystem => shr_pio_getiosys(atm_id)
    ! Hack to get around no pio_geterrorhandling
    curr_peh = pio_subsystem%error_handling
    call pio_seterrorhandling(pio_subsystem, PIO_BCAST_ERROR)
    do ie = 1, nelemd
      do i = 1, np
        do j = 1, np
          array(j, i, ie) = sarray(j, i, ie)%lat * conv
        end do
      end do
    end do
    call pio_write_darray(diagfile, lat_id, dyn_iodesc, array, ierr)
    if (ierr /= 0) call fatal_error(ierr, 'writing 3d lat coordinate variable')
    do ie = 1, nelemd
      do i = 1, np
        do j = 1, np
          array(j, i, ie) = sarray(j, i, ie)%lon * conv
        end do
      end do
    end do
    call pio_write_darray(diagfile, lon_id, dyn_iodesc, array, ierr)
    if (ierr /= 0) call fatal_error(ierr, 'writing 3d lon coordinate variable')
    call pio_seterrorhandling(pio_subsystem, curr_peh)
  end subroutine write_grid_file_coord_3d

  subroutine dptest_close_grid_file()
    use shr_pio_mod,     only: shr_pio_getiosys
    use cam_instance,    only: atm_id
    use pio_types,       only: iosystem_desc_t
    use pio,             only: pio_seterrorhandling, pio_closefile,           &
                               pio_freedecomp, PIO_INTERNAL_ERROR

    type(iosystem_desc_t), pointer             :: pio_subsystem

    pio_subsystem => shr_pio_getiosys(atm_id)
    call pio_seterrorhandling(pio_subsystem, PIO_INTERNAL_ERROR)
    call pio_freedecomp(pio_subsystem, phys_iodesc)
    call pio_freedecomp(pio_subsystem, dyn_iodesc)
    call pio_closefile(diagfile)
  end subroutine dptest_close_grid_file

  subroutine dptest_bndry_exchange(elem, array)
    use element_mod,     only: element_t
    use dimensions_mod,  only: nelemd
    use spmd_utils,      only: iam
    use parallel_mod,    only: par
    use bndry_mod,       only: bndry_exchangeV
    use edge_mod,        only: initEdgeBuffer, edgeVpack, edgeVunpack
    use edgetype_mod,    only: EdgeBuffer_t

    type(element_t), dimension(:), intent(in)     :: elem
    real(r8),    dimension(:,:,:), intent(inout)  :: array

    type (EdgeBuffer_t), save                     :: edgebuf    ! edge buffer
    logical,             save                     :: initialized = .false.
    integer                                       :: ie

    if (.not. initialized) then
      if(iam < par%nprocs) then
        ! Only going to exchange the one test field
        call initEdgeBuffer(par, edgebuf, elem, 1,nthreads=1)
      end if
      initialized = .true.
    endif

    do ie = 1, nelemd
      array(:,:,ie) = array(:,:,ie)*elem(ie)%spheremp(:,:)
      call edgeVpack(edgebuf,array(:,:,ie),1,0,ie)
    enddo

    call bndry_exchangeV(par, edgebuf,location='dptest_bndry_exchange')
    do ie = 1, nelemd
      call edgeVunpack(edgebuf,array(:,:,ie),1,0,ie)
      array(:,:,ie) = array(:,:,ie)*elem(ie)%rspheremp(:,:)
    enddo

  end subroutine dptest_bndry_exchange

  function get_latlon(elem, i, j, coords) result(sphere)
    use element_mod,            only: element_t
    use coordinate_systems_mod, only: spherical_polar_t
    use cube_mod,               only: ref2sphere

    type(element_t),         intent(in)  :: elem
    integer,                 intent(in)  :: i
    integer,                 intent(in)  :: j
    real(r8),                intent(in)  :: coords(:)
    type(spherical_polar_t)              :: sphere

    sphere = ref2sphere(coords(i), coords(j), elem)
!%corners, elem%vertex%face_number)

  end function get_latlon
#endif
!! End obsolete?
#endif
!! End TEST_GRID_COUPLING

end module test_dp_grids
