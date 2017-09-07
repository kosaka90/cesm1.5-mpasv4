module dyn_inic_held_suarez

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set Held-Suarez initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  implicit none
  private

  ! Public interface
  public :: hs94_dyn_set_inic

!==============================================================================
CONTAINS
!==============================================================================

  subroutine hs94_dyn_set_inic(latvals, lonvals, U, V, T, PS, PHIS,           &
       Q, m_cnst, mask, verbose)
    use constituents,  only: cnst_name

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set Held-Suarez initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical           , intent(in)    :: verbose    ! For internal use

    ! Local variables
    logical, allocatable              :: mask_use(:)
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'HS94_DYN_SET_INIC'

    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun('cnst_init_default: input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    ncol = size(latvals, 1)
    nlev = -1
    if (present(U)) then
      nlev = size(U, 2)
      do k = 1, nlev
        where(mask_use)
          U(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose) then
        write(iulog,*) '          U initialized by "',subname,'"'
      end if
    end if

    if (present(V)) then
      nlev = size(V, 2)
      do k = 1, nlev
        where(mask_use)
          V(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose) then
        write(iulog,*) '          V initialized by "',subname,'"'
      end if
    end if

    if (present(T)) then
      nlev = size(T, 2)
      do k = 1, nlev
        where(mask_use)
          T(:,k) = 250.0_r8
        end where
      end do
      if(masterproc .and. verbose) then
        write(iulog,*) '          T initialized by "',subname,'"'
      end if
    end if

    if (present(PS)) then
      where(mask_use)
        PS = 100000.0_r8
      end where
      if(masterproc .and. verbose) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if

    if (present(PHIS)) then
      where(mask_use)
        PHIS = 0.0_r8
      end where
      if(masterproc .and. verbose) then
        write(iulog,*) '          V initialized by "',subname,'"'
      end if
    end if

    if (present(Q)) then
      nlev = size(Q, 2)
      ncnst = size(m_cnst, 1)
      do m = 1, ncnst
        do k = 1, nlev
          where(mask_use)
            Q(:,k,m_cnst(m)) = 0.0_r8
          end where
        end do
        if(masterproc .and. verbose) then
          write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by "',subname,'"'
        end if
      end do
    end if

    deallocate(mask_use)

  end subroutine hs94_dyn_set_inic

end module dyn_inic_held_suarez
