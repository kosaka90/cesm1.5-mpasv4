!
!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling
  use constituents,   only: pcnst, cnst_name
  use cam_history,    only: outfld, write_inithist, hist_fld_active
  use dimensions_mod, only: np, npsq, nelemd, nlev, nc, qsize, ntrac
  use dof_mod,        only: UniquePoints, PutUniquePoints
  use dyn_comp,       only: dyn_export_t, dyn_import_t
  use dyn_grid,       only: get_gcol_block_d, TimeLevel
  use element_mod,    only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use kinds,          only: real_kind, int_kind
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physics_types,  only: physics_state, physics_tend

  use phys_grid,      only: get_ncols_p, get_gcol_all_p, block_to_chunk_send_pters, transpose_block_to_chunk, &
       block_to_chunk_recv_pters, chunk_to_block_send_pters, transpose_chunk_to_block, chunk_to_block_recv_pters
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
  use control_mod,    only: smooth_phis_numcycle
  use cam_logfile,    only : iulog
  use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use spmd_utils,   only: mpicom, iam
  use perf_mod,    only : t_startf, t_stopf, t_barrierf
  use parallel_mod, only : par
  private

  public :: d_p_coupling, p_d_coupling

  real (kind=r8), allocatable :: q_prev(:,:,:,:) ! Previous Q for computing tendencies

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(phys_state, phys_tend,  pbuf2d, dyn_out)
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, &
                                     pbuf_get_field
    use shr_vmath_mod,          only: shr_vmath_exp
    use time_manager,           only: is_first_step
    use viscosity_mod,          only: compute_zeta_C0
    use cam_abortutils,         only: endrun
    use gravity_waves_sources,  only: gws_src_fnct
    use dyn_comp,               only: frontgf_idx, frontga_idx
    use phys_control,           only: use_gw_front, use_gw_front_igw
    use dp_grids,               only: nphys_pts
    use dimensions_mod,         only: fv_nphys
    use hycoef,                 only: hyai, hybi, ps0
    use fvm_control_volume_mod, only: n0_fvm
    use dimensions_mod,         only: ldry_mass_vertical_coordinates
    use derivative_mod, only : subcell_integration
    use thread_mod, only : horz_num_threads
    use fvm_mapping     , only: dyn2phys_vector, dyn2phys_all_vars
    implicit none
!-----------------------------------------------------------------------
! !INPUT PARAMETERS:
!
    type(dyn_export_t), intent(inout) :: dyn_out    ! dynamics export
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
! !OUTPUT PARAMETERS:

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend


! LOCAL VARIABLES
    type(element_t), pointer :: elem(:)               ! pointer to dyn_out element array
    integer (kind=int_kind)  :: ie               ! indices over elements
    integer (kind=int_kind)  :: lchnk, icol, ilyr      ! indices over chunks, columns, layers

    real (kind=r8), allocatable :: ps_tmp(:,:)        ! temporary array to hold ps
    real (kind=r8), allocatable :: dp3d_tmp(:,:,:)    ! temporary array to hold dp3d
    real (kind=r8), allocatable :: dp3d_tmp_tmp(:,:)
    real (kind=r8), allocatable :: phis_tmp(:,:)      ! temporary array to hold phis
    real (kind=r8), allocatable :: T_tmp(:,:,:)       ! temporary array to hold T
    real (kind=r8), allocatable :: uv_tmp(:,:,:,:)    ! temporary array to hold u and v
    real (kind=r8), allocatable :: q_tmp(:,:,:,:)     ! temporary to hold advected constituents
    real (kind=r8), allocatable :: omega_tmp(:,:,:)   ! temporary array to hold omega
    real (kind=r8), allocatable :: inv_darea_dp_fvm(:,:)
!!XXgoldyXX: v Need to add GW to physgrid
    ! Frontogenesis
    real (kind=r8), allocatable :: frontgf(:,:,:) ! temporary arrays to hold frontogenesis
    real (kind=r8), allocatable :: frontga(:,:,:) !   function (frontgf) and angle (frontga)
    ! Pointers to pbuf
    real (kind=r8),        pointer     :: pbuf_frontgf(:,:)
    real (kind=r8),        pointer     :: pbuf_frontga(:,:)
!!XXgoldyXX: ^ Need to add GW to physgrid

    integer :: ncols,i,j,ierr,k,iv

    integer (kind=int_kind)  :: ioff, m, m_cnst
    integer :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer, allocatable :: bpter(:,:)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real (kind=r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers
    integer :: tl_f

    logical :: lmono

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    !----------------------------------------------------------------------

    nullify(pbuf_chnk)
    nullify(pbuf_frontgf)
    nullify(pbuf_frontga)

    ! Allocate temporary arrays to hold data for physics decomposition
    allocate(ps_tmp(nphys_pts,nelemd))
    allocate(dp3d_tmp(nphys_pts,pver,nelemd))
    allocate(dp3d_tmp_tmp(nphys_pts,pver))
    allocate(phis_tmp(nphys_pts,nelemd))
    allocate(T_tmp(nphys_pts,pver,nelemd))
    allocate(uv_tmp(nphys_pts,2,pver,nelemd))
    allocate(q_tmp(nphys_pts,pver,MAX(qsize,ntrac),nelemd))
    allocate(omega_tmp(nphys_pts,pver,nelemd))

    if (use_gw_front .or. use_gw_front_igw) then

      allocate(frontgf(npsq,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of frontgf failed.")

      allocate(frontga(npsq,pver,nelemd), stat=ierr)
      if (ierr /= 0) call endrun("dp_coupling: Allocate of frontga failed.")

    end if
    if( iam < par%nprocs) then

      elem => dyn_out%elem

      tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)

#ifdef debug_coupling
!      call dptest_write_dyn_state_history(elem)
#endif
      if (fv_nphys>0) then
         !
         !******************************************************************
         !
         ! physics grid: map GLL vars to physics grid
         !
         !******************************************************************
         !
         call t_startf('dyn2phys')
!         allocate(inv_darea_dp_fvm(fv_nphys,fv_nphys))
!         allocate(inv_darea(fv_nphys,fv_nphys))
!         allocate(dp_phys(fv_nphys,fv_nphys,nlev))
!         allocate(q_phys(fv_nphys,fv_nphys,nlev,qsize))
         do ie = 1,nelemd
           !
           ! note that the fvm halo has been filled in prim_run_subcycle
           ! of physics grid resolution is not equal to fvm resolution
           !
           call dyn2phys_all_vars(ie,&
                 !
                 ! spectral element state
                 !
                 elem(ie)%state%dp3d(:,:,:,tl_f),elem(ie)%state%ps(:,:,tl_f),&
                 elem(ie)%state%q(:,:,:,1:qsize),elem(ie)%state%T(:,:,:,tl_f),&
                 elem(ie)%derived%omega(:,:,:),elem(ie)%state%phis(:,:),&
                 !
                 ! fvm state
                 !
                 dyn_out%fvm(ie)%dp_fvm(:,:,:,n0_fvm),dyn_out%fvm(ie)%psC(:,:),&
                 dyn_out%fvm(ie)%c(:,:,:,1:qsize,n0_fvm),&
                 dyn_out%fvm(ie)%inv_area_sphere(:,:),&
                 qsize,elem(ie)%metdet,ntrac>0,&
                 dyn_out%fvm(ie),&
                 !
                 hyai(1)*ps0,&
                 !
                 ! output
                 !
                 dp3d_tmp(:,:,ie),ps_tmp(:,ie),q_tmp(:,:,:,ie),T_tmp(:,:,ie),omega_tmp(:,:,ie),&
                 phis_tmp(:,ie)&
                 )
               uv_tmp(:,:,:,ie) = &
                    dyn2phys_vector(elem(ie)%state%v(:,:,:,:,tl_f),elem(ie))
         end do
         call t_stopf('dyn2phys')
         !!XXgoldy!!: v How do we want to do this for the physics grid?
         if (use_gw_front .or. use_gw_front_igw) then
            call endrun("d_p_coupling: Gravity waves not supported on physics grid")
            !          call gws_src_fnct(elem, tl_f, frontgf, frontga)
         end if
      else
        !
        !******************************************************************
        !
        ! No physics grid
        !
        !******************************************************************
        !
        if (qsize < pcnst) then
          call endrun('d_p_coupling: Fewer GLL tracers advected than required')
        end if
        call t_startf('UniquePoints')
        do ie=1,nelemd
          ncols = elem(ie)%idxP%NumUniquePts
          call UniquePoints(elem(ie)%idxP, elem(ie)%state%ps(:,:,tl_f), ps_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%dp3d(:,:,:,tl_f), dp3d_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxP, nlev, elem(ie)%state%T(:,:,:,tl_f), T_tmp(1:ncols,:,ie))
          call UniquePoints(elem(ie)%idxV, 2, nlev, elem(ie)%state%V(:,:,:,:,tl_f), uv_tmp(1:ncols,:,:,ie))
          call UniquePoints(elem(ie)%idxV, nlev, elem(ie)%derived%omega, omega_tmp(1:ncols,:,ie))

          call UniquePoints(elem(ie)%idxP, elem(ie)%state%phis, phis_tmp(1:ncols,ie))
          call UniquePoints(elem(ie)%idxP, nlev,pcnst, elem(ie)%state%Q(:,:,:,:), Q_tmp(1:ncols,:,:,ie))
        end do
        call t_stopf('UniquePoints')
      end if!physgrid
      if (use_gw_front .or. use_gw_front_igw) then
        call gws_src_fnct(elem, tl_f, frontgf, frontga)
      end if
    else
      ps_tmp(:,:) = 0._r8
      T_tmp(:,:,:) = 0._r8
      uv_tmp(:,:,:,:) = 0._r8
      omega_tmp(:,:,:) = 0._r8
      if (use_gw_front .or. use_gw_front_igw) then
        frontgf(:,:,:) = 0._r8
        frontga(:,:,:) = 0._r8
      end if
      phis_tmp(:,:) = 0._r8
      Q_tmp(:,:,:,:) = 0._r8
    endif !! iam .lt. par%nprocs

    ! q_prev is for saving the tracer fields for calculating tendencies
    if (.not. allocated(q_prev)) then
      allocate(q_prev(pcols, pver, pcnst, begchunk:endchunk))
    end if
    q_prev = 0.0_R8

    call t_startf('dpcopy')
    if (local_dp_map) then

      !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
      do lchnk=begchunk,endchunk
        ncols=get_ncols_p(lchnk)
        call get_gcol_all_p(lchnk,pcols,pgcols)

        pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

        if (use_gw_front .or. use_gw_front_igw) then
          call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
          call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
        end if

        do icol=1,ncols
          call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
          ie = idmb3(1)
          ioff=idmb2(1)
          phys_state(lchnk)%ps(icol)=ps_tmp(ioff,ie)
          phys_state(lchnk)%phis(icol)=phis_tmp(ioff,ie)
          do ilyr=1,pver
            phys_state(lchnk)%pdel(icol,ilyr)=dp3d_tmp(ioff,ilyr,ie)
            phys_state(lchnk)%t(icol,ilyr)=T_tmp(ioff,ilyr,ie)
            phys_state(lchnk)%u(icol,ilyr)=uv_tmp(ioff,1,ilyr,ie)
            phys_state(lchnk)%v(icol,ilyr)=uv_tmp(ioff,2,ilyr,ie)
            phys_state(lchnk)%omega(icol,ilyr)=omega_tmp(ioff,ilyr,ie)

            if (use_gw_front .or. use_gw_front_igw) then
              pbuf_frontgf(icol,ilyr) = frontgf(ioff,ilyr,ie)
              pbuf_frontga(icol,ilyr) = frontga(ioff,ilyr,ie)
            endif
          end do

          do m=1,pcnst
            do ilyr=1,pver
              phys_state(lchnk)%q(icol,ilyr,m)=Q_tmp(ioff,ilyr,m,ie)
            end do
          end do
        end do
      end do

    else  ! .not. local_dp_map

      tsize = 5 + pcnst
      if (use_gw_front .or. use_gw_front_igw) tsize = tsize + 2

      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))
      if (fv_nphys>0) then
        allocate(bpter(fv_nphys*fv_nphys,0:pver))
      else
        allocate(bpter(npsq,0:pver))
      end if
      if(iam .lt. par%nprocs) then
        !$omp parallel do num_threads(horz_num_threads) private (ie, bpter, icol, ilyr, m, ncols, ioff)
        do ie=1,nelemd

          if (fv_nphys>0) then
            call block_to_chunk_send_pters(elem(ie)%GlobalID, fv_nphys*fv_nphys,  &
                 pver+1, tsize, bpter)
            ncols = fv_nphys*fv_nphys
          else
            call block_to_chunk_send_pters(elem(ie)%GlobalID, npsq,       &
                 pver+1, tsize, bpter)
            ncols = elem(ie)%idxP%NumUniquePts
          end if

          do icol=1,ncols
            bbuffer(bpter(icol,0)+2:bpter(icol,0)+tsize-1) = 0.0_r8
            bbuffer(bpter(icol,0))   = ps_tmp(icol,ie)
            bbuffer(bpter(icol,0)+1) = phis_tmp(icol,ie)

            do ilyr=1,pver
              ioff = 0
              bbuffer(bpter(icol,ilyr)+ioff) = T_tmp(icol,ilyr,ie)
              ioff = ioff + 1
              bbuffer(bpter(icol,ilyr)+ioff) = uv_tmp(icol,1,ilyr,ie)
              ioff = ioff + 1
              bbuffer(bpter(icol,ilyr)+ioff) = uv_tmp(icol,2,ilyr,ie)
              ioff = ioff + 1
              bbuffer(bpter(icol,ilyr)+ioff) = omega_tmp(icol,ilyr,ie)
              ioff = ioff + 1
              bbuffer(bpter(icol,ilyr)+ioff) = dp3d_tmp(icol,ilyr,ie)
              if (use_gw_front .or. use_gw_front_igw) then
                ioff = ioff + 1
                bbuffer(bpter(icol,ilyr)+ioff) = frontgf(icol,ilyr,ie)
                ioff = ioff + 1
                bbuffer(bpter(icol,ilyr)+ioff) = frontga(icol,ilyr,ie)
              end if

              do m=1,pcnst
                bbuffer(bpter(icol,ilyr)+tsize-pcnst-1+m) = Q_tmp(icol,ilyr,m,ie)
              end do
            end do

          end do

        end do
      else
        bbuffer(:) = 0._r8
      end if

      call t_barrierf ('sync_blk_to_chk', mpicom)
      call t_startf ('block_to_chunk')
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      call t_stopf  ('block_to_chunk')
      !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, icol, ilyr, m, pbuf_chnk, pbuf_frontgf, pbuf_frontga, ioff)
      do lchnk = begchunk,endchunk
        ncols = phys_state(lchnk)%ncol

        pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

        if (use_gw_front .or. use_gw_front_igw) then
          call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
          call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
        end if

        call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

        do icol=1,ncols
          phys_state(lchnk)%ps(icol)   = cbuffer(cpter(icol,0))
          phys_state(lchnk)%phis(icol) = cbuffer(cpter(icol,0)+1)

          do ilyr=1,pver
            ioff = 0
            phys_state(lchnk)%t    (icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
            ioff = ioff + 1
            phys_state(lchnk)%u    (icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
            ioff = ioff + 1
            phys_state(lchnk)%v    (icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
            ioff = ioff + 1
            phys_state(lchnk)%omega(icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
            ioff = ioff + 1
            phys_state(lchnk)%pdel (icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)

            if (use_gw_front .or. use_gw_front_igw) then
              ioff = ioff + 1
              pbuf_frontgf(icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
              ioff = ioff + 1
              pbuf_frontga(icol,ilyr) = cbuffer(cpter(icol,ilyr)+ioff)
            endif

            do m=1,pcnst
              phys_state(lchnk)%q  (icol,ilyr,m) = cbuffer(cpter(icol,ilyr)+tsize-pcnst-1+m)
            end do

          end do

        end do

      end do

      deallocate( bbuffer )
      deallocate( cbuffer )

    end if
    call t_stopf('dpcopy')

    ! Save the tracer fields for calculating tendencies
    do lchnk = begchunk, endchunk
      ncols = phys_state(lchnk)%ncol
      q_prev(1:ncols, 1:pver, 1:pcnst, lchnk) = phys_state(lchnk)%q(1:ncols, 1:pver, 1:pcnst)
    end do


    ! Deallocate the temporary arrays
    deallocate(ps_tmp)
    deallocate(dp3d_tmp)
    deallocate(phis_tmp)
    deallocate(T_tmp)
    deallocate(uv_tmp)
    deallocate(q_tmp)
    deallocate(omega_tmp)

    call t_startf('derived_phys')
    if (ldry_mass_vertical_coordinates) then
      call derived_phys_dry(phys_state,phys_tend,pbuf2d)
    else
      call derived_phys_wet(phys_state,phys_tend,pbuf2d)
    end if
    call t_stopf('derived_phys')

    if (.not.ldry_mass_vertical_coordinates) then
      !
      ! phys_state%omega hold omega/p for moist vertical coordinates
      !
    !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, ilyr, icol)
      do lchnk=begchunk,endchunk
        ncols=get_ncols_p(lchnk)
        do ilyr=1,pver
          do icol=1,ncols
            phys_state(lchnk)%omega(icol,ilyr)=phys_state(lchnk)%omega(icol,ilyr)*phys_state(lchnk)%pmid(icol,ilyr)
          end do
        end do
      end do
    end if

    if (write_inithist() ) then
!      if ((fv_nphys) .and. (iam < par%nprocs)) then
!        ! Maybe write IC data
!        allocate(T_tmp(npsq,nlev,1))
!        allocate(uv_tmp(npsq,2,nlev,1))
!        allocate(ps_tmp(npsq,1))
!        elem => dyn_out%elem
!        do ie = 1, nelemd
!          do j = 1, np
!            do i = 1, np
!              T_tmp(i+(j-1)*np,:,1) = elem(ie)%state%t(i,j,:,tl_f)
!              uv_tmp(i+(j-1)*np,1,:,1) = elem(ie)%state%v(i,j,1,:,tl_f)
!              uv_tmp(i+(j-1)*np,2,:,1) = elem(ie)%state%v(i,j,2,:,tl_f)
!              ps_tmp(i+(j-1)*np,1) = elem(ie)%state%psdry(i,j,tl_f)!phl ??? *100._r8
!            end do
!          end do
!          call outfld('Tgll&IC', T_tmp(:,:,1), npsq, ie)
!          call outfld('Ugll&IC', uv_tmp(:,1,:,1), npsq, ie)
!          call outfld('Vgll&IC', uv_tmp(:,2,:,1), npsq, ie)
!          call outfld('PSgll&IC', ps_tmp(:,1), npsq, ie) !this is currently dry pressure - should be full
!        end do
!        deallocate(T_tmp)
!        deallocate(uv_tmp)
!        deallocate(ps_tmp)
!      else
      do lchnk=begchunk,endchunk
        call outfld('T&IC',phys_state(lchnk)%t,pcols,lchnk)
        call outfld('U&IC',phys_state(lchnk)%u,pcols,lchnk)
        call outfld('V&IC',phys_state(lchnk)%v,pcols,lchnk)
        call outfld('PS&IC',phys_state(lchnk)%ps,pcols,lchnk)!this is currently full pressure (see derived_phys)
        do m=1,pcnst
          call outfld(trim(cnst_name(m))//'&IC',phys_state(lchnk)%q(1,1,m), pcols,lchnk)!mixing ratios are wet for water variables
        end do
      end do
      if (ntrac>0) then
        write(*,*) "write_inithist not coded for fvm tracers yet"
!        stop
      end if
      !      end if
    endif
  end subroutine d_p_coupling

  subroutine p_d_coupling(phys_state, phys_tend,  dyn_in, tl_qdp)
    use cam_control_mod, only : adiabatic
    use dp_grids,        only: nphys_pts
    use hycoef,          only: hyai, hybi, ps0
    use cam_abortutils,  only: endrun
    use constituents,    only: cnst_type
    use dimensions_mod  ,only: fv_nphys
    use dimensions_mod  ,only: ldry_mass_vertical_coordinates
    use fvm_control_volume_mod, only: n0_fvm
    use hybrid_mod,       only: config_thread_region, get_loop_ranges
    use hybrid_mod,       only: hybrid_t
    use fvm_mapping     , only: phys2dyn_forcings_fvm
#ifdef debug_coupling
    use constituents,   only: cnst_get_ind
    use dyn_inic_baroclinic, only: test_func
#endif

    use thread_mod,       only: horz_num_threads
    implicit none

    ! !INPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend


    ! !OUTPUT PARAMETERS:
    type(dyn_import_t),  intent(inout)   :: dyn_in
    type(hybrid_t)                       :: hybrid

    integer, intent(in) :: tl_qdp

    ! LOCAL VARIABLES
    integer :: ic , ncols                            ! index
    type(element_t), pointer :: elem(:)               ! pointer to dyn_in element array
    integer (kind=int_kind)  :: ie, iep               ! indices over elements
    integer (kind=int_kind)  :: lchnk, icol, ilyr      ! indices over chunks, columns, layers

    real (kind=r8), allocatable :: dp_phys(:,:,:)     ! temporary array to hold dp on physics grid
    real (kind=r8), allocatable :: T_tmp(:,:,:)     ! temporary array to hold T
    real (kind=r8), allocatable :: dq_tmp(:,:,:,:)     ! temporary array to hold q
    real (kind=r8), allocatable :: uv_tmp(:,:,:,:)  ! temporary array to hold uv
    !   real (kind=r8), allocatable :: omega_tmp(:,:,:) ! temporary array to hold omega
    integer (kind=int_kind)  :: ioff, m, i, j, k, dtrac
    integer(kind=int_kind)   :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)

    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data
    integer, allocatable :: bpter(:,:)    ! offsets into block buffer for unpacking data

    real (kind=r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

    real (kind=r8) :: factor

    integer :: tl_f, num_trac

    integer :: nets,nete
#ifdef debug_coupling
    integer :: ixtt
    call cnst_get_ind('TT_UN',ixtt, abort=.false.)
#endif
    tl_f = TimeLevel%n0  ! time split physics (with forward-in-time RK)

    if (iam .lt. par%nprocs) then
      elem => dyn_in%elem
    else
      nullify(elem)
    end if

    allocate(T_tmp(nphys_pts,pver,nelemd))
    allocate(uv_tmp(nphys_pts,2,pver,nelemd))
    allocate(dq_tmp(nphys_pts,pver,pcnst,nelemd))
    !    allocate(omega_tmp(nphys_pts,pver,nelemd))
    allocate(dp_phys(nphys_pts,pver,nelemd))

    T_tmp=0.0_r8
    uv_tmp=0.0_r8
    dq_tmp=0.0_r8

    if (.not. allocated(q_prev)) then
      call endrun('p_d_coupling: q_prev not allocated')
    end if

    dtrac = MIN(qsize, pcnst)

    call t_startf('pd_copy')
    if(local_dp_map) then

       !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, ilyr, m, factor)
       do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)
          call get_gcol_all_p(lchnk,pcols,pgcols)

          do icol=1,ncols
             call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
             ie = idmb3(1)
             ioff=idmb2(1)

             do ilyr=1,pver
                !
                ! phys_state%ps units is Pa
                !
                if (ldry_mass_vertical_coordinates) then
                   dp_phys(ioff,ilyr,ie) = phys_state(lchnk)%pdeldry(icol,ilyr)
                else
                   dp_phys(ioff,ilyr,ie) = ((hyai(ilyr+1) - hyai(ilyr)) * ps0) +     &
                        ((hybi(ilyr+1) - hybi(ilyr)) * phys_state(lchnk)%ps(icol))
                   dp_phys(ioff,ilyr,ie) = phys_state(lchnk)%pdel(icol,ilyr)
                end if

#ifdef debug_coupling
                phys_tend(lchnk)%dtdt(icol,ilyr)   = test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,9)
                phys_tend(lchnk)%dudt(icol,ilyr)   = test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,12)
                phys_tend(lchnk)%dvdt(icol,ilyr)   = test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,13)
                q_prev(icol, ilyr, 2:pcnst, lchnk) = 0.0D0
                do m=2,pcnst
                   phys_state(lchnk)%Q(icol,ilyr,m)=test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,m)
                end do
#endif


                T_tmp(ioff,ilyr,ie)      = phys_tend(lchnk)%dtdt(icol,ilyr)
                uv_tmp(ioff,1,ilyr,ie)   = phys_tend(lchnk)%dudt(icol,ilyr)
                uv_tmp(ioff,2,ilyr,ie)   = phys_tend(lchnk)%dvdt(icol,ilyr)
                !
                ! convert wet mixing ratios to dry
                !
                if (ldry_mass_vertical_coordinates) then
                  !
                  ! Note that the outcommented formula for "factor" is wrong since physics has updated q
                  !
                  ! factor = 1/(1-SUM(phys_state(lchnk)%q(icol,ilyr,1:qsize_condensate_loading)))
                  !
                  ! this code is equivalent to dme_adjust
                  !
                   factor = phys_state(lchnk)%pdel(icol,ilyr)/phys_state(lchnk)%pdeldry(icol,ilyr)
                else
                   factor = 1.0_r8
                end if
#ifdef debug_coupling
                factor = 1.0_r8
#endif
                do m=1,pcnst
                   if (cnst_type(m).eq.'wet') phys_state(lchnk)%q(icol,ilyr,m) = factor*phys_state(lchnk)%q(icol,ilyr,m)
                   dq_tmp(ioff,ilyr,m,ie) = (phys_state(lchnk)%q(icol,ilyr,m) - q_prev(icol, ilyr, m, lchnk))!*dp_phys(ioff,ilyr,ie)
                end do
             end do
          end do
       end do

    else

       tsize = 4 + pcnst

       allocate( bbuffer(tsize*block_buf_nrecs) )
       allocate( cbuffer(tsize*chunk_buf_nrecs) )

       !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, i, icol, ilyr, m, factor)
       do lchnk = begchunk,endchunk
          ncols = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do i=1,ncols
             cbuffer(cpter(i,0):cpter(i,0)+2+pcnst) = 0.0_r8
          end do

          do icol=1,ncols
             do ilyr=1,pver
#ifdef debug_coupling
                phys_tend(lchnk)%dtdt(icol,ilyr)   = test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,9)
                phys_tend(lchnk)%dudt(icol,ilyr)   = test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,12)
                phys_tend(lchnk)%dvdt(icol,ilyr)   = test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,13)
                q_prev(icol, ilyr, 2:pcnst, lchnk) = 0.0D0
                do m=2,pcnst
                   phys_state(lchnk)%Q(icol,ilyr,m)=test_func(phys_state(lchnk)%lat(icol),phys_state(lchnk)%lon(icol),ilyr,m)
                end do
#endif

                cbuffer   (cpter(icol,ilyr))     = phys_tend(lchnk)%dtdt(icol,ilyr)
                cbuffer   (cpter(icol,ilyr)+1)   = phys_tend(lchnk)%dudt(icol,ilyr)
                cbuffer   (cpter(icol,ilyr)+2)   = phys_tend(lchnk)%dvdt(icol,ilyr)

                if (ldry_mass_vertical_coordinates) then
                   cbuffer   (cpter(icol,ilyr)+3) = phys_state(lchnk)%pdeldry(icol,ilyr)
                else
                   cbuffer   (cpter(icol,ilyr)+3) = phys_state(lchnk)%pdel(icol,ilyr)
                end if

                if (ldry_mass_vertical_coordinates) then
                  !
                  ! Note that the outcommented formula for "factor" is wrong since physics has updated q
                  !
                  ! factor = 1/(1-SUM(phys_state(lchnk)%q(icol,ilyr,1:qsize_condensate_loading)))
                  !
                  ! this code is equivalent to dme_adjust
                  !
                   factor = phys_state(lchnk)%pdel(icol,ilyr)/phys_state(lchnk)%pdeldry(icol,ilyr)
                else
                   factor = 1.0_r8
                end if
#ifdef debug_coupling
                factor=1.0_r8
#endif
                do m=1,pcnst
                  if (cnst_type(m).eq.'wet') &
                    phys_state(lchnk)%q(icol,ilyr,m) = factor*phys_state(lchnk)%q(icol,ilyr,m)

                    cbuffer(cpter(icol,ilyr)+3+m) =  (&!(cbuffer(cpter(icol,ilyr)+3)*(&
                         phys_state(lchnk)%q(icol,ilyr,m) - q_prev(icol, ilyr, m, lchnk))
                end do
             end do

          end do
       end do

       call t_barrierf('sync_chk_to_blk', mpicom)
       call t_startf ('chunk_to_block')
       call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
       call t_stopf  ('chunk_to_block')

       if(iam < par%nprocs) then
          if (fv_nphys>0) then
             allocate(bpter(fv_nphys*fv_nphys,0:pver))
          else
             allocate(bpter(npsq,0:pver))
          end if
          !$omp parallel do num_threads(horz_num_threads) private (ie, bpter, icol, ilyr, m, ncols)
          do ie=1,nelemd

             if (fv_nphys>0) then
                call chunk_to_block_recv_pters(elem(ie)%GlobalID, fv_nphys*fv_nphys,  &
                     pver+1, tsize, bpter)
                ncols = fv_nphys*fv_nphys
             else
                call chunk_to_block_recv_pters(elem(ie)%GlobalID, npsq,       &
                     pver+1, tsize, bpter)
                ncols = elem(ie)%idxP%NumUniquePts
             end if

             do icol=1,ncols

                do ilyr=1,pver

                   T_tmp   (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr))
                   uv_tmp  (icol,1,ilyr,ie) = bbuffer(bpter(icol,ilyr)+1)
                   uv_tmp  (icol,2,ilyr,ie) = bbuffer(bpter(icol,ilyr)+2)
                   dp_phys (icol,ilyr,ie)   = bbuffer(bpter(icol,ilyr)+3)

                   do m=1,pcnst
                      dq_tmp(icol,ilyr,m,ie) = bbuffer(bpter(icol,ilyr)+3+m)
                   end do

                end do

             end do

          end do
          deallocate(bpter)
       end if
       deallocate( bbuffer )
       deallocate( cbuffer )

    end if
    call t_stopf('pd_copy')
    if(iam < par%nprocs) then
       if (fv_nphys>0) then
          !
          ! put forcings into fvm structure
          !
          num_trac = max(qsize,ntrac)
          do ie = 1, nelemd
             dyn_in%fvm(ie)%ft(1:fv_nphys,1:fv_nphys,:)           = RESHAPE(T_tmp(:,:,ie)   ,(/fv_nphys,fv_nphys,pver/))
             dyn_in%fvm(ie)%fm(1:fv_nphys,1:fv_nphys,:,:)         = RESHAPE(uv_tmp(:,:,:,ie),(/fv_nphys,fv_nphys,2,pver/))
             dyn_in%fvm(ie)%fc_phys(1:fv_nphys,1:fv_nphys,:,1:num_trac)= &
                  RESHAPE(dq_tmp(:,:,1:num_trac,ie),(/fv_nphys,fv_nphys,pver,num_trac/))
             dyn_in%fvm(ie)%dp_phys(1:fv_nphys,1:fv_nphys,:)      = RESHAPE(dp_phys(:,:,ie),(/fv_nphys,fv_nphys,pver/))
          end do

          !JMD $OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n)
          !JMD        hybrid = config_thread_region(par,'horizontal')
          hybrid = config_thread_region(par,'serial')
          call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
        !
          ! high-order mapping of ft and fm (and fq if no cslam) using fvm technology
          !
          call phys2dyn_forcings_fvm(elem, dyn_in%fvm, hybrid,nets,nete,ntrac==0,tl_f,tl_qdp)
       else
          call t_startf('putUniquePoints')
          !$omp parallel do num_threads(horz_num_threads) private(ie,ncols)
          do ie=1,nelemd
             ncols = elem(ie)%idxP%NumUniquePts
             call putUniquePoints(elem(ie)%idxP, nlev, T_tmp(1:ncols,:,ie), elem(ie)%derived%fT(:,:,:,1))
             call putUniquePoints(elem(ie)%idxV, 2, nlev, uv_tmp(1:ncols,:,:,ie), &
                  elem(ie)%derived%fM(:,:,:,:,1))
             call putUniquePoints(elem(ie)%idxV, nlev,pcnst, dq_tmp(1:ncols,:,:,ie), &
                  elem(ie)%derived%fQ(:,:,:,:,1))
          end do
          call t_stopf('putUniquePoints')
       end if
    end if
    ! Deallocate the temporary arrays
    deallocate(T_tmp)
    deallocate(uv_tmp)
    deallocate(dq_tmp)
    !    deallocate(omega_tmp)
  end subroutine p_d_coupling

  subroutine derived_phys_dry(phys_state, phys_tend, pbuf2d)
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use constituents,  only: qmin
    use physconst,     only: cpair, gravit, rair, zvir, cappa, rairv,rh2o,rair
    use spmd_utils,    only: masterproc
    use ppgrid,        only: pver
    use geopotential,  only: geopotential_t
    use physics_types, only: set_state_pdry, set_wet_to_dry
    use check_energy,  only: check_energy_timestep_init
    use hycoef,   only : hyam, hybm, hyai, hybi, ps0
    use shr_vmath_mod, only: shr_vmath_log
    use phys_gmean,    only: gmean
    use constituents,  only: cnst_type
    use dimensions_mod,only: qsize_condensate_loading
    use thread_mod, only: horz_num_threads

    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),      pointer     :: pbuf2d(:,:)

    integer :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1

    real(r8) :: ke(pcols,begchunk:endchunk)
    real(r8) :: se(pcols,begchunk:endchunk)
    real(r8) :: ke_glob(1),se_glob(1)
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
    real(r8) :: factor_array(pcols,nlev)

    integer :: m, i, k, ncol

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    !
    ! Evaluate derived quantities
    !
    !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncol, k, i, m , zvirv, pbuf_chnk, factor_array)
    !
    ! dry pressure variables
    !
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do i=1,ncol
        phys_state(lchnk)%psdry(i)=hyai(1)*ps0+sum(phys_state(lchnk)%pdel(i,:))
      end do
      do i=1,ncol
        phys_state(lchnk)%pintdry(i,1)=hyai(1)*ps0
      end do
      call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,1),phys_state(lchnk)%lnpintdry(1:ncol,1),ncol)
      do k=1,nlev
        do i=1,ncol
          phys_state(lchnk)%pintdry(i,k+1)=phys_state(lchnk)%pintdry(i,k)+&
               phys_state(lchnk)%pdel(i,k)
        end do
        call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,k+1),phys_state(lchnk)%lnpintdry(1:ncol,k+1),ncol)
      end do

      do k=1,nlev
        do i=1,ncol
          phys_state(lchnk)%pdeldry (i,k) = phys_state(lchnk)%pdel(i,k)
          phys_state(lchnk)%rpdeldry(i,k) = 1._r8/phys_state(lchnk)%pdeldry(i,k)
          phys_state(lchnk)%pmiddry (i,k) = 0.5D0*(phys_state(lchnk)%pintdry(i,k+1) + phys_state(lchnk)%pintdry(i,k))
        end do
        call shr_vmath_log(phys_state(lchnk)%pmiddry(1:ncol,k),phys_state(lchnk)%lnpmiddry(1:ncol,k),ncol)
      end do

      !
      ! wet pressure variables (should be removed from physics!)
      !
      do k=1,nlev
        do i=1,ncol
!          factor_array(i,k) = 1+SUM(phys_state(lchnk)%q(i,k,1:qsize_condensate_loading))
          !
          ! to be consistent with total energy formula in physic's check_energy module only
          ! include water vapor in in moist dp
          !
          factor_array(i,k) = 1+phys_state(lchnk)%q(i,k,1)
        end do
      end do

      do k=1,nlev
        do i=1,ncol
          phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pdeldry(i,k)*factor_array(i,k)
        end do
      end do
      !
      ! initialize vertical loop - model top pressure
      !
      do i=1,ncol
        phys_state(lchnk)%ps(i)     = phys_state(lchnk)%pintdry(i,1)
        phys_state(lchnk)%pint(i,1) = phys_state(lchnk)%pintdry(i,1)
      end do
      do k = 1, nlev
        do i=1,ncol
          phys_state(lchnk)%pint(i,k+1) =  phys_state(lchnk)%pint(i,k)+phys_state(lchnk)%pdel(i,k)
          phys_state(lchnk)%pmid(i,k)   = (phys_state(lchnk)%pint(i,k+1)+phys_state(lchnk)%pint(i,k))/2._r8
          phys_state(lchnk)%ps  (i)     =  phys_state(lchnk)%ps(i) + phys_state(lchnk)%pdel(i,k)
        end do
        call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
        call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
      end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

      do k=1,nlev
        do i=1,ncol
          phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
          phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
               / phys_state(lchnk)%pmid(i,k))**cappa
        end do
      end do
      !
      ! all tracers (including moisture) are in dry mixing ratio units - physics expect water variables moist
      !
      factor_array(1:ncol,1:nlev) = 1/factor_array(1:ncol,1:nlev)

      do m = 1,pcnst
        if (cnst_type(m).eq.'wet') then
          do k=1,nlev
            do i=1,ncol
              phys_state(lchnk)%q(i,k,m) = factor_array(i,k)*phys_state(lchnk)%q(i,k,m)
            end do
          end do
        endif
      end do


      !-----------------------------------------------------------------------------------
      !  Need to fill zvirv 2D variables to be compatible with geopotential_t interface
      !-----------------------------------------------------------------------------------
      zvirv(:,:) = zvir
      !
      ! Compute initial geopotential heights - based on full pressure
      !
      call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
           phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
           phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk),  gravit,  zvirv       , &
           phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
          do i=1,ncol
#if FIX_TOTE
            ! general formula:  E = CV_air T + phis + gravit*zi )
            ! hydrostatic case: integrate zi term by parts, use CP=CV+R to get:
            ! E = CP_air T + phis   (Holton Section 8.3)
            ! to use this, update geopotential.F90, and other not-yet-found physics routines:
            ! (check boundary layer code, others which have gravit and zi() or zm()
            phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                 + phys_state(lchnk)%phis(i)
#else
            phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                 + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
#endif
          end do
        end do

        !
        ! Ensure tracers are all positive
        !
        call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
             1, pcnst, qmin  ,phys_state(lchnk)%q)

        ! Compute energy and water integrals of input state
        pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
        call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)


#if 0
        ke(:,lchnk) = 0._r8
        se(:,lchnk) = 0._r8
        !       wv = 0._r8
        !       wl = 0._r8
        !       wi = 0._r8
        do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + &
                  phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
             !             wv = wv + phys_state(lchnk)%q(i,k,1       )*phys_state(lchnk)%pdel(i,k)
             !             wl = wl + phys_state(lchnk)%q(i,k,ixcldliq)*phys_state(lchnk)%pdel(i,k)
             !             wi = wi + phys_state(lchnk)%q(i,k,ixcldice)*phys_state(lchnk)%pdel(i,k)
           end do
         end do
#endif
       end do!lchnk

#if 0
       ! This wont match SE exactly.  SE computes KE at half levels
       ! SE includes cp_star (physics SE uses cp )
       ! CAM stdout of total energy also includes latent energy of Q,Q1,Q2
       ! making it a little larger
       call gmean(ke,ke_glob,1)
       call gmean(se,se_glob,1)
       if (masterproc) then
         write(iulog,'(a,e20.8,a,e20.8)') 'KE = ',ke_glob(1),' SE = ',se_glob(1)
         write(iulog,'(a,e20.8)') 'TOTE = ',ke_glob(1)+se_glob(1)
       endif
#endif
     end subroutine derived_phys_dry

     subroutine derived_phys_wet(phys_state, phys_tend, pbuf2d)
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use constituents,  only: qmin
    use physconst,     only: cpair, gravit, rair, zvir, cappa, rairv
    use spmd_utils,    only: masterproc
    use ppgrid,        only: pver
    use geopotential,  only: geopotential_t
    use physics_types, only: set_state_pdry, set_wet_to_dry
    use check_energy,  only: check_energy_timestep_init
    use hycoef,   only : hyam, hybm, hyai, hybi, ps0
    use shr_vmath_mod, only: shr_vmath_log
    use phys_gmean,      only: gmean
    use thread_mod,      only: horz_num_threads


    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),      pointer     :: pbuf2d(:,:)

    integer :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1

    real(r8) :: ke(pcols,begchunk:endchunk)
    real(r8) :: se(pcols,begchunk:endchunk)
    real(r8) :: ke_glob(1),se_glob(1)
    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer

    integer :: m, i, k, ncol

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

!
! Evaluate derived quantities
!
    !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncol, k, i, zvirv, pbuf_chnk)
    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)
       do k=1,nlev
          do i=1,ncol
             phys_state(lchnk)%pint(i,k)=hyai(k)*ps0+hybi(k)*phys_state(lchnk)%ps(i)
             phys_state(lchnk)%pmid(i,k)=hyam(k)*ps0+hybm(k)*phys_state(lchnk)%ps(i)
          end do
          call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
          call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
       end do
       do i=1,ncol
          phys_state(lchnk)%pint(i,pverp)=hyai(pverp)*ps0+hybi(pverp)*phys_state(lchnk)%ps(i)
       end do
       call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)


       do k=1,nlev
          do i=1,ncol
            phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                                             / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do

!-----------------------------------------------------------------------------------
!  Need to fill zvirv 2D variables to be compatible with geopotential_t interface
!-----------------------------------------------------------------------------------
       zvirv(:,:) = zvir
!
! Compute initial geopotential heights
!

       call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
            phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
            phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk),  gravit,  zvirv       , &
            phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

! Compute initial dry static energy, include surface geopotential
       do k = 1, pver
          do i=1,ncol
#if FIX_TOTE
             ! general formula:  E = CV_air T + phis + gravit*zi )
             ! hydrostatic case: integrate zi term by parts, use CP=CV+R to get:
             ! E = CP_air T + phis   (Holton Section 8.3)
             ! to use this, update geopotential.F90, and other not-yet-found physics routines:
             ! (check boundary layer code, others which have gravit and zi() or zm()
             phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                                      + phys_state(lchnk)%phis(i)
#else
             phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                  + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
#endif
          end do
       end do

! NOTE:  if a tracer is marked "dry", that means physics wants it dry
!        if dycore advects it wet, it should be converted here
!        FV dycore does this, and in physics/cam/tphysac.F90 it will
!        be converted back to wet, BUT ONLY FOR FV dycore
!
!        EUL: advects all tracers (except q1) as dry.  so it never
!        calls this.
!
!        SE:  we should follow FV and advect all tracers wet (especially
!        since we will be switching to conservation form of advection).
!        So this is broken since dry tracers will never get converted
!        back to wet. (in APE, all tracers are wet, so it is ok for now)
!
! Convert dry type constituents from moist to dry mixing ratio
       call set_state_pdry(phys_state(lchnk))    ! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.
!
! Ensure tracers are all positive
!
       call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
                  1, pcnst, qmin  ,phys_state(lchnk)%q)

! Compute energy and water integrals of input state
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)


#if 0
       ke(:,lchnk) = 0._r8
       se(:,lchnk) = 0._r8
!       wv = 0._r8
!       wl = 0._r8
!       wi = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + &
                  phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
!             wv = wv + phys_state(lchnk)%q(i,k,1       )*phys_state(lchnk)%pdel(i,k)
!             wl = wl + phys_state(lchnk)%q(i,k,ixcldliq)*phys_state(lchnk)%pdel(i,k)
!             wi = wi + phys_state(lchnk)%q(i,k,ixcldice)*phys_state(lchnk)%pdel(i,k)
          end do
       end do
#endif
    end do

#if 0
! This wont match SE exactly.  SE computes KE at half levels
! SE includes cp_star (physics SE uses cp )
! CAM stdout of total energy also includes latent energy of Q,Q1,Q2
! making it a little larger
    call gmean(ke,ke_glob,1)
    call gmean(se,se_glob,1)
    if (masterproc) then
       write(iulog,'(a,e20.8,a,e20.8)') 'KE = ',ke_glob(1),' SE = ',se_glob(1)
       write(iulog,'(a,e20.8)') 'TOTE = ',ke_glob(1)+se_glob(1)
    endif
#endif

  end subroutine derived_phys_wet


#ifdef xxx
#ifdef debug_coupling
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
        call outfld('xUgll', ftmp(:,:,1), npsq, ie)
        call outfld('xVgll', ftmp(:,:,2), npsq, ie)
        do j = 1, np
          do i = 1, np
            do k = 1, nlev
              ftmp(i+(j-1)*np,k,1) = elem(ie)%state%T(i,j,k,tl_f)
              ftmp(i+(j-1)*np,k,2) = elem(ie)%state%dp3d(i,j,k,tl_f)
            end do
          end do
        end do
        call outfld('xTgll', ftmp(:,:,1), npsq, ie)
        call outfld('xpdeldrygll', ftmp(:,:,2), npsq, ie)
        do j = 1, np
          do i = 1, np
            ftmp(i+(j-1)*np,1,1) = elem(ie)%state%ps(i,j,tl_f)
          end do
        end do
        call outfld('xpsdrygll', ftmp(:,1,1), npsq, ie)

        do j = 1, np
          do i = 1, np
            ftmp(i+(j-1)*np,1,1) = elem(ie)%state%pswet(i,j)
          end do
        end do
        call outfld('xpswetgll', ftmp(:,1,1), npsq, ie)

        !    end if
      end do

  end subroutine dptest_write_dyn_state_history

#endif
#endif
end module dp_coupling
