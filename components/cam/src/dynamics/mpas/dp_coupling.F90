
! Uncomment to assume traditional CAM hybrid vertical coordinate and
!  that dynamics surface pressure is total rather than dry

!
!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------
module dp_coupling
  use cam_abortutils, only: endrun
  use constituents,   only: pcnst, cnst_name
  use cam_history,    only: outfld, write_inithist 
  use dyn_comp,       only: dyn_export_t, dyn_import_t
  use dyn_grid,       only: get_gcol_block_d, max_col_per_block, nblocks_per_pe, num_col_per_block
  use dyn_grid,       only: col_indices_in_block
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use physics_types,  only: physics_state, physics_tend
  use phys_grid
  use dyn_grid,       only: global_to_local_cell
  use pmgrid
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
  use cam_abortutils, only: endrun
  use cam_logfile,    only : iulog
  use spmd_utils,     only: iam
#ifdef SPMD
  use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use mpishorthand,   only: mpicom
#else
  use spmd_dyn,       only: local_dp_map
#endif
  use perf_mod
  use mpas_cam_interface,  only: mpas_to_cam, cam_to_mpas
  private
  public :: d_p_coupling, p_d_coupling
!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
  subroutine d_p_coupling(phys_state, phys_tend,  pbuf2d, dyn_out)

    use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field
    use shr_vmath_mod,  only: shr_vmath_exp
    use time_manager,   only: is_first_step
    use cam_abortutils, only: endrun

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

! Variables from dynamics export container
    real(r8), pointer :: phis(:)
    real(r8), pointer :: psd(:)
    real(r8), pointer :: pint(:,:)
    real(r8), pointer :: pmid(:,:)
    real(r8), pointer :: zint(:,:)
    real(r8), pointer :: zmid(:,:)
    real(r8), pointer :: ux(:,:)
    real(r8), pointer :: uy(:,:)
    real(r8), pointer :: temp(:,:)
    real(r8), pointer :: omega(:,:)
    real(r8), pointer :: tracer(:,:,:)
!++KSA
    real(r8), pointer :: div(:,:)
    real(r8), pointer :: vor(:,:)
!--KSA

    integer :: lchnk, icol, k      ! indices over chunks, columns, layers
    integer :: ncols, ig, nblk, nb, m, i

    integer :: pgcols(pcols)
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: bpter(max_col_per_block,0:pver)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

#if (! defined SPMD)
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
#endif

    phis => dyn_out%phis
    psd => dyn_out%psd
    pint => dyn_out%pint
    pmid => dyn_out%pmid
    zint => dyn_out%zint
    zmid => dyn_out%zmid
    ux => dyn_out%ux
    uy => dyn_out%uy
    temp => dyn_out%t
    omega => dyn_out%omega
    tracer => dyn_out%tracer
!++KSA
    div => dyn_out%div
    vor => dyn_out%vor
!--KSA

! Note: Omega is dimensioned plev+1, but this interface transmits only the first plev levels.
    call mpas_to_cam(numcols, pver, pcnst, psd, phis, pint(:,pverp:1:-1), pmid(:,pver:1:-1), &
                     zint(:,pverp:1:-1), zmid(:,pver:1:-1), &
                     temp(:,plev:1:-1), ux(:,plev:1:-1), uy(:,plev:1:-1), &
                     omega(:,plev:1:-1), tracer(:,plev:1:-1,:), div(:,plev:1:-1),vor(:,plev:1:-1)) 
    !KSA div(:,plev:1:-1),vor(:,plev:1:-1) added
    omega(:,plev+1) = 0._r8


    call t_startf('dpcopy')
    if (local_dp_map) then

!$omp parallel do private (lchnk, ncols, icol, i, k, m, pgcols)
       do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)                         ! number of columns in this chunk
          call get_gcol_all_p(lchnk,pcols,pgcols)          ! global column indices
          do icol=1,ncols
             i = global_to_local_cell(pgcols(icol))        ! local (to process) column index

             !phys_state(lchnk)%ps(icol)=psd(i)
             phys_state(lchnk)%psdry(icol)=psd(i)

             phys_state(lchnk)%phis(icol)=phis(i)
             do k=1,pver
                phys_state(lchnk)%t(icol,k)=temp(i,k)	   
                phys_state(lchnk)%u(icol,k)=ux(i,k)
                phys_state(lchnk)%v(icol,k)=uy(i,k)
                phys_state(lchnk)%omega(icol,k)=omega(i,k)
                phys_state(lchnk)%pintdry(icol,k)=pint(i,k)
                phys_state(lchnk)%pmiddry(icol,k)=pmid(i,k)
                !phys_state(lchnk)%pint(icol,k)=pint(i,k)
                !phys_state(lchnk)%pmid(icol,k)=pmid(i,k)
                phys_state(lchnk)%zi(icol,k)=zint(i,k)
                phys_state(lchnk)%zm(icol,k)=zmid(i,k)
!++KSA
                phys_state(lchnk)%div(icol,k)=div(i,k)
                phys_state(lchnk)%vor(icol,k)=vor(i,k)
!--KSA
             end do
             phys_state(lchnk)%pintdry(icol,pverp)=pint(i,pverp)
             !phys_state(lchnk)%pint(icol,pverp)=pint(i,pverp)
             phys_state(lchnk)%zi(icol,pverp)=zint(i,pverp)

             do m=1,pcnst
                do k=1,pver
                   phys_state(lchnk)%q(icol,k,m)=tracer(i,k,m)
                end do
             end do
          end do
       end do


    else  ! .not. local_dp_map

! MGD what are we to do here?

       tsize = 4 + pcnst
       allocate(bbuffer(tsize*block_buf_nrecs))
       allocate(cbuffer(tsize*chunk_buf_nrecs))

!$omp parallel do private (nb, nblk, ncols, icol, ig, i, k, m, bpter)
       do nb = 1, nblocks_per_pe
          nblk = iam * nblocks_per_pe + nb   !  global block index
          ncols = num_col_per_block(nblk)    !  number of columns in this block

          call block_to_chunk_send_pters(nblk,max_col_per_block,pver+1,tsize,bpter)

          do icol=1,ncols
             ig = col_indices_in_block(icol,nblk)   !  global column index
             i = global_to_local_cell(ig)           !  local (to process) column index

             bbuffer(bpter(icol,0)+2:bpter(icol,0)+3+pcnst) = 0.0_r8

             bbuffer(bpter(icol,0))   = psd(i)
             bbuffer(bpter(icol,0)+1) = phis(i)

             do k=1,pver
                bbuffer(bpter(icol,k))   = temp(i,k)
                bbuffer(bpter(icol,k)+1) = ux(i,k)
                bbuffer(bpter(icol,k)+2) = uy(i,k)
                bbuffer(bpter(icol,k)+3) = omega(i,k)

                do m=1,pcnst
                   bbuffer(bpter(icol,k)+3+m) = tracer(i,k,m)
                end do
             end do

          end do

       end do

       call t_barrierf ('sync_blk_to_chk', mpicom)
       call t_startf ('block_to_chunk')
       call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
       call t_stopf  ('block_to_chunk')

!$omp parallel do private (lchnk, ncols, icol, k, m, cpter)
       do lchnk = begchunk,endchunk
          ncols = phys_state(lchnk)%ncol

          call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

          do icol=1,ncols

             !phys_state(lchnk)%ps   (icol)     = cbuffer(cpter(icol,0))
             phys_state(lchnk)%psdry (icol)     = cbuffer(cpter(icol,0))

             phys_state(lchnk)%phis (icol)     = cbuffer(cpter(icol,0)+1)

             do k=1,pver

                phys_state(lchnk)%t     (icol,k)   = cbuffer(cpter(icol,k))
                phys_state(lchnk)%u     (icol,k)   = cbuffer(cpter(icol,k)+1)
                phys_state(lchnk)%v     (icol,k)   = cbuffer(cpter(icol,k)+2)
                phys_state(lchnk)%omega (icol,k)   = cbuffer(cpter(icol,k)+3)

                do m=1,pcnst
                   phys_state(lchnk)%q  (icol,k,m) = cbuffer(cpter(icol,k)+3+m)
                end do

             end do

          end do

       end do

       deallocate( bbuffer )
       deallocate( cbuffer )

    end if
    call t_stopf('dpcopy')

    call t_startf('derived_phys')
    call derived_phys(phys_state,phys_tend,pbuf2d)
    call t_stopf('derived_phys')

!++KSA
   do lchnk=begchunk,endchunk
        call outfld('DIV',phys_state(lchnk)%div,pcols,lchnk)
        call outfld('VOR',phys_state(lchnk)%vor,pcols,lchnk)
   end do
!--KSA

    if (write_inithist() ) then
     do lchnk=begchunk,endchunk
        call outfld('T&IC',phys_state(lchnk)%t,pcols,lchnk)
        call outfld('U&IC',phys_state(lchnk)%u,pcols,lchnk)
        call outfld('V&IC',phys_state(lchnk)%v,pcols,lchnk)
        call outfld('PS&IC',phys_state(lchnk)%ps,pcols,lchnk)
!jjang        call outfld('PHIS&IC',phys_state(lchnk)%phis,pcols,lchnk)
        do m=1,pcnst
           call outfld(trim(cnst_name(m))//'&IC',phys_state(lchnk)%q(1,1,m), pcols,lchnk)
        end do
     end do
    endif
       
  end subroutine d_p_coupling

  subroutine p_d_coupling(phys_state, phys_tend,  dyn_in)

    use shr_vmath_mod, only: shr_vmath_log
    use cam_control_mod, only : adiabatic
    use cam_history, only : outfld

    implicit none

! !INPUT PARAMETERS:
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ),   intent(inout), dimension(begchunk:endchunk) :: phys_tend


! !OUTPUT PARAMETERS:
    type(dyn_import_t),  intent(inout)   :: dyn_in

    ! LOCAL VARIABLES
    integer :: ic , ncols                            ! index
    integer :: lchnk, icol, k      ! indices over chunks, columns, layers

    integer :: m, i, j, nb, nblk, ig
    integer :: pgcols(pcols)

    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data
    integer :: bpter(max_col_per_block,0:pver)    ! offsets into block buffer for unpacking data

    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

#if (! defined SPMD)
    integer  :: block_buf_nrecs = 0
    integer  :: chunk_buf_nrecs = 0
    integer  :: mpicom = 0
#endif

    real(r8), pointer :: temp(:,:)
    real(r8), pointer :: tracer(:,:,:)
    real(r8), pointer :: ux_tend(:,:)
    real(r8), pointer :: uy_tend(:,:)
    real(r8), pointer :: temp_tend(:,:)

    temp => dyn_in%t
    tracer => dyn_in%tracer
    ux_tend => dyn_in%ux_tend
    uy_tend => dyn_in%uy_tend
    temp_tend => dyn_in%t_tend

!   do lchnk=begchunk,endchunk
!       call outfld('FU',phys_tend(lchnk)%dudt,pcols,lchnk)
!       call outfld('FV',phys_tend(lchnk)%dvdt,pcols,lchnk)
!       call outfld('FT',phys_tend(lchnk)%dtdt,pcols,lchnk)
!   end do

    call t_startf('pd_copy')
    if(local_dp_map) then

!$omp parallel do private (lchnk, ncols, icol, i, k, m, pgcols)
       do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)                         ! number of columns in this chunk
          call get_gcol_all_p(lchnk,pcols,pgcols)          ! global column indices
          do icol=1,ncols
             i = global_to_local_cell(pgcols(icol))        ! local (to process) column index

             do k=1,pver
                temp_tend(i,k) = phys_tend(lchnk)%dtdt(icol,k)
                ux_tend(i,k)   = phys_tend(lchnk)%dudt(icol,k)
                uy_tend(i,k)   = phys_tend(lchnk)%dvdt(icol,k)
                do m=1,pcnst
                   tracer(i,k,m) = phys_state(lchnk)%q(icol,k,m)
                end do
             end do

          end do

       end do

    else

       tsize = 3 + pcnst

       allocate( bbuffer(tsize*block_buf_nrecs) )
       allocate( cbuffer(tsize*chunk_buf_nrecs) )

!$omp parallel do private (lchnk, ncols, icol, k, m, cpter)
       do lchnk = begchunk,endchunk
          ncols = get_ncols_p(lchnk)

          call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

          do icol=1,ncols
             cbuffer(cpter(icol,0):cpter(icol,0)+2+pcnst) = 0.0_r8
          end do

          do icol=1,ncols

             do k=1,pver
                cbuffer   (cpter(icol,k))     = phys_tend(lchnk)%dtdt(icol,k)
                cbuffer   (cpter(icol,k)+1)   = phys_tend(lchnk)%dudt(icol,k)
                cbuffer   (cpter(icol,k)+2)   = phys_tend(lchnk)%dvdt(icol,k)

                do m=1,pcnst
                   cbuffer(cpter(icol,k)+2+m) = phys_state(lchnk)%q(icol,k,m)
                end do
             end do

          end do

       end do

       call t_barrierf('sync_chk_to_blk', mpicom)
       call t_startf ('chunk_to_block')
       call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
       call t_stopf  ('chunk_to_block')

!$omp parallel do private (nb, nblk, ncols, icol, ig, i, k, m, bpter)
       do nb = 1, nblocks_per_pe
          nblk = iam * nblocks_per_pe + nb   !  global block index
          ncols = num_col_per_block(nblk)    !  number of columns in this block

          call chunk_to_block_recv_pters(nblk,max_col_per_block,pver+1,tsize,bpter)

          do icol=1,ncols
             ig = col_indices_in_block(icol,nblk)   !  global column index
             i = global_to_local_cell(ig)           !  local (to process) column index

             do k=1,pver

                temp_tend(i,k) = bbuffer(bpter(icol,k))
                ux_tend  (i,k) = bbuffer(bpter(icol,k)+1)
                uy_tend  (i,k) = bbuffer(bpter(icol,k)+2)

                do m=1,pcnst
                   tracer(i,k,m) = bbuffer(bpter(icol,k)+2+m)
                end do

             end do

          end do

       end do

       deallocate( bbuffer )
       deallocate( cbuffer )

    end if

    call cam_to_mpas(numcols, pver, pcnst, &
                     temp_tend(:,plev:1:-1), ux_tend(:,plev:1:-1),   &
                     uy_tend(:,plev:1:-1), tracer(:,plev:1:-1,:))

    call t_stopf('pd_copy')

  end subroutine p_d_coupling

  subroutine derived_phys(phys_state, phys_tend, pbuf2d)
    use constituents,  only: qmin
    use physconst,     only: cpair, gravit, rair, zvir, cappa
    use spmd_utils,    only: masterproc
    use ppgrid,        only: pver
    use physics_types, only: set_state_pdry, set_wet_to_dry
    use check_energy,  only: check_energy_timestep_init
    use shr_vmath_mod, only: shr_vmath_log
    use infnan,        only: inf
    use phys_gmean,      only: gmean
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk


    implicit none
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ),   intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: lchnk
    real(r8) :: qbot                 ! bottom level q before change
    real(r8) :: qbotm1               ! bottom-1 level q before change
    real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
    real(r8) :: qmavl                ! available q at level pver-1

    real(r8) :: ke(pcols,begchunk:endchunk)   
    real(r8) :: se(pcols,begchunk:endchunk)   
    real(r8) :: ke_glob(1),se_glob(1)

    integer :: m, i, k, ncol
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    real(r8):: p00

    p00   = 1.e5_r8
!
! Evaluate derived quantities
!
!$omp parallel do private (lchnk, ncol, k, i)
    do lchnk = begchunk,endchunk
       ncol = get_ncols_p(lchnk)

!       do k=1,pver
!          call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
!          call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
!       end do
!       call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)
!
!       do k=1,pver
!          do i=1,ncol
!             phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
!             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
!             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
!                                             / phys_state(lchnk)%pmid(i,k))**cappa
!          end do
!       end do

       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%pdeldry(i,k) = phys_state(lchnk)%pintdry(i,k+1) - phys_state(lchnk)%pintdry(i,k)
          end do
       end do
       phys_state(lchnk)%rpdeldry(:ncol,:) = 1._r8/phys_state(lchnk)%pdeldry(:ncol,:)
       phys_state(lchnk)%lnpmiddry(:ncol,:) = log(phys_state(lchnk)%pmiddry(:ncol,:))
       phys_state(lchnk)%lnpintdry(:ncol,:) = log(phys_state(lchnk)%pintdry(:ncol,:))

       do k=1,pver
          phys_state(lchnk)%pdel(:ncol,k) = phys_state(lchnk)%pdeldry(:ncol,k)/(1._r8-phys_state(lchnk)%q(:ncol,k,1))

          !SHP-PS
          !phys_state(lchnk)%pdel(:ncol,k) = phys_state(lchnk)%pdeldry(:ncol,k)*(1._r8 + &
          !                                  phys_state(lchnk)%q(:ncol,k,1)/(rair/rh2o*(1._r8-phys_state(lchnk)%q(:ncol,k,1))))
       end do

       phys_state(lchnk)%ps(:ncol) = phys_state(lchnk)%pintdry(:ncol,1)
       phys_state(lchnk)%pint(:ncol,1) = phys_state(lchnk)%pintdry(:ncol,1)
       do k=1,pver
          phys_state(lchnk)%pint(:ncol,k+1) = phys_state(lchnk)%pint(:ncol,k)+phys_state(lchnk)%pdel(:ncol,k)
          phys_state(lchnk)%pmid(:ncol,k) = (phys_state(lchnk)%pint(:ncol,k+1)+phys_state(lchnk)%pint(:ncol,k))/2._r8
          phys_state(lchnk)%ps(:ncol) = phys_state(lchnk)%ps(:ncol) + phys_state(lchnk)%pdel(:ncol,k)
       end do

       do k=1,pver
          call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
          call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
       end do
       call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

       do k=1,pver
          do i=1,ncol
             phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
!             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) / phys_state(lchnk)%pmid(i,k))**cappa 
             phys_state(lchnk)%exner (i,k) = (p00 / phys_state(lchnk)%pmid(i,k))**cappa
          end do
       end do


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
!        HOMME:  we should follow FV and advect all tracers wet (especially
!        since we will be switching to conservation form of advection).  
!        So this is broken since dry tracers will never get converted 
!        back to wet. (in APE, all tracers are wet, so it is ok for now)  
!
!        MPAS: advects dry mixing ratios, but coupler did converted whole dry mixing ratio to wet mixing ratio
!

! Convert dry type constituents from moist to dry mixing ratio
       call set_state_pdry(phys_state(lchnk))	 ! First get dry pressure to use for this timestep
       call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.

! Compute energy and water integrals of input state
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

	
#if 0
       ke(:,lchnk) = 0._r8
       se(:,lchnk) = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
             se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
          end do
       end do
#endif

    end do

  end subroutine derived_phys


end module dp_coupling
