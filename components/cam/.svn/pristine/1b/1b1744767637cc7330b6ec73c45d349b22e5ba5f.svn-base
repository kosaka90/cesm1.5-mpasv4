module fvm_mapping
  use kinds                 , only: real_kind, int_kind
  use dimensions_mod        , only: ngpc, irecons_tracer
  use element_mod           , only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use perf_mod             , only : t_startf, t_stopf ! _EXTERNAL
  use cam_abortutils,         only: endrun

  implicit none
  private

  logical :: ltmp !xxx

  save

  real (kind=real_kind), dimension(ngpc), private :: gsweights, gspts
  real (kind=real_kind),parameter       , private :: eps=2.0e-14_real_kind

  public :: phys2dyn_forcings_fvm, dyn2phys, dyn2phys_vector, dyn2phys_all_vars,dyn2fvm_mass_vars
contains


  !
  ! map all mass variables from gll to fvm
  !
  subroutine phys2dyn_forcings_fvm(elem, fvm, hybrid,nets,nete,no_cslam,tl_f,tl_qdp)
    use dimensions_mod, only : np, nc,nlev, nhe, qsize, nht, phys2dyn_interp_method
    use dimensions_mod, only : qsize_condensate_loading,fv_nphys, nhc_phys,ntrac,nhc
    use fvm_control_volume_mod, only: n0_fvm
    use hybrid_mod, only : hybrid_t
    use cam_abortutils,         only: endrun
#ifdef debug_coupling
   use dyn_inic_baroclinic, only: test_func
#endif


    type (element_t), intent(inout):: elem(:)
    type(fvm_struct), intent(inout):: fvm(:)

    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    logical, intent(in) :: no_cslam
    integer, intent(in) :: tl_f,tl_qdp

    integer, intent(in)                     :: nets, nete

    integer                           :: i, j, ie, k, qsize_local, m_cnst

    real (kind=real_kind), dimension(:,:,:,:,:)  , allocatable :: fld_phys, fld_gll, fld_fvm
    real (kind=real_kind), dimension(:,:,:,:)    , allocatable :: dp_fvm_tmp
    real (kind=real_kind), dimension(:,:,:)      , allocatable :: inv_dp_phys
    !
    ! for tensor product Lagrange interpolation
    !
    integer :: iwidth, nflds
    real (kind=real_kind):: v(2)
    logical, allocatable :: llimiter(:)

    if (no_cslam) then
      !
      !*************************************
      !
      ! no cslam case
      !
      !*************************************
      !
      nflds = qsize+3
      allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,nflds,nets:nete))
      allocate(fld_gll(np,np,nlev,nflds,nets:nete))
      allocate(llimiter(nflds))
      llimiter          = .false.
      llimiter(4:nflds) = .true.
      do ie=nets,nete
        !
        ! pack fields that need to be interpolated
        !
        fld_phys(1:nhc_phys,1:fv_nphys,:,1,ie)       = fvm(ie)%ft(1:fv_nphys,1:fv_nphys,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,2,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,1,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,3,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,2,:)
        !
        ! convert forcing from mass per unit area to mixing ratio per unit area
        !
        do m_cnst=1,qsize
          do k=1,nlev
            fld_phys(1:fv_nphys,1:fv_nphys,k,m_cnst+3,ie) = &
                 fvm(ie)%fc_phys(1:fv_nphys,1:fv_nphys,k,m_cnst)!/fvm(ie)%dp_phys(1:fv_nphys,1:fv_nphys,k)
          end do
        end do
      end do
      !
      ! do mapping
      !
      call fill_halo_phys(elem,fld_phys,fld_gll,hybrid,nets,nete,nflds,fvm)
      call phys2dyn(elem,fld_phys,fld_gll,nets,nete,nflds,fvm,llimiter,2)
      do ie=nets,nete
        elem(ie)%derived%fT(:,:,:,1)   = fld_gll(:,:,:,1,ie)
        elem(ie)%derived%fM(:,:,1,:,1) = fld_gll(:,:,:,2,ie)
        elem(ie)%derived%fM(:,:,2,:,1) = fld_gll(:,:,:,3,ie)
      end do
      do ie=nets,nete
        do m_cnst=1,qsize
          !
          ! convert to mass tendency
          !
          elem(ie)%derived%fq(:,:,:,m_cnst,1) = fld_gll(:,:,:,m_cnst+3,ie)!*elem(ie)%state%dp3d(:,:,:,tl_f)
        end do
      end do
    else if (nc.ne.fv_nphys) then
      !
      !***********************************************************
      !
      ! using cslam and different resolution physics grid
      !
      !***********************************************************
      !
      nflds = 3+ntrac+1
      allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,nflds,nets:nete))
      allocate(fld_gll(np,np,nlev,3,nets:nete))
      allocate(llimiter(nflds))
      allocate(inv_dp_phys(fv_nphys,fv_nphys,nlev))

      llimiter          = .false.
      do ie=nets,nete
         inv_dp_phys = 1.0D0/fvm(ie)%dp_phys(1:fv_nphys,1:fv_nphys,:)
         !
         ! pack fields that need to be interpolated
         !
         fld_phys(1:nhc_phys,1:fv_nphys,:,1,ie)       = fvm(ie)%ft(1:fv_nphys,1:fv_nphys,:)
         fld_phys(1:fv_nphys,1:fv_nphys,:,2,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,1,:)
         fld_phys(1:fv_nphys,1:fv_nphys,:,3,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,2,:)
         fld_phys(1:fv_nphys,1:fv_nphys,:,4,ie)       = fvm(ie)%dp_phys(1:fv_nphys,1:fv_nphys,:)
         do m_cnst=1,ntrac
            fld_phys(1:fv_nphys,1:fv_nphys,:,4+m_cnst,ie) = &
                 fvm(ie)%fc_phys(1:fv_nphys,1:fv_nphys,:,m_cnst)!*inv_dp_phys(:,:,:)
         end do
      end do

      call fill_halo_phys(elem,fld_phys,fld_gll,hybrid,nets,nete,nflds,fvm)
      !
      ! do mapping of fu,fv,ft
      !
      call phys2dyn(elem,fld_phys(:,:,:,1:3,:),fld_gll,nets,nete,3,fvm,llimiter(1:3),2)
      do ie=nets,nete
        elem(ie)%derived%fT(:,:,:,1)   = fld_gll(:,:,:,1,ie)
        elem(ie)%derived%fM(:,:,1,:,1) = fld_gll(:,:,:,2,ie)
        elem(ie)%derived%fM(:,:,2,:,1) = fld_gll(:,:,:,3,ie)
      end do
      !
      ! map fq from phys to fvm
      !
      allocate(dp_fvm_tmp(nc,nc,nlev,nets:nete))
      do ie=nets,nete
         do k=1,nlev
            fld_phys(:,:,k,4,ie) = 1.0D0!xxx
            call phys2fvm(ie,fvm(ie),fld_phys(:,:,k,4,ie),dp_fvm_tmp(:,:,k,ie),&
                 fld_phys(:,:,k,5:4+ntrac,ie),fvm(ie)%fc(:,:,k,1:ntrac),ntrac)
         end do
      end do
      !
      ! overwrite SE Q with cslam Q
      !
      allocate(fld_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,qsize_condensate_loading,nets:nete))
      do ie=nets,nete
        !
        ! compute cslam updated Q value
        !
        do m_cnst=1,qsize_condensate_loading
          fld_fvm(1:nc,1:nc,:,m_cnst,ie) = fvm(ie)%c(1:nc,1:nc,:,m_cnst,n0_fvm)+&
               fvm(ie)%fc(1:nc,1:nc,:,m_cnst)/fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm)
        enddo
        llimiter(1:qsize_condensate_loading) = .false.
      end do
      call fvm2dyn(elem,fld_fvm,fld_gll(:,:,:,1:qsize_condensate_loading,:),hybrid,nets,nete,&
           qsize_condensate_loading,fvm,llimiter(1:qsize_condensate_loading))
      !
      ! fld_gll now holds q cslam value on gll grid
      !
      ! convert fld_gll to increment (q_new-q_old)*dp
      !
      do ie=nets,nete
        do m_cnst=1,qsize_condensate_loading
          elem(ie)%derived%fq(:,:,:,m_cnst,1)   = &
!               fld_gll(:,:,:,m_cnst,ie)*elem(ie)%state%dp3d(:,:,:,tl_f)-&
!               elem(ie)%state%Qdp(:,:,:,m_cnst,tl_qdp)
               elem(ie)%state%dp3d(:,:,:,tl_f)*(&
               fld_gll(:,:,:,m_cnst,ie)-elem(ie)%state%Q(:,:,:,m_cnst))
        end do
      end do
      deallocate(fld_fvm)
    else
      !
      !
      !*****************************************************************************************
      !
      ! using cslam with same physics grid resolution as cslam resolution
      !
      !*****************************************************************************************
      !
      nflds = 3+qsize_condensate_loading
      allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,nflds,nets:nete))
      allocate(fld_gll(np,np,nlev,nflds,nets:nete))
      allocate(llimiter(nflds))
      llimiter          = .false.
      llimiter(4:nflds) = .true.
      do ie=nets,nete
        !
        ! pack fields that need to be interpolated
        !
        fld_phys(1:nhc_phys,1:fv_nphys,:,1,ie)       = fvm(ie)%ft(1:fv_nphys,1:fv_nphys,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,2,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,1,:)
        fld_phys(1:fv_nphys,1:fv_nphys,:,3,ie)       = fvm(ie)%fm(1:fv_nphys,1:fv_nphys,2,:)
        !
        ! compute cslam mixing ratio with physics update
        !
        do m_cnst=1,qsize_condensate_loading
          do k=1,nlev
            fld_phys(1:fv_nphys,1:fv_nphys,k,m_cnst+3,ie) = &
                 fvm(ie)%c(1:fv_nphys,1:fv_nphys,k,m_cnst,n0_fvm)+&
                 fvm(ie)%fc_phys(1:fv_nphys,1:fv_nphys,k,m_cnst)!/fvm(ie)%dp_phys(1:fv_nphys,1:fv_nphys,k)
          end do
        end do
      end do
      !
      ! do mapping
      !
      call fill_halo_phys(elem,fld_phys,fld_gll,hybrid,nets,nete,nflds,fvm)
      call phys2dyn(elem,fld_phys,fld_gll,nets,nete,nflds,fvm,llimiter,2)
      do ie=nets,nete
        elem(ie)%derived%fT(:,:,:,1)   = fld_gll(:,:,:,1,ie)
        elem(ie)%derived%fM(:,:,1,:,1) = fld_gll(:,:,:,2,ie)
        elem(ie)%derived%fM(:,:,2,:,1) = fld_gll(:,:,:,3,ie)
      end do
      do ie=nets,nete
        do m_cnst=1,qsize_condensate_loading
          !
          ! convert fq so that is will effectively overwrite SE q with CSLAM q
          !
!          elem(ie)%derived%fq(:,:,:,m_cnst,1) = fld_gll(:,:,:,m_cnst+3,ie)*elem(ie)%state%dp3d(:,:,:,tl_f)-&
!               elem(ie)%state%Qdp(:,:,:,m_cnst,tl_qdp)
          elem(ie)%derived%fq(:,:,:,m_cnst,1) = fld_gll(:,:,:,m_cnst+3,ie)-&
               elem(ie)%state%Q(:,:,:,m_cnst)


        end do
        fvm(ie)%fc(1:nc,1:nc,:,1:ntrac) = fvm(ie)%fc_phys(1:nc,1:nc,:,1:ntrac)
     end do
  end if


  deallocate(fld_phys,llimiter,fld_gll)

#ifdef debug_coupling
    !
    ! map fc tendencies to gll grid
    !
    allocate(fld_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,ntrac,nets:nete))
    allocate(fld_gll(np,np,nlev,ntrac,nets:nete))
    allocate(llimiter(ntrac))
    do ie=nets,nete
      !
      ! compute cslam updated Q value
      !
      do m_cnst=1,ntrac
        fld_fvm(1:nc,1:nc,:,m_cnst,ie) = fvm(ie)%fc(1:nc,1:nc,:,m_cnst)!/fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm)
      end do
   end do
   do ie=nets,nete
      do m_cnst=2,ntrac
         do k=1,nlev
            do j=1,nc
               do i=1,nc
                  fld_fvm(i,j,k,m_cnst,ie)=fld_fvm(i,j,k,m_cnst,ie)-&
                       test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon,k,m_cnst)

!                 write(*,*) "xxx ",i,j,k,ie,m_cnst,test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon,k,m_cnst)-&
!                      fvm(ie)%fc_phys(i,j,k,m_cnst)/fvm(ie)%dp_fvm(i,j,k,n0_fvm)
!                 write(*,*) "zzz",i,j,k,ie,m_cnst,test_func(fvm(ie)%center_cart(i,j)%lat,fvm(ie)%center_cart(i,j)%lon,k,m_cnst),&
!                      fvm(ie)%fc_phys(i,j,k,m_cnst)/fvm(ie)%dp_fvm(i,j,k,n0_fvm),fvm(ie)%dp_fvm(i,j,k,n0_fvm)
!                 write(*,*) "yyy",i,j,k,ie,m_cnst,fvm(ie)%fc_phys(i,j,k,m_cnst)-&
!                      fvm(ie)%fc(i,j,k,m_cnst),fvm(ie)%fc_phys(i,j,k,m_cnst),&
!                      fvm(ie)%fc(i,j,k,m_cnst)
               end do
            end do
         end do
      end do
   end do
   llimiter = .false.

!    call phys2dyn(elem,fld_fvm,fld_gll(:,:,:,1:ntrac,:),hybrid,nets,nete,&
    call fvm2dyn(elem,fld_fvm,fld_gll(:,:,:,1:ntrac,:),hybrid,nets,nete,&
         ntrac,fvm,llimiter(1:ntrac))
    !
    ! fld_gll now holds q cslam value on gll grid
    !
    ! convert fld_gll to increment (q_new-q_old)*dp
    !
    do ie=nets,nete
       do m_cnst=1,ntrac
          elem(ie)%derived%fq(:,:,:,m_cnst,1)   = fld_gll(:,:,:,m_cnst,ie)!*&
          !             elem(ie)%state%dp3d(:,:,:,tl_f)
       end do
    end do
    deallocate(fld_fvm,llimiter,fld_gll)
#endif
  end subroutine phys2dyn_forcings_fvm

  subroutine fvm2dyn(elem,fld_fvm,fld_gll,hybrid,nets,nete,num_flds,fvm,llimiter)
    use dimensions_mod, only: np, nlev, nhc, nc
    use hybrid_mod    , only: hybrid_t
    use bndry_mod     , only: ghost_exchangeV
    use edge_mod      , only: initghostbufferTR, freeghostbuffertr, ghostVpack, ghostVunpack
    use edgetype_mod  , only: ghostbuffertr_t
    use cam_abortutils, only: endrun

    integer              , intent(in)    :: nets,nete,num_flds
    real (kind=real_kind), intent(inout) :: fld_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,num_flds,nets:nete)
    real (kind=real_kind), intent(out)   :: fld_gll(np,np,nlev,num_flds,nets:nete)
    type (hybrid_t)      , intent(in)    :: hybrid  ! distributed parallel structure (shared)
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(in)    :: fvm(:)
    logical              , intent(in)    :: llimiter(num_flds)

    integer                 :: i, j, ie, k, iwidth
    type (ghostBuffertr_t)  :: cellghostbuf
    !
    !*********************************************
    !
    ! halo exchange
    !
    !*********************************************
    !
    call initghostbufferTR(cellghostbuf,nlev,num_flds,nhc,nc)
    do ie=nets,nete
       call ghostVpack(cellghostbuf, fld_fvm(:,:,:,:,ie),nhc,nc,nlev,num_flds,0,elem(ie)%desc)
    end do
    call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc,num_flds)
    do ie=nets,nete
       call ghostVunpack(cellghostbuf, fld_fvm(:,:,:,:,ie),nhc,nc,nlev,num_flds,0,elem(ie)%desc)
    end do
    call freeghostbuffertr(cellghostbuf)
    !
    ! mapping
    !
    iwidth=2
    do ie=nets,nete
      call tensor_lagrange_interp(fvm(ie)%cubeboundary,np,nc,nhc,nlev,num_flds,fld_fvm(:,:,:,:,ie),&
           fld_gll(:,:,:,:,ie),llimiter,iwidth,fvm(ie)%norm_elem_coord)
    end do
  end subroutine fvm2dyn


  subroutine fill_halo_phys(elem,fld_phys,fld_gll,hybrid,nets,nete,num_flds,fvm)
    use dimensions_mod, only: np, nlev, nhc_phys, fv_nphys
    use hybrid_mod    , only: hybrid_t
    use bndry_mod     , only: ghost_exchangeV
    use edge_mod      , only: initghostbufferTR, freeghostbuffertr, ghostVpack, ghostVunpack
    use edgetype_mod  , only: ghostbuffertr_t

    integer              , intent(in)    :: nets,nete,num_flds
    real (kind=real_kind), intent(inout) :: fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,num_flds, &
         nets:nete)
    real (kind=real_kind), intent(out)   :: fld_gll(np,np,nlev,num_flds,nets:nete)
    type (hybrid_t)      , intent(in)    :: hybrid  ! distributed parallel structure (shared)
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(in)    :: fvm(:)

    integer                 :: ie
    type (ghostBuffertr_t)  :: cellghostbuf
    !
    !*********************************************
    !
    ! halo exchange
    !
    !*********************************************
    !
    call initghostbufferTR(cellghostbuf,nlev,num_flds,nhc_phys,fv_nphys)
    do ie=nets,nete
       call ghostVpack(cellghostbuf, fld_phys(:,:,:,:,ie),nhc_phys,fv_nphys,nlev,num_flds,0,elem(ie)%desc)
    end do
    call ghost_exchangeV(hybrid,cellghostbuf,nhc_phys,fv_nphys,num_flds)
    do ie=nets,nete
       call ghostVunpack(cellghostbuf, fld_phys(:,:,:,:,ie),nhc_phys,fv_nphys,nlev,num_flds,0,elem(ie)%desc)
    end do
    call freeghostbuffertr(cellghostbuf)
  end subroutine fill_halo_phys

  !
  ! must call fill_halo_phys before calling this subroutine
  !
  subroutine phys2dyn(elem,fld_phys,fld_gll,nets,nete,num_flds,fvm,llimiter,istart_vector)
    use dimensions_mod, only: np, nlev, nhc_phys, fv_nphys

    integer              , intent(in)    :: nets,nete,num_flds
    real (kind=real_kind), intent(inout) :: fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,num_flds, &
         nets:nete)
    real (kind=real_kind), intent(out)   :: fld_gll(np,np,nlev,num_flds,nets:nete)
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(in)    :: fvm(:)
    integer, optional    , intent(in)    :: istart_vector
    logical              , intent(in)    :: llimiter(num_flds)

    integer                 :: i, j, ie, k, iwidth
    real (kind=real_kind)   :: v(2)

    if (present(istart_vector)) then
      do ie=nets,nete
        do k=1,nlev
          do j=1-nhc_phys,fv_nphys+nhc_phys
            do i=1-nhc_phys,fv_nphys+nhc_phys
              !
              ! convert lat-lon vectors to contra-variant gnomonic
              !
              v(1) = fld_phys(i,j,k,istart_vector  ,ie)
              v(2) = fld_phys(i,j,k,istart_vector+1,ie)
              fld_phys(i,j,k,istart_vector  ,ie)=fvm(ie)%Dinv_physgrid(i,j,1,1)*v(1) + fvm(ie)%Dinv_physgrid(i,j,1,2)*v(2)
              fld_phys(i,j,k,istart_vector+1,ie)=fvm(ie)%Dinv_physgrid(i,j,2,1)*v(1) + fvm(ie)%Dinv_physgrid(i,j,2,2)*v(2)
            end do
          end do
        end do
      end do
    end if

    !
    ! mapping
    !
    iwidth=2
    if (fv_nphys==1) iwidth=1
    do ie=nets,nete
      call tensor_lagrange_interp(fvm(ie)%cubeboundary,np,fv_nphys,nhc_phys,nlev,num_flds,fld_phys(:,:,:,:,ie),&
           fld_gll(:,:,:,:,ie),llimiter,iwidth,fvm(ie)%norm_elem_coord_physgrid)
    end do

    if (present(istart_vector)) then
      !
      ! convert contra-variant to lat-lon
      !
      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              v(1) = fld_gll(i,j,k,istart_vector  ,ie)
              v(2) = fld_gll(i,j,k,istart_vector+1,ie)
              fld_gll(i,j,k,istart_vector  ,ie) = elem(ie)%D(i,j,1,1)*v(1) + elem(ie)%D(i,j,1,2)*v(2)
              fld_gll(i,j,k,istart_vector+1,ie) = elem(ie)%D(i,j,2,1)*v(1) + elem(ie)%D(i,j,2,2)*v(2)
            end do
          end do
        end do
      end do
    end if

  end subroutine phys2dyn


  !
  ! map all mass variables from gll to fvm
  !
  subroutine dyn2fvm_mass_vars(dp_gll,ps_gll,q_gll,&
       dp_fvm,ps_fvm,q_fvm,num_trac,metdet,inv_area,lcslam)
    use dimensions_mod, only: np, nc,nlev
    integer, intent(in) :: num_trac
    real (kind=real_kind), dimension(np,np,nlev)         , intent(in) :: dp_gll
    real (kind=real_kind), dimension(np,np,nlev,num_trac), intent(in) :: q_gll
    real (kind=real_kind), dimension(np,np)              , intent(in) :: ps_gll


    real (kind=real_kind), dimension(nc,nc,nlev)         , intent(inout) :: dp_fvm
    real (kind=real_kind), dimension(nc,nc,nlev,num_trac), intent(inout) :: q_fvm
    real (kind=real_kind), dimension(nc,nc)     , intent(inout)        :: ps_fvm
    real (kind=real_kind), dimension(nc,nc)     , intent(out)          :: inv_area

    real (kind=real_kind), intent(in)           :: metdet(np,np)
    logical              , intent(in)           :: lcslam

    real (kind=real_kind) :: se_area_sphere(nc,nc), tmp(np,np)
    real (kind=real_kind) :: inv_darea_dp_fvm(nc,nc)
    integer :: i,j,k,m_cnst

    tmp = 1.0D0
    se_area_sphere = dyn2fvm(tmp,metdet)
    inv_area = 1.0D0/se_area_sphere
    if (.not.lcslam) then
      ps_fvm(:,:) = dyn2fvm(ps_gll,metdet,inv_area)
      do k=1,nlev
        dp_fvm(:,:,k) = dyn2fvm(dp_gll(:,:,k),metdet,inv_area)
        inv_darea_dp_fvm = inv_area/dp_fvm(:,:,k)
        do m_cnst=1,num_trac
          q_fvm(:,:,k,m_cnst) = &
               dyn2fvm(q_gll(:,:,k,m_cnst)*dp_gll(:,:,k),metdet,&
               inv_darea_dp_fvm,q_gll(:,:,k,m_cnst))
        end do
!        dp_fvm(:,:,k) = dp_fvm(:,:,k)*inv_area(:,:)
      end do
    end if
  end subroutine dyn2fvm_mass_vars
  !
  ! this subroutine assumes that the fvm halo has already been filled
  ! (if nc/=fv_nphys)
  !
  subroutine dyn2phys_all_vars(ie,dp_gll,ps_gll,q_gll,T_gll,omega_gll,phis_gll,&
       dp_fvm,ps_fvm,q_fvm,inv_area_fvm,num_trac,metdet,lcslam,fvm,ptop,&
       dp3d_phys,ps_phys,q_phys,T_phys,omega_phys,phis_phys)
    use dimensions_mod, only: np, nc,nlev,fv_nphys,nhc
    use dp_grids      , only: nphys_pts
    integer, intent(in) :: ie,num_trac
    real (kind=real_kind), dimension(np,np,nlev)         , intent(in) :: dp_gll,T_gll,omega_gll
    real (kind=real_kind), dimension(np,np,nlev,num_trac), intent(in) :: q_gll
    real (kind=real_kind), dimension(np,np)              , intent(in) :: ps_gll,phis_gll


    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev)         , intent(inout) :: dp_fvm
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,num_trac), intent(inout) :: q_fvm
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)              , intent(inout) :: ps_fvm
    real (kind=real_kind), dimension(nc,nc)              , intent(in)    :: inv_area_fvm
    type(fvm_struct)                                                         , intent(in)    :: fvm

    real (kind=real_kind), intent(in)           :: metdet(np,np)
    logical              , intent(in)           :: lcslam
    real (kind=real_kind), intent(in)           :: ptop

    real (kind=real_kind), dimension(nphys_pts)               , intent(out) :: ps_phys,phis_phys
    real (kind=real_kind), dimension(nphys_pts,nlev)          , intent(out) :: dp3d_phys,T_phys,omega_phys
    real (kind=real_kind), dimension(nphys_pts,nlev,num_trac) , intent(out) :: q_phys


    real (kind=real_kind) :: tmp(np,np)
    real (kind=real_kind), dimension(fv_nphys,fv_nphys)          :: inv_area,inv_darea_dp_phys,se_area_sphere,dp3d_tmp
    real (kind=real_kind), dimension(fv_nphys,fv_nphys)          :: dp_phys_tmp
    real (kind=real_kind), dimension(fv_nphys,fv_nphys,num_trac) :: q_phys_tmp

    integer :: i,j,k,m_cnst

    tmp = 1.0D0
    se_area_sphere = dyn2phys(tmp,metdet)
    inv_area = 1.0D0/se_area_sphere
    phis_phys(:) = RESHAPE(dyn2phys(phis_gll,metdet,inv_area),SHAPE(phis_phys(:)))

    ps_phys = ptop
    do k=1,nlev
      dp3d_tmp       = dyn2phys(dp_gll(:,:,k),metdet,inv_area)
      inv_darea_dp_phys = inv_area/dp3d_tmp
      T_phys(:,k) = RESHAPE(dyn2phys(T_gll(:,:,k)*dp_gll(:,:,k),metdet,&
           inv_darea_dp_phys),SHAPE(T_phys(:,k)))
      Omega_phys(:,k) = RESHAPE(dyn2phys(Omega_gll(:,:,k),metdet,inv_area),SHAPE(Omega_phys(:,k)))

      if (lcslam) then
        call fvm2phys(ie,fvm,dp_fvm(:,:,k),dp_phys_tmp,q_fvm(:,:,k,:),q_phys_tmp,num_trac)
        dp3d_phys(:,k) = RESHAPE(dp_phys_tmp,SHAPE(dp3d_phys(:,k)))
        ps_phys(:) = ps_phys(:)+RESHAPE(dp_phys_tmp,SHAPE(ps_phys(:)))
        do m_cnst=1,num_trac
          q_phys(:,k,m_cnst) = RESHAPE(q_phys_tmp(:,:,m_cnst),SHAPE(q_phys(:,k,m_cnst)))
        end do
      else
        dp3d_phys(:,k) = RESHAPE(dp3d_tmp,SHAPE(dp3d_phys(:,k)))
        ps_phys(:) = ps_phys(:)+RESHAPE(dp_phys_tmp,SHAPE(ps_phys(:)))
        do m_cnst=1,num_trac
          q_phys(:,k,m_cnst) = RESHAPE( &
               dyn2phys(q_gll(:,:,k,m_cnst)*dp_gll(:,:,k),metdet,&
               inv_darea_dp_phys,q_gll(:,:,k,m_cnst)),&
               SHAPE(q_phys(:,k,m_cnst)))
        end do
      end if
    end do
  end subroutine dyn2phys_all_vars

  function dyn2phys(qdp_gll,metdet,inv_dp_darea_phys,q_gll) result(qdp_phys)
    use dimensions_mod, only: np, nc, fv_nphys
    use derivative_mod, only: subcell_integration
    real (kind=real_kind), intent(in)           :: qdp_gll(np,np)
    real (kind=real_kind)                       :: qdp_phys(fv_nphys,fv_nphys)
    real (kind=real_kind), intent(in)           :: metdet(np,np)
    real (kind=real_kind), intent(in), optional :: inv_dp_darea_phys(fv_nphys,fv_nphys)
    real (kind=real_kind), intent(in), optional :: q_gll(np,np)
    integer :: i,j
    real (kind=real_kind) :: min_val, max_val

    qdp_phys = subcell_integration(qdp_gll(:,:), np, fv_nphys, metdet,nc.ne.fv_nphys)
    if (present(inv_dp_darea_phys)) then
      !
      ! convert qdp to q
      !
      qdp_phys = qdp_phys*inv_dp_darea_phys
      !
      ! simple limiter
      !
      if (present(q_gll)) then
        min_val = minval(q_gll)
        max_val = maxval(q_gll)
        do j = 1, fv_nphys
          do i = 1, fv_nphys
            !
            ! simple limiter: only coded for nc=3 and np4
            !
!            min_val = minval(q_gll(i:i+1,j:j+1))
!            max_val = maxval(q_gll(i:i+1,j:j+1))
            qdp_phys(i,j) = max(min_val,min(max_val,qdp_phys(i,j)))
          end do
        end do
      end if
    end if
  end function dyn2phys


  function dyn2fvm(qdp_gll,metdet,inv_dp_darea_phys,q_gll) result(qdp_phys)
    use dimensions_mod, only: np, nc
    use derivative_mod, only: subcell_integration
    real (kind=real_kind), intent(in)           :: qdp_gll(np,np)
    real (kind=real_kind)                       :: qdp_phys(nc,nc)
    real (kind=real_kind), intent(in)           :: metdet(np,np)
    real (kind=real_kind), intent(in), optional :: inv_dp_darea_phys(nc,nc)
    real (kind=real_kind), intent(in), optional :: q_gll(np,np)
    integer :: i,j
    real (kind=real_kind) :: min_val, max_val

    qdp_phys = subcell_integration(qdp_gll(:,:), np, nc, metdet)
    if (present(inv_dp_darea_phys)) then
      !
      ! convert qdp to q
      !
      qdp_phys = qdp_phys*inv_dp_darea_phys
      !
      ! simple limiter
      !
      if (present(q_gll)) then
        do j = 1, nc
          do i = 1, nc
            !
            ! simple limiter: only coded for nc=3 and np4
            !
            min_val = minval(q_gll(i:i+1,j:j+1))
            max_val = maxval(q_gll(i:i+1,j:j+1))
            qdp_phys(i,j) = max(min_val,min(max_val,qdp_phys(i,j)))
          end do
        end do
      end if
    end if
  end function dyn2fvm

  function dyn2phys_vector(v_gll,elem) result(v_phys)
    use dimensions_mod, only: np, nc, nlev,  fv_nphys
    use interpolate_mod,only: interpolate_vector, interpdata_t,interpolate_2d,interpolate_t
    use cube_mod       ,only: dmap
    use control_mod    ,only: cubed_sphere_map

    type (interpdata_t):: interpdata
    type (element_t), intent(in)   :: elem
    type (interpolate_t) , target :: interp_p
    real (kind=real_kind), intent(in)           :: v_gll(np,np,2,nlev)
    real (kind=real_kind)                       :: v_phys(fv_nphys*fv_nphys,2,nlev)

    integer :: i,j,k

    ! Local variables
    real (kind=real_kind)    ::  fld_contra(np,np,2,nlev) ! vector field

    real (kind=real_kind)    ::  v1,v2
    real (kind=real_kind)    ::  D(2,2,fv_nphys*fv_nphys)   ! derivative of gnomonic mapping
!    real (kind=real_kind)    ::  JJ(2,2), tmpD(2,2)   ! derivative of gnomonic mapping
    !
    ! this could be done at initialization and does not need to be repeated
    !
    call setup_interpdata_for_gll_to_phys_vec_mapping(elem,interpdata, interp_p)
    ! convert to contra
    do k=1,nlev
      do j=1,np
        do i=1,np
          ! latlon->contra
          fld_contra(i,j,1,k) = elem%Dinv(i,j,1,1)*v_gll(i,j,1,k) + elem%Dinv(i,j,1,2)*v_gll(i,j,2,k)
          fld_contra(i,j,2,k) = elem%Dinv(i,j,2,1)*v_gll(i,j,1,k) + elem%Dinv(i,j,2,2)*v_gll(i,j,2,k)
        enddo
      enddo
    end do

    do k=1,nlev
      do i=1,interpdata%n_interp
        v_phys(i,1,k)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp_p,np)
        v_phys(i,2,k)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp_p,np)
!      v_phys(i,1)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,np)
        !      v_phys(i,2)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,np)
      end do
    end do
    do i=1,interpdata%n_interp
      ! convert fld from contra->latlon
      call dmap(D(:,:,i),interpdata%interp_xy(i)%x,interpdata%interp_xy(i)%y,&
           elem%corners3D,cubed_sphere_map,elem%corners,elem%u2qmap,elem%facenum)
    end do
    do k=1,nlev
      do i=1,interpdata%n_interp
        ! convert fld from contra->latlon
        v1 = v_phys(i,1,k)
        v2 = v_phys(i,2,k)

        v_phys(i,1,k)=D(1,1,i)*v1 + D(1,2,i)*v2
        v_phys(i,2,k)=D(2,1,i)*v1 + D(2,2,i)*v2
      end do
    end do
  end function dyn2phys_vector

  subroutine setup_interpdata_for_gll_to_phys_vec_mapping(elem,interpdata,interp_p)
    !
    ! initialize interpolation data structures to interpolate to phys grid
    ! using interpolate_mod subroutines
    !
    use interpolate_mod, only: interpolate_t, interpdata_t, interpolate_create
    use dimensions_mod, only : np
    use quadrature_mod, only : quadrature_t, gausslobatto
    use dimensions_mod, only : fv_nphys
    type (element_t)     , intent(in) , target :: elem
    type (interpdata_t)  , intent(out)         :: interpdata
    type (interpolate_t) , intent(out), target :: interp_p

    ! local
    type (quadrature_t)   :: gp_quadrature
    integer i,j,ioff,ngrid
    real (kind=real_kind) ::  dx

    ngrid = fv_nphys*fv_nphys
    interpdata%n_interp=ngrid
    !
    ! initialize interpolation stuff related to basis functions
    !
    gp_quadrature = gausslobatto(np)
    call interpolate_create(gp_quadrature,interp_p)
    allocate(interpdata%interp_xy(ngrid))
    allocate(interpdata%ilat(ngrid) )
    allocate(interpdata%ilon(ngrid) )
    !
    !WARNING: THIS CODE INTERFERES WITH LAT-LON OUTPUT
    !         OF REGULAR SE IF nc>0
    !
    ioff=1
    dx = 2.0D0/dble(fv_nphys)
    do j=1,fv_nphys
      do i=1,fv_nphys
        interpdata%interp_xy(ioff)%x = -1D0+(i-0.5D0)*dx
        interpdata%interp_xy(ioff)%y = -1D0+(j-0.5D0)*dx
        interpdata%ilon(ioff) = i
        interpdata%ilat(ioff) = j
        ioff=ioff+1
      enddo
    enddo
  end subroutine setup_interpdata_for_gll_to_phys_vec_mapping


  function lagrange_1d(src_grid,src_val,ngrid,dst_point,iwidth) result(val)
    integer              , intent(in)  :: ngrid,iwidth
    real (kind=real_kind), intent(in)  :: src_grid(ngrid), src_val(ngrid)
    real (kind=real_kind)              :: val

    real (kind=real_kind), intent(in)  :: dst_point

    integer :: iref, j,k
    real (kind=real_kind)              :: w(ngrid)

    if (dst_point.LE.src_grid(1)) then
      iref=1
    else
      iref=1
      do while (dst_point>src_grid(iref))
        iref = iref + 1
        if (iref>ngrid) then
!          write(*,*) "extrapolating",dst_point,src_grid,iref
          !        call endrun("search out of bounds")
          exit
        end if
      end do
      iref=iref-1
    end if

!    if (iref<2.or.iref>ngrid-2) then
!      call endrun("src_grid not wide enought for cubic interpolation")
!      write(*,*) "xxx min/max involved",iref,dst_point,src_grid
    iref=MIN(MAX(iref,iwidth),ngrid-iwidth)
!    endif

    w = 1.0D0
    do j=iref-(iwidth-1),iref+iwidth
      do k=iref-(iwidth-1),iref+iwidth
        if (k.ne.j) then
          w(j)=w(j)*(dst_point-src_grid(k))/(src_grid(j)-src_grid(k))
        end if
      end do
    end do

    val=0.0D0
    do j=iref-(iwidth-1),iref+iwidth
      val=val+w(j)*src_val(j)
    end do
  end function lagrange_1d

  subroutine tensor_lagrange_interp(cubeboundary,np,nc,nhc,nlev,nflds,psi,interp_value,llimiter,iwidth,norm_elem_coord)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use cam_abortutils,         only: endrun
    implicit none

    integer              , intent(in)    :: cubeboundary,nc, np, iwidth,nhc,nlev,nflds
    logical              , intent(in)    :: llimiter(nflds)                   !apply limiter
    real (kind=real_kind), intent(inout) :: psi(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,nflds) !fvm grid values with filled halo
    real (kind=real_kind), intent(out)   :: interp_value(np,np,nlev,nflds)            !interpolated field
    real (kind=real_kind), intent(in)    :: norm_elem_coord(2,1-nhc:nc+nhc,1-nhc:nc+nhc)
    integer :: which_nc_cell(np)

    real (kind=real_kind):: dx,gll_points(np)
    real (kind=real_kind):: nc_points(1-nc:nc+nc)

    real (kind=real_kind):: value(1-iwidth:nc+iwidth)
    real (kind=real_kind):: val_tmp(1-nhc:nc+nhc,1-nhc:nc+nhc)

    real (kind=real_kind):: min_value(np,np,nlev,nflds), max_value(np,np,nlev,nflds)

    integer :: imin(1-nhc:nc+nhc), imax(1-nhc:nc+nhc)
    integer :: k,i,j,isearch,igll,jgll,jrow,ilev,h,irow,itr

    gll_points(1) = -1.0D0
    gll_points(2) = -sqrt(1D0/5D0)
    gll_points(3) =  sqrt(1D0/5D0)
    gll_points(4) =  1.0D0

    dx = 2D0/dble(nc)
    do k=1-nc,2*nc
      nc_points(k) = -1D0+dx*0.5D0+dble(k-1)*dx
    end do
    !
    ! find fvm point surrounding gll points for simple limiter
    !
    do k=1,np
      do isearch=0,nc+1
        if (nc_points(isearch)<gll_points(k).and.nc_points(isearch+1).ge.gll_points(k)) exit
      end do
      which_nc_cell(k)=isearch
    end do
    do itr=1,nflds
      if (llimiter(itr)) then
        !
        ! fill non-existent halo cells for limiter
        !
        if (cubeboundary>4) then
          h=1
          select case(cubeboundary)
          case (nwest)
            psi(0,nc+h  ,:,itr) = psi(1-h,nc  ,:,itr)
            psi(1-h,nc+1,:,itr) = psi(1  ,nc+h,:,itr)
          case (swest)
            psi(1-h,0,:,itr) = psi(1,1-h,:,itr)
            psi(0,1-h,:,itr) = psi(1-h,1,:,itr)
          case (seast)
            psi(nc+h,0,:,itr) = psi(nc,1-h,:,itr)
            psi(nc+1,1-h,:,itr) = psi(nc+h,1,:,itr)
          case (neast)
            psi(nc+h,nc+1,:,itr) = psi(nc,nc+h,:,itr)
            psi(nc+1,nc+h,:,itr) = psi(nc+h,nc,:,itr)
          end select
        end if
        do k=1,nlev
          do j=1,np
            do i=1,np
              max_value(i,j,k,itr) = max(&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)+1,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)+1,k,itr) &
                   )
              min_value(i,j,k,itr) = min(&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)  ,k,itr),&
                   psi(which_nc_cell(i)  ,which_nc_cell(j)+1,k,itr),&
                   psi(which_nc_cell(i)+1,which_nc_cell(j)+1,k,itr) &
                   )
            end do
          end do
        end do
      end if
    end do

    imin=1-nhc
    imax=nc+nhc
    !
    ! special corner treatment
    !
    if (cubeboundary==swest) then
      do itr=1,nflds
        do k=1,nlev
          do jrow=1,nc+iwidth
            !
            ! cubic along constant x (i=irow) in west halo to fvm points in halo
            !
            do irow=1-iwidth,0
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1:nc+nhc),psi(irow,1:nc+nhc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(1-iwidth:0,1:nc+iwidth,k,itr) = val_tmp(1-iwidth:0,1:nc+iwidth)
        enddo
      end do
      imin(1-nhc:0) = 1
    end if
    if (cubeboundary==nwest) then
      do itr=1,nflds
        do k=1,nlev
          do jrow=1-iwidth,nc
            !
            ! cubic along constant x (i=irow) in west halo to fvm points in halo
            !
            do irow=1-iwidth,0
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1-nhc:nc),psi(irow,1-nhc:nc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(1-iwidth:0,1-iwidth:nc,k,itr) = val_tmp(1-iwidth:0,1-iwidth:nc)
        end do
      end do
      imin(nc+1:nc+nhc) = 1
    end if

    if (cubeboundary==seast) then
      do itr=1,nflds
        do k=1,nlev
          do jrow=1,nc+iwidth
            value=0.0D0
            !
            ! cubic along constant y in ease halo to fvm points in halo
            !
            do irow=nc+1,nc+iwidth
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1:nc+nhc),psi(irow,1:nc+nhc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(nc+1:nc+iwidth,1:nc+iwidth,k,itr) = val_tmp(nc+1:nc+iwidth,1:nc+iwidth)
        end do
      end do
      imax(1-nhc:0) = nc
    end if

    if (cubeboundary==neast) then
      do itr=1,nflds
        do k=1,nlev
          do jrow=1-iwidth,nc
            !
            ! cubic along constant y in ease halo to fvm points in halo
            !
            do irow=nc+1,nc+iwidth
              val_tmp(irow,jrow) = lagrange_1d(norm_elem_coord(2,irow,1-nhc:nc),psi(irow,1-nhc:nc,k,itr),nc+nhc,&
                   norm_elem_coord(2,1,jrow),iwidth)
            end do
          end do
          psi(nc+1:nc+iwidth,1-iwidth:nc,k,itr) = val_tmp(nc+1:nc+iwidth,1-iwidth:nc)
        end do
      end do
      imax(nc+1:nc+nhc) = nc
    end if
    !
    ! mapping
    !
    !
    if (cubeboundary==0.or.cubeboundary==north.or.cubeboundary==south.or.&
         cubeboundary==swest.or.cubeboundary==nwest.or.&
         cubeboundary==seast.or.cubeboundary==neast) then
      do itr=1,nflds
        do k=1,nlev
          do igll=1,np
            !
            ! cubic along constant y (j=jrow)
            !
            do jrow=1-iwidth,nc+iwidth
              value(jrow) = lagrange_1d(norm_elem_coord(1,imin(jrow):imax(jrow),jrow),psi(imin(jrow):imax(jrow),jrow,k,itr),&
                   imax(jrow)-imin(jrow)+1,gll_points(igll),iwidth)
            end do
            do jgll=1,np
              interp_value(igll,jgll,k,itr) = lagrange_1d(norm_elem_coord(2,1,1-iwidth:nc+iwidth),value,nc+2*iwidth,&
                   gll_points(jgll),iwidth)
            end do
          end do
        end do
      end do
    else if (cubeboundary==east.or.cubeboundary==west) then
      do itr=1,nflds
        do k=1,nlev
          do jgll=1,np
            !
            ! cubic along constant x (i=irow)
            !
            do irow=1-iwidth,nc+iwidth
              value(irow) = lagrange_1d(norm_elem_coord(2,irow,1-nhc:nc+nhc),psi(irow,1-nhc:nc+nhc,k,itr),nc+2*nhc,&
                   gll_points(jgll),iwidth)
            end do
            do igll=1,np
              interp_value(igll,jgll,k,itr) = lagrange_1d(norm_elem_coord(1,1-iwidth:nc+iwidth,1),value,nc+2*iwidth,&
                   gll_points(igll),iwidth)
            end do
          end do
        end do
      end do
    end if
    do itr=1,nflds
      if (llimiter(itr)) then
        do k=1,nlev
          do j=1,np
            do i=1,np
              interp_value(i,j,k,itr)=max(min_value(i,j,k,itr),min(max_value(i,j,k,itr),interp_value(i,j,k,itr)))
            end do
          enddo
        end do
      end if
    end do
  end subroutine tensor_lagrange_interp


  subroutine fvm2phys(ie,fvm,dp_fvm,dp_phys,q_fvm,q_phys,num_trac)
    use dimensions_mod, only: nc,nhr,nhc,ns,irecons_tracer,fv_nphys
    use fvm_reconstruction_mod, only: reconstruction
    !
    ! weights must be initialized in fvm2phys_init before using this function
    !
    use dp_grids      , only: weights_all_fvm2phys,weights_eul_index_all_fvm2phys,weights_lgr_index_all_fvm2phys
    use dp_grids      , only: jall_fvm2phys,num_weights_fvm2phys,num_weights_fvm2phys

    !
    ! setting nhe=0 because we do not need reconstruction outside of element
    !
    integer, parameter :: nhe_local=0
    integer, parameter :: nh = nhr!+(nhe-1) ! = 2 (nhr=2; nhe_local=1),! = 3 (nhr=2; nhe_local=2)

    !
    ! xxx change nhe=0
    !
    type(fvm_struct)     , intent(in)           :: fvm
    integer              , intent(in)           :: ie
    integer              , intent(in)           :: num_trac
    real (kind=real_kind), intent(inout)        :: dp_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,1)
    real (kind=real_kind), intent(out)          :: dp_phys(fv_nphys,fv_nphys)

    real (kind=real_kind), intent(inout)        :: q_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,num_trac)
    real (kind=real_kind), intent(out)          :: q_phys(fv_nphys,fv_nphys,num_trac)

    real (kind=real_kind)                       :: recons    (1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,irecons_tracer,1)
    real (kind=real_kind)                       :: recons_q  (1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,irecons_tracer, &
         num_trac)
    real (kind=real_kind)                       :: recons_tmp(irecons_tracer)

    logical                                     :: llimiter(1),llimiter_q(num_trac)
    integer :: h,jx,jy,jdx,jdy,jall,m_cnst
    integer :: jx_min_local(3), jx_max_local(3), jy_min_local(3), jy_max_local(3)
    real (kind=real_kind)                       :: dp_phys_inv(fv_nphys,fv_nphys)

    real (kind=real_kind)                       :: dp_tmp
    integer :: nht_local

    llimiter=.false.
    nht_local=nhe_local+nhr     !total halo width where reconstruction is needed (nht<=nc) - phl
    !
    ! to accomodate nhe=0 make sure nothing is done for neighboring panel recontructions
    !
    jx_min_local(1) = 1            ; jx_max_local(1) = nc+1
    jy_min_local(1) = 1            ; jy_max_local(1) = nc+1
    jx_min_local(2) = 0            ; jx_max_local(2) = -1
    jy_min_local(2) = 0            ; jy_max_local(2) = -1
    jx_min_local(3) = 0            ; jx_max_local(3) = -1
    jy_min_local(3) = 0            ; jy_max_local(3) = -1

    call reconstruction(dp_fvm,recons,irecons_tracer,llimiter,1,&
       nc,nhe_local,nhr,nhc,nht_local,ns,nh,&
       jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
       fvm%cubeboundary,fvm%halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,:),fvm%ibase(1-nh:nc+nh,1:nhr,:),&
       fvm%spherecentroid(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,:),&
       fvm%recons_metrics(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,:),&
       fvm%recons_metrics_integral(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,:)    ,&
       fvm%rot_matrix,fvm%centroid_stretch(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,1:7),&
       fvm%vertex_recons_weights(:,1:irecons_tracer-1,1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local),&
       fvm%vtx_cart(1-nhc:nc+nhc,1-nhc:nc+nhc,:,:))


    dp_phys = 0.0D0
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie)
       jy  = weights_lgr_index_all_fvm2phys(h,2,ie)
       jdx = weights_eul_index_all_fvm2phys(h,1,ie)
       jdy = weights_eul_index_all_fvm2phys(h,2,ie)
       dp_phys(jx,jy) = dp_phys(jx,jy) + SUM(weights_all_fvm2phys(h,:,ie)*recons(jdx,jdy,:,1))
    end do

    llimiter_q=.true.
    call reconstruction(q_fvm,recons_q,irecons_tracer,llimiter_q,num_trac,&
       nc,nhe_local,nhr,nhc,nht_local,ns,nh,&
       jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
       fvm%cubeboundary,fvm%halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,:),fvm%ibase(1-nh:nc+nh,1:nhr,:),&
       fvm%spherecentroid(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,:),&
       fvm%recons_metrics(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,:),&
       fvm%recons_metrics_integral(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,:)    ,&
       fvm%rot_matrix,fvm%centroid_stretch(1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local,1:7),&
       fvm%vertex_recons_weights(:,1:irecons_tracer-1,1-nhe_local:nc+nhe_local,1-nhe_local:nc+nhe_local),&
       fvm%vtx_cart(1-nhc:nc+nhc,1-nhc:nc+nhc,:,:))
    !
    ! q-dp coupling as described in equation (55) in Appendinx B of
    ! Nair and Lauritzen, 2010: A Class of Deformational Flow Test Cases for Linear Transport Problems on the Sphere.
    ! J. Comput. Phys.: Vol. 229, Issue 23, pp. 8868–8887, DOI:10.1016/j.jcp.2010.08.014.
    !
    q_phys = 0.0D0
    do h=1,jall_fvm2phys(ie)
       jx  = weights_lgr_index_all_fvm2phys(h,1,ie)
       jy  = weights_lgr_index_all_fvm2phys(h,2,ie)
       jdx = weights_eul_index_all_fvm2phys(h,1,ie)
       jdy = weights_eul_index_all_fvm2phys(h,2,ie)
       recons_tmp    = recons(jdx,jdy,:,1)
       recons_tmp(1) = recons(jdx,jdy,1,1)-dp_fvm(jdx,jdy,1)
       dp_tmp = SUM(weights_all_fvm2phys(h,:,ie)*recons_tmp(:))
       do m_cnst=1,num_trac
         q_phys(jx,jy,m_cnst) = q_phys(jx,jy,m_cnst) + &
              dp_fvm(jdx,jdy,1)*SUM(weights_all_fvm2phys(h,:,ie)*recons_q(jdx,jdy,:,m_cnst))+&
              q_fvm(jdx,jdy,m_cnst) *dp_tmp
       end do
    end do
    !
    ! convert to mixing ratio
    !
    dp_phys_inv=1.0D0/dp_phys
    do m_cnst=1,num_trac
      q_phys(:,:,m_cnst) = q_phys(:,:,m_cnst)*dp_phys_inv(:,:)
    end do
  end subroutine fvm2phys

#ifdef xxx
  subroutine phys2fvm(ie,fvm,q_fvm,q_phys,num_trac)
    use dimensions_mod, only: nc,nhr_phys,nhc_phys,ns_phys,irecons_tracer,fv_nphys,nhe_phys,nlev
    use fvm_reconstruction_mod, only: reconstruction
    !
    ! weights must be initialized in phys2fvm_init before using this function
    !
    use dp_grids      , only: weights_all_phys2fvm,weights_eul_index_all_phys2fvm,weights_lgr_index_all_phys2fvm
    use dp_grids      , only: jall_phys2fvm,num_weights_phys2fvm,num_weights_phys2fvm

    !
    ! setting nhe=0 because we do not need reconstruction outside of element
    !
    integer, parameter :: nh_phys = nhr_phys!+(nhe-1) ! = 2 (nhr=2; nhe_phys=1),! = 3 (nhr=2; nhe_phys=2)
    integer, parameter :: nht_local=nhe_phys+nhr_phys     !total halo width where reconstruction is needed (nht<=nc) - phl
    !
    ! xxx change nhe=0
    !
    type(fvm_struct)     , intent(in)           :: fvm
    integer              , intent(in)           :: ie

    real (kind=real_kind), intent(inout)        :: q_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,num_trac)
    real (kind=real_kind), intent(out)          :: q_fvm(nc,nc,nlev,num_trac)
    integer              , intent(in)           :: num_trac

    real (kind=real_kind)                       :: recons_q  (1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys, &
         irecons_tracer,num_trac)

    logical                                     :: llimiter_q(num_trac)
    integer :: h,jx,jy,jdx,jdy,jall,m_cnst,k
    integer :: jx_min_local(3), jx_max_local(3), jy_min_local(3), jy_max_local(3)

    llimiter_q=.true.
    !
    ! to accomodate nhe=0 make sure nothing is done for neighboring panel recontructions
    !
    jx_min_local(1) = 1            ; jx_max_local(1) = fv_nphys+1
    jy_min_local(1) = 1            ; jy_max_local(1) = fv_nphys+1
    jx_min_local(2) = 0            ; jx_max_local(2) = -1
    jy_min_local(2) = 0            ; jy_max_local(2) = -1
    jx_min_local(3) = 0            ; jx_max_local(3) = -1
    jy_min_local(3) = 0            ; jy_max_local(3) = -1

    q_fvm = 0.0D0
    do k=1,nlev
      call reconstruction(q_phys(:,:,k,:),recons_q,irecons_tracer,llimiter_q,num_trac,&
           fv_nphys,nhe_phys,nhr_phys,nhc_phys,nht_local,ns_phys,nh_phys,&
           jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
           fvm%cubeboundary,fvm%halo_interp_weight_physgrid(1:ns_phys,1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
           fvm%ibase_physgrid(1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
           fvm%spherecentroid_physgrid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:),&
           fvm%recons_metrics_physgrid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:),&
           fvm%recons_metrics_integral_physgrid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:)    ,&
           fvm%rot_matrix,fvm%centroid_stretch_physgrid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,1:7),&
           fvm%vertex_recons_weights_physgrid(:,1:irecons_tracer-1,1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys),&
           fvm%vtx_cart_physgrid(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,:,:))

      do m_cnst=1,num_trac
        do h=1,jall_phys2fvm(ie)
          jx  = weights_lgr_index_all_phys2fvm(h,1,ie)
          jy  = weights_lgr_index_all_phys2fvm(h,2,ie)
          jdx = weights_eul_index_all_phys2fvm(h,1,ie)
          jdy = weights_eul_index_all_phys2fvm(h,2,ie)
          q_fvm(jx,jy,k,m_cnst) = q_fvm(jx,jy,k,m_cnst) + SUM(weights_all_phys2fvm(h,:,ie)*recons_q(jdx,jdy,:,m_cnst))
        end do
      end do
    end do
  end subroutine phys2fvm
#endif

  subroutine phys2fvm(ie,fvm,dp_phys,dp_fvm,q_phys,q_fvm,num_trac)
    use dimensions_mod, only: fv_nphys,nhr_phys,nhc_phys,ns_phys,irecons_tracer,fv_nphys,nhe_phys,nc
    use fvm_reconstruction_mod, only: reconstruction
    !
    ! weights must be initialized in phys2fvm_init before using this function
    !
    use dp_grids      , only: weights_all_phys2fvm,weights_eul_index_all_phys2fvm,weights_lgr_index_all_phys2fvm
    use dp_grids      , only: jall_phys2fvm,num_weights_phys2fvm,num_weights_phys2fvm

    !
    ! xxx change nhe=0
    !
    type(fvm_struct)     , intent(in)           :: fvm
    integer              , intent(in)           :: ie
    integer              , intent(in)           :: num_trac
    real (kind=real_kind), intent(inout)        :: dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,1)
    real (kind=real_kind), intent(out)          :: dp_fvm(nc,nc)

    real (kind=real_kind), intent(inout)        :: q_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,num_trac)
    real (kind=real_kind), intent(out)          :: q_fvm (nc,nc,num_trac)

    real (kind=real_kind)                       :: recons    (1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys, &
         irecons_tracer,1)
    real (kind=real_kind)                       :: recons_q  (1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys, &
         irecons_tracer,num_trac)
    real (kind=real_kind)                       :: recons_tmp(irecons_tracer)

    logical                                     :: llimiter(1),llimiter_q(num_trac)
    integer :: h,jx,jy,jdx,jdy,jall,m_cnst
    integer :: jx_min_local(3), jx_max_local(3), jy_min_local(3), jy_max_local(3)
    real (kind=real_kind)                       :: dp_fvm_inv(nc,nc)

    real (kind=real_kind)                       :: dp_tmp

    !
    ! setting nhe=0 because we do not need reconstruction outside of element
    !
    integer:: nh_phys
    integer:: nht_local

    nh_phys = nhr_phys!-1!+(nhe-1) ! = 2 (nhr=2; nhe_phys=1),! = 3 (nhr=2; nhe_phys=2)
    nht_local=nhe_phys+nhr_phys     !total halo width where reconstruction is needed (nht<=nc) - phl


    llimiter=.false.
    !
    ! to accomodate nhe=0 make sure nothing is done for neighboring panel recontructions
    !
    jx_min_local(1) = 1            ; jx_max_local(1) = fv_nphys+1
    jy_min_local(1) = 1            ; jy_max_local(1) = fv_nphys+1
    jx_min_local(2) = 0            ; jx_max_local(2) = -1
    jy_min_local(2) = 0            ; jy_max_local(2) = -1
    jx_min_local(3) = 0            ; jx_max_local(3) = -1
    jy_min_local(3) = 0            ; jy_max_local(3) = -1

    call reconstruction(dp_phys,recons,irecons_tracer,llimiter,1,&
       fv_nphys,nhe_phys,nhr_phys,nhc_phys,nht_local,ns_phys,nh_phys,&
       jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
       fvm%cubeboundary,fvm%halo_interp_weight(1:ns_phys,1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
       fvm%ibase(1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
       fvm%spherecentroid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:),&
       fvm%recons_metrics(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:),&
       fvm%recons_metrics_integral(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:)    ,&
       fvm%rot_matrix,fvm%centroid_stretch(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,1:7),&
       fvm%vertex_recons_weights(:,1:irecons_tracer-1,1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys),&
       fvm%vtx_cart(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,:,:))


    dp_fvm = 0.0D0
    do h=1,jall_phys2fvm(ie)
       jx  = weights_lgr_index_all_phys2fvm(h,1,ie)
       jy  = weights_lgr_index_all_phys2fvm(h,2,ie)
       jdx = weights_eul_index_all_phys2fvm(h,1,ie)
       jdy = weights_eul_index_all_phys2fvm(h,2,ie)
       dp_fvm(jx,jy) = dp_fvm(jx,jy) + SUM(weights_all_phys2fvm(h,:,ie)*recons(jdx,jdy,:,1))
    end do

    llimiter_q=.true.
    call reconstruction(q_phys,recons_q,irecons_tracer,llimiter_q,num_trac,&
       fv_nphys,nhe_phys,nhr_phys,nhc_phys,nht_local,ns_phys,nh_phys,&
       jx_min_local,jx_max_local,jy_min_local,jy_max_local,&
       fvm%cubeboundary,fvm%halo_interp_weight(1:ns_phys,1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
       fvm%ibase(1-nh_phys:fv_nphys+nh_phys,1:nhr_phys,:),&
       fvm%spherecentroid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:),&
       fvm%recons_metrics(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:),&
       fvm%recons_metrics_integral(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,:)    ,&
       fvm%rot_matrix,fvm%centroid_stretch(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,1:7),&
       fvm%vertex_recons_weights(:,1:irecons_tracer-1,1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys),&
       fvm%vtx_cart(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,:,:))
    !
    ! q-dp coupling as described in equation (55) in Appendinx B of
    ! Nair and Lauritzen, 2010: A Class of Deformational Flow Test Cases for Linear Transport Problems on the Sphere.
    ! J. Comput. Phys.: Vol. 229, Issue 23, pp. 8868–8887, DOI:10.1016/j.jcp.2010.08.014.
    !
    q_fvm = 0.0D0
    do h=1,jall_phys2fvm(ie)
       jx  = weights_lgr_index_all_phys2fvm(h,1,ie)
       jy  = weights_lgr_index_all_phys2fvm(h,2,ie)
       jdx = weights_eul_index_all_phys2fvm(h,1,ie)
       jdy = weights_eul_index_all_phys2fvm(h,2,ie)
       recons_tmp    = recons(jdx,jdy,:,1)
       recons_tmp(1) = recons(jdx,jdy,1,1)-dp_fvm(jdx,jdy)
       dp_tmp = SUM(weights_all_phys2fvm(h,:,ie)*recons_tmp(:))
       do m_cnst=1,num_trac
         q_fvm(jx,jy,m_cnst) = q_fvm(jx,jy,m_cnst) + &
              dp_fvm(jdx,jdy)*SUM(weights_all_phys2fvm(h,:,ie)*recons_q(jdx,jdy,:,m_cnst))+&
              q_fvm(jdx,jdy,m_cnst) *dp_tmp
       end do
    end do
    !
    ! convert to mixing ratio
    !
    dp_fvm_inv=1.0D0/dp_fvm
    do m_cnst=1,num_trac
      q_fvm(:,:,m_cnst) = q_fvm(:,:,m_cnst)*dp_fvm_inv(:,:)
    end do
  end subroutine phys2fvm


end module fvm_mapping
