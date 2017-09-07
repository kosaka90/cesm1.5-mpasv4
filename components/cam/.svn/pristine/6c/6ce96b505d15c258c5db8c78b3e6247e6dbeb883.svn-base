!-----------------------------------------------------------------------------------!
!MODULE FVM_MOD-----------------------------------------------------------CE-for FVM!
! FVM_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
! 7.Februar 2012: cslam_run and cslam_runair                                        !
!-----------------------------------------------------------------------------------!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_mod    
  use kinds, only : real_kind
  use edge_mod, only : initghostbufferTR, freeghostbuffertr, &
       ghostVpack, ghostVunpack!,  initEdgebuffer
  use edgetype_mod, only : ghostbuffertr_t, edgebuffer_t
  use bndry_mod, only: ghost_exchangeV                     
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod, only : hybrid_t
  
  implicit none
  private
  save
  
  type (EdgeBuffer_t)                         :: edgeveloc

  public :: edgeveloc, fvm_init1,fvm_init2, fill_halo_fvm, physgrid_init2,fvm_init3
contains

  subroutine fill_halo_fvm(elem,fvm,hybrid,nets,nete,tnp0,ndepth,kmin,kmax)
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    use dimensions_mod, only: irecons_tracer,nc, nlev, ntrac
    implicit none
    type (element_t),intent(inout)            :: elem(:)
    type (fvm_struct),intent(inout)           :: fvm(:)
    type (hybrid_t),intent(in)                :: hybrid

    type (ghostBuffertr_t)                    :: cellghostbuf


    integer,intent(in)                        :: nets,nete,tnp0,ndepth,kmin,kmax
    integer                                   :: ie,i1,i2,num_levels
    !
    ! note "call initghostbufferTR(cellghostbuf,nlev,ntrac+1,nhc,nc)" in fvm_init1.
    ! should initghostbuffer be called here?
    !

    i1=1-ndepth
    i2=nc+ndepth
    num_levels = kmax-kmin+1
    call initghostbufferTR(cellghostbuf,num_levels,ntrac+1,ndepth,nc)
    call t_startf('FVM pack')
    do ie=nets,nete
       call ghostVpack(cellghostbuf, fvm(ie)%dp_fvm(i1:i2,i1:i2,kmin:kmax,tnp0),ndepth,nc,num_levels,1,    0,elem(ie)%desc)
       call ghostVpack(cellghostbuf, fvm(ie)%c(i1:i2,i1:i2,kmin:kmax,:,tnp0)   ,ndepth,nc,num_levels,ntrac,1,elem(ie)%desc)
    end do
    call t_stopf('FVM pack')
    call t_startf('FVM Communication')
    call ghost_exchangeV(hybrid,cellghostbuf,ndepth,nc,ntrac+1)
    call t_stopf('FVM Communication')
    !-----------------------------------------------------------------------------------!                        
    call t_startf('FVM Unpack')
    do ie=nets,nete
       call ghostVunpack(cellghostbuf, fvm(ie)%dp_fvm(i1:i2,i1:i2,kmin:kmax,tnp0), ndepth, nc,num_levels,1    ,0,elem(ie)%desc)
       call ghostVunpack(cellghostbuf, fvm(ie)%c(i1:i2,i1:i2,kmin:kmax,:,tnp0),    ndepth, nc,num_levels,ntrac,1,elem(ie)%desc)
    enddo
    call t_stopf('FVM Unpack')
    call freeghostbuffertr(cellghostbuf)
  end subroutine fill_halo_fvm

  
  ! initialize global buffers shared by all threads
  subroutine fvm_init1(par,elem)
    use parallel_mod, only : parallel_t, haltmp
    use control_mod, only : tracer_transport_type, tracer_grid_type, rsplit
    use control_mod, only : TRACERTRANSPORT_CONSISTENT_SE_FVM
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm, fvm_supercycling
    use dimensions_mod, only: qsize, qsize_d, irecons_tracer
    use dimensions_mod, only: nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc

    type (parallel_t) :: par
    type (element_t)  :: elem(:)
    
    !
    ! initialize fvm time-levels
    !
    n0_fvm  = 1
    np1_fvm = 2
    !
    if (ntrac>0) then
      if (par%masterproc) then 
        print *, "                                           "
        print *, "|-----------------------------------------|"
        print *, "| FVM tracer transport scheme information |"
        print *, "|-----------------------------------------|"
        print *, "                                           "
      end if
      if (tracer_transport_type == TRACERTRANSPORT_CONSISTENT_SE_FVM) then
        if (par%masterproc) then 
          print *, "Running consistent SE-CSLAM, Lauritzen et al. (2015, in prep)."
          print *, "Air flux prescribed by SE (Taylor, Overfelt, Ullrich)"
          print *, "SE = Spectral Element"
          print *, "CSLAM = Conservative Semi-LAgrangian Multi-tracer scheme"
          print *, "Lauritzen et al., (2010), J. Comput. Phys."
          print *, "  "
        end if
      else
!        call haltmp("going into fvm_init1 with inconsistent tracer_transport_type")
      end if

      if (ntrac>ntrac_d) &
           call haltmp("PARAMETER ERROR for fvm: ntrac > ntrac_d")

            if (qsize>0.and.mod(rsplit,fvm_supercycling).ne.0) then
        if (par%masterproc) then
          print *,'cannot supercycle fvm tracers with respect to se tracers'
          print *,'with this choice of rsplit =',rsplit
          print *,'rsplit must be a multiple of fvm_supercycling=',fvm_supercycling
          call haltmp("PARAMETER ERROR for fvm: mod(rsplit,fvm_supercycling<>0")
        end if
      endif
      
      
      if (par%masterproc) then 
        print *, "                                            "
        print *, "Done Tracer transport scheme information    "
        print *, "                                            "
      end if
    end if

      
    if (par%masterproc) print *, "fvm resolution is nc*nc in each element: nc = ",nc
    if (par%masterproc) print *,'ntrac,ntrac_d=',ntrac,ntrac_d      
    if (par%masterproc) print *,'ntrac,ntrac_d=',qsize,qsize_d

    
    if (nc<3) then
      if (par%masterproc) then 
        print *, "NUMBER OF CELLS ERROR for fvm: Number of cells parameter"
        print *, "parameter nc at least 3 (nc>=3), nc*nc cells per element. This is"
        print *, "needed for the cubic reconstruction, which is only implemented yet! STOP"
      endif
      call haltmp("stopping")
    end if
    
    if (par%masterproc) then
      print *, "  "
      if (ns==1) then
        print *, "ns==1: using no interpolation for mapping cell averages values across edges"
        print *, "Note: this is not a recommended setting - large errors at panel edges!"
      else if (ns==2) then
        print *, "ns==2: using linear interpolation for mapping cell averages values across edges"
        print *, "Note that ns=4 is default CSLAM setting used in Lauritzen et al. (2010)"
        print *, "so this option is slightly less accurate (but the stencil is smaller near panel edges!)"
        
      else if (ns==3) then
        print *, "ns==3: using quadratic interpolation for mapping cell averages values across edges"
        print *, "Note that ns=4 is default CSLAM setting used in Lauritzen et al. (2010)"
        print *, "so this option is slightly less accurate (but the stencil is smaller near panel edges!)"
      else if (ns==4) then
        print *, "ns==4: using cubic interpolation for mapping cell averages values across edges"
        print *, "This is default CSLAM setting used in Lauritzen et al. (2010)"
      else 
        print *, "Not a tested value for ns but it should work! You choose ns = ",ns
      end if
      
      !       if (ns.NE.3) then
      !         write(*,*) "In fvm_reconstruction_mod function matmul_w has been hard-coded for ns=3 for performance"
      !         write(*,*) "Revert to general code - outcommented above"
      !         call haltmp("stopping")
      !       end if
    end if
    
    if (MOD(ns,2)==0.and.nhr+(nhe-1)+ns/2>nc+nc) then
      print *, "to run this combination of ns and nhr you need to increase nc to ",nhr+ns/2+nhe-1
      print *, "You choose (ns,nhr,nc,nhe)=",ns,nhr,nc,nhe
      call haltmp("stopping")
    end if
    if (MOD(ns,2)==1.and.nhr+(ns-1)/2+(nhe-1)>nc+nc) then
      print *, "to run this combination of ns and nhr you need to increase nc to ",nhr+(ns-1)/2+nhe-1
      print *, "You choose (ns,nhr,nc,nhe)=",ns,nhr,nc,nhe
      call haltmp("stopping")
    end if
    
    if (nc==3.and.ns.ne.3) then
      if (par%masterproc) then
        print *, "Recommended setting for nc=3 is ns=3 (linear interpolation in halo)"
        print *, "You choose ns=",ns
        print *, "Goto dimensions_mod to change value of ns"
        print *, "or outcomment call haltmop below (i.e. you know what you are doing!)"
      endif
      call haltmp("stopping")
    end if
    
    if (nc==4.and.ns.ne.4) then
      if (par%masterproc) then
        print *, "Recommended setting for nc=4 is ns=4 (cubic interpolation in halo)"
        print *, "You choose ns=",ns
        print *, "Goto dimensions_mod to change value of ns"
        print *, "or outcomment call haltmop below (i.e. you know what you are doing!)"
      endif
      call haltmp("stopping")
    end if
    
    if (nhe .ne. 1) then
      if (par%masterproc) then
        print *, "PARAMETER ERROR for fvm: Number of halo zone for the extended"
        print *,"element nhe has to be 1, only this is available now! STOP!"
      endif
      call haltmp("stopping")
    end if
    
    
!    call initghostbufferTR(cellghostbuf,nlev,ntrac+1,nhc,nc)
!    call initEdgebuffer(par,edgeveloc,elem,2*nlev)
  end subroutine fvm_init1
  
  
  
  
  
  ! initialization that can be done in threaded regions
  subroutine fvm_init2(elem,fvm,hybrid,nets,nete)
    use fvm_control_volume_mod, only: fvm_mesh,fvm_set_cubeboundary
    use bndry_mod, only: compute_ghost_corner_orientation
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm
    use dimensions_mod, only: irecons_tracer
    use dimensions_mod, only: nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc

    !  use edge_mod, only : ghostvpackfull, ghostvunpackfull
    
    type (fvm_struct) :: fvm(:)
    type (element_t) :: elem(:)
    type (hybrid_t)                             :: hybrid
    integer :: ie,nets,nete

    n0_fvm  = 1 !in case no cslam but physgrid
    np1_fvm = 2
    
    call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
    ! run some tests:
    !    call test_ghost(hybrid,elem,nets,nete)
    
    do ie=nets,nete
      call fvm_set_cubeboundary(elem(ie),fvm(ie))
      call fvm_mesh(elem(ie),fvm(ie))
      fvm(ie)%inv_area_sphere    = 1.0D0/fvm(ie)%area_sphere
      fvm(ie)%inv_se_area_sphere = fvm(ie)%inv_area_sphere

      fvm(ie)%fc(:,:,:,:) = 0.0D0
      fvm(ie)%fm(:,:,:,:) = 0.0D0
      fvm(ie)%ft(:,:,:  ) = 0.0D0
    enddo
  end subroutine fvm_init2

  
  subroutine fvm_init3(elem,fvm,hybrid,nets,nete,irecons)
    use control_mod     , only: north, south, east, west, neast, nwest, seast, swest
    use fvm_analytic_mod, only: compute_reconstruct_matrix
    use dimensions_mod  , only: fv_nphys
    use dimensions_mod, only: nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc
    use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
    use fvm_control_volume_mod, only: n0_fvm
    use coordinate_systems_mod, only : cubedsphere2cart, cart2cubedsphere
    implicit none
    type (element_t) ,intent(inout)  :: elem(:)
    type (fvm_struct),intent(inout)  :: fvm(:) 
    type (hybrid_t)  ,intent(in)     :: hybrid                      
    integer          ,intent(in)     :: nets,nete,irecons
    !
    type (ghostBuffertr_t)  :: cellghostbuf
    integer                 :: ie, ixy, ivertex, i, j,istart,itot,ishft,imin,imax
    integer, dimension(2,4) :: unit_vec
    integer                 :: rot90_matrix(2,2), iside

    type (cartesian2D_t)                :: tmpgnom
    type (cartesian2D_t)                :: gnom
    type(cartesian3D_t)                 :: tmpcart3d


    if (ntrac>0.and.nc.ne.fv_nphys) then
      !
      ! fill the fvm halo for mapping in d_p_coupling if
      ! physics grid resolution is different than fvm resolution
      !
      call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm,nhc,1,nlev)
    end if


    imin=1-nhc
    imax=nc+nhc
    !
    ! fill halo start
    !
    itot=9+irecons-1+2
    call initghostbufferTR(cellghostbuf,1,itot,nhc,nc)
    do ie=nets,nete
      istart = 0
      call ghostVpack(cellghostbuf, fvm(ie)%norm_elem_coord(1,:,:),nhc,nc,1,1,istart,elem(ie)%desc)
      istart = istart+1
      call ghostVpack(cellghostbuf, fvm(ie)%norm_elem_coord(2,:,:),nhc,nc,1,1,istart,elem(ie)%desc)
      istart = istart+1        
      do ixy=1,2
        do ivertex=1,4
          call ghostVpack(cellghostbuf, fvm(ie)%vtx_cart(:,:,ixy,ivertex) ,nhc,nc,1,1,istart,elem(ie)%desc)
          istart = istart+1
        end do
      end do
      call ghostVpack(cellghostbuf, fvm(ie)%flux_orient(1,:,:) ,nhc,nc,1,1,istart,elem(ie)%desc)
      do ixy=1,irecons-1
        istart=istart+1
        call ghostVpack(cellghostbuf, fvm(ie)%spherecentroid(:,:,ixy) ,nhc,nc,1,1,istart,elem(ie)%desc)
      end do
    end do
    call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc,itot)
    do ie=nets,nete
      istart = 0
      call ghostVunpack(cellghostbuf, fvm(ie)%norm_elem_coord(1,:,:),nhc,nc,1,1,istart,elem(ie)%desc)
      istart = istart+1
      call ghostVunpack(cellghostbuf, fvm(ie)%norm_elem_coord(2,:,:),nhc,nc,1,1,istart,elem(ie)%desc)
      istart = istart+1        
      do ixy=1,2
        do ivertex=1,4
          call ghostVunpack(cellghostbuf, fvm(ie)%vtx_cart(:,:,ixy,ivertex) ,nhc,nc,1,1,istart,elem(ie)%desc)
          istart = istart+1
        end do
      end do
      call ghostVunpack(cellghostbuf, fvm(ie)%flux_orient(1,:,:) ,nhc,nc,1,1,istart,elem(ie)%desc)
      do ixy=1,irecons-1
        istart=istart+1
        call ghostVunpack(cellghostbuf, fvm(ie)%spherecentroid(:,:,ixy) ,nhc,nc,1,1,istart,elem(ie)%desc)
      end do
    enddo
    call freeghostbuffertr(cellghostbuf)    
    !
    ! indicator for non-existing cells 
    ! set vtx_cart to corner value in non-existent cells
    !
    do ie=nets,nete
       if (fvm(ie)%cubeboundary==nwest) then
         fvm(ie)%flux_orient     (:  ,1-nhc      :0     ,nc      +1 :nc      +nhc      ) = -1
         fvm(ie)%spherecentroid  (    1-nhc      :0     ,nc      +1 :nc      +nhc    ,:) = -1D5
         fvm(ie)%vtx_cart(1-nhc:0     ,nc+1 :nc+nhc,1,:) = fvm(ie)%vtx_cart(1,nc,1,4)
         fvm(ie)%vtx_cart(1-nhc:0     ,nc+1 :nc+nhc,2,:) = fvm(ie)%vtx_cart(1,nc,2,4)
       else if (fvm(ie)%cubeboundary==swest) then
         fvm(ie)%flux_orient     (:,1-nhc      :0     ,1-nhc      :0     ) = -1
         fvm(ie)%spherecentroid  (  1-nhc      :0     ,1-nhc      :0   ,:) = -1D5
         fvm(ie)%vtx_cart(1-nhc:0     ,1-nhc:0     ,1,:) = fvm(ie)%vtx_cart(1,1,1,1)
         fvm(ie)%vtx_cart(1-nhc:0     ,1-nhc:0     ,2,:) = fvm(ie)%vtx_cart(1,1,2,1)
       else if (fvm(ie)%cubeboundary==neast) then
         fvm(ie)%flux_orient     (:,nc      +1 :nc      +nhc      ,nc      +1 :nc      +nhc      ) = -1
         fvm(ie)%spherecentroid  (  nc      +1 :nc      +nhc      ,nc      +1 :nc      +nhc    ,:) = -1D5
         fvm(ie)%vtx_cart(nc+1 :nc+nhc,nc+1 :nc+nhc,1,:) = fvm(ie)%vtx_cart(nc,nc,1,3)
         fvm(ie)%vtx_cart(nc+1 :nc+nhc,nc+1 :nc+nhc,2,:) = fvm(ie)%vtx_cart(nc,nc,2,3)
       else if (fvm(ie)%cubeboundary==seast) then
         fvm(ie)%flux_orient     (:,nc      +1 :nc      +nhc      ,1-nhc      :0     ) = -1
         fvm(ie)%spherecentroid  (  nc      +1 :nc      +nhc      ,1-nhc      :0   ,:) = -1D5
         fvm(ie)%vtx_cart(nc+1 :nc+nhc,1-nhc:0     ,1,:) = fvm(ie)%vtx_cart(nc,1,1,2)
         fvm(ie)%vtx_cart(nc+1 :nc+nhc,1-nhc:0     ,2,:) = fvm(ie)%vtx_cart(nc,1,2,2)
       end if
     end do
     
     !
     ! set vectors for perpendicular flux vector
     !
     rot90_matrix(1,1) = 0; rot90_matrix(2,1) =  1 !counter-clockwise rotation matrix
     rot90_matrix(1,2) =-1; rot90_matrix(2,2) =  0 !counter-clockwise rotation matrix 
     
     iside = 1
     unit_vec(1,iside) = 0 !x-component of displacement vector for side 1
     unit_vec(2,iside) = 1 !y-component of displacement vector for side 1
     
     do iside=2,4
       unit_vec(:,iside) = MATMUL(rot90_matrix(:,:),unit_vec(:,iside-1))
     end do
     
     !
     ! fill halo done
     !
     !-------------------------------
     
     do ie=nets,nete
       fvm(ie)%displ_max = 0.0D0
       do j=imin,imax
         do i=imin,imax
           !
           ! rotate gnomonic coordinate vector
           !
           !           fvm(ie)%norm_elem_coord(:,i,j) = MATMUL(fvm(ie)%rot_matrix(:,:,i,j),fvm(ie)%norm_elem_coord(:,i,j))
           !    
           ishft = NINT(fvm(ie)%flux_orient(2,i,j))
           do ixy=1,2
             !
             ! rotate coordinates if needed through permutation
             !
             fvm(ie)%vtx_cart(i,j,ixy,1:4) = cshift(fvm(ie)%vtx_cart(i,j,ixy,1:4),shift=ishft)
             fvm(ie)%flux_vec        (ixy,i,j,1:4) = cshift(unit_vec                (ixy,1:4    ),shift=ishft)
             !
             ! set flux vector to zero in non-existent cells (corner halo) 
             !
             fvm(ie)%flux_vec        (ixy,i,j,1:4) = fvm(ie)%ifct(i,j)*fvm(ie)%flux_vec(ixy,i,j,1:4)
             
             iside=1
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(i,j,ixy,4)-fvm(ie)%vtx_cart(i,j,ixy,1))
             iside=2
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(i,j,ixy,1)-fvm(ie)%vtx_cart(i,j,ixy,2))
             iside=3
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(i,j,ixy,2)-fvm(ie)%vtx_cart(i,j,ixy,3))
             iside=4
             fvm(ie)%displ_max(i,j,iside) = fvm(ie)%displ_max(i,j,iside)+&
                  ABS(fvm(ie)%vtx_cart(i,j,ixy,2)-fvm(ie)%vtx_cart(i,j,ixy,1))
           end do
         end do
       end do
     end do
     !
     ! pre-compute derived metric terms used for integration, polynomial
     ! evaluation at fvm cell vertices, etc.
     !    
     do ie=nets,nete
       call compute_reconstruct_matrix(nc,nhe,nhc,irecons,fvm(ie)%dalpha,fvm(ie)%dbeta,&
           fvm(ie)%spherecentroid,fvm(ie)%vtx_cart,fvm(ie)%centroid_stretch,&
           fvm(ie)%vertex_recons_weights,fvm(ie)%recons_metrics,fvm(ie)%recons_metrics_integral)
     end do
      !
      ! create a normalized element coordinate system with a halo
      !    
     do ie=nets,nete
       do j=1-nhc,nc+nhc
         do i=1-nhc,nc+nhc
           !
           ! only compute for physically existent cells
           !
           if (fvm(ie)%ifct(i,j)>0) then
             gnom%x = fvm(ie)%norm_elem_coord(1,i,j)
             gnom%y = fvm(ie)%norm_elem_coord(2,i,j)
             !
             ! coordinate transform only necessary for points on another panel
             !
             if (NINT(fvm(ie)%flux_orient(1,1,1)).NE.NINT(fvm(ie)%flux_orient(1,i,j))) then
               tmpcart3d=cubedsphere2cart(gnom,NINT(fvm(ie)%flux_orient(1,i,j)))
               tmpgnom=cart2cubedsphere(tmpcart3d,NINT(fvm(ie)%flux_orient(1,1,1)))
             else
               tmpgnom%x = fvm(ie)%norm_elem_coord(1,i,j)
               tmpgnom%y = fvm(ie)%norm_elem_coord(2,i,j)
             end if
             !
             ! convert to element normalized coordinates
             !
             fvm(ie)%norm_elem_coord(1,i,j) =(tmpgnom%x-elem(ie)%corners(1)%x)/&
                  (0.5D0*dble(nc)*fvm(ie)%dalpha)-1.0D0
             fvm(ie)%norm_elem_coord(2,i,j) =(tmpgnom%y-elem(ie)%corners(1)%y)/&
                  (0.5D0*dble(nc)*fvm(ie)%dalpha)-1.0D0
           else
             fvm(ie)%norm_elem_coord(1,i,j) = 1D9
             fvm(ie)%norm_elem_coord(2,i,j) = 1D9
           end if
         end do
       end do
     end do

   end subroutine fvm_init3
  

  subroutine physgrid_init2(elem, fvm, hybrid, nets, nete,irecons)
    use edge_mod, only : initghostbufferTR, freeghostbuffertr, ghostVpack, ghostVunpack
    use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use coordinate_systems_mod, only : cubedsphere2cart, cart2cubedsphere
    use dimensions_mod, only: fv_nphys, nhe_phys,ns_phys,nhr_phys,nhc_phys
    use dimensions_mod, only: nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc
    use cube_mod       ,only: dmap
    use control_mod    ,only: cubed_sphere_map
    use fvm_analytic_mod, only: compute_reconstruct_matrix

    type (element_t) , intent(in)    :: elem(:)
    type (fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t)  , intent(in)    :: hybrid

    type (cartesian2D_t)                :: gnom
    type(cartesian3D_t)                 :: tmpcart3d
    type (cartesian2D_t)                :: tmpgnom


    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete,irecons  ! ending thread element number   (private)

    ! ==================================
    ! Local variables
    ! ==================================

    integer                 :: ie, ixy, ivertex, i, j,istart,itot,ishft,imin,imax
    integer, dimension(2,4) :: unit_vec
    integer                 :: rot90_matrix(2,2), iside
    real (kind=real_kind)   :: cartx,carty    

    type (ghostBuffertr_t)                    :: cellghostbuf

    ! D is derivative of gnomonic mapping
    real (kind=real_kind)              :: D(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,2)
    real (kind=real_kind)              :: detD,dx,x1,x2

    if (fv_nphys>0) then
      !
      ! do the same as fvm_init3 for the metric terms of physgrid
      !
      imin=1-nhc_phys
      imax=fv_nphys+nhc_phys
      !
      ! fill halo start
      !
      itot=9+irecons-1+2
      call initghostbufferTR(cellghostbuf,1,itot,nhc_phys,fv_nphys)
      do ie=nets,nete
        istart = 0
        call ghostVpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(1,:,:),nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        istart = istart+1
        call ghostVpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(2,:,:),nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        istart = istart+1
        do ixy=1,2
          do ivertex=1,4
            call ghostVpack(cellghostbuf, fvm(ie)%vtx_cart_physgrid(:,:,ixy,ivertex) ,nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
            istart = istart+1
          end do
        end do
        call ghostVpack(cellghostbuf, fvm(ie)%flux_orient_physgrid(1,:,:) ,nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        do ixy=1,irecons-1
          istart=istart+1
          call ghostVpack(cellghostbuf, fvm(ie)%spherecentroid_physgrid(:,:,ixy) ,nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        end do
      end do
      call ghost_exchangeV(hybrid,cellghostbuf,nhc_phys,fv_nphys,itot)
      do ie=nets,nete
        istart = 0
        call ghostVunpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(1,:,:),nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        istart = istart+1
        call ghostVunpack(cellghostbuf, fvm(ie)%norm_elem_coord_physgrid(2,:,:),nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        istart = istart+1
        do ixy=1,2
          do ivertex=1,4
            call ghostVunpack(cellghostbuf, fvm(ie)%vtx_cart_physgrid(:,:,ixy,ivertex) ,nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
            istart = istart+1
          end do
        end do
        call ghostVunpack(cellghostbuf, fvm(ie)%flux_orient_physgrid(1,:,:) ,nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        do ixy=1,irecons-1
          istart=istart+1
          call ghostVunpack(cellghostbuf, fvm(ie)%spherecentroid_physgrid(:,:,ixy) ,nhc_phys,fv_nphys,1,1,istart,elem(ie)%desc)
        end do
      enddo
      call freeghostbuffertr(cellghostbuf)    
      !
      ! indicator for non-existing cells 
      ! set vtx_cart to corner value in non-existent cells
      !
      do ie=nets,nete
        if (fvm(ie)%cubeboundary==nwest) then
          fvm(ie)%flux_orient_physgrid   (:  ,1-nhc_phys      :0     ,fv_nphys      +1 :fv_nphys      +nhc_phys      ) = -1
          fvm(ie)%spherecentroid_physgrid(    1-nhc_phys      :0     ,fv_nphys      +1 :fv_nphys      +nhc_phys    ,:) = -1D5
          fvm(ie)%vtx_cart_physgrid(1-nhc_phys:0     ,fv_nphys+1 :fv_nphys+nhc_phys,1,:) = &
               fvm(ie)%vtx_cart_physgrid(1,fv_nphys,1,4)
          fvm(ie)%vtx_cart_physgrid(1-nhc_phys:0     ,fv_nphys+1 :fv_nphys+nhc_phys,2,:) = &
               fvm(ie)%vtx_cart_physgrid(1,fv_nphys,2,4)
        else if (fvm(ie)%cubeboundary==swest) then
          fvm(ie)%flux_orient_physgrid   (:,1-nhc_phys      :0     ,1-nhc_phys      :0     ) = -1
          fvm(ie)%spherecentroid_physgrid(  1-nhc_phys      :0     ,1-nhc_phys      :0   ,:) = -1D5
          fvm(ie)%vtx_cart_physgrid(1-nhc_phys:0     ,1-nhc_phys:0     ,1,:) = fvm(ie)%vtx_cart_physgrid(1,1,1,1)
          fvm(ie)%vtx_cart_physgrid(1-nhc_phys:0     ,1-nhc_phys:0     ,2,:) = fvm(ie)%vtx_cart_physgrid(1,1,2,1)
        else if (fvm(ie)%cubeboundary==neast) then
          fvm(ie)%flux_orient_physgrid   (:,fv_nphys      +1 :fv_nphys      +nhc_phys      , &
               fv_nphys      +1 :fv_nphys      +nhc_phys      ) = -1
          fvm(ie)%spherecentroid_physgrid(  fv_nphys      +1 :fv_nphys      +nhc_phys      , &
               fv_nphys      +1 :fv_nphys      +nhc_phys    ,:) = -1D5
          fvm(ie)%vtx_cart_physgrid(fv_nphys+1 :fv_nphys+nhc_phys,fv_nphys+1 :fv_nphys+nhc_phys,1,:) = &
               fvm(ie)%vtx_cart_physgrid(fv_nphys,fv_nphys,1,3)
          fvm(ie)%vtx_cart_physgrid(fv_nphys+1 :fv_nphys+nhc_phys,fv_nphys+1 :fv_nphys+nhc_phys,2,:) = &
               fvm(ie)%vtx_cart_physgrid(fv_nphys,fv_nphys,2,3)
        else if (fvm(ie)%cubeboundary==seast) then
          fvm(ie)%flux_orient_physgrid   (:,fv_nphys      +1 :fv_nphys      +nhc_phys      ,1-nhc_phys      :0     ) = -1
          fvm(ie)%spherecentroid_physgrid(  fv_nphys      +1 :fv_nphys      +nhc_phys      ,1-nhc_phys      :0   ,:) = -1D5
          fvm(ie)%vtx_cart_physgrid(fv_nphys+1 :fv_nphys+nhc_phys,1-nhc_phys:0     ,1,:) = &
               fvm(ie)%vtx_cart_physgrid(fv_nphys,1,1,2)
          fvm(ie)%vtx_cart_physgrid(fv_nphys+1 :fv_nphys+nhc_phys,1-nhc_phys:0     ,2,:) = &
               fvm(ie)%vtx_cart_physgrid(fv_nphys,1,2,2)
        end if
      end do
      
      !
      ! set vectors for perpendicular flux vector
      !
      rot90_matrix(1,1) = 0; rot90_matrix(2,1) =  1 !counter-clockwise rotation matrix
      rot90_matrix(1,2) =-1; rot90_matrix(2,2) =  0 !counter-clockwise rotation matrix 
      
      iside = 1
      unit_vec(1,iside) = 0 !x-component of displacement vector for side 1
      unit_vec(2,iside) = 1 !y-component of displacement vector for side 1
      
      do iside=2,4
        unit_vec(:,iside) = MATMUL(rot90_matrix(:,:),unit_vec(:,iside-1))
      end do
      
      !
      ! fill halo done
      !
      !-------------------------------
      
      do ie=nets,nete
        do j=imin,imax
          do i=imin,imax
            !
            ! rotate gnomonic coordinate vector
            !    
            ishft = NINT(fvm(ie)%flux_orient_physgrid(2,i,j))
            do ixy=1,2
              !
              ! rotate coordinates if needed through permutation
              !
              fvm(ie)%vtx_cart_physgrid(i,j,ixy,1:4) = cshift(fvm(ie)%vtx_cart_physgrid(i,j,ixy,1:4),shift=ishft)
            end do
          end do
        end do
      end do
      !
      ! pre-compute derived metric terms used for integration, polynomial
      ! evaluation at fvm cell vertices, etc.
      !    
      do ie=nets,nete
        call compute_reconstruct_matrix(fv_nphys,nhe_phys,nhc_phys,irecons,fvm(ie)%dalpha_physgrid,fvm(ie)%dbeta_physgrid,&
             fvm(ie)%spherecentroid_physgrid,fvm(ie)%vtx_cart_physgrid,fvm(ie)%centroid_stretch_physgrid,&
             fvm(ie)%vertex_recons_weights_physgrid,fvm(ie)%recons_metrics_physgrid,fvm(ie)%recons_metrics_integral_physgrid)
      end do      
      !
      ! code specific for physgrid
      !
      !
      ! create a normalized element coordinate system with a halo
      !    
      do ie=nets,nete
        do j=1-nhc_phys,fv_nphys+nhc_phys
          do i=1-nhc_phys,fv_nphys+nhc_phys
            !
            ! only compute for physically existent cells
            !
            if (fvm(ie)%ifct_physgrid(i,j)>0) then
              gnom%x = fvm(ie)%norm_elem_coord_physgrid(1,i,j)
              gnom%y = fvm(ie)%norm_elem_coord_physgrid(2,i,j)
              !
              ! coordinate transform only necessary for points on another panel
              !
              if (NINT(fvm(ie)%flux_orient_physgrid(1,1,1)).NE.NINT(fvm(ie)%flux_orient_physgrid(1,i,j))) then
                tmpcart3d=cubedsphere2cart(gnom,NINT(fvm(ie)%flux_orient_physgrid(1,i,j)))
                tmpgnom=cart2cubedsphere(tmpcart3d,NINT(fvm(ie)%flux_orient_physgrid(1,1,1)))
              else
                tmpgnom%x = fvm(ie)%norm_elem_coord_physgrid(1,i,j)
                tmpgnom%y = fvm(ie)%norm_elem_coord_physgrid(2,i,j)
              end if
              !
              ! convert to element normalized coordinates
              !
              fvm(ie)%norm_elem_coord_physgrid(1,i,j) =(tmpgnom%x-elem(ie)%corners(1)%x)/&
                   (0.5D0*dble(fv_nphys)*fvm(ie)%dalpha_physgrid)-1.0D0
              fvm(ie)%norm_elem_coord_physgrid(2,i,j) =(tmpgnom%y-elem(ie)%corners(1)%y)/&
                   (0.5D0*dble(fv_nphys)*fvm(ie)%dalpha_physgrid)-1.0D0
            else
              fvm(ie)%norm_elem_coord_physgrid(1,i,j) = 1D9
              fvm(ie)%norm_elem_coord_physgrid(2,i,j) = 1D9
            end if
          end do
        end do
      end do
      !
      ! compute Dinv
      !
      do ie=nets,nete
        do j=1-nhc_phys,fv_nphys+nhc_phys
          do i=1-nhc_phys,fv_nphys+nhc_phys
            x1 = fvm(ie)%norm_elem_coord_physgrid(1,i,j)
            x2 = fvm(ie)%norm_elem_coord_physgrid(2,i,j)
            !        x1=-1D0+(i-0.5D0)*dx
            !        x2=-1D0+(j-0.5D0)*dx
            call Dmap(D(i,j,:,:),x1,x2,elem(ie)%corners3D,cubed_sphere_map,elem(ie)%corners,elem(ie)%u2qmap,elem(ie)%facenum)
            detD = D(i,j,1,1)*D(i,j,2,2) - D(i,j,1,2)*D(i,j,2,1)      
            
            fvm(ie)%Dinv_physgrid(i,j,1,1) =  D(i,j,2,2)/detD
            fvm(ie)%Dinv_physgrid(i,j,1,2) = -D(i,j,1,2)/detD
            fvm(ie)%Dinv_physgrid(i,j,2,1) = -D(i,j,2,1)/detD
            fvm(ie)%Dinv_physgrid(i,j,2,2) =  D(i,j,1,1)/detD
          end do
        end do
      end do            
    end if

  end subroutine physgrid_init2


end module fvm_mod
