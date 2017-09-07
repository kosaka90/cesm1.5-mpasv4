#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_CONTROL_VOLUME_MOD---------------------------------------------CE-for FVM
! AUTHOR: Christoph Erath, 11.June 2011                                             !
! This module contains everything to initialize the arrival. It also provides the   !
! interpolation points for the reconstruction (projection from one face to another  !
! when the element is on the cube edge)                                             !
! It also intialize the start values, see also fvm_analytic                         !
!-----------------------------------------------------------------------------------! 
module fvm_control_volume_mod
  ! ---------------------------------------------------------------------------------
  use kinds, only : real_kind, int_kind
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2D_t
  ! ---------------------------------------------------------------------------------
  use element_mod, only: element_t
  ! ---------------------------------------------------------------------------------
  use dimensions_mod, only: nc, nhe, nlev, ntrac, ntrac_d, qsize_d,ne, np, nhr, ns, nhc
  use dimensions_mod, only: fv_nphys, nhe_phys, nhr_phys, ns_phys, nhc_phys,fv_nphys
  use dimensions_mod, only: irecons_tracer, qsize_condensate_loading,qsize
  ! ---------------------------------------------------------------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  ! ---------------------------------------------------------------------------------

  use parallel_mod, only : abortmp


  implicit none
  private
  integer, parameter, private:: nh = nhr+(nhe-1) ! = 2 (nhr=2; nhe=1)
                                                 ! = 3 (nhr=2; nhe=2)

  type, public :: fvm_struct
    ! fvm tracer mixing ratio: (kg/kg)
    real (kind=real_kind) :: c(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev,ntrac_d,2)
    real (kind=real_kind) :: se_flux  (1-nhe:nc+nhe,1-nhe:nc+nhe,   4,nlev) 
    real (kind=real_kind) :: mass(ntrac_d+3) !total tracer mass, fvm air mass (psc,dp_dvm) and se air mass - for diagnostics only

    real (kind=real_kind) :: dp_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,nlev        ,2)
    real (kind=real_kind) :: dp_ref(nlev)
    real (kind=real_kind) :: dp_ref_inverse(nlev)
    real (kind=real_kind) :: psc(1-nhc:nc+nhc,1-nc:nc+nhc)

    real (kind=real_kind)    :: inv_area_sphere(nc,nc)    ! inverse area_sphere    
    real (kind=real_kind)    :: inv_se_area_sphere(nc,nc) ! inverse area_sphere    

    integer                  :: faceno         !face number
    ! number of south,....,swest and 0 for interior element 
    integer                  :: cubeboundary                                                 


    real (kind=real_kind) :: displ_max(1-nhc:nc+nhc,1-nhc:nc+nhc,4)
    integer               :: flux_vec   (2,1-nhc:nc+nhc,1-nhc:nc+nhc,4) 

    real (kind=real_kind)    :: maxcfl(2,nlev)


    !
    ! cartesian location of vertices for flux sides
    !
    !  x-coordinate of vertex 1: vtx_cart(i,j,1,1) = fvm%acartx(i)
    !  y-coordinate of vertex 1: vtx_cart(i,j,2,1) = fvm%acarty(j)
    !
    !  x-coordinate of vertex 2: vtx_cart(i,j,1,2) = fvm%acartx(i+1)
    !  y-coordinate of vertex 2: vtx_cart(i,j,2,2) = fvm%acarty(j  )
    !
    !  x-coordinate of vertex 3: vtx_cart(i,j,1,3) = fvm%acartx(i+1)
    !  y-coordinate of vertex 3: vtx_cart(i,j,2,3) = fvm%acarty(j+1)
    !
    !  x-coordinate of vertex 4: vtx_cart(i,j,1,4) = fvm%acartx(i  )
    !  y-coordinate of vertex 4: vtx_cart(i,j,2,4) = fvm%acarty(j+1)
    !
    real (kind=real_kind) :: vtx_cart (1-nhc:nc+nhc,1-nhc:nc+nhc,2,4)
    !
    ! flux_orient(1,i,j) = panel on which control volume (i,j) is located
    ! flux_orient(2,i,j) = cshift value for vertex permutation
    !
    real (kind=real_kind) :: flux_orient(2  ,1-nhc:nc+nhc,1-nhc:nc+nhc) 
    !
    ! i,j: indicator function for non-existent cells (0 for corner halo and 1 elsewhere)
    !
    integer                  :: ifct   (1-nhc:nc+nhc,1-nhc:nc+nhc) 
    integer                  :: rot_matrix(2,2,1-nhc:nc+nhc,1-nhc:nc+nhc)
    !    
    real (kind=real_kind)    :: dalpha, dbeta             ! central-angle for gnomonic coordinates
    type (spherical_polar_t) :: center_cart(nc,nc)        ! center of fvm cell in gnomonic coordinates
    real (kind=real_kind)    :: area_sphere(nc,nc)        ! spherical area of fvm cell
    real (kind=real_kind)    :: spherecentroid(1-nhc:nc+nhc,1-nhc:nc+nhc,irecons_tracer-1) ! centroids
    !
    ! pre-computed metric terms (for efficiency)
    !
    ! recons_metrics(:,:,1) = spherecentroid(:,:,1)**2 -spherecentroid(:,:,3)
    ! recons_metrics(:,:,2) = spherecentroid(:,:,2)**2 -spherecentroid(:,:,4)
    ! recons_metrics(:,:,3) = spherecentroid(:,:,1)*spherecentroid(:,:,2)-spherecentroid(:,:,5)
    !
    real (kind=real_kind)    :: recons_metrics(1-nhe:nc+nhe,1-nhe:nc+nhe,3)    
    !
    ! recons_metrics_integral(:,:,1) = 2.0D0*spherecentroid(:,:,1)**2 -spherecentroid(:,:,3)
    ! recons_metrics_integral(:,:,2) = 2.0D0*spherecentroid(:,:,2)**2 -spherecentroid(:,:,4)
    ! recons_metrics_integral(:,:,3) = 2.0D0*spherecentroid(:,:,1)*spherecentroid(:,:,2)-spherecentroid(:,:,5)
    !
    real (kind=real_kind)    :: recons_metrics_integral(1-nhe:nc+nhe,1-nhe:nc+nhe,3)    
    integer                  :: jx_min(3), jx_max(3), jy_min(3), jy_max(3) !bounds for computation

    ! provide fixed interpolation points with respect to the arrival grid for 
    ! reconstruction   
    integer                  :: ibase(1-nh:nc+nh,1:nhr,2)  
    real (kind=real_kind)    :: halo_interp_weight(1:ns,1-nh:nc+nh,1:nhr,2)
    real (kind=real_kind)    :: centroid_stretch(1-nhe:nc+nhe,1-nhe:nc+nhe,1:7) !for finite-difference reconstruction
    !
    ! pre-compute weights for reconstruction at cell vertices
    !
    !  ! Evaluate constant order terms
    !  value = fcube(a,b) + &
    !  ! Evaluate linear order terms
    !          recons(a,b,1) * (cartx - centroid(a,b,1)) + &
    !          recons(a,b,2) * (carty - centroid(a,b,2)) + &
    !  ! Evaluate second order terms
    !          recons(a,b,3) * (centroid(a,b,1)**2 - centroid(a,b,3)) + &
    !          recons(a,b,4) * (centroid(a,b,2)**2 - centroid(a,b,4)) + &
    !          recons(a,b,5) * (centroid(a,b,1) * centroid(a,b,2) - centroid(a,b,5)) + &
    !
    !          recons(a,b,3) * (cartx - centroid(a,b,1))**2 + &
    !          recons(a,b,4) * (carty - centroid(a,b,2))**2 + &
    !          recons(a,b,5) * (cartx - centroid(a,b,1)) * (carty - centroid(a,b,2))
    !   
    real (kind=real_kind)    :: vertex_recons_weights(4,1:irecons_tracer-1,1-nhe:nc+nhe,1-nhe:nc+nhe)
    !
    ! for mapping fvm2dyn
    !
    real (kind=real_kind)    :: norm_elem_coord(2,1-nhc:nc+nhc,1-nhc:nc+nhc)
    !
    !******************************************
    !
    ! separate physics grid variables
    !
    !******************************************
    !
    real (kind=real_kind)    , allocatable :: vtx_cart_physgrid(:,:,:,:)
    real (kind=real_kind)    , allocatable :: flux_orient_physgrid(:,:,:)
    integer                  , allocatable :: ifct_physgrid(:,:)
    integer                  , allocatable :: rot_matrix_physgrid(:,:,:,:)
    real (kind=real_kind)    , allocatable :: spherecentroid_physgrid(:,:,:)
    real (kind=real_kind)    , allocatable :: recons_metrics_physgrid(:,:,:)
    real (kind=real_kind)    , allocatable :: recons_metrics_integral_physgrid(:,:,:)
    ! centroid_stretch_physgrid for finite-difference reconstruction
    real (kind=real_kind)    , allocatable :: centroid_stretch_physgrid       (:,:,:)
    real (kind=real_kind)    :: dalpha_physgrid, dbeta_physgrid             ! central-angle for gnomonic coordinates
    type (spherical_polar_t) , allocatable :: center_cart_physgrid(:,:)        ! center of fvm cell in gnomonic coordinates
    real (kind=real_kind)    , allocatable :: area_sphere_physgrid(:,:)        ! spherical area of fvm cell
    integer                  :: jx_min_physgrid(3), jx_max_physgrid(3) !bounds for computation
    integer                  :: jy_min_physgrid(3), jy_max_physgrid(3) !bounds for computation
    integer                  , allocatable :: ibase_physgrid(:,:,:)
    real (kind=real_kind)    , allocatable :: halo_interp_weight_physgrid(:,:,:,:)
    real (kind=real_kind)    , allocatable :: vertex_recons_weights_physgrid(:,:,:,:)

    real (kind=real_kind), allocatable :: norm_elem_coord_physgrid(:,:,:)
    real (kind=real_kind)    , allocatable :: Dinv_physgrid(:,:,:,:)

    real (kind=real_kind) , allocatable :: fc(:,:,:,:)
    real (kind=real_kind) , allocatable :: fc_phys(:,:,:,:)
    real (kind=real_kind) , allocatable :: ft(:,:,:)
    real (kind=real_kind) , allocatable :: fm(:,:,:,:)
    real (kind=real_kind) , allocatable :: dp_phys(:,:,:)
  end type fvm_struct

  public :: fvm_mesh, fvm_set_cubeboundary, allocate_physgrid_vars

  
  real (kind=real_kind),parameter, public   :: bignum = 1.0D20

  integer, public :: n0_fvm, np1_fvm !fvm time-levels
  integer, parameter, public :: fvm_supercycling = 3

contains
  subroutine fvm_set_cubeboundary(elem, fvm)
    implicit none
    type (element_t) , intent(in)      :: elem
    type (fvm_struct), intent(inout)   :: fvm
    
    integer                              :: i,j
    logical                              :: corner
    integer                              :: mynbr_cnt, cnt, mystart, start
    integer                              :: nbrsface(8)    ! store the neighbours in north, south 
        
    fvm%faceno=elem%FaceNum
    ! write the neighbors in the structure
    fvm%cubeboundary=0
    corner=.FALSE.
    do j=1,8
       mynbr_cnt = elem%vertex%nbrs_ptr(j+1) - elem%vertex%nbrs_ptr(j) !length of neighbor location  
       mystart = elem%vertex%nbrs_ptr(j) 
       !NOTE: assuming that we do not have multiple corner neighbors (so not a refined mesh)
       if (mynbr_cnt > 0 ) then
          nbrsface(j)=elem%vertex%nbrs_face(mystart)
          ! note that if the element lies on a corner, it will be at j=5,6,7,8
          if ((nbrsface(j) /= fvm%faceno) .AND. (j<5)) then
             fvm%cubeboundary=j
          endif
       else   ! corner on the cube
          if (.NOT. corner) then
             nbrsface(j)=-1
             fvm%cubeboundary=j
             corner=.TRUE.
          else
             if ( ne == 0 ) then
                ! dont check this condition.  note that we call these code
                ! generate phys grid template files, so we need to be able
                ! to call create_ari() to create the subcells even though
                ! cslam cant run with the unstructed ne=0 case
             else
                print *,'Error in fvm_CONTROL_VOLUME_MOD - Subroutine fvm_MESH_ARI: '
                call abortmp('Do not allow one element per face for fvm, please increase ne!')
             endif
          endif
       end if
    end do
  end subroutine fvm_set_cubeboundary

  subroutine fvm_mesh(elem, fvm)
    use fvm_analytic_mod, only : compute_halo_vars
    use fvm_analytic_mod, only : create_interpolation_points
    use derivative_mod  , only : subcell_integration

    implicit none
    type (element_t), intent(in)     :: elem
    type (fvm_struct), intent(inout) :: fvm
    integer :: i,j
    real (kind=real_kind)            :: tmp(np,np)
    !
    ! initialize metric and related terms on panel
    !    
    call compute_halo_vars(elem,&    !input
         fvm%faceno,fvm%cubeboundary,nc,nhc,nhe,np,&  !input
         fvm%jx_min,fvm%jx_max,fvm%jy_min,fvm%jy_max,&!output
         fvm%flux_orient,fvm%ifct,fvm%rot_matrix)     !output
    do j=1,nc
      do i=1,nc
        fvm%norm_elem_coord(1,i,j) = elem%corners(1)%x+(i-0.5D0)*fvm%dalpha
        fvm%norm_elem_coord(2,i,j) = elem%corners(1)%y+(j-0.5D0)*fvm%dalpha
      end do
    end do

    !
    ! overwrite areas for consistency with SE areas (that are O(10E-5) incorrect)
    !
    tmp = 1.0D0
    fvm%area_sphere = subcell_integration(tmp, np, nc, elem%metdet)
    !
    ! do the same for physics grid
    !
    call compute_halo_vars(elem,&
         fvm%faceno,fvm%cubeboundary,fv_nphys,nhc_phys,nhe_phys,np,&
         fvm%jx_min_physgrid,fvm%jx_max_physgrid,fvm%jy_min_physgrid,fvm%jy_max_physgrid,&
         fvm%flux_orient_physgrid,fvm%ifct_physgrid,fvm%rot_matrix_physgrid)
    do j=1,fv_nphys
      do i=1,fv_nphys
        fvm%norm_elem_coord_physgrid(1,i,j) = elem%corners(1)%x+(i-0.5D0)*fvm%dalpha_physgrid
        fvm%norm_elem_coord_physgrid(2,i,j) = elem%corners(1)%y+(j-0.5D0)*fvm%dalpha_physgrid
      end do
    end do
    !
    ! initialize halo interpolation variables
    !
    call create_interpolation_points(elem,&
         nc,nhc,nhr,ns,nh,fvm%cubeboundary,&
         fvm%dalpha,fvm%dbeta,fvm%ibase,fvm%halo_interp_weight)
    call create_interpolation_points(elem,&
         fv_nphys,nhc_phys,nhr_phys,ns_phys,nhr_phys,fvm%cubeboundary,&
         fvm%dalpha_physgrid,fvm%dbeta_physgrid,fvm%ibase_physgrid,fvm%halo_interp_weight_physgrid)
  end subroutine fvm_mesh


  subroutine allocate_physgrid_vars(fvm,par)
    use cam_logfile   , only : iulog
    use parallel_mod  , only : parallel_t
    use dimensions_mod, only : nelemd
    type (fvm_struct), intent(inout) :: fvm(:)
    type (parallel_t), intent(in)    :: par
    integer :: ie

    nhc_phys = fv_nphys
    nhe_phys = 0
    nhr_phys = 2
    ns_phys  = MAX(fv_nphys,2)

    if(par%masterproc) then
      write(iulog,*)"allocating physgrid grid vars"
      write(iulog,*)"fv_nphys,nhc_phys,nhe_phys,nhr_phys,ns_phys = ",&
           fv_nphys,nhc_phys,nhe_phys,nhr_phys,ns_phys
    end if

    do ie=1,nelemd
      allocate(fvm(ie)%vtx_cart_physgrid      (1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,4))
      allocate(fvm(ie)%flux_orient_physgrid   (2,1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys))
      allocate(fvm(ie)%ifct_physgrid         (1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys))
      allocate(fvm(ie)%rot_matrix_physgrid   (2,2,1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys))
      allocate(fvm(ie)%spherecentroid_physgrid(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys, &
           irecons_tracer-1)) ! centroids
      allocate(fvm(ie)%recons_metrics_physgrid         (1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,3))
      allocate(fvm(ie)%recons_metrics_integral_physgrid(1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,3))
      ! centroid_stretch_physgrid for finite-difference reconstruction
      allocate(fvm(ie)%centroid_stretch_physgrid       (1-nhe_phys:fv_nphys+nhe_phys,1-nhe_phys:fv_nphys+nhe_phys,1:7))
      allocate(fvm(ie)%center_cart_physgrid(fv_nphys,fv_nphys))
      allocate(fvm(ie)%area_sphere_physgrid(fv_nphys,fv_nphys))       ! spherical area of fvm cell
      allocate(fvm(ie)%ibase_physgrid(1-nhr_phys:fv_nphys+nhr_phys,1:nhr_phys,2))
      allocate(fvm(ie)%halo_interp_weight_physgrid(1:ns_phys,1-nhr_phys:fv_nphys+nhr_phys,1:nhr_phys,2))
      allocate(fvm(ie)%vertex_recons_weights_physgrid(4,1:irecons_tracer-1,1-nhe_phys:fv_nphys+nhe_phys,&
           1-nhe_phys:fv_nphys+nhe_phys))
      
      allocate(fvm(ie)%norm_elem_coord_physgrid(2,1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys    ))
      allocate(fvm(ie)%Dinv_physgrid           (  1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,2))
      
      allocate(fvm(ie)%fc(nc,nc,nlev,max(ntrac_d,qsize_d)))
      allocate(fvm(ie)%fc_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev,max(ntrac_d,qsize_d)))
      allocate(fvm(ie)%ft(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev))
      allocate(fvm(ie)%fm(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,2,nlev))
      allocate(fvm(ie)%dp_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,nlev))
    end do
  end subroutine allocate_physgrid_vars
end module fvm_control_volume_mod
