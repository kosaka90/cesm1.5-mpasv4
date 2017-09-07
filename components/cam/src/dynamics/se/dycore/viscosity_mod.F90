#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
use thread_mod, only : max_num_threads, omp_get_num_threads
use kinds, only : real_kind, iulog
use dimensions_mod, only : np, nc, nlev,qsize,nelemd
use hybrid_mod, only : hybrid_t, get_loop_ranges, config_thread_region
use parallel_mod, only : parallel_t,iam
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit, divergence_sphere
use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
use edge_mod, only : edgevpack, edgerotate, edgevunpack, edgevunpackmin, &
    edgevunpackmax, initEdgeBuffer, FreeEdgeBuffer, edgeSunpackmax, edgeSunpackmin,edgeSpack

use bndry_mod, only : bndry_exchangev, bndry_exchangeS, bndry_exchangeS_start,bndry_exchangeS_finish
use control_mod, only : hypervis_scaling, nu, nu_div
use thread_mod, only : vert_num_threads

implicit none
save

public :: biharmonic_wk
public :: biharmonic_wk_scalar
public :: biharmonic_wk_scalar_minmax
public :: neighbor_minmax, neighbor_minmax_start,neighbor_minmax_finish

!
! compute vorticity/divergence and then project to make continious
! high-level routines uses only for I/O
public :: compute_zeta_C0
public :: compute_div_C0
interface compute_zeta_C0
    module procedure compute_zeta_C0_hybrid       ! hybrid version
    module procedure compute_zeta_C0_par          ! single threaded
end interface
interface compute_div_C0
    module procedure compute_div_C0_hybrid
    module procedure compute_div_C0_par
end interface

public :: compute_zeta_C0_contra    ! for older versions of sweq which carry
public :: compute_div_C0_contra     ! velocity around in contra-coordinates

type (EdgeBuffer_t)          :: edge1

contains

subroutine biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,kbeg,kend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nt,nets,nete,kbeg,kend
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind), dimension(np,np,nets:nete) :: pstens

! local
integer :: k,kptr,i,j,ie,ic,kblk
real (kind=real_kind), dimension(:,:), pointer :: rspheremv
real (kind=real_kind), dimension(np,np) :: lap_ps
real (kind=real_kind), dimension(np,np,nlev) :: T
real (kind=real_kind), dimension(np,np,2) :: v
real (kind=real_kind) ::  nu_ratio1,nu_ratio2
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk

   kblk = kend - kbeg + 1

   var_coef1 = .true.
   if(hypervis_scaling > 0)  var_coef1= .false.

   ! note: there is a scaling bug in the treatment of nu_div
   ! nu_ratio is applied twice, once in each laplace operator
   ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
   ! We should fix this, but it requires adjusting all CAM defaults
   nu_ratio1=1
   nu_ratio2=1
   if (nu_div/=nu) then
      if(hypervis_scaling /= 0) then
         ! we have a problem with the tensor in that we cant seperate
         ! div and curl components.  So we do, with tensor V:
         ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
         nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
         nu_ratio2=1
      else
         nu_ratio1=nu_div/nu
         nu_ratio2=nu_div/nu
      endif
   endif


   do ie=nets,nete
      
      ! should filter lnps + PHI_s/RT?
      pstens(:,:,ie)=laplace_sphere_wk(elem(ie)%state%ps(:,:,nt),deriv,elem(ie),var_coef=var_coef1)
      
      do k=kbeg,kend
         do j=1,np
            do i=1,np
               T(i,j,k)=elem(ie)%state%T(i,j,k,nt) 
            enddo
         enddo
        
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=var_coef1)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
              elem(ie),var_coef=var_coef1,nu_ratio=nu_ratio1)

      enddo

      kptr = kbeg - 1
      call edgeVpack(edge3,ptens(:,:,kbeg:kend,ie),kblk,kptr,ie)

      kptr = (kbeg - 1) + nlev
      call edgeVpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)

      kptr = (kbeg - 1) + 2*nlev
      call edgeVpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)

      kptr = (kbeg - 1) + 3*nlev
      call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie) ! need logic for surface field
   enddo
   
   call bndry_exchangeV(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv => elem(ie)%rspheremp(:,:)
      
      kptr = kbeg - 1
      call edgeVunpack(edge3,ptens(:,:,kbeg:kend,ie),kblk,kptr,ie)

      kptr = (kbeg - 1) + nlev
      call edgeVunpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)
      
      kptr = (kbeg - 1) + 2*nlev
      call edgeVunpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)

      ! apply inverse mass matrix, then apply laplace again
      do k=kbeg,kend
         do j=1,np
            do i=1,np
               T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
               v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
               v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
            enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.true.)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),var_coef=.true.,&
              nu_ratio=nu_ratio2)
      enddo
         
      kptr = (kbeg - 1) + 3*nlev 
      call edgeVunpack(edge3,pstens(:,:,ie),1,kptr,ie) ! need logic for surface field

      ! apply inverse mass matrix, then apply laplace again
      lap_ps(:,:)=rspheremv(:,:)*pstens(:,:,ie)
      pstens(:,:,ie)=laplace_sphere_wk(lap_ps,deriv,elem(ie),var_coef=.true.)

   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


subroutine biharmonic_wk_dp3d(elem,dptens,dpflux,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,kbeg,kend)
  use derivative_mod, only : subcell_Laplace_fluxes
  use dimensions_mod, only : ntrac
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute weak biharmonic operator
  !    input:  h,v (stored in elem()%, in lat-lon coordinates
  !    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer              , intent(in)  :: nt,nets,nete
  integer              , intent(in)  :: kbeg, kend
  real (kind=real_kind), intent(out), dimension(nc,nc,4,nlev,nets:nete) :: dpflux
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens,dptens
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  
  ! local
  integer :: i,j,k,kptr,ie,kblk
  real (kind=real_kind), dimension(:,:), pointer :: rspheremv
  real (kind=real_kind), dimension(np,np) :: tmp
  real (kind=real_kind), dimension(np,np) :: tmp2
  real (kind=real_kind), dimension(np,np,2) :: v
  real (kind=real_kind) :: nu_ratio1, nu_ratio2
  logical var_coef1
  
  kblk = kend - kbeg + 1
  
  if (ntrac>0) dpflux = 0
  !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
  !so tensor is only used on second call to laplace_sphere_wk
  var_coef1 = .true.
  if(hypervis_scaling > 0)    var_coef1 = .false.
  
  ! note: there is a scaling bug in the treatment of nu_div
  ! nu_ratio is applied twice, once in each laplace operator
  ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
  ! We should fix this, but it requires adjusting all CAM defaults
  nu_ratio1=1
  nu_ratio2=1
  if (nu_div/=nu) then
    if(hypervis_scaling /= 0) then
      ! we have a problem with the tensor in that we cant seperate
      ! div and curl components.  So we do, with tensor V:
      ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
      nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
      nu_ratio2=1
    else
      nu_ratio1=nu_div/nu
      nu_ratio2=nu_div/nu
    endif
  endif
  
  
  do ie=nets,nete
    
!$omp parallel do num_threads(vert_num_threads) private(tmp)
    do k=kbeg,kend
      tmp=elem(ie)%state%T(:,:,k,nt) 
      ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
      tmp=elem(ie)%state%dp3d(:,:,k,nt) 
      dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
      vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),&
           var_coef=var_coef1,nu_ratio=nu_ratio1)
    enddo
    
    kptr = kbeg - 1
    call edgeVpack(edge3,ptens(:,:,kbeg:kend,ie),kblk,kptr,ie)
    
    kptr = kbeg - 1 + nlev 
    call edgeVpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)
    
    kptr = kbeg - 1 + 2*nlev 
    call edgeVpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)
    
    kptr = kbeg - 1 + 3*nlev 
    call edgeVpack(edge3,dptens(:,:,kbeg:kend,ie),kblk,kptr,ie)
    
  enddo
  
  call bndry_exchangeV(hybrid,edge3)
  
  do ie=nets,nete
    rspheremv     => elem(ie)%rspheremp(:,:)
    
    kptr = kbeg - 1
    call edgeVunpack(edge3,ptens(:,:,kbeg:kend,ie),kblk,kptr,ie)
    
    kptr = kbeg - 1 + nlev 
    call edgeVunpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)
    
    kptr = kbeg - 1 + 2*nlev 
    call edgeVunpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)
    
    kptr = kbeg - 1 + 3*nlev 
    call edgeVunpack(edge3,dptens(:,:,kbeg:kend,ie),kblk,kptr,ie)
    
    if (ntrac>0) then
      do k=1,nlev
        tmp(:,:)=rspheremv(:,:)*dptens(:,:,k,ie)
        dpflux(:,:,:,k,ie) = subcell_Laplace_fluxes(tmp, deriv, elem(ie), np, nc) 
      enddo
    endif
    
    ! apply inverse mass matrix, then apply laplace again
!$omp parallel do num_threads(vert_num_threads) private(v,tmp,tmp2)
    do k=kbeg,kend
      tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
      ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)
      tmp2(:,:)=rspheremv(:,:)*dptens(:,:,k,ie)
      dptens(:,:,k,ie)=laplace_sphere_wk(tmp2,deriv,elem(ie),var_coef=.true.)
      v(:,:,1)=rspheremv(:,:)*vtens(:,:,1,k,ie)
      v(:,:,2)=rspheremv(:,:)*vtens(:,:,2,k,ie)
      vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),&
           var_coef=.true.,nu_ratio=nu_ratio2)
      
    enddo
  enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk_dp3d


subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,i,j,ie,ic,q
integer :: kbeg,kend,qbeg,qend 
real (kind=real_kind), dimension(np,np) :: lap_p
logical var_coef1
integer :: kblk,qblk   ! The per thead size of the vertical and tracers

  call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.


   kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
   qblk = qend - qbeg + 1   ! calculate size of the block of tracers

   do ie=nets,nete
      do q=qbeg,qend    
         do k=kbeg,kend  !  Potential loop inversion (AAM)
            lap_p(:,:)=qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
         enddo
         kptr = nlev*(q-1) + kbeg - 1
         call edgeVpack(edgeq, qtens(:,:,kbeg:kend,q,ie),kblk,kptr,ie)
      enddo
   enddo


   call bndry_exchangeV(hybrid,edgeq,location='biharmonic_wk_scalar')
   
   do ie=nets,nete

      ! apply inverse mass matrix, then apply laplace again
      do q=qbeg,qend      
        kptr = nlev*(q-1) + kbeg - 1
!JMD        call edgeVunpack(edgeq, qtens(:,:,kbeg:kend,q,ie),kblk,kptr,ie,rank=hybrid%par%rank)
        call edgeVunpack(edgeq, qtens(:,:,kbeg:kend,q,ie),kblk,kptr,ie)
        do k=kbeg,kend    !  Potential loop inversion (AAM)
           lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
        enddo
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk_scalar

subroutine biharmonic_wk_scalar_minmax(elem,qtens,deriv,edgeq,edgeminmax,hybrid,nets,nete,emin,emax)
! NOT yet worked on
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q and Q element min/max
!
!    note: emin/emax must be initialized with Q element min/max.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (EdgeBuffer_t)  , intent(inout) :: edgeminmax
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind), intent(inout), dimension(nlev,qsize,nets:nete) :: emin,emax

! local
integer :: k,kptr,i,j,ie,ic,q
integer :: kbeg,kend,qbeg,qend
real (kind=real_kind), dimension(np,np) :: lap_p
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.

   call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

   do ie=nets,nete
      do q=qbeg,qend    
      do k=kend,kend
         Qmin(:,:,k,q)=emin(k,q,ie)  ! need to set all values in element for
         Qmax(:,:,k,q)=emax(k,q,ie)  ! edgeVpack routine below
         lap_p(:,:) = qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
      enddo
      enddo
      call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,ie)
      call edgeVpack(edgeq,Qmin,nlev*qsize,nlev*qsize,ie)
      call edgeVpack(edgeq,Qmax,nlev*qsize,2*nlev*qsize,ie)
   enddo
   
   call bndry_exchangeV(hybrid,edgeq)
   
   do ie=nets,nete
      do q=qbeg,qend    
      do k=kbeg,kend
         Qmin(:,:,k,q)=emin(k,q,ie)  ! restore element data.  we could avoid
         Qmax(:,:,k,q)=emax(k,q,ie)  ! this by adding a "ie" index to Qmin/max
      enddo
      enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
      call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,ie)
      call edgeVunpackMin(edgeq, Qmin,qsize*nlev,qsize*nlev,ie)
      call edgeVunpackMax(edgeq, Qmax,qsize*nlev,2*qsize*nlev,ie)

      ! apply inverse mass matrix, then apply laplace again
      do q=qbeg,qend    
        do k=kbeg,kend
           lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
           ! note: only need to consider the corners, since the data we packed was
           ! constant within each element
           emin(k,q,ie)=min(qmin(1,1,k,q),qmin(1,np,k,q),qmin(np,1,k,q),qmin(np,np,k,q))
! dont add threshold in this routine - it should be done by calling routine, if needed
!           emin(k,q,ie)=max(emin(k,q,ie),0d0)
           emax(k,q,ie)=max(qmax(1,1,k,q),qmax(1,np,k,q),qmax(np,1,k,q),qmax(np,np,k,q))
        enddo
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine biharmonic_wk_scalar_minmax

subroutine make_C0(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic,kptr,nthread_save


    call initEdgeBuffer(hybrid%par,edge1,elem,nlev)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge1, zeta(1,1,1,ie),nlev,kptr,ie)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,1,ie),nlev,kptr, ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

call FreeEdgeBuffer(edge1) 

end subroutine


subroutine make_C0_vector(v,elem,hybrid,nets,nete)
#if 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS to a velocity vector
! this is a low-performance routine used for I/O and analysis.
! no need to optimize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: v

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge2
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: v1

v1(:,:,:,:) = v(:,:,1,:,:)
call make_C0(v1,elem,hybrid,nets,nete)
v(:,:,1,:,:) = v1(:,:,:,:)

v1(:,:,:,:) = v(:,:,2,:,:)
call make_C0(v1,elem,hybrid,nets,nete)
v(:,:,2,:,:) = v1(:,:,:,:)
#else
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: v

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge2
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: v1



    call initEdgeBuffer(hybrid%par,edge2,elem,2*nlev)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,ie)
enddo
call bndry_exchangeV(hybrid,edge2)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,ie)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

call FreeEdgeBuffer(edge2) 
#endif
end subroutine






subroutine compute_zeta_C0_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine



subroutine compute_div_C0_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_zeta_C0_par(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (parallel_t) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
integer :: nt

! local
type (hybrid_t)              :: hybrid
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

! single thread
hybrid = config_thread_region(par,'serial')
!hybrid = hybrid_create(par,0,1)

call compute_zeta_C0_hybrid(zeta,elem,hybrid,1,nelemd,nt)

end subroutine


subroutine compute_div_C0_par(zeta,elem,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t) :: par
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: zeta
integer :: nt

! local
type (hybrid_t)              :: hybrid
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

! single thread
!hybrid = hybrid_create(par,0,1)
hybrid = config_thread_region(par,'serial')

call compute_div_C0_hybrid(zeta,elem,hybrid,1,nelemd,nt)

end subroutine



subroutine compute_zeta_C0_hybrid(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0_hybrid(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine








subroutine neighbor_minmax(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
 
   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)
   integer :: kblk, qblk
   ! local 
   integer:: ie, q, k, kptr
   integer:: kbeg, kend, qbeg, qend

   call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)
   
   kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
   qblk = qend - qbeg + 1   ! calculate size of the block of tracers
   
   do ie=nets,nete
      do q = qbeg, qend
         kptr = nlev*(q - 1) + kbeg - 1
         call  edgeSpack(edgeMinMax,min_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) + kbeg - 1
         call  edgeSpack(edgeMinMax,max_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
      enddo
   enddo
   
   call bndry_exchangeS(hybrid,edgeMinMax)

   do ie=nets,nete
      do q=qbeg,qend
         kptr = nlev*(q - 1) + kbeg - 1
         call  edgeSunpackMIN(edgeMinMax,min_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) + kbeg - 1
         call  edgeSunpackMAX(edgeMinMax,max_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         do k=kbeg,kend
            min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
         enddo
      enddo
   enddo

end subroutine neighbor_minmax
  

subroutine neighbor_minmax_start(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)

   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)
   integer :: kblk, qblk
   integer :: kbeg, kend, qbeg, qend

   ! local 
   integer :: ie,q, k,kptr

   call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

   kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
   qblk = qend - qbeg + 1   ! calculate size of the block of tracers

   do ie=nets,nete
      do q=qbeg, qend
         kptr = nlev*(q - 1) + kbeg - 1
         call  edgeSpack(edgeMinMax,min_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) + kbeg - 1
         call  edgeSpack(edgeMinMax,max_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
      enddo
   enddo

   call bndry_exchangeS_start(hybrid,edgeMinMax)

end subroutine neighbor_minmax_start

subroutine neighbor_minmax_finish(hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)

   type (hybrid_t)      , intent(in) :: hybrid
   type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
   integer :: nets,nete
   real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
   real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)
   integer :: kblk, qblk
   integer :: ie,q, k,kptr
   integer :: kbeg, kend, qbeg, qend

   call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend,qbeg=qbeg,qend=qend)

   kblk = kend - kbeg + 1   ! calculate size of the block of vertical levels
   qblk = qend - qbeg + 1   ! calculate size of the block of tracers

   call bndry_exchangeS_finish(hybrid,edgeMinMax)

   do ie=nets,nete
      do q=qbeg, qend
         kptr = nlev*(q - 1) + kbeg - 1
         call  edgeSunpackMIN(edgeMinMax,min_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         kptr = qsize*nlev + nlev*(q - 1) + kbeg - 1
         call  edgeSunpackMAX(edgeMinMax,max_neigh(kbeg:kend,q,ie),kblk,kptr,ie)
         do k=kbeg,kend
            min_neigh(k,q,ie) = max(min_neigh(k,q,ie),0d0)
         enddo
      enddo
   enddo

end subroutine neighbor_minmax_finish

end module viscosity_mod
