module spmd_dyn

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: SPMD implementation of CAM mpas dynamics.
   ! 
   !-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use spmd_utils, only : iam, masterproc, npes
   use cam_abortutils, only : endrun

   implicit none
   private

   integer, allocatable, public :: proc(:)

   public compute_gsfactors, spmdinit_dyn, spmdbuf

   logical,public :: local_dp_map=.true.    ! flag indicates that mapping between dynamics 
   !  and physics decompositions does not require 
   !  interprocess communication
#ifdef SPMD
   integer, allocatable, public :: proc_smp_map(:)
   integer, allocatable, public :: nlat_p(:)
   integer,public :: nsmps
   
   integer, public :: block_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                      !  in dynamics decomposition (including level 0)
   integer, public :: chunk_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                      !  in physics decomposition (including level 0)
#endif

CONTAINS

!========================================================================

subroutine spmdinit_dyn ()

   return

end subroutine spmdinit_dyn

!========================================================================

subroutine spmdbuf

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: placeholder for buffer allocation routine 
   ! 
   ! Method: 
   ! 
   ! Author: CCM Core Group
   ! 
   !-----------------------------------------------------------------------
   implicit none

   return

end subroutine spmdbuf

!========================================================================

subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: Compute arguments for gatherv, scatterv
   ! 
   ! Author: CCM Core Group
   ! 
   !-----------------------------------------------------------------------
   !
   ! Input arguments
   !
   integer, intent(in) :: numperlat    ! number of elements per latitude
   !
   ! Output arguments
   !
   integer, intent(out) :: numtot               ! total number of elements (to send or recv)
   integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
   integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements

end subroutine compute_gsfactors

!========================================================================

end module spmd_dyn
