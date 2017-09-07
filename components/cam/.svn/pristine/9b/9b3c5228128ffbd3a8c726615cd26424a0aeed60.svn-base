#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dimensions_mod
!#ifdef FVM_TRACERS
  use constituents, only : ntrac_d=>pcnst ! _EXTERNAL
!!XXgoldyXX:
  use constituents, only : qsize_d=>pcnst ! _EXTERNAL
!!XXgoldyXX:
!#else
!  use constituents, only : qsize_d=>pcnst ! _EXTERNAL
!#endif

  implicit none
  private

! set MAX number of tracers.  actual number of tracers is a run time argument  
!#ifdef FVM_TRACERS
!!XXgoldyXX:
!  integer, parameter         :: qsize_d = 1        ! SE tracers  
!!XXgoldyXX:
!#else
!  integer, parameter         :: ntrac_d = 0        ! fvm tracers
!#endif



  integer, parameter, public :: nvar = 4 ! FI # dependent variables 
  integer, public            :: qsize_condensate_loading = 1 !how many water variables to include in full density
  logical, parameter, public :: ldry_mass_vertical_coordinates = .true.
  !
  !moist cp in energy conversion term
  !
  ! .false.: force dycore to use cpd (cp dry) instead of moist cp
  ! .true. : use moist cp in dycore
  !
  logical           , public :: lcp_moist = .true. 
  !
  ! phys2dyn_interp_method = 0: bilinear tensor product interpolation
  ! phys2dyn_interp_method = 1: cubic tensor product interpolation (recommended)
  ! phys2dyn_interp_method = 2: cslam reconstruction function evaluation
  ! phys2dyn_interp_method = 3: cslam integration over control volumes constructed over gll points
  !
  integer, parameter, public :: phys2dyn_interp_method = 1
 
  integer, parameter, public :: np = NP
  integer, parameter, public :: nc = 3       !cslam resolution
  integer           , public :: fv_nphys !physics-grid resolution - the "MAX" is so that the code compiles with NC=0

  integer         :: ntrac = 0 !ntrac is set in dyn_comp
  integer         :: qsize = 0 !qsize is set in dyn_comp

  ! fvm dimensions:
  logical, public :: lprint!for debugging
  integer, parameter, public :: ngpc=3          !number of Gausspoints for the fvm integral approximation   !phl change from 4
  integer, parameter, public :: irecons_tracer=6!=1 is PCoM, =3 is PLM, =6 is PPM for tracer reconstruction
  integer, parameter, public :: nhe=1           !Max. Courant number
  integer, parameter, public :: nhr=2           !halo width needed for reconstruction - phl
  integer, parameter, public :: nht=nhe+nhr     !total halo width where reconstruction is needed (nht<=nc) - phl
  integer, parameter, public :: ns=3!quadratic halo interpolation - recommended setting for nc=3
  !nhc determines width of halo exchanged with neighboring elements
  integer, parameter, public :: nhc = nhr+(nhe-1)+(ns-MOD(ns,2))/2
                                                !(different from halo needed for elements on edges and corners


  integer, parameter, public :: kmin_jet=10,kmax_jet=20 !min and max level index for the jet

  integer,  public :: nhc_phys 
  integer,  public :: nhe_phys 
  integer,  public :: nhr_phys 
  integer,  public :: ns_phys  

  integer, public :: npdg = 0  ! dg degree for hybrid cg/dg element  0=disabled 

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1


!  params for a mesh 
!  integer, public, parameter :: max_elements_attached_to_node = 7
!  integer, public, parameter :: s_nv = 2*max_elements_attached_to_node 

  !default for non-refined mesh (note that these are *not* parameters now)
  integer, public  :: max_elements_attached_to_node = 4
  integer, public  :: s_nv = 6
  integer, public  :: max_corner_elem               = 1 !max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8 !4 + 4*max_corner_elem

  public :: qsize,qsize_d,ntrac_d,ntrac

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols



  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    ! new "params"
    max_elements_attached_to_node = 7  ! variable resolution
    s_nv = 2*max_elements_attached_to_node 

    !recalculate these
    max_corner_elem               = max_elements_attached_to_node-3
    max_neigh_edges               = 4 + 4*max_corner_elem


  end subroutine set_mesh_dimensions


end module dimensions_mod

