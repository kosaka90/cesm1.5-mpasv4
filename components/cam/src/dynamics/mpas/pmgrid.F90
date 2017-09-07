module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: Parameters and variables related to the dynamics grid
! 
! Author: Art Mirin, Michael Duda
! 
! PLON and PLAT do not correspond to the number of latitudes and longitudes in
! this version of dynamics. 
! 
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8

   integer, parameter :: plon   = 1                     ! number of longitudes
   integer, parameter :: plev   = PLEV                  ! number of vertical levels
   integer, parameter :: plat   = 1                     ! number of latitudes

   integer :: nEdges                                    ! number of edges
   integer :: nVertices                                 ! number of vertices
   integer :: numcols                                   ! local number of cells (columns)
   integer :: maxEdges                                  ! maximum number of edges per cell

   integer, parameter :: plevp  = plev + 1              ! plev + 1

   real(kind=r8), parameter :: MPAS_SPHERE_RAD = 6371229.0_r8   ! Radius of spherical earth in MPAS-A

   !logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set

end module pmgrid
