module dycore
!
! Data and utility routines related to the dycore
!
   implicit none

PRIVATE

   public :: dycore_is, get_resolution

CONTAINS

   logical function dycore_is (name)
!
! Input arguments
!
      character(len=*) :: name
      
      dycore_is = .false.
      if (name == 'unstructured' .or. name == 'UNSTRUCTURED' .or. &
           name == 'mpas' .or. name == 'MPAS') then
         dycore_is = .true.
      end if
      
      return
   end function dycore_is

   character(len=7) function get_resolution()

      ! Deprecated -- remove as soon as external references have been removed
      get_resolution = 'UNKNOWN'

   end function get_resolution

end module dycore


