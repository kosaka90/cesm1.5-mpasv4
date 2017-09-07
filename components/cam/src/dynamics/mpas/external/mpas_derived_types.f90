! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_derived_types
!
!> \brief   MPAS Derived data types
!> \author  Doug Jacobsen, Michael Duda
!> \date    03/10/2015
!> \details 
!> This module defines derived data types related to fields, and variable structures.
!> It also includes routines for allocating and deallocating these types.
!
!-----------------------------------------------------------------------
module mpas_derived_types

   use mpas_kind_types

   use pio
   use pio_types

   use ESMF
   use ESMF_BaseMod
   use ESMF_Stubs
   use ESMF_CalendarMod
   use ESMF_ClockMod
   use ESMF_TimeMod
   use ESMF_TimeIntervalMod

   integer, parameter :: MPAS_ATT_INT   = 1
   integer, parameter :: MPAS_ATT_INTA  = 2
   integer, parameter :: MPAS_ATT_REAL  = 3
   integer, parameter :: MPAS_ATT_REALA = 4
   integer, parameter :: MPAS_ATT_TEXT  = 5

   ! Derived type for holding field attributes
   type att_list_type
      character (len=StrKIND) :: attName
      integer :: attType
      integer :: attValueInt
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=StrKIND) :: attValueText
      type (att_list_type), pointer :: next => null()
   end type att_list_type

   integer, parameter :: TABLESIZE=27183     !< Number of spaces in the table (the number of linked lists)

   type hashnode
      integer :: key
      type (hashnode), pointer :: next
   end type hashnode
 
   type hashnode_ptr
      type (hashnode), pointer :: p        !< Pointer to a list of entries
   end type hashnode_ptr
 
   type hashtable
      integer :: size
      type (hashnode_ptr), dimension(TABLESIZE) :: table !< The hashtable array
   end type hashtable

   type dm_info
     integer :: nprocs, my_proc_id, comm, info
     logical :: using_external_comm
   end type dm_info


   type mpas_exchange_list
     integer :: endPointID
     integer :: nlist
     integer, dimension(:), pointer :: srcList
     integer, dimension(:), pointer :: destList
     type (mpas_exchange_list), pointer :: next
   end type mpas_exchange_list


   type mpas_exchange_list_pointer
     type (mpas_exchange_list), pointer :: exchList
   end type mpas_exchange_list_pointer


   type mpas_multihalo_exchange_list
     type (mpas_exchange_list_pointer), dimension(:), pointer :: halos
     ! Pointers to the mulithalo exchange lists for this variable on the prev and next blocks on this processor
     type (mpas_multihalo_exchange_list), pointer :: prev, next
   end type mpas_multihalo_exchange_list


   type mpas_communication_list
     integer :: procID
     integer :: nlist
     real (kind=RKIND), dimension(:), pointer :: rbuffer
     integer, dimension(:), pointer :: ibuffer
     integer :: reqID
     type (mpas_communication_list), pointer :: next
   end type mpas_communication_list

   integer, parameter :: MPAS_MISSING_DIM = -999


   integer, parameter :: MPAS_DECOMP_NONDECOMP = 1013, &
                         MPAS_DECOMP_CELLS     = 1014, &
                         MPAS_DECOMP_EDGES     = 1015, &
                         MPAS_DECOMP_VERTICES  = 1016

   ! Derived type for storing fields
   type field5DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(5) :: dimNames
      integer, dimension(5) :: dimSizes
      real (kind=RKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field5DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field5DReal


   ! Derived type for storing fields
   type field4DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(4) :: dimNames
      integer, dimension(4) :: dimSizes
      real (kind=RKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field4DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field4DReal



   ! Derived type for storing fields
   type field3DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer, dimension(3) :: dimSizes
      real (kind=RKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DReal


   ! Derived type for storing fields
   type field2DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer, dimension(2) :: dimSizes
      real (kind=RKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DReal


   ! Derived type for storing fields
   type field1DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      real (kind=RKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DReal


   ! Derived type for storing fields
   type field0DReal
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      real (kind=RKIND) :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      real (kind=RKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DReal), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DReal


   ! Derived type for storing fields
   type field3DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(3) :: dimNames
      integer :: defaultValue
      integer, dimension(3) :: dimSizes
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field3DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field3DInteger


   ! Derived type for storing fields
   type field2DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:,:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(2) :: dimNames
      integer :: defaultValue
      integer, dimension(2) :: dimSizes
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field2DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field2DInteger


   ! Derived type for storing fields
   type field1DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer, dimension(:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer :: defaultValue
      integer, dimension(1) :: dimSizes
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DInteger


   ! Derived type for storing fields
   type field0DInteger
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      integer :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      integer :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DInteger), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DInteger


   ! Derived type for storing fields
   type field1DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND), dimension(:), pointer :: array

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND), dimension(1) :: dimNames
      integer, dimension(1) :: dimSizes
      character (len=StrKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      logical :: isPersistent
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field1DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field1DChar


   ! Derived type for storing fields
   type field0DChar
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      character (len=StrKIND) :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      character (len=StrKIND) :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DChar), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DChar


   ! Derived type for storing fields
   type field0DLogical
  
      ! Back-pointer to the containing block
      type (block_type), pointer :: block

      ! Raw array holding field data on this block
      logical :: scalar

      ! Information used by the I/O layer
      character (len=StrKIND) :: fieldName
      character (len=StrKIND), dimension(:), pointer :: constituentNames => null()
      logical :: defaultValue
      logical :: isDecomposed
      logical :: hasTimeDimension
      logical :: isActive
      logical :: isVarArray
      type (att_list_type), pointer :: attList => null()     

      ! Pointers to the prev and next blocks for this field on this task
      type (field0DLogical), pointer :: prev, next

      ! Halo communication lists
      type (mpas_multihalo_exchange_list), pointer :: sendList
      type (mpas_multihalo_exchange_list), pointer :: recvList
      type (mpas_multihalo_exchange_list), pointer :: copyList
   end type field0DLogical


   integer, parameter :: MPAS_POOL_TABLE_SIZE = 128

   integer, parameter :: MPAS_POOL_SILENT = 1001, &
                         MPAS_POOL_WARN   = 1002, &
                         MPAS_POOL_FATAL  = 1003

   integer, parameter :: MPAS_POOL_FIELD     = 1004, &
                         MPAS_POOL_CONFIG    = 1005, &
                         MPAS_POOL_DIMENSION = 1006, &
                         MPAS_POOL_SUBPOOL   = 1007, &
                         MPAS_POOL_PACKAGE   = 1008

   integer, parameter :: MPAS_POOL_REAL      = 1009, &
                         MPAS_POOL_INTEGER   = 1010, &
                         MPAS_POOL_LOGICAL   = 1011, &
                         MPAS_POOL_CHARACTER = 1012

   type mpas_pool_data_type
      integer :: contentsType
      integer :: contentsDims
      integer :: contentsTimeLevs

      ! For storing fields
      type (field0DReal), pointer :: r0 => null()
      type (field1DReal), pointer :: r1 => null()
      type (field2DReal), pointer :: r2 => null()
      type (field3DReal), pointer :: r3 => null()
      type (field4DReal), pointer :: r4 => null()
      type (field5DReal), pointer :: r5 => null()
      type (field0DReal), dimension(:), pointer :: r0a => null()
      type (field1DReal), dimension(:), pointer :: r1a => null()
      type (field2DReal), dimension(:), pointer :: r2a => null()
      type (field3DReal), dimension(:), pointer :: r3a => null()
      type (field4DReal), dimension(:), pointer :: r4a => null()
      type (field5DReal), dimension(:), pointer :: r5a => null()
      type (field0DInteger), pointer :: i0 => null()
      type (field1DInteger), pointer :: i1 => null()
      type (field2DInteger), pointer :: i2 => null()
      type (field3DInteger), pointer :: i3 => null()
      type (field0DInteger), dimension(:), pointer :: i0a => null()
      type (field1DInteger), dimension(:), pointer :: i1a => null()
      type (field2DInteger), dimension(:), pointer :: i2a => null()
      type (field3DInteger), dimension(:), pointer :: i3a => null()
      type (field0DChar), pointer :: c0 => null()
      type (field1DChar), pointer :: c1 => null()
      type (field0DChar), dimension(:), pointer :: c0a => null()
      type (field1DChar), dimension(:), pointer :: c1a => null()
      type (field0DLogical), pointer :: l0 => null()
      type (field0DLogical), dimension(:), pointer :: l0a => null()
      type (mpas_pool_type), pointer :: p => null()
 
      ! For storing config options, dimensions, and packages
      integer, pointer :: simple_int => null()
      integer, dimension(:), pointer :: simple_int_arr => null()
      real(kind=RKIND), pointer :: simple_real => null()
      logical, pointer :: simple_logical => null()
      character(len=StrKIND), pointer :: simple_char => null()
   end type mpas_pool_data_type

   type mpas_pool_member_type
      character (len=StrKIND) :: key
      integer :: keyLen
      integer :: contentsType
      type (mpas_pool_data_type), pointer :: data => null()
      type (mpas_pool_member_type), pointer :: next => null()
      type (mpas_pool_member_type), pointer :: iteration_next => null()
      type (mpas_pool_member_type), pointer :: iteration_prev => null()
   end type mpas_pool_member_type

   type mpas_pool_head_type
      type (mpas_pool_member_type), pointer :: head => null()
   end type mpas_pool_head_type

   type mpas_pool_type
      integer :: size
      type (mpas_pool_head_type), dimension(:), pointer :: table => null()
      type (mpas_pool_member_type), pointer :: iterator => null()
      type (mpas_pool_member_type), pointer :: iteration_head => null()
      type (mpas_pool_member_type), pointer :: iteration_tail => null()
   end type mpas_pool_type

   type mpas_pool_iterator_type
      character (len=StrKIND) :: memberName
      integer :: memberType
      integer :: dataType
      integer :: nDims
      integer :: nTimeLevels
   end type mpas_pool_iterator_type

   type mpas_pool_field_info_type
      integer :: fieldType
      integer :: nDims
      integer :: nTimeLevels
      logical :: isActive
   end type mpas_pool_field_info_type

   integer, parameter :: MPAS_IO_OFFSET_KIND = PIO_OFFSET_KIND

   ! File access modes
   integer, parameter :: MPAS_IO_READ  = 1, &
                         MPAS_IO_WRITE = 2

   ! I/O formats
   integer, parameter :: MPAS_IO_NETCDF   = 3, &
                         MPAS_IO_PNETCDF  = 4, &
                         MPAS_IO_NETCDF4  = 5, &
                         MPAS_IO_PNETCDF5 = 6

   ! Field and attribute types
   integer, parameter :: MPAS_IO_REAL     =  7,  &
                         MPAS_IO_DOUBLE   =  8,  &
                         MPAS_IO_INT      =  9,  &
                         MPAS_IO_LOGICAL  = 10,  &
                         MPAS_IO_CHAR     = 11

   ! Field precision
   integer, parameter :: MPAS_IO_SINGLE_PRECISION = 12, &
                         MPAS_IO_DOUBLE_PRECISION = 13, &
                         MPAS_IO_NATIVE_PRECISION = MPAS_IO_DOUBLE_PRECISION

   ! Unlimited / record dimension
   integer, parameter :: MPAS_IO_UNLIMITED_DIM = -123456

   ! Error codes
   integer, parameter :: MPAS_IO_NOERR              =  0, &
                         MPAS_IO_ERR_INVALID_MODE   = -1, &
                         MPAS_IO_ERR_INVALID_FORMAT = -2, &
                         MPAS_IO_ERR_LONG_FILENAME  = -3, &
                         MPAS_IO_ERR_UNINIT_HANDLE  = -4, &
                         MPAS_IO_ERR_PIO            = -5, &
                         MPAS_IO_ERR_DATA_MODE      = -6, &
                         MPAS_IO_ERR_NOWRITE        = -7, &
                         MPAS_IO_ERR_REDEF_DIM      = -8, &
                         MPAS_IO_ERR_REDEF_VAR      = -9, &
                         MPAS_IO_ERR_UNDEFINED_DIM  = -10, &
                         MPAS_IO_ERR_UNDEFINED_VAR  = -11, &
                         MPAS_IO_ERR_REDEF_ATT      = -12, &
                         MPAS_IO_ERR_WRONG_ATT_TYPE = -13, &
                         MPAS_IO_ERR_NO_DECOMP      = -14, &
                         MPAS_IO_ERR_TWO_UNLIMITED_DIMS = -15, &
                         MPAS_IO_ERR_WRONG_MODE         = -16, &
                         MPAS_IO_ERR_NO_UNLIMITED_DIM   = -17, &
                         MPAS_IO_ERR_UNIMPLEMENTED      = -18, &
                         MPAS_IO_ERR_WOULD_CLOBBER      = -19, &
                         MPAS_IO_ERR_NOEXIST_READ       = -20, &
                         MPAS_IO_ERR_MISSING_DIM        = -21

   type MPAS_IO_Handle_type
      logical :: initialized = .false.
      logical :: preexisting_file = .false.
      logical :: data_mode = .false.
      type (file_desc_t) :: pio_file
      character (len=StrKIND) :: filename
      integer :: iomode
      integer :: ioformat
      integer :: pio_unlimited_dimid
      integer :: preexisting_records = 0
      integer (kind=MPAS_IO_OFFSET_KIND) :: frame_number = 1
      type (dimlist_type), pointer :: dimlist_head => null()
      type (dimlist_type), pointer :: dimlist_tail => null()
      type (fieldlist_type), pointer :: fieldlist_head => null()
      type (fieldlist_type), pointer :: fieldlist_tail => null()
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      !SHP-MPAS2
      logical :: external_file_desc = .false.
      logical :: restart = .false.
   end type MPAS_IO_Handle_type

   type decomphandle_type
      integer :: field_type
      integer, dimension(:), pointer :: dims
      integer, dimension(:), pointer :: indices
      type (io_desc_t) :: pio_iodesc
   end type decomphandle_type

   type atthandle_type
      character (len=StrKIND) :: attName
      integer :: attType
      integer :: attValueInt
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=StrKIND) :: attValueText
   end type atthandle_type

   type dimhandle_type
      character (len=StrKIND) :: dimname
      logical :: is_unlimited_dim = .false.
      integer :: dimsize
      integer :: dimid
   end type dimhandle_type

   type fieldhandle_type
      character (len=StrKIND) :: fieldname
      integer :: fieldid
      type (Var_desc_t) :: field_desc
      integer :: field_type
      logical :: has_unlimited_dim = .false.
      integer :: ndims
      integer :: precision
      type (dimhandle_type), pointer, dimension(:) :: dims
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      type (decomphandle_type), pointer :: decomp => null()
   end type fieldhandle_type

   type decomplist_type
      type (decomphandle_type), pointer :: decomphandle
      type (decomplist_type), pointer :: next => null()
   end type decomplist_type

   type attlist_type
      type (atthandle_type), pointer :: atthandle
      type (attlist_type), pointer :: next => null()
   end type attlist_type

   type dimlist_type
      type (dimhandle_type), pointer :: dimhandle
      type (dimlist_type), pointer :: next => null()
   end type dimlist_type

   type fieldlist_type
      type (fieldhandle_type), pointer :: fieldhandle
      type (fieldlist_type), pointer :: next => null()
   end type fieldlist_type


   integer, parameter :: MPAS_STREAM_EXACT_TIME              = 100, &
                         MPAS_STREAM_NEAREST                 = 101, &
                         MPAS_STREAM_LATEST_BEFORE           = 102, &
                         MPAS_STREAM_EARLIEST_AFTER          = 103, &
                         MPAS_STREAM_LATEST_STRICTLY_BEFORE  = 104, &
                         MPAS_STREAM_EARLIEST_STRICTLY_AFTER = 105
   
   integer, parameter :: MPAS_STREAM_NOERR              =  0, &
                         MPAS_STREAM_NOT_INITIALIZED    = -1, &
                         MPAS_STREAM_FIELD_NOT_FOUND    = -2, &
                         MPAS_STREAM_CLOBBER_FILE       = -3, &
                         MPAS_STREAM_CLOBBER_RECORD     = -4, &
                         MPAS_IO_ERR                    = -5

   integer, parameter :: FIELD_0D_INT   =  1, &
                         FIELD_1D_INT   =  2, &
                         FIELD_2D_INT   =  3, &
                         FIELD_3D_INT   =  4, &
                         FIELD_0D_REAL  =  5, &
                         FIELD_1D_REAL  =  6, &
                         FIELD_2D_REAL  =  7, &
                         FIELD_3D_REAL  =  8, &
                         FIELD_4D_REAL  =  9, &
                         FIELD_5D_REAL  =  10, &
                         FIELD_0D_CHAR  =  11, &
                         FIELD_1D_CHAR  =  12

   type field_list_type
      integer :: field_type = -999
      logical :: isDecomposed
      integer :: totalDimSize     ! Total size of outer dimension across all blocks for decomposed fields
      logical, dimension(:), pointer :: isAvailable => null() ! Used for reading var-arrays where one or more 
                                                              !   constitutent arrays may not be present in the input file
      type (field0dInteger), pointer :: int0dField => null()
      type (field1dInteger), pointer :: int1dField => null()
      type (field2dInteger), pointer :: int2dField => null()
      type (field3dInteger), pointer :: int3dField => null()
      type (field0dReal), pointer :: real0dField => null()
      type (field1dReal), pointer :: real1dField => null()
      type (field2dReal), pointer :: real2dField => null()
      type (field3dReal), pointer :: real3dField => null()
      type (field4dReal), pointer :: real4dField => null()
      type (field5dReal), pointer :: real5dField => null()
      type (field0dChar), pointer :: char0dField => null()
      type (field1dChar), pointer :: char1dField => null()
      type (field_list_type), pointer :: next => null()
   end type field_list_type

   type MPAS_Stream_type
      logical :: isInitialized = .false.
      integer :: ioFormat
      integer :: ioDirection
      integer :: defaultPrecision = MPAS_IO_NATIVE_PRECISION
      logical :: clobberRecords = .false.
      character(len=StrKIND) :: filename
      type (MPAS_IO_Handle_type) :: fileHandle
      type (att_list_type), pointer :: attList => null()
      type (field_list_type), pointer :: fieldList => null()
      !SHP-MPAS2
      integer :: framesPerFile
   end type MPAS_Stream_type

  !SHP-MPAS2
  type io_output_object
     character (len=StrKIND) :: filename
     integer :: stream

     integer :: time

     type (MPAS_Stream_type) :: io_stream
  end type io_output_object


    integer, parameter :: MPAS_STREAM_LIST_NOERR     = 0, &
                          MPAS_STREAM_LIST_DUPLICATE = 1, &
                          MPAS_STREAM_LIST_NOT_FOUND = 2

    type MPAS_stream_list_type

        ! Used by list head
        integer :: nItems = 0
        type (MPAS_stream_list_type), pointer :: head => null()

        ! Used by streams
        integer :: direction
        logical :: valid = .false.
        logical :: immutable = .false.
        logical :: active_stream = .true.
        character(len=StrKIND) :: filename
        character(len=StrKIND) :: filename_template
        character(len=StrKIND) :: filename_interval
        type (MPAS_Stream_type), pointer :: stream => null()
        integer :: timeLevel = 0
        integer :: nRecords
        integer :: precision = MPAS_IO_NATIVE_PRECISION
        integer :: clobber_mode
        integer :: io_type
        type (MPAS_TimeInterval_type), pointer :: recordInterval => null()
        type (MPAS_stream_list_type), pointer :: alarmList_in => null()
        type (MPAS_stream_list_type), pointer :: alarmList_out => null()
        type (mpas_pool_type), pointer :: att_pool => null()
        type (mpas_pool_type), pointer :: field_pool => null()
        type (mpas_pool_type), pointer :: field_pkg_pool => null()
        type (mpas_pool_type), pointer :: pkg_pool => null()
        type (MPAS_Time_type), pointer :: referenceTime => null()

        ! Used by alarms
        type (MPAS_stream_list_type), pointer :: streamList => null()

        ! Used by streams and alarms
        character(len=StrKIND) :: name
        type (MPAS_stream_list_type), pointer :: xref => null()
        type (MPAS_stream_list_type), pointer :: next => null()

    end type MPAS_stream_list_type


    integer, public, parameter :: MPAS_STREAM_ERR_FATAL  = 1000, &
                                  MPAS_STREAM_ERR_WARN   = 1001, &
                                  MPAS_STREAM_ERR_SILENT = 1002

    integer, public, parameter :: MPAS_STREAM_MGR_NOERR            =  0, &
                                  MPAS_STREAM_MGR_ERR_CLOBBER_FILE = -1, &
                                  MPAS_STREAM_MGR_ERR_CLOBBER_REC  = -2, &
                                  MPAS_STREAM_MGR_ERROR            = -3

    integer, public, parameter :: MPAS_STREAM_INPUT        = 1, &
                                  MPAS_STREAM_OUTPUT       = 2, &
                                  MPAS_STREAM_INPUT_OUTPUT = 3, &
                                  MPAS_STREAM_NONE         = 4

    integer, public, parameter :: MPAS_STREAM_PROPERTY_ACTIVE        =  5, &
                                  MPAS_STREAM_PROPERTY_IMMUTABLE     =  6, &
                                  MPAS_STREAM_PROPERTY_FILENAME      =  7, &
                                  MPAS_STREAM_PROPERTY_REF_TIME      =  8, &
                                  MPAS_STREAM_PROPERTY_RECORD_INTV   =  9, &
                                  MPAS_STREAM_PROPERTY_PRECISION     = 10, &
                                  MPAS_STREAM_PROPERTY_FILENAME_INTV = 11, &
                                  MPAS_STREAM_PROPERTY_CLOBBER       = 12, &
                                  MPAS_STREAM_PROPERTY_IOTYPE        = 13

    integer, public, parameter :: MPAS_STREAM_CLOBBER_NEVER     = 100, &
                                  MPAS_STREAM_CLOBBER_APPEND    = 101, &
                                  MPAS_STREAM_CLOBBER_TRUNCATE  = 102, &
                                  MPAS_STREAM_CLOBBER_OVERWRITE = 103

    type MPAS_streamManager_type

        integer :: numStreams = 0
        integer :: errorLevel

        type (MPAS_Clock_type), pointer :: streamClock
        type (MPAS_Pool_type), pointer :: allFields
        type (MPAS_Pool_type), pointer :: allPackages
        type (MPAS_Pool_type), pointer :: allStructs
        type (MPAS_Pool_type), pointer :: defaultAtts

        type (MPAS_stream_list_type), pointer :: currentStream => NULL()
        type (MPAS_stream_list_type), pointer :: streams
        type (MPAS_stream_list_type), pointer :: alarms_in
        type (MPAS_stream_list_type), pointer :: alarms_out

    end type MPAS_streamManager_type



   integer, parameter :: MPAS_MAX_ALARMS = 40

   integer, parameter :: MPAS_NOW = 0, &
                         MPAS_START_TIME = 1, &
                         MPAS_STOP_TIME = 2

   integer, parameter :: MPAS_FORWARD = 1, &
                         MPAS_BACKWARD = -1

   integer, parameter :: MPAS_GREGORIAN = 0, &
                         MPAS_GREGORIAN_NOLEAP = 1, &
                         MPAS_360DAY = 2

   type MPAS_Time_type
      type (ESMF_Time) :: t
   end type

   type MPAS_TimeInterval_type
      type (ESMF_TimeInterval) :: ti
   end type

   type MPAS_Alarm_type
      character (len=ShortStrKIND) :: alarmID
      logical :: isRecurring
      logical :: isSet
      type (MPAS_Time_type) :: ringTime
      type (MPAS_Time_type) :: prevRingTime
      type (MPAS_TimeInterval_type) :: ringTimeInterval
      type (MPAS_Alarm_type), pointer :: next
   end type
   
   type MPAS_Clock_type
      integer :: direction
      integer :: nAlarms
      type (ESMF_Clock) :: c
      type (MPAS_Alarm_type), pointer :: alarmListHead
   end type


        type timer_node
          character (len=StrKIND) :: timer_name
          logical :: running, printable
          integer :: levels, calls, nlen
          real (kind=R8KIND) :: start_time, end_time, total_time
          real (kind=RKIND) :: max_time, min_time, avg_time
          real (kind=RKIND) :: efficiency
          type (timer_node), pointer :: next
        end type timer_node



   ! Type for storing (possibly architecture specific) information concerning to parallelism
   type parallel_info
      type (mpas_multihalo_exchange_list), pointer :: cellsToSend            ! List of types describing which cells to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToRecv            ! List of types describing which cells to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToCopy            ! List of types describing which cells to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: edgesToSend            ! List of types describing which edges to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToRecv            ! List of types describing which edges to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToCopy            ! List of types describing which edges to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: verticesToSend         ! List of types describing which vertices to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToRecv         ! List of types describing which vertices to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToCopy         ! List of types describing which vertices to copy from other blocks
   end type parallel_info


   ! Derived type for storing part of a domain; used as a basic unit of work for a process
   type block_type

      integer :: blockID   ! Unique global ID number for this block
      integer :: localBlockID  ! Unique local ID number for this block

      type (domain_type), pointer :: domain

      type (parallel_info), pointer :: parinfo

      type (block_type), pointer :: prev => null()
      type (block_type), pointer :: next => null()

      type (mpas_pool_type), pointer :: structs, dimensions, configs, packages
      type (mpas_pool_type), pointer :: allFields, allStructs
   end type block_type

   integer, parameter :: MPAS_DECOMP_NOERR = 1000, &
                         MPAS_DECOMP_ERROR = 1001

   abstract interface
      function mpas_decomp_function(block, manager, globalDimSize, numBlocks, ownedIndices) result(iErr)
         import block_type
         import mpas_streamManager_type

         type (block_type), intent(in) :: block
         type (mpas_streamManager_type), intent(inout) :: manager
         integer, intent(in) :: globalDimSize
         integer, intent(in) :: numBlocks
         integer, dimension(:), pointer :: ownedIndices
         integer :: iErr
      end function
   end interface

   type mpas_decomp_list
      integer :: nameLen
      character (len=StrKIND) :: decompName
      procedure (mpas_decomp_function), pointer, nopass :: decompFunc => null()
      type (mpas_decomp_list), pointer :: next => null()
   end type mpas_decomp_list

   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
      type (mpas_pool_type), pointer :: configs, packages

      type (MPAS_Clock_type), pointer :: clock
      type (MPAS_streamManager_type), pointer :: streamManager
      type (mpas_decomp_list), pointer :: decompositions

      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo

      ! Domain specific constants
      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius
      character (len=StrKIND) :: namelist_filename !< Constant: Name of namelist file
      character (len=StrKIND) :: streams_filename !< Constant: Name of stream configuration file
      character (len=StrKIND) :: mesh_spec !< mesh_spec attribute, read in from input file.
      character (len=StrKIND) :: parent_id !< parent_id attribute, read in from input file.

      ! Back pointer to core
      type (core_type), pointer :: core => null()

      ! Domain_type is a linked list
      type (domain_type), pointer :: next => null()

      !SHP-MPAS2
      character (len=StrKIND*2) :: history !< History attribute, read in from input file. 
   end type domain_type

   abstract interface
      function mpas_setup_namelist_function(configs, namelistFilename, dminfo) result(iErr)
         use mpas_kind_types
         import mpas_pool_type
         import dm_info

         type (mpas_pool_type), intent(inout) :: configs
         character (len=*), intent(in) :: namelistFilename
         type (dm_info), intent(in) :: dminfo
         integer :: iErr
      end function mpas_setup_namelist_function
   end interface

   abstract interface
      function mpas_define_packages_function(packages) result(iErr)
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: packages
         integer :: iErr
      end function mpas_define_packages_function
   end interface

   abstract interface
      function mpas_setup_packages_function(configs, packages) result(iErr)
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: configs
         type (mpas_pool_type), intent(inout) :: packages
         integer :: iErr
      end function mpas_setup_packages_function
   end interface

   abstract interface
      function mpas_setup_decompositions_function(decompList) result(iErr)
         import mpas_decomp_list

         type (mpas_decomp_list), pointer :: decompList
         integer :: iErr
      end function mpas_setup_decompositions_function
   end interface

   abstract interface
      function mpas_get_mesh_stream_function(configs, stream) result(iErr)
         use mpas_kind_types
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: configs
         character (len=StrKIND), intent(out) :: stream
         integer :: iErr
      end function mpas_get_mesh_stream_function
   end interface

   abstract interface
      function mpas_setup_clock_function(clock, configs) result(iErr)
         import mpas_clock_type
         import mpas_pool_type

         type (mpas_clock_type), intent(inout) :: clock
         type (mpas_pool_type), intent(inout) :: configs
         integer :: iErr
      end function mpas_setup_clock_function
   end interface

   abstract interface
      function mpas_setup_immutable_streams_function(manager) result(iErr)
         import mpas_streamManager_type

         type (mpas_streamManager_type), pointer :: manager
         integer :: iErr
      end function mpas_setup_immutable_streams_function
   end interface

   abstract interface
      function mpas_setup_block_function(block) result(iErr)
         import block_type

         type (block_type), pointer :: block
         integer :: iErr
      end function mpas_setup_block_function
   end interface

   abstract interface
      function mpas_setup_setup_derived_dimensions_function(readDimensions, dimensionPool, configPool) result(iErr)
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: readDimensions
         type (mpas_pool_type), intent(inout) :: dimensionPool
         type (mpas_pool_type), intent(inout) :: configPool
         integer :: iErr
      end function mpas_setup_setup_derived_dimensions_function
   end interface


   abstract interface
      function mpas_core_init_function(domain, timeStamp) result(iErr)
         import domain_type
         type (domain_type), intent(inout) :: domain
         character (len=*), intent(out) :: timeStamp
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_core_run_function(domain) result(iErr)
         import domain_type
         type (domain_type), intent(inout) :: domain
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_core_finalize_function(domain) result(iErr)
         import domain_type
         type (domain_type), intent(inout) :: domain
         integer :: iErr
      end function
   end interface

   type core_type
      type (domain_type), pointer :: domainlist => null()

      character (len=StrKIND) :: modelName !< Constant: Name of model
      character (len=StrKIND) :: coreName !< Constant: Name of core
      character (len=StrKIND) :: modelVersion !< Constant: Version number
      character (len=StrKIND) :: executableName !< Constant: Name of executable generated at build time.
      character (len=StrKIND) :: git_version !< Constant: Version string from git-describe.
      character (len=StrKIND*2) :: history !< History attribute, read in from input file.
      character (len=StrKIND) :: Conventions !< Conventions attribute, read in from input file.
      character (len=StrKIND) :: source !< source attribute, read in from input file.

      ! Core init, run, and finalize function pointers
      procedure (mpas_core_init_function), pointer, nopass :: core_init => null()
      procedure (mpas_core_run_function), pointer, nopass :: core_run => null()
      procedure (mpas_core_finalize_function), pointer, nopass :: core_finalize => null()

      ! Core framework function pointers
      procedure (mpas_setup_namelist_function), pointer, nopass :: setup_namelist => null()
      procedure (mpas_define_packages_function), pointer, nopass :: define_packages => null()
      procedure (mpas_setup_packages_function), pointer, nopass :: setup_packages => null()
      procedure (mpas_setup_decompositions_function), pointer, nopass :: setup_decompositions => null()
      procedure (mpas_get_mesh_stream_function), pointer, nopass :: get_mesh_stream => null()
      procedure (mpas_setup_clock_function), pointer, nopass :: setup_clock => null()
      procedure (mpas_setup_block_function), pointer, nopass :: setup_block => null()
      procedure (mpas_setup_immutable_streams_function), pointer, nopass :: setup_immutable_streams => null()
      procedure (mpas_setup_setup_derived_dimensions_function), pointer, nopass :: setup_derived_dimensions => null()

      ! core_type is a linked list
      type (core_type), pointer :: next => null()
   end type core_type

  ! abstract interface for variable interface function
  abstract interface
     function variable_interval(currentTime) result(variableInterval)
       import MPAS_Time_type
       import MPAS_TimeInterval_type
       type(MPAS_Time_type), intent(in) :: currentTime
       type(MPAS_TimeInterval_type) :: variableInterval
     end function variable_interval
  end interface

  ! individual forcing field associated with a forcing stream
  type :: mpas_forcing_field_type

     ! variable identificaion
     character(len=strKIND) :: forcingName

     character(len=strKIND) :: poolname  ! pool name of the input field
     character(len=strKIND) :: fieldname ! field name of the output field

     ! linked list next pointer
     type (mpas_forcing_field_type), pointer :: next => null()

  end type mpas_forcing_field_type

  ! forcing stream associated with a forcing group
  type :: mpas_forcing_stream_type

     ! stream identification
     character(len=strKIND) :: forcingStreamID ! the stream ID for this forcing stream

     ! alarm for reading in the next forcing time
     character(len=ShortStrKind) :: forcingAlarmID ! the alarm ID for this forcing stream

     ! interpolation
     character(len=strKIND) :: interpolationType ! the interpolation type (e.g. 'linear', 'constant')
     integer :: nTimeStencil ! the number of forcing data times
     integer :: nTimeStencilLower ! the number of forcing data times less than the current forcing time
     integer :: nTimeStencilUpper ! the number of forcing data times more than the current forcing time
     type(MPAS_time_type), dimension(:), allocatable :: forcingTimes ! the forcing data times
     
     ! forcing times definition
     type(MPAS_TimeInterval_type) :: forcingIntervalConstant ! the forcing interval
     type(MPAS_Time_type) :: forcingReferenceTime ! the forcing reference time

     ! forcing initialization type
     character(len=strKIND) :: forcingInitializationType

     ! optional functions to calculate variable intervals
     procedure (variable_interval), pointer, nopass :: variable_interval_forward_ptr => null ()
     procedure (variable_interval), pointer, nopass :: variable_interval_backward_ptr => null ()

     ! function pointer to new functions to test same as stream ones
     procedure (variable_interval), pointer, nopass :: variable_interval_forward_test_ptr => null ()
     procedure (variable_interval), pointer, nopass :: variable_interval_backward_test_ptr => null ()

     ! linked list of individual forcing fields
     type (mpas_forcing_field_type), pointer :: field => null()

     ! linked list next pointer
     type (mpas_forcing_stream_type), pointer :: next => null()

  end type mpas_forcing_stream_type

  ! collection of forcing steams and fields using the same forcing clock
  type :: mpas_forcing_group_type

     ! the forcings name
     character(len=strKIND) :: forcingGroupName ! the forcing group identifying name

     ! pointer to the forcing group domain
     type(domain_type), pointer :: domain_ptr

     ! forcing clock and cycling times
     type (MPAS_Clock_type) :: forcingClock ! the forcing clock

     logical                       :: forcingCycleUse            ! whether we cycle the forcing clock
     type (MPAS_Time_Type)         :: forcingCycleStart          ! the start year of the forcing cycle
     type (MPAS_Time_Type)         :: forcingCycleEnd            ! the end year of the forcing cycle
     type (MPAS_TimeInterval_Type) :: forcingCycleDuration       ! the duration of the forcing cycle
     logical                       :: forcingCycleStartInclusive ! cycle start time is inclusive to the cycle

      ! the alarm ID for cycling the clock
     character(len=ShortStrKIND) :: forcingCycleAlarmID = "forcingCycleAlarmID"

     ! linked list of individual streams
     type (mpas_forcing_stream_type), pointer :: stream => null()

     ! linked list next pointer
     type (mpas_forcing_group_type), pointer :: next => null()

  end type mpas_forcing_group_type


   contains

end module mpas_derived_types
