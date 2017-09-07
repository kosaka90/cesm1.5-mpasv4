module mpas_io_restart

   use mpas_attlist
   use mpas_dmpar
   use mpas_io_units, only : stdoutUnit, stderrUnit

   use pio
   use piolib_mod
   use pionfatt_mod

   use mpas_io
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_framework
   use mpas_kind_types


   type (iosystem_desc_t), pointer, private, save :: pio_iosystem
   type (decomplist_type), pointer, private :: decomp_list => null()
   type (dm_info), private :: local_dminfo
   integer, private:: master_pio_iotype = -999

   contains
   
   subroutine mpas_io_output_init_MPAS2( domain, output_obj, &!{{{
                               dminfo, &
                               mesh, &
                               PIOFileDesc &
                             )

      use mpas_io_streams
      use mpas_stream_manager
  
      implicit none

      type (domain_type), intent(in) :: domain
      type (io_output_object), intent(inout) :: output_obj
      type (dm_info), intent(in) :: dminfo
      type (mpas_pool_type), intent(in) :: mesh
      type (file_desc_t), intent(inout), optional :: PIOFileDesc

      type (mpas_pool_field_info_type) :: fieldInfo

      character(len=StrKIND) :: fieldName
      logical :: fieldActive

!      real(kind=RKIND), dimension(:,:,:), pointer :: real3d
!      real(kind=RKIND), dimension(:,:), pointer :: real2d
!      real(kind=RKIND), dimension(:), pointer :: real1d
!      real(kind=RKIND), pointer :: real0d
!
!      integer, dimension(:,:), pointer :: int2d
!      integer, dimension(:), pointer :: int1d
!      integer, pointer :: int0d
!
!      character(len=StrKIND), dimension(:), pointer :: char1d
!      character(len=StrKIND), pointer :: char0d

        type (field3DReal), pointer :: real3d
        type (field2DReal), pointer :: real2d
        type (field1DReal), pointer :: real1d
        type (field0DReal), pointer :: real0d

        type (field3DInteger), pointer :: int3d
        type (field2DInteger), pointer :: int2d
        type (field1DInteger), pointer :: int1d
        type (field0DInteger), pointer :: int0d

        type (field1DChar), pointer :: char1d
        type (field0DChar), pointer :: char0d


      integer :: nferr, ierr
      integer, dimension(10) :: dimlist
      character (len=StrKIND*4) :: runCmd

      if(len(trim(domain % history)) > 0) then
          write(runCmd,'(a,a,i0,a,a,a)') trim(domain % history),' mpirun -n ',domain % dminfo % nProcs, ' ', trim(domain % core %coreName), '_model; '
      else
          write(runCmd,'(a,i0,a,a,a)') 'mpirun -n ',domain % dminfo % nProcs, ' ', trim(domain % core % coreName), '_model; '
      end if

      !if (trim(config_pio_format) == 'pnetcdf') then
      !   call MPAS_createStream_MPAS2(output_obj % io_stream, trim(output_obj % filename), MPAS_IO_PNETCDF, MPAS_IO_WRITE, 1, PIOFileDesc, ierr=nferr)
      !else
         call MPAS_createStream_MPAS2(output_obj % io_stream, trim(output_obj % filename), MPAS_IO_NETCDF, MPAS_IO_WRITE, 1, PIOFileDesc, ierr=nferr)
      !end if

      call MPAS_stream_mgr_begin_iteration(domain % streamManager, streamID='restart')

      do while ( MPAS_stream_mgr_get_next_field(domain % streamManager, 'restart', fieldName, fieldActive) )
          if (fieldActive) then
              call mpas_pool_get_field_info(domain % blocklist % allFields, trim(fieldName), fieldInfo)
              if (fieldInfo % fieldType == MPAS_POOL_REAL) then
                  if (fieldInfo % nDims == 3) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), real3d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, real3d)
                  else if (fieldInfo % nDims == 2) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), real2d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, real2d)
                  else if (fieldInfo % nDims == 1) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), real1d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, real1d)
                  else if (fieldInfo % nDims == 0) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), real0d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, real0d)
                  end if
              else if (fieldInfo % fieldType == MPAS_POOL_INTEGER) then
                  if (fieldInfo % nDims == 3) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), int3d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, int3d)
                  else if (fieldInfo % nDims == 2) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), int2d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, int2d)
                  else if (fieldInfo % nDims == 1) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), int1d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, int1d)
                  else if (fieldInfo % nDims == 0) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), int0d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, int0d)
                  end if
              else if (fieldInfo % fieldType == MPAS_POOL_CHARACTER) then
                  if (fieldInfo % nDims == 0) then
                      call mpas_pool_get_field(domain % blocklist % allFields, trim(fieldName), char0d, 1)
                      call MPAS_streamAddField(output_obj % io_stream, char0d)
                  end if
              end if

          end if
      end do

      if (domain % on_a_sphere) then
         call MPAS_writeStreamAtt(output_obj % io_stream, 'on_a_sphere', 'YES             ', nferr)
      else
         call MPAS_writeStreamAtt(output_obj % io_stream, 'on_a_sphere', 'NO              ', nferr)
      end if

      call MPAS_writeStreamAtt(output_obj % io_stream, 'sphere_radius', domain % sphere_radius, nferr)

   end subroutine mpas_io_output_init_MPAS2!}}}

   subroutine MPAS_createStream_MPAS2(stream, fileName, ioFormat, ioDirection, framesPerFile, PIOFileDesc, ierr)

      use pio_types

      implicit none

      type (MPAS_Stream_type), intent(out) :: stream
      character (len=*), intent(in) :: fileName
      integer, intent(in) :: ioFormat
      integer, intent(in) :: ioDirection
      integer, intent(in) :: framesPerFile
      type (file_desc_t), intent(in), optional :: PIOFileDesc
      integer, intent(out), optional :: ierr

      integer :: io_err

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      stream % fileHandle = MPAS_io_open_MPAS2(fileName, ioDirection, ioFormat, file_desc=PIOFileDesc, ierr=io_err)
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

      stream % ioDirection = ioDirection
      stream % ioFormat = ioFormat
      !stream % framesPerFile = framesPerFile

      stream % isInitialized = .true.

   end subroutine MPAS_createStream_MPAS2

   type (MPAS_IO_Handle_type) function MPAS_io_open_MPAS2(filename, mode, ioformat, file_desc, ierr)

      implicit none

      character (len=*), intent(in) :: filename
      integer, intent(in) :: mode
      integer, intent(in) :: ioformat
      type (file_desc_t), intent(in), optional :: file_desc
      integer, intent(out), optional :: ierr

      integer :: pio_ierr, pio_iotype

!      write(stderrUnit,*) 'Called MPAS_io_open()'
      if (present(ierr)) ierr = MPAS_IO_NOERR


      ! Sanity checks
      if (mode /= MPAS_IO_READ .and. &
          mode /= MPAS_IO_WRITE) then
         if (present(ierr)) ierr = MPAS_IO_ERR_INVALID_MODE
         return
      end if
      if (ioformat /= MPAS_IO_NETCDF .and. &
          ioformat /= MPAS_IO_PNETCDF) then
         if (present(ierr)) ierr = MPAS_IO_ERR_INVALID_FORMAT
         return
      end if
      if (len(filename) > 1024) then
         if (present(ierr)) ierr = MPAS_IO_ERR_LONG_FILENAME
         return
      end if

      MPAS_io_open_MPAS2 % filename = filename
      MPAS_io_open_MPAS2 % iomode   = mode
      MPAS_io_open_MPAS2 % ioformat = ioformat

      if (.not. present(file_desc)) then
         if (master_pio_iotype /= -999) then
            pio_iotype = master_pio_iotype
         else
            if (ioformat == MPAS_IO_PNETCDF) then
               pio_iotype = PIO_iotype_pnetcdf
            else
               pio_iotype = PIO_iotype_netcdf
            end if
         end if

         if (mode == MPAS_IO_WRITE) then
!write(stderrUnit,*) 'MGD PIO_createfile'
            pio_ierr = PIO_createfile(pio_iosystem, MPAS_io_open_MPAS2 % pio_file, pio_iotype, trim(filename), PIO_64BIT_OFFSET)
         else
!write(stderrUnit,*) 'MGD PIO_openfile'
            pio_ierr = PIO_openfile(pio_iosystem, MPAS_io_open_MPAS2 % pio_file, pio_iotype, trim(filename), PIO_nowrite)
         endif
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            return
         end if

         MPAS_io_open_MPAS2 % external_file_desc = .false.
      else
         MPAS_io_open_MPAS2 % pio_file = file_desc
         MPAS_io_open_MPAS2 % external_file_desc = .true.
         !MPAS_io_open_MPAS2 % preexisting_file = .true.
         MPAS_io_open_MPAS2 % restart = .true.
      end if

      if (mode == MPAS_IO_READ) then
!MPAS_io_open % pio_unlimited_dimid = 44
         pio_ierr = PIO_inquire(MPAS_io_open_MPAS2 % pio_file, unlimitedDimID=MPAS_io_open_MPAS2 % pio_unlimited_dimid)
!write(stderrUnit,*) 'Found unlimited dim ', MPAS_io_open_MPAS2 % pio_unlimited_dimid
         if (pio_ierr /= PIO_noerr) then
            if (present(ierr)) ierr = MPAS_IO_ERR_PIO
            return
         end if
      end if

      MPAS_io_open_MPAS2 % initialized = .true.

      return

   end function MPAS_io_open_MPAS2


   subroutine mpas_output_state_for_domain_MPAS2(output_obj, domain, itime)!{{{
   
      implicit none
   
      type (io_output_object), intent(inout) :: output_obj
      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itime

      type (mpas_pool_type), pointer :: mesh
      type(block_type), pointer :: block_ptr

      integer, pointer :: nCells, nEdges, nVertices, vertexDegree
      integer, pointer :: maxEdges, maxEdges2, nEdgesSolve, nCellsSolve, nVerticesSolve
      integer :: ierr
      integer :: i, j

      type (field1dInteger), pointer :: nEdgesOnCell, nEdgesOnEdge, indexToCellID, indexToEdgeID, indexToVertexID

      type (field2dInteger), pointer :: cellsOnCell_save, edgesOnCell_save, verticesOnCell_save, &
                               cellsOnEdge_save, verticesOnEdge_save, edgesOnEdge_save, &
                               cellsOnVertex_save, edgesOnVertex_save

      type (field2dInteger), pointer :: cellsOnCell_ptr, edgesOnCell_ptr, verticesOnCell_ptr, &
                               cellsOnEdge_ptr, verticesOnEdge_ptr, edgesOnEdge_ptr, &
                               cellsOnVertex_ptr, edgesOnVertex_ptr

      type (field2dInteger), pointer :: cellsOnCell, edgesOnCell, verticesOnCell, &
                               cellsOnEdge, verticesOnEdge, edgesOnEdge, &
                               cellsOnVertex, edgesOnVertex

      output_obj % time = itime

      !
      ! Convert connectivity information from local to global indices
      ! Needs to be done block by block
      !
      ! Also, backup local indices to be copied back into blocks after output is complete.
      !
      allocate(cellsOnCell_save)
      allocate(edgesOnCell_save) 
      allocate(verticesOnCell_save)
      allocate(cellsOnEdge_save)
      allocate(verticesOnEdge_save)
      allocate(edgesOnEdge_save)
      allocate(cellsOnVertex_save)
      allocate(edgesOnVertex_save)

      cellsOnCell_ptr => cellsOnCell_save
      edgesOnCell_ptr => edgesOnCell_save 
      verticesOnCell_ptr => verticesOnCell_save
      cellsOnEdge_ptr => cellsOnEdge_save 
      verticesOnEdge_ptr => verticesOnEdge_save 
      edgesOnEdge_ptr => edgesOnEdge_save
      cellsOnVertex_ptr => cellsOnVertex_save 
      edgesOnVertex_ptr => edgesOnVertex_save

      block_ptr => domain % blocklist
      do while(associated(block_ptr))

        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)

        call mpas_pool_get_dimension(mesh,'maxEdges',maxEdges)
        call mpas_pool_get_dimension(mesh,'maxEdges2',maxEdges2)
        call mpas_pool_get_dimension(mesh,'vertexDegree',vertexDegree)
        call mpas_pool_get_dimension(mesh,'nCells',nCells)
        call mpas_pool_get_dimension(mesh,'nEdges',nEdges)
        call mpas_pool_get_dimension(mesh,'nVertices',nVertices)
        call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
        call mpas_pool_get_dimension(mesh,'nEdgesSolve',nEdgesSolve)
        call mpas_pool_get_dimension(mesh,'nVerticesSolve',nVerticesSolve)

        call mpas_pool_get_field(mesh, 'cellsOnCell', cellsOnCell)
        !nullify(cellsOncell_ptr % ioinfo)
        cellsOnCell_ptr % array => cellsOnCell % array
        allocate(cellsOnCell % array(maxEdges, nCells+1))

        call mpas_pool_get_field(mesh, 'edgesOnCell', edgesOnCell)
        !nullify(edgesOnCell_ptr % ioinfo)
        edgesOnCell_ptr % array => edgesOnCell % array
        allocate(edgesOnCell % array(maxEdges, nCells+1))

        call mpas_pool_get_field(mesh, 'verticesOnCell', verticesOnCell)
        !nullify(verticesOnCell_ptr % ioinfo)
        verticesOnCell_ptr % array => verticesOnCell % array
        allocate(verticesOnCell % array(maxEdges, nCells+1))

        call mpas_pool_get_field(mesh, 'cellsOnEdge', cellsOnEdge)
        !nullify(cellsOnEdge_ptr % ioinfo)
        cellsOnEdge_ptr % array => cellsOnEdge % array
        allocate(cellsOnEdge % array(2, nEdges+1))

        call mpas_pool_get_field(mesh, 'verticesOnEdge', verticesOnEdge)
        !nullify(verticesOnEdge_ptr % ioinfo)
        verticesOnEdge_ptr % array => verticesOnEdge % array
        allocate(verticesOnEdge % array(2, nEdges+1))

        call mpas_pool_get_field(mesh, 'edgesOnEdge', edgesOnEdge)
        !nullify(edgesOnEdge_ptr % ioinfo)
        edgesOnEdge_ptr % array => edgesOnEdge % array
        allocate(edgesOnEdge % array(maxEdges2, nEdges+1))

        call mpas_pool_get_field(mesh, 'cellsOnVertex', cellsOnVertex)
        !nullify(cellsOnVertex_ptr % ioinfo)
        cellsOnVertex_ptr % array => cellsOnVertex % array
        allocate(cellsOnVertex % array(vertexDegree, nVertices+1))

        call mpas_pool_get_field(mesh, 'edgesOnVertex', edgesOnVertex)
        !nullify(edgesOnVertex_ptr % ioinfo)
        edgesOnVertex_ptr % array => edgesOnVertex % array
        allocate(edgesOnVertex % array(vertexDegree, nVertices+1))

        call mpas_pool_get_field(mesh, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_field(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
        call mpas_pool_get_field(mesh, 'indexToCellID', indexToCellID)
        call mpas_pool_get_field(mesh, 'indexToEdgeID', indexToEdgeID)
        call mpas_pool_get_field(mesh, 'indexToVertexID', indexToVertexID)

        do i = 1, nCellsSolve
          do j = 1, nEdgesOnCell % array(i)
            cellsOnCell % array(j, i) = indexToCellID % array(cellsOnCell_ptr % array(j, i))
            edgesOnCell % array(j, i) = indexToEdgeID % array(edgesOnCell_ptr % array(j, i))
            verticesOnCell % array(j, i) = indexToVertexID % array(verticesOnCell_ptr % array(j, i))
          end do

          cellsOnCell % array(nEdgesOnCell % array(i) + 1:maxEdges, i) = nCells+1
          edgesOnCell % array(nEdgesOnCell % array(i) + 1:maxEdges, i) = nEdges+1
          verticesOnCell % array(nEdgesOnCell % array(i) + 1:maxEdges, i) = nVertices+1
        end do

        do i = 1, nEdgesSolve
          cellsOnEdge % array(1, i) = indexToCellID % array(cellsOnEdge_ptr % array(1, i))
          cellsOnEdge % array(2, i) = indexToCellID % array(cellsOnEdge_ptr % array(2, i))

          verticesOnedge % array(1, i) = indexToVertexID % array(verticesOnEdge_ptr % array(1,i))
          verticesOnedge % array(2, i) = indexToVertexID % array(verticesOnEdge_ptr % array(2,i))

          do j = 1, nEdgesOnEdge % array(i)
            edgesOnEdge % array(j, i) = indexToEdgeID % array(edgesOnEdge_ptr % array(j, i))
          end do

          edgesOnEdge % array(nEdgesOnEdge % array(i)+1:maxEdges2, i) = nEdges+1
        end do

        do i = 1, nVerticesSolve
          do j = 1, vertexDegree
            cellsOnVertex % array(j, i) = indexToCellID % array(cellsOnVertex_ptr % array(j, i))
            edgesOnVertex % array(j, i) = indexToEdgeID % array(edgesOnVertex_ptr % array(j, i))
          end do
        end do

        block_ptr => block_ptr % next
        if(associated(block_ptr)) then
          allocate(cellsOnCell_ptr % next)
          allocate(edgesOnCell_ptr % next)
          allocate(verticesOnCell_ptr % next)
          allocate(cellsOnEdge_ptr % next)
          allocate(verticesOnEdge_ptr % next)
          allocate(edgesOnEdge_ptr % next)
          allocate(cellsOnVertex_ptr % next)
          allocate(edgesOnVertex_ptr % next)

          cellsOnCell_ptr => cellsOnCell_ptr % next
          edgesOnCell_ptr => edgesOnCell_ptr % next
          verticesOnCell_ptr => verticesOnCell_ptr % next
          cellsOnEdge_ptr => cellsOnEdge_ptr % next
          verticesOnEdge_ptr => verticesOnEdge_ptr % next
          edgesOnEdge_ptr => edgesOnEdge_ptr % next
          cellsOnVertex_ptr => cellsOnVertex_ptr % next
          edgesOnVertex_ptr => edgesOnVertex_ptr % next
        end if

        nullify(cellsOnCell_ptr % next)
        nullify(edgesOnCell_ptr % next)
        nullify(verticesOnCell_ptr % next)
        nullify(cellsOnEdge_ptr % next)
        nullify(verticesOnEdge_ptr % next)
        nullify(edgesOnEdge_ptr % next)
        nullify(cellsOnVertex_ptr % next)
        nullify(edgesOnVertex_ptr % next)
      end do

      ! Write output file
      call MPAS_writeStream_MPAS2(output_obj % io_stream, output_obj % time, ierr)

      ! Converge indices back to local indices, and deallocate all temporary arrays.
      #if 0
      cellsOnCell_ptr => cellsOnCell_save
      edgesOnCell_ptr => edgesOnCell_save 
      verticesOnCell_ptr => verticesOnCell_save
      cellsOnEdge_ptr => cellsOnEdge_save 
      verticesOnEdge_ptr => verticesOnEdge_save 
      edgesOnEdge_ptr => edgesOnEdge_save
      cellsOnVertex_ptr => cellsOnVertex_save 
      edgesOnVertex_ptr => edgesOnVertex_save

      block_ptr => domain % blocklist
      do while(associated(block_ptr))

        deallocate(cellsOnCell % array)
        deallocate(edgesOnCell % array)
        deallocate(verticesOnCell % array)
        deallocate(cellsOnEdge % array)
        deallocate(verticesOnEdge % array)
        deallocate(edgesOnEdge % array)
        deallocate(cellsOnVertex % array)
        deallocate(edgesOnVertex % array)

        cellsOncell % array => cellsOnCell_ptr % array
        edgesOnCell % array => edgesOnCell_ptr % array
        verticesOnCell % array => verticesOnCell_ptr % array
        cellsOnEdge % array => cellsOnEdge_ptr % array
        verticesOnEdge % array => verticesOnEdge_ptr % array
        edgesOnEdge % array => edgesOnEdge_ptr % array
        cellsOnVertex % array => cellsOnVertex_ptr % array
        edgesOnVertex % array => edgesOnVertex_ptr % array

        nullify(cellsOnCell_ptr % array)
        nullify(edgesOnCell_ptr % array)
        nullify(verticesOnCell_ptr % array)
        nullify(cellsOnEdge_ptr % array)
        nullify(verticesOnEdge_ptr % array)
        nullify(edgesOnEdge_ptr % array)
        nullify(cellsOnVertex_ptr % array)
        nullify(edgesOnVertex_ptr % array)

        block_ptr => block_ptr % next
        cellsOnCell_ptr => cellsOnCell_ptr % next
        edgesOnCell_ptr => edgesOnCell_ptr % next
        verticesOnCell_ptr => verticesOnCell_ptr % next
        cellsOnEdge_ptr => cellsOnEdge_ptr % next
        verticesOnEdge_ptr => verticesOnEdge_ptr % next
        edgesOnEdge_ptr => edgesOnEdge_ptr % next
        cellsOnVertex_ptr => cellsOnVertex_ptr % next
        edgesOnVertex_ptr => edgesOnVertex_ptr % next
      end do

      call mpas_deallocate_field(cellsOnCell_save)
      call mpas_deallocate_field(edgesOnCell_save) 
      call mpas_deallocate_field(verticesOnCell_save)
      call mpas_deallocate_field(cellsOnEdge_save)
      call mpas_deallocate_field(verticesOnEdge_save)
      call mpas_deallocate_field(edgesOnEdge_save)
      call mpas_deallocate_field(cellsOnVertex_save)
      call mpas_deallocate_field(edgesOnVertex_save)
      #endif

   end subroutine mpas_output_state_for_domain_MPAS2!}}}


   subroutine MPAS_writeStream_MPAS2(stream, frame, ierr)

      implicit none

      type (MPAS_Stream_type), intent(inout) :: stream
      integer, intent(in) :: frame
      integer, intent(out), optional :: ierr

      integer :: io_err
      integer :: i, j
      integer :: ncons
      integer, pointer :: ownedSize
      type (field0dInteger), pointer :: field_0dint_ptr
      type (field1dInteger), pointer :: field_1dint_ptr
      type (field2dInteger), pointer :: field_2dint_ptr
      type (field3dInteger), pointer :: field_3dint_ptr
      type (field0dReal), pointer :: field_0dreal_ptr
      type (field1dReal), pointer :: field_1dreal_ptr
      type (field2dReal), pointer :: field_2dreal_ptr
      type (field3dReal), pointer :: field_3dreal_ptr
      type (field4dReal), pointer :: field_4dreal_ptr
      type (field5dReal), pointer :: field_5dreal_ptr
      type (field0dChar), pointer :: field_0dchar_ptr
      type (field1dChar), pointer :: field_1dchar_ptr
      type (field_list_type), pointer :: field_cursor
      integer                            :: int0d_temp
      integer, dimension(:),     pointer :: int1d_temp
      integer, dimension(:,:),   pointer :: int2d_temp
      integer, dimension(:,:,:), pointer :: int3d_temp
      real (kind=RKIND)                                :: real0d_temp
      real (kind=RKIND), dimension(:),         pointer :: real1d_temp
      real (kind=RKIND), dimension(:,:),       pointer :: real2d_temp
      real (kind=RKIND), dimension(:,:,:),     pointer :: real3d_temp
      real (kind=RKIND), dimension(:,:,:,:),   pointer :: real4d_temp
      real (kind=RKIND), dimension(:,:,:,:,:), pointer :: real5d_temp

      if (present(ierr)) ierr = MPAS_STREAM_NOERR

      !
      ! Sanity checks
      !
      if (.not. stream % isInitialized) then
         if (present(ierr)) ierr = MPAS_STREAM_NOT_INITIALIZED
         return
      end if

      !
      ! Set time frame to write
      !
      call MPAS_io_set_frame(stream % fileHandle, frame, io_err) 
      call MPAS_io_err_mesg(io_err, .false.)
      if (io_err /= MPAS_IO_NOERR) then
         if (present(ierr)) ierr = MPAS_IO_ERR
         return
      end if

!      !
!      ! Check whether we will clobber any records
!      !
!      if (MPAS_io_would_clobber_records(stream % fileHandle, io_err)) then
!         if (.not. stream % clobberRecords) then
!            if (present(ierr)) ierr = MPAS_STREAM_CLOBBER_RECORD
!            return
!         else
!            write(stderrUnit,'(a,i4,a)') 'MPAS I/O: Overwriting existing record ', frame, &
!                                         ' in output file '//trim(stream % filename)
!         end if
!      end if

      !
      ! Loop over fields in the stream
      !
      field_cursor => stream % fieldList
      do while (associated(field_cursor))

         if (field_cursor % field_type == FIELD_0D_INT) then

!write(stderrUnit,*) 'Writing out field '//trim(field_cursor % int0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'Copying field from first block'
            int0d_temp = field_cursor % int0dField % scalar

!write(stderrUnit,*) 'MGD calling MPAS_io_put_var now...'
            call MPAS_io_put_var(stream % fileHandle, field_cursor % int0dField % fieldName, int0d_temp, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_INT) then

            if (field_cursor % int1dField % isVarArray) then
               ncons = size(field_cursor % int1dField % constituentNames)
            else
               ncons = 1
               allocate(int1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dint_ptr => field_cursor % int1dField
                  i = 1
                  do while (associated(field_1dint_ptr))
                     if (trim(field_1dint_ptr % dimNames(1)) == 'nCells') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dint_ptr % dimNames(1)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_1dint_ptr % block % dimensions, field_1dint_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % int1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        int0d_temp = field_1dint_ptr % array(j)
                     else
                        int1d_temp(i:i+ownedSize-1) = field_1dint_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_1dint_ptr => field_1dint_ptr % next
                  end do
               else
                  if (field_cursor % int1dField % isVarArray) then
                     int0d_temp = field_cursor % int1dField % array(j)
                  else
                     int1d_temp(:) = field_cursor % int1dField % array(:)
                  end if
               end if

               if (field_cursor % int1dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int1dField % constituentNames(j), int0d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int1dField % fieldName, int1d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (.not. field_cursor % int1dField % isVarArray) then
               deallocate(int1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_INT) then

            if (field_cursor % int2dField % isVarArray) then
               ncons = size(field_cursor % int2dField % constituentNames)
               allocate(int1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int2d_temp(field_cursor % int2dField % dimSizes(1), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dint_ptr => field_cursor % int2dField
                  i = 1
                  do while (associated(field_2dint_ptr))
                     if (trim(field_2dint_ptr % dimNames(2)) == 'nCells') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dint_ptr % dimNames(2)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_2dint_ptr % block % dimensions, field_2dint_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % int2dField % isVarArray) then
                        int1d_temp(i:i+ownedSize-1) = field_2dint_ptr % array(j,1:ownedSize)
                     else
                        int2d_temp(:,i:i+ownedSize-1) = field_2dint_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_2dint_ptr => field_2dint_ptr % next
                  end do
               else
                  if (field_cursor % int2dField % isVarArray) then
                     int1d_temp(:) = field_cursor % int2dField % array(j,:)
                  else
                     int2d_temp(:,:) = field_cursor % int2dField % array(:,:)
                  end if
               end if

               if (field_cursor % int2dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int2dField % constituentNames(j), int1d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int2dField % fieldName, int2d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % int2dField % isVarArray) then
               deallocate(int1d_temp)
            else
               deallocate(int2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_INT) then

            if (field_cursor % int3dField % isVarArray) then
               ncons = size(field_cursor % int3dField % constituentNames)
               allocate(int2d_temp(field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(int3d_temp(field_cursor % int3dField % dimSizes(1), &
                                   field_cursor % int3dField % dimSizes(2), &
                                   field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dint_ptr => field_cursor % int3dField
                  i = 1
                  do while (associated(field_3dint_ptr))
                     if (trim(field_3dint_ptr % dimNames(3)) == 'nCells') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dint_ptr % dimNames(3)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_3dint_ptr % block % dimensions, field_3dint_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % int3dField % isVarArray) then
                        int2d_temp(:,i:i+ownedSize-1) = field_3dint_ptr % array(j,:,1:ownedSize)
                     else
                        int3d_temp(:,:,i:i+ownedSize-1) = field_3dint_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_3dint_ptr => field_3dint_ptr % next
                  end do
               else
                  if (field_cursor % int3dField % isVarArray) then
                     int2d_temp(:,:) = field_cursor % int3dField % array(j,:,:)
                  else
                     int3d_temp(:,:,:) = field_cursor % int3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % int3dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int3dField % constituentNames(j), int2d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % int3dField % fieldName, int3d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % int3dField % isVarArray) then
               deallocate(int2d_temp)
            else
               deallocate(int3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_REAL) then

!write(stderrUnit,*) 'Writing out field '//trim(field_cursor % real0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'Copying field from first block'
            real0d_temp = field_cursor % real0dField % scalar

!write(stderrUnit,*) 'MGD calling MPAS_io_put_var now...'
            call MPAS_io_put_var(stream % fileHandle, field_cursor % real0dField % fieldName, real0d_temp, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_REAL) then

            if (field_cursor % real1dField % isVarArray) then
               ncons = size(field_cursor % real1dField % constituentNames)
            else
               ncons = 1
               allocate(real1d_temp(field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_1dreal_ptr => field_cursor % real1dField
                  i = 1
                  do while (associated(field_1dreal_ptr))
                     if (trim(field_1dreal_ptr % dimNames(1)) == 'nCells') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_1dreal_ptr % dimNames(1)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_1dreal_ptr % block % dimensions, field_1dreal_ptr % dimNames(1), ownedSize)
                     end if

                     if (field_cursor % real1dField % isVarArray) then
! I suspect we will never hit this code, as it doesn't make sense, really
                        real0d_temp = field_1dreal_ptr % array(j)
                     else
                        real1d_temp(i:i+ownedSize-1) = field_1dreal_ptr % array(1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_1dreal_ptr => field_1dreal_ptr % next
                  end do
               else
                  if (field_cursor % real1dField % isVarArray) then
                     real0d_temp = field_cursor % real1dField % array(j)
                  else
                     real1d_temp(:) = field_cursor % real1dField % array(:)
                  end if
               end if

               if (field_cursor % real1dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real1dField % constituentNames(j), real0d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real1dField % fieldName, real1d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (.not. field_cursor % real1dField % isVarArray) then
               deallocate(real1d_temp)
            end if

         else if (field_cursor % field_type == FIELD_2D_REAL) then

            if (field_cursor % real2dField % isVarArray) then
               ncons = size(field_cursor % real2dField % constituentNames)
               allocate(real1d_temp(field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real2d_temp(field_cursor % real2dField % dimSizes(1), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_2dreal_ptr => field_cursor % real2dField
                  i = 1
                  do while (associated(field_2dreal_ptr))
                     if (trim(field_2dreal_ptr % dimNames(2)) == 'nCells') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_2dreal_ptr % dimNames(2)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_2dreal_ptr % block % dimensions, field_2dreal_ptr % dimNames(2), ownedSize)
                     end if

                     if (field_cursor % real2dField % isVarArray) then
                        real1d_temp(i:i+ownedSize-1) = field_2dreal_ptr % array(j,1:ownedSize)
                     else
                        real2d_temp(:,i:i+ownedSize-1) = field_2dreal_ptr % array(:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_2dreal_ptr => field_2dreal_ptr % next
                  end do
               else
                  if (field_cursor % real2dField % isVarArray) then
                     real1d_temp(:) = field_cursor % real2dField % array(j,:)
                  else
                     real2d_temp(:,:) = field_cursor % real2dField % array(:,:)
                  end if
               end if

               if (field_cursor % real2dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real2dField % constituentNames(j), real1d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real2dField % fieldName, real2d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real2dField % isVarArray) then
               deallocate(real1d_temp)
            else
               deallocate(real2d_temp)
            end if

         else if (field_cursor % field_type == FIELD_3D_REAL) then

            if (field_cursor % real3dField % isVarArray) then
               ncons = size(field_cursor % real3dField % constituentNames)
               allocate(real2d_temp(field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real3d_temp(field_cursor % real3dField % dimSizes(1), &
                                    field_cursor % real3dField % dimSizes(2), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_3dreal_ptr => field_cursor % real3dField
                  i = 1
                  do while (associated(field_3dreal_ptr))
                     if (trim(field_3dreal_ptr % dimNames(3)) == 'nCells') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_3dreal_ptr % dimNames(3)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_3dreal_ptr % block % dimensions, field_3dreal_ptr % dimNames(3), ownedSize)
                     end if

                     if (field_cursor % real3dField % isVarArray) then
                        real2d_temp(:,i:i+ownedSize-1) = field_3dreal_ptr % array(j,:,1:ownedSize)
                     else
                        real3d_temp(:,:,i:i+ownedSize-1) = field_3dreal_ptr % array(:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_3dreal_ptr => field_3dreal_ptr % next
                  end do
               else
                  if (field_cursor % real3dField % isVarArray) then
                     real2d_temp(:,:) = field_cursor % real3dField % array(j,:,:)
                  else
                     real3d_temp(:,:,:) = field_cursor % real3dField % array(:,:,:)
                  end if
               end if

               if (field_cursor % real3dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real3dField % constituentNames(j), real2d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real3dField % fieldName, real3d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real3dField % isVarArray) then
               deallocate(real2d_temp)
            else
               deallocate(real3d_temp)
            end if

         else if (field_cursor % field_type == FIELD_4D_REAL) then

            if (field_cursor % real4dField % isVarArray) then
               ncons = size(field_cursor % real4dField % constituentNames)
               allocate(real3d_temp(field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real4d_temp(field_cursor % real4dField % dimSizes(1), &
                                    field_cursor % real4dField % dimSizes(2), &
                                    field_cursor % real4dField % dimSizes(3), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_4dreal_ptr => field_cursor % real4dField
                  i = 1
                  do while (associated(field_4dreal_ptr))
                     if (trim(field_4dreal_ptr % dimNames(4)) == 'nCells') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_4dreal_ptr % dimNames(4)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_4dreal_ptr % block % dimensions, field_4dreal_ptr % dimNames(4), ownedSize)
                     end if

                     if (field_cursor % real4dField % isVarArray) then
                        real3d_temp(:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(j,:,:,1:ownedSize)
                     else
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_4dreal_ptr % array(:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_4dreal_ptr => field_4dreal_ptr % next
                  end do
               else
                  if (field_cursor % real4dField % isVarArray) then
                     real3d_temp(:,:,:) = field_cursor % real4dField % array(j,:,:,:)
                  else
                     real4d_temp(:,:,:,:) = field_cursor % real4dField % array(:,:,:,:)
                  end if
               end if

               if (field_cursor % real4dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real4dField % constituentNames(j), real3d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real4dField % fieldName, real4d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real4dField % isVarArray) then
               deallocate(real3d_temp)
            else
               deallocate(real4d_temp)
            end if

         else if (field_cursor % field_type == FIELD_5D_REAL) then

            if (field_cursor % real5dField % isVarArray) then
               ncons = size(field_cursor % real5dField % constituentNames)
               allocate(real4d_temp(field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            else
               ncons = 1
               allocate(real5d_temp(field_cursor % real5dField % dimSizes(1), &
                                    field_cursor % real5dField % dimSizes(2), &
                                    field_cursor % real5dField % dimSizes(3), &
                                    field_cursor % real5dField % dimSizes(4), &
                                    field_cursor % totalDimSize))
            end if

            do j=1,ncons
               if (field_cursor % isDecomposed) then
                  ! Gather field from across multiple blocks
                  field_5dreal_ptr => field_cursor % real5dField
                  i = 1
                  do while (associated(field_5dreal_ptr))
                     if (trim(field_5dreal_ptr % dimNames(5)) == 'nCells') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nCellsSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nEdges') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nEdgesSolve', ownedSize)
                     else if (trim(field_5dreal_ptr % dimNames(5)) == 'nVertices') then
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, 'nVerticesSolve', ownedSize)
                     else
                        call mpas_pool_get_dimension(field_5dreal_ptr % block % dimensions, field_5dreal_ptr % dimNames(5), ownedSize)
                     end if

                     if (field_cursor % real5dField % isVarArray) then
                        real4d_temp(:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(j,:,:,:,1:ownedSize)
                     else
                        real5d_temp(:,:,:,:,i:i+ownedSize-1) = field_5dreal_ptr % array(:,:,:,:,1:ownedSize)
                     end if
                     i = i + ownedSize
                     field_5dreal_ptr => field_5dreal_ptr % next
                  end do
               else
                  if (field_cursor % real5dField % isVarArray) then
                     real4d_temp(:,:,:,:) = field_cursor % real5dField % array(j,:,:,:,:)
                  else
                     real5d_temp(:,:,:,:,:) = field_cursor % real5dField % array(:,:,:,:,:)
                  end if
               end if

               if (field_cursor % real5dField % isVarArray) then
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real5dField % constituentNames(j), real4d_temp, io_err)
               else
                  call MPAS_io_put_var(stream % fileHandle, field_cursor % real5dField % fieldName, real5d_temp, io_err)
               end if
               call MPAS_io_err_mesg(io_err, .false.)
               if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR
            end do

            if (field_cursor % real5dField % isVarArray) then
               deallocate(real4d_temp)
            else
               deallocate(real5d_temp)
            end if

         else if (field_cursor % field_type == FIELD_0D_CHAR) then

!write(stderrUnit,*) 'Writing out field '//trim(field_cursor % char0dField % fieldName)
!write(stderrUnit,*) '   > is the field decomposed? ', field_cursor % isDecomposed
!write(stderrUnit,*) '   > outer dimension size ', field_cursor % totalDimSize

!write(stderrUnit,*) 'Copying field from first block'
!write(stderrUnit,*) 'MGD calling MPAS_io_put_var now...'
            call MPAS_io_put_var(stream % fileHandle, field_cursor % char0dField % fieldName, field_cursor % char0dField % scalar, io_err)
            call MPAS_io_err_mesg(io_err, .false.)
            if (io_err /= MPAS_IO_NOERR .and. present(ierr)) ierr = MPAS_IO_ERR

         else if (field_cursor % field_type == FIELD_1D_CHAR) then
         end if
         field_cursor => field_cursor % next
      end do


      !
      ! Sync all fields with disk
      !
      call MPAS_io_sync(stream % fileHandle)

   end subroutine MPAS_writeStream_MPAS2

end module mpas_io_restart
