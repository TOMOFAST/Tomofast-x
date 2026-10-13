
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!           Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.
!
!               (c) 2021 The University of Western Australia.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!================================================================================================
! A single-process replacement for the subset of MPI that Tomofast-x uses.
! It is compiled only when MPI is not available (make USE_MPI=0, which defines NO_MPI).
! The parallelism is then provided by OpenMP threads (make USE_OPENMP=1).
!
! The routines are external procedures (not module procedures), because the MPI buffer
! arguments accept any type and rank; this needs -fallow-argument-mismatch with GCC 10+.
!================================================================================================

!================================================================================================
! The MPI constants. This module is used by global_typedefs.
!================================================================================================
module mpi_stub_defs

  use, intrinsic :: iso_c_binding, only: c_int8_t

  implicit none

  private

  public :: MPI_ADDRESS_KIND
  public :: MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, MPI_INFO_NULL
  public :: MPI_INTEGER, MPI_INT, MPI_INTEGER8, MPI_REAL, MPI_DOUBLE_PRECISION
  public :: MPI_CHARACTER, MPI_CHAR, MPI_LOGICAL
  public :: MPI_SUM, MPI_MAX
  public :: MPI_IN_PLACE

  integer, parameter :: MPI_ADDRESS_KIND = selected_int_kind(18)

  integer, parameter :: MPI_COMM_WORLD = 1
  integer, parameter :: MPI_COMM_TYPE_SHARED = 1
  integer, parameter :: MPI_INFO_NULL = 0

  ! The data type codes (see type_size() in mpi_stub_state for their sizes).
  integer, parameter :: MPI_INTEGER = 1
  integer, parameter :: MPI_INT = 1
  integer, parameter :: MPI_INTEGER8 = 2
  integer, parameter :: MPI_REAL = 3
  integer, parameter :: MPI_DOUBLE_PRECISION = 4
  integer, parameter :: MPI_CHARACTER = 5
  integer, parameter :: MPI_CHAR = 5
  integer, parameter :: MPI_LOGICAL = 6

  ! The reduction operations (not used with one process).
  integer, parameter :: MPI_SUM = 1
  integer, parameter :: MPI_MAX = 2

  ! Only the address of this variable matters: it marks the "in place" buffer argument.
  integer(kind=c_int8_t), target, save :: MPI_IN_PLACE = 0_c_int8_t

end module mpi_stub_defs

!================================================================================================
! Helpers for the stub routines.
!================================================================================================
module mpi_stub_state

  use, intrinsic :: iso_c_binding
  use mpi_stub_defs

  implicit none

  public

  integer, parameter :: MAX_WINDOWS = 4096

  ! Memory of the "shared memory windows".
  type(c_ptr), save :: win_ptr(MAX_WINDOWS)
  logical, save :: win_used(MAX_WINDOWS) = .false.

  interface
    function c_calloc(n, sz) bind(C, name='calloc') result(p)
      import :: c_ptr, c_size_t
      integer(c_size_t), value :: n, sz
      type(c_ptr) :: p
    end function c_calloc

    subroutine c_free(p) bind(C, name='free')
      import :: c_ptr
      type(c_ptr), value :: p
    end subroutine c_free
  end interface

contains

!> Size in bytes of one element of an MPI data type.
function type_size(dtype) result(res)
  integer, intent(in) :: dtype
  integer(c_int64_t) :: res

  select case (dtype)
  case (MPI_INTEGER, MPI_REAL, MPI_LOGICAL)
    res = 4
  case (MPI_INTEGER8, MPI_DOUBLE_PRECISION)
    res = 8
  case default
    res = 1
  end select
end function type_size

!> Address of a buffer.
function addr(x) result(res)
  integer(c_int8_t), target, intent(in) :: x(*)
  integer(c_intptr_t) :: res

  res = transfer(c_loc(x(1)), res)
end function addr

!> Address of the MPI_IN_PLACE marker.
function in_place_addr() result(res)
  integer(c_intptr_t) :: res

  res = transfer(c_loc(MPI_IN_PLACE), res)
end function in_place_addr

!> Copy nbytes bytes from src (starting at byte src_off) to dst (starting at byte dst_off).
subroutine copy_bytes(src, dst, nbytes, src_off, dst_off)
  integer(c_int8_t), intent(in) :: src(*)
  integer(c_int8_t), intent(inout) :: dst(*)
  integer(c_int64_t), intent(in) :: nbytes, src_off, dst_off

  integer(c_int64_t) :: i

  do i = 1, nbytes
    dst(dst_off + i) = src(src_off + i)
  enddo
end subroutine copy_bytes

end module mpi_stub_state

!================================================================================================
! Environment.
!================================================================================================
subroutine MPI_Init(ierror)
  implicit none
  integer :: ierror
  ierror = 0
end subroutine MPI_Init

subroutine MPI_Finalize(ierror)
  implicit none
  integer :: ierror
  ierror = 0
end subroutine MPI_Finalize

subroutine MPI_Abort(comm, errorcode, ierror)
  implicit none
  integer :: comm, errorcode, ierror
  ierror = 0
  error stop "MPI_Abort called (build without MPI)."
end subroutine MPI_Abort

subroutine MPI_Comm_rank(comm, rank, ierror)
  implicit none
  integer :: comm, rank, ierror
  rank = 0
  ierror = 0
end subroutine MPI_Comm_rank

subroutine MPI_Comm_size(comm, size, ierror)
  implicit none
  integer :: comm, size, ierror
  size = 1
  ierror = 0
end subroutine MPI_Comm_size

subroutine MPI_Comm_split_type(comm, split_type, key, info, newcomm, ierror)
  implicit none
  integer :: comm, split_type, key, info, newcomm, ierror
  newcomm = comm
  ierror = 0
end subroutine MPI_Comm_split_type

subroutine MPI_Barrier(comm, ierror)
  implicit none
  integer :: comm, ierror
  ierror = 0
end subroutine MPI_Barrier

!================================================================================================
! Collective operations. With one process, the result is only a copy of the send buffer.
!================================================================================================
subroutine MPI_Bcast(buffer, count, datatype, root, comm, ierror)
  implicit none
  integer :: buffer(*)
  integer :: count, datatype, root, comm, ierror
  ierror = 0
end subroutine MPI_Bcast

subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_state
  implicit none
  integer(c_int8_t), target :: sendbuf(*), recvbuf(*)
  integer :: count, datatype, op, comm, ierror

  ierror = 0
  if (addr(sendbuf) == in_place_addr() .or. addr(sendbuf) == addr(recvbuf)) return

  call copy_bytes(sendbuf, recvbuf, int(count, c_int64_t) * type_size(datatype), 0_c_int64_t, 0_c_int64_t)
end subroutine MPI_Allreduce

subroutine MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_state
  implicit none
  integer(c_int8_t), target :: sendbuf(*), recvbuf(*)
  integer :: sendcount, sendtype, recvcount, recvtype, comm, ierror

  ierror = 0
  if (addr(sendbuf) == in_place_addr() .or. addr(sendbuf) == addr(recvbuf)) return

  call copy_bytes(sendbuf, recvbuf, int(sendcount, c_int64_t) * type_size(sendtype), 0_c_int64_t, 0_c_int64_t)
end subroutine MPI_Allgather

subroutine MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_state
  implicit none
  integer(c_int8_t), target :: sendbuf(*), recvbuf(*)
  integer :: sendcount, sendtype, recvcounts(*), displs(*), recvtype, root, comm, ierror

  integer(c_int64_t) :: off

  ierror = 0
  off = int(displs(1), c_int64_t) * type_size(recvtype)
  if (addr(sendbuf) == in_place_addr() .or. addr(sendbuf) == addr(recvbuf) + off) return

  call copy_bytes(sendbuf, recvbuf, int(sendcount, c_int64_t) * type_size(sendtype), 0_c_int64_t, off)
end subroutine MPI_Gatherv

subroutine MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_state
  implicit none
  integer(c_int8_t), target :: sendbuf(*), recvbuf(*)
  integer :: sendcount, sendtype, recvcounts(*), displs(*), recvtype, comm, ierror

  integer(c_int64_t) :: off

  ierror = 0
  off = int(displs(1), c_int64_t) * type_size(recvtype)
  if (addr(sendbuf) == in_place_addr() .or. addr(sendbuf) == addr(recvbuf) + off) return

  call copy_bytes(sendbuf, recvbuf, int(sendcount, c_int64_t) * type_size(sendtype), 0_c_int64_t, off)
end subroutine MPI_Allgatherv

subroutine MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_state
  implicit none
  integer(c_int8_t), target :: sendbuf(*), recvbuf(*)
  integer :: sendcounts(*), displs(*), sendtype, recvcount, recvtype, root, comm, ierror

  integer(c_int64_t) :: off

  ierror = 0
  off = int(displs(1), c_int64_t) * type_size(sendtype)
  if (addr(recvbuf) == in_place_addr() .or. addr(recvbuf) == addr(sendbuf) + off) return

  call copy_bytes(sendbuf, recvbuf, int(recvcount, c_int64_t) * type_size(recvtype), off, 0_c_int64_t)
end subroutine MPI_Scatterv

!================================================================================================
! Shared memory windows (one process owns all the memory).
!================================================================================================
subroutine MPI_Win_allocate_shared(wsize, disp_unit, info, comm, baseptr, win, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_defs, only: MPI_ADDRESS_KIND
  use mpi_stub_state
  implicit none
  integer(MPI_ADDRESS_KIND) :: wsize
  integer :: disp_unit, info, comm, win, ierror
  type(c_ptr) :: baseptr

  integer :: i

  win = 0
  ierror = 1

  do i = 1, MAX_WINDOWS
    if (.not. win_used(i)) then
      ! The memory is zeroed, as an MPI shared window usually is.
      win_ptr(i) = c_calloc(1_c_size_t, int(max(wsize, 1_MPI_ADDRESS_KIND), c_size_t))
      if (c_associated(win_ptr(i))) then
        win_used(i) = .true.
        baseptr = win_ptr(i)
        win = i
        ierror = 0
      endif
      return
    endif
  enddo
end subroutine MPI_Win_allocate_shared

subroutine MPI_Win_shared_query(win, rank, wsize, disp_unit, baseptr, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_defs, only: MPI_ADDRESS_KIND
  use mpi_stub_state
  implicit none
  integer :: win, rank, disp_unit, ierror
  integer(MPI_ADDRESS_KIND) :: wsize
  type(c_ptr) :: baseptr

  ierror = 1
  if (win >= 1 .and. win <= MAX_WINDOWS) then
    if (win_used(win)) then
      baseptr = win_ptr(win)
      ierror = 0
    endif
  endif
end subroutine MPI_Win_shared_query

subroutine MPI_Win_free(win, ierror)
  use, intrinsic :: iso_c_binding
  use mpi_stub_state
  implicit none
  integer :: win, ierror

  ierror = 1
  if (win >= 1 .and. win <= MAX_WINDOWS) then
    if (win_used(win)) then
      call c_free(win_ptr(win))
      win_ptr(win) = c_null_ptr
      win_used(win) = .false.
      win = 0
      ierror = 0
    endif
  endif
end subroutine MPI_Win_free
