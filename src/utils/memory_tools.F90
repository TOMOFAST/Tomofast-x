
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

!========================================================================
! Functions to profile memory usage.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!========================================================================
module memory_tools

  use global_typedefs

  implicit none

  private

  public :: get_max_mem_usage
  private :: get_max_mem_usage_proc

contains

!===========================================================================
! Returns the maximum RAM used by all ranks.
!===========================================================================
function get_max_mem_usage() result(memory)
  real(kind=CUSTOM_REAL) :: memory
  real(kind=CUSTOM_REAL) :: memory_loc
  integer :: ierr

  memory_loc = get_max_mem_usage_proc()

  call mpi_allreduce(memory_loc, memory, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Convert kB to GB.
  memory = memory / 1024**2

end function get_max_mem_usage

!===========================================================================
! Returns the maximum RAM used by a current rank.
!===========================================================================
function get_max_mem_usage_proc() result(value)
  character(len=80) :: line
  integer :: ios, fu, value
  value = -1

  open(newunit=fu, file='/proc/self/status', action='read')
  do
    read(fu, '(a)', iostat=ios) line
    if (ios /= 0) exit
    ! Use VmHWM for the maximum used RAM memory. Use VmRSS for the current memory.
    if(line(1:6) == 'VmHWM:') then
      read(line(7:), *) value
      exit
    endif
  enddo
  close(fu)
end function get_max_mem_usage_proc

end module memory_tools

