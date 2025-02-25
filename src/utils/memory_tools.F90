
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

  public :: get_current_mem_usage
  private :: get_current_mem_usage_proc

contains

!===========================================================================
! Returns the currernt PSS memory used by all ranks.
! Use PSS memory to correctly handle shared memory usage.
!===========================================================================
function get_current_mem_usage() result(memory)
  real(kind=CUSTOM_REAL) :: memory
  real(kind=CUSTOM_REAL) :: memory_loc
  integer :: ierr

  memory_loc = get_current_mem_usage_proc()

  call mpi_allreduce(memory_loc, memory, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Convert kB to GB.
  memory = memory / 1024**2

end function get_current_mem_usage

!===========================================================================
! Returns the current PSS memory used by a current rank.
!===========================================================================
function get_current_mem_usage_proc() result(value)
  character(len=80) :: line
  integer :: ios, fu, value
  character(len=80) :: filename
  logical :: exists

  value = 0

  filename = '/proc/self/smaps_rollup'

  inquire(file=filename, exist=exists)

  if (exists) then
    open(newunit=fu, file=filename, action='read')
    do
      read(fu, '(a)', iostat=ios) line
      if (ios /= 0) exit
      ! Use Pss to account for shared memory usage.
      if(line(1:4) == 'Pss:') then
        read(line(5:), *) value
        exit
      endif
    enddo
    close(fu)
  endif
end function get_current_mem_usage_proc

end module memory_tools

