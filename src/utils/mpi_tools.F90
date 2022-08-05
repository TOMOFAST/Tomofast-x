
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

module mpi_tools

  use global_typedefs

  implicit none

  private

  public :: exit_MPI

contains

!===========================================================================
! Subroutine to exit MPI in case of a problem.
!===========================================================================
subroutine exit_MPI(message, myrank, ierror_code_from_sender)
  ! Error message.
  character(len=*), intent(in) :: message
  ! Rank of the calling process; error code.
  integer, intent(in) :: myrank, ierror_code_from_sender

  integer :: ierr

  print *
  print *
  print *, "***************************************************************"
  print *, "***************************************************************"
  print *, "Error number ", ierror_code_from_sender, " from proc ", myrank, ": ", trim(message)
  print *, "***************************************************************"
  print *, "***************************************************************"
  print *
  print *

#ifdef USE_FLUSH6
  flush(6)
#endif

  call MPI_Abort(MPI_COMM_WORLD, ierror_code_from_sender, ierr)

end subroutine exit_MPI

end module mpi_tools

