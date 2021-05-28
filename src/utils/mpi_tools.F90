
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module mpi_tools

  use global_typedefs

  implicit none

  private

  public :: exit_MPI
  public :: mpisendrecv

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
  call flush(6)
#endif

  call MPI_Abort(MPI_COMM_WORLD, ierror_code_from_sender, ierr)

end subroutine exit_MPI

!=========================================================================================================
! Exchanges the boundary values, for e.g., sparse matrix to vector multiplication in Conjugate Gradients.
!
! Algorithm to do the exchange: classical "even ones send, odd ones receive,
! and vise versa", with some extra care for the first and last subdomain,
! recall that we are doing MPI among entire (x,y) chunks along the z axis only
!
! The data that needs to be send is always nr*ntheta values, stored contiguously in memory.
!=========================================================================================================
subroutine mpisendrecv(x, nr, ntheta, nz, myrank, nbproc, ierr)

  ! problem dimensions
  integer , intent(in):: nr,ntheta,nz
  ! local part of the distributed vector to perform neighbour exchanges on
  real(kind=CUSTOM_REAL), intent(inout) :: x(0:nr+1,0:ntheta+1,0:nz+1)
  ! MPI rank and total number of processes
  integer, intent(in) :: myrank,nbproc
  ! MPI error code
  integer, intent(in) :: ierr

  integer :: nb_of_elems_to_send
  ! MPI status variable (needed for calls to mpi_recv
  integer :: mpi_stat(MPI_STATUS_SIZE)


  ! only do actual work if there is more than one rank total
  ! otherwise, this routine is empty
  ! DG DG this check is also performed at each call of the subroutine (where it
  ! DG DG is actually cheaper, unless this routine gets inlined anyway)
  if(nbproc > 1) then

    nb_of_elems_to_send = (nr+2)*(ntheta+2)

    !! DK DK and DG DG: there might be room for improvement here by switching to non-blocking MPI
    !! DK DK and DG DG: (if we can overlap communications with calculations; let us see if and how)

    if ( mod(myrank,2) == 0 ) then
      if (myrank < nbproc-1) then
        call mpi_send(x(0,0,nz),   nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 42, MPI_COMM_WORLD, ierr)
        call mpi_recv(x(0,0,nz+1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 43, MPI_COMM_WORLD, mpi_stat, ierr)
      endif

      if (myrank > 0) then
        call mpi_send(x(0,0,1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 44, MPI_COMM_WORLD, ierr)
        call mpi_recv(x(0,0,0), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 45, MPI_COMM_WORLD, mpi_stat, ierr)
      endif

    else

      if (myrank > 0) then
        call mpi_recv(x(0,0,0), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 42, MPI_COMM_WORLD, mpi_stat, ierr)
        call mpi_send(x(0,0,1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 43, MPI_COMM_WORLD, ierr)
      endif

      if (myrank < nbproc-1) then
        call mpi_recv(x(0,0,nz+1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 44, MPI_COMM_WORLD, mpi_stat, ierr)
        call mpi_send(x(0,0,nz),   nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 45, MPI_COMM_WORLD, ierr)
      endif

    endif

  endif

end subroutine mpisendrecv

end module mpi_tools

