
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
!
! (c) CNRS, France, and University of Western Australia. January 2018
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module sanity_check

  use mpi_tools, only: exit_MPI

  implicit none

  private

  public :: sanity_nz
  public :: sanity_ntheta_nel

contains

!===========================================================================
! checks if given nz (number of partitions in z direction) is a multiple
! of the number of MPI processes
!===========================================================================
subroutine sanity_nz(nz, nbproc, myrank)

  ! number of slices in Z direction
  integer, intent(in) :: nz
  ! number of processes
  integer, intent(in) :: nbproc
  ! rank of current process
  integer, intent(in) :: myrank

  if (nbproc == 0) &
    call exit_MPI("nbproc=0; exiting...", myrank, 0)

  ! simple check
  if(mod(nz, nbproc) /= 0) then
    call exit_MPI("nz is not a multiple of nbproc; exiting...", myrank, 0)
  endif

end subroutine sanity_nz

!===========================================================================
! checks if given ntheta (number of finite volumes in theta direction) is
! a multiple of the number of electrodes.
!
! see module geometry on the various limitations that apply
!===========================================================================
subroutine sanity_ntheta_nel(ntheta, nel, myrank)

  ! number of elements in theta direction
  integer, intent(in) :: ntheta
  ! number of electrodes per ring
  integer, intent(in) :: nel
  ! rank of current process
  integer, intent(in) :: myrank

  ! simple check:
  if(mod(ntheta, nel) /= 0) then
    call exit_MPI("ntheta is not a multiple of nel; exiting...", myrank, 0)
  endif

end subroutine sanity_ntheta_nel

end module sanity_check

