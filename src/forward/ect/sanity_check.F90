
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

