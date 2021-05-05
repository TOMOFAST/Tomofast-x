
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

!==========================================================================
! Contains various helper functions to work with timers.
!==========================================================================
module timers

  use global_typedefs
  use stopwatch

  implicit none

  private

  public :: ReadWatchTime
  public :: get_CPU_times

contains

!==========================================================================
! Extract times from watchtype object.
!==========================================================================
subroutine ReadWatchTime(w, e, c, u)

  type(watchtype), intent(in) :: w
  real, intent(out) :: e, c, u

  integer :: ierr
  integer :: err
  real :: readval
  real :: t
  !real :: hours, minutes, sec

  call read_watch(readval, w, "wall", err=err)
  t = readval
  call mpi_allreduce(t, e, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  call read_watch(readval, w, "cpu", err=err)
  t = readval
  call mpi_allreduce(t, c, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  !hours = int(t / 3600.)
  !minutes = int((t - hours * 3600.) / 60.)
  !sec = t - hours * 3600. - minutes * 60.
  !write(83,*) hours,minutes,sec

  call read_watch(readval, w, "user", err=err)
  t = readval
  call mpi_allreduce(t, u, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine ReadWatchTime

!==========================================================================
! Extract accumulated and maximum CPU times from watchtype object.
!==========================================================================
subroutine get_CPU_times(watch, time_accum, time_min, time_max)

  type(watchtype), intent(in) :: watch
  real, intent(out) :: time_accum, time_min, time_max

  integer :: ierr
  real :: readval

  call read_watch(readval, watch, "cpu", err=ierr)

  call mpi_allreduce(readval, time_accum, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(readval, time_min, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(readval, time_max, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)

end subroutine get_CPU_times

end module timers

