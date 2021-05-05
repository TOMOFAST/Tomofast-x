
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

!========================================================================================
! A class that stores gravity specific parameters,
! and inherits common grav/mag parameters from the base class.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015.
!========================================================================================
module parameters_grav

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parameters_gravmag

  implicit none

  ! Main gravity input parameters read from the Parfile.
  type, extends(t_parameters_base), public :: t_parameters_grav
    ! Add here gravity-only specific parameters.

  contains
    private

    procedure, public, pass :: broadcast => parameters_grav_broadcast

  end type t_parameters_grav

contains

!=========================================================================
! MPI broadcast parameters that are read from a Parfile.
!=========================================================================
subroutine parameters_grav_broadcast(this, myrank)
  class(t_parameters_grav), intent(in) :: this
  integer, intent(in) :: myrank

  ! Broadcast the base parameters.
  call this%t_parameters_base%broadcast(myrank)

end subroutine parameters_grav_broadcast

end module parameters_grav

