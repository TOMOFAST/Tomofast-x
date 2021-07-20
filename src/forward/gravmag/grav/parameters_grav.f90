
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

