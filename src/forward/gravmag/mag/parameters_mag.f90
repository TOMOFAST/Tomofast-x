
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
! A class that stores magnetism specific parameters,
! and inherits common grav/mag parameters from the base class.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015.
!========================================================================================
module parameters_mag

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use parameters_gravmag

  implicit none

  ! Main magnetism input parameters read from the Parfile.
  type, extends(t_parameters_base), public :: t_parameters_mag
    ! Angles are in degrees, with inclinations positive below
    ! horizontal and declinations positive east of true north.

    ! Magnetic field inclination.
    real(kind=CUSTOM_REAL) :: mi
    ! Magnetic field declination.
    real(kind=CUSTOM_REAL) :: md
    ! Ambient field inclination.
    real(kind=CUSTOM_REAL) :: fi
    ! Ambient field declination.
    real(kind=CUSTOM_REAL) :: fd
    ! X axis declination
    real(kind=CUSTOM_REAL) :: theta
    ! Ambient field intensity.
    real(kind=CUSTOM_REAL) :: intensity

  contains
    private

    procedure, public, pass :: broadcast => parameters_mag_broadcast

  end type t_parameters_mag

contains

!=========================================================================
! MPI broadcast parameters that are read from a Parfile.
!=========================================================================
subroutine parameters_mag_broadcast(this, myrank)
  class(t_parameters_mag), intent(in) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ! Broadcast the base parameters.
  call this%t_parameters_base%broadcast(myrank)

  ierr = 0

  call MPI_Bcast(this%mi, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%md, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%fi, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%fd, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%theta, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%intensity, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI_Bcast error in parameters_mag_broadcast!", myrank, ierr)

end subroutine parameters_mag_broadcast

end module parameters_mag
