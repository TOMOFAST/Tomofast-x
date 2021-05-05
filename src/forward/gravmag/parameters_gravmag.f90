
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
! Parameters needed for forward problems that are common for gravity and magnetism.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!========================================================================================
module parameters_gravmag

  use global_typedefs
  use mpi_tools, only: exit_MPI

  implicit none

  private

  type, public :: t_parameters_base

    ! Problem dimensions.
    integer :: nx, ny, nz
    ! The number of elements (on current CPU).
    integer :: nelements
    ! Number of data.
    integer :: ndata
    ! Number of data components (x, y, z, zz, yy, zz, etc).
    integer :: ncomponents
    ! File name for the data.
    character(len=256) :: data_file
    ! File name for the data grid.
    character(len=256) :: data_grid_file

    ! File names for the model (with grid).
    ! (1): Read model.
    ! (2): Prior model.
    ! (3): Starting model.
    character(len=256) :: model_files(3)

    ! Type of prior model.
    integer :: prior_model_type
    ! Number of prior models (for model type 2).
    integer :: number_prior_models
    ! Set prior model to this value.
    real(kind=CUSTOM_REAL) :: prior_model_val

    ! Type of starting model.
    integer :: start_model_type
    ! Set starting model to this value.
    real(kind=CUSTOM_REAL) :: start_model_val

    ! Flag to calculate data from prior model directly
    ! without storing sensitivity matrix.
    ! Then program stops after writing data to a file.
    integer :: calc_data_directly

    ! Type of the depth weighting (1-power law, 2-sensitivity column below the data, 3-integrated sensitivity).
    integer :: depth_weighting_type
    ! Power constant for depth weighting (type=1).
    real(kind=CUSTOM_REAL) :: beta
    ! Empirical constant for depth weighting.
    real(kind=CUSTOM_REAL) :: Z0

    ! ------ Matrix compression ---------------------------------------------
    ! Parameters for reduction of the memory requirements (to store the sensitivity matrix).
    ! Need to choose the compression rate accordingly to the distance_threshold.
    real(kind=CUSTOM_REAL) :: distance_threshold
    real(kind=CUSTOM_REAL) :: compression_rate

  contains
    private

    procedure, public, pass :: broadcast => parameters_base_broadcast
    procedure, public, pass :: get_nnz_compressed => parameters_base_get_nnz_compressed
  end type t_parameters_base

contains

!==============================================================================
! Returns the number of nonzero elements in the compressed sensitivity matrix.
!==============================================================================
function parameters_base_get_nnz_compressed(this) result(nnz)
  class(t_parameters_base), intent(in) :: this
  integer :: nnz

  nnz = int(this%compression_rate * this%nelements * this%ndata) + 1

  ! Sanity check.
  if (nnz <= 1) then
    print *, "Bad nnz value in get_nnz_compressed!"
    stop
  endif
end function

!=========================================================================
! MPI broadcast parameters that are read from a Parfile.
!=========================================================================
subroutine parameters_base_broadcast(this, myrank)
  class(t_parameters_base), intent(in) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  ierr = 0

  call MPI_Bcast(this%nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%ndata, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%prior_model_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%number_prior_models, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%prior_model_val, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%start_model_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%start_model_val, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%calc_data_directly, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%depth_weighting_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%Z0, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%beta, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(this%distance_threshold, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(this%compression_rate, 1, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI_Bcast error in parameters_base_broadcast!", myrank, ierr)

end subroutine parameters_base_broadcast

end module parameters_gravmag
