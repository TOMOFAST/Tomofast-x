
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

!================================================================================================
! A class that extends t_model_IO to work with parallel models (split between CPUs) for inversion.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!================================================================================================
module model

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use paraview
  use grid
  use string, only: str
  use parallel_tools
  use string
  use model_base
  use model_IO
  use sparse_matrix

  implicit none

  private

  type, extends(t_model_IO), public :: t_model

  contains
    private

    procedure, public, pass :: update => model_update
    procedure, public, pass :: calculate_data => model_calculate_data

  end type t_model

  public :: rescale_model

contains

!======================================================================================================
! Update model after inversion.
!======================================================================================================
subroutine model_update(this, delta_model)
  class(t_model), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: delta_model(:)

  integer :: i

  do i = 1, this%nelements
    this%val(i) = this%val(i) + delta_model(i)
  enddo

end subroutine model_update

!======================================================================================================
! Calculate the (linear) data using original (not scaled)
! sensitivity kernel (S) and model (m) as d = S * m.
!======================================================================================================
subroutine model_calculate_data(this, ndata, matrix_sensit, data, myrank)
  class(t_model), intent(in) :: this
  integer, intent(in) :: ndata, myrank
  type(t_sparse_matrix), intent(in) :: matrix_sensit

  real(kind=CUSTOM_REAL), intent(out) :: data(:)

  integer :: ierr

  call matrix_sensit%mult_vector(this%val, data)

  ! NOTE: Not sure this is correct way of calling the function (manual says MPI_IN_PLACE should be used on the root only).
  !       But this works somehow.
  call MPI_Allreduce(MPI_IN_PLACE, data, ndata, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(data, ndata, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in model_calculate_data!", myrank, ierr)

end subroutine model_calculate_data

!================================================================================================
! Weights the model parameters.
!================================================================================================
subroutine rescale_model(model, weight, nc)
  real(kind=CUSTOM_REAL), intent(inout) :: model(:)
  integer, intent(in) :: nc
  real(kind=CUSTOM_REAL), intent(in) :: weight(:)

  integer :: i

  do i = 1, nc
    model(i) = model(i) * weight(i)
  enddo
end subroutine rescale_model

end module model

