
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

!==============================================================================================
! A class to add data misfit contribution (sensitivity matrix and right hand side)
! to the System of Linear Algebraic Equations (SLAE) that is stored
! using Compressed Sparse Row (CSR) format.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015.
!==============================================================================================
! TODO: rename to misfit
module sensitivity_matrix

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use model
  use sparse_matrix
  use sparse_matrix
  !use wavelet_transforms (DISABLED)

  implicit none

  private

  type, public :: t_sensitivity_matrix
    private

    ! Number of data.
    integer :: ndata

    ! Number of model parameters..
    integer :: nelements

    ! Weight of the whole problem (damping + misfit) in joint inversion.
    real(kind=CUSTOM_REAL) :: problem_weight

    ! Cost of the misfit objective functional phi_d.
    real(kind=CUSTOM_REAL) :: cost

  contains
    private

    procedure, public, pass :: initialize => sensitivity_matrix_initialize
    procedure, public, pass :: add => sensitivity_matrix_add
    procedure, public, pass :: get_cost => sensitivity_get_cost

    procedure, private, pass :: scale => sensitivity_matrix_scale
  end type t_sensitivity_matrix

contains

!===================================================================================================
! Initialization.
!===================================================================================================
subroutine sensitivity_matrix_initialize(this, ndata, nelements, problem_weight, myrank)
  class(t_sensitivity_matrix), intent(inout) :: this
  integer, intent(in) :: ndata, nelements
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  integer, intent(in) :: myrank

  this%ndata = ndata
  this%nelements = nelements
  this%problem_weight = problem_weight

end subroutine sensitivity_matrix_initialize

!==============================================================================================
! (DISABLED currently) 1. Compresses the sensitivity matrix using wavelets.
! 2. Adds the sensitivity matrix in the sparse format to the big matrix 'matrix'.
! 3. Initializes the corresponding rows of the right hand side (of the SLAE).
!
! Use legacy sensitivity matrix format, when the flag USE_LEGACY_SENSIT_MATRIX is true,
! to support the ECT inversion.
!==============================================================================================
subroutine sensitivity_matrix_add(this, matrix, b_RHS, param_shift, sensitivity, sensit_sparse, &
                                  model, column_weight, residuals, USE_LEGACY_SENSIT_MATRIX, myrank)
  class(t_sensitivity_matrix), intent(inout) :: this
  type(t_model), intent(in) :: model
  real(kind=CUSTOM_REAL), intent(in) :: sensitivity(:, :)
  type(t_sparse_matrix), intent(in) :: sensit_sparse
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: residuals(:)
  logical, intent(in) :: USE_LEGACY_SENSIT_MATRIX
  integer, intent(in) :: param_shift
  integer, intent(in) :: myrank

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  ! Local variables.
  integer :: i, p, ierr
  integer :: row_beg, row_end
  ! Sensitivity matrix row.
  real(kind=CUSTOM_REAL), allocatable :: line(:)

  allocate(line(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in sensitivity_matrix_add!", myrank, ierr)

  ! First matrix row (in the big matrix) of the sensitivity matrix that will be added.
  row_beg = matrix%get_current_row_number() + 1

  ! Loop on all the data lines.
  do i = 1, this%ndata

    ! Adding new matrix row.
    call matrix%new_row(myrank)

    ! Sensitivity matrix row.
    ! Copy to do not modify the original sensitivity matrix (by scaling),
    ! as we will need it to calculate the new (linear) data.
    if (USE_LEGACY_SENSIT_MATRIX) then
      line = sensitivity(:, i)

    else
      ! Extracting the i-th line from the sparse sensitivity matrix.
      call sensit_sparse%get_line(i, line)
    endif

    line = this%problem_weight * line

    ! Sensitivity matrix columns scaling.
    call this%scale(line, column_weight, myrank)

! DISABLED.
!    if (par%compress_matrix) then
!      ! Wavelet transform.
!      call FWT(line, par%nx, par%ny, par%nz)
!    endif

    do p = 1, this%nelements
      call matrix%add(line(p), param_shift + p, myrank)
    enddo
  enddo

  deallocate(line)

  ! Last matrix row (in the big matrix) of the added sensitivity matrix.
  row_end = matrix%get_current_row_number()

  ! Add corresponding contribution to the right-hand-side.
  b_RHS(row_beg:row_end) = this%problem_weight * residuals

  ! Calculate misfit function cost.
  this%cost = sum(b_RHS(row_beg:row_end)**2)

#ifdef USE_FLUSH6
  call flush(6)
#endif

end subroutine sensitivity_matrix_add

!================================================================================
! Scales the columns of the sensitivity matrix (excluding damping).
!================================================================================
subroutine sensitivity_matrix_scale(this, sensit_line, column_weight, myrank)
  class(t_sensitivity_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(inout) :: sensit_line(:)

  integer :: i

  do i = 1, this%nelements
    sensit_line(i) = sensit_line(i) * column_weight(i)
  enddo

end subroutine sensitivity_matrix_scale

!===========================================================================================
! Returns misfit objective functional cost (norm).
!===========================================================================================
pure function sensitivity_get_cost(this) result(res)
  class(t_sensitivity_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = this%cost

end function sensitivity_get_cost

end module sensitivity_matrix
