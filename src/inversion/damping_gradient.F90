
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

!===========================================================================================
! A class to (vertical) gradient minimization to inversion.
!
! Vitaliy Ogarko and Jeremie Giraud, UWA, CET, Australia, 2016-2017.
!===========================================================================================
module damping_gradient

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use model
  use parallel_tools
  use gradient
  use vector

  implicit none

  private

  type, public :: t_damping_gradient
    private

    ! Model dimensions.
    integer :: nx, ny, nz
    ! Number of model parameters (local per CPU and total).
    integer :: nelements, nelements_total
    ! Weight.
    real(kind=CUSTOM_REAL) :: beta
    ! Weight of the whole problem (damping + misfit) in joint inversion.
    real(kind=CUSTOM_REAL) :: problem_weight

    ! Cost.
    real(kind=CUSTOM_REAL) :: cost

    ! Gradient weight per cell.
    real(kind=CUSTOM_REAL), public, allocatable :: grad_weight(:)

  contains
    private

    procedure, public, pass :: initialize => damping_gradient_initialize

    procedure, public, pass :: add => damping_gradient_add
    procedure, public, pass :: get_cost => damping_gradient_get_cost

  end type t_damping_gradient
contains

!===========================================================================================
! Initialization.
!===========================================================================================
subroutine damping_gradient_initialize(this, beta, problem_weight, nx, ny, nz, nelements, myrank)
  class(t_damping_gradient), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: beta, problem_weight
  integer, intent(in) :: nx, ny, nz, nelements
  integer, intent(in) :: myrank
  integer :: ierr

  this%beta = beta
  this%problem_weight = problem_weight
  this%nx = nx
  this%ny = ny
  this%nz = nz
  this%nelements = nelements

  this%nelements_total = nx * ny * nz
  this%cost = 0.d0

  if (.not. allocated(this%grad_weight)) &
    allocate(this%grad_weight(this%nelements_total), source=1._CUSTOM_REAL, stat=ierr)

end subroutine damping_gradient_initialize

!==================================================================================================
! Returns the cost for every component (this is what we want to minimize).
!==================================================================================================
pure function damping_gradient_get_cost(this) result(res)
  class(t_damping_gradient), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = this%cost

end function damping_gradient_get_cost

!==================================================================================================================
! Adding gradient constraints to the matrix and right-hand-side.
!
! grad_weight(:) is an array of gradient weights, one per cell.
!==================================================================================================================
subroutine damping_gradient_add(this, model, grad_weight, column_weight, matrix, b_RHS, param_shift, direction, myrank, nbproc)
  class(t_damping_gradient), intent(inout) :: this
  type(t_model), intent(in) :: model
  real(kind=CUSTOM_REAL), intent(in) :: grad_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  integer, intent(in) :: param_shift
  integer, intent(in) :: direction
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(:)

  integer :: row_beg, row_end, nsmaller
  integer :: i, j, k, p
  integer :: ind(2)
  real(kind=CUSTOM_REAL) :: val(2), delta, gradient_val
  type(t_parallel_tools) :: pt

  type(t_gradient) :: grad
  type(t_vector) :: gradient_fwd, gradient_bwd

  ! Number of parameters on ranks smaller than current one.
  nsmaller = pt%get_nsmaller(this%nelements, myrank, nbproc)

  ! First matrix row (in the big matrix).
  row_beg = matrix%get_current_row_number() + 1

  ! Add matrix lines (Jacobian).
  do p = 1, this%nelements_total

    ! Grid (i,j,k)-index in the full grid.
    i = model%grid_full%i_(p)
    j = model%grid_full%j_(p)
    k = model%grid_full%k_(p)

    call matrix%new_row(myrank)


    gradient_fwd = grad%get_grad(model%val_full, model%grid_full, i, j, k, FWD_TYPE)
    gradient_bwd = grad%get_grad(model%val_full, model%grid_full, i, j, k, BWD_TYPE)

    ! NOTE: Use only gradient in one direction per time.

    ! Derivatives for FWD_TYPE, on the right boundary use BWD_TYPE.
    if (direction == 1) then
      delta = model%grid_full%get_hx()

      if (i /= this%nx) then
        ind(1) = model%grid_full%get_ind(i + 1, j, k) ! f(i + 1, j, k)
        ind(2) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
        gradient_val = gradient_fwd%x
      else
        ind(1) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
        ind(2) = model%grid_full%get_ind(i - 1, j, k) ! f(i - 1, j, k)
        gradient_val = gradient_bwd%x
      endif

    else if (direction == 2) then
      delta = model%grid_full%get_hy()

      if (j /= this%ny) then
        ind(1) = model%grid_full%get_ind(i, j + 1, k) ! f(i, j + 1, k)
        ind(2) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
        gradient_val = gradient_fwd%y
      else
        ind(1) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
        ind(2) = model%grid_full%get_ind(i, j - 1, k) ! f(i, j - 1, k)
        gradient_val = gradient_bwd%y
      endif

    else if (direction == 3) then
      delta = model%grid_full%get_hz()

      if (k /= this%nz) then
        ind(1) = model%grid_full%get_ind(i, j, k + 1) ! f(i, j, k + 1)
        ind(2) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
        gradient_val = gradient_fwd%z
      else
        ind(1) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
        ind(2) = model%grid_full%get_ind(i, j, k - 1) ! f(i, j, k - 1)
        gradient_val = gradient_bwd%z
      endif

    else
      call exit_MPI("Wrong direction in damping_gradient_add!", myrank, 0)
    endif

    val(1) = 1.d0 / delta
    val(2) = - val(1)

    ! Add elements to the matrix.
    do i = 1, 2
      if (ind(i) > nsmaller .and. ind(i) <= nsmaller + this%nelements) then
        ind(i) = ind(i) - nsmaller
        val(i) = val(i) * this%problem_weight * this%beta * column_weight(ind(i)) * grad_weight(p)

        call matrix%add(val(i), param_shift + ind(i), myrank)
      endif
    enddo

    ! Setting the right-hand side.
    b_RHS(matrix%get_current_row_number()) = - this%problem_weight * this%beta * gradient_val * grad_weight(p)
  enddo

  ! Last matrix row (in the big matrix).
  row_end = matrix%get_current_row_number()

  ! Calculate the cost.
  this%cost = sum(b_RHS(row_beg:row_end)**2)

  !if (myrank == 0) print *, 'Damping_gradient added lines: ', row_end - row_beg + 1

end subroutine damping_gradient_add

end module damping_gradient
