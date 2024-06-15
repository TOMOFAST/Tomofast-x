
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

!===========================================================================================
! A class to (vertical) gradient minimization to inversion.
!
! Vitaliy Ogarko and Jeremie Giraud, UWA, CET, Australia.
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

  contains
    private

    procedure, public, pass :: initialize => damping_gradient_initialize

    procedure, public, pass :: add => damping_gradient_add
    procedure, public, pass :: get_cost => damping_gradient_get_cost

  end type t_damping_gradient
contains

!=================================================================================================
! Initialization.
!=================================================================================================
subroutine damping_gradient_initialize(this, beta, problem_weight, nx, ny, nz, nelements)
  class(t_damping_gradient), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: beta, problem_weight
  integer, intent(in) :: nx, ny, nz, nelements

  this%beta = beta
  this%problem_weight = problem_weight
  this%nx = nx
  this%ny = ny
  this%nz = nz
  this%nelements = nelements

  this%nelements_total = nx * ny * nz
  this%cost = 0.d0
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
!==================================================================================================================
subroutine damping_gradient_add(this, model, column_weight, local_weight, matrix, nrows, &
                                b_RHS, param_shift, direction, icomp, myrank, nbproc)
  class(t_damping_gradient), intent(inout) :: this
  type(t_model), intent(inout) :: model
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)
  real(kind=CUSTOM_REAL), intent(in) :: local_weight(:)
  integer, intent(in) :: param_shift, nrows
  integer, intent(in) :: direction, icomp
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: b_RHS(nrows)

  integer :: nsmaller
  integer :: i, j, k, l, p
  integer :: ind(2)
  real(kind=CUSTOM_REAL) :: val(2), delta, gradient_val
  type(t_vector) :: gradient_fwd
  !type(t_vector) :: gradient_bwd

  ! Number of parameters on ranks smaller than current one.
  nsmaller = get_nsmaller(this%nelements, myrank, nbproc)

  this%cost = 0.d0

  ! Add matrix lines (Jacobian).
  p = 0
  do k = 1, this%nz
    do j = 1, this%ny
      do i = 1, this%nx
        p = p + 1

        gradient_fwd = get_grad(model%val_full(:, icomp), model%grid_full, i, j, k, FWD_TYPE)
        !gradient_bwd = get_grad(model%val_full(:, icomp), model%grid_full, i, j, k, BWD_TYPE)

        ! NOTE: Use only gradient in one direction per time.

        if (direction == 1) then
          delta = model%grid_full%get_hx(p)

          if (i /= this%nx) then
            ind(1) = model%grid_full%get_ind(i + 1, j, k) ! f(i + 1, j, k)
            ind(2) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
            gradient_val = gradient_fwd%x
          else
            !ind(1) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
            !ind(2) = model%grid_full%get_ind(i - 1, j, k) ! f(i - 1, j, k)
            !gradient_val = gradient_bwd%x
            call matrix%new_row(myrank)
            cycle
          endif

        else if (direction == 2) then
          delta = model%grid_full%get_hy(p)

          if (j /= this%ny) then
            ind(1) = model%grid_full%get_ind(i, j + 1, k) ! f(i, j + 1, k)
            ind(2) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
            gradient_val = gradient_fwd%y
          else
            !ind(1) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
            !ind(2) = model%grid_full%get_ind(i, j - 1, k) ! f(i, j - 1, k)
            !gradient_val = gradient_bwd%y
            call matrix%new_row(myrank)
            cycle
          endif

        else if (direction == 3) then
          delta = model%grid_full%get_hz(p)

          if (k /= this%nz) then
            ind(1) = model%grid_full%get_ind(i, j, k + 1) ! f(i, j, k + 1)
            ind(2) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
            gradient_val = gradient_fwd%z
          else
            !ind(1) = model%grid_full%get_ind(i, j, k)     ! f(i, j, k)
            !ind(2) = model%grid_full%get_ind(i, j, k - 1) ! f(i, j, k - 1)
            !gradient_val = gradient_bwd%z
            call matrix%new_row(myrank)
            cycle
          endif

        else
          call exit_MPI("Wrong direction in damping_gradient_add!", myrank, 0)
        endif

        val(1) = 1.d0 / delta
        val(2) = - val(1)

        ! Add elements to the matrix.
        do l = 1, 2
          if (ind(l) > nsmaller .and. ind(l) <= nsmaller + this%nelements) then
            ind(l) = ind(l) - nsmaller
            val(l) = val(l) * this%problem_weight * this%beta * column_weight(ind(l)) * local_weight(p)

            call matrix%add(val(l), param_shift + ind(l), myrank)
          endif
        enddo

        call matrix%new_row(myrank)

        ! Setting the right-hand side.
        b_RHS(matrix%get_current_row_number()) = - this%problem_weight * this%beta * gradient_val * local_weight(p)

        ! Calculate the cost.
        this%cost = this%cost + gradient_val**2
      enddo
    enddo
  enddo

end subroutine damping_gradient_add

end module damping_gradient
