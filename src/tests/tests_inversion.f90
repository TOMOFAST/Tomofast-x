
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

!===============================================================================================
! Unit tests for inversion part.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===============================================================================================
module tests_inversion

  use global_typedefs
  use ftnunit

  use parameters_inversion
  use inversion_arrays
  use sparse_matrix
  use model
  use damping
  use cross_gradient

  implicit none

  private

  public :: test_add_damping_identity_matrix

  ! Cross-gradient tests.
  public :: test_cross_gradient_calculate_all

  private :: test_cross_gradient_calculate

contains

!=============================================================================================
! Testing add_damping_identity_matrix()
! by multiplying the identity damping matrix by a vector (1, 2, 3, ..., nelements_total).
!=============================================================================================
subroutine test_add_damping_identity_matrix(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix) :: isensit
  type(t_damping) :: damping
  type(t_parameters_inversion) :: par
  type(t_inversion_arrays) :: arr
  real(kind=CUSTOM_REAL), allocatable :: x(:)
  real(kind=CUSTOM_REAL), allocatable :: b_loc(:)
  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  integer :: i, ierr

  par%ndata = 0
  par%nx = 10
  par%ny = 72
  par%nz = 4
  par%nelements_total = par%nx * par%ny * par%nz
  par%nelements = par%nelements_total / nbproc
  par%alpha = 1._CUSTOM_REAL
  par%problem_weight = 1._CUSTOM_REAL
  par%compression_type = 0

  if (mod(par%nelements_total, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: nelements_total mod nbproc /= 0, skipping the test."
    return
  endif

  call arr%model%initialize(par%nelements, myrank, nbproc)

  allocate(arr%model_prior(par%nelements), source=0._CUSTOM_REAL)
  allocate(arr%column_weight(par%nelements), source=0._CUSTOM_REAL)
  allocate(arr%damping_weight(par%nelements), source=0._CUSTOM_REAL)

  allocate(x(par%nelements), source=0._CUSTOM_REAL)
  allocate(b_loc(par%nelements_total), source=0._CUSTOM_REAL)
  allocate(b(par%nelements_total), source=0._CUSTOM_REAL)
  allocate(b_RHS(par%nelements_total), source=0._CUSTOM_REAL)

  arr%column_weight   = 1._CUSTOM_REAL
  arr%damping_weight  = 1._CUSTOM_REAL
  par%norm_power      = 2._CUSTOM_REAL

  do i = 1, par%nelements
    ! Form a vector (1, 2, 3,..., nelements_total) split between CPUs.
    x(i) = (par%nelements * myrank) + dble(i)
  enddo

  call isensit%initialize(par%ndata(1) + par%nelements_total, &
                          par%ndata(1) * par%nelements + par%nelements, myrank)

  call damping%initialize(par%nelements, par%alpha(1), par%problem_weight(1), par%norm_power, &
                          par%compression_type, par%nx, par%ny, par%nz)

  ! Create an identity matrix.
  call damping%add(isensit, b_RHS, arr%column_weight, arr%damping_weight, &
                   arr%model, arr%model_prior, 0, myrank, nbproc)

  ! Store the index of last element.
  call isensit%finalize(par%nelements, myrank)

  ! Multiply the identity matrix by x.
  call isensit%mult_vector(x, b_loc)

  call mpi_allreduce(b_loc, b, par%nelements_total, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Test the solution b.
  if (myrank == 0) then
    do i = 1, par%nelements_total
      call assert_comparable_real(b(i), dble(i), tol, "I * x /= b in test_add_damping_identity_matrix.")
    enddo
  endif

  deallocate(x)
  deallocate(b_loc)
  deallocate(b)
  deallocate(b_RHS)

  deallocate(arr%column_weight)
  deallocate(arr%damping_weight)
  deallocate(arr%model_prior)

end subroutine test_add_damping_identity_matrix

!=============================================================================================
! Testing cross_gradient_calculate().
!=============================================================================================
subroutine test_cross_gradient_calculate_all(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  call test_cross_gradient_calculate(myrank, nbproc, 1)
  call test_cross_gradient_calculate(myrank, nbproc, 2)

end subroutine test_cross_gradient_calculate_all

!=============================================================================================
! Testing cross_gradient_calculate().
!=============================================================================================
subroutine test_cross_gradient_calculate(myrank, nbproc, derivative_type)
  integer, intent(in) :: myrank, nbproc
  integer, intent(in) :: derivative_type

  type(t_cross_gradient) :: cross_grad
  type(t_sparse_matrix) :: matrix
  type(t_model) :: model1
  type(t_model) :: model2
  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: column_weight1(:)
  real(kind=CUSTOM_REAL), allocatable :: column_weight2(:)
  integer :: nx, ny, nz
  integer :: nelements, nelements_total, ierr
  integer :: matrix_nel_loc, matrix_nel_glob
  integer :: i, j, k, p

  ! Note: changing these dimensions will affect the test result (matrix_nel_glob).
  nx = 20
  ny = 20
  nz = 144

  nelements_total = nx * ny * nz
  nelements = nelements_total / nbproc

  if (mod(nelements_total, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: nelements_total mod nbproc /= 0, skipping the test."
    return
  endif

  call model1%initialize(nelements_total, myrank, nbproc)
  call model1%init_grid(nx, ny, nz, myrank)

  call model2%initialize(nelements_total, myrank, nbproc)
  call model2%init_grid(nx, ny, nz, myrank)

  allocate(b_RHS(3 * nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  allocate(column_weight1(nelements), source=1._CUSTOM_REAL, stat=ierr)
  allocate(column_weight2(nelements), source=1._CUSTOM_REAL, stat=ierr)

  p = 0

  ! Initialize models and grids.
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        p = p + 1

        model1%val_full(p) = dble(i)
        model2%val_full(p) = dble(i + 1)

        model1%grid_full%X1(p) = dble(i)
        model1%grid_full%X2(p) = dble(i + 1)
        model1%grid_full%Y1(p) = dble(j)
        model1%grid_full%Y2(p) = dble(j + 1)
        model1%grid_full%Z1(p) = dble(k)
        model1%grid_full%Z2(p) = dble(k + 1)

        model1%grid_full%i_(p) = i
        model1%grid_full%j_(p) = j
        model1%grid_full%k_(p) = k

        model1%grid_full%ind(i, j, k) = p
      enddo
    enddo
  enddo

  model2%grid_full = model1%grid_full

  call cross_grad%initialize(nx, ny, nz, nelements, myrank)

  call matrix%initialize(3 * nelements_total, cross_grad%get_num_elements(derivative_type), myrank)

  call cross_grad%calculate(model1, model2, column_weight1, column_weight2, &
                            matrix, b_RHS, .true., derivative_type, myrank, nbproc)

  matrix_nel_loc = matrix%get_number_elements()

  call mpi_allreduce(matrix_nel_loc, matrix_nel_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Note, the results were obtained by running the test using single CPU.
  ! So we are testing if in parallel case a big matrix has the same number of elements.
  if (derivative_type == 1) then
    call assert_equal_int(matrix_nel_glob, 461592, "matrix_nel_glob /= 461592 in test_cross_gradient_calculate.")
  else
    call assert_equal_int(matrix_nel_glob, 457904, "matrix_nel_glob /= 457904 in test_cross_gradient_calculate.")
  endif

  deallocate(b_RHS)
  deallocate(column_weight1)
  deallocate(column_weight2)

end subroutine test_cross_gradient_calculate

end module tests_inversion
