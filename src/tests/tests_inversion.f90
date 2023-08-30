
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
! Author: Vitaliy Ogarko, UWA, CET, Australia.
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
  type(t_model) :: model
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
  par%nmodel_components = 1

  if (mod(par%nelements_total, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: nelements_total mod nbproc /= 0, skipping the test."
    return
  endif

  call model%initialize(par%nelements, par%nmodel_components, .true., myrank, nbproc)

  allocate(arr%column_weight(par%nelements), source=0._CUSTOM_REAL)
  allocate(x(par%nelements), source=0._CUSTOM_REAL)
  allocate(b_loc(par%nelements_total), source=0._CUSTOM_REAL)
  allocate(b(par%nelements_total), source=0._CUSTOM_REAL)
  allocate(b_RHS(par%nelements_total), source=0._CUSTOM_REAL)

  arr%column_weight   = 1._CUSTOM_REAL
  par%norm_power      = 2._CUSTOM_REAL

  do i = 1, par%nelements
    ! Form a vector (1, 2, 3,..., nelements_total) split between CPUs.
    x(i) = (par%nelements * myrank) + dble(i)
  enddo

  call isensit%initialize(par%ndata(1) + par%nelements_total, par%nelements, &
                          int(par%ndata(1) * par%nelements + par%nelements, 8), myrank)

  call damping%initialize(par%nelements, par%alpha(1), par%problem_weight(1), par%norm_power, &
                          par%compression_type, par%nx, par%ny, par%nz)

  ! Create an identity matrix.
  call damping%add(isensit, isensit%get_total_row_number(), b_RHS, arr%column_weight, &
                   model%val(:, 1), model%val_prior(:, 1), 0, myrank, nbproc)

  ! Store the index of last element.
  call isensit%finalize(myrank)

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
  real(kind=CUSTOM_REAL) :: glob_weight
  integer :: nx, ny, nz, ncomponents
  integer :: nelements, nelements_total, ierr
  integer :: matrix_nel_loc, matrix_nel_glob
  integer :: i, j, k, p

  ! Note: changing these dimensions will affect the test result (matrix_nel_glob).
  nx = 20
  ny = 20
  nz = 144
  ncomponents = 1

  nelements_total = nx * ny * nz
  nelements = nelements_total / nbproc

  if (mod(nelements_total, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: nelements_total mod nbproc /= 0, skipping the test."
    return
  endif

  call model1%initialize(nelements_total, ncomponents, .true., myrank, nbproc)
  call model1%grid_full%allocate(nx, ny, nz, myrank)

  call model2%initialize(nelements_total, ncomponents, .true., myrank, nbproc)
  call model2%grid_full%allocate(nx, ny, nz, myrank)

  allocate(b_RHS(3 * nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  allocate(column_weight1(nelements), source=1._CUSTOM_REAL, stat=ierr)
  allocate(column_weight2(nelements), source=1._CUSTOM_REAL, stat=ierr)

  p = 0

  ! Initialize models and grids.
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        p = p + 1

        model1%val_full(p, 1) = dble(i)
        model2%val_full(p, 1) = dble(i + 1)

        model1%grid_full%X1(p) = real(i, CUSTOM_REAL)
        model1%grid_full%X2(p) = real(i + 1, CUSTOM_REAL)
        model1%grid_full%Y1(p) = real(j, CUSTOM_REAL)
        model1%grid_full%Y2(p) = real(j + 1, CUSTOM_REAL)
        model1%grid_full%Z1(p) = real(k, CUSTOM_REAL)
        model1%grid_full%Z2(p) = real(k + 1, CUSTOM_REAL)

        model1%grid_full%i_(p) = i
        model1%grid_full%j_(p) = j
        model1%grid_full%k_(p) = k

        model1%grid_full%ind(i, j, k) = p
      enddo
    enddo
  enddo

  model2%grid_full = model1%grid_full

  ! Scatter the full model as we update the full model inside cross_gradient_calculate() from its local parts.
  call model1%distribute(myrank, nbproc)
  call model2%distribute(myrank, nbproc)

  call cross_grad%initialize(nx, ny, nz, nelements, myrank)

  call matrix%initialize(3 * nelements_total, 2 * nelements, &
                         int(cross_grad%get_num_elements(derivative_type), 8), myrank)

  glob_weight = 1.d0

  call cross_grad%calculate(model1, model2, column_weight1, column_weight2, &
                            matrix, b_RHS, .true., derivative_type, glob_weight, myrank, nbproc)

  call matrix%finalize(myrank)

  matrix_nel_loc = int(matrix%get_number_elements())

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
