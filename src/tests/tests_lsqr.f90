
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
! Unit tests for the Least Squares (LSQR) ans SCA solvers, parallelized by model (parameters).
! The matrix is stored using Compressed Sparse Row (CSR) format.
!
! Note: The LSQR solver is also tested in tests_method_of_weights module.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===============================================================================================
module tests_lsqr

  use global_typedefs
  use ftnunit

  use sparse_matrix
  use lsqr_solver
  use sca_solver

  implicit none

  private

  ! Testing determined system (# equations = # unknowns).
  public :: test_lsqr_determined

  ! Testing overdetermined system (# equations > # unknowns).
  public :: test_lsqr_overdetermined_1

  ! Testing overdetermined system (# equations > # unknowns).
  public :: test_lsqr_overdetermined_2

  ! Testing underdetermined system (# equations < # unknowns).
  public :: test_lsqr_underdetermined_1

  ! Testing underdetermined system (# equations < # unknowns).
  public :: test_lsqr_underdetermined_2

  ! Testing underdetermined system (# equations < # unknowns).
  public :: test_lsqr_underdetermined_3

contains

!=============================================================================================
! Perform test by building a matrix N x N:
!
!   1 1 1 ... 1
!   2 2 2 ... 2
!   3 3 3 ... 3
!   ...
!   N N N ... N
!
! And the right hand side vector b = (1N, 2N, 3N, ..., N^2).
! Solve LSQR for x (i.e., A * x = b).
! Compare the LSQR result with known x = (1, 1, 1, ..., 1), i.e., all ones.
!
! NOTE: this is a 'bad' matrix and it has many solutions, e.g. (N, 0, 0, ...., 0).
! The SCA method finds the (N, 0, 0, ...., 0) solution for this system.
!=============================================================================================
subroutine test_lsqr_determined(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:)
  type(t_sparse_matrix) :: matrix
  real(kind=CUSTOM_REAL) :: rmin, gamma
  integer :: nelements, nelements_total, niter
  integer :: i, j, ierr
  integer :: nrows

  nelements_total = 1440
  nrows = nelements_total
  rmin = 1.e-13_CUSTOM_REAL
  niter = 100
  gamma = 0._CUSTOM_REAL

  nelements = nelements_total / nbproc

  if (mod(nelements_total, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: nelements_total mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix%initialize(nelements_total, int8(nelements * nrows), myrank)

  allocate(b_RHS(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(delta_model(nelements), source=0._CUSTOM_REAL, stat=ierr)

  ! Building the matrix with right hand side.
  do j = 1, nrows
    call matrix%new_row(myrank)

    do i = 1, nelements
      call matrix%add(dble(j), i, myrank)
    enddo

    b_RHS(j) = dble(j * nelements_total)
  enddo

  call matrix%finalize(nelements, myrank)

  call lsqr_solve(niter, rmin, gamma, matrix, b_RHS, delta_model, myrank)

  ! Check the result.
  do i = 1, nelements
    call assert_comparable_real(delta_model(i), 1._CUSTOM_REAL, tol, "delta_model(i) /= 1.0 in test_lsqr_determined.")
  enddo

  deallocate(b_RHS)
  deallocate(delta_model)

end subroutine test_lsqr_determined

!=============================================================================================
! Consider a regression with constant, linear and quadratic terms:
! f(x) = b1 + b2 * x + b3 * x^2
!
! We build a matrix 3 x N:
!
!   1  x1  x1^2
!   1  x2  x2^2
!   ...
!   1  xn  xn^2,
!
! with x_i = i/N, i = 1, ..., N.
!
! Apply LSQR method and compare to the solution x = (b1, b2, b3).
! An example from [1].
! [1] Least Squares Estimation, Sara A. van de Geer, Volume 2, pp. 1041-1045,
!     in Encyclopedia of Statistics in Behavioral Science, 2005.
!
! Note: this test can only be run on 1 and 3 CPUs.
!=============================================================================================
subroutine test_lsqr_overdetermined_1(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:)
  type(t_sparse_matrix) :: matrix
  real(kind=CUSTOM_REAL) :: rmin, gamma
  integer :: nelements, nelements_total, niter
  integer :: i, ierr
  integer :: nrows
  real(kind=CUSTOM_REAL) :: xi
  real(kind=CUSTOM_REAL) :: b(3)

  nelements_total = 3
  nrows = 1000
  rmin = 1.e-14_CUSTOM_REAL
  niter = 100
  gamma = 0._CUSTOM_REAL

  nelements = nelements_total / nbproc

  if (mod(nelements_total, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: nelements_total mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix%initialize(nrows, int8(nelements * nrows), myrank)

  allocate(b_RHS(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(delta_model(nelements), source=0._CUSTOM_REAL, stat=ierr)

  b(1) = 1._CUSTOM_REAL
  b(2) = - 3._CUSTOM_REAL
  b(3) = 0._CUSTOM_REAL

  ! Building the matrix.
  do i = 1, nrows
    call matrix%new_row(myrank)

    xi = dble(i) / dble(nrows)

    if (nbproc == 1) then
      call matrix%add(xi**0, 1, myrank)
      call matrix%add(xi**1, 2, myrank)
      call matrix%add(xi**2, 3, myrank)

    else if (nbproc == 3) then
      call matrix%add(xi**myrank, 1, myrank)
    endif

    b_RHS(i) = b(1) + b(2) * xi + b(3) * xi**2
  enddo

  call matrix%finalize(nelements, myrank)

  call lsqr_solve(niter, rmin, gamma, matrix, b_RHS, delta_model, myrank)

  ! Check the result.
  if (nbproc == 1) then
    call assert_comparable_real(delta_model(1), b(1), tol, "delta_model(i) /= b(i) in test_lsqr_overdetermined.")
    call assert_comparable_real(delta_model(2), b(2), tol, "delta_model(i) /= b(i) in test_lsqr_overdetermined.")
    call assert_true(abs(delta_model(3)) < tol, "delta_model(3) /= 0 in test_lsqr_overdetermined.")

  else if (nbproc == 3) then
    if (myrank == 0) &
      call assert_comparable_real(delta_model(1), b(1), tol, "delta_model(i) /= b(i) in test_lsqr_overdetermined.")
    if (myrank == 1) &
      call assert_comparable_real(delta_model(1), b(2), tol, "delta_model(i) /= b(i) in test_lsqr_overdetermined.")
    if (myrank == 2) &
      call assert_true(abs(delta_model(1)) < tol, "delta_model(3) /= 0 in test_lsqr_overdetermined.")
  endif

  deallocate(b_RHS)
  deallocate(delta_model)

end subroutine test_lsqr_overdetermined_1

!=============================================================================================
! A test from [1], Chapter 1.2, page 12, Eq.(1.2.13).
!
! [1] Carl Wunsch, The ocean circulation inverse problem, 1996.
!
! Note: this test can only be run on 1 and 3 CPUs.
!=============================================================================================
subroutine test_lsqr_overdetermined_2(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: Ax(:)
  real(kind=CUSTOM_REAL), allocatable :: x(:)
  type(t_sparse_matrix) :: matrix
  real(kind=CUSTOM_REAL) :: rmin, gamma
  integer :: niter
  integer :: i, j, ierr
  integer :: nrows, ncols, ncols_loc
  real(kind=CUSTOM_REAL) :: a(3, 5)
  real(kind=CUSTOM_REAL) :: tolerance_local

  rmin = 1.e-13_CUSTOM_REAL
  niter = 100
  gamma = 0._CUSTOM_REAL

  ncols = 3
  nrows = 5

  ncols_loc = ncols / nbproc

  if (mod(ncols, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: ncols mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix%initialize(nrows, int8(ncols_loc * nrows), myrank)

  allocate(b_RHS(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(x(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(Ax(nrows), source=0._CUSTOM_REAL, stat=ierr)

  a(1, 1) = 1.2550d0
  a(2, 1) = 1.6731d0
  a(3, 1) = - 1.3927d0

  a(1, 2) = 0.4891d0
  a(2, 2) = 0.0943d0
  a(3, 2) = - 0.7829d0

  a(1, 3) = - 0.1755d0
  a(2, 3) = 1.8612d0
  a(3, 3) = 1.0972d0

  a(1, 4) = 0.4189d0
  a(2, 4) = 0.2469d0
  a(3, 4) = - 0.5990d0

  a(1, 5) = - 0.2900d0
  a(2, 5) = 0.7677d0
  a(3, 5) = 0.8188d0

  ! Right hand side.
  b_RHS(1) = 0.3511d0
  b_RHS(2) = - 1.6710d0
  b_RHS(3) = 6.838d0
  b_RHS(4) = - 0.8843d0
  b_RHS(5) = 3.7018d0

  ! Building the sparse matrix.
  do j = 1, nrows
    call matrix%new_row(myrank)

    if (nbproc == 1) then
      do i = 1, ncols
        call matrix%add(a(i, j), i, myrank)
      enddo

    else if (nbproc == 3) then
      call matrix%add(a(myrank + 1, j), 1, myrank)
    endif
  enddo

  call matrix%finalize(ncols_loc, myrank)

  ! Solve the least squares problem.
  call lsqr_solve(niter, rmin, gamma, matrix, b_RHS, x, myrank)
  !call sca_solve(niter, rmin, matrix, b_RHS, x, myrank, nbproc)

  ! Note the solution given in the Wunsch's book seems to be not correct one.
  ! As it leads to a larger residual norm.
  ! Besides, my solution coincides with the one obtained using mathematics software using the following script:
  ! LeastSquares[{{1.2550, 1.6731, - 1.3927}, {0.4891, 0.0943, - 0.7829}, {- 0.1755, 1.8612, 1.0972}, {0.4189, 0.2469, - 0.5990}, {- 0.2900, 0.7677, 0.8188}}, {0.3511, - 1.6710, 6.838, - 0.8843, 3.7018}]
  ! Mathematica output is: {157.611, -38.0747, 96.0291}

  if (MATRIX_PRECISION == SIZE_DOUBLE) then
  ! Double precision.
    tolerance_local = 1.e-3_CUSTOM_REAL
  else
  ! Single precision.
    tolerance_local = 1.e-2_CUSTOM_REAL
  endif

  if (nbproc == 1) then
    print *, 'x1 =', x(1)
    print *, 'x2 =', x(2)
    print *, 'x3 =', x(3)

    call assert_true(abs(x(1) - 157.611) < tolerance_local, "x(1) is wrong in test_lsqr_overdetermined_2.")
    call assert_true(abs(x(2) + 38.0747) < tolerance_local, "x(2) is wrong in test_lsqr_overdetermined_2.")
    call assert_true(abs(x(3) - 96.0291) < tolerance_local, "x(3) is wrong in test_lsqr_overdetermined_2.")

    call matrix%mult_vector(x, Ax)
    print *, 'norm 1 =', norm2(Ax - b_RHS)

    ! Testing the norm according to the C. Wunsch's book solution.
    x(1) = 1.d0
    x(2) = 2.d0
    x(3) = 3.d0

    call matrix%mult_vector(x, Ax)
    print *, 'norm 2 =', norm2(Ax - b_RHS)

  else if (nbproc == 3) then
    if (myrank == 0) call assert_true(abs(x(1) - 157.611) < tolerance_local, "x(1) is wrong in test_lsqr_overdetermined_2.")
    if (myrank == 1) call assert_true(abs(x(1) + 38.0747) < tolerance_local, "x(2) is wrong in test_lsqr_overdetermined_2.")
    if (myrank == 2) call assert_true(abs(x(1) - 96.0291) < tolerance_local, "x(3) is wrong in test_lsqr_overdetermined_2.")
  endif

  deallocate(b_RHS)
  deallocate(x)
  deallocate(Ax)

end subroutine test_lsqr_overdetermined_2

!=============================================================================================
! Solving the following system:
!
! x1 + x2 = 1,
! 2x1 + x2 - q = 0.
!
! Which has the minimum norm underdetermined solution x1 = 0, x2 = 1, q = 1.
! (See Carl Wunsch, The ocean circulation inverse problem, Eq.(3.4.120).)
!
! (!)NOTE: The SCA solver finds another solution (0.5, 0.5, 1.5), which is not minimum norm though.
!
! Note: this test can only be run on 1 and 3 CPUs.
!=============================================================================================
subroutine test_lsqr_underdetermined_1(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:)
  type(t_sparse_matrix) :: matrix
  real(kind=CUSTOM_REAL) :: rmin, gamma
  integer :: niter
  integer :: i, j, ierr
  integer :: nrows, ncols, ncols_loc
  real(kind=CUSTOM_REAL) :: a(3, 2)

  rmin = 1.e-13_CUSTOM_REAL
  niter = 100
  gamma = 0._CUSTOM_REAL

  ncols = 3
  nrows = 2

  ncols_loc = ncols / nbproc

  if (mod(ncols, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: ncols mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix%initialize(nrows, int8(ncols_loc * nrows), myrank)

  allocate(b_RHS(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(delta_model(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)

  a(1, 1) = 1._CUSTOM_REAL
  a(2, 1) = 1._CUSTOM_REAL
  a(3, 1) = 0._CUSTOM_REAL

  a(1, 2) = 2._CUSTOM_REAL
  a(2, 2) = 1._CUSTOM_REAL
  a(3, 2) = - 1._CUSTOM_REAL

  ! Right hand side.
  b_RHS(1) = 1._CUSTOM_REAL
  b_RHS(2) = 0._CUSTOM_REAL

  ! Building the matrix.
  do j = 1, nrows
    call matrix%new_row(myrank)

    if (nbproc == 1) then
      do i = 1, ncols
        call matrix%add(a(i, j), i, myrank)
      enddo

    else if (nbproc == 3) then
      call matrix%add(a(myrank + 1, j), 1, myrank)
    endif
  enddo

  call matrix%finalize(ncols_loc, myrank)

  !------------------------------------------------------------------
  ! Testing the LSQR solver.
  !------------------------------------------------------------------
  call lsqr_solve(niter, rmin, gamma, matrix, b_RHS, delta_model, myrank)

  ! Check the result.
  if (nbproc == 1) then
    call assert_true(abs(delta_model(1)) < 1.e-15_CUSTOM_REAL, "delta_model(1) /= 0 in test_lsqr_underdetermined.")
    call assert_comparable_real(delta_model(2), 1._CUSTOM_REAL, tol, "delta_model(i) /= x(i) in test_lsqr_underdetermined.")
    call assert_comparable_real(delta_model(3), 1._CUSTOM_REAL, tol, "delta_model(i) /= x(i) in test_lsqr_underdetermined.")

  else if (nbproc == 3) then
    if (myrank == 0) &
      call assert_true(abs(delta_model(1)) < 1.e-15_CUSTOM_REAL, "delta_model(1) /= 0 in test_lsqr_underdetermined.")
    if (myrank == 1) &
      call assert_comparable_real(delta_model(1), 1._CUSTOM_REAL, tol, "delta_model(i) /= x(i) in test_lsqr_underdetermined.")
    if (myrank == 2) &
      call assert_comparable_real(delta_model(1), 1._CUSTOM_REAL, tol, "delta_model(i) /= x(i) in test_lsqr_underdetermined.")
  endif

  !------------------------------------------------------------------
  ! Testing the SCA solver.
  !------------------------------------------------------------------
  ! (!) NOTE: It finds another correct solution: (0.5, 0.5, 1.5)
  if (nbproc == 1) then
    call sca_solve(100, rmin, matrix, b_RHS, delta_model, myrank, nbproc)

    call assert_comparable_real(delta_model(1), 0.5_CUSTOM_REAL, tol, "delta_model(2) /= x(1) in test_lsqr_underdetermined.")
    call assert_comparable_real(delta_model(2), 0.5_CUSTOM_REAL, tol, "delta_model(2) /= x(2) in test_lsqr_underdetermined.")
    call assert_comparable_real(delta_model(3), 1.5_CUSTOM_REAL, tol, "delta_model(3) /= x(3) in test_lsqr_underdetermined.")
  endif

  deallocate(b_RHS)
  deallocate(delta_model)

end subroutine test_lsqr_underdetermined_1

!=============================================================================================
! A test from the book of W. Menke, Geophysical Data Analysis: Descrete Inverse Theory, 1989, page 102.
!
!                          [m1]
! Gm = [1/4, 1/4, 1/4, 1/4][m2] = [d1]
!                          [m3]
!                          [m4]
!
! The minimum length solution is m = [d1, d1, d1, d1]'.
!
! Note: this test can only be run on 1, 2 or 4 CPUs.
!=============================================================================================
subroutine test_lsqr_underdetermined_2(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: x(:)
  type(t_sparse_matrix) :: matrix
  real(kind=CUSTOM_REAL) :: rmin, gamma
  real(kind=CUSTOM_REAL) :: val, d1
  integer :: niter
  integer :: i, j, ierr
  integer :: nrows, ncols, ncols_loc

  rmin = 1.e-14_CUSTOM_REAL
  niter = 100
  gamma = 0._CUSTOM_REAL

  nrows = 1
  ncols = 4

  ncols_loc = ncols / nbproc

  if (mod(ncols, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: ncols mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix%initialize(nrows, int8(ncols_loc * nrows), myrank)

  allocate(b_RHS(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(x(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)

  ! Building the matrix.
  val = 1.d0 / 4.d0
  do j = 1, nrows
    call matrix%new_row(myrank)

    do i = 1, ncols_loc
      call matrix%add(val, i, myrank)
    enddo
  enddo

  call matrix%finalize(ncols_loc, myrank)

  d1 = 1.d0
  b_RHS(1) = d1

  !------------------------------------------------------------------
  ! Testing the LSQR solver.
  !------------------------------------------------------------------
  call lsqr_solve(niter, rmin, gamma, matrix, b_RHS, x, myrank)

  do i = 1, ncols_loc
    call assert_comparable_real(x(i), d1, tol, "x(i) /= d1 in test_lsqr_underdetermined_2.")
  enddo

  !------------------------------------------------------------------
  ! Testing the SCA solver.
  !------------------------------------------------------------------
  ! (!) NOTE: It finds another correct solution: (4 * d1, 0, 0, 0)
  if (nbproc == 1) then
    call sca_solve(niter, rmin, matrix, b_RHS, x, myrank, nbproc)

    call assert_comparable_real(x(1), 4.d0 * d1, tol, "x(1) /= 4 * d1 in test_lsqr_underdetermined_2.")
    do i = 2, ncols
      call assert_comparable_real(x(i), 0.d0, tol, "x(i) /= 0 in test_lsqr_underdetermined_2.")
    enddo
  endif

  deallocate(b_RHS)
  deallocate(x)

end subroutine test_lsqr_underdetermined_2


!=============================================================================================
! A test from the book of Carl Wunsch, The ocean circulation inverse problem,
! an example after Eq.(3.6.26) on page 187.
!
! (1  1  1  1) x + n = ( 1 ),
! (1 -1 -1  1)         (-1 )
!
! where n is the noise.
!
! Note: this test can only be run on 1, 2 or 4 CPUs.
!=============================================================================================
subroutine test_lsqr_underdetermined_3(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b_RHS(:)
  real(kind=CUSTOM_REAL), allocatable :: x(:)
  real(kind=CUSTOM_REAL), allocatable :: n(:)
  type(t_sparse_matrix) :: matrix
  real(kind=CUSTOM_REAL) :: rmin, gamma
  real(kind=CUSTOM_REAL) :: a(4, 2)
  integer :: niter
  integer :: i, j, ierr
  integer :: nrows, ncols, ncols_loc

  rmin = 1.e-14_CUSTOM_REAL
  niter = 100
  gamma = 0._CUSTOM_REAL

  nrows = 2
  ncols = 4

  ncols_loc = ncols / nbproc

  if (mod(ncols, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: ncols mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix%initialize(nrows, int8(ncols_loc * nrows), myrank)

  allocate(b_RHS(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(x(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(n(nrows), source=0._CUSTOM_REAL, stat=ierr)

  a(:, 1) = 1.d0
  a(1, 2) = 1.d0
  a(2, 2) = - 1.d0
  a(3, 2) = - 1.d0
  a(4, 2) = 1.d0

  ! Building the matrix.
  do j = 1, nrows
    call matrix%new_row(myrank)

    do i = 1, ncols_loc
      call matrix%add(a(myrank * ncols_loc + i, j), i, myrank)
    enddo
  enddo

  call matrix%finalize(ncols_loc, myrank)

  b_RHS(1) = 1
  b_RHS(2) = - 1

  n = 0.0d0

  b_RHS = b_RHS - n

  !------------------------------------------------------------------
  ! Testing the LSQR solver.
  !------------------------------------------------------------------
  call lsqr_solve(niter, rmin, gamma, matrix, b_RHS, x, myrank)

  if (nbproc == 1) then
    call assert_comparable_real(x(1), 0.0d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
    call assert_comparable_real(x(2), 0.5d0, tol, "x(2) is wrong in test_lsqr_underdetermined_3.")
    call assert_comparable_real(x(3), 0.5d0, tol, "x(3) is wrong in test_lsqr_underdetermined_3.")
    call assert_comparable_real(x(4), 0.0d0, tol, "x(4) is wrong in test_lsqr_underdetermined_3.")

  else if (nbproc == 2) then
    if (myrank == 0) then
      call assert_comparable_real(x(1), 0.0d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
      call assert_comparable_real(x(2), 0.5d0, tol, "x(2) is wrong in test_lsqr_underdetermined_3.")
    else if (myrank == 1) then
      call assert_comparable_real(x(1), 0.5d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
      call assert_comparable_real(x(2), 0.0d0, tol, "x(2) is wrong in test_lsqr_underdetermined_3.")
    endif

  else if (nbproc == 4) then
    if (myrank == 0) then
      call assert_comparable_real(x(1), 0.0d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
    else if (myrank == 1) then
      call assert_comparable_real(x(1), 0.5d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
    else if (myrank == 2) then
      call assert_comparable_real(x(1), 0.5d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
    else if (myrank == 3) then
      call assert_comparable_real(x(1), 0.0d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
    endif
  endif

  !------------------------------------------------------------------
  ! Testing the SCA solver.
  !------------------------------------------------------------------
  ! (!) NOTE: It finds another correct solution: (0, 1, 0, 0)
  if (nbproc == 1) then
    x = 0.d0
    call sca_solve(niter, rmin, matrix, b_RHS, x, myrank, nbproc)

    call assert_comparable_real(x(1), 0.0d0, tol, "x(1) is wrong in test_lsqr_underdetermined_3.")
    call assert_comparable_real(x(2), 1.0d0, tol, "x(2) is wrong in test_lsqr_underdetermined_3.")
    call assert_comparable_real(x(3), 0.0d0, tol, "x(3) is wrong in test_lsqr_underdetermined_3.")
    call assert_comparable_real(x(4), 0.0d0, tol, "x(4) is wrong in test_lsqr_underdetermined_3.")
  endif

  deallocate(b_RHS)
  deallocate(x)
  deallocate(n)

end subroutine test_lsqr_underdetermined_3

end module tests_lsqr
