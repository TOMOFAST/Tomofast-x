
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
! Unit tests for methods of weights (equality-constrained LSQR).
!
! Author: Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module tests_method_of_weights

  use global_typedefs
  use ftnunit

  use sparse_matrix
  use lsqr_solver
  use method_of_weights

  implicit none

  private

  ! Testing a well-conditioned system, whose matrix rows vary in norm.
  public :: test_method_of_weights_1

  ! Testing a well-conditioned system.
  public :: test_method_of_weights_2

contains

!=============================================================================================
! Use an example (small system) from Ref. [1], page 856 (top).
! [1] Charles Van Loan, On the Method of Weighting for Equality-Constrained Least-Squares Problems,
!     SIAM J. Numer. Anal., Vo. 22, No. 5, 1985.
!
! A = ( 1 2 ), b = ( 1 ), B = (1, - 1), d = (2).
!     ( 3 4 )      ( 1 )
!
!=============================================================================================
subroutine test_method_of_weights_1(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: y(:)
  real(kind=CUSTOM_REAL), allocatable :: y_LSE(:)
  real(kind=CUSTOM_REAL), allocatable :: g(:)
  type(t_sparse_matrix) :: matrix_A, matrix_C
  type(t_parameters_lsqr) :: par
  integer :: ncols, ncols_loc
  integer :: ierr
  integer :: nrows_A, nrows_C, niter
  real(kind=CUSTOM_REAL) :: tau

  if (nbproc > 1) then
    if (myrank == 0) print *, "WARNING: test_method_of_weights_1 is not yet adjusted for parallel runs, skipping the test."
    return
  endif

  tau = 1.d3
  niter = 3

  ncols = 2
  nrows_A = 3
  nrows_C = 1

  ncols_loc = ncols / nbproc

  par%rmin = 1.e-13_CUSTOM_REAL
  par%niter = 20
  par%gamma = 0._CUSTOM_REAL

  if (mod(ncols, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: ncols mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix_A%initialize(nrows_A, ncols_loc, int(ncols_loc * nrows_A, 8), myrank)
  call matrix_C%initialize(nrows_C, ncols_loc, int(ncols_loc * nrows_C, 8), myrank)

  allocate(b(nrows_A), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y_LSE(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(g(nrows_C), source=0._CUSTOM_REAL, stat=ierr)

  ! Build matrix C. ---------------
  call matrix_C%new_row(myrank)

  call matrix_C%add(1.d0, 1, myrank)
  call matrix_C%add(- 1.d0, 2, myrank)

  call matrix_C%finalize(myrank)

  ! Build matrix A. --------------
  call matrix_A%new_row(myrank)

  call matrix_A%add(1.d0, 1, myrank)
  call matrix_A%add(2.d0, 2, myrank)

  call matrix_A%new_row(myrank)

  call matrix_A%add(3.d0, 1, myrank)
  call matrix_A%add(4.d0, 2, myrank)

  call matrix_A%add_matrix(matrix_C, tau, myrank)

  call matrix_A%finalize(myrank)

  ! Build right-hand sides. ---------
  b(1) = 1.d0
  b(2) = 1.d0

  g(1) = 2.d0

  b(3:3) = tau * g

  !----------------------------------------------------------
  ! Solve: find y = y(tau).
  call lsqr_solve(size(b), size(y), par%niter, par%rmin, par%gamma, matrix_A, b, y, myrank)

  ! Apply method of weights to correct y.
  call apply_method_of_weights(par, niter, matrix_A, matrix_C, &
                               matrix_A%get_ncolumns(), matrix_A%get_total_row_number(), &
                               y, b, g, tau, 3, myrank)

  !----------------------------------------------------------
  ! Compare results with analytical solution y_LSE.

  y_LSE(1) =  (1.d0 / 29.d0) * 39.d0
  y_LSE(2) =  (1.d0 / 29.d0) * (- 19.d0)

  call assert_true(abs(y(1) - y_LSE(1)) / abs(y_LSE(1)) < 1.e-14_CUSTOM_REAL, "y(1) is not correct in test_method_of_weights_1.")
  call assert_true(abs(y(2) - y_LSE(2)) / abs(y_LSE(2)) < 1.e-14_CUSTOM_REAL, "y(2) is not correct in test_method_of_weights_1.")

  deallocate(b)
  deallocate(y)
  deallocate(y_LSE)
  deallocate(g)

end subroutine test_method_of_weights_1

!=============================================================================================
! Use an example of a well-conditioned system in Ref. [1], by Per-Ake Wedin, page 856.
! [1] Charles Van Loan, SIAM J. Numer. Anal., Vo. 22, No. 5, 1985.
!=============================================================================================
subroutine test_method_of_weights_2(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: y(:)
  real(kind=CUSTOM_REAL), allocatable :: y_LSE(:)
  real(kind=CUSTOM_REAL), allocatable :: g(:)
  type(t_sparse_matrix) :: matrix_A, matrix_C
  type(t_parameters_lsqr) :: par
  integer :: ncols, ncols_loc
  integer :: ierr
  integer :: nrows_A, nrows_C, niter
  real(kind=CUSTOM_REAL) :: tau

  if (nbproc > 1) then
    if (myrank == 0) print *, "WARNING: test_method_of_weights_2 is not yet adjusted for parallel runs, skipping the test."
    return
  endif

  tau = 1.d3
  niter = 3

  ncols = 3
  nrows_A = 6
  nrows_C = 2

  ncols_loc = ncols / nbproc

  par%rmin = 1.e-13_CUSTOM_REAL
  par%niter = 20
  par%gamma = 0._CUSTOM_REAL

  if (mod(ncols, nbproc) /= 0) then
    if (myrank == 0) print *, "WARNING: ncols mod nbproc /= 0, skipping the test."
    return
  endif

  call matrix_A%initialize(nrows_A, ncols_loc, int(ncols_loc * nrows_A, 8), myrank)
  call matrix_C%initialize(nrows_C, ncols_loc, int(ncols_loc * nrows_C, 8), myrank)

  allocate(b(nrows_A), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y_LSE(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(g(nrows_C), source=0._CUSTOM_REAL, stat=ierr)

  ! Build matrix C. ---------------
  call matrix_C%new_row(myrank)

  call matrix_C%add(1.d0, 1, myrank)
  call matrix_C%add(1.d0, 2, myrank)
  call matrix_C%add(1.d0, 3, myrank)

  call matrix_C%new_row(myrank)

  call matrix_C%add(1.d0, 1, myrank)
  call matrix_C%add(1.d0, 2, myrank)
  call matrix_C%add(- 1.d0, 3, myrank)

  call matrix_C%finalize(myrank)

  ! Build matrix A. --------------
  call matrix_A%new_row(myrank)

  call matrix_A%add(1.d0, 1, myrank)
  call matrix_A%add(1.d0, 2, myrank)
  call matrix_A%add(1.d0, 3, myrank)

  call matrix_A%new_row(myrank)

  call matrix_A%add(1.d0, 1, myrank)
  call matrix_A%add(3.d0, 2, myrank)
  call matrix_A%add(1.d0, 3, myrank)

  call matrix_A%new_row(myrank)

  call matrix_A%add(1.d0, 1, myrank)
  call matrix_A%add(- 1.d0, 2, myrank)
  call matrix_A%add(1.d0, 3, myrank)

  call matrix_A%new_row(myrank)

  call matrix_A%add(1.d0, 1, myrank)
  call matrix_A%add(1.d0, 2, myrank)
  call matrix_A%add(1.d0, 3, myrank)

  call matrix_A%add_matrix(matrix_C, tau, myrank)

  call matrix_A%finalize(myrank)

  ! Build right-hand sides. ---------
  b(1) = 1.d0
  b(2) = 2.d0
  b(3) = 3.d0
  b(4) = 4.d0

  g(1) = 7.d0
  g(2) = 4.d0

  b(5:6) = tau * g

  !----------------------------------------------------------
  ! Solve: find y = y(tau).
  call lsqr_solve(size(b), size(y), par%niter, par%rmin, par%gamma, matrix_A, b, y, myrank)

  ! Apply method of weights to correct y.
  call apply_method_of_weights(par, niter, matrix_A, matrix_C, &
                               matrix_A%get_ncolumns(), matrix_A%get_total_row_number(), &
                               y, b, g, tau, 3, myrank)

  !----------------------------------------------------------
  ! Compare results with analytical solution y_LSE.

  y_LSE(1) =  (1.d0 / 8.d0) * 46.d0
  y_LSE(2) =  (1.d0 / 8.d0) * (- 2.d0)
  y_LSE(3) =  (1.d0 / 8.d0) * 12.d0

  call assert_true(abs(y(1) - y_LSE(1)) / abs(y_LSE(1)) < 1.e-14_CUSTOM_REAL, "y(1) is not correct in test_method_of_weights_2.")
  call assert_true(abs(y(2) - y_LSE(2)) / abs(y_LSE(2)) < 1.e-14_CUSTOM_REAL, "y(2) is not correct in test_method_of_weights_2.")
  call assert_true(abs(y(3) - y_LSE(3)) / abs(y_LSE(3)) < 1.e-14_CUSTOM_REAL, "y(3) is not correct in test_method_of_weights_2.")

  deallocate(b)
  deallocate(y)
  deallocate(y_LSE)
  deallocate(g)

end subroutine test_method_of_weights_2

end module tests_method_of_weights
