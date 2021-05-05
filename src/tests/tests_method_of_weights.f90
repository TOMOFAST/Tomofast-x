
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

!===============================================================================================
! Unit tests for methods of weights (equality-constrained LSQR).
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2016.
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

  call matrix_A%initialize(nrows_A, ncols_loc * nrows_A, myrank)
  call matrix_C%initialize(nrows_C, ncols_loc * nrows_C, myrank)

  allocate(b(nrows_A), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(y_LSE(ncols_loc), source=0._CUSTOM_REAL, stat=ierr)
  allocate(g(nrows_C), source=0._CUSTOM_REAL, stat=ierr)

  ! Build matrix C. ---------------
  call matrix_C%new_row(myrank)

  call matrix_C%add(1.d0, 1, myrank)
  call matrix_C%add(- 1.d0, 2, myrank)

  call matrix_C%finalize(ncols_loc, myrank)

  ! Build matrix A. --------------
  call matrix_A%new_row(myrank)

  call matrix_A%add(1.d0, 1, myrank)
  call matrix_A%add(2.d0, 2, myrank)

  call matrix_A%new_row(myrank)

  call matrix_A%add(3.d0, 1, myrank)
  call matrix_A%add(4.d0, 2, myrank)

  call matrix_A%add_matrix(matrix_C, tau, myrank)

  call matrix_A%finalize(ncols_loc, myrank)

  ! Build right-hand sides. ---------
  b(1) = 1.d0
  b(2) = 1.d0

  g(1) = 2.d0

  b(3:3) = tau * g

  !----------------------------------------------------------
  ! Solve: find y = y(tau).
  call lsqr_solve(par%niter, par%rmin, par%gamma, matrix_A, b, y, myrank)

  ! Apply method of weights to correct y.
  call apply_method_of_weights(par, niter, matrix_A, matrix_C, y, b, g, tau, 3, myrank)

  !----------------------------------------------------------
  ! Compare results with analytical solution y_LSE.

  y_LSE(1) =  (1.d0 / 29.d0) * 39.d0
  y_LSE(2) =  (1.d0 / 29.d0) * (- 19.d0)

!  call assert_true(abs(y(1) - y_LSE(1)) / abs(y_LSE(1)) < 1.e-15_CUSTOM_REAL, "y(1) is not correct in test_method_of_weights_1.")
!  call assert_true(abs(y(2) - y_LSE(2)) / abs(y_LSE(2)) < 1.e-15_CUSTOM_REAL, "y(2) is not correct in test_method_of_weights_1.")

  ! Somehow the exact solution is found (with gcc 4.9.4). Not sure if this is true for other machines/compilers...
  ! If this fails, use the above version.
  call assert_true(abs(y(1) - y_LSE(1)) / abs(y_LSE(1)) == 0._CUSTOM_REAL, "y(1) is not correct in test_method_of_weights_1.")
  call assert_true(abs(y(2) - y_LSE(2)) / abs(y_LSE(2)) == 0._CUSTOM_REAL, "y(2) is not correct in test_method_of_weights_1.")

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

  call matrix_A%initialize(nrows_A, ncols_loc * nrows_A, myrank)
  call matrix_C%initialize(nrows_C, ncols_loc * nrows_C, myrank)

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

  call matrix_C%finalize(ncols_loc, myrank)

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

  call matrix_A%finalize(ncols_loc, myrank)

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
  call lsqr_solve(par%niter, par%rmin, par%gamma, matrix_A, b, y, myrank)

  ! Apply method of weights to correct y.
  call apply_method_of_weights(par, niter, matrix_A, matrix_C, y, b, g, tau, 3, myrank)

  !----------------------------------------------------------
  ! Compare results with analytical solution y_LSE.

  y_LSE(1) =  (1.d0 / 8.d0) * 46.d0
  y_LSE(2) =  (1.d0 / 8.d0) * (- 2.d0)
  y_LSE(3) =  (1.d0 / 8.d0) * 12.d0

  call assert_true(abs(y(1) - y_LSE(1)) / abs(y_LSE(1)) < 1.e-15_CUSTOM_REAL, "y(1) is not correct in test_method_of_weights_2.")
  call assert_true(abs(y(2) - y_LSE(2)) / abs(y_LSE(2)) < 1.e-15_CUSTOM_REAL, "y(2) is not correct in test_method_of_weights_2.")
  call assert_true(abs(y(3) - y_LSE(3)) / abs(y_LSE(3)) < 1.e-15_CUSTOM_REAL, "y(3) is not correct in test_method_of_weights_2.")

  deallocate(b)
  deallocate(y)
  deallocate(y_LSE)
  deallocate(g)

end subroutine test_method_of_weights_2

end module tests_method_of_weights
