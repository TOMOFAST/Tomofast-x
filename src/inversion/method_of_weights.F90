
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
! Implementation of methods of weights that is applied to a SLAE,
! to enforce equality-constrained LSQR.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2016.
!===============================================================================================
module method_of_weights

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix
  use lsqr_solver

  implicit none

  private

  public :: apply_method_of_weights

contains

!===============================================================================================
! method = 1: Van Loan's method of weights [Charles Van Loan, SIAM J. Numer. Anal., Vo. 22, No. 5, 1985].
! method = 2: Modified Van Loan's method of weights [1].
! method = 3: A more robust version of method #2 [2].
!
! [1] Barlow, Jesse L., Error Analysis and Implementation Aspects of Deferred Correction
! for Equality Constrained Least Squares Problems,
! SIAM Journal on Numerical Analysis, 1988, Vol.25(6), pp.1340-1358.
!
! [2] J. L. Barlow and U. B. Vemulapati, A Note on Deferred Correction for Equality Constrained Least Squares Problems,
! SIAM Journal on Numerical Analysis, Vol. 29, No. 1 (Feb., 1992), pp. 249-256.
!
! The notation used here is the same as in Refs. [1, 2].
!
! Solve min|f - E y|, with equality constraints: C y = g.
!
! Do this by solving min|b - A y|, where
!
! A = (   E   )
!     ( tau C )
!
! b = (   f   )
!     ( tau g )
!
! This subroutine is unit-tested in tests_method_of_weights.f90.
!
!===============================================================================================
subroutine apply_method_of_weights(par, num_iter, matrix_A, matrix_C, y, b, g, tau, method, myrank)
  type(t_parameters_lsqr), intent(in) :: par
  integer, intent(in) :: num_iter
  type(t_sparse_matrix), intent(inout) :: matrix_A
  type(t_sparse_matrix), intent(in) :: matrix_C
  real(kind=CUSTOM_REAL), intent(inout) :: y(:)
  real(kind=CUSTOM_REAL), intent(inout) :: b(:)
  real(kind=CUSTOM_REAL), intent(in) :: g(:)
  real(kind=CUSTOM_REAL), intent(in) :: tau
  integer, intent(in) :: method, myrank

  ! Residual.
  real(kind=CUSTOM_REAL), allocatable :: w1(:)
  ! Perturbation.
  real(kind=CUSTOM_REAL), allocatable :: delta_y(:)
  ! Store here original b.
  real(kind=CUSTOM_REAL), allocatable :: b0(:)
  ! Store here C * delta_y.
  real(kind=CUSTOM_REAL), allocatable :: C_dy(:)
  ! A vector of Lagrange multipliers.
  real(kind=CUSTOM_REAL), allocatable :: lambda(:)

  integer :: k, ierr
  integer :: nrows, ncols
  integer :: b_size, ibeg, iend

  nrows = matrix_C%get_total_row_number()
  ncols = size(y)
  b_size = ubound(b, 1)

  ierr = 0

  ! Memory allocation.
  allocate(w1(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(delta_y(ncols), source=0._CUSTOM_REAL, stat=ierr)

  if (method == 2 .or. method == 3) then
    allocate(b0(b_size), source=0._CUSTOM_REAL, stat=ierr)
    allocate(lambda(nrows), source=0._CUSTOM_REAL, stat=ierr)
  endif

  if (method == 3) then
    allocate(C_dy(nrows), source=0._CUSTOM_REAL, stat=ierr)
  endif

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in method_of_weights_solve!", myrank, ierr)
  !---------------------------------------------------------------------------

  ibeg = b_size - nrows + 1
  iend = b_size

  if (method == 1) then
    b = 0._CUSTOM_REAL
  else
    b0 = b
    lambda = 0._CUSTOM_REAL
  endif

  ! Main loop.
  do k = 1, num_iter

    if (method == 3 .and. k > 1) then
      ! w1[k] = w1[k-1] - C delta_y[k-1].
      call matrix_C%mult_vector(delta_y, C_dy)
      w1 = w1 - C_dy
    else
      ! w1[k] = g - C y[k].
      call matrix_C%mult_vector(y, w1)
      w1 = g - w1
    endif

    if (myrank == 0) print *, '|w1|^2 = ', sum(w1**2)

    ! Construct the right hand side.
    if (method == 1) then
      !
      ! b = (     0     )
      !     ( tau w1[k] )
      !
      b(ibeg:iend) = tau * w1
    else
      lambda = lambda + tau**2 * w1

      !
      ! b = (             r[k]               )
      !     ( tau w1[k] + tau^{-1} lambda[k] )
      !
      ! Note: use b-array and matrix A, to calculate r[k].
      if (method == 3 .and. k > 1) then
        ! r[k] = r[k-1] - E delta_y[k-1].
        b0 = b
        call matrix_A%mult_vector(delta_y, b)
        b = b0 - b
      else
        ! r[k] = f - E y[k].
        call matrix_A%mult_vector(y, b)
        b = b0 - b
      endif

      if (myrank == 0) print *, '|r|^2 = ', sum(b(1:ibeg - 1)**2)

      b(ibeg:iend) = tau * w1 + (1.d0 / tau) * lambda
    end if

    ! Solve for new delta_y.
    delta_y = 0._CUSTOM_REAL
    call lsqr_solve(par%niter, par%rmin, par%gamma, matrix_A, b, delta_y, myrank)

    ! Update solution.
    y = y + delta_y

    print *, 'k, |dy[k]| / |y[k]| =', k, norm2(delta_y) / norm2(y)
  enddo

  if (allocated(w1)) deallocate(w1)
  if (allocated(delta_y)) deallocate(delta_y)
  if (allocated(b0)) deallocate(b0)
  if (allocated(lambda)) deallocate(lambda)
  if (allocated(C_dy)) deallocate(C_dy)

end subroutine apply_method_of_weights

end module method_of_weights
