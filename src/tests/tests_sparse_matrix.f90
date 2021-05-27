
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
! Unit tests for the t_sparse_matrix class.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia, 2016.
!===============================================================================================
module tests_sparse_matrix

  use global_typedefs
  use ftnunit

  use sparse_matrix

  implicit none

  private

  ! Testing retrieving the matrix element.
  public :: test_get_value

  ! Testing matrix columns normalization.
  public :: test_normalize_columns

  ! Testing sparse_matrix_trans_mult_matrix.
  public :: test_trans_mult_matrix

contains

!=============================================================================================
! Perform test for a (ncolumns x nrows) matrix with some zero columns.
!=============================================================================================
subroutine test_get_value(myrank)
  integer, intent(in) :: myrank

  type(t_sparse_matrix) :: matrix
  integer :: nrows, ncolumns
  integer :: i, j, counter

  ! Set matrix size.
  ncolumns = 10
  nrows = 30

  call matrix%initialize(nrows, ncolumns * nrows, myrank)

  ! Building the matrix.
  counter = 0
  do j = 1, nrows
    call matrix%new_row(myrank)
    do i = 1, ncolumns
      counter = counter + 1

      if (i <= ncolumns / 2) then
        call matrix%add(dble(counter), i, myrank)
      else
        ! These columns get zero values.
      endif
    enddo
  enddo

  call matrix%finalize(ncolumns, myrank)

  ! Testing matrix values.
  counter = 0
  do j = 1, nrows
    do i = 1, ncolumns
      counter = counter + 1

      if (i <= ncolumns / 2) then
        call assert_comparable_real(matrix%get_value(i, j), dble(counter), tol, "Wrong value in test_get_value!")
      else
        call assert_comparable_real(matrix%get_value(i, j), 0.d0, tol, "Wrong value in test_get_value!")
      endif
    enddo
  enddo

end subroutine test_get_value

!=============================================================================================
! Perform test for a (ncolumns x nrows) matrix with some zero columns.
!=============================================================================================
subroutine test_normalize_columns(myrank)
  integer, intent(in) :: myrank

  type(t_sparse_matrix) :: matrix
  integer :: nrows, ncolumns
  integer :: i, j, counter

  real(kind=CUSTOM_REAL), allocatable :: A(:, :)
  real(kind=CUSTOM_REAL), allocatable :: column_norm(:)
  real(kind=CUSTOM_REAL), allocatable :: ej(:, :), A_ej(:, :)

  ! Set matrix size.
  ncolumns = 10
  nrows = 30

  allocate(A(ncolumns, nrows))
  allocate(column_norm(ncolumns))
  allocate(ej(ncolumns, ncolumns))
  allocate(A_ej(ncolumns, nrows))

  ! Building the matrix.
  counter = 0
  do j = 1, nrows
    do i = 1, ncolumns
      counter = counter + 1

      if (i <= ncolumns / 2) then
        A(i, j) = dble(counter)
      else
        ! Set zero columns (all column elements are equal to zero).
        A(i, j) = 0.d0
      endif
    enddo
  enddo

  call matrix%initialize(nrows, ncolumns * nrows, myrank)

  ! Copy values from A to the sparse matrix.
  do j = 1, nrows
    call matrix%new_row(myrank)
    do i = 1, ncolumns
      call matrix%add(A(i, j), i, myrank)
    enddo
  enddo

  call matrix%finalize(ncolumns, myrank)

  call matrix%normalize_columns(column_norm)

  ! Basis vectors.
  do i = 1, ncolumns
    ej(i, :) = 0.d0
    ej(i, i) = 1.d0
  enddo

  do i = 1, ncolumns
    ! Compare matrix norms.
    call assert_comparable_real(column_norm(i), norm2(A(i, :)), tol, "Wrong column_norm array in test_normalize_columns!")

    ! Put the columns of (sparse) matrix into A_ej.
    call matrix%mult_vector(ej(i, :), A_ej(i, :))

    ! Test the norm of the (sparse) matrix columns.
    if (norm2(A(i, :)) /= 0.d0) then
      call assert_comparable_real(norm2(A_ej(i, :)), 1.d0, tol, "Wrong matrix column norm in test_normalize_columns!")
    else
      call assert_comparable_real(norm2(A_ej(i, :)), 0.d0, tol, "Wrong matrix column norm in test_normalize_columns!")
    endif
  enddo

  deallocate(A)
  deallocate(column_norm)
  deallocate(ej)
  deallocate(A_ej)

end subroutine test_normalize_columns

!=============================================================================================
! Perform test of A'A calculation for a matrix:
!     (1 2)
! A = (3 4)
!     (5 6)
!=============================================================================================
subroutine test_trans_mult_matrix(myrank)
  integer, intent(in) :: myrank

  type(t_sparse_matrix) :: matrix
  integer :: nrows, ncolumns
  integer :: i, j

  real(kind=CUSTOM_REAL), allocatable :: vec(:)
  real(kind=CUSTOM_REAL), allocatable :: vec2(:)
  real(kind=CUSTOM_REAL), allocatable :: Avi(:)
  real(kind=CUSTOM_REAL), allocatable :: Avj(:)
  real(kind=CUSTOM_REAL), allocatable :: H(:, :)

  ! Set matrix size.
  ncolumns = 2
  nrows = 3

  ! Allocate auxiliary arrays.
  allocate(vec(ncolumns))
  allocate(vec2(nrows))
  allocate(Avi(nrows))
  allocate(Avj(nrows))
  allocate(H(ncolumns, ncolumns))

  ! Building the matrix.
  call matrix%initialize(nrows, ncolumns * nrows, myrank)

  call matrix%new_row(myrank)
  call matrix%add(1.d0, 1, myrank)
  call matrix%add(2.d0, 2, myrank)

  call matrix%new_row(myrank)
  call matrix%add(3.d0, 1, myrank)
  call matrix%add(4.d0, 2, myrank)

  call matrix%new_row(myrank)
  call matrix%add(5.d0, 1, myrank)
  call matrix%add(6.d0, 2, myrank)

  call matrix%finalize(ncolumns, myrank)

  ! Calculate the resulting matrix H = A'A.
  do i = 1, ncolumns
    do j = 1, ncolumns
      H(i, j) = matrix%trans_mult_matrix(i, j, vec, Avi, Avj)
    enddo
  enddo

  call assert_comparable_real(H(1, 1), 35.d0, tol, "Wrong H[1, 1] in test_trans_mult_matrix!")
  call assert_comparable_real(H(2, 1), 44.d0, tol, "Wrong H[2, 1] in test_trans_mult_matrix!")
  call assert_comparable_real(H(1, 2), 44.d0, tol, "Wrong H[1, 2] in test_trans_mult_matrix!")
  call assert_comparable_real(H(2, 2), 56.d0, tol, "Wrong H[2, 2] in test_trans_mult_matrix!")

  !------------------------------------------------------------------------------------
  ! Test extraction of the matrix row.
  ! TODO: Move to a separate routine.
  call matrix%get_line(1, vec)

  call assert_comparable_real(vec(1), 1.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec(2), 2.d0, tol, "Wrong value in test_trans_mult_matrix!")

  call matrix%get_line(2, vec)

  call assert_comparable_real(vec(1), 3.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec(2), 4.d0, tol, "Wrong value in test_trans_mult_matrix!")

  call matrix%get_line(3, vec)

  call assert_comparable_real(vec(1), 5.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec(2), 6.d0, tol, "Wrong value in test_trans_mult_matrix!")

  !------------------------------------------------------------------------------------
  ! Test extraction of the matrix column.
  ! TODO: Move to a separate routine.
  call matrix%get_column(1, vec2)

  call assert_comparable_real(vec2(1), 1.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec2(2), 3.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec2(3), 5.d0, tol, "Wrong value in test_trans_mult_matrix!")

  call matrix%get_column(2, vec2)

  call assert_comparable_real(vec2(1), 2.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec2(2), 4.d0, tol, "Wrong value in test_trans_mult_matrix!")
  call assert_comparable_real(vec2(3), 6.d0, tol, "Wrong value in test_trans_mult_matrix!")

  deallocate(vec)
  deallocate(vec2)
  deallocate(Avi)
  deallocate(Avj)
  deallocate(H)

end subroutine test_trans_mult_matrix

end module tests_sparse_matrix
