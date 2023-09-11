
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
! Unit tests for the t_sparse_matrix class.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module tests_sparse_matrix

  use global_typedefs
  use ftnunit

  use sparse_matrix

  implicit none

  private

  ! Testing matrix columns normalization.
  public :: test_normalize_columns

contains

!=============================================================================================
! Perform test for a (ncolumns x nrows) matrix with some zero columns.
!=============================================================================================
subroutine test_normalize_columns(myrank, nbproc)
  integer, intent(in) :: myrank, nbproc

  type(t_sparse_matrix) :: matrix
  integer :: nrows, ncolumns
  integer :: i, j, counter

  real(kind=CUSTOM_REAL), allocatable :: A(:, :)
  real(kind=CUSTOM_REAL), allocatable :: column_norm(:)
  real(kind=CUSTOM_REAL), allocatable :: A_column(:)
  real(kind=CUSTOM_REAL), allocatable :: vi(:)

  if (nbproc > 0) continue

  ! Set matrix size.
  ncolumns = 10
  nrows = 30

  allocate(A(ncolumns, nrows))
  allocate(column_norm(ncolumns))
  allocate(A_column(nrows))
  allocate(vi(ncolumns))

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

  call matrix%initialize(nrows, ncolumns, int(ncolumns * nrows, 8), myrank)

  ! Copy values from A to the sparse matrix.
  do j = 1, nrows
    do i = 1, ncolumns
      call matrix%add(A(i, j), i, myrank)
    enddo
    call matrix%new_row(myrank)
  enddo

  call matrix%finalize(myrank)

  call matrix%normalize_columns(column_norm)

  do i = 1, ncolumns
    ! Compare matrix norms.
    call assert_comparable_real(column_norm(i), norm2(A(i, :)), tol, "Wrong column_norm array in test_normalize_columns!")

    ! Extract the sparse matrix column.
    vi = 0.d0
    vi(i) = 1.d0
    call matrix%mult_vector(vi, A_column)

    ! Test the norm of the sparse matrix column.
    if (norm2(A(i, :)) /= 0.d0) then
      call assert_comparable_real(norm2(A_column), 1.d0, tol, "Wrong matrix column norm in test_normalize_columns!")
    else
      call assert_comparable_real(norm2(A_column), 0.d0, tol, "Wrong matrix column norm in test_normalize_columns!")
    endif
  enddo

  deallocate(A)
  deallocate(column_norm)
  deallocate(A_column)
  deallocate(vi)

end subroutine test_normalize_columns

end module tests_sparse_matrix
