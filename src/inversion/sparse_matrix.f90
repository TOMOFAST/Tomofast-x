
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

!========================================================================================
! A class to work with sparse matrices that are
! stored using Compressed Sparse Row (CSR) format.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!========================================================================================
module sparse_matrix

  use global_typedefs, only: CUSTOM_REAL
  use mpi_tools, only: exit_MPI
  use string, only: str

  implicit none

  private

  type, public :: t_sparse_matrix
    !private

    ! Predicted number of nonzero elements.
    integer, private :: nnz
    ! Actual number of nonzero elements.
    integer, private :: nel
    ! Total number of rows in the matrix.
    integer, private :: nl
    ! Current number of the added lines (rows) in the matrix.
    integer, private :: nl_current

    ! An array of the (left-to-right, then top-to-bottom) non-zero values of the matrix.
    real(kind=CUSTOM_REAL), allocatable, private :: sa(:)
    ! The column indexes corresponding to the values.
    integer, allocatable, private :: ija(:)
    ! The list of 'sa' indexes where each row starts.
    integer, allocatable, private :: ijl(:)

    ! Auxilary array to calculate the solution variance (var = diag(D.D')).
    ! [1] E. Kostina, M. A. Saunders, and I. Schierle, Computation of covariance matrices 
    !     for constrained parameter estimation problems using lsqr, 2009.
    real(kind=CUSTOM_REAL), allocatable, public :: lsqr_var(:)
    ! For passing some info.
    integer, public :: tag

  contains
    private

    procedure, public, pass :: initialize => sparse_matrix_initialize
    procedure, public, pass :: reset => sparse_matrix_reset
    procedure, public, pass :: finalize => sparse_matrix_finalize

    procedure, public, pass :: add => sparse_matrix_add
    procedure, public, pass :: new_row => sparse_matrix_new_row
    procedure, public, pass :: add_empty_rows => sparse_matrix_add_empty_rows
    procedure, public, pass :: add_matrix => sparse_matrix_add_matrix

    procedure, public, pass :: mult_vector => sparse_matrix_mult_vector
    procedure, public, pass :: trans_mult_vector => sparse_matrix_trans_mult_vector
    procedure, public, pass :: trans_mult_matrix => sparse_matrix_trans_mult_matrix

    procedure, public, pass :: normalize_columns => sparse_matrix_normalize_columns

    procedure, public, pass :: get_total_row_number => sparse_matrix_get_total_row_number
    procedure, public, pass :: get_current_row_number => sparse_matrix_get_current_row_number
    procedure, public, pass :: get_number_elements => sparse_matrix_get_number_elements
    procedure, public, pass :: get_nnz => sparse_matrix_get_nnz
    procedure, public, pass :: get_value => sparse_matrix_get_value
    procedure, public, pass :: get_line => sparse_matrix_get_line
    procedure, public, pass :: get_column => sparse_matrix_get_column
    procedure, public, pass :: get_integrated_sensit => sparse_matrix_get_integrated_sensit

    procedure, public, pass :: allocate_variance_array => sparse_matrix_allocate_variance_array

    procedure, private, pass :: validate => sparse_matrix_validate

    procedure, private, pass :: allocate_arrays => sparse_matrix_allocate_arrays

  end type t_sparse_matrix

contains

!=========================================================================
! Returns the total number of rows in the matrix.
!=========================================================================
pure function sparse_matrix_get_total_row_number(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer :: res

  res = this%nl
end function sparse_matrix_get_total_row_number

!=========================================================================
! Returns the current number of rows in the matrix.
!=========================================================================
pure function sparse_matrix_get_current_row_number(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer :: res

  res = this%nl_current
end function sparse_matrix_get_current_row_number

!=========================================================================
! Returns the current number of elements stored in the matrix.
!=========================================================================
pure function sparse_matrix_get_number_elements(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer :: res

  res = this%nel
end function sparse_matrix_get_number_elements

!===========================================================================
! Returns the predicted number of elements.
!===========================================================================
pure function sparse_matrix_get_nnz(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer :: res

  res = this%nnz
end function sparse_matrix_get_nnz

!===========================================================================
! Returns the (i, j)-element's value of the matrix.
! i - column number, j - row number.
!===========================================================================
pure function sparse_matrix_get_value(this, i, j) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer, intent(in) :: i, j
  real(kind=CUSTOM_REAL) :: res
  integer :: k

  res = 0.d0
  do k = this%ijl(j), this%ijl(j + 1) - 1
    if (this%ija(k) == i) then
    ! Found a non-zero element at the column i.
      res = this%sa(k)
      exit
    endif
  enddo

end function sparse_matrix_get_value

!=========================================================================
! Initializes the sparse matrix.
!=========================================================================
subroutine sparse_matrix_initialize(this, nl, nnz, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: nl, nnz, myrank

  this%nl_current = 0
  this%nel = 0
  this%nl = nl
  this%nnz = nnz

  call this%allocate_arrays(myrank)

end subroutine sparse_matrix_initialize

!=========================================================================
! Resets the sparse matrix.
!=========================================================================
pure subroutine sparse_matrix_reset(this)
  class(t_sparse_matrix), intent(inout) :: this

  this%nl_current = 0
  this%nel = 0

  this%sa = 0._CUSTOM_REAL
  this%ija = 0
  this%ijl = 0

end subroutine sparse_matrix_reset

!=========================================================================
! (1) Stores the index of last element.
! (2) Validates the matrix indexes.
!=========================================================================
subroutine sparse_matrix_finalize(this, ncolumns, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: ncolumns, myrank

  ! Sanity check.
  if (this%nl_current /= this%nl) &
    call exit_MPI("Error in total number of rows in sparse_matrix_finalize!"//new_line('a') &
                  //"nl_current="//str(this%nl_current)//new_line('a') &
                  //"nl="//str(this%nl), myrank, 0)

  this%ijl(this%nl + 1) = this%nel + 1

  if (ncolumns > 0) call this%validate(ncolumns, myrank)

end subroutine sparse_matrix_finalize

!============================================================================
! Validates the boundaries of column indexes.
!============================================================================
subroutine sparse_matrix_validate(this, ncolumns, myrank)
  class(t_sparse_matrix), intent(in) :: this
  integer, intent(in) :: ncolumns, myrank
  integer :: i, j, k

  ! Use the same loop as in sparse_matrix_trans_mult_vector().
  do i = 1, this%nl
    do k = this%ijl(i), this%ijl(i + 1) - 1
      if (k < 1 .or. k > this%nnz) &
        call exit_MPI("Sparse matrix validation failed (k)!", myrank, k)

      j = this%ija(k)

      if (j < 1 .or. j > ncolumns) &
        call exit_MPI("Sparse matrix validation failed (j)!", myrank, j)
    enddo
  enddo

end subroutine sparse_matrix_validate

!=========================================================================
! Adds one element at specified position (column index).
!=========================================================================
subroutine sparse_matrix_add(this, value, column, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: value
  integer, intent(in) :: column, myrank

  ! Do not add zero values to a sparse matrix.
  if (value == 0._CUSTOM_REAL) return

  ! Sanity check.
  if (this%nel >= this%nnz) &
    call exit_MPI("Error in total number of elements in sparse_matrix_add!", myrank, this%nnz)

  this%nel = this%nel + 1
  this%sa(this%nel) = value
  this%ija(this%nel) = column

end subroutine sparse_matrix_add

!=========================================================================
! Adds a new row.
!=========================================================================
subroutine sparse_matrix_new_row(this, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: myrank

  ! Sanity check.
  if (this%nl_current >= this%nl) &
    call exit_MPI("Error in number of rows in sparse_matrix_new_row!"//new_line('a') &
                  //"nl_current="//str(this%nl_current)//new_line('a') &
                  //"nl="//str(this%nl), myrank, 0)

  this%nl_current = this%nl_current + 1
  this%ijl(this%nl_current) = this%nel + 1

end subroutine sparse_matrix_new_row

!=========================================================================
! Add empty rows.
!=========================================================================
subroutine sparse_matrix_add_empty_rows(this, nrows, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: nrows, myrank
  integer :: i

  do i = 1, nrows
    call this%new_row(myrank)
  enddo
end subroutine sparse_matrix_add_empty_rows

!=================================================================================
! Add a matrix B in the current matrix A below, so that a resulting matrix C is
!
! C = (    A )
!     ( mu B ), where mu is a weighting factor.
!=================================================================================
subroutine sparse_matrix_add_matrix(this, matrix_B, mu, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  type(t_sparse_matrix), intent(in) :: matrix_B
  real(kind=CUSTOM_REAL), intent(in) :: mu
  integer, intent(in) :: myrank
  integer :: i, k

  do i = 1, matrix_B%nl
    call this%new_row(myrank)

    do k = matrix_B%ijl(i), matrix_B%ijl(i + 1) - 1
      call this%add(mu * matrix_B%sa(k), matrix_B%ija(k), myrank)
    enddo
  enddo

end subroutine sparse_matrix_add_matrix

!=========================================================================
! Computes the product between the sparse matrix and vector x.
!=========================================================================
pure subroutine sparse_matrix_mult_vector(this, x, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: x(:)
  real(kind=CUSTOM_REAL), intent(out) :: b(:)
  integer :: i, k

  b = 0._CUSTOM_REAL

  do i = 1, this%nl
    do k = this%ijl(i), this%ijl(i + 1) - 1
      b(i) = b(i) + this%sa(k) * x(this%ija(k))
    enddo
  enddo

end subroutine sparse_matrix_mult_vector

!=========================================================================
! Extract the j-th line from the matrix and stores in b-vector.
!=========================================================================
pure subroutine sparse_matrix_get_line(this, j, b)
  class(t_sparse_matrix), intent(in) :: this
  integer, intent(in) :: j
  real(kind=CUSTOM_REAL), intent(out) :: b(:)
  integer :: k

  b = 0._CUSTOM_REAL

  do k = this%ijl(j), this%ijl(j + 1) - 1
    b(this%ija(k)) = this%sa(k)
  enddo

end subroutine sparse_matrix_get_line

!============================================================================
! Computes the product between the transpose of sparse matrix and vector x.
!============================================================================
pure subroutine sparse_matrix_trans_mult_vector(this, x, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: x(:)
  real(kind=CUSTOM_REAL), intent(out) :: b(:)
  integer :: i, j, k

  b = 0._CUSTOM_REAL

  do i = 1, this%nl
!IBM* ASSERT (NODEPS,ITERCNT(1000))
!DIR$ IVDEP
    do k = this%ijl(i), this%ijl(i + 1) - 1
      j = this%ija(k)
      b(j) = b(j) + this%sa(k) * x(i)
    enddo
  enddo

end subroutine sparse_matrix_trans_mult_vector

!============================================================================
! Extract a column from the matrix and stores in b-vector.
!============================================================================
pure subroutine sparse_matrix_get_column(this, column, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(out) :: b(:)
  integer, intent(in) :: column
  integer :: i, j, k

  b = 0._CUSTOM_REAL

  do i = 1, this%nl
!IBM* ASSERT (NODEPS,ITERCNT(1000))
!DIR$ IVDEP
    do k = this%ijl(i), this%ijl(i + 1) - 1
      j = this%ija(k)
      if (j == column) then
        b(i) = this%sa(k)
      endif
    enddo
  enddo

end subroutine sparse_matrix_get_column

!============================================================================
! Calculates the integrated sensitivity (norm of columns).
! The output dimension = the number of columns.
!============================================================================
pure subroutine sparse_matrix_get_integrated_sensit(this, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(out) :: b(:)
  integer :: i, j, k

  b = 0._CUSTOM_REAL

  do i = 1, this%nl
!IBM* ASSERT (NODEPS,ITERCNT(1000))
!DIR$ IVDEP
    do k = this%ijl(i), this%ijl(i + 1) - 1
      ! Column number.
      j = this%ija(k)
      ! The sum of squared column elements.
      b(j) = b(j) + this%sa(k) * this%sa(k)
    enddo
  enddo
  b = sqrt(b)

end subroutine sparse_matrix_get_integrated_sensit

!================================================================================
! Computes the (i,j)-element of the product between
! the transpose of (this) matrix and (this) matrix, i.e, H = A'A, return H[i,j].
! SLOW VERSION (use auxiliary arrays), but simple to implement.
! vec has size on ncolumns, and Avi, Avj have size of nrows.
!
! This function is unit tested in tests_sparse_matrix.f90.
!================================================================================
function sparse_matrix_trans_mult_matrix(this, i, j, vec, Avi, Avj) result (Hij)
  class(t_sparse_matrix), intent(in) :: this
  integer, intent(in) :: i, j
  ! Use for storage only, to do not allocate memory here.
  real(kind=CUSTOM_REAL), intent(inout) :: vec(:)
  real(kind=CUSTOM_REAL), intent(inout) :: Avi(:)
  real(kind=CUSTOM_REAL), intent(inout) :: Avj(:)
  integer :: k

  real(kind=CUSTOM_REAL) :: Hij

  vec = 0.d0
  vec(i) = 1.d0

  call this%mult_vector(vec, Avi)

  vec(i) = 0.d0
  vec(j) = 1.d0

  call this%mult_vector(vec, Avj)

  Hij = 0.d0

  do k = 1, this%nl
    Hij = Hij + Avi(k) * Avj(k)
  enddo

end function sparse_matrix_trans_mult_matrix


!================================================================================
! Scale all matrix columns to have unit length and returns the column original norm.
! According to [Paige and Saunders, 1982] this usually removes
! some unnecessary ill-conditioning from the problem.
!
! This function is unit tested in tests_sparse_matrix.f90.
!================================================================================
subroutine sparse_matrix_normalize_columns(this, column_norm)
  class(t_sparse_matrix), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(out) :: column_norm(:)
  integer :: i, j, k

  column_norm = 0._CUSTOM_REAL

  ! Calculate the column norm.
  do i = 1, this%nl
    do k = this%ijl(i), this%ijl(i + 1) - 1
      j = this%ija(k)
      column_norm(j) = column_norm(j) + this%sa(k)**2
    enddo
  enddo

  column_norm = sqrt(column_norm)

  ! Normalize matrix columns.
  do i = 1, this%nl
    do k = this%ijl(i), this%ijl(i + 1) - 1
      j = this%ija(k)

      if (column_norm(j) /= 0.d0) then
        this%sa(k) = this%sa(k) / column_norm(j)
      endif
    enddo
  enddo

end subroutine sparse_matrix_normalize_columns

!=========================================================================
! Allocates dynamic arrays.
!=========================================================================
subroutine sparse_matrix_allocate_arrays(this, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: myrank
  integer :: ierr

  if (this%nnz <= 0 .or. this%nl <= 0) &
    call exit_MPI("Wrong sizes in sparse_matrix_allocate_arrays!", myrank, 0)

  ierr = 0

  if (.not. allocated(this%sa)) allocate(this%sa(this%nnz), source=0._CUSTOM_REAL, stat=ierr)
  if (.not. allocated(this%ijl)) allocate(this%ijl(this%nl + 1), source=0, stat=ierr)
  if (.not. allocated(this%ija)) allocate(this%ija(this%nnz), source=0, stat=ierr)

  if (ierr /= 0) &
    call exit_MPI("Dynamic memory allocation error in sparse_matrix_allocate_arrays!", myrank, ierr)

end subroutine sparse_matrix_allocate_arrays

!=========================================================================
! Allocates dynamic array to store solution variance.
!=========================================================================
subroutine sparse_matrix_allocate_variance_array(this, isize, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: isize, myrank
  integer :: ierr

  ierr = 0

  if (.not. allocated(this%lsqr_var)) allocate(this%lsqr_var(isize), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) &
    call exit_MPI("Dynamic memory allocation error in sparse_matrix_allocate_variance_array!", myrank, ierr)

end subroutine sparse_matrix_allocate_variance_array

end module sparse_matrix
