
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

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use string, only: str

  implicit none

  private

  type, public :: t_sparse_matrix
    private

    ! Predicted number of nonzero elements.
    integer(kind=8), private :: nnz
    ! Actual number of nonzero elements.
    integer(kind=8), private :: nel
    ! Last number of elements (when adding a new row).
    integer(kind=8), private :: nel_last
    ! Total number of rows.
    integer, private :: nl
    ! Total actual number of non-empty rows.
    integer, private :: nl_nonempty
    ! The expected number of non-empty rows.
    integer, private :: nl_nonempty_allocated
    ! Current number of the added rows.
    integer, private :: nl_current ! Excludes empty rows.
    integer, private :: nl_current_all
    ! Total number of matrix columns.
    integer, private :: ncolumns

    ! An array of the (left-to-right, then top-to-bottom) non-zero values of the matrix.
    real(kind=MATRIX_PRECISION), allocatable, private :: sa(:)
    ! The column indexes corresponding to the values.
    integer, allocatable, private :: ija(:)
    ! The list of 'sa' indexes where each row starts.
    integer(kind=8), allocatable, private :: ijl(:)

    ! An array to identify non-empty rows.
    integer, allocatable, private :: rowptr(:)

    ! Auxilary array to calculate the solution variance.
    real(kind=CUSTOM_REAL), allocatable, public :: lsqr_var(:)

    ! For passing some info.
    integer, public :: tag

  contains
    private

    procedure, public, pass :: initialize => sparse_matrix_initialize
    procedure, public, pass :: reset => sparse_matrix_reset
    procedure, public, pass :: finalize => sparse_matrix_finalize

    procedure, public, pass :: add => sparse_matrix_add
    procedure, public, pass :: add_row => sparse_matrix_add_row
    procedure, public, pass :: new_row => sparse_matrix_new_row
    procedure, public, pass :: add_empty_rows => sparse_matrix_add_empty_rows

    procedure, public, pass :: mult_vector => sparse_matrix_mult_vector
    procedure, public, pass :: add_mult_vector => sparse_matrix_add_mult_vector
    procedure, public, pass :: part_mult_vector => sparse_matrix_part_mult_vector
    procedure, public, pass :: trans_mult_vector => sparse_matrix_trans_mult_vector
    procedure, public, pass :: add_trans_mult_vector => sparse_matrix_add_trans_mult_vector

    procedure, public, pass :: normalize_columns => sparse_matrix_normalize_columns

    procedure, public, pass :: get_total_row_number => sparse_matrix_get_total_row_number
    procedure, public, pass :: get_current_row_number => sparse_matrix_get_current_row_number
    procedure, public, pass :: get_ncolumns => sparse_matrix_get_ncolumns
    procedure, public, pass :: get_number_elements => sparse_matrix_get_number_elements
    procedure, public, pass :: get_nnz => sparse_matrix_get_nnz

    procedure, private, pass :: validate => sparse_matrix_validate

    procedure, private, pass :: allocate_arrays => sparse_matrix_allocate_arrays

  end type t_sparse_matrix

contains

!=========================================================================
! Initializes the sparse matrix.
!=========================================================================
subroutine sparse_matrix_initialize(this, nl, ncolumns, nnz, myrank, nl_empty)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: nl, ncolumns, myrank
  integer, intent(in), optional :: nl_empty
  integer(kind=8), intent(in) :: nnz

  this%nl_current = 0
  this%nl_current_all = 0
  this%nel = 0
  this%nel_last = 0
  this%nl_nonempty = 0

  this%nl = nl
  this%ncolumns = ncolumns
  this%nnz = nnz

  ! Define the expected number of non-empty rows.
  if (present(nl_empty)) then
    this%nl_nonempty_allocated = this%nl - nl_empty
  else
    this%nl_nonempty_allocated = this%nl
  endif

  call this%allocate_arrays(myrank)

end subroutine sparse_matrix_initialize

!=========================================================================
! Resets the sparse matrix.
!=========================================================================
pure subroutine sparse_matrix_reset(this)
  class(t_sparse_matrix), intent(inout) :: this

  this%nl_current = 0
  this%nl_current_all = 0
  this%nel = 0
  this%nel_last = 0
  this%nl_nonempty = 0

  if (allocated(this%sa)) then
    this%sa = 0._MATRIX_PRECISION
    this%ija = 0
    this%ijl = 0
    this%rowptr = 0
  endif

end subroutine sparse_matrix_reset

!=========================================================================
! (1) Stores the index of last element.
! (2) Validates the matrix indexes.
!=========================================================================
subroutine sparse_matrix_finalize(this, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: myrank

  ! Sanity check.
  if (this%nl_current_all /= this%nl) &
    call exit_MPI("Error in total number of rows in sparse_matrix_finalize!"//new_line('a') &
                  //"nl_current="//str(this%nl_current)//new_line('a') &
                  //"nl="//str(this%nl), myrank, 0)

  ! Sanity check.
  if (this%nel_last /= this%nel) then
    call exit_MPI("Elements were added to the matrix after calling new_row() and before calling finalize()!", myrank, 0)
  endif

  this%ijl(this%nl_current + 1) = this%nel + 1

  this%nl_nonempty = this%nl_current

  if (this%nl_nonempty_allocated > this%nl_nonempty) then
    if (myrank == 0) print *, "nl_nonempty, nl_nonempty_allocated, nl =", &
                              this%nl_nonempty, this%nl_nonempty_allocated, this%nl
  endif

  call this%validate(myrank)

end subroutine sparse_matrix_finalize

!============================================================================
! Validates the boundaries of column indexes.
!============================================================================
subroutine sparse_matrix_validate(this, myrank)
  class(t_sparse_matrix), intent(in) :: this
  integer, intent(in) :: myrank
  integer :: i, j
  integer(kind=8) :: k

  ! Use the same loop as in sparse_matrix_trans_mult_vector().
  do i = 1, this%nl_nonempty
    do k = this%ijl(i), this%ijl(i + 1) - 1
      if (k < 1 .or. k > this%nnz) &
        call exit_MPI("Sparse matrix element-index validation failed!", myrank, 0)

      ! Column index.
      j = this%ija(k)

      if (j < 1 .or. j > this%ncolumns) &
        call exit_MPI("Sparse matrix column-index validation failed!", myrank, j)
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
    call exit_MPI("Error in total number of elements in sparse_matrix_add!", myrank, 0)

  this%nel = this%nel + 1
  this%sa(this%nel) = real(value, MATRIX_PRECISION)
  this%ija(this%nel) = column

end subroutine sparse_matrix_add

!=========================================================================
! Adds one matrix row.
!=========================================================================
subroutine sparse_matrix_add_row(this, nel_add, values, columns, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: nel_add
  real(kind=MATRIX_PRECISION), intent(in) :: values(nel_add)
  integer, intent(in) :: columns(nel_add)
  integer, intent(in) :: myrank

  ! Sanity check.
 if (this%nel + nel_add > this%nnz) &
    call exit_MPI("Error in total number of elements in sparse_matrix_add_row!", myrank, 0)

  this%sa(this%nel + 1 : this%nel + nel_add) = values
  this%ija(this%nel + 1 : this%nel + nel_add) = columns
  this%nel = this%nel + nel_add

end subroutine sparse_matrix_add_row

!=========================================================================
! Adds a new row.
!=========================================================================
subroutine sparse_matrix_new_row(this, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: myrank

  ! Sanity check.
  if (this%nl_current >= this%nl_nonempty_allocated) &
    call exit_MPI("Error in number of rows in sparse_matrix_new_row!"//new_line('a') &
                  //"nl_current="//str(this%nl_current)//new_line('a') &
                  //"nl="//str(this%nl), myrank, 0)

  this%nl_current_all = this%nl_current_all + 1

  ! Store only non-empty rows.
  if (this%nel > this%nel_last) then
    this%nl_current = this%nl_current + 1
    this%ijl(this%nl_current) = this%nel_last + 1

    this%rowptr(this%nl_current) = this%nl_current_all

    this%nel_last = this%nel
  endif

end subroutine sparse_matrix_new_row

!=========================================================================
! Add empty rows.
!=========================================================================
subroutine sparse_matrix_add_empty_rows(this, nrows, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: nrows, myrank
  integer :: i

  if (myrank > 0) continue

  do i = 1, nrows
    ! We exclude empty rows, so just increase the row counter here.
    this%nl_current_all = this%nl_current_all + 1
  enddo

end subroutine sparse_matrix_add_empty_rows

!=========================================================================
! Computes the product between the sparse matrix and vector x: b = Ax.
!=========================================================================
pure subroutine sparse_matrix_mult_vector(this, x, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: x(this%ncolumns)

  real(kind=CUSTOM_REAL), intent(out) :: b(this%nl)

  b = 0._CUSTOM_REAL
  call this%add_mult_vector(x, b)

end subroutine sparse_matrix_mult_vector

!=========================================================================
! Computes the product between the sparse matrix and vector x,
! and adds the result as: b = b + Ax.
!=========================================================================
pure subroutine sparse_matrix_add_mult_vector(this, x, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: x(this%ncolumns)

  real(kind=CUSTOM_REAL), intent(inout) :: b(this%nl)

  integer :: i, i_all
  integer(kind=8) :: k

  do i = 1, this%nl_nonempty
    i_all = this%rowptr(i)
    do k = this%ijl(i), this%ijl(i + 1) - 1
      b(i_all) = b(i_all) + this%sa(k) * x(this%ija(k))
    enddo
  enddo

end subroutine sparse_matrix_add_mult_vector

!===============================================================================================
! Computes the product between the part of sparse matrix and vector x.
! The required part of the matrix is specified via: line_start and param_shift.
!===============================================================================================
subroutine sparse_matrix_part_mult_vector(this, nelements, x, ndata, b, line_start, param_shift, myrank)
  class(t_sparse_matrix), intent(in) :: this
  integer, intent(in) :: nelements, ndata
  real(kind=CUSTOM_REAL), intent(in) :: x(nelements)
  integer, intent(in) :: line_start, param_shift
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: b(ndata)

  integer :: i, l, line_end, i_all
  integer(kind=8) :: k

  line_end = line_start + ndata - 1

  ! Sanity check.
  if (line_start < 1 .or. line_start > this%nl_current_all .or. &
      line_end < 1 .or. line_end > this%nl_current_all) then
    call exit_MPI("Wrong line index in sparse_matrix_part_mult_vector!", myrank, 0)
  endif

  b = 0._CUSTOM_REAL

  do i = 1, this%nl_nonempty
    i_all = this%rowptr(i)
    if (i_all >= line_start .and. i_all <= line_end) then
      l = i_all - line_start + 1
      do k = this%ijl(i), this%ijl(i + 1) - 1
        b(l) = b(l) + this%sa(k) * x(this%ija(k) - param_shift)
      enddo
    endif
  enddo

end subroutine sparse_matrix_part_mult_vector

!============================================================================
! Computes the product between the transpose of sparse matrix and vector x:
! b = A'x.
!============================================================================
pure subroutine sparse_matrix_trans_mult_vector(this, x, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: x(this%nl)

  real(kind=CUSTOM_REAL), intent(out) :: b(this%ncolumns)

  b = 0._CUSTOM_REAL
  call this%add_trans_mult_vector(x, b)

end subroutine sparse_matrix_trans_mult_vector

!============================================================================
! Computes the product between the transpose of sparse matrix and vector x,
! and adds the result as: b = b + A'x.
!============================================================================
pure subroutine sparse_matrix_add_trans_mult_vector(this, x, b)
  class(t_sparse_matrix), intent(in) :: this
  real(kind=CUSTOM_REAL), intent(in) :: x(this%nl)

  real(kind=CUSTOM_REAL), intent(inout) :: b(this%ncolumns)

  integer :: i, j, i_all
  integer(kind=8) :: k

  do i = 1, this%nl_nonempty
    i_all = this%rowptr(i)
    do k = this%ijl(i), this%ijl(i + 1) - 1
      j = this%ija(k)
      b(j) = b(j) + this%sa(k) * x(i_all)
    enddo
  enddo

end subroutine sparse_matrix_add_trans_mult_vector

!================================================================================
! Scale all matrix columns to have unit length and returns the column original norm.
! According to [Paige and Saunders, 1982] this usually removes
! some unnecessary ill-conditioning from the problem.
!
! This function is unit tested in tests_sparse_matrix.f90.
!================================================================================
subroutine sparse_matrix_normalize_columns(this, column_norm)
  class(t_sparse_matrix), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(out) :: column_norm(this%ncolumns)
  integer :: i, j
  integer(kind=8) :: k

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
        this%sa(k) = real(this%sa(k) / column_norm(j), MATRIX_PRECISION)
      endif
    enddo
  enddo

end subroutine sparse_matrix_normalize_columns

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

  res = this%nl_current_all
end function sparse_matrix_get_current_row_number

!=========================================================================
! Returns the total number of columnts in the matrix.
!=========================================================================
pure function sparse_matrix_get_ncolumns(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer :: res

  res = this%ncolumns
end function sparse_matrix_get_ncolumns

!=========================================================================
! Returns the current number of elements stored in the matrix.
!=========================================================================
pure function sparse_matrix_get_number_elements(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer(kind=8) :: res

  res = this%nel
end function sparse_matrix_get_number_elements

!===========================================================================
! Returns the predicted number of elements.
!===========================================================================
pure function sparse_matrix_get_nnz(this) result(res)
  class(t_sparse_matrix), intent(in) :: this
  integer(kind=8) :: res

  res = this%nnz
end function sparse_matrix_get_nnz

!=========================================================================
! Allocates dynamic arrays.
!=========================================================================
subroutine sparse_matrix_allocate_arrays(this, myrank)
  class(t_sparse_matrix), intent(inout) :: this
  integer, intent(in) :: myrank
  integer :: ierr
  real(kind=CUSTOM_REAL) :: mem, mem_loc

  if (this%nnz < 0 .or. this%nl < 0) then
    call exit_MPI("Wrong sizes in sparse_matrix_allocate_arrays!", myrank, 0)
  endif

  mem_loc = kind(this%sa) * this%nnz + &
            kind(this%ijl) * (this%nl_nonempty_allocated + 1) + &
            kind(this%ija) * this%nnz + &
            kind(this%rowptr) * this%nl_nonempty_allocated

  call mpi_allreduce(mem_loc, mem, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (myrank == 0) print *, "Allocating the matrix, memory (GB) =", mem / 1024**3

  allocate(this%sa(this%nnz), source=0._MATRIX_PRECISION, stat=ierr)
  allocate(this%ijl(this%nl_nonempty_allocated + 1), source=int(0, 8), stat=ierr)
  allocate(this%ija(this%nnz), source=0, stat=ierr)
  allocate(this%rowptr(this%nl_nonempty_allocated), source=0, stat=ierr)

  if (ierr /= 0) then
    call exit_MPI("Dynamic memory allocation error in sparse_matrix_allocate_arrays!", myrank, ierr)
  endif

end subroutine sparse_matrix_allocate_arrays

end module sparse_matrix
