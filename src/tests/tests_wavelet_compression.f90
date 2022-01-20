
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
! Unit tests for the matrix wavelet compression.
!
! Author: Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module tests_wavelet_compression

  use global_typedefs
  use ftnunit

  use wavelet_transform

  implicit none

  private

  ! Testing data calculation in the wavelet domain.
  public :: test_wavelet_calculate_data

  ! Testing the application of wavelet transform to diagonal matrix.
  public :: test_wavelet_diagonal_matrix

  ! Testing the norm preserving for wavelet transform.
  public :: test_wavelet_norm_preserving
  private :: test_wavelet_norm_preserving_kind

  ! Testing the inverse wavelet transform.
  public :: test_wavelet_inverse
  private :: test_wavelet_inverse_kind

  ! Returns the matrix-vector product.
  private :: matvecmul

contains

!=============================================================================================
! Returns the matrix-vector product: b = A * x.
!=============================================================================================
subroutine matvecmul(A, x, b)
  real(kind=CUSTOM_REAL), intent(in) :: A(:, :)
  real(kind=CUSTOM_REAL), intent(in) :: x(:)
  real(kind=CUSTOM_REAL), intent(out) :: b(:)

  integer :: nrows, j

  nrows = size(b)
  do j = 1, nrows
    b(j) = dot_product(A(:, j), x)
  enddo

end subroutine matvecmul

!=============================================================================================
! Perform test the data calculation in the wavelet domain.
!=============================================================================================
subroutine test_wavelet_calculate_data()

  real(kind=CUSTOM_REAL), allocatable :: A(:, :)
  real(kind=CUSTOM_REAL), allocatable :: x(:)
  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: b2(:)
  integer :: nrows, ncolumns
  integer :: i, j
  integer :: nx, ny, nz

  nx = 3
  ny = 4
  nz = 5

  ! Set matrix size.
  ncolumns = nx * ny * nz
  nrows = 5

  allocate(A(ncolumns, nrows))
  allocate(x(ncolumns))
  allocate(b(nrows))
  allocate(b2(nrows))

  ! Define the matrix.
  do j = 1, nrows
    do i = 1, ncolumns
        A(i, j) = dble(2 * i - j) / dble(i + j)
    enddo
  enddo

  ! Define the vector x.
  do i = 1, ncolumns
    x(i) = dble(2 * j + 1)
  enddo

  ! Calculate b = A * x.
  call matvecmul(A, x, b)

  print *, "b (normal) = ", b

  ! Wavelet transform the matrix rows, transforming matrix to the wavelet domain: A --> A_w
  do j = 1, nrows
    call Haar3D(A(:, j), nx, ny, nz)
  enddo

  ! Calculate A * x in the wavelet domain: b2 = A_w * x_w.
  ! To perform this, we first transform x to the wavelet domain: x --> x_w
  call Haar3D(x, nx, ny, nz)
  call matvecmul(A, x, b2)

  print *, "b (wavelet) = ", b2

  do j = 1, nrows
    call assert_comparable_real(b(j), b2(j), tol, "Wrong result in test_wavelet_calculate_data!")
  enddo

  deallocate(A)
  deallocate(x)
  deallocate(b)
  deallocate(b2)

end subroutine test_wavelet_calculate_data

!=============================================================================================
! Testing the application of wavelet transform to diagonal matrix.
!=============================================================================================
subroutine test_wavelet_diagonal_matrix()

  real(kind=CUSTOM_REAL), allocatable :: A(:, :)
  integer :: nrows, ncolumns
  integer :: i, j
  integer :: nx, ny, nz, nnz

  nx = 10
  ny = 10
  nz = 10

  ! Set matrix size.
  ncolumns = nx * ny * nz

  ! A square matrix.
  nrows = ncolumns

  allocate(A(ncolumns, nrows))

  ! Define the diagonal matrix.
  A = 0.d0
  do i = 1, ncolumns
    A(i, i) = 1.d0
  enddo

  ! Wavelet transform the matrix rows, transforming matrix to the wavelet domain: A --> A_w
  do j = 1, nrows
    call Haar3D(A(:, j), nx, ny, nz)
  enddo

  ! The number of non-zero elements in the matrix.
  nnz = count(A /= 0.d0)

  print *, 'nnz =', nnz

  ! The test result was taken from execution of this test.
  call assert_equal_int(nnz, 46656, "nnz /= 46656 in test_wavelet_diagonal_matrix!")

  deallocate(A)
end subroutine test_wavelet_diagonal_matrix

!=============================================================================================
! Testing that the wavelet transform is L2 norm preserving.
!=============================================================================================
subroutine test_wavelet_norm_preserving()

  ! Haar wavelet.
  call test_wavelet_norm_preserving_kind(1)

  ! TODO: DaubD43D wavelet not passing the test! Should it preserve the norm?
  ! Daubechies D4 wavelet.
  !call test_wavelet_norm_preserving_kind(2)

end subroutine test_wavelet_norm_preserving
!=============================================================================================

subroutine test_wavelet_norm_preserving_kind(waveletType)
  integer, intent(in) :: waveletType

  real(kind=CUSTOM_REAL), allocatable :: x(:)
  integer :: N, i
  integer :: nx, ny, nz
  real(kind=CUSTOM_REAL) :: norm, norm_w

  nx = 10
  ny = 11
  nz = 12

  N = nx * ny * nz

  allocate(x(N))

  ! Forming the image to compress.
  do i = 1, N
    x(i) = dble(i)
  enddo

  norm = norm2(x)
  print *, 'norm =', norm

  if (waveletType == 1) then
    call Haar3D(x, nx, ny, nz)
  else
    call DaubD43D(x, nx, ny, nz)
  endif

  norm_w = norm2(x)
  print *, 'norm_w =', norm_w

  call assert_comparable_real(norm, norm_w, tol, "Wrong result in test_wavelet_norm_preserving_kind!")

  deallocate(x)

end subroutine test_wavelet_norm_preserving_kind

!=============================================================================================
! Testing the inverse wavelet transform.
!=============================================================================================
subroutine test_wavelet_inverse()

  ! Haar wavelet.
  call test_wavelet_inverse_kind(1)

  ! Daubechies D4 wavelet.
  call test_wavelet_inverse_kind(2)

end subroutine test_wavelet_inverse

!=============================================================================================
! Testing the inverse wavelet transform.
!=============================================================================================
subroutine test_wavelet_inverse_kind(waveletType)
  integer, intent(in) :: waveletType

  real(kind=CUSTOM_REAL), allocatable :: A(:, :)
  integer :: nrows, ncolumns
  integer :: i, j
  integer :: nx, ny, nz, nnz

  nx = 10
  ny = 11
  nz = 12

  ! Set matrix size.
  ncolumns = nx * ny * nz

  ! A square matrix.
  nrows = ncolumns

  allocate(A(ncolumns, nrows))

  ! Define the identity matrix.
  A = 0.d0
  do i = 1, ncolumns
    A(i, i) = 1.d0
  enddo

  print *, 'nnz (orig) = ', count(A > 1.d-15)

  ! Wavelet transform the matrix rows.
  do j = 1, nrows
    if (waveletType == 1) then
      call Haar3D(A(:, j), nx, ny, nz)
    else
      call DaubD43D(A(:, j), nx, ny, nz)
    endif
  enddo

  print *, 'nnz (wavelet) = ', count(A > 1.d-15)

  ! Applying the inverse wavelet transform.
  do j = 1, nrows
    if (waveletType == 1) then
      call iHaar3D(A(:, j), nx, ny, nz)
    else
      call iDaubD43D(A(:, j), nx, ny, nz)
    endif
  enddo

  nnz = count(A > 1.d-15)
  print *, 'nnz (inverse) = ', nnz

  call assert_equal_int(nnz, ncolumns, "Wrong nnz in test_wavelet_inverse_kind!")

  ! Testing that the inverse wavelet transform leads to original identity matrix.
  do i = 1, ncolumns
    do j = 1, nrows
      if (i == j) then
        call assert_comparable_real(A(i, j), 1.d0, tol, "Wrong result in test_wavelet_inverse_kind!")
      else
        call assert_true(abs(A(i, j)) < 1.d-15, "Wrong result in test_wavelet_inverse_kind!")
      endif
    enddo
  enddo

  deallocate(A)
end subroutine test_wavelet_inverse_kind

end module tests_wavelet_compression
