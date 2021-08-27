
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
  use sensitivity_gravmag

  implicit none

  private

  ! Testing data calculation in the compressed domain.
  public :: test_wavelet_calculate_data
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
! Perform test the data calculation in the compressed variables.
!=============================================================================================
subroutine test_wavelet_calculate_data(myrank)
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), allocatable :: A(:, :)
  real(kind=CUSTOM_REAL), allocatable :: x(:)
  real(kind=CUSTOM_REAL), allocatable :: b(:)
  real(kind=CUSTOM_REAL), allocatable :: b2(:)
  integer :: nrows, ncolumns
  integer :: i, j, counter
  integer :: nx, ny, nz
  real(kind=CUSTOM_REAL) :: threshold, comp_rate
  type(t_sensitivity_gravmag) :: sens

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
  threshold = 0.d0
  do j = 1, nrows
    call sens%compress_matrix_line_wavelet(nx, ny, nz, A(:, j), threshold, comp_rate)
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

end module tests_wavelet_compression
