
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

!=============================================================================
! Sequential Coordinate-wise Algorithm for Non-negative Lease Squares (NNLS) Problem.
!
! V. Franc, V. Hlavac, and M. Navara,
! Sequential coordinate-wise algorithm for nonnegative least squares problem,
! Computer Analysis of Images and Patterns,
! Volume 3691 of the series Lecture Notes in Computer Science pp. 407-414.
! http://cmp.felk.cvut.cz/ftp/articles/franc/Franc-Hlavac-Navara-CAIP05.pdf
!
! Vitaliy Ogarko, UWA, CET, Australia, 2016.
!=============================================================================
module sca_solver

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix

  implicit none

  private

  public :: sca_solve

contains

!================================================================================
! NNLS solver for a sparse matrix.
! Solve min||Ax - b||, where A is stored in sparse format in 'matrix' variable.
!
! Niter - maximum number of iterations.
! rmin - stopping criterion (relative residual).
!
! NOTE: it is not parallelized!
!
!================================================================================
subroutine sca_solve(niter, rmin, matrix, b, x, myrank, nbproc)
  integer, intent(in) :: myrank, nbproc
  integer, intent(in) :: niter
  real(kind=CUSTOM_REAL), intent(in) :: rmin
  real(kind=CUSTOM_REAL), intent(in) :: b(:)
  type(t_sparse_matrix), intent(in) :: matrix

  real(kind=CUSTOM_REAL), intent(out) :: x(:)

  ! Local variables.
  integer :: iter, k, j, ierr
  integer :: nrows, nelements
  real(kind=CUSTOM_REAL) :: r, xk, b_norm
  real(kind=CUSTOM_REAL), allocatable :: mu(:)
  !real(kind=CUSTOM_REAL), allocatable :: h(:)
  real(kind=CUSTOM_REAL), allocatable :: H(:, :)
  ! Auxiliary arrays.
  real(kind=CUSTOM_REAL), allocatable :: vec(:)
  real(kind=CUSTOM_REAL), allocatable :: Avi(:)
  real(kind=CUSTOM_REAL), allocatable :: Avj(:)

  if (myrank == 0) print *, 'Entered subroutine sca_solve.'

  if (nbproc > 1) then
    if (myrank == 0) print *, "WARNING: SCA solver is not yet adjusted for parallel runs! Exiting."
    return
  endif

  nrows = matrix%get_total_row_number()
  nelements = size(x)

  print *, 'nrows, nelements =', nrows, nelements

  ! Allocate memory.
  allocate(mu(nelements), source=0._CUSTOM_REAL, stat=ierr)
  !allocate(h(nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(vec(nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(Avi(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(Avj(nrows), source=0._CUSTOM_REAL, stat=ierr)
  allocate(H(nelements, nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in sca_solve!", myrank, ierr)

  ! Sanity check.
  if (norm2(b) == 0.d0) then
    call exit_MPI("|b| = 0, prior model is exact? Exiting.", myrank, 0)
  end if

  ! Initialization.
  call matrix%trans_mult_vector(b, mu)
  mu = - mu
  x = 0.d0

  if (myrank == 0) print *, "Calculating the matrix H = A'A ..."
  do k = 1, nelements
    do j = 1, nelements
      H(k, j) = matrix%trans_mult_matrix(k, j, vec, Avi, Avj)
    enddo
  enddo
  if (myrank == 0) print *, "Finished calculating the matrix H = A'A."

  iter = 1
  r = 1._CUSTOM_REAL
  b_norm = norm2(b)

  ! Use an exit (threshold) criterion in case we can exit the loop before reaching the max iteration count.
  do while (iter <= niter .and. r > rmin)

    do k = 1, nelements
      ! Calculate h - the kth column of the matrix H = A'A.
      ! (!) Note: H is a huge matrix and in general cannot be pre-calculated due to memory consumption.
      ! (!) Note: this below loop is too slow as we use sparse row format,
      !           while the matrix columns are needed for h(j) calculation.
      !           Maybe using Compressed Sparse Column (CSC) format can speed-up this code substantially.
!      do j = 1, nelements
!        h(j) = matrix%trans_mult_matrix(k, j, vec, Avi, Avj)
!      enddo

      xk = x(k)
      !x(k) = max(0.d0, x(k) - mu(k) / h(k))

      if (H(k, k) /= 0.d0) then
        x(k) = max(0.d0, x(k) - mu(k) / H(k, k))
      else
        x(k) = 0.d0
      endif

      !mu = mu + (x(k) - xk) * h
      mu = mu + (x(k) - xk) * H(k, :)

      !if (myrank == 0) print *, 'k =', k
    enddo

    ! Calculate the residual. Use Avi just for storage.
    call matrix%mult_vector(x, Avi)

    ! Norm of the relative residual.
    r = norm2(Avi - b) / b_norm

    if (mod(iter, 10) == 0) then
      if (myrank == 0) print *, 'it, r =', iter, r
    endif

    iter = iter + 1
  enddo

  if (myrank == 0) print *, 'End of subroutine sca_solve, r =', r, ' iter =', iter - 1

#ifdef USE_FLUSH6
  call flush(6)
#endif

  deallocate(mu)
  !deallocate(h)
  deallocate(H)
  deallocate(vec)
  deallocate(Avi)
  deallocate(Avj)

end subroutine sca_solve

end module sca_solver
