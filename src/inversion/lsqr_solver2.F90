
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
! Least Square (LSQR) solver parallelized by data.
!
! This solver is unit tested in tests_lsqr.f90 and in tests_method_of_weights.f90.
! The unit tests are available for serial and parallel versions.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!=============================================================================
module lsqr_solver

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix

  implicit none

  private

  public :: lsqr_solve

  private :: normalize
  private :: apply_soft_thresholding

contains

!================================================================================
! LSQR solver for a sparse matrix.
! Solve min||Ax - b||, where A is stored in sparse format in 'matrix' variable.
!
! nlines - local number of matrix rows.
! nelements - number of parameters.
! Niter - maximum number of iterations.
! rmin - stopping criterion (relative residual).
! gamma - soft thresholding parameter (see ISTA, proximate L1 norm). Use gamma=0 for pure L2 norm.
!================================================================================
subroutine lsqr_solve(nlines, nelements, niter, rmin, gamma, matrix, b, x, myrank)
  integer, intent(in) :: nlines, nelements, niter
  real(kind=CUSTOM_REAL), intent(in) :: rmin, gamma
  real(kind=CUSTOM_REAL), intent(in) :: b(nlines)
  integer, intent(in) :: myrank

  type(t_sparse_matrix), intent(inout) :: matrix
  real(kind=CUSTOM_REAL), intent(inout) :: x(nelements)

  ! Local variables.
  integer :: iter, ierr
  real(kind=CUSTOM_REAL) :: alpha, beta, rho, rhobar, phi, phibar, theta
  real(kind=CUSTOM_REAL) :: b1, c, r, s, t1, t2
  real(kind=CUSTOM_REAL) :: rho_inv
  ! Flag for calculating variance.
  logical :: calculateVariance

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: v, w, u

  if (myrank == 0) print *, 'Entered subroutine lsqr_solve, gamma =', gamma

  ! Sanity check.
  if (matrix%get_total_row_number() /= nlines) then
    call exit_MPI("Wrong matrix size in lsqr_solve! Exiting.", myrank, 0)
  endif

  if (allocated(matrix%lsqr_var) .and. size(matrix%lsqr_var) == nelements) then
    calculateVariance = .true.

    ! Initialize variance array.
    matrix%lsqr_var = 0.d0
  else
    calculateVariance = .false.
  endif

  ! Allocate memory.
  allocate(u(nlines))
  allocate(v(nelements))
  allocate(w(nelements))

  ! Sanity check.
  if (norm2(b) == 0.d0) then
    call exit_MPI("|b| = 0, starting model is exact? Exiting.", myrank, 0)
  end if

  ! Initialization.
  u = b

  ! Normalize u and initialize beta.
  call normalize(nlines, u, beta, .true., ierr)
  if (ierr /= 0) then
    call exit_MPI("Could not normalize initial u, zero denominator!", myrank, 0)
  endif

  b1 = beta

  ! Required by the algorithm.
  x = 0._CUSTOM_REAL

  ! Compute v = Ht.u.
  call matrix%trans_mult_vector(u, v)

  call MPI_Allreduce(MPI_IN_PLACE, v, nelements, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Normalize v and initialize alpha.
  call normalize(nelements, v, alpha, .false., ierr)
  if (ierr /= 0) then
    call exit_MPI("Could not normalize initial v, zero denominator!", myrank, 0)
  endif

  rhobar = alpha
  phibar = beta
  w = v

  iter = 1
  r = 1._CUSTOM_REAL

  ! Use an exit (threshold) criterion in case we can exit the loop before reaching the max iteration count.
  do while (iter <= niter .and. r > rmin)

    ! Scale u.
    u = - alpha * u

    ! Compute u = u + H.v.
    call matrix%add_mult_vector(v, u)

    ! Normalize u and update beta.
    call normalize(nlines, u, beta, .true., ierr)
    if (ierr /= 0) then
      ! Found an exact solution. Happens in the unit test test_lsqr_underdetermined_2.
      if (myrank == 0) print *, 'WARNING: u = 0. Possibly found an exact solution in the LSQR solver!'
    endif

    !-----------------------------------------------------------------
    ! Compute v = - beta.v + Ht.u in parallel.
    !-----------------------------------------------------------------
    if (myrank == 0) then
      v = - beta * v
    else
      v = 0._CUSTOM_REAL
    endif

    ! Compute v = v + Ht_loc.u.
    call matrix%add_trans_mult_vector(u, v)

    ! Sum partial results from all ranks: v = (- beta.v + Htu_loc1) + Htu_loc2 + ... + Htu_locN = - beta.v + Ht.u.
    call MPI_Allreduce(MPI_IN_PLACE, v, nelements, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    !---------------------------------------------------------------

    ! Normalize v and update alpha.
    call normalize(nelements, v, alpha, .false., ierr)
    if (ierr /= 0) then
      ! Found an exact solution. Happens in the unit test test_method_of_weights_1.
      if (myrank == 0) print *, 'WARNING: v = 0. Possibly found an exact solution in the LSQR solver!'
    endif

    ! Compute scalars for updating the solution.
    rho = sqrt(rhobar * rhobar + beta * beta)

    ! Sanity check (avoid zero division).
    if (rho == 0._CUSTOM_REAL) then
      print *, 'WARNING: rho = 0. Exiting.'
      exit
    endif

    ! Compute scalars for updating the solution.
    rho_inv = 1.d0 / rho

    c       = rhobar * rho_inv
    s       = beta * rho_inv
    theta   = s * alpha
    rhobar  = - c * alpha
    phi     = c * phibar
    phibar  = s * phibar
    t1      = phi * rho_inv
    t2      = - theta * rho_inv

    if (calculateVariance) then
      ! Calculate solution variance (see Paige & Saunders 1982, page 53).
      matrix%lsqr_var = matrix%lsqr_var + (rho_inv * w)**2
    endif

    ! Update the current solution x (w is an auxiliary array in order to compute the solution).
    x = t1 * w + x
    w = t2 * w + v

    if (gamma /= 0._CUSTOM_REAL) then
    ! Soft thresholding.
      call apply_soft_thresholding(x, nelements, gamma)

      ! Calculate the residual for thresholded solution.
      !call matrix%mult_vector(x, Hv)
      !call MPI_Allreduce(MPI_IN_PLACE, Hv, nlines, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Norm of the relative residual.
      !r = norm2(Hv - b) / b1

      ! Approximate residual.
      r = phibar / b1
    else
      ! Norm of the relative residual.
      r = phibar / b1
    endif

#ifndef SUPPRESS_OUTPUT
    ! Commented as this badly affects the performance.
    if (mod(iter, 10) == 0) then
      if (myrank == 0) print *, 'it, r =', iter, r
    endif
#endif

    ! To avoid floating point exception of denormal value.
    if (abs(rhobar) < 1.e-30) then
      if (myrank == 0) print *, 'WARNING: Small rhobar! Possibly algorithm has converged, exiting the loop.'
      exit
    endif

    iter = iter + 1
  enddo

  if (calculateVariance) then
    ! Scale with data residual (see the bottom of the page 53 in Paige & Saunders 1982).
    ! Note, we are using that phibar = ||Ax - b||
    if (matrix%tag > 1) then
      matrix%lsqr_var = phibar**2 * matrix%lsqr_var
    else
    ! First major iteration (stored in matrix%tag). Scale with original residual to calculate the prior variance.
      matrix%lsqr_var = b1**2 * matrix%lsqr_var
    endif
    ! Calculate the standard error (s_i).
    matrix%lsqr_var = sqrt(matrix%lsqr_var)
  endif

  if (myrank == 0) print *, 'End of subroutine lsqr_solve, r =', r, ' iter =', iter - 1

  deallocate(u)
  deallocate(v)
  deallocate(w)

end subroutine lsqr_solve

!==========================================================================
! Applies soft thresholding (ISTA, proximate L1 norm).
!==========================================================================
pure subroutine apply_soft_thresholding(x, nelements, gamma)
  integer, intent(in) :: nelements
  real(kind=CUSTOM_REAL), intent(in) :: gamma
  real(kind=CUSTOM_REAL), intent(inout) :: x(nelements)

  integer :: i

  do i = 1, nelements
    if (abs(x(i)) <= gamma) then
      x(i) = 0._CUSTOM_REAL
    else if (x(i) <= - gamma) then
      x(i) = x(i) + gamma
    else if (x(i) >= gamma) then
      x(i) = x(i) - gamma
    endif
  enddo
end subroutine apply_soft_thresholding

!============================================================================
! Normalizes vector x by its L2 norm.
! Set in_parallel=true, if x is split between CPUs.
! Also returns the norm (s).
!============================================================================
subroutine normalize(n, x, s, in_parallel, ierr)
  integer, intent(in) :: n
  logical, intent(in) :: in_parallel

  real(kind=CUSTOM_REAL), intent(out) :: s
  real(kind=CUSTOM_REAL), intent(inout) :: x(n)
  integer, intent(out) :: ierr

  real(kind=CUSTOM_REAL) :: ss, s_glob

  if (in_parallel) then
    ! Calculate the L2 norm of a vector x that is split between CPUs.
    s = sum(x**2)
    call MPI_Allreduce(s, s_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    s = sqrt(s_glob)
  else
    s = norm2(x)
  endif

  ierr = 0
  if (s /= 0._CUSTOM_REAL) then
    ss = 1._CUSTOM_REAL / s
  else
    ierr = -1
    return
  endif

  x = ss * x

end subroutine normalize

end module lsqr_solver
