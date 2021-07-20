
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
! Least Square (LSQR) solver parallelized by the model parameters.
!
! This solver is unit tested in tests_lsqr.f90 and in tests_method_of_weights.f90.
! The unit tests are available for serial and parallel versions.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!=============================================================================
module lsqr_solver

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sparse_matrix

  implicit none

  private

  public :: lsqr_solve

  private :: get_norm_parallel
  private :: normalize
  private :: apply_soft_thresholding

  ! Input scalar parameters used in LSQR solver.
  type, public :: t_parameters_lsqr
    integer :: niter
    real(kind=CUSTOM_REAL) :: rmin
    real(kind=CUSTOM_REAL) :: gamma
  end type t_parameters_lsqr

contains

!================================================================================
! LSQR solver for a sparse matrix.
! Solve min||Ax - b||, where A is stored in sparse format in 'matrix' variable.
!
! Niter - maximum number of iterations.
! rmin - stopping criterion (relative residual).
! gamma - soft thresholding parameter (see ISTA, proximate L1 norm). Use gamma=0 for pure L2 norm.
!================================================================================
subroutine lsqr_solve(niter, rmin, gamma, matrix, b, x, myrank)
  integer, intent(in) :: myrank
  integer, intent(in) :: niter
  real(kind=CUSTOM_REAL), intent(in) :: rmin, gamma
  real(kind=CUSTOM_REAL), intent(in) :: b(:)
  type(t_sparse_matrix), intent(inout) :: matrix

  real(kind=CUSTOM_REAL), intent(inout) :: x(:)

  ! Local variables.
  integer :: iter, ierr
  integer :: N_lines
  ! Local (at current CPU) number of parameters.
  integer :: nelements
  real(kind=CUSTOM_REAL) :: alpha, beta, rho, rhobar, phi, phibar, theta
  real(kind=CUSTOM_REAL) :: b1, c, r, g, s, t1, t2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: v, w, u, v0
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: Hv, Hv_loc
  real(kind=CUSTOM_REAL) :: rho_inv
  ! Flag for calculating variance.
  logical :: calculateVariance

  if (myrank == 0) print *, 'Entered subroutine lsqr_solve, gamma =', gamma

  N_lines = matrix%get_total_row_number()
  nelements = size(x)
  
  if (allocated(matrix%lsqr_var) .and. size(matrix%lsqr_var) == nelements) then
    calculateVariance = .true.
  else
    calculateVariance = .false.
  endif

  ! Allocate memory.
  allocate(u(N_lines))
  allocate(Hv(N_lines))
  allocate(Hv_loc(N_lines))
  allocate(v0(nelements))
  allocate(v(nelements))
  allocate(w(nelements))

  ! Sanity check.
  if (norm2(b) == 0.d0) then
    call exit_MPI("|b| = 0, prior model is exact? Exiting.", myrank, 0)
  end if

  ! Initialization.
  u = b

  ! Normalize u and initialize beta.
  if (.not. normalize(u, beta, .false.)) then
    call exit_MPI("Could not normalize initial u, zero denominator!", myrank, 0)
  endif

  b1 = beta
  ! Required by the algorithm.
  x = 0._CUSTOM_REAL

  ! Compute v = Ht.u.
  call matrix%trans_mult_vector(u, v)

  ! Normalize v and initialize alpha.
  if (.not. normalize(v, alpha, .true.)) then
    call exit_MPI("Could not normalize initial v, zero denominator!", myrank, 0)
  endif

  rhobar = alpha
  phibar = beta
  w = v

  !print *, 'myrank, alpha0=', myrank, alpha

  iter = 1
  r = 1._CUSTOM_REAL

  ! Use an exit (threshold) criterion in case we can exit the loop before reaching the max iteration count.
  do while (iter <= niter .and. r > rmin)

    ! Scale u.
    u = - alpha * u

    ! Compute u = u + H.v parallel.
    call matrix%mult_vector(v, Hv_loc)
    call mpi_allreduce(Hv_loc, Hv, N_lines, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

    u = u + Hv

    ! Normalize u and update beta.
    if (.not. normalize(u, beta, .false.)) then
      ! Achieved an exact solution. Happens in the unit test test_lsqr_underdetermined_2.
      if (myrank == 0) print *, 'WARNING: u = 0. Possibly found an exact solution in the LSQR solver!'
    endif

    ! Scale v.
    v = - beta * v

    ! Compute v = v + Ht.u.
    call matrix%trans_mult_vector(u, v0)

    v = v + v0

    ! Normalize v and update alpha.
    if (.not. normalize(v, alpha, .true.)) then
      ! Achieved an exact solution. Happens in the unit test test_method_of_weights_1.
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
      call matrix%mult_vector(x, Hv_loc)
      call mpi_allreduce(Hv_loc, Hv, N_lines, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Norm of the relative residual.
      r = norm2(Hv - b) / b1
    else
      ! Norm of the relative residual.
      r = phibar / b1
    endif

#ifndef SUPPRESS_OUTPUT
    if (mod(iter, 10) == 0) then
      ! Calculate the gradient: 2A'(Ax - b).
      call matrix%mult_vector(x, Hv_loc)
      call mpi_allreduce(Hv_loc, Hv, N_lines, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      Hv = Hv - b
      call matrix%trans_mult_vector(Hv, v0)

      ! Norm of the gradient.
      g = 2.d0 * get_norm_parallel(v0)

      if (myrank == 0) print *, 'it, r, g =', iter, r, g
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
  endif

  if (myrank == 0) print *, 'End of subroutine lsqr_solve, r =', r, ' iter =', iter - 1

#ifdef USE_FLUSH6
  call flush(6)
#endif

  deallocate(u)
  deallocate(Hv)
  deallocate(Hv_loc)
  deallocate(v0)
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
! Returns L2 norm of a vector x that is split between CPUs.
!============================================================================
function get_norm_parallel(x) result(s)
  real(kind=CUSTOM_REAL), intent(in) :: x(:)
  real(kind=CUSTOM_REAL) :: s

  integer :: ierr
  real(kind=CUSTOM_REAL) :: s0

  s = sum(x**2)

  call mpi_allreduce(s, s0, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
  s = sqrt(s0)

end function get_norm_parallel

!============================================================================
! Normalizes vector x by its L2 norm.
! Set in_parallel=true, if x is split between CPUs.
! Also returns the norm (s).
!============================================================================
function normalize(x, s, in_parallel) result(res)
  logical, intent(in) :: in_parallel
  real(kind=CUSTOM_REAL), intent(out) :: s
  real(kind=CUSTOM_REAL), intent(inout) :: x(:)
  logical :: res

  real(kind=CUSTOM_REAL) :: ss

  res = .true.

  if (in_parallel) then
    s = get_norm_parallel(x)
  else
    s = norm2(x)
  endif

  if (s /= 0._CUSTOM_REAL) then
    ss = 1._CUSTOM_REAL / s
  else
    res = .false.
    return
  endif

  x = ss * x

end function normalize

end module lsqr_solver
