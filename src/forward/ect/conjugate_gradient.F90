
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
!
! (c) CNRS, France, and University of Western Australia. January 2018
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

module conjugate_gradient

  use global_typedefs
  use m_matmul
  use utils, only: compute_snrm

  implicit none

  private

  ! Container for all auxiliary arrays needed by the conjugate gradient solver.
  type, public :: t_pcg_auxarrays
    ! descent direction
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: p
    ! residual vector
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: r
    ! inverse of the matrix diagonal, for Jacobi preconditioning
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: inv_diagA
    ! preconditioned residual
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: z
    ! temporary array, for Jacobi preconditioning
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: tmp
  end type t_pcg_auxarrays

  public :: solver_pcg
  public :: mpisendrecv

  private :: calculate_precon_residual

contains

!=========================================================================================================
! Exchanges the boundary values, for e.g., sparse matrix to vector multiplication in Conjugate Gradients.
!
! Algorithm to do the exchange: classical "even ones send, odd ones receive,
! and vise versa", with some extra care for the first and last subdomain,
! recall that we are doing MPI among entire (x,y) chunks along the z axis only
!
! The data that needs to be send is always nr*ntheta values, stored contiguously in memory.
!=========================================================================================================
subroutine mpisendrecv(x, nr, ntheta, nz, myrank, nbproc)
  ! Problem dimensions.
  integer , intent(in):: nr, ntheta, nz
  ! Local part of the distributed vector to perform neighbour exchanges on.
  real(kind=CUSTOM_REAL), intent(inout) :: x(0:nr+1, 0:ntheta+1, 0:nz+1)
  ! MPI rank and total number of processes.
  integer, intent(in) :: myrank,nbproc
  ! MPI error code.
  integer :: ierr

  integer :: nb_of_elems_to_send
  ! MPI status variable (needed for calls to mpi_recv).
  integer :: mpi_stat(MPI_STATUS_SIZE)

  ! Only do actual work if there is more than one rank total
  ! otherwise, this routine is empty.
  if(nbproc > 1) then

    nb_of_elems_to_send = (nr+2)*(ntheta+2)

    !! There might be room for improvement here by switching to non-blocking MPI
    !! (if we can overlap communications with calculations; let us see if and how).

    if (mod(myrank, 2) == 0) then
      if (myrank < nbproc-1) then
        call mpi_send(x(0,0,nz),   nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 42, MPI_COMM_WORLD, ierr)
        call mpi_recv(x(0,0,nz+1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 43, MPI_COMM_WORLD, mpi_stat, ierr)
      endif

      if (myrank > 0) then
        call mpi_send(x(0,0,1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 44, MPI_COMM_WORLD, ierr)
        call mpi_recv(x(0,0,0), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 45, MPI_COMM_WORLD, mpi_stat, ierr)
      endif

    else

      if (myrank > 0) then
        call mpi_recv(x(0,0,0), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 42, MPI_COMM_WORLD, mpi_stat, ierr)
        call mpi_send(x(0,0,1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank-1, 43, MPI_COMM_WORLD, ierr)
      endif

      if (myrank < nbproc-1) then
        call mpi_recv(x(0,0,nz+1), nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 44, MPI_COMM_WORLD, mpi_stat, ierr)
        call mpi_send(x(0,0,nz),   nb_of_elems_to_send, CUSTOM_MPI_TYPE, myrank+1, 45, MPI_COMM_WORLD, ierr)
      endif

    endif

  endif

end subroutine mpisendrecv

!======================================================================================================
! Calculate preconditioned residual.
! z = omega * r /diagA followed by precond arr = M^{-1}r.
!======================================================================================================
subroutine calculate_precon_residual(auxarrays, A, iprecond, omega1, nr, ntheta, nz, myrank, nbproc)
  integer, intent(in) :: nr, ntheta, nz
  real(kind=CUSTOM_REAL), intent(in) :: A(nr, ntheta, nz, 7)
  real(kind=CUSTOM_REAL), intent(in) :: omega1
  integer, intent(in) :: iprecond
  integer, intent(in) :: myrank, nbproc

  type (t_pcg_auxarrays), intent(inout) :: auxarrays

  integer :: i, j, k
  integer :: itt

  if (iprecond == 0) then
  ! No preconditioner: just divide by matrix diagonal.
    do k = 1, nz
      do j = 1, ntheta
        do i = 1, nr
          auxarrays%z(i,j,k) = auxarrays%r(i,j,k) * auxarrays%inv_diagA(i,j,k)
        enddo
      enddo
    enddo

  else
  ! Apply Jacobi preconditioner.
    ! Initialization (z = omega * r /diagA).
    do k = 1, nz
      do j = 1, ntheta
        do i = 1, nr
          auxarrays%z(i,j,k) = omega1 * auxarrays%r(i,j,k) * auxarrays%inv_diagA(i,j,k)
        enddo
      enddo
    enddo

    ! Preconditioning (z = M^{-1}r).
    if (.NOT. allocated(auxarrays%tmp)) then
      print *, "Error! Array auxarrays%tmp is not allocated in calculate_precon_residual!"
      stop
    endif
    do itt = 1, 1
      ! Exchange halo values with neighbouring z-slices
      if (nbproc > 1) call mpisendrecv(auxarrays%z, nr, ntheta, nz, myrank, nbproc)
      call periodic_and_matmul(A, auxarrays%z, auxarrays%tmp, nr, ntheta, nz)

      do k = 1, nz
        do j = 1, ntheta
          do i = 1, nr
            auxarrays%z(i,j,k) = auxarrays%z(i,j,k) + &
              omega1 * (auxarrays%r(i,j,k) - auxarrays%tmp(i,j,k)) * auxarrays%inv_diagA(i,j,k)
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine calculate_precon_residual

!====================================================================================
! (Jacobi) Preconditioned conjugate gradient solver.
!====================================================================================
subroutine solver_pcg(A, b, x, itypenorm, iprecond, omega1, tol, itmax, iter, &
                      output_frequency, suppress_output, &
                      nr, ntheta, nz, myrank, nbproc, auxarrays)

  ! Number of unknowns in each dimension.
  integer, intent(in) :: nr, ntheta, nz
  ! Matrix (7-diagonal).
  real(kind=CUSTOM_REAL), intent(in) :: A(nr, ntheta, nz, 7)
  ! RHS (right-hand side).
  real(kind=CUSTOM_REAL), intent(in) :: b(0:nr+1, 0:ntheta+1, 0:nz+1)
  ! Type of norm to use (1=L2, 2=max)
  integer, intent(in) :: itypenorm
  ! Preconditioning (1 = Jacobi, 0 = diagonal division = no preconditioning).
  integer, intent(in) :: iprecond
  ! Relaxation factor for preconditioning.
  real(kind=CUSTOM_REAL), intent(in) :: omega1
  ! Tolerance to reach (reduction of residual by this factor).
  real(kind=CUSTOM_REAL), intent(in) :: tol
  ! Maximum number of iterations.
  integer, intent(in) :: itmax
  ! frequency at which convergence is checked and printed out
  integer, intent(in) :: output_frequency
  ! suppress output except in emergencies?
  logical, intent(in) :: suppress_output
  ! MPI variables
  integer, intent(in) :: myrank, nbproc

  ! auxiliary arrays
  type (t_pcg_auxarrays), intent(inout) :: auxarrays
  ! iteration (solution) vector
  real(kind=CUSTOM_REAL), intent(inout) :: x(0:nr+1, 0:ntheta+1, 0:nz+1)
  ! actual number of iterations (reported back to calling routine)
  integer, intent(out) :: iter

  ! Local variables:

  ! some counters and misc variables
  integer :: i, j, k, ierr
  ! result of norm computations (local part)
  real(kind=CUSTOM_REAL) :: err
  ! steplength parameter, dot(r,p), both locally (per process intermediate value)
  ! and also fully mpi_reduced globally
  real(kind=CUSTOM_REAL) :: bkold, bknew, bkold_glob, bknew_glob
  ! DG DG a few more dots and norms that need documentation
  real(kind=CUSTOM_REAL) :: ak, ak1, akden, bk, bnrm
  real(kind=CUSTOM_REAL) :: akden_glob, err_glob
  real(kind=CUSTOM_REAL) :: bnrm_glob

#ifdef USE_TIMERS_OLD
  ! timers
  ! DG DG TODO (later): switch to the same stopwatch timer we use in part II
  real(kind=CUSTOM_REAL) time_new,time_old, time_init
  real(kind=CUSTOM_REAL) time_Ax_new, time_Ax_old
  real(kind=CUSTOM_REAL) time_glob_min, time_glob_max
  real(kind=CUSTOM_REAL) time_Ax_tot, time_Ax_tot_max, time_Ax_tot_min
#endif

  if (.not. suppress_output .and. myrank == 0) print *, 'STARTING CG'

  !
  ! Compute initial residual: r = b - Ax.
  !
  ! exchange halo values with neighbouring z-slices
  if (nbproc > 1) call mpisendrecv(x, nr, ntheta, nz, myrank, nbproc)
  ! r = Ax (temporary)
  call periodic_and_matmul(A, x, auxarrays%r, nr, ntheta, nz)
  ! r = b - Ax
  do k = 1, nz
    do j = 1, ntheta
      do i = 1, nr
        auxarrays%r(i,j,k) = b(i,j,k) - auxarrays%r(i,j,k)
      enddo
    enddo
  enddo

  !
  ! Compute initial preconditioned residual z.
  ! z0 = r/diagA = (b-Ax)/diagA
  !
  auxarrays%z = 0.d0
  call calculate_precon_residual(auxarrays, A, iprecond, omega1, nr, ntheta, nz, myrank, nbproc)

  !
  ! Compute initial descent direction, in preconditioned CG.
  ! This is simply the preconditioned initial residual.
  !
  auxarrays%p = auxarrays%z

  !
  ! Compute initial steplength as the dot product of r and p.
  ! As for all reductions, do this first locally and then globally.
  !
  ! initializations
  bkold = 0._CUSTOM_REAL
  bkold_glob = 0._CUSTOM_REAL
  ! local computations
  do k = 1, nz
    do j = 1, ntheta
      do i = 1, nr
        ! DG DG -> DK DK: possible BLAS call
        bkold = bkold + auxarrays%p(i,j,k) * auxarrays%r(i,j,k)
      enddo
    enddo
  enddo

  ! global computations (with the usual switch for serial mode)
  ! we store the global result in the local variable explicitly when we are done,
  ! to avoid problems with serial mode.
  if (nbproc > 1) then
    call mpi_allreduce(bkold, bkold_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    bkold = bkold_glob
  endif

  call compute_snrm(auxarrays%z, itypenorm, nr, ntheta, nz, 1, 1, 1, nr, ntheta, nz, bnrm)

  ! bnrm = inverse of L2 norm of z0 = inverse of preconditioned initial residual
  bnrm_glob = 0._CUSTOM_REAL
  if (nbproc > 1) then
    call mpi_allreduce(bnrm, bnrm_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
    bnrm = bnrm_glob
  endif

  if (bnrm /= 0._CUSTOM_REAL) then
    bnrm = 1._CUSTOM_REAL/sqrt(bnrm)
  else
    print *, "Error! bnrm = 0."
    stop
  endif

  if (.not.suppress_output) then
    if (myrank == 0) print *, 'CG iter=', 0, ' residual error=', 1._CUSTOM_REAL / bnrm
  endif

  !
  ! update timers
  !
#ifdef USE_TIMERS_OLD
  call cpu_time(time_init)
  call cpu_time(time_old)

  time_Ax_tot=0._CUSTOM_REAL
  time_Ax_tot_max=0._CUSTOM_REAL
  time_Ax_tot_min=0._CUSTOM_REAL
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     CONJUGATE GRADIENT MAIN LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! here we enter the iterative process
  do iter = 1, itmax

    !
    ! z=Ap
    !
#ifdef USE_TIMERS_OLD
    call cpu_time(time_Ax_old)
#endif
    if (nbproc > 1) call mpisendrecv(auxarrays%p, nr, ntheta, nz, myrank, nbproc)
    call periodic_and_matmul(A, auxarrays%p, auxarrays%z, nr, ntheta, nz)
#ifdef USE_TIMERS_OLD
    call cpu_time(time_Ax_new)
    time_Ax_tot = time_Ax_tot + time_Ax_new - time_Ax_old
#endif

    !
    ! (p,z)=(p,Ap)
    !
    akden = 0._CUSTOM_REAL
    akden_glob = 0._CUSTOM_REAL
    do k = 1, nz
      do j = 1, ntheta
        do i = 1, nr
          akden = akden + auxarrays%z(i,j,k) * auxarrays%p(i,j,k)
        enddo
      enddo
    enddo

    if (nbproc > 1) then
      call mpi_allreduce(akden, akden_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      akden = akden_glob
    endif

    !
    ! alpha=(r,M^{-1}r)/(p,Ap)
    ! note that z=M^{-1}r
    !
    ak = bkold / akden
    ak1 = - ak

    !
    ! x=x+alpha*p
    ! r=r-alpha*Ap
    !
    do k = 1, nz
      do j = 1, ntheta
        do i = 1, nr
          x(i,j,k) = x(i,j,k) + ak * auxarrays%p(i,j,k)
          auxarrays%r(i,j,k) = auxarrays%r(i,j,k) + ak1 * auxarrays%z(i,j,k)
        enddo
      enddo
    enddo


    !
    ! Compute again the new preconditioned residual.
    !
    !auxarrays%z = 0.d0
    call calculate_precon_residual(auxarrays, A, iprecond, omega1, nr, ntheta, nz, myrank, nbproc)

    !
    ! (r,z)
    !
    bknew = 0._CUSTOM_REAL
    bknew_glob = 0._CUSTOM_REAL
    do k = 1, nz
      do j = 1, ntheta
        do i = 1, nr
          bknew = bknew + auxarrays%r(i,j,k) * auxarrays%z(i,j,k)
        enddo
      enddo
    enddo

    if (nbproc > 1) then
      call mpi_allreduce(bknew, bknew_glob, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      bknew = bknew_glob
    endif

    bk = bknew / bkold

    !
    ! p = (r,z)*p + z
    !
    do k = 1, nz
      do j = 1, ntheta
        do i = 1, nr
          auxarrays%p(i,j,k) = bk * auxarrays%p(i,j,k) + auxarrays%z(i,j,k)
        enddo
      enddo
    enddo

    bkold = bknew

    !
    ! compute norm of residual for convergence check
    ! this is done only every output_frequency iteration to avoid global reductions in every iteration.
    ! we compute normL2 of z again = err=err_glob.
    !
    if (iter == 1 .or. (mod(iter, output_frequency) == 0 .and. iter > 1) .or. iter == itmax) then

      call compute_snrm(auxarrays%z, itypenorm, nr, ntheta, nz, 1, 1, 1, nr, ntheta, nz, err)
      err_glob = 0._CUSTOM_REAL
      if (nbproc > 1) then
        call mpi_allreduce(err, err_glob, 1, CUSTOM_MPI_TYPE,MPI_SUM, MPI_COMM_WORLD, ierr)
        err = err_glob
      endif
      ! we compute the RELATIVE error err=sqrt(err)*bnrm=normL2(z)/normL2(zo)
      ! This is a relative error from Roland's point of view.
      ! Now, if you want we can change the criterion: using relative errors on residuals r and not on preconditioned residuals z.
      ! The idea of using preconditionned residuals z is to have a criterion close to relative errors on solution x.
      err = sqrt(err) * bnrm

      if (.not.suppress_output) then
        if (myrank == 0) then
          print *, 'CG iter=', iter, ' residual error=', err
        endif
      endif

#ifdef USE_TIMERS_OLD
      call cpu_time(time_new)
      time_old=time_new
      if(nbproc>1) then
        call mpi_allreduce(time_Ax_tot, time_Ax_tot_max, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(time_Ax_tot, time_Ax_tot_min, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(time_new - time_init, time_glob_max, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(time_new - time_init, time_glob_min, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD, ierr)
      else
        time_Ax_tot_max = time_Ax_tot
        time_Ax_tot_min = time_Ax_tot
        time_glob_max = time_new - time_init
        time_glob_min = time_new - time_init
      endif

      if (myrank == 0) then
        print *,'time glob =',time_new-time_init, 'time glob min=',time_glob_min, &
                'time glob max =',time_glob_max, 'time Ax glob =',time_Ax_tot, 'time Ax glob min =',time_Ax_tot_min, &
                'time Ax glob max =',time_Ax_tot_max
        print *
      endif
#endif

#ifdef USE_FLUSH6
      flush(6)
#endif

      !
      ! final convergence check
      ! note that this is actually a RELATIVE convergence criterion, because err is multiplied
      ! by bnrm directly before the timing stuff, which is in turn the reciprocal of the initial residual
      !
      if (err <= tol) exit

    endif
  enddo ! main conjugate gradient loop

  !
  ! final output (note that in case we actually suppress output, we do not suppress if something went haywire)
  !
  if (myrank == 0) then
    if (suppress_output .and. iter >= itmax) print *, "CG: max. iters reached!", err

#ifdef USE_TIMERS_OLD
    if (.not.suppress_output) print *, 'CG end of iterations: iter=',iter,' residual error=',err,&
                                       'time =',time_new-time_init,'time Ax glob =',time_Ax_tot
#else
    if (.not.suppress_output) print *, 'CG end of iterations: iter=',iter,' residual error=',err
#endif

#ifdef USE_FLUSH6
    flush(6)
#endif
  endif

  ! Exchange halo values with neighbouring z-slices.
  ! Need them to calculate capacitance in sensitivity() function.
  if (nbproc > 1) call mpisendrecv(x, nr, ntheta, nz, myrank, nbproc)

end subroutine solver_pcg

end module conjugate_gradient

