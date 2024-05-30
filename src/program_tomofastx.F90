
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
! If you use this code, please cite the following two papers:
!
! [1] V. Ogarko, K. Frankcombe, T. Liu, J. Giraud, R. Martin, and M. Jessell (2024),
!     "Tomofast-x 2.0: an open-source parallel code for inversion of potential field data with topography using wavelet compression",
!     Geosci. Model Dev., 17, 2325–2345, https://doi.org/10.5194/gmd-17-2325-2024
!
! [2] J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell (2021),
!     "Structural, petrophysical, and geological constraints in potential field inversion using the Tomofast-x v1.0 open-source code",
!     Geosci. Model Dev., 14, 6681–6709, https://doi.org/10.5194/gmd-14-6681-2021
!
!========================================================================

program program_tomofastx

  use global_typedefs
  use ftnunit, only: runtests_init, runtests, runtests_final
  use unit_tests, only: test_all
  use parameters_grav
  use parameters_mag
  use parameters_inversion
  use init_parameters
  use problem_joint_gravmag
  use memory_tools

  implicit none

  ! MPI variables (error code, rank of this process, total number of ranks).
  integer :: ierr, myrank, nbproc

  ! Gravity (forward) problem parameters.
  type(t_parameters_grav) :: gpar
  ! Magnetism (forward) problem parameters.
  type(t_parameters_mag) :: mpar
  ! Inversion parameters.
  type(t_parameters_inversion) :: ipar
  ! Type of problem to solve.
  integer :: problem_type
  ! Memory usage.
  real(kind=CUSTOM_REAL) :: memory

  !----------------------------------------------------------------------------
  ! These initializations will work for both serial (no MPI) and parallel runs,
  ! because correctness is guaranteed by the MPI standard.
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierr)

  memory = get_max_mem_usage()
  if (myrank == 0) print *, "MEMORY USED (mpi init) [GB] =", memory

  !----------------------------------------------------------------------------
  ! UNIT TESTING.
  ! The routine runtests will check if unit tests are requested (if a file ftnunit.run exists).
  ! If not, it will return immediately.
  ! This way we make sure the unit tests remain part of the program.
  ! The routine test_all runs all unit tests (see the unit_tests module).
  call runtests_init
  call runtests(test_all, myrank, nbproc)
  call runtests_final(myrank)

  !----------------------------------------------------------------------------
  ! INITIALIZATION.
  if (myrank == 0) print *, "Started Tomofast-x, version >= 2.0.3"

  if (command_argument_count() /= 2) then
    if (myrank == 0) print *, "Usage: tomofastx -p <Parfile_path>"
    stop
  endif

  ! Get the problem type from the command line.
  call get_problem_type(problem_type, myrank)

  ! Read the Parfile and initialize forward and inverse problem parameters.
  call initialize_parameters(problem_type, gpar, mpar, ipar, myrank, nbproc)

  ! Create the output directory. If it already exists there is no problem.
  if (myrank == 0) call execute_command_line('mkdir -p "'//trim(path_output)//'"')

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !----------------------------------------------------------------------------
  ! MAIN CALCULATIONS.
  if (problem_type == 1) then
    ! Gravity and Magnetism problem.
    call solve_problem_joint_gravmag(gpar, mpar, ipar, myrank, nbproc)
  endif

  memory = get_max_mem_usage()
  if (myrank == 0) print *, "MEMORY USED (end) [GB] =", memory

  if (myrank == 0) print *, "THE END."

  !----------------------------------------------------------------------------
  ! ALL DONE.
  call MPI_Finalize(ierr)

end program program_tomofastx
