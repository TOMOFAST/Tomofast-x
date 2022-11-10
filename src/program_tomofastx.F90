
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
! If you use this code for your own research, please cite these articles:
!
! [1] V. Ogarko, J. Giraud, R. Martin, and M. Jessell (2021),
!     "Disjoint interval bound constraints using the alternating direction method of multipliers
!     for geologically constrained inversion: Application to gravity data,"
!     GEOPHYSICS 86: G1-G11. https://doi.org/10.1190/geo2019-0633.1
!
! [2] J. Giraud, V. Ogarko, R. Martin, M. Lindsay, M. Jessell (2021),
!     "Structural, petrophysical and geological constraints in potential field inversion 
!     using the Tomofast-x open-source code", Geoscientific Model Development Discussions, 
!     https://doi.org/10.5194/gmd-2021-14
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

  implicit none

  ! MPI variables (error code, rank of this process, total number of ranks).
  integer :: ierr, myrank, nbproc

  ! Gravity (forward) problem parameters.
  type(t_parameters_grav) :: gpar
  ! Magnetism (forward) problem parameters.
  type(t_parameters_mag) :: mpar
  ! Inversion parameters.
  type(t_parameters_inversion) :: ipar
  ! Type of problem (ECT/Grav/Mag).
  integer :: problem_type

  !----------------------------------------------------------------------------
  ! These initializations will work for both serial (no MPI) and parallel runs,
  ! because correctness is guaranteed by the MPI standard.
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierr)

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
  if (myrank == 0) print *, "Started Tomofast-x, version >= v.1.5"

  ! Get the problem type from the command line (ECT = 1, Grav/Mag = 2).
  call get_problem_type(problem_type, myrank)

  ! Read the Parfile and initialize forward and inverse problem parameters.
  call initialize_parameters(problem_type, gpar, mpar, ipar, myrank, nbproc)

  ! Create the output directory. If it already exists there is no problem.
  if (myrank == 0) call execute_command_line('mkdir -p '//path_output)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !----------------------------------------------------------------------------
  ! MAIN CALCULATIONS.
  if (problem_type == 1) then

  else if (problem_type == 2) then
  ! Gravity and Magnetism problem.
    call solve_problem_joint_gravmag(gpar, mpar, ipar, myrank, nbproc)
  endif

  if (myrank == 0) print *, "THE END."

  !----------------------------------------------------------------------------
  ! ALL DONE.
  call MPI_Finalize(ierr)

end program program_tomofastx

