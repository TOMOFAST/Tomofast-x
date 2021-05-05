
!========================================================================
!
!                    T O M O F A S T X  Version 1.0
!                  ----------------------------------
!
!              Main authors: Vitaliy Ogarko, Roland Martin,
!                   Jeremie Giraud, Dimitri Komatitsch.
! CNRS, France, and University of Western Australia.
! (c) CNRS, France, and University of Western Australia. January 2018
!
!
! If you use this code for your own research, please cite some (or all) of
! these articles:
!
!  @article{MaMoKoPeJeBoLi2013,
!  author = {Martin, Roland and Monteiller, Vadim and Komatitsch, Dimitri
!  and Perrouty, Stephane and Jessell, Mark and Bonvalot, Sylvain and Lindsay, Mark},
!  title = {Gravity inversion using wavelet-based compression on parallel hybrid
!  {CPU/GPU} systems: application to southwest {G}hana},
!  volume = {195},
!  number = {3},
!  pages = {1594-1619},
!  year = {2013},
!  doi = {10.1093/gji/ggt334},
!  journal = {Geophysical Journal International}}
!
!
! This software is a computer program whose purpose is to perform
! capacitance, gravity, magnetic, or joint gravity and magnetic tomography.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================


program program_tomofast3D

  use stopwatch
  use timers, only: get_CPU_times
  use global_typedefs
  use ftnunit, only: runtests_init, runtests, runtests_final
  use unit_tests, only: test_all
  use parameters_ect
  use parameters_grav
  use parameters_mag
  use parameters_inversion
  use init_parameters
  use problem_ect
  use problem_joint_gravmag

  implicit none

  ! MPI variables (error code, rank of this process, total number of ranks).
  integer :: ierr, myrank, nbproc

  ! Timer variables.
  type(watchtype) :: watch_prog
  real :: time_accum, time_min, time_max

  ! ECT (forward) problem parameters.
  type(t_parameters_ect) :: epar
  ! Gravity (forward) problem parameters.
  type(t_parameters_grav) :: gpar
  ! Magnetism (forward) problem parameters.
  type(t_parameters_mag) :: mpar
  ! Inversion parameters.
  type(t_parameters_inversion) :: ipar
  ! Type of problem (ECT/Grav/Mag).
  integer :: problem_type
  ! Joint problem of gravity and magnetism.
  type(t_problem_joint_gravmag) :: prj

  !----------------------------------------------------------------------------
  ! These initializations will work for both serial (no MPI) and parallel runs,
  ! because correctness is guaranteed by the MPI standard.
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierr)

  ! Start timers.
  call create_watch(watch_prog)
  call start_watch(watch_prog)

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
  ! Get the problem type from the command line (ECT = 1, Gravity = 2, Magnetism = 3).
  call get_problem_type(problem_type, myrank)

  ! Read the Parfile and initialize forward and inverse problem parameters.
  call initialize_parameters(problem_type, epar, gpar, mpar, ipar, myrank, nbproc)

  ! Create output directory. If it already exists there is no problem.
  ! TODO: replace with a proper use of execute_command_line().
  call system('mkdir -p '//path_output)

  !----------------------------------------------------------------------------
  ! MAIN CALCULATIONS.
  if (problem_type == 1) then
  ! ECT problem.
    call solve_problem_ect(epar, ipar, myrank, nbproc)

  else if (problem_type == 2 .or. problem_type == 3) then
  ! Gravity and Magnetism problems.
    !call solve_problem_gravmag(gpar, mpar, ipar, problem_type, myrank, nbproc)

  else if (problem_type == 4) then
  ! Joint Gravity and Magnetism problem.

    call prj%initialize(ipar%nelements, myrank)

    call prj%solve_problem_joint_gravmag(gpar, mpar, ipar, myrank, nbproc)

  endif

  if (myrank == 0) print *, 'THE END.'

  call stop_watch(watch_prog)
  call get_CPU_times(watch_prog, time_accum, time_min, time_max)
  call destroy_watch(watch_prog)

  if (myrank == 0) print *, 'Total CPU accumulated time: ', time_accum
  if (myrank == 0) print *, 'Total CPU minimum time:     ', time_min
  if (myrank == 0) print *, 'Total CPU maximum time:     ', time_max

  !----------------------------------------------------------------------------
  ! ALL DONE.
#ifdef USE_FLUSH6
  call flush(6)
#endif

  call MPI_Finalize(ierr)

end program program_tomofast3D

