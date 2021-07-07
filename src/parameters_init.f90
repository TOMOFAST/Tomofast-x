
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

module init_parameters

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use sanity_check
  use parameters_ect
  use parameters_grav
  use parameters_mag
  use parameters_inversion
  use geometry, only: get_refined_ntheta
  use parallel_tools

  implicit none

  private

  public :: initialize_parameters
  public :: get_problem_type

  private :: read_parfile
  private :: read_parfile2
  private :: set_default_parameters

  private :: print_arg_int
  private :: print_arg_dbl
  private :: read_filename

  ! For overloading the 'print_arg' function.
  interface print_arg
    module procedure print_arg_int, print_arg_dbl, print_arg_str
  end interface

contains

!==========================================================================
! Functions for printing the parameter values.
!==========================================================================
subroutine print_arg_int(myrank, name, value)
  character(len=128), intent(in) :: name
  integer, intent(in) :: myrank
  integer, intent(in) :: value

  if (myrank == 0) print *, trim(name)//" =", value
end subroutine print_arg_int
!-----------------------------------------------
subroutine print_arg_dbl(myrank, name, value)
  character(len=128), intent(in) :: name
  real(kind=CUSTOM_REAL), intent(in) :: value
  integer, intent(in) :: myrank

  if (myrank == 0) print *, trim(name)//" =", value
end subroutine print_arg_dbl
!-----------------------------------------------
subroutine print_arg_str(myrank, name, value)
  character(len=128), intent(in) :: name
  character(len=256), intent(in) :: value
  integer, intent(in) :: myrank

  if (myrank == 0) print *, trim(name)//" =", trim(value)
end subroutine print_arg_str

!==========================================================================
! Read the filename from file.
!==========================================================================
subroutine read_filename(file_id, val)
  integer, intent(in) :: file_id
  character(len=256), intent(out) :: val

  read(file_id, '(a)') val
  val = adjustl(val)
end subroutine read_filename

!==========================================================================
! Get problem type (ECT / Gravity).
!==========================================================================
subroutine get_problem_type(problem_type, myrank)
  integer, intent(in) :: myrank
  integer, intent(out) :: problem_type
  character(len=256) :: arg

  ! ECT problem is set by default.
  arg = '-e'
  problem_type = 1

  if (command_argument_count() > 0) call get_command_argument(1, arg)

  if (arg == '-e') then
    problem_type = 1
    if (myrank == 0) print *, '===== START ECT PROBLEM ====='
  else if (arg == '-g') then
    problem_type = 2
    if (myrank == 0) print *, '===== START GRAVITY PROBLEM ====='
  else if (arg == '-m') then
    problem_type = 3
    if (myrank == 0) print *, '===== START MAGNETISM PROBLEM ====='
  else if (arg == '-j') then
    problem_type = 4
    if (myrank == 0) print *, '===== START JOINT GRAV/MAG PROBLEM ====='
  else
    call exit_MPI("UNKNOWN PROBLEM TYPE! arg ="//arg, myrank, 0)
  endif

end subroutine get_problem_type

!=======================================================================================
! TODO: Split to ect, grav & mag routines.
! Initialize parameters for forward and inverse problems.
!=======================================================================================
subroutine initialize_parameters(problem_type, epar, gpar, mpar, ipar, myrank, nbproc)
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank,nbproc

  type(t_parameters_ect), intent(out) :: epar
  type(t_parameters_grav), intent(out) :: gpar
  type(t_parameters_mag), intent(out) :: mpar
  type(t_parameters_inversion), intent(out) :: ipar

  type(t_parameters_base) :: gmpar
  type(t_parallel_tools) :: pt
  integer :: nelements, ierr

  if (myrank == 0) then
    ! Setting default parameter values.
    call set_default_parameters(epar, gpar, mpar, ipar)

    ! Read Parfile data, only the master does this,
    ! and then broadcasts all the information to the other processes.
    call read_parfile2(epar, gpar, mpar, ipar, myrank)

    stop

    ! Global sanity checks.
    if (problem_type == 1) then
      ! ntheta is a multiple of the number of electrodes.
      call sanity_ntheta_nel(epar%dims%ntheta, epar%nel / epar%nrings, myrank)
      ! nz is a multiple of the number of processes.
      call sanity_nz(epar%dims%nz, nbproc, myrank)
    endif
  endif

  ! Use barrier to do not mix the parameters output from the master CPU with other log messages.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! All sanity checks passed so far, so we may broadcast the Parfile.

  ! MPI broadcast parameters.
  call MPI_Bcast(path_output, len(path_output), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  if (problem_type == 1) then
    call epar%broadcast()

  else if (problem_type == 2) then
    call gpar%broadcast(myrank)

  else if (problem_type == 3) then
    call mpar%broadcast(myrank)

  else if (problem_type == 4) then
    call gpar%broadcast(myrank)
    call mpar%broadcast(myrank)
  endif

  call ipar%broadcast(myrank)

  !--------------------------------
  ! Some extra initializations.
  !--------------------------------
  if (problem_type == 1) then
  ! Electrical capacitance tomography (ECT) problem.

    epar%read_guess_from_file = .false.

    !-----------------------
    ! PARTITIONING FOR MPI

    ! Cut the model into vertical slices for MPI,
    ! always divides evenly because of the previous sanity check.
    epar%dims%nzlocal = epar%dims%nz / nbproc

    ! Initialize dimensions of the sensitivity matrix in inverse problem (nelements x ndata+nelements).
    ipar%nelements_total = epar%dims%nr * epar%dims%ntheta * (epar%dims%nz + 1)

    ipar%nx = epar%dims%nr
    ipar%ny = epar%dims%ntheta

    ! We have nzlocal+1 k-elements in the last processor,
    ! since the model is defined in the middle of the potential (phi) grid nodes, which run from 0 to nz+1.
    ! So for nz=4, phi-nodes are:   0   1   2   3   4   5, and
    !              model-nodes are:   1   2   3   4   5, i.e., not even number.
    ipar%nz = epar%dims%nzlocal + (myrank + 1) / nbproc

    ipar%nelements = ipar%nx * ipar%ny * ipar%nz

    ipar%ndata = epar%get_ndata()

    !-----------------------------------
    ! CHANGE NTHETA FOR MESH REFINEMENT

    epar%dims%ntheta0 = epar%dims%ntheta
    ! Increase theta-dimension for mesh refinement.
    if (epar%irefine == 1) then
      ! Sanity checks.
      if (epar%linear_solver == LINSOLV_MG) then
        call exit_MPI("Mesh refinement is not implemented for the MG solver!", myrank, 0)
      endif
      if (epar%sens%space_electrodes == 0._CUSTOM_REAL) then
        call exit_MPI("Mesh refinement is implemented only for the case with gaps between electrodes!", myrank, 0)
      endif

      epar%dims%ntheta = get_refined_ntheta(epar%dims%ntheta, epar%nel)
      if (myrank == 0) print *, 'ntheta_read0, ntheta(new)', epar%dims%ntheta0, epar%dims%ntheta
    endif

  else if (problem_type == 2 .or. problem_type == 3 .or. problem_type == 4) then
  ! Gravity and magnetism problems.

    if (problem_type == 2 .or. problem_type == 4) then
    ! Gravity.

      gpar%ncomponents = 1

      gmpar = gpar%t_parameters_base
    endif

    if (problem_type == 3 .or. problem_type == 4) then
    ! Magnetism.

      mpar%ncomponents = 1

      gmpar = mpar%t_parameters_base
    endif

    ! Inverse problem parameters. -------------------------------
    ipar%nx = gmpar%nx
    ipar%ny = gmpar%ny
    ipar%nz = gmpar%nz

    ipar%ndata(1) = gpar%ndata
    ipar%ndata(2) = mpar%ndata

    ipar%nelements_total = ipar%nx * ipar%ny * ipar%nz

    ! Define model splitting for parallelization.
    nelements = pt%calculate_nelements_at_cpu(ipar%nelements_total, myrank, nbproc)

    ipar%nelements = nelements
    gpar%nelements = nelements
    mpar%nelements = nelements

  endif

  !---------------------------------------------------------------------------
  ! Print out some useful debug information.
  if (myrank == 0 .or. myrank == nbproc - 1) &
    print *, 'myrank=', myrank, ' nbproc=', nbproc, ' nelements_total=', ipar%nelements_total, &
             'nelements=', ipar%nelements, 'ndata=', ipar%ndata

end subroutine initialize_parameters

!===================================================================================
! Set default parameters.
!===================================================================================
subroutine set_default_parameters(epar, gpar, mpar, ipar)
  type(t_parameters_ect), intent(out) :: epar
  type(t_parameters_grav), intent(out) :: gpar
  type(t_parameters_mag), intent(out) :: mpar
  type(t_parameters_inversion), intent(out) :: ipar

  !-----------------------------------------------
  ! Define here the DEFAULT parameter values:
  !-----------------------------------------------

  ! DATA parameters.
  gpar%ndata = 0
  mpar%ndata = 0
  gpar%data_grid_file = "NILL"
  mpar%data_grid_file = "NILL"
  gpar%data_file = "NILL"
  mpar%data_file = "NILL"
  gpar%calc_data_directly = 0
  mpar%calc_data_directly = 0

  ! PRIOR MODEL parameters.
  gpar%prior_model_type = 1 ! 1-set value, 2-read from file.
  mpar%prior_model_type = 1
  gpar%number_prior_models = 1 ! Number of prior models, for the model type = 2
  mpar%number_prior_models = 1
  gpar%prior_model_val = 0.d0
  mpar%prior_model_val = 1.d-9
  gpar%model_files(2) = "NILL"
  mpar%model_files(2) = "NILL"

  ! STARTING MODEL parameters.
  gpar%start_model_type = 1 ! 1-set value, 2-read from file.
  mpar%start_model_type = 1
  gpar%start_model_val = 0.d0
  mpar%start_model_val = 1.d-9
  gpar%model_files(3) = "NILL"
  mpar%model_files(3) = "NILL"

  ! MAGNETIC FIELD constants.
  mpar%mi = 75.d0
  mpar%md = 25.d0
  mpar%fi = 75.d0
  mpar%fd = 25.d0
  mpar%intensity = 50000.d0
  mpar%theta = 0.d0

  ! DEPTH WEIGHTING parameters.
  gpar%depth_weighting_type = 3 ! 1-power, 2-sens, 3-isens
  ! Power weight function: W(Z) = 1 / (Z + Z0)**(beta / 2)
  gpar%beta = 1.4d0   ! TODO: update beta-power to a number from tests. Also to update the default value in ParametersDefault.md file
  gpar%Z0 = 0.d0
  mpar%beta = 1.4d0   ! TODO: update beta-power to a number from tests. Also to update the default value in ParametersDefault.md file
  mpar%Z0 = 0.d0

  ! MATRIX COMPRESSION parameters.
  gpar%distance_threshold = 1.d+10 ! Source to the cell distance (m).
  gpar%compression_rate = 1.d0     ! 1.0 = full matrix
  mpar%distance_threshold = gpar%distance_threshold
  mpar%compression_rate = gpar%compression_rate

  ! INVERSION parameters.
  ipar%ninversions = 10
  ipar%niter = 100
  ipar%rmin = 1.d-13
  ipar%method = 1 ! LSQR = 1
  ipar%gamma = 0. ! soft threshold ("L1-norm", no=0.)

  ! MODEL DAMPING (m - m_prior).
  ipar%alpha(1) = 1.d-9
  ipar%alpha(2) = 1.d-11
  ipar%norm_power = 2.d0

  ! JOINT INVERSION parameters.
  ipar%problem_weight(1) = 1.d0
  ipar%problem_weight(2) = 0.d0
  ipar%column_weight_multiplier(1) = 4.d+3
  ipar%column_weight_multiplier(2) = 1.d0
  ipar%niter_single(1) = 0
  ipar%niter_single(2) = 0

  ! DAMPING-GRADIENT constraints.
  ipar%damp_grad_weight_type = 1 ! 1-global, 2-local
  ipar%beta(1) = 0.d0
  ipar%beta(2) = 0.d0

  ! CROSS-GRADIENT constraints.
  ipar%cross_grad_weight = 0.d0
  ipar%method_of_weights_niter = 0
  ipar%derivative_type = 1 ! 1-fwd, 2-cent, 3-mixed

  ! CLUSTERING constraints.
  ipar%clustering_weight_glob(1) = 0.d0
  ipar%clustering_weight_glob(2) = 0.d0
  ipar%nclusters = 4
  ipar%mixture_file = "NILL"
  ipar%cell_weights_file = "NILL"
  ipar%clustering_opt_type = 2 ! 1-normal, 2-log
  ipar%clustering_constraints_type = 2 ! 1-global, 2-local

  ! ADMM constraints.
  ipar%admm_type = 0 ! 0-no admm, 1-enable admm
  ipar%nlithos = 5
  ipar%bounds_ADMM_file(1) = "NILL"
  ipar%bounds_ADMM_file(2) = "NILL"
  ipar%rho_ADMM(1) = 1.d-7
  ipar%rho_ADMM(2) = 1.d+5

end subroutine set_default_parameters

!===================================================================================
! Read input parameters from Parfile.
!===================================================================================
subroutine read_parfile2(epar, gpar, mpar, ipar, myrank)
  integer, intent(in) :: myrank

  type(t_parameters_ect), intent(inout) :: epar
  type(t_parameters_grav), intent(inout) :: gpar
  type(t_parameters_mag), intent(inout) :: mpar
  type(t_parameters_inversion), intent(inout) :: ipar

  integer :: itmp

  ! This is junk in order to ignore the variable name at the beginning of the line.
  ! This ignores exactly 40 characters.
  character(len=1) :: ch
  character(len=40) :: junk
  character(len=256) :: parfile_name
  character(len=128) :: parname
  character(len=256) :: line
  integer :: symbol_index, i
  integer :: ios

  ! The name of the Parfile can be passed in via the command line,
  ! if no argument is given, the default value is used.
  call get_command_argument(2,parfile_name)
  if (len_trim(parfile_name) == 0) parfile_name = "parfiles/Parfile_MASTER.txt"

  open(unit=10,file=parfile_name,status='old',iostat=itmp,action='read')
  if (itmp /= 0) call exit_MPI("Parfile """ // trim(parfile_name) // """ cannot be opened!",myrank,15)

  !---------------------------------------------------------------------------------
  ! Reading parameter values from Parfile.
  !---------------------------------------------------------------------------------
  do
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit

    symbol_index = index(line, '=')
    parname = line(:symbol_index - 1)
    backspace(10)

    ! Reading file line until the end of '=' symbol.
    do i = 1, symbol_index
      read(10, '(a)', advance='NO') ch
    enddo

    select case(trim(parname))

      ! DATA parameters -------------------------------------

      case("forward.gravmag.data.grav.nData")
        read(10, 2) gpar%ndata
        call print_arg(myrank, parname, gpar%ndata)

      case("forward.gravmag.data.mag.nData")
        read(10, 2) mpar%ndata
        call print_arg(myrank, parname, mpar%ndata)

      case("forward.gravmag.data.grav.dataGridFile")
        call read_filename(10, gpar%data_grid_file)
        call print_arg(myrank, parname, gpar%data_grid_file)

      case("forward.gravmag.data.mag.dataGridFile")
        call read_filename(10, mpar%data_grid_file)
        call print_arg(myrank, parname, mpar%data_grid_file)

      case("forward.gravmag.data.grav.dataValuesFile")
        call read_filename(10, gpar%data_file)
        call print_arg(myrank, parname, gpar%data_file)

      case("forward.gravmag.data.mag.dataValuesFile")
        call read_filename(10, mpar%data_file)
        call print_arg(myrank, parname, mpar%data_file)

      case("forward.gravmag.data.calcDataWithoutSensit")
        read(10, 2) gpar%calc_data_directly
        call print_arg(myrank, parname, gpar%calc_data_directly)
        mpar%calc_data_directly = gpar%calc_data_directly

      ! PRIOR MODEL -----------------------------------------

      case("forward.gravmag.priorModel.type")
        read(10, 2) gpar%prior_model_type
        call print_arg(myrank, parname, gpar%prior_model_type)
        mpar%prior_model_type = gpar%prior_model_type

      case("forward.gravmag.priorModel.nModels")
        read(10, 2) gpar%number_prior_models
        call print_arg(myrank, parname, gpar%number_prior_models)
        mpar%number_prior_models = gpar%number_prior_models

      case("forward.gravmag.priorModel.grav.value")
        read(10, 1) gpar%prior_model_val
        call print_arg(myrank, parname, gpar%prior_model_val)

      case("forward.gravmag.priorModel.mag.value")
        read(10, 1) mpar%prior_model_val
        call print_arg(myrank, parname, mpar%prior_model_val)

      case("forward.gravmag.priorModel.grav.file")
        call read_filename(10, gpar%model_files(2))
        call print_arg(myrank, parname, gpar%model_files(2))

      case("forward.gravmag.priorModel.mag.file")
        call read_filename(10, mpar%model_files(2))
        call print_arg(myrank, parname, mpar%model_files(2))

      ! STARTING MODEL -------------------------------------

      case("forward.gravmag.startingModel.type")
        read(10, 2) gpar%start_model_type
        call print_arg(myrank, parname, gpar%start_model_type)
        mpar%start_model_type = gpar%start_model_type

      case("forward.gravmag.startingModel.grav.value")
        read(10, 1) gpar%start_model_val
        call print_arg(myrank, parname, gpar%start_model_val)

      case("forward.gravmag.startingModel.mag.value")
        read(10, 1) mpar%start_model_val
        call print_arg(myrank, parname, mpar%start_model_val)

      case("forward.gravmag.startingModel.grav.file")
        call read_filename(10, gpar%model_files(3))
        call print_arg(myrank, parname, gpar%model_files(3))

      case("forward.gravmag.startingModel.mag.file")
        call read_filename(10, mpar%model_files(3))
        call print_arg(myrank, parname, mpar%model_files(3))

      ! MAGNETIC FIELD constants ---------------------------

      case("forward.magneticField.inclination")
        read(10, 1) mpar%mi
        call print_arg(myrank, parname, mpar%mi)

      case("forward.magneticField.declination")
        read(10, 1) mpar%md
        call print_arg(myrank, parname, mpar%md)

      case("forward.magneticField.ambient.inclination")
        read(10, 1) mpar%fi
        call print_arg(myrank, parname, mpar%fi)

      case("forward.magneticField.ambient.declination")
        read(10, 1) mpar%fd
        call print_arg(myrank, parname, mpar%fd)

      case("forward.magneticField.ambient.intensity_nT")
        read(10, 1) mpar%intensity
        call print_arg(myrank, parname, mpar%intensity)

      case("forward.magneticField.XaxisDeclination")
        read(10, 1) mpar%theta
        call print_arg(myrank, parname, mpar%theta)

      ! DEPTH WEIGHTING parameters -------------------------

      case("forward.gravmag.depthWeighting.type")
        read(10, 2) gpar%depth_weighting_type
        call print_arg(myrank, parname, gpar%depth_weighting_type)
        mpar%depth_weighting_type = gpar%depth_weighting_type

      case("forward.gravmag.depthWeighting.powerWeight.grav.beta")
        read(10, 1) gpar%beta
        call print_arg(myrank, parname, gpar%beta)

      case("forward.gravmag.depthWeighting.powerWeight.grav.Z0")
        read(10, 1) gpar%Z0
        call print_arg(myrank, parname, gpar%Z0)

      case("forward.gravmag.depthWeighting.powerWeight.mag.beta")
        read(10, 1) mpar%beta
        call print_arg(myrank, parname, mpar%beta)

      case("forward.gravmag.depthWeighting.powerWeight.mag.Z0")
        read(10, 1) mpar%Z0
        call print_arg(myrank, parname, mpar%Z0)

      ! MATRIX COMPRESSION parameters ----------------------

      case("forward.gravmag.matrixCompression.distanceThreshold")
        read(10, 1) gpar%distance_threshold
        call print_arg(myrank, parname, gpar%distance_threshold)
        mpar%distance_threshold = gpar%distance_threshold

      case("forward.gravmag.matrixCompression.compressionRate")
        read(10, 1) gpar%compression_rate
        call print_arg(myrank, parname, gpar%compression_rate)
        mpar%compression_rate = gpar%compression_rate

      ! INVERSION parameters -------------------------------

      case("inversion.nMajorIterations")
        read(10, 2) ipar%ninversions
        call print_arg(myrank, parname, ipar%ninversions)

      case("inversion.nMinorIterations")
        read(10, 2) ipar%niter
        call print_arg(myrank, parname, ipar%niter)

      case("inversion.minResidual")
        read(10, 1) ipar%rmin
        call print_arg(myrank, parname, ipar%rmin)

      case("inversion.solver")
        read(10, 2) ipar%method
        call print_arg(myrank, parname, ipar%method)

      case("inversion.softThresholdL1")
        read(10, 1) ipar%gamma
        call print_arg(myrank, parname, ipar%gamma)

      ! MODEL DAMPING (m - m_prior) ------------------------

      case("inversion.modelDamping.weightProblem1")
        read(10, 1) ipar%alpha(1)
        call print_arg(myrank, parname, ipar%alpha(1))

      case("inversion.modelDamping.weightProblem2")
        read(10, 1) ipar%alpha(2)
        call print_arg(myrank, parname, ipar%alpha(2))

      case("inversion.modelDamping.normPower")
        read(10, 1) ipar%norm_power
        call print_arg(myrank, parname, ipar%norm_power)

      ! JOINT INVERSION parameters -------------------------------

      case("inversion.joint.problemWeight1")
        read(10, 1) ipar%problem_weight(1)
        call print_arg(myrank, parname, ipar%problem_weight(1))

      case("inversion.joint.problemWeight2")
        read(10, 1) ipar%problem_weight(2)
        call print_arg(myrank, parname, ipar%problem_weight(2))

      case("inversion.joint.columnWeightMultiplier1")
        read(10, 1) ipar%column_weight_multiplier(1)
        call print_arg(myrank, parname, ipar%column_weight_multiplier(1))

      case("inversion.joint.columnWeightMultiplier2")
        read(10, 1) ipar%column_weight_multiplier(2)
        call print_arg(myrank, parname, ipar%column_weight_multiplier(2))

      case("inversion.joint.nIterSingle1")
        read(10, 2) ipar%niter_single(1)
        call print_arg(myrank, parname, ipar%niter_single(1))

      case("inversion.joint.nIterSingle2")
        read(10, 2) ipar%niter_single(2)
        call print_arg(myrank, parname, ipar%niter_single(2))

      ! DAMPING-GRADIENT constraints -------------------------------

      case("inversion.dampingGradient.weightType")
        read(10, 2) ipar%damp_grad_weight_type
        call print_arg(myrank, parname, ipar%damp_grad_weight_type)

      case("inversion.dampingGradient.weightProblem1")
        read(10, 1) ipar%beta(1)
        call print_arg(myrank, parname, ipar%beta(1))

      case("inversion.dampingGradient.weightProblem2")
        read(10, 1) ipar%beta(2)
        call print_arg(myrank, parname, ipar%beta(2))

      ! CROSS-GRADIENT constraints ---------------------------------

      case("inversion.crossGradient.weight")
        read(10, 1) ipar%cross_grad_weight
        call print_arg(myrank, parname, ipar%cross_grad_weight)

      case("inversion.crossGradient.nIterMethodOfWeight")
        read(10, 2) ipar%method_of_weights_niter
        call print_arg(myrank, parname, ipar%method_of_weights_niter)

      case("inversion.crossGradient.derivativeType")
        read(10, 2) ipar%derivative_type
        call print_arg(myrank, parname, ipar%derivative_type)

      ! CLUSTERING constraints ---------------------------------------

      case("inversion.clustering.weightProblem1")
        read(10, 1) ipar%clustering_weight_glob(1)
        call print_arg(myrank, parname, ipar%clustering_weight_glob(1))

      case("inversion.clustering.weightProblem2")
        read(10, 1) ipar%clustering_weight_glob(2)
        call print_arg(myrank, parname, ipar%clustering_weight_glob(2))

      case("inversion.clustering.nClusters")
        read(10, 2) ipar%nclusters
        call print_arg(myrank, parname, ipar%nclusters)

      case("inversion.clustering.mixtureFile")
        call read_filename(10, ipar%mixture_file)
        call print_arg(myrank, parname, ipar%mixture_file)

      case("inversion.clustering.cellWeightsFile")
        call read_filename(10, ipar%cell_weights_file)
        call print_arg(myrank, parname, ipar%cell_weights_file)

      case("inversion.clustering.optimizationType")
        read(10, 2) ipar%clustering_opt_type
        call print_arg(myrank, parname, ipar%clustering_opt_type)

      case("inversion.clustering.constraintsType")
        read(10, 2) ipar%clustering_constraints_type
        call print_arg(myrank, parname, ipar%clustering_constraints_type)

      ! ADMM constraints --------------------------------------------

      case("inversion.admm.enableADMM")
        read(10, 2) ipar%admm_type
        call print_arg(myrank, parname, ipar%admm_type)

      case("inversion.admm.nLithologies")
        read(10, 2) ipar%nlithos
        call print_arg(myrank, parname, ipar%nlithos)

      case("inversion.admm.localBoundsFile1")
        call read_filename(10, ipar%bounds_ADMM_file(1))
        call print_arg(myrank, parname, ipar%bounds_ADMM_file(1))

      case("inversion.admm.localBoundsFile2")
        call read_filename(10, ipar%bounds_ADMM_file(2))
        call print_arg(myrank, parname, ipar%bounds_ADMM_file(2))

      case("inversion.admm.weightProblem1")
        read(10, 1) ipar%rho_ADMM(1)
        call print_arg(myrank, parname, ipar%rho_ADMM(1))

      case("inversion.admm.weightProblem2")
        read(10, 1) ipar%rho_ADMM(2)
        call print_arg(myrank, parname, ipar%rho_ADMM(2))

      case default
        read(10, 3, iostat=ios) line
    end select
  enddo

  print *, "Finished reading the file."

! Format to read a floating-point value.
 1 format(f16.8)

! Format to read an integer value.
 2 format(i8)

! Format to read a string.
 3 format(a)

 end subroutine read_parfile2

!===================================================================================
! Read input parameters from Parfile.
!===================================================================================
subroutine read_parfile(epar, gpar, mpar, ipar, myrank)
  integer, intent(in) :: myrank

  type(t_parameters_ect), intent(out) :: epar
  type(t_parameters_grav), intent(out) :: gpar
  type(t_parameters_mag), intent(out) :: mpar
  type(t_parameters_inversion), intent(out) :: ipar

  integer :: itmp

  ! This is junk in order to ignore the variable name at the beginning of the line.
  ! This ignores exactly 40 characters.
  character(len=40) :: junk
  character(len=256) :: parfile_name
  character(len=128) :: dum

  ! The name of the Parfile can be passed in via the command line,
  ! if no argument is given, the default value is used.
  call get_command_argument(2,parfile_name)
  if (len_trim(parfile_name) == 0) parfile_name = "parfiles/Parfile_MASTER.txt"

  open(unit=10,file=parfile_name,status='old',iostat=itmp,action='read')
  if (itmp /= 0) call exit_MPI("Parfile """ // trim(parfile_name) // """ cannot be opened!",myrank,15)

  ! GLOBAL -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,3) junk, path_output
  if (myrank == 0) print *,junk,trim(path_output)

  ! DIMENSIONS -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk,epar%dims%nr
  if (myrank == 0) print *,junk,epar%dims%nr

  read(10,2) junk,epar%dims%ntheta
  if (epar%dims%ntheta == 0) epar%dims%ntheta = epar%dims%nr
  if (myrank == 0) print *,junk,epar%dims%ntheta

  read(10,2) junk,epar%dims%nz
  if (epar%dims%nz == 0) epar%dims%nz = epar%dims%nr
  if (myrank == 0) print *,junk,epar%dims%nz

  ! GEOMETRY -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk,epar%nel
  if (myrank == 0) print *,junk,epar%nel

  read(10,2) junk,epar%nrings
  if (myrank == 0) print *,junk,epar%nrings

  read(10,2) junk,epar%dims%kguards
  if (epar%dims%kguards == 0) epar%dims%kguards = epar%dims%nz/4
  if (myrank == 0) print *,junk,epar%dims%kguards

  read(10,2) junk,epar%ifixed_elecgeo
  if (myrank == 0) print *,junk,epar%ifixed_elecgeo

  read(10,2) junk,epar%irefine
  if (myrank == 0) print *,junk,epar%irefine

  read(10,1) junk,epar%sens%radiusin
  if (myrank == 0) print *,junk,epar%sens%radiusin

  read(10,1) junk,epar%sens%radiusout
  if (myrank == 0) print *,junk,epar%sens%radiusout

  read(10,1) junk,epar%sens%radiusoutout
  if (myrank == 0) print *,junk,epar%sens%radiusoutout

  read(10,1) junk,epar%sens%heicyl
  if (myrank == 0) print *,junk,epar%sens%heicyl

  read(10,1) junk,epar%sens%space_elec_guards
  if (myrank == 0) print *,junk,epar%sens%space_elec_guards

  read(10,1) junk,epar%sens%space_electrodes
  if (myrank == 0) print *,junk,epar%sens%space_electrodes

  ! MODEL -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk,epar%num_bubbles
  if (myrank == 0) print *,junk,epar%num_bubbles

  read(10,3) junk, epar%filename_bubbles
  if (myrank == 0) print *, junk, trim(epar%filename_bubbles)

  read(10,1) junk,epar%permit0
  if (myrank == 0) print *,junk,epar%permit0

  read(10,1) junk,epar%permit_air
  if (myrank == 0) print *,junk,epar%permit_air

  read(10,1) junk,epar%permit_isolated_tube
  if (myrank == 0) print *,junk,epar%permit_isolated_tube

  read(10,1) junk,epar%permit_oil
  if (myrank == 0) print *,junk,epar%permit_oil

  ! SOLVER parameters -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  epar%linear_solver = LINSOLV_PCG

  read(10,2) junk,epar%iprecond
  if (myrank == 0) print *,junk,epar%iprecond

  read(10,1) junk,epar%omega1
  if (myrank == 0) print *,junk, epar%omega1

  read(10,2) junk,epar%itypenorm
  if (epar%itypenorm == 1) then
    epar%itypenorm = NORM_L2
  else if (epar%itypenorm == 2) then
    epar%itypenorm = NORM_MAX
  else
   call exit_MPI("wrong setting for type of norm, must be 1 for L2- or 2 for max-norm; exiting...",myrank,17)
  endif
  if (myrank == 0) print *,junk,epar%itypenorm

  read(10,2) junk,epar%itmax
  if (myrank == 0) print *,junk,epar%itmax

  read(10,2) junk,epar%output_frequency
  if (myrank == 0) print *,junk,epar%output_frequency

  read(10,1) junk,epar%tol
  if (myrank == 0) print *,junk,epar%tol

  ! MULTIGRID parameters -------------------------------

  ! Removed multigrid, keep this not to change a lot of code.
  epar%ilevel_coarse = 1
  epar%coarse_solver = LINSOLV_PCG

  ! GRAVITY / MAGNETISM parameters -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,'(a)', advance='NO') junk
  read(10,*) gpar%nx, gpar%ny, gpar%nz
  if (myrank == 0) print *, junk, gpar%nx, gpar%ny, gpar%nz

  mpar%nx = gpar%nx
  mpar%ny = gpar%ny
  mpar%nz = gpar%nz

  read(10,3) junk, gpar%model_files(1)
  if (myrank == 0) print *, junk, trim(gpar%model_files(1))

  read(10,3) junk, mpar%model_files(1)
  if (myrank == 0) print *, junk, trim(mpar%model_files(1))

  read(10,2) junk, gpar%depth_weighting_type
  if (myrank == 0) print *, junk, gpar%depth_weighting_type

  mpar%depth_weighting_type = gpar%depth_weighting_type

  ! GRAV / MAG DATA parameters -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk, gpar%ndata
  if (myrank == 0) print *, junk, gpar%ndata

  read(10,2) junk, mpar%ndata
  if (myrank == 0) print *, junk, mpar%ndata

  read(10,3) junk, gpar%data_grid_file
  if (myrank == 0) print *, junk, trim(gpar%data_grid_file)

  read(10,3) junk, mpar%data_grid_file
  if (myrank == 0) print *, junk, trim(mpar%data_grid_file)

  read(10,3) junk, gpar%data_file
  if (myrank == 0) print *, junk, trim(gpar%data_file)

  read(10,3) junk, mpar%data_file
  if (myrank == 0) print *, junk, trim(mpar%data_file)

  read(10,2) junk, gpar%calc_data_directly
  if (myrank == 0) print *, junk, gpar%calc_data_directly

  mpar%calc_data_directly = gpar%calc_data_directly

  ! PRIOR MODEL -----------------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk, gpar%prior_model_type
  if (myrank == 0) print *, junk, gpar%prior_model_type
  
  mpar%prior_model_type = gpar%prior_model_type
  
  read(10,2) junk, gpar%number_prior_models
  if (myrank == 0) print *, junk, gpar%number_prior_models

  mpar%number_prior_models = gpar%number_prior_models

  read(10,1) junk, gpar%prior_model_val
  if (myrank == 0) print *, junk, gpar%prior_model_val

  read(10,1) junk, mpar%prior_model_val
  if (myrank == 0) print *, junk, mpar%prior_model_val

  read(10,3) junk, gpar%model_files(2)
  if (myrank == 0) print *, junk, trim(gpar%model_files(2))

  read(10,3) junk, mpar%model_files(2)
  if (myrank == 0) print *, junk, trim(mpar%model_files(2))

  ! STARTING MODEL ---------------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk, gpar%start_model_type
  if (myrank == 0) print *, junk, gpar%start_model_type

  mpar%start_model_type = gpar%start_model_type

  read(10,1) junk, gpar%start_model_val
  if (myrank == 0) print *, junk, gpar%start_model_val

  read(10,1) junk, mpar%start_model_val
  if (myrank == 0) print *, junk, mpar%start_model_val

  read(10,3) junk, gpar%model_files(3)
  if (myrank == 0) print *, junk, trim(gpar%model_files(3))

  read(10,3) junk, mpar%model_files(3)
  if (myrank == 0) print *, junk, trim(mpar%model_files(3))

  ! MAGNETIC constants -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, mpar%mi
  if (myrank == 0) print *, junk, mpar%mi

  read(10,1) junk, mpar%md
  if (myrank == 0) print *, junk, mpar%md

  read(10,1) junk, mpar%fi
  if (myrank == 0) print *, junk, mpar%fi

  read(10,1) junk, mpar%fd
  if (myrank == 0) print *, junk, mpar%fd

  read(10,1) junk, mpar%intensity
  if (myrank == 0) print *, junk, mpar%intensity

  read(10,1) junk, mpar%theta
  if (myrank == 0) print *, junk, mpar%theta

  read(10,1) junk, mpar%beta
  if (myrank == 0) print *, junk, mpar%beta

  read(10,1) junk, mpar%Z0
  if (myrank == 0) print *, junk, mpar%Z0

  ! GRAVITY constants -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, gpar%beta
  if (myrank == 0) print *, junk, gpar%beta

  read(10,1) junk, gpar%Z0
  if (myrank == 0) print *, junk, gpar%Z0

  ! MATRIX COMPRESSION parameters -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, gpar%distance_threshold
  if (myrank == 0) print *, junk, gpar%distance_threshold

  read(10,1) junk, gpar%compression_rate
  if (myrank == 0) print *, junk, gpar%compression_rate

  mpar%distance_threshold = gpar%distance_threshold
  mpar%compression_rate = gpar%compression_rate

  ! INVERSION parameters -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk, ipar%ninversions
  if (myrank == 0) print *, junk, ipar%ninversions

  read(10,2) junk, ipar%niter
  if (myrank == 0) print *, junk, ipar%niter

  read(10,1) junk, ipar%rmin
  if (myrank == 0) print *, junk, ipar%rmin

  read(10,2) junk, ipar%method
  if (myrank == 0) print *, junk, ipar%method

  read(10,1) junk, ipar%gamma
  if (myrank == 0) print *, junk, ipar%gamma

  ! MODEL DAMPING (m - m_prior) -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, ipar%alpha(1)
  if (myrank == 0) print *, junk, ipar%alpha(1)

  read(10,1) junk, ipar%alpha(2)
  if (myrank == 0) print *, junk, ipar%alpha(2)

  read(10,1) junk, ipar%norm_power
  if (myrank == 0) print *, junk, ipar%norm_power

  ! JOINT INVERSION parameters -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, ipar%problem_weight(1)
  if (myrank == 0) print *, junk, ipar%problem_weight(1)

  read(10,1) junk, ipar%problem_weight(2)
  if (myrank == 0) print *, junk, ipar%problem_weight(2)

  read(10,1) junk, ipar%column_weight_multiplier(1)
  if (myrank == 0) print *, junk, ipar%column_weight_multiplier(1)

  read(10,1) junk, ipar%column_weight_multiplier(2)
  if (myrank == 0) print *, junk, ipar%column_weight_multiplier(2)

  read(10,2) junk, ipar%niter_single(1)
  if (myrank == 0) print *, junk, ipar%niter_single(1)

  read(10,2) junk, ipar%niter_single(2)
  if (myrank == 0) print *, junk, ipar%niter_single(2)

  ! Damping-gradient constraints -------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk, ipar%damp_grad_weight_type
  if (myrank == 0) print *, junk, ipar%damp_grad_weight_type

  read(10,1) junk, ipar%beta(1)
  if (myrank == 0) print *, junk, ipar%beta(1)

  read(10,1) junk, ipar%beta(2)
  if (myrank == 0) print *, junk, ipar%beta(2)

  ! Cross-gradient constraints ---------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, ipar%cross_grad_weight
  if (myrank == 0) print *, junk, ipar%cross_grad_weight

  read(10,2) junk, ipar%method_of_weights_niter
  if (myrank == 0) print *, junk, ipar%method_of_weights_niter

  read(10,2) junk, ipar%derivative_type
  if (myrank == 0) print *, junk, ipar%derivative_type

  ! Clustering constraints -------------------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,1) junk, ipar%clustering_weight_glob(1)
  if (myrank == 0) print *, junk, ipar%clustering_weight_glob(1)

  read(10,1) junk, ipar%clustering_weight_glob(2)
  if (myrank == 0) print *, junk, ipar%clustering_weight_glob(2)

  read(10,2) junk, ipar%nclusters
  if (myrank == 0) print *, junk, ipar%nclusters

  read(10,3) junk, ipar%mixture_file
  if (myrank == 0) print *, junk, trim(ipar%mixture_file)

  read(10,3) junk, ipar%cell_weights_file
  if (myrank == 0) print *, junk, trim(ipar%cell_weights_file)

  read(10,2) junk, ipar%clustering_opt_type
  if (myrank == 0) print *, junk, ipar%clustering_opt_type

  read(10,2) junk, ipar%clustering_constraints_type
  if (myrank == 0) print *, junk, ipar%clustering_constraints_type

  ! ADMM constraints ------------------------------------------------

  read(10,'(a)') dum
  if (myrank == 0) print *, trim(dum)

  read(10,2) junk, ipar%admm_type
  if (myrank == 0) print *, junk, ipar%admm_type
  
  read(10,2) junk, ipar%nlithos
  if (myrank == 0) print *, junk, ipar%nlithos

  read(10,3) junk, ipar%bounds_ADMM_file(1)
  if (myrank == 0) print *, junk, trim(ipar%bounds_ADMM_file(1))
  
  read(10,3) junk, ipar%bounds_ADMM_file(2)
  if (myrank == 0) print *, junk, trim(ipar%bounds_ADMM_file(2))

  read(10,'(a)', advance='NO') junk
  read(10,*) ipar%rho_ADMM(1)
  if (myrank == 0) print *, junk, ipar%rho_ADMM(1)

  read(10,'(a)', advance='NO') junk
  read(10,*) ipar%rho_ADMM(2)
  if (myrank == 0) print *, junk, ipar%rho_ADMM(2)

  close(10)

  if (myrank == 0) print *, '**********************************************'

  ! Print out if we do this in double or single precision.
  if (myrank == 0) then
    if (CUSTOM_REAL == SIZE_DOUBLE) then
      print *, "precision                    = DOUBLE"
    else
      print *, "precision                    = SINGLE"
    endif
  endif

! Format to read a floating-point value.
 1 format(a, f16.8)

! Format to read an integer value.
 2 format(a, i8)

! Format to read a string.
 3 format(a, a)

end subroutine read_parfile

end module init_parameters
