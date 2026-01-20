
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

!========================================================================================
! Read the input code parameters from the Parfile, and initialize the code parameters.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!========================================================================================
module init_parameters

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use file_utils, only: create_directory
  use parameters_grav
  use parameters_mag
  use parameters_inversion
  use parallel_tools

  implicit none

  private

  public :: initialize_parameters
  public :: get_problem_type

  private :: get_parfile_path
  private :: read_parfile
  private :: set_default_parameters

  private :: print_arg_int
  private :: print_arg_dbl
  private :: print_arg_dblarr
  private :: print_arg_str
  private :: read_filename
  private :: copy_file

  ! For overloading the 'print_arg' function.
  interface print_arg
    module procedure print_arg_int, print_arg_dbl, print_arg_dblarr, print_arg_str
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
subroutine print_arg_dblarr(myrank, name, value)
  character(len=128), intent(in) :: name
  real(kind=CUSTOM_REAL), intent(in) :: value(:)
  integer, intent(in) :: myrank

  if (myrank == 0) print *, trim(name)//" =", value
end subroutine print_arg_dblarr
!-----------------------------------------------
subroutine print_arg_str(myrank, name, value)
  character(len=128), intent(in) :: name
  character(len=*), intent(in) :: value
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
! Get problem type.
!==========================================================================
subroutine get_problem_type(problem_type, myrank)
  integer, intent(in) :: myrank
  integer, intent(out) :: problem_type
  character(len=256) :: arg

  arg = ''

  if (command_argument_count() > 0) call get_command_argument(1, arg)

  if (arg == '-j' .or. arg == '-p') then
    problem_type = 1
  else
    call exit_MPI("Unknown problem type! arg ="//arg, myrank, 0)
  endif

end subroutine get_problem_type

!=======================================================================================
! Initialize parameters for forward and inverse problems.
!=======================================================================================
subroutine initialize_parameters(problem_type, gpar, mpar, ipar, myrank, nbproc)
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc
  type(t_parameters_grav), intent(out) :: gpar
  type(t_parameters_mag), intent(out) :: mpar
  type(t_parameters_inversion), intent(out) :: ipar

  character(len=256) :: parfile_path
  integer :: nelements, ierr

  if (myrank == 0) then
    ! Setting default parameter values.
    call set_default_parameters(gpar, mpar, ipar)

    ! Get the Parfile path.
    call get_parfile_path(parfile_path, myrank)

    ! Read the Parfile parameters.
    call read_parfile(parfile_path, gpar, mpar, ipar, myrank)

    ! Create the output directory. If it already exists there is no problem.
    call create_directory(trim(path_output))

    ! Copy the Parfile to the output folder.
    call copy_file(parfile_path, trim(path_output)//'/Parfile_copy.txt')

    ! Print out if we use double or single precision.
    if (CUSTOM_REAL == SIZE_DOUBLE) then
      print *, "precision = DOUBLE"
    else
      print *, "precision = SINGLE"
    endif
  endif

  ! Use barrier to do not mix the parameters output from the master CPU with other log messages.
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! All sanity checks passed so far, so we may broadcast the Parfile.

  ! MPI broadcast parameters.
  call MPI_Bcast(path_output, len(path_output), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  if (problem_type == 1) then
    call gpar%broadcast(myrank)
    call mpar%broadcast(myrank)
  endif

  call ipar%broadcast(myrank)

  !--------------------------------
  ! Some extra initializations.
  !--------------------------------
  if (problem_type == 1) then
    ! Gravity and magnetism problems.

    ! Inverse problem parameters.
    ipar%nx = gpar%nx
    ipar%ny = gpar%ny
    ipar%nz = gpar%nz

    ipar%ndata(1) = gpar%ndata
    ipar%ndata(2) = mpar%ndata

    ipar%ndata_components(1) = gpar%ndata_components
    ipar%ndata_components(2) = mpar%ndata_components

    ipar%nmodel_components = 1

    if (mpar%nmodel_components > 1) then
      if (ipar%problem_weight(1) /= 0.d0) then
        call exit_MPI("For the magnetisation inversion the gravity problem should be disabled!", myrank, 0)
      endif
      ipar%nmodel_components = mpar%nmodel_components
    endif

    ipar%nelements_total = ipar%nx * ipar%ny * ipar%nz

    ipar%compression_type = gpar%compression_type

    ! Define model splitting for parallelization.
    nelements = calculate_nelements_at_cpu(ipar%nelements_total, myrank, nbproc)

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
subroutine set_default_parameters(gpar, mpar, ipar)
  type(t_parameters_grav), intent(out) :: gpar
  type(t_parameters_mag), intent(out) :: mpar
  type(t_parameters_inversion), intent(out) :: ipar

  !-----------------------------------------------
  ! Define here the DEFAULT parameter values:
  !-----------------------------------------------

  ! GLOBAL parameters.
  path_output = "output/test/"
  gpar%data_units_mult = 1.d0
  mpar%data_units_mult = 1.d0
  gpar%model_units_mult = 1.d0
  mpar%model_units_mult = 1.d0
  gpar%z_axis_dir = 1
  mpar%z_axis_dir = 1
  gpar%vtk_model_label = "rho"
  mpar%vtk_model_label = "k"

  gpar%model_files = "NILL"
  mpar%model_files = "NILL"

  ! MODEL GRID parameters.
  gpar%nx = 0
  gpar%ny = 0
  gpar%nz = 0
  mpar%nx = 0
  mpar%ny = 0
  mpar%nz = 0
  gpar%model_grid_file = "NILL"
  mpar%model_grid_file = "NILL"
  gpar%nmodel_components = 1
  mpar%nmodel_components = 1

  ! DATA parameters.
  gpar%ndata = 0
  mpar%ndata = 0
  gpar%data_grid_file = "NILL"
  mpar%data_grid_file = "NILL"
  gpar%data_file = "NILL"
  mpar%data_file = "NILL"
  gpar%ndata_components = 1
  mpar%ndata_components = 1
  gpar%data_type = 1
  mpar%data_type = 1

  gpar%useSyntheticModelForDataValues = .false.
  mpar%useSyntheticModelForDataValues = .false.

  ! Data error.
  gpar%use_data_error = 0
  mpar%use_data_error = 0
  gpar%data_error_file = "NILL"
  mpar%data_error_file = "NILL"

  ! MAGNETIC FIELD constants.
  mpar%mi = 90.d0
  mpar%md = 0.d0
  mpar%intensity = 50000.d0
  mpar%theta = 0.d0

  ! DEPTH WEIGHTING parameters.
  gpar%depth_weighting_type = 2 ! 1-depth weighting, 2-distance weighting
  mpar%depth_weighting_type = 2
  gpar%depth_weighting_power = 2.0d0
  mpar%depth_weighting_power = 3.0d0
  gpar%depth_weighting_beta = 1.0d0
  mpar%depth_weighting_beta = 1.0d0
  gpar%Z0 = 0.d0
  mpar%Z0 = 0.d0

  gpar%apply_local_weight = 0
  mpar%apply_local_weight = 0
  gpar%local_weight_file = "NILL"
  mpar%local_weight_file = "NILL"

  ! SENSITIVITY KERNEL parameters.
  gpar%sensit_read = 0
  gpar%sensit_path = "SENSIT/"
  mpar%sensit_read = gpar%sensit_read
  mpar%sensit_path = gpar%sensit_path

  ! MATRIX COMPRESSION parameters.
  gpar%compression_type = 0
  mpar%compression_type = 0
  gpar%compression_rate = 1.d-1
  mpar%compression_rate = 1.d-1

  ! PRIOR MODEL parameters.
  gpar%prior_model_type = 1 ! 1-set value, 2-read from file.
  mpar%prior_model_type = 1
  gpar%number_prior_models = 1 ! Number of prior models, for the model type = 2
  mpar%number_prior_models = 1
  gpar%prior_model_val = 0.d0
  mpar%prior_model_val = 0.d0

  ! STARTING MODEL parameters.
  gpar%start_model_type = 1 ! 1-set value, 2-read from file.
  mpar%start_model_type = 1
  gpar%start_model_val = 0.d0
  mpar%start_model_val = 0.d0

  ! INVERSION parameters.
  ipar%ninversions = 10
  ipar%niter = 100
  ipar%target_misfit = 0.d0
  ipar%write_model_niter = 0
  ipar%rmin = 1.d-13
  ipar%method = 1 ! LSQR = 1
  ipar%gamma = 0. ! soft threshold ("L1-norm", no=0.)

  ! MODEL DAMPING (m - m_prior).
  ipar%alpha(1) = 1.d-11
  ipar%alpha(2) = 1.d-8
  ipar%norm_power = 2.d0

  ipar%apply_local_damping_weight = 0
  ipar%damping_weight_file = "NILL"

  ! JOINT INVERSION parameters.
  ipar%problem_weight(1) = 1.d0
  ipar%problem_weight(2) = 0.d0
  ipar%column_weight_multiplier(1) = 4.d+3
  ipar%column_weight_multiplier(2) = 1.d0

  ! ADMM constraints.
  ipar%admm_type = 0 ! 0-no admm, 1-enable admm
  ipar%admm_bound_type = 1 ! 1-global, 2-local from file
  ipar%nlithos = 1
  ipar%bounds_ADMM_file(1) = "NILL"
  ipar%bounds_ADMM_file(2) = "NILL"
  ipar%rho_ADMM(1) = 1.d-7
  ipar%rho_ADMM(2) = 1.d+5

  ipar%data_cost_threshold_ADMM = 1.e-4
  ipar%weight_multiplier_ADMM = 1.d0
  ipar%max_weight_ADMM = 1.e+10

  ! DAMPING-GRADIENT constraints.
  ipar%damp_grad_weight_type = 1 ! 1-global, 2-local
  ipar%beta(1) = 0.d0
  ipar%beta(2) = 0.d0
  ipar%damping_gradient_file = "NILL"

  ! CROSS-GRADIENT constraints.
  ipar%cross_grad_weight = 0.d0
  ipar%derivative_type = 1 ! 1-fwd, 2-cent
  ipar%keep_model_constant = 0
  ipar%vec_field_type = 0
  ipar%vec_field_file = "NILL"

  ! CLUSTERING constraints.
  ipar%clustering_weight_glob(1) = 0.d0
  ipar%clustering_weight_glob(2) = 0.d0
  ipar%nclusters = 4
  ipar%mixture_file = "NILL"
  ipar%cell_weights_file = "NILL"
  ipar%clustering_opt_type = 2 ! 1-normal, 2-log
  ipar%clustering_constraints_type = 2 ! 1-global, 2-local

end subroutine set_default_parameters

!========================================================================================
! Retrieve the Parfile path from the command line.
!========================================================================================
subroutine get_parfile_path(parfile_path, myrank)
  integer, intent(in) :: myrank
  character(len=*), intent(out) :: parfile_path

  ! The name of the Parfile should be passed in via the command line.
  call get_command_argument(2, parfile_path)
  if (len_trim(parfile_path) == 0) then
    call exit_MPI("No Parfile supplied!", myrank, 0)
    stop
  endif
end subroutine get_parfile_path

!===================================================================================
! Copy a file.
!===================================================================================
subroutine copy_file(path_from, path_to)
  character(len=*) :: path_from, path_to

  call execute_command_line('cp "' // trim(path_from) // '" "' // trim(path_to) // '"')
end subroutine copy_file

!===================================================================================
! Read input parameters from Parfile.
!===================================================================================
subroutine read_parfile(parfile_path, gpar, mpar, ipar, myrank)
  character(len=*), intent(in) :: parfile_path
  integer, intent(in) :: myrank
  type(t_parameters_grav), intent(inout) :: gpar
  type(t_parameters_mag), intent(inout) :: mpar
  type(t_parameters_inversion), intent(inout) :: ipar

  character(len=1) :: ch
  character(len=128) :: parname
  character(len=256) :: line
  character(len=256) :: parfile_description
  integer :: itmp
  integer :: symbol_index, i
  integer :: ios
  logical :: global_bounds_defined
  integer :: useSyntheticModelForDataValuesInt

  open(unit=10, file=trim(parfile_path), status='old', iostat=itmp, action='read')
  if (itmp /= 0) call exit_MPI("Parfile """ // trim(parfile_path) // """ cannot be opened!", myrank, 0)

  global_bounds_defined = .false.

  !---------------------------------------------------------------------------------
  ! Reading parameter values from Parfile.
  !---------------------------------------------------------------------------------
  do
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit

    ! Skip a comment line.
    if (line(1:1) == '#') cycle

    symbol_index = index(line, '=')
    parname = line(:symbol_index - 1)

    ! Skip this line.
    if (trim(parname) == '') cycle

    backspace(10)

    ! Reading file line until the end of '=' symbol.
    do i = 1, symbol_index
      read(10, '(a)', advance='NO') ch
    enddo

    select case(trim(parname))

      ! GLOBAL parameters -----------------------------------

      case("global.outputFolderPath")
        call read_filename(10, path_output)
        call print_arg(myrank, parname, path_output)

      case("global.description")
        read(10, '(a)') parfile_description
        call print_arg(myrank, parname, parfile_description)

      case("global.grav.dataUnitsMultiplier")
        read(10, *) gpar%data_units_mult
        call print_arg(myrank, parname, gpar%data_units_mult)

      case("global.magn.dataUnitsMultiplier")
        read(10, *) mpar%data_units_mult
        call print_arg(myrank, parname, mpar%data_units_mult)

      case("global.grav.modelUnitsMultiplier")
        read(10, *) gpar%model_units_mult
        call print_arg(myrank, parname, gpar%model_units_mult)

      case("global.magn.modelUnitsMultiplier")
        read(10, *) mpar%model_units_mult
        call print_arg(myrank, parname, mpar%model_units_mult)

      case("global.zAxisDirection")
        read(10, *) gpar%z_axis_dir
        call print_arg(myrank, parname, gpar%z_axis_dir)
        mpar%z_axis_dir = gpar%z_axis_dir

      ! MODEL GRID parameters -------------------------------

      case("modelGrid.size")
        read(10, *) gpar%nx, gpar%ny, gpar%nz
        if (myrank == 0) print *, trim(parname)//" =", gpar%nx, gpar%ny, gpar%nz
        mpar%nx = gpar%nx
        mpar%ny = gpar%ny
        mpar%nz = gpar%nz

      case("modelGrid.grav.file")
        call read_filename(10, gpar%model_grid_file)
        call print_arg(myrank, parname, gpar%model_grid_file)

      case("modelGrid.magn.file")
        call read_filename(10, mpar%model_grid_file)
        call print_arg(myrank, parname, mpar%model_grid_file)

      case("modelGrid.magn.nModelComponents")
        read(10, *) mpar%nmodel_components
        call print_arg(myrank, parname, mpar%nmodel_components)

      ! DATA parameters -------------------------------------

      case("forward.data.grav.nData")
        read(10, *) gpar%ndata
        call print_arg(myrank, parname, gpar%ndata)

      case("forward.data.magn.nData")
        read(10, *) mpar%ndata
        call print_arg(myrank, parname, mpar%ndata)

      case("forward.data.grav.dataGridFile")
        call read_filename(10, gpar%data_grid_file)
        call print_arg(myrank, parname, gpar%data_grid_file)

      case("forward.data.magn.dataGridFile")
        call read_filename(10, mpar%data_grid_file)
        call print_arg(myrank, parname, mpar%data_grid_file)

      case("forward.data.grav.nDataComponents")
        read(10, *) gpar%ndata_components
        call print_arg(myrank, parname, gpar%ndata_components)

      case("forward.data.magn.nDataComponents")
        read(10, *) mpar%ndata_components
        call print_arg(myrank, parname, mpar%ndata_components)

      case("forward.data.grav.type")
        read(10, *) gpar%data_type
        call print_arg(myrank, parname, gpar%data_type)

      case("forward.data.grav.useError")
        read(10, *) gpar%use_data_error
        call print_arg(myrank, parname, gpar%use_data_error)

      case("forward.data.magn.useError")
        read(10, *) mpar%use_data_error
        call print_arg(myrank, parname, mpar%use_data_error)

      case("forward.data.grav.errorFile")
        call read_filename(10, gpar%data_error_file)
        call print_arg(myrank, parname, gpar%data_error_file)

      case("forward.data.magn.errorFile")
        call read_filename(10, mpar%data_error_file)
        call print_arg(myrank, parname, mpar%data_error_file)

      case("forward.data.grav.useSyntheticModelForDataValues")
        read(10, *) useSyntheticModelForDataValuesInt
        call print_arg(myrank, parname, useSyntheticModelForDataValuesInt)

        if (useSyntheticModelForDataValuesInt == 0) then
          gpar%useSyntheticModelForDataValues = .false.
        else
          gpar%useSyntheticModelForDataValues = .true.
        endif

      case("forward.data.magn.useSyntheticModelForDataValues")
        read(10, *) useSyntheticModelForDataValuesInt
        call print_arg(myrank, parname, useSyntheticModelForDataValuesInt)

        if (useSyntheticModelForDataValuesInt == 0) then
          mpar%useSyntheticModelForDataValues = .false.
        else
          mpar%useSyntheticModelForDataValues = .true.
        endif

      case("forward.data.grav.syntheticModelFile")
        call read_filename(10, gpar%model_files(1))
        call print_arg(myrank, parname, gpar%model_files(1))

      case("forward.data.magn.syntheticModelFile")
        call read_filename(10, mpar%model_files(1))
        call print_arg(myrank, parname, mpar%model_files(1))

      ! MAGNETIC FIELD constants ---------------------------

      case("forward.magneticField.inclination")
        read(10, *) mpar%mi
        call print_arg(myrank, parname, mpar%mi)

      case("forward.magneticField.declination")
        read(10, *) mpar%md
        call print_arg(myrank, parname, mpar%md)

      case("forward.magneticField.intensity_nT")
        read(10, *) mpar%intensity
        call print_arg(myrank, parname, mpar%intensity)

      case("forward.magneticField.XaxisDeclination")
        read(10, *) mpar%theta
        call print_arg(myrank, parname, mpar%theta)

      ! DEPTH WEIGHTING parameters -------------------------

      case("forward.depthWeighting.type")
        read(10, *) gpar%depth_weighting_type
        call print_arg(myrank, parname, gpar%depth_weighting_type)
        mpar%depth_weighting_type = gpar%depth_weighting_type

      case("forward.depthWeighting.grav.power")
        read(10, *) gpar%depth_weighting_power
        call print_arg(myrank, parname, gpar%depth_weighting_power)

      case("forward.depthWeighting.grav.beta")
        read(10, *) gpar%depth_weighting_beta
        call print_arg(myrank, parname, gpar%depth_weighting_beta)

      case("forward.depthWeighting.grav.Z0")
        read(10, *) gpar%Z0
        call print_arg(myrank, parname, gpar%Z0)

      case("forward.depthWeighting.magn.power")
        read(10, *) mpar%depth_weighting_power
        call print_arg(myrank, parname, mpar%depth_weighting_power)

      case("forward.depthWeighting.magn.beta")
        read(10, *) mpar%depth_weighting_beta
        call print_arg(myrank, parname, mpar%depth_weighting_beta)

      case("forward.depthWeighting.magn.Z0")
        read(10, *) mpar%Z0
        call print_arg(myrank, parname, mpar%Z0)

      case("forward.depthWeighting.applyLocalWeight")
        read(10, *) gpar%apply_local_weight
        call print_arg(myrank, parname, gpar%apply_local_weight)
        mpar%apply_local_weight = gpar%apply_local_weight

      case("forward.depthWeighting.grav.file")
        call read_filename(10, gpar%local_weight_file)
        call print_arg(myrank, parname, gpar%local_weight_file)

      case("forward.depthWeighting.magn.file")
        call read_filename(10, mpar%local_weight_file)
        call print_arg(myrank, parname, mpar%local_weight_file)

      ! SENSITIVITY KERNEL parameters -----------------------

      case("sensit.readFromFiles")
        read(10, *) gpar%sensit_read
        call print_arg(myrank, parname, gpar%sensit_read)
        mpar%sensit_read = gpar%sensit_read

      case("sensit.folderPath")
        call read_filename(10, gpar%sensit_path)
        call print_arg(myrank, parname, gpar%sensit_path)
        mpar%sensit_path = gpar%sensit_path

      ! MATRIX COMPRESSION parameters ----------------------

      case("forward.matrixCompression.type")
        read(10, *) gpar%compression_type
        call print_arg(myrank, parname, gpar%compression_type)
        mpar%compression_type = gpar%compression_type

      case("forward.matrixCompression.rate")
        read(10, *) gpar%compression_rate
        call print_arg(myrank, parname, gpar%compression_rate)
        mpar%compression_rate = gpar%compression_rate

      ! PRIOR MODEL -----------------------------------------

      case("inversion.priorModel.type")
        read(10, *) gpar%prior_model_type
        call print_arg(myrank, parname, gpar%prior_model_type)
        mpar%prior_model_type = gpar%prior_model_type

      case("inversion.priorModel.nModels")
        read(10, *) gpar%number_prior_models
        call print_arg(myrank, parname, gpar%number_prior_models)
        mpar%number_prior_models = gpar%number_prior_models

      case("inversion.priorModel.grav.value")
        read(10, *) gpar%prior_model_val
        call print_arg(myrank, parname, gpar%prior_model_val)

      case("inversion.priorModel.magn.value")
        read(10, *) mpar%prior_model_val
        call print_arg(myrank, parname, mpar%prior_model_val)

      case("inversion.priorModel.grav.file")
        call read_filename(10, gpar%model_files(2))
        call print_arg(myrank, parname, gpar%model_files(2))

      case("inversion.priorModel.magn.file")
        call read_filename(10, mpar%model_files(2))
        call print_arg(myrank, parname, mpar%model_files(2))

      ! STARTING MODEL -------------------------------------

      case("inversion.startingModel.type")
        read(10, *) gpar%start_model_type
        call print_arg(myrank, parname, gpar%start_model_type)
        mpar%start_model_type = gpar%start_model_type

      case("inversion.startingModel.grav.value")
        read(10, *) gpar%start_model_val
        call print_arg(myrank, parname, gpar%start_model_val)

      case("inversion.startingModel.magn.value")
        read(10, *) mpar%start_model_val
        call print_arg(myrank, parname, mpar%start_model_val)

      case("inversion.startingModel.grav.file")
        call read_filename(10, gpar%model_files(3))
        call print_arg(myrank, parname, gpar%model_files(3))

      case("inversion.startingModel.magn.file")
        call read_filename(10, mpar%model_files(3))
        call print_arg(myrank, parname, mpar%model_files(3))

      ! INVERSION parameters -------------------------------

      case("inversion.nMajorIterations")
        read(10, *) ipar%ninversions
        call print_arg(myrank, parname, ipar%ninversions)

      case("inversion.nMinorIterations")
        read(10, *) ipar%niter
        call print_arg(myrank, parname, ipar%niter)

      case("inversion.targetMisfit")
        read(10, *) ipar%target_misfit
        call print_arg(myrank, parname, ipar%target_misfit)

      case("inversion.writeModelEveryNiter")
        read(10, *) ipar%write_model_niter
        call print_arg(myrank, parname, ipar%write_model_niter)

      case("inversion.minResidual")
        read(10, *) ipar%rmin
        call print_arg(myrank, parname, ipar%rmin)

      case("inversion.solver")
        read(10, *) ipar%method
        call print_arg(myrank, parname, ipar%method)

      case("inversion.softThresholdL1")
        read(10, *) ipar%gamma
        call print_arg(myrank, parname, ipar%gamma)

      ! MODEL DAMPING (m - m_prior) ------------------------

      case("inversion.modelDamping.grav.weight")
        read(10, *) ipar%alpha(1)
        call print_arg(myrank, parname, ipar%alpha(1))

      case("inversion.modelDamping.magn.weight")
        read(10, *) ipar%alpha(2)
        call print_arg(myrank, parname, ipar%alpha(2))

      case("inversion.modelDamping.normPower")
        read(10, *) ipar%norm_power
        call print_arg(myrank, parname, ipar%norm_power)

      case("inversion.modelDamping.applyLocalWeight")
        read(10, *) ipar%apply_local_damping_weight
        call print_arg(myrank, parname, ipar%apply_local_damping_weight)

      case("inversion.modelDamping.grav.file")
        call read_filename(10, ipar%damping_weight_file(1))
        call print_arg(myrank, parname, ipar%damping_weight_file(1))

      case("inversion.modelDamping.magn.file")
        call read_filename(10, ipar%damping_weight_file(2))
        call print_arg(myrank, parname, ipar%damping_weight_file(2))

      ! JOINT INVERSION parameters -------------------------------

      case("inversion.joint.grav.problemWeight")
        read(10, *) ipar%problem_weight(1)
        call print_arg(myrank, parname, ipar%problem_weight(1))

      case("inversion.joint.magn.problemWeight")
        read(10, *) ipar%problem_weight(2)
        call print_arg(myrank, parname, ipar%problem_weight(2))

      case("inversion.joint.grav.columnWeightMultiplier")
        read(10, *) ipar%column_weight_multiplier(1)
        call print_arg(myrank, parname, ipar%column_weight_multiplier(1))

      case("inversion.joint.magn.columnWeightMultiplier")
        read(10, *) ipar%column_weight_multiplier(2)
        call print_arg(myrank, parname, ipar%column_weight_multiplier(2))

      ! ADMM constraints --------------------------------------------

      case("inversion.admm.enableADMM")
        read(10, *) ipar%admm_type
        call print_arg(myrank, parname, ipar%admm_type)

      case("inversion.admm.boundType")
        read(10, *) ipar%admm_bound_type
        call print_arg(myrank, parname, ipar%admm_bound_type)

      case("inversion.admm.nLithologies")
        read(10, *) ipar%nlithos
        call print_arg(myrank, parname, ipar%nlithos)

      case("inversion.admm.grav.bounds")
        if (ipar%admm_type > 0 .and. ipar%admm_bound_type == 1) then
        ! Read these only when the global bound type is used.
          allocate(ipar%admm_bounds(1)%val(2 * ipar%nlithos), source=0._CUSTOM_REAL)
          read(10, *) ipar%admm_bounds(1)%val
          call print_arg(myrank, parname, ipar%admm_bounds(1)%val)
          global_bounds_defined = .true.
        endif

      case("inversion.admm.magn.bounds")
        if (ipar%admm_type > 0 .and. ipar%admm_bound_type == 1) then
        ! Read these only when the global bound type is used.
          allocate(ipar%admm_bounds(2)%val(2 * ipar%nlithos), source=0._CUSTOM_REAL)
          read(10, *) ipar%admm_bounds(2)%val
          call print_arg(myrank, parname, ipar%admm_bounds(2)%val)
          global_bounds_defined = .true.
        endif

      case("inversion.admm.grav.boundsFile")
        call read_filename(10, ipar%bounds_ADMM_file(1))
        call print_arg(myrank, parname, ipar%bounds_ADMM_file(1))

      case("inversion.admm.magn.boundsFile")
        call read_filename(10, ipar%bounds_ADMM_file(2))
        call print_arg(myrank, parname, ipar%bounds_ADMM_file(2))

      case("inversion.admm.grav.weight")
        read(10, *) ipar%rho_ADMM(1)
        call print_arg(myrank, parname, ipar%rho_ADMM(1))

      case("inversion.admm.magn.weight")
        read(10, *) ipar%rho_ADMM(2)
        call print_arg(myrank, parname, ipar%rho_ADMM(2))

      case("inversion.admm.dataCostThreshold")
        read(10, *) ipar%data_cost_threshold_ADMM
        call print_arg(myrank, parname, ipar%data_cost_threshold_ADMM)

      case("inversion.admm.weightMultiplier")
        read(10, *) ipar%weight_multiplier_ADMM
        call print_arg(myrank, parname, ipar%weight_multiplier_ADMM)

      case("inversion.admm.maxWeight")
        read(10, *) ipar%max_weight_ADMM
        call print_arg(myrank, parname, ipar%max_weight_ADMM)

      ! DAMPING-GRADIENT constraints -------------------------------

      case("inversion.dampingGradient.weightType")
        read(10, *) ipar%damp_grad_weight_type
        call print_arg(myrank, parname, ipar%damp_grad_weight_type)

      case("inversion.dampingGradient.grav.weight")
        read(10, *) ipar%beta(1)
        call print_arg(myrank, parname, ipar%beta(1))

      case("inversion.dampingGradient.magn.weight")
        read(10, *) ipar%beta(2)
        call print_arg(myrank, parname, ipar%beta(2))

      case("inversion.dampingGradient.grav.weightsFile")
        call read_filename(10, ipar%damping_gradient_file(1))
        call print_arg(myrank, parname, ipar%damping_gradient_file(1))

      case("inversion.dampingGradient.magn.weightsFile")
        call read_filename(10, ipar%damping_gradient_file(2))
        call print_arg(myrank, parname, ipar%damping_gradient_file(2))

      ! CROSS-GRADIENT constraints ---------------------------------

      case("inversion.crossGradient.weight")
        read(10, *) ipar%cross_grad_weight
        call print_arg(myrank, parname, ipar%cross_grad_weight)

      case("inversion.crossGradient.derivativeType")
        read(10, *) ipar%derivative_type
        call print_arg(myrank, parname, ipar%derivative_type)

      case("inversion.crossGradient.grav.keepModelConstant")
        read(10, *) ipar%keep_model_constant(1)
        call print_arg(myrank, parname, ipar%keep_model_constant(1))

      case("inversion.crossGradient.magn.keepModelConstant")
        read(10, *) ipar%keep_model_constant(2)
        call print_arg(myrank, parname, ipar%keep_model_constant(2))

      case("inversion.crossGradient.vectorFieldType")
        read(10, *) ipar%vec_field_type
        call print_arg(myrank, parname, ipar%vec_field_type)

      case("inversion.crossGradient.vectorFieldFile")
        call read_filename(10, ipar%vec_field_file)
        call print_arg(myrank, parname, ipar%vec_field_file)

      ! CLUSTERING constraints ---------------------------------------

      case("inversion.clustering.grav.weight")
        read(10, *) ipar%clustering_weight_glob(1)
        call print_arg(myrank, parname, ipar%clustering_weight_glob(1))

      case("inversion.clustering.magn.weight")
        read(10, *) ipar%clustering_weight_glob(2)
        call print_arg(myrank, parname, ipar%clustering_weight_glob(2))

      case("inversion.clustering.nClusters")
        read(10, *) ipar%nclusters
        call print_arg(myrank, parname, ipar%nclusters)

      case("inversion.clustering.mixtureFile")
        call read_filename(10, ipar%mixture_file)
        call print_arg(myrank, parname, ipar%mixture_file)

      case("inversion.clustering.cellWeightsFile")
        call read_filename(10, ipar%cell_weights_file)
        call print_arg(myrank, parname, ipar%cell_weights_file)

      case("inversion.clustering.optimizationType")
        read(10, *) ipar%clustering_opt_type
        call print_arg(myrank, parname, ipar%clustering_opt_type)

      case("inversion.clustering.constraintsType")
        read(10, *) ipar%clustering_constraints_type
        call print_arg(myrank, parname, ipar%clustering_constraints_type)

      ! OUTPUT parameters -------------------------------------------------------

      case("output.paraview.grav.modelLabel")
        read(10, '(a)') gpar%vtk_model_label
        call print_arg(myrank, parname, gpar%vtk_model_label)

      case("output.paraview.magn.modelLabel")
        read(10, '(a)') mpar%vtk_model_label
        call print_arg(myrank, parname, mpar%vtk_model_label)

      case default
        print *, "WARNING: Unknown parameter name! Name =", parname
        read(10, *, iostat=ios)

    end select
  enddo

  ! Sanity check.
  if (ipar%admm_type > 0 .and. ipar%admm_bound_type == 1) then
    if (.not. global_bounds_defined) then
      call exit_MPI("Global admm bounds are not defined! They must be defined in the Parfile.", myrank, 0)
    endif
  endif

  ! Sanity check.
  if (index(trim(adjustl(gpar%vtk_model_label)), " ") > 0 .or. &
      index(trim(adjustl(mpar%vtk_model_label)), " ") > 0) then
    call exit_MPI("The vtk model label cannot contain spaces!", myrank, 0)
  endif

  print *, "Finished reading the parameter file."

end subroutine read_parfile

end module init_parameters
