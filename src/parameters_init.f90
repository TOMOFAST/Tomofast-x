
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
    call read_parfile(epar, gpar, mpar, ipar, myrank)

    ! Print out if we do this in double or single precision.
    if (myrank == 0) then
      if (CUSTOM_REAL == SIZE_DOUBLE) then
        print *, "precision = DOUBLE"
      else
        print *, "precision = SINGLE"
      endif
    endif

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

    ipar%compression_type = gpar%compression_type
    ipar%wavelet_threshold = gpar%wavelet_threshold

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

  ! GLOBAL parameters.
  path_output = "output/test/"

  ! MODEL GRID parameters.
  gpar%nx = 0
  gpar%ny = 0
  gpar%nz = 0
  mpar%nx = 0
  mpar%ny = 0
  mpar%nz = 0
  gpar%model_files(1) = "NILL"
  mpar%model_files(1) = "NILL"

  ! DATA parameters.
  gpar%ndata = 0
  mpar%ndata = 0
  gpar%data_grid_file = "NILL"
  mpar%data_grid_file = "NILL"
  gpar%data_file = "NILL"
  mpar%data_file = "NILL"

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
  gpar%compression_type = 0
  gpar%wavelet_threshold = 1.d-7
  gpar%distance_threshold = 1.d+10 ! Source to the cell distance (m).
  gpar%compression_rate = 1.d0     ! 1.0 = full matrix

  mpar%compression_type = gpar%compression_type
  mpar%wavelet_threshold = gpar%wavelet_threshold
  mpar%distance_threshold = gpar%distance_threshold
  mpar%compression_rate = gpar%compression_rate

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

  ! INVERSION parameters.
  ipar%ninversions = 10
  ipar%niter = 100
  ipar%rmin = 1.d-13
  ipar%method = 1 ! LSQR = 1
  ipar%gamma = 0. ! soft threshold ("L1-norm", no=0.)

  ! MODEL DAMPING (m - m_prior).
  ipar%alpha(1) = 1.d-11
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

  !***********************************************
  ! ECT parameters:
  !***********************************************
  ! ECT GRID parameters.
  epar%dims%nr = 36
  epar%dims%ntheta = 36
  epar%dims%nz = 36

  ! ECT GEOMETRY parameters.
  epar%nel = 36
  epar%nrings = 3
  epar%dims%kguards = 0
  epar%ifixed_elecgeo = 0 ! NO=0, YES=1
  epar%irefine = 0        ! NO=0, YES=1
  epar%sens%radiusin = 0.045d0
  epar%sens%radiusout = 0.06d0
  epar%sens%radiusoutout = 0.07d0
  epar%sens%heicyl = 0.2d0
  epar%sens%space_elec_guards = 0.d0
  epar%sens%space_electrodes = 0.d0

  ! ECT MODEL parameters.
  epar%num_bubbles = 4
  epar%filename_bubbles = "data/ECT/bubble_4vert.dat"
  epar%permit0 = 1.d0
  epar%permit_air = 1.d0
  epar%permit_isolated_tube = 3.5d0
  epar%permit_oil = 2.d0

  ! ECT SOLVER parameters.
  epar%linear_solver = LINSOLV_PCG ! Not exposed to Parfile.
  epar%iprecond = 1 ! 0=NO, YES>0
  epar%omega1 = 0.8d0
  epar%itypenorm = 1 ! 1=L2, 2=max
  epar%itmax = 1000
  epar%output_frequency = 20
  epar%tol = 1.d-12

  ! MULTIGRID parameters.
  ! Removed multigrid, keep this not to change a lot of code.
  epar%ilevel_coarse = 1
  epar%coarse_solver = LINSOLV_PCG

end subroutine set_default_parameters

!===================================================================================
! Read input parameters from Parfile.
!===================================================================================
subroutine read_parfile(epar, gpar, mpar, ipar, myrank)
  integer, intent(in) :: myrank

  type(t_parameters_ect), intent(inout) :: epar
  type(t_parameters_grav), intent(inout) :: gpar
  type(t_parameters_mag), intent(inout) :: mpar
  type(t_parameters_inversion), intent(inout) :: ipar

  character(len=1) :: ch
  character(len=256) :: parfile_name
  character(len=128) :: parname
  character(len=256) :: line
  character(len=256) :: parfile_description
  integer :: itmp
  integer :: symbol_index, i
  integer :: ios

  ! The name of the Parfile should be passed in via the command line.
  call get_command_argument(2, parfile_name)
  if (len_trim(parfile_name) == 0) then
    call exit_MPI("No Parfile supplied!", myrank, 0)
    stop
  endif

  open(unit=10, file=parfile_name, status='old', iostat=itmp, action='read')
  if (itmp /= 0) call exit_MPI("Parfile """ // trim(parfile_name) // """ cannot be opened!", myrank, 0)

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

      ! GLOBAL parameters -----------------------------------
      case("global.outputFolderPath")
        call read_filename(10, path_output)
        call print_arg(myrank, parname, path_output)

      case("global.description")
        read(10, 3) parfile_description
        call print_arg(myrank, parname, parfile_description)

      ! ECT GRID parameters ---------------------------------

      case("forward.ect.grid.nr")
        read(10, 2) epar%dims%nr
        call print_arg(myrank, parname, epar%dims%nr)

      case("forward.ect.grid.ntheta")
        read(10, 2) epar%dims%ntheta
        call print_arg(myrank, parname, epar%dims%ntheta)

      case("forward.ect.grid.nz")
        read(10, 2) epar%dims%nz
        call print_arg(myrank, parname, epar%dims%nz)

      ! ECT GEOMETRY parameters -----------------------------

      case("forward.ect.geometry.nElectrodes")
        read(10, 2) epar%nel
        call print_arg(myrank, parname, epar%nel)

      case("forward.ect.geometry.nRings")
        read(10, 2) epar%nrings
        call print_arg(myrank, parname, epar%nrings)

      case("forward.ect.geometry.kguards")
        read(10, 2) epar%dims%kguards
        if (epar%dims%kguards == 0) epar%dims%kguards = epar%dims%nz / 4
        call print_arg(myrank, parname, epar%dims%kguards)

      case("forward.ect.geometry.fixedElectrodes")
        read(10, 2) epar%ifixed_elecgeo
        call print_arg(myrank, parname, epar%ifixed_elecgeo)

      case("forward.ect.geometry.meshRefinement")
        read(10, 2) epar%irefine
        call print_arg(myrank, parname, epar%irefine)

      case("forward.ect.geometry.locationR1")
        read(10, 1) epar%sens%radiusin
        call print_arg(myrank, parname, epar%sens%radiusin)

      case("forward.ect.geometry.locationR2")
        read(10, 1) epar%sens%radiusout
        call print_arg(myrank, parname, epar%sens%radiusout)

      case("forward.ect.geometry.locationR3")
        read(10, 1) epar%sens%radiusoutout
        call print_arg(myrank, parname, epar%sens%radiusoutout)

      case("forward.ect.geometry.sensorHeight")
        read(10, 1) epar%sens%heicyl
        call print_arg(myrank, parname, epar%sens%heicyl)

      case("forward.ect.geometry.spaceBetweenGuards")
        read(10, 1) epar%sens%space_elec_guards
        call print_arg(myrank, parname, epar%sens%space_elec_guards)

      case("forward.ect.geometry.spaceBetweenElectrodes")
        read(10, 1) epar%sens%space_electrodes
        call print_arg(myrank, parname, epar%sens%space_electrodes)

      ! ECT MODEL parameters --------------------------------

      case("forward.ect.model.nBubbles")
        read(10, 2) epar%num_bubbles
        call print_arg(myrank, parname, epar%num_bubbles)

      case("forward.ect.model.bubblesLocationFile")
        call read_filename(10, epar%filename_bubbles)
        call print_arg(myrank, parname, epar%filename_bubbles)

      case("forward.ect.model.absolutePermittivity")
        read(10, 1) epar%permit0
        call print_arg(myrank, parname, epar%permit0)

      case("forward.ect.model.permittivityAir")
        read(10, 1) epar%permit_air
        call print_arg(myrank, parname, epar%permit_air)

      case("forward.ect.model.permittivityIsolatedTube")
        read(10, 1) epar%permit_isolated_tube
        call print_arg(myrank, parname, epar%permit_isolated_tube)

      case("forward.ect.model.permittivityOil")
        read(10, 1) epar%permit_oil
        call print_arg(myrank, parname, epar%permit_oil)

      ! ECT SOLVER parameters -------------------------------

      case("forward.ect.solver.precond")
        read(10, 2) epar%iprecond
        call print_arg(myrank, parname, epar%iprecond)

      case("forward.ect.solver.precond.relaxOmega")
        read(10, 1) epar%omega1
        call print_arg(myrank, parname, epar%omega1)

      case("forward.ect.solver.normType")
        read(10, 2) epar%itypenorm
        call print_arg(myrank, parname, epar%itypenorm)

      case("forward.ect.solver.nMaxIterations")
        read(10, 2) epar%itmax
        call print_arg(myrank, parname, epar%itmax)

      case("forward.ect.solver.outputFrequencyIter")
        read(10, 2) epar%output_frequency
        call print_arg(myrank, parname, epar%output_frequency)

      case("forward.ect.solver.tolerance")
        read(10, 1) epar%tol
        call print_arg(myrank, parname, epar%tol)

      ! MODEL GRID parameters -------------------------------

      case("modelGrid.size")
        read(10, *) gpar%nx, gpar%ny, gpar%nz
        if (myrank == 0) print *, trim(parname)//" =", gpar%nx, gpar%ny, gpar%nz
        mpar%nx = gpar%nx
        mpar%ny = gpar%ny
        mpar%nz = gpar%nz

      case("modelGrid.grav.file")
        call read_filename(10, gpar%model_files(1))
        call print_arg(myrank, parname, gpar%model_files(1))

      case("modelGrid.magn.file")
        call read_filename(10, mpar%model_files(1))
        call print_arg(myrank, parname, mpar%model_files(1))

      ! DATA parameters -------------------------------------

      case("forward.data.grav.nData")
        read(10, 2) gpar%ndata
        call print_arg(myrank, parname, gpar%ndata)

      case("forward.data.magn.nData")
        read(10, 2) mpar%ndata
        call print_arg(myrank, parname, mpar%ndata)

      case("forward.data.grav.dataGridFile")
        call read_filename(10, gpar%data_grid_file)
        call print_arg(myrank, parname, gpar%data_grid_file)

      case("forward.data.magn.dataGridFile")
        call read_filename(10, mpar%data_grid_file)
        call print_arg(myrank, parname, mpar%data_grid_file)

      case("forward.data.grav.dataValuesFile")
        call read_filename(10, gpar%data_file)
        call print_arg(myrank, parname, gpar%data_file)

      case("forward.data.magn.dataValuesFile")
        call read_filename(10, mpar%data_file)
        call print_arg(myrank, parname, mpar%data_file)

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

      case("forward.depthWeighting.type")
        read(10, 2) gpar%depth_weighting_type
        call print_arg(myrank, parname, gpar%depth_weighting_type)
        mpar%depth_weighting_type = gpar%depth_weighting_type

      case("forward.depthWeighting.powerWeight.grav.beta")
        read(10, 1) gpar%beta
        call print_arg(myrank, parname, gpar%beta)

      case("forward.depthWeighting.powerWeight.grav.Z0")
        read(10, 1) gpar%Z0
        call print_arg(myrank, parname, gpar%Z0)

      case("forward.depthWeighting.powerWeight.magn.beta")
        read(10, 1) mpar%beta
        call print_arg(myrank, parname, mpar%beta)

      case("forward.depthWeighting.powerWeight.magn.Z0")
        read(10, 1) mpar%Z0
        call print_arg(myrank, parname, mpar%Z0)

      ! MATRIX COMPRESSION parameters ----------------------

      case("forward.matrixCompression.type")
        read(10, 2) gpar%compression_type
        call print_arg(myrank, parname, gpar%compression_type)
        mpar%compression_type = gpar%compression_type

      case("forward.matrixCompression.waveletThreshold")
        read(10, 1) gpar%wavelet_threshold
        call print_arg(myrank, parname, gpar%wavelet_threshold)
        mpar%wavelet_threshold = gpar%wavelet_threshold

      case("forward.matrixCompression.distanceThreshold")
        read(10, 1) gpar%distance_threshold
        call print_arg(myrank, parname, gpar%distance_threshold)
        mpar%distance_threshold = gpar%distance_threshold

      case("forward.matrixCompression.compressionRate")
        read(10, 1) gpar%compression_rate
        call print_arg(myrank, parname, gpar%compression_rate)
        mpar%compression_rate = gpar%compression_rate

      ! PRIOR MODEL -----------------------------------------

      case("inversion.priorModel.type")
        read(10, 2) gpar%prior_model_type
        call print_arg(myrank, parname, gpar%prior_model_type)
        mpar%prior_model_type = gpar%prior_model_type

      case("inversion.priorModel.nModels")
        read(10, 2) gpar%number_prior_models
        call print_arg(myrank, parname, gpar%number_prior_models)
        mpar%number_prior_models = gpar%number_prior_models

      case("inversion.priorModel.grav.value")
        read(10, 1) gpar%prior_model_val
        call print_arg(myrank, parname, gpar%prior_model_val)

      case("inversion.priorModel.magn.value")
        read(10, 1) mpar%prior_model_val
        call print_arg(myrank, parname, mpar%prior_model_val)

      case("inversion.priorModel.grav.file")
        call read_filename(10, gpar%model_files(2))
        call print_arg(myrank, parname, gpar%model_files(2))

      case("inversion.priorModel.magn.file")
        call read_filename(10, mpar%model_files(2))
        call print_arg(myrank, parname, mpar%model_files(2))

      ! STARTING MODEL -------------------------------------

      case("inversion.startingModel.type")
        read(10, 2) gpar%start_model_type
        call print_arg(myrank, parname, gpar%start_model_type)
        mpar%start_model_type = gpar%start_model_type

      case("inversion.startingModel.grav.value")
        read(10, 1) gpar%start_model_val
        call print_arg(myrank, parname, gpar%start_model_val)

      case("inversion.startingModel.magn.value")
        read(10, 1) mpar%start_model_val
        call print_arg(myrank, parname, mpar%start_model_val)

      case("inversion.startingModel.grav.file")
        call read_filename(10, gpar%model_files(3))
        call print_arg(myrank, parname, gpar%model_files(3))

      case("inversion.startingModel.magn.file")
        call read_filename(10, mpar%model_files(3))
        call print_arg(myrank, parname, mpar%model_files(3))

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

      case("inversion.modelDamping.grav.weight")
        read(10, 1) ipar%alpha(1)
        call print_arg(myrank, parname, ipar%alpha(1))

      case("inversion.modelDamping.magn.weight")
        read(10, 1) ipar%alpha(2)
        call print_arg(myrank, parname, ipar%alpha(2))

      case("inversion.modelDamping.ect.weight")
        read(10, 1) ipar%alpha(1)
        call print_arg(myrank, parname, ipar%alpha(1))

      case("inversion.modelDamping.normPower")
        read(10, 1) ipar%norm_power
        call print_arg(myrank, parname, ipar%norm_power)

      ! JOINT INVERSION parameters -------------------------------

      case("inversion.joint.grav.problemWeight")
        read(10, 1) ipar%problem_weight(1)
        call print_arg(myrank, parname, ipar%problem_weight(1))

      case("inversion.joint.magn.problemWeight")
        read(10, 1) ipar%problem_weight(2)
        call print_arg(myrank, parname, ipar%problem_weight(2))

      case("inversion.joint.grav.columnWeightMultiplier")
        read(10, 1) ipar%column_weight_multiplier(1)
        call print_arg(myrank, parname, ipar%column_weight_multiplier(1))

      case("inversion.joint.magn.columnWeightMultiplier")
        read(10, 1) ipar%column_weight_multiplier(2)
        call print_arg(myrank, parname, ipar%column_weight_multiplier(2))

      case("inversion.joint.grav.nIterSingle")
        read(10, 2) ipar%niter_single(1)
        call print_arg(myrank, parname, ipar%niter_single(1))

      case("inversion.joint.magn.nIterSingle")
        read(10, 2) ipar%niter_single(2)
        call print_arg(myrank, parname, ipar%niter_single(2))

      ! DAMPING-GRADIENT constraints -------------------------------

      case("inversion.dampingGradient.weightType")
        read(10, 2) ipar%damp_grad_weight_type
        call print_arg(myrank, parname, ipar%damp_grad_weight_type)

      case("inversion.dampingGradient.grav.weight")
        read(10, 1) ipar%beta(1)
        call print_arg(myrank, parname, ipar%beta(1))

      case("inversion.dampingGradient.magn.weight")
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

      case("inversion.clustering.grav.weight")
        read(10, 1) ipar%clustering_weight_glob(1)
        call print_arg(myrank, parname, ipar%clustering_weight_glob(1))

      case("inversion.clustering.magn.weight")
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

      case("inversion.admm.grav.boundsFile")
        call read_filename(10, ipar%bounds_ADMM_file(1))
        call print_arg(myrank, parname, ipar%bounds_ADMM_file(1))

      case("inversion.admm.magn.boundsFile")
        call read_filename(10, ipar%bounds_ADMM_file(2))
        call print_arg(myrank, parname, ipar%bounds_ADMM_file(2))

      case("inversion.admm.grav.weight")
        read(10, 1) ipar%rho_ADMM(1)
        call print_arg(myrank, parname, ipar%rho_ADMM(1))

      case("inversion.admm.magn.weight")
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

end subroutine read_parfile

end module init_parameters
