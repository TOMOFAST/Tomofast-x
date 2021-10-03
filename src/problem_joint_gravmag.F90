
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

!===============================================================================================
! A class for solving joint gravity & magnetism problem.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!===============================================================================================
module problem_joint_gravmag

  use global_typedefs
  use parameters_grav
  use parameters_mag
  use parameters_gravmag
  use parameters_inversion
  use inversion_arrays
  use forward_problem_gravmag
  use weights_gravmag
  use joint_inverse_problem
  use costs
  use data_gravmag
  use parallel_tools
  use filter
  use noise
  use compare_models
  use string, only: str

  implicit none

  private

  ! Flags for writing data for visualization.
  logical, parameter :: WRITE_SENSITIVITY     = .true.
  logical, parameter :: WRITE_DAMPING_WEIGHT  = .false.

  ! Unit number for cost file handle.
  integer, parameter :: FILE_COSTS = 1234567

  type, public :: t_problem_joint_gravmag
    private

    ! Local number of elements of one problem.
    integer :: nelements

    ! Model change (update) at inversion iteration.
    real(kind=CUSTOM_REAL), allocatable :: delta_model(:)

  contains
    private

    procedure, public, pass :: initialize => problem_joint_gravmag_initialize
    procedure, public, pass :: solve_problem_joint_gravmag

    procedure, private, nopass :: calculate_model_costs
    procedure, private, nopass :: write_sensitivity_matrix
    procedure, private, nopass :: read_model

  end type t_problem_joint_gravmag

contains

!================================================================================================
! Initialization.
!================================================================================================
subroutine problem_joint_gravmag_initialize(this, nelements, myrank)
  class(t_problem_joint_gravmag), intent(inout) :: this
  integer, intent(in) :: nelements, myrank

  integer :: ierr

  this%nelements = nelements

  ierr = 0

  if (.not. allocated(this%delta_model)) &
    allocate(this%delta_model(2 * this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in problem_joint_gravmag_initialize!", myrank, ierr)

end subroutine problem_joint_gravmag_initialize

!===================================================================================
! Solves gravity AND magnetism joint problem (forward + inversion).
!===================================================================================
subroutine solve_problem_joint_gravmag(this, gpar, mpar, ipar, myrank, nbproc)
  class(t_problem_joint_gravmag), intent(inout) :: this
  type(t_parameters_grav), intent(in) :: gpar
  type(t_parameters_mag), intent(in) :: mpar
  type(t_parameters_inversion), intent(in) :: ipar
  integer, intent(in) :: myrank, nbproc

  type(t_joint_inversion) :: joint_inversion
  ! Data arrays for (1) gravity and (2) magnetism inversions.
  type(t_inversion_arrays) :: iarr(2)
  type(t_data) :: data(2)
  type(t_weights) :: weights
  type(t_parallel_tools) :: pt
  type(t_compare_models) :: comp
  real(kind=CUSTOM_REAL) :: compres(2)
  real(kind=CUSTOM_REAL) :: cost_data(2)
  real(kind=CUSTOM_REAL) :: cost_model(2)
  integer :: it, i, m, number_prior_models, ierr
  character(len=256) :: path_output_parfile
  character(len=256) :: grav_prior_model_filename, mag_prior_model_filename

  logical :: SOLVE_PROBLEM(2)
  integer :: nnz(2)

  if (myrank == 0) print *, "Solving problem joint grav/mag."

  do i = 1, 2
    SOLVE_PROBLEM(i) = (ipar%problem_weight(i) /= 0.d0)
    if (myrank == 0) print *, "SOLVE_PROBLEM("//trim(str(i))//") = ", SOLVE_PROBLEM(i)
  enddo

  ! (I) MODEL ALLOCATION.  ---------------------------------------------------------------

  if (myrank == 0) print *, "(I) MODEL ALLOCATION."

  ! Initialize inversion arrays dimensions.
  if (SOLVE_PROBLEM(1)) call iarr(1)%initialize(ipar%nelements, ipar%ndata(1), ipar%nx, ipar%ny, ipar%nz)
  if (SOLVE_PROBLEM(2)) call iarr(2)%initialize(ipar%nelements, ipar%ndata(2), ipar%nx, ipar%ny, ipar%nz)

  ! Allocate memory for model (with grid) objects.
  if (SOLVE_PROBLEM(1)) call iarr(1)%init_model(myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call iarr(2)%init_model(myrank, nbproc)

  ! Reading the full grid and model.
  if (SOLVE_PROBLEM(1)) call iarr(1)%model%read(gpar%model_files(1), .true., myrank)
  if (SOLVE_PROBLEM(2)) call iarr(2)%model%read(mpar%model_files(1), .true., myrank)

#ifndef SUPPRESS_OUTPUT
  ! Write the model read to a file for Paraview visualization.
  if (SOLVE_PROBLEM(1)) call iarr(1)%model%write('grav_read_', .false., myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call iarr(2)%model%write('mag_read_', .false., myrank, nbproc)
#endif

  ! Distribute the model and grid among CPUs.
  if (SOLVE_PROBLEM(1)) call iarr(1)%model%distribute(myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call iarr(2)%model%distribute(myrank, nbproc)

  ! (I2) SETTING ADMM BOUNDS --------------------------------------------------------------

  if (ipar%admm_type > 0) then
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) then
        ! Reading min/max ADMM bounds from file.
        call iarr(i)%model%allocate_bound_arrays(ipar%nlithos, myrank)
        call iarr(i)%model%read_bound_constraints(ipar%bounds_ADMM_file(i), myrank, nbproc)
      endif
    enddo
  endif

  ! (II) DATA ALLOCATION. -----------------------------------------------------------------

  if (myrank == 0) print *, "(II) DATA ALLOCATION."

  ! Allocate memory for data objects.
  if (SOLVE_PROBLEM(1)) call data(1)%initialize(gpar%ndata, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%initialize(mpar%ndata, myrank)

  ! Reading the GRID ONLY for data point (needed to generate sensitivity matrix).
  if (SOLVE_PROBLEM(1)) call data(1)%read(gpar%data_grid_file, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%read(mpar%data_grid_file, myrank)
 
  ! (III) SENSITIVITY MATRIX ALLOCATION  ---------------------------------------------------

  if (myrank == 0) print *, "(III) SENSITIVITY MATRIX ALLOCATION."

  ! Memory allocation for auxiliarily inversion arrays.
  if (SOLVE_PROBLEM(1)) call iarr(1)%allocate_aux(myrank)
  if (SOLVE_PROBLEM(2)) call iarr(2)%allocate_aux(myrank)

  !-----------------------------------------------------------------------------------------
  ! Calculates weights.
  if (SOLVE_PROBLEM(1)) call weights%calculate(gpar, iarr(1), data(1), 1, myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call weights%calculate(mpar, iarr(2), data(2), 2, myrank, nbproc)

  ! Precondition the column weights.
  if (SOLVE_PROBLEM(1)) iarr(1)%column_weight = ipar%column_weight_multiplier(1) * iarr(1)%column_weight
  if (SOLVE_PROBLEM(2)) iarr(2)%column_weight = ipar%column_weight_multiplier(2) * iarr(2)%column_weight

  !-------------------------------------------------------------------------------------------------------
  ! Calculate nnz for the sensitivity kernel.
  nnz = 0
  if (ipar%compression_type > 0) then
  ! Wavelet compression. Calculate the number of nonzero elements in the compressed kernel (on every rank).
  ! Also calculate the sensitivity-based depth-weight.
    if (SOLVE_PROBLEM(1)) &
      call calculate_kernel_size_and_weight(gpar, iarr(1), data(1), nnz(1), ipar%column_weight_multiplier(1), myrank, nbproc)

    if (SOLVE_PROBLEM(2)) &
      call calculate_kernel_size_and_weight(mpar, iarr(2), data(2), nnz(2), ipar%column_weight_multiplier(2), myrank, nbproc)
  else
    nnz(1) = ipar%nelements * ipar%ndata(1)
    nnz(2) = ipar%nelements * ipar%ndata(2)
  endif

  ! Allocate the sensitivity kernel.
  if (SOLVE_PROBLEM(1)) call iarr(1)%allocate_sensit(.false., nnz(1), myrank)
  if (SOLVE_PROBLEM(2)) call iarr(2)%allocate_sensit(.false., nnz(2), myrank)

  !-----------------------------------------------------------------------------------------
  ! Solve forward problems for gravity and magnetism.
  if (SOLVE_PROBLEM(1)) call solve_forward_problem(gpar, iarr(1), data(1), myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call solve_forward_problem(mpar, iarr(2), data(2), myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
  ! Write data calculated from the model read.
  if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_read_', 2, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_read_', 2, myrank)
#endif

  ! Reading the data. Read here to allow the use of the above calculated data from the (original) model read.
  if (SOLVE_PROBLEM(1)) call data(1)%read(gpar%data_file, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%read(mpar%data_file, myrank)

#ifndef SUPPRESS_OUTPUT
  ! Write the observed (measured) data.
  if (SOLVE_PROBLEM(1)) call data(1)%write('grav_observed_', 1, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%write('mag_observed_', 1, myrank)
#endif

  !-----------------------------------------------------------------------------------------
  number_prior_models = gpar%number_prior_models
  path_output_parfile = path_output

  ! Initialize joint inversion object.
  ! Read the clustering parameters from file inside.
  call joint_inversion%initialize(ipar, nnz(1) + nnz(2), myrank)

  !******************************************************************************************
  ! Loop over different prior models.
  !******************************************************************************************
  do m = 1, number_prior_models

    ! Define prior model and output folder names.
    if (m > 1) then
      path_output = trim(path_output_parfile)//'_'//trim(str(m))//'/'
      grav_prior_model_filename = trim(gpar%model_files(2))//'_'//trim(str(m))
      mag_prior_model_filename = trim(mpar%model_files(2))//'_'//trim(str(m))
    else
      grav_prior_model_filename = gpar%model_files(2)
      mag_prior_model_filename = mpar%model_files(2)
    endif

    if (myrank == 0) then
      print *
      print *, '********************************************************************************'
      print *, 'Solve problem for prior model #', m, ', output folder = ', trim(path_output)
      print *, '********************************************************************************'
    endif

    if (m > 1) call joint_inversion%reset()
    
    ! SETTING PRIOR MODEL FOR INVERSION  -----------------------------------------------------
    if (SOLVE_PROBLEM(1)) &
      call read_model(iarr(1), gpar%prior_model_type, gpar%prior_model_val, grav_prior_model_filename, myrank, nbproc)
    if (SOLVE_PROBLEM(2)) &
      call read_model(iarr(2), mpar%prior_model_type, mpar%prior_model_val, mag_prior_model_filename, myrank, nbproc)

    ! Set the prior model.
    if (SOLVE_PROBLEM(1)) iarr(1)%model_prior = iarr(1)%model%val
    if (SOLVE_PROBLEM(2)) iarr(2)%model_prior = iarr(2)%model%val

#ifndef SUPPRESS_OUTPUT
    ! Write the prior model to a file for visualization.
    if (SOLVE_PROBLEM(1)) call iarr(1)%model%write('grav_prior_', .false., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call iarr(2)%model%write('mag_prior_', .false., myrank, nbproc)
#endif

    !-----------------------------------------------------------------------------------------
    ! Calculate data from the prior model.
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) call iarr(i)%model%calculate_data(ipar%ndata(i), iarr(i)%matrix_sensit, &
        iarr(i)%column_weight, data(i)%val_calc, ipar%compression_type, myrank, nbproc)
    enddo

#ifndef SUPPRESS_OUTPUT
    ! Write data calculated from the prior model.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_prior_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_prior_', 2, myrank)
#endif

    ! SETTING STARTING MODEL FOR INVERSION  -----------------------------------------------------
    if (SOLVE_PROBLEM(1)) call read_model(iarr(1), gpar%start_model_type, gpar%start_model_val, gpar%model_files(3), myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call read_model(iarr(2), mpar%start_model_type, mpar%start_model_val, mpar%model_files(3), myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
    ! Write the starting model to a file for visualization.
    if (SOLVE_PROBLEM(1)) call iarr(1)%model%write('grav_starting_', .true., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call iarr(2)%model%write('mag_starting_', .true., myrank, nbproc)
#endif

    !-----------------------------------------------------------------------------------------
    ! Calculate data from the starting model.
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) call iarr(i)%model%calculate_data(ipar%ndata(i), iarr(i)%matrix_sensit, &
        iarr(i)%column_weight, data(i)%val_calc, ipar%compression_type, myrank, nbproc)
    enddo

#ifndef SUPPRESS_OUTPUT
    ! Write data calculated from the starting model.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_starting_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_starting_', 2, myrank)
#endif

    !-----------------------------------------------------------------------------------------
    ! Calculate initial cost (misfit).
    if (SOLVE_PROBLEM(1)) call calculate_cost(ipar%ndata(1), data(1)%val_meas, data(1)%val_calc, cost_data(1), myrank)
    if (SOLVE_PROBLEM(2)) call calculate_cost(ipar%ndata(2), data(2)%val_meas, data(2)%val_calc, cost_data(2), myrank)

    ! Calculate costs for the models (damping term in the cost function).
    call calculate_model_costs(ipar, iarr, cost_model, SOLVE_PROBLEM, myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
    ! Stores costs.
    if (myrank == 0) &
      open(FILE_COSTS, file=trim(path_output)//'costs', access='stream', form='formatted', status='unknown', action='write')
#endif

    ! Non-linear inversion loop.
    do it = 1, ipar%ninversions

      if (myrank == 0) then
        print *, '======================================================='
        print *, 'Iteration, prior model =', it, ' ', m
        print *, '======================================================='
      endif

      ! Calculate data residuals.
      if (SOLVE_PROBLEM(1)) iarr(1)%residuals = data(1)%val_meas - data(1)%val_calc
      if (SOLVE_PROBLEM(2)) iarr(2)%residuals = data(2)%val_meas - data(2)%val_calc

      ! Resets the joint inversion.
      if (it > 1) call joint_inversion%reset()

      ! Solve joint inverse problem.
      call joint_inversion%solve(ipar, iarr, this%delta_model, myrank, nbproc)

      ! Update the local models.
      if (SOLVE_PROBLEM(1)) call iarr(1)%model%update(this%delta_model(1:this%nelements))
      if (SOLVE_PROBLEM(2)) call iarr(2)%model%update(this%delta_model(this%nelements + 1:))

      ! Update the full models (needed e.g. for cross-gradient right-hand-side).
      do i = 1, 2
        if (SOLVE_PROBLEM(i)) &
          call pt%get_full_array(iarr(i)%model%val, ipar%nelements, iarr(i)%model%val_full, .true., myrank, nbproc)
      enddo

      ! Calculate data based on the new model from inversion.
      do i = 1, 2
        if (SOLVE_PROBLEM(i)) call iarr(i)%model%calculate_data(ipar%ndata(i), iarr(i)%matrix_sensit, &
          iarr(i)%column_weight, data(i)%val_calc, ipar%compression_type, myrank, nbproc)
      enddo

#ifndef SUPPRESS_OUTPUT
      ! Write costs (for the previous iteration).
      if (myrank == 0) write(FILE_COSTS, *) it - 1, cost_data(1), cost_data(2), cost_model(1), cost_model(2), &
                                            joint_inversion%get_cross_grad_cost(), &
                                            joint_inversion%get_clustering_cost(1), joint_inversion%get_clustering_cost(2)
#endif

      ! Calculate new costs for data misfits.
      if (SOLVE_PROBLEM(1)) call calculate_cost(ipar%ndata(1), data(1)%val_meas, data(1)%val_calc, cost_data(1), myrank)
      if (SOLVE_PROBLEM(2)) call calculate_cost(ipar%ndata(2), data(2)%val_meas, data(2)%val_calc, cost_data(2), myrank)

      ! Calculate new costs for the models (damping term in the cost function).
      call calculate_model_costs(ipar, iarr, cost_model, SOLVE_PROBLEM, myrank, nbproc)

      ! Store final models from the single (non-joint) inversions.
      do i = 1, 2
        if (SOLVE_PROBLEM(i)) then
          if (ipar%single_problem_complete(i, it) .and. .not. ipar%single_problem_complete(i, it - 1)) then
            iarr(i)%model%val_final0 = iarr(i)%model%val
          endif
        endif
      enddo

    enddo

#ifndef SUPPRESS_OUTPUT
    ! Write final costs (excluding cross-gradient cost, as it is being calculated only during solution).
    if (myrank == 0) write(FILE_COSTS, *) ipar%ninversions, cost_data(1), cost_data(2), cost_model(1), cost_model(2)
    if (myrank == 0) close(FILE_COSTS)
#endif

#ifndef SUPPRESS_OUTPUT
    if (myrank == 0) then
      ! Print model value bounds.
      do i = 1, 2
        if (SOLVE_PROBLEM(i)) &
          print *, 'Model', i , 'min/max values =', minval(iarr(i)%model%val_full), maxval(iarr(i)%model%val_full)
      enddo
    endif
#endif

#ifndef SUPPRESS_OUTPUT
    ! Compare final models of single and joint inversions.
    compres = 0.d0
    if (SOLVE_PROBLEM(1)) call comp%compare(iarr(1)%model, ipar%derivative_type, compres(1), myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call comp%compare(iarr(2)%model, ipar%derivative_type, compres(2), myrank, nbproc)

    if (myrank == 0) print *, 'Model comparison:', ipar%column_weight_multiplier(1), ipar%column_weight_multiplier(2), &
                              compres(1), compres(2), ipar%cross_grad_weight, joint_inversion%get_cross_grad_cost()
#endif

#ifndef SUPPRESS_OUTPUT
    ! Write the final model to a file.
    if (SOLVE_PROBLEM(1)) call iarr(1)%model%write('grav_final_', .false., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call iarr(2)%model%write('mag_final_', .false., myrank, nbproc)
#endif

#ifndef SUPPRESS_OUTPUT
    ! Write data calculated from final model.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_final_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_final_', 2, myrank)
#endif

#ifndef SUPPRESS_OUTPUT
    if (joint_inversion%add_cross_grad) then
      ! Write final cross-gradient vector magnitude to a file.
      iarr(1)%model%val_full = joint_inversion%get_cross_grad()

      call iarr(1)%model%write('cross_grad_final_', .false., myrank, nbproc)
    endif
#endif

#ifndef SUPPRESS_OUTPUT
    if (joint_inversion%add_clustering) then
      ! Write final clustering probabilities, i.e., P(m) per cell.
      iarr(1)%model%val_full = joint_inversion%get_clustering()

      call iarr(1)%model%write('clustering_final_', .false., myrank, nbproc)

      call joint_inversion%clustering%write_data('clustering_data.txt', iarr(1)%model%grid_full, myrank)
    endif
#endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  enddo ! loop over prior models
  !******************************

#ifndef SUPPRESS_OUTPUT
  if (WRITE_DAMPING_WEIGHT) then
    ! Write the damping weight.
    if (SOLVE_PROBLEM(1)) iarr(1)%model%val = iarr(1)%damping_weight
    if (SOLVE_PROBLEM(2)) iarr(2)%model%val = iarr(2)%damping_weight

    if (SOLVE_PROBLEM(1)) call iarr(1)%model%write('damping_weight_grav_', .true., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call iarr(2)%model%write('damping_weight_mag_', .true., myrank, nbproc)
  endif
#endif

#ifndef SUPPRESS_OUTPUT
  if (WRITE_SENSITIVITY) then
    ! Write a root mean square sensitivity (integrated sensitivity).
    call write_sensitivity_matrix(iarr, SOLVE_PROBLEM, myrank, nbproc)
  endif
#endif

end subroutine solve_problem_joint_gravmag

!========================================================================================
! Write a root mean square sensitivity (integrated sensitivity) to a file.
!========================================================================================
subroutine write_sensitivity_matrix(iarr, solve_problem, myrank, nbproc)
  type(t_inversion_arrays), intent(inout) :: iarr(2)
  logical, intent(in) :: solve_problem(2)
  integer, intent(in) :: myrank, nbproc

  integer :: i

  ! Loop over problems.
  do i = 1, 2
    if (solve_problem(i)) call iarr(i)%matrix_sensit%get_integrated_sensit(iarr(i)%model%val)
  enddo

   ! Write sensitivity to files.
  if (solve_problem(1)) call iarr(1)%model%write('sensit_grav_', .true., myrank, nbproc)
  if (solve_problem(2)) call iarr(2)%model%write('sensit_mag_', .true., myrank, nbproc)

end subroutine write_sensitivity_matrix

!========================================================================================
! Computes and prints norm Lp of the difference between inverted and prior models.
!========================================================================================
subroutine calculate_model_costs(ipar, iarr, cost_model, solve_problem, myrank, nbproc)
  type(t_parameters_inversion), intent(in) :: ipar
  type(t_inversion_arrays), intent(in) :: iarr(2)
  logical, intent(in) :: solve_problem(2)
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(out) :: cost_model(2)

  integer :: i

  cost_model = 0.d0

  do i = 1, 2
    if (solve_problem(i)) then
      call calculate_cost_model(ipar%nelements, ipar%norm_power, iarr(i)%model%val, iarr(i)%model_prior, &
                                iarr(i)%damping_weight, cost_model(i), nbproc)

      if (myrank == 0) print *, 'model cost =', cost_model(i)
    endif
  enddo
end subroutine calculate_model_costs

!========================================================================================
! Computes and prints norm Lp of the difference between inverted and prior models.
!========================================================================================
subroutine read_model(iarr, model_type, model_val, model_file, myrank, nbproc)
  integer, intent(in) :: model_type, myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: model_val
  character(len=256), intent(in) :: model_file
  type(t_inversion_arrays), intent(inout) :: iarr

  if (model_type == 1) then
    ! Setting homogeneous starting value.
    iarr%model%val_full = model_val

  else if (model_type == 2) then
    ! Reading from file.
    call iarr%model%read(model_file, .false., myrank)

  else
    print *, "Unknown model type!"
    stop
  endif

  ! Distribute the model and grid among CPUs.
  call iarr%model%distribute(myrank, nbproc)
end subroutine read_model

end module problem_joint_gravmag
