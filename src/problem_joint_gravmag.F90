
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
  use model
  use sensitivity_gravmag
  use weights_gravmag
  use joint_inverse_problem
  use costs
  use data_gravmag
  use string, only: str
  use model_IO
  use memory_tools

  implicit none

  private

  public :: solve_problem_joint_gravmag

  private :: calculate_model_costs
  private :: set_model
  private :: set_model_bounds
  private :: adjust_admm_weight
  private :: calculate_residual
  private :: exit_loop

contains

!===================================================================================
! Solves gravity AND magnetism joint problem (forward + inversion).
!===================================================================================
subroutine solve_problem_joint_gravmag(gpar, mpar, ipar, myrank, nbproc)
  type(t_parameters_grav), intent(inout) :: gpar
  type(t_parameters_mag), intent(inout) :: mpar
  type(t_parameters_inversion), intent(inout) :: ipar
  integer, intent(in) :: myrank, nbproc

  type(t_joint_inversion) :: jinv
  ! Data arrays for (1) gravity and (2) magnetism inversions.
  type(t_inversion_arrays) :: iarr(2)
  ! Data with its grid.
  type(t_data) :: data(2)
  ! Solution model (with grid).
  type(t_model) :: model(2)

  real(kind=CUSTOM_REAL) :: cost_data(2)
  real(kind=CUSTOM_REAL) :: cost_model(2)
  real(kind=CUSTOM_REAL) :: damping_gradient_cost(6)

  integer :: it, i, m, number_prior_models, ierr
  integer :: line_start(2), line_end(2), param_shift(2)

  character(len=256) :: path_output_parfile
  character(len=256) :: grav_prior_model_filename
  character(len=256) :: mag_prior_model_filename

  logical :: SOLVE_PROBLEM(2)
  logical :: allocate_full_model_on_all_cpus(2)
  integer(kind=8) :: nnz(2)
  integer :: nelements_new

  ! Unit number for cost file handle.
  integer, parameter :: FILE_COSTS = 1234567

  ! Model change (update) at inversion iteration.
  real(kind=CUSTOM_REAL), allocatable :: delta_model(:, :, :)
  ! Data change (update) at inversion iteration.
  type(t_real2d) :: delta_data(2)
  ! Memory usage.
  real(kind=CUSTOM_REAL) :: memory

  if (myrank == 0) print *, "Solving problem joint grav/mag."

  do i = 1, 2
    SOLVE_PROBLEM(i) = (ipar%problem_weight(i) /= 0.d0)
    if (myrank == 0) print *, "SOLVE_PROBLEM("//trim(str(i))//") = ", SOLVE_PROBLEM(i)
  enddo

  ! Define if we need to allocatye the full model array on all CPUs.
  ! Allocate it only when cross-gradient, damping gradient, or clustering constraints are used.
  allocate_full_model_on_all_cpus = .false.
  do i = 1, 2
    if (SOLVE_PROBLEM(i)) then
      if (ipar%cross_grad_weight /= 0.d0 .or. ipar%beta(i) /= 0.d0 .or. ipar%clustering_weight_glob(i) /= 0.d0) then
        allocate_full_model_on_all_cpus(i) = .true.
      endif
    endif
  enddo

  if (myrank == 0) print *, "allocate_full_model_on_all_cpus =", allocate_full_model_on_all_cpus

  nnz = 0
  cost_data = 0._CUSTOM_REAL
  cost_model = 0._CUSTOM_REAL

  ! (I) MODEL GRID ALLOCATION.  -----------------------------------------------------------

  if (myrank == 0) print *, "(I) MODEL GRID ALLOCATION."

  ! Allocate the model grid.
  if (SOLVE_PROBLEM(1)) call model(1)%grid_full%allocate(ipar%nx, ipar%ny, ipar%nz, myrank)
  if (SOLVE_PROBLEM(2)) call model(2)%grid_full%allocate(ipar%nx, ipar%ny, ipar%nz, myrank)

  ! Reading the full model grid.
  if (SOLVE_PROBLEM(1)) call read_model_grid(model(1)%grid_full, ipar%nmodel_components, gpar%model_files(1), myrank)
  if (SOLVE_PROBLEM(2)) call read_model_grid(model(2)%grid_full, ipar%nmodel_components, mpar%model_files(1), myrank)

  memory = get_max_mem_usage()
  if (myrank == 0) print *, "MEMORY USED (model grid) [GB] =", memory

  ! (II) DATA ALLOCATION. -----------------------------------------------------------------

  if (myrank == 0) print *, "(II) DATA ALLOCATION."

  ! Allocate memory for data objects.
  if (SOLVE_PROBLEM(1)) call data(1)%initialize(gpar%ndata, gpar%ndata_components, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%initialize(mpar%ndata, mpar%ndata_components, myrank)

  ! Reading the GRID ONLY for data points (needed to generate sensitivity matrix).
  if (SOLVE_PROBLEM(1)) call data(1)%read_grid(gpar%data_grid_file, myrank)
  if (SOLVE_PROBLEM(2)) call data(2)%read_grid(mpar%data_grid_file, myrank)

  memory = get_max_mem_usage()
  if (myrank == 0) print *, "MEMORY USED (data grid) [GB] =", memory

  ! (III) SENSITIVITY MATRIX CALCULATION  ---------------------------------------------------

  if (myrank == 0) print *, "(III) SENSITIVITY MATRIX CALCULATION."

  ! Memory allocation for auxiliarily inversion arrays.
  if (SOLVE_PROBLEM(1)) call iarr(1)%allocate_aux(ipar%nelements, ipar%ndata(1), ipar%ndata_components(1), myrank)
  if (SOLVE_PROBLEM(2)) call iarr(2)%allocate_aux(ipar%nelements, ipar%ndata(2), ipar%ndata_components(2), myrank)

  !-------------------------------------------------------------------------------------------------------
  if (gpar%sensit_read == 0) then
    ! Calculates the depth weights.
    if (SOLVE_PROBLEM(1)) call calculate_depth_weight(gpar, iarr(1), model(1)%grid_full, data(1), myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call calculate_depth_weight(mpar, iarr(2), model(2)%grid_full, data(2), myrank, nbproc)

    ! Precondition the column weights (to balance the columns in joint inversion).
    if (SOLVE_PROBLEM(1)) iarr(1)%column_weight = ipar%column_weight_multiplier(1) * iarr(1)%column_weight
    if (SOLVE_PROBLEM(2)) iarr(2)%column_weight = ipar%column_weight_multiplier(2) * iarr(2)%column_weight

    ! Calculate and write the sensitivity kernel to files.
    if (SOLVE_PROBLEM(1)) call calculate_and_write_sensit(gpar, model(1)%grid_full, data(1), iarr(1)%column_weight, &
                                                          myrank, nbproc)

    if (SOLVE_PROBLEM(2)) call calculate_and_write_sensit(mpar, model(2)%grid_full, data(2), iarr(2)%column_weight, &
                                                          myrank, nbproc)
  endif

  ! Calculate new partitioning for the load balancing.
  if (SOLVE_PROBLEM(1)) call calculate_new_partitioning(gpar, nnz(1), nelements_new, 1, myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call calculate_new_partitioning(mpar, nnz(2), nelements_new, 2, myrank, nbproc)

  ! Update the nelements for the nnz load balancing.
  if (SOLVE_PROBLEM(1)) gpar%nelements = nelements_new
  if (SOLVE_PROBLEM(2)) mpar%nelements = nelements_new
  ipar%nelements = nelements_new

  ! Reallocate the inversion arrays using the updated nelements value (for the nnz load balancing).
  if (SOLVE_PROBLEM(1)) call iarr(1)%reallocate_aux(ipar%nelements, ipar%ndata(1), ipar%ndata_components(1), myrank)
  if (SOLVE_PROBLEM(2)) call iarr(2)%reallocate_aux(ipar%nelements, ipar%ndata(2), ipar%ndata_components(2), myrank)

  !-------------------------------------------------------------------------------------------------------
  ! Deallocate the model grid.
  ! Keep the grid only on rank 0 for writing the models.
  ! Keep the grid on all ranks if we use gradient-based constraints (cross-gradient or damping gradient).
  !-------------------------------------------------------------------------------------------------------
  if (myrank /= 0) then
    if (ipar%cross_grad_weight == 0.d0) then
      do i = 1, 2
        if (SOLVE_PROBLEM(i) .and. ipar%beta(i) == 0.d0) then ! beta is the damping gradient weight.
          call model(i)%grid_full%deallocate()
        endif
      enddo
    endif
  endif

  memory = get_max_mem_usage()
  if (myrank == 0) print *, "MEMORY USED (sensit calc) [GB] =", memory

  ! (IV) MATRIX ALLOCATION  ------------------------------------------------------------------------------

  if (myrank == 0) print *, "(IV) MATRIX ALLOCATION."

  ! Allocate the sensitivity matrix.
  call jinv%initialize(ipar, nnz(1) + nnz(2), myrank)

  ! READING THE SENSITIVITY KERNEL ----------------------------------------------------------------------

  ! Reading the sensitivity kernel and depth weight from files.
  if (SOLVE_PROBLEM(1)) &
    call read_sensitivity_kernel(gpar, jinv%matrix, iarr(1)%column_weight, ipar%problem_weight(1), 1, myrank, nbproc)
  if (SOLVE_PROBLEM(2)) &
    call read_sensitivity_kernel(mpar, jinv%matrix, iarr(2)%column_weight, ipar%problem_weight(2), 2, myrank, nbproc)

  memory = get_max_mem_usage()
  if (myrank == 0) print *, "MEMORY USED (sensit read) [GB] =", memory

  ! RHS ALLOCATION -----------------------------------------------------------------------------------

  ! Allocate the right-hand side array.
  call jinv%initialize2(myrank)

  ! MODEL ALLOCATION -----------------------------------------------------------------------------------

  ! Allocate memory for the model.
  if (SOLVE_PROBLEM(1)) &
    call model(1)%initialize(ipar%nelements, ipar%nmodel_components, allocate_full_model_on_all_cpus(1), myrank, nbproc)
  if (SOLVE_PROBLEM(2)) &
    call model(2)%initialize(ipar%nelements, ipar%nmodel_components, allocate_full_model_on_all_cpus(2), myrank, nbproc)

  !-----------------------------------------------------------------------------------------------------
  ! Writing the column weight for Paraview visualisation.
  do i = 1, 2
    if (SOLVE_PROBLEM(i)) then
      model(i)%val(:, 1) = iarr(i)%column_weight
      call model_write(model(i), merge('grav_weight_', 'magn_weight_', i == 1), .true., .false., myrank, nbproc)
    endif
  enddo

  ! READING THE READ MODEL (SYNTHETIC) -----------------------------------------------------------------

  ! Reading the read model - that is stored in the model grid file.
  if (SOLVE_PROBLEM(1)) call model_read(model(1), gpar%model_files(1), myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call model_read(model(2), mpar%model_files(1), myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
  ! Write the model read to a file for Paraview visualization.
  if (SOLVE_PROBLEM(1)) call model_write(model(1), 'grav_read_', .false., .false., myrank, nbproc)
  if (SOLVE_PROBLEM(2)) call model_write(model(2), 'mag_read_', .false., .false., myrank, nbproc)
#endif

  ! SETTING THE ADMM BOUNDS -----------------------------------------------------------------------------

  if (ipar%admm_type > 0) then
    if (SOLVE_PROBLEM(1)) call set_model_bounds(ipar, model(1), 1, myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call set_model_bounds(ipar, model(2), 2, myrank, nbproc)
  endif

  ! SETTING damping gradient weights --------------------------------------------------------------------

  do i = 1, 2
    if (SOLVE_PROBLEM(i)) then
      if (ipar%beta(i) /= 0.) then
        call model(i)%allocate_damping_gradient_arrays(myrank)
        if (ipar%damp_grad_weight_type > 1) then
          call read_damping_gradient_weights(model(i), ipar%damping_gradient_file(i), myrank)
        endif
      endif
    endif
  enddo

  !-------------------------------------------------------------------------------------------------------
  ! Calculate parameters for calculating the data using the big (joint inversion) parallel sparse matrix.
  call jinv%calculate_matrix_partitioning(ipar, line_start, line_end, param_shift)

  !-------------------------------------------------------------------------------------------------------
  ! Calculate the data from the read model.
  do i = 1, 2
    if (SOLVE_PROBLEM(i)) call model(i)%calculate_data(ipar%ndata(i), ipar%ndata_components(i), jinv%matrix, &
      ipar%problem_weight(i), iarr(i)%column_weight, data(i)%val_calc, ipar%compression_type, &
      line_start(i), param_shift(i), myrank, nbproc)
  enddo

#ifndef SUPPRESS_OUTPUT
  ! Write data calculated from the read model.
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

  ! Allocate memory.
  allocate(delta_model(ipar%nelements, ipar%nmodel_components, 2), source=0._CUSTOM_REAL, stat=ierr)

  do i = 1, 2
    if (SOLVE_PROBLEM(i)) then
      allocate(delta_data(i)%val(ipar%ndata_components(i), ipar%ndata(i)), source=0._CUSTOM_REAL, stat=ierr)
    else
      ! Allocate 1 element not to have unallocated arrays.
      allocate(delta_data(i)%val(ipar%ndata_components(i), 1), source=0._CUSTOM_REAL, stat=ierr)
    endif
  enddo

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
      print *, '******************************************************************************'
      print *, 'Solve problem for prior model #', m, ', output folder = ', trim(path_output)
      print *, '******************************************************************************'
    endif

    if (m > 1) call jinv%reset(myrank)

    ! SETTING PRIOR MODEL FOR INVERSION -----------------------------------------------------
    if (SOLVE_PROBLEM(1)) &
      call set_model(model(1), gpar%prior_model_type, gpar%prior_model_val, grav_prior_model_filename, myrank, nbproc)
    if (SOLVE_PROBLEM(2)) &
      call set_model(model(2), mpar%prior_model_type, mpar%prior_model_val, mag_prior_model_filename, myrank, nbproc)

    ! TODO: Read values directly to val_prior in set_model(), by setting the flag for the model type.
    ! Set the prior model.
    if (SOLVE_PROBLEM(1)) model(1)%val_prior = model(1)%val
    if (SOLVE_PROBLEM(2)) model(2)%val_prior = model(2)%val

#ifndef SUPPRESS_OUTPUT
    ! Write the prior model to a file for visualization.
    if (SOLVE_PROBLEM(1)) call model_write(model(1), 'grav_prior_', .false., .false., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call model_write(model(2), 'mag_prior_', .false., .false., myrank, nbproc)
#endif

    !-----------------------------------------------------------------------------------------
    ! Calculate data from the prior model.
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) call model(i)%calculate_data(ipar%ndata(i), ipar%ndata_components(i), jinv%matrix, &
        ipar%problem_weight(i), iarr(i)%column_weight, data(i)%val_calc, ipar%compression_type, &
        line_start(i), param_shift(i), myrank, nbproc)
    enddo

#ifndef SUPPRESS_OUTPUT
    ! Write data calculated from the prior model.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_prior_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_prior_', 2, myrank)
#endif

    ! SETTING STARTING MODEL FOR INVERSION -----------------------------------------------------
    if (SOLVE_PROBLEM(1)) call set_model(model(1), gpar%start_model_type, gpar%start_model_val, gpar%model_files(3), myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call set_model(model(2), mpar%start_model_type, mpar%start_model_val, mpar%model_files(3), myrank, nbproc)

#ifndef SUPPRESS_OUTPUT
    ! Write the starting model to a file for visualization.
    if (SOLVE_PROBLEM(1)) call model_write(model(1), 'grav_starting_', .false., .false., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call model_write(model(2), 'mag_starting_', .false., .false., myrank, nbproc)
#endif

    !-----------------------------------------------------------------------------------------
    ! Calculate data from the starting model.
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) call model(i)%calculate_data(ipar%ndata(i), ipar%ndata_components(i), jinv%matrix, &
        ipar%problem_weight(i), iarr(i)%column_weight, data(i)%val_calc, ipar%compression_type, &
        line_start(i), param_shift(i), myrank, nbproc)
    enddo

#ifndef SUPPRESS_OUTPUT
    ! Write data calculated from the starting model.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_starting_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_starting_', 2, myrank)
#endif

    !-----------------------------------------------------------------------------------------
    ! Calculate costs for the models (damping term in the cost function).
    call calculate_model_costs(ipar, iarr, model, cost_model, SOLVE_PROBLEM, myrank, nbproc)

    ! Calculate initial cost (misfit).
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) then
        call calculate_cost(size(data(i)%val_meas), data(i)%val_meas, data(i)%val_calc, cost_data(i), .false., nbproc)
        if (myrank == 0) print *, 'data cost =', cost_data(i)
      endif
    enddo

#ifndef SUPPRESS_OUTPUT
    ! Stores costs.
    if (myrank == 0) &
      open(FILE_COSTS, file=trim(path_output)//'/costs.txt', access='stream', form='formatted', status='replace', action='write')
#endif

    memory = get_max_mem_usage()
    if (myrank == 0) print *, "MEMORY USED (major loop start) [GB] =", memory

    ! Major inversion loop.
    do it = 1, ipar%ninversions

      ! Exiting the loop when a stop file is found.
      if (exit_loop(myrank)) then
        if (myrank == 0) print *, 'Stop file found! Exiting the loop.'
        exit
      endif

      if (myrank == 0) then
        print *, '======================================================='
        print *, 'Iteration, prior model =', it, ' ', m
        print *, '======================================================='
      endif

      ! Calculate data residuals.
      if (SOLVE_PROBLEM(1)) call calculate_residual(size(data(1)%val_meas), data(1)%val_meas, data(1)%val_calc, iarr(1)%residuals)
      if (SOLVE_PROBLEM(2)) call calculate_residual(size(data(2)%val_meas), data(2)%val_meas, data(2)%val_calc, iarr(2)%residuals)

      ! Resets the joint inversion.
      if (it > 1) call jinv%reset(myrank)

      ! Solve joint inverse problem.
      call jinv%solve(ipar, iarr, model, delta_model, delta_data, myrank, nbproc)

      ! Update the local models.
      if (SOLVE_PROBLEM(1)) call model(1)%update(delta_model(:, :, 1))
      if (SOLVE_PROBLEM(2)) call model(2)%update(delta_model(:, :, 2))

      ! Write intermediate models to file.
      if (ipar%write_model_niter > 0) then
        if (mod(it, ipar%write_model_niter) == 0) then
          if (SOLVE_PROBLEM(1)) call model_write(model(1), 'grav_inter_'//trim(str(it))//'_', .true., .false., myrank, nbproc)
          if (SOLVE_PROBLEM(2)) call model_write(model(2), 'mag_inter_'//trim(str(it))//'_', .true., .false., myrank, nbproc)
        endif
      endif

      ! Calculate new data. Using the data update as the grav/mag problems are linear.
      if (SOLVE_PROBLEM(1)) data(1)%val_calc = data(1)%val_calc + delta_data(1)%val
      if (SOLVE_PROBLEM(2)) data(2)%val_calc = data(2)%val_calc + delta_data(2)%val

#ifndef SUPPRESS_OUTPUT
      ! Write costs (for the previous iteration).
      if (myrank == 0) then
        damping_gradient_cost = jinv%get_damping_gradient_cost()
        write(FILE_COSTS, *) it - 1, cost_data(1), cost_data(2), cost_model(1), cost_model(2), &
                             jinv%get_admm_cost(), &
                             damping_gradient_cost, &
                             jinv%get_cross_grad_cost(), &
                             jinv%get_clustering_cost(1), jinv%get_clustering_cost(2)
        flush(FILE_COSTS)
      endif
#endif

      ! Calculate new costs for the models (damping term in the cost function).
      call calculate_model_costs(ipar, iarr, model, cost_model, SOLVE_PROBLEM, myrank, nbproc)

      ! Calculate new costs for data misfits.
      do i = 1, 2
        if (SOLVE_PROBLEM(i)) then
          call calculate_cost(size(data(i)%val_meas), data(i)%val_meas, data(i)%val_calc, cost_data(i), .false., nbproc)
          if (myrank == 0) print *, 'data cost =', cost_data(i)
        endif
      enddo

      ! Adjust the ADMM weight dynamically.
      if (ipar%admm_type > 0) then
        call adjust_admm_weight(ipar, SOLVE_PROBLEM, cost_data, myrank)
      endif

    enddo ! Major inversion loop.

#ifndef SUPPRESS_OUTPUT
    ! Write final costs.
    if (myrank == 0) write(FILE_COSTS, *) ipar%ninversions, cost_data(1), cost_data(2), cost_model(1), cost_model(2)
    if (myrank == 0) close(FILE_COSTS)
#endif

#ifndef SUPPRESS_OUTPUT
    ! Write the final model to a file.
    if (SOLVE_PROBLEM(1)) call model_write(model(1), 'grav_final_', .true., .true., myrank, nbproc)
    if (SOLVE_PROBLEM(2)) call model_write(model(2), 'mag_final_', .true., .true., myrank, nbproc)
#endif

#ifndef SUPPRESS_OUTPUT
    if (myrank == 0) then
      ! Print model value bounds.
      do i = 1, 2
        if (SOLVE_PROBLEM(i)) &
          print *, 'Model', i , 'min/max values =', minval(model(i)%val_full), maxval(model(i)%val_full)
      enddo
    endif
#endif

#ifndef SUPPRESS_OUTPUT
    ! Write data calculated from final model.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_calc_final_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_calc_final_', 2, myrank)

    ! Calculate final data residual.
    do i = 1, 2
      if (SOLVE_PROBLEM(i)) then
        data(i)%val_calc = data(i)%val_meas - data(i)%val_calc
      endif
    enddo

    ! Write final data residual.
    if (SOLVE_PROBLEM(1)) call data(1)%write('grav_misfit_final_', 2, myrank)
    if (SOLVE_PROBLEM(2)) call data(2)%write('mag_misfit_final_', 2, myrank)
#endif

#ifndef SUPPRESS_OUTPUT
    if (jinv%add_cross_grad) then
      ! Write final cross-gradient vector magnitude to a file.
      call jinv%get_cross_grad(model(1)%val_full)
      call model_write(model(1), 'cross_grad_final_', .false., .false., myrank, nbproc)
    endif
#endif

#ifndef SUPPRESS_OUTPUT
    if (jinv%add_clustering) then
      ! Write final clustering probabilities, i.e., P(m) per cell.
      call jinv%get_clustering(model(1)%val_full)
      call model_write(model(1), 'clustering_final_', .false., .false., myrank, nbproc)

      call jinv%clustering%write_data('clustering_data.txt', model(1)%grid_full, myrank)
    endif
#endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  enddo ! loop over prior models
  !******************************

  deallocate(delta_model)
  deallocate(delta_data(1)%val)
  deallocate(delta_data(2)%val)

end subroutine solve_problem_joint_gravmag

!==============================================================================================
! Adjusts the ADMM weight dynamically.
!==============================================================================================
subroutine adjust_admm_weight(ipar, SOLVE_PROBLEM, cost_data, myrank)
  type(t_parameters_inversion), intent(inout) :: ipar
  logical, intent(in) :: SOLVE_PROBLEM(2)
  real(kind=CUSTOM_REAL), intent(in) :: cost_data(2)
  integer, intent(in) :: myrank

  integer :: i

  do i = 1, 2
    if (SOLVE_PROBLEM(i)) then
      if (cost_data(i) < ipar%data_cost_threshold_ADMM &
          .and. ipar%rho_ADMM(i) < ipar%max_weight_ADMM &
          .and. ipar%weight_multiplier_ADMM /= 1.d0) then
        ! Adjust the weight.
        ipar%rho_ADMM(i) = ipar%weight_multiplier_ADMM * ipar%rho_ADMM(i)
        if (myrank == 0) print *, 'Increased the ADMM weight to:', ipar%rho_ADMM(i)
      endif
    endif
  enddo

end subroutine adjust_admm_weight

!==============================================================================================
! Computes and prints norm Lp of the difference between inverted and prior models.
!==============================================================================================
subroutine calculate_model_costs(ipar, iarr, model, cost_model, solve_problem, myrank, nbproc)
  type(t_parameters_inversion), intent(in) :: ipar
  type(t_inversion_arrays), intent(in) :: iarr(2)
  type(t_model), intent(in) :: model(2)
  logical, intent(in) :: solve_problem(2)
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(out) :: cost_model(2)

  integer :: i

  cost_model = 0.d0

  do i = 1, 2
    if (solve_problem(i)) then
      call calculate_cost_model(ipar%nelements, ipar%norm_power, model(i)%val(:, 1), model(i)%val_prior(:, 1), &
                                iarr(i)%column_weight, cost_model(i), nbproc)

      if (myrank == 0) print *, 'model cost =', cost_model(i)
    endif
  enddo
end subroutine calculate_model_costs

!========================================================================================
! Sets the model values: via constant from Parfile or via reading it from a file.
!========================================================================================
subroutine set_model(model, model_type, model_val, model_file, myrank, nbproc)
  integer, intent(in) :: model_type, myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: model_val
  character(len=256), intent(in) :: model_file
  type(t_model), intent(inout) :: model

  if (model_type == 1) then
    ! Setting homogeneous starting value.
    model%val = model_val
    if (myrank == 0) model%val_full = model_val

  else if (model_type == 2) then
    ! Reading from file.
    call model_read(model, model_file, myrank, nbproc)

  else
    call exit_MPI("Unknown model type in set_model!", myrank, model_type)
  endif
end subroutine set_model

!========================================================================================
! Sets the model bounds.
!========================================================================================
subroutine set_model_bounds(ipar, model, problem_type, myrank, nbproc)
  type(t_parameters_inversion), intent(in) :: ipar
  integer, intent(in) :: problem_type
  integer, intent(in) :: myrank, nbproc
  type(t_model), intent(inout) :: model

  integer :: i

  ! Allocate bound arrays.
  call model%allocate_bound_arrays(ipar%nlithos, myrank)

  if (ipar%admm_bound_type == 1) then
    ! Global bounds - define from Parfile parameters.
    do i = 1, model%nelements
      model%min_local_bound(:, i) = ipar%admm_bounds(problem_type)%val(1 : ipar%nlithos)
      model%max_local_bound(:, i) = ipar%admm_bounds(problem_type)%val(ipar%nlithos + 1 : 2 * ipar%nlithos)
      model%local_bound_constraints_weight(i) = ipar%admm_bounds(problem_type)%val(2 * ipar%nlithos + 1)
    enddo
  else
    ! Local bounds - read from file.
    call read_bound_constraints(model, ipar%bounds_ADMM_file(problem_type), myrank, nbproc)
  endif

end subroutine set_model_bounds
!========================================================================================
! Calcualte data residual.
! Note: utilize conversion 2D to 1D array via subroutine interface.
!========================================================================================
pure subroutine calculate_residual(n, d_obs, d_calc, residual)
  integer, intent(in) :: n
  real(kind=CUSTOM_REAL), intent(in) :: d_obs(n)
  real(kind=CUSTOM_REAL), intent(in) :: d_calc(n)
  real(kind=CUSTOM_REAL), intent(out) :: residual(n)

  residual = d_obs - d_calc

end subroutine calculate_residual

!========================================================================================
! Check if need to exit the inversion loop.
!========================================================================================
function exit_loop(myrank) result(res)
  integer, intent(in) :: myrank
  logical :: res, file_exists
  integer :: exist_loc, exist_glob, ierr

  exist_loc = 0

  if (myrank == 0) then
    inquire(file="stop", exist=file_exists)
    if (file_exists) exist_loc = 1
  endif

  ! Communicate the result of inquire() to other ranks.
  call mpi_allreduce(exist_loc, exist_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (exist_glob > 0) then
    res = .true.
  else
    res = .false.
  endif
end function exit_loop

end module problem_joint_gravmag
