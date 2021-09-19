
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

!===============================================================================================
! Electrical capacitance tomography (ECT) problem.
!
! Vitaliy Ogarko, UWA, CET, Australia, 2015-2016.
!===============================================================================================
module problem_ect

  use global_typedefs, only: CUSTOM_REAL
  use mpi_tools, only: exit_MPI
  use parameters_ect
  use parameters_inversion
  use inversion_arrays
  use forward_problem_ect
  use data_ect
  use inverse_problem
  use costs

  implicit none

  private

  public :: solve_problem_ect

contains

!=============================================================
! Solve electrical capacitance tomography (ECT) problem.
!=============================================================
subroutine solve_problem_ect(epar, ipar, myrank, nbproc)

  type(t_parameters_ect), intent(inout) :: epar
  type(t_parameters_inversion), intent(inout) :: ipar
  integer, intent(in) :: myrank, nbproc

  ! Local variables.

  ! Inversion allocatable arrays.
  type(t_inversion_arrays) :: iarr
  ! Measured data (or calculated using original model).
  real(kind=CUSTOM_REAL), allocatable :: data_measured(:)
  ! Predicted data using model from inversion.
  real(kind=CUSTOM_REAL), allocatable :: data_calculated(:)
  ! Data calculated using model low permittivity material (empty sensor), used for normalization.
  real(kind=CUSTOM_REAL), allocatable :: data_low(:)
  ! Data calculated using model high permittivity material (full sensor), used for normalization.
  real(kind=CUSTOM_REAL), allocatable :: data_high(:)
  ! Real model.
  real(kind=CUSTOM_REAL), allocatable :: model_real(:)

  type(t_inversion) :: inversion

  integer :: it, ierr
  ! Costs (misfit).
  real(kind=CUSTOM_REAL) :: cost_data
  real(kind=CUSTOM_REAL) :: cost_model

  ! Memory allocation.
  call iarr%initialize(ipar%nelements, ipar%ndata(1), ipar%nx, ipar%ny, ipar%nz)
  call iarr%allocate_aux(myrank)
  call iarr%allocate_sensit(.true., 0, myrank)

  ! Allocate memory for model (with grid) objects.
  call iarr%init_model(myrank, nbproc)

  allocate(data_measured(ipar%ndata(1)), source=0._CUSTOM_REAL, stat=ierr)
  allocate(data_calculated(ipar%ndata(1)), source=0._CUSTOM_REAL, stat=ierr)
  allocate(data_low(ipar%ndata(1)), source=0._CUSTOM_REAL, stat=ierr)
  allocate(data_high(ipar%ndata(1)), source=0._CUSTOM_REAL, stat=ierr)
  allocate(model_real(ipar%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error, you probably ran out of memory; exiting...", myrank, ierr)

  ! Stores misfit.
  if (myrank == 0) open(10,file=trim(path_output)//'misfit',status='unknown',action='write')

  ! Initially, the synthetic model with bubbles is used.
  epar%read_model_from_inversion = 0

  ! Use real model (with bubbles), and calculate the capacitance data used for inversion.
  call solve_forward_problem_ect(epar, iarr%sensitivity, data_measured, model_real, &
                                 iarr%column_weight, iarr%damping_weight, myrank, nbproc)

  ! Remove bubbles.
  epar%num_bubbles = 0

  ! Calculated data using an empty tube (for normalization).
  epar%permit_matrix = epar%permit_air
  call solve_forward_problem_ect(epar, iarr%sensitivity, data_low, iarr%model%val, &
                                 iarr%column_weight, iarr%damping_weight, myrank, nbproc)

  ! Calculated data using a full tube (for normalization).
  epar%permit_matrix = epar%permit_isolated_tube
  call solve_forward_problem_ect(epar, iarr%sensitivity, data_high, iarr%model%val, &
                                 iarr%column_weight, iarr%damping_weight, myrank, nbproc)

  if (myrank == 0) then
    ! Check that data_high > data_low.
    call check_data(epar, ipar%ndata(1), data_low, data_high, myrank)
  endif

  ! Use oil for the prior model.
  epar%permit_matrix = epar%permit_oil

  ! Solve forward problem, and also generate the initial prior model and weights needed for inversion.
  call solve_forward_problem_ect(epar, iarr%sensitivity, data_calculated, iarr%model_prior, &
                                 iarr%column_weight, iarr%damping_weight, myrank, nbproc)

  if (myrank == 0) then
    ! Calculate initial misfit.
    call calculate_cost(ipar%ndata(1), data_measured, data_calculated, cost_data, myrank)

    ! Write initial misfit to a file.
    write(10,*) 0, cost_data
  endif

  ! The model from inversion is used in next forward problems.
  epar%read_model_from_inversion = 1

  ! Initial model for inversion model[0] = prior model (sets the damping).
  iarr%model%val = iarr%model_prior

  ! Initialize inversion.
  call inversion%initialize(ipar, myrank)

  ! Nonlinear inversion loop.
  do it = 1, ipar%ninversions
    if (myrank == 0) print *, 'it = ', it

    ! Calculate data residuals.
    iarr%residuals = data_measured - data_calculated

    ! Apply normalization typical in ECT: data_new = (data - data_low) / (data_high - data_low),
    call ect_normalization(ipar%ndata(1), iarr%sensitivity, iarr%residuals, data_high, data_low)
    ! Solve inverse problem.
    call inversion%solve(ipar, iarr, myrank, nbproc)

    ! Update the model.
    iarr%model%val = iarr%model%val + inversion%get_model_change()

    ! Compute norm Lp of the difference between inverted and prior model.
    call calculate_cost_model(ipar%nelements, ipar%norm_power, iarr%model%val, iarr%model_prior, &
                              iarr%damping_weight, cost_model, nbproc)

    ! Reset the inversion object.
    call inversion%reset()

    ! Solve forward problem.
    call solve_forward_problem_ect(epar, iarr%sensitivity, data_calculated, iarr%model%val, &
                                   iarr%column_weight, iarr%damping_weight, myrank, nbproc)

    if (myrank == 0) then
      ! Calculate cost (misfit).
      call calculate_cost(ipar%ndata(1), data_measured, data_calculated, cost_data, myrank)

      ! Write misfit to a file.
      write(10, *) it, cost_data, cost_model
    endif
  enddo

  if (myrank == 0) close(10)
  !----------------------------------------------------------------

  if (allocated(model_real)) deallocate(model_real)
  if (allocated(data_high)) deallocate(data_high)
  if (allocated(data_low)) deallocate(data_low)
  if (allocated(data_calculated)) deallocate(data_calculated)
  if (allocated(data_measured)) deallocate(data_measured)

end subroutine solve_problem_ect

end module problem_ect
