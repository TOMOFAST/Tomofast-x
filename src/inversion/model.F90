
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

!================================================================================================
! This class contains the model data and basic operations with the model.
!
! Vitaliy Ogarko, UWA, CET, Australia.
!================================================================================================
module model

  use global_typedefs
  use mpi_tools, only: exit_MPI
  use grid
  use string, only: str
  use parallel_tools
  use string
  use sparse_matrix
  use wavelet_transform
  use wavelet_utils

  implicit none

  private

  type, public :: t_model

    ! Local model parameters.
    real(kind=CUSTOM_REAL), allocatable :: val(:, :)

    ! Prior model (local).
    real(kind=CUSTOM_REAL), allocatable :: val_prior(:, :)

    ! Full model parameters.
    ! (Read initial model here, write final model from here, and use also for cross-gradient.)
    real(kind=CUSTOM_REAL), allocatable :: val_full(:, :)

    ! Data arrays for the ADMM constraints.
    integer :: nlithos
    real(kind=CUSTOM_REAL), allocatable :: min_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: max_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: bound_weight(:)

    ! Local weights for damping gradient constraints, for every direction (X, Y, Z).
    real(kind=CUSTOM_REAL), allocatable :: damping_grad_weight(:, :)

    ! Full grid.
    type(t_grid) :: grid_full

    ! Total number of model parameters.
    integer :: nelements_total
    ! Local number of model parameters (on current CPU).
    integer :: nelements
    ! The number of model components.
    integer :: ncomponents

    ! Units multiplier.
    real(kind=CUSTOM_REAL) :: units_mult
    ! The model output label in the vtk.
    character(len=16) :: vtk_label

  contains
    private

    procedure, public, pass :: initialize => model_initialize
    procedure, public, pass :: allocate_bound_arrays => model_allocate_bound_arrays
    procedure, public, pass :: allocate_damping_gradient_arrays => model_allocate_damping_gradient_arrays

    procedure, public, pass :: distribute => model_distribute

    procedure, public, pass :: update => model_update
    procedure, public, pass :: update_full => model_update_full
    procedure, public, pass :: calculate_data => model_calculate_data

  end type t_model

  public :: rescale_model

contains

!================================================================================================
! Initialization.
!================================================================================================
subroutine model_initialize(this, nelements, ncomponents, alloc_full_on_all_cpus, &
                            units_mult, vtk_label, myrank, nbproc)
  class(t_model), intent(inout) :: this
  logical, intent(in) :: alloc_full_on_all_cpus
  integer, intent(in) :: nelements, ncomponents, myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: units_mult
  character(len=*), intent(in) :: vtk_label

  integer :: ierr

  ! Sanity check.
  if (nelements <= 0) then
    call exit_MPI("Bad nelements value in model_initialize!", myrank, nelements)
  endif

  this%nelements = nelements
  this%nelements_total = get_total_number_elements(nelements, myrank, nbproc)
  this%ncomponents = ncomponents

  this%units_mult = units_mult
  this%vtk_label = vtk_label

  ierr = 0

  if (alloc_full_on_all_cpus) then
  ! Allocate the full model array on all cpus.
    allocate(this%val_full(this%nelements_total, this%ncomponents), source=0._CUSTOM_REAL, stat=ierr)
  else
  ! Allocate the full model array on the master cpu only.
    if (myrank == 0) then
      allocate(this%val_full(this%nelements_total, this%ncomponents), source=0._CUSTOM_REAL, stat=ierr)
    else
      allocate(this%val_full(1, 1), source=0._CUSTOM_REAL, stat=ierr)
    endif
  endif

  allocate(this%val(this%nelements, this%ncomponents), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%val_prior(this%nelements, this%ncomponents), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_initialize!", myrank, ierr)

end subroutine model_initialize

!================================================================================================
! Allocate bound arrays (used for ADMM).
!================================================================================================
subroutine model_allocate_bound_arrays(this, nlithos, myrank)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: nlithos, myrank

  integer :: ierr
  ierr = 0

  this%nlithos = nlithos

  allocate(this%min_bound(this%nlithos, this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%max_bound(this%nlithos, this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%bound_weight(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_allocate_bound_arrays!", myrank, ierr)

end subroutine model_allocate_bound_arrays

!================================================================================================
! Allocate damping gradient arrays.
!================================================================================================
subroutine model_allocate_damping_gradient_arrays(this, myrank)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: myrank

  integer :: ierr
  ierr = 0

  ! Note: allocate the full array as some local weights will be needed on other procs due to gradient involves several elements.
  allocate(this%damping_grad_weight(this%nelements_total, 3), source=1._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_allocate_damping_gradient_arrays!", myrank, ierr)

end subroutine model_allocate_damping_gradient_arrays

!=================================================================================
! Distribute the full model among CPUs.
!=================================================================================
subroutine model_distribute(this, myrank, nbproc)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: myrank, nbproc
  integer :: k

  do k = 1, this%ncomponents
    ! Scatter the local array parts.
    call scatter_full_array(this%nelements, this%val_full(:, k), this%val(:, k), myrank, nbproc)
  enddo
end subroutine model_distribute

!======================================================================================================
! Update the local model after inversion.
!======================================================================================================
subroutine model_update(this, delta_model)
  class(t_model), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: delta_model(this%nelements, this%ncomponents)

  this%val = this%val + delta_model

end subroutine model_update

!======================================================================================================
! Update the full model from its local parts (split between CPUs).
!======================================================================================================
subroutine model_update_full(this, broadcast, myrank, nbproc)
  class(t_model), intent(inout) :: this
  logical, intent(in) :: broadcast
  integer, intent(in) :: myrank, nbproc
  integer :: k

  do k = 1, this%ncomponents
    call get_full_array(this%val(:, k), this%nelements, this%val_full(:, k), broadcast, myrank, nbproc)
  enddo
end subroutine model_update_full

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start, line_end, param_shift to calculate the data using part of the big (joint) matrix.
!======================================================================================================
subroutine model_calculate_data(this, ndata, ndata_components, matrix_sensit, problem_weight, column_weight, data_cov, &
                                data_calc, compression_type, line_start, param_shift, myrank, nbproc)
  class(t_model), intent(in) :: this
  integer, intent(in) :: ndata, ndata_components, compression_type
  integer, intent(in) :: line_start, param_shift
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(this%nelements)
  real(kind=CUSTOM_REAL), intent(in) :: data_cov(ndata)

  real(kind=CUSTOM_REAL), intent(out) :: data_calc(ndata_components, ndata)

  real(kind=CUSTOM_REAL), allocatable :: model_scaled(:, :)
  real(kind=CUSTOM_REAL), allocatable :: model_scaled_full(:)
  integer :: i, k, ierr
  logical :: SOLVE_PROBLEM(1)

  allocate(model_scaled(this%nelements, this%ncomponents), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

  ! Rescale the model to calculate data, as we store the depth-weighted sensitivity kernel.
  do k = 1, this%ncomponents
    do i = 1, this%nelements
      model_scaled(i, k) = this%val(i, k) / column_weight(i)
    enddo
  enddo

  if (compression_type > 0) then
  ! Apply wavelet transform to the model to calculate data using compressed sensitivity kernel.

    if (nbproc > 1) then
    ! Parallel version.

      ! Allocate memory for the full model on the master rank only.
      if (myrank == 0) then
        allocate(model_scaled_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      else
        ! Fortran standard requires that allocatable array is allocated when passing by argument.
        allocate(model_scaled_full(1), source=0._CUSTOM_REAL, stat=ierr)
      endif

      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

      SOLVE_PROBLEM = .true.

      ! Transform the model_scaled to wavelet domain.
      call apply_wavelet_transform(this%nelements, this%grid_full%nx, this%grid_full%ny, this%grid_full%nz, this%ncomponents, &
                                   model_scaled, model_scaled_full, &
                                   .true., compression_type, 1, SOLVE_PROBLEM, myrank, nbproc)

      deallocate(model_scaled_full)

    else
    ! Serial version.
      do k = 1, this%ncomponents
        call forward_wavelet(model_scaled(:, k), this%grid_full%nx, this%grid_full%ny, this%grid_full%nz, compression_type)
      enddo
    endif
  endif

  ! Calculate data: d = S * m
  ! Note, the 2D-array model_scaled is remapped to 1D on the input of part_mult_vector().
  call matrix_sensit%part_mult_vector(size(model_scaled), model_scaled, &
                                      size(data_calc), data_calc, line_start, param_shift, myrank)

  deallocate(model_scaled)

  call MPI_Allreduce(MPI_IN_PLACE, data_calc, size(data_calc), CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in model_calculate_data!", myrank, ierr)

  ! Apply the problem weight.
  if (problem_weight /= 0.d0) then
    data_calc = data_calc / problem_weight
  else
    call exit_MPI("Zero problem weight in model_calculate_data!", myrank, 0)
  endif

  ! Apply data covariance transform (as the sensitivity kernel is scaled by covariance).
  do i = 1, ndata
    data_calc(:, i) = data_calc(:, i) / data_cov(i)
  enddo

end subroutine model_calculate_data

!================================================================================================
! Weights the model parameters.
!================================================================================================
subroutine rescale_model(nelements, ncomponents, model, weight)
  integer, intent(in) :: nelements, ncomponents
  real(kind=CUSTOM_REAL), intent(in) :: weight(nelements)
  real(kind=CUSTOM_REAL), intent(inout) :: model(nelements, ncomponents)

  integer :: i, k

  do k = 1, ncomponents
    do i = 1, nelements
      model(i, k) = model(i, k) * weight(i)
    enddo
  enddo
end subroutine rescale_model

end module model
