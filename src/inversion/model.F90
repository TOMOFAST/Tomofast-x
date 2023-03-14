
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
    real(kind=CUSTOM_REAL), allocatable :: min_local_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: max_local_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: local_bound_constraints_weight(:)

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

    ! Flag to check if the full model has been updated (to the state of the local model).
    logical :: full_model_updated

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
  public  :: calculate_data_unscaled

contains

!================================================================================================
! Initialization.
!================================================================================================
subroutine model_initialize(this, nelements, ncomponents, myrank, nbproc)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: nelements, ncomponents, myrank, nbproc

  integer :: ierr

  ! Sanity check.
  if (nelements <= 0) then
    call exit_MPI("Bad nelements value in model_initialize!", myrank, nelements)
  endif

  this%nelements = nelements
  this%nelements_total = get_total_number_elements(nelements, myrank, nbproc)

  this%ncomponents = ncomponents

  this%full_model_updated = .true.

  ierr = 0

  allocate(this%val_full(this%nelements_total, this%ncomponents), source=0._CUSTOM_REAL, stat=ierr)
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

  allocate(this%min_local_bound(this%nlithos, this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%max_local_bound(this%nlithos, this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%local_bound_constraints_weight(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

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
! Distribute the grid and the prior model among CPUs.
!=================================================================================
subroutine model_distribute(this, myrank, nbproc)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: myrank, nbproc

  ! Displacement for mpi_scatterv.
  integer :: displs(nbproc)
  ! The number of elements on every CPU for mpi_scatterv.
  integer :: nelements_at_cpu(nbproc)
  integer :: ierr, k

  ! Partitioning for MPI_Scatterv.
  call get_mpi_partitioning(this%nelements, displs, nelements_at_cpu, myrank, nbproc)

  do k = 1, this%ncomponents
    call MPI_Scatterv(this%val_full(:, k), nelements_at_cpu, displs, CUSTOM_MPI_TYPE, &
                      this%val(:, k), this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)
  enddo

  if (ierr /= 0) call exit_MPI("Error in MPI_Scatterv in model_distribute!", myrank, ierr)

end subroutine model_distribute

!======================================================================================================
! Update the local model after inversion.
!======================================================================================================
subroutine model_update(this, delta_model)
  class(t_model), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: delta_model(this%nelements, this%ncomponents)

  this%val = this%val + delta_model

  this%full_model_updated = .false.

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

  if (broadcast) then
    this%full_model_updated = .true.
  endif

end subroutine model_update_full

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start, line_end, param_shift to calculate the data using part of the big (joint) matrix.
!======================================================================================================
subroutine model_calculate_data(this, ndata, ndata_components, matrix_sensit, problem_weight, column_weight, data_calc, &
                                compression_type, line_start, param_shift, myrank, nbproc)
  class(t_model), intent(in) :: this
  integer, intent(in) :: ndata, ndata_components, compression_type
  integer, intent(in) :: line_start, param_shift
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(this%nelements)

  real(kind=CUSTOM_REAL), intent(out) :: data_calc(ndata_components, ndata)

  real(kind=CUSTOM_REAL), allocatable :: model_scaled(:, :)
  real(kind=CUSTOM_REAL), allocatable :: model_scaled_full(:)
  integer :: i, k, ierr
  integer :: nsmaller

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

      ! Allocate memory for the full model.
      allocate(model_scaled_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

      nsmaller = get_nsmaller(this%nelements, myrank, nbproc)

      do k = 1, this%ncomponents
        ! Gather the full model from all processors.
        call get_full_array(model_scaled(:, k), this%nelements, model_scaled_full, .true., myrank, nbproc)

        ! Compress the full model.
        call Haar3D(model_scaled_full, this%grid_full%nx, this%grid_full%ny, this%grid_full%nz)

        ! Extract the local model part.
        model_scaled(:, k) = model_scaled_full(nsmaller + 1 : nsmaller + this%nelements)
      enddo
      deallocate(model_scaled_full)

    else
    ! Serial version.
      do k = 1, this%ncomponents
        call Haar3D(model_scaled(:, k), this%grid_full%nx, this%grid_full%ny, this%grid_full%nz)
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

end subroutine model_calculate_data

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start, line_end, param_shift to calculate the data using part of the big (joint) matrix.
! This version uses unscaled model (in wavelet domain).
!======================================================================================================
subroutine calculate_data_unscaled(nelements, nmodel_components, model, matrix_sensit, problem_weight, &
                                   ndata, ndata_components, data_calc, line_start, param_shift, myrank)
  integer, intent(in) :: nelements, nmodel_components
  real(kind=CUSTOM_REAL), intent(in) :: model(nelements, nmodel_components)
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  integer, intent(in) :: ndata, ndata_components, line_start, param_shift
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: data_calc(ndata_components, ndata)

  integer :: ierr

  ! Calculate data: d = S' * m'
  ! Assume that both the kernel and the model are unscaled (depth weighted and in wavelet domain if compression is activated).
  call matrix_sensit%part_mult_vector(size(model), model, &
                                      size(data_calc), data_calc, line_start, param_shift, myrank)

  call MPI_Allreduce(MPI_IN_PLACE, data_calc, size(data_calc), CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in calculate_data_unscaled!", myrank, ierr)

  ! Apply the problem weight.
  if (problem_weight /= 0.d0) then
    data_calc = data_calc / problem_weight
  else
    call exit_MPI("Zero problem weight in calculate_data_unscaled!", myrank, 0)
  endif

end subroutine calculate_data_unscaled

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
