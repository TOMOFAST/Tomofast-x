
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
    real(kind=CUSTOM_REAL), allocatable :: val(:)

    ! Prior model (local).
    real(kind=CUSTOM_REAL), allocatable :: val_prior(:)

    ! Full model parameters.
    ! (Read initial model here, write final model from here, and use also for cross-gradient.)
    real(kind=CUSTOM_REAL), allocatable :: val_full(:)

    ! Data arrays for local ADMM constraints.
    integer :: nlithos
    real(kind=CUSTOM_REAL), allocatable :: min_local_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: max_local_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: local_bound_constraints_weight(:)

    ! Full grid.
    type(t_grid) :: grid_full

    ! Total number of model parameters.
    integer :: nelements_total
    ! Local number of model parameters (on current CPU).
    integer :: nelements

    ! Flag to check if the full model has been updated (to the state of the local model).
    logical :: full_model_updated

  contains
    private

    procedure, public, pass :: initialize => model_initialize
    procedure, public, pass :: allocate_bound_arrays => model_allocate_bound_arrays

    procedure, public, pass :: distribute => model_distribute

    procedure, public, pass :: get_Xmin => model_get_Xmin
    procedure, public, pass :: get_Xmax => model_get_Xmax
    procedure, public, pass :: get_Ymin => model_get_Ymin
    procedure, public, pass :: get_Ymax => model_get_Ymax

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
subroutine model_initialize(this, nelements, myrank, nbproc)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: nelements, myrank, nbproc

  integer :: ierr

  ! Sanity check.
  if (nelements <= 0) then
    call exit_MPI("Bad nelements value in model_initialize!", myrank, nelements)
  endif

  this%nelements = nelements
  this%nelements_total = get_total_number_elements(nelements, myrank, nbproc)

  this%full_model_updated = .true.

  ierr = 0

  allocate(this%val_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%val(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%val_prior(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

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

  allocate(this%min_local_bound(this%nelements, this%nlithos), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%max_local_bound(this%nelements, this%nlithos), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%local_bound_constraints_weight(this%nelements), source=0._CUSTOM_REAL, stat=ierr)

  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_allocate_bound_arrays!", myrank, ierr)

end subroutine model_allocate_bound_arrays

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
  integer :: ierr

  ! Partitioning for MPI_Scatterv.
  call get_mpi_partitioning(this%nelements, displs, nelements_at_cpu, myrank, nbproc)

  call MPI_Scatterv(this%val_full, nelements_at_cpu, displs, CUSTOM_MPI_TYPE, &
                    this%val, this%nelements, CUSTOM_MPI_TYPE, 0, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("Error in MPI_Scatterv in model_distribute!", myrank, ierr)

end subroutine model_distribute

!================================================================================================
! Get the minimum X-coordinate of the model grid.
!================================================================================================
pure function model_get_Xmin(this) result(res)
  class(t_model), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = minval(this%grid_full%X1)
end function model_get_Xmin

!================================================================================================
! Get the maximum X-coordinate of the model grid.
!================================================================================================
pure function model_get_Xmax(this) result(res)
  class(t_model), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = maxval(this%grid_full%X2)
end function model_get_Xmax

!================================================================================================
! Get the minimum Y-coordinate of the model grid.
!================================================================================================
pure function model_get_Ymin(this) result(res)
  class(t_model), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = minval(this%grid_full%Y1)
end function model_get_Ymin

!================================================================================================
! Get the maximum Y-coordinate of the model grid.
!================================================================================================
pure function model_get_Ymax(this) result(res)
  class(t_model), intent(in) :: this
  real(kind=CUSTOM_REAL) :: res

  res = maxval(this%grid_full%Y2)
end function model_get_Ymax

!======================================================================================================
! Update the local model after inversion.
!======================================================================================================
subroutine model_update(this, delta_model)
  class(t_model), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: delta_model(:)

  integer :: i

  do i = 1, this%nelements
    this%val(i) = this%val(i) + delta_model(i)
  enddo

  this%full_model_updated = .false.

end subroutine model_update

!======================================================================================================
! Update the full model from its local parts (split between CPUs).
!======================================================================================================
subroutine model_update_full(this, myrank, nbproc)
  class(t_model), intent(inout) :: this
  integer, intent(in) :: myrank, nbproc

  call get_full_array(this%val, this%nelements, this%val_full, .true., myrank, nbproc)

  this%full_model_updated = .true.

end subroutine model_update_full

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start, line_end, param_shift to calculate the data using part of the big (joint) matrix.
!======================================================================================================
subroutine model_calculate_data(this, ndata, matrix_sensit, problem_weight, column_weight, data, compression_type, &
                                line_start, line_end, param_shift, myrank, nbproc)
  class(t_model), intent(in) :: this
  integer, intent(in) :: ndata, compression_type
  integer, intent(in) :: line_start, line_end, param_shift
  integer, intent(in) :: myrank, nbproc
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(:)

  real(kind=CUSTOM_REAL), intent(out) :: data(ndata)

  real(kind=CUSTOM_REAL), allocatable :: model_scaled(:)
  real(kind=CUSTOM_REAL), allocatable :: model_scaled_full(:)
  integer :: i, ierr
  integer :: nsmaller

  allocate(model_scaled(this%nelements), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

  ! Rescale the model to calculate data, as we store the depth-weighted sensitivity kernel.
  do i = 1, this%nelements
    model_scaled(i) = this%val(i) / column_weight(i)
  enddo

  if (compression_type > 0) then
  ! Apply wavelet transform to the model to calculate data using compressed sensitivity kernel.

    if (nbproc > 1) then
    ! Parallel version.

      ! Allocate memory for the full model.
      allocate(model_scaled_full(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
      if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

      ! Gather the full model from all processors.
      call get_full_array(model_scaled, this%nelements, model_scaled_full, .true., myrank, nbproc)

      ! Compress the full model.
      call Haar3D(model_scaled_full, this%grid_full%nx, this%grid_full%ny, this%grid_full%nz)

      ! Extract the local model part.
      nsmaller = get_nsmaller(this%nelements, myrank, nbproc)
      model_scaled = model_scaled_full(nsmaller + 1 : nsmaller + this%nelements)

      deallocate(model_scaled_full)

    else
    ! Serial version.
      call Haar3D(model_scaled, this%grid_full%nx, this%grid_full%ny, this%grid_full%nz)
    endif
  endif

  ! Calculate data: d = S * m
  call matrix_sensit%part_mult_vector(model_scaled, data, line_start, line_end, param_shift, myrank)

  deallocate(model_scaled)

  call MPI_Allreduce(MPI_IN_PLACE, data, ndata, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in model_calculate_data!", myrank, ierr)

  ! Apply the problem weight.
  if (problem_weight /= 0.d0) then
    data = data / problem_weight
  else
    call exit_MPI("Zero problem weight in model_calculate_data!", myrank, 0)
  endif

end subroutine model_calculate_data

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start, line_end, param_shift to calculate the data using part of the big (joint) matrix.
! This version uses unscaled model (in wavelet domain).
!======================================================================================================
subroutine calculate_data_unscaled(model, matrix_sensit, problem_weight, ndata, data, &
                                   line_start, line_end, param_shift, myrank)
  real(kind=CUSTOM_REAL), intent(in) :: model(:)
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  integer, intent(in) :: ndata, line_start, line_end, param_shift
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: data(ndata)

  integer :: ierr

  ! Calculate data: d = S' * m'
  ! Assume that both the kernel and the model are unscaled (in wavelet domain).
  call matrix_sensit%part_mult_vector(model, data, line_start, line_end, param_shift, myrank)

  call MPI_Allreduce(MPI_IN_PLACE, data, ndata, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (ierr /= 0) call exit_MPI("MPI error in calculate_data_unscaled!", myrank, ierr)

  ! Apply the problem weight.
  if (problem_weight /= 0.d0) then
    data = data / problem_weight
  else
    call exit_MPI("Zero problem weight in calculate_data_unscaled!", myrank, 0)
  endif

end subroutine calculate_data_unscaled

!================================================================================================
! Weights the model parameters.
!================================================================================================
subroutine rescale_model(model, weight, nc)
  real(kind=CUSTOM_REAL), intent(inout) :: model(:)
  integer, intent(in) :: nc
  real(kind=CUSTOM_REAL), intent(in) :: weight(:)

  integer :: i

  do i = 1, nc
    model(i) = model(i) * weight(i)
  enddo
end subroutine rescale_model

end module model
