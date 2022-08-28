
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

    ! Model parameters.
    real(kind=CUSTOM_REAL), allocatable :: val(:)

    ! Prior model.
    real(kind=CUSTOM_REAL), allocatable :: val_prior(:)

    ! Data arrays for local ADMM constraints.
    integer :: nlithos
    real(kind=CUSTOM_REAL), allocatable :: min_local_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: max_local_bound(:, :)
    real(kind=CUSTOM_REAL), allocatable :: local_bound_constraints_weight(:)

    ! Full grid.
    type(t_grid) :: grid_full

    ! Total number of model parameters.
    integer :: nelements_total
    ! Local number of model parameters.
    integer :: nelements

  contains
    private

    procedure, public, pass :: initialize => model_initialize
    procedure, public, pass :: allocate_bound_arrays => model_allocate_bound_arrays

    procedure, public, pass :: get_Xmin => model_get_Xmin
    procedure, public, pass :: get_Xmax => model_get_Xmax
    procedure, public, pass :: get_Ymin => model_get_Ymin
    procedure, public, pass :: get_Ymax => model_get_Ymax

    procedure, public, pass :: update => model_update
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

  ierr = 0

  allocate(this%val(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  allocate(this%val_prior(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)

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
  real(kind=CUSTOM_REAL), intent(in) :: delta_model(this%nelements_total)

  this%val = this%val + delta_model

end subroutine model_update

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start and param_shift to calculate the data using part of the big (joint) matrix.
!======================================================================================================
subroutine model_calculate_data(this, ndata_loc, matrix_sensit, problem_weight, column_weight, data, &
                                compression_type, line_start, param_shift, myrank)
  class(t_model), intent(in) :: this
  integer, intent(in) :: ndata_loc, compression_type
  integer, intent(in) :: line_start, param_shift
  integer, intent(in) :: myrank
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: column_weight(this%nelements_total)

  real(kind=CUSTOM_REAL), intent(out) :: data(ndata_loc)

  real(kind=CUSTOM_REAL), allocatable :: model_scaled(:)
  integer :: i, ierr

  allocate(model_scaled(this%nelements_total), source=0._CUSTOM_REAL, stat=ierr)
  if (ierr /= 0) call exit_MPI("Dynamic memory allocation error in model_calculate_data!", myrank, ierr)

  ! Rescale the model to calculate data, as we store the depth-weighted sensitivity kernel.
  do i = 1, this%nelements_total
    model_scaled(i) = this%val(i) / column_weight(i)
  enddo

  if (compression_type > 0) then
  ! Apply wavelet transform to the model to calculate data using compressed sensitivity kernel.
    call Haar3D(model_scaled, this%grid_full%nx, this%grid_full%ny, this%grid_full%nz)
  endif

  ! Calculate data: d = S * m
  call matrix_sensit%part_mult_vector(this%nelements_total, model_scaled, ndata_loc, data, line_start, param_shift, myrank)

  deallocate(model_scaled)

  ! Apply the problem weight.
  if (problem_weight /= 0.d0) then
    data = data / problem_weight
  else
    call exit_MPI("Zero problem weight in model_calculate_data!", myrank, 0)
  endif

end subroutine model_calculate_data

!======================================================================================================
! Calculate the linear data using the sensitivity kernel (S) and model (m) as d = S * m.
! Use line_start and param_shift to calculate the data using part of the big (joint) matrix.
! This version uses unscaled model (in wavelet domain).
!======================================================================================================
subroutine calculate_data_unscaled(nelements, model, matrix_sensit, problem_weight, ndata, data, &
                                   line_start, param_shift, myrank)
  integer, intent(in) :: nelements
  real(kind=CUSTOM_REAL), intent(in) :: model(nelements)
  type(t_sparse_matrix), intent(in) :: matrix_sensit
  real(kind=CUSTOM_REAL), intent(in) :: problem_weight
  integer, intent(in) :: ndata, line_start, param_shift
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL), intent(out) :: data(ndata)

  ! Calculate data: d = S' * m'
  ! Assume that both the kernel and the model are unscaled (in wavelet domain).
  call matrix_sensit%part_mult_vector(nelements, model, ndata, data, line_start, param_shift, myrank)

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
subroutine rescale_model(nelements, model, weight)
  integer, intent(in) :: nelements
  real(kind=CUSTOM_REAL), intent(in) :: weight(nelements)
  real(kind=CUSTOM_REAL), intent(inout) :: model(nelements)

  integer :: i

  do i = 1, nelements
    model(i) = model(i) * weight(i)
  enddo
end subroutine rescale_model

end module model
